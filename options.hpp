#ifndef DYNEARTHSOL3D_OPTIONS_HPP
#define DYNEARTHSOL3D_OPTIONS_HPP

/* A thin registry over CLI11 for the config-file parameters.
 *
 * Every parameter is named "section.key". Each section becomes a configurable
 * CLI11 subcommand, which is how CLI11 maps the [section] headers of an INI
 * file onto options. The declaration syntax is
 *
 *     reg.add_options()
 *         ("section.key", opt::value<T>(&p.section.key)->default_value(d), "help")
 *         ("section.list", opt::list()->default_value("[1, 2]"), "help")
 *         ;
 *
 * where opt::value<T> binds a scalar (or string) to a variable and opt::list
 * holds a bracketed list that the caller reads back with str() and parses.
 *
 * The file is read in two passes. scan_lines() is the only authority on line
 * structure and applies the grammar the Boost parser had: "[section]",
 * "key = value", blank lines, and "#" comments; anything else is an invalid
 * line, so a value cannot continue on the next line. It also rejects a repeated
 * key and an unbracketed list, and keeps each raw value. CLI11 then converts
 * the values. It would read "010" as octal, "0x10" as a number, and split or
 * unquote a string, so numbers and bools are validated by a transform and a
 * string takes the raw text from the first pass.
 */

#include <any>
#include <cctype>
#include <fstream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <sys/stat.h>

#include <CLI11.hpp>

namespace opt {

// Thrown for any failure to read or interpret the config file. `bad_value` is
// true when the file is readable but its content is wrong (unknown key,
// duplicate key, unparsable value), false for a missing file or key.
struct ConfigError : public std::runtime_error {
    bool bad_value;
    ConfigError(const std::string &msg, bool bad) : std::runtime_error(msg), bad_value(bad) {}
};


class Registry;

namespace detail {

struct Entry {
    CLI::Option *opt = nullptr;
    bool has_default = false;
    bool required = false;
    std::any ref;                           // T*, for typed access
    std::string (*to_str)(const Entry &) = nullptr;
    const std::vector<std::string> *tokens = nullptr;  // list options only
    std::string *text = nullptr;                       // string options only
};

class SpecBase {
public:
    virtual ~SpecBase() {}
    virtual void attach(CLI::App *section, const std::string &key,
                        const std::string &help, Entry &e) = 0;
};

// Validators that give the strictness lexical_cast had. CLI11 alone reads
// "010" as octal 8 and accepts "0x10" for integers and floats, and any of
// y/n/t/f for a bool.
inline std::string normalize_integer(std::string &s)
{
    static const std::regex re("[+-]?[0-9]+");
    if (! std::regex_match(s, re)) return "not a decimal integer";
    const std::size_t sign = (s[0] == '+' || s[0] == '-') ? 1 : 0;
    std::size_t nz = sign;
    while (nz + 1 < s.size() && s[nz] == '0') ++nz;   // keep a final "0"
    s.erase(sign, nz - sign);
    return "";
}

inline std::string check_float(std::string &s)
{
    static const std::regex re(
        "[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?"
        "|[+-]?([iI][nN][fF]([iI][nN][iI][tT][yY])?|[nN][aA][nN])");
    return std::regex_match(s, re) ? "" : "not a decimal number";
}

inline std::string normalize_bool(std::string &s)
{
    std::string t;
    for (char c : s) t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (t == "1" || t == "true" || t == "yes" || t == "on")  { s = "1"; return ""; }
    if (t == "0" || t == "false" || t == "no" || t == "off") { s = "0"; return ""; }
    return "not one of true/false, yes/no, on/off, 1/0";
}

inline std::string trim(const std::string &s)
{
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string::npos) return "";
    return s.substr(b, s.find_last_not_of(" \t\r") - b + 1);
}

// True if `v` is a list as the old parser read it: "[" elements "]", the elements
// separated by commas, one trailing comma allowed, "[]" empty. An element is plain
// text with no brackets, quotes or blanks, so the reader that converts it and the
// text CLI11 would hand over agree. Anything else, such as "[[1]]", "[\"1\"]",
// "[1,,2]" or "[1 2]", was an error before and stays one.
inline bool is_plain_list(const std::string &v)
{
    if (v.size() < 2 || v.front() != '[' || v.back() != ']') return false;
    const std::string inner = trim(v.substr(1, v.size() - 2));
    if (inner.empty()) return true;

    std::size_t begin = 0;
    for (;;) {
        const std::size_t comma = inner.find(',', begin);
        const bool last = comma == std::string::npos;
        const std::string item = trim(inner.substr(begin, last ? std::string::npos : comma - begin));
        if (item.empty()) {
            // Only a trailing comma leaves the last item empty, and then not alone.
            if (! (last && begin > 0)) return false;
        }
        else if (item.find_first_not_of("0123456789+-.abcdefghijklmnopqrstuvwxyz"
                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ") != std::string::npos)
            return false;
        if (last) return true;
        begin = comma + 1;
    }
}

// A name is letters, digits and _, and does not start with a digit.
inline bool is_identifier(const std::string &s)
{
    if (s.empty() || std::isdigit(static_cast<unsigned char>(s[0]))) return false;
    for (char c : s)
        if (! (std::isalnum(static_cast<unsigned char>(c)) || c == '_')) return false;
    return true;
}

// "section.name" with both parts an identifier, and exactly one dot.
inline bool is_option_name(const std::string &key)
{
    const std::size_t dot = key.find('.');
    return dot != std::string::npos
        && is_identifier(key.substr(0, dot)) && is_identifier(key.substr(dot + 1));
}

inline std::string join_tokens(const std::vector<std::string> &v)
{
    std::string s = "[";
    for (std::size_t i=0; i<v.size(); ++i) {
        if (i) s += ", ";
        s += v[i];
    }
    return s + "]";
}

// Split "[a, b, c]" (brackets optional) into tokens.
inline std::vector<std::string> split_tokens(const std::string &str)
{
    std::string s;
    for (char c : str) if (c != '[' && c != ']') s += c;
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        const auto b = tok.find_first_not_of(" \t");
        if (b == std::string::npos) continue;
        out.push_back(tok.substr(b, tok.find_last_not_of(" \t") - b + 1));
    }
    return out;
}

}  // namespace detail


// A scalar or string parameter, bound to a variable or (when ref is null) owned.
template<class T>
class Value : public detail::SpecBase {
public:
    explicit Value(T *ref) : ref_(ref ? ref : &own_) {}
    Value* default_value(const T &v) { has_def_ = true; def_ = v; return this; }
    Value* required() { req_ = true; return this; }

    void attach(CLI::App *section, const std::string &key,
                const std::string &help, detail::Entry &e) override
    {
        if (has_def_) *ref_ = def_;
        CLI::Option *o = section->add_option("--" + key, *ref_, help);
        if constexpr (std::is_same<T, std::string>::value) {
            // CLI11 splits a value on whitespace and commas and strips quotes.
            // Take it from the raw text instead (see Registry::read_file); this
            // only keeps CLI11 from rejecting the line.
            o->multi_option_policy(CLI::MultiOptionPolicy::Join)->delimiter(' ');
            e.text = ref_;
        }
        else if constexpr (std::is_same<T, bool>::value)
            o->transform(CLI::Validator(detail::normalize_bool, "BOOLEAN"));
        else if constexpr (std::is_integral<T>::value)
            o->transform(CLI::Validator(detail::normalize_integer, "INT"));
        else if constexpr (std::is_floating_point<T>::value)
            o->transform(CLI::Validator(detail::check_float, "FLOAT"));
        if (has_def_) o->capture_default_str();
        e.opt = o;
        e.has_default = has_def_;
        e.required = req_;
        e.ref = ref_;
        e.to_str = [](const detail::Entry &x) {
            std::ostringstream os;
            os << *std::any_cast<T*>(x.ref);
            return os.str();
        };
    }

private:
    T own_{};
    T *ref_;
    bool has_def_ = false, req_ = false;
    T def_{};
};


// A bracketed list of numbers, e.g. "[1e3, 2e3]". Elements are kept as text so
// the caller can parse them as int or double.
class List : public detail::SpecBase {
public:
    List* default_value(const std::string &s) {
        has_def_ = true; def_ = s; tokens_ = detail::split_tokens(s); return this;
    }
    List* required() { req_ = true; return this; }

    void attach(CLI::App *section, const std::string &key,
                const std::string &help, detail::Entry &e) override
    {
        CLI::Option *o = section->add_option("--" + key, tokens_, help);
        o->expected(0, -1);  // allow "[]"
        if (has_def_) o->default_str(def_);
        e.opt = o;
        e.has_default = has_def_;
        e.required = req_;
        e.tokens = &tokens_;
        e.to_str = [](const detail::Entry &x) { return detail::join_tokens(*x.tokens); };
    }

private:
    std::vector<std::string> tokens_;
    std::string def_;
    bool has_def_ = false, req_ = false;
};


template<class T>
inline Value<T>* value(T *ref = nullptr) { return new Value<T>(ref); }

inline List* list() { return new List(); }


class Registry {
public:
    Registry() {
        app_.allow_config_extras(CLI::config_extras_mode::error);
        app_.set_help_flag();   // -h/--help is handled by the caller
    }

    class Adder {
    public:
        explicit Adder(Registry &r) : reg_(r) {}
        Adder& operator()(const char *name, detail::SpecBase *spec, const char *help)
        {
            reg_.add(name, std::unique_ptr<detail::SpecBase>(spec), help);
            return *this;
        }
    private:
        Registry &reg_;
    };

    Adder add_options() { return Adder(*this); }

    // Read the INI file and fill every bound variable. Throws ConfigError.
    void read_file(const char *filename)
    {
        // stat(), not std::filesystem: g++ 8 needs -lstdc++fs for the latter.
        struct stat st;
        std::ifstream in(filename);
        if (! in.good() || stat(filename, &st) != 0 || ! S_ISREG(st.st_mode))
            throw ConfigError(std::string("cannot open file '") + filename + "'", false);
        const std::map<std::string, std::string> raw = scan_lines(in);

        app_.set_config("--des-config-file", filename, "", true);
        try {
            app_.parse(std::vector<std::string>{});
        }
        catch (const CLI::ConfigError &e) {
            throw ConfigError(std::string(e.what()) + " (a misspelled or retired key)", true);
        }
        catch (const CLI::ParseError &e) {
            throw ConfigError(e.what(), true);
        }

        for (const auto &kv : entries_) {
            if (kv.second.required && kv.second.opt->count() == 0)
                throw ConfigError("the option '" + kv.first + "' is required but missing", false);
            // A string is the raw text: the whole value, quotes and commas kept.
            if (kv.second.text && kv.second.opt->count() > 0) {
                const auto it = raw.find(kv.first);
                if (it == raw.end())   // CLI11 set an option the line scan never saw
                    throw ConfigError("option '" + kv.first + "' was read from the "
                                      "file, but not by name; check how it is written", true);
                *kv.second.text = it->second;
            }
        }
    }

    // First pass over the file; see the comment at the top. Returns the raw
    // value of every key, by "section.key".
    std::map<std::string, std::string> scan_lines(std::istream &in) const
    {
        std::map<std::string, std::string> raw;
        std::string text, section;
        bool in_section = false;
        // An invalid line stops the scan at once. A repeated key or an unbracketed
        // list is held back until the whole file has been read, so that a
        // malformed line further down is reported first, as the Boost parser did.
        std::unique_ptr<ConfigError> deferred;
        const auto defer = [&](const std::string &msg) {
            if (! deferred) deferred.reset(new ConfigError(msg, true));
        };
        // getline writes only to text, never to a string this code assigned. With
        // g++ 8 headers in C++17, nvc++ -gpu=mem:managed allocates a string made
        // here with managed new, but the out-of-line getline in libstdc++.so frees
        // with glibc: growing such a string aborts (munmap_chunk(): invalid pointer).
        while (std::getline(in, text)) {
            const std::string shown = detail::trim(text);
            const std::string line = detail::trim(text.substr(0, text.find('#')));
            if (line.empty()) continue;

            const auto bad_line = [&]() {
                return ConfigError("the options configuration file contains an invalid line '"
                                   + shown + "'", false);
            };
            if (line[0] == '[') {
                if (line.back() != ']') throw bad_line();
                section = detail::trim(line.substr(1, line.size() - 2));
                in_section = true;
                continue;
            }
            const auto eq = line.find('=');
            if (eq == std::string::npos) throw bad_line();
            const std::string name = detail::trim(line.substr(0, eq));

            // A key at the top of the file may carry its own section, "mat.rho0".
            const std::string key = in_section ? section + "." + name : name;
            const std::string value = detail::trim(line.substr(eq + 1));

            // From here on a key is "section.name", each part letters, digits and _
            // and not starting with a digit, and it names a registered option. CLI11
            // reads names more freely: it strips quotes and backticks, folds "[[a]]"
            // into the section a, reads "[default]" as no section, and takes a line
            // that starts with ";" as a comment. Refusing every other name here keeps
            // the name CLI11 finds the one this scan saw. Boost found no such option
            // for any of them either.
            if (! detail::is_option_name(key)) {
                defer("unrecognised option '" + key + "' (a misspelled or retired key; a "
                      "name may only use letters, digits and _, and cannot start with a digit)");
                continue;
            }
            const detail::Entry *e = find(key);
            if (! e) {
                defer("unrecognised option '" + key + "' (a misspelled or retired key)");
                continue;
            }
            if (! raw.emplace(key, value).second)
                defer("option '" + key + "' cannot be specified more than once");

            // CLI11 reads on past the end of the line for a value that opens with
            // three quotes (a multi-line string) or with an unclosed "[" (an array),
            // and the keys after it vanish without a word. DES takes neither.
            if (value.compare(0, 3, "\"\"\"") == 0 || value.compare(0, 3, "'''") == 0)
                defer("option '" + key + "' cannot start with three quotes (multi-line "
                      "strings are not supported), got '" + value + "'");
            else if (e->tokens) {
                if (! detail::is_plain_list(value))
                    defer("option '" + key + "' must be a list written as "
                          "'[d0, d1, ...]' of plain numbers, got '" + value + "'");
            }
            else if (! value.empty() && value.front() == '[')
                defer("option '" + key + "' is not a list, but its value starts with '[': '"
                      + value + "'");
        }
        if (deferred) throw *deferred;
        return raw;
    }

    // True if the key was set in the file, or has a default.
    std::size_t count(const std::string &name) const
    {
        const detail::Entry *e = find(name);
        return (e && (e->opt->count() > 0 || e->has_default)) ? 1 : 0;
    }

    // True if the key was not set in the file and its default is in effect.
    bool defaulted(const std::string &name) const
    {
        const detail::Entry *e = find(name);
        return e && e->has_default && e->opt->count() == 0;
    }

    // The value as text. A list reads back as "[a, b, c]".
    std::string str(const std::string &name) const
    {
        const detail::Entry *e = find(name);
        if (! e) throw std::out_of_range("unknown option " + name);
        return e->to_str(*e);
    }

    template<class T>
    T as(const std::string &name) const
    {
        const detail::Entry *e = find(name);
        if (! e) throw std::out_of_range("unknown option " + name);
        return *std::any_cast<T*>(e->ref);
    }

    // One block per parameter, in declaration order:
    //   --section.key TYPE [=default]
    //       help text
    std::string help() const
    {
        std::ostringstream os;
        for (const auto &name : order_) {
            const detail::Entry &e = entries_.at(name);
            os << "  --" << name << ' ' << e.opt->get_type_name();
            if (e.required) os << " (required)";
            else if (e.has_default) os << " (=" << e.opt->get_default_str() << ')';
            os << '\n';
            std::istringstream lines(e.opt->get_description());
            std::string line;
            while (std::getline(lines, line))
                os << "      " << line << '\n';
            os << '\n';
        }
        return os.str();
    }

private:
    void add(const std::string &name, std::unique_ptr<detail::SpecBase> spec,
             const std::string &help)
    {
        const auto dot = name.find('.');
        if (dot == std::string::npos || entries_.count(name))
            throw std::logic_error("bad or duplicate option name " + name);
        const std::string sec = name.substr(0, dot);

        CLI::App *section = sections_.count(sec) ? sections_[sec]
                          : (sections_[sec] = app_.add_subcommand(sec)->configurable());
        order_.push_back(name);
        detail::Entry &e = entries_[name];
        spec->attach(section, name.substr(dot+1), help, e);
        specs_.push_back(std::move(spec));
    }

    const detail::Entry* find(const std::string &name) const
    {
        auto it = entries_.find(name);
        return it == entries_.end() ? nullptr : &it->second;
    }

    CLI::App app_{"Config file options"};
    std::map<std::string, CLI::App*> sections_;
    std::map<std::string, detail::Entry> entries_;
    std::vector<std::string> order_;
    std::vector<std::unique_ptr<detail::SpecBase>> specs_;
};

}  // namespace opt

#endif
