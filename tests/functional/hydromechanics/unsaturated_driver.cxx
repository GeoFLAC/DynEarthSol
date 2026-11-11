// Command driver for check_unsaturated_header.py; reads stdin, one command
// per line, and prints results with 17 significant digits.
//   eval model alpha n sr p        -> se dse sat kr bishop
//   pinv model alpha n sr se       -> p
//   advance hmass p dm count {model alpha n sr weight}*count -> x change(x) invalid
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "unsaturated.hpp"

namespace {

struct Patch
{
    std::vector<Retention> laws;
    std::vector<double> weights;
    int size() const { return static_cast<int>(laws.size()); }
    const Retention &law(int k) const { return laws[k]; }
    double weight(int k) const { return weights[k]; }
};

// strtod also accepts nan and inf, which operator>> rejects.
double num(std::istream &in)
{
    std::string s;
    in >> s;
    return std::strtod(s.c_str(), nullptr);
}

Retention read_law(std::istream &in)
{
    Retention r;
    r.model = static_cast<int>(num(in));
    r.alpha = num(in);
    r.n = num(in);
    r.sr = num(in);
    return r;
}

} // namespace

int main()
{
    std::string cmd;
    while (std::cin >> cmd) {
        if (cmd == "eval") {
            const Retention r = read_law(std::cin);
            const double p = num(std::cin);
            double dse;
            const double se = retention_se(r, p, &dse);
            std::printf("%.17g %.17g %.17g %.17g %.17g\n", se, dse,
                        retention_sat(r, p), retention_kr(r, p), bishop_pressure(r, p));
        } else if (cmd == "pinv") {
            const Retention r = read_law(std::cin);
            std::printf("%.17g\n", retention_p_of_se(r, num(std::cin)));
        } else if (cmd == "advance") {
            const double hmass = num(std::cin), p = num(std::cin), dm = num(std::cin);
            const int count = static_cast<int>(num(std::cin));
            Patch patch;
            for (int k = 0; k < count; ++k) {
                patch.laws.push_back(read_law(std::cin));
                patch.weights.push_back(num(std::cin));
            }
            int invalid = 0;
            const double x = advance_content(patch, hmass, p, dm, invalid);
            std::printf("%.17g %.17g %d\n", x, content_change(patch, hmass, p, x), invalid);
        } else {
            std::fprintf(stderr, "unknown command %s\n", cmd.c_str());
            return 1;
        }
    }
    return 0;
}
