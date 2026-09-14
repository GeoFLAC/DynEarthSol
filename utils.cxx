#include <cstdlib>
#include <iostream>

#include "utils.hpp"

namespace {

const char* exit_category(int code)
{
    switch (code / 10) {
    case 0: return "normal exit";
    case 1: return "user input";
    case 2: return "I/O";
    case 3: return "unsupported in this build";
    case 4: return "mesh generation";
    case 5: return "runtime";
    case 6: return "internal error";
    default: return "unknown";
    }
}

} // namespace

void die(ExitCode code)
{
    std::cerr << "[DES exit " << int(code) << "] " << exit_category(code) << '\n';
    std::exit(code);
}

void die(ExitCode code, const char* msg)
{
    std::cerr << "[DES exit " << int(code) << "] " << exit_category(code)
              << ": " << msg << '\n';
    std::exit(code);
}
