#include <cstdlib>
#include <iostream>
#ifdef USEMPI
#include <mpi.h>
#endif

#include "utils.hpp"

namespace {

void abort_parallel_run(ExitCode code)
{
#ifdef USEMPI
    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Finalized(&finalized);
    // A worker may be waiting for the next force call when rank zero fails in MMG or I/O.
    if (initialized && !finalized && code != EXIT_OK) MPI_Abort(MPI_COMM_WORLD, code);
#endif
}

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
    abort_parallel_run(code);
    std::exit(code);
}

void die(ExitCode code, const char* msg)
{
    std::cerr << "[DES exit " << int(code) << "] " << exit_category(code)
              << ": " << msg << '\n';
    abort_parallel_run(code);
    std::exit(code);
}

void die_together(ExitCode code, const char* msg)
{
#ifdef USEMPI
    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Finalized(&finalized);
    int rank = 0;
    if (initialized && !finalized) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
        // Rank zero never enters this barrier: its MPI_Abort ends the wait.
        MPI_Barrier(MPI_COMM_WORLD);
        std::exit(code);
    }
#endif
    if (msg) die(code, msg);
    die(code);
}
