#ifndef DES_EXIT_CODES_HPP
#define DES_EXIT_CODES_HPP

// Process exit codes: first digit the category, second the cause, so a bare status
// says where to look -- 1x the user to fix, 2x the environment, 3x-6x ours. Keep in
// step with CONTRIBUTING.md. Nothing above 99; the shell reserves 126+.
enum ExitCode {
    EXIT_OK                   =  0,

    EXIT_CONFIG               = 10,  // generic config error
    EXIT_CONFIG_VALUE         = 11,  // value out of range / unknown option
    EXIT_CONFIG_DATA          = 12,  // malformed input data file (.poly, .exo)

    EXIT_IO_OPEN              = 20,  // cannot open file
    EXIT_IO_RW                = 21,  // read/write failed (incl. HDF5)
    EXIT_IO_RESTART           = 22,  // restart/checkpoint mismatch

    EXIT_UNSUPPORTED_DIM      = 30,  // not implemented for this NDIMS
    EXIT_UNSUPPORTED_LIB      = 31,  // optional library not compiled in

    EXIT_MESH_TETGEN          = 40,  // Triangle / TetGen
    EXIT_MESH_MMG             = 41,  // MMG
    EXIT_MESH_QUALITY         = 42,  // mesh quality / topology

    EXIT_RUNTIME_NAN          = 50,  // NaN or non-finite state
    EXIT_RUNTIME_LOOKUP       = 51,  // marker/geometry lookup failed
    EXIT_RUNTIME_RESOURCE     = 52,  // resource exhausted
    EXIT_RUNTIME_GOSPL        = 53,  // GoSPL step failed
    EXIT_RUNTIME_NONCONVERGENCE = 54, // mechanical solve did not converge

    EXIT_INTERNAL_ASSERT      = 60,  // assertion / invariant violated
    EXIT_INTERNAL_UNREACHABLE = 61   // unreachable branch
};

#endif
