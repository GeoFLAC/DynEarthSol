#ifndef DYNEARTHSOL3D_MMG_UTILS_HPP
#define DYNEARTHSOL3D_MMG_UTILS_HPP

#ifdef USEMMG

#include <vector>
#include "parameters.hpp"

// Mesh + metric (+ optional freeze flags) handed to a single MMG adaptation pass.
// All index arrays are 0-indexed; mmg_adapt makes its own 1-indexed copies internally.
struct MMGInput {
    int nnode, nelem, nseg;
    const double *coord;          // AoS, nnode*NDIMS
    const int    *conn;           // nelem*NODES_PER_ELEM  (tetra 3D / triangles 2D)
    const int    *seg;            // nseg*NODES_PER_FACET   (boundary triangles 3D / edges 2D)
    const int    *segflag;        // nseg boundary references
    const double *metric;         // nnode isotropic target edge length
    const char   *required_node;  // nnode; 1 = frozen (MMG-required). nullptr => none
    const char   *required_elem;  // nelem; 1 = frozen (MMG-required). nullptr => none
    bool tolerate_low_failure;    // true: keep the saved mesh on MMG5_LOWFAILURE (init); false: abort (remesh)
};

// Adapted mesh returned by mmg_adapt (0-indexed, flat buffers ready for load_from_buffer).
struct MMGOutput {
    int nnode, nelem, nseg;
    std::vector<double> coord;    // nnode*NDIMS
    std::vector<int>    conn;     // nelem*NODES_PER_ELEM
    std::vector<int>    seg;      // nseg*NODES_PER_FACET
    std::vector<int>    segflag;  // nseg
};

// Run one MMG (2D or 3D, per build) adaptation of `in` into `out`. The size band and solver
// parameters come from `mesh` (largest_size/smallest_size/mmg_hausd_factor/resolution/
// mmg_verbose/mmg_debug): the single owner of the MMG lifecycle, shared by mesh init and
// remeshing. Dies with EXIT_MESH_MMG on any MMG failure.
void mmg_adapt(const Mesh &mesh, const MMGInput &in, MMGOutput &out);

#endif // USEMMG
#endif
