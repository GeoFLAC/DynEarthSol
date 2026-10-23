#ifndef DYNEARTHSOL3D_MMG_UTILS_HPP
#define DYNEARTHSOL3D_MMG_UTILS_HPP

#ifdef USEMMG

#include <vector>
#include "parameters.hpp"

// The input of one mmg_adapt pass. Index arrays are 0-based; mmg_adapt never writes them.
struct MMGInput {
    int nnode, nelem, nseg;
    const double *coord;          // AoS, nnode*NDIMS
    const int    *conn;           // nelem*NODES_PER_ELEM  (tetra 3D / triangles 2D)
    const int    *seg;            // nseg*NODES_PER_FACET   (boundary triangles 3D / edges 2D)
    const int    *segflag;        // nseg boundary references
    const double *metric;         // nnode isotropic target edge length
    const char   *required_node;  // nnode; 1 = pinned (MMG-required vertex). nullptr => none
    const char   *required_elem;  // nelem; 1 = frozen (MMG-required). nullptr => none
    bool tolerate_low_failure;    // true: keep the saved mesh on MMG5_LOWFAILURE (init); false: abort (remesh)
    // Optional (nullptr => scalar `metric`): nnode symmetric tensors, row-major per node (2D:
    // m11,m12,m22; 3D: m11,m12,m13,m22,m23,m33), used INSTEAD of `metric`.
    const double *metric_aniso;
    // Optional, 2D only: n_req_edges node pairs kept as required interface edges (passed with a
    // sentinel ref, dropped from out.seg). nullptr/0 => none.
    const int *req_edges;
    int n_req_edges;
    // Optional, nnode: 1 = MMG starts from the metric at this node even on a required edge.
    const char *keep_metric;
    // Optional, nelem, 2D only: element references; each new element gets its parent's
    // (out.elem_ref), and an edge between two references stays as an internal boundary.
    const int *elem_ref;
    // Optional, nnode, 2D only: corners, never moved or removed; a boundary or reference curve
    // runs straight into one.
    const char *corner_node;
};

// Adapted mesh returned by mmg_adapt (0-indexed, flat buffers ready for load_from_buffer).
struct MMGOutput {
    int nnode, nelem, nseg;
    std::vector<double> coord;    // nnode*NDIMS
    std::vector<int>    conn;     // nelem*NODES_PER_ELEM
    std::vector<int>    seg;      // nseg*NODES_PER_FACET
    std::vector<int>    segflag;  // nseg
    std::vector<int>    elem_ref; // nelem when MMGInput::elem_ref is set, else empty
};

// Run one MMG adaptation (2D or 3D, per build) of `in` into `out`, for mesh init and remeshing.
// Size band from `mesh`: hmax = largest_size^(1/NDIMS) x resolution, hmin = (smallest_size x
// mmg_remesh_defensive_tiny_ratio)^(1/NDIMS) x resolution (MMG's floor above the tiny-element
// trigger), hausd = mmg_hausd_factor x resolution; optim 0. Every node of a required element must
// be required (asserted). At a node on a required edge (a required element's or req_edges) MMG
// starts from the mean length of the required edges there, unless keep_metric marks the node
// (apply_required_edge_size). out.seg holds boundary facets only (drop_unflagged_facets). Dies
// with EXIT_MESH_MMG on any MMG failure, except MMG5_LOWFAILURE under tolerate_low_failure.
void mmg_adapt(const Mesh &mesh, const MMGInput &in, MMGOutput &out);

#endif // USEMMG
#endif
