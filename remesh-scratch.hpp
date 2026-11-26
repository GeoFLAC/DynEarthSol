#ifndef DYNEARTHSOL3D_REMESH_SCRATCH_HPP
#define DYNEARTHSOL3D_REMESH_SCRATCH_HPP

#include <memory>
#include <vector>

#include "parameters.hpp"

// Per NEW element: if its part inside the old mesh lies in exactly one old element whose markers
// are all one mattype, that mattype (mat), else -1; time: that old element's lowest-index marker's
// time. Filled by nearest_neighbor_interpolation, read by remap_markers' replenishment.
struct SplitParents {
    std::vector<int> mat;
    std::vector<double> time;
};

// What one remesh() passes between its steps and nothing outside it reads. remesh() owns it,
// so all of it is freed when the remesh returns.
struct RemeshScratch {
    // The stress-remap transients, on the heap: OpenACC kernels read and write them through
    // pointers, which reach managed memory but not the host stack. stress_n and spr_blend_weight
    // double as the switch for the SPR recovery: remesh() leaves them null when the rheology has
    // no viscous component, and each consumer skips its part of the chain on that.
    // - stress_n, stressyy_n: SPR-recovered nodal stress
    // - spr_blend_weight: Deborah weight toward the NN-remapped stress (1 = NN, 0 = SPR)
    // - spr_p_ref_old: reference pressure added at pressure-centering, per element
    // - remesh_is_changed: the NN pass's per-NEW-element is_changed mapping
    //   (0 = geometry identical to an old element, 1 = remapped, -1 = ACM failed)
    std::unique_ptr<tensor_t> stress_n;
    std::unique_ptr<double_vec> stressyy_n;
    std::unique_ptr<double_vec> spr_blend_weight;
    std::unique_ptr<double_vec> spr_p_ref_old;
    std::unique_ptr<int_vec> remesh_is_changed;

    // Per OLD node: 1 if this remesh moved (flatten_*), deleted or collapsed it
    // (record_affected_old_nodes); read by barycentric_node_interpolation's "interior node not
    // found" warning.
    std::vector<char> affected_old_node;

    // Per NEW node: 1 if barycentric_node_interpolation found it outside the old mesh and copied
    // the nearest old node's fields. restore_side_fields (option 13) and remesh() overwrite some
    // of those copies (see there); the rest stay.
    std::vector<char> node_outside;

    SplitParents split_parent;
};

#endif
