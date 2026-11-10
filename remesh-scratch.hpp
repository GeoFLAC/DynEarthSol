#ifndef DYNEARTHSOL3D_REMESH_SCRATCH_HPP
#define DYNEARTHSOL3D_REMESH_SCRATCH_HPP

#include <vector>

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
