#ifndef DYNEARTHSOL3D_NN_INTERPOLATION_HPP
#define DYNEARTHSOL3D_NN_INTERPOLATION_HPP

typedef Array2D<int, 32> nn_t;
typedef Array2D<double, 32> ratio_t;

struct RemeshScratch;

// The element pass (is_surface false) fills rs.split_parent and rs.remesh_is_changed and carries
// rs.spr_blend_weight and rs.spr_p_ref_old to the new elements (remesh-scratch.hpp).
void nearest_neighbor_interpolation(const Param& Param, Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity,
                                    RemeshScratch &rs,
                                    const bool is_surface = false);

#endif
