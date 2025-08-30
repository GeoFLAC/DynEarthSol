#ifndef DYNEARTHSOL3D_NN_INTERPOLATION_HPP
#define DYNEARTHSOL3D_NN_INTERPOLATION_HPP

void nearest_neighbor_interpolation(const Param& Param, Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity,
                                    const bool is_surface = false,
                                    const int_pair_vec &old_bfacets_surface = int_pair_vec());

#endif
