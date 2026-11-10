#ifndef DYNEARTHSOL3D_BRC_INTERPOLATION_HPP
#define DYNEARTHSOL3D_BRC_INTERPOLATION_HPP

struct RemeshScratch;

// Reads rs.affected_old_node and fills rs.node_outside (remesh-scratch.hpp).
void barycentric_node_interpolation(const Param& param, Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity,
                                    RemeshScratch &rs);

void barycentric_node_interpolation_forT(const Param& param, const Variables &var,
                                         const Barycentric_transformation &bary,
                                         const array_t &input_coord,
                                         const conn_t &input_connectivity,
                                         const Support &input_support,
                                         const double_vec &inputtemperature,
                                         double_vec &outputtemperature);

#endif
