#ifndef DYNEARTHSOL3D_BRC_INTERPOLATION_HPP
#define DYNEARTHSOL3D_BRC_INTERPOLATION_HPP

void barycentric_node_interpolation(Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity);

void barycentric_node_interpolation_forT(const Variables &var,
                                         const Barycentric_transformation &bary,
                                         const array_t &input_coord,
                                         const conn_t &input_connectivity,
                                         const int_vec2D &input_support,
                                         const double_vec &inputtemperature,
                                         double_vec &outputtemperature);

#endif
