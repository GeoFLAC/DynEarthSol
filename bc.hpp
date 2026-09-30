#ifndef DYNEARTHSOL3D_BC_HPP
#define DYNEARTHSOL3D_BC_HPP

bool is_on_boundary(const Variables &var, int node);
double find_max_vbc(const BC &bc);
void create_boundary_normals(const Variables &var, array_t &bnormals,
                             double_vec& edge_vec, int* edge_slot);
// Homogeneous/normalized modes are used to construct the corrected PT subspace.
// Default calls retain the legacy boundary map.
void apply_vbcs(const Param &param, const Variables &var, array_t &vel,
                bool homogeneous = false, bool normalized = false);

struct VelocityConstraints {
    Array2D<double, NDIMS*NDIMS> projector;
    array_t prescribed;
    int_vec free_rank;
    explicit VelocityConstraints(int n) : projector(n), prescribed(n), free_rank(n) {}
};
void build_velocity_constraints(const Param& param, const Variables& var,
                                VelocityConstraints& constraints, bool homogeneous = false);
void project_free_vectors(const Variables& var, const VelocityConstraints& constraints,
                          array_t& vectors, bool include_prescribed = false);
void apply_stress_bcs(const Param& param, const Variables& var, array_t& force,
                      double displacement_dt = 0);
void apply_stress_bcs_neumann(const Param& param, const Variables& var, array_t& force);
void surface_plstrain_diffusion(const Param &param, const Variables& var, double_vec& plstrain);
void correct_surface_element(const Variables& var, double_vec& volume, double_vec& volume_n, tensor_t& stress, \
                              tensor_t& strain, tensor_t& strain_rate, double_vec& plstrain);
void surface_processes(const Param& param, const Variables& var, array_t& coord, tensor_t& stress, tensor_t& strain, \
                       tensor_t& strain_rate, double_vec& plstrain, double_vec& volume, double_vec& volume_n, \
                       SurfaceInfo& surfinfo, std::vector<MarkerSet*> &markersets, \
                       int_vec2D& elemmarkers, int_vec2D& markers_in_elem);

#endif
