#ifndef DYNEARTHSOL3D_REMESHING_HPP
#define DYNEARTHSOL3D_REMESHING_HPP

int bad_mesh_quality(const Param&, const Variables&, int&, double&);
void remesh(const Param&, Variables&, int);
void initialize_elem_size_n(const Variables&, double_vec &);
void detect_side_profile(const Param&, Variables&, uint side_bits);

// Depth lookups in a recorded SideProfile (parameters.hpp); depths relative to the side's
// top point. Used by the post-remesh field restore (remeshing.cxx) and the side-wall
// marker replenishment (markerset.cxx). node_value interpolates any (depth, value) table
// sorted by depth (node_reldepth vs node_temperature, a coord0 component, ...).
double side_profile_node_value(const double_vec &depths, const double_vec &values,
                               double rel, double tol);
void side_profile_elems_at(const SideProfile&, double rel, double tol, int_vec &out);

#endif
