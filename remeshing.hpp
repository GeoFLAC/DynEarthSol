#ifndef DYNEARTHSOL3D_REMESHING_HPP
#define DYNEARTHSOL3D_REMESHING_HPP

int bad_mesh_quality(const Param&, const Variables&, int&, double&);
void remesh(const Param&, Variables&, int);
void initialize_elem_size_n(const Variables&, double_vec &);
void detect_side_profile(const Param&, Variables&, uint side_bits);

// Depth lookups in a recorded SideProfile (parameters.hpp); depths relative to the side's
// top point. Used by the post-remesh field restore (remeshing.cxx); side_profile_elems_at also
// by the side-wall marker replenishment (markerset.cxx). node_value interpolates any (depth, value) table
// sorted by depth (node_reldepth vs node_temperature, node_coord0_z, ...).
// side_profile_elems_at requires a recorded column (asserted) and always returns at least one
// element.
double side_profile_node_value(const double_vec &depths, const double_vec &values,
                               double rel, double tol);
void side_profile_elems_at(const SideProfile&, double rel, double tol, int_vec &out);
int side_profile_mattype_at(const SideProfile&, double rel, double tol, int_vec &cand);

// Same-depth tolerance of the side-profile lookups (tie runs, division levels, the division
// snap's "already there" test). A degeneracy guard, untuned: 1e-6 x resolution lies far above
// coordinate round-off (~1e-16 x domain size) and far below any node spacing (about hmin, a
// fraction of resolution), and the lookups do not depend on it in that range. In 3D, wall
// nodes tie into a row only where their depths agree to round-off: every row of a structured
// wall under a level top, and on any wall a flat bottom edge and a level top edge.
inline double side_profile_depth_tol(const Param &param) { return 1e-6 * param.mesh.resolution; }

#endif
