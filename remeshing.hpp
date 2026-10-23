#ifndef DYNEARTHSOL3D_REMESHING_HPP
#define DYNEARTHSOL3D_REMESHING_HPP

int bad_mesh_quality(const Param&, const Variables&, int&, double&);
void remesh(const Param&, Variables&, int);
void initialize_elem_size_n(const Variables&, double_vec &);
void detect_side_profile(const Param&, Variables&, uint side_bits);
double side_wall_top_z(const Variables&, int t);

// Along-strike coordinate of point x (indexed by axis) on restored side wall t: y on the x walls,
// x on the y walls. A 2D wall is a line: 0.
template <class X>
inline double side_wall_strike(const X &x, int t)
{
#ifdef THREED
    return x[1 - SIDEWALL_AXIS[t]];
#else
    (void)x; (void)t;
    return 0.;
#endif
}

// SideProfile lookups at (s, rel): s along strike (side_wall_strike), rel below the column top
// (SideProfile). Both locate the point as side_profile_locate does and assert a recorded column.
// The recorded material of the profile element at (s, rel).
int side_profile_mattype_at(const SideProfile&, double s, double rel);
// Marks in in_win (nmat entries) every material the profile assigns at s within depths [lo, hi].
void side_profile_mattypes_within(const SideProfile&, double s, double lo, double hi,
                                  std::vector<char> &in_win);
// Builds the derived lookup index; call after a profile is recorded or read back.
void side_profile_index_elems(SideProfile&);

#ifdef USEMMG
struct MMGInput;
struct MMGOutput;
// The floor of the size targets in the metric handed to MMG.
double mmg_metric_floor(const Mesh&);
// Element sizes and per-node plastic-strain size targets of any mesh (R5, the initial weak zone).
void plastic_size_targets(double coeff, const conn_t &connectivity, const double_vec &volume,
                          const double_vec &plstrain, const double_vec &init_size,
                          int nnode, int nelem, double_vec &esize, double_vec &target);
// R5's oversized test: an element at least twice its own size target.
inline bool r5_oversized(double esize, double target) { return esize >= 2.0 * target; }
// mmg_adapt under the post-output quality gate (FREEZE / REQUIRED POLICY), for remesh() and the
// initial refinements (mesh.cxx). req_node/req_elem must be the vectors in.required_node /
// required_elem point into (asserted); the gate edits them in place.
void mmg_adapt_quality_gated(const Param&, const MMGInput &in, MMGOutput &out,
                             std::vector<char> &req_node, std::vector<char> &req_elem);
#ifndef THREED
// Bisects the over-long required interface edges (iface) and their triangles before MMG; grows
// the arrays in place. Called by remesh() (pull_down) and refine_initial_side_walls.
void presubdivide_interface_edges(std::vector<double> &coord, std::vector<int> &conn,
                                  std::vector<double> &metric, std::vector<char> &reqn,
                                  std::vector<char> &reqe, std::vector<int> &iface,
                                  int &nn, int &ne, std::vector<int> *elem_ref, bool pull_down);
#endif
#endif

#endif
