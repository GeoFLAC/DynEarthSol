#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_map>
#include <iomanip>
#include <math.h>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "markerset.hpp"
#include "barycentric-fn.hpp"
#include "geometry.hpp"
#include "utils.hpp"
#ifdef HAS_GOSPL_CPP_INTERFACE
#include "gospl_driver/gospl-driver.hpp"
#endif


#include "bc.hpp"


namespace {

#pragma acc routine seq
void normal_vector_of_facet(ConstArrayIndirectAccessor facet_coord, double *normal, double &zcenter)
{
#ifdef THREED
    // vectors: n0-n1 and n0-n2
    double v01[NDIMS], v02[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = facet_coord[1][i] - facet_coord[0][i];
        v02[i] = facet_coord[2][i] - facet_coord[0][i];
    }

    // the outward normal vector is parallel to the cross product
    // of v01 and v02, with magnitude = area of the triangle
    normal[0] = (v01[1] * v02[2] - v01[2] * v02[1]) / 2;
    normal[1] = (v01[2] * v02[0] - v01[0] * v02[2]) / 2;
    normal[2] = (v01[0] * v02[1] - v01[1] * v02[0]) / 2;

    zcenter = (facet_coord[0][2] + facet_coord[1][2] + facet_coord[2][2]) / NODES_PER_FACET;
#else
    // the edge vector
    double v01[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = facet_coord[1][i] - facet_coord[0][i];
    }

    // the normal vector to the edge, pointing outward
    normal[0] = v01[1];
    normal[1] = -v01[0];

    zcenter = (facet_coord[0][1] + facet_coord[1][1]) / NODES_PER_FACET;
#endif
}

}

#pragma acc routine seq
bool is_on_boundary(const uint_vec &bcflag, int node)
{
    uint flag = bcflag[node];
    return flag & BOUND_ANY;
}


namespace {
// Code 0 on a side: no velocity component held (tractions aside).
bool all_free(const SideBC &s)
{
    for (int d=0; d<NDIMS; ++d)
        if (s.type[d] == 'v') return false;
    return true;
}

// Codes 1 and 11: the normal velocity held, the other components not.
bool normal_only(const SideBC &s)
{
    if (s.type[0] != 'v') return false;
    for (int d=1; d<NDIMS; ++d)
        if (s.type[d] == 'v') return false;
    return true;
}
} // anonymous namespace


// Read the velocity and stress codes into one table per side, in the side's frame
// (see SideBC). Zeros on a negative axis are stored as -0.0, so that turning a value
// back into a global component gives +0.0.
void convert_bc_codes(Param &p)
{
    // What each code holds, per global component.
    // N is the normal value, which follows the period table; V the value; L the lateral
    // value; Z holds the component at 0; F leaves it free. Side walls list (normal,
    // lateral, z) in 3D and (x, z) in 2D; z sides list (x, y, z) or (x, z).
    enum { F, Z, N, V, L };
#ifdef THREED
    static const int wall_rows[8][NDIMS] = {{F, F, F}, {N, F, F}, {F, Z, Z}, {N, Z, Z},
                                            {F, V, Z}, {Z, V, Z}, {N, L, F}, {N, Z, F}};
    static const int z_rows[4][NDIMS] = {{F, F, F}, {F, F, V}, {Z, Z, F}, {Z, Z, V}};
#else
    static const int wall_rows[8][NDIMS] = {{F, F}, {N, F}, {F, Z}, {N, Z}, {Z, V}, {F, F}, {N, L}, {F, F}};
    static const int z_rows[5][NDIMS] = {{F, F}, {F, V}, {Z, F}, {Z, V}, {V, Z}};
#endif

    BC &bc = p.bc;
    for (int i=0; i<nbdrytypes; ++i) {
        SideBC &s = bc.side[i];
        s = SideBC{};
        for (int d=0; d<NDIMS; ++d) s.type[d] = 'f';
        s.release_time = (i == iboundz1) ? bc.vbc_val_z1_loading_period : std::numeric_limits<double>::max();
        s.stress_axis = -1;
        const int code = bc.*side_vbc_code[i];
        const double v = bc.*side_vbc_value[i];

        const std::string &text = bc.vbc_text[i];
        if (!text.empty() && std::isalpha(static_cast<unsigned char>(text[0]))) {
            // Given by components: the values are already in the side's frame. A normal not
            // held carries the support; only a held normal value follows the period table, and
            // in 2D the depth profile and the bottom shear zone of x0.
            for (int d=0; d<NDIMS; ++d) {
                if (text[d] == 'f') continue;
                s.type[d] = text[d];
                s.value[d] = bc.vbc_component_values[i][d];
                if (text[d] == 'v') s.max_speed = std::max(s.max_speed, std::fabs(s.value[d]));
            }
            const bool free_normal = (text[0] != 'v');
            s.support = (i == iboundz0) ? free_normal && bc.has_winkler_foundation
                      : (i == iboundz1) ? free_normal && bc.has_water_loading : free_normal;
            if (i < iboundz0) {
                s.normal_follows_period = !free_normal;
                s.taper = (NDIMS == 2 && i == iboundx0 && !free_normal);
            }
            continue;
        }

        if (i >= iboundn0) {
            // n0 to n3 accept 1, 3, 11 and 13 only (input.cxx checks). The frame is the
            // boundary's normal, or its horizontal part for 11 and 13.
            s.horizontal_normal = (code == 11 || code == 13);
            const bool all = (code == 3 || code == 13);
            for (int d=0; d<NDIMS; ++d) {
                if (d > 0 && !all) break;
                s.type[d] = 'v';
                s.value[d] = d ? 0. : v;
            }
            s.max_speed = std::fabs(v);
            continue;
        }

        // stress_bc_*: the global stress component k, applied on each facet as S n_k with the
        // facet's area-weighted normal n (input.cxx refuses what adds nothing).
        const int stress_code = bc.*side_stress_code[i];
        if (stress_code >= 1 && stress_code <= 3) {
#ifdef THREED
            s.stress_axis = stress_code - 1;
#else
            if (stress_code != 2) s.stress_axis = (stress_code == 1) ? 0 : 1;
#endif
            s.stress_value = bc.*side_stress_value[i];
        }
#ifndef THREED
        if (i == iboundy0 || i == iboundy1) continue;
#endif

        const bool z_side = (i == iboundz0 || i == iboundz1);
        const int *row = nullptr;
        int global[NDIMS];  // the global axis of each row entry
        if (z_side) {
            if (code >= 0 && code < int(sizeof(z_rows) / sizeof(z_rows[0]))) row = z_rows[code];
            for (int d=0; d<NDIMS; ++d) global[d] = d;
        }
        else {
            // input.cxx refuses the codes that set nothing: anything else than 0 to 7, and 5
            // and 7 in 2D.
            const bool sets_nothing = (code < 0 || code > 7)
#ifndef THREED
                || code == 5 || code == 7
#endif
                ;
            if (!sets_nothing) row = wall_rows[code];
#ifdef THREED
            const int normal = (i < iboundy0) ? 0 : 1;
            global[0] = normal; global[1] = 1 - normal; global[2] = 2;
#else
            global[0] = 0; global[1] = 1;
#endif
            s.taper = (NDIMS == 2 && i == iboundx0 && code == 3);
            s.normal_follows_period = (row && row[0] == N);
        }
        if (!row) continue;

        for (int d=0; d<NDIMS; ++d) {
            const int axis = frame_axis(i, d);
            const double sign = frame_sign(i, d);
            int k = 0;
            while (global[k] != axis) ++k;
            if (row[k] == F) continue;
            const double held = (row[k] == Z) ? 0. : (row[k] == L) ? bc.*side_vbc_lateral[i] : v;
            s.type[d] = 'v';
            s.value[d] = sign * held;
            s.max_speed = std::max(s.max_speed, std::fabs(held));
        }
        // A free normal carries the support, as with letters.
        const bool free_normal = (s.type[0] == 'f');
        s.support = (i == iboundz0) ? free_normal && bc.has_winkler_foundation
                  : (i == iboundz1) ? free_normal && bc.has_water_loading : free_normal;
    }
}


double find_max_vbc(const BC &bc)
{
    double max_vbc_val = 1e-12; // equal to 0.03 mm/yr
    for (int i=0; i<nbdrytypes; ++i)
        max_vbc_val = std::max(max_vbc_val, bc.side[i].max_speed);
    return max_vbc_val;
}


void create_boundary_normals(const Param &param, const Variables &var, array_t &bnormals,
                             double_vec& edge_vec, int* edge_slot)
{
    /* This subroutine finds the outward normal unit vectors of boundaries.
     * There are two types of boundaries: ibound{x,y,z}? and iboundn?.
     * The first type ibound{x,y,z}? has well defined "normal" direction
     * as their name implied. Thus they can be curved. (Ex. the top boundary
     * can have topography.)
     * The second type iboundn?, on the other hand, must be a plane so that
     * their vbc are well defined.
     *
     * The normals are those of the starting mesh: the run keeps them, and so do
     * its checkpoints, even where a boundary bends later. A side without
     * facets gets 0.
     */

    #pragma acc wait

    for (int i=0; i<nbdrytypes; i++) {
        double normal[NDIMS] = {0};
        if (var.bfacets[i]->size() == 0) {
            for (int d=0; d<NDIMS; d++) bnormals[i][d] = 0;
            continue;
        }

        for (auto j=var.bfacets[i]->begin(); j<var.bfacets[i]->end(); ++j) {
            int e = j->first;
            int f = j->second;
            double tmp;

            int idx[NODES_PER_FACET];
            for (int k=0; k<NODES_PER_FACET; ++k)
                idx[k] = (*var.connectivity)[e][NODE_OF_FACET[f][k]];

            ConstArrayIndirectAccessor facet_coord = var.coord->view_const(idx);

            normal_vector_of_facet(facet_coord, normal, tmp);

            // make an unit vector
            double len = 0;
            for(int d=0; d<NDIMS; d++)
                len += normal[d] * normal[d];

            len = std::sqrt(len);
            for(int d=0; d<NDIMS; d++)
                normal[d] = normal[d] / len;

            if (j == var.bfacets[i]->begin()) {
                for(int d=0; d<NDIMS; d++)
                    bnormals[i][d] = normal[d];

                if (i < iboundn0) break; // slant boundaries start at iboundn0, other boundary can be curved
            }
            else {
                // Make sure the boundary is a plane, ie. all facets have the same normal vector.
                const double eps2 = 1e-12;
                double diff2 = 0;
                for(int d=0; d<NDIMS; d++)
                    diff2 += (bnormals[i][d] - normal[d]) * (bnormals[i][d] - normal[d]);
                if (diff2 > eps2) {
                    std::cerr << "Error: boundary " << i << " is curved.\n";
                    std::cerr << "Got  -  Expected\n";
                    std::cerr << bnormals[i][0] << " - " << normal[0] << '\n';
                    std::cerr << bnormals[i][1] << " - " << normal[1] << '\n';
#ifdef THREED
                    std::cerr << bnormals[i][2] << " - " << normal[2] << '\n';
#endif
                    die(EXIT_MESH_QUALITY);
                }
            }
        }
    }

    create_boundary_edges(param, var, edge_vec, edge_slot);
}


void create_boundary_edges(const Param &param, const Variables &var, double_vec& edge_vec,
                           int* edge_slot)
{
    // The edge direction of each pair of boundaries, from var.bnormals.
    // -1 = this pair of boundaries shares no edge. Only pairs with facets on both sides
    // get an entry, so the table is sparse and the absent case has to be representable.
    std::fill_n(edge_slot, nbdrytypes * nbdrytypes, -1);
    edge_vec.clear();

    for (int i=0; i<nbdrytypes; i++) {
        if (var.bfacets[i]->size() == 0) continue;

        const double eps = 1e-15;
        for (int j=i+1; j<nbdrytypes; j++) {
            if (var.bfacets[j]->size() == 0) continue;
            double_vec s (NDIMS);  // intersection of two boundaries
                                            // whole-application lifetime, no need to delete manually
#ifdef THREED
            // quick path: both walls are vertical
            if (std::abs((*var.bnormals)[i][NDIMS-1]) < eps &&
                std::abs((*var.bnormals)[j][NDIMS-1]) < eps) {
                s[0] = s[1] = 0;
                s[NDIMS-1] = 1;
            }
            else {
                // cross product of 2 normal vectors
                s[0] = (*var.bnormals)[i][1]*(*var.bnormals)[j][2] - (*var.bnormals)[i][2]*(*var.bnormals)[j][1];
                s[1] = (*var.bnormals)[i][2]*(*var.bnormals)[j][0] - (*var.bnormals)[i][0]*(*var.bnormals)[j][2];
                s[2] = (*var.bnormals)[i][0]*(*var.bnormals)[j][1] - (*var.bnormals)[i][1]*(*var.bnormals)[j][0];
            }
#else
            // 2D
            s[0] = 0;
            s[1] = 1;
#endif
            edge_slot[i*nbdrytypes + j] = edge_vec.size() / NDIMS;
            edge_vec.push_back(s[0]);
            edge_vec.push_back(s[1]);
#ifdef THREED
            edge_vec.push_back(s[2]);
#endif
        }
    }

    // apply_vbcs projects the velocity onto this edge for any node flagged on both
    // boundaries of a pair, from inside a device kernel that can neither report nor
    // abort. A pair can be flagged on a node while the loop above skipped it -- either
    // boundary having no facets is enough -- and the constraint is then dropped in
    // silence. Report it here, the last place with somewhere to say it, for the pairs
    // the kernel can actually ask about.
    for (int ib=iboundn0; ib<=iboundn3; ib++) {
        if (!normal_only(param.bc.side[ib])) continue;
        for (int ic=iboundx0; ic<ib; ic++) {
            const SideBC &s = param.bc.side[ic];
            if (!normal_only(s) || s.horizontal_normal) continue;
            if (edge_slot[ic*nbdrytypes + ib] >= 0) continue;
            for (int n=0; n<var.nnode; ++n) {
                uint flag = (*var.bcflag)[n];
                if (!(flag & (1U << ib)) || !(flag & (1U << ic))) continue;
                std::cerr << "Warning: boundaries " << ic << " and " << ib << " meet at node "
                          << n << " but share no edge; apply_vbcs leaves the edge-parallel "
                          << "velocity constraint unapplied there\n";
                break;
            }
        }
    }
}


namespace {
#pragma acc routine seq
void project_boundary_intersection(const Variables& var, int ic, int ib,
                                   ConstArrayAccessor n, ArrayAccessor v, bool normalized,
                                   double prescribed_normal, bool horizontal_earlier,
                                   bool horizontal_current)
{
    if (normalized) {
        // Add the new affine normal constraint while preserving the normal
        // component fixed by the earlier boundary. This retains imposed motion
        // at edges; projecting onto an edge through the origin would erase it.
        double earlier[NDIMS] = {}, current[NDIMS] = {};
        if (ic < iboundn0) {
            const int direction = ic <= iboundx1 ? 0 : (ic <= iboundy1 ? 1 : NDIMS-1);
            earlier[direction] = 1;
        } else {
            for (int d=0; d<NDIMS; ++d) earlier[d] = (*var.bnormals)[ic][d];
            if (horizontal_earlier) earlier[NDIMS-1] = 0;
        }
        double aa = 0, bb = 0;
        for (int d=0; d<NDIMS; ++d) {
            current[d] = (horizontal_current && d == NDIMS-1) ? 0 : n[d];
            aa += earlier[d]*earlier[d];
            bb += current[d]*current[d];
        }
        if (!(aa > 0 && bb > 0)) return; // invalid normals rejected before kernel
        double ab = 0, vn = 0;
        for (int d=0; d<NDIMS; ++d) {
            earlier[d] /= std::sqrt(aa);
            current[d] /= std::sqrt(bb);
            ab += earlier[d]*current[d];
            vn += v[d]*current[d];
        }
        const double determinant = 1-ab*ab;
        const double delta = prescribed_normal-vn;
        for (int d=0; d<NDIMS; ++d)
            v[d] += determinant > 1e-12
                ? delta*(current[d]-ab*earlier[d])/determinant
                : delta*current[d]; // parallel normals: later boundary wins
        return;
    }
    // Preserve the original edge map exactly for non-PT callers.
    const int slot = var.edge_slot[ic*nbdrytypes + ib];
    if (slot < 0) return;
    const double* edge = &var.edge_vec[slot*NDIMS];
    double ve = 0;
    for (int d=0; d<NDIMS; ++d) ve += v[d]*edge[d];
    for (int d=0; d<NDIMS; ++d) v[d] = ve*edge[d];
}
}

// Normal velocity apply_vbcs prescribes on side wall t (0 = x0, 1 = x1, 2 = y0, 3 = y1) at
// time_in_yr, along the global axis, before the 2D depth profile and homogeneous constraints:
// the wall's normal value times its period table; 0 unless that value follows the table.
namespace {
double side_wall_normal_vbc(const BC &bc, int t, double time_in_yr)
{
    const bool lo = (t % 2 == 0), y = (t >= 2);
    const SideBC &s = bc.side[t];
    if (!s.normal_follows_period) return 0.;
    const double v = frame_sign(t, 0) * s.value[0];
    const double_vec &time = y ? (lo ? bc.vbc_period_y0_time_in_yr : bc.vbc_period_y1_time_in_yr)
                               : (lo ? bc.vbc_period_x0_time_in_yr : bc.vbc_period_x1_time_in_yr);
    const double_vec &ratio = y ? (lo ? bc.vbc_period_y0_ratio : bc.vbc_period_y1_ratio)
                                : (lo ? bc.vbc_period_x0_ratio : bc.vbc_period_x1_ratio);
    return v * interp1(time, ratio, time_in_yr);
}
} // anonymous namespace

// True if side_wall_normal_vbc points into the domain (lower wall: > 0; upper wall: < 0).
bool side_wall_vbc_is_inflow(const BC &bc, int t, double time_in_yr)
{
    const double v = side_wall_normal_vbc(bc, t, time_in_yr);
    return (t % 2 == 0) ? (v > 0.) : (v < 0.);
}

// True if side_wall_normal_vbc points out of the domain (lower wall: < 0; upper wall: > 0).
bool side_wall_vbc_is_outflow(const BC &bc, int t, double time_in_yr)
{
    const double v = side_wall_normal_vbc(bc, t, time_in_yr);
    return (t % 2 == 0) ? (v < 0.) : (v > 0.);
}

namespace {
// What apply_vbcs does to one global component of a side wall or z side.
enum { FREE, CONSTANT, PROFILED };  // PROFILED: a 2D normal value times the depth profile

// A side wall or z side in global components, for one apply_vbcs call.
struct WallBC {
    int mode[NDIMS];
    double value[NDIMS];
    bool taper;
};

// The z sides of the node set the components they hold.
#pragma acc routine seq
void apply_z_sides(const WallBC *walls, uint flag, ArrayAccessor v)
{
    for (int b=iboundz0; b<=iboundz1; ++b) {
        if (!(flag & (1U << b))) continue;
        for (int d=0; d<NDIMS; ++d)
            if (walls[b].mode[d] != FREE) v[d] = walls[b].value[d];
    }
}

// What a side is to an n boundary of code 1 or 11 that meets it: a side that holds no
// velocity, or no longer does, lets it set its normal; a side holding only its normal (11
// only when normalized) projects onto the edge. Where a side at the node holds other
// components, the n boundary instead sets its normal while keeping what the sides there
// hold, and does nothing where its normal lies in that.
enum { MEETS_FREE, MEETS_NORMAL, MEETS_HORIZONTAL_NORMAL, MEETS_HELD };

// The part of the n-boundary normal n (its horizontal part, for 11 and 13) in r.
#pragma acc routine seq
void normal_part(ConstArrayAccessor n, bool horizontal, double r[NDIMS])
{
    for (int d=0; d<NDIMS; ++d) r[d] = (horizontal && d == NDIMS-1) ? 0 : n[d];
}

// Removes from r its parts along the first k vectors of the orthonormal basis, and puts r.r
// in left; returns whether more than rounding is left.
#pragma acc routine seq
bool remove_parts(const double basis[NDIMS][NDIMS], int k, double r[NDIMS], double &left)
{
    double before = 0;
    for (int d=0; d<NDIMS; ++d) before += r[d] * r[d];
    for (int j=0; j<k; ++j) {
        double rb = 0;
        for (int d=0; d<NDIMS; ++d) rb += r[d] * basis[j][d];
        for (int d=0; d<NDIMS; ++d) r[d] -= rb * basis[j][d];
    }
    left = 0;
    for (int d=0; d<NDIMS; ++d) left += r[d] * r[d];
    return left > 1e-12 * before;
}

// Sets v.n of n boundary ib (over the horizontal components, for 11) to target by the
// smallest change that keeps the global components in axes and, before ib at the node, the
// normals (horizontal normals, for 11) of the n boundaries; nothing where n lies in what they
// keep. The z sides of the node, which apply their components last, give them first.
#pragma acc routine seq
void set_normal_outside(const Variables &var, int ib, uint flag, int axes, const int *meets,
                        const bool *horizontal, const WallBC *walls, ArrayAccessor v,
                        double target)
{
    if (horizontal[ib]) axes |= 1 << (NDIMS-1);
    // An orthonormal basis of what is kept, the axes first.
    double basis[NDIMS][NDIMS], left;
    int k = 0;
    for (int d=0; d<NDIMS; ++d) {
        if (!(axes & (1 << d))) continue;
        for (int e=0; e<NDIMS; ++e) basis[k][e] = (e == d);
        ++k;
    }
    for (int ic=iboundn0; ic<ib && k<NDIMS; ++ic) {
        if (!(flag & (1U << ic)) || meets[ic] == MEETS_FREE) continue;
        double r[NDIMS];
        normal_part((*var.bnormals)[ic], horizontal[ic], r);
        if (!remove_parts(basis, k, r, left)) continue;  // already kept
        for (int d=0; d<NDIMS; ++d) basis[k][d] = r[d] / std::sqrt(left);
        ++k;
    }
    double n[NDIMS], p[NDIMS];
    normal_part((*var.bnormals)[ib], horizontal[ib], n);
    for (int d=0; d<NDIMS; ++d) p[d] = n[d];
    if (!remove_parts(basis, k, p, left)) return;  // p: the part of n free to change
    apply_z_sides(walls, flag, v);
    double vn = 0;
    for (int d=0; d<NDIMS; ++d) vn += v[d] * n[d];
    for (int d=0; d<NDIMS; ++d)
        if (!(axes & (1 << d))) v[d] += (target - vn) * p[d] / left;
}
} // anonymous namespace


// The velocity conditions apply_vbcs gives at the current time, as read from the input, so
// that a restart can tell whether its velocity boundary conditions are the run's at its
// checkpoint: which components the sides hold and at what values, a normal value that follows
// its period table at the value the table gives then, the frame of an n boundary, and in 2D
// the depth profile and bottom shear zone of such a normal value. Tractions and free components
// do not count, and a side holds nothing once released or if no node lies on it.
double_vec vbc_controls(const Param &param, const Variables &var)
{
    const BC &bc = param.bc;
    const double t_now = var.time / YEAR2SEC;
    uint present = 0;  // the sides some node lies on
    for (int n=0; n<var.nnode; ++n) present |= (*var.bcflag)[n];
    double_vec c;
    for (int i=0; i<nbdrytypes; ++i) {
        const SideBC &s = bc.side[i];
        const bool in_force = ((present >> i) & 1) && !(var.time > s.release_time) &&
                              !all_free(s);
        c.push_back(in_force);
        if (!in_force) continue;
        for (int d=0; d<NDIMS; ++d) {
            c.push_back(s.type[d] == 'v');
            if (s.type[d] != 'v') continue;
            c.push_back((d == 0 && s.normal_follows_period) ? side_wall_normal_vbc(bc, i, t_now)
                                                            : s.value[d]);
        }
        c.push_back(s.horizontal_normal);
        c.push_back(s.normal_follows_period);
#ifndef THREED
        if (!s.normal_follows_period) continue;
        auto add = [&c](const double_vec &table) {
            c.push_back(table.size());
            c.insert(c.end(), table.begin(), table.end());
        };
        add(i == iboundx0 ? var.vbc_vertical_div_x0 : var.vbc_vertical_div_x1);
        add(i == iboundx0 ? var.vbc_vertical_ratio_x0 : var.vbc_vertical_ratio_x1);
        const bool taper = s.taper && bc.bottom_shear_zone_thickness > 0;
        c.push_back(taper);
        if (taper) c.push_back(bc.bottom_shear_zone_thickness);
#endif
    }
    return c;
}


void apply_vbcs(const Param &param, const Variables &var, array_t &vel, bool homogeneous, bool normalized)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    // Each side's table (param.bc.side, in the side's frame) is turned back into global
    // components: the normal value of a side wall follows its period table, and a
    // homogeneous call holds every component at 0.
    const BC &bc = param.bc;
    const double t_now = var.time / YEAR2SEC;

    WallBC walls[nbdrytypes_hydro];
    bool any_z = false;
    for (int i=0; i<nbdrytypes_hydro; ++i) {
        const SideBC &s = bc.side[i];
        WallBC &w = walls[i];
        w.taper = s.taper;
        for (int d=0; d<NDIMS; ++d) {
            w.mode[d] = FREE;
            w.value[d] = 0;
        }
        if (var.time > s.release_time) continue;
        for (int d=0; d<NDIMS; ++d) {
            if (s.type[d] != 'v') continue;
            const int axis = frame_axis(i, d);
            const double v = frame_sign(i, d) * s.value[d];
            w.mode[axis] = CONSTANT;
            if (d == 0 && s.normal_follows_period) {
                const double applied = homogeneous ? 0.0 : side_wall_normal_vbc(bc, i, t_now);
                w.value[axis] = applied;
#ifndef THREED
                w.mode[axis] = PROFILED;
#endif
            }
            else
                w.value[axis] = homogeneous ? 0.0 : v;
            if (i == iboundz0 || i == iboundz1) any_z = true;
        }
    }

    int meets[nbdrytypes], held_axes[nbdrytypes];
    bool held[nbdrytypes], horizontal[nbdrytypes], full[nbdrytypes];
    double n_value[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        const SideBC &s = bc.side[i];
        held[i] = (s.type[0] == 'v');
        horizontal[i] = s.horizontal_normal;
        full[i] = true;
        for (int d=0; d<NDIMS; ++d)
            if (s.type[d] != 'v') full[i] = false;
        n_value[i] = homogeneous ? 0.0 : s.value[0];
        const bool released = var.time > s.release_time;
        meets[i] = (released || all_free(s)) ? MEETS_FREE : !normal_only(s) ? MEETS_HELD
                 : s.horizontal_normal ? MEETS_HORIZONTAL_NORMAL : MEETS_NORMAL;
        // The global axes a side wall or z side holds, and every axis for an n boundary that
        // holds every component; one holding its normal keeps that instead.
        held_axes[i] = (i >= iboundn0 && meets[i] == MEETS_HELD) ? (1 << NDIMS) - 1 : 0;
        for (int d=0; d<NDIMS && i < iboundn0; ++d)
            if (s.type[d] == 'v') held_axes[i] |= 1 << frame_axis(i, d);
    }

#ifndef THREED
    // find max and min coordinate for BOUNDX0
    // Side-wall extents as device reductions (a host scan would migrate coord/bcflag every step).
    // max/min are exact, so the result is order-independent; a wall with no node gets extent 0.
    const double big = std::numeric_limits<double>::max();
    double BOUNDX0_max = -big, BOUNDX0_min = big;
    double BOUNDX1_max = -big, BOUNDX1_min = big;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var) \
        reduction(max:BOUNDX0_max, BOUNDX1_max) reduction(min:BOUNDX0_min, BOUNDX1_min)
#endif
    #pragma acc parallel loop gang vector reduction(max:BOUNDX0_max, BOUNDX1_max) reduction(min:BOUNDX0_min, BOUNDX1_min) async
    for (int i=0; i<var.nnode; ++i) {
        if (! is_on_boundary(*var.bcflag, i)) continue;
        uint flag = (*var.bcflag)[i];
        if (flag & BOUNDX0) {
            const double z = (*var.coord)[i][1];
            if (z > BOUNDX0_max) BOUNDX0_max = z;
            if (z < BOUNDX0_min) BOUNDX0_min = z;
        } else if (flag & BOUNDX1) {
            const double z = (*var.coord)[i][1];
            if (z > BOUNDX1_max) BOUNDX1_max = z;
            if (z < BOUNDX1_min) BOUNDX1_min = z;
        }
    }
    #pragma acc wait
    if (BOUNDX0_max < BOUNDX0_min) { BOUNDX0_max = 0.; BOUNDX0_min = 0.; }
    if (BOUNDX1_max < BOUNDX1_min) { BOUNDX1_max = 0.; BOUNDX1_min = 0.; }
    double BOUNDX0_width = BOUNDX0_max - BOUNDX0_min;
    double BOUNDX1_width = BOUNDX1_max - BOUNDX1_min;

    const double_vec& vbc_vertical_ratios_x0 = var.vbc_vertical_ratio_x0;
    const double_vec& vbc_vertical_ratios_x1 = var.vbc_vertical_ratio_x1;

    double_vec vbc_vertical_divisions_x0(4,0.);
    double_vec vbc_vertical_divisions_x1(4,0.);
    for (int i=0;i<4;i++) {
        vbc_vertical_divisions_x0[i] = - (BOUNDX0_max - var.vbc_vertical_div_x0[i] * BOUNDX0_width);
        vbc_vertical_divisions_x1[i] = - (BOUNDX1_max - var.vbc_vertical_div_x1[i] * BOUNDX1_width);
    }
    const double thickness = bc.bottom_shear_zone_thickness;

    double zmin = 0;
    #pragma acc enter data copyin(zmin) async
#ifndef ACC
    #pragma omp parallel for default(none) shared(var) reduction(min:zmin)
#endif
    #pragma acc parallel loop gang vector async reduction(min:zmin)
    for (int k=0; k<var.nnode; ++k) {
        ConstArrayAccessor tmpx = (*var.coord)[k];
        if ( tmpx[NDIMS-1] < zmin )
            zmin = tmpx[NDIMS-1];
    }
#endif
#ifndef ACC
#ifdef THREED
    #pragma omp parallel for default(none) \
        shared(var, vel, normalized, walls, any_z, meets, held_axes, held, horizontal, full, n_value)
#else
    #pragma omp parallel for default(none) \
        shared(var, vel, normalized, walls, any_z, meets, held_axes, held, horizontal, full, n_value, zmin, thickness, \
        vbc_vertical_divisions_x0, vbc_vertical_divisions_x1, \
        vbc_vertical_ratios_x0, vbc_vertical_ratios_x1)
#endif
#endif
#ifdef THREED
    #pragma acc parallel loop gang vector async \
        copyin(walls[0:nbdrytypes_hydro], meets[0:nbdrytypes], held_axes[0:nbdrytypes], \
               held[0:nbdrytypes], horizontal[0:nbdrytypes], full[0:nbdrytypes], n_value[0:nbdrytypes])
#else
    #pragma acc parallel loop gang vector async present(zmin) \
        copyin(walls[0:nbdrytypes_hydro], meets[0:nbdrytypes], held_axes[0:nbdrytypes], \
               held[0:nbdrytypes], horizontal[0:nbdrytypes], full[0:nbdrytypes], n_value[0:nbdrytypes])
#endif
    for (int i=0; i<var.nnode; ++i) {

        // fast path: skip nodes not on boundary
        if (! is_on_boundary(*var.bcflag, i)) continue;

        uint flag = (*var.bcflag)[i];
        ArrayAccessor v = vel[i];

        //
        // Side walls, in the order x0, x1, y0, y1
        //
        for (int b=iboundx0; b<=iboundy1; ++b) {
            if (!(flag & (1U << b))) continue;
            const WallBC &w = walls[b];
            for (int d=0; d<NDIMS; ++d) {
                if (w.mode[d] == FREE) continue;
                v[d] = w.value[d];
#ifndef THREED
                if (w.mode[d] == PROFILED) {
                    ConstArrayAccessor x = (*var.coord)[i];
                    v[d] = w.value[d] * interp1(b == iboundx0 ? vbc_vertical_divisions_x0 : vbc_vertical_divisions_x1,
                                                b == iboundx0 ? vbc_vertical_ratios_x0 : vbc_vertical_ratios_x1, -x[1]);
                    if (w.taper && thickness > 0.) {
                        double dz = x[NDIMS-1] - zmin;
                        if (dz < thickness)
                            v[d] = v[d] * dz / thickness;
                    }
                }
#endif
            }
        }

        //
        // N
        //
        for (int ib=iboundn0; ib<=iboundn3; ib++) {
            ConstArrayAccessor n = (*var.bnormals)[ib]; // unit normal vector

            if ((flag & (1 << ib)) && held[ib]) {
                double fac = 0;
                int axes = 0;  // the global axes held by the sides at this node
                bool held_other = false;  // one of them holds more than its normal
                for (int ic=iboundx0; ic<ib; ic++) {
                    if (!(flag & (1 << ic)) || meets[ic] == MEETS_FREE) continue;
                    axes |= held_axes[ic];
                    held_other |= meets[ic] == MEETS_HELD;
                }
                if (!full[ib] && !horizontal[ib]) {  // the normal held, as code 1
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS; d++)
                            vn += v[d] * n[d];  // normal velocity

                        for (int d=0; d<NDIMS; d++)
                            v[d] += (n_value[ib] - vn) * n[d];  // setting normal velocity
                    }
                    else if (held_other) {  // the sides at this node together
                        set_normal_outside(var, ib, flag, axes, meets, horizontal, walls, v, n_value[ib]);
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (meets[ic] == MEETS_FREE) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        vn += v[d] * n[d];  // normal velocity

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] += (n_value[ib] - vn) * n[d];  // setting normal velocity
                                }
                                else if (meets[ic] == MEETS_NORMAL || (normalized && meets[ic] == MEETS_HORIZONTAL_NORMAL)) {
                                    project_boundary_intersection(var, ic, ib, n, v, normalized, n_value[ib],
                                                                  horizontal[ic], horizontal[ib]);
                                }
                            }
                        }
                    }
                }
                else if (!horizontal[ib]) {  // every component held, as code 3
                    for (int d=0; d<NDIMS; d++)
                        v[d] = n_value[ib] * n[d];  // v must be normal to n
                }
                else if (!full[ib]) {  // the horizontal normal held, as code 11
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS-1; d++)
                            vn += v[d] * n[d];  // normal velocity

                        for (int d=0; d<NDIMS-1; d++)
                            v[d] += (n_value[ib] * (normalized ? 1/fac : fac) - vn) * n[d] / (normalized ? (1 - n[NDIMS-1]*n[NDIMS-1]) : 1.0);  // setting normal velocity
                    }
                    else if (held_other) {  // the sides at this node together
                        set_normal_outside(var, ib, flag, axes, meets, horizontal, walls, v, n_value[ib] / fac);
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (meets[ic] == MEETS_FREE) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS-1; d++)
                                        vn += v[d] * n[d];  // normal velocity

                                    for (int d=0; d<NDIMS-1; d++)
                                        v[d] += (n_value[ib] * (normalized ? 1/fac : fac) - vn) * n[d] / (normalized ? (1 - n[NDIMS-1]*n[NDIMS-1]) : 1.0);  // setting normal velocity
                                }
                                else if (meets[ic] == MEETS_NORMAL || (normalized && meets[ic] == MEETS_HORIZONTAL_NORMAL)) {
                                    project_boundary_intersection(var, ic, ib, n, v, normalized, n_value[ib],
                                                                  horizontal[ic], horizontal[ib]);
                                }
                            }
                        }
                    }
                }
                else {  // the horizontal velocity along the horizontal normal, as code 13
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    for (int d=0; d<NDIMS-1; d++)
                        v[d] = n_value[ib] * fac * n[d];
                    v[NDIMS-1] = 0;
                }
            }
        }

        //
        // Z, must be dealt last
        //

        // fast path: vz is usually free in the models
        if (!any_z) continue;

        apply_z_sides(walls, flag, v);
    }
#ifndef THREED
    #pragma acc exit data copyout(zmin) async
#endif
    #pragma acc wait
#ifdef NPROF
    nvtxRangePop();
#endif
}

void apply_stress_bcs(const Param& param, const Variables& var, array_t& force,
                      double displacement_dt)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    if (param.control.gravity == 0) return;
    //
    // Gravity-induced (hydrostatic and lithostatic) stress BCs
    //

#ifndef ACC
    #pragma omp parallel for default(none) shared(var)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<var.nelem; ++e)
        (*var.etmp_int)[e] = -1;

    #pragma acc wait

    for (int i=0; i<nbdrytypes; i++) {
        // a free normal on a side wall, or on z0 with a Winkler foundation and z1 with water loading
        if (!param.bc.side[i].support) continue;

        // Plain host reads of two vector sizes: an `acc kernels` block here becomes a serial device
        // kernel plus a synchronization per boundary per step (nvc++).
        int bound = static_cast<int>(var.bfacets[i]->size());
        int nbdry_nodes = static_cast<int>(var.bnodes[i]->size());

#ifndef ACC
        #pragma omp parallel default(none) \
            shared(param, var, force, i, NODE_OF_FACET, bound, nbdry_nodes, displacement_dt)
#endif
        {
            // loops over all bdry facets
#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop gang vector async
            for (int n=0; n<bound; ++n) {
                // this facet belongs to element e
                int e = (*var.bfacets[i])[n].first;
                // this facet is the f-th facet of e
                int f = (*var.bfacets[i])[n].second;
                ConstConnAccessor conn = (*var.connectivity)[e];

                // the outward-normal vector
                double normal[NDIMS];
                // the z-coordinate of the facet center
                double zcenter;

                int idx[NODES_PER_FACET];
                for (int j=0; j<NODES_PER_FACET; ++j)
                    idx[j] = conn[NODE_OF_FACET[f][j]];

                ConstArrayIndirectAccessor facet_coord = var.coord->view_const(idx);

                normal_vector_of_facet(facet_coord, normal, zcenter);

                // Incremental PT uses reference facets but evaluates the support
                // pressure at the candidate height. Never move the actual mesh
                // or advance surface processes while testing a mechanical trial.
                if (displacement_dt != 0) {
                    double vz = 0;
                    for (int j=0; j<NODES_PER_FACET; ++j)
                        vz += (*var.vel)[idx[j]][NDIMS-1];
                    zcenter += displacement_dt * vz / NODES_PER_FACET;
                }

                double p;
                if (i==iboundz0 && param.bc.has_winkler_foundation) {
                    double rho_effective = var.mat->rho(e);  // Base density of the solid material

                    // If hydraulic diffusion is active, modify the density to account for porosity and fluid content
                    if (param.control.has_hydraulic_diffusion) {
                        rho_effective = (var.mat->rho(e) * (1 - var.mat->phi(e)) + 1000.0 * var.mat->phi(e));
                        // 1000.0 is the fluid density (e.g., water, in kg/m³)
                        // Winkler foundation for the bottom boundary with adjusted effective density
                        p = var.compensation_pressure - 
                            (rho_effective + param.bc.winkler_delta_rho) *  // Effective density with hydraulic diffusion
                            param.control.gravity * (zcenter + param.mesh.zlength);  // Adjust for depth from base
                    }
                    else
                    {
                        // Winkler foundation for the bottom boundary with adjusted effective density
                        p = var.compensation_pressure - 
                            (rho_effective + param.bc.winkler_delta_rho) *  // Effective density with hydraulic diffusion
                            param.control.gravity * (zcenter + param.mesh.zlength);  // Adjust for depth from base
                    }
                }
                else if (i==iboundz1 && param.bc.has_water_loading) {
                    // hydrostatic water loading for the surface boundary
                    p = 0;
                    if (zcenter < param.control.surf_base_level) {
                        // below sea level
                        p = param.bc.sea_water_density * param.control.gravity * (param.control.surf_base_level - zcenter);
                    }
                }
                else {
                    // sidewalls: lithostatic support pressure (Archimedes) from the exterior column.
                    // Floor at zero: for a wall facet ABOVE the datum (zcenter > 0, wall-edge
                    // topography) ref_pressure_option 0/4 extrapolate linearly NEGATIVE, which would
                    // flip the traction into an unphysical OUTWARD suction that actively pulls the
                    // topographic strip out of the domain. The true exterior there is air/vacuum
                    // (zero traction), so clamp. Options 1/2/3 already return 0 above the datum, so
                    // this only changes options 0/4; below the datum (the usual case) p > 0 and the
                    // clamp is inert -> normal runs are bit-identical.
                    p = ref_pressure(param, zcenter);
                    if (p < 0.0) p = 0.0;
                }

                (*var.etmp_int)[e] = n;
                ElemCacheAccessor tmp = (*var.tmp_result)[n];
                // lithostatc support - Archimed force (normal to the surface)
                for (int j=0; j<NODES_PER_FACET; ++j) {
                    int nn = conn[NODE_OF_FACET[f][j]];
                    for (int d=0; d<NDIMS; ++d)
                        tmp[j*NDIMS + d] = p * normal[d] / NODES_PER_FACET;
                }
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop gang vector async
            for (int j=0; j<nbdry_nodes; ++j) {
                const int n = (*var.bnodes[i])[j];
                const int npatch = var.support.size(n);
                const int* patch = var.support.patch(n);
                for (int k=0; k<npatch; ++k) {
                    int e = patch[k];
                    int ibound = (*var.etmp_int)[e];
                    if (ibound < 0) continue;  // not a boundary element

                    int f = (*var.bfacets[i])[ibound].second;  // facet index
                    ConstConnAccessor conn = (*var.connectivity)[e];
                    for (int l=0; l<NODES_PER_FACET; ++l) {
                        if (n == conn[NODE_OF_FACET[f][l]]) {
                            for (int d=0; d<NDIMS; ++d)
                                force[n][d] -= (*var.tmp_result)[ibound][l*NDIMS + d];  // subtract the force from the facet
                            break;  // found the node in the facet
                        }
                    }
                }
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop gang vector async
            for (int n=0; n<bound; ++n) {
                const int e = (*var.bfacets[i])[n].first;
                (*var.etmp_int)[e] = -1;
            }
        }
    }

    #pragma acc wait

    if (param.bc.has_elastic_foundation) {
        /* A restoration force on the bottom nodes proportional to total vertical displacement */
        for (auto i=var.bnodes[iboundz0]->begin(); i<var.bnodes[iboundz0]->end(); ++i) {
            int n = *i;
            force[n][NDIMS-1] -= param.bc.elastic_foundation_constant * ((*var.coord)[n][NDIMS-1] - (*var.coord0)[n][NDIMS-1]);
        }
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

namespace {
// The frame axes of side i in global components: a side wall's or z side's from side_frame,
// an n boundary's from its normal n as t_h = up x n / |up x n| (+x for a level n) and
// t_v = n x t_h; in 2D the in-plane t is t_v, or +x for a level n.
void side_axes(const Variables &var, int i, double axes[NDIMS][NDIMS])
{
    for (int c=0; c<NDIMS; ++c)
        for (int d=0; d<NDIMS; ++d) axes[c][d] = 0;
    if (i < iboundn0) {
        for (int c=0; c<NDIMS; ++c)
            axes[c][frame_axis(i, c)] = frame_sign(i, c);
        return;
    }
    ConstArrayAccessor n = (*var.bnormals)[i];
    for (int d=0; d<NDIMS; ++d) axes[0][d] = n[d];
#ifdef THREED
    const double h = std::sqrt(n[0]*n[0] + n[1]*n[1]);
    if (h < 1e-12) {
        axes[1][0] = 1;
    } else {
        axes[1][0] = -n[1] / h;
        axes[1][1] = n[0] / h;
    }
    axes[2][0] = n[1]*axes[1][2] - n[2]*axes[1][1];
    axes[2][1] = n[2]*axes[1][0] - n[0]*axes[1][2];
    axes[2][2] = n[0]*axes[1][1] - n[1]*axes[1][0];
#else
    if (std::fabs(n[0]) < 1e-12) {
        axes[1][0] = 1;
    } else {
        const double s = (n[0] < 0) ? -1. : 1.;
        axes[1][0] = -s * n[1];
        axes[1][1] = s * n[0];
    }
#endif
}
} // anonymous namespace


void apply_stress_bcs_neumann(const Param& param, const Variables& var, array_t& force)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    #pragma acc wait

    // stress_bc_* (codes only) loads each facet with S n_k on component k, n being the
    // facet's area-weighted normal. A traction t (letters only) loads it with T times
    // the facet area along the frame axis. The facet's nodes share the load.
    for (int i = 0; i < nbdrytypes; ++i) {
    const SideBC &side = param.bc.side[i];
    bool any = false;
    for (int c = 0; c < NDIMS; ++c) any = any || side.type[c] == 't';
    if ((!any && side.stress_axis < 0) || !var.bfacets[i] || var.bfacets[i]->empty()) continue;
    double stress[NDIMS] = {0.0}, traction[NDIMS] = {0.0};
    if (side.stress_axis >= 0) stress[side.stress_axis] = side.stress_value;
    if (any) {
        double axes[NDIMS][NDIMS];
        side_axes(var, i, axes);
        for (int c = 0; c < NDIMS; ++c)
            if (side.type[c] == 't')
                for (int d = 0; d < NDIMS; ++d) traction[d] += side.value[c] * axes[c][d];
    }

    const auto& bdry = *(var.bfacets[i]);
    const int bound = static_cast<int>(bdry.size());

    // Loop over all boundary facets
    for (int n = 0; n < bound; ++n) {
        int e = bdry[n].first;      // This facet belongs to element e
        int f = bdry[n].second;     // This facet is the f-th facet of element e
        ConstConnAccessor conn = (*var.connectivity)[e];  // Element connectivity

        // Outward-normal vector, whose length is the facet area, and z-center for the facet
        double normal[NDIMS] = {0.0}; 
        double zcenter = 0.0;

        int idx[NODES_PER_FACET];

        for (int j=0; j<NODES_PER_FACET; ++j)
            idx[j] = conn[NODE_OF_FACET[f][j]];

        ConstArrayIndirectAccessor facet_coord = var.coord->view_const(idx);

        normal_vector_of_facet(facet_coord, normal, zcenter);
        double area = 0;
        for (int d = 0; d < NDIMS; ++d) area += normal[d] * normal[d];
        area = std::sqrt(area);
        double load[NDIMS];
        for (int d = 0; d < NDIMS; ++d)
            load[d] = (side.stress_axis >= 0) ? stress[d] * normal[d] : traction[d] * area;

        // Apply the load to the boundary facet nodes
        for (int j = 0; j < NODES_PER_FACET; ++j) {
            int nn = conn[NODE_OF_FACET[f][j]];  // Get global node index for this facet
            ArrayAccessor force_node = force[nn];      // Pointer to the force array for this node
            for (int d = 0; d < NDIMS; ++d)
                force_node[d] += load[d] / NODES_PER_FACET;
        }
    }
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

namespace {

    void simple_diffusion(const Variables& var)
    {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif
        /* Diffusing surface topography to simulate the effect of erosion and
         * sedimentation.
         */

#ifndef ACC
        #pragma omp parallel default(none) shared(var,NODE_OF_FACET)
#endif
        {
            const array_t& coord = *var.coord;
            const int_vec& top_nodes = *var.surfinfo.top_nodes;

            const auto& top = *(var.bfacets[iboundz1]);

            const int ntop = var.surfinfo.ntop;

            double_vec& total_dx = *var.surfinfo.total_dx;
            double_vec& total_slope = *var.surfinfo.total_slope;
            double_vec& dh = *var.surfinfo.dh;

            // loops over all top facets
#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop gang vector async
            for (int i=0; i<var.nnode; i++) {
                total_dx[i] = 0.;
                total_slope[i] = 0.;
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop gang vector async
            for (int i=0; i<var.surfinfo.etop; ++i) {
#ifdef THREED
                // this facet belongs to element e
                int e = top[i].first;
                // this facet is the f-th facet of e
                int f = top[i].second;

                ConstConnAccessor conn = (*var.connectivity)[e];
                int idx[NODES_PER_FACET];
                for (int j=0; j<NODES_PER_FACET; ++j)
                    idx[j] = conn[NODE_OF_FACET[f][j]];

                ConstArrayIndirectAccessor coord_f = var.coord->view_const(idx);

                double projected_area;
                {
                    // normal vector of this facet
                    double normal[NDIMS];

                    // two vectors n0-n1 and n0-n2
                    // n is the cross product of these two vectors
                    // the length of n is 2 * triangle area
                    double x01, y01, z01, x02, y02, z02;
                    x01 = coord_f[1][0] - coord_f[0][0];
                    y01 = coord_f[1][1] - coord_f[0][1];
                    z01 = coord_f[1][2] - coord_f[0][2];
                    x02 = coord_f[2][0] - coord_f[0][0];
                    y02 = coord_f[2][1] - coord_f[0][1];
                    z02 = coord_f[2][2] - coord_f[0][2];

                    normal[0] = y01*z02 - z01*y02;
                    normal[1] = z01*x02 - x01*z02;
                    normal[2] = x01*y02 - y01*x02;

                    /* the area of this facet is:
                    *   0.5 * std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
                    *
                    * theta is the angle between this facet and horizontal
                    *   tan_theta = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]) / normal[2]
                    *   cos_theta = normal[2] / (2 * area)
                    *
                    * the area projected on the horizontal plane is:
                    *   projected_area = area * cos_theta = 0.5 * normal[2]
                    */
                    projected_area = 0.5 * normal[2];
                }

                (*var.etmp)[i] = projected_area;

                double shp2dx[NODES_PER_FACET], shp2dy[NODES_PER_FACET];
                double iv = 1 / (2 * projected_area);
                shp2dx[0] = iv * (coord_f[1][1] - coord_f[2][1]);
                shp2dx[1] = iv * (coord_f[2][1] - coord_f[0][1]);
                shp2dx[2] = iv * (coord_f[0][1] - coord_f[1][1]);
                shp2dy[0] = iv * (coord_f[2][0] - coord_f[1][0]);
                shp2dy[1] = iv * (coord_f[0][0] - coord_f[2][0]);
                shp2dy[2] = iv * (coord_f[1][0] - coord_f[0][0]);

                double D[NODES_PER_FACET][NODES_PER_FACET];
                for (int j=0; j<NODES_PER_FACET; j++) {
                    for (int k=0; k<NODES_PER_FACET; k++) {
                        D[j][k] = (shp2dx[j] * shp2dx[k] +
                                shp2dy[j] * shp2dy[k]);
                    }
                }

                for (int j=0; j<NODES_PER_FACET; j++) {
                    double slope = 0;
                    for (int k=0; k<NODES_PER_FACET; k++)
                        slope += D[j][k] * coord_f[k][2];
                    (*var.tmp_result)[i][j] = slope * projected_area;
                }
#else
                /* The 1D diffusion operation is implemented ad hoc, not using FEM
                * formulation (e.g. computing shape function derivation on the edges).
                */
                int n0 = top_nodes[i];
                int n1 = top_nodes[i+1];

                double dx = std::fabs(coord[n1][0] - coord[n0][0]);
                (*var.etmp)[i] = dx;

                (*var.tmp_result)[i][0] = -(coord[n1][1] - coord[n0][1]) / dx;
                (*var.tmp_result)[i][1] = (coord[n1][1] - coord[n0][1]) / dx;
#endif
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop gang vector async
            for (int i=0; i<ntop; ++i) {
                int n = top_nodes[i];
#ifdef THREED
                const int nsurf_sup = var.surfinfo.support_surf.size(i);
                const int* surf_sup = var.surfinfo.support_surf.patch(i);

                for (int j=0; j<nsurf_sup; ++j) {
                    int k = surf_sup[j];
                    total_dx[n] += (*var.etmp)[k];

                    // this facet belongs to element e
                    int e = top[k].first;
                    // this facet is the f-th facet of e
                    int f = top[k].second;

                    ConstConnAccessor conn = (*var.connectivity)[e];
                    int n0 = conn[NODE_OF_FACET[f][0]];
                    int n1 = conn[NODE_OF_FACET[f][1]];
                    int n2 = conn[NODE_OF_FACET[f][2]];

                    const int nn[NODES_PER_FACET] = {n0, n1, n2};

                    for (int m=0; m<NODES_PER_FACET; ++m) {
                        if (nn[m] == n) {
                            total_slope[n] += (*var.tmp_result)[k][m];
                            break;
                        }
                    }
                }
#else
                if (i == 0) {
                    total_dx[n] = (*var.etmp)[i];
                    total_slope[n] = (*var.tmp_result)[i][0];
                } else if (i == ntop-1) {
                    total_dx[n] = (*var.etmp)[i-1];
                    total_slope[n] = (*var.tmp_result)[i-1][1];
                } else {
                    total_dx[n] = (*var.etmp)[i-1] + (*var.etmp)[i];
                    total_slope[n] = (*var.tmp_result)[i-1][1] + (*var.tmp_result)[i][0];
                }
#endif
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop gang vector async
            for (int i=0; i<ntop; ++i) {
                // we don't treat edge nodes specially, i.e. reflecting bc is used for erosion.
                int n = top_nodes[i];
                double conv =  var.surfinfo.surf_diff * var.dt * total_slope[n] / total_dx[n];
    #ifdef THREED
                dh[i] -= conv;
    #else
                if ( coord[n][1] >  var.surfinfo.base_level && conv > 0.) {
                    dh[i] -= var.surfinfo.diff_ratio_terrig * conv;
                } else if ( coord[n][1] <= var.surfinfo.base_level && conv < 0. ) {
                    dh[i] -= var.surfinfo.diff_ratio_marine * conv;
                } else {
                    dh[i] -= conv;
                }
    #endif
            }
        }
#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }


    void custom_surface_processes(const Variables& var, array_t& coord)
    {
        #pragma acc wait // here is not ACC parallelized yet

        const int top_bdry = iboundz1;
        const int_vec& top_nodes = *var.bnodes[top_bdry];
        const std::size_t ntop = top_nodes.size();

        // loops over all top nodes
        for (std::size_t i=0; i<ntop; ++i) {
            int n = top_nodes[i];

            // coordinate of this node
            double x = coord[n][0];
            double z = coord[n][NDIMS-1];

            // compute topography changes here ...
            // positive dh means sedimentation, negative dh means erosion
            double dh = 0;
            {
                // dh == ....
            }

            coord[n][NDIMS-1] += dh;
        }
    }

#ifdef HAS_GOSPL_CPP_INTERFACE
    void use_gospl(const Param& param, const Variables& var)
    {
        if (!var.gospl_driver || !var.gospl_driver->is_initialized()) {
            std::cerr << "Error: GoSPL driver is not initialized" << std::endl;
            return;
        }

        // Accumulate time and steps since last coupling
        const double SEC_PER_YR = 3.1536e7;
        double dt_years = var.dt / SEC_PER_YR;
        var.gospl_driver->accumulated_dt += dt_years;
        var.gospl_driver->step_counter++;

        // Check whether it's time to couple
        bool do_couple;
        if (var.gospl_driver->coupling_by_time)
            do_couple = (var.gospl_driver->accumulated_dt >= var.gospl_driver->coupling_interval_in_yr);
        else
            do_couple = (var.gospl_driver->step_counter >= var.gospl_driver->coupling_frequency);

        if (!do_couple)
            return;  // GoSPL idle this step

        // ===== COUPLING STEP =====
        double total_dt = var.gospl_driver->accumulated_dt;
        var.gospl_driver->accumulated_dt = 0.0;
        var.gospl_driver->step_counter   = 0;

        const int_vec& top_nodes = *var.bnodes[iboundz1];
        const std::size_t ntop = top_nodes.size();

        // Only the first event is logged in full, as a check on the setup; later
        // ones go into the summary printed with DES's progress line.
        const bool log_event = (var.gospl_driver->total_events == 0);

        if (log_event)
            std::cout << "GoSPL coupling: dt=" << total_dt << " yr, "
                      << ntop << " surface nodes, vel_coupling="
                      << (var.gospl_driver->velocity_coupling ? "on" : "off")
                      << std::endl;

        // Build coordinate array for current DES surface nodes
        std::vector<double> coords(ntop * 3);
        for (std::size_t i = 0; i < ntop; ++i) {
            int n = top_nodes[i];
            for (int d = 0; d < NDIMS; d++)
                coords[i * 3 + d] = (*var.coord)[n][d];
        }

        // Static storage for time-averaged velocity reference.
        // Stored as static (not in GoSPLDriver) to avoid any allocation in
        // GoSPLDriver's memory space, which is adjacent to PETSc-managed memory
        // that can interfere with the C++ heap.
        static std::vector<double> surface_coords_start;
        static bool has_surface_start = false;
        // var.nremesh when surface_coords_start was recorded. A remesh renumbers the nodes
        // (renumbering_mesh), so top_nodes[i] may name a different node afterwards even
        // when the surface node count is unchanged.
        static int nremesh_at_start = -1;

        // Re-initialize GoSPL topography on the first coupling event only.
        // After that, GoSPL owns its topography (ASPECT-FastScape scheme).
        if (var.gospl_driver->needs_elevation_reset) {
            std::vector<double> des_elev(ntop);
            for (std::size_t i = 0; i < ntop; ++i)
                des_elev[i] = (*var.coord)[top_nodes[i]][NDIMS-1];
            int r = var.gospl_driver->apply_elevation_data(
                coords.data(), des_elev.data(), ntop, 3, 1.0);
            if (r != 0)
                std::cerr << "Warning: Failed to reset GoSPL elevation" << std::endl;
            var.gospl_driver->needs_elevation_reset = false;
            has_surface_start = false;  // force instantaneous vel this event
        }

        // Pass all three DES surface velocity components to GoSPL.
        // Use time-averaged velocities (Δcoord/Δt over the coupling interval) when
        // available to filter out inertial (quasi-dynamic) oscillations in the
        // instantaneous DES velocity.  Fall back to instantaneous on the first
        // coupling event or after remeshing.
        if (var.gospl_driver->velocity_coupling) {
            std::vector<double> vx_yr(ntop), vy_yr(ntop), vz_yr(ntop);
            bool use_time_avg = has_surface_start && nremesh_at_start == var.nremesh;
            for (std::size_t i = 0; i < ntop; ++i) {
                int n = top_nodes[i];
                if (use_time_avg) {
                    vx_yr[i] = (coords[i*3+0]       - surface_coords_start[i*3+0])       / total_dt;
#ifdef THREED
                    vy_yr[i] = (coords[i*3+1]       - surface_coords_start[i*3+1])       / total_dt;
#else
                    vy_yr[i] = 0.0;
#endif
                    vz_yr[i] = (coords[i*3+NDIMS-1] - surface_coords_start[i*3+NDIMS-1]) / total_dt;
                } else {
                    vx_yr[i] = (*var.vel)[n][0]       * SEC_PER_YR;
#ifdef THREED
                    vy_yr[i] = (*var.vel)[n][1]        * SEC_PER_YR;
#else
                    vy_yr[i] = 0.0;
#endif
                    vz_yr[i] = (*var.vel)[n][NDIMS-1]  * SEC_PER_YR;
                }
            }
            int r = var.gospl_driver->set_surface_velocity(
                coords.data(), vx_yr.data(), vy_yr.data(), vz_yr.data(), ntop, 3, 1.0);
            if (r != 0)
                std::cerr << "Warning: Failed to set surface velocity in GoSPL" << std::endl;

            // Velocity transfer diagnostics
            double vx_min = *std::min_element(vx_yr.begin(), vx_yr.end());
            double vx_max = *std::max_element(vx_yr.begin(), vx_yr.end());
            double vy_min = *std::min_element(vy_yr.begin(), vy_yr.end());
            double vy_max = *std::max_element(vy_yr.begin(), vy_yr.end());
            double vz_min = *std::min_element(vz_yr.begin(), vz_yr.end());
            double vz_max = *std::max_element(vz_yr.begin(), vz_yr.end());
            const double vmin[3] = {vx_min, vy_min, vz_min};
            const double vmax[3] = {vx_max, vy_max, vz_max};
            var.gospl_driver->record_velocity(vmin, vmax, !use_time_avg);
            if (log_event)
                std::cout << "  vel DES->GoSPL (" << (use_time_avg ? "time-avg" : "instantaneous") << ", m/yr):"
                          << " vx=[" << vx_min << ", " << vx_max << "]"
                          << " vy=[" << vy_min << ", " << vy_max << "]"
                          << " vz=[" << vz_min << ", " << vz_max << "]"
                          << std::endl;
            else if (!use_time_avg) {
                // The fallback after the first event means a remesh. "same count" marks
                // the case a node-count check cannot see.
                const std::size_t nstart = surface_coords_start.size() / 3;
                std::cout << "GoSPL: instantaneous velocity at t=" << var.time / YEAR2SEC
                          << " yr: remeshed since the last coupling (surface nodes "
                          << nstart << " -> " << ntop
                          << (nstart == ntop ? ", same count" : "") << ")" << std::endl;
            }
        }

        // Run GoSPL and get net dh in one call.
        // Internally: horiz advection (if vel set) + uplift + erosion/deposition,
        // differenced on the native GoSPL mesh, IDW-interpolated to DES query points.
        std::vector<double> net_erosion(ntop, 0.0);
        int run_result = var.gospl_driver->run_and_get_erosion(
            total_dt, coords.data(), ntop, net_erosion.data(), 3, 1.0);
        if (run_result != 0)
            // A failed call can leave GoSPL's surface partly advanced, so the
            // coupling cannot continue from it.
            die(EXIT_RUNTIME_GOSPL, "GoSPL run_and_get_erosion failed");

        // Apply dh to all DES surface nodes.
        // net_erosion contains only erosion+diffusion (uplift stripped in Python);
        // tectonic uplift is already present in DES coords from the mechanical solver.
        // The GoSPL mesh is padded beyond the DES domain (gospl_mesh_padding),
        // so all DES nodes query interior GoSPL nodes, free of boundary artifacts.
        double min_dh = 0, max_dh = 0, sum_dh = 0;
        for (std::size_t i = 0; i < ntop; ++i) {
            int n = top_nodes[i];
            double dh = net_erosion[i];
            min_dh = std::min(min_dh, dh);
            max_dh = std::max(max_dh, dh);
            sum_dh += dh;
            (*var.coord)[n][NDIMS-1] += dh;
            coords[i*3+NDIMS-1]      += dh;  // keep coords in sync for recording below
        }

        // Record end-of-coupling surface coords for next interval's time-averaged velocity.
        // Done after GoSPL returns; surface_coords_start is a static vector allocated in
        // DES's own memory space (not adjacent to PETSc memory), so assignment is safe.
        surface_coords_start = coords;
        nremesh_at_start = var.nremesh;
        has_surface_start = true;

        var.gospl_driver->record_coupling(total_dt, min_dh, max_dh, sum_dh / ntop);
        if (log_event)
            std::cout << "GoSPL: dh [" << min_dh << ", " << max_dh << "]"
                      << " mean=" << sum_dh / ntop << std::endl;
    }
#endif

    void get_surface_info(const Variables& var, \
        double_vec& top_base, double_vec& top_depth) {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const std::size_t ntop = top_nodes.size();

        for (std::size_t i=0; i<ntop-1;i++) {
            int n0 = top_nodes[i];
            int n1 = top_nodes[i+1];

            double dx = std::fabs(coord[n1][0] - coord[n0][0]);
            top_base[i] += dx;
            top_base[i+1] += dx;
        }

        for (std::size_t i=0; i<ntop; i++) {
            top_depth[i] = surfinfo.base_level - coord[top_nodes[i]][1];
            top_base[i] *= 0.5;
        }
#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }

    void out_basin_info(const double_vec& top_depth) {
        const int chart_width = 60;
        const int ntop = top_depth.size();
        const double ratio = double(chart_width) / ntop;

        int sum = 0, count_sign = 0;
        int_vec signs(chart_width);

        for (int i=0;i<ntop;i++) {
            double prec = i * ratio;
            if(prec > sum) {
                int diff = ceil(prec) - sum;
                sum += diff;
                for (int j=0;j<diff;j++)
                    if (top_depth[i] > 0.)
                        signs[count_sign++] = 0;
                    else
                        signs[count_sign++] = 1;
            }
        }
        std::cout << "\t";

        for (int i=0;i<chart_width;i++)
            if (signs[i] == 0)
                std::cout << " ";
            else
                if (i == 0 || i == (chart_width-1))
                    std::cout << " ";
                else if (signs[i+1] == 0)
                    std::cout << " ";
                else if (signs[i-1] == 0)
                    std::cout << " ";
                else
                    std::cout << "_";

        std::cout << std::endl << "\t";

        for (int i=0;i<chart_width;i++)
            if (signs[i] == 0)
                std::cout << "_";
            else
                if (i == 0 || i == (chart_width-1))
                    std::cout << "|";
                else if (signs[i+1] == 0)
                    std::cout << "\\";
                else if (signs[i-1] == 0)
                    std::cout << "/";
                else
                    std::cout << " ";

        std::cout << " (ntop: " << std::setw(3) << ntop << ")\n";
    }


    void terrigenous_diffusion(const Param& param,const Variables& var, const double_vec& basin_x, const double_vec& basin_dx, const double_vec& basin_depth, \
            const int nbasin,const int option, double_vec& dh_terrig, double dt_cycle) {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif
        const double S0 = param.control.terrig_sediment_area;
        const double C0 = param.control.terrig_sediment_diffusivity;
        const double C1 = param.control.terrig_depth_coefficient;
        const double coeff = dt_cycle * C0;
        double_vec basin_slope(nbasin+1);

        // calculate basin slope
        for (int i=0; i<nbasin+1;i++)
            basin_slope[i] = -(basin_depth[i+1] - basin_depth[i]) / (basin_x[i+1] - basin_x[i]);

        // set boundary condition of basin_slope
        if (option == 0)
            basin_slope[0] = - S0 / C0;
        else
            basin_slope[nbasin] = S0 / C0;

        // calculate dh for depth-dependent terrigenous diffusion
        for (int i=0; i<nbasin;i++)
            dh_terrig[i] = coeff * exp(-C1*basin_depth[i+1]) * (basin_slope[i+1] - basin_slope[i]) / basin_dx[i];

        // set boundary condition if basin is larger than 1
        if (nbasin > 1) {
            if (option == 0)
                dh_terrig[nbasin-1] = 0.;
            else
                dh_terrig[0] = 0.;
        }

        for (int i=0;i<nbasin;i++)
            // make sure the sedimentation is positive
            if (dh_terrig[i] < 0.)
                dh_terrig[i] = 0.;
            // avoid deposition above the base level
            else if (dh_terrig[i] > basin_depth[i+1])
                dh_terrig[i] = basin_depth[i+1]+1e-2;

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
        return;
    }

    int* find_basin(const double_vec& depth_tmp,const int ntop, const int option) {

        int *basin = new int[2];
        basin[0] = -1;
        basin[1] = -1;

        double boundary[ntop-1];
        for (int i=0;i<(ntop-1);i++)
            boundary[i] = depth_tmp[i] * depth_tmp[i+1];

        bool flag = false;

        if (option == 0) {
            for (int i=0;i<(ntop-1);i++) {
                if (boundary[i] <= 0. && depth_tmp[i] < 0.) {
                    basin[0] = i+1;
                    for (int j=i+1;j<ntop;j++) {
                        if (boundary[j] <= 0. && depth_tmp[j+1] < 0.) {
                            basin[1] = j;
                            flag = true;
                            break;
                        } else if (boundary[j] <= 0.) {
                            break;
                        }
                    }
                    if (flag) break;
                }
            }
        } else {
            for (int i=ntop-2;i>-1;i--) {
                if (boundary[i] <= 0. && depth_tmp[i+1] < 0.) {
                    basin[1] = i;
                    for (int j=i-1;j>=0;j--) {
                        if (boundary[j] <= 0. && depth_tmp[j] < 0.) {
                            basin[0] = j+1;
                            flag = true;
                            break;
                        } else if (boundary[j] <= 0.) {
                            break;
                        }
                    }
                    if (flag) break;
                }
            }
        }
        return basin;
    }

    void hemipelagic_deposition(const Param& param,const Variables& var) {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif
        #pragma acc wait

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const int ntop = surfinfo.ntop;
        double_vec& dh = *var.surfinfo.dh;

        const bool is_report = param.control.is_reporting_terrigenous_info && var.steps%10000 == 0;
        // hemipelagic sedimentation (Emmerich et al., 2009)
        // maximum accumulation rate (m/s)
        const double dh_hemipelagic = param.control.hemipelagic_sedimentation_rate * var.dt;
        // pelagic sedimentation (Emmerich et al., 2009)
        // maximum accumulation rate (m/s)
        const double dh_pelagic = param.control.pelagic_sedimentation_rate * var.dt;

        double max_dhi = 0.;

        for (int i=0;i<ntop;i++) {
            double depth = surfinfo.base_level - coord[top_nodes[i]][1];
            if (depth > 0. ) {
                double dhi = dh_hemipelagic * exp(-pow((depth - param.control.hemipelagic_max_depth)/param.control.hemipelagic_width,2));
                dhi += dh_pelagic * (1. - exp(-pow(depth/param.control.pelagic_increasing_width,2)));
                dh[i] += dhi;
                if (dhi > max_dhi)
                    max_dhi = dhi;
            }
        }
        if (is_report)
            std::cout << "\tMax hemi/pelagic sedimentation rate (km/Myr): " << std::fixed << std::setprecision(3) << \
                max_dhi*YEAR2SEC*1e3/var.dt << std::endl;

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }


    void terrigenous_deposition(const Param& param,const Variables& var) {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif
        #pragma acc wait // here is not ACC parallelized yet

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const int ntop = surfinfo.ntop;
        double_vec& dh = *var.surfinfo.dh;

        double_vec top_depth(ntop), dh_tmp(ntop);

        double max_depth = 0., min_depth = 0.;
        const bool is_report = param.control.is_reporting_terrigenous_info && var.steps%10000 == 0;
        int iter = 10;

        for (int i=0;i<ntop;i++) {
            dh_tmp[i] = 0.;
            top_depth[i] = surfinfo.base_level - coord[top_nodes[i]][1];
            if (top_depth[i] > max_depth)
                max_depth = top_depth[i];
            if (top_depth[i] < min_depth)
                min_depth = top_depth[i];
        }

        if (is_report) {
            std::cout << "\t** Sedimentation report: \n";
            out_basin_info(top_depth);
        }

        if (max_depth*min_depth >= 0.) {
            if (is_report)
                std::cout << "\tNo basin exist.\n";
#ifdef NPROF_DETAIL
            nvtxRangePop();
#endif
            return;
        }

        // run source from both side
        for (int iside=0;iside<2;iside++){
            double dt_next = 0.;
            int itry = 0;

            // decreae dt of terrigenous sedimentation to increase the stability
            do {
                double_vec depth_tmp(ntop);
                for (int j=0;j<ntop;j++)
                    depth_tmp[j] = top_depth[j] - dh_tmp[j];

                int *basin = find_basin(depth_tmp,ntop,iside);

                // check if basin exist
                if (basin[0] == -1 || basin[1] == -1) {
                    if (is_report)
                        std::cout << "\tSource " << iside << " : No accommodation space.\n";
                    break;
                } else {
                    if (is_report && itry == 0)
                        std::cout << "\tSource " << iside << " : Basin at " << basin[0] << " " << basin[1] << "\n";
                }

                itry++;
                int nbasin = basin[1] - basin[0] + 1;
                double_vec basin_x(nbasin+2), basin_dx(nbasin);
                double_vec basin_depth(nbasin+2), dh_basin(nbasin);
                double basin_area = 0.;
                double dt_cycle = var.dt / iter;
                if (dt_next > 0.){
                    dt_cycle = dt_next;
                    dt_next = 0.;
                }
                double area_ref = param.control.terrig_sediment_area * dt_cycle;

                for (int j=0;j<nbasin+2;j++) {
                    basin_x[j] = coord[top_nodes[basin[0]+j-1]][0];
                    basin_depth[j] = depth_tmp[basin[0]+j-1];
                }
                // loops over all top facets
                for (int j=0;j<nbasin;j++) {
                    basin_dx[j] = std::fabs(basin_x[j+2] - basin_x[j]) / 2.;
                    // calculate the area of the basin
                    basin_area += basin_dx[j] * (basin_depth[j+1] + 1e-2);
                }

                // if the area of the basin is smaller than the reference area,
                // use basin depth as the dh
                if (basin_area <= area_ref) {
                    itry -= 1;
                    dt_next = dt_cycle * (1. - basin_area/area_ref);
                    printf("dt_next %e dt_cycle %e\n",dt_next,dt_cycle);
                    for (int j=0; j<nbasin; j++)
                        dh_basin[j] = basin_depth[j+1] + 1e-2;
                } else
                    // calculate the terrigenous sedimentation
                    terrigenous_diffusion(param,var,basin_x,basin_dx,basin_depth,nbasin,iside,dh_basin,dt_cycle);

                for (int j=0; j<nbasin; j++)
                    dh_tmp[basin[0]+j] += dh_basin[j];

                delete [] basin;
            } while (itry < iter);
            if (itry == 0) break;
        }

        for (int i=0; i<ntop; i++)
            dh[i] += dh_tmp[i];

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }
}


void surface_plstrain_diffusion(const Variables& var, double_vec& plstrain)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    double half_life = 1.e2 * YEAR2SEC;
    double lambha = 0.69314718056 / half_life; // ln2
    // Top elements only, but a host loop here would drag plstrain off the device on every call.
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, plstrain, lambha)
#endif
    #pragma acc parallel loop gang vector async
    for (int ti = 0; ti < var.ntop_elems; ++ti) {
        const int e = (*var.top_elems)[ti];
        plstrain[e] -= plstrain[e] * lambha * var.dt;
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void correct_surface_element(const Variables& var, double_vec& volume, double_vec& volume_n, \
    tensor_t& stress, tensor_t& strain, tensor_t& strain_rate, double_vec& plstrain)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
#ifndef ACC
    #pragma omp parallel default(none) \
        shared(var, volume, volume_n, stress, strain, strain_rate, plstrain)
#endif
    {
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop gang vector async
        for (int i=0; i<var.ntop_elems; i++) {
            const int e = (*var.top_elems)[i];

            ConstArrayIndirectAccessor coord1 = var.coord->view_const((*var.connectivity)[e]);
            double new_volumes = compute_volume(coord1);
            double rdv =  new_volumes / volume[e];
            volume[e] = new_volumes;
            if (rdv < 1.0) continue;

            // correct stress and strain
            plstrain[e] /= rdv;
            // correct the plastic strain overestimation of surface element caused by sedimentation.
            for (int j=0;j<NSTR;j++) {
                stress[e][j] /= rdv;
                strain[e][j] /= rdv;
                strain_rate[e][j] /= rdv;
            }
        }

#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop gang vector async
        for (int n=0;n<var.surfinfo.ntop;n++) {
            int nt = (*var.surfinfo.top_nodes)[n];
            const int npatch = var.support.size(nt);
            const int* patch = var.support.patch(nt);
            double acc = 0.;
            for (int i=0;i<npatch;i++)
                acc += volume[patch[i]];
            volume_n[nt] = acc;
        }
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void surface_processes(const Param& param, const Variables& var, array_t& coord, tensor_t& stress, tensor_t& strain, \
                       tensor_t& strain_rate, double_vec& plstrain, double_vec& volume, double_vec& volume_n, \
                       SurfaceInfo& surfinfo, std::vector<MarkerSet*> &markersets, \
                       int_vec2D& elemmarkers, int_vec2D& markers_in_elem)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0;i<var.surfinfo.ntop;i++)
        (*var.surfinfo.dh)[i] = 0.;

    switch (param.control.surface_process_option) {
    case 0:
        // no surface process
        break;
    case 1:
        simple_diffusion(var);
        break;
#ifdef HAS_GOSPL_CPP_INTERFACE
    case 11:
        if (var.gospl_driver != nullptr) {
            use_gospl(param, var);
        } else {
            std::cerr << "Error: GoSPL driver not initialized for surface process option 11" << std::endl;
        }
        break;
#endif
    case 101:
        custom_surface_processes(var, coord);
        break;
    case 102:
#ifdef THREED
        die(EXIT_UNSUPPORTED_DIM, "3D deposition of sediment processes is not ready yet.");
#else
        simple_diffusion(var);
        if (var.steps != 0) {
            terrigenous_deposition(param, var);
            if (param.control.hemipelagic_sedimentation_rate > 0. || param.control.pelagic_sedimentation_rate > 0.)
                hemipelagic_deposition(param, var);
        }
#endif
        break;
    default:
        std::cout << "Error: unknown surface process option: " << param.control.surface_process_option << '\n';
        die(EXIT_CONFIG_VALUE);
    }

#ifndef ACC
    #pragma omp parallel default(none) shared(param, var, coord)
#endif
    {
        // go through all surface nodes and abject all surface node by dh
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop gang vector async
        for (int i=0; i<var.surfinfo.ntop; i++) {
            // get global index of node
            // update coordinate via dh
            // update dhacc for marker correction
            int nt = (*var.surfinfo.top_nodes)[i];
            coord[nt][NDIMS-1] += (*var.surfinfo.dh)[i];
            (*var.surfinfo.dhacc)[nt] += (*var.surfinfo.dh)[i];
        }

        // update edvacc_surf of connected elements
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop gang vector async
        for (int i=0;i<var.surfinfo.etop;i++) {
            double dh_e = 0.;

            for (int j=0;j<NDIMS;j++)
                dh_e += (*var.surfinfo.dh)[(*var.surfinfo.elem_and_nodes)[i][j]];

            ConstArrayIndirectAccessor coord_surf = var.coord->view_const((*var.connectivity_surface)[i]);
            double base = compute_area_facet(coord_surf);

            (*var.surfinfo.edvacc_surf)[i] += dh_e * base / NDIMS;
        }
    }

//     if (var.steps%1000 == 0) {
//         double sum = 0.;
// #ifndef ACC
//         #pragma omp parallel default(none) shared(var) reduction(+:sum)
// #endif
//         #pragma acc parallel loop gang vector async reduction(+:sum)
//         for (int i=0; i<var.surfinfo.etop; i++) {
//             sum += (*var.surfinfo.edvacc_surf)[i];
//         }
//         #pragma acc wait
//         printf("Sum of accu. surface deposition volume: %.2e\n", sum);
//     }

    // find max surface velocity
    double maxdh = 0.;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var) reduction(max:maxdh)
#endif
    #pragma acc parallel loop gang vector reduction(max:maxdh) async
    for (int i=0; i<var.surfinfo.ntop; ++i) {
        double tmp = fabs((*var.surfinfo.dh)[i]);

        if (maxdh < tmp)
            maxdh = tmp;
    }

    #pragma acc wait

    surfinfo.max_surf_vel = maxdh / var.dt;

    correct_surface_element(var, volume, volume_n, stress, strain, strain_rate, plstrain);

    if (var.steps != 0) {
        // The main loop rebuilds the serpentinization grid under this same step
        // test, because it must index every marker; change the two together.
        if ( var.steps % param.mesh.quality_check_step_interval == 0) {
            markersets[0]->correct_surface_marker(param, var, *var.surfinfo.dhacc, elemmarkers, markers_in_elem);

#ifndef ACC
            #pragma omp parallel for default(none) shared(var)
#endif
            #pragma acc parallel loop gang vector
            for (int i=0;i<var.surfinfo.ntop;i++)
                (*var.surfinfo.dhacc)[(*var.surfinfo.top_nodes)[i]] = 0.;

            // set marker of sediment.
            markersets[0]->set_surface_marker(param, var, param.mesh.smallest_size, param.mat.mattype_sediment, *var.surfinfo.edvacc_surf, elemmarkers, markers_in_elem);
        }
    }

#ifdef THREED
    // todo
#else
    if (!(var.steps % param.mesh.quality_check_step_interval &&  var.steps != 0))
        // correct the plastic strain of urface element for preventing surface landslide.
        surface_plstrain_diffusion(var,plstrain);

    if ( param.control.is_reporting_terrigenous_info && var.steps%10000 == 0 &&  var.steps != 0  ) {
        double max_dh = 0., min_dh = 0.;
        for (int i=0;i<var.surfinfo.ntop;i++) {
            max_dh = std::max(max_dh, (*var.surfinfo.dh)[i]);
            min_dh = std::min(min_dh, (*var.surfinfo.dh)[i]);
            // std::cout << n << "  dh:  " << dh << '\n';
        }

        std::cout << "\tMax erosion / sedimentation rate (mm/yr):  "
                    << std::fixed << std::setprecision(3) << min_dh / var.dt * 1000. * YEAR2SEC << " / "
                    << std::fixed << std::setprecision(3) << max_dh / var.dt * 1000. * YEAR2SEC << '\n';

    }
#endif

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif


}

void build_velocity_constraints(const Param& param, const Variables& var,
                                VelocityConstraints& constraints, bool homogeneous)
{
    // F(v)=Jv+b is the existing ordered BC map. Its admissible velocities solve
    // (I-J)v=b. Homogeneous calls give J without subtracting prescribed velocities,
    // avoiding cancellation when a load is large relative to a perturbation.
    #pragma acc wait
    const int nn = var.nnode;
    std::unique_ptr<array_t> basis_owner(new array_t(nn, 0.0));
    array_t& basis = *basis_owner;
    auto& p = constraints.projector;
    auto& offset = constraints.prescribed;
    for (int ib=iboundn0; ib<=iboundn3; ++ib) {
        if (all_free(param.bc.side[ib])) continue;
        for (int n=0; n<nn; ++n) {
            if (!((*var.bcflag)[n] & (1U << ib))) continue;
            double norm2 = 0;
            for (int d=0; d<NDIMS; ++d) norm2 += (*var.bnormals)[ib][d]*(*var.bnormals)[ib][d];
            if (!std::isfinite(norm2) || std::abs(norm2-1) > 1e-10)
                die(EXIT_CONFIG_VALUE, "PT oblique boundary requires a finite unit normal.");
            const double nz = (*var.bnormals)[ib][NDIMS-1];
            if (param.bc.side[ib].horizontal_normal && !(1-nz*nz > 1e-14))
                die(EXIT_CONFIG_VALUE, "PT horizontal oblique BC requires a nonzero horizontal normal.");
        }
    }
    for (int n=0; n<nn; ++n)
        for (int d=0; d<NDIMS; ++d) offset[n][d] = 0;
    apply_vbcs(param, var, offset, homogeneous, true);
    for (int column=0; column<NDIMS; ++column) {
        for (int n=0; n<nn; ++n)
            for (int d=0; d<NDIMS; ++d) basis[n][d] = d == column ? 1 : 0;
        apply_vbcs(param, var, basis, true, true);
        for (int n=0; n<nn; ++n)
            for (int d=0; d<NDIMS; ++d)
                p[n][d*NDIMS+column] = (d == column ? 1.0 : 0.0) - basis[n][d];
    }
    // At most NDIMS rows: reorthogonalized elimination supplies both the minimum
    // norm affine solution and an orthogonal free-space projector, including edges.
    for (int n=0; n<nn; ++n) {
        double q[NDIMS][NDIMS] = {}, rhs[NDIMS] = {};
        int rank = 0;
        double scale = 0;
        for (int d=0; d<NDIMS; ++d) scale = std::max(scale, std::abs(offset[n][d]));
        for (int row=0; row<NDIMS; ++row) {
            double a[NDIMS], b = offset[n][row];
            for (int d=0; d<NDIMS; ++d) a[d] = p[n][row*NDIMS+d];
            for (int pass=0; pass<2; ++pass)
                for (int k=0; k<rank; ++k) {
                    double dot = 0;
                    for (int d=0; d<NDIMS; ++d) dot += a[d]*q[k][d];
                    for (int d=0; d<NDIMS; ++d) a[d] -= dot*q[k][d];
                    b -= dot*rhs[k];
                }
            double norm2 = 0;
            for (int d=0; d<NDIMS; ++d) norm2 += a[d]*a[d];
            if (norm2 <= 1e-24) {
                if (std::abs(b) > 1e-12*std::max(scale, std::numeric_limits<double>::min()))
                    die(EXIT_CONFIG_VALUE, "PT velocity boundary constraints have no consistent fixed point.");
                continue;
            }
            const double norm = std::sqrt(norm2);
            for (int d=0; d<NDIMS; ++d) q[rank][d] = a[d]/norm;
            rhs[rank++] = b/norm;
        }
        constraints.free_rank[n] = NDIMS-rank;
        for (int d=0; d<NDIMS; ++d) {
            offset[n][d] = 0;
            for (int k=0; k<rank; ++k) offset[n][d] += rhs[k]*q[k][d];
            for (int j=0; j<NDIMS; ++j) {
                p[n][d*NDIMS+j] = d == j ? 1 : 0;
                for (int k=0; k<rank; ++k) p[n][d*NDIMS+j] -= q[k][d]*q[k][j];
            }
        }
    }
}

void project_free_vectors(const Variables& var, const VelocityConstraints& constraints,
                          array_t& vectors, bool include_prescribed)
{
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, constraints, vectors) firstprivate(include_prescribed)
#endif
    #pragma acc parallel loop gang vector async
    for (int n=0; n<var.nnode; ++n) {
        double x[NDIMS];
        for (int d=0; d<NDIMS; ++d) x[d] = vectors[n][d];
        for (int d=0; d<NDIMS; ++d) {
            double value = include_prescribed ? constraints.prescribed[n][d] : 0;
            for (int j=0; j<NDIMS; ++j) value += constraints.projector[n][d*NDIMS+j]*x[j];
            vectors[n][d] = value;
        }
    }
    #pragma acc wait
}
