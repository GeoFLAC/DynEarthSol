#include <algorithm>
#include <array>
#include <climits>
#include <cstring>
#include <functional>
#include <iostream>
#include <numeric>
#include <limits>
#include <map>
#include <unordered_map>

#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
#include "brc-interpolation.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "mmg_utils.hpp"
#include "nn-interpolation.hpp"
#include "utils.hpp"
#include "markerset.hpp"
#include "remeshing.hpp"

#ifdef USEMMG
#ifdef THREED
#include "mmg/mmg3d/libmmg3d.h"
#else
#include "mmg/mmg2d/libmmg2d.h"
#endif
#endif

#ifdef ACC
#include "knn_bvh.hpp"
#endif

namespace { // anonymous namespace


// equilateral triangle area = 0.433*s^2, equilateral tetrahedron volume = 0.118*s^3
#ifdef THREED
const double sizefactor = 0.118;
#else
const double sizefactor = 0.433;
#endif

// Shared refine-policy size floors, derived in ONE place so the triangle decimation, the MMG
// size band, the tiny-element remesh trigger, and the boundary-merge distance cannot drift apart.
// (These were previously spelled out at six sites with subtly different sizefactor placement -- the
// exact class of bug that caused earlier over/under-refinement.) All are in the mesh's own units.
struct RefineFloors {
    double hmin;          // min edge length      = smallest_size^(1/NDIMS) * resolution
    double smallest_vol;  // min element measure  = smallest_size * sizefactor * resolution^NDIMS
    double min_dist;      // boundary-merge dist  = (smallest_size*sizefactor)^(1/NDIMS) * resolution
};
RefineFloors refine_floors(const Mesh &m)
{
    RefineFloors f;
    f.hmin         = std::pow(m.smallest_size, 1.0/NDIMS) * m.resolution;
    f.smallest_vol = m.smallest_size * sizefactor * std::pow(m.resolution, NDIMS);
    f.min_dist     = std::pow(m.smallest_size * sizefactor, 1.0/NDIMS) * m.resolution;
    return f;
}

const int DELETED_FACET = -1;
const int DEBUG = 0;

bool is_boundary(uint flag)
{
    return flag & BOUND_ANY;
}


#pragma acc routine seq
bool is_bottom(uint flag)
{
    return flag & BOUNDZ0;
}


#pragma acc routine seq
bool is_x0(uint flag)
{
    return flag & BOUNDX0;
}


#pragma acc routine seq
bool is_x1(uint flag)
{
    return flag & BOUNDX1;
}


#pragma acc routine seq
bool is_y0(uint flag)
{
    return flag & BOUNDY0;
}


#pragma acc routine seq
bool is_y1(uint flag)
{
    return flag & BOUNDY1;
}

bool is_corner(uint flag)
{
    uint f = flag & BOUND_ANY;
    if (!f) return 0;

    // A corner node will have multiple bits (2 in 2D; 3 or more in 3D) set in its flag.
    int nbits = 0;
    for (int j=0; j<nbdrytypes; j++) {
        // counting how many bits are set
        if (f & (1<<j)) nbits++;
    }

#ifdef THREED
    return (nbits >= NDIMS);
#else
    return (nbits == NDIMS);
#endif
}


bool is_bottom_corner(uint flag)
{
    if ((flag & BOUNDZ0) && is_corner(flag)) return 1;
    return 0;
}

void flatten_bottom(const uint_vec &old_bcflag, double *qcoord,
                    double bottom, int_vec &points_to_delete, double min_dist)
{
    // find old nodes that are on or close to the bottom boundary

    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_bottom(flag)) {
            // restore edge nodes to initial depth
            qcoord[i*NDIMS + NDIMS-1] = bottom;
        }
        else if (flag & BOUNDZ1 && qcoord[i*NDIMS + NDIMS-1] < bottom + min_dist) {
            // A TOP node within min_dist of (or below) the restored bottom: the domain has necked
            // through. Collapsing it would weld the free surface onto the floor (garbage mesh, NaN);
            // stop instead.
            std::cerr << "Error: top-surface node " << i << " (z = "
                      << qcoord[i*NDIMS + NDIMS-1] << ") has reached the bottom boundary (z = "
                      << bottom << ").\n       The model domain has necked through; the mesh "
                      "can no longer represent the model. Stopping.\n";
            die(EXIT_MESH_QUALITY);   // terminal physical state reached in the mesh-repair path
        }
        else if (qcoord[i*NDIMS + NDIMS-1] < bottom + min_dist) {
            // Delete every non-bottom node within min_dist of (or below) the bottom: an interior
            // one inverts its elements once the bottom snaps back up, a side-wall one below the
            // corner degenerates the corner elements. The MMG collapse or
            // delete_points_and_merge_boundary (see its keep rule) removes it.
            points_to_delete.push_back(i);
        }
    }
}

void flatten_x0(const uint_vec &old_bcflag, double *qcoord,
                int_vec &points_to_delete, double min_dist)
{
    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_x0(flag))
            qcoord[i*NDIMS] = 0.;
        // Delete any non-x0 node within min_dist of (or past) the plane, a top or bottom node
        // carried past the corner too (skipping boundary nodes leaves an orphan wedge there).
        // flatten_x1/y0/y1 do the same.
        else if (! is_x0(flag))
            if (qcoord[i*NDIMS] < min_dist)
                points_to_delete.push_back(i);
    }
}

void flatten_x1(const uint_vec &old_bcflag, double *qcoord,
                double side, int_vec &points_to_delete, double min_dist)
{
    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_x1(flag))
            qcoord[i*NDIMS] = side;
        else if (! is_x1(flag))
            if (qcoord[i*NDIMS] > (side - min_dist))
                points_to_delete.push_back(i);
    }
}

void flatten_y0(const uint_vec &old_bcflag, double *qcoord,
                    int_vec &points_to_delete, double min_dist)
{
    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_y0(flag))
            qcoord[i*NDIMS + 1] = 0.;
        else if (! is_y0(flag))
            if ( qcoord[i*NDIMS + 1] < min_dist)
                points_to_delete.push_back(i);
    }
}

void flatten_x0_corner(const uint_vec &old_bcflag, double *qcoord,
                int_vec &points_to_delete)
{

    // for cutting x0 boundary to (x0 & z1) in X
    double x1x, bx0;

    int j=0;
    while (j>-1) {
        uint flag = old_bcflag[j];
        if ((flag & BOUNDX0 ) && (flag & BOUNDZ1)) {
            x1x = qcoord[j*NDIMS]; break;
        }
        j++;
    }

    // interpolation or extrapolation 4 bx0
    j=0;
    int l=0;
    double bo_X[2], bo_depth[2];
    do {
        uint flag = old_bcflag[j];
        if (is_bottom(flag)) {
           bo_depth[l]  = qcoord[j*NDIMS + NDIMS-1];
           bo_X[l] = qcoord[j*NDIMS];
           l++;
        }
        j++;
    } while (l<2);

    bx0 = bo_depth[0] + (bo_depth[1]-bo_depth[0])*(x1x-bo_X[0])/(bo_X[1]-bo_X[0]);

    std::cout << "x1x: " << x1x << "; bx0: "<<bx0<< '\n';
    double x0_zmin = std::numeric_limits<double>::max();
    double bot_xmax = std::numeric_limits<double>::lowest();
    double b0_exceed_xmax = std::numeric_limits<double>::lowest();

    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_x0(flag)) {
            if ((x0_zmin > qcoord[i*NDIMS + NDIMS-1]) && !(qcoord[i*NDIMS] > x1x))
                x0_zmin = qcoord[i*NDIMS + NDIMS-1];
            else if (is_bottom(flag))
                b0_exceed_xmax = qcoord[i*NDIMS];
            // set all x0 to x1x in X
            qcoord[i*NDIMS] = x1x;
        }
    }

    double z0_xmax_ref = 2*x1x - b0_exceed_xmax;

    for (std::size_t i=0; i<old_bcflag.size(); ++i)
        if (is_bottom(old_bcflag[i]))
            if ((bot_xmax < qcoord[i*NDIMS]) && (qcoord[i*NDIMS] < z0_xmax_ref))
                bot_xmax = qcoord[i*NDIMS];

    double v = (x1x - bot_xmax) / (x0_zmin - bx0);

#ifdef USEMMG
    double shrink_x = (x1x - bot_xmax) / (b0_exceed_xmax - bot_xmax);
#else
    double b0_clean_zmin_ref = x0_zmin - (x0_zmin - bx0)/2.;
#endif

    for (std::size_t i=0; i < old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        double x = qcoord[i*NDIMS];
        double y = qcoord[i*NDIMS + NDIMS-1];
        if (!is_x0(flag)) {
            double fx = v * (y - x0_zmin) + x1x;
            if (fx < x) {
#ifdef USEMMG
                // adjust all nodes at the right side of edge line
                double dx = (x - fx) * shrink_x;
                qcoord[i*NDIMS] = fx + dx;
#else
                // remove all nodes at the right side of edge line
                points_to_delete.push_back(i);
#endif
            }
#ifdef USEMMG
#else
        } else if (!is_bottom(flag)) {
            // remove b0 nodes close to the z0-b0 corner
            if (y < b0_clean_zmin_ref)
                points_to_delete.push_back(i);
#endif
        }
    }
}

void flatten_y1(const uint_vec &old_bcflag, double *qcoord,
                double side, int_vec &points_to_delete, double min_dist)
{
    // Mirror of flatten_x1 for y-direction
    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_y1(flag))
            qcoord[i*NDIMS + 1] = side;
        else if (! is_y1(flag))
            if (qcoord[i*NDIMS + 1] > (side - min_dist))
                points_to_delete.push_back(i);
    }
}


void new_bottom(const uint_vec &old_bcflag, double *qcoord,
                double bottom_depth, int_vec &points_to_delete, double min_dist,
                int *segment, int *segflag, int nseg)
{
    /* deleting nodes that are on or close to the bottom boundary,
     * excluding nodes on the side walls
     */

#ifdef THREED
    /* In 3D, if bottom nodes were deleted, the facets on the side boundaries
     * near the bottom are affected as well. The code does not deal with this
     * complexity and can only work in 2D.
     */
    die(EXIT_UNSUPPORTED_DIM, "new_bottom() does not work in 3D.");
#endif

    int_vec bottom_corners;
    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_bottom(flag)) {
            if(is_bottom_corner(flag))
                bottom_corners.push_back(i);
            else
                points_to_delete.push_back(i);
        }
        else if (! is_boundary(flag) &&
                 std::fabs(qcoord[i*NDIMS + NDIMS-1] - bottom_depth) < min_dist) {
            points_to_delete.push_back(i);
        }
    }

    if (DEBUG) {
        std::cout << "bottom points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
        std::cout << "segment before delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag before delete: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }

    // must have 2 bottom corners in 2D
    if (bottom_corners.size() != 2) {
        std::cerr << "Error: cannot find all bottom corners before remeshing. n_bottom_corners = "
                  << bottom_corners.size() << " (2 expected).\n";
        std::cout << "bottom corners: ";
        print(std::cout, bottom_corners);
        std::cout << '\n';
        die(EXIT_MESH_QUALITY);
    }

    // move the corners to the same depth
    for (std::size_t i=0; i<bottom_corners.size(); i++) {
        int n = bottom_corners[i];
        qcoord[n*NDIMS + NDIMS-1] = bottom_depth;
    }

    // mark all bottom nodes as deleted
    for (int i=0; i<nseg; ++i) {
        if (static_cast<uint>(segflag[i]) == BOUNDZ0) {
            for (int j=0; j<NODES_PER_FACET; j++)
                segment[i*NODES_PER_FACET + j] = DELETED_FACET;
        }
    }

    // create new bottom segments from corner nodes
    for (int i=0; i<nseg; ++i) {
        if (static_cast<uint>(segflag[i]) == BOUNDZ0) {
            segment[i*NODES_PER_FACET + 0] = bottom_corners[0];
            segment[i*NODES_PER_FACET + 1] = bottom_corners[1];
            break;
        }
    }

    if (DEBUG) {
        std::cout << "bottom corners: ";
        print(std::cout, bottom_corners);
        std::cout << '\n';
        std::cout << "segment with new bottom: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag with new bottom: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }
}


// local cmp functor
struct cmp {
    const array_t &coord;
    const int d;
    cmp (const array_t &coord_, int dim) : coord(coord_), d(dim) {};
    bool operator()(const int &a, const int &b) {return coord[a][d] < coord[b][d];}
};


typedef std::pair<int,int> edge_t;
struct equal_to1
{
    bool operator()(const edge_t &lhs, const edge_t &rhs) const {
        return (lhs.first == rhs.first && lhs.second == rhs.second);
    }
};


struct hash1
{
    std::size_t operator()(const edge_t &k) const {
        // first as the upper half
        const int halfbytes = 4;
        return (static_cast<std::size_t>(k.first) << sizeof(std::size_t)*halfbytes) | k.second;
    }
};


void find_tiny_element(const Param &param, const double_vec &volume,
                       int_vec &tiny_elems)
{
    const double smallest_vol = refine_floors(param.mesh).smallest_vol;

    for (std::size_t e=0; e<volume.size(); e++) {
        if (volume[e] < smallest_vol)
            tiny_elems.push_back(e);
    }

    if (DEBUG) {
        std::cout << "tiny elements: ";
        print(std::cout, tiny_elems);
        std::cout << '\n';
    }
}


void find_points_of_tiny_elem(const array_t &coord, const conn_t &connectivity,
                              const double_vec &volume, const int_vec &tiny_elems,
                              int npoints, const array_t &old_points,
                              const uint_vec &old_bcflag, int_vec &points_to_delete,
                              bool excl_func(uint))
{
    // collecting the nodes of tiny_elems
    int tiny_nelem = tiny_elems.size();
    array_t tiny_coord(tiny_nelem * NODES_PER_ELEM);
    conn_t tiny_conn(tiny_nelem);
    double_vec tiny_vol(tiny_nelem);
    int ii = 0;
    for (int ee=0; ee<tiny_nelem; ++ee) {
        int e = tiny_elems[ee];

        tiny_vol[ee] = volume[e];

        ConstArrayIndirectAccessor coord_e = coord.view_const(connectivity[e]);
        for (int j=0; j<NODES_PER_ELEM; ++j) {
            tiny_conn[ee][j] = ii;

            for (int d=0; d<NDIMS; ++d) {
                tiny_coord[ii][d] = coord_e[j][d];
            }
            ii ++;
        }
    }

    Barycentric_transformation bary(tiny_coord, tiny_conn, tiny_vol);

    #pragma acc wait

    // find old nodes that are connected to tiny elements and are not excluded
    // (most of the nodes of tiny elements are newly inserted by the remeshing library)
    for (int i=0; i<npoints; ++i) {
        // excluded nodes
        if (excl_func(old_bcflag[i])) continue;

        ConstArrayAccessor p = old_points[i];
        for (int ee=0; ee<tiny_nelem; ++ee) {
            if (bary.is_inside_elem(p, ee)) {
                points_to_delete.push_back(i);
                break;
            }
        }
    }

    if (DEBUG) {
        std::cout << "points of tiny elements: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
    }
}


// Swap-delete points_to_delete. A deleted point must be in no segment (its segments would end on
// the node moved into its slot): option 2's bottom nodes qualify, as new_bottom replaced their
// segments; other boundary points go through delete_points_and_merge_boundary.
void delete_points(const int_vec &points_to_delete, int &npoints,
                   int nseg, double *points, int *segment)
{
    if (points_to_delete.size() == 0) return;

    if (DEBUG) {
        std::cout << "old points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
    }

    int *endsegment = segment + nseg * NODES_PER_FACET;

    int end = npoints - 1;

    // delete points from the end
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
        // when a point is deleted, replace it with the last point
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        std::replace(segment, endsegment, end, *i);
        // std::cout << *i << " <- " << end << "\n";

        end --;
    }
    npoints -= points_to_delete.size();
}


void delete_facets(int &nseg, int *segment, int *segflag)
{
    // delete facets from the end
    for (int i=nseg-1; i>=0; i--) {
        if (segment[i*NODES_PER_FACET] == DELETED_FACET) {
            // safety check
            if (segment[i*NODES_PER_FACET + 1] != DELETED_FACET
#ifdef THREED
                || segment[i*NODES_PER_FACET + 2] != DELETED_FACET
#endif
                ) {
                std::cerr << "Error: segment array is corrupted before delete_facets()!\n";
                print(std::cerr, segment, nseg*NODES_PER_FACET);
                die(EXIT_MESH_QUALITY);
            }

            // replace deleted segment with the last segment
            for (int j=0; j<NODES_PER_FACET; ++j) {
                segment[i*NODES_PER_FACET + j] = segment[(nseg-1)*NODES_PER_FACET + j];
            }
            segflag[i] = segflag[nseg-1];
            nseg --;
        }
    }

    if (DEBUG) {
        std::cout << "segment: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }
}


void delete_points_and_merge_segments(const int_vec &points_to_delete, int &npoints,
                                      int nseg, double *points, int *segment,
                                      uint_vec &bcflag, double min_length)
{
#ifdef THREED
    die(EXIT_INTERNAL_ASSERT, "delete_points_and_merge_segments() doesn't work in 3D!");
#endif

    int *endsegment = segment + nseg * NODES_PER_FACET;

    int end = npoints - 1;

    // delete points from the end
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
        if (DEBUG) {
            std::cout << " deleting point " << *i << " replaced by point " << end << '\n';
        }

        // if the deleted point is a segment point, merge the neighboring two segments
        uint flag = bcflag[*i];
        if (is_boundary(flag)) {
            // segment and points:  aa------a = b-------bb

            int *a = std::find(segment, endsegment, *i);
            int *b = std::find(a+1, endsegment, *i);
            if (b == endsegment) {
                die(EXIT_MESH_QUALITY, "segment array is corrupted when merging segment!");
            }

            // a could be the either first or second node of the segment,
            // which will affect the location of aa in the segment array.
            bool not_first;
            not_first = (a - segment) % NODES_PER_FACET;
            int *aa = not_first? (a - 1) : (a + 1);
            not_first = (b - segment) % NODES_PER_FACET;
            int *bb = not_first? (b - 1) : (b + 1);

            // if the length of the two segments are too long
            // do not delete this point
            double la2, lb2;
            la2 = dist2(points + (*a)*NDIMS, points + (*aa)*NDIMS);
            lb2 = dist2(points + (*b)*NDIMS, points + (*bb)*NDIMS);
            if (la2 > min_length*min_length && lb2 > min_length*min_length) {
                if (DEBUG) {
                    std::cout << " the segments of point " << *i << " have length^2 "
                              << la2 << ", " << lb2 << " -- skip deletion."<< '\n';
                }
                continue;
            }

            // merge the two segments
            *a = *bb;
            *bb = DELETED_FACET;
            *b = DELETED_FACET;

            if (DEBUG) {
                std::cout << "a: " << (a-segment) << "  b: " << (b-segment)
                          << " not_first? " << not_first << " bb: " << (bb-segment)
                          << " = " << *bb << '\n';
                std::cout << "segment after merging: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        // when a point is deleted, replace it with the last point
        flag = bcflag[end];
        bcflag[*i] = flag;
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index needs to be
        // updated as well. Not gated on is_boundary(flag): a node of a kept internal segment
        // (is_discarding_internal_segments = no) has bcflag 0.
        std::replace(segment, endsegment, end, *i);
        if (DEBUG) {
            std::cout << "segment after replace: ";
            print(std::cout, segment, nseg*NODES_PER_FACET);
            std::cout << '\n';
        }

        end --;
        npoints --;
    }

    /* PS: We don't check whether the merged segments belong to the same boundary or not,
     *     e.g., one segment is on the top boundary while the other segment is on the left
     *     boundary. The check is performed in delete_facets() instead.
     */
}


#ifdef THREED
// Lawson edge flips over one boundary's facets (node triples). An edge flips only if its two
// diagonals cross in the projection plane, the flip raises the smaller 3D angle of its two facets,
// and the surface stays within max_change of the input along the projection axis: the move is
// measured where the diagonals cross and summed per facet, so flips cannot add up past it.
// Outline edges (one facet) never flip. Each flip strictly raises an angle, so the loop ends. The
// output is not oriented (TetGen orients the faces it meshes).
static void flip_boundary_facets(int_vec &facets, const double *points, int major_direction,
                                 double max_change)
{
    auto P = [&](int n, int d) { return points[n*NDIMS + d]; };
    const int u = (major_direction + 1) % NDIMS, v = (major_direction + 2) % NDIMS;
    auto orient = [&](int a, int b, int c) {   // twice the signed area in the projection plane
        return (P(b,u) - P(a,u)) * (P(c,v) - P(a,v)) - (P(b,v) - P(a,v)) * (P(c,u) - P(a,u));
    };
    auto min_angle = [&](int a, int b, int c) {
        const int n[3] = {a, b, c};
        double m = 4.0;   // above pi
        for (int k=0; k<3; ++k) {
            double d1 = 0, d2 = 0, dot = 0;
            for (int d=0; d<NDIMS; ++d) {
                const double e1 = P(n[(k+1) % 3], d) - P(n[k], d), e2 = P(n[(k+2) % 3], d) - P(n[k], d);
                d1 += e1*e1; d2 += e2*e2; dot += e1*e2;
            }
            m = std::min(m, std::acos(std::max(-1.0, std::min(1.0, dot / std::sqrt(d1*d2)))));
        }
        return m;
    };
    auto key = [](int a, int b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
    std::unordered_map<edge_t, std::array<int, 2>, hash1, equal_to1> efacets;   // edge -> its facets
    const int nf = facets.size() / NODES_PER_FACET;
    std::vector<double> moved(nf, 0.0);   // how far each facet's area has moved, along the axis
    for (int f=0; f<nf; ++f)
        for (int k=0; k<NODES_PER_FACET; ++k) {
            auto it = efacets.emplace(key(facets[3*f + k], facets[3*f + (k+1) % 3]),
                                      std::array<int, 2>{{f, -1}});
            if (!it.second) it.first->second[1] = f;
        }
    std::vector<edge_t> stack;
    for (const auto &e : efacets)
        if (e.second[1] >= 0) stack.push_back(e.first);
    auto third = [&](int f, int a, int b) {
        for (int k=0; k<NODES_PER_FACET; ++k)
            if (facets[3*f + k] != a && facets[3*f + k] != b) return facets[3*f + k];
        return -1;
    };
    while (!stack.empty()) {
        const edge_t e = stack.back();
        stack.pop_back();
        auto it = efacets.find(e);
        if (it == efacets.end() || it->second[1] < 0) continue;   // flipped away, or outline
        const int f = it->second[0], g = it->second[1], a = e.first, b = e.second;
        const int c = third(f, a, b), d = third(g, a, b);
        const double oa = orient(c, d, a), ob = orient(c, d, b), oc = orient(a, b, c), od = orient(a, b, d);
        if (oa * ob >= 0 || oc * od >= 0) continue;     // the diagonals do not cross in projection
        if (efacets.count(key(c, d))) continue;          // cd is an edge already (the projection folds)
        if (!(std::min(min_angle(c, d, a), min_angle(c, d, b)) >
              std::min(min_angle(a, b, c), min_angle(a, b, d)))) continue;
        const double s = oa / (oa - ob), t = oc / (oc - od);   // the crossing on ab and on cd
        const int M = major_direction;
        const double gap = std::abs(P(a,M) + s*(P(b,M) - P(a,M)) - P(c,M) - t*(P(d,M) - P(c,M)));
        const double moved_fg = std::max(moved[f], moved[g]) + gap;
        if (moved_fg > max_change) continue;
        moved[f] = moved[g] = moved_fg;
        const int nf3[3] = {a, c, d}, ng3[3] = {b, d, c};
        std::copy(nf3, nf3 + 3, &facets[3*f]);
        std::copy(ng3, ng3 + 3, &facets[3*g]);
        efacets.erase(it);
        efacets[key(c, d)] = {{f, g}};
        for (auto &x : efacets[key(b, c)]) if (x == f) x = g;   // (b,c) moved from f to g
        for (auto &x : efacets[key(a, d)]) if (x == g) x = f;   // (a,d) moved from g to f
        stack.push_back(key(a, c)); stack.push_back(key(b, c));
        stack.push_back(key(a, d)); stack.push_back(key(b, d));
    }
}
#endif


// Delete points_to_delete and re-triangulate only the cavity of each boundary that loses a point:
// the facets touching a deleted point, plus the facets they enclose (the fill has no holes). The
// fill is a constrained triangulation of the cavity's remaining points in the boundary's
// projection plane; every other facet and boundary is kept. A deleted point on the outline (a
// ridge with the next boundary) drops out of it, so its two outline edges become one chord on both
// boundaries and the surface stays closed. The fill is not bounded; flip_boundary_facets then
// improves the boundary's facets, moving its surface by at most max_change in all.
void delete_points_and_merge_facets(const int_vec &points_to_delete,
                                    uint deleting_bdry,
                                    const array_t &bnormals,
                                    int &npoints,
                                    int &nseg, double *points,
                                    int *segment, int *segflag,
                                    uint_vec &bcflag, double max_change)
{
#ifndef THREED
    die(EXIT_INTERNAL_ASSERT, "delete_points_and_merge_facets() doesn't work in 2D!");
#else
    std::vector<char> deleting(npoints, 0);
    for (int p : points_to_delete) deleting[p] = 1;

    std::vector<char> replaced(nseg, 0);   // facet replaced by its boundary's cavity triangulation
    int_vec new_facets, new_flags;         // NODES_PER_FACET nodes per new facet

    for (int i=0; i<nbdrytypes; ++i) {
        const int flag = 1 << i;
        if (!(deleting_bdry & flag)) continue;

        // this boundary's facets, and the facet owning each directed edge; across() below needs a
        // boundary's facets consistently oriented, as TetGen outputs the facets of a mesh
        int_vec bfacets;
        std::unordered_map<edge_t, int, hash1, equal_to1> owner;
        for (int f=0; f<nseg; ++f) {
            if (segflag[f] != flag) continue;
            bfacets.push_back(f);
            const int *s = segment + f*NODES_PER_FACET;
            for (int k=0; k<NODES_PER_FACET; ++k)
                if (!owner.emplace(std::make_pair(s[k], s[(k+1) % NODES_PER_FACET]), f).second)
                    die(EXIT_INTERNAL_ASSERT, "delete_points_and_merge_facets: two facets of a "
                        "boundary share a directed edge (inconsistent facet orientation)");
        }
        auto across = [&](int a, int b) {   // the facet on the other side of edge a->b, or -1
            auto it = owner.find(std::make_pair(b, a));
            return it == owner.end() ? -1 : it->second;
        };

        // the cavity, then every non-cavity facet not linked through non-cavity facets to the
        // boundary's outline (a hole) joins it
        for (int f : bfacets) {
            const int *s = segment + f*NODES_PER_FACET;
            if (deleting[s[0]] || deleting[s[1]] || deleting[s[2]]) replaced[f] = 1;
        }
        std::vector<char> outside(nseg, 0);
        int_vec stack;
        for (int f : bfacets) {
            if (replaced[f]) continue;
            const int *s = segment + f*NODES_PER_FACET;
            for (int k=0; k<NODES_PER_FACET; ++k)
                if (across(s[k], s[(k+1) % NODES_PER_FACET]) < 0 && !outside[f]) {
                    outside[f] = 1;
                    stack.push_back(f);
                }
        }
        while (!stack.empty()) {
            const int f = stack.back();
            stack.pop_back();
            const int *s = segment + f*NODES_PER_FACET;
            for (int k=0; k<NODES_PER_FACET; ++k) {
                const int g = across(s[k], s[(k+1) % NODES_PER_FACET]);
                if (g >= 0 && !replaced[g] && !outside[g]) {
                    outside[g] = 1;
                    stack.push_back(g);
                }
            }
        }
        for (int f : bfacets)
            if (!outside[f]) replaced[f] = 1;

        // outline successor of each point on the boundary's outline (edges with no facet across)
        std::unordered_map<int, int> outline_next;
        for (int f : bfacets) {
            const int *s = segment + f*NODES_PER_FACET;
            for (int k=0; k<NODES_PER_FACET; ++k) {
                const int a = s[k], b = s[(k+1) % NODES_PER_FACET];
                if (across(a, b) < 0 && !outline_next.emplace(a, b).second)
                    die(EXIT_MESH_QUALITY, "delete_points_and_merge_facets: the boundary outline "
                        "is not a simple polygon");
            }
        }

        // the cavity's remaining points (local numbering) and its outline as segments, a deleted
        // outline point bridged by the chord between its remaining outline neighbours
        int_vec inv;                          // local -> global point
        std::unordered_map<int, int> local;   // global -> local point
        auto local_of = [&](int n) {
            auto it = local.emplace(n, (int)inv.size());
            if (it.second) inv.push_back(n);
            return it.first->second;
        };
        int_vec surf_segment;
        for (int f : bfacets) {
            if (!replaced[f]) continue;
            const int *s = segment + f*NODES_PER_FACET;
            for (int k=0; k<NODES_PER_FACET; ++k) {
                const int a = s[k];
                int b = s[(k+1) % NODES_PER_FACET];
                if (deleting[a]) continue;
                local_of(a);
                const int g = across(a, b);
                if (g >= 0 && replaced[g]) continue;          // inside the cavity
                while (deleting[b]) b = outline_next.at(b);   // only outline points can be deleted here
                surf_segment.push_back(local_of(a));
                surf_segment.push_back(local_of(b));
            }
        }
        // a cavity left with fewer than 3 points is a facet whose bridged chord is one of its own
        // edges: dropping the facet closes it
        if (inv.size() < 3) continue;
        const int nsurfseg = surf_segment.size() / 2;
        int_vec surf_segflag(nsurfseg);   // all 0, its value does not matter

        // 3d coordinate -> 2d, dropping the normal's largest component
        int major_direction = i / 2;
        if (i >= iboundn0) { // slant boundaries start at iboundn0
            double max = 0;
            for (int d=0; d<NDIMS; d++) {
                if (std::abs(bnormals[i][d]) > max) {
                    max = std::abs(bnormals[i][d]);
                    major_direction = d;
                }
            }
        }
        double_vec coord2d(2 * inv.size());
        for (std::size_t j=0; j<inv.size(); ++j)
            for (int d=0, dd=0; d<NDIMS; d++)
                if (d != major_direction) coord2d[2*j + dd++] = points[inv[j]*NDIMS + d];

        // re-triangulate the cavity (connect only, not adding new points)
        int new_nnode, new_nelem, new_nseg;
        double *pcoord, *pregattr;
        int *pconnectivity, *psegment, *psegflag;
        Mesh mesh;
        mesh.min_angle = 0;
        mesh.meshing_verbosity = 0;
        points_to_new_surface(mesh, inv.size(), coord2d.data(),
                              nsurfseg, surf_segment.data(), surf_segflag.data(),
                              0, NULL,
                              0, 3,
                              new_nnode, new_nelem, new_nseg,
                              pcoord, pconnectivity, psegment, psegflag, pregattr);
        if (static_cast<std::size_t>(new_nnode) != inv.size())
            die(EXIT_MESH_QUALITY, "delete_points_and_merge_facets: the cavity triangulation "
                "added or dropped points");
        if (new_nseg != nsurfseg)
            die(EXIT_MESH_QUALITY, "delete_points_and_merge_facets: the cavity triangulation "
                "split a segment");
        // the bridged outline must keep the replaced facets' orientation and enclose the fill
        // exactly: a chord passing beyond a remaining point inverts part of it, and the fill then
        // overlaps the kept facets (areas doubled, in coord2d's plane)
        {
            int d0 = -1, d1 = -1;
            for (int d=0; d<NDIMS; d++)
                if (d != major_direction) (d0 < 0 ? d0 : d1) = d;
            double old_area = 0, outline_area = 0, fill_area = 0;
            for (int f : bfacets) {
                if (!replaced[f]) continue;
                const int *n = segment + f*NODES_PER_FACET;
                auto X = [&](int k, int d) { return points[n[k]*NDIMS + d]; };
                old_area += (X(1,d0) - X(0,d0)) * (X(2,d1) - X(0,d1)) - (X(1,d1) - X(0,d1)) * (X(2,d0) - X(0,d0));
            }
            for (int j=0; j<nsurfseg; ++j) {
                const double *p = &coord2d[2*surf_segment[2*j]], *q = &coord2d[2*surf_segment[2*j + 1]];
                outline_area += p[0]*q[1] - p[1]*q[0];
            }
            for (int j=0; j<new_nelem; ++j) {
                const double *p = &coord2d[2*pconnectivity[3*j]], *q = &coord2d[2*pconnectivity[3*j + 1]],
                             *r = &coord2d[2*pconnectivity[3*j + 2]];
                fill_area += std::abs((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]));
            }
            // 1e-9: far above the round-off of these sums (~1e-13 relative), far below one facet
            if (outline_area * old_area <= 0 || std::abs(std::abs(outline_area) - fill_area) > 1e-9 * fill_area)
                die(EXIT_MESH_QUALITY, "delete_points_and_merge_facets: the cavity fill overlaps the "
                    "boundary (a deleted outline point's chord passes beyond a remaining point)");
        }
        // the boundary's facets after the fill, re-emitted with the new ones after bounded flips
        int_vec bnew;
        for (int f : bfacets) {
            if (replaced[f]) continue;
            bnew.insert(bnew.end(), segment + f*NODES_PER_FACET, segment + (f+1)*NODES_PER_FACET);
            replaced[f] = 1;
        }
        for (int j=0; j<new_nelem*NODES_PER_FACET; ++j)
            bnew.push_back(inv[pconnectivity[j]]);
        flip_boundary_facets(bnew, points, major_direction, max_change);
        new_facets.insert(new_facets.end(), bnew.begin(), bnew.end());
        new_flags.insert(new_flags.end(), bnew.size() / NODES_PER_FACET, flag);
        delete [] pcoord;
        delete [] pconnectivity;
        delete [] psegment;
        delete [] psegflag;
        delete [] pregattr;
    }

    // kept facets first, then the cavities' new ones; the rest marked deleted. Internal facets
    // (no boundary bit) are dropped: the deletion below renumbers only boundary points.
    int_vec kept_facets, kept_flags;
    for (int f=0; f<nseg; ++f) {
        if (replaced[f] || !(segflag[f] & BOUND_ANY)) continue;
        kept_facets.insert(kept_facets.end(), segment + f*NODES_PER_FACET, segment + (f+1)*NODES_PER_FACET);
        kept_flags.push_back(segflag[f]);
    }
    kept_facets.insert(kept_facets.end(), new_facets.begin(), new_facets.end());
    kept_flags.insert(kept_flags.end(), new_flags.begin(), new_flags.end());
    const int nseg2 = kept_flags.size();
    if (nseg2 > nseg)
        die(EXIT_INTERNAL_ASSERT, "delete_points_and_merge_facets: more facets than before");
    for (int n : kept_facets)
        if (deleting[n])
            die(EXIT_INTERNAL_ASSERT, "delete_points_and_merge_facets: a facet keeps a deleted point");
    std::copy(kept_facets.begin(), kept_facets.end(), segment);
    std::copy(kept_flags.begin(), kept_flags.end(), segflag);
    for (int i=nseg2; i<nseg; ++i) {
        for (int j=0; j<NODES_PER_FACET; ++j)
            segment[i*NODES_PER_FACET + j] = DELETED_FACET;
        segflag[i] = 0;
    }
    nseg = nseg2;

    // delete points from the end
    int *endsegment = segment + nseg * NODES_PER_FACET;  // last segment
    int end = npoints - 1;  // last point
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
        if (DEBUG) {
            std::cout << " deleting point " << *i << " replaced by point " << end << '\n';
        }

        // when a point is deleted, replace it with the last point
        uint flag = bcflag[end];
        bcflag[*i] = flag;
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        if (is_boundary(flag)) {
            std::replace(segment, endsegment, end, *i);
            if (DEBUG > 1) {
                std::cout << *i << " <- " << end << "\n";
                std::cout << "segment after replace: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        end --;
        npoints --;
    }

#endif
}


// Delete the points in points_to_delete, merging the neighbouring boundary
// segments (2D) / facets (3D) for the points that lie on a boundary so the
// boundary stays closed; in 2D a boundary point whose two segments are both longer
// than min_size is kept, in 3D the facet flips move the surface by at most min_size in all.
// Non-boundary points are simply removed. This is the
// deletion path for remeshing_option 1 and 10-12 (boundaries may be modified); it
// fully handles the removal, so the caller must NOT also call delete_points().
void delete_points_and_merge_boundary(int_vec &points_to_delete,
                               const array_t &bnormals,
                               int &npoints,
                               int &nseg, double *points,
                               int *segment, int *segflag,
                               uint_vec &bcflag, double min_size)
{
    if (DEBUG > 1) {
        std::cout << "old points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
        std::cout << "segment before delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
    }

#ifdef THREED
    // the boundaries that lose a point
    uint deleting_bdry = 0;
    for (int p : points_to_delete) deleting_bdry |= bcflag[p] & BOUND_ANY;

    // non-boundary points changed,
    if (! deleting_bdry) {
        delete_points(points_to_delete, npoints, nseg,
                      points, segment);
        delete_facets(nseg, segment, segflag);
        return;
    }

    delete_points_and_merge_facets(points_to_delete, deleting_bdry, bnormals, npoints, nseg,
                                   points, segment, segflag, bcflag, min_size);
    delete_facets(nseg, segment, segflag);
#else
    delete_points_and_merge_segments(points_to_delete, npoints, nseg,
                                     points, segment, bcflag, min_size);
    delete_facets(nseg, segment, segflag);
#endif

    if (DEBUG > 1) {
        std::cout << "segment after  delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
    }
}


void refine_surface_elem(const Param &param, const Variables &var,
                         const array_t &old_coord, const conn_t &old_connectivity,
                         const double_vec &old_volume, int &old_nnode, double *qcoord)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    const double surface_vol = param.mesh.sediment_size * sizefactor * std::pow(param.mesh.resolution, NDIMS);

    std::cout << "    Checking surface element volume.\n";

//    #pragma omp parallel for default(none) shared(param, var, old_coord, old_connectivity, old_volume, old_nnode, qcoord)

    for (size_t i=0; i<(*var.surfinfo.top_facet_elems).size(); i++) {
        int e = (*var.surfinfo.top_facet_elems)[i];

        int_vec &a = (*var.elemmarkers)[e];
        if (a[param.mat.mattype_sed] == 0) continue;
//        int mat = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
//        if (mat != param.mat.mattype_sed) continue;

        if (old_volume[e] < surface_vol) continue;

        ConstConnAccessor conn = old_connectivity[e];
        int_vec n(NDIMS);

//        if (DEBUG)
            std::printf("      Surface node added (%4d %.1e %.1e)\n",e, old_volume[e], surface_vol);
        // get the nodes of the element on surface
        for (int j=0; j<NDIMS; j++)
            n[j] = (*var.connectivity_surface)[e][j];

        int nsub_node = -1;
        for (int j=0;j<NODES_PER_ELEM; j++)
            if (std::find(n.begin(),n.end(),conn[j]) == n.end()) {
                nsub_node = conn[j];
                break;
            }

        if (nsub_node >= 0) {
            // for nodes on surface
            for (int j=0; j<NDIMS; j++) {
                double mcoord[NDIMS]    ;

                for (int d=0;d<NDIMS; d++) {
                    mcoord[d] = old_coord[ n[j] ][d];
                    mcoord[d] += old_coord[nsub_node][d];
                    mcoord[d] /= 2.;
                }

//                #pragma omp critical(refine_surface_elem)
                {
                    for (int d=0;d<NDIMS; d++)
                        qcoord[old_nnode*NDIMS + d] = mcoord[d];
                    old_nnode++;
                }
            }
        }
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

// Node-spacing floor of the Triangle/TetGen path (new_mesh): drop INTERIOR nodes within hmin of a
// retained node, so re-triangulating the advected points rarely forms a tiny element (spacing
// bounds edges, not measure: new_mesh's smallest_vol retry catches the rest). Boundary nodes stay
// (flatten_* owns their spacing), and so do kept internal segments' nodes (deleting one misroutes
// its segment). Uniform grid of cell hmin: O(nnode).
// Do not exempt shear-band nodes: that lowered the peak plastic strain and added tiny-element
// remeshes.
void decimate_below_hmin(const Param &param, const array_t &coord, int nnode,
                         const uint_vec &bcflag, int nseg, const int *segment,
                         int_vec &points_to_delete)
{
    const double hmin = refine_floors(param.mesh).hmin;
    if (hmin <= 0.0 || nnode <= 0) return;
    const double inv_h = 1.0/hmin, h2 = hmin*hmin;
    const long long OFF = 1LL<<20, BITS = 21;   // cell index range ~[-1M, 1M) -> fits domains up to ~1M*hmin

    auto cell = [&](int i, long long c[NDIMS]) {
        for (int d=0; d<NDIMS; ++d) c[d] = (long long)std::floor(coord[i][d]*inv_h);
    };
    auto pack = [&](const long long c[NDIMS]) -> long long {
        long long k = c[0] + OFF;
        for (int d=1; d<NDIMS; ++d) k = (k<<BITS) | (c[d] + OFF);
        return k;
    };
    std::unordered_map<long long, int_vec> grid;   // cell -> retained node ids

    auto probe = [&](const long long nb[NDIMS], int i) -> bool {
        auto it = grid.find(pack(nb));
        if (it == grid.end()) return false;
        for (int j : it->second) {
            double d2 = 0;
            for (int d=0; d<NDIMS; ++d) { double dd = coord[i][d]-coord[j][d]; d2 += dd*dd; }
            if (d2 < h2) return true;
        }
        return false;
    };
    auto near_retained = [&](int i) -> bool {
        long long c[NDIMS]; cell(i, c);
        long long nb[NDIMS];
#ifdef THREED
        for (int a=-1;a<=1;++a){ nb[0]=c[0]+a;
          for (int b=-1;b<=1;++b){ nb[1]=c[1]+b;
            for (int e=-1;e<=1;++e){ nb[2]=c[2]+e;
              if (probe(nb, i)) return true;
            }}}
#else
        for (int a=-1;a<=1;++a){ nb[0]=c[0]+a;
          for (int b=-1;b<=1;++b){ nb[1]=c[1]+b;
            if (probe(nb, i)) return true;
          }}
#endif
        return false;
    };
    auto retain = [&](int i) { long long c[NDIMS]; cell(i, c); grid[pack(c)].push_back(i); };
    // Pass 1: retain all boundary and segment nodes (never decimated) so interior nodes crowding
    // them are removed. With discarded internal segments every segment node is a boundary node.
    std::vector<char> kept(nnode, 0);
    for (int i=0; i<nnode; ++i) kept[i] = is_boundary(bcflag[i]);
    for (int k=0; k<nseg*NODES_PER_FACET; ++k)
        if (segment[k] != DELETED_FACET) kept[segment[k]] = 1;   // new_bottom (option 2) deletes
    for (int i=0; i<nnode; ++i) if (kept[i]) retain(i);
    // Pass 2: greedily thin the other nodes to >= hmin spacing.
    for (int i=0; i<nnode; ++i) {
        if (kept[i]) continue;
        if (near_retained(i)) points_to_delete.push_back(i);
        else retain(i);
    }
}


void new_mesh(const Param &param, Variables &var, int bad_quality,
              const array_t &original_coord, const conn_t &original_connectivity,
              const segment_t &original_segment, const segflag_t &original_segflag)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // create a copy of original mesh
    array_t old_coord(original_coord);
    conn_t old_connectivity(original_connectivity);
    segment_t old_segment(original_segment);
    segflag_t old_segflag(original_segflag);

    double_vec qcoord_vec;
    int_vec qconn_vec, qsegment_vec, qsegflag_vec;

    old_coord.pack_to(qcoord_vec);
    old_connectivity.pack_to(qconn_vec);
    old_segment.pack_to(qsegment_vec);
    old_segflag.pack_to(qsegflag_vec);

    // raw pointers to old mesh
    double *qcoord = qcoord_vec.data();
    int *qconn = qconn_vec.data();
    int *qsegment = qsegment_vec.data();
    int *qsegflag = qsegflag_vec.data();

    // size of old mesh
    int old_nnode = old_coord.size();
    int old_nelem = old_connectivity.size();
    int old_nseg = old_segment.size();

    // copying useful arrays of old mesh
    double_vec old_volume(*var.volume);
    uint_vec old_bcflag(*var.bcflag);
    bool (*excl_func)(uint) = NULL; // function pointer indicating which point cannot be deleted
    switch (param.mesh.remeshing_option) {
    case 0:
    case 1:
    case 2:
        // DO NOT change the boundary
        excl_func = &is_boundary;
        break;
    case 10:
    case 11:
    case 12:
        // DO NOT change the corners
        excl_func = &is_corner;
        break;
    case 13:
        // new_mesh is the non-MMG remesher, and validate_parameters rejects option 13 without MMG.
        die(EXIT_INTERNAL_UNREACHABLE, "new_mesh: remeshing_option 13 requires MMG remeshing");
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        die(EXIT_CONFIG_VALUE);
    }

    /* choosing which way to remesh the boundary */
    int_vec points_to_delete;
    const double min_dist = refine_floors(param.mesh).min_dist;
    switch (param.mesh.remeshing_option) {
    case 0:
    case 10:
        // Nothing
        break;
    case 1:
    case 11:
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 2:
        new_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
        break;
    case 12:
        flatten_x0_corner(old_bcflag, qcoord, points_to_delete);
        break;
    }

    array_t array_t_qcoord(qcoord, old_nnode);

    if (bad_quality == 3) { // there is a tiny element
        // Marking (non-boundary) points of small elements, which will be deleted later
        int_vec tiny_elems;
        find_tiny_element(param, old_volume, tiny_elems);

        if (tiny_elems.size() > 0) {
            find_points_of_tiny_elem(old_coord, old_connectivity, old_volume,
                                     tiny_elems, old_nnode, array_t_qcoord, old_bcflag, points_to_delete,
                                     excl_func);
        }
    }

    decimate_below_hmin(param, array_t_qcoord, old_nnode, old_bcflag, old_nseg, qsegment,
                        points_to_delete);

    // sort points_to_delete and remove duplicates
    {
        std::sort(points_to_delete.begin(), points_to_delete.end());
        auto last = std::unique(points_to_delete.begin(), points_to_delete.end());
        points_to_delete.resize(last - points_to_delete.begin());
    }

    // Record every ORIGINAL node deleted here (interior ones included) or moved by flatten_*:
    // barycentric_node_interpolation silences its "interior node not found" warning for a new
    // node whose nearest old node is one of them. Before delete_points() renumbers qcoord.
    if ((int)var.remesh_affected_old_node.size() == old_nnode) {
        for (int n : points_to_delete)
            if (n >= 0 && n < old_nnode) var.remesh_affected_old_node[n] = 1;
        for (int n = 0; n < old_nnode; ++n)
            for (int d = 0; d < NDIMS; ++d)
                if (qcoord[n*NDIMS + d] != original_coord[n][d]) {
                    var.remesh_affected_old_node[n] = 1;
                    break;
                }
    }

    // delete points
    switch (param.mesh.remeshing_option) {
    case 0:
    case 2:
        // deleting points that are in no segment (option 2: new_bottom replaced the bottom ones)
        delete_points(points_to_delete, old_nnode, old_nseg,
                      qcoord, qsegment);
        delete_facets(old_nseg, qsegment, qsegflag);
        break;
    case 1:
    case 10:
    case 11:
    case 12:
        // deleting points, some of them might be on the boundary (option 1: the side-wall nodes
        // flatten_bottom marks).
        // delete_points_and_merge_boundary() already removes the points in
        // points_to_delete (merging segments for the boundary ones) and reduces
        // old_nnode accordingly. Calling delete_points() again on the same list
        // would be a *second* deletion pass over an already-compacted array with
        // now-stale indices: its swap-delete + std::replace(segment, end, *i)
        // rewrites the boundary segments of the highest-index nodes (typically the
        // most recently added surface/corner nodes), which detaches the top surface
        // from the side wall and collapses the corner into a spurious cliff. Do NOT
        // call delete_points() here.
        delete_points_and_merge_boundary(points_to_delete, *var.bnormals,
                                  old_nnode, old_nseg,
                                  qcoord, qsegment, qsegflag, old_bcflag, min_dist);
        break;
    }

    // refine surface element where volume is too large
#ifdef THREED
    // todo
#else
    if (param.mesh.meshing_sediment) {
        double *nqcoord = new double[(old_nnode + var.surfinfo.top_nodes->size() * 2 ) * NDIMS];
        std::memcpy(nqcoord, qcoord, sizeof(double) * old_nnode * NDIMS);
        qcoord = nqcoord;
        refine_surface_elem(param, var, old_coord, old_connectivity, old_volume, old_nnode, qcoord);
    }
#endif
    int new_nnode, new_nelem, new_nseg;
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

    int nloops = 0;
    Mesh mesh = param.mesh;  // temporary copy
    mesh.poly_filename = ""; // not used, but still need to be init'd
    while (1) {

        if (bad_quality == 3) {
            // lessen the quality constraint so that less new points got inserted to the mesh
            // and less chance of having tiny elements
            mesh.min_angle *= 0.9;
            mesh.max_ratio *= 1.1;
            mesh.min_tet_angle *= 0.9;
        }
#ifdef THREED
        if (nloops != 0 && bad_quality == 1) {
            // enable tetgen optimization for mesh with higher quality
            // tetgen might consume too much memory and crash!
            mesh.tetgen_optlevel = 3;
        }
#endif
        pregattr = NULL;

        //
        // new mesh
        //

        // We don't want to refine large elements during remeshing,
        // so using negative size as the max area
        const double max_elem_size = -1;
        const int vertex_per_polygon = 3;
        points_to_new_mesh(mesh, old_nnode, qcoord,
                           old_nseg, qsegment, qsegflag,
                           0, pregattr,
                           max_elem_size, vertex_per_polygon,
                           new_nnode, new_nelem, new_nseg,
                           pcoord, pconnectivity, psegment, psegflag, pregattr);

        array_t new_coord(pcoord, new_nnode);
        conn_t new_connectivity(pconnectivity, new_nelem);

        // deleting (non-boundary) nodes to avoid having tiny elements
        double_vec new_volume(new_nelem);
        compute_volume(new_coord, new_connectivity, new_volume);

        #pragma acc wait

        const double smallest_vol = refine_floors(param.mesh).smallest_vol;
        bad_quality = 0;
        for (int e=0; e<new_nelem; e++) {
            if (new_volume[e] < smallest_vol) {
                bad_quality = 3;
                break;
            }
        }
        int worst_elem;
        double q = worst_elem_quality(new_coord, new_connectivity,
                                      new_volume, worst_elem);
#ifdef THREED
        // normalizing q so that its magnitude is about the same in 2D and 3D
        q = std::pow(q, 1.0/3);
#endif
        if (q < param.mesh.min_quality) {
            bad_quality = 1;
        }

        // new_coord.nullify();
        // new_connectivity.nullify();
        if (! bad_quality) break;

        nloops ++;
        if (nloops > 5) {
            std::cout << "Warning: exceeding loop limit in remeshing. Proceeding with risks.\n";
            break;
        }

        delete [] pcoord;
        delete [] pconnectivity;
        delete [] psegment;
        delete [] psegflag;
        delete [] pregattr;
    }

    var.nnode = new_nnode;
    var.nelem = new_nelem;
    var.nseg = new_nseg;
    var.coord->load_from_buffer(pcoord, new_nnode);
    var.connectivity->load_from_buffer(pconnectivity, new_nelem);
    var.segment->load_from_buffer(psegment, var.nseg);
    var.segflag->load_from_buffer(psegflag, var.nseg);

    delete [] pcoord;
    delete [] pconnectivity;
    delete [] psegment;
    delete [] psegflag;

    if (param.mesh.meshing_sediment)
        delete [] qcoord;
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif

}


// use linear interpolation to calculate uniformly distributed 2D curve points
void interpolate_curve(const Param &param,const array_t &input_points, 
                               array_t &output_points, int primary_index) {
    int np = input_points.size();
    int np_out = output_points.size();

    int secondary_index = 1 - primary_index;

    double min_val = input_points[0][primary_index];
    double max_val = input_points[(np - 1)][primary_index];

    bool reverse = false;
    if (min_val > max_val) {
        std::swap(min_val, max_val);
        reverse = true;
    }

    if (primary_index == 1) {
        switch (param.mesh.remeshing_option)
        {
        case 1:
            min_val = -param.mesh.zlength;
            break;
        case 11:
            min_val = -param.mesh.zlength;
            break;
        default:
            break;
        }
    }

    double step = (max_val - min_val) / (np - 1);

    int ind;
    for (int i = 0, j = 0; i < np_out; ++i) {
        if (reverse)
            ind = np_out - 1 - i;
        else
            ind = i;
        ;
        double target = output_points[ind][primary_index];

        while (j < np_out - 2 && target > input_points[(j + 1)][primary_index]) {
            j++;
        }

        double t = (target - input_points[j][primary_index]) /
                   (input_points[j + 1][primary_index] - input_points[j][primary_index]);

        output_points[ind][secondary_index] = 
            (1 - t) * input_points[j][secondary_index] + 
            t * input_points[(j + 1)][secondary_index];
    }
}


// use linear interpolation to calculate uniformly distributed 2D curve points
void interpolate_uniform_curve(const Param &param,const array_t &input_points, 
                               array_t &output_points, int primary_index) {
    int np = input_points.size();

    double min_val = input_points[0][primary_index];
    double max_val = input_points[(np - 1)][primary_index];

    bool reverse = false;
    if (min_val > max_val) {
        std::swap(min_val, max_val);
        reverse = true;
    }
    switch (param.mesh.remeshing_option)
    {
    case 1:
    case 11:
        if (primary_index == (NDIMS - 1))
            min_val = -param.mesh.zlength;
        break;
    case 13:
        if (primary_index == (NDIMS - 1))
            min_val = -param.mesh.zlength;
        else if (primary_index == 0) {
            min_val = 0;
            max_val = param.mesh.xlength;
#ifdef THREED
        } else if (primary_index == 1) {
            min_val = 0;
            max_val = param.mesh.ylength;
#endif
        }
        break;
    default:
        break;
    }

    double step = (max_val - min_val) / (np - 1);

    int ind;
    for (int i = 0, j = 0; i < np; ++i) {
        if (reverse)
            ind = np - 1 - i;
        else
            ind = i;
        double target = min_val + i * step;
        output_points[ind][primary_index] = target;

        while (j < np - 2 && target > input_points[(j + 1)][primary_index]) {
            j++;
        }

        double t = (target - input_points[j][primary_index]) /
                   (input_points[j + 1][primary_index] - input_points[j][primary_index]);

        for (int k = 0; k < NDIMS; k++) {
            if (k != primary_index) {
                output_points[ind][k] = 
                    (1 - t) * input_points[j][k] + 
                    t * input_points[(j + 1)][k];
            }
        }
    }
}


double quadraticInterpolation(double x, 
                              double x0, double x1, double x2, 
                              double y0, double y1, double y2) {
    return y0 * ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2)) +
           y1 * ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) +
           y2 * ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1));
}


void new_equilateral_info(const Param& param, const Variables& var, double xlength, int *nx, int *nz, int *nnode, int *nelem, int *nseg) {
    
    double sqrt3_to_2 = 2./std::sqrt(3.0);

    double x_mid = xlength / 2;
    *nx = int((x_mid-0.5*param.mesh.resolution) / param.mesh.resolution)*2 + 2;
    *nz = int(param.mesh.zlength*sqrt3_to_2 / param.mesh.resolution) + 1;
    *nnode = (*nx)*((int)(((*nz)-1)/2)+1) + ((*nx)+1)*(int(((*nz)-1)/2)+(1-(*nz)%2));
    *nelem = (2*(*nx)-1) * ((*nz)-1);
    *nseg = 2*((*nz)+(*nx)-2) + 1 - (*nz)%2;
}

void get_side_nodes(const Variables& var, const segflag_t& old_segflag, const segment_t &old_segment, int side, int* side_tips) {
    int ind = 0;
    for (int i = 0; i < var.nseg; ++i) {
        if(old_segflag[i][0] == side) {
            if (ind==0) {
                side_tips[ind] = old_segment[i][0];
                ind++;
                side_tips[ind] = old_segment[i][1];
                ind++;
            } else {
                side_tips[ind] = old_segment[i][1];
                ind++;
            }
        }
    }
}


void new_uniformed_equilateral_mesh(const Param &param, Variables &var,
              const array_t &old_coord, const conn_t &old_conn,
              const segment_t &old_segment, const segflag_t &old_segflag)
{
    int nx_new, nz_new, nnode_new, nelem_new, nseg_new;

    // find sides    
    int side_top[var.nx], side_bottom[var.nx+1-var.nz%2], side_left[var.nz], side_right[var.nz];

    get_side_nodes(var, old_segflag, old_segment, BOUNDZ1, side_top);
    get_side_nodes(var, old_segflag, old_segment, BOUNDZ0, side_bottom);
    get_side_nodes(var, old_segflag, old_segment, BOUNDX0, side_left);
    get_side_nodes(var, old_segflag, old_segment, BOUNDX1, side_right);

    double xlength = old_coord[side_top[var.nx-1]][0] - old_coord[side_top[0]][0];

    new_equilateral_info(param, var, xlength, &nx_new, &nz_new, &nnode_new, &nelem_new, &nseg_new);

    double *qcoord = new double[nnode_new* NDIMS];

    array_t inz(var.nz);
    array_t outz(nz_new);

    double dx = param.mesh.resolution;
    double dz = -param.mesh.resolution * std::sqrt(3.0) / 2.;

    var.coord->reset(qcoord, nnode_new);

    int istart_n2 = nx_new * (int(nz_new/2)+nz_new%2);

    // interpolate left side
    for (int j = 0; j < var.nz; ++j) {
        ConstArrayAccessor p = old_coord[side_left[j]];
        inz[j][0] = p[0];
        inz[j][1] = p[1];
    }
    // give assign z value to the new mesh
    double ddz = (old_coord[side_left[0]][1] - 0.)/(nz_new-1);

    for (int j=0; j<nz_new; ++j) {
        if (j == nz_new-1) {
            outz[j][1] = -param.mesh.zlength;
        } else {
            outz[j][1] = j * dz + old_coord[side_left[0]][1] + ddz*j;
        }
    }

    interpolate_curve(param, inz, outz, 1);
    for (int j = 0; j < nz_new; ++j) {
        int inc = j%2;
        int idx;
        if (inc) {
            idx = istart_n2 + int(j/2)*(nx_new+1);
        } else {
            idx = int(j/2)*nx_new;
        }
        ArrayAccessor p = (*var.coord)[idx];
        p[0] = outz[j][0];
        p[1] = outz[j][1];
    }

    // interpolate right side
    for (int j = 0; j < var.nz; ++j) {
        ConstArrayAccessor p = old_coord[side_right[j]];
        inz[j][0] = p[0];
        inz[j][1] = p[1];
    }

    // give assign z value to the new mesh
    ddz = (old_coord[side_right[0]][1] - 0.)/(nz_new-1);
    for (int j=0; j<nz_new; ++j) {
        if (j == nz_new-1) {
            outz[j][1] = -param.mesh.zlength;
        } else {
            outz[j][1] = j * dz + old_coord[side_right[0]][1] + ddz*j;;
        }
    }
    interpolate_curve(param, inz, outz, 1);
    for (int j = 0; j < nz_new; ++j) {
        int inc = j%2;
        int idx;
        if (inc) {
            idx = istart_n2 + (int(j/2)+1)*(nx_new+1)-1;
        } else {
            idx = (int(j/2)+1)*nx_new-1;
        }
        ArrayAccessor p = (*var.coord)[idx];
        p[0] = outz[j][0];
        p[1] = outz[j][1];
    }

    array_t inx_top(var.nx);
    array_t outx_top(nx_new);

    // interpolate top side
    for (int i = 0; i < var.nx; ++i) {
        ConstArrayAccessor p = old_coord[side_top[i]];
        inx_top[i][0] = p[0];
        inx_top[i][1] = p[1];
    }

    double bdy_dx = (xlength - (nx_new-1)*dx) / 2.;
    double bdy_dz = param.mesh.zlength - (nz_new-1)*dz;

    outx_top[0][0] = old_coord[side_top[0]][0];
    for (int i=1; i <nx_new-1; ++i)
        outx_top[i][0] = i * dx + bdy_dx + old_coord[side_top[0]][0];
    outx_top[nx_new-1][0] = old_coord[side_top[var.nx-1]][0];

    interpolate_curve(param, inx_top, outx_top, 0);
    for (int i = 0; i <nx_new; ++i) {
        ArrayAccessor p = (*var.coord)[i];
        p[0] = outx_top[i][0];
        p[1] = outx_top[i][1];
    }

    array_t inx_bot(var.nx+1-var.nz%2);
    array_t outx_bot(nx_new+1-nz_new%2);

    // interpolate top side
    for (int i = 0; i < var.nx+1-var.nz%2; ++i) {
        ConstArrayAccessor p = old_coord[side_bottom[i]];
        inx_bot[i][0] = p[0];
        switch (param.mesh.remeshing_option) {
        case 1:
            inx_bot[i][1] = -param.mesh.zlength;
            break;
        case 11:
            inx_bot[i][1] = -param.mesh.zlength;
            break;
        default:
            inx_bot[i][1] = p[1];
        }
    }
    int nbot = nx_new + 1-nz_new%2;

    outx_bot[0][0] = old_coord[side_bottom[0]][0];
    for (int i=1; i <nx_new-nz_new%2; ++i)
        outx_bot[i][0] = (i+(nz_new%2-1)+0.5*(1-nz_new%2)) * dx + bdy_dx + old_coord[side_bottom[0]][0];
    outx_bot[nx_new-nz_new%2][0] = old_coord[side_bottom[var.nx-var.nz%2]][0];

    interpolate_curve(param, inx_bot, outx_bot, 0);
    int bot_node0;
    if (nz_new % 2 == 0) {
        // last odd row
        bot_node0 = nnode_new - nbot;
    } else {
        // last even row
        bot_node0 = int(nz_new/2) * nx_new;
    }
    for (int i = 0; i < nbot; ++i) {
        ArrayAccessor p = (*var.coord)[bot_node0 + i];
        p[0] = outx_bot[i][0];
        p[1] = outx_bot[i][1];
    }

    // interpolate the x with left and right side nodes
    for (int j = 1; j < int(nz_new/2)+nz_new%2; ++j) {
        double xi = (*var.coord)[j*nx_new][0];
        double xe = (*var.coord)[(j+1)*nx_new-1][0];
        for (int i = 1; i < nx_new-1; ++i)
            (*var.coord)[i + j*nx_new][0] = xi + bdy_dx + i * dx;
    }
    for (int j = 0; j < int(nz_new/2)-1+nz_new%2; ++j) {
        double xi = (*var.coord)[istart_n2+j*(nx_new+1)][0];
        double xe = (*var.coord)[istart_n2+(j+1)*(nx_new+1)-1][0];
        for (int i = 1; i < nx_new; ++i)
            (*var.coord)[istart_n2 + i + j*(nx_new+1)][0] = xi + bdy_dx + (i-0.5) * dx;
    }

    array_t outx_top_fine(nx_new*2-1);
    array_t outx_bot_fine(nx_new*2-1);

    outx_top_fine[0][0] = old_coord[side_top[0]][0];
    for (int i=1; i <nx_new*2-2; ++i)
        outx_top_fine[i][0] = i * dx/2 + bdy_dx + old_coord[side_top[0]][0];
    outx_top_fine[nx_new*2-2][0] = old_coord[side_top[var.nx-1]][0];

    interpolate_curve(param, inx_top, outx_top_fine, 0);

    outx_bot_fine[0][0] = old_coord[side_bottom[0]][0];
    for (int i=1; i <nx_new*2-2; ++i)
        outx_bot_fine[i][0] = i * dx/2. + bdy_dx + old_coord[side_bottom[0]][0];
    outx_bot_fine[nx_new*2-2][0] = old_coord[side_bottom[var.nx-var.nz%2]][0];

    interpolate_curve(param, inx_bot, outx_bot_fine, 0);

    bdy_dz = param.mesh.zlength - (nz_new-2)*dz;
    dz = -param.mesh.resolution * std::sqrt(3.0) / 2.;

    double zi, ze;
    // interpolate the z with bottom and top side nodes
    for (int i = 1; i < nx_new-1; i++) {
        zi = outx_top_fine[i*2][1];
        ze = outx_bot_fine[i*2][1];
        ddz = ((zi - ze) - dz*(nz_new-2)-bdy_dz)/(nz_new-1);
        for (int j = 1; j < int(nz_new/2); j++) {
            (*var.coord)[i + j*nx_new][1] = zi + (j*2)*dz - (j*2)*ddz;
        }
    }

    for (int i = 1; i < nx_new; i++) {
        zi = outx_top_fine[i*2-1][1];
        ze = outx_bot_fine[i*2-1][1];
        ddz = ((zi - ze) - dz*(nz_new-2)-bdy_dz)/(nz_new-1);
        for (int j = 0; j < int(nz_new/2)-1+nz_new%2; j++) {
            (*var.coord)[istart_n2 + i + j*(nx_new+1)][1] = zi + (j*2+1)*dz - (j*2+1)*ddz;
        }
    }

    var.nx = nx_new;
    var.nz = nz_new;
    var.nnode = nnode_new;
    var.nelem = nelem_new;
    var.nseg = nseg_new;

    int *qconn, *qsegment, *qsegflag;

    create_equilateral_elem(var, qconn);
    create_equilateral_segments(var, qsegment, qsegflag);

    var.connectivity->load_from_buffer(qconn, nelem_new);
    var.segment->load_from_buffer(qsegment, nseg_new);
    delete [] qconn;
    delete [] qsegment;

    var.segflag->reset(qsegflag, nseg_new);
}


void new_uniformed_regular_mesh(const Param &param, Variables &var,
              const array_t &old_coord, const conn_t &old_conn,
              const segment_t &old_segment, const segflag_t &old_segflag)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    double *qcoord = new double[var.nnode * NDIMS];
    int *qconn = new int[var.nelem * NODES_PER_ELEM];
    int *qsegment = new int[var.nseg * NODES_PER_FACET];
    int *qsegflag = new int[var.nseg];

    var.coord->reset(qcoord, var.nnode);
    var.connectivity->reset(qconn, var.nelem);
    var.segment->reset(qsegment, var.nseg);
    var.segflag->reset(qsegflag, var.nseg);

    #pragma omp parallel for default(none) shared(var,old_conn)
    for (int i = 0; i < var.nelem; ++i) {
        ConnAccessor  p = (*var.connectivity)[i];
        p[0] = old_conn[i][0];
        p[1] = old_conn[i][1];
        p[2] = old_conn[i][2];
#ifdef THREED
        p[3] = old_conn[i][3];
#endif
    }
    #pragma omp parallel for default(none) shared(var,old_segment,old_segflag)
    for (int i = 0; i < var.nseg; ++i) {
        SegmentAccessor  p = (*var.segment)[i];
        p[0] = old_segment[i][0];
        p[1] = old_segment[i][1];
#ifdef THREED
        p[2] = old_segment[i][2];
#endif
        (*var.segflag)[i][0] = old_segflag[i][0];
    }

    // interpolate coordinates
#ifdef THREED
    int_vec nxyz = {var.nx, var.ny, var.nz};

    // interpolate edges
    for (int n0=0;n0<NDIMS;n0++) {
        for (int n1=n0+1;n1<NDIMS;n1++) {
            if (n0 >= n1) continue;
            int n2 = 3 - n0 - n1;

            #pragma omp parallel for default(none) shared(param,var,old_coord,nxyz,n0,n1,n2) collapse(2)
            for (int ii=0; ii<2; ii++) {
                for (int jj=0; jj<2; jj++) {
                    array_t in(nxyz[n2]);
                    array_t out(nxyz[n2]);
                    int_vec idx(3);
                    idx[n0] = ii*(nxyz[n0]-1);
                    idx[n1] = jj*(nxyz[n1]-1);

                    for (int kk=0; kk<nxyz[n2]; kk++) {
                        idx[n2] = kk;
                        ConstArrayAccessor p = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        in[kk][0] = p[0];
                        in[kk][1] = p[1];
                        in[kk][2] = p[2];
                    }
                    switch (param.mesh.remeshing_option)
                    {
                    case 1:
                    case 11:
                        if (idx[2] == 0 && n2 != 2)
                            for (int kk=0; kk<nxyz[n2]; ++kk)
                                in[kk][2] = -param.mesh.zlength;
                        break;
                    case 13:
                        if (idx[2] == 0 && n2 != 2)
                            for (int kk=0; kk<nxyz[n2]; ++kk)
                                in[kk][2] = -param.mesh.zlength;

                        if (idx[0] == 0 && n2 != 0)
                            for (int kk=0; kk<nxyz[n2]; ++kk)
                                in[kk][0] = 0;

                        if (idx[0] == nxyz[0]-1 && n2 != 0)
                            for (int kk=0; kk<nxyz[n2]; ++kk)
                                in[kk][0] = param.mesh.xlength;

                        if (idx[1] == 0 && n2 != 1)
                            for (int kk=0; kk<nxyz[n2]; ++kk)
                                in[kk][1] = 0;

                        if (idx[1] == nxyz[1]-1 && n2 != 1)
                            for (int kk=0; kk<nxyz[n2]; ++kk)
                                in[kk][1] = param.mesh.ylength;

                        break;
                    default:
                        break;
                    }
    
                    interpolate_uniform_curve(param, in, out, n2);

                    for (int kk=0; kk<nxyz[n2]; kk++) {
                        idx[n2] = kk;
                        ArrayAccessor p = (*var.coord)[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        p[0] = out[kk][0];
                        p[1] = out[kk][1];
                        p[2] = out[kk][2];
                    }
                }
            }
        }
    }

    // interpolation x, y, z in each plane
    for (int n0=0; n0<NDIMS; n0++) {
        #pragma omp parallel for default(none) shared(param,var,old_coord,nxyz,n0)
        for (int ii=1; ii<nxyz[n0]-1; ii++) {
            for (int n1=0; n1<NDIMS; n1++) {
                if (n0 == n1) continue;
                int n2 = 3 - n0 - n1;
                int i0 = 0;
                int i1 = nxyz[n1]-1;
                int_vec idx(3);
                idx[n0] = ii;
                for (int jj=0; jj<2; jj++) {
                    idx[n2] = jj * (nxyz[n2]-1);
                    idx[n1] = i0;
                    int ind0 = idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2];
                    idx[n1] = i1;
                    int ind1 = idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2];
                    double delta = ((*var.coord)[ind1][n1] - (*var.coord)[ind0][n1]) / i1;
                    for (int kk=1; kk<nxyz[n1]-1; kk++) {
                        idx[n1] = kk;
                        int ind = idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2];
                        (*var.coord)[ind][n1] = (*var.coord)[ind0][n1] + kk * delta;
                    }
                }
            }
        }   
    }

    // 2D interpolation of x-y, y-z, x-z plane
    for (int n0=0; n0<NDIMS; n0++) {
        for (int n1=0; n1<NDIMS; n1++) {
            if (n0 >= n1) continue;
            int n2 = 3 - n0 - n1;

            #pragma omp parallel for default(none) shared(param,var,old_coord,nxyz,n0,n1,n2) collapse(2)
            for (int ii=1; ii<nxyz[n0]-1;ii++) {
                for (int jj=1; jj<nxyz[n1]-1;jj++) {
                    for (int kk=0; kk<2; kk++) {
                        int_vec idx(3);
                        idx[n0] = ii;
                        idx[n1] = jj;
                        idx[n2] = kk * (nxyz[n2]-1);
                        int ind = idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2];


                        switch (param.mesh.remeshing_option)
                        {
                        case 1:
                        case 11:
                            if (n2==2 && kk==0) {
                                (*var.coord)[ind][2] = -param.mesh.zlength;
                                continue;
                            }
                            break;
                        case 13:
                            if (n2==2 && kk==0) {
                                (*var.coord)[ind][2] = -param.mesh.zlength;
                                continue;
                            } else if (n2==0 && kk==0) {
                                (*var.coord)[ind][0] = 0;
                                continue;
                            } else if (n2==0 && kk==1) {
                                (*var.coord)[ind][0] = param.mesh.xlength;
                                continue;
                            } else if (n2==1 && kk==0) {
                                (*var.coord)[ind][1] = 0;
                                continue;
                            } else if (n2==1 && kk==1) {
                                (*var.coord)[ind][1] = param.mesh.ylength;
                                continue;
                            }
                            break;
                        default:
                            break;
                        }

                        // Grid indices of the 3-point stencil bracketing p0 and p1. All -1
                        // until a search below finds the first grid coordinate above its
                        // target; a target at or above the top of the grid never matches,
                        // and the stencil is then the last three points.
                        int ind_x0=-1, ind_x1=-1, ind_x2=-1, ind_y0=-1, ind_y1=-1, ind_y2=-1;
                        double p0 = (*var.coord)[ind][n0];
                        double p1 = (*var.coord)[ind][n1];

                        for (int iii=0; iii<nxyz[n0]; iii++) {
                            idx[n0] = iii;
                            if ((*var.coord)[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]][n0] > p0) {
                                ind_x0 = iii - 1;
                                ind_x1 = iii;
                                ind_x2 = iii + 1;
                                break;
                            }
                        }
                        // ind_x2 < 0 is the no-match case; it takes the same last-three
                        // stencil as a match whose upper point ran off the end.
                        if (ind_x2 < 0 || ind_x2 > nxyz[n0] - 1) {
                            ind_x2 = nxyz[n0] - 1;
                            ind_x1 = nxyz[n0] - 2;
                            ind_x0 = nxyz[n0] - 3;
                        }
                        if (ind_x0 < 0) {
                            ind_x2 = 2;
                            ind_x1 = 1;
                            ind_x0 = 0;
                        }
                        if (ind_x2 > nxyz[n0] - 1) {
                            ind_x2 = nxyz[n0] - 1;
                        }

                        idx[n0] = ii;
                        for (int jjj=0; jjj<nxyz[n1]; jjj++) {
                            idx[n1] = jjj;
                            if ((*var.coord)[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]][n1] > p1) {
                                ind_y0 = jjj - 1;
                                ind_y1 = jjj;
                                ind_y2 = jjj + 1;
                                break;
                            }
                        }
                        if (ind_y2 < 0 || ind_y2 > nxyz[n1] - 1) {
                            ind_y2 = nxyz[n1] - 1;
                            ind_y1 = nxyz[n1] - 2;
                            ind_y0 = nxyz[n1] - 3;
                        }
                        if (ind_y0 < 0) {
                            ind_y2 = 2;
                            ind_y1 = 1;
                            ind_y0 = 0;
                        }
                        if (ind_y2 > nxyz[n1] - 1) {
                            ind_y2 = nxyz[n1] - 1;
                        }
                        idx[n0] = ind_x0;
                        idx[n1] = ind_y0;
                        ConstArrayAccessor x0y0 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n1] = ind_y1;
                        ConstArrayAccessor x0y1 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n1] = ind_y2;
                        ConstArrayAccessor x0y2 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n0] = ind_x1;
                        idx[n1] = ind_y0;
                        ConstArrayAccessor x1y0 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n1] = ind_y1;
                        ConstArrayAccessor x1y1 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n1] = ind_y2;
                        ConstArrayAccessor x1y2 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n0] = ind_x2;
                        idx[n1] = ind_y0;
                        ConstArrayAccessor x2y0 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n1] = ind_y1;
                        ConstArrayAccessor x2y1 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];
                        idx[n1] = ind_y2;
                        ConstArrayAccessor x2y2 = old_coord[idx[0]*nxyz[1]*nxyz[2] + idx[1]*nxyz[2] + idx[2]];

                        double interp0 = quadraticInterpolation(p1,
                            x0y0[n1], x0y1[n1], x0y2[n1], x0y0[n2], x0y1[n2], x0y2[n2]);
                        double interp1 = quadraticInterpolation(p1, 
                            x1y0[n1], x1y1[n1], x1y2[n1], x1y0[n2], x1y1[n2], x1y2[n2]);
                        double interp2 = quadraticInterpolation(p1, 
                            x2y0[n1], x2y1[n1], x2y2[n1], x2y0[n2], x2y1[n2], x2y2[n2]);
                        (*var.coord)[ind][n2] = quadraticInterpolation(p0, 
                            x0y1[n0], x1y1[n0], x2y1[n0], interp0, interp1, interp2);
                    }
                }
            }
        }
    }

    // interpolate x,y,z inside the mesh
    for (int n0=0; n0<NDIMS;n0++) {
        for (int n1=0; n1<NDIMS;n1++) {
            if (n0 >= n1) continue;
            int n2 = 3 - n0 - n1;

            #pragma omp parallel for default(none) shared(var,nxyz,n0,n1,n2) collapse(2)
            for (int ii=1; ii<nxyz[n0]-1;ii++) {
                for (int jj=1; jj<nxyz[n1]-1; jj++) {
                    int_vec idx0(3);
                    idx0[n0] = ii;
                    idx0[n1] = jj;
                    idx0[n2] = 0;
                    double zs = (*var.coord)[idx0[0] * var.ny * var.nz + idx0[1] * var.nz + idx0[2]][n2];
                    idx0[n2] = nxyz[n2] - 1;
                    double ze = (*var.coord)[idx0[0] * var.ny * var.nz + idx0[1] * var.nz + idx0[2]][n2];
                    double delta = (ze - zs) / (nxyz[n2] - 1);
                    for (int kk=1; kk<nxyz[n2]-1; kk++) {
                        idx0[n2] = kk;
                        (*var.coord)[idx0[0] * var.ny * var.nz + idx0[1] * var.nz + idx0[2]][n2] = zs + kk * delta;
                    }
                }
            }
        }
    }

#else
    array_t inz(var.nz);
    array_t outz(var.nz);
    array_t inx(var.nx);
    array_t outx(var.nx);
    // interpolate left side
    for (int j = 0; j < var.nz; ++j) {
        ConstArrayAccessor p = old_coord[j];
        inz[j][1] = p[1];

        switch (param.mesh.remeshing_option) {
        case 13:
            inz[j][0] = 0;
            break;
        default:
            inz[j][0] = p[0];
        }
    }
    interpolate_uniform_curve(param, inz, outz, 1);
    for (int j = 0; j < var.nz; ++j) {
        ArrayAccessor p = (*var.coord)[j];
        p[0] = outz[j][0];
        p[1] = outz[j][1];
    }

    // interpolate right side
    for (int j = 0; j < var.nz; ++j) {
        ConstArrayAccessor p = old_coord[(var.nx - 1) * var.nz + j];
        inz[j][1] = p[1];

        switch (param.mesh.remeshing_option) {
        case 13:
            inz[j][0] = param.mesh.xlength;
            break;
        default:
            inz[j][0] = p[0];
        }
    }
    interpolate_uniform_curve(param, inz, outz, 1);
    for (int j = 0; j < var.nz; ++j) {
        ArrayAccessor p = (*var.coord)[(var.nx - 1) * var.nz + j];
        p[0] = outz[j][0];
        p[1] = outz[j][1];
    }
    // interpolate botton side
    for (int i = 0; i < var.nx; ++i) {
        ConstArrayAccessor p = old_coord[var.nz * i];
        inx[i][0] = p[0];
        switch (param.mesh.remeshing_option) {
        case 1:
        case 11:
        case 13:
            inx[i][1] = -param.mesh.zlength;
            break;
        default:
            inx[i][1] = p[1];
        }
    }
    interpolate_uniform_curve(param, inx, outx, 0);
    for (int i = 0; i < var.nx; ++i) {
        ArrayAccessor p = (*var.coord)[var.nz * i];
        p[0] = outx[i][0];
        p[1] = outx[i][1];
    }
    // interpolate top side
    for (int i = 0; i < var.nx; ++i) {
        ConstArrayAccessor p = old_coord[var.nz * (i + 1) - 1];
        inx[i][0] = p[0];
        inx[i][1] = p[1];
    }
    interpolate_uniform_curve(param, inx, outx, 0);
    for (int i = 0; i < var.nx; ++i) {
        ArrayAccessor p = (*var.coord)[var.nz * (i + 1) - 1];
        p[0] = outx[i][0];
        p[1] = outx[i][1];
    }
    // interpolate the x with left and right side nodes
    for (int j = 1; j < var.nz - 1; ++j) {
        double xi = (*var.coord)[j][0];
        double xe = (*var.coord)[(var.nx - 1) * var.nz - 1 + j][0];
        double dx = (xe - xi) / (var.nx - 1);
        for (int i = 1; i < var.nx - 1; ++i)
            (*var.coord)[i * var.nz + j][0] = xi + i * dx;
    }
    // interpolate the z with bottom and top side nodes
    for (int i = 1; i < var.nx - 1; ++i) {
        double zi = (*var.coord)[i * var.nz][1];
        double ze = (*var.coord)[(i + 1) * var.nz - 1][1];
        double dz = (ze - zi) / (var.nz - 1);
        for (int j = 1; j < var.nz - 1; ++j)
            (*var.coord)[i * var.nz + j][1] = zi + j * dz; 
    }
#endif

        // create_uniform_interpolated_mesh(param, var, old_coord, var.coord);
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

#ifdef USEMMG
void compute_metric_field(const Param &param, const Variables &var, double_vec &metric, double_vec &etmp)
{
    /* Compute the desired element size (nodal metric = edge length) for MMG remeshing.
     *
     * The base size is init_elem_size_n (the frozen initial nodal element size), so away
     * from plastic strain element sizes are MAINTAINED across remeshing (metric == base).
     *
     * Where plastic strain has accumulated, the target is reduced to refine the mesh,
     * following Triangle's area/volume-constraint convention: the target element VOLUME is
     * scaled by 1/(1 + coeff*plstrain), so the target EDGE LENGTH that MMG consumes is that
     * volume ratio raised to 1/NDIMS. (The earlier form scaled the edge length by the
     * volume ratio directly -- i.e. treated a volume ratio as a length ratio -- which was
     * ~NDIMS times too aggressive and refined the whole domain to the floor whenever
     * yielding was broad, ratcheting the element count up at every remesh.)
     *
     * coeff = mmg_metric_refine_coeff makes the refinement sensitivity easy to manage.
     */
    // Not const-qualified on purpose: a const scalar is predetermined-shared in OpenMP and
    // may not appear in a shared() clause on some toolchains (gcc<=10, icpc) -- see the
    // #ifdef GPP1X dance used for sizefactor elsewhere. A plain local sidesteps that.
    double coeff = param.mesh.mmg_metric_refine_coeff;
    std::fill_n(metric.begin(), var.nnode, 0);

    // volume-based refinement target where plastic strain is present (Triangle max_area style)
    #pragma omp parallel for default(none) shared(var, etmp, coeff)
    for (int e = 0; e < var.nelem; e++)
        etmp[e] = (*var.volume)[e] / (1.0 + coeff * (*var.plstrain)[e]);

    #pragma omp parallel for default(none) shared(var, metric, etmp)
    for (int n = 0; n < var.nnode; n++) {
        double vol_ratio = 0.0;
        const int npatch = var.support.size(n);
        const int* patch = var.support.patch(n);
        for (int i=0; i<npatch; ++i)
            vol_ratio += etmp[patch[i]];
        vol_ratio /= (*var.volume_n)[n];
#ifdef THREED
        metric[n] = (*var.init_elem_size_n)[n] * std::cbrt(vol_ratio);
#else
        metric[n] = (*var.init_elem_size_n)[n] * std::sqrt(vol_ratio);
#endif
    }
}


// Signed measure of an element (2x area in 2D, 6x volume in 3D) from packed node coords.
double signed_elem_measure(const double *coord, const int *elem)
{
#ifdef THREED
    const double *a = coord + elem[0]*3, *b = coord + elem[1]*3;
    const double *c = coord + elem[2]*3, *d = coord + elem[3]*3;
    double bx=b[0]-a[0], by=b[1]-a[1], bz=b[2]-a[2];
    double cx=c[0]-a[0], cy=c[1]-a[1], cz=c[2]-a[2];
    double dx=d[0]-a[0], dy=d[1]-a[1], dz=d[2]-a[2];
    return bx*(cy*dz-cz*dy) - by*(cx*dz-cz*dx) + bz*(cx*dy-cy*dx);
#else
    const double *a = coord + elem[0]*2, *b = coord + elem[1]*2, *c = coord + elem[2]*2;
    return (b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1]);
#endif
}

// Orientation + scale reference for a flat/inverted-element test: the measure of the first element
// `pristine` accepts that is wound like the majority. Not simply the first: flatten_* can turn an
// untouched element inside out, which would invert the whole test.
template <typename Pristine, typename Measure>
static double majority_reference(int nelem, const Pristine &pristine, const Measure &measure)
{
    double first_pos = 0.0, first_neg = 0.0;
    int npos = 0, nneg = 0;
    for (int e = 0; e < nelem; ++e) {
        if (!pristine(e)) continue;
        const double m = measure(e);
        if (m > 0.0) { if (npos++ == 0) first_pos = m; }
        else if (m < 0.0) { if (nneg++ == 0) first_neg = m; }
    }
    return (npos >= nneg) ? first_pos : first_neg;
}

// Rebuild a collapsed mesh's boundary segments to match its surviving elements (nconn, new ids):
// drop a segment collapsed to a point and a flagged one that is no longer the face of exactly one
// survivor. A face a dropped flat/inverted element exposed is new boundary, on the one boundary its
// nodes share (left unflagged, MMG re-adds it with ref 0: a hole in the flagged boundary). A face
// on no single boundary means the dropped element was a fold, not a sliver on a wall: stop.
template <typename Remap>
static void rebuild_collapsed_segments(int nseg, const int *segment, const int *segflag,
                                       const Remap &remap, int ne, const int_vec &nconn,
                                       const int_vec &new_to_old_node, const uint_vec &bcflag,
                                       const double_vec &ncoord,
                                       int_vec &nsegment, int_vec &nsegflag, int &ns)
{
    // faces of the surviving elements: sorted node set -> {owner count, the owner's node order}
    typedef std::array<int, NODES_PER_FACET> Facet;
    std::map<Facet, std::pair<int, Facet>> faces;
    for (int e = 0; e < ne; ++e)
        for (int f = 0; f < FACETS_PER_ELEM; ++f) {
            Facet v, key;
            for (int k = 0; k < NODES_PER_FACET; ++k) v[k] = nconn[e*NODES_PER_ELEM + NODE_OF_FACET[f][k]];
            key = v;
            std::sort(key.begin(), key.end());
            auto &fc = faces[key];
            ++fc.first;
            fc.second = v;
        }

    nsegment.clear(); nsegflag.clear();
    ns = 0;
    int n_orphan = 0, n_exposed = 0;
    for (int s = 0; s < nseg; ++s) {
        int v[NODES_PER_FACET];
        bool degenerate = false;
        for (int k = 0; k < NODES_PER_FACET; ++k) {
            v[k] = remap(segment[s*NODES_PER_FACET + k]);
            for (int j = 0; j < k; ++j) if (v[j] == v[k]) degenerate = true;
        }
        if (degenerate) continue;
        if (segflag[s] != 0) {
            Facet key;
            std::copy(v, v + NODES_PER_FACET, key.begin());
            std::sort(key.begin(), key.end());
            auto it = faces.find(key);
            if (it == faces.end() || it->second.first != 1) { ++n_orphan; continue; }
            it->second.first = 0;   // covered
        }
        for (int k = 0; k < NODES_PER_FACET; ++k) nsegment.push_back(v[k]);
        nsegflag.push_back(segflag[s]);
        ++ns;
    }

    for (const auto &kv : faces) {
        if (kv.second.first != 1) continue;
        uint flag = BOUND_ANY;
        for (int n : kv.first) flag &= bcflag[new_to_old_node[n]];
        const uint wall = flag & ~BOUNDZ1;   // every node on a wall's top edge: the face lies in the
        if ((flag & BOUNDZ1) && wall && !(wall & (wall - 1))) flag = wall;   // wall (the top is no plane)
        if (flag == 0 || (flag & (flag - 1)) != 0) {
            std::cerr << "Error: a collapse before MMG left the mesh folded at (";
            for (int d = 0; d < NDIMS; ++d) std::cerr << (d ? ", " : "") << ncoord[kv.first[0]*NDIMS + d];
            std::cerr << "): a face exposed by a dropped flat/inverted element lies on no single"
                         " boundary.\n";
            die(EXIT_MESH_QUALITY);
        }
        for (int n : kv.second.second) nsegment.push_back(n);
        nsegflag.push_back(flag);
        ++ns;
        ++n_exposed;
    }
    if (n_orphan || n_exposed)
        std::cout << "      [collapse] boundary: dropped " << n_orphan << " orphan segment(s), flagged "
                  << n_exposed << " exposed face(s).\n";
}

// MMG-native handling of material that moved OUTSIDE a restored boundary. flatten_* snapped the
// boundary nodes back onto their planes and collected the nodes past, or within min_dist of, a
// restored plane (`pts`, called outside nodes here); those past it would invert the adjoining
// elements. Remove them by edge collapse instead of falling back to a whole-domain
// Triangle/TetGen remesh:
//   * CORNER-AWARE: each node flatten_* marked records the planes it is past or within min_dist
//     of (BOUND* bitmask `crossed`: the restored plane it was marked for and any other plane near
//     it) and only collapses onto a boundary node of one of them, never onto a node of the top
//     surface alone.
//   * ITERATIVE: each pass merges only the outside nodes currently touching such a boundary node,
//     peeling deep stacks one layer at a time (mapping a deep node straight to the boundary tangled).
//   * STRANDED: a last pass merges any node left without a surviving element, outside or not,
//     onto its nearest well-connected node (a boundary node: one on a shared plane), whatever
//     its mask.
// Builds fresh 0-indexed output arrays; inputs untouched. Survivor node order is preserved so the
// metric maps by renumbering; segments are rebuilt to match the surviving elements.
void collapse_outside_nodes(const int_vec &pts, int nnode, int nelem, int nseg,
                            const double *coord, const uint_vec &bcflag,
                            const int *conn, const int *segment, const int *segflag,
                            const double_vec &metric,
                            double zlen, double xlen, double ylen, double min_dist,
                            double_vec &ncoord, int_vec &nconn,
                            int_vec &nsegment, int_vec &nsegflag, double_vec &nmetric,
                            int &nn, int &ne, int &ns)
{
    // crossed-plane bitmask per node (BOUND* convention, so a boundary node m lies on a plane
    // that s crossed iff (bcflag[m] & crossed[s]) != 0; a corner node matches on either bit).
    std::vector<uint> crossed(nnode, 0);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        int s = pts[i];
        const double *p = coord + s*NDIMS;
        uint c = 0;
        if (p[0]       <  0.0  + min_dist) c |= BOUNDX0;
        if (p[0]       >  xlen - min_dist) c |= BOUNDX1;
#ifdef THREED
        if (p[1]       <  0.0  + min_dist) c |= BOUNDY0;
        if (p[1]       >  ylen - min_dist) c |= BOUNDY1;
#endif
        if (p[NDIMS-1] < -zlen + min_dist) c |= BOUNDZ0;
        if (!c)   // flatten_* marks a node only past a restored plane, whose test is one of these
            die(EXIT_INTERNAL_ASSERT, "collapse_outside_nodes: marked node crossed no plane");
        crossed[s] = c;
    }

    // merged_to[n] = the node n was collapsed onto, or -1 if alive. Peel and stalemate never merge
    // onto an outside node (merges apply after the scan: a cycle could close); the stranded pass can
    // merge a node others were merged onto, so resolve() follows the chain to its alive end.
    int_vec merged_to(nnode, -1);
    auto alive   = [&](int n) { return merged_to[n] < 0; };
    auto resolve = [&](int n) {
        for (int hops = 0; merged_to[n] >= 0; n = merged_to[n])
            if (++hops > nnode) die(EXIT_INTERNAL_ASSERT, "collapse_outside_nodes: merge cycle");
        return n;
    };

    auto dist2 = [&](int s, int m) {
        double d2 = 0.0;
        for (int d = 0; d < NDIMS; ++d) { double dd = coord[s*NDIMS+d] - coord[m*NDIMS+d]; d2 += dd*dd; }
        return d2;
    };
    // an element no outside and no merged node touches: its shape is the pre-collapse one
    auto pristine = [&](int e) {
        const int *el = conn + e*NODES_PER_ELEM;
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            if (crossed[el[k]] || !alive(el[k])) return false;
        return true;
    };

    // A merge must not invert an element: the nearest crossed-plane node can lie beyond the element's
    // far face (a wall corner above the subsided free surface), and the rebuild would drop it and
    // leave a fold. An element the snap already inverted does not count. An outside node's elements
    // are its original ones (no node merges onto an outside node).
    int_vec slot(nnode, -1);
    for (std::size_t i = 0; i < pts.size(); ++i) slot[pts[i]] = (int)i;
    std::vector<int_vec> elems_of(pts.size());
    for (int e = 0; e < nelem; ++e)
        for (int k = 0; k < NODES_PER_ELEM; ++k) {
            const int i = slot[conn[e*NODES_PER_ELEM + k]];
            if (i >= 0) elems_of[i].push_back(e);
        }
    const double refp = majority_reference(nelem, pristine, [&](int e) {
        return signed_elem_measure(coord, conn + e*NODES_PER_ELEM); });
    std::vector<char> passed_over(pts.size(), 0);   // the nearest target would have inverted
    auto keeps_orientation = [&](int s, int m) {   // merging s onto m inverts none of s's elements
        for (int e : elems_of[slot[s]]) {
            int r[NODES_PER_ELEM], q[NODES_PER_ELEM];
            bool degen = false;
            for (int k = 0; k < NODES_PER_ELEM; ++k) {
                r[k] = resolve(conn[e*NODES_PER_ELEM + k]);
                q[k] = (r[k] == s) ? m : r[k];
                for (int j = 0; j < k; ++j) if (q[j] == q[k]) degen = true;
            }
            if (degen) continue;   // collapses with the merged edge, as intended
            // refuse a real inversion only: a merge that leaves the element flat passes, and the
            // rebuild's flat cut (m/ref <= 1e-6) then drops it, as the collapse intends
            if (signed_elem_measure(coord, r) * refp > 0.0 &&
                signed_elem_measure(coord, q) * refp < -1e-6 * refp * refp)
                return false;
        }
        return true;
    };

    // PEEL: each pass merges every outside node touching a boundary node of a crossed plane onto the
    // nearest such node that keeps its elements oriented (one the re-check below rejects waits a
    // pass), one layer at a time: mapping a deep node straight to the boundary tangles the mesh.
    for (int pass = 0; pass < nnode; ++pass) {
        int_vec    tgt(nnode, -1);
        double_vec best(nnode, std::numeric_limits<double>::max()), nearest(best);
        for (int e = 0; e < nelem; ++e) {
            const int *el = conn + e*NODES_PER_ELEM;
            for (int a = 0; a < NODES_PER_ELEM; ++a) {
                int s = resolve(el[a]);
                if (!crossed[s]) continue;                       // only still-outside nodes
                for (int b = 0; b < NODES_PER_ELEM; ++b) {
                    int m = resolve(el[b]);
                    if (m == s || crossed[m]) continue;          // never onto an outside node
                    if (is_boundary(bcflag[m]) && (bcflag[m] & crossed[s])) {  // on a crossed plane
                        double d2 = dist2(s, m);
                        nearest[s] = std::min(nearest[s], d2);
                        if (d2 < best[s] && keeps_orientation(s, m)) { best[s] = d2; tgt[s] = m; }
                    }
                }
            }
        }
        // applied one by one, each re-checked against this pass's earlier merges: two outside nodes
        // of one element can each keep it oriented alone and invert it together
        bool progress = false;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            int s = pts[i];
            if (!alive(s) || !crossed[s]) continue;
            if (tgt[s] >= 0 && !keeps_orientation(s, tgt[s])) { tgt[s] = -1; passed_over[slot[s]] = 1; }
            if (nearest[s] < best[s]) passed_over[slot[s]] = 1;
            if (tgt[s] >= 0) { merged_to[s] = tgt[s]; progress = true; }
        }
        if (!progress) break;
    }

    // STALEMATE: an outside node the peel could not reach merges onto the nearest orientation-keeping
    // boundary node of a crossed plane, searched globally (rare); sweep again while a sweep merges (a
    // merge can make a target valid for a node passed before). With none it stays (logged).
    for (bool again = true; again; ) {
        again = false;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            int s = pts[i];
            if (!alive(s) || !crossed[s]) continue;
            double best = std::numeric_limits<double>::max(), nearest = best; int tg = -1;
            for (int m = 0; m < nnode; ++m) {
                if (!alive(m) || crossed[m] || !is_boundary(bcflag[m]) || !(bcflag[m] & crossed[s]))
                    continue;
                double d2 = dist2(s, m);
                nearest = std::min(nearest, d2);
                if (d2 < best && keeps_orientation(s, m)) { best = d2; tg = m; }
            }
            if (nearest < best) passed_over[slot[s]] = 1;
            if (tg >= 0) { merged_to[s] = tg; again = true; }
        }
    }
    int n_passed = 0, n_left = 0;
    for (char c : passed_over) n_passed += c;
    for (int n = 0; n < nnode; ++n) n_left += (crossed[n] && alive(n));
    if (n_passed)
        std::cout << "      [collapse] " << n_passed << " outside node(s) deferred or redirected a merge "
                     "that would invert an element; " << n_left
                  << " left unmerged by the crossed-plane passes.\n";

    // STRANDED: a node alive but with no surviving element, outside or not. A sunk blob collapsing
    // onto its plane can leave a boundary node whose elements all degenerate while its segments
    // survive; MMG would keep it as an isolated 0-mass vertex -> NaN. Merge each onto its nearest
    // well-connected node, ignoring `crossed` (the dangling segment collapses with it); iterate,
    // since absorbing one node can strand the next.
    {
        // Reference measure (orientation + scale), computed on `coord` like the connectivity
        // rebuild's so the flat/inverted test matches exactly.
        const double refm = majority_reference(nelem, pristine, [&](int e) {
            return signed_elem_measure(coord, conn + e*NODES_PER_ELEM); });
        for (int pass = 0; pass < nnode; ++pass) {
            int_vec live_deg(nnode, 0);
            std::vector<char> in_mesh(nnode, 0);
            for (int e = 0; e < nelem; ++e) {
                const int *el = conn + e*NODES_PER_ELEM;
                int r[NODES_PER_ELEM];
                bool degen = false;
                for (int k = 0; k < NODES_PER_ELEM; ++k) {
                    r[k] = resolve(el[k]);
                    in_mesh[el[k]] = 1;
                    for (int j = 0; j < k; ++j) if (r[j] == r[k]) degen = true;
                }
                if (degen) continue;
                double m = signed_elem_measure(coord, r);
                // same flat/inverted cut as the rebuild below (reason there)
                if (refm != 0.0 && m*refm <= 1e-6*refm*refm) continue;   // flat/inverted -> dropped
                for (int k = 0; k < NODES_PER_ELEM; ++k) live_deg[r[k]]++;
            }
            bool progress = false;
            for (int s = 0; s < nnode; ++s) {
                if (!alive(s) || !in_mesh[s] || live_deg[s] > 0) continue;   // stranded node
                const uint sflag = bcflag[s];
                double best = std::numeric_limits<double>::max(); int tg = -1;
                for (int m = 0; m < nnode; ++m) {
                    if (m == s || !alive(m) || live_deg[m] == 0) continue;   // target must be well-connected
                    // keep boundary continuity: a boundary node merges only onto a node sharing a plane.
                    if (is_boundary(sflag) && !(is_boundary(bcflag[m]) && (bcflag[m] & sflag))) continue;
                    double d2 = dist2(s, m);
                    if (d2 < best) { best = d2; tg = m; }
                }
                if (tg >= 0) { merged_to[s] = tg; progress = true; }
            }
            if (!progress) break;
        }
    }

    // renumber survivors, rebuild coord + metric.
    int_vec node_map(nnode, -1);
    ncoord.clear(); nmetric.clear();
    nn = 0;
    for (int n = 0; n < nnode; ++n) {
        if (!alive(n)) continue;
        node_map[n] = nn++;
        for (int d = 0; d < NDIMS; ++d) ncoord.push_back(coord[n*NDIMS + d]);
        nmetric.push_back(metric[n]);
    }
    auto remap = [&](int n) { return node_map[resolve(n)]; };

    // reference orientation/scale from the pristine elements, after the stranded pass
    const double ref = majority_reference(nelem, pristine, [&](int e) {
        int v[NODES_PER_ELEM];
        for (int k = 0; k < NODES_PER_ELEM; ++k) v[k] = node_map[conn[e*NODES_PER_ELEM + k]];
        return signed_elem_measure(ncoord.data(), v); });

    // rebuild connectivity, dropping elements that collapse onto a boundary: a repeated node,
    // OR near-null / inverted measure (all nodes land on the flat boundary, or orientation flip).
    nconn.clear();
    ne = 0;
    for (int e = 0; e < nelem; ++e) {
        const int *el = conn + e*NODES_PER_ELEM;
        int v[NODES_PER_ELEM];
        bool degenerate = false;
        for (int k = 0; k < NODES_PER_ELEM; ++k) {
            v[k] = remap(el[k]);
            for (int j = 0; j < k; ++j) if (v[j] == v[k]) degenerate = true;
        }
        if (degenerate) continue;
        double m = signed_elem_measure(ncoord.data(), v);
        // Flat/inverted: m/ref <= 1e-6 (ref = a pristine element's measure). An element flattened
        // onto an axis-aligned restored plane measures exactly 0, and 1e-6 stays below the
        // tiny-element trigger up to largest_size/smallest_size ~1e5: no valid element is dropped.
        if (ref != 0.0 && m * ref <= 1e-6 * ref * ref) continue;
        for (int k = 0; k < NODES_PER_ELEM; ++k) nconn.push_back(v[k]);
        ++ne;
    }

    int_vec old_id(nn);   // a survivor's old node id
    for (int n = 0; n < nnode; ++n)
        if (node_map[n] >= 0) old_id[node_map[n]] = n;
    rebuild_collapsed_segments(nseg, segment, segflag, remap, ne, nconn, old_id, bcflag, ncoord,
                               nsegment, nsegflag, ns);
}

// True if any node lies beyond a boundary plane this remeshing option snaps back (bottom for
// 1/11/13; all sides for 13): material outside the fixed domain, to be collapsed before MMG.
// Mirrors what flatten_* marks for deletion: non-bottom nodes below the base, and for a side every
// node NOT on that side's own plane that crossed it (adjacent-boundary nodes included). Excluding a
// plane's own wall nodes keeps genuine extension inert. Option 12 is excluded; option 2 never
// reaches MMG.
bool has_outside_material(const Param &param, const array_t &coord, const uint_vec &bcflag,
                          int nnode, double min_dist)
{
    const int opt = param.mesh.remeshing_option;
    const bool restore_bottom = (opt==1 || opt==11 || opt==13);
    const bool restore_sides  = (opt==13);
    if (!restore_bottom && !restore_sides) return false;
    const double zb = -param.mesh.zlength;
    for (int i = 0; i < nnode; ++i) {
        const uint f = bcflag[i];
        // Thresholds MUST match the flatten_* deletion tests exactly (`z < bottom + min_dist`,
        // `x < min_dist`, ...): a stricter test here left a marked node in the mesh, and the
        // element inverted when flatten snapped a neighbour past it.
        if (restore_bottom && !is_bottom(f) && coord[i][NDIMS-1] < zb + min_dist) return true;
        if (restore_sides) {
            if (!is_x0(f) && coord[i][0] < min_dist) return true;
            if (!is_x1(f) && coord[i][0] > param.mesh.xlength - min_dist) return true;
#ifdef THREED
            if (!is_y0(f) && coord[i][1] < min_dist) return true;
            if (!is_y1(f) && coord[i][1] > param.mesh.ylength - min_dist) return true;
#endif
        }
    }
    return false;
}


#ifdef THREED
void optimize_mesh(const Param &param, Variables &var, int bad_quality,
              const array_t &original_coord, const conn_t &original_connectivity,
              const segment_t &original_segment, const segflag_t &original_segflag)
{
    // We don't want to refine large elements during remeshing,
    // so using negative size as the max area
    const double max_elem_size = -1;
    const int vertex_per_polygon = 3;
    const double min_dist = refine_floors(param.mesh).min_dist;
    Mesh mesh_param = param.mesh;
    mesh_param.poly_filename = "";

    // create a copy of original_coord and original_segment
    double_vec qcoord_vec;
    int_vec qconn_vec, qsegment_vec, qconn_to_mmg_vec, qsegment_to_mmg_vec, qsegflag_vec;

    original_coord.pack_to(qcoord_vec);
    original_connectivity.pack_to(qconn_vec);
    original_segment.pack_to(qsegment_vec);
    original_connectivity.pack_to(qconn_to_mmg_vec);
    original_segment.pack_to(qsegment_to_mmg_vec);
    original_segflag.pack_to(qsegflag_vec);

    double *qcoord = qcoord_vec.data();
    int *qconn = qconn_vec.data();
    int *qsegment = qsegment_vec.data();
    int *qconn_to_mmg = qconn_to_mmg_vec.data();
    int *qsegment_to_mmg = qsegment_to_mmg_vec.data();
    int *qsegflag = qsegflag_vec.data();

    int old_nnode = original_coord.size();
    int old_nelem = original_connectivity.size();
    int old_nseg = original_segment.size();

    // copy
    double_vec old_volume(*var.volume);
    uint_vec old_bcflag(*var.bcflag);
    int_vec old_bnodes[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        old_bnodes[i] = *(var.bnodes[i]);
    }

    int_vec points_to_delete;
    bool (*excl_func)(uint) = NULL; // function pointer indicating which point cannot be deleted

    /* choosing which way to remesh the boundary */
    switch (param.mesh.remeshing_option) {
    case 0:
        // DO NOT change the boundary
        excl_func = &is_boundary;
        break;
    case 1:
        excl_func = &is_boundary;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 10:
        excl_func = &is_corner;
        break;
    case 11:
        excl_func = &is_corner;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 13:
        excl_func = &is_corner;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist);
        flatten_x0(old_bcflag, qcoord, points_to_delete, min_dist);
        flatten_x1(old_bcflag, qcoord, param.mesh.xlength, points_to_delete, min_dist);
        flatten_y0(old_bcflag, qcoord, points_to_delete, min_dist);
        flatten_y1(old_bcflag, qcoord, param.mesh.ylength, points_to_delete, min_dist);
        break;
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        die(EXIT_CONFIG_VALUE);
    }

    // --- Prepare the mesh + metric handed to MMG ---------------------------------
    // Compute the nodal metric (target element size) on the current mesh.
    compute_metric_field(param, var, *var.ntmp, *var.etmp);

    // By default MMG adapts the current mesh in place.
    int   mnode = old_nnode, melem = old_nelem, mseg = old_nseg;
    double *mcoord   = qcoord;
    int    *mconn1   = qconn_to_mmg;      // 0-indexed; mmg_adapt 1-indexes internally
    int    *mseg1    = qsegment_to_mmg;   // 0-indexed; mmg_adapt 1-indexes internally
    int    *msegflag = qsegflag;
    double *mmetric  = (*var.ntmp).data();

    // Material outside a restored boundary (sunk below the fixed bottom; with option 13 also past
    // a side wall): flatten_* snaps the boundary back and the elements beyond invert, which MMG
    // cannot adapt. Collapse the marked nodes onto their boundary plane and let MMG adapt
    // the result; no fallback to a whole-domain Triangle/TetGen remesh (has_outside_material).
    double_vec c_coord, c_metric;
    int_vec c_conn, c_seg, c_segflag;
    bool outside_material = has_outside_material(param, original_coord, old_bcflag,
                                                 old_nnode, min_dist);
    if (outside_material) {
        std::cerr << "  Material moved outside a restored boundary; collapsing it back onto "
                     "the boundary before MMG.\n";
        int cn = 0, ce = 0, cs = 0;
        collapse_outside_nodes(points_to_delete, old_nnode, old_nelem, old_nseg,
                               qcoord, old_bcflag, qconn, qsegment, qsegflag, *var.ntmp,
                               param.mesh.zlength, param.mesh.xlength, param.mesh.ylength, min_dist,
                               c_coord, c_conn, c_seg, c_segflag, c_metric, cn, ce, cs);
        mnode = cn; melem = ce; mseg = cs;
        mcoord = c_coord.data();
        mconn1 = c_conn.data();
        mseg1 = c_seg.data();
        msegflag = c_segflag.data();
        mmetric = c_metric.data();
    }

    // --- Adapt the mesh with the shared MMG driver (mmg_utils.cxx) ----------------

    MMGInput  mmg_in = { mnode, melem, mseg, mcoord, mconn1, mseg1, msegflag, mmetric, nullptr, nullptr, false };
    MMGOutput mmg_out;
    mmg_adapt(param.mesh, mmg_in, mmg_out);

    var.nnode = mmg_out.nnode;
    var.nelem = mmg_out.nelem;
    var.nseg  = mmg_out.nseg;
    std::cerr << "New number of vertices:" << var.nnode
              << ", elements:" << var.nelem << ", segments:" << var.nseg << std::endl;

    array_t   new_coord(var.nnode);
    conn_t    new_connectivity(var.nelem);
    segment_t new_segment(var.nseg);
    segflag_t new_segflag(var.nseg);
    new_coord.load_from_buffer(mmg_out.coord.data(), var.nnode);
    new_connectivity.load_from_buffer(mmg_out.conn.data(), var.nelem);
    new_segment.load_from_buffer(mmg_out.seg.data(), var.nseg);
    new_segflag.load_from_buffer(mmg_out.segflag.data(), var.nseg);
    var.coord->steal_ref( new_coord );
    var.connectivity->steal_ref( new_connectivity );
    var.segment->steal_ref( new_segment );
    var.segflag->steal_ref( new_segflag );

    std::cerr << "\nMesh optimization done" <<std::endl;
}

#else

void optimize_mesh_2d(const Param &param, Variables &var, int bad_quality,
              const array_t &original_coord, const conn_t &original_connectivity,
              const segment_t &original_segment, const segflag_t &original_segflag)
{
    // We don't want to refine large elements during remeshing,
    // so using negative size as the max area
    const double max_elem_size = -1;
    const int vertex_per_polygon = 3;
    const double min_dist = refine_floors(param.mesh).min_dist;
    Mesh mesh_param = param.mesh;
    mesh_param.poly_filename = "";

    // create a copy of original_coord and original_segment
    double_vec qcoord_vec;
    int_vec qconn_vec, qsegment_vec, qconn_to_mmg_vec, qsegment_to_mmg_vec, qsegflag_vec;

    original_coord.pack_to(qcoord_vec);
    original_connectivity.pack_to(qconn_vec);
    original_segment.pack_to(qsegment_vec);
    original_connectivity.pack_to(qconn_to_mmg_vec);
    original_segment.pack_to(qsegment_to_mmg_vec);
    original_segflag.pack_to(qsegflag_vec);

    double *qcoord = qcoord_vec.data();
    int *qconn = qconn_vec.data();
    int *qsegment = qsegment_vec.data();
    int *qconn_to_mmg = qconn_to_mmg_vec.data();
    int *qsegment_to_mmg = qsegment_to_mmg_vec.data();
    int *qsegflag = qsegflag_vec.data();

    int old_nnode = original_coord.size();
    int old_nelem = original_connectivity.size();
    int old_nseg = original_segment.size();

    // copy
    double_vec old_volume(*var.volume);
    uint_vec old_bcflag(*var.bcflag);
    int_vec old_bnodes[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        old_bnodes[i] = *(var.bnodes[i]);
    }

    int_vec points_to_delete;
    bool (*excl_func)(uint) = NULL; // function pointer indicating which point cannot be deleted

    /* choosing which way to remesh the boundary */
    switch (param.mesh.remeshing_option) {
    case 0:
        // DO NOT change the boundary
        excl_func = &is_boundary;
        break;
    case 1:
        excl_func = &is_boundary;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 10:
        excl_func = &is_corner;
        break;
    case 11:
        excl_func = &is_corner;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    case 12:
        flatten_x0_corner(old_bcflag, qcoord, points_to_delete);
        break;
    case 13:
        // restore the bottom AND both side walls to their initial planes
        excl_func = &is_corner;
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        flatten_x0(old_bcflag, qcoord, points_to_delete, min_dist);
        flatten_x1(old_bcflag, qcoord, param.mesh.xlength, points_to_delete, min_dist);
        break;
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        die(EXIT_CONFIG_VALUE);
    }

    // --- Prepare the mesh + metric handed to MMG ---------------------------------
    // Compute the nodal metric (target element size) on the current mesh.
    compute_metric_field(param, var, *var.ntmp, *var.etmp);

    // By default MMG adapts the current mesh in place.
    int   mnode = old_nnode, melem = old_nelem, mseg = old_nseg;
    double *mcoord   = qcoord;
    int    *mconn1   = qconn_to_mmg;      // 0-indexed; mmg_adapt 1-indexes internally
    int    *mseg1    = qsegment_to_mmg;   // 0-indexed; mmg_adapt 1-indexes internally
    int    *msegflag = qsegflag;
    double *mmetric  = (*var.ntmp).data();

    // Material outside a restored boundary (sunk below the fixed bottom; with option 13 also past
    // a side wall): flatten_* snaps the boundary back and the elements beyond invert, which MMG
    // cannot adapt. Collapse the marked nodes onto their boundary plane and let MMG adapt
    // the result; no fallback to a whole-domain Triangle/TetGen remesh (has_outside_material).
    double_vec c_coord, c_metric;
    int_vec c_conn, c_seg, c_segflag;
    bool outside_material = has_outside_material(param, original_coord, old_bcflag,
                                                 old_nnode, min_dist);
    if (outside_material) {
        std::cerr << "  Material moved outside a restored boundary; collapsing it back onto "
                     "the boundary before MMG.\n";
        int cn = 0, ce = 0, cs = 0;
        collapse_outside_nodes(points_to_delete, old_nnode, old_nelem, old_nseg,
                               qcoord, old_bcflag, qconn, qsegment, qsegflag, *var.ntmp,
                               param.mesh.zlength, param.mesh.xlength, param.mesh.ylength, min_dist,
                               c_coord, c_conn, c_seg, c_segflag, c_metric, cn, ce, cs);
        mnode = cn; melem = ce; mseg = cs;
        mcoord = c_coord.data();
        mconn1 = c_conn.data();
        mseg1 = c_seg.data();
        msegflag = c_segflag.data();
        mmetric = c_metric.data();
    }

    // --- Adapt the mesh with the shared MMG driver (mmg_utils.cxx) ----------------

    MMGInput  mmg_in = { mnode, melem, mseg, mcoord, mconn1, mseg1, msegflag, mmetric, nullptr, nullptr, false };
    MMGOutput mmg_out;
    mmg_adapt(param.mesh, mmg_in, mmg_out);

    var.nnode = mmg_out.nnode;
    var.nelem = mmg_out.nelem;
    var.nseg  = mmg_out.nseg;
    std::cerr << "New number of vertices:" << var.nnode
              << ", elements:" << var.nelem << ", segments:" << var.nseg << std::endl;

    array_t   new_coord(var.nnode);
    conn_t    new_connectivity(var.nelem);
    segment_t new_segment(var.nseg);
    segflag_t new_segflag(var.nseg);
    new_coord.load_from_buffer(mmg_out.coord.data(), var.nnode);
    new_connectivity.load_from_buffer(mmg_out.conn.data(), var.nelem);
    new_segment.load_from_buffer(mmg_out.seg.data(), var.nseg);
    new_segflag.load_from_buffer(mmg_out.segflag.data(), var.nseg);
    var.coord->steal_ref( new_coord );
    var.connectivity->steal_ref( new_connectivity );
    var.segment->steal_ref( new_segment );
    var.segflag->steal_ref( new_segflag );

    std::cerr << "\nMesh optimization done" <<std::endl;
}
#endif  // end of if THREED
#endif // end of if USEMMG

} // anonymous namespace


void initialize_elem_size_n(const Variables &var, double_vec &init_elem_size_n)
{
    /* Compute and freeze the initial nodal element size distribution.
     * Called once at step 0 to capture the mesh refinement zones defined
     * in the input file. This frozen field is subsequently interpolated
     * to new nodes during remeshing to prevent refinement zones from
     * diffusing away.
     */
    if (init_elem_size_n.size() > 0) return;

#ifndef ACC
#ifdef GPP1X
    // for Apple clang version 17.0.0
    // future version of clang might fix this problem
    #pragma omp parallel for default(none) shared(var, sizefactor)
#else
    #pragma omp parallel for default(none) shared(var)
#endif
#endif
    #pragma acc parallel loop gang vector async
    for (int e = 0; e < var.nelem; e++) {
#ifdef THREED
        double elem_size = std::cbrt((*var.volume)[e] / sizefactor);
#else
        double elem_size = std::sqrt((*var.volume)[e] / sizefactor);
#endif
        (*var.etmp)[e] = elem_size * (*var.volume)[e];
    }

    // The node loop below divides by var.volume_n, which compute_mass fills first and
    // correct_surface_element rewrites at top nodes (a call before compute_mass gives +inf sizes;
    // the size test also catches a fill from a smaller nnode). compute_mass leaves its node kernel
    // async: without the wait this host read sees the pre-kernel zeros.
    #pragma acc wait
    if (var.volume_n->size() != std::size_t(var.nnode) || (*var.volume_n)[0] <= 0.)
        die(EXIT_INTERNAL_ASSERT, "initialize_elem_size_n ran before compute_mass filled volume_n.");

    init_elem_size_n.resize(var.nnode);
    std::fill_n(init_elem_size_n.begin(), var.nnode, 0);

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, init_elem_size_n)
#endif
    #pragma acc parallel loop gang vector async
    for (int n = 0; n < var.nnode; n++) {
        const int npatch = var.support.size(n);
        const int* patch = var.support.patch(n);
        for (int i=0; i<npatch; ++i)
            init_elem_size_n[n] += (*var.etmp)[patch[i]];
        init_elem_size_n[n] /= (*var.volume_n)[n];
    }
}


int bad_mesh_quality(const Param &param, const Variables &var, int &index, double &min_quality)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    /* Check the quality of the mesh, return 0 if the mesh quality (by several
     * measures) is good. Non-zero returned values indicate --
     * 1: an element has bad quality (too acute / narrow / flat).
     * 2: a bottom node has moved too far away from the flat bottom.
     * 3: an element is smaller than (mesh.smallest_size * [volume of a equilateral triangle/tetrahedron
     *    of side = mesh.resolution]).
     */

    // check tiny elements
    double smallest_vol = refine_floors(param.mesh).smallest_vol;
    // Device scans (a host loop would pull volume/coord/bcflag off the device on every quality
    // check); each reports the FIRST offender = the minimum index.
    int tiny = INT_MAX;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, smallest_vol) reduction(min:tiny)
#endif
    #pragma acc parallel loop gang vector reduction(min:tiny) async
    for (int e=0; e<var.nelem; e++)
        if ((*var.volume)[e] < smallest_vol && e < tiny) tiny = e;
    #pragma acc wait
    if (tiny < INT_MAX) {
        {
            const int e = tiny;
            index = e;
            // report location + the nodes' boundary flags, so a chronic re-trigger spot
            // is identifiable (interior vs which boundary it is pinned to)
            double cent[NDIMS] = {0.};
            ConstConnAccessor conn = (*var.connectivity)[e];
            for (int k = 0; k < NODES_PER_ELEM; ++k)
                for (int d = 0; d < NDIMS; ++d)
                    cent[d] += (*var.coord)[conn[k]][d] / NODES_PER_ELEM;
            std::cout << "    The size of element #" << index << " is too small at (";
            for (int d = 0; d < NDIMS; ++d)
                std::cout << (d ? ", " : "") << cent[d];
            std::cout << "), node bcflags:";
            for (int k = 0; k < NODES_PER_ELEM; ++k)
                std::cout << ' ' << (*var.bcflag)[conn[k]];
            std::cout << ".\n";
#ifdef NPROF
            nvtxRangePop();
#endif
            return 3;
        }
    }

    // check if any bottom node is too far away from the bottom depth
    if (param.mesh.remeshing_option == 1 ||
        param.mesh.remeshing_option == 2 ||
        param.mesh.remeshing_option == 11 ||
        param.mesh.remeshing_option == 13) {
        double bottom = - param.mesh.zlength;
        double dist = param.mesh.max_boundary_distortion * param.mesh.resolution;
        int far = INT_MAX;
#ifndef ACC
        #pragma omp parallel for default(none) shared(var, bottom, dist) reduction(min:far)
#endif
        #pragma acc parallel loop gang vector reduction(min:far) async
        for (int i=0; i<var.nnode; ++i) {
            if (is_bottom((*var.bcflag)[i]) && i < far) {
                double z = (*var.coord)[i][NDIMS-1];
                if (std::fabs(z - bottom) > dist) far = i;
            }
        }
        #pragma acc wait
        if (far < INT_MAX) {
            {
                const int i = far;
                {
                    double z = (*var.coord)[i][NDIMS-1];
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the bottom: z = " << z
                              << " (plane z = " << bottom << ", limit +-" << dist << ")\n";
#ifdef NPROF
                    nvtxRangePop();
#endif
                    return 2;
                }
            }
        }
    }
    // check if any side node is too far away from the side
    if (param.mesh.remeshing_option == 13) {
        index = -1;
        double dist = param.mesh.max_boundary_distortion * param.mesh.resolution;
        double xlen = param.mesh.xlength;
#ifdef THREED
        double ylen = param.mesh.ylength;
#endif
        int side = INT_MAX;
#ifndef ACC
#ifdef THREED
        #pragma omp parallel for default(none) shared(var, dist, xlen, ylen) reduction(min:side)
#else
        #pragma omp parallel for default(none) shared(var, dist, xlen) reduction(min:side)
#endif
#endif
        #pragma acc parallel loop gang vector reduction(min:side) async
        for (int i=0; i<var.nnode; ++i) {
            if (i >= side) continue;
            const uint f = (*var.bcflag)[i];
            bool bad = false;
            if (is_x0(f))      bad = std::fabs((*var.coord)[i][0]) > dist;
            else if (is_x1(f)) bad = std::fabs((*var.coord)[i][0] - xlen) > dist;
#ifdef THREED
            else if (is_y0(f)) bad = std::fabs((*var.coord)[i][1]) > dist;
            else if (is_y1(f)) bad = std::fabs((*var.coord)[i][1] - ylen) > dist;
#endif
            if (bad) side = i;
        }
        #pragma acc wait
        for (int i = side; i < var.nnode && i == side; ++i) {   // at most the one offender
            if (is_x0((*var.bcflag)[i])) {
                double x = (*var.coord)[i][0];
                if (std::fabs(x) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the x0 side: x = " << x
                              << " (plane x = 0, limit +-" << dist << ")\n";
                }
            } else if (is_x1((*var.bcflag)[i])) {
                double x = (*var.coord)[i][0];
                if (std::fabs(x - param.mesh.xlength) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the x1 side: x = " << x
                              << " (plane x = " << param.mesh.xlength << ", limit +-" << dist << ")\n";
                }
#ifdef THREED
            } else if (is_y0((*var.bcflag)[i])) {
                double y = (*var.coord)[i][1];
                if (std::fabs(y) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the y0 side: y = " << y
                              << " (plane y = 0, limit +-" << dist << ")\n";
                }
            } else if (is_y1((*var.bcflag)[i])) {
                double y = (*var.coord)[i][1];
                if (std::fabs(y - param.mesh.ylength) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the y1 side: y = " << y
                              << " (plane y = " << param.mesh.ylength << ", limit +-" << dist << ")\n";
                }
#endif
            }
            if (index >= 0) break;
        }
        if (index >= 0) {
#ifdef NPROF
            nvtxRangePop();
#endif
            return 2;
        }
    }

    // check element distortion
    int worst_elem;
    double q = worst_elem_quality(*var.coord, *var.connectivity,
                                  *var.volume, worst_elem);
#ifdef THREED
    // normalizing q so that its magnitude is about the same in 2D and 3D
    q = std::pow(q, 1.0/3);
#endif
    min_quality = q;
    if (q < param.mesh.min_quality) {
        index = worst_elem;
        // same format as the too-small trigger above: centroid + the nodes' boundary flags
        double cent[NDIMS] = {0.};
        ConstConnAccessor conn = (*var.connectivity)[worst_elem];
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            for (int d = 0; d < NDIMS; ++d)
                cent[d] += (*var.coord)[conn[k]][d] / NODES_PER_ELEM;
        std::cout << "    The quality of element #" << worst_elem << " is too low (" << q
                  << ") at (";
        for (int d = 0; d < NDIMS; ++d)
            std::cout << (d ? ", " : "") << cent[d];
        std::cout << "), node bcflags:";
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            std::cout << ' ' << (*var.bcflag)[conn[k]];
        std::cout << ".\n";
#ifdef NPROF
        nvtxRangePop();
#endif
        return 1;
    }
#ifdef NPROF
    nvtxRangePop();
#endif
    return 0;
}

// Fail-fast guard against a TANGLED mesh (a self-collapsed free surface that no remesh can repair):
//   (1) FOLDED element: signed measure opposite to the mesh's dominant sign (not hard-coded CCW;
//       the 2D solver is orientation-agnostic). It can fire only after a uniform mesher
//       (meshing_elem_shape 1 or 2): Triangle and TetGen orient every element alike, and MMG
//       reorients an inverted input element.
//   (2) SELF-INTERSECTING boundary (2D only): two boundary segments sharing no node cross. MMG
//       does not split such a crossing; Triangle splits it at a vertex the pieces share.
// So it can stop only a uniform-mesher or a 2D MMG remesh, never a Triangle, TetGen or 3D MMG one.
// Prints the offending geometry and dies with the mesh error code. `when` labels the call site.
static double tangle_orient2d(const double *a, const double *b, const double *c)
{
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}
void check_mesh_tangle(const char *when, const array_t &coord, const conn_t &connectivity,
                       int nelem, const segment_t &segment, const segflag_t &segflag, int nseg)
{
    // Signed measure of element e (2x area in 2D / 6x volume in 3D) and the tolerance
    // REL * maxedge^NDIMS: round-off in m is ~1e-15 of maxedge^NDIMS, a fold is O(1), so 1e-9
    // (not tuned) ignores noise and leaves anything flatter to the tiny/quality triggers.
    const double REL = 1e-9;
    auto elem_measure = [&](int e, double &m, double &tol) {
        const auto &cn = connectivity[e];
        double maxedge2 = 0.0;
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            for (int j = 0; j < i; ++j) {
                double d2 = 0.0;
                for (int k = 0; k < NDIMS; ++k) { double dd = coord[cn[i]][k] - coord[cn[j]][k]; d2 += dd*dd; }
                if (d2 > maxedge2) maxedge2 = d2;
            }
#ifdef THREED
        double bx=coord[cn[1]][0]-coord[cn[0]][0], by=coord[cn[1]][1]-coord[cn[0]][1], bz=coord[cn[1]][2]-coord[cn[0]][2];
        double cx=coord[cn[2]][0]-coord[cn[0]][0], cy=coord[cn[2]][1]-coord[cn[0]][1], cz=coord[cn[2]][2]-coord[cn[0]][2];
        double dx=coord[cn[3]][0]-coord[cn[0]][0], dy=coord[cn[3]][1]-coord[cn[0]][1], dz=coord[cn[3]][2]-coord[cn[0]][2];
        m = bx*(cy*dz-cz*dy) - by*(cx*dz-cz*dx) + bz*(cx*dy-cy*dx);
        tol = REL * maxedge2 * std::sqrt(maxedge2);   // measure ~ length^3
#else
        double ax=coord[cn[0]][0], ay=coord[cn[0]][1];
        m = (coord[cn[1]][0]-ax)*(coord[cn[2]][1]-ay) - (coord[cn[2]][0]-ax)*(coord[cn[1]][1]-ay);
        tol = REL * maxedge2;                          // measure ~ length^2
#endif
    };

    // (1) folded element -- opposite orientation to the mesh's dominant sign.
    // Pass 1: dominant orientation sign (ignoring near-degenerate elements).
    int npos = 0, nneg = 0;
    for (int e = 0; e < nelem; ++e) {
        double m, tol; elem_measure(e, m, tol);
        if (m >  tol) ++npos;
        else if (m < -tol) ++nneg;
    }
    const double ref = (npos >= nneg) ? 1.0 : -1.0;   // majority winding is the "correct" one
    // Pass 2: flag the first element wound against the majority.
    for (int e = 0; e < nelem; ++e) {
        double m, tol; elem_measure(e, m, tol);
        if (std::fabs(m) > tol && m * ref < 0.0) {
            const auto &cn = connectivity[e];
            std::cerr << "Error: tangled mesh (" << when << "): element #" << e
                      << " is folded (signed measure = " << m << ", opposite the mesh's orientation). ";
            std::cerr << "nodes";
            for (int k=0;k<NODES_PER_ELEM;++k) {
                std::cerr << " (";
                for (int dd=0; dd<NDIMS; ++dd) std::cerr << coord[cn[k]][dd] << (dd<NDIMS-1?",":"");
                std::cerr << ")";
            }
            std::cerr << "\n       The mesh has folded over itself (typically a self-collapsed free "
                         "surface) and can no longer be repaired by remeshing. Stopping.\n";
            die(EXIT_MESH_QUALITY);
        }
    }

#ifndef THREED
    // (2) self-intersecting boundary (2D), O(nseg^2), remesh time only. Read coords through the
    // Array2D accessor (SoA: components are not adjacent). A per-segment tolerance requires each
    // endpoint clearly on one side, so exact collinearity is never a crossing.
    auto pt = [&](int n, double p[2]) { p[0] = coord[n][0]; p[1] = coord[n][1]; };
    auto seg_len2 = [&](const double p[2], const double q[2]) {
        double dx=q[0]-p[0], dy=q[1]-p[1]; return dx*dx+dy*dy;
    };
    for (int s = 0; s < nseg; ++s) {
        const int a0 = segment[s][0], a1 = segment[s][1];
        double p1[2], p2[2]; pt(a0, p1); pt(a1, p2);
        const double tolA = REL * seg_len2(p1, p2);   // scale of orient() about segment A
        for (int t = s + 1; t < nseg; ++t) {
            const int b0 = segment[t][0], b1 = segment[t][1];
            if (a0==b0 || a0==b1 || a1==b0 || a1==b1) continue;   // share a node -> not a crossing
            double p3[2], p4[2]; pt(b0, p3); pt(b1, p4);
            const double tolB = REL * seg_len2(p3, p4);
            double d1 = tangle_orient2d(p3, p4, p1);   // p1,p2 sides of segment B
            double d2 = tangle_orient2d(p3, p4, p2);
            double d3 = tangle_orient2d(p1, p2, p3);   // p3,p4 sides of segment A
            double d4 = tangle_orient2d(p1, p2, p4);
            bool oppB = (d1 >  tolB && d2 < -tolB) || (d1 < -tolB && d2 >  tolB);
            bool oppA = (d3 >  tolA && d4 < -tolA) || (d3 < -tolA && d4 >  tolA);
            if (oppA && oppB) {
                std::cerr << "Error: tangled mesh (" << when << "): boundary segments #" << s
                          << " [flag " << segflag[s][0] << "] (" << p1[0] << "," << p1[1] << ")-("
                          << p2[0] << "," << p2[1] << ") and #" << t << " [flag " << segflag[t][0]
                          << "] (" << p3[0] << "," << p3[1] << ")-(" << p4[0] << "," << p4[1]
                          << ") cross.\n       The free surface has self-collapsed (self-intersecting "
                             "boundary); the mesh can no longer be repaired by remeshing. Stopping.\n";
                die(EXIT_MESH_QUALITY);
            }
        }
    }
#endif
}


void remesh(const Param &param, Variables &var, int bad_quality)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    int64_t time_tmp = get_nanoseconds();

    std::cout << "  Remeshing starts (#" << var.nremesh_prior + var.nremesh + 1
              << ", step " << var.steps << ", t = " << var.time / YEAR2SEC << " yr)...\n";
#ifdef ACC
    {
        size_t free_bytes, total_bytes;
        knn_bvh_mem_info(&free_bytes, &total_bytes);
        std::cout << "      [GPU mem] remesh start: free="
                  << free_bytes/(1<<20) << " MB / total=" << total_bytes/(1<<20) << " MB\n";
    }
#endif

    double_vec old_surface_area(var.surfinfo.etop);

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,old_surface_area)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0; i<var.surfinfo.etop; ++i) {
        ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity_surface)[i]);
        old_surface_area[i] = compute_area_facet(coord);
    }

    // convert value field to average field
#ifndef ACC
    #pragma omp parallel for default(none) shared(var,old_surface_area)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0; i<var.surfinfo.etop; i++) {
        double inv_volume = 1.0 / old_surface_area[i];
        (*var.surfinfo.edvacc_surf)[i] *= inv_volume;
    }

    // convert volume_old to dv = volume/volume_old - 1 for NN interpolation.
#ifndef ACC
    #pragma omp parallel for default(none) shared(var)
#endif
    #pragma acc parallel loop gang vector async
    for (int e = 0; e < var.nelem; ++e)
        (*var.volume_old)[e] = (*var.volume)[e] / (*var.volume_old)[e] - 1.0;

    // Superconvergent patch recovery for stress before remeshing -- but only where
    // an element can actually take the SPR average. A rheology with no viscous
    // component never relaxes stress, so its Maxwell time is infinite, De = inf and
    // every element keeps the NN stress: the recovery, its nodal transfer across the
    // remesh, the free-surface pin and the Deborah blend are then all dead work.
    // Leaving these pointers null is what switches them off; each consumer tests the
    // one it uses. Note visc() would not have revealed this -- it answers for every
    // rheology, clamping the viscous law into [min_viscosity, max_viscosity], so
    // those two knobs would have picked the remap operator for a run that otherwise
    // never touches them.
    const bool spr_stress = (var.mat->rheol_type & MatProps::rh_viscous) != 0;
    if (spr_stress) {
        var.stress_n = new tensor_t(var.nnode);
        if (param.mat.is_plane_strain)
            var.stressyy_n = new double_vec(var.nnode);
        var.spr_blend_weight = new double_vec(var.nelem);
    }
    var.spr_p_ref_old = new double_vec(var.nelem);

    {
        SurfaceTopo topo_old;
        topo_old.build(param, var);   // old mesh: coord/bnodes still pre-remesh here
        if (spr_stress)
            compute_spr_blend_weight(param, var);   // before centering: visc() reads the trace
        center_stress_to_ref(param, var, topo_old);
        if (spr_stress)
            spr_elem_to_node(param, var, var.stress_n, var.stressyy_n);
    }

    {
        // creating a "copy" of mesh pointer so that they are not deleted
        array_t old_coord;
        conn_t old_connectivity;
        conn_t old_connectivity_surface;
        segment_t old_segment;
        segflag_t old_segflag;
        old_coord.steal_ref(*var.coord);
        old_connectivity.steal_ref(*var.connectivity);
        old_connectivity_surface.steal_ref(*var.connectivity_surface);
        old_segment.steal_ref(*var.segment);
        old_segflag.steal_ref(*var.segflag);

        const int old_nnode = old_coord.size();
        const int old_nelem = old_connectivity.size();
        const int old_nseg  = old_segment.size();
        const char *remesh_engine = "?";

        // Reset the per-remesh record of old nodes the remesh moved or deleted; new_mesh
        // (Triangle/TetGen) fills it for barycentric_node_interpolation.
        var.remesh_affected_old_node.assign(old_nnode, 0);

#ifdef THREED
        if (param.mesh.meshing_elem_shape == 0) {
#if defined USEMMG
            remesh_engine = "MMG";
            optimize_mesh(param, var, bad_quality, old_coord, old_connectivity,
                    old_segment, old_segflag);
#else
            remesh_engine = "TetGen";
            new_mesh(param, var, bad_quality, old_coord, old_connectivity,
                    old_segment, old_segflag);
#endif
        } else if (param.mesh.meshing_elem_shape == 1) {
            remesh_engine = "uniform regular";
            new_uniformed_regular_mesh(param, var, old_coord, old_connectivity,
                    old_segment, old_segflag);
        } else {
            std::cerr << "Error: unknown meshing_elem_shape: " << param.mesh.meshing_elem_shape << '\n';
            die(EXIT_CONFIG_VALUE);
        }
#else  // if 2d
        if (param.mesh.meshing_elem_shape == 0) {
#if defined USEMMG
            remesh_engine = "MMG";
            optimize_mesh_2d(param, var, bad_quality, old_coord, old_connectivity,
                old_segment, old_segflag);
#else
            remesh_engine = "Triangle";
            new_mesh(param, var, bad_quality, old_coord, old_connectivity,
                old_segment, old_segflag);
#endif
        } else if (param.mesh.meshing_elem_shape == 1) {
            remesh_engine = "uniform regular";
            new_uniformed_regular_mesh(param, var, old_coord, old_connectivity,
                old_segment, old_segflag);
        } else if (param.mesh.meshing_elem_shape == 2) {
            remesh_engine = "uniform equilateral";
            new_uniformed_equilateral_mesh(param, var, old_coord, old_connectivity,
                old_segment, old_segflag);
        } else {
            std::cerr << "Error: unknown meshing_elem_shape: " << param.mesh.meshing_elem_shape << '\n';
            die(EXIT_CONFIG_VALUE);
        }        
#endif

        // Per-engine remesh report, with the free-surface (BOUNDZ1) facet count.
        int old_ntop_facets = 0;
        for (int i=0; i<old_nseg; ++i)
            if (static_cast<uint>(old_segflag[i][0]) & BOUNDZ1) ++old_ntop_facets;
        int new_ntop_facets = 0;
        for (int i=0; i<var.nseg; ++i)
            if (static_cast<uint>((*var.segflag)[i][0]) & BOUNDZ1) ++new_ntop_facets;
        std::cout << "    New mesh (" << remesh_engine << "): nodes " << old_nnode << " -> " << var.nnode
                  << ", elements " << old_nelem << " -> " << var.nelem
                  << ", top surface facets " << old_ntop_facets << " -> " << new_ntop_facets << "\n";
        // Fail fast on a tangled mesh before interpolating onto it.
        check_mesh_tangle("post-remesh", *var.coord, *var.connectivity, var.nelem,
                          *var.segment, *var.segflag, var.nseg);
        // Drain all async GPU work before freeing/reallocating temporary arrays.
        // This gives CUDA a sync point to reclaim migrated managed-memory pages.
        #pragma acc wait
#ifdef ACC
        {
            size_t free_bytes, total_bytes;
            knn_bvh_mem_info(&free_bytes, &total_bytes);
            std::cout << "      [GPU mem] before reallocate_tmp: free="
                      << free_bytes/(1<<20) << " MB / total=" << total_bytes/(1<<20) << " MB\n";
        }
#endif
        reallocate_tmp(param, var);

        // Per-NEW-element is_changed mapping, filled by the element NN pass
        // and consumed by spr_node_to_elem; freed with the other stress-remap
        // transients below.
        var.remesh_is_changed = new int_vec(var.nelem);

        if (param.mesh.meshing_elem_shape == 0) {
            // renumbering mesh
            renumbering_mesh(param, *var.coord, *var.connectivity, *var.segment, nullptr);
        }

        create_boundary_flags(var);
        for (int i=0; i<nbdrytypes; ++i)
            var.bfacets[i]->clear();
        delete var.connectivity_surface;
        create_boundary_facets(var);

        {
            Barycentric_transformation bary(old_coord, old_connectivity, *var.volume);

            // interpolating fields defined on elements
            nearest_neighbor_interpolation(param, var, bary, old_coord, old_connectivity);

            Barycentric_transformation bary_surface(old_coord, old_connectivity_surface, old_surface_area, true);

            // interpolating fields defined on surface elements
            nearest_neighbor_interpolation(param, var, bary_surface, old_coord, old_connectivity_surface, true);

            // interpolating fields defined on nodes
            barycentric_node_interpolation(param, var, bary, old_coord, old_connectivity);
        }

        // done with the old mesh; drop the affected-old-node record so it is never
        // consulted against a mismatched mesh later (e.g. temperature-file interpolation).
        var.remesh_affected_old_node.clear();
        create_support(var);
        // delete var.neighbor;
        // delete var.contact;
        // delete var.ctmp;
        // create_neighbor(var);

        // remap markers. elemmarkers and markers_in_elem are updated here, too.
        remap_markers(param, var, old_coord, old_connectivity);
  
        // old_coord et al. are destroyed before exiting this block
    }

    // Drain GPU queue before freeing old field arrays and reallocating new ones.
    // Interpolation functions above may have launched async GPU work; syncing here
    // gives CUDA a chance to reclaim managed-memory pages before new allocations.
    #pragma acc wait
#ifdef ACC
    {
        size_t free_bytes, total_bytes;
        knn_bvh_mem_info(&free_bytes, &total_bytes);
        std::cout << "      [GPU mem] before reallocate_variables: free="
                  << free_bytes/(1<<20) << " MB / total=" << total_bytes/(1<<20) << " MB\n";
    }
#endif

    // memory for new fields
    reallocate_variables(param, var);

    // updating other arrays
    for (int i=0; i<nbdrytypes; ++i)
        var.bnodes[i]->clear();
    create_boundary_nodes(var);

    delete var.top_elems;
    create_top_elems(var);

    update_surface_info(var, var.surfinfo);

    {
        SurfaceTopo topo_new;
        topo_new.build(param, var);   // new mesh: create_boundary_nodes already ran
        if (spr_stress)
            spr_node_to_elem(param, var, topo_new, var.stress, var.stressyy);
        restore_stress_from_ref(param, var, topo_new, var.stress, var.stressyy);
    }

    delete var.stress_n;
    var.stress_n = nullptr;
    delete var.stressyy_n;
    var.stressyy_n = nullptr;
    delete var.spr_blend_weight;
    var.spr_blend_weight = nullptr;
    delete var.spr_p_ref_old;
    var.spr_p_ref_old = nullptr;
    delete var.remesh_is_changed;
    var.remesh_is_changed = nullptr;

    // Timescale reference for the next remesh's Deborah-number stress blend.
    var.last_remesh_time = var.time;

    compute_volume(*var.coord, *var.connectivity, *var.volume);

    double_vec surface_area(var.surfinfo.etop);

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,surface_area)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0; i<var.surfinfo.etop; ++i) {
        ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity_surface)[i]);
        surface_area[i] = compute_area_facet(coord);
    }

    // convert value field back to portional field for nn interpolation
#ifndef ACC
    #pragma omp parallel for default(none) shared(var,surface_area)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0; i<var.surfinfo.etop; i++) {
        (*var.surfinfo.edvacc_surf)[i] *= surface_area[i];
    }

    // convert dv back to actual old volume using the new mesh volumes.
#ifndef ACC
    #pragma omp parallel for default(none) shared(var)
#endif
    #pragma acc parallel loop gang vector async
    for (int e=0; e<var.nelem; ++e)
        (*var.volume_old)[e] = (*var.volume)[e] / (1.0 + (*var.volume_old)[e]);

    // Refresh dt on the NEW mesh unconditionally. dt is otherwise only recomputed every
    // slow_updates_interval (10) steps in the main loop, so a remesh that refines the mesh
    // (smaller minl -> smaller stable dt) would run up to 10 steps on the stale, too-large
    // dt -- enough for an explicit CFL runaway to invert elements (then the step-10
    // compute_dt hits negative volumes and dies with "dt <= 0"). compute_mass below
    // consumes var.dt, so the refresh must come first.
    var.dt = compute_dt(param, var);
    compute_mass(param, var, var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass, *var.hmass, *var.ymass, *var.tmp_result);


    #pragma acc wait

#ifdef NPROF_DETAIL
    nvtxRangePush("reset bounrdary condition");
#endif

    if (param.mesh.remeshing_option==1 ||
        param.mesh.remeshing_option==2 ||
        param.mesh.remeshing_option==11 ||
        param.mesh.remeshing_option==13) {
        /* Reset coord0 of the bottom nodes */
        int nbot = static_cast<int>(var.bnodes[iboundz0]->size());
#ifndef ACC
        #pragma omp parallel for default(none) shared(param, var, nbot)
#endif
        #pragma acc parallel loop gang vector async
        for (int i=0; i<nbot; ++i) {
            int n = (*var.bnodes[iboundz0])[i];
            (*var.coord0)[n][NDIMS-1] = -param.mesh.zlength;
            // Reest temperature of the bottom nodes to mantle temperature
            (*var.temperature)[n] = var.bottom_temperature;
        }
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
    if (param.sim.has_output_during_remeshing) {
        // the following variables need to be re-computed only when we are
        // outputing right after remeshing
        update_strain_rate(var, *var.strain_rate);
        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
    }

    #pragma acc wait

#ifdef ACC
    {
        size_t free_bytes, total_bytes;
        knn_bvh_mem_info(&free_bytes, &total_bytes);
        std::cout << "      [GPU mem] remesh end:   free="
                  << free_bytes/(1<<20) << " MB / total=" << total_bytes/(1<<20) << " MB\n";
    }
#endif

    std::cout << "  Remeshing finished.\n";

    var.nremesh += 1;
    var.func_time.remesh_time += get_nanoseconds() - time_tmp;
#ifdef NPROF
    nvtxRangePop();
#endif
}
