#include <algorithm>
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
// ADAPT-based optimization and related VTK/Adaptivity usage removed.

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

// Refine-policy size floors, derived once so the Triangle decimation, the MMG size band, the
// tiny-element trigger and the boundary-merge distance share one definition (mesh units).
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


bool is_bottom(uint flag)
{
    return flag & BOUNDZ0;
}


bool is_x0(uint flag)
{
    return flag & BOUNDX0;
}


bool is_x1(uint flag)
{
    return flag & BOUNDX1;
}


bool is_y0(uint flag)
{
    return flag & BOUNDY0;
}


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
            // A TOP node at/below the restored bottom: the domain has necked through. Collapsing it
            // would weld the free surface onto the floor (garbage mesh, NaN); stop instead.
            std::cerr << "Error: top-surface node " << i << " (z = "
                      << qcoord[i*NDIMS + NDIMS-1] << ") has reached the bottom boundary (z = "
                      << bottom << ").\n       The model domain has necked through -- "
                      "this is a terminal model state, not a remeshing problem. Stopping.\n";
            die(EXIT_MESH_QUALITY);   // terminal physical state reached in the mesh-repair path
        }
        else if (qcoord[i*NDIMS + NDIMS-1] < bottom + min_dist) {
            // Mark every NON-bottom node at/below the flattened bottom plane for deletion: interior
            // nodes (they invert the adjoining elements once the boundary snaps back up) and
            // side-wall nodes that sank below the corner (the corner snaps up, its wall neighbour
            // would stay below and degenerate the corner elements). The collapse (MMG) /
            // delete_points (Triangle) then removes them.
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
        // Mark for deletion any NON-x0 node left of the restored plane, including nodes of an
        // ADJACENT boundary (top/bottom) carried past the corner: guarding on `!is_boundary`
        // left an orphan wedge there. Inert on healthy runs.
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


void assemble_bdry_polygons(const Variables &var, const array_t &old_coord,
                             const conn_t &old_connectivity,
                             int_vec (&bdry_polygons)[nbdrytypes])
{
    /* bdry_polygons[i] contains a polygon, ie. a list of vertex, enclosing the i-th boundary */

#ifdef THREED
    const int nodes_per_edge = 2;  // an edge has 2 nodes
    const int edges_per_facet = 3;  // a facet has 3 edges
    const int edgenodes[edges_per_facet][nodes_per_edge] = { {1, 2}, {2, 0}, {0, 1} };

    for (int ibound=0; ibound<nbdrytypes; ibound++) {
        if (var.bnodes[ibound]->size() == 0) continue;  // skip empty boundary

        //
        // Collecting edges on each boundary.
        //

        std::unordered_map<edge_t, int, hash1, equal_to1> edges;
        for (std::size_t i=0; i<var.bfacets[ibound]->size(); i++) {
            const auto &facet = (*(var.bfacets[ibound]))[i];
            int e = facet.first;
            int f = facet.second;
            ConstConnAccessor conn = old_connectivity[e];

            for (int j=0; j<edges_per_facet; j++) {
                int n0 = conn[ NODE_OF_FACET[f][ edgenodes[j][0] ] ];
                int n1 = conn[ NODE_OF_FACET[f][ edgenodes[j][1] ] ];
                if (n0 > n1) {  // ensure n0 < n1
                    int tmp;
                    tmp = n0;
                    n0 = n1;
                    n1 = tmp;
                }
                edge_t g = std::make_pair(n0, n1);
                auto search = edges.find(g);
                if (search == edges.end()) {
                    edges[g] = 1;
                }
                else {
                    search->second ++;
                }
            }
        }
        if (DEBUG > 1) {
            std::cout << ibound << "-th edge:\n";
            for (auto kk=edges.begin(); kk!=edges.end(); ++kk) {
                std::cout << kk->first.first << ",\t" << kk->first.second << "\t: " << kk->second << '\n';
            }
            std::cout << '\n';
        }

        //
        // Collecting edges enclosing the boundary
        //
        std::vector<const edge_t*> enclosing_edges;
        for (auto kk=edges.begin(); kk!=edges.end(); ++kk) {
            int count = kk->second;
            if (count == 2) {
                // this is an "internal" edge
                continue;
            }
            else if (count == 1) {
                // this edge encloses the boundary
                enclosing_edges.push_back(&(kk->first));
            }
            else {
                // not possible
                die(EXIT_MESH_QUALITY, "an edge is belonged to more than 2 facets. The mesh is corrupted.");
            }
        }

        //
        // Connecting edges to form a polygon
        //

        int_vec &polygon = bdry_polygons[ibound];
        const auto g0 = enclosing_edges.begin();
        int head = (*g0)->first;
        int tail = (*g0)->second;
        polygon.push_back(head);
        polygon.push_back(tail);
        auto g1 = g0+1;
        while (head != tail) {
            for (auto g=g1; g!=enclosing_edges.end(); ++g) {
                if ((*g)->first == tail) {
                    tail = (*g)->second;
                    polygon.push_back(tail);
                }
                else if ((*g)->second == tail) {
                    tail = (*g)->first;
                    polygon.push_back(tail);
                }
            }
        }
        // the starting point and end point must be the same
        if (polygon.front() != polygon.back()) {
            die(EXIT_MESH_QUALITY, "boundary polygon is not closed. The mesh is corrupted.");
        }
        polygon.pop_back();  // removed the duplicating end point

        if (DEBUG > 1) {
            std::cout << "nodes for " << ibound << "-th boundary polygon:\n";
            print(std::cout, polygon);
            std::cout << '\n';
        }
    }

#endif
}


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

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        if (is_boundary(flag)) {
            std::replace(segment, endsegment, end, *i);
            // std::cout << *i << " <- " << end << "\n";
            if (DEBUG) {
                std::cout << "segment after replace: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        end --;
        npoints --;
    }

    /* PS: We don't check whether the merged segments belong to the same boundary or not,
     *     e.g., one segment is on the top boundary while the other segment is on the left
     *     boundary. The check is performed in delete_facets() instead.
     */
}


void delete_points_and_merge_facets(const int_vec &points_to_delete,
                                    const int_vec (&bnodes)[nbdrytypes],
                                    const int_vec (&bdry_polygons)[nbdrytypes],
                                    const int_vec (&bdrynode_deleting)[nbdrytypes],
                                    const array_t &bnormals,
                                    int &npoints,
                                    int &nseg, double *points,
                                    int *segment, int *segflag,
                                    uint_vec &bcflag, double min_length)
{
#ifndef THREED
    die(EXIT_INTERNAL_ASSERT, "delete_points_and_merge_facets() doesn't work in 2D!");
#else

    int_vec inverse[nbdrytypes], nfacets;
    std::vector<int*> facet;

    // before deleting boundary points, create a new triangulation
    // TODO: skip boundary bdrynode_deleting.size()==0, but retaining its facets
    for (int i=0; i<nbdrytypes; ++i) {  // looping over all possible boundaries
        const int_vec& bdeleting = bdrynode_deleting[i];  // nodes to be deleted on this boundary, sorted
        const int_vec& bdry_nodes = bnodes[i];  // all nodes on this boundary, sorted
        if (bdry_nodes.size() == 0) continue;

        if (DEBUG) {
            std::cout << i << "-th boundary to be merged: ";
            print(std::cout, bdeleting);
            std::cout << '\n';
        }

        // 3d coordinate -> 2d
        Array2D<double,2> coord2d(bdry_nodes.size() - bdeleting.size());
        int_vec& inv = inverse[i];
        {
            if (DEBUG > 1) {
                std::cout << "bdry nodes:\n";
                print(std::cout, bdry_nodes);
                std::cout << "\n";
                std::cout << bdry_nodes.size() << ' ' << bdeleting.size() << '\n';
            }

            int major_direction = i / 2;  // largest component of normal vector
            if (i >= iboundn0) { // slant boundaries start at iboundn0
                double max = 0;
                for (int d=0; d<NDIMS; d++) {
                    if (std::abs(bnormals[i][d]) > max) {
                        max = std::abs(bnormals[i][d]);
                        major_direction = d;
                    }
                }
            }

            for (std::size_t j=0, k=0, n=0; j<bdry_nodes.size(); j++) {

                // is j belonged to bdeleting[]?
                if (k < bdeleting.size() && bdry_nodes[j] == bdeleting[k]) {
                    k++;
                    continue;
                }

                // record the node # for inverse mapping
                inv.push_back(bdry_nodes[j]);

                int dd = 0;
                for (int d=0; d<NDIMS; d++) {
                    if (d == major_direction) continue;
                    coord2d[n][dd] = points[bdry_nodes[j]*NDIMS + d];
                    dd++;
                }
                n++;
            }

            if (DEBUG > 1) {
                std::cout << i << "-th boundary to be remeshed: ";
                print(std::cout, coord2d);
                std::cout << '\n';
            }
        }

        // re-triangulate the affected boundary (connect only, not adding new points)
        {
            // converting polygon vertex to segment array
            const int_vec& polygon = bdry_polygons[i];
            int_vec surf_segflag(polygon.size()); // all 0, its value does not matter
            int_vec surf_segment(2 * polygon.size());
            std::size_t first = 0;
            int new_polygon_size = 0;
            for (std::size_t j=0; j<polygon.size(); j++) {
                auto search = std::find(inv.begin(), inv.end(), polygon[j]);
                if (search != inv.end()) {  // the vertex is not deleted
                    std::size_t ia = search - inv.begin();
                    if (new_polygon_size == 0) {
                        // start of the polygon
                        surf_segment[0] = ia;
                        first = ia;
                    }
                    else {
                        surf_segment[2*new_polygon_size-1] = surf_segment[2*new_polygon_size] = ia;
                    }
                    ++new_polygon_size;
                }
            }
            // end of the polygon
            surf_segment[2*new_polygon_size-1] = first;

            if (DEBUG) {
                std::cout << "inverse: \n";
                print(std::cout, inv);
                std::cout << '\n';
                std::cout << "polygon segments: ";
                print(std::cout, surf_segment, new_polygon_size*2);
                std::cout << '\n';
            }

            // temporary arrays, some will be allocated inside Triangle library
            int new_nnode, new_nelem, new_nseg;
            double *pcoord, *pregattr;
            int *pconnectivity, *psegment, *psegflag;

            Mesh mesh;
            mesh.min_angle = 0;
            mesh.meshing_verbosity = 0;

            double_vec coord2d_vec;
            coord2d.pack_to(coord2d_vec);
            double* coord2d_ptr = coord2d_vec.data();

            points_to_new_surface(mesh, coord2d.size(), coord2d_ptr,
                                  new_polygon_size, surf_segment.data(), surf_segflag.data(),
                                  0, NULL,
                                  0, 3,
                                  new_nnode, new_nelem, new_nseg,
                                  pcoord, pconnectivity, psegment, psegflag, pregattr);

            if (static_cast<std::size_t>(new_nnode) != coord2d.size()) {
                std::cerr << "Error: ponits_to_new_surface is adding new points!\n";
                std::cout << new_nnode << ' ' << coord2d.size() << '\n';
                std::cout << "old points: ";
                print(std::cout, coord2d);
                std::cout << '\n';
                std::cout << "new points: ";
                print(std::cout, pcoord, new_nnode*2);
                std::cout << '\n';
                std::cout << "new conn: ";
                print(std::cout, pconnectivity, new_nelem*3);
                std::cout << '\n';
                die(EXIT_INTERNAL_ASSERT);
            }
            if (new_nseg != new_polygon_size) {
                die(EXIT_INTERNAL_ASSERT, "points_to_new_surface is adding new segments!");
            }

            delete [] pcoord;
            delete [] psegment;
            delete [] psegflag;
            delete [] pregattr;
            // remember to free pconnectivity later

            // translating index of local (per-boundary) node # to global (whole-mesh) node #
            for (int j=0; j<new_nelem; ++j) {
                for (int k=0; k<NODES_PER_FACET; ++k) {
                    int n = pconnectivity[NODES_PER_FACET*j + k];
                    pconnectivity[NODES_PER_FACET*j + k] = inv[n];
                }
            }

            // storing the new boundary facets for later
            // ownership of "pconnectivity" array is transferred to "facet"
            facet.push_back(pconnectivity);
            nfacets.push_back(new_nelem);
        }
    }

    int nseg2 = std::accumulate(nfacets.begin(), nfacets.end(), 0);
    if (nseg2 > nseg) {
        die(EXIT_INTERNAL_ASSERT, "ponits_to_new_surface too many segments!");
    }

    // appending facets of all boundaries into segment array
    for (int i=0, n=0; i<nbdrytypes; ++i) {
        if (bnodes[i].size() == 0) continue;

        for (int k=0; k<nfacets[i]; ++k, ++n) {
            for (int j=0; j<NODES_PER_FACET; ++j)
                segment[n*NODES_PER_FACET + j] = facet[i][k*NODES_PER_FACET + j];
            segflag[n] = 1 << i;
        }
        delete [] facet[i];
    }

    // mark deleted facets
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


// Delete every point in points_to_delete, merging the neighbouring boundary
// segments (2D) / facets (3D) for the points that lie on a boundary so the
// boundary stays closed. Non-boundary points are simply removed. This is the
// deletion path for remeshing_option >= 10 (boundaries may be modified); it fully
// handles the removal, so the caller must NOT also call delete_points().
void delete_points_and_merge_boundary(int_vec &points_to_delete,
                               const int_vec (&bnodes)[nbdrytypes],
                               const int_vec (&bdry_polygons)[nbdrytypes],
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
    // are there any points in points_to_delete on the boundary? and on which boundary?
    // if the deleted point is a boundary point, store its index
    bool changed = 0;
    int_vec bdrynode_deleting[nbdrytypes];
    for (auto i=points_to_delete.begin(); i<points_to_delete.end(); ++i) {
        uint flag = bcflag[*i];
        for (int j=0; j<nbdrytypes; ++j) {
            uint bc = 1 << j;
            if (flag & bc) {
                bdrynode_deleting[j].push_back(*i);
                changed = 1;
            }
        }
    }

    // non-boundary points changed,
    if (! changed) {
        delete_points(points_to_delete, npoints, nseg,
                      points, segment);
        delete_facets(nseg, segment, segflag);
        return;
    }

    delete_points_and_merge_facets(points_to_delete, bnodes, bdry_polygons,
                                   bdrynode_deleting, bnormals, npoints, nseg,
                                   points, segment, segflag, bcflag, min_size);
    delete_facets(nseg, segment, segflag);
#else
    delete_points_and_merge_segments(points_to_delete, npoints, nseg,
                                     points, segment, bcflag, min_size);
    delete_facets(nseg, segment, segflag);
    // silence 'unused variable' complier warning
    (void) bnodes;
    (void) bdry_polygons;
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

// Per-side INFLOW state of the restored side walls (remeshing_option 13), measured on the OLD
// mesh: a side receives material when its wall moved inward off its restored plane since the
// last remesh. Static or outflow walls get none of the side treatments.
struct SideInflowState {
    double x0_in, x1_in;    // inward extent of the old wall (x0: max x; x1: min x)
    bool in_x0, in_x1;
#ifdef THREED
    double y0_in, y1_in;
    bool in_y0, in_y1;
#endif
    SideInflowState() : x0_in(0.), x1_in(0.), in_x0(false), in_x1(false)
#ifdef THREED
        , y0_in(0.), y1_in(0.), in_y0(false), in_y1(false)
#endif
    {}
};

SideInflowState side_inflow_state(const Param &param, int nnode,
                                  const array_t &coord, const uint_vec &bcflag)
{
    SideInflowState s;
    s.x0_in = 0.0; s.x1_in = param.mesh.xlength;
#ifdef THREED
    s.y0_in = 0.0; s.y1_in = param.mesh.ylength;
#endif
    for (int n = 0; n < nnode; ++n) {
        const uint f = bcflag[n];
        if (f & BOUNDX0) s.x0_in = std::max(s.x0_in, coord[n][0]);
        if (f & BOUNDX1) s.x1_in = std::min(s.x1_in, coord[n][0]);
#ifdef THREED
        if (f & BOUNDY0) s.y0_in = std::max(s.y0_in, coord[n][1]);
        if (f & BOUNDY1) s.y1_in = std::min(s.y1_in, coord[n][1]);
#endif
    }
    // above numeric jitter of a truly static wall, far below one remesh interval of any
    // geological inflow rate
    const double tol = 1e-3 * param.mesh.resolution;
    s.in_x0 = s.x0_in > tol;
    s.in_x1 = s.x1_in < param.mesh.xlength - tol;
#ifdef THREED
    s.in_y0 = s.y0_in > tol;
    s.in_y1 = s.y1_in < param.mesh.ylength - tol;
#endif
    return s;
}

// Give the Triangle path MMG's size floor hmin = resolution * smallest_size^(1/NDIMS): drop
// INTERIOR nodes within hmin of a retained node so re-triangulation of the advected point cloud
// cannot form a sub-hmin element and thrash the tiny-element trigger. Boundary nodes are kept
// (their spacing is the flatten logic's job). Uniform grid of cell hmin -> O(nnode).
// Do not exempt "active" (shear-band) nodes: it lowered the peak plastic strain and added
// elements and tiny-element remeshes.
void decimate_below_hmin(const Param &param, const array_t &coord, int nnode,
                         const uint_vec &bcflag, int_vec &points_to_delete)
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
    // Pass 1: retain all boundary nodes (never decimated) so interior nodes crowding them are removed.
    for (int i=0; i<nnode; ++i) if (is_boundary(bcflag[i])) retain(i);
    // Pass 2: greedily thin interior nodes to >= hmin spacing.
    for (int i=0; i<nnode; ++i) {
        if (is_boundary(bcflag[i])) continue;
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
    int_vec bdry_polygons[nbdrytypes];
    assemble_bdry_polygons(var, original_coord, original_connectivity, bdry_polygons);

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
    int_vec old_bnodes[nbdrytypes];
    for (int i=0; i<nbdrytypes; ++i) {
        old_bnodes[i] = *(var.bnodes[i]);  // copying whole vector
    }

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
    case 13:
        // DO NOT change the corners
        excl_func = &is_corner;
        break;
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
    case 13:
        flatten_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        flatten_x0(old_bcflag, qcoord, points_to_delete, min_dist);
        flatten_x1(old_bcflag, qcoord, param.mesh.xlength, points_to_delete, min_dist);
#ifdef THREED
        flatten_y0(old_bcflag, qcoord, points_to_delete, min_dist);
        flatten_y1(old_bcflag, qcoord, param.mesh.ylength, points_to_delete, min_dist);
#endif
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

    // Triangle path: enforce MMG's hmin floor so no sub-hmin (tiny) element is emitted.
    decimate_below_hmin(param, array_t_qcoord, old_nnode, old_bcflag, points_to_delete);

    // sort points_to_delete and remove duplicates
    {
        std::sort(points_to_delete.begin(), points_to_delete.end());
        auto last = std::unique(points_to_delete.begin(), points_to_delete.end());
        points_to_delete.resize(last - points_to_delete.begin());
    }

    // Record the ORIGINAL nodes boundary remeshing reshaped (deleted here + moved by flatten_*):
    // barycentric_node_interpolation silences its "interior node not found" warning where a new
    // node legitimately maps outside the old outline. Before delete_points() renumbers qcoord.
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
    case 1:
    case 2:
        // deleting non-boundary points
        delete_points(points_to_delete, old_nnode, old_nseg,
                      qcoord, qsegment);
        delete_facets(old_nseg, qsegment, qsegflag);
        break;
    case 10:
    case 11:
    case 12:
    case 13:
        // deleting points, some of them might be on the boundary.
        // delete_points_and_merge_boundary() already removes every point in
        // points_to_delete (merging segments for the boundary ones) and reduces
        // old_nnode accordingly. Calling delete_points() again on the same list
        // would be a *second* deletion pass over an already-compacted array with
        // now-stale indices: its swap-delete + std::replace(segment, end, *i)
        // rewrites the boundary segments of the highest-index nodes (typically the
        // most recently added surface/corner nodes), which detaches the top surface
        // from the side wall and collapses the corner into a spurious cliff. Do NOT
        // call delete_points() here.
        delete_points_and_merge_boundary(points_to_delete, old_bnodes, bdry_polygons, *var.bnormals,
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
void compute_metric_field(const Param &param, const Variables &var, double_vec &metric, double_vec &etmp,
                          const array_t &old_coord, const conn_t &old_connectivity)
{
    /* Desired element size (nodal metric = edge length) for MMG remeshing.
     *
     * Base = init_elem_size_n (frozen initial nodal size), so sizes are MAINTAINED away from
     * any trigger. The target VOLUME is scaled by 1/(1 + coeff*plstrain) and the edge length
     * MMG consumes is that ratio^(1/NDIMS) (aggregation at the end); scaling the edge by the
     * volume ratio was ~NDIMS times too aggressive. Do not re-anchor the base to the CURRENT
     * size: at coeff=0 it leaves the band coarse, at coeff>0 the element count ratchets.
     *
     * conn must be the PRE-remesh connectivity: var.connectivity is steal_ref'd null here.
     */
    // Not const-qualified on purpose: a const scalar is predetermined-shared in OpenMP and
    // may not appear in a shared() clause on some toolchains -- a plain local sidesteps that.
    double coeff = param.mesh.mmg_metric_refine_coeff;

    // Do not inflate the target by the convergence rate exp(-div v * tau) (no measurable effect)
    // nor treat sub-floor REFINE elements as REPAIR (it frees the nodes that preserve the band).

    // Reuse metric (= *var.ntmp) as scratch for the nodal aggregation below.
    std::fill_n(metric.begin(), var.nnode, 0.0);

    // etmp[e] = vol[e] / (1 + coeff*plstrain[e])  (plastic-strain-reduced element volume, aggregated to nodes)
    #pragma omp parallel for default(none) shared(var, etmp, coeff)
    for (int e = 0; e < var.nelem; e++)
        etmp[e] = (*var.volume)[e] / (1.0 + coeff * (*var.plstrain)[e]);

#ifndef THREED
    // SURFACE-MOTION ANTICIPATION for the top elements: inflate a top element's target volume by
    // the deformation expected over one remesh interval, so it is rebuilt with room to be crushed.
    // Movers: the surface-process dh (bc.cxx) and the advective nodal velocity (a fast slide
    // moves the surface by km per interval; dh alone is blind to it).
    //   dA2/dt  = shoelace area rate from (vx, vz + dh/dt), exact for the linear triangle;
    //             uniform advection cancels, only DEFORMING motion inflates.
    //   inflate = exp(-rate * tau), tau = time since the last remesh, capped at 8x volume.
    // 2D only (dh is the 2D top-chain bookkeeping). *var.coord / *var.connectivity are EMPTY here.
    if (var.top_elems && var.surfinfo.dh && var.dt > 0.0) {
        // global node id -> last step's surface-process rate dh/dt (m/s)
        std::vector<double> dhdt(var.nnode, 0.0);
        for (int i = 0; i < var.surfinfo.ntop; ++i)
            dhdt[(*var.surfinfo.top_nodes)[i]] = (*var.surfinfo.dh)[i] / var.dt;
        const double tau = var.time - var.last_remesh_time;
        int n_inflated = 0;
        for (int k = 0; k < var.ntop_elems; ++k) {
            const int e = (*var.top_elems)[k];
            ConstConnAccessor conn = old_connectivity[e];
            const double xa = old_coord[conn[0]][0], za = old_coord[conn[0]][1];
            const double xb = old_coord[conn[1]][0], zb = old_coord[conn[1]][1];
            const double xc = old_coord[conn[2]][0], zc = old_coord[conn[2]][1];
            const double area2 = (xb - xa) * (zc - za) - (xc - xa) * (zb - za);
            // nodal motion rates: advection + surface process (vertical)
            double ux[3], uz[3];
            for (int i = 0; i < 3; ++i) {
                const int n = conn[i];
                ux[i] = (*var.vel)[n][0];
                uz[i] = (*var.vel)[n][NDIMS-1] + dhdt[n];
            }
            // d(2A)/dt of the shoelace with moving vertices (exact for linear motion)
            const double dA2dt = ux[0]*(zb - zc) + ux[1]*(zc - za) + ux[2]*(za - zb)
                               + uz[0]*(xc - xb) + uz[1]*(xa - xc) + uz[2]*(xb - xa);
            const double rate = dA2dt / area2;          // signed d(lnV)/dt
            if (rate < 0.0) {
                // cap 8x volume (~2.8x edge in 2D): heuristic headroom for a fast slide, untuned
                etmp[e] *= std::min(8.0, std::exp(-rate * tau));
                ++n_inflated;
            }
        }
        if (n_inflated)
            std::cout << "    Surface-motion anticipation: inflated the metric of "
                      << n_inflated << " top element(s).\n";
    }
#endif

    #pragma omp parallel for default(none) shared(var, metric, etmp)
    for (int n = 0; n < var.nnode; n++) {
        // vol_ratio = (sum of reduced element volumes) / (nodal volume); == 1 where pls == 0
        double vol_ratio = 0.0;
        const int npatch = var.support.size(n);
        const int* patch = var.support.patch(n);
        for (int i=0; i<npatch; ++i)
            vol_ratio += etmp[patch[i]];
        vol_ratio /= (*var.volume_n)[n];
        // edge-length metric = frozen edge length * (volume ratio)^(1/NDIMS)
#ifdef THREED
        metric[n] = (*var.init_elem_size_n)[n] * std::cbrt(vol_ratio);
#else
        metric[n] = (*var.init_elem_size_n)[n] * std::sqrt(vol_ratio);
#endif
    }
    // FLOOR the metric one hysteresis step above MMG's hmin: a target at hmin puts the band on the
    // size floor, where MMG's output scatter lands inside the tiny-element re-trigger buffer and
    // every remesh re-fires. hmin * margin^(1/NDIMS) keeps the scatter's low tail a full margin
    // (in volume) above the trigger.
    {
        double metric_floor = refine_floors(param.mesh).hmin
                              * std::pow(param.mesh.remesh_tiny_margin, 1.0 / NDIMS);
        #pragma omp parallel for default(none) shared(var, metric, metric_floor)
        for (int n = 0; n < var.nnode; n++)
            metric[n] = std::max(metric[n], metric_floor);
    }

    // RESTORED side walls (remeshing_option 13): incoming material is restored from per-side
    // profiles sampled at wall-node spacing, so wall-faceted elements stay at ~resolution. Clamp
    // only the METRIC (init_elem_size_n keeps the cfg size, so off-wall elements coarsen back).
    // Applies to ALL restored sides; only the unfreeze and the field restore are inflow-gated.
    // *var.coord is EMPTY here: never feed side_inflow_state from var in this function.
    if (param.mesh.remeshing_option == 13) {
        uint side_bits = BOUNDX0 | BOUNDX1;
#ifdef THREED
        side_bits |= BOUNDY0 | BOUNDY1;
#endif
        double res = param.mesh.resolution;
        #pragma omp parallel for default(none) shared(var, metric, side_bits, res)
        for (int n = 0; n < var.nnode; n++)
            if ((*var.bcflag)[n] & side_bits)
                metric[n] = std::min(metric[n], res);
    }
}

// ============================================================================
//  FREEZE / REQUIRED POLICY
// ============================================================================
// Remeshing must repair bad elements while PRESERVING the fields, above all the localized
// plastic strain: every rebuilt entity is re-interpolated, which diffuses it. Triangle
// re-triangulates the existing points, so preservation is implicit (its only policy is the
// shared size floors, RefineFloors). MMG rebuilds everything not forbidden, so preservation is
// explicit through MMG's required VERTEX (never moved/deleted) and required ELEMENT (untouched).
//
// Per-entity freedom; the complement is handed to MMG as constraints:
//   node : MOVABLE, or PINNED (required vertex, nodal fields carried verbatim).
//   elem : FREE   -- split/collapse/swap AND its nodes MOVABLE (full repair);
//          SPLIT  -- split/swap only, nodes NOT freed (refinement inserts nodes without moving
//                    the old ones, so the band's plastic strain is carried verbatim);
//          FROZEN -- MMG-required IF all its nodes ended up PINNED (DERIVED REQUIRED).
//
// Rules, applied in order by mark_quiet_required; each may only RAISE freedom:
//   R1 REPAIR    quality < min_quality * mmg_remesh_defensive_quality_ratio, or volume <
//                smallest_vol -> FREE (a node where disjoint material bodies meet stays pinned).
//   R2 REFINE    plstrain - plstrain_remesh > mmg_remesh_active_plstrain -> SPLIT.
//   R3 BOUNDARY  only an ACTIVE free-surface (TOP) node is MOVABLE (surface-process displacement
//                > 1% of resolution per interval). Quiet top nodes, every pure wall node of a
//                restored side (option 13) and every BOTTOM node stay pinned, carried verbatim.
//                A fully pinned top lid starves MMG of slack and death-spirals.
//   R4 SIDE STRIP (option 13, INFLOW sides) wall-faceted elements + 1 ring -> FREE with movable
//                interior nodes, so the ~resolution wall band coarsens back once off the wall.
//   R5 SIZE RECOVERY a quiet element FINER than the gradation envelope (the smallest size MMG's
//                hgrad could legally assign, grown from every intentionally-fine target) by more
//                than mmg_remesh_size_recovery_ratio -> FREE to coarsen; >= 2x its own target -> FREE
//                to split. Skips top-connected elements, material-interface elements (for
//                coarsening) and the wall band + 4 rings of a FIXED side; never unpins a wall or
//                bottom node.
//   R6 FLATTEN-BROKEN elements flatten_* broke (post-flatten quality < min_quality) -> FREE.
//   R7 COLLAPSE-REGION (collapse path) elements reshaped by collapse_outside_nodes + 1 ring -> FREE.
//   DERIVED REQUIRED: required_node = NOT movable; required_elem = FROZEN with all nodes pinned.
//   An element with any movable node is never required: the built-in one-element hinge layer.
//
// POST-OUTPUT QUALITY GATE (mmg_adapt_quality_gated): if MMG's output contains an element that
// would re-trigger remeshing, unfreeze a graded region around it (CORE: full freedom, grown to
// frozen contact + `attempt` rings; HINGE: bordering elements lose required-element status but
// keep pinned nodes) and re-run MMG. Bounded attempts.
//
// INIT-TIME (mesh.cxx refine_initial_side_walls_2d, option 13): the same idea at t=0 --
// everything required except the wall band + 2 rings.
//
// KNOBS ([mesh]): min_quality, mmg_remesh_defensive_quality_ratio (R1), mmg_remesh_active_plstrain (R2),
// smallest_size/resolution (floors), remesh_tiny_margin (trigger buffer, metric floor, gate),
// mmg_remesh_size_recovery_ratio (R5, 0 = off).
//
// Measured dead ends, do not retry: an unbounded undersized flood-fill unfreeze (re-frees the
// hgrad transition forever); freeing wider than REPAIR around flatten damage (NaN); dissolving
// small required islands (more retries); raising R1's tiny floor to the metric floor (diffuses
// the band); an explicit extra hinge ring (no-op); SPLIT instead of FREE for R1's quality term
// (more remeshes, no preservation gain).
// ============================================================================


// ----------------------------------------------------------------------------
//  Packed-coords element measure / quality (policy tests + collapse machinery)
// ----------------------------------------------------------------------------

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
// elem_quality on packed coords (0..1, equilateral = 1), measured from THESE coords instead of
// var.volume so it can judge an element after flatten_* moved its nodes; <= 0 = degenerate/inverted.
double packed_elem_quality(const double *coord, const int *elem)
{
    double m = signed_elem_measure(coord, elem);
#ifdef THREED
    if (m <= 0.0) return m;
    const double vol = m / 6.0;
    double area_sum = 0.0;
    static const int faces[4][3] = {{0,1,2},{0,1,3},{2,3,0},{2,3,1}};
    for (int f = 0; f < 4; ++f) {
        const double *a = coord + elem[faces[f][0]]*3;
        const double *b = coord + elem[faces[f][1]]*3;
        const double *c = coord + elem[faces[f][2]]*3;
        double ux=b[0]-a[0], uy=b[1]-a[1], uz=b[2]-a[2];
        double vx=c[0]-a[0], vy=c[1]-a[1], vz=c[2]-a[2];
        double nx=uy*vz-uz*vy, ny=uz*vx-ux*vz, nz=ux*vy-uy*vx;
        area_sum += 0.5 * std::sqrt(nx*nx + ny*ny + nz*nz);
    }
    return 216.0 * std::sqrt(3.0) * vol * vol / (area_sum * area_sum * area_sum);
#else
    if (m <= 0.0) return m;
    const double area = 0.5 * m;
    double d2_sum = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double *a = coord + elem[k]*2, *b = coord + elem[(k+1)%3]*2;
        double dx=b[0]-a[0], dy=b[1]-a[1];
        d2_sum += dx*dx + dy*dy;
    }
    return 4.0 * std::sqrt(3.0) * area / d2_sum;
#endif
}


// ----------------------------------------------------------------------------
//  Freedom masks + rules R1-R5 (the OLD-mesh part of the policy)
// ----------------------------------------------------------------------------

// Element freedom grades (section VOCABULARY). Only the FROZEN/not-FROZEN distinction
// reaches MMG directly; FREE vs SPLIT differ through the nodes each rule frees.
enum ElemGrade : char { ELEM_FROZEN = 0, ELEM_SPLIT = 1, ELEM_FREE = 2 };

struct FreezeMask {
    std::vector<char> node_movable;   // 1 = MOVABLE, 0 = PINNED
    std::vector<char> elem_grade;     // ElemGrade
};

// Wall-band ring growth shared by R4, R5's fixed-side band and mesh.cxx's init-time refinement:
// pass 0 marks the elements faceting a seed node, passes 1..nrings free one connectivity ring
// each. band_node is in/out (seed with the wall nodes); band_elem receives the freed elements;
// elem_gate skips ineligible elements. Returns the number freed.
// Array2D is SoA: never row-index .data() (conn.data() + e*NODES_PER_ELEM reads a transposed
// garbage triangle of valid node ids); always go through the accessors.
int free_wall_band(int nnode, int nelem, const conn_t &conn, int nrings,
                   const std::vector<char> *elem_gate,
                   std::vector<char> &band_node, std::vector<char> &band_elem)
{
    int n_freed = 0;
    for (int r = 0; r <= nrings; ++r) {
        std::vector<char> next(nnode, 0);
        for (int e = 0; e < nelem; ++e) {
            if (band_elem[e]) continue;
            if (elem_gate && !(*elem_gate)[e]) continue;
            ConstConnAccessor el = conn[e];
            bool touch = false;
            for (int k = 0; k < NODES_PER_ELEM && !touch; ++k) touch = band_node[el[k]];
            if (!touch) continue;
            band_elem[e] = 1; ++n_freed;
            for (int k = 0; k < NODES_PER_ELEM; ++k) next[el[k]] = 1;
        }
        for (int n = 0; n < nnode; ++n)
            if (next[n]) band_node[n] = 1;
    }
    return n_freed;
}

// R1 REPAIR + R2 REFINE (one pass):
//   REPAIR = distorted (quality < min_quality) or tiny (measure < smallest_vol) -> FREE. The tiny
//            term is mandatory (the return-3 trigger is a hard failure). The quality term keeps a
//            node pinned where disjoint material bodies meet (moving it drags the boundary).
//   REFINE = yielded since the last remesh -> SPLIT: MMG adds resolution without moving the
//            existing nodes, so their plastic strain is carried verbatim.
// Do not give the quality term SPLIT instead of FREE: MMG cannot fix distortion without node
// motion, the badness survives to the gate's blunter unfreeze, and preservation does not improve.
void rule_repair_and_refine(const Param &param, const Variables &var,
                            const array_t &coord, const conn_t &connectivity,
                            int nelem, FreezeMask &m)
{
    // DEFENSIVE band: repair everything below min_quality * ratio, so a marginal element is
    // rebuilt with headroom now instead of re-triggering a remesh a few steps later.
    const double q_thr = param.mesh.min_quality * param.mesh.mmg_remesh_defensive_quality_ratio;
    // Keep the raw floor: raising it to the metric floor frees REFINE-pinned band nodes.
    const double repair_vol = refine_floors(param.mesh).smallest_vol;
    const double pls_thr = param.mesh.mmg_remesh_active_plstrain;

    // Marker-material bitmask per element (mattypes >= 64 lumped into the top bit): two elements
    // are disjoint bodies iff the masks do not intersect; a markerless element counts as disjoint.
    std::vector<uint64_t> mat_mask(nelem, 0);
    for (int e = 0; e < nelem; ++e) {
        const int_vec &em = (*var.elemmarkers)[e];
        for (int t = 0; t < (int)em.size(); ++t)
            if (em[t] > 0) mat_mask[e] |= 1ULL << std::min(t, 63);
    }

    int n_iface_pinned = 0;
    for (int e = 0; e < nelem; ++e) {
        // Match bad_mesh_quality's quality convention: normalize by 1/NDIMS in 3D before the
        // min_quality comparison, so the freeze test and the remesh trigger use the same scale.
        double q = elem_quality(coord, connectivity, *var.volume, e);
#ifdef THREED
        q = std::pow(q, 1.0 / 3);
#endif
        const bool low_q = q < q_thr;
        const bool tiny  = (*var.volume)[e] < repair_vol;
        const bool refine = ((*var.plstrain)[e] - (*var.plstrain_remesh)[e]) > pls_thr; // TODO: change to pls rate
        m.elem_grade[e] = (low_q || tiny) ? ELEM_FREE
                                          : (refine ? ELEM_SPLIT : ELEM_FROZEN);
        if (tiny) {                        // tiny-only: unguarded (hard failure)
            ConstConnAccessor conn = connectivity[e];
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                m.node_movable[conn[i]] = 1;
        } else if (low_q) {
            // LOW QUALITY: free the nodes, except one that also touches an element sharing NO
            // mattype with this one (disjoint bodies meet there; moving it drags the boundary).
            ConstConnAccessor conn = connectivity[e];
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                const int n = conn[i];
                bool disjoint_nbr = false;
                const int npatch = var.support.size(n);
                const int* patch = var.support.patch(n);
                for (int k = 0; k < npatch; ++k) {
                    const int ne = patch[k];
                    if (ne != e && (mat_mask[ne] & mat_mask[e]) == 0) {
                        disjoint_nbr = true;
                        break;
                    }
                }
                if (disjoint_nbr) { ++n_iface_pinned; continue; }
                m.node_movable[n] = 1;
            }
        }
    }
    if (n_iface_pinned)
        std::cout << "    Repair: kept " << n_iface_pinned
                  << " material-interface node ref(s) of low-quality elements required.\n";
}

// R3 BOUNDARY: only an ACTIVE free-surface (TOP) node becomes movable (flatten_* moves it).
// Quiet top nodes, pure wall nodes of restored sides (clamp_side_bits) and bottom nodes stay
// pinned so their row/column is carried verbatim.
void rule_boundary_nodes(const Param &param, const Variables &var, int nnode,
                         uint clamp_side_bits, FreezeMask &m)
{
    // Surface-process ACTIVITY per top node: expected displacement over one remesh interval,
    // |dh|/dt * tau. Below a small fraction of the resolution the node is static and carried
    // verbatim. The gate is mandatory: pinning the WHOLE top row starves MMG of the slack to
    // rebuild a converging wedge underneath and escalates bad outputs until an element folds.
    std::vector<char> top_active(nnode, 0);
#ifndef THREED
    if (var.surfinfo.dh && var.dt > 0.0) {
        const double tau = var.time - var.last_remesh_time;
        // 1% of resolution per remesh interval: heuristic static/active split, untuned
        const double tol = 0.01 * param.mesh.resolution;
        for (int i = 0; i < var.surfinfo.ntop; ++i) {
            const int n = (*var.surfinfo.top_nodes)[i];
            if (std::abs((*var.surfinfo.dh)[i]) / var.dt * tau > tol)
                top_active[n] = 1;
        }
    }
#endif

    for (int n = 0; n < nnode; ++n) {
        const uint f = (*var.bcflag)[n];
        if (!f) continue;
        // QUIET TOP nodes stay required: nothing in the pipeline moves the free surface, and R1
        // already frees the nodes of a broken surface element. Freeing them wholesale
        // re-tessellated the whole top row every remesh (field diffusion, topography noise within
        // hausd). Includes a quiet top corner shared with a restored wall: left movable, MMG
        // diagonal-flipped the corner elements every remesh.
        if ((f & BOUNDZ1) && !(f & (BOUND_ANY & ~(BOUNDZ1 | clamp_side_bits)))
            && !top_active[n]) continue;
        // Pure wall nodes of ANY restored side stay required (column carried verbatim); an R1
        // REPAIR element still frees its wall nodes, so a broken wall element is always fixable.
        if ((f & clamp_side_bits) && !(f & ~clamp_side_bits)) continue;

        // BOTTOM nodes (incl. bottom corners) stay required: flatten_bottom snaps the row in the
        // boundary outline and remesh() resets its coord0, so nothing needs them movable.
        if (f & BOUNDZ0) continue;
        m.node_movable[n] = 1;              // active free-surface nodes: movable for flattening
    }
}

// R4 SIDE STRIP (option 13, INFLOW sides): wall-faceted elements are held at ~resolution while
// their material advects inland; off the wall they would stay frozen fine forever and the band
// would grow by the advection distance per remesh. Free exactly the wall-faceted elements + ONE
// connectivity rings and make their INTERIOR nodes movable (boundary nodes stay R3-pinned), so
// MMG merges the just-off-wall elements into their coarse neighbours. Nothing wider: an
// unbounded undersized flood-fill re-freed the hgrad transition zone at every remesh.
// Only sides actually RECEIVING material (side_inflow_state) are treated.
void rule_inflow_side_strip(const Variables &var, const conn_t &connectivity,
                            int nnode, int nelem, uint inflow_side_bits, FreezeMask &m)
{
    std::vector<char> band_node(nnode, 0), band_elem(nelem, 0);
    for (int n = 0; n < nnode; ++n)
        if ((*var.bcflag)[n] & inflow_side_bits)
            band_node[n] = 1;                       // ring 0 grows from the wall nodes
    const int nrings = 1;
    const int n_freed = free_wall_band(nnode, nelem, connectivity, nrings,
                                       NULL, band_node, band_elem);
    for (int e = 0; e < nelem; ++e) {
        if (!band_elem[e]) continue;
        m.elem_grade[e] = ELEM_FREE;
        ConstConnAccessor conn = connectivity[e];

        // Find the most abundant marker mattype in this element
        int_vec &a = (*var.elemmarkers)[e];
        int mattype = std::distance(a.begin(), std::max_element(a.begin(), a.end()));

        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            const int npatch = var.support.size(conn[i]);
            const int* patch = var.support.patch(conn[i]);
            bool is_iface = false;
            for (int j = 0; j < npatch; ++j) {
                int e_sup = patch[j];
                int_vec &b = (*var.elemmarkers)[e_sup];
                int mat_sup = std::distance(b.begin(), std::max_element(b.begin(), b.end()));
                if (mat_sup != mattype) {
                    is_iface = true;  // mark as interface node
                    break;
                }
            }
            if (is_iface) continue;

            // Free the band ELEMENTS; only interior band nodes become movable (R3 keeps the
            // boundary pinned), so MMG re-meshes the band without moving the boundary.
            if (!is_boundary((*var.bcflag)[conn[i]]))
                m.node_movable[conn[i]] = 1;
        }
    }
    if (n_freed)
        std::cout << "    Side strip: freed " << n_freed
                  << " wall-faceted element(s) + " << nrings
                  << " connectivity rings to resize back.\n";
}

// Required MATERIAL-INTERFACE edges of the inflow band (option 13, 2D): an edge shared by two
// band elements of different dominant mattype lies on the interface. Marked MMG-required, the
// remesh keeps the boundary conforming and cannot coarsen a straddling element across a thin
// layer. Re-derived every remesh from the markers. Returns 0-based OLD-mesh node pairs; the
// collapse path translates them through its remap. 3D returns empty (an interface there is a
// face -> required triangles, not edges).
std::vector<int> build_interface_req_edges(const Variables &var, const conn_t &connectivity,
                                           int nnode, int nelem, uint inflow_side_bits)
{
#ifdef THREED
    (void)var; (void)connectivity; (void)nnode; (void)nelem; (void)inflow_side_bits;
    return {};
#else
    if (inflow_side_bits == 0) return {};
    std::vector<char> band_node(nnode, 0), band_elem(nelem, 0);
    for (int n = 0; n < nnode; ++n)
        if ((*var.bcflag)[n] & inflow_side_bits) band_node[n] = 1;
    free_wall_band(nnode, nelem, connectivity, 1, NULL, band_node, band_elem);

    auto dominant = [&](int e) -> int {
        const int_vec &a = (*var.elemmarkers)[e];
        return (int)std::distance(a.begin(), std::max_element(a.begin(), a.end()));
    };
    // An interior edge is shared by exactly two triangles; emit it once, on the second share,
    // when the two dominant mattypes differ.
    std::map<std::pair<int,int>, int> first_mat;      // edge -> dominant mattype of its first elem
    std::vector<int> out;
    for (int e = 0; e < nelem; ++e) {
        if (!band_elem[e]) continue;
        const int mt = dominant(e);
        ConstConnAccessor el = connectivity[e];
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            for (int j = i + 1; j < NODES_PER_ELEM; ++j) {
                std::pair<int,int> key(std::min(el[i], el[j]), std::max(el[i], el[j]));
                auto res = first_mat.emplace(key, mt);
                if (!res.second && res.first->second != mt) {
                    out.push_back(key.first);
                    out.push_back(key.second);
                    res.first->second = mt;           // guard against a 3rd share re-emitting
                }
            }
    }
    return out;
#endif
}

// Pre-subdivide the required interface edges to ~metric BEFORE MMG: required edges (the only tag
// that reliably KEEPS an internal line; opnbdy/ref edges let MMG curve or cross it) cannot be
// split by MMG, so an over-long edge is bisected here at its midpoint (collinear -> exact line),
// the two adjacent triangles split, and the midpoint + half-edges marked required. Converges
// over remeshes (each pass halves an over-long edge). 2D only. Grows coord/conn/metric/reqn/reqe
// in place and rewrites `iface`; a triangle is split at most once per pass (tri_consumed).
void presubdivide_interface_edges(std::vector<double> &coord, std::vector<int> &conn,
                                  std::vector<double> &metric, std::vector<char> &reqn,
                                  std::vector<char> &reqe, std::vector<int> &iface,
                                  int &nn, int &ne, double resolution)
{
#ifndef THREED
    if (iface.empty()) return;
    const double FACTOR = 1.41;   // MMG splits an edge longer than ~sqrt(2)*target; pre-do it here
    auto keyof = [](int a, int b){ return std::make_pair(std::min(a,b), std::max(a,b)); };
    std::map<std::pair<int,int>, int> ifmap;                             // iface-edge membership
    for (std::size_t k = 0; k + 1 < iface.size(); k += 2) ifmap[keyof(iface[k], iface[k+1])] = 1;
    std::map<std::pair<int,int>, std::vector<std::pair<int,int>>> adj;   // iface edge -> [(tri, opp node)]
    for (int e = 0; e < ne; ++e) {
        const int n0 = conn[e*3], n1 = conn[e*3+1], n2 = conn[e*3+2];
        const int ed[3][2] = {{n0,n1},{n1,n2},{n2,n0}}; const int opp[3] = {n2,n0,n1};
        for (int s = 0; s < 3; ++s) {
            auto k = keyof(ed[s][0], ed[s][1]);
            if (ifmap.count(k)) adj[k].push_back({e, opp[s]});
        }
    }
    std::vector<char> tri_consumed(ne, 0);
    std::vector<int> out; out.reserve(iface.size() * 2);
    for (std::size_t k = 0; k + 1 < iface.size(); k += 2) {
        const int a = iface[k], b = iface[k+1];
        const double dx = coord[a*2] - coord[b*2], dz = coord[a*2+1] - coord[b*2+1];
        const double L = std::sqrt(dx*dx + dz*dz);
        double tgt = std::min(metric[a], metric[b]); if (tgt <= 0.0) tgt = resolution;
        auto &tris = adj[keyof(a,b)];
        const bool ok = (L > FACTOR * tgt) && tris.size() == 2
                        && !tri_consumed[tris[0].first] && !tri_consumed[tris[1].first];
        if (!ok) { out.push_back(a); out.push_back(b); continue; }
        const int m = nn++;
        coord.push_back(0.5*(coord[a*2]   + coord[b*2]));
        coord.push_back(0.5*(coord[a*2+1] + coord[b*2+1]));
        // Size the divided line to the SMALLEST target of its surrounding nodes (endpoints + the
        // two apexes) and pull those nodes down to it: the fine wall size then propagates along
        // the in-band interface as it advects inland instead of creeping upward. Monotone (min
        // only, floored where compute_metric_field floored it).
        const int c0 = tris[0].second, c1 = tris[1].second;
        const double h = std::min({metric[a], metric[b], metric[c0], metric[c1]});
        metric.push_back(h);
        metric[a] = h; metric[b] = h; metric[c0] = h; metric[c1] = h;
        reqn.push_back(1);
        for (auto &tc : tris) {
            const int e = tc.first, c = tc.second;
            tri_consumed[e] = 1;
            conn[e*3] = a; conn[e*3+1] = m; conn[e*3+2] = c;              // reuse slot e -> (a,m,c)
            conn.push_back(m); conn.push_back(b); conn.push_back(c);      // append (m,b,c)
            reqe.push_back(reqe[e]); ++ne;                               // child inherits requiredness
        }
        out.push_back(a); out.push_back(m); out.push_back(m); out.push_back(b);
    }
    iface.swap(out);
#else
    (void)coord;(void)conn;(void)metric;(void)reqn;(void)reqe;(void)iface;(void)nn;(void)ne;(void)resolution;
#endif
}

// R5 SIZE RECOVERY: free fossil-fine quiet elements so MMG can coarsen them back. The wall clamp
// (option 13) and the plstrain refinement create elements finer than init_elem_size_n; once quiet
// and inland they would be carried verbatim forever (the element count roughly doubled over a
// long inflow run). The criterion is the GRADATION ENVELOPE, not the hint: freeing everything
// undersized vs the hint re-frees the hgrad transition that MMG itself forces.
//   1. per-node fine target t(n), mirroring compute_metric_field (hint reduced by the plstrain
//      refinement, floored, clamped at restored walls);
//   2. achievable(n) = min over nodes s of t(s) * hgrad^rings(n,s), by fixpoint relaxation;
//   3. FOSSIL EXCESS iff esize * ratio < achievable at all its nodes -> FREE (nodes deletable).
//      ratio (mmg_remesh_size_recovery_ratio > 1) is the hysteresis; 0 disables the rule.
// Band and transition elements sit on the envelope and are never freed. Also skipped: R2
// elements, top-surface-connected elements (surface guard), material-interface elements for
// coarsening (interface guard), and the wall band + 4 rings of a FIXED side. A freed element's
// wall (clamp_side_bits) and bottom nodes are never unpinned.
void rule_size_recovery(const Param &param, const Variables &var,
                        const conn_t &connectivity,
                        int nnode, int nelem, uint clamp_side_bits,
                        uint fixed_side_bits, FreezeMask &m)
{
    const double ratio = param.mesh.mmg_remesh_size_recovery_ratio;
    if (ratio <= 0.0) return;                                    // disabled
    if ((int)var.init_elem_size_n->size() != nnode) return;      // no hint yet (pre-first-remesh)

    // Fixed-side band: the FIXED restored wall column is carried verbatim, so mark its
    // wall-faceted elements + 4 rings and skip them below. Empty without a fixed side.
    std::vector<char> fixed_band_elem(nelem, 0);
    if (fixed_side_bits) {
        std::vector<char> fixed_band_node(nnode, 0);
        for (int n = 0; n < nnode; ++n)
            if ((*var.bcflag)[n] & fixed_side_bits)
                fixed_band_node[n] = 1;
        free_wall_band(nnode, nelem, connectivity, 4, NULL, fixed_band_node, fixed_band_elem);
    }

    const double hgrad = 1.3;   // MMG's default gradation; mmg_adapt does not override it
    const double coeff = param.mesh.mmg_metric_refine_coeff;
    const double metric_floor = refine_floors(param.mesh).hmin
                                * std::pow(param.mesh.remesh_tiny_margin, 1.0 / NDIMS);

    // (1) per-node fine target
    std::vector<double> achievable(nnode, std::numeric_limits<double>::max());
    std::vector<double> esize(nelem);
    for (int e = 0; e < nelem; ++e) {
        ConstConnAccessor conn = connectivity[e];
#ifdef THREED
        esize[e] = std::cbrt((*var.volume)[e] / sizefactor);
        const double red = std::cbrt(1.0 / (1.0 + coeff * (*var.plstrain)[e]));
#else
        esize[e] = std::sqrt((*var.volume)[e] / sizefactor);
        const double red = std::sqrt(1.0 / (1.0 + coeff * (*var.plstrain)[e]));
#endif
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            const int n = conn[i];
            const double t = (*var.init_elem_size_n)[n] * red;
            if (t < achievable[n]) achievable[n] = t;
        }
    }
    for (int n = 0; n < nnode; ++n) {
        if ((*var.bcflag)[n] & clamp_side_bits)
            achievable[n] = std::min(achievable[n], param.mesh.resolution);  // wall clamp
        achievable[n] = std::max(achievable[n], metric_floor);
    }

    // Per-node target BEFORE the gradation sweep: the oversized rule compares an element against
    // its OWN target, never against an envelope a fine neighbour graded down.
    std::vector<double> target_raw = achievable;

    // (2) gradation envelope: fixpoint of achievable(n) <= min(elem neighbours) * hgrad.
    //     Sweeps needed = mesh diameter in rings; 512 is a heuristic cap on that diameter.
    bool converged = false;
    for (int sweep = 0; sweep < 512 && !converged; ++sweep) {
        bool changed = false;
        for (int e = 0; e < nelem; ++e) {
            ConstConnAccessor conn = connectivity[e];
            double amin = achievable[conn[0]];
            for (int i = 1; i < NODES_PER_ELEM; ++i)
                amin = std::min(amin, achievable[conn[i]]);
            const double lim = amin * hgrad;
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                if (achievable[conn[i]] > lim) { achievable[conn[i]] = lim; changed = true; }
        }
        converged = !changed;
    }
    if (!converged)
        std::cout << "      [freeze] WARNING: R5 gradation envelope not converged in 512 sweeps; "
                     "size recovery may free transition elements.\n";

    // Surface guard: never free an element connected to a top-surface node (var.top_elems, old
    // mesh). R5's freeing is additive and would unpin QUIET surface nodes that R3 pinned; such an
    // element stays a required hinge. R3's activity gate still frees active surface spots.
    std::vector<char> is_top_elem(nelem, 0);
    for (int i = 0; i < var.ntop_elems; ++i)
        is_top_elem[(*var.top_elems)[i]] = 1;

    // Material-interface guard: never free an element with mixed markers, or sharing a node with
    // an element of a different mattype: coarsening remaps the patch's markers and diffuses the
    // boundary. elem_mat[e] = the sole mattype, or -1 if mixed/markerless.
    std::vector<int> elem_mat(nelem, -1);
    for (int e = 0; e < nelem; ++e) {
        const int_vec &em = (*var.elemmarkers)[e];
        int mt = -1;
        for (int t = 0; t < (int)em.size(); ++t) {
            if (em[t] <= 0) continue;
            if (mt >= 0) { mt = -1; break; }   // second mattype -> mixed
            mt = t;
        }
        elem_mat[e] = mt;
    }
    // (3) free the size excess in BOTH directions: fossil-fine (finer than the envelope by more
    //     than `ratio`) -> free to COARSEN; oversized (>= 2x its own target, e.g. frozen coarse
    //     while the metric refined under it) -> free to SPLIT (MMG cannot split a required element).
    int n_freed = 0, n_oversized = 0, n_top_skipped = 0, n_mat_skipped = 0, n_fixed_skipped = 0;
    for (int e = 0; e < nelem; ++e) {
        if (m.elem_grade[e] == ELEM_SPLIT) continue;   // R2 growing band: nodes stay pinned
        ConstConnAccessor conn = connectivity[e];
        double ach = achievable[conn[0]], tgt = target_raw[conn[0]];
        for (int i = 1; i < NODES_PER_ELEM; ++i) {
            ach = std::min(ach, achievable[conn[i]]);
            tgt = std::max(tgt, target_raw[conn[i]]);
        }
        const bool fossil_fine = (esize[e] * ratio < ach);   // finer than envelope -> coarsen
        const bool oversized   = (esize[e] >= 2.0 * tgt);     // >= 2x its target -> split
        if (!fossil_fine && !oversized) continue;             // within the legal size band
        if (is_top_elem[e]) { ++n_top_skipped; continue; }
        if (fixed_band_elem[e]) { ++n_fixed_skipped; continue; }  // wall + 4 rings of a fixed side
        // The interface guard applies to COARSENING only: splitting keeps the (required) interface
        // edge and only adds resolution. Freeable to coarsen only if this element and all its node
        // neighbours carry the same single mattype.
        if (fossil_fine) {
            bool on_interface = (elem_mat[e] < 0);
            for (int i = 0; i < NODES_PER_ELEM && !on_interface; ++i) {
                const int npatch = var.support.size(conn[i]);
                const int* patch = var.support.patch(conn[i]);
                for (int k = 0; k < npatch; ++k)
                    if (elem_mat[patch[k]] != elem_mat[e]) { on_interface = true; break; }
            }
            if (on_interface) { ++n_mat_skipped; continue; }
        }
        if (m.elem_grade[e] != ELEM_FREE) { if (fossil_fine) ++n_freed; else ++n_oversized; }
        m.elem_grade[e] = ELEM_FREE;
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            // Never unpin a wall (clamp_side_bits) or bottom node: only the element's off-wall,
            // off-bottom nodes become movable.
            if (!((*var.bcflag)[conn[i]] & clamp_side_bits) && !((*var.bcflag)[conn[i]] & BOUNDZ0))
                m.node_movable[conn[i]] = 1;
    }
    if (n_freed || n_oversized || n_top_skipped || n_mat_skipped || n_fixed_skipped)
        std::cout << "    Size recovery: freed " << n_freed
                  << " fossil-fine + " << n_oversized << " oversized(>=2x) element(s) (skipped "
                  << n_top_skipped << " top-surface-connected, "
                  << n_mat_skipped << " on a material interface, "
                  << n_fixed_skipped << " in a fixed-side band).\n";
}


// Rules R1-R5 on the OLD mesh. The masks are consumed by mark_quiet_required, which applies
// R6/R7 and derives MMG's required arrays.
FreezeMask compute_active_mask(const Param &param, const Variables &var,
                               const array_t &coord, const conn_t &connectivity,
                               int nnode, int nelem)
{
    FreezeMask m;
    m.node_movable.assign(nnode, 0);
    m.elem_grade.assign(nelem, ELEM_FROZEN);

    rule_repair_and_refine(param, var, coord, connectivity, nelem, m);   // R1 + R2

    // DIAG: per-rule snapshots of the node mask, consumed by the top-row report below
    std::vector<char> diag_mv_r1 = m.node_movable;   // movable after R1/R2

    // Side-wall masks of the restored walls (option 13):
    //   clamp_side_bits  = ALL restored sides (R3 pins, R5 protects; matches the metric clamp)
    //   fixed_side_bits  = the NON-inflow sides (R5 leaves their wall band verbatim)
    //   inflow_side_bits = the sides receiving material (gates R4)
    SideInflowState sis;
    uint  clamp_side_bits = 0, fixed_side_bits = 0, inflow_side_bits = 0;
    if (param.mesh.remeshing_option == 13) {
        clamp_side_bits = BOUNDX0 | BOUNDX1;
        sis = side_inflow_state(param, nnode, coord, *var.bcflag);
        if (!sis.in_x0) fixed_side_bits |= BOUNDX0; else inflow_side_bits |= BOUNDX0;
        if (!sis.in_x1) fixed_side_bits |= BOUNDX1; else inflow_side_bits |= BOUNDX1;
#ifdef THREED
        clamp_side_bits |= BOUNDY0 | BOUNDY1;
        if (!sis.in_y0) fixed_side_bits |= BOUNDY0; else inflow_side_bits |= BOUNDY0;
        if (!sis.in_y1) fixed_side_bits |= BOUNDY1; else inflow_side_bits |= BOUNDY1;
#endif
    }

    rule_boundary_nodes(param, var, nnode, clamp_side_bits, m);          // R3

    // DIAG: movable after R3
    std::vector<char> diag_mv_r3 = m.node_movable;

    if (param.mesh.remeshing_option == 13 && inflow_side_bits)
        rule_inflow_side_strip(var, connectivity, nnode, nelem, inflow_side_bits, m);  // R4

    // DIAG: movable after R4 (for the DES_DEBUG_MASK attribution below)
    std::vector<char> diag_mv_r4_snapshot = m.node_movable;
    const char *diag_mv_r4_set = diag_mv_r4_snapshot.data();

    // DIAG top-row report: attribute, for the TOP-row entities, which rule removed their
    // required status. This is how the whole-top-row re-tessellation was pinned on R3
    // (surface_churn_trace.py / compare_surface_churn_fix.py in the subd-serp benchmark);
    // it stays as the cheap per-remesh telemetry for the surface freeze. Element
    // attribution follows derive_required: FROZEN + all nodes pinned = required; anything
    // else is rebuildable by MMG.
    {
        int tn = 0, tmv = 0, t_r1 = 0, t_r3 = 0, t_r4 = 0;
        for (int n = 0; n < nnode; ++n) {
            if (!((*var.bcflag)[n] & BOUNDZ1)) continue;
            ++tn;
            if (!m.node_movable[n]) continue;
            ++tmv;
            if (diag_mv_r1[n]) ++t_r1;
            else if (diag_mv_r3[n]) ++t_r3;
            else ++t_r4;
        }
        int te = 0, n_free = 0, n_split = 0, n_frozen = 0, n_req = 0;
        int blk_r1 = 0, blk_r3 = 0, blk_r4 = 0;
        for (int e = 0; e < nelem; ++e) {
            ConstConnAccessor conn = connectivity[e];
            bool top = false;
            for (int i = 0; i < NODES_PER_ELEM && !top; ++i)
                top = ((*var.bcflag)[conn[i]] & BOUNDZ1) != 0;
            if (!top) continue;
            ++te;
            if (m.elem_grade[e] == ELEM_FREE) { ++n_free; continue; }
            if (m.elem_grade[e] == ELEM_SPLIT) { ++n_split; continue; }
            ++n_frozen;
            bool all_fixed = true;
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                if (m.node_movable[conn[i]]) all_fixed = false;
            if (all_fixed) { ++n_req; continue; }
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                const int n = conn[i];
                if (!m.node_movable[n]) continue;
                if (diag_mv_r1[n]) ++blk_r1;
                else if (diag_mv_r3[n]) ++blk_r3;
                else ++blk_r4;
            }
        }
        std::cout << "    [DIAG top-bdry nodes] total=" << tn << " movable=" << tmv
                  << " (R1=" << t_r1 << " R3=" << t_r3 << " R4=" << t_r4 << ")\n"
                  << "    [DIAG top-row elems]  total=" << te
                  << " free(R1)=" << n_free << " split(R2)=" << n_split
                  << " frozen=" << n_frozen << " -> required=" << n_req
                  << "; blocking movable-node refs: R1=" << blk_r1
                  << " R3=" << blk_r3 << " R4=" << blk_r4 << "\n";
    }

    // R5 wall-clamp targets: ALL restored sides (the metric clamp is not inflow-gated). Runs after
    // the DIAG report so the R1-R4 attribution stays pure.
    rule_size_recovery(param, var, connectivity, nnode, nelem, clamp_side_bits,
                       fixed_side_bits, m);  // R5

    // DES_DEBUG_MASK: per-rule attribution of every movable node (1=R1/R2, 3=R3, 4=R4,
    // 5=R5; env-gated diagnostic, appended per remesh).
    if (const char *fn = std::getenv("DES_DEBUG_MASK")) {
        std::FILE *fp = std::fopen(fn, "a");
        std::fprintf(fp, "# mask nnode=%d nelem=%d\n", nnode, nelem);
        for (int n = 0; n < nnode; ++n) {
            if (!m.node_movable[n]) continue;
            int rule = 5;
            if (diag_mv_r1[n]) rule = 1;
            else if (diag_mv_r3[n]) rule = 3;
            else if (diag_mv_r4_set && diag_mv_r4_set[n]) rule = 4;
            std::fprintf(fp, "N %d %.3f %.3f %u %d\n", n, coord[n][0], coord[n][NDIMS-1],
                         (*var.bcflag)[n], rule);
        }
        std::fclose(fp);
    }
    return m;
}

// DERIVED REQUIRED: required_node = pinned; required_elem = FROZEN grade with all nodes
// pinned. An element with any movable node is left non-required (the one-element hinge/
// transition layer), so MMG can adjust it when the neighbour moves.
void derive_required(const FreezeMask &m, int mnode, int melem, const int *conn,
                     std::vector<char> &required_node, std::vector<char> &required_elem,
                     int &n_req_node, int &n_req_elem)
{
    required_node.assign(mnode, 0);
    required_elem.assign(melem, 0);
    n_req_node = n_req_elem = 0;
    for (int n = 0; n < mnode; ++n)
        if (!m.node_movable[n]) { required_node[n] = 1; ++n_req_node; }
    for (int e = 0; e < melem; ++e) {
        if (m.elem_grade[e] != ELEM_FROZEN) continue;   // repair/refine element: MMG may split/remesh it
        const int *el = conn + (std::size_t)e*NODES_PER_ELEM;
        bool all_fixed = true;
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            if (m.node_movable[el[i]]) { all_fixed = false; break; }
        if (all_fixed) { required_elem[e] = 1; ++n_req_elem; }
    }
}


// ----------------------------------------------------------------------------
//  R6 input: what did flatten_* break?
// ----------------------------------------------------------------------------
// Elements flatten_* broke: a boundary node was snapped onto its restored plane (qcoord !=
// original_coord) and the post-flatten quality fell below min_quality. Also records the
// boundary-reshaped old nodes in var.remesh_affected_old_node. OLD element ids.
std::vector<char> detect_flatten_broken(const Param &param, Variables &var,
                                        const array_t &original_coord,
                                        const double *qcoord, const int *qconn,
                                        int old_nnode, int old_nelem,
                                        const int_vec &points_to_delete)
{
    std::vector<char> displaced(old_nnode, 0);
    for (int n = 0; n < old_nnode; ++n)
        for (int d = 0; d < NDIMS; ++d)
            if (qcoord[n*NDIMS + d] != original_coord[n][d]) { displaced[n] = 1; break; }
    // Boundary-reshaped original nodes: moved by flatten_* plus deleted/collapsed
    // (points_to_delete); barycentric_node_interpolation silences its warning there.
    if ((int)var.remesh_affected_old_node.size() == old_nnode) {
        for (int n = 0; n < old_nnode; ++n)
            if (displaced[n]) var.remesh_affected_old_node[n] = 1;
        for (int n : points_to_delete)
            if (n >= 0 && n < old_nnode) var.remesh_affected_old_node[n] = 1;
    }
    std::vector<char> flatten_broken(old_nelem, 0);
    for (int e = 0; e < old_nelem; ++e) {
        const int *el = qconn + e*NODES_PER_ELEM;
        bool touched = false;
        for (int i = 0; i < NODES_PER_ELEM; ++i) if (displaced[el[i]]) { touched = true; break; }
        if (!touched) continue;
        double q = packed_elem_quality(qcoord, el);
#ifdef THREED
        if (q > 0.0) q = std::cbrt(q);   // match rule_repair_and_refine's normalization
#endif
        if (q < param.mesh.min_quality) flatten_broken[e] = 1;
    }
    return flatten_broken;
}


// ----------------------------------------------------------------------------
//  R6/R7 + derived required: the required arrays handed to MMG
// ----------------------------------------------------------------------------

// Renumbering info for the collapse path (collapse_outside_nodes /
// collapse_short_boundary_segments_2d rebuilt the mesh handed to MMG, so old and MMG entity
// ids no longer match 1:1). All arrays describe the NEW (collapsed) mesh except
// collapsed_pts, which lists the OLD node ids that were merged away.
struct CollapseRemap {
    const int_vec *new_to_old_node;   // [mnode] -> old node id
    const int_vec *new_to_old_elem;   // [melem] -> old elem id
    const int     *new_conn;          // packed melem * NODES_PER_ELEM
    const int_vec *collapsed_pts;     // old node ids removed by the collapse
};

// Conservative remeshing: mark the "quiet" part of the mesh as MMG-required so MMG leaves it
// untouched and only remeshes the active region. This mirrors the triangle path, which
// re-triangulates existing points and does not disturb elements without plastic strain or
// deformation -- avoiding needless node motion and field interpolation in quiet regions.
//
// The criterion is plastic strain that is still INCREASING, not accumulated plastic strain.
// Every MMG remesh re-interpolates the elements it re-adapts (inject_field's weighted average),
// which numerically diffuses their fields; keying on accumulated strain would re-adapt the whole
// shear band every remesh and smear the localization (peak plstrain decays vs the triangle path).
// Instead we compare against plstrain_remesh (the snapshot taken at the previous remesh): an
// element is ACTIVE only if it yielded MORE since then (R2). A fossil band that has stopped
// growing is therefore frozen and its sharp plstrain is carried across the remesh verbatim
// (is_changed==0 -> direct copy), exactly like the triangle path preserves untouched points.
//
// MMG entity ids are 1-based and, on the direct path (remap == NULL), match the old-mesh
// node/element order 1:1 (vertices = packed old_coord, triangles/tets built from old
// connectivity in order). On the collapse path the mesh handed to MMG has been RENUMBERED
// (sunk material merged out, elements dropped), so the old<->MMG id map is no longer 1:1:
//   (1) compute the freeze policy on the OLD mesh (its var fields match old ids exactly),
//   (2) TRANSLATE the masks through the collapse renumbering (remap->new_to_old_*),
//   (3) additionally FREE the collapse region (R6): any new element whose source old element
//       referenced a collapsed (pts) node was reshaped by the collapse; mark it and its nodes
//       movable, plus one ring of connected elements. Those freed nodes carry the
//       init_elem_size_n-based metric (compute_metric_field's base), so MMG re-refines the
//       collapsed region back to the initial element size instead of leaving it coarse,
//   (4) derive required on the NEW connectivity.
//
// R6 on both paths: flatten_* may have snapped boundary nodes onto their restored plane
// (up to max_boundary_distortion), squashing or inverting elements AFTER the R1-R5 masks were
// computed on the PRE-flatten coords; their interior nodes would stay pinned as required
// vertices and MMG could never repair the flat band it would otherwise emit along the restored
// boundary (cmp_bottom remesh 176/177). flatten_broken[e] applies the same REPAIR rule on the
// post-flatten geometry; free those elements exactly like R1 frees a repair element.
// Deliberately NO wider freeing: every extra freed element enlarges the re-interpolated region
// and its post-remesh disequilibrium shock (freeing the whole flatten band NaN'd cmp_bottom at
// step 118600; dies at: whole displaced band 118600 / quality-halved band 99200 / REPAIR-only
// 177400 = the model's physical neck-through end).
//
// NOTE (2026-07-10): dissolving small REQUIRED-element islands (components < 6 elements
// losing required-element status, nodes kept) was tried here as a preventive against the
// post-remesh quality retries and REVERTED after an A/B/C/D comparison on the
// subd-serp-remesh frame-30 restart window: with dissolution 14 remeshes / 17 MMG re-runs,
// without it 10 / 8 -- re-tessellating the dissolved islands at every remesh CREATED more
// bad output elements (and mesh churn) than the pinned islands ever did. Do not re-add.
void mark_quiet_required(const Param &param, const Variables &var,
                         const array_t &old_coord, const conn_t &old_connectivity,
                         int mnode, int melem,
                         const std::vector<char> &flatten_broken,
                         const CollapseRemap *remap,
                         std::vector<char> &required_node, std::vector<char> &required_elem)
{
    const int old_nnode = old_coord.size();
    const int old_nelem = old_connectivity.size();

    // R1-R5 on the OLD mesh (var's fields match old ids exactly).
    FreezeMask m = compute_active_mask(param, var, old_coord, old_connectivity,
                                       old_nnode, old_nelem);

    int n_req_node = 0, n_req_elem = 0;
    if (!remap) {
        // Direct path: old and MMG ids match 1:1 (mnode == old_nnode, melem == old_nelem).
        // R6: free what flatten broke.
        for (int e = 0; e < melem; ++e) {
            if (!flatten_broken[e]) continue;
            m.elem_grade[e] = ELEM_FREE;
            ConstConnAccessor conn = old_connectivity[e];
            for (int i = 0; i < NODES_PER_ELEM; ++i) m.node_movable[conn[i]] = 1;
        }
        // Array2D is SoA -- .data() must not be row-indexed (see free_wall_band's note).
        // Materialize a row-major copy for derive_required, which shares its signature
        // with the collapse path's genuinely row-major plain buffer.
        std::vector<int> conn_rows((std::size_t)melem * NODES_PER_ELEM);
        for (int e = 0; e < melem; ++e) {
            ConstConnAccessor cc = old_connectivity[e];
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                conn_rows[(std::size_t)e*NODES_PER_ELEM + i] = cc[i];
        }
        derive_required(m, mnode, melem, conn_rows.data(),
                        required_node, required_elem, n_req_node, n_req_elem);
        std::cout << "    Conservative remesh: froze " << n_req_elem << "/" << melem
                  << " elements (" << (melem - n_req_elem) << " modifiable), "
                  << n_req_node << "/" << mnode << " nodes required.\n";
        // DES_DEBUG_MASK=<file>: dump this remesh's freeze mask for offline attribution
        // (node id, x, z, bcflag, movable, required; element grade, required). Env-gated
        // diagnostic, direct path only.
        if (const char *fn = std::getenv("DES_DEBUG_MASK")) {
            std::FILE *fp = std::fopen(fn, "a");
            std::fprintf(fp, "# remesh nnode=%d nelem=%d\n", mnode, melem);
            for (int n = 0; n < mnode; ++n)
                std::fprintf(fp, "N %d %.6f %.6f %u %d %d\n", n,
                             old_coord[n][0], old_coord[n][NDIMS-1], (*var.bcflag)[n],
                             (int)m.node_movable[n], (int)required_node[n]);
            for (int e = 0; e < melem; ++e)
                std::fprintf(fp, "E %d %d %d\n", e, (int)m.elem_grade[e], (int)required_elem[e]);
            std::fclose(fp);
        }
    }
    else {
        const int *c_conn = remap->new_conn;

        // (2) translate OLD masks -> NEW (collapsed) masks
        FreezeMask t;
        t.node_movable.assign(mnode, 0);
        t.elem_grade.assign(melem, ELEM_FROZEN);
        for (int n = 0; n < mnode; ++n) t.node_movable[n] = m.node_movable[(*remap->new_to_old_node)[n]];
        for (int e = 0; e < melem; ++e) t.elem_grade[e]   = m.elem_grade[(*remap->new_to_old_elem)[e]];

        // (3) R7 + R6: free the collapse region (elements whose source old element touched a
        //     collapsed node, + their nodes, + one ring) and the elements flatten BROKE (looked up
        //     through the renumbering): otherwise their pinned interior nodes leave MMG unable to
        //     repair the flat band along the restored boundary.
        std::vector<char> is_pts(old_nnode, 0);
        for (std::size_t i = 0; i < remap->collapsed_pts->size(); ++i)
            is_pts[(*remap->collapsed_pts)[i]] = 1;
        std::vector<char> collapse_node(mnode, 0);
        for (int e = 0; e < melem; ++e) {
            const int old_e = (*remap->new_to_old_elem)[e];
            ConstConnAccessor oc = old_connectivity[old_e];
            bool touched = flatten_broken[old_e] != 0;
            for (int i = 0; i < NODES_PER_ELEM && !touched; ++i)
                if (is_pts[oc[i]]) touched = true;
            if (!touched) continue;
            t.elem_grade[e] = ELEM_FREE;
            for (int i = 0; i < NODES_PER_ELEM; ++i) { int nn = c_conn[e*NODES_PER_ELEM + i]; t.node_movable[nn] = 1; collapse_node[nn] = 1; }
        }
        for (int e = 0; e < melem; ++e) {
            bool touch = false;
            for (int i = 0; i < NODES_PER_ELEM; ++i) if (collapse_node[c_conn[e*NODES_PER_ELEM + i]]) { touch = true; break; }
            if (!touch) continue;
            t.elem_grade[e] = ELEM_FREE;
            for (int i = 0; i < NODES_PER_ELEM; ++i) t.node_movable[c_conn[e*NODES_PER_ELEM + i]] = 1;
        }
        // (4) derive required on the NEW connectivity
        derive_required(t, mnode, melem, c_conn,
                        required_node, required_elem, n_req_node, n_req_elem);
        std::cout << "    Conservative remesh (collapse path): froze " << n_req_elem << "/" << melem
                  << " elements (" << (melem - n_req_elem) << " modifiable), "
                  << n_req_node << "/" << mnode << " nodes required.\n";
    }
}


// ----------------------------------------------------------------------------
//  Post-output quality gate (retry with graded unfreeze)
// ----------------------------------------------------------------------------

// Centroid + longest-edge length of a packed-coords element (locality scale for the
// post-remesh unfreeze below).
void elem_centroid_scale(const double *coord, const int *elem, double *centroid, double &scale)
{
    for (int d = 0; d < NDIMS; ++d) centroid[d] = 0.0;
    scale = 0.0;
    for (int k = 0; k < NODES_PER_ELEM; ++k) {
        const double *a = coord + (std::size_t)elem[k]*NDIMS;
        for (int d = 0; d < NDIMS; ++d) centroid[d] += a[d];
        for (int j = 0; j < k; ++j) {
            const double *b = coord + (std::size_t)elem[j]*NDIMS;
            double d2 = 0.0;
            for (int d = 0; d < NDIMS; ++d) { double dd = a[d]-b[d]; d2 += dd*dd; }
            if (d2 > scale) scale = d2;
        }
    }
    for (int d = 0; d < NDIMS; ++d) centroid[d] /= NODES_PER_ELEM;
    scale = std::sqrt(scale);
}

// Collect the below-min_quality / tiny elements of an MMG output (the same criteria as
// bad_mesh_quality's element checks -- an output that would immediately re-trigger remeshing).
// Returns the count; when bad_c/bad_r are non-null the bad elements' centroids and locality
// radii (longest edge) are appended for the unfreeze region seeding below.
int collect_bad_output(const Param &param, const MMGOutput &out,
                       std::vector<double> *bad_c, std::vector<double> *bad_r)
{
    const double q_thr = param.mesh.min_quality;
    // A fresh output element must clear the tiny-element trigger by HALF the hysteresis gap
    // (> smallest_vol / sqrt(margin)): at the trigger itself a constrained spot re-crosses within
    // a few thousand steps; at the full floor MMG's output scatter routinely fails.
    const double tiny_vol = refine_floors(param.mesh).smallest_vol
                            / std::sqrt(param.mesh.remesh_tiny_margin);
    int nbad = 0;
    for (int e = 0; e < out.nelem; ++e) {
        const int *el = &out.conn[(std::size_t)e*NODES_PER_ELEM];
        const double m = signed_elem_measure(out.coord.data(), el);
#ifdef THREED
        const double vol = m / 6.0;
#else
        const double vol = m / 2.0;
#endif
        double q = packed_elem_quality(out.coord.data(), el);
#ifdef THREED
        if (q > 0.0) q = std::cbrt(q);   // bad_mesh_quality's normalization
#endif
        if (q >= q_thr && vol >= tiny_vol) continue;
        ++nbad;
        if (!bad_c || !bad_r) continue;
        double c[NDIMS], r;
        elem_centroid_scale(out.coord.data(), el, c, r);
        for (int d = 0; d < NDIMS; ++d) bad_c->push_back(c[d]);
        bad_r->push_back(r);
    }
    return nbad;
}

// Post-remesh quality gate: bad_mesh_quality's element checks applied to the MMG OUTPUT, so a
// remesh cannot emit an element that immediately re-triggers remeshing (MMG repairs the trigger
// element, frozen neighbours block the fix, the badness moves next door). On failure, unfreeze
// the still-required input entities around each bad output element with GRADED freedom:
//   * CORE  -- the bad element's input pre-image grown to frozen contact (capped) + `attempt`
//     extra rings, so a persistent failure gets more room instead of being chased one element
//     per retry: free req_elem AND req_node.
//   * HINGE -- one more ring: free req_elem only, keep req_node (split/swap allowed, nodal fields
//     still carried verbatim).
// Returns the number unfreed; 0 = clean output, or nothing frozen in reach.
int unfreeze_near_bad_output(const Param &param, const MMGOutput &out,
                             int mnode, int melem, const double *mcoord, const int *mconn,
                             int attempt,
                             std::vector<char> &req_node, std::vector<char> &req_elem)
{
    std::vector<double> bad_c;    // bad output element centroids (NDIMS each)
    std::vector<double> bad_r;    // matching locality radii (longest edge)
    const int nbad = collect_bad_output(param, out, &bad_c, &bad_r);
    if (nbad == 0) return 0;

    // Node -> element support of the MMG INPUT mesh, built locally: on the collapse path the mesh
    // is renumbered and var.support does not match mconn. Centroids/radii are cached alongside.
    std::vector<std::vector<int>> node_elems(mnode);
    std::vector<double> ic((std::size_t)melem * NDIMS);   // input-element centroids
    std::vector<double> ir(melem);                        // input-element radii (longest edge)
    for (int e = 0; e < melem; ++e) {
        const int *el = mconn + (std::size_t)e*NODES_PER_ELEM;
        elem_centroid_scale(mcoord, el, &ic[(std::size_t)e*NDIMS], ir[e]);
        for (int k = 0; k < NODES_PER_ELEM; ++k) node_elems[el[k]].push_back(e);
    }

    // Seed: input elements overlapping a bad output element (centroids within the summed radii),
    // plus each bad element's single nearest input element so the seed is never empty.
    std::vector<char> region(melem, 0);
    std::vector<int> frontier;
    std::vector<int> nearest(nbad, -1);
    std::vector<double> nearest_d2(nbad, 1e300);
    for (int e = 0; e < melem; ++e) {
        const double *c = &ic[(std::size_t)e*NDIMS];
        for (int j = 0; j < nbad; ++j) {
            double d2 = 0.0;
            for (int d = 0; d < NDIMS; ++d) { double dd = c[d] - bad_c[(std::size_t)j*NDIMS + d]; d2 += dd*dd; }
            if (d2 < (bad_r[j] + ir[e]) * (bad_r[j] + ir[e]) && !region[e]) { region[e] = 1; frontier.push_back(e); }
            if (d2 < nearest_d2[j]) { nearest_d2[j] = d2; nearest[j] = e; }
        }
    }
    for (int j = 0; j < nbad; ++j)
        if (nearest[j] >= 0 && !region[nearest[j]]) { region[nearest[j]] = 1; frontier.push_back(nearest[j]); }

    // Does the current region touch any frozen (still-required) entity?
    auto region_has_frozen = [&]() {
        for (int e = 0; e < melem; ++e) {
            if (!region[e]) continue;
            if (req_elem[e]) return true;
            const int *el = mconn + (std::size_t)e*NODES_PER_ELEM;
            for (int k = 0; k < NODES_PER_ELEM; ++k) if (req_node[el[k]]) return true;
        }
        return false;
    };

    // Grow one node-support ring: every element sharing a node with the current frontier.
    auto grow_ring = [&](std::vector<int> &fr) {
        std::vector<int> next;
        for (int e : fr) {
            const int *el = mconn + (std::size_t)e*NODES_PER_ELEM;
            for (int k = 0; k < NODES_PER_ELEM; ++k)
                for (int nb : node_elems[el[k]])
                    if (!region[nb]) { region[nb] = 1; next.push_back(nb); }
        }
        fr.swap(next);
    };

    // CORE growth. Phase A: grow to the first frozen entity (capped), so a blocking frozen band a
    // few elements away is reached. Phase B: `attempt` extra rings per retry.
    const int max_contact_rings = 3;   // matches brc-interpolation's BFS layer cap
    int rings = 0;
    while (!frontier.empty() && rings < max_contact_rings && !region_has_frozen()) {
        grow_ring(frontier);
        ++rings;
    }
    for (int r = 0; r < attempt && !frontier.empty(); ++r) {
        grow_ring(frontier);
        ++rings;
    }

    // Free the CORE: both entity kinds -- full repair freedom.
    int nfreed_core = 0;
    for (int e = 0; e < melem; ++e) {
        if (!region[e]) continue;
        const int *el = mconn + (std::size_t)e*NODES_PER_ELEM;
        if (req_elem[e]) { req_elem[e] = 0; ++nfreed_core; }
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            if (req_node[el[k]]) { req_node[el[k]] = 0; ++nfreed_core; }
    }

    // Required-element invariant: an element with any movable node cannot be required. This IS
    // the hinge ring: required elements touching the core lose element status, nodes stay pinned.
    // (An extra explicitly-freed ring beyond it measured as an exact no-op.)
    int nfreed_hinge = 0;
    for (int e = 0; e < melem; ++e) {
        if (!req_elem[e]) continue;
        const int *el = mconn + (std::size_t)e*NODES_PER_ELEM;
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            if (!req_node[el[k]]) { req_elem[e] = 0; ++nfreed_hinge; break; }
    }
    const int nfreed = nfreed_core + nfreed_hinge;
    std::cout << "    Post-remesh quality check: " << nbad
              << " output element(s) below min_quality/tiny; unfroze " << nfreed_core
              << " core entities within " << rings << " connectivity ring(s) (+"
              << nfreed_hinge << " hinge element(s), nodes kept)"
              << (nfreed ? "." : " -- nothing frozen in reach, keeping this mesh.")
              << "\n";
    return nfreed;
}

// Quality-gated adaptation: if the output would immediately re-trigger remeshing
// (below-min_quality or tiny element), unfreeze the frozen entities around the bad
// spots and re-run MMG. req_node/req_elem MUST be the vectors mmg_in.required_node/
// required_elem point into: the retry unfreeze mutates them in place and re-runs MMG
// on the same input. Bounded retries; a failed final attempt keeps the last mesh.
void mmg_adapt_quality_gated(const Param &param, const MMGInput &in, MMGOutput &out,
                             std::vector<char> &req_node, std::vector<char> &req_elem)
{
    for (int attempt = 0; ; ++attempt) {
        mmg_adapt(param.mesh, in, out);
        if (attempt >= 3) {   // 4 MMG passes: heuristic cap, untuned
            // Retries exhausted: keep the last mesh, but say so when it is still bad --
            // a silent bad element re-triggers remeshing a few steps later.
            const int nbad = collect_bad_output(param, out, NULL, NULL);
            if (nbad)
                std::cout << "    Warning: post-remesh quality retries exhausted; keeping a mesh "
                             "with " << nbad << " below-min_quality/tiny element(s).\n";
            break;
        }
        if (unfreeze_near_bad_output(param, out, in.nnode, in.nelem, in.coord, in.conn,
                                     attempt, req_node, req_elem) == 0) break;
        std::cout << "    Re-running MMG with the unfrozen neighbourhood (attempt "
                  << attempt + 2 << ").\n";
    }
}

// ======================= end FREEZE / REQUIRED POLICY =======================

// MMG-native handling of material that moved OUTSIDE a restored boundary. flatten_* snapped the
// boundary nodes back onto their planes and collected the nodes left on the far side (`pts`);
// those would invert the adjoining elements. Remove them by edge collapse instead of falling
// back to a whole-domain Triangle/Tetgen remesh:
//   * CORNER-AWARE: each outside node records which plane(s) it crossed (BOUND* bitmask) and
//     only collapses onto a boundary node of a plane it crossed, never the adjacent top/bottom.
//   * ITERATIVE: each pass merges only the outside nodes currently touching such a boundary node,
//     peeling deep stacks one layer at a time (mapping a deep node straight to the boundary tangled).
// Builds fresh 0-indexed output arrays; inputs untouched. Survivor node order is preserved so the
// metric maps by renumbering; outside nodes are interior, so segments only need renumbering.
void collapse_outside_nodes(const int_vec &pts, int nnode, int nelem, int nseg,
                            const double *coord, const uint_vec &bcflag,
                            const int *conn, const int *segment, const int *segflag,
                            const double_vec &metric,
                            double zlen, double xlen, double ylen, double min_dist,
                            double_vec &ncoord, int_vec &nconn,
                            int_vec &nsegment, int_vec &nsegflag, double_vec &nmetric,
                            int_vec &new_to_old_node, int_vec &new_to_old_elem,
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
        crossed[s] = c ? c : BOUND_ANY;   // safety: allow any boundary if none matched
    }

    // merged_to[n] = the (alive, boundary) node n was collapsed onto, or -1 if still alive.
    // Targets are boundary nodes (never in pts -> never merge), so chains have depth 1.
    int_vec merged_to(nnode, -1);
    auto alive   = [&](int n) { return merged_to[n] < 0; };
    auto resolve = [&](int n) { return merged_to[n] >= 0 ? merged_to[n] : n; };

    auto dist2 = [&](int s, int m) {
        double d2 = 0.0;
        for (int d = 0; d < NDIMS; ++d) { double dd = coord[s*NDIMS+d] - coord[m*NDIMS+d]; d2 += dd*dd; }
        return d2;
    };

    // ITERATIVE peel: each pass merges every outside node touching a boundary node of a crossed
    // plane onto the nearest such node; converges in ~(number of outside layers) passes.
    for (int pass = 0; pass < nnode; ++pass) {
        int_vec    tgt(nnode, -1);
        double_vec best(nnode, std::numeric_limits<double>::max());
        for (int e = 0; e < nelem; ++e) {
            const int *el = conn + e*NODES_PER_ELEM;
            for (int a = 0; a < NODES_PER_ELEM; ++a) {
                int s = resolve(el[a]);
                if (!crossed[s] || !alive(s)) continue;          // only still-outside nodes
                for (int b = 0; b < NODES_PER_ELEM; ++b) {
                    int m = resolve(el[b]);
                    if (m == s || !alive(m)) continue;
                    if (is_boundary(bcflag[m]) && (bcflag[m] & crossed[s])) {  // on a crossed plane
                        double d2 = dist2(s, m);
                        if (d2 < best[s]) { best[s] = d2; tgt[s] = m; }
                    }
                }
            }
        }
        bool progress = false;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            int s = pts[i];
            if (alive(s) && crossed[s] && tgt[s] >= 0) { merged_to[s] = tgt[s]; progress = true; }
        }
        if (!progress) break;
    }

    // Stalemate: any outside node still not reachable from its crossed plane through the mesh
    // graph -> merge onto the nearest boundary node on a crossed plane, globally (rare).
    for (std::size_t i = 0; i < pts.size(); ++i) {
        int s = pts[i];
        if (!alive(s) || !crossed[s]) continue;
        double best = std::numeric_limits<double>::max(); int tg = -1;
        for (int m = 0; m < nnode; ++m) {
            if (!alive(m) || !is_boundary(bcflag[m]) || !(bcflag[m] & crossed[s])) continue;
            double d2 = dist2(s, m);
            if (d2 < best) { best = d2; tg = m; }
        }
        if (tg >= 0) merged_to[s] = tg;
    }

#ifdef DEBUG_COLLAPSE
    {
        int merged = 0, unmerged = 0;
        std::fprintf(stderr, "[collapse] pts=%zu\n", pts.size());
        for (std::size_t i = 0; i < pts.size(); ++i) {
            int s = pts[i];
            if (!alive(s)) { ++merged; continue; }
            ++unmerged;
            if (unmerged <= 12)
                std::fprintf(stderr, "[collapse]  UNMERGED node %d: x=%.1f z=%.1f crossed=%u bc=%u\n",
                             s, coord[s*NDIMS], coord[s*NDIMS+NDIMS-1], crossed[s], (uint)bcflag[s]);
        }
        std::fprintf(stderr, "[collapse] merged=%d unmerged=%d\n", merged, unmerged);
    }
#endif

    // Finish merging STRANDED nodes: alive but with zero surviving elements. A sunk blob collapsing
    // onto its plane can leave a boundary node whose triangles all degenerate while its boundary
    // segments survive; MMG would keep it as an isolated 0-mass vertex -> NaN. Merge each onto its
    // nearest alive neighbour on the shared plane (the dangling segment collapses with it);
    // iterate, since absorbing one node can strand the next.
    {
        // Reference measure (orientation + scale) from a pristine element -- same value the
        // connectivity rebuild uses, computed on `coord` so the flat/inverted test matches exactly.
        double refm = 0.0;
        for (int e = 0; e < nelem && refm == 0.0; ++e) {
            const int *el = conn + e*NODES_PER_ELEM;
            bool pristine = true;
            for (int k = 0; k < NODES_PER_ELEM; ++k)
                if (crossed[el[k]] || !alive(el[k])) { pristine = false; break; }
            if (pristine) refm = signed_elem_measure(coord, el);
        }
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

    // renumber survivors, rebuild coord + metric. new_to_old_node[new] = the survivor's old id,
    // so the caller can look up per-node old-mesh state (bcflag, freeze masks) after renumbering.
    int_vec node_map(nnode, -1);
    ncoord.clear(); nmetric.clear(); new_to_old_node.clear();
    nn = 0;
    for (int n = 0; n < nnode; ++n) {
        if (!alive(n)) continue;
        node_map[n] = nn++;
        new_to_old_node.push_back(n);
        for (int d = 0; d < NDIMS; ++d) ncoord.push_back(coord[n*NDIMS + d]);
        nmetric.push_back(metric[n]);
    }
    auto remap = [&](int n) { return node_map[resolve(n)]; };

    // reference orientation/scale from a pristine element (no crossed node).
    double ref = 0.0;
    for (int e = 0; e < nelem && ref == 0.0; ++e) {
        const int *el = conn + e*NODES_PER_ELEM;
        bool pristine = true;
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            if (crossed[el[k]] || !alive(el[k])) { pristine = false; break; }
        if (!pristine) continue;
        int v[NODES_PER_ELEM];
        for (int k = 0; k < NODES_PER_ELEM; ++k) v[k] = node_map[el[k]];
        ref = signed_elem_measure(ncoord.data(), v);
    }

    // rebuild connectivity, dropping elements that collapse onto a boundary: a repeated node,
    // OR near-null / inverted measure (all nodes land on the flat boundary, or orientation flip).
    // new_to_old_elem[new] = the surviving element's old id (same map purpose as new_to_old_node).
    nconn.clear(); new_to_old_elem.clear();
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
        if (ref != 0.0 && m * ref <= 1e-6 * ref * ref) continue;
        for (int k = 0; k < NODES_PER_ELEM; ++k) nconn.push_back(v[k]);
        new_to_old_elem.push_back(e);
        ++ne;
    }

    // Rebuild boundary segments, dropping any collapsed to a point (an adjacent-boundary node
    // merged onto the corner); the neighbouring segment re-links to the merged endpoint.
    nsegment.clear(); nsegflag.clear();
    ns = 0;
    for (int s = 0; s < nseg; ++s) {
        int v[NODES_PER_FACET];
        bool degenerate = false;
        for (int k = 0; k < NODES_PER_FACET; ++k) {
            v[k] = remap(segment[s*NODES_PER_FACET + k]);
            for (int j = 0; j < k; ++j) if (v[j] == v[k]) degenerate = true;
        }
        if (degenerate) continue;
        for (int k = 0; k < NODES_PER_FACET; ++k) nsegment.push_back(v[k]);
        nsegflag.push_back(segflag[s]);
        ++ns;
    }
}

// True if any node lies beyond a boundary plane this remeshing option snaps back (bottom for
// 1/2/11/13; all sides for 13): material outside the fixed domain, to be collapsed before MMG.
// Mirrors exactly what flatten_* marks for deletion: non-bottom nodes below the base, and for a
// side every node NOT on that side's own plane that crossed it (adjacent-boundary nodes included).
// Excluding a plane's own wall nodes keeps genuine extension inert. Option 12 is excluded.
bool has_outside_material(const Param &param, const array_t &coord, const uint_vec &bcflag,
                          int nnode, double min_dist)
{
    const int opt = param.mesh.remeshing_option;
    const bool restore_bottom = (opt==1 || opt==2 || opt==11 || opt==13);
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

// Under strong convergence the free surface piles boundary nodes into a short near-vertical
// "cliff" (a boundary segment far below hmin). MMG will not collapse a boundary edge (its corner
// detection pins the ends; turning it off wrecks the box corners), so the sliver re-fires the
// tiny-element trigger forever. Collapse such segments BEFORE MMG by merging the removable
// endpoint onto its segment PARTNER (deterministic), like the Triangle path's boundary decimation.
// Rules: never remove a domain corner; each node joins at most one collapse per pass (depth-1);
// the connectivity rebuild drops degenerate AND near-null/inverted elements.
// AREA CONSERVATION: a plain snap sweeps the triangle (p, rem, keep) out of / into the domain.
// For pure TOP collapses the target is repositioned along n = rot90(q - p):
//     keep' = keep + (2 * signed_area(p, rem, keep) / |q - p|^2) * n
// so the shoelace contribution of p->keep'->q equals that of p->rem->keep->q. Falls back to the
// plain snap at corners / non-top segments / degenerate chords. 2D only. Returns the number of
// nodes collapsed (0 => outputs untouched). new_to_old_{node,elem}[new] = old id.
#ifndef THREED
int collapse_short_boundary_segments_2d(int nnode, int nelem, int nseg,
        const double *coord, const uint_vec &bcflag,
        const int *conn, const int *segment, const int *segflag, const double_vec &metric,
        double min_len,
        double_vec &ncoord, int_vec &nconn, int_vec &nsegment, int_vec &nsegflag, double_vec &nmetric,
        int_vec &new_to_old_node, int_vec &new_to_old_elem, int_vec &removed_old_nodes,
        int &nn, int &ne, int &ns)
{
    // node -> boundary-segment adjacency (a simple boundary chain has exactly 2 per node),
    // for locating the collapse neighbours p and q of the area-preserving reposition.
    int_vec nadj_of(nnode, 0), adj(2*nnode, -1);
    for (int s = 0; s < nseg; ++s)
        for (int k = 0; k < NODES_PER_FACET; ++k) {
            const int n = segment[s*NODES_PER_FACET + k];
            if (nadj_of[n] < 2) adj[2*n + nadj_of[n]] = s;
            ++nadj_of[n];
        }
    // other endpoint of boundary node n's OTHER segment (not segment s); -1 if not a chain
    auto chain_nbr = [&](int n, int s) {
        if (nadj_of[n] != 2) return -1;
        const int o = (adj[2*n] == s) ? adj[2*n + 1] : adj[2*n];
        const int e0 = segment[o*NODES_PER_FACET], e1 = segment[o*NODES_PER_FACET + 1];
        return e0 == n ? e1 : e0;
    };

    int_vec merged_to(nnode, -1);
    std::vector<char> touched(nnode, 0), moved(nnode, 0);
    double_vec moved_pos(2*nnode, 0.0);
    const double thr2 = min_len * min_len;
    int ncoll = 0, nadjusted = 0;
    double a_conserved = 0.0;
    for (int s = 0; s < nseg; ++s) {
        const int a = segment[s*NODES_PER_FACET], b = segment[s*NODES_PER_FACET + 1];
        if (touched[a] || touched[b]) continue;                 // keep every collapse independent
        double d2 = 0.0;
        for (int d = 0; d < NDIMS; ++d) { double dd = coord[a*NDIMS+d] - coord[b*NDIMS+d]; d2 += dd*dd; }
        if (d2 >= thr2) continue;
        const bool ca = is_corner(bcflag[a]), cb = is_corner(bcflag[b]);
        if (ca && cb) continue;                                 // never collapse between two corners
        const int rem = ca ? b : a, keep = ca ? a : b;          // keep the corner if there is one

        // area-preserving reposition of the merge target (pure top-surface chains only)
        const bool keep_pure_top =
            (bcflag[keep] & BOUNDZ1) && !(bcflag[keep] & (BOUND_ANY & ~BOUNDZ1));
        const uint sflag = static_cast<uint>(segflag[s]);
        if (keep_pure_top && (sflag & BOUNDZ1) && !(sflag & (BOUND_ANY & ~BOUNDZ1))) {
            const int p = chain_nbr(rem, s), q = chain_nbr(keep, s);
            if (p >= 0 && q >= 0 && p != q && p != keep && q != rem &&
                !touched[p] && !touched[q]) {
                const double *P = coord + p*NDIMS,   *R = coord + rem*NDIMS;
                const double *K = coord + keep*NDIMS, *Q = coord + q*NDIMS;
                // 2 * signed area of the triangle a plain snap would sweep
                const double A2 = (P[0]*R[1] - P[1]*R[0])
                                + (R[0]*K[1] - R[1]*K[0])
                                + (K[0]*P[1] - K[1]*P[0]);
                const double nx = Q[1] - P[1], ny = P[0] - Q[0];   // rot90(q - p)
                const double L2 = nx*nx + ny*ny;
                if (L2 > 1e-6 * thr2) {
                    const double t = A2 / L2;
                    moved[keep] = 1;
                    moved_pos[2*keep]     = K[0] + t*nx;
                    moved_pos[2*keep + 1] = K[1] + t*ny;
                    touched[p] = touched[q] = 1;   // a chained collapse would shift the chord
                    a_conserved += std::fabs(0.5 * A2);
                    ++nadjusted;
                }
            }
        }

        merged_to[rem] = keep;
        touched[a] = touched[b] = 1;
        removed_old_nodes.push_back(rem);
        ++ncoll;
    }
    if (ncoll == 0) { nn = ne = ns = 0; return 0; }
    if (nadjusted)
        std::cout << "    Area-preserving collapse: repositioned " << nadjusted
                  << " merge target(s), conserving " << a_conserved << " m^2.\n";
    auto resolve = [&](int n) { return merged_to[n] >= 0 ? merged_to[n] : n; };

    // renumber survivors, rebuild coord + metric (repositioned merge targets get their
    // area-conserving coordinates)
    int_vec node_map(nnode, -1);
    ncoord.clear(); nmetric.clear(); new_to_old_node.clear(); nn = 0;
    for (int n = 0; n < nnode; ++n) {
        if (merged_to[n] >= 0) continue;                        // removed -> merged onto its partner
        node_map[n] = nn++;
        new_to_old_node.push_back(n);
        for (int d = 0; d < NDIMS; ++d)
            ncoord.push_back(moved[n] ? moved_pos[n*NDIMS + d] : coord[n*NDIMS + d]);
        nmetric.push_back(metric[n]);
    }
    auto remap = [&](int n) { return node_map[resolve(n)]; };

    // reference orientation/scale from a pristine element (no merged or repositioned node)
    // to detect flips.
    double ref = 0.0;
    for (int e = 0; e < nelem && ref == 0.0; ++e) {
        const int *el = conn + e*NODES_PER_ELEM;
        bool pristine = true;
        for (int k = 0; k < NODES_PER_ELEM; ++k)
            if (merged_to[el[k]] >= 0 || moved[el[k]]) { pristine = false; break; }
        if (!pristine) continue;
        int v[NODES_PER_ELEM];
        for (int k = 0; k < NODES_PER_ELEM; ++k) v[k] = node_map[el[k]];
        ref = signed_elem_measure(ncoord.data(), v);
    }

    // rebuild connectivity, dropping repeated-node and near-null/inverted elements.
    nconn.clear(); new_to_old_elem.clear(); ne = 0;
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
        if (ref != 0.0 && m * ref <= 1e-6 * ref * ref) continue;
        for (int k = 0; k < NODES_PER_ELEM; ++k) nconn.push_back(v[k]);
        new_to_old_elem.push_back(e);
        ++ne;
    }

    // rebuild boundary segments, dropping any that collapsed to a point.
    nsegment.clear(); nsegflag.clear(); ns = 0;
    for (int s = 0; s < nseg; ++s) {
        int v[NODES_PER_FACET];
        bool degenerate = false;
        for (int k = 0; k < NODES_PER_FACET; ++k) {
            v[k] = remap(segment[s*NODES_PER_FACET + k]);
            for (int j = 0; j < k; ++j) if (v[j] == v[k]) degenerate = true;
        }
        if (degenerate) continue;
        for (int k = 0; k < NODES_PER_FACET; ++k) nsegment.push_back(v[k]);
        nsegflag.push_back(segflag[s]);
        ++ns;
    }
    return ncoll;
}
#endif
// Per-node ANISOTROPIC (tensor) metric for MMG from the scalar target: isotropic diag(1/h^2)
// everywhere except nodes on an INFLOW restored wall (option 13, sides in `sis`): fine (h)
// tangentially, coarse (h*ratio) wall-perpendicular, so the incoming column is carried by thin
// wide elements. Wall membership is by COORDINATE (collapse-safe). `aniso` is nnode*NTENSOR
// (2D: m11,m12,m22; 3D: m11,m12,m13,m22,m23,m33).
static void build_wall_anisotropic_metric(const Param &param, const SideInflowState &sis,
                                          int nnode, const double *coord, const double *scalar,
                                          double tol, double_vec &aniso)
{
    const int NT = (NDIMS == 2) ? 3 : 6;
    const double ratio = param.mesh.mmg_aniso_wall_ratio;    // > 1; coarsening factor ⟂ the wall
    const double xlen = param.mesh.xlength;
#ifdef THREED
    const double ylen = param.mesh.ylength;
#endif
    // Coarsen the wall-PERPENDICULAR axis only: R3 pins the wall nodes and R4 frees the interior
    // band, so inland is the only direction with movable nodes; wide wall elements also absorb the
    // incoming advection slowly. Coarsening ALONG the wall instead lengthens the transition and
    // degrades the profile's depth sampling. Limitation: a wide wall element's centroid sits
    // inland of the band, so restore_side_fields does not re-impose its ELEMENT fields (nodal
    // temperature still is; deep-wall element fields are ~0). 3D also coarsens the y walls.
    aniso.assign((std::size_t)nnode * NT, 0.0);
    for (int n = 0; n < nnode; ++n) {
        const double *p = coord + (std::size_t)n * NDIMS;
        const double h = scalar[n];
        double hx = h, hz = h;             // per-axis target; default isotropic (non-wall nodes)
#ifdef THREED
        double hy = h;
#endif
        // Coarsen the inflow-PERPENDICULAR axis at an inflow wall; tangential axes stay = h.
        if (sis.in_x0 && p[0] < tol)          hx = h * ratio;
        if (sis.in_x1 && p[0] > xlen - tol)   hx = h * ratio;
#ifdef THREED
        if (sis.in_y0 && p[1] < tol)          hy = h * ratio;
        if (sis.in_y1 && p[1] > ylen - tol)   hy = h * ratio;
#endif
        double *m = &aniso[(std::size_t)n * NT];
#ifdef THREED
        // (m11,m12,m13,m22,m23,m33) = diag(1/hx^2, 1/hy^2, 1/hz^2)
        m[0] = 1.0/(hx*hx); m[1] = 0.0; m[2] = 0.0;
        m[3] = 1.0/(hy*hy); m[4] = 0.0;
        m[5] = 1.0/(hz*hz);
#else
        // (m11,m12,m22) = diag(1/hx^2, 1/hz^2)  (2D coords are (x, z))
        m[0] = 1.0/(hx*hx); m[1] = 0.0; m[2] = 1.0/(hz*hz);
#endif
    }
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

    int_vec bdry_polygons[nbdrytypes];
    assemble_bdry_polygons(var, original_coord, original_connectivity,
                           bdry_polygons);

    // create a copy of original_coord and original_segment
    double_vec qcoord_vec;
    int_vec qconn_vec, qsegment_vec, qconn_from_1_vec, qsegment_from_1_vec, qsegflag_vec;

    original_coord.pack_to(qcoord_vec);
    original_connectivity.pack_to(qconn_vec);
    original_segment.pack_to(qsegment_vec);
    original_connectivity.pack_to(qconn_from_1_vec);
    original_segment.pack_to(qsegment_from_1_vec);
    original_segflag.pack_to(qsegflag_vec);

    double *qcoord = qcoord_vec.data();
    int *qconn = qconn_vec.data();
    int *qsegment = qsegment_vec.data();
    int *qconn_from_1 = qconn_from_1_vec.data();
    int *qsegment_from_1 = qsegment_from_1_vec.data();
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
    case 2:
        excl_func = &is_boundary;
        new_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
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
    compute_metric_field(param, var, *var.ntmp, *var.etmp, original_coord, original_connectivity);

    // By default MMG adapts the current mesh in place.
    int   mnode = old_nnode, melem = old_nelem, mseg = old_nseg;
    double *mcoord   = qcoord;
    int    *mconn1   = qconn_from_1;      // 0-indexed; mmg_adapt 1-indexes internally
    int    *mseg1    = qsegment_from_1;   // 0-indexed; mmg_adapt 1-indexes internally
    int    *msegflag = qsegflag;
    double *mmetric  = (*var.ntmp).data();

    // Material sunk below the fixed bottom (large max_boundary_distortion): flatten_bottom snaps
    // the boundary back up and the elements below invert, which MMG cannot adapt. Collapse the
    // outside nodes onto the boundary they crossed and let MMG adapt the result; no fallback to a
    // whole-domain Triangle/Tetgen remesh. Only fires in this pathological case.
    double_vec c_coord, c_metric;
    int_vec c_conn, c_seg, c_segflag;
    int_vec c_new_to_old_node, c_new_to_old_elem;   // collapse renumbering maps (new id -> old id)
    bool outside_material = has_outside_material(param, original_coord, old_bcflag,
                                                 old_nnode, min_dist);
    bool collapse_needed = outside_material;
    if (outside_material) {
        std::cerr << "  Material moved outside a restored boundary; collapsing it back onto "
                     "the boundary before MMG.\n";
        int cn = 0, ce = 0, cs = 0;
        collapse_outside_nodes(points_to_delete, old_nnode, old_nelem, old_nseg,
                               qcoord, old_bcflag, qconn, qsegment, qsegflag, *var.ntmp,
                               param.mesh.zlength, param.mesh.xlength, param.mesh.ylength, min_dist,
                               c_coord, c_conn, c_seg, c_segflag, c_metric,
                               c_new_to_old_node, c_new_to_old_elem, cn, ce, cs);
        mnode = cn; melem = ce; mseg = cs;
        mcoord = c_coord.data();
        mconn1 = c_conn.data();
        mseg1 = c_seg.data();
        msegflag = c_segflag.data();
        mmetric = c_metric.data();
    }
#ifndef THREED
    // Otherwise collapse crowded free-surface cliffs MMG cannot fix (short near-vertical boundary
    // segments that persist as tiny slivers), via the same renumber-aware collapse path.
    else {
        int cn = 0, ce = 0, cs = 0;
        int_vec removed;
        // Collapse at hmin, not at the metric floor: the more aggressive collapse perturbed the
        // hausd-constrained surface and produced more bad outputs.
        int n_surf = collapse_short_boundary_segments_2d(old_nnode, old_nelem, old_nseg,
                         qcoord, old_bcflag, qconn, qsegment, qsegflag, *var.ntmp,
                         refine_floors(param.mesh).hmin,
                         c_coord, c_conn, c_seg, c_segflag, c_metric,
                         c_new_to_old_node, c_new_to_old_elem, removed, cn, ce, cs);
        if (n_surf > 0) {
            std::cout << "    Collapsing " << n_surf
                      << " crowded free-surface node(s) below hmin before MMG.\n";
            mnode = cn; melem = ce; mseg = cs;
            mcoord = c_coord.data();
            mconn1 = c_conn.data();
            mseg1 = c_seg.data();
            msegflag = c_segflag.data();
            mmetric = c_metric.data();
            for (int n : removed) points_to_delete.push_back(n);   // record for remesh_affected + freeze region
            collapse_needed = true;
        }
    }
#endif

    // Adapt with the shared MMG driver: detect_flatten_broken supplies R6, mark_quiet_required
    // applies R1-R7 and derives the required arrays, mmg_adapt_quality_gated runs MMG under the
    // post-output gate. See the FREEZE / REQUIRED POLICY section.
    std::vector<char> flatten_broken =
        detect_flatten_broken(param, var, original_coord, qcoord, qconn,
                              old_nnode, old_nelem, points_to_delete);
    std::vector<char> req_node, req_elem;
    CollapseRemap remap = { &c_new_to_old_node, &c_new_to_old_elem, mconn1, &points_to_delete };
    mark_quiet_required(param, var, original_coord, original_connectivity, mnode, melem,
                        flatten_broken, collapse_needed ? &remap : nullptr,
                        req_node, req_elem);
    const char *rn = req_node.data();
    const char *re = req_elem.data();

    // Optional anisotropic wall metric (nullptr => scalar path). Built from the post-collapse
    // metric/coords handed to MMG; the gate never rescales the metric, so once is valid for all retries.
    double_vec wall_aniso;
    const double *mmetric_aniso = nullptr;
    if (param.mesh.remeshing_option == 13 && param.mesh.mmg_aniso_wall_ratio > 0.0) {
        SideInflowState sis = side_inflow_state(param, old_nnode, original_coord, old_bcflag);
        build_wall_anisotropic_metric(param, sis, mnode, mcoord, mmetric, min_dist, wall_aniso);
        mmetric_aniso = wall_aniso.data();
    }

    // Front-tracking (option 13): keep the incoming material boundary an MMG-required edge. On the
    // collapse path translate the old-mesh node ids through the remap (dropping edges whose
    // endpoint was collapsed away); on the direct path old ids == MMG ids.
    std::vector<int> iface_edges;
    const int *ie_ptr = nullptr; int n_ie = 0;
    std::vector<double> aug_coord, aug_metric;   // augmented mmg-input mesh (interface pre-subdivision)
    std::vector<int>    aug_conn;                // kept alive through the quality-gate retries below
    if (param.mesh.remeshing_option == 13) {
        SideInflowState sis_ie = side_inflow_state(param, old_nnode, original_coord, old_bcflag);
        uint inflow_bits = 0;
        if (sis_ie.in_x0) inflow_bits |= BOUNDX0;
        if (sis_ie.in_x1) inflow_bits |= BOUNDX1;
        iface_edges = build_interface_req_edges(var, original_connectivity, old_nnode, old_nelem, inflow_bits);
        if (collapse_needed && !iface_edges.empty()) {
            std::vector<int> old_to_new(old_nnode, -1);
            for (int n = 0; n < mnode; ++n) old_to_new[(*remap.new_to_old_node)[n]] = n;
            std::vector<int> tr; tr.reserve(iface_edges.size());
            for (std::size_t k = 0; k + 1 < iface_edges.size(); k += 2) {
                const int a = old_to_new[iface_edges[k]], b = old_to_new[iface_edges[k+1]];
                if (a >= 0 && b >= 0) { tr.push_back(a); tr.push_back(b); }
            }
            iface_edges.swap(tr);
        }
        // BOTH paths: pre-subdivide over-long interface edges (MMG cannot split a required edge).
        // Grows the local mesh buffers, the required masks and the aniso tensor, and repoints the
        // mesh handed to MMG. Gating this to the direct path let a wall-attached edge ratchet
        // longer for hundreds of kyr while every remesh took the collapse path.
        if (!iface_edges.empty()) {
            const int mnode0 = mnode;
            aug_coord.assign(mcoord, mcoord + (std::size_t)mnode*NDIMS);
            aug_conn.assign(mconn1, mconn1 + (std::size_t)melem*NODES_PER_ELEM);
            aug_metric.assign(mmetric, mmetric + mnode);
            int nn = mnode, ne = melem;
            presubdivide_interface_edges(aug_coord, aug_conn, aug_metric, req_node, req_elem,
                                         iface_edges, nn, ne, param.mesh.resolution);
            if (nn != mnode0) {
                if (mmetric_aniso) {   // extend the aniso tensor for the new interior nodes (isotropic)
                    for (int i = mnode0; i < nn; ++i) {
                        const double h = aug_metric[i]; const double inv = (h > 0.0) ? 1.0/(h*h) : 0.0;
                        wall_aniso.push_back(inv); wall_aniso.push_back(0.0); wall_aniso.push_back(inv);
                    }
                    mmetric_aniso = wall_aniso.data();
                }
                mnode = nn; melem = ne;
                mcoord = aug_coord.data(); mconn1 = aug_conn.data(); mmetric = aug_metric.data();
                rn = req_node.data(); re = req_elem.data();
            }
        }
        if (!iface_edges.empty()) { ie_ptr = iface_edges.data(); n_ie = (int)iface_edges.size() / 2; }
    }

    MMGInput  mmg_in = { mnode, melem, mseg, mcoord, mconn1, mseg1, msegflag, mmetric, rn, re, false, mmetric_aniso, ie_ptr, n_ie };
    MMGOutput mmg_out;
    mmg_adapt_quality_gated(param, mmg_in, mmg_out, req_node, req_elem);

    var.nnode = mmg_out.nnode;
    var.nelem = mmg_out.nelem;
    var.nseg  = mmg_out.nseg;

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

    int_vec bdry_polygons[nbdrytypes];
    assemble_bdry_polygons(var, original_coord, original_connectivity,
                           bdry_polygons);

    // create a copy of original_coord and original_segment
    double_vec qcoord_vec;
    int_vec qconn_vec, qsegment_vec, qconn_from_1_vec, qsegment_from_1_vec, qsegflag_vec;

    original_coord.pack_to(qcoord_vec);
    original_connectivity.pack_to(qconn_vec);
    original_segment.pack_to(qsegment_vec);
    original_connectivity.pack_to(qconn_from_1_vec);
    original_segment.pack_to(qsegment_from_1_vec);
    original_segflag.pack_to(qsegflag_vec);

    double *qcoord = qcoord_vec.data();
    int *qconn = qconn_vec.data();
    int *qsegment = qsegment_vec.data();
    int *qconn_from_1 = qconn_from_1_vec.data();
    int *qsegment_from_1 = qsegment_from_1_vec.data();
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
    case 2:
        excl_func = &is_boundary;
        new_bottom(old_bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
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
    compute_metric_field(param, var, *var.ntmp, *var.etmp, original_coord, original_connectivity);

    // By default MMG adapts the current mesh in place.
    int   mnode = old_nnode, melem = old_nelem, mseg = old_nseg;
    double *mcoord   = qcoord;
    int    *mconn1   = qconn_from_1;      // 0-indexed; mmg_adapt 1-indexes internally
    int    *mseg1    = qsegment_from_1;   // 0-indexed; mmg_adapt 1-indexes internally
    int    *msegflag = qsegflag;
    double *mmetric  = (*var.ntmp).data();

    // Material sunk below the fixed bottom (large max_boundary_distortion): flatten_bottom snaps
    // the boundary back up and the elements below invert, which MMG cannot adapt. Collapse the
    // outside nodes onto the boundary they crossed and let MMG adapt the result; no fallback to a
    // whole-domain Triangle/Tetgen remesh. Only fires in this pathological case.
    double_vec c_coord, c_metric;
    int_vec c_conn, c_seg, c_segflag;
    int_vec c_new_to_old_node, c_new_to_old_elem;   // collapse renumbering maps (new id -> old id)
    bool outside_material = has_outside_material(param, original_coord, old_bcflag,
                                                 old_nnode, min_dist);
    bool collapse_needed = outside_material;
    if (outside_material) {
        std::cerr << "  Material moved outside a restored boundary; collapsing it back onto "
                     "the boundary before MMG.\n";
        int cn = 0, ce = 0, cs = 0;
        collapse_outside_nodes(points_to_delete, old_nnode, old_nelem, old_nseg,
                               qcoord, old_bcflag, qconn, qsegment, qsegflag, *var.ntmp,
                               param.mesh.zlength, param.mesh.xlength, param.mesh.ylength, min_dist,
                               c_coord, c_conn, c_seg, c_segflag, c_metric,
                               c_new_to_old_node, c_new_to_old_elem, cn, ce, cs);
        mnode = cn; melem = ce; mseg = cs;
        mcoord = c_coord.data();
        mconn1 = c_conn.data();
        mseg1 = c_seg.data();
        msegflag = c_segflag.data();
        mmetric = c_metric.data();
    }
#ifndef THREED
    // Otherwise collapse crowded free-surface cliffs MMG cannot fix (short near-vertical boundary
    // segments that persist as tiny slivers), via the same renumber-aware collapse path.
    else {
        int cn = 0, ce = 0, cs = 0;
        int_vec removed;
        // Collapse at hmin, not at the metric floor: the more aggressive collapse perturbed the
        // hausd-constrained surface and produced more bad outputs.
        int n_surf = collapse_short_boundary_segments_2d(old_nnode, old_nelem, old_nseg,
                         qcoord, old_bcflag, qconn, qsegment, qsegflag, *var.ntmp,
                         refine_floors(param.mesh).hmin,
                         c_coord, c_conn, c_seg, c_segflag, c_metric,
                         c_new_to_old_node, c_new_to_old_elem, removed, cn, ce, cs);
        if (n_surf > 0) {
            std::cout << "    Collapsing " << n_surf
                      << " crowded free-surface node(s) below hmin before MMG.\n";
            mnode = cn; melem = ce; mseg = cs;
            mcoord = c_coord.data();
            mconn1 = c_conn.data();
            mseg1 = c_seg.data();
            msegflag = c_segflag.data();
            mmetric = c_metric.data();
            for (int n : removed) points_to_delete.push_back(n);   // record for remesh_affected + freeze region
            collapse_needed = true;
        }
    }
#endif

    // Adapt with the shared MMG driver: detect_flatten_broken supplies R6, mark_quiet_required
    // applies R1-R7 and derives the required arrays, mmg_adapt_quality_gated runs MMG under the
    // post-output gate. See the FREEZE / REQUIRED POLICY section.
    std::vector<char> flatten_broken =
        detect_flatten_broken(param, var, original_coord, qcoord, qconn,
                              old_nnode, old_nelem, points_to_delete);
    std::vector<char> req_node, req_elem;
    CollapseRemap remap = { &c_new_to_old_node, &c_new_to_old_elem, mconn1, &points_to_delete };
    mark_quiet_required(param, var, original_coord, original_connectivity, mnode, melem,
                        flatten_broken, collapse_needed ? &remap : nullptr,
                        req_node, req_elem);
    const char *rn = req_node.data();
    const char *re = req_elem.data();

    // DES_DEBUG_BSEG: diagnose bottom-BOUNDZ0 loss across MMG (env-gated, off by default).
    const bool dbg_bseg = std::getenv("DES_DEBUG_BSEG") != NULL;
    if (dbg_bseg) {
        const double zb = -param.mesh.zlength;
        std::map<std::pair<int,int>, std::vector<int>> edges;
        int nflag0 = 0;
        for (int s = 0; s < mseg; ++s) {
            int a = mseg1[s*NODES_PER_FACET], b = mseg1[s*NODES_PER_FACET+1];
            if (a > b) std::swap(a, b);
            edges[{a,b}].push_back(msegflag[s]);
            if (msegflag[s] == 0) ++nflag0;
        }
        int ndup = 0, ndup_mixed = 0;
        for (auto &kv : edges) {
            if (kv.second.size() < 2) continue;
            ++ndup;
            bool mixed = false;
            for (std::size_t i = 1; i < kv.second.size(); ++i)
                if (kv.second[i] != kv.second[0]) mixed = true;
            if (mixed) {
                ++ndup_mixed;
                int a = kv.first.first, b = kv.first.second;
                std::fprintf(stderr, "[bseg-in] DUP-MIXED edge %d-%d (%.0f,%.0f)-(%.0f,%.0f) flags:",
                             a, b, mcoord[a*NDIMS], mcoord[a*NDIMS+1], mcoord[b*NDIMS], mcoord[b*NDIMS+1]);
                for (int fl : kv.second) std::fprintf(stderr, " %d", fl);
                std::fprintf(stderr, "\n");
            }
        }
        std::fprintf(stderr, "[bseg-in] outside=%d nseg=%d flag0=%d dup_pairs=%d dup_mixed=%d\n",
                     (int)outside_material, mseg, nflag0, ndup, ndup_mixed);
        // bottom-line nodes of the MMG INPUT without a BOUNDZ0 segment endpoint
        std::vector<uint> nodeflag(mnode, 0);
        for (int s = 0; s < mseg; ++s) {
            nodeflag[mseg1[s*NODES_PER_FACET]]   |= (uint)msegflag[s];
            nodeflag[mseg1[s*NODES_PER_FACET+1]] |= (uint)msegflag[s];
        }
        for (int n = 0; n < mnode; ++n)
            if (std::fabs(mcoord[n*NDIMS+NDIMS-1] - zb) < 1.0 && !(nodeflag[n] & BOUNDZ0))
                std::fprintf(stderr, "[bseg-in] node %d ON LINE no-BOUNDZ0 x=%.0f segsum=%u req=%d\n",
                             n, mcoord[n*NDIMS], nodeflag[n], (int)rn[n]);
    }

    // Optional anisotropic wall metric (nullptr => scalar path). Built from the post-collapse
    // metric/coords handed to MMG; the gate never rescales the metric, so once is valid for all retries.
    double_vec wall_aniso;
    const double *mmetric_aniso = nullptr;
    if (param.mesh.remeshing_option == 13 && param.mesh.mmg_aniso_wall_ratio > 0.0) {
        SideInflowState sis = side_inflow_state(param, old_nnode, original_coord, old_bcflag);
        build_wall_anisotropic_metric(param, sis, mnode, mcoord, mmetric, min_dist, wall_aniso);
        mmetric_aniso = wall_aniso.data();
    }

    // Front-tracking (option 13): keep the incoming material boundary an MMG-required edge. On the
    // collapse path translate the old-mesh node ids through the remap (dropping edges whose
    // endpoint was collapsed away); on the direct path old ids == MMG ids.
    std::vector<int> iface_edges;
    const int *ie_ptr = nullptr; int n_ie = 0;
    std::vector<double> aug_coord, aug_metric;   // augmented mmg-input mesh (interface pre-subdivision)
    std::vector<int>    aug_conn;                // kept alive through the quality-gate retries below
    if (param.mesh.remeshing_option == 13) {
        SideInflowState sis_ie = side_inflow_state(param, old_nnode, original_coord, old_bcflag);
        uint inflow_bits = 0;
        if (sis_ie.in_x0) inflow_bits |= BOUNDX0;
        if (sis_ie.in_x1) inflow_bits |= BOUNDX1;
        iface_edges = build_interface_req_edges(var, original_connectivity, old_nnode, old_nelem, inflow_bits);
        if (collapse_needed && !iface_edges.empty()) {
            std::vector<int> old_to_new(old_nnode, -1);
            for (int n = 0; n < mnode; ++n) old_to_new[(*remap.new_to_old_node)[n]] = n;
            std::vector<int> tr; tr.reserve(iface_edges.size());
            for (std::size_t k = 0; k + 1 < iface_edges.size(); k += 2) {
                const int a = old_to_new[iface_edges[k]], b = old_to_new[iface_edges[k+1]];
                if (a >= 0 && b >= 0) { tr.push_back(a); tr.push_back(b); }
            }
            iface_edges.swap(tr);
        }
        // BOTH paths: pre-subdivide over-long interface edges (MMG cannot split a required edge).
        // Grows the local mesh buffers, the required masks and the aniso tensor, and repoints the
        // mesh handed to MMG. Gating this to the direct path let a wall-attached edge ratchet
        // longer for hundreds of kyr while every remesh took the collapse path.
        if (!iface_edges.empty()) {
            const int mnode0 = mnode;
            aug_coord.assign(mcoord, mcoord + (std::size_t)mnode*NDIMS);
            aug_conn.assign(mconn1, mconn1 + (std::size_t)melem*NODES_PER_ELEM);
            aug_metric.assign(mmetric, mmetric + mnode);
            int nn = mnode, ne = melem;
            presubdivide_interface_edges(aug_coord, aug_conn, aug_metric, req_node, req_elem,
                                         iface_edges, nn, ne, param.mesh.resolution);
            if (nn != mnode0) {
                if (mmetric_aniso) {   // extend the aniso tensor for the new interior nodes (isotropic)
                    for (int i = mnode0; i < nn; ++i) {
                        const double h = aug_metric[i]; const double inv = (h > 0.0) ? 1.0/(h*h) : 0.0;
                        wall_aniso.push_back(inv); wall_aniso.push_back(0.0); wall_aniso.push_back(inv);
                    }
                    mmetric_aniso = wall_aniso.data();
                }
                mnode = nn; melem = ne;
                mcoord = aug_coord.data(); mconn1 = aug_conn.data(); mmetric = aug_metric.data();
                rn = req_node.data(); re = req_elem.data();
            }
        }
        if (!iface_edges.empty()) { ie_ptr = iface_edges.data(); n_ie = (int)iface_edges.size() / 2; }
    }

    MMGInput  mmg_in = { mnode, melem, mseg, mcoord, mconn1, mseg1, msegflag, mmetric, rn, re, false, mmetric_aniso, ie_ptr, n_ie };
    MMGOutput mmg_out;
    mmg_adapt_quality_gated(param, mmg_in, mmg_out, req_node, req_elem);

    if (dbg_bseg) {
        const double zb = -param.mesh.zlength;
        int nflag0 = 0;
        for (int s = 0; s < mmg_out.nseg; ++s) if (mmg_out.segflag[s] == 0) ++nflag0;
        std::vector<uint> nodeflag(mmg_out.nnode, 0);
        std::vector<std::vector<int>> nodesegs(mmg_out.nnode);
        for (int s = 0; s < mmg_out.nseg; ++s) {
            for (int k = 0; k < NODES_PER_FACET; ++k) {
                int n = mmg_out.seg[s*NODES_PER_FACET+k];
                nodeflag[n] |= (uint)mmg_out.segflag[s];
                nodesegs[n].push_back(s);
            }
        }
        std::fprintf(stderr, "[bseg-out] nseg=%d flag0=%d\n", mmg_out.nseg, nflag0);
        for (int n = 0; n < mmg_out.nnode; ++n) {
            if (std::fabs(mmg_out.coord[n*NDIMS+NDIMS-1] - zb) < 1.0 && !(nodeflag[n] & BOUNDZ0)) {
                std::fprintf(stderr, "[bseg-out] node %d ON LINE no-BOUNDZ0 x=%.0f nsegs=%zu:",
                             n, mmg_out.coord[n*NDIMS], nodesegs[n].size());
                for (int s : nodesegs[n]) {
                    int a = mmg_out.seg[s*NODES_PER_FACET], b = mmg_out.seg[s*NODES_PER_FACET+1];
                    std::fprintf(stderr, " [s%d %d-%d f%d (%.0f,%.0f)-(%.0f,%.0f)]",
                                 s, a, b, mmg_out.segflag[s],
                                 mmg_out.coord[a*NDIMS], mmg_out.coord[a*NDIMS+1],
                                 mmg_out.coord[b*NDIMS], mmg_out.coord[b*NDIMS+1]);
                }
                std::fprintf(stderr, "\n");
            }
        }
    }

    var.nnode = mmg_out.nnode;
    var.nelem = mmg_out.nelem;
    var.nseg  = mmg_out.nseg;

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
}
#endif  // end of if THREED
#endif // end of if USEMMG

/* ADAPT (libadaptivity) optimize_mesh implementation removed. */
#if 0
// ADAPT implementation removed
#endif

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

    // compute_mass is the ONLY writer of var.volume_n, which the node loop below divides by;
    // running this before it makes every nodal size +inf. The size test also catches a fill left
    // over from a smaller nnode. compute_mass leaves its node kernel async, so wait before the
    // host read or this guard reads the pre-kernel zeros and fires on correct code.
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

    // Tiny-element trigger at smallest_vol / remesh_tiny_margin, BELOW the size the remesher floors
    // elements at: margin == 1 puts the floor on the trigger and remeshing thrashes; margin > 1
    // opens a hysteresis gap.
    const double smallest_vol = refine_floors(param.mesh).smallest_vol / param.mesh.remesh_tiny_margin;
    for (int e=0; e<var.nelem; e++) {
        if ((*var.volume)[e] < smallest_vol) {
            index = e;
            // report location + the nodes' boundary flags, so a chronic re-trigger spot
            // is identifiable (interior vs which boundary it is pinned to)
            double cent[NDIMS] = {0.};
            ConstConnAccessor conn = (*var.connectivity)[e];
            for (int k = 0; k < NODES_PER_ELEM; ++k)
                for (int d = 0; d < NDIMS; ++d)
                    cent[d] += (*var.coord)[conn[k]][d] / NODES_PER_ELEM;
            std::cout << "    The size of element #" << index << " is too small at ("
                      << cent[0] << ", " << cent[NDIMS-1] << "), node bcflags:";
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
        const double dist = param.mesh.max_boundary_distortion * param.mesh.resolution;
        for (int i=0; i<var.nnode; ++i) {
            if (is_bottom((*var.bcflag)[i])) {
                double z = (*var.coord)[i][NDIMS-1];
                if (std::fabs(z - bottom) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the bottm: z = " << z << "\n";
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
        const double dist = param.mesh.max_boundary_distortion * param.mesh.resolution;
        for (int i=0; i<var.nnode; ++i) {
            if (is_x0((*var.bcflag)[i])) {
                double x = (*var.coord)[i][0];
                if (std::fabs(x) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the x0 side: x = " << x << "\n";
                }
            } else if (is_x1((*var.bcflag)[i])) {
                double x = (*var.coord)[i][0];
                if (std::fabs(x - param.mesh.xlength) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the x1 side: x = " << x << "\n";
                }
#ifdef THREED
            } else if (is_y0((*var.bcflag)[i])) {
                double y = (*var.coord)[i][1];
                if (std::fabs(y) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the y0 side: y = " << y << "\n";
                }
            } else if (is_y1((*var.bcflag)[i])) {
                double y = (*var.coord)[i][1];
                if (std::fabs(y - param.mesh.ylength) > dist) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the y1 side: y = " << y << "\n";
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
                  << ") at (" << cent[0] << ", " << cent[NDIMS-1] << "), node bcflags:";
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
//   (1) FOLDED element: signed measure of the opposite sign to the mesh's dominant sign (not
//       hard-coded CCW; the 2D solver is orientation-agnostic). A size-relative tolerance ignores
//       near-degenerate tiny elements (the tiny-element trigger's job).
//   (2) SELF-INTERSECTING boundary (2D): two boundary segments sharing no node cross. Not tested
//       in 3D.
// Prints the offending geometry and dies with the mesh error code. `when` labels the call site.
static double tangle_orient2d(const double *a, const double *b, const double *c)
{
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}
void check_mesh_tangle(const char *when, const array_t &coord, const conn_t &connectivity,
                       int nelem, const segment_t &segment, const segflag_t &segflag, int nseg)
{
    // Signed measure of element e (2x area in 2D / 6x volume in 3D) and a size-relative
    // tolerance = REL * (characteristic element measure), so a genuine flip is caught but
    // floating-point noise on a near-degenerate element is not mistaken for one.
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


// Auto-detect the incoming-material profile at each side wall. For every x-side boundary,
// snapshot the initial wall column into var.side_profile[b] (see SideProfile in
// parameters.hpp): the wall nodes' locations and temperatures (sorted top->bottom), their
// node->element support, and the field values AND dominant marker material of the elements
// touching the wall. Depths are measured relative to the side's top point, pinning the
// profile to the surface as it evolves. Called once at run start (fresh or restart) when
// the sides are restored (remeshing_option 13); restore_side_fields() re-imposes the field
// profile on material entering through a restored side wall after every remesh, and the
// side-wall marker replenishment (markerset.cxx) assigns new markers the elem_mattype of
// their depth. The profile resolution is the wall-node spacing, which is why restored side
// walls are kept at ~mesh.resolution (compute_metric_field).
void detect_side_profile(const Param& param, Variables& var)
{
    #pragma acc wait   // init_elem_size_n, temperature, coord0 are async kernel outputs read here on the host
    const uint side_flag[2] = { BOUNDX0, BOUNDX1 };
    const int  side_idx [2] = { iboundx0, iboundx1 };
    const char *side_name[2] = { "x0", "x1" };
    for (int si = 0; si < 2; ++si) {
        const uint F = side_flag[si];
        SideProfile &prof = var.side_profile[side_idx[si]];
        prof.clear();

        // top point of this side = highest node on the wall
        double top_z = -std::numeric_limits<double>::max();
        for (int n = 0; n < var.nnode; ++n)
            if ((*var.bcflag)[n] & F)
                top_z = std::max(top_z, (*var.coord)[n][NDIMS-1]);
        if (top_z == -std::numeric_limits<double>::max()) continue;

        // wall nodes sorted top -> bottom (ties -- a constant-depth row in 3D -- broken by
        // the remaining coordinates for a deterministic order)
        std::vector<int> wnode;
        for (int n = 0; n < var.nnode; ++n)
            if ((*var.bcflag)[n] & F)
                wnode.push_back(n);
        if (wnode.size() < 2) continue;
        std::sort(wnode.begin(), wnode.end(), [&](int a, int b) {
            double da = top_z - (*var.coord)[a][NDIMS-1];
            double db = top_z - (*var.coord)[b][NDIMS-1];
            if (da != db) return da < db;
            for (int d = 0; d < NDIMS-1; ++d)
                if ((*var.coord)[a][d] != (*var.coord)[b][d])
                    return (*var.coord)[a][d] < (*var.coord)[b][d];
            return a < b;
        });

        // wall elements = union of the wall nodes' supports (every element with at least
        // one node on the wall), ascending; e2l maps global -> local profile index
        std::vector<int> welem;
        for (int n : wnode) {
            const int npatch = var.support.size(n);
            const int* patch = var.support.patch(n);
            for (int k = 0; k < npatch; ++k)
                welem.push_back(patch[k]);
        }
        std::sort(welem.begin(), welem.end());
        welem.erase(std::unique(welem.begin(), welem.end()), welem.end());
        std::map<int,int> e2l;
        for (std::size_t le = 0; le < welem.size(); ++le)
            e2l[welem[le]] = static_cast<int>(le);

        // node arrays + CSR node -> element support
        prof.node_support_idx.push_back(0);
        for (int n : wnode) {
            prof.node_reldepth.push_back(top_z - (*var.coord)[n][NDIMS-1]);
            for (int d = 0; d < NDIMS; ++d)
                prof.node_coord.push_back((*var.coord)[n][d]);
            for (int d = 0; d < NDIMS; ++d)
                prof.node_coord0.push_back((*var.coord0)[n][d]);
            prof.node_temperature.push_back((*var.temperature)[n]);
            const int npatch = var.support.size(n);
            const int* patch = var.support.patch(n);
            for (int k = 0; k < npatch; ++k)
                prof.node_support_arr.push_back(e2l[patch[k]]);
            prof.node_support_idx.push_back(static_cast<int>(prof.node_support_arr.size()));
        }
        // per wall element: dominant marker material (the side's layering, read by the side-wall
        // replenishment) and radiogenic source. Strain/stress are not recorded: incoming elements
        // enter pristine.
        for (int e : welem) {
            ConstConnAccessor conn = (*var.connectivity)[e];
            double cz = 0;
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                cz += (*var.coord)[conn[i]][NDIMS-1];
            cz /= NODES_PER_ELEM;
            prof.elem_reldepth.push_back(top_z - cz);
            const int_vec &em = (*var.elemmarkers)[e];
            int mt = 0, best = -1;
            for (int m = 0; m < (int)em.size(); ++m)
                if (em[m] > best) { best = em[m]; mt = m; }
            prof.elem_mattype.push_back(mt);
            prof.elem_radiogenic_source.push_back((*var.radiogenic_source)[e]);
        }

        // layering summary: dominant materials top->bottom, consecutive repeats collapsed
        std::vector<std::pair<double,int>> col;
        for (int le = 0; le < prof.nelem(); ++le)
            col.emplace_back(prof.elem_reldepth[le], prof.elem_mattype[le]);
        std::sort(col.begin(), col.end());
        std::cout << "  Incoming profile at side " << side_name[si]
                  << " (top z=" << top_z << "): " << prof.nnode() << " nodes, "
                  << prof.nelem() << " elements, T "
                  << prof.node_temperature.front() << " -> "
                  << prof.node_temperature.back() << ", mattype";
        for (std::size_t i = 0; i < col.size(); ++i)
            if (i == 0 || col[i].second != col[i-1].second)
                std::cout << " " << col[i].second;
        std::cout << "\n";
    }
}


namespace {

// Average of the profile values over the tie-run of nodes at (numerically) the same depth
// as node k (3D: a whole constant-depth row of wall nodes; 2D: exact duplicates only).
double side_profile_run_avg(const double_vec &d, const double_vec &v,
                            std::size_t k, double tol)
{
    std::size_t lo = k, hi = k;
    while (lo > 0 && d[k] - d[lo-1] < tol) --lo;
    while (hi + 1 < d.size() && d[hi+1] - d[k] < tol) ++hi;
    double s = 0;
    for (std::size_t m = lo; m <= hi; ++m) s += v[m];
    return s / (hi - lo + 1);
}

} // anonymous namespace

// Nodal profile value at relative depth `rel` in a depth-sorted (depth, value) table: linear
// interpolation, clamped at the ends, same-depth rows averaged. Shared with markerset.cxx.
double side_profile_node_value(const double_vec &d, const double_vec &v,
                               double rel, double tol)
{
    if (rel <= d.front()) return side_profile_run_avg(d, v, 0, tol);
    if (rel >= d.back())  return side_profile_run_avg(d, v, d.size()-1, tol);
    std::size_t i = std::upper_bound(d.begin(), d.end(), rel) - d.begin();  // d[i-1] <= rel < d[i]
    if (d[i] - d[i-1] < tol)   // rel sits inside one same-depth row
        return side_profile_run_avg(d, v, i-1, tol);
    double vl = side_profile_run_avg(d, v, i-1, tol);
    double vh = side_profile_run_avg(d, v, i, tol);
    double t = (rel - d[i-1]) / (d[i] - d[i-1]);
    return vl + t * (vh - vl);
}

// Profile elements spanning relative depth `rel`, via the stored node->element support (2D: the
// 1-2 wall elements of that interval; 3D: the element row); falls back to the nearest centroid.
// Shared with markerset.cxx.
void side_profile_elems_at(const SideProfile &p, double rel, double tol, int_vec &out)
{
    out.clear();
    if (p.nelem() == 0) return;
    const double_vec &d = p.node_reldepth;
    double r = std::min(std::max(rel, d.front()), d.back());
    std::size_t i = std::upper_bound(d.begin(), d.end(), r) - d.begin();
    if (i == 0) ++i;                        // r == d.front()
    if (i == d.size()) --i;                 // r == d.back()
    // union of the support over each bracketing tie-group
    auto group_elems = [&](std::size_t k, int_vec &g) {
        std::size_t lo = k, hi = k;
        while (lo > 0 && d[k] - d[lo-1] < tol) --lo;
        while (hi + 1 < d.size() && d[hi+1] - d[k] < tol) ++hi;
        for (std::size_t m = lo; m <= hi; ++m)
            for (int j = p.node_support_idx[m]; j < p.node_support_idx[m+1]; ++j)
                g.push_back(p.node_support_arr[j]);
        std::sort(g.begin(), g.end());
        g.erase(std::unique(g.begin(), g.end()), g.end());
    };
    int_vec below, above;
    group_elems(i-1, below);
    group_elems(i, above);
    std::set_intersection(below.begin(), below.end(), above.begin(), above.end(),
                          std::back_inserter(out));
    if (!out.empty()) return;
    std::size_t best = 0;
    for (std::size_t e = 1; e < p.elem_reldepth.size(); ++e)
        if (std::fabs(p.elem_reldepth[e] - rel) < std::fabs(p.elem_reldepth[best] - rel))
            best = e;
    out.push_back(static_cast<int>(best));
}

namespace {

// Give the material that entered through a RESTORED side wall its original field structure,
// just like the marker replenishment keeps the side's original material layering
// (replenish_markers_by_side_profile). Between remeshes an inflow BC carries the wall inward; the
// remesh restores the wall to its initial plane, and the new nodes in the band between the
// restored plane and the old wall -- the incoming material -- fall OUTSIDE the old mesh, where
// barycentric_node_interpolation used a nearest-node fallback that smears the field
// structure. Those nodes (recorded in var.remesh_node_outside) inside the band of a side with
// a detected profile get the profile temperature for their depth below the side's CURRENT top
// point, and the profile reference coordinates (coord0): the depth (z) component interpolated
// from the recorded wall coord0, the wall-normal component set to plane + wall_shift, where
// wall_shift accumulates the displacement of every wall restore (the remesher snaps the
// drifted wall back to x=0 / x=xlength each remesh; the increments sum in
// SideProfile::wall_shift, updated here). That makes coord0[x] a material (Lagrangian) entry
// coordinate OUTSIDE the domain -- each inflow generation gets a distinct reference and
// coord - coord0 stays a consistent displacement field (3D: the y component is left to the
// barycentric interpolation -- an x-side profile carries no y structure);
// the band's elements (centroid in the band, at least one node outside the old mesh) enter
// PRISTINE: plstrain, delta_plstrain and strain are reset to zero, and the radiogenic source
// is restored from the recorded wall elements spanning their centroid depth (node-support
// lookup). Stress/stressyy and the RSF fields are NOT touched -- they keep the regular
// NN/SPR-remapped values. Gated per side on ACTUAL inflow (sis, measured on the
// old mesh) on top of the natural self-gating (a static or outflow wall produces no outside
// nodes in its band), so genuine evolution at e.g. a fixed backstop is never overwritten.
// x-sides only (like detect_side_profile). MUST run before the plstrain_remesh snapshot,
// so the zeroed plstrain seeds the conservative-freeze baseline.
void restore_side_fields(const Param &param, Variables &var,
                         const SideInflowState &sis)
{
    if (param.mesh.remeshing_option != 13 || var.remesh_node_outside.empty()) return;

    const uint side_flag[2] = { BOUNDX0, BOUNDX1 };
    const int  side_idx [2] = { iboundx0, iboundx1 };
    const bool   side_in[2] = { sis.in_x0, sis.in_x1 };
    const double old_x0_in = sis.x0_in, old_x1_in = sis.x1_in;
    // the planes the remesher snaps the restored walls back to (cf. has_outside_material /
    // flatten_x0/x1): the wall-normal reference coordinate of incoming material
    const double wall_plane[2] = { 0.0, param.mesh.xlength };
    const double tol = refine_floors(param.mesh).min_dist;
    const double depth_tol = 1e-6 * param.mesh.resolution;  // same-depth node merge tolerance
    int count_n = 0, count_e = 0;
    for (int t = 0; t < 2; ++t) {
        if (!side_in[t]) continue;   // no incoming material at this side
        SideProfile &prof = var.side_profile[side_idx[t]];
        if (prof.empty()) continue;

        // Record this remesh's wall restore: the signed snap-back displacement accumulates in
        // wall_shift so consecutive inflow generations get distinct entry references.
        prof.wall_shift += (t == 0) ? (wall_plane[0] - old_x0_in)
                                    : (wall_plane[1] - old_x1_in);
        const double entry_x = wall_plane[t] + prof.wall_shift;

        // current top point of this side; the profile depths are measured below it
        double top_z = -std::numeric_limits<double>::max();
        for (int n = 0; n < var.nnode; ++n)
            if ((*var.bcflag)[n] & side_flag[t])
                top_z = std::max(top_z, (*var.coord)[n][NDIMS-1]);
        if (top_z == -std::numeric_limits<double>::max()) continue;

        // de-interleave the recorded wall coord0 depth component for the table lookup
        double_vec c0z(prof.nnode());
        for (int i = 0; i < prof.nnode(); ++i)
            c0z[i] = prof.node_coord0[i*NDIMS + NDIMS-1];

        for (int n = 0; n < var.nnode; ++n) {
            if (!var.remesh_node_outside[n]) continue;
            double x = (*var.coord)[n][0];
            bool in_band = (t == 0) ? (x <= old_x0_in + tol) : (x >= old_x1_in - tol);
            if (!in_band) continue;
            double rel = top_z - (*var.coord)[n][NDIMS-1];
            (*var.temperature)[n] =
                side_profile_node_value(prof.node_reldepth, prof.node_temperature, rel, depth_tol);
            (*var.coord0)[n][0] = entry_x;
            (*var.coord0)[n][NDIMS-1] =
                side_profile_node_value(prof.node_reldepth, c0z, rel, depth_tol);
            ++count_n;
        }

        if (prof.nelem() == 0) continue;
        int_vec cand;
        for (int e = 0; e < var.nelem; ++e) {
            ConstConnAccessor conn = (*var.connectivity)[e];
            bool any_outside = false;
            double cx = 0, cz = 0;
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                int n = conn[i];
                if (var.remesh_node_outside[n]) any_outside = true;
                cx += (*var.coord)[n][0];
                cz += (*var.coord)[n][NDIMS-1];
            }
            if (!any_outside) continue;
            cx /= NODES_PER_ELEM;
            cz /= NODES_PER_ELEM;
            bool in_band = (t == 0) ? (cx <= old_x0_in + tol) : (cx >= old_x1_in - tol);
            if (!in_band) continue;
            // pristine incoming material: no accumulated strain. Stress/stressyy and the
            // RSF fields are NOT touched -- they keep the regular NN/SPR-remapped values.
            (*var.plstrain)[e] = 0.;
            (*var.delta_plstrain)[e] = 0.;
            for (int d = 0; d < NSTR; ++d)
                (*var.strain)[e][d] = 0.;
            // radiogenic source of the recorded wall elements spanning this depth
            side_profile_elems_at(prof, top_z - cz, depth_tol, cand);
            if (!cand.empty()) {
                double rad = 0;
                for (int le : cand)
                    rad += prof.elem_radiogenic_source[le];
                (*var.radiogenic_source)[e] = rad / cand.size();
            }
            ++count_e;
        }
    }
    if (count_n || count_e)
        std::cout << "    Restored the side profile on " << count_n
                  << " incoming node(s) and " << count_e << " incoming element(s)"
                  << " (cumulative wall shift x0/x1: "
                  << var.side_profile[iboundx0].wall_shift << "/"
                  << var.side_profile[iboundx1].wall_shift << " m).\n";
}

} // anonymous namespace


void remesh(const Param &param, Variables &var, int bad_quality)
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    int64_t time_tmp = get_nanoseconds();

    std::cout << "  Remeshing starts...\n";
#ifdef ACC
    {
        size_t free_bytes, total_bytes;
        knn_bvh_mem_info(&free_bytes, &total_bytes);
        std::cout << "  [GPU mem] remesh start: free="
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

    // Per-side inflow state of the restored side walls on the OLD mesh (remeshing_option 13),
    // filled inside the block below while var.bcflag still describes the old mesh; consumed by
    // restore_side_fields() after spr_node_to_elem, when the new-mesh element stress is final.
    SideInflowState side_inflow;

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

        // Reset the per-remesh record of boundary-reshaped old nodes; the meshing routine fills it
        // for barycentric_node_interpolation.
        var.remesh_affected_old_node.assign(old_nnode, 0);

        // Inflow state taken while var.bcflag still describes the old mesh: a new node outside the
        // old mesh AND inside a receiving side's band is incoming material (restore_side_fields).
        if (param.mesh.remeshing_option == 13)
            side_inflow = side_inflow_state(param, old_nnode, old_coord, *var.bcflag);

#ifdef THREED
        if (param.mesh.meshing_elem_shape == 0) {
#if defined USEMMG
            remesh_engine = "MMG";
            optimize_mesh(param, var, bad_quality, old_coord, old_connectivity,
                    old_segment, old_segflag);
#else
            remesh_engine = "Tetgen";
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

        // Per-engine remesh report. Reports the free-surface (BOUNDZ1) facet count rather than raw
        // nseg: the initial mesh holds only input segments while remeshers emit the full edge set.
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
            std::cout << "  [GPU mem] before reallocate_tmp: free="
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

        // Done with the old mesh: drop the boundary-reshaped-node record. remesh_node_outside stays
        // alive for restore_side_fields() below.
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
        std::cout << "  [GPU mem] before reallocate_variables: free="
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

    // Incoming side material gets the recorded side profile (nodal temperature and coord0, zeroed
    // strain and the profile radiogenic source on the band's elements). Before the plstrain_remesh
    // snapshot, so the zeroed plstrain seeds the freeze baseline.
    restore_side_fields(param, var, side_inflow);
    var.remesh_node_outside.clear();

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

    // Snapshot plastic strain on the new mesh: the baseline the NEXT remesh's R2 test compares
    // against, so a fossil band that stops growing is frozen. After the interpolation set the
    // new-mesh var.plstrain; sized to the new var.nelem.
    delete var.plstrain_remesh;
    var.plstrain_remesh = new double_vec(*var.plstrain);

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

    // Option 13: re-seat the HORIZONTAL coord0 of pure-bottom nodes so the incoming reference is
    // continuous down to the floor. The flat bottom edge is re-discretized every remesh, and the
    // barycentric interpolation of coord0 there re-samples the OLD bottom edge, a 1-D chain that
    // never receives the wall shift: a spurious coord0[x] jump at the inflow bottom corner. Copy
    // the displacement (coord - coord0) of the nearest connected node ABOVE the floor. Side-wall
    // bottom corners are skipped (the wall column owns their coord0).
    if (param.mesh.remeshing_option==13) {
        int nbot13 = static_cast<int>(var.bnodes[iboundz0]->size());   // not const: named in shared(), gcc 8 rejects a const there
#ifndef ACC
        #pragma omp parallel for default(none) shared(var, nbot13)
#endif
        for (int i=0; i<nbot13; ++i) {
            const int n = (*var.bnodes[iboundz0])[i];
            const uint f = (*var.bcflag)[n];
            if (is_x0(f) || is_x1(f)
#ifdef THREED
                || is_y0(f) || is_y1(f)
#endif
               ) continue;                       // side-wall corner: owned by the wall logic
            // nearest connected node that is NOT itself on the floor (i.e. lies above it)
            int m = -1;
            double best = std::numeric_limits<double>::max();
            const int npatch = var.support.size(n);
            const int* patch = var.support.patch(n);
            for (int k = 0; k < npatch; ++k) {
                const int e = patch[k];
                ConstConnAccessor conn = (*var.connectivity)[e];
                for (int j=0; j<NODES_PER_ELEM; ++j) {
                    const int p = conn[j];
                    if (p == n || ((*var.bcflag)[p] & BOUNDZ0)) continue;
                    double d2 = 0;
                    for (int d=0; d<NDIMS-1; ++d) {
                        const double dd = (*var.coord)[p][d] - (*var.coord)[n][d];
                        d2 += dd*dd;
                    }
                    if (d2 < best) { best = d2; m = p; }
                }
            }
            if (m < 0) continue;
            for (int d=0; d<NDIMS-1; ++d)
                (*var.coord0)[n][d] =
                    (*var.coord)[n][d] - ((*var.coord)[m][d] - (*var.coord0)[m][d]);
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
        std::cout << "  [GPU mem] remesh end:   free="
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
