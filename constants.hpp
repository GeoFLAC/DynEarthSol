#ifndef DYNEARTHSOL3D_CONSTANTS_HPP
#define DYNEARTHSOL3D_CONSTANTS_HPP

# include <cmath>
# include <map>

#if !defined(M_PI)
#define M_PI 3.1415926535897932384626433832795
#endif

// # of spatial dimension, either 2 or 3
#ifdef THREED
const int NDIMS = 3;
#else
const int NDIMS = 2;
#endif

// triangles (3 node) in 2D and tetrahedra (4 nodes) in 3D
const int NODES_PER_ELEM = NDIMS + 1;
const int NODES_PER_REGULAR = (NDIMS-1) * 4;

// # of indep. components of a symmetric tensor
// 2D -> 3;  3D -> 6
const int NSTR = NDIMS * (NDIMS + 1) / 2;

// Flags for boundary
const int iboundx0 = 0;
const int iboundx1 = 1;
const int iboundy0 = 2;
const int iboundy1 = 3;
const int iboundz0 = 4;
const int iboundz1 = 5;
const int iboundn0 = 6;
const int iboundn1 = 7;
const int iboundn2 = 8;
const int iboundn3 = 9;
const int nbdrytypes = 10;
const int nbdrytypes_hydro = 6;

typedef unsigned int uint;
const uint BOUNDX0 = 1 << iboundx0;  //    1, western (left in 2D)
const uint BOUNDX1 = 1 << iboundx1;  //    2, eastern (right in 2D)
const uint BOUNDY0 = 1 << iboundy0;  //    4, southern
const uint BOUNDY1 = 1 << iboundy1;  //    8, northern
const uint BOUNDZ0 = 1 << iboundz0;  //   16, bottom
const uint BOUNDZ1 = 1 << iboundz1;  //   32, top
const uint BOUNDN0 = 1 << iboundn0;  //   64, arbitrary (not parallel to x,y,z)
const uint BOUNDN1 = 1 << iboundn1;  //  128, arbitrary (not parallel to x,y,z)
const uint BOUNDN2 = 1 << iboundn2;  //  256, arbitrary (not parallel to x,y,z)
const uint BOUNDN3 = 1 << iboundn3;  //  512, arbitrary (not parallel to x,y,z)

const uint BOUND_ANY = BOUNDX0 | BOUNDX1 | BOUNDY0 | BOUNDY1 | BOUNDZ0 | BOUNDZ1
                          | BOUNDN0 | BOUNDN1 | BOUNDN2 | BOUNDN3;

// # of facets (edges) per element, 3 for 2D, 4 for 3D
const int FACETS_PER_ELEM = NDIMS + 1;
// # of nodes per facets (edges), 2 for 2D, 3 for 3D
const int NODES_PER_FACET = NDIMS;

// local node # of a facet
#ifdef THREED
// the nodes are ordered counter-clockwise, if viewed from outside
const int NODE_OF_FACET[FACETS_PER_ELEM][NODES_PER_FACET] =
    {{1,2,3},
     {0,3,2},
     {0,1,3},
     {0,2,1}};
#else
// the nodes are ordered counter-clockwise, if viewed from above
const int NODE_OF_FACET[FACETS_PER_ELEM][NODES_PER_FACET] =
    {{1,2},
     {2,0},
     {0,1}};
#endif

const double YEAR2SEC = 365.2422 * 86400;
const double DEG2RAD = M_PI / 180;

#endif
