#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <sstream>
#include <set>

#ifdef THREED

#define TETLIBRARY
#include "tetgen/tetgen.h"
#undef TETLIBRARY

#ifdef USEMMG
#include "mmg/mmg3d/libmmg3d.h"
#endif

#else // !THREED

#ifdef USEMMG
#include "mmg/mmg2d/libmmg2d.h"
#endif

#endif // THREED

#ifdef USEEXODUS
#include "netcdf.h"
#include "exodusII.h"
#endif // USEEXODUS

#define REAL double
#define VOID void
#define ANSI_DECLARATORS
#include "triangle/triangle.h"
#undef REAL
#undef VOID
#undef ANSI_DECLARATORS


#include "constants.hpp"
#include "parameters.hpp"
#include "sortindex.hpp"
#include "utils.hpp"
#include "mesh.hpp"
#include "mmg_utils.hpp"
#include "markerset.hpp"
#include "remeshing.hpp"
#include "ic.hpp"

#ifdef WIN32
#ifndef _MSC_VER
namespace std {
  static std::string to_string(long double t)
  {
    char temp[32];
    sprintf(temp,"%f",double(t));
    return std::string(temp);
  }
}
#endif //_MSC_VER
#endif // WIN32

// equilateral triangle area = 0.433*s^2, equilateral tetrahedron volume = 0.118*s^3
#ifdef THREED
const double sizefactor = 0.118;
#else
const double sizefactor = 0.433;
#endif


namespace { // anonymous namespace

void set_verbosity_str(std::string &verbosity, int meshing_verbosity)
{
    switch (meshing_verbosity) {
    case -1:
        verbosity = "Q";
        break;
    case 0:
        verbosity = "";
        break;
    case 1:
        verbosity = "V";
        break;
    case 2:
        verbosity = "VV";
        break;
    case 3:
        verbosity = "VVV";
        break;
    default:
        verbosity = "";
        break;
    }
}


void set_volume_str(std::string &vol, double max_volume)
{
    vol.clear();
    if (max_volume > 0) {
        vol += 'a';
        vol += std::to_string((long double)max_volume);
    }
    else if (max_volume == 0) {
        // max_volume is set inside regional attributes (usu. read from poly files)
        vol += 'a';
    }
}


void set_2d_quality_str(std::string &quality, double min_angle)
{
    quality.clear();
    if (min_angle > 0) {
        quality += 'q';
        quality += std::to_string((long double)min_angle);
    }
}


void set_steiner_str(std::string &steiner, double max_steiner_factor, int npoints,
                     double max_area)
{
    // Triangle's -S switch caps the number of Steiner points (vertices not in the
    // input) it may insert. Without it (steinerleft == -1), Triangle's encroached-
    // subsegment splitting can fail to terminate when two input segments meet at a
    // very small angle -- each split encroaches the other, forever (see the loop at
    // triangle.c splitencsegs()). A cap set well above any healthy remesh bounds this
    // runaway without changing the mesh in the normal case, because steinerleft never
    // reaches 0.
    //
    // The cap is factor*npoints, so it is only meaningful when the output size stays
    // close to the input size -- i.e. during remeshing (max_area < 0, no area target).
    // When an area constraint is active (initial meshing), the legitimate point count
    // is driven by max_area rather than npoints, so the cap must NOT be applied or it
    // would truncate a perfectly good mesh.
    steiner.clear();
    if (max_steiner_factor > 0 && npoints > 0 && max_area < 0) {
        long cap = (long)(max_steiner_factor * npoints);
        if (cap < 1) cap = 1;
        if (cap > 100000000L) cap = 100000000L; // keep within int, far above any real need
        steiner += 'S';
        steiner += std::to_string(cap);
    }
}


void create_quadrilateral_cells(Variables &var, int *&cells) {
#ifndef THREED
    cells = new int[(var.nx-1)*(var.nz-1)*4];
#else
    cells = new int[(var.nx-1)*(var.ny-1)*(var.nz-1)*8];
#endif
    int cell_idx = 0;

#ifndef THREED
    for (int i = 0; i < var.nx - 1; ++i) {
        for (int j = 0; j < var.nz - 1; ++j) {
            int idx0 = i * var.nz + j;
            int idx1 = idx0 + var.nz;
            cells[cell_idx] = idx0;
            cells[cell_idx+1] = idx1;
            cells[cell_idx+2] = idx1 + 1;
            cells[cell_idx+3] = idx0 + 1;
            cell_idx += 4;
        }
    }
#else
    for (int i = 0; i < var.nx - 1; ++i) {
        for (int j = 0; j < var.ny - 1; ++j) {
            for (int k = 0; k < var.nz - 1; ++k) {
                // 8 vertices index
                int idx0 = i * var.ny * var.nz + j * var.nz + k;
                int idx1 = idx0 + var.nz;
                int idx2 = idx1 + var.ny * var.nz;
                int idx3 = idx2 - var.nz;
                cells[cell_idx] = idx0;
                cells[cell_idx+1] = idx1;
                cells[cell_idx+2] = idx2;
                cells[cell_idx+3] = idx3;
                cells[cell_idx+4] = idx0 + 1;
                cells[cell_idx+5] = idx1 + 1;
                cells[cell_idx+6] = idx2 + 1;
                cells[cell_idx+7] = idx3 + 1;
                cell_idx += 8;
            }
        }
    }
#endif
}

void divide_hexahedron_to_tetrahedra_index(ConstRegularAccessor cell, int order, int n, int *conn) {

    if (order == 0) {
        switch (n)
        {
        case 0:
            conn[0] = cell[0];
            conn[1] = cell[1];
            conn[2] = cell[2];
            conn[3] = cell[5];
            break;
        case 1:
            conn[0] = cell[0];
            conn[1] = cell[2];
            conn[2] = cell[3];
            conn[3] = cell[7];
            break;
        case 2:
            conn[0] = cell[0];
            conn[1] = cell[4];
            conn[2] = cell[5];
            conn[3] = cell[7];
            break;
        case 3:
            conn[0] = cell[2];
            conn[1] = cell[5];
            conn[2] = cell[6];
            conn[3] = cell[7];
            break;
        case 4:
            conn[0] = cell[0];
            conn[1] = cell[5];
            conn[2] = cell[2];
            conn[3] = cell[7];
            break;
        default:
            die(EXIT_INTERNAL_UNREACHABLE);
            break;
        }
    } else {
        switch (n)
        {
        case 0:
            conn[0] = cell[1];
            conn[1] = cell[2];
            conn[2] = cell[3];
            conn[3] = cell[6];
            break;
        case 1:
            conn[0] = cell[0];
            conn[1] = cell[1];
            conn[2] = cell[3];
            conn[3] = cell[4];
            break;
        case 2:
            conn[0] = cell[1];
            conn[1] = cell[4];
            conn[2] = cell[5];
            conn[3] = cell[6];
            break;
        case 3:
            conn[0] = cell[3];
            conn[1] = cell[4];
            conn[2] = cell[6];
            conn[3] = cell[7];
            break;
        case 4:
            conn[0] = cell[1];
            conn[1] = cell[3];
            conn[2] = cell[4];
            conn[3] = cell[6];
            break;
        default:
            die(EXIT_INTERNAL_UNREACHABLE);
            break;
        }
    }
}


void create_elem_from_cell(const Variables& var, int *&connectivity) {
    connectivity = new int[var.nelem*NODES_PER_ELEM];
#ifndef THREED

    for (int i = 0; i < var.nx - 1; ++i) {
        for (int j = 0; j < var.nz - 1; ++j) {
            int idx = i * (var.nz - 1) + j;
            // should be counter-clockwise
            if ((i+j) % 2 == 0){
                connectivity[idx*6] = (*var.cell)[idx][0];
                connectivity[idx*6+1] = (*var.cell)[idx][2];
                connectivity[idx*6+2] = (*var.cell)[idx][1];

                connectivity[idx*6+3] = (*var.cell)[idx][0];
                connectivity[idx*6+4] = (*var.cell)[idx][3];
                connectivity[idx*6+5] = (*var.cell)[idx][2];                
            } else {
                connectivity[idx*6] = (*var.cell)[idx][0];
                connectivity[idx*6+1] = (*var.cell)[idx][3];
                connectivity[idx*6+2] = (*var.cell)[idx][1];

                connectivity[idx*6+3] = (*var.cell)[idx][1];
                connectivity[idx*6+4] = (*var.cell)[idx][3];
                connectivity[idx*6+5] = (*var.cell)[idx][2];
            }
        }
    }
#else
    int conn[NODES_PER_ELEM];
    for (int i = 0; i<var.nx - 1; ++i) {
        for (int j = 0; j<var.ny -1; ++j) {
            for (int k = 0; k<var.nz - 1; ++k) {
                int idx = i * (var.ny - 1) * (var.nz - 1) + j * (var.nz - 1) + k;
                int order = (i+j+k)%2;
                for (int n = 0; n < 5; ++n) {
                    int idx_elem = idx*5+n;
                    divide_hexahedron_to_tetrahedra_index((*var.cell)[idx], order, n, conn);
                    connectivity[idx_elem*4] = conn[0];
                    connectivity[idx_elem*4+1] = conn[1];
                    connectivity[idx_elem*4+2] = conn[2];
                    connectivity[idx_elem*4+3] = conn[3];
                }
            }
        }
    }
#endif
}


void create_rect_node(const Param& param, const Variables& var, double *&points) {
    points = new double[var.nnode*NDIMS];
    double dx = param.mesh.xlength / (var.nx - 1);
    double dz = -param.mesh.zlength / (var.nz - 1);
#ifdef THREED
    double dy = param.mesh.ylength / (var.ny - 1);
#endif

#ifndef THREED
    for (int i = 0; i < var.nx; ++i) {
        for (int j = 0; j < var.nz; ++j) {
            points[(j + i * var.nz)*2] = i * dx;
            points[(j + i * var.nz)*2+1] = j * dz;
        }
    }
#else
    for (int i = 0; i < var.nx; ++i) {
        for (int j = 0; j < var.ny; ++j) {
            for (int k = 0; k < var.nz; ++k) {
                int idx = k + j * var.nz + i * var.ny * var.nz;
                points[idx*3] = i * dx;
                points[idx*3+1] = j * dy;
                points[idx*3+2] = k * dz;
            }
        }
    }
#endif
}


void create_regular_segments(const Variables& var, int *&segments, int *&segflags) {
    segments = new int[var.nseg*NODES_PER_FACET];
    segflags = new int[var.nseg];
    int seg_idx = 0;
    int flag_idx = 0;
#ifndef THREED
    for (int i = 0; i < var.nx - 1; ++i) {
        for (int j = 0; j < var.nz - 1; ++j) {
            if (i == 0) {
                segments[seg_idx] = j + i * var.nz;
                seg_idx++;
                segments[seg_idx] = j + i * var.nz + 1;
                seg_idx++;
                segflags[flag_idx] = 1;
                flag_idx++;
            }
            if (j == 0) {
                segments[seg_idx] = j + i * var.nz;
                seg_idx++;
                segments[seg_idx] = j + i * var.nz + var.nz;
                seg_idx++;
                segflags[flag_idx] = 32;
                flag_idx++;
            }
        }
    }
    for (int i = 1; i < var.nx; ++i) {
        segments[seg_idx] = i * var.nz - 1;
        seg_idx++;
        segments[seg_idx] = i * var.nz + var.nz - 1;
        seg_idx++;
        segflags[flag_idx] = 16;
        flag_idx++;
    }
    for (int j = 0; j < var.nz - 1; ++j) {
        segments[seg_idx] = var.nz * (var.nx - 1) + j;
        seg_idx++;
        segments[seg_idx] = var.nz * (var.nx - 1) + j + 1;
        seg_idx++;
        segflags[flag_idx] = 2;
        flag_idx++;
    }
#else
    int idx_elem, idx_cell;
    // top and bottom
    for (int i = 0; i < var.nx - 1; ++i) {
        for (int j = 0; j < var.ny - 1; ++j) {
            // top
            int k = 0;
            idx_cell = i * (var.ny - 1) * (var.nz - 1) + j * (var.nz - 1) + k;
            {
                idx_elem = idx_cell*5;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[2];
                segments[seg_idx*3+2] = conn[1];
                segflags[seg_idx] = BOUNDZ1;
                seg_idx++;
            }
            {
                idx_elem = idx_cell*5 + 1;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[2];
                segments[seg_idx*3+2] = conn[1];
                segflags[seg_idx] = BOUNDZ1;
                seg_idx++;
            }
            // bottom
            {
                k = var.nz - 2;
                idx_cell = i * (var.ny - 1) * (var.nz - 1) + j * (var.nz - 1) + k;
                idx_elem = idx_cell*5 + 2;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[1];
                segments[seg_idx*3+1] = conn[2];
                segments[seg_idx*3+2] = conn[3];
                segflags[seg_idx] = BOUNDZ0;
                seg_idx++;
            }
            {
                idx_elem = idx_cell*5 + 3;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[1];
                segments[seg_idx*3+1] = conn[2];
                segments[seg_idx*3+2] = conn[3];
                segflags[seg_idx] = BOUNDZ0;
                seg_idx++;
            }
        }
    }

    // left and right
    for (int j = 0; j < var.ny - 1; ++j) {
        for (int k = 0; k < var.nz - 1; ++k) {
            // left
            int i = 0;
            idx_cell = i * (var.ny - 1) * (var.nz - 1) + j * (var.nz - 1) + k;
            int order = (i+j+k)%2;

            if (order)
                idx_elem = idx_cell*5 + 1;
            else
                idx_elem = idx_cell*5;
            {
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[1];
                segments[seg_idx*3+2] = conn[3];
                segflags[seg_idx] = BOUNDX0;
                seg_idx++;
            }
            {
                idx_elem = idx_cell*5 + 2;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[2];
                segments[seg_idx*3+2] = conn[1];
                segflags[seg_idx] = BOUNDX0;
                seg_idx++;
            }

            // right
            i = var.nx - 2;
            order = (i+j+k)%2;
            idx_cell = i * (var.ny - 1) * (var.nz - 1) + j * (var.nz - 1) + k;
            if (order)
                idx_elem = idx_cell*5;
            else
                idx_elem = idx_cell*5 + 1;
            {
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[1];
                segments[seg_idx*3+1] = conn[2];
                segments[seg_idx*3+2] = conn[3];
                segflags[seg_idx] = BOUNDX1;
                seg_idx++;
            }
            {
                idx_elem = idx_cell*5 + 3;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[3];
                segments[seg_idx*3+2] = conn[2];
                segflags[seg_idx] = BOUNDX1;
                seg_idx++;
            }
        }
    }

    // front and back
    for (int i = 0; i < var.nx - 1; ++i) {
        for (int k = 0; k < var.nz - 1; ++k) {
            // front
            int j = 0;
            idx_cell = i * (var.ny - 1) * (var.nz - 1) + j * (var.nz - 1) + k;
            int order = (i+j+k)%2;
            idx_elem = idx_cell*5 + 1;
            {
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[3];
                segments[seg_idx*3+2] = conn[2];
                segflags[seg_idx] = BOUNDY0;
                seg_idx++;
            }

            if (order)
                idx_elem = idx_cell*5 + 3;
            else
                idx_elem = idx_cell*5 + 2;
            {
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[1];
                segments[seg_idx*3+2] = conn[3];
                segflags[seg_idx] = BOUNDY0;
                seg_idx++;
            }
 
            // back
            j = var.ny - 2;
            order = (i+j+k)%2;
            idx_cell = i * (var.ny - 1) * (var.nz - 1) + j * (var.nz - 1) + k;
            idx_elem = idx_cell*5;
            {
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                if (order) {
                    segments[seg_idx*3] = conn[0];
                    segments[seg_idx*3+1] = conn[1];
                    segments[seg_idx*3+2] = conn[3];
                } else {
                    segments[seg_idx*3] = conn[1];
                    segments[seg_idx*3+1] = conn[2];
                    segments[seg_idx*3+2] = conn[3];
                }
                segflags[seg_idx] = BOUNDY1;
                seg_idx++;
            }

            if (order) {
                idx_elem = idx_cell*5 + 2;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[3];
                segments[seg_idx*3+2] = conn[2];
            } else {
                idx_elem = idx_cell*5 + 3;
                ConstConnAccessor conn = (*var.connectivity)[idx_elem];
                segments[seg_idx*3] = conn[0];
                segments[seg_idx*3+1] = conn[2];
                segments[seg_idx*3+2] = conn[1];
            }
            segflags[seg_idx] = BOUNDY1;
            seg_idx++;
        }
    }
#endif
}

void create_regular_regattr(const Variables& var, double *&regattr) {
    regattr = new double[var.nelem];

    for (int i = 0; i < var.nelem; i++)
        regattr[i] = 0.;
}


void create_equilateral_node(const Param& param, const Variables& var, double *&points) {
    /*
    Create equilateral mesh nodes.
    
                     nx = 7
    0------1-----2-----3-----4-----5-----6
    | \  /  \  /  \  /  \  /  \  /  \  / |
    14-15----16----17----18----19----20-21
    | /  \  /  \  /  \  /  \  /  \  /  \ |   nz = 4
    7-----8-----9-----10----11----12----13
    | \  /  \  /  \  /  \  /  \  /  \  / |
    22-23----24----25----26----27----28-29

    var.nnode = 30
    var.nelem = 42
    */
    points = new double[var.nnode*2];
    double dx = param.mesh.resolution;
    double dz = -param.mesh.resolution * std::sqrt(3.0) / 2.;

    double bdy_dx = (param.mesh.xlength - (var.nx-1)*dx) / 2.;

    int ind = 0;
    for (int j=0; j<var.nz; j+=2) {
        int ind_tmp = ind;
        if (j == var.nz-1) {
            for (int i=0; i<var.nx; ++i) {
                points[ind_tmp*2+1] = -param.mesh.zlength;
                ind_tmp++;
            }
        } else {
            for (int i=0; i<var.nx; ++i) {
                points[ind_tmp*2+1] = j * dz;
                ind_tmp++;
            }
        }
        points[ind*2] = 0.;
        ind++;

        for (int i=1; i <var.nx-1; ++i) {
            points[ind*2] = i * dx + bdy_dx;
            ind++;
        }

        points[ind*2] = param.mesh.xlength;
        ind++;
    }

    for (int j=1; j<var.nz; j+=2) {
        int ind_tmp = ind;
        if (j == var.nz-1) {
            ind_tmp = ind;
            for (int i = 0; i < var.nx+1; ++i) {
                points[ind_tmp*2+1] = -param.mesh.zlength;
                ind_tmp++;
            }
        } else {
            for (int i = 0; i<var.nx+1; ++i) {
                points[ind_tmp*2+1] = j * dz;
                ind_tmp++;
            }
        }
        points[ind*2] = 0.;
        ind++;

        for (int i=0; i<var.nx-1; ++i) {
            points[ind*2] = (i+0.5) * dx + bdy_dx;
            ind++;
        }

        points[ind*2] = param.mesh.xlength;
        ind++;
    }
}


void new_mesh_regular_equilateral(const Param& param, Variables& var)
{
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

    double sqrt3_to_2 = 2./std::sqrt(3.0);

    double x_mid = param.mesh.xlength / 2;
    var.nx = int((x_mid-0.5*param.mesh.resolution) / param.mesh.resolution)*2 + 2;
    var.nz = int(param.mesh.zlength*sqrt3_to_2 / param.mesh.resolution) + 1;
    var.nnode = var.nx*((int)((var.nz-1)/2)+1) + (var.nx+1)*((int)((var.nz-1)/2)+(1-var.nz%2));
    var.nelem = (2*var.nx-1) * (var.nz-1);
    var.nseg = 2*(var.nz+var.nx-2) + 1 - var.nz%2;

    create_equilateral_node(param, var, pcoord);
    create_equilateral_elem(var, pconnectivity);
    create_equilateral_segments(var, psegment, psegflag);

    var.coord = new array_t(pcoord, var.nnode);
    var.connectivity = new conn_t(pconnectivity, var.nelem);
    var.segment = new segment_t(psegment, var.nseg);
    var.segflag = new segflag_t(psegflag, var.nseg);

    create_regular_regattr(var, pregattr);
    var.regattr = new regattr_t(pregattr, var.nelem);

    delete [] pcoord;
    delete [] pconnectivity;
    delete [] psegment;
    delete [] psegflag;
    delete [] pregattr;
}


void triangulate_polygon
(double min_angle, double max_area,
 int meshing_verbosity, double max_steiner_factor,
 int npoints, int nsegments,
 const double *points, const int *segments, const int *segflags,
 const int nregions, const double *regionattributes,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags, double **outregattr)
{
    char options[255];
    triangulateio in, out;

    std::string verbosity, vol, quality, steiner;
    set_verbosity_str(verbosity, meshing_verbosity);
    set_volume_str(vol, max_area);
    set_2d_quality_str(quality, min_angle);
    set_steiner_str(steiner, max_steiner_factor, npoints, max_area);

    if( nregions > 0 )
        std::sprintf(options, "%s%spjz%s%sA", verbosity.c_str(), quality.c_str(), vol.c_str(), steiner.c_str());
    else
        std::sprintf(options, "%s%spjz%s%s", verbosity.c_str(), quality.c_str(), vol.c_str(), steiner.c_str());

    if( meshing_verbosity >= 0 )
        std::cout << "The meshing option is: " << options << '\n';

    in.pointlist = const_cast<double*>(points);
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    in.numberofpoints = npoints;
    in.numberofpointattributes = 0;

    in.trianglelist = NULL;
    in.triangleattributelist = NULL;
    in.trianglearealist = NULL;
    in.numberoftriangles = 0;
    in.numberofcorners = 3;
    in.numberoftriangleattributes = 0;

    in.segmentlist = const_cast<int*>(segments);
    in.segmentmarkerlist = const_cast<int*>(segflags);
    in.numberofsegments = nsegments;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.numberofregions = nregions;
    if( nregions > 0 )
        in.regionlist = const_cast<double*>(regionattributes);
    else
        in.regionlist = NULL;

    out.pointlist = NULL;
    out.pointattributelist = NULL;
    out.pointmarkerlist = NULL;
    out.trianglelist = NULL;
    out.triangleattributelist = NULL;
    out.neighborlist = NULL;
    out.segmentlist = NULL;
    out.segmentmarkerlist = NULL;
    out.edgelist = NULL;
    out.edgemarkerlist = NULL;

    /*******************************/
    triangulate(options, &in, &out, NULL);
    /*******************************/

    *noutpoints = out.numberofpoints;
    *outpoints = out.pointlist;

    *ntriangles = out.numberoftriangles;
    *triangles = out.trianglelist;

    *noutsegments = out.numberofsegments;
    *outsegments = out.segmentlist;
    *outsegflags = out.segmentmarkerlist;
    *outregattr = out.triangleattributelist;

    trifree(out.pointmarkerlist);
}


void set_3d_quality_str(std::string &quality, double max_ratio,
                        double min_dihedral_angle, double max_dihedral_angle)
{
    quality.clear();
    if (max_ratio > 0) {
        quality += 'q';
        quality += std::to_string((long double)max_ratio);
        quality += "qq";
        quality += std::to_string((long double)min_dihedral_angle);
        quality += "qqq";
        quality += std::to_string((long double)max_dihedral_angle);
    }
}

#ifdef USEMMG
struct CentroidCloud {
    const double *pts;  // AoS layout: x0,y0,z0, x1,y1,z1, ...
    int n;
    CentroidCloud(const double *p, int n_) : pts(p), n(n_) {}
    inline size_t kdtree_get_point_count() const { return (size_t)n; }
    inline double kdtree_get_pt(size_t i, size_t d) const {
        return pts[i * NDIMS + d];
    }
    template<class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};
using CentroidKDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, CentroidCloud>,
    CentroidCloud, NDIMS>;


static void propagate_regattr_nearest(
    // Coarse mesh (TetGen output, AoS layout)
    int coarse_nelem,
    const double *coarse_coord,       // [coarse_nnode × NDIMS]
    const int    *coarse_conn,        // [coarse_nelem × NODES_PER_ELEM], 0-indexed
    const double *coarse_regattr,     // [coarse_nelem], per-element material
    // Fine mesh (MMG output, AoS layout)
    int fine_nelem,
    const double *fine_coord,         // [fine_nnode × NDIMS]
    const int    *fine_conn,          // [fine_nelem × NODES_PER_ELEM], 0-indexed
    // Output
    double *fine_regattr)             // [fine_nelem], caller-allocated
{
    // 1. Build coarse element centroid array (AoS)
    std::vector<double> coarse_cents(coarse_nelem * NDIMS, 0.0);
    for (int e = 0; e < coarse_nelem; ++e) {
        for (int n = 0; n < NODES_PER_ELEM; ++n)
            for (int d = 0; d < NDIMS; ++d)
                coarse_cents[e*NDIMS+d] +=
                    coarse_coord[coarse_conn[e*NODES_PER_ELEM+n]*NDIMS+d];
        for (int d = 0; d < NDIMS; ++d)
            coarse_cents[e*NDIMS+d] /= NODES_PER_ELEM;
    }

    // 2. Build KD-tree over coarse centroids
    CentroidCloud cloud(coarse_cents.data(), coarse_nelem);
    CentroidKDTree kdtree(NDIMS, cloud,
        nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdtree.buildIndex();

    // 3. For each fine element: find nearest coarse centroid and copy regattr
    std::vector<double> fine_cent(NDIMS);
    for (int e = 0; e < fine_nelem; ++e) {
        fine_cent.assign(NDIMS, 0.0);
        for (int n = 0; n < NODES_PER_ELEM; ++n)
            for (int d = 0; d < NDIMS; ++d)
                fine_cent[d] += fine_coord[fine_conn[e*NODES_PER_ELEM+n]*NDIMS+d];
        for (int d = 0; d < NDIMS; ++d)
            fine_cent[d] /= NODES_PER_ELEM;

        std::size_t nearest_idx;
        double dist_sq;
        nanoflann::KNNResultSet<double> result(1);
        result.init(&nearest_idx, &dist_sq);
        kdtree.findNeighbors(result, fine_cent.data(),
            nanoflann::SearchParameters());

        fine_regattr[e] = coarse_regattr[nearest_idx];
    }
}


// Measure (area in 2D, volume in 3D) of element e of a raw mesh.
static double elem_measure(const double *coord, const int *conn, int e)
{
    const int *c = conn + e * NODES_PER_ELEM;
    const double *a = coord + c[0] * NDIMS;
    const double *b = coord + c[1] * NDIMS;
    const double *d = coord + c[2] * NDIMS;
#ifdef THREED
    const double *f = coord + c[3] * NDIMS;
    double ab[3], ad[3], af[3];
    for (int i = 0; i < 3; ++i) {
        ab[i] = b[i] - a[i];
        ad[i] = d[i] - a[i];
        af[i] = f[i] - a[i];
    }
    return std::fabs(ab[0] * (ad[1]*af[2] - ad[2]*af[1])
                   - ab[1] * (ad[0]*af[2] - ad[2]*af[0])
                   + ab[2] * (ad[0]*af[1] - ad[1]*af[0])) / 6.0;
#else
    return 0.5 * std::fabs((b[0] - a[0]) * (d[1] - a[1])
                         - (d[0] - a[0]) * (b[1] - a[1]));
#endif
}

static void compute_init_metric(
    int nnode, int nelem,
    const double *coord,           // node coordinates [nnode × NDIMS]
    const int *conn,               // coarse connectivity [nelem × NODES_PER_ELEM]
    const double *elem_region,     // coarse pregattr [nelem]: region indices when per-region, ignored when uniform
    int n_regions, const double *regattr,  // original region attrs (NDIMS+2 per region)
    double max_elem_size,          // global volume constraint (>0 if uniform, else 0)
    double resolution,
    std::vector<double> &metric)   // output [nnode]
{
    metric.resize(nnode);

    if (max_elem_size > 0) {
        // Uniform target: convert volume constraint to edge length
        // Convention in DES: max_elem_size ≈ sizefactor * resolution^NDIMS
        double target = std::pow(max_elem_size / sizefactor, 1.0 / NDIMS);
        std::fill(metric.begin(), metric.end(), target);
    }
    else {
        // Per-region: connectivity-based projection from coarse elements to nodes.
        // elem_region[e] is the region index (0..n_regions-1) for element e,
        // as stamped by TetGen's -A flag using the temporarily encoded region index.
        // Nodes shared by elements of different regions get an averaged metric,
        // providing a smooth transition at zone boundaries.
        const int attr_stride = NDIMS + 2;
        std::vector<int> count(nnode, 0);
        std::fill(metric.begin(), metric.end(), 0.0);

        for (int e = 0; e < nelem; ++e) {
            int r = (int)elem_region[e];
            if (r < 0 || r >= n_regions) r = 0;
            double vol = regattr[r * attr_stride + NDIMS + 1];
            // A region without a size constraint (.poly size <= 0) keeps each element's as-meshed
            // size, as on the non-MMG path (a mesh.resolution fallback would refine the region).
            double target = (vol > 0)
                ? std::pow(vol / sizefactor, 1.0 / NDIMS)
                : std::pow(elem_measure(coord, conn, e) / sizefactor, 1.0 / NDIMS);

            for (int k = 0; k < NODES_PER_ELEM; ++k) {
                int ni = conn[e * NODES_PER_ELEM + k];
                metric[ni] += target;
                count[ni]++;
            }
        }

        for (int i = 0; i < nnode; ++i)
            metric[i] = (count[i] > 0) ? metric[i] / count[i] : resolution;
    }
}

#ifdef THREED
static void mmg_refine_init_mesh_3d(
    const Mesh &mesh,
    double max_elem_size,        // original TetGen target (for metric)
    int n_regions, const double *regattr,   // original region attrs (for metric)
    // In/out mesh arrays (replaced with fine arrays on return)
    int &nnode, int &nelem, int &nseg,
    double *&pcoord,    // [nnode × 3], AoS
    int    *&pconn,     // [nelem × 4], 0-indexed
    int    *&pseg,      // [nseg × 3], 0-indexed
    int    *&psegflag,  // [nseg]
    double *&pregattr,  // [nelem]
    double_vec *fine_init_metric_n)  // optional: per-node desired edge length on fine mesh
{
    // --- Save coarse mesh for regattr propagation ---
    const int cn = nnode, ce = nelem, cs = nseg;
    double_vec cc(pcoord, pcoord + cn * NDIMS);
    int_vec    ck(pconn,  pconn  + ce * NODES_PER_ELEM);
    double_vec cra(pregattr ? pregattr : nullptr,
                            pregattr ? pregattr + ce : nullptr);

    // --- Metric on the coarse mesh (uniform or per-region) ---
    double_vec metric;
    compute_init_metric(cn, ce, cc.data(), ck.data(), cra.data(),
                        n_regions, regattr, max_elem_size, mesh.resolution, metric);

    // --- Adapt the coarse mesh (mmg_adapt) ---
    MMGInput  mmg_in = { cn, ce, cs, pcoord, pconn, pseg, psegflag, metric.data(), nullptr, nullptr, true, nullptr };
    MMGOutput mmg_out;
    mmg_adapt(mesh, mmg_in, mmg_out);

    // --- Free old coarse arrays; reallocate outputs from the adapted mesh ---
    delete [] pcoord;   delete [] pconn;   delete [] pseg;   delete [] psegflag;
    if (pregattr) delete [] pregattr;
    nnode = mmg_out.nnode;  nelem = mmg_out.nelem;  nseg = mmg_out.nseg;
    pcoord   = new double[nnode * NDIMS];
    pconn    = new int   [nelem * NODES_PER_ELEM];
    pseg     = new int   [nseg  * NODES_PER_FACET];
    psegflag = new int   [nseg];
    pregattr = new double[nelem];
    std::memcpy(pcoord,   mmg_out.coord.data(),   mmg_out.coord.size()   * sizeof(double));
    std::memcpy(pconn,    mmg_out.conn.data(),    mmg_out.conn.size()    * sizeof(int));
    std::memcpy(pseg,     mmg_out.seg.data(),     mmg_out.seg.size()     * sizeof(int));
    std::memcpy(psegflag, mmg_out.segflag.data(), mmg_out.segflag.size() * sizeof(int));

    // --- Compute fine-mesh init metric hint (for init_elem_size_n) ---
    // Propagate coarse REGION INDICES (from cra) to fine elements, then project to nodes.
    // This gives a zone-aware per-node desired edge length free of MMG's transition smoothing.
    if (!cra.empty()) {
        if (max_elem_size <= 0 && n_regions > 0) {
            // Per-region case: propagate region indices to fine elements, then compute metric
            std::vector<double> fine_region_idx(nelem);
            propagate_regattr_nearest(ce, cc.data(), ck.data(), cra.data(),
                                        nelem, pcoord, pconn, fine_region_idx.data());
            compute_init_metric(nnode, nelem, pcoord, pconn, fine_region_idx.data(),
                                n_regions, regattr, max_elem_size, mesh.resolution, *fine_init_metric_n);
        } else {
            // Uniform case: exact constant metric
            const double target = std::pow(max_elem_size / sizefactor, 1.0 / NDIMS);
            fine_init_metric_n->assign(nnode, target);
        }
    }

    // --- Propagate regattr via nearest-centroid from coarse ---
    // mmg_adapt takes element references in 2D only, so a zone boundary here lands on the bisector
    // of two coarse centroids, up to about half a coarse element off, as the hint's above does.
    // When per-region volumes were used (max_elem_size<=0 and n_regions>0),
    // cra contains region indices (temporarily encoded by the use_mmg_init coarse pass of
    // points_to_mesh or of the 3D .poly path).
    // Decode back to actual mattypes before propagation.
    if (!cra.empty()) {
        const double *mattype_src = cra.data();
        std::vector<double> coarse_actual_mattype;
        if (max_elem_size <= 0 && n_regions > 0) {
            const int attr_stride = NDIMS + 2;
            coarse_actual_mattype.resize(ce);
            for (int e = 0; e < ce; ++e)
                coarse_actual_mattype[e] = regattr[(int)cra[e] * attr_stride + NDIMS];
            mattype_src = coarse_actual_mattype.data();
        }
        propagate_regattr_nearest(ce, cc.data(), ck.data(), mattype_src,
                                    nelem, pcoord, pconn, pregattr);
    }
    else
        std::fill(pregattr, pregattr + nelem, 0.0);

    std::cerr << "MMG init refinement done: "
                << nnode << " nodes, " << nelem << " elements\n";
}
#endif // THREED

#ifndef THREED
static void mmg_refine_init_mesh_2d(
    const Mesh &mesh,
    bool keep_regions,           // mattype_option 0: the coarse mattype boundaries stay put
    int npoints, const double *points,   // the input vertices; those on a mattype
                                         // boundary become MMG corners
    double max_elem_size,        // original TetGen target (for metric)
    int n_regions, const double *regattr,   // original region attrs (for metric)
    // In/out mesh arrays (replaced with fine arrays on return)
    int &nnode, int &nelem, int &nseg,
    double *&pcoord,    // [nnode × 2], AoS
    int    *&pconn,     // [nelem × 3], 0-indexed
    int    *&pseg,      // [nseg × 2], 0-indexed
    int    *&psegflag,  // [nseg]
    double *&pregattr,  // [nelem]
    double_vec *fine_init_metric_n)  // optional: per-node desired edge length on fine mesh
{
    // --- Save coarse mesh for regattr propagation ---
    const int cn = nnode, ce = nelem, cs = nseg;
    double_vec cc(pcoord, pcoord + cn * NDIMS);
    int_vec    ck(pconn,  pconn  + ce * NODES_PER_ELEM);
    double_vec cra(pregattr ? pregattr : nullptr,
                            pregattr ? pregattr + ce : nullptr);

    // --- Metric on the coarse mesh (uniform or per-region) ---
    double_vec metric;
    compute_init_metric(cn, ce, cc.data(), ck.data(), cra.data(),
                        n_regions, regattr, max_elem_size, mesh.resolution, metric);

    // --- Coarse mattypes: per-region volumes stamped region indices, decoded here ---
    double_vec cmat(cra);
    if (!cra.empty() && max_elem_size <= 0 && n_regions > 0)
        for (int e = 0; e < ce; ++e) cmat[e] = regattr[(int)cra[e] * (NDIMS + 2) + NDIMS];
    // Under mattype_option 0, regattr's only reader, the mattype (truncated, as markerset reads it)
    // rides through MMG as the element reference, and each input vertex on a mattype boundary is
    // an MMG corner, so MMG keeps the boundary on its straight segments between them.
    std::vector<int> eref;
    std::vector<char> corner;
    if (keep_regions && !cmat.empty()) {
        eref.resize(ce);
        for (int e = 0; e < ce; ++e) eref[e] = (int)cmat[e];
        std::vector<int> node_ref(cn, 0);   // the reference of one element through the node
        for (int e = 0; e < ce; ++e)
            for (int k = 0; k < NODES_PER_ELEM; ++k) node_ref[ck[e*NODES_PER_ELEM + k]] = eref[e];
        std::vector<char> on_iface(cn, 0);  // another element through it has another reference
        for (int e = 0; e < ce; ++e)
            for (int k = 0; k < NODES_PER_ELEM; ++k) {
                const int n = ck[e*NODES_PER_ELEM + k];
                if (node_ref[n] != eref[e]) on_iface[n] = 1;
            }
        typedef std::array<double, NDIMS> Point;
        auto point_at = [](const double *x) -> Point {
            Point p; std::copy(x, x + NDIMS, p.begin()); return p; };
        std::set<Point> vertex;   // exact: Triangle copies the input coordinates
        for (int i = 0; i < npoints; ++i) vertex.insert(point_at(points + i*NDIMS));
        corner.assign(cn, 0);
        for (int n = 0; n < cn; ++n)
            corner[n] = on_iface[n] && vertex.count(point_at(&cc[n*NDIMS]));
        const long n_corner = std::count(corner.begin(), corner.end(), 1);
        if (n_corner)
            std::cout << "      [corner] " << n_corner << " input vertex(es) on a material "
                         "interface of the initial mesh handed to MMG as corners.\n";
    }

    // --- Adapt the coarse mesh (mmg_adapt) ---
    MMGInput  mmg_in = { cn, ce, cs, pcoord, pconn, pseg, psegflag, metric.data(), nullptr, nullptr,
                         true, nullptr, nullptr, 0, nullptr, eref.empty() ? nullptr : eref.data(),
                         corner.empty() ? nullptr : corner.data() };
    MMGOutput mmg_out;
    mmg_adapt(mesh, mmg_in, mmg_out);

    // --- Free old coarse arrays; reallocate outputs from the adapted mesh ---
    delete [] pcoord;   delete [] pconn;   delete [] pseg;   delete [] psegflag;
    if (pregattr) delete [] pregattr;
    nnode = mmg_out.nnode;  nelem = mmg_out.nelem;  nseg = mmg_out.nseg;
    pcoord   = new double[nnode * NDIMS];
    pconn    = new int   [nelem * NODES_PER_ELEM];
    pseg     = new int   [nseg  * NODES_PER_FACET];
    psegflag = new int   [nseg];
    pregattr = new double[nelem];
    std::memcpy(pcoord,   mmg_out.coord.data(),   mmg_out.coord.size()   * sizeof(double));
    std::memcpy(pconn,    mmg_out.conn.data(),    mmg_out.conn.size()    * sizeof(int));
    std::memcpy(pseg,     mmg_out.seg.data(),     mmg_out.seg.size()     * sizeof(int));
    std::memcpy(psegflag, mmg_out.segflag.data(), mmg_out.segflag.size() * sizeof(int));

    // --- Compute fine-mesh init metric hint (for init_elem_size_n) ---
    if (fine_init_metric_n && !cra.empty()) {
        if (max_elem_size <= 0 && n_regions > 0) {
            // The hint's zone boundaries sit on the bisector of two coarse centroids (up to about
            // half a coarse element off); regattr's are exact under mattype_option 0 (eref).
            std::vector<double> fine_region_idx(nelem);
            propagate_regattr_nearest(ce, cc.data(), ck.data(), cra.data(),
                                        nelem, pcoord, pconn, fine_region_idx.data());
            compute_init_metric(nnode, nelem, pcoord, pconn, fine_region_idx.data(),
                                n_regions, regattr, 0.0, mesh.resolution, *fine_init_metric_n);
        } else {
            const double target = std::pow(max_elem_size / sizefactor, 1.0 / NDIMS);
            fine_init_metric_n->assign(nnode, target);
        }
    }

    // --- Propagate regattr: the MMG element references when set, else nearest coarse centroid ---
    if (!mmg_out.elem_ref.empty())
        for (int e = 0; e < nelem; ++e) pregattr[e] = mmg_out.elem_ref[e];
    else if (!cra.empty())
        propagate_regattr_nearest(ce, cc.data(), ck.data(), cmat.data(),
                                    nelem, pcoord, pconn, pregattr);
    else
        std::fill(pregattr, pregattr + nelem, 0.0);

    std::cerr << "MMG2D init refinement done: "
                << nnode << " nodes, " << nelem << " elements\n";
}

#endif // !THREED

// Per-node element size of a raw mesh: esize[e] = the edge of the equilateral element of e's measure;
// size_n = its measure-weighted average over each node's elements (initialize_elem_size_n's
// convention).
static void measured_size_n(int nnode, int nelem, const double *pcoord, const int *pconn,
                            std::vector<double> &esize, std::vector<double> &size_n)
{
    std::vector<double> wsum(nnode, 0.0);
    size_n.assign(nnode, 0.0);
    esize.assign(nelem, 0.0);
    for (int e = 0; e < nelem; ++e) {
        const int *c = pconn + e*NODES_PER_ELEM;
        const double measure = elem_measure(pcoord, pconn, e);
        esize[e] = std::pow(measure / sizefactor, 1.0 / NDIMS);
        for (int k = 0; k < NODES_PER_ELEM; ++k) {
            size_n[c[k]] += esize[e] * measure;
            wsum [c[k]] += measure;
        }
    }
    for (int n = 0; n < nnode; ++n)
        if (wsum[n] > 0) size_n[n] /= wsum[n];
}

// out[n] = old_val at the old node nearest to new node n.
static void transfer_nearest_node(int old_nnode, const double *old_coord,
                                  const std::vector<double> &old_val,
                                  int new_nnode, const double *new_coord, double_vec &out)
{
    CentroidCloud cloud(old_coord, old_nnode);   // plain AoS point cloud; works for nodes too
    CentroidKDTree kdtree(NDIMS, cloud,
        nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdtree.buildIndex();
    out.assign(new_nnode, 0.0);
    for (int n = 0; n < new_nnode; ++n) {
        std::size_t nearest_idx;
        double dist_sq;
        nanoflann::KNNResultSet<double> result(1);
        result.init(&nearest_idx, &dist_sq);
        kdtree.findNeighbors(result, &new_coord[(std::size_t)n*NDIMS], nanoflann::SearchParameters());
        out[n] = old_val[nearest_idx];
    }
}

// Adopt the output of an MMG pass over the raw initial mesh: var.init_elem_size_n = the
// nearest-old-node transfer of hint; the regions from the MMG element references when set (2D
// only), else from the nearest old element centroid; the mesh buffers are replaced.
static void adopt_initial_mmg_output(Variables &var, const std::vector<double> &hint,
        const MMGOutput &mmg_out, int &nnode, int &nelem, int &nseg,
        double *&pcoord, int *&pconn, int *&pseg, int *&psegflag, double *&pregattr)
{
    transfer_nearest_node(nnode, pcoord, hint, mmg_out.nnode, mmg_out.coord.data(),
                          *var.init_elem_size_n);

    double *fine_regattr = new double[mmg_out.nelem];
    if (!mmg_out.elem_ref.empty())
        for (int e = 0; e < mmg_out.nelem; ++e) fine_regattr[e] = mmg_out.elem_ref[e];
    else if (pregattr)
        propagate_regattr_nearest(nelem, pcoord, pconn, pregattr,
                                  mmg_out.nelem, mmg_out.coord.data(), mmg_out.conn.data(),
                                  fine_regattr);
    else
        std::fill(fine_regattr, fine_regattr + mmg_out.nelem, 0.0);

    delete [] pcoord;   delete [] pconn;   delete [] pseg;   delete [] psegflag;
    if (pregattr) delete [] pregattr;
    nnode = mmg_out.nnode;  nelem = mmg_out.nelem;  nseg = mmg_out.nseg;
    pcoord   = new double[(std::size_t)nnode * NDIMS];
    pconn    = new int   [(std::size_t)nelem * NODES_PER_ELEM];
    pseg     = new int   [(std::size_t)nseg  * NODES_PER_FACET];
    psegflag = new int   [nseg];
    pregattr = fine_regattr;
    std::memcpy(pcoord,   mmg_out.coord.data(),   mmg_out.coord.size()   * sizeof(double));
    std::memcpy(pconn,    mmg_out.conn.data(),    mmg_out.conn.size()    * sizeof(int));
    std::memcpy(pseg,     mmg_out.seg.data(),     mmg_out.seg.size()     * sizeof(int));
    std::memcpy(psegflag, mmg_out.segflag.data(), mmg_out.segflag.size() * sizeof(int));
}

// Weak-zone plastic strain of a raw mesh's elements (initial_weak_zone), zero outside the zone. Heap
// copies of the mesh and the result: initial_weak_zone's loop is an OpenACC kernel.
static double_vec *raw_weak_zone(const Param &param, int nnode, int nelem, const double *pcoord,
                                 const int *pconn)
{
    array_t *coord = new array_t(pcoord, nnode);
    conn_t *conn = new conn_t(pconn, nelem);
    double_vec *pls = new double_vec(nelem, 0.0);
    initial_weak_zone(param, nelem, *coord, *conn, *pls);
    delete coord;
    delete conn;
    return pls;
}

#ifdef THREED
// Six times the signed volume of tetrahedron (a, b, c, d).
static double signed_vol6(const double *a, const double *b, const double *c, const double *d)
{
    const double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
    const double v[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
    const double w[3] = { d[0] - a[0], d[1] - a[1], d[2] - a[2] };
    return u[0] * (v[1]*w[2] - v[2]*w[1]) - u[1] * (v[0]*w[2] - v[2]*w[0])
         + u[2] * (v[0]*w[1] - v[1]*w[0]);
}

// 3D regions after an initial MMG pass: the nearest old centroid (adopt_initial_mmg_output) can
// lie across a region boundary beside a flat element, so each element takes the region of the old
// element containing its centroid, searched among the 64 nearest old centroids (heuristic: the
// containing one ranked at most 26 in the meshes measured); none found keeps the nearest one's.
static void regions_by_containing_element(const std::vector<double> &old_coord,
                                          const std::vector<int> &old_conn,
                                          const std::vector<double> &old_reg,
                                          int nelem, const double *pcoord, const int *pconn,
                                          double *pregattr)
{
    const int old_nelem = (int)(old_conn.size() / NODES_PER_ELEM);
    std::vector<double> cents((std::size_t)old_nelem * NDIMS, 0.0);
    for (int e = 0; e < old_nelem; ++e)
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            for (int d = 0; d < NDIMS; ++d)
                cents[e*NDIMS + d] += old_coord[old_conn[e*NODES_PER_ELEM + i]*NDIMS + d]
                                      / NODES_PER_ELEM;
    CentroidCloud cloud(cents.data(), old_nelem);
    CentroidKDTree kdtree(NDIMS, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdtree.buildIndex();
    const int k = std::min(64, old_nelem);
    std::vector<std::size_t> idx(k);
    std::vector<double> dist_sq(k);
    for (int e = 0; e < nelem; ++e) {
        const int *c = pconn + e*NODES_PER_ELEM;
        double p[NDIMS] = { 0 };
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            for (int d = 0; d < NDIMS; ++d) p[d] += pcoord[c[i]*NDIMS + d] / NODES_PER_ELEM;
        nanoflann::KNNResultSet<double> result(k);
        result.init(idx.data(), dist_sq.data());
        kdtree.findNeighbors(result, p, nanoflann::SearchParameters());
        for (std::size_t j = 0; j < result.size(); ++j) {
            const int o = (int)idx[j];
            const double *x[4];
            for (int i = 0; i < 4; ++i) x[i] = &old_coord[old_conn[o*NODES_PER_ELEM + i]*NDIMS];
            const double v = signed_vol6(x[0], x[1], x[2], x[3]);
            // inside (or on a face) when replacing any vertex by p keeps the volume's sign
            bool inside = true;
            for (int i = 0; i < 4 && inside; ++i) {
                const double *y[4] = { x[0], x[1], x[2], x[3] };
                y[i] = p;
                // a barycentric coordinate; its round-off is ~1e-16 x the element's aspect, so
                // -1e-12 admits points on a face (both sides then match; the first wins)
                inside = signed_vol6(y[0], y[1], y[2], y[3]) / v >= -1e-12;
            }
            if (inside) { pregattr[e] = old_reg[o]; break; }
        }
    }
}
#endif

// Refine the raw initial mesh's weak zone as the first MMG remesh's R5 would, before any field
// exists, so that remesh does not remap the prescribed plastic strain. Free: a zone element that
// is r5_oversized against the largest of its nodes' floored plastic_size_targets, with no
// top-surface node and no node between two regions (R5's other skips do not apply). Its nodes
// move, except those on the bottom, the top, a restored side wall (option 13) or a domain edge or
// corner, and those of a region-boundary element. MMG rebuilds the free elements and their hinges
// at the first remesh's metric, and a region the quality gate frees around a bad output element,
// but no region-boundary element: a pass whose gate frees one is discarded. So the regions come
// back exactly (2D element references, 3D regions_by_containing_element). A pass that leaves no
// fewer free elements is undone. The size base (the use_mmg_init hint, else the measured size) stays
// var.init_elem_size_n, so the refinement does not count twice.
static void refine_initial_weak_zone(const Param &param, Variables &var,
        int &nnode, int &nelem, int &nseg,
        double *&pcoord, int *&pconn, int *&pseg, int *&psegflag, double *&pregattr)
{
    if (!param.control.has_moving_mesh || param.mesh.mmg_remesh_size_recovery_ratio <= 0) return;
    if (!pregattr)   // Triangle/TetGen stamp a region attribute (nregions >= 1), MMG copies it
        die(EXIT_INTERNAL_ASSERT, "refine_initial_weak_zone: the initial mesh has no regions");
    const double floor = mmg_metric_floor(param.mesh);
    const int max_passes = 2;   // heuristic: the second takes what MMG's gradation left
    // The mesh before the last pass, for the undo (MMG's coordinate scaling rounds every node).
    int prev_nfree = -1, prev_ntop = 0, prev_nmixed = 0;
    std::vector<double> prev_coord, prev_regattr, prev_size;
    std::vector<int> prev_conn, prev_seg, prev_segflag;

    for (int pass = 1; ; ++pass) {
        std::vector<double> size_n, esize;
        measured_size_n(nnode, nelem, pcoord, pconn, esize, size_n);
        const double_vec hint = ((int)var.init_elem_size_n->size() == nnode)
                                ? *var.init_elem_size_n : size_n;

        double_vec *pls = raw_weak_zone(param, nnode, nelem, pcoord, pconn);
        double_vec volume(nelem), target;
        for (int e = 0; e < nelem; ++e) volume[e] = elem_measure(pcoord, pconn, e);
        const conn_t conn(pconn, nelem);
        plastic_size_targets(param.mesh.mmg_metric_refine_coeff, conn, volume, *pls, hint,
                             nnode, nelem, esize, target);
        for (int n = 0; n < nnode; ++n) target[n] = std::max(target[n], floor);

        std::vector<uint> bflag(nnode, 0);
        for (int s = 0; s < nseg; ++s)
            for (int i = 0; i < NODES_PER_FACET; ++i)
                bflag[pseg[s*NODES_PER_FACET + i]] |= psegflag[s];
        std::vector<char> mixed(nnode, 0), seen(nnode, 0), on_iface(nelem, 0);
        std::vector<double> node_reg(nnode, 0.0);
        for (int e = 0; e < nelem; ++e)
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                const int n = pconn[e*NODES_PER_ELEM + i];
                if (!seen[n]) { seen[n] = 1; node_reg[n] = pregattr[e]; }
                else if (node_reg[n] != pregattr[e]) mixed[n] = 1;
            }
        for (int e = 0; e < nelem; ++e)
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                if (mixed[pconn[e*NODES_PER_ELEM + i]]) on_iface[e] = 1;

        std::vector<char> free_elem(nelem, 0);
        int nfree = 0, ntop = 0, nmixed = 0;
        for (int e = 0; e < nelem; ++e) {
            if ((*pls)[e] <= 0) continue;
            const int *c = pconn + e*NODES_PER_ELEM;
            double tmax = 0;
            bool on_top = false;
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                tmax = std::max(tmax, target[c[i]]);
                on_top = on_top || (bflag[c[i]] & BOUNDZ1);
            }
            if (!r5_oversized(esize[e], tmax)) continue;
            if (on_top) { ++ntop; continue; }
            if (on_iface[e]) { ++nmixed; continue; }
            free_elem[e] = 1;
            ++nfree;
        }
        // what stays due at t = 0, logged on every return
        auto report = [](const char *why, int f, int t, int m) {
            if (f + t + m)
                std::cout << "  [weakzone] " << why << "due at t = 0: " << f << " free, " << t
                          << " on a top-surface node, " << m << " on a region boundary.\n";
        };
        if (prev_nfree >= 0 && nfree >= prev_nfree) {
            delete pls;
            delete [] pcoord;   delete [] pconn;   delete [] pseg;   delete [] psegflag;
            delete [] pregattr;
            nnode = (int)(prev_coord.size() / NDIMS);
            nelem = (int)(prev_conn.size() / NODES_PER_ELEM);
            nseg = (int)prev_segflag.size();
            pcoord = new double[prev_coord.size()];
            pconn = new int[prev_conn.size()];
            pseg = new int[prev_seg.size()];
            psegflag = new int[nseg];
            pregattr = new double[nelem];
            std::copy(prev_coord.begin(), prev_coord.end(), pcoord);
            std::copy(prev_conn.begin(), prev_conn.end(), pconn);
            std::copy(prev_seg.begin(), prev_seg.end(), pseg);
            std::copy(prev_segflag.begin(), prev_segflag.end(), psegflag);
            std::copy(prev_regattr.begin(), prev_regattr.end(), pregattr);
            *var.init_elem_size_n = prev_size;
            std::cout << "  [weakzone] pass " << pass - 1 << " undone: it left " << nfree
                      << " free element(s).\n";
            report("", prev_nfree, prev_ntop, prev_nmixed);
            return;
        }
        if (!nfree || pass > max_passes) {
            delete pls;
            report(nfree ? "after the passes, " : "", nfree, ntop, nmixed);
            return;
        }

        // The first remesh's metric: plastic_size_target (remeshing.cxx), floored; keep them equal.
        std::vector<double> metric(nnode, 0.0), vol_n(nnode, 0.0);
        for (int e = 0; e < nelem; ++e)
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                const int n = pconn[e*NODES_PER_ELEM + i];
                metric[n] += volume[e] / (1.0 + param.mesh.mmg_metric_refine_coeff * (*pls)[e]);
                vol_n[n] += volume[e];
            }
        delete pls;
        for (int n = 0; n < nnode; ++n)
            metric[n] = std::max(hint[n] * std::pow(metric[n] / vol_n[n], 1.0 / NDIMS), floor);

        uint pinned = BOUNDZ0 | BOUNDZ1;
        if (param.mesh.remeshing_option == 13)
            for (int t = 0; t < NSIDEWALL; ++t) pinned |= SIDEWALL_FLAG[t];
        std::vector<char> req_node(nnode, 1), req_elem(nelem, 1);
        for (int e = 0; e < nelem; ++e)
            if (free_elem[e])
                for (int i = 0; i < NODES_PER_ELEM; ++i) {
                    const int n = pconn[e*NODES_PER_ELEM + i];
                    const uint b = bflag[n] & BOUND_ANY;
                    if (!(b & pinned) && !(b & (b - 1))) req_node[n] = 0;
                }
        for (int e = 0; e < nelem; ++e)
            if (on_iface[e])
                for (int i = 0; i < NODES_PER_ELEM; ++i) req_node[pconn[e*NODES_PER_ELEM + i]] = 1;
        for (int e = 0; e < nelem; ++e)
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                if (!req_node[pconn[e*NODES_PER_ELEM + i]]) req_elem[e] = 0;
        const long nrebuild = std::count(req_elem.begin(), req_elem.end(), 0);
#ifndef THREED
        std::vector<int> eref(nelem);
        for (int e = 0; e < nelem; ++e) eref[e] = (int)pregattr[e];
        MMGInput  mmg_in = { nnode, nelem, nseg, pcoord, pconn, pseg, psegflag,
                             metric.data(), req_node.data(), req_elem.data(), true, nullptr,
                             nullptr, 0, nullptr, eref.data() };
#else
        MMGInput  mmg_in = { nnode, nelem, nseg, pcoord, pconn, pseg, psegflag,
                             metric.data(), req_node.data(), req_elem.data(), true, nullptr };
#endif
        MMGOutput mmg_out;
        // Gated as refine_initial_side_walls is; the gate edits req_node/req_elem in place.
        mmg_adapt_quality_gated(param, mmg_in, mmg_out, req_node, req_elem);
        bool iface_freed = false;   // the gate drops req_elem wherever it frees a node
        for (int e = 0; e < nelem && !iface_freed; ++e)
            iface_freed = on_iface[e] && !req_elem[e];
        if (iface_freed) {
            std::cout << "  [weakzone] pass " << pass << " discarded: the quality gate freed an "
                         "element on a region boundary.\n";
            report("", nfree, ntop, nmixed);
            return;
        }

        prev_nfree = nfree;
        prev_ntop = ntop;
        prev_nmixed = nmixed;
        prev_coord.assign(pcoord, pcoord + (std::size_t)nnode*NDIMS);
        prev_conn.assign(pconn, pconn + (std::size_t)nelem*NODES_PER_ELEM);
        prev_seg.assign(pseg, pseg + (std::size_t)nseg*NODES_PER_FACET);
        prev_segflag.assign(psegflag, psegflag + nseg);
        prev_regattr.assign(pregattr, pregattr + nelem);
        prev_size = *var.init_elem_size_n;
        const int nnode0 = nnode, nelem0 = nelem;
        adopt_initial_mmg_output(var, hint, mmg_out, nnode, nelem, nseg,
                                 pcoord, pconn, pseg, psegflag, pregattr);
#ifdef THREED
        regions_by_containing_element(prev_coord, prev_conn, prev_regattr, nelem, pcoord, pconn,
                                      pregattr);
#endif
        std::cout << "  [weakzone] pass " << pass << ": " << nfree << " free weak-zone element(s) "
                  << "(" << ntop << " on a top-surface node and " << nmixed << " on a region "
                  << "boundary stay), " << nrebuild << " element(s) MMG may rebuild at the first "
                  << "remesh's metric: " << nnode0 << " -> " << nnode << " nodes, " << nelem0
                  << " -> " << nelem << " elements.\n";
    }
}

// Refine the restored side walls of the initial mesh to ~mesh.resolution (remeshing_option 13):
// the incoming-material profile is sampled at wall-node spacing, so a coarse wall would freeze a
// kinked profile for the whole run. MMG rebuilds the wall band (each element touching a wall node
// + 2 rings, unless already at resolution) at a metric clamped at the wall nodes;
// var.init_elem_size_n keeps the UNCLAMPED size, so off-wall elements coarsen back at later
// remeshes. Regions: in 2D under mattype_option 0, their only reader, the band's region boundaries
// are required interface edges and each new element keeps its region (element references);
// otherwise a new element takes the nearest old centroid's.
static void refine_initial_side_walls(const Param &param, Variables &var,
        int &nnode, int &nelem, int &nseg,
        double *&pcoord, int *&pconn, int *&pseg, int *&psegflag, double *&pregattr)
{
    const double res = param.mesh.resolution;

    // Per-node actual element size.
    std::vector<double> size_n, esize;
    measured_size_n(nnode, nelem, pcoord, pconn, esize, size_n);

    // Clamp at every restored wall, whatever its BC: a static wall keeps a full-resolution profile
    // should the flow turn inward.
    auto on_wall = [&](int n) {
        for (int t = 0; t < NSIDEWALL; ++t) {
            const double plane = side_wall_plane(param.mesh, t);
            const double tol = 1e-9 * (SIDEWALL_AXIS[t] == 0 ? param.mesh.xlength : param.mesh.ylength);
            if (std::abs(pcoord[n*NDIMS + SIDEWALL_AXIS[t]] - plane) < tol) return true;
        }
        return false;
    };
    std::vector<double> metric(size_n);
    bool needs_refining = false;
    for (int n = 0; n < nnode; ++n) {
        if (on_wall(n) && metric[n] > res) {
            metric[n] = res;
            needs_refining = true;
        }
    }
    if (!needs_refining) return;

    // The cfg-intended size: the one use_mmg_init or refine_initial_weak_zone set, else the
    // measured one. A copy: init_elem_size_n is reassigned below.
    const double_vec hint = ((int)var.init_elem_size_n->size() == nnode)
                            ? *var.init_elem_size_n : size_n;

    // Free the wall band; the rest of the mesh is required.
    std::vector<char> free_elem(nelem, 0), free_node(nnode, 0);
    for (int n = 0; n < nnode; ++n)
        if (on_wall(n))
            free_node[n] = 1;                        // ring 0 grows from the wall nodes
    for (int r = 0; r <= 2; ++r) {
        std::vector<char> next(nnode, 0);
        for (int e = 0; e < nelem; ++e) {
            if (free_elem[e]) continue;
            if (esize[e] <= res) continue;           // already at/below resolution
            const int *c = pconn + e*NODES_PER_ELEM;
            bool touch = false;
            for (int k = 0; k < NODES_PER_ELEM && !touch; ++k) touch = free_node[c[k]];
            if (!touch) continue;
            free_elem[e] = 1;
            for (int k = 0; k < NODES_PER_ELEM; ++k) next[c[k]] = 1;
        }
        for (int n = 0; n < nnode; ++n)
            if (next[n]) free_node[n] = 1;
    }
    std::vector<char> req_elem(nelem), req_node(nnode);
    for (int e = 0; e < nelem; ++e) req_elem[e] = !free_elem[e];
    for (int n = 0; n < nnode; ++n) req_node[n] = !free_node[n];
    // mmg_adapt's contract: a required element's nodes are required.
    for (int e = 0; e < nelem; ++e)
        if (req_elem[e])
            for (int k = 0; k < NODES_PER_ELEM; ++k) req_node[pconn[e*NODES_PER_ELEM + k]] = 1;

#ifndef THREED
    // An edge of two free elements of different regattr (truncated, as markerset reads it) is a
    // required interface edge; a frozen element's edges are required anyway.
    const bool keep_regions = (param.ic.mattype_option == 0);
    std::vector<int> eref, iface;
    if (keep_regions) {
        eref.assign(nelem, 0);
        if (pregattr)
            for (int e = 0; e < nelem; ++e) eref[e] = (int)pregattr[e];
        std::map<std::pair<int,int>, int> first_ref;   // edge -> its first free element's reference
        for (int e = 0; e < nelem; ++e) {
            if (!free_elem[e]) continue;
            const int *c = pconn + e*NODES_PER_ELEM;
            for (int i = 0; i < NODES_PER_ELEM; ++i)
                for (int j = i + 1; j < NODES_PER_ELEM; ++j) {
                    const std::pair<int,int> key(std::min(c[i], c[j]), std::max(c[i], c[j]));
                    auto ins = first_ref.emplace(key, eref[e]);
                    if (!ins.second && ins.first->second != eref[e]) {
                        iface.push_back(key.first);
                        iface.push_back(key.second);
                    }
                }
        }
    }
    // MMG cannot split a required edge, so bisect each to the size MMG grades up to from the wall,
    // resolution + ln(1.3) x the wall distance (MMG5_gradsiz_iso; mmg_adapt keeps hgrad 1.3); a
    // long one leaves poor elements beside it. The targets are scratch and never pulled down (a
    // midpoint gets its own), so the wall size does not spread along an edge.
    auto graded = [&](const double *x) {
        double d = std::numeric_limits<double>::max();
        for (int t = 0; t < NSIDEWALL; ++t)
            d = std::min(d, std::abs(x[SIDEWALL_AXIS[t]] - side_wall_plane(param.mesh, t)));
        return res + 0.26236426446 * d;   // MMG5_HGRAD = ln(1.3), mmgcommon_private.h
    };
    std::vector<double> in_coord(pcoord, pcoord + (std::size_t)nnode*NDIMS);
    std::vector<int> in_conn(pconn, pconn + (std::size_t)nelem*NODES_PER_ELEM);
    int in_nnode = nnode, in_nelem = nelem;
    std::vector<double> tgt(nnode);
    for (int n = 0; n < nnode; ++n) tgt[n] = graded(&in_coord[(std::size_t)n*NDIMS]);
    for (int ne0 = -1; ne0 != in_nelem; ) {
        ne0 = in_nelem;
        const int nn0 = in_nnode;
        presubdivide_interface_edges(in_coord, in_conn, tgt, req_node, req_elem, iface,
                                     in_nnode, in_nelem, &eref, false);
        // `metric` (MMG's sizes) takes each midpoint's h, as at a remesh; `tgt` its graded size.
        metric.insert(metric.end(), tgt.begin() + nn0, tgt.end());
        for (int n = nn0; n < in_nnode; ++n) tgt[n] = graded(&in_coord[(std::size_t)n*NDIMS]);
    }
    if (!iface.empty())
        std::cout << "      [require] " << iface.size() / 2 << " material-interface edge(s) of the "
                     "initial wall band handed to MMG as required.\n";
    MMGInput  mmg_in = { in_nnode, in_nelem, nseg, in_coord.data(), in_conn.data(), pseg, psegflag,
                         metric.data(), req_node.data(), req_elem.data(), true, nullptr,
                         iface.empty() ? nullptr : iface.data(), (int)iface.size() / 2, nullptr,
                         keep_regions ? eref.data() : nullptr };
#else
    MMGInput  mmg_in = { nnode, nelem, nseg, pcoord, pconn, pseg, psegflag,
                         metric.data(), req_node.data(), req_elem.data(), true, nullptr };
#endif
    MMGOutput mmg_out;
    // Gated as at a remesh: a bad output element would make the first quality check remesh. The
    // gate scans the whole output, so it also repairs a bad element left from the initial mesh.
    mmg_adapt_quality_gated(param, mmg_in, mmg_out, req_node, req_elem);

    adopt_initial_mmg_output(var, hint, mmg_out, nnode, nelem, nseg,
                             pcoord, pconn, pseg, psegflag, pregattr);

    std::cout << "  Refined restored side walls of the initial mesh to ~resolution: "
              << nnode << " nodes, " << nelem << " elements\n";
}
#endif // USEMMG


#ifdef THREED
void tetrahedralize_polyhedron
(double max_ratio, double min_dihedral_angle, double max_volume,
 int vertex_per_polygon, int meshing_verbosity, int optlevel,
 int npoints, int nsegments,
 const double *points, const int *segments, const int *segflags,
 const tetgenio::facet *facets,
 const int nregions, const double *regionattributes,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags, double **outregattr)
{
    //
    // Setting Tetgen options.
    //
    char options[255];
    double max_dihedral_angle = 180 - 3 * min_dihedral_angle;

    std::string verbosity, vol, quality;
    set_verbosity_str(verbosity, meshing_verbosity);
    set_volume_str(vol, max_volume);
    set_3d_quality_str(quality, max_ratio, min_dihedral_angle, max_dihedral_angle);

    if( nregions > 0 )
        std::sprintf(options, "%s%s%spzs%dA", verbosity.c_str(), quality.c_str(), vol.c_str(), optlevel);
    else
        std::sprintf(options, "%s%s%spzs%d", verbosity.c_str(), quality.c_str(), vol.c_str(), optlevel);

    if( meshing_verbosity >= 0 )
        std::cout << "The meshing option is: " << options << '\n';

    //
    // Setting input arrays to tetgen
    //

    // NOTE: all tetgenio pointers will need to be
    // reset to NULL to prevent double-free error
    tetgenio in;
    in.pointlist = const_cast<double*>(points);
    in.numberofpoints = npoints;

    tetgenio::polygon *polys;
    tetgenio::facet *fl;
    if (facets == NULL) {
        polys = new tetgenio::polygon[nsegments];
        for (int i=0; i<nsegments; ++i) {
            polys[i].vertexlist = const_cast<int*>(&segments[i*vertex_per_polygon]);
            polys[i].numberofvertices = vertex_per_polygon;
        }

        fl = new tetgenio::facet[nsegments];
        for (int i=0; i<nsegments; ++i) {
            fl[i].polygonlist = &polys[i];
            fl[i].numberofpolygons = 1;
            fl[i].holelist = NULL;
            fl[i].numberofholes = 0;
        }
    }
    else {
        fl = const_cast<tetgenio::facet*>(facets);
    }

    in.facetlist = fl;
    in.facetmarkerlist = const_cast<int*>(segflags);
    in.numberoffacets = nsegments;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.numberofregions = nregions;
    if( nregions > 0 )
        in.regionlist = const_cast<double*>(regionattributes);
    else
        in.regionlist = NULL;

    tetgenio out;
    /*******************************/
    tetrahedralize(options, &in, &out, NULL, NULL);
    /*******************************/

    // the destructor of tetgenio will free any non-NULL pointer
    // set in.pointers to NULL to prevent double-free
    in.pointlist = NULL;
    in.facetmarkerlist = NULL;
    in.facetlist = NULL;
    in.regionlist = NULL;
    if (facets == NULL) {
        delete [] polys;
        delete [] fl;
    }

    *noutpoints = out.numberofpoints;
    *outpoints = out.pointlist;
    out.pointlist = NULL;

    *ntriangles = out.numberoftetrahedra;
    *triangles = out.tetrahedronlist;
    out.tetrahedronlist = NULL;

    *noutsegments = out.numberoftrifaces;
    *outsegments = out.trifacelist;
    *outsegflags = out.trifacemarkerlist;
    *outregattr = out.tetrahedronattributelist;
    out.trifacelist = NULL;
    out.trifacemarkerlist = NULL;
    out.tetrahedronattributelist = NULL;

}
#endif


// Shared tail of every Triangle/TetGen-built initial mesh (points_to_mesh, the 3D .poly path):
// refine_initial_weak_zone, refine_initial_side_walls (option 13), then var copies the buffers and
// this function frees them. An imported .exo mesh (meshing_option 95) skips it.
void adopt_initial_mesh(const Param &param, Variables &var,
                        double *pcoord, int *pconnectivity, int *psegment, int *psegflag,
                        double *pregattr)
{
#ifdef USEMMG
    refine_initial_weak_zone(param, var, var.nnode, var.nelem, var.nseg,
                             pcoord, pconnectivity, psegment, psegflag, pregattr);
    if (param.mesh.remeshing_option == 13 && param.mesh.meshing_elem_shape == 0)
        refine_initial_side_walls(param, var, var.nnode, var.nelem, var.nseg,
                                  pcoord, pconnectivity, psegment, psegflag, pregattr);
#endif

    var.coord = new array_t(pcoord, var.nnode);
    var.connectivity = new conn_t(pconnectivity, var.nelem);
    var.segment = new segment_t(psegment, var.nseg);
    var.segflag = new segflag_t(psegflag, var.nseg);
    var.regattr = new regattr_t(pregattr, var.nelem);

    delete [] pcoord;
    delete [] pconnectivity;
    delete [] psegment;
    delete [] psegflag;
    delete [] pregattr;
}


void points_to_mesh(const Param &param, Variables &var,
                    int npoints, const double *points,
                    int n_init_segments, const int *init_segments, const int *init_segflags,
                    int nregions, const double *regattr,
                    double max_elem_size, int vertex_per_polygon)
{
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

    double coarse_max_elem_size = max_elem_size;
    std::vector<double> coarse_regattr_vec;
    const double *coarse_regattr_ptr = regattr;

#ifdef USEMMG
    if (param.mesh.use_mmg_init) {
        // Step 1: scale volume for coarse TetGen mesh
        if (max_elem_size > 0) {
            coarse_max_elem_size = max_elem_size *
                std::pow(param.mesh.mmg_init_coarsening_factor, NDIMS);
        } else {
            // per-region: encode region index as temporary mattype so TetGen -A
            // stamps coarse elements with their zone index (0, 1, ...).
            // This allows compute_init_metric() to use connectivity-based projection,
            // correctly preserving per-zone resolution (e.g. meshing_option=2).
            const int stride = NDIMS + 2;
            coarse_regattr_vec.assign(regattr, regattr + nregions * stride);
            // Find the fine zone: region with the smallest volume constraint.
            // Keep it at target density so TetGen places elements near the fine zone
            // resolution directly, avoiding a large MMG transition zone that would
            // cause over-refinement at the first remesh.
            double min_vol = std::numeric_limits<double>::max();
            for (int r = 0; r < nregions; ++r) {
                double vol = coarse_regattr_vec[r * stride + NDIMS + 1];
                if (vol > 0 && vol < min_vol) min_vol = vol;
            }
            for (int r = 0; r < nregions; ++r) {
                coarse_regattr_vec[r * stride + NDIMS] = (double)r;  // temporary: region index as mattype
                double &vol = coarse_regattr_vec[r * stride + NDIMS + 1];
                // Only coarsen non-fine zones; fine zone stays at target density.
                if (vol > 0 && vol > min_vol) vol *= std::pow(param.mesh.mmg_init_coarsening_factor, NDIMS);
            }
            coarse_regattr_ptr = coarse_regattr_vec.data();
        }

        // For the coarse TetGen mesh, hardcode optlevel=0 to save time
        Mesh coarse_mesh = param.mesh;
        coarse_mesh.tetgen_optlevel = 0;

        points_to_new_mesh(coarse_mesh, npoints, points,
                           n_init_segments, init_segments, init_segflags,
                           nregions, coarse_regattr_ptr,
                           coarse_max_elem_size, vertex_per_polygon,
                           var.nnode, var.nelem, var.nseg,
                           pcoord, pconnectivity, psegment, psegflag, pregattr);

        // Step 2: MMG refinement — also compute fine-mesh init size hint
#ifdef THREED
        mmg_refine_init_mesh_3d(param.mesh, max_elem_size, nregions, regattr,
                                var.nnode, var.nelem, var.nseg,
                                pcoord, pconnectivity, psegment, psegflag, pregattr,
                                var.init_elem_size_n);
#else
        mmg_refine_init_mesh_2d(param.mesh, param.ic.mattype_option == 0, npoints, points,
                                max_elem_size, nregions, regattr,
                                var.nnode, var.nelem, var.nseg,
                                pcoord, pconnectivity, psegment, psegflag, pregattr,
                                var.init_elem_size_n);
#endif
    } else {
        // Plain triangle/tetgen path — init_elem_size_n left empty;
        // initialize_elem_size_n() will compute it from the mesh at step 0, unless
        // adopt_initial_mesh's refinements fill it.
        points_to_new_mesh(param.mesh, npoints, points,
                           n_init_segments, init_segments, init_segflags,
                           nregions, regattr,
                           max_elem_size, vertex_per_polygon,
                           var.nnode, var.nelem, var.nseg,
                           pcoord, pconnectivity, psegment, psegflag, pregattr);
    }
#else // !USEMMG
    points_to_new_mesh(param.mesh, npoints, points,
                       n_init_segments, init_segments, init_segflags,
                       nregions, regattr,
                       max_elem_size, vertex_per_polygon,
                       var.nnode, var.nelem, var.nseg,
                       pcoord, pconnectivity, psegment, psegflag, pregattr);
#endif

    adopt_initial_mesh(param, var, pcoord, pconnectivity, psegment, psegflag, pregattr);
}

void new_mesh_regular(const Param& param, Variables& var)
{
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag, *pcell;

    create_quadrilateral_cells(var, pcell);
    var.cell = new regular_t(pcell, var.ncell);

    create_rect_node(param, var, pcoord);
    create_elem_from_cell(var, pconnectivity);

    var.coord = new array_t(pcoord, var.nnode);
    var.connectivity = new conn_t(pconnectivity, var.nelem);

    create_regular_segments(var, psegment, psegflag);
    var.segment = new segment_t(psegment, var.nseg);
    var.segflag = new segflag_t(psegflag, var.nseg);

    create_regular_regattr(var, pregattr);
    var.regattr = new regattr_t(pregattr, var.nelem);

    delete [] pcell;
    delete [] pcoord;
    delete [] pconnectivity;
    delete [] psegment;
    delete [] psegflag;
    delete [] pregattr;
}


void new_mesh_uniform_resolution(const Param& param, Variables& var)
{
    int npoints = 4 * (NDIMS - 1); // 2D:4;  3D:8
    double *points = new double[npoints*NDIMS];

    int n_init_segments = 2 * NDIMS; // 2D:4;  3D:6
    int n_segment_nodes = 2 * (NDIMS - 1); // 2D:2; 3D:4
    int *init_segments = new int[n_init_segments*n_segment_nodes];
    int *init_segflags = new int[n_init_segments];

    const int attr_ndata = NDIMS+2;
    const int nregions = (param.ic.mattype_option == 0) ? param.mat.nmat : 1;
    double *regattr = new double[nregions*attr_ndata]; // each region has (NDIMS+2) data fields: x, (y,) z, region marker (mat type), and volume.

    double elem_size;  // size of a typical element
    int vertex_per_polygon = 4;

#ifndef THREED
    {
	/* Define 4 corner points of the rectangle, with this order:
         *            BOUNDZ1
         *          0 ------- 3
         *  BOUNDX0 |         | BOUNDX1
         *          1 ------- 2
         *            BOUNDZ0
         */
        // corner 0
        points[0] = 0;
	points[1] = 0;
        // corner 1
	points[2] = 0;
	points[3] = -param.mesh.zlength;
        // corner 2
	points[4] = param.mesh.xlength;
	points[5] = -param.mesh.zlength;
        // corner 3
	points[6] = param.mesh.xlength;
	points[7] = 0;

	for (int i=0; i<n_init_segments; ++i) {
            // 0th node of the i-th segment
	    init_segments[2*i] = i;
            // 1st node of the i-th segment
	    init_segments[2*i+1] = i+1;
	}
	// the 1st node of the last segment is connected back to the 0th node
	init_segments[2*n_init_segments-1] = 0;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDZ0;
        init_segflags[2] = BOUNDX1;
        init_segflags[3] = BOUNDZ1;

        elem_size = 1.5 * param.mesh.resolution * param.mesh.resolution;

        // Need to modify when nregions > 1.
        for (int i = 0; i < nregions; i++) {
            regattr[i * attr_ndata] = 0.5*param.mesh.xlength;
            regattr[i * attr_ndata + 1] = -0.5*param.mesh.zlength;
            regattr[i * attr_ndata + 2] = 0;
            regattr[i * attr_ndata + 3] = -1;
        }
    }
#else
    {
	/* Define 8 corner points of the box, with this order:
         *         4 ------- 7
         *        /         /|
         *       /         / 6
         *      0 ------- 3 /
         *      |         |/
         *      1 ------- 2
         *
         * Cut-out diagram with boundary flag:
         *             4 ------- 7
         *             | BOUNDZ1 |
         *   4 ------- 0 ------- 3 ------- 7 ------- 4
         *   | BOUNDX0 | BOUNDY0 | BOUNDX1 | BOUNDY1 |
         *   5 ------- 1 ------- 2 ------- 6 ------- 5
         *             | BOUNDZ0 |
         *             5 ------- 6
         */

        // corner 0
        points[0] = 0;
	points[1] = 0;
	points[2] = 0;
        // corner 1
	points[3] = 0;
	points[4] = 0;
	points[5] = -param.mesh.zlength;
        // corner 2
	points[6] = param.mesh.xlength;
	points[7] = 0;
	points[8] = -param.mesh.zlength;
        // corner 3
	points[9] = param.mesh.xlength;
	points[10] = 0;
	points[11] = 0;
        // corner 4
	points[12] = 0;
	points[13] = param.mesh.ylength;
	points[14] = 0;
        // corner 5
	points[15] = 0;
	points[16] = param.mesh.ylength;
	points[17] = -param.mesh.zlength;
        // corner 6
	points[18] = param.mesh.xlength;
	points[19] = param.mesh.ylength;
	points[20] = -param.mesh.zlength;
        // corner 7
	points[21] = param.mesh.xlength;
	points[22] = param.mesh.ylength;
	points[23] = 0;

        // BOUNDX0
        init_segments[0] = 0;
        init_segments[1] = 1;
        init_segments[2] = 5;
        init_segments[3] = 4;
        // BOUNDY0
        init_segments[4] = 0;
        init_segments[5] = 3;
        init_segments[6] = 2;
        init_segments[7] = 1;
        // BOUNDZ0
        init_segments[8] = 1;
        init_segments[9] = 2;
        init_segments[10] = 6;
        init_segments[11] = 5;
        // BOUNDX1
        init_segments[12] = 3;
        init_segments[13] = 7;
        init_segments[14] = 6;
        init_segments[15] = 2;
        // BOUNDY1
        init_segments[16] = 7;
        init_segments[17] = 4;
        init_segments[18] = 5;
        init_segments[19] = 6;
        // BOUNDZ1
        init_segments[20] = 0;
        init_segments[21] = 4;
        init_segments[22] = 7;
        init_segments[23] = 3;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDY0;
        init_segflags[2] = BOUNDZ0;
        init_segflags[3] = BOUNDX1;
        init_segflags[4] = BOUNDY1;
        init_segflags[5] = BOUNDZ1;

        elem_size = 0.7 * param.mesh.resolution
            * param.mesh.resolution * param.mesh.resolution;

        // Need to modify when nregions > 1. Using .poly file might be better.
        for (int i = 0; i < nregions; i++) {
            regattr[i * attr_ndata] = 0.5*param.mesh.xlength;
            regattr[i * attr_ndata + 1] = 0.5*param.mesh.ylength;
            regattr[i * attr_ndata + 2] = -0.5*param.mesh.zlength;
            regattr[i * attr_ndata + 3] = 0;
            regattr[i * attr_ndata + 4] = -1;
        }
    }
#endif

    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags, nregions, regattr,
                   elem_size, vertex_per_polygon);

    delete [] points;
    delete [] init_segments;
    delete [] init_segflags;
    delete [] regattr;
}


void new_mesh_refined_zone(const Param& param, Variables& var)
{
    const Mesh& m = param.mesh;

    // bounding box
    const double Lx = m.xlength;
#ifdef THREED
    const double Ly = m.ylength;
#endif
    const double Lz = m.zlength;

    // typical distance between nodes in the refined zone
    const double d = m.resolution;

    // adjust the bounds of the refined zone so that nodes are not on the boundary
    const double x0 = std::max(m.refined_zonex.first, d / Lx);
    const double x1 = std::min(m.refined_zonex.second, 1 - d / Lx);
#ifdef THREED
    const double y0 = std::max(m.refined_zoney.first, d / Ly);
    const double y1 = std::min(m.refined_zoney.second, 1 - d / Ly);
#endif
    const double z0 = std::max(m.refined_zonez.first, d / Lz);
    const double z1 = std::min(m.refined_zonez.second, 1 - d / Lz);

    const int npoints = 2 * 4 * (NDIMS - 1);
    double *points = new double[npoints*NDIMS];

    const int nsegments = 2 * (2 * NDIMS); // 2 x 2D:4;  3D:6
    const int nnodes = 2 * (NDIMS - 1);    // 2D:2; 3D:4
    const int vertex_per_polygon = 4;
    int *segments = new int[nsegments*nnodes];
    int *segflags = new int[nsegments];

#ifndef THREED
    {
      /* Define 4 corner points of the rectangle, with this order:
       *            BOUNDZ1
       *          0 ------- 3
       *  BOUNDX0 |         | BOUNDX1
       *          1 ------- 2
       *            BOUNDZ0
       */

      // outer box            // inner box
      // corner 0             // corner 0
      points[0] = 0.0;        points[8+0] = +x0*Lx;
      points[1] = 0.0;        points[8+1] = -z0*Lz;
      // corner 1             // corner 1
      points[2] = 0.0;        points[8+2] = +x0*Lx;
      points[3] = -Lz;        points[8+3] = -z1*Lz;
      // corner 2             // corner 2
      points[4] = +Lx;        points[8+4] = +x1*Lx;
      points[5] = -Lz;        points[8+5] = -z1*Lz;
      // corner 3             // corner 3
      points[6] = +Lx;        points[8+6] = +x1*Lx;
      points[7] = 0.0;        points[8+7] = -z0*Lz;

      // outer segments       // inner segments
      // BOUNDX0              // BOUNDX0
      segments[0] = 0;        segments[8+0] = 4+0;
      segments[1] = 1;        segments[8+1] = 4+1;
      // BOUNDZ0              // BOUNDZ0
      segments[2] = 1;        segments[8+2] = 4+1;
      segments[3] = 2;        segments[8+3] = 4+2;
      // BOUNDX1              // BOUNDX1
      segments[4] = 2;        segments[8+4] = 4+2;
      segments[5] = 3;        segments[8+5] = 4+3;
      // BOUNDZ1              // BOUNDZ1
      segments[6] = 3;        segments[8+6] = 4+3;
      segments[7] = 0;        segments[8+7] = 4+0;

      // outer boundary       // inner boundary
      segflags[0] = BOUNDX0;  segflags[4+0] = 0;
      segflags[1] = BOUNDZ0;  segflags[4+1] = 0;
      segflags[2] = BOUNDX1;  segflags[4+2] = 0;
      segflags[3] = BOUNDZ1;  segflags[4+3] = 0;
    }
#else
    {
      /* Define 8 corner points of the box, with this order:
       *         4 ------- 7
       *        /         /|
       *       /         / 6
       *      0 ------- 3 /
       *      |         |/
       *      1 ------- 2
       *
       * Cut-out diagram with boundary flag:
       *             4 ------- 7
       *             | BOUNDZ1 |
       *   4 ------- 0 ------- 3 ------- 7 ------- 4
       *   | BOUNDX0 | BOUNDY0 | BOUNDX1 | BOUNDY1 |
       *   5 ------- 1 ------- 2 ------- 6 ------- 5
       *             | BOUNDZ0 |
       *             5 ------- 6
       */

      // outer box            // inner box
      // corner 0             // corner 0
      points[ 0] = 0.0;       points[24+ 0] = +x0*Lx;
      points[ 1] = 0.0;       points[24+ 1] = +y0*Ly;
      points[ 2] = 0.0;       points[24+ 2] = -z0*Lz;
      // corner 1             // corner 1
      points[ 3] = 0.0;       points[24+ 3] = +x0*Lx;
      points[ 4] = 0.0;       points[24+ 4] = +y0*Ly;
      points[ 5] = -Lz;       points[24+ 5] = -z1*Lz;
      // corner 2             // corner 2
      points[ 6] = +Lx;       points[24+ 6] = +x1*Lx;
      points[ 7] = 0.0;       points[24+ 7] = +y0*Ly;
      points[ 8] = -Lz;       points[24+ 8] = -z1*Lz;
      // corner 3             // corner 3
      points[ 9] = +Lx;       points[24+ 9] = +x1*Lx;
      points[10] = 0.0;       points[24+10] = +y0*Ly;
      points[11] = 0.0;       points[24+11] = -z0*Lz;
      // corner 4             // corner 4
      points[12] = 0.0;       points[24+12] = +x0*Lx;
      points[13] = +Ly;       points[24+13] = +y1*Ly;
      points[14] = 0.0;       points[24+14] = -z0*Lz;
      // corner 5             // corner 5
      points[15] = 0.0;       points[24+15] = +x0*Lx;
      points[16] = +Ly;       points[24+16] = +y1*Ly;
      points[17] = -Lz;       points[24+17] = -z1*Lz;
      // corner 6             // corner 6
      points[18] = +Lx;       points[24+18] = +x1*Lx;
      points[19] = +Ly;       points[24+19] = +y1*Ly;
      points[20] = -Lz;       points[24+20] = -z1*Lz;
      // corner 7             // corner 7
      points[21] = +Lx;       points[24+21] = +x1*Lx;
      points[22] = +Ly;       points[24+22] = +y1*Ly;
      points[23] = 0.0;       points[24+23] = -z0*Lz;

      // outer segments       // inner segments
      // BOUNDX0              // BOUNDX0
      segments[ 0] = 0;       segments[24+ 0] = 8+0;
      segments[ 1] = 1;       segments[24+ 1] = 8+1;
      segments[ 2] = 5;       segments[24+ 2] = 8+5;
      segments[ 3] = 4;       segments[24+ 3] = 8+4;
      // BOUNDY0              // BOUNDY0
      segments[ 4] = 0;       segments[24+ 4] = 8+0;
      segments[ 5] = 3;       segments[24+ 5] = 8+3;
      segments[ 6] = 2;       segments[24+ 6] = 8+2;
      segments[ 7] = 1;       segments[24+ 7] = 8+1;
      // BOUNDZ0              // BOUNDZ0
      segments[ 8] = 1;       segments[24+ 8] = 8+1;
      segments[ 9] = 2;       segments[24+ 9] = 8+2;
      segments[10] = 6;       segments[24+10] = 8+6;
      segments[11] = 5;       segments[24+11] = 8+5;
      // BOUNDX1              // BOUNDX1
      segments[12] = 3;       segments[24+12] = 8+3;
      segments[13] = 7;       segments[24+13] = 8+7;
      segments[14] = 6;       segments[24+14] = 8+6;
      segments[15] = 2;       segments[24+15] = 8+2;
      // BOUNDY1              // BOUNDY1
      segments[16] = 7;       segments[24+16] = 8+7;
      segments[17] = 4;       segments[24+17] = 8+4;
      segments[18] = 5;       segments[24+18] = 8+5;
      segments[19] = 6;       segments[24+19] = 8+6;
      // BOUNDZ1              // BOUNDZ1
      segments[20] = 0;       segments[24+20] = 8+0;
      segments[21] = 4;       segments[24+21] = 8+4;
      segments[22] = 7;       segments[24+22] = 8+7;
      segments[23] = 3;       segments[24+23] = 8+3;

      // outer boundary       // inner boundary
      segflags[0] = BOUNDX0;  segflags[6+0] = 0;
      segflags[1] = BOUNDY0;  segflags[6+1] = 0;
      segflags[2] = BOUNDZ0;  segflags[6+2] = 0;
      segflags[3] = BOUNDX1;  segflags[6+3] = 0;
      segflags[4] = BOUNDY1;  segflags[6+4] = 0;
      segflags[5] = BOUNDZ1;  segflags[6+5] = 0;
    }
#endif

    const int nregions = 2;
    const int ndata = NDIMS+2; // each region has (NDIMS+2) data fields: x, (y,) z, region marker (mat type), area/volume.
    double *regattr = new double[nregions*ndata];
    double *region0 = regattr;
    double *region1 = regattr + ndata;
#ifndef THREED
    const double area1 = 1.5 * (d*d);
#else
    const double area1 = 0.7 * (d*d*d);
#endif
    // Coarse-region max area in largest_size's relative-volume convention (that of smallest_size,
    // RefineFloors); area1 above keeps its own fixed factor.
    const double area0 = m.largest_size * sizefactor * std::pow(d, NDIMS);
    int pos = 0;
    region0[pos] = +d/2;  region1[pos++] = +x0*Lx+d/2; // x
#if THREED
    region0[pos] = +d/2;  region1[pos++] = +y0*Ly+d/2; // y
#endif
    region0[pos] = -d/2;  region1[pos++] = -z0*Lz-d/2; // z
    region0[pos] = 0.0;   region1[pos++] = 0.0;        // attribute
    region0[pos] = area0; region1[pos++] = area1;      // area/volume

    const double max_elem_size = 0;

    points_to_mesh(param, var, npoints, points,
                   nsegments, segments, segflags, nregions, regattr,
                   max_elem_size, vertex_per_polygon);

    delete [] points;
    delete [] segments;
    delete [] segflags;
    delete [] regattr;
}


void my_fgets(char *buffer, std::size_t size, std::FILE *fp,
              int &lineno, const std::string &filename)
{
    char *s;
    while (1) {
        ++ lineno;
        s = std::fgets(buffer, size, fp);
        if (! s) {
            std::cerr << "Error: reading line " << lineno
                      << " of '" << filename << "'\n";
            die(EXIT_IO_RW);
        }
        if (std::strlen(buffer) == size-1 && buffer[size-2] != '\n') {
            std::cerr << "Error: reading line " << lineno
                      << " of '" << filename << "', line is too long.\n";
            die(EXIT_IO_RW);
        }

        // check for blank lines and comments
        if (buffer[0] != '\n' && buffer[0] != '#') break;
    }
}


void new_mesh_from_polyfile(const Param& param, Variables& var)
{
    /* The format specifiction for the poly file can be found in:
     *   2D:  http://www.cs.cmu.edu/~quake/triangle.poly.html
     *   3D:  http://wias-berlin.de/software/tetgen/fformats.poly.html
     *
     * Note:
     * 1. the poly file in 3D has a more complicated format.
     * 2. The last row in the regional attributes is the regional constraint
     *    on the maximum triangle area or tetrahedra volume, which is the max
     *    size of an element.
     *    when meshing_option=90, the unit is meter^dim.
     *    when meshing_option=91, the unit is mesh.resolution^dim.
     */

#ifdef THREED
    const double std_elem_size = 0.7 * param.mesh.resolution
        * param.mesh.resolution * param.mesh.resolution;
#else
    const double std_elem_size = 1.5 * param.mesh.resolution * param.mesh.resolution;
#endif

    std::FILE *fp = std::fopen(param.mesh.poly_filename.c_str(), "r");
    if (! fp) {
        std::cerr << "Error: Cannot open poly_filename '" << param.mesh.poly_filename << "'\n";
        die(EXIT_IO_OPEN);
    }

    int lineno = 0;
    int n;
    char buffer[2550];

    // get header of node list
    int npoints;
    {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        int dim, nattr, nbdrym;
        n = std::sscanf(buffer, "%d %d %d %d", &npoints, &dim, &nattr, &nbdrym);
        if (n != 4) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }

        if (dim != NDIMS ||
            nattr != 0 ||
            nbdrym != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
    }

    // get node list
    double *points = new double[npoints * NDIMS];
    for (int i=0; i<npoints; i++) {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        int k;
        double *x = &points[i*NDIMS];
#ifdef THREED
        n = std::sscanf(buffer, "%d %lf %lf %lf", &k, x, x+1, x+2);
#else
        n = std::sscanf(buffer, "%d %lf %lf", &k, x, x+1);
#endif
        if (n != NDIMS+1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
        if (k != i) {
            std::cerr << "Error: node number is continuous from 0 at line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
    }

    // get header of segment (facet) list
    int n_init_segments;
    {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        int has_bdryflag;
        n = std::sscanf(buffer, "%d %d", &n_init_segments, &has_bdryflag);
        if (n != 2) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }

        if (has_bdryflag != 1) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
    }

    // get segment (facet) list
#ifdef THREED
    auto facets = new tetgenio::facet[n_init_segments];
    int *init_segflags = new int[n_init_segments];
    for (int i=0; i<n_init_segments; i++) {
        my_fgets(buffer, 2550, fp, lineno, param.mesh.poly_filename);

        auto &f = facets[i];
        int npolygons, nholes, bdryflag;
        n = std::sscanf(buffer, "%d %d %d", &npolygons, &nholes, &bdryflag);
        if (n != 3) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
        if (npolygons <= 0 ||
            nholes != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
        if (bdryflag == 0) goto flag_ok;
        for (int j=0; j<nbdrytypes; j++) {
            if (bdryflag == 1 << j) goto flag_ok;
        }
        std::cerr << "Error: bdry_flag has multiple bits set in line " << lineno
                  << " of '" << param.mesh.poly_filename << "'\n";
        die(EXIT_CONFIG_DATA);
    flag_ok:
        init_segflags[i] = bdryflag;

        f.polygonlist = new tetgenio::polygon[npolygons];
        f.numberofpolygons = npolygons;
        f.holelist = NULL;
        f.numberofholes = 0;

        for (int j=0; j<npolygons; j++) {
            my_fgets(buffer, 4096, fp, lineno, param.mesh.poly_filename);

            std::istringstream inbuf(std::string(buffer, 4096));
            int nvertex;
            inbuf >> nvertex;
            if (nvertex < NODES_PER_FACET || nvertex > 9999) {
                std::cerr << "Error: unsupported number of polygon points in line " << lineno
                          << " of '" << param.mesh.poly_filename << "'\n";
                die(EXIT_CONFIG_DATA);
            }

            f.polygonlist[j].vertexlist = new int[nvertex];
            f.polygonlist[j].numberofvertices = nvertex;
            for (int k=0; k<nvertex; k++) {
                inbuf >> f.polygonlist[j].vertexlist[k];
                if (f.polygonlist[j].vertexlist[k] < 0 ||
                    f.polygonlist[j].vertexlist[k] >= npoints) {
                    std::cerr << "Error: segment contains out-of-range node # [0-" << npoints
                              <<"] in line " << lineno << " of '"
                              << param.mesh.poly_filename << "'\n";
                    die(EXIT_CONFIG_DATA);
                }
            }
        }
    }
#else
    int *init_segments = new int[n_init_segments * NODES_PER_FACET];
    int *init_segflags = new int[n_init_segments];
    for (int i=0; i<n_init_segments; i++) {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int *x = &init_segments[i*NODES_PER_FACET];
        int junk, bdryflag;
        n = std::sscanf(buffer, "%d %d %d %d", &junk, x, x+1, &bdryflag);
        if (n != NODES_PER_FACET+2) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
        if (bdryflag == 0) goto flag_ok;
        for (int j=0; j<nbdrytypes; j++) {
            if (bdryflag == 1 << j) goto flag_ok;
        }
        std::cerr << "Error: bdry_flag has multiple bits set in line " << lineno
                  << " of '" << param.mesh.poly_filename << "'\n";
        die(EXIT_CONFIG_DATA);
    flag_ok:
        init_segflags[i] = bdryflag;
    }
    for (int i=0; i<n_init_segments; i++) {
        int *x = &init_segments[i*NODES_PER_FACET];
        for (int j=0; j<NODES_PER_FACET; j++) {
            if (x[j] < 0 || x[j] >= npoints) {
                std::cerr << "Error: segment contains out-of-range node # [0-" << npoints
                          <<"] in line " << lineno << " of '"
                          << param.mesh.poly_filename << "'\n";
                die(EXIT_CONFIG_DATA);
            }
        }
    }
#endif

    // get header of hole list
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int nholes;
        n = std::sscanf(buffer, "%d", &nholes);
        if (n != 1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }

        if (nholes != 0) {
            std::cerr << "Error: unsupported value in line " << lineno
                      << " of '" << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
    }

    // get header of region list
    int nregions;
    {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        n = std::sscanf(buffer, "%d", &nregions);
        if (n != 1) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
        if (nregions <= 0) {
            std::cerr << "Error: nregions <= 0, at line " << lineno << " of '"
                      << param.mesh.poly_filename << "'\n";
            die(EXIT_CONFIG_DATA);
        }
    }

    // get region list
    double *regattr = new double[nregions * (NDIMS+2)]; // each region has these data fields: x, (y,) z, region marker (mattype), and volume.
    bool has_max_size = false;
    for (int i=0; i<nregions; i++) {
        my_fgets(buffer, 255, fp, lineno, param.mesh.poly_filename);

        int junk;
        double *x = &regattr[i*(NDIMS+2)];
#ifdef THREED
        n = std::sscanf(buffer, "%d %lf %lf %lf %lf %lf", &junk, x, x+1, x+2, x+3, x+4);
#else
        n = std::sscanf(buffer, "%d %lf %lf %lf %lf", &junk, x, x+1, x+2, x+3);
#endif
        if (n != NDIMS+3) {
            std::cerr << "Error: parsing line " << lineno << " of '"
                      << param.mesh.poly_filename << "'. "<<NDIMS+3<<" values should be given but only "<<n<<" found.\n";
            die(EXIT_CONFIG_DATA);
        }

        if ( x[NDIMS] < 0 || x[NDIMS] >= param.mat.nmat ) {
            std::cerr << "Error: "<<NDIMS+2<<"-th value in line "<<lineno<<" should be >=0 and < "<<param.mat.nmat<<" (=mat.num_materials) but is "<<x[NDIMS]<<"\n";
            die(EXIT_CONFIG_DATA, "Note that this parameter is directly used as the index of mat. prop. arrays.");
        }

        if ( x[NDIMS+1] > 0 ) {
            has_max_size = true; // max area is set for this region

            // normalize the element size w.r.t mesh.resolution
            if (param.mesh.meshing_option == 91) {
                x[NDIMS+1] *= std_elem_size;
            }
        }
    }
    double max_elem_size = std_elem_size;
    if ( has_max_size ) max_elem_size = 0; // special value, see set_volume_str() above.

#ifdef THREED
    // shortcutting points_to_mesh() for polygon boundary; the tail is shared (adopt_initial_mesh)
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;

#ifdef USEMMG
    if (param.mesh.use_mmg_init) {
        // Step 1: scale volume for coarse TetGen mesh
        double coarse_max_elem_size = max_elem_size;
        std::vector<double> coarse_regattr_vec;
        const double *coarse_regattr_ptr = regattr;

        if (max_elem_size > 0) {
            coarse_max_elem_size = max_elem_size *
                std::pow(param.mesh.mmg_init_coarsening_factor, NDIMS);
        } else {
            // per-region: encode region index as temporary mattype so TetGen -A
            // stamps coarse elements with their zone index (0, 1, ...).
            // This allows compute_init_metric() to use connectivity-based projection,
            // correctly preserving per-zone resolution.
            const int stride = NDIMS + 2;
            coarse_regattr_vec.assign(regattr, regattr + nregions * stride);
            for (int r = 0; r < nregions; ++r) {
                coarse_regattr_vec[r * stride + NDIMS] = (double)r;  // temporary: region index as mattype
                double &vol = coarse_regattr_vec[r * stride + NDIMS + 1];
                if (vol > 0) vol *= std::pow(param.mesh.mmg_init_coarsening_factor, NDIMS);
            }
            coarse_regattr_ptr = coarse_regattr_vec.data();
        }

        tetrahedralize_polyhedron(param.mesh.max_ratio,
                                  param.mesh.min_tet_angle, coarse_max_elem_size,
                                  0,
                                  param.mesh.meshing_verbosity,
                                  0, // HARDCODED optlevel=0 for coarse mesh
                                  npoints, n_init_segments, points,
                                  NULL, init_segflags,
                                  facets,
                                  nregions, coarse_regattr_ptr,
                                  &var.nnode, &var.nelem, &var.nseg,
                                  &pcoord, &pconnectivity,
                                  &psegment, &psegflag, &pregattr);

        // Step 2: MMG refinement — also compute fine-mesh init size hint
        mmg_refine_init_mesh_3d(param.mesh, max_elem_size, nregions, regattr,
                                var.nnode, var.nelem, var.nseg,
                                pcoord, pconnectivity, psegment, psegflag, pregattr,
                                var.init_elem_size_n);
    } else {
        tetrahedralize_polyhedron(param.mesh.max_ratio,
                                  param.mesh.min_tet_angle, max_elem_size,
                                  0,
                                  param.mesh.meshing_verbosity,
                                  param.mesh.tetgen_optlevel,
                                  npoints, n_init_segments, points,
                                  NULL, init_segflags,
                                  facets,
                                  nregions, regattr,
                                  &var.nnode, &var.nelem, &var.nseg,
                                  &pcoord, &pconnectivity,
                                  &psegment, &psegflag, &pregattr);
    }
#else // !USEMMG
    tetrahedralize_polyhedron(param.mesh.max_ratio,
                              param.mesh.min_tet_angle, max_elem_size,
                              0,
                              param.mesh.meshing_verbosity,
                              param.mesh.tetgen_optlevel,
                              npoints, n_init_segments, points,
                              NULL, init_segflags,
                              facets,
                              nregions, regattr,
                              &var.nnode, &var.nelem, &var.nseg,
                              &pcoord, &pconnectivity,
                              &psegment, &psegflag, &pregattr);
#endif

    adopt_initial_mesh(param, var, pcoord, pconnectivity, psegment, psegflag, pregattr);
#else
    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags, nregions, regattr,
                   max_elem_size, NODES_PER_FACET);
#endif

    delete [] points;
#ifdef THREED
    for (int i=0; i<n_init_segments; i++) {
        auto f = facets[i];
        for (int j=0; j<f.numberofpolygons; j++) {
            delete [] f.polygonlist[j].vertexlist;
        }
        delete [] f.polygonlist;
    }
    delete [] facets;
#else
    delete [] init_segments;
#endif
    delete [] init_segflags;
    delete [] regattr;
}

#ifdef USEEXODUS
void new_mesh_from_exofile(const Param& param, Variables& var)
{

#ifndef THREED
    die(EXIT_UNSUPPORTED_DIM, "Importing an exofile currently works in 3D only.");
#endif

    // Open an .exo file.
    int CPU_word_size = 0;
    int IO_word_size = 0;
    float version = 0.0;
    int exoid = ex_open( param.mesh.exo_filename.c_str(), /* filename path */
                         EX_READ, /* access mode = READ */
                         &CPU_word_size, /* CPU word size */
                         &IO_word_size, /* IO word size */
                         &version); /* ExodusII library version */
    if (exoid < 0) {
        std::cerr << "Error: Cannot open exo_filename '" << param.mesh.exo_filename << "'\n";
        die(EXIT_IO_OPEN);
    }

    // Read database parameters.
    int error;
    char title[MAX_LINE_LENGTH+1];
    int num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;

    error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem,
                         &num_elem_blk, &num_node_sets, &num_side_sets);
    if( error != 0 ) {
        std::cerr << "Error: Unable to read database parameters from '" << param.mesh.exo_filename << std::endl;
        die(EXIT_IO_RW);
    }
    var.nnode = num_nodes;
    var.nelem = num_elem;
    std::cerr <<" Reading " << param.mesh.exo_filename <<"."<<std::endl;
    std::cerr <<" Numbers of nodes and elements are " << var.nnode <<" and "<< var.nelem <<"."<<std::endl;
    std::cerr <<" Number of element blocks is " << num_elem_blk <<"."<<std::endl;
    std::cerr <<" Number of node sets is " << num_node_sets <<"."<<std::endl;
    std::cerr <<" Number of side sets is " << num_side_sets <<"."<<std::endl;
    if( param.mat.nmat != num_elem_blk) {
        std::cerr <<"param.mat.nmat is not equal to # of element blocks in this exo file!"<<std::endl;
        std::cerr <<"Check if your material parameters are properly set!!"<<std::endl;
        die(EXIT_IO_RW);
    }

    // Assign node coordinates.
    float *x = (float *) calloc(var.nnode, sizeof(float));
    float *y = (float *) calloc(var.nnode, sizeof(float));
    float *z = (float *) calloc(var.nnode, sizeof(float));

    error = ex_get_coord (exoid, x, y, z);
    if( error != 0 ) {
        std::cerr << "Error: Unable to read coordinates from '" << param.mesh.exo_filename << "'\n";
        die(EXIT_IO_RW);
    }

    // assign node coordinates to var.coord.
    var.coord = new array_t(var.nnode);
    for(int i=0; i<var.nnode; i++) {
        (*var.coord)[i][0] = static_cast<double>(x[i]);
        (*var.coord)[i][1] = static_cast<double>(y[i]);
        (*var.coord)[i][2] = static_cast<double>(z[i]);
    }
    free(x);
    free(y);
    free(z);

    // Process element blocks.
    //
    // Note that blocks are the basic units for
    // organizing connectivity and material id.
    //
    int *ids = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_elem_in_block = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_nodes_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_edges_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_faces_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    int *num_attr = (int *) calloc(num_elem_blk, sizeof(int));
    char elem_type[MAX_STR_LENGTH+1];

    // - Read element block ids.
    error = ex_get_ids (exoid, EX_ELEM_BLOCK, ids);
    if( error != 0 ) {
        std::cerr << "Error: Unable to get element block ids." << std::endl;
        die(EXIT_IO_RW);
    }
    // - Read element block parameters.
    for (int i=0; i<num_elem_blk; i++) {
        error = ex_get_block (exoid, EX_ELEM_BLOCK, ids[i], elem_type,
                                &(num_elem_in_block[i]), &(num_nodes_per_elem[i]),
                                &(num_edges_per_elem[i]), &(num_faces_per_elem[i]), 
                                &(num_attr[i]));
        if( error != 0 ) {
            std::cerr << "Error: Unable to element block " << ids[i] ;
            std::cerr << " out of " << num_elem_blk << " blocks." << std::endl;
            die(EXIT_IO_RW);
        }
        if( NODES_PER_ELEM != num_nodes_per_elem[i] ) {
            std::cerr << "Error: Element has " << num_nodes_per_elem[i] << " nodes per element but should have "<<NODES_PER_ELEM<<" because element type should be uniformly tetrahedral."<< std::endl;
            die(EXIT_IO_RW);
        }
    }

    // Create the array of region attributes.
    // Element blocks must be numbered to be consistent
    // with the following simplistic mapping between block id and mat. id:
    //
    // block id: 1 2 3 4 ...
    // mat.  id: 0 1 2 3 ...
    //
    // Will need a more general mapping from block id to mat id.
    // One way is to provide a map from block names to mat id. TBD.
    //
    var.regattr = new regattr_t(var.nelem);
    int start = 0;
    for (int i=0; i<num_elem_blk; i++) {
        for(int j=0; j<num_elem_in_block[i]; j++)
            (*var.regattr)[start+j][0] = static_cast<double>(ids[i]-1);
        start += num_elem_in_block[i];
    }

    // - Read and append element connectivity
    int* connect[num_elem_blk];
    for (int i=0; i<num_elem_blk; i++) {
        connect[i] = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), sizeof(int));
        error = ex_get_conn (exoid, EX_ELEM_BLOCK, ids[i], connect[i], 0, 0);
        if( error != 0 ) {
            std::cerr << "Error: Unable to connectivity for element block " << ids[i] ;
            std::cerr << " out of " << num_elem_blk << " blocks." << std::endl;
            die(EXIT_IO_RW);
        }
    }

    { // To render 'int *conn' local to this block.
        var.connectivity = new conn_t(var.nelem);
        start = 0;
        for (int i=0; i<num_elem_blk; i++) {
            for(int j=0; j<num_elem_in_block[i]; j++) {
                const int elem_num = start + j;
                for(int k=0; k < NODES_PER_ELEM; k++ )
                    (*var.connectivity)[elem_num][k] = connect[i][ NODES_PER_ELEM*j + k ]-1;
            }
            // can the j and k loops be replaced with memcpy?
            // std::memcpy( &(conn[NODES_PER_ELEM*start]), &(connect[i]), num_elem_in_block[i]*NODES_PER_ELEM*sizeof(int) );
            start += num_elem_in_block[i];
        }
    } // code block to localize conn

    for (int i=0; i<num_elem_blk; i++) free(connect[i]);
    free (ids);
    free (num_elem_in_block);
    free (num_nodes_per_elem);
    free (num_edges_per_elem);
    free (num_faces_per_elem);
    free (num_attr);

    //
    // Read side sets (i.e., boundaries) one by one.
    //
    ids = (int *) calloc(num_side_sets, sizeof(int));
    int *num_sides_in_set = (int *) calloc(num_side_sets, sizeof(int));
    int *num_df_in_set = (int *) calloc(num_side_sets, sizeof(int));

    // - Read side set ids.
    error = ex_get_ids (exoid, EX_SIDE_SET, ids);
    if( error != 0 ) {
            std::cerr << "Error: Unable to get side set ids." << std::endl;
            die(EXIT_IO_RW);
    }
    var.nseg = 0;
    for(int i=0; i<num_side_sets; i++) {
        error = ex_get_set_param (exoid, EX_SIDE_SET, ids[i], &(num_sides_in_set[i]),
                                       &(num_df_in_set[i]));
        if( error != 0 ) {
            std::cerr << "Error: Unable to read "<<i<<"-th side set parameters." << std::endl;
            die(EXIT_IO_RW);
        }
        var.nseg += num_sides_in_set[i];
    }
    std::cerr <<" Numbers of segments are " << var.nseg <<"."<<std::endl;

    // list of elements in the side set
    int *elem_list[num_side_sets];
    // list of facets (sides) of a corresponding element in the side set
    int *side_list[num_side_sets];
    // list of the number of nodes of a facet. Uniformly 3 in DES.
    int *node_cnt_list[num_side_sets];
    // list of node IDs of the facet-composing nodes
    int *node_list[num_side_sets];
    // list of distribution factors
    float *dist_fact[num_side_sets];

    for(int i=0; i<num_side_sets; i++) {
        elem_list[i]     = (int *) calloc(num_sides_in_set[i], sizeof(int));
        side_list[i]     = (int *) calloc(num_sides_in_set[i], sizeof(int));
        node_cnt_list[i] = (int *) calloc(num_sides_in_set[i], sizeof(int));
        node_list[i]     = (int *) calloc(num_sides_in_set[i]*FACETS_PER_ELEM, sizeof(int));
        dist_fact[i]     = (float *) calloc(num_df_in_set[i], sizeof(float));

        // list of elements of which facet belongs to a side set
        // Note: The # of elements is same as # of sides!
        error = ex_get_set (exoid, EX_SIDE_SET, ids[i], elem_list[i], side_list[i]);
        if( error != 0 ) {
            std::cerr << "Error: Unable to read "<< i <<"-th side set." << std::endl;
            die(EXIT_IO_RW);
        }
        if (num_df_in_set != NULL) {
            error = ex_get_set_dist_fact (exoid, EX_SIDE_SET, ids[i], dist_fact[i]);
            if( error != 0 ) {
                std::cerr << "Error: Unable to read "<< i;
                std::cerr <<"-th side set's distribution factor list." << std::endl;
                die(EXIT_IO_RW);
            }
        }
    }

    // Assign memory to segments (boundary facets)
    var.segment = new segment_t(var.nseg, 0);

    // Assign memory to segment flags (which boundary a segments belong to).
    var.segflag = new segflag_t(var.nseg, 0);

    // Populate segment and segflag arrays
    // For sideset node ordering,
    // refer Table 4.2 on p. 30, "Exodus: A finite element data model" by Sjaardema et al.
    // Retrieved on 2019/09/28 from gsjaardema.github.io/seacas/exodusII-new.pdf.
    int_vec2D local_node_list{{1,2,4},{2,3,4},{1,4,3},{1,3,2}};
    start = 0;
    for (int i=0; i<num_side_sets; i++) {
        for (int j=0; j<num_sides_in_set[i]; j++) {
            const int elem_num = elem_list[i][j] - 1;
            const int side_num = side_list[i][j] - 1;
            for (int k=0; k<NODES_PER_FACET; k++) { // 3 nodes per triangular facet
                const int local_node_number = local_node_list[side_num][k] - 1;
                (*var.segment)[ start + j][k] = (*var.connectivity)[elem_num][local_node_number] - 1;
            }
            (*var.segflag)[start + j][0] = ids[i];
        }
        start += num_sides_in_set[i];
    }

    // free up memory.
    for (int i=0; i<num_side_sets; i++) {
        free (elem_list[i]);
        free (side_list[i]);
        free (node_cnt_list[i]);
        free (node_list[i]);
        free (dist_fact[i]);
    }
    free (num_sides_in_set);
    free (num_df_in_set);
    free (ids);

} // end of function new_mesh_from_exofile
#endif // end of USEEXODUS

} // end of anonymous namespace


void create_equilateral_elem(const Variables& var, int *&connectivity) {
    connectivity = new int[var.nelem*3];

    int istart_n1 = var.nx * (int(var.nz/2)+var.nz%2);

    int ind = 0;
    for (int j=0; j<var.nz-1;j++) {
        int jnc = j%2;
        int in0 = var.nx * int((j+1)/2);
        int in1 = istart_n1 + (var.nx+1)*int(j/2);
        int istarts[2] = {in0, in1};

        for (int k=0; k<2;k++) {
            int inc = (k+j)%2;
            int istart0 = istarts[inc];
            int istart1 = istarts[1-inc];

            if (inc==0) {
                for (int i=0; i<var.nx-1;i++) {
                    connectivity[ind*3] = istart0+i;
                    connectivity[ind*3+1+jnc] = istart1+i+1;
                    connectivity[ind*3+2-jnc] = istart0+i+1;
                    ind++;
                }
            } else {
                for (int i=0; i<var.nx;i++) {
                    connectivity[ind*3] = istart0+i;
                    connectivity[ind*3+1+jnc] = istart0+i+1;
                    connectivity[ind*3+2-jnc] = istart1+i;
                    ind++;
                }
            }
        }
    }
}

void create_equilateral_segments(const Variables& var, int *&segments, int *&segflags) {
    segments = new int[var.nseg*2];
    segflags = new int[var.nseg];
    int seg_idx = 0;
    int flag_idx = 0;

    // top
    for (int i = 0; i<var.nx-1; ++i) {
        segments[seg_idx*2] = i;
        segments[seg_idx*2+1] = i+1;
        seg_idx++;
        segflags[flag_idx] = 32;
        flag_idx++;
    }
    // bottom
    int nbot = var.nx - var.nz%2;
    int bot_node0;
    if (var.nz % 2 == 0) {
        // last odd row: tail of node array
        bot_node0 = var.nnode - nbot - 1;
    } else {
        // last even row: middle of node array
        bot_node0 = int(var.nz/2) * var.nx;
    }
    for (int i = 0; i < nbot; ++i) {
        segments[seg_idx*2]   = bot_node0 + i;
        segments[seg_idx*2+1] = bot_node0 + i + 1;
        seg_idx++;
        segflags[flag_idx] = 16;
        flag_idx++;
    }

    // left and right
    int istart_n2 = var.nx * (int(var.nz/2)+var.nz%2);
    for (int j=0; j<var.nz-1;j++) {
        int jnc = j%2;
        segments[seg_idx*2+jnc] = (int(j/2)+jnc)*var.nx;
        segments[seg_idx*2+1-jnc] = istart_n2 + (int(j/2))*(var.nx+1);
        seg_idx++;
        segflags[flag_idx] = 1;
        flag_idx++;

        segments[seg_idx*2+jnc] = (int(j/2)+jnc+1)*var.nx -1;
        segments[seg_idx*2+1-jnc] = istart_n2 + (int(j/2)+1)*(var.nx+1) -1;
        seg_idx++;
        segflags[flag_idx] = 2;
        flag_idx++;
    }
}

void points_to_new_mesh(const Mesh &mesh, int npoints, const double *points,
                        int n_init_segments, const int *init_segments, const int *init_segflags,
                        int n_regions, const double *regattr,
                        double max_elem_size, int vertex_per_polygon,
                        int &nnode, int &nelem, int &nseg, double *&pcoord,
                        int *&pconnectivity, int *&psegment, int *&psegflag, double *&pregattr)
{
#ifdef THREED

    tetrahedralize_polyhedron(mesh.max_ratio,
                              mesh.min_tet_angle, max_elem_size,
                              vertex_per_polygon,
                              mesh.meshing_verbosity,
                              mesh.tetgen_optlevel,
                              npoints, n_init_segments, points,
                              init_segments, init_segflags,
                              NULL,
                              n_regions, regattr,
                              &nnode, &nelem, &nseg,
                              &pcoord, &pconnectivity,
                              &psegment, &psegflag, &pregattr);

#else

    triangulate_polygon(mesh.min_angle, max_elem_size,
                        mesh.meshing_verbosity, mesh.max_steiner_factor,
                        npoints, n_init_segments, points,
                        init_segments, init_segflags,
                        n_regions, regattr,
                        &nnode, &nelem, &nseg,
                        &pcoord, &pconnectivity,
                        &psegment, &psegflag, &pregattr);

#endif

    if (nelem <= 0) {
#ifdef THREED
        std::cerr << "Error: tetrahedralization failed\n";
#else
        std::cerr << "Error: triangulation failed\n";
#endif
        die(EXIT_MESH_TETGEN);
    }

}


void points_to_new_surface(const Mesh &mesh, int npoints, const double *points,
                           int n_init_segments, const int *init_segments, const int *init_segflags,
                           int n_regions, const double *regattr,
                           double max_elem_size, int vertex_per_polygon,
                           int &nnode, int &nelem, int &nseg, double *&pcoord,
                           int *&pconnectivity, int *&psegment, int *&psegflag, double *&pregattr)
{
#ifdef THREED
    /* For triangulation of boundary surfaces in 3D */

    triangulate_polygon(mesh.min_angle, max_elem_size,
                        mesh.meshing_verbosity, mesh.max_steiner_factor,
                        npoints, n_init_segments, points,
                        init_segments, init_segflags,
                        n_regions, regattr,
                        &nnode, &nelem, &nseg,
                        &pcoord, &pconnectivity,
                        &psegment, &psegflag, &pregattr);

    if (nelem <= 0) {
        die(EXIT_MESH_TETGEN, "surface triangulation failed");
    }
#endif
}


void discard_internal_segments(int &nseg, segment_t &segment, segflag_t &segflag)
{
    int n = 0;
    while (n < nseg) {
        if (segflag[n][0] & BOUND_ANY) {
            // This is a boundary segment, go to next.
            n++;
        }
        else {
            // This is an internal boundary segment, replace it with the last segment.
            nseg--;
            segflag[n][0] = segflag[nseg][0];
            for (int i=0; i<NDIMS; i++) {
                segment[n][i] = segment[nseg][i];
            }
        }
    }

    segment.resize(nseg);
    segflag.resize(nseg);

}


void renumbering_mesh(const Param& param, array_t &coord, conn_t &connectivity,
                      segment_t &segment, regattr_t *regattr, double_vec *init_elem_size_n)
/* Note: the last argument is a pointer to regional attribute. If the pointer
 * is not NULL (when the initial mesh is created, before creating markers),
 * regattr is renumbered.
 */
{
    /* Renumbering nodes and elements to enhance cache coherance and better parallel performace. */

    const int nnode = coord.size();
    const int nelem = connectivity.size();
    const int nseg = segment.size();

    //
    // find the longest/shortest dimension
    //
    double_vec lengths = {param.mesh.xlength,
#ifdef THREED
                          param.mesh.ylength,
#endif
                          param.mesh.zlength};
    size_t_vec idx(NDIMS);
    sortindex(lengths, idx);
    int dmin, dmax;
#ifdef THREED
    int dmid;   // the middle axis only enters the 3D sort weight below
#endif

    if (param.mesh.meshing_elem_shape == 0) {
        dmin = idx[0];
#ifdef THREED
        dmid = idx[1];
#endif
        dmax = idx[NDIMS-1];
    } else {
#ifdef THREED
        dmax = 0;
        dmid = NDIMS-2;
#else
        dmax = NDIMS-2;
#endif
        dmin = NDIMS-1;
    }
    //
    // sort coordinate of nodes and element centers
    //
    double_vec wn(nnode);
    const double f = 1e-3;
    for(int i=0; i<nnode; i++) {
        wn[i] = coord[i][dmax]
#ifdef THREED
            + f * coord[i][dmid]
#endif
            + f * f * coord[i][dmin];
    }

    double_vec we(nelem);
    for(int i=0; i<nelem; i++) {
        ConstConnAccessor conn = connectivity[i];
        we[i] = wn[conn[0]] + wn[conn[1]]
#ifdef THREED
            + wn[conn[2]]
#endif
            + wn[conn[NODES_PER_ELEM-1]];
    }

    // arrays to store the result of sorting
    int_vec nd_idx(nnode);
    int_vec el_idx(nelem);
    sortindex(wn, nd_idx);
    sortindex(we, el_idx);

    // inverse permutation
    int_vec nd_inv(nnode);
    for(int i=0; i<nnode; i++)
      nd_inv[nd_idx[i]] = i;

    //
    // renumbering
    //
    array_t coord2(nnode);
    for(int i=0; i<nnode; i++) {
        int n = nd_idx[i];
        for(int j=0; j<NDIMS; j++)
            coord2[i][j] = coord[n][j];
    }
    coord.steal_ref(coord2);

    conn_t conn2(nelem);
    for(int i=0; i<nelem; i++) {
        int n = el_idx[i];
        for(int j=0; j<NODES_PER_ELEM; j++) {
            int k = connectivity[n][j];
            conn2[i][j] = nd_inv[k];
        }
    }
    connectivity.steal_ref(conn2);

    segment_t seg2(nseg);
    for(int i=0; i<nseg; i++) {
        for(int j=0; j<NDIMS; j++) {
            int k = segment[i][j];
            seg2[i][j] = nd_inv[k];
        }
    }
    segment.steal_ref(seg2);

    if (regattr != NULL) {
        regattr_t regattr2(nelem);
        for(int i=0; i<nelem; i++) {
            int n = el_idx[i];
            regattr2[i][0] = (*regattr)[n][0];
        }
        regattr->steal_ref(regattr2);
    }

    if (init_elem_size_n != nullptr) {
        if (init_elem_size_n->size() > 0) {
            double_vec hint2(nnode);
            for (int i=0; i<nnode; ++i) {
                hint2[i] = (*init_elem_size_n)[nd_idx[i]];
            }
            *init_elem_size_n = std::move(hint2);
        }
    }
}


void create_boundary_flags2(uint_vec &bcflag, int nseg,
                            const segment_t &segment, const segflag_t& segflag)
{
    for (int i=0; i<nseg; ++i) {
        uint flag = static_cast<uint>(segflag[i][0]);
        ConstSegmentAccessor n = segment[i];
        for (int j=0; j<NODES_PER_FACET; ++j) {
            bcflag[n[j]] |= flag;
        }
    }
}


void create_boundary_flags(Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // allocate and init to 0
    if (var.bcflag) delete var.bcflag;
    var.bcflag = new uint_vec(var.nnode);

    create_boundary_flags2(*var.bcflag, var.segment->size(),
                           *var.segment, *var.segflag);
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void create_boundary_nodes(Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    /* var.bnodes[i] contains a list of nodes on the i-th boundary.
     * (See constants.hpp for the order of boundaries.)
     */
    for (std::size_t i=0; i<var.bcflag->size(); ++i) {
        uint f = (*var.bcflag)[i];
        for (int j=0; j<nbdrytypes; ++j) {
            if (f & (1<<j)) {
                // this node belongs to a boundary
                (var.bnodes[j])->push_back(i);
            }
        }
    }

    // for (int j=0; j<nbdrytypes; ++j) {
    //     std::cout << "boundary " << j << '\n';
    //     print(std::cout, var.bnodes[j]);
    //     std::cout << '\n';
    // }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void create_top_elems(Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    const int top_bdry = iboundz1;
    const int_vec& top_tmp = *var.bnodes[top_bdry];
    const std::size_t ntop = top_tmp.size();
    int_vec top_ind(ntop,0);
    int_vec top_nodes(ntop,0);

    for (size_t i=0;i<ntop;i++)
        top_ind[i] = i;

    sort(top_ind.begin(), top_ind.end(), [&](const int& a, const int& b) {
            return ((*var.coord)[top_tmp[a]][0] < (*var.coord)[top_tmp[b]][0]);
    });

    for (size_t i=0;i<ntop;i++)
        top_nodes[i] = top_tmp[top_ind[i]];
            


    int_vec telems;
    telems.reserve(100);
    std::set<int> elem_set;

    //  get top elements
    for (std::size_t i=0; i<ntop; i++)
    {
        int n = top_nodes[i];
        // grep the element connected surface
        const int npatch = var.support.size(n);
        const int* patch = var.support.patch(n);
        for (int j=0; j<npatch; j++)
            //use set to avoid duplicated elements.
            elem_set.insert(patch[j]);
    }
    // turn set to vector
    for (auto it = elem_set.begin(); it != elem_set.end();it++)
        telems.push_back(*it);

    var.ntop_elems = telems.size();

    var.top_elems = new int_vec(var.ntop_elems);
    for (int i=0; i<var.ntop_elems; i++)
        (*var.top_elems)[i] = telems[i];

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

// Inverse of elem_and_nodes, which the caller must have filled. Same two passes
// as create_support(); lidx is left empty because only facet ids are read.
static void create_support_surf(SurfaceInfo& surfinfo, int ntop, size_t etop)
{
    Support& sup = surfinfo.support_surf;
    // assign, not resize: pass 1 accumulates, so idx must start zeroed.
    sup.idx_data.assign(ntop + 1, 0);

    // Pass 1: count the facets touching each surface node.
    for (size_t i=0; i<etop; i++)
        for (int k=0; k<NDIMS; k++)
            sup.idx_data[(*surfinfo.elem_and_nodes)[i][k] + 1]++;
    // prefix-sum counts into start offsets
    for (int n=1; n<=ntop; ++n)
        sup.idx_data[n] += sup.idx_data[n-1];

    sup.arr_data.resize(sup.idx_data[ntop]);

    // Pass 2: fill, in the same facet order as pass 1, so each node's patch comes
    // out ascending in facet id -- the order the total_dx gather sums in.
    {
        int_vec cursor(sup.idx_data.begin(), sup.idx_data.end() - 1);
        for (size_t i=0; i<etop; i++)
            for (int k=0; k<NDIMS; k++)
                sup.arr_data[cursor[(*surfinfo.elem_and_nodes)[i][k]]++] = i;
    }
    sup.rebind();
}


void update_surface_info(const Variables& var, SurfaceInfo& surfinfo)
{

    const size_t etop = var.bfacets[iboundz1]->size();
    const size_t ntop = var.bnodes[iboundz1]->size();

    const int top_bdry = iboundz1;
    const int_vec& top_tmp = *var.bnodes[top_bdry];
    int_vec top_ind(ntop,0);
    int_vec top_nodes(ntop,0);
    double_vec top_x(ntop,0.);

    surfinfo.ntop = ntop;
    surfinfo.etop = etop;
    for (size_t i=0;i<ntop;i++)
        top_x[i] = (*var.coord)[top_tmp[i]][0];

    for (size_t i=0;i<ntop;i++)
        top_ind[i] = i;
    sort(top_ind.begin(), top_ind.end(),
        [&](const int& a, const int& b) {
            return (top_x[a] < top_x[b]);
        }
    );
    for (size_t i=0;i<ntop;i++)
        top_nodes[i] = top_tmp[top_ind[i]];
    delete surfinfo.top_nodes;
    surfinfo.top_nodes =  new int_vec(top_nodes.begin(),top_nodes.end());

    // create globle --> local map
    surfinfo.arctop_nodes.clear();
    for (size_t i=0; i<ntop; i++)
        surfinfo.arctop_nodes[(*surfinfo.top_nodes)[i]] = i;
    
    delete surfinfo.dh;
    surfinfo.dh = new double_vec(ntop,0.);
    delete surfinfo.dhacc;
    surfinfo.dhacc = new double_vec(var.nnode,0);
    delete surfinfo.top_facet_elems;
    surfinfo.top_facet_elems = new int_vec(etop,0);
    delete surfinfo.elem_and_nodes;
    surfinfo.elem_and_nodes = new segment_t(etop);

    delete surfinfo.total_dx;
    surfinfo.total_dx = new double_vec(var.nnode,0.);
    delete surfinfo.total_slope;
    surfinfo.total_slope = new double_vec(var.nnode,0.);


    for (size_t i=0; i<etop; i++) {
        (*surfinfo.top_facet_elems)[i] = (*(var.bfacets[iboundz1]))[i].first;

        // Same facets in the same order, so connectivity_surface has already gathered
        // this facet's global node ids: elem_and_nodes is just those mapped into the
        // surface numbering -- one contiguous read, not a second random gather.
        for (int k=0; k<NDIMS; k++)
            (*surfinfo.elem_and_nodes)[i][k] =
                surfinfo.arctop_nodes[(*var.connectivity_surface)[i][k]];
    }
    create_support_surf(surfinfo, ntop, etop);

}

void create_surface_info(const Param& param, const Variables& var, SurfaceInfo& surfinfo)
{

    const size_t etop = var.bfacets[iboundz1]->size();
    const size_t ntop = var.bnodes[iboundz1]->size();

    const int top_bdry = iboundz1;
    const int_vec& top_tmp = *var.bnodes[top_bdry];
    int_vec top_ind(ntop,0);
    int_vec top_nodes(ntop,0);

    double_vec top_x(ntop,0.);

    for (size_t i=0;i<ntop;i++)
        top_x[i] = (*var.coord)[top_tmp[i]][0];

    for (size_t i=0;i<ntop;i++)
        top_ind[i] = i;
    sort(top_ind.begin(), top_ind.end(),
        [&](const int& a, const int& b) {
            return (top_x[a] < top_x[b]);
        }
    );
    for (size_t i=0;i<ntop;i++)
        top_nodes[i] = top_tmp[top_ind[i]];
    surfinfo.top_nodes =  new int_vec(top_nodes.begin(),top_nodes.end());

    surfinfo.ntop = ntop;
    surfinfo.etop = etop;
    surfinfo.max_surf_vel = 0;
    surfinfo.base_level = param.control.surf_base_level;
    surfinfo.surf_diff = param.control.surface_diffusivity;
    surfinfo.diff_ratio_terrig = param.control.surf_diff_ratio_terrig;
    surfinfo.diff_ratio_marine = param.control.surf_diff_ratio_marine;
    surfinfo.depo_universal = param.control.surf_depo_universal;
    surfinfo.terrig_diffusivity = param.control.terrig_sediment_diffusivity;
    surfinfo.terrig_dpeth_coeff = param.control.terrig_depth_coefficient;

    // create globle --> local map
    surfinfo.arctop_nodes.clear();
    for (size_t i=0; i<ntop; i++)
        surfinfo.arctop_nodes[(*surfinfo.top_nodes)[i]] = i;
    surfinfo.dh = new double_vec(ntop,0.);
    surfinfo.dhacc = new double_vec(var.nnode,0);
    surfinfo.edvacc_surf = new double_vec(etop,0);
    surfinfo.top_facet_elems = new int_vec(etop,0);
    surfinfo.elem_and_nodes = new segment_t(etop);

    surfinfo.src_locs = new double_vec(2,0.);

    surfinfo.total_dx = new double_vec(var.nnode,0.);
    surfinfo.total_slope = new double_vec(var.nnode,0.);

    for (size_t i=0; i<etop; i++) {
        (*surfinfo.top_facet_elems)[i] = (*(var.bfacets[iboundz1]))[i].first;

        // Same facets in the same order, so connectivity_surface has already gathered
        // this facet's global node ids: elem_and_nodes is just those mapped into the
        // surface numbering -- one contiguous read, not a second random gather.
        for (int k=0; k<NDIMS; k++)
            (*surfinfo.elem_and_nodes)[i][k] =
                surfinfo.arctop_nodes[(*var.connectivity_surface)[i][k]];
    }
    create_support_surf(surfinfo, ntop, etop);

    //***** to do *****
//    surface_edhacc_geometry_interpolation(var,info);

}
namespace {

    struct OrderedInt
    {
#ifdef THREED
        int a, b, c;
        OrderedInt(int x, int y, int z)
        {
            // ensure that a < b < c
            if (x < y) {
                if (y < z)
                    a = x, b = y, c = z;
                else {
                    if (x < z)
                        a = x, b = z, c = y;
                    else
                        a = z, b = x, c = y;
                }
            }
            else {
                if (x < z)
                    a = y, b = x, c = z;
                else {
                    if (y < z)
                        a = y, b = z, c = x;
                    else
                        a = z, b = y, c = x;
                }
            }
        }

        bool operator==(OrderedInt &rhs)
        {
            return a==rhs.a && b==rhs.b && c==rhs.c;
        }

#else

        int a, b;
        OrderedInt(int x, int y)
        {
            // ensure that a < b
            if (x < y)
                a = x, b = y;
            else
                a = y, b = x;
        }

        bool operator==(OrderedInt &rhs)
        {
            return a==rhs.a && b==rhs.b;
        }
#endif
    };
}


void create_boundary_facets(Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    /* var.bfacets[i] contains a list of facets (or segments in 2D)
     * on the i-th boundary. (See constants.hpp for the order of boundaries.)
     */

    // Looping through var.segment
    #pragma omp parallel default(none) shared(var, NODE_OF_FACET, std::cerr) // TODO: fix no deterministic behavior
    {
        // local storage for each thread
        int_pair_vec bfacet_local[nbdrytypes];

        #pragma omp for nowait
        for (int i=0; i<var.nseg; ++i) {
            uint flag = static_cast<uint>((*var.segflag)[i][0]);
            if ((flag & BOUND_ANY) == 0) continue; // not a boundary facet
            // Nodes of this facet
            OrderedInt af((*var.segment)[i][0], (*var.segment)[i][1]
#ifdef THREED
                          , (*var.segment)[i][2]
#endif
                          );

            // Finding the corresponding element and facet #
            bool found = false;
            for (int e=0; e<var.nelem && !found; ++e) {
                ConstConnAccessor conn = (*var.connectivity)[e];
                for (int f=0; f<FACETS_PER_ELEM && !found; ++f) {
                    uint facet_flag = (*var.bcflag)[conn[NODE_OF_FACET[f][0]]] &
                                      (*var.bcflag)[conn[NODE_OF_FACET[f][1]]]
#ifdef THREED
                                      & (*var.bcflag)[conn[NODE_OF_FACET[f][2]]]
#endif
                                      ;

                    if ((flag & facet_flag) == 0U) continue; // skip

                    OrderedInt bf(conn[NODE_OF_FACET[f][0]], conn[NODE_OF_FACET[f][1]]
#ifdef THREED
                                  , conn[NODE_OF_FACET[f][2]]
#endif
                                  );
                    if (af == bf) {
                        for (int k=0; k<nbdrytypes; ++k) {
                            if (flag == (1U << k)) {
                                bfacet_local[k].emplace_back(e, f);
                                found = true; // break out of nested loops
                                break;
                            }
                        }
                    }
                }
            }
            // not found
            if (!found) {
                #pragma omp critical
                {
                    std::cerr << "Error: " << i << "-th segment is not on any element\n";
                    die(EXIT_INTERNAL_UNREACHABLE);
                }
            }
        }

        // Merge local results into global results
        #pragma omp critical
        {
            for (int k = 0; k < nbdrytypes; ++k) {
                var.bfacets[k]->insert(
                    var.bfacets[k]->end(),
                    bfacet_local[k].begin(),
                    bfacet_local[k].end());
            }
        }

        #pragma omp barrier

        #pragma omp for
        for (int n=0; n<nbdrytypes; ++n) {
            std::sort(var.bfacets[n]->begin(), var.bfacets[n]->end(),
                    [](const int_pair &a, const int_pair &b) {
                        return a.first < b.first;
                    });
        }
    }

    int etop = var.bfacets[iboundz1]->size();
    var.connectivity_surface = new conn_t(etop);

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,NODE_OF_FACET,etop)
#endif
    #pragma acc parallel loop gang vector
    for (int i=0; i<etop; ++i) {
        int e = (*var.bfacets[iboundz1])[i].first;
        int f = (*var.bfacets[iboundz1])[i].second;
        for (int j=0; j<NODES_PER_FACET; ++j) {
            (*var.connectivity_surface)[i][j] = (*var.connectivity)[e][ NODE_OF_FACET[f][j] ];
        }
    }

    // for (int n=0; n<nbdrytypes; ++n) {
    //     std::cout << "boundary facet " << n << ":\n";
    //     print(std::cout, var.bfacets[n]);
    //     std::cout << '\n';
    //     for (int i=0; i<var.bfacets[n].size(); ++i) {
    //         int e = var.bfacets[n][i].first;
    //         int f = var.bfacets[n][i].second;
    //         const int *conn = (*var.connectivity)[e];
    //         std::cout << i << ", " << e << ":";
    //         for (int j=0; j<NODES_PER_FACET; ++j) {
    //             std::cout << " " << conn[NODE_OF_FACET[f][j]];
    //         }
    //         std::cout << '\n';
    //     }
    //     std::cout << '\n';
    // }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}



void create_support(Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // assign, not resize: pass 1 accumulates, so idx must start zeroed.
    var.support.idx_data.assign(var.nnode + 1, 0);

    // Pass 1: count the elements supporting each node.
    for (int e=0; e<var.nelem; ++e) {
        ConstConnAccessor conn = (*var.connectivity)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i)
            var.support.idx_data[conn[i] + 1]++;
    }
    // prefix-sum counts into start offsets
    for (int n=1; n<=var.nnode; ++n)
        var.support.idx_data[n] += var.support.idx_data[n-1];

    int nsup = var.support.idx_data[var.nnode];

    var.support.arr_data.resize(nsup);
    var.support.lidx_data.resize(nsup);

    // Pass 2: fill, in the same element order as pass 1, so each node's patch
    // comes out ascending in element id -- the order every gather sums in.
    {
        int_vec cursor(var.support.idx_data.begin(), var.support.idx_data.end() - 1);
        for (int e=0; e<var.nelem; ++e) {
            ConstConnAccessor conn = (*var.connectivity)[e];
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                const int n = conn[i];
                const int slot = cursor[n]++;
                var.support.arr_data[slot] = e;
                var.support.lidx_data[slot] = i;   // free here: it is the loop counter
            }
        }
    }
    var.support.rebind();

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

// Currently unused: all three call sites are commented out, so var.neighbor /
// var.contact / var.ctmp are written and read by nobody.
void create_neighbor(Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    var.neighbor = new conn_t(var.nelem, -1);
    var.contact = new int_pair_vec(var.nelem*4);
    int ncontact = 0;

    // create the inverse mapping of connectivity
#ifndef ACC
#ifdef GPP1X
    // for Apple clang version 17.0.0
    // future version of clang might fix this problem
    #pragma omp parallel for default(none) \
        shared(var, NODE_OF_FACET, ncontact, NODES_PER_ELEM) collapse(2)
#else
    #pragma omp parallel for default(none) \
        shared(var, NODE_OF_FACET, ncontact) collapse(2)
#endif
#endif
    #pragma acc parallel loop collapse(2) copy(ncontact) async
    for (int e=0; e<var.nelem; ++e) {
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            // No "already set" early-exit: it would race the cross-row write below, and the
            // `neigh > e` filter already registers each internal facet exactly once.
            int n[NDIMS], n2[NDIMS];
            for (int j=0; j<NDIMS; ++j)
                n[j] = (*var.connectivity)[e][NODE_OF_FACET[i][j]];

            // sort the nodes in n
            for (int j=0; j<NDIMS-1; ++j) {
                for (int k=j+1; k<NDIMS; ++k) {
                    if (n[j] > n[k])
                        std::swap(n[j], n[k]);
                }
            }

            // CSR row, not a copy of it: copying the patch heap-allocated per iteration.
            const int nsup = var.support.size(n[0]);
            const int* sup = var.support.patch(n[0]);
            bool found = false;
            for (int j=0; j<nsup && !found; ++j) {
                int neigh = sup[j];
                if (neigh > e) {
                    ConstConnAccessor conn2 = (*var.connectivity)[neigh];
                    for (int k=0; k<NODES_PER_ELEM; ++k) {
                        bool match = true;
                        for (int l=0; l<NDIMS; ++l)
                            n2[l] = conn2[NODE_OF_FACET[k][l]];

                        // sort the nodes in n2
                        for (int l=0; l<NDIMS-1; ++l) {
                            for (int m=l+1; m<NDIMS; ++m) {
                                if (n2[l] > n2[m])
                                    std::swap(n2[l], n2[m]);
                            }
                        }
                        for (int l=0; l<NDIMS; ++l) {
                            if (n[l] != n2[l]) {
                                match = false;
                                break;
                            }
                        }
                        if (match) {
                            // found a neighbor element
                            int pos;
#ifndef ACC
                            #pragma omp atomic capture
#endif
                            #pragma acc atomic capture
                            pos = ncontact++;

                            (*var.contact)[pos*2] = {e, i};
                            (*var.contact)[pos*2+1] = {neigh, k};

                            (*var.neighbor)[e][i] = pos*2+1;
                            (*var.neighbor)[neigh][k] = pos*2;
                            found = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    #pragma acc wait

    var.contact->resize(ncontact*2);
    var.ncontact = ncontact;
    var.ctmp = new double_vec(ncontact);

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

void create_elemmarkers(const Param& param, Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    var.elemmarkers = new int_vec2D( var.nelem, int_vec(param.mat.nmat, 0) );
    var.markers_in_elem = new int_vec2D(var.nelem, int_vec(0));
    if (param.control.has_hydration_processes) {
        var.hydrous_elemmarkers = new Array2D<int,1>( var.nelem, 0 );
        var.hydrous_markers_in_elem = new int_vec2D(var.nelem, int_vec(0));   
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void create_markers(const Param& param, Variables& var)
{
    var.markersets.push_back(new MarkerSet(param, var, std::string("markerset")));
    if (param.control.has_hydration_processes) {
        var.hydrous_marker_index = var.markersets.size();
        var.markersets.push_back(new MarkerSet(std::string("hydrous-markerset")));
    }
}


void create_new_mesh(const Param& param, Variables& var)
{
    switch (param.mesh.meshing_option) {
    case 1:
        if (param.mesh.meshing_elem_shape == 0) {
            new_mesh_uniform_resolution(param, var);
        } else if (param.mesh.meshing_elem_shape == 1) {
            new_mesh_regular(param, var);
        } else if (param.mesh.meshing_elem_shape == 2) {
            new_mesh_regular_equilateral(param, var);
        } else {
            std::cout << "Error: unknown meshing_elem_shape: " << param.mesh.meshing_elem_shape << '\n';
            die(EXIT_CONFIG_VALUE);
        }
        break;
    case 2:
        new_mesh_refined_zone(param, var);
        break;
    case 90:
    case 91:
        new_mesh_from_polyfile(param, var);
        break;
    case 95:
#ifdef USEEXODUS
        new_mesh_from_exofile(param, var);
#else
        std::cout << "Error: Install Exodus library and rebuild with 'useexo' turned on in Makefile." << std::endl;
        die(EXIT_UNSUPPORTED_LIB);
#endif
        break;
    default:
        std::cout << "Error: unknown meshing option: " << param.mesh.meshing_option << '\n';
        die(EXIT_CONFIG_VALUE);
    }

    if (param.mesh.is_discarding_internal_segments)
        discard_internal_segments(var.nseg, *var.segment, *var.segflag);

    renumbering_mesh(param, *var.coord, *var.connectivity, *var.segment, var.regattr, var.init_elem_size_n);

    // std::cout << "segment:\n";
    // print(std::cout, *var.segment);
    // std::cout << '\n';
    // std::cout << "segflag:\n";
    // print(std::cout, *var.segflag);
    // std::cout << '\n';
}

// Average per-node scalar values to element centroids (simple arithmetic mean over vertices).
template<typename T_in, typename T_out>
void average_nodal_to_elem(T_in nodal, const conn_t &connectivity,
                                   int nelem, T_out elem, bool is_surface)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    int nnodes_per_elem = is_surface ? NODES_PER_FACET : NODES_PER_ELEM;

#ifndef ACC
    #pragma omp parallel for default(none) \
            shared(nodal, connectivity, nelem, elem, nnodes_per_elem)
#endif
    #pragma acc parallel loop gang vector async
    for (int e = 0; e < nelem; ++e) {
        ConstConnAccessor conn = connectivity[e];
        double avg = 0.0;
        for (int k = 0; k < nnodes_per_elem; ++k)
            avg += nodal[conn[k]];
        elem[e] = avg / nnodes_per_elem;
    }
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

template
void average_nodal_to_elem<array_t::ConstComponentAccessor, array_t::ComponentAccessor>(
    array_t::ConstComponentAccessor nodal, 
    const conn_t &connectivity, int nelem, array_t::ComponentAccessor elem, bool);
template
void average_nodal_to_elem<tensor_t::ConstComponentAccessor, tensor_t::ComponentAccessor>(
    tensor_t::ConstComponentAccessor nodal, 
    const conn_t &connectivity, int nelem, tensor_t::ComponentAccessor elem, bool);
template
void average_nodal_to_elem<const double*, double*>(const double* nodal,
    const conn_t &connectivity, int nelem, double* elem, bool);


//
// VirtualTopology -- background Cartesian grid of per-cell marker-id
// buckets.  See mesh.hpp for the interface and the algorithm's contract.
//

VirtualTopology::VirtualTopology(const double* domain_lo, const double* domain_hi,
                                 double resolution)
    : dx_(resolution)
{
    for (int d=0; d<NDIMS; d++) {
        lo_[d] = domain_lo[d];
        dims_[d] = std::max(1, static_cast<int>(std::ceil((domain_hi[d] - domain_lo[d]) / dx_)));
    }
    nz_ = dims_[NDIMS-1];
    ncol_ = 1;
    for (int d=0; d<NDIMS-1; d++)
        ncol_ *= dims_[d];

    cell_start_.assign(static_cast<std::size_t>(ncol_) * nz_ + 1, 0);
    dominant_.assign(static_cast<std::size_t>(ncol_) * nz_, 0);
}


void VirtualTopology::serialize(int_vec& ints, double_vec& dbls) const
{
    const std::size_t ncells = static_cast<std::size_t>(ncol_) * nz_;

    ints.clear();
    ints.push_back(layout_tag);
    for (int d=0; d<NDIMS; d++) ints.push_back(dims_[d]);
    ints.push_back(static_cast<int>(label_.size()));
    ints.insert(ints.end(), label_.begin(), label_.end());
    ints.push_back(built_last_id_);

    ints.push_back(xext_built_ ? 1 : 0);
    ints.push_back(static_cast<int>(
        std::count_if(dominant_.begin(), dominant_.end(), [](char f) { return f != 0; })));

    dbls.clear();
    for (int d=0; d<NDIMS; d++) dbls.push_back(lo_[d]);
    dbls.push_back(dx_);

    // One loop for ids AND extents so the two streams stay aligned by
    // position. With xext_built_ every dominant cell has an extent:
    // build_serpentinization_grid() records one for each.
    if (xext_built_ && xext_.size() != ncells * xext_stride)
        die(EXIT_INTERNAL_ASSERT, "VirtualTopology: the grid has no wet-crust extents");
    for (std::size_t c=0; c<ncells; c++) {
        if (!dominant_[c]) continue;
        ints.push_back(static_cast<int>(c));
        if (xext_built_) {
            for (int k=0; k<xext_stride; k++) {
                if (std::isnan(xext_[c*xext_stride + k]))
                    die(EXIT_INTERNAL_ASSERT, "VirtualTopology: a dominant cell has no wet-crust extent");
                dbls.push_back(xext_[c*xext_stride + k]);
            }
        }
    }
}


VirtualTopology* VirtualTopology::deserialize(const int_vec& ints, const double_vec& dbls)
{
    if (ints.empty() || ints[0] != layout_tag)
        die(EXIT_IO_RESTART, "vtopo_int has another VirtualTopology layout: restart from a checkpoint this build wrote");
    VirtualTopology* vt = new VirtualTopology();   // private default ctor
    // Every count and id below comes from the file, so every read is bounded:
    // a corrupt checkpoint dies here instead of reading past either array.
    std::size_t p = 1;
    auto next = [&]() {
        if (p >= ints.size())
            die(EXIT_IO_RESTART, "vtopo_int ends before its counts do");
        return ints[p++];
    };
    auto check = [](bool ok) {
        if (!ok)
            die(EXIT_IO_RESTART, "vtopo_int holds a count or id outside its range");
    };
    std::size_t ncells = 1;
    for (int d=0; d<NDIMS; d++) {
        vt->dims_[d] = next();
        check(vt->dims_[d] > 0);
        ncells *= vt->dims_[d];
        check(ncells <= std::numeric_limits<int>::max());   // cell ids are ints
    }
    vt->nz_ = vt->dims_[NDIMS-1];
    vt->ncol_ = 1;
    for (int d=0; d<NDIMS-1; d++) vt->ncol_ *= vt->dims_[d];

    const int nlabels = next();
    check(nlabels >= 0 && p + nlabels <= ints.size());
    vt->label_.assign(ints.begin()+p, ints.begin()+p+nlabels);
    p += nlabels;
    for (int c : vt->label_)
        check(c >= 0 && static_cast<std::size_t>(c) < ncells);
    vt->built_last_id_ = next();

    // Wet-crust extents ride in dbls, one per dominant id in the same order (see
    // serialize()), exactly when the stored xext_built_ says so.
    std::size_t r = NDIMS + 1;   // dbls cursor: lo_ and dx_ come first
    if (dbls.size() < r)
        die(EXIT_IO_RESTART, "vtopo_dbl ends before its counts do");
    auto next_dbl = [&]() {
        if (r >= dbls.size())
            die(EXIT_IO_RESTART, "vtopo_dbl ends before its counts do");
        return dbls[r++];
    };
    {
        vt->dominant_.assign(ncells, 0);
        vt->xext_built_ = (next() != 0);
        if (vt->xext_built_)
            vt->xext_.assign(ncells * xext_stride, std::numeric_limits<double>::quiet_NaN());
        const int ndominant = next();
        check(ndominant >= 0);
        for (int i=0; i<ndominant; i++) {
            const int c = next();
            check(c >= 0 && static_cast<std::size_t>(c) < ncells);
            vt->dominant_[c] = 1;
            if (vt->xext_built_)
                for (int k=0; k<xext_stride; k++) {
                    // A stored extent is a wet-crust-marker coordinate, so finite.
                    const double v = next_dbl();
                    if (!std::isfinite(v))
                        die(EXIT_IO_RESTART, "vtopo_dbl holds a wet-crust extent that is not finite");
                    vt->xext_[c*xext_stride + k] = v;
                }
        }
    }

    vt->cell_start_.assign(ncells + 1, 0);
    vt->fill_buckets();

    // The reads stop at either array's end; an array longer than its counts
    // is caught here.
    if (p != ints.size() || r != dbls.size())
        die(EXIT_IO_RESTART, "vtopo_int / vtopo_dbl do not decode to their end");

    for (int d=0; d<NDIMS; d++) vt->lo_[d] = dbls[d];
    vt->dx_ = dbls[NDIMS];
    bool finite = std::isfinite(vt->dx_) && vt->dx_ > 0;
    for (int d=0; d<NDIMS; d++) finite = finite && std::isfinite(vt->lo_[d]);
    if (!finite)
        die(EXIT_IO_RESTART, "vtopo_dbl holds a grid origin or cell size that is not finite and positive");
    return vt;
}


int VirtualTopology::column_index(const double* p) const
{
    // Fold the NDIMS-1 horizontal axes: i in 2D, i*dims_[1]+j in 3D.
    int col = 0;
    for (int d=0; d<NDIMS-1; d++) {
        int i = static_cast<int>((p[d] - lo_[d]) / dx_);
        i = std::min(std::max(i, 0), dims_[d] - 1);
        col = col * dims_[d] + i;
    }
    return col;
}


int VirtualTopology::vertical_index(const double* p) const
{
    int k = static_cast<int>((p[NDIMS-1] - lo_[NDIMS-1]) / dx_);
    return std::min(std::max(k, 0), nz_ - 1);
}


void VirtualTopology::build_index(const array_t& coords, int npoints, int last_id)
{
    built_last_id_ = last_id;
    label_.assign(npoints, -1);
    for (int m=0; m<npoints; m++) {
        double p[NDIMS];
        for (int d=0; d<NDIMS; d++)
            p[d] = coords[m][d];
        label_[m] = cell_of(p);
    }
    fill_buckets();
}


void VirtualTopology::fill_buckets()
{
    // Serial: it runs once per rebuild, not per phase_changes() call.
    const std::size_t ncells = cell_start_.size() - 1;
    std::fill(cell_start_.begin(), cell_start_.end(), 0);
    for (int c : label_)
        ++cell_start_[c + 1];
    for (std::size_t c = 0; c < ncells; ++c)
        cell_start_[c + 1] += cell_start_[c];
    ids_.resize(label_.size());
    int_vec next(cell_start_.begin(), cell_start_.end() - 1);
    for (int m = 0; m < static_cast<int>(label_.size()); ++m)
        ids_[next[label_[m]]++] = m;
}


int VirtualTopology::neighbor_cell(int cell, int axis, int dir) const
{
    int idx[NDIMS];
    int rem = cell;
    for (int d=NDIMS-1; d>=0; d--) {
        idx[d] = rem % dims_[d];
        rem /= dims_[d];
    }
    idx[axis] += dir;
    if (idx[axis] < 0 || idx[axis] >= dims_[axis])
        return -1;
    int flat = 0;
    for (int d=0; d<NDIMS; d++)
        flat = flat * dims_[d] + idx[d];
    return flat;
}


VirtualTopology::IdRange VirtualTopology::cell_markers(int cell) const
{
    if (cell < 0 || cell >= num_cells())
        die(EXIT_INTERNAL_ASSERT, "VirtualTopology::cell_markers: cell outside [0, num_cells())");
    const int* base = ids_.data();
    return IdRange{base + cell_start_[cell], base + cell_start_[cell + 1]};
}


void VirtualTopology::clear_dominant()
{
    dominant_.assign(dominant_.size(), 0);
    xext_.clear();
    xext_built_ = false;
}


void VirtualTopology::set_dominant(int cell)
{
    dominant_[cell] = 1;
}


void VirtualTopology::set_dominant_xext(int cell, const double* xlo, const double* xhi)
{
    if (xext_.empty())
        xext_.assign(dominant_.size() * xext_stride, std::numeric_limits<double>::quiet_NaN());
    for (int d=0; d<NDIMS-1; d++) {
        xext_[cell*xext_stride + d] = xlo[d];
        xext_[cell*xext_stride + NDIMS-1 + d] = xhi[d];
    }
}


void VirtualTopology::mark_xext_built()
{
    // Built extents have their array even when no cell is dominant.
    if (xext_.empty())
        xext_.assign(dominant_.size() * xext_stride, std::numeric_limits<double>::quiet_NaN());
    xext_built_ = true;
}


bool VirtualTopology::dominant_xext_ok(int cell, const double* p, double margin) const
{
    // Called only with the margin on, when xext_built_ holds and so every
    // dominant cell has an extent (see serialize()). Written so a NaN extent
    // (a cell without one) rejects: device code cannot die.
    const double* e = &xext_[static_cast<std::size_t>(cell) * xext_stride];
    for (int d=0; d<NDIMS-1; d++)
        if (!(p[d] >= e[d] - margin && p[d] <= e[NDIMS-1 + d] + margin))
            return false;
    return true;
}


void VirtualTopology::cells_with_dominant_below(int nsteps, std::vector<char>& below) const
{
    below.assign(num_cells(), 0);
    for (int dom = 0; dom < num_cells(); ++dom) {
        if (!dominant_[dom]) continue;
        int c = neighbor_cell(dom, NDIMS-1, +1);
        for (int k = 0; c >= 0 && k < nsteps; ++k) {
            below[c] = 1;
            c = neighbor_cell(c, NDIMS-1, +1);
        }
    }
}


void VirtualTopology::invalidate_dominant_cell(int cell)
{
    dominant_[cell] = 0;
}


