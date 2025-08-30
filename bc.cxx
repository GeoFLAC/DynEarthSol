#include <iostream>
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


#include "bc.hpp"


namespace {

#pragma acc routine seq
void normal_vector_of_facet(int f, const int *conn, const array_t &coord,
                            double *normal, double &zcenter)
{
    int n0 = conn[NODE_OF_FACET[f][0]];
    int n1 = conn[NODE_OF_FACET[f][1]];
#ifdef THREED
    int n2 = conn[NODE_OF_FACET[f][2]];

    // vectors: n0-n1 and n0-n2
    double v01[NDIMS], v02[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = coord[n1][i] - coord[n0][i];
        v02[i] = coord[n2][i] - coord[n0][i];
    }

    // the outward normal vector is parallel to the cross product
    // of v01 and v02, with magnitude = area of the triangle
    normal[0] = (v01[1] * v02[2] - v01[2] * v02[1]) / 2;
    normal[1] = (v01[2] * v02[0] - v01[0] * v02[2]) / 2;
    normal[2] = (v01[0] * v02[1] - v01[1] * v02[0]) / 2;

    zcenter = (coord[n0][2] + coord[n1][2] + coord[n2][2]) / NODES_PER_FACET;
#else
    // the edge vector
    double v01[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = coord[n1][i] - coord[n0][i];
    }

    // the normal vector to the edge, pointing outward
    normal[0] = v01[1];
    normal[1] = -v01[0];

    zcenter = (coord[n0][1] + coord[n1][1]) / NODES_PER_FACET;
#endif
}

}

#pragma acc routine seq
bool is_on_boundary(const uint_vec &bcflag, int node)
{
    uint flag = bcflag[node];
    return flag & BOUND_ANY;
}


double find_max_vbc(const BC &bc)
{
    double max_vbc_val = 1e-12; // equal to 0.03 mm/yr
    if (bc.vbc_x0 % 2 == 1) // odd number indicates fixed velocity component
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x0));
    if (bc.vbc_x1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x1));
    if (bc.vbc_y0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y0));
    if (bc.vbc_y1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y1));
    if (bc.vbc_z0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z0));
    if (bc.vbc_z1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z1));
    if (bc.vbc_n0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n0));
    if (bc.vbc_n1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n1));
    if (bc.vbc_n2 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n2));
    if (bc.vbc_n3 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n3));

    return max_vbc_val;
}


void create_boundary_normals(const Variables &var, array_t &bnormals,
                             std::map<std::pair<int,int>, double*>  &edge_vectors, double_vec& edge_vec, int_vec &edge_vec_idx)
{
    /* This subroutine finds the outward normal unit vectors of boundaries.
     * There are two types of boundaries: ibound{x,y,z}? and iboundn?.
     * The first type ibound{x,y,z}? has well defined "normal" direction
     * as their name implied. Thus they can be curved. (Ex. the top boundary
     * can have topography.)
     * The second type iboundn?, on the other hand, must be a plane so that
     * their vbc are well defined.
     *
     * If the normal component of the boundary velocity is fixed, the boundary
     * normal will not change with time.
     */

    #pragma acc wait

    for (int i=0; i<nbdrytypes; i++) {
        double normal[NDIMS] = {0};
        if (var.bfacets[i]->size() == 0) continue;

        for (auto j=var.bfacets[i]->begin(); j<var.bfacets[i]->end(); ++j) {
            int e = j->first;
            int f = j->second;
            double tmp;
            normal_vector_of_facet(f, (*var.connectivity)[e], *var.coord,
                                   normal, tmp);
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
                    std::exit(1);
                }
            }
        }
    }

    for (int i=0; i<nbdrytypes; i++) {
        if (var.bfacets[i]->size() == 0) continue;

        const double eps = 1e-15;
        for (int j=i+1; j<nbdrytypes; j++) {
            if (var.bfacets[j]->size() == 0) continue;
            double *s = new double[NDIMS];  // intersection of two boundaries
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
            edge_vectors[std::make_pair(i, j)] = s;
            edge_vec_idx.push_back(i*nbdrytypes + j);
            edge_vec.push_back(s[0]);
            edge_vec.push_back(s[1]);
#ifdef THREED
            edge_vec.push_back(s[2]);
#endif
        }
    }
}


void apply_vbcs(const Param &param, const Variables &var, array_t &vel)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // meaning of vbc flags (odd: free; even: fixed) --
    // 0: all components free
    // 1: normal component fixed, shear components free
    // 2: normal component free, shear components fixed at 0
    // 3: normal component fixed, shear components fixed at 0
    // 4: normal component free, shear component (not z) fixed, only in 3D
    // 5: normal component fixed at 0, shear component (not z) fixed, only in 3D

    const BC &bc = param.bc;
#ifdef THREED
    // vel period is not ready for 3D

#else

    // new way to find vbc_applied_x
    double t_now = var.time / YEAR2SEC;
    double vbc_applied_x0 = bc.vbc_val_x0 * interp1(bc.vbc_period_x0_time_in_yr,bc.vbc_period_x0_ratio, t_now);
    double vbc_applied_x1 = bc.vbc_val_x1 * interp1(bc.vbc_period_x1_time_in_yr,bc.vbc_period_x1_ratio, t_now);

    // find max and min coordinate for BOUNDX0
    double BOUNDX0_max = 0.;
    double BOUNDX0_min = 0.;
    double BOUNDX1_max = 0.;
    double BOUNDX1_min = 0.;
    bool if_init0 = false;
    bool if_init1 = false;
    for (int i=0; i<var.nnode; ++i) {
        if (! is_on_boundary(*var.bcflag, i)) continue;

        uint flag = (*var.bcflag)[i];
        if (flag & BOUNDX0) {
            const double *x = (*var.coord)[i];
            if (!if_init0) {
               BOUNDX0_max = x[1];
               BOUNDX0_min = x[1];
               if_init0 = true;
            } else {
                if (x[1]>BOUNDX0_max) {
                    BOUNDX0_max = x[1];
                }
                if (x[1]<BOUNDX0_min) {
                    BOUNDX0_min = x[1];
                }
            }
        } else if (flag & BOUNDX1) {
            const double *x = (*var.coord)[i];
            if (!if_init1) {
               BOUNDX1_max = x[1];
               BOUNDX1_min = x[1];
               if_init1 = true;
            } else {
                if (x[1]>BOUNDX1_max) {
                    BOUNDX1_max = x[1];
                }
                if (x[1]<BOUNDX1_min) {
                    BOUNDX1_min = x[1];
                }
            }            
        }
    }
    double BOUNDX0_width = BOUNDX0_max - BOUNDX0_min;
    double BOUNDX1_width = BOUNDX1_max - BOUNDX1_min;

    const double_vec& vbc_vertical_ratios_x0 = var.vbc_vertical_ratio_x0;
    const double_vec& vbc_vertical_ratios_x1 = var.vbc_vertical_ratio_x1;

    double_vec vbc_vertical_divisions_x0(4,0.);
    double_vec vbc_vertical_divisions_x1(4,0.);
    for (int i=0;i<4;i++) {
        vbc_vertical_divisions_x0[i] = - (BOUNDX0_max - var.vbc_vertical_div_x0[i] * BOUNDX0_width);
        vbc_vertical_divisions_x1[i] = - (BOUNDX0_max - var.vbc_vertical_div_x1[i] * BOUNDX0_width);
    }
#endif

    // diverging x-boundary
    // const std::map<std::pair<int,int>, double*>  *edgevec = &(var.edge_vectors);


    int bc_x0 = bc.vbc_x0;
    int bc_x1 = bc.vbc_x1;
    int bc_y0 = bc.vbc_y0;
    int bc_y1 = bc.vbc_y1;
    int bc_z0 = bc.vbc_z0;
    int bc_z1 = bc.vbc_z1;

    double bc_vx0 = bc.vbc_val_x0;
    double bc_vx1 = bc.vbc_val_x1;
    double bc_vy0 = bc.vbc_val_y0;
    double bc_vy1 = bc.vbc_val_y1;
    double bc_vz0 = bc.vbc_val_z0;
    double bc_vz1 = bc.vbc_val_z1;

    if (var.time > var.vbc_val_z1_loading_period) {bc_z1 = 0;}

    double bc_vx0min = bc.vbc_val_division_x0_min;
    double bc_vx0max = bc.vbc_val_division_x0_max;

    double bc_vx0r0 = bc.vbc_val_x0_ratio0;
    double bc_vx0r1 = bc.vbc_val_x0_ratio1;
    double bc_vx0r2 = bc.vbc_val_x0_ratio2;

#ifdef THREED
    #pragma acc parallel loop async
#else
    double zmin = 0;
    for (int k=0; k<var.nnode; ++k) {
        double* tmpx = (*var.coord)[k];
        if ( tmpx[NDIMS-1] < zmin )
            zmin = tmpx[NDIMS-1];
    }
    #pragma acc parallel loop async
#endif
    for (int i=0; i<var.nnode; ++i) {

        // fast path: skip nodes not on boundary
        if (! is_on_boundary(*var.bcflag, i)) continue;

        uint flag = (*var.bcflag)[i];
        double *v = vel[i];

#ifdef THREED
#else
        const double *x = (*var.coord)[i];
        double ratio, rr, dvr;
        double vbc_exact_x0 = vbc_applied_x0 * interp1(vbc_vertical_divisions_x0, vbc_vertical_ratios_x0,-x[1]);
        double vbc_exact_x1 = vbc_applied_x1 * interp1(vbc_vertical_divisions_x1, vbc_vertical_ratios_x1,-x[1]);
#endif
        //
        // X
        //
        if (flag & BOUNDX0) {
            switch (bc_x0) {
            case 0:
                break;
            case 1:
#ifdef THREED
                v[0] = bc_vx0;
#else
                v[0] = vbc_exact_x0;
#endif
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
#ifdef THREED
                v[0] = bc_vx0;
#else
                v[0] = vbc_exact_x0;
                if (bc.bottom_shear_zone_thickness > 0.) {
                    double dz = x[NDIMS-1] - zmin;
                    if (dz < bc.bottom_shear_zone_thickness)
                        v[0] = v[0] * dz / bc.bottom_shear_zone_thickness;
                }
#endif
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc_vx0;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc_vx0;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc_vx0;
                v[1] = 0;
                break;
#endif
            }
        }
        if (flag & BOUNDX1) {
            switch (bc_x1) {
            case 0:
                break;
            case 1:
#ifdef THREED
                v[0] = bc_vx1;
#else
                v[0] = vbc_exact_x1;
#endif
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
#ifdef THREED
                v[0] = bc_vx1;
#else
                v[0] = vbc_exact_x1;
#endif
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc_vx1;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc_vx1;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc_vx1;
                v[1] = 0;
                break;
#endif
            }
        }
#ifdef THREED
        //
        // Y
        //
        if (flag & BOUNDY0) {
            switch (bc_y0) {
            case 0:
                break;
            case 1:
                v[1] = bc_vy0;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = 0;
                v[1] = bc_vy0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc_vy0;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc_vy0;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc_vy0;
                break;
            }
        }
        if (flag & BOUNDY1) {
            switch (bc_y1) {
            case 0:
                break;
            case 1:
                v[1] = bc_vy1;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = bc_vy1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc_vy1;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc_vy1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc_vy1;
                break;
            }
        }
#endif

        //
        // N
        //
        for (int ib=iboundn0; ib<=iboundn3; ib++) {
            const double eps = 1e-15;
            const double *n = (*var.bnormals)[ib]; // unit normal vector

            if (flag & (1 << ib)) {
                double fac = 0;
                switch (var.vbc_types[ib]) {
                case 1:
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS; d++)
                            vn += v[d] * n[d];  // normal velocity

			for (int d=0; d<NDIMS; d++)
                            v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var.vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        vn += v[d] * n[d];  // normal velocity

				    for (int d=0; d<NDIMS; d++)
                                        v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                                }
                                else if (var.vbc_types[ic] == 1) {
                                    const double *edge;
                                    for (int j=0; j<var.edge_vec_idx.size();j++) {
                                        int ei = var.edge_vec_idx[j]/nbdrytypes;
                                        int ej = var.edge_vec_idx[j]%nbdrytypes;
                                        if (ei == ic && ej == ib)
                                            edge = &var.edge_vec[j*NDIMS];
                                    }
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 3:
                    for (int d=0; d<NDIMS; d++)
                        v[d] = var.vbc_values[ib] * n[d];  // v must be normal to n
                    break;
                case 11:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS-1; d++)
                            vn += v[d] * n[d];  // normal velocity

			for (int d=0; d<NDIMS-1; d++)
                            v[d] += (var.vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var.vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS-1; d++)
                                        vn += v[d] * n[d];  // normal velocity

				    for (int d=0; d<NDIMS-1; d++)
                                        v[d] += (var.vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                                }
                                else if (var.vbc_types[ic] == 1) {
                                    const double *edge;
                                    for (int j=0; j<var.edge_vec_idx.size();j++) {
                                        int ei = var.edge_vec_idx[j]/nbdrytypes;
                                        int ej = var.edge_vec_idx[j]%nbdrytypes;
                                        if (ei == ic && ej == ib)
                                            edge = &var.edge_vec[j*NDIMS];
                                    }
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 13:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    for (int d=0; d<NDIMS-1; d++)
                        v[d] = var.vbc_values[ib] * fac * n[d];
                    v[NDIMS-1] = 0;
                    break;
                }
            }
        }

        //
        // Z, must be dealt last
        //

        // fast path: vz is usually free in the models
        if (bc_z0==0 && bc_z1==0) continue;

        if (flag & BOUNDZ0) {
            switch (bc_z0) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc_vz0;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc_vz0;
                break;
            }
        }
        if (flag & BOUNDZ1) {
            switch (bc_z1) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc_vz1;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0.0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc_vz1;
                break;
            case 4:
                v[0] = bc_vz1;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = 0;
                break;
            }
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void apply_vbcs_PT(const Param &param, const Variables &var, array_t &vel)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // meaning of vbc flags (odd: free; even: fixed) --
    // 0: all components free
    // 1: normal component fixed, shear components free
    // 2: normal component free, shear components fixed at 0
    // 3: normal component fixed, shear components fixed at 0
    // 4: normal component free, shear component (not z) fixed, only in 3D
    // 5: normal component fixed at 0, shear component (not z) fixed, only in 3D

    const BC &bc = param.bc;
#ifdef THREED
    // vel period is not ready for 3D

#else

    // new way to find vbc_applied_x
    double t_now = var.time / YEAR2SEC;
    double vbc_applied_x0 = bc.vbc_val_x0 * interp1(bc.vbc_period_x0_time_in_yr,bc.vbc_period_x0_ratio, t_now);
    double vbc_applied_x1 = bc.vbc_val_x1 * interp1(bc.vbc_period_x1_time_in_yr,bc.vbc_period_x1_ratio, t_now);

    // find max and min coordinate for BOUNDX0
    double BOUNDX0_max = 0.;
    double BOUNDX0_min = 0.;
    double BOUNDX1_max = 0.;
    double BOUNDX1_min = 0.;
    bool if_init0 = false;
    bool if_init1 = false;
    for (int i=0; i<var.nnode; ++i) {
        if (! is_on_boundary(*var.bcflag, i)) continue;

        uint flag = (*var.bcflag)[i];
        if (flag & BOUNDX0) {
            const double *x = (*var.coord)[i];
            if (!if_init0) {
               BOUNDX0_max = x[1];
               BOUNDX0_min = x[1];
               if_init0 = true;
            } else {
                if (x[1]>BOUNDX0_max) {
                    BOUNDX0_max = x[1];
                }
                if (x[1]<BOUNDX0_min) {
                    BOUNDX0_min = x[1];
                }
            }
        } else if (flag & BOUNDX1) {
            const double *x = (*var.coord)[i];
            if (!if_init1) {
               BOUNDX1_max = x[1];
               BOUNDX1_min = x[1];
               if_init1 = true;
            } else {
                if (x[1]>BOUNDX1_max) {
                    BOUNDX1_max = x[1];
                }
                if (x[1]<BOUNDX1_min) {
                    BOUNDX1_min = x[1];
                }
            }            
        }
    }
    double BOUNDX0_width = BOUNDX0_max - BOUNDX0_min;
    double BOUNDX1_width = BOUNDX1_max - BOUNDX1_min;

    const double_vec& vbc_vertical_ratios_x0 = var.vbc_vertical_ratio_x0;
    const double_vec& vbc_vertical_ratios_x1 = var.vbc_vertical_ratio_x1;

    double_vec vbc_vertical_divisions_x0(4,0.);
    double_vec vbc_vertical_divisions_x1(4,0.);
    for (int i=0;i<4;i++) {
        vbc_vertical_divisions_x0[i] = - (BOUNDX0_max - var.vbc_vertical_div_x0[i] * BOUNDX0_width);
        vbc_vertical_divisions_x1[i] = - (BOUNDX0_max - var.vbc_vertical_div_x1[i] * BOUNDX0_width);
    }
#endif

    // diverging x-boundary
    // const std::map<std::pair<int,int>, double*>  *edgevec = &(var.edge_vectors);


    int bc_x0 = bc.vbc_x0;
    int bc_x1 = bc.vbc_x1;
    int bc_y0 = bc.vbc_y0;
    int bc_y1 = bc.vbc_y1;
    int bc_z0 = bc.vbc_z0;
    int bc_z1 = bc.vbc_z1;

    double bc_vx0 = 0.0;
    double bc_vx1 = 0.0;
    double bc_vy0 = 0.0;
    double bc_vy1 = 0.0;
    double bc_vz0 = 0.0;
    double bc_vz1 = 0.0;

    if (var.time > var.vbc_val_z1_loading_period) {bc_z1 = 0;}

    double bc_vx0min = bc.vbc_val_division_x0_min;
    double bc_vx0max = bc.vbc_val_division_x0_max;

    double bc_vx0r0 = bc.vbc_val_x0_ratio0;
    double bc_vx0r1 = bc.vbc_val_x0_ratio1;
    double bc_vx0r2 = bc.vbc_val_x0_ratio2;

#ifdef THREED
    #pragma acc parallel loop async
#else
    double zmin = 0;
    for (int k=0; k<var.nnode; ++k) {
        double* tmpx = (*var.coord)[k];
        if ( tmpx[NDIMS-1] < zmin )
            zmin = tmpx[NDIMS-1];
    }
    #pragma acc parallel loop async
#endif
    for (int i=0; i<var.nnode; ++i) {

        // fast path: skip nodes not on boundary
        if (! is_on_boundary(*var.bcflag, i)) continue;

        uint flag = (*var.bcflag)[i];
        double *v = vel[i];

#ifdef THREED
#else
        const double *x = (*var.coord)[i];
        double ratio, rr, dvr;
        double vbc_exact_x0 = vbc_applied_x0 * interp1(vbc_vertical_divisions_x0, vbc_vertical_ratios_x0,-x[1]) * 0.0;
        double vbc_exact_x1 = vbc_applied_x1 * interp1(vbc_vertical_divisions_x1, vbc_vertical_ratios_x1,-x[1]) * 0.0;

#endif
        //
        // X
        //
        if (flag & BOUNDX0) {
            switch (bc_x0) {
            case 0:
                break;
            case 1:
#ifdef THREED
                v[0] = bc_vx0;
#else
                v[0] = vbc_exact_x0;
#endif
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
#ifdef THREED
                v[0] = bc_vx0;
#else
                v[0] = vbc_exact_x0;
                if (bc.bottom_shear_zone_thickness > 0.) {
                    double dz = x[NDIMS-1] - zmin;
                    if (dz < bc.bottom_shear_zone_thickness)
                        v[0] = v[0] * dz / bc.bottom_shear_zone_thickness;
                }
#endif
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc_vx0;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc_vx0;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc_vx0;
                v[1] = 0;
                break;
#endif
            }
        }
        if (flag & BOUNDX1) {
            switch (bc_x1) {
            case 0:
                break;
            case 1:
#ifdef THREED
                v[0] = bc_vx1;
#else
                v[0] = vbc_exact_x1;
#endif
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
#ifdef THREED
                v[0] = bc_vx1;
#else
                v[0] = vbc_exact_x1;
#endif
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc_vx1;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc_vx1;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc_vx1;
                v[1] = 0;
                break;
#endif
            }
        }
#ifdef THREED
        //
        // Y
        //
        if (flag & BOUNDY0) {
            switch (bc_y0) {
            case 0:
                break;
            case 1:
                v[1] = bc_vy0;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = 0;
                v[1] = bc_vy0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc_vy0;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc_vy0;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc_vy0;
                break;
            }
        }
        if (flag & BOUNDY1) {
            switch (bc_y1) {
            case 0:
                break;
            case 1:
                v[1] = bc_vy1;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = bc_vy1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc_vy1;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc_vy1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc_vy1;
                break;
            }
        }
#endif

        //
        // N
        //
        for (int ib=iboundn0; ib<=iboundn3; ib++) {
            const double eps = 1e-15;
            const double *n = (*var.bnormals)[ib]; // unit normal vector

            if (flag & (1 << ib)) {
                double fac = 0;
                switch (var.vbc_types[ib]) {
                case 1:
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS; d++)
                            vn += v[d] * n[d];  // normal velocity

			for (int d=0; d<NDIMS; d++)
                            v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var.vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        vn += v[d] * n[d];  // normal velocity

				    for (int d=0; d<NDIMS; d++)
                                        v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                                }
                                else if (var.vbc_types[ic] == 1) {
                                    const double *edge;
                                    for (int j=0; j<var.edge_vec_idx.size();j++) {
                                        int ei = var.edge_vec_idx[j]/nbdrytypes;
                                        int ej = var.edge_vec_idx[j]%nbdrytypes;
                                        if (ei == ic && ej == ib)
                                            edge = &var.edge_vec[j*NDIMS];
                                    }
                                    // auto edge = edgevec->at(std::make_pair(ic, ib));
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 3:
                    for (int d=0; d<NDIMS; d++)
                        v[d] = var.vbc_values[ib] * n[d];  // v must be normal to n
                    break;
                case 11:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS-1; d++)
                            vn += v[d] * n[d];  // normal velocity

			for (int d=0; d<NDIMS-1; d++)
                            v[d] += (var.vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var.vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS-1; d++)
                                        vn += v[d] * n[d];  // normal velocity

				    for (int d=0; d<NDIMS-1; d++)
                                        v[d] += (var.vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                                }
                                else if (var.vbc_types[ic] == 1) {
                                    const double *edge;
                                    for (int j=0; j<var.edge_vec_idx.size();j++) {
                                        int ei = var.edge_vec_idx[j]/nbdrytypes;
                                        int ej = var.edge_vec_idx[j]%nbdrytypes;
                                        if (ei == ic && ej == ib)
                                            edge = &var.edge_vec[j*NDIMS];
                                    }
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 13:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    for (int d=0; d<NDIMS-1; d++)
                        v[d] = var.vbc_values[ib] * fac * n[d];
                    v[NDIMS-1] = 0;
                    break;
                }
            }
        }

        //
        // Z, must be dealt last
        //

        // fast path: vz is usually free in the models
        if (bc_z0==0 && bc_z1==0) continue;

        if (flag & BOUNDZ0) {
            switch (bc_z0) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc_vz0;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc_vz0;
                break;
            }
        }
        if (flag & BOUNDZ1) {
            switch (bc_z1) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc_vz1;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc_vz1;
                break;
            }
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void apply_stress_bcs(const Param& param, const Variables& var, array_t& force)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    if (param.control.gravity == 0) return;
    //
    // Gravity-induced (hydrostatic and lithostatic) stress BCs
    //

#ifndef ACC
    #pragma omp parallel for
#endif
    #pragma acc parallel loop async
    for (int e=0; e<var.nelem; ++e)
        (*var.etmp_int)[e] = -1;

    for (int i=0; i<nbdrytypes; i++) {
        if (var.vbc_types[i] != 0 &&
            var.vbc_types[i] != 2 &&
            var.vbc_types[i] != 4) continue;

        if (i==iboundz0 && !param.bc.has_winkler_foundation) continue;
        if (i==iboundz1 && !param.bc.has_water_loading) continue;

        int bound, nbdry_nodes;

        #pragma acc kernels async
        {
            bound = static_cast<int>(var.bfacets[i]->size());
            nbdry_nodes = static_cast<int>(var.bnodes[i]->size());
        }

#ifndef ACC
        #pragma omp parallel default(none) \
            shared(param, var, force, i, NODE_OF_FACET, bound, nbdry_nodes)
#endif
        {
            // loops over all bdry facets
#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop async
            for (int n=0; n<bound; ++n) {
                // this facet belongs to element e
                int e = (*var.bfacets[i])[n].first;
                // this facet is the f-th facet of e
                int f = (*var.bfacets[i])[n].second;
                const int *conn = (*var.connectivity)[e];

                // the outward-normal vector
                double normal[NDIMS];
                // the z-coordinate of the facet center
                double zcenter;

                normal_vector_of_facet(f, conn, *var.coord, normal, zcenter);

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
                        const double sea_water_density = 1030;
                        p = sea_water_density * param.control.gravity * (param.control.surf_base_level - zcenter);
                    }
                }
                else {
                    // sidewalls
                    p = ref_pressure(param, zcenter);
                }

                (*var.etmp_int)[e] = n;
                double *tmp = (*var.tmp_result)[n];
                // lithostatc support - Archimed force (normal to the surface)
                for (int j=0; j<NODES_PER_FACET; ++j) {
                    int nn = conn[NODE_OF_FACET[f][j]];
                    double *f = force[nn];
                    for (int d=0; d<NDIMS; ++d) {
                            tmp[j*NDIMS + d] = p * normal[d] / NODES_PER_FACET;
                    }
                }
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop async
            for (int j=0; j<nbdry_nodes; ++j) {
                const int n = (*var.bnodes[i])[j];
                const int_vec& sup = (*var.support)[n];
                for (int k=0; k<sup.size(); ++k) {
                    int e = sup[k];
                    int ibound = (*var.etmp_int)[e];
                    if (ibound < 0) continue;  // not a boundary element

                    int f = (*var.bfacets[i])[ibound].second;  // facet index
                    const int *conn = (*var.connectivity)[e];
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
            #pragma acc parallel loop async
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
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void apply_stress_bcs_neumann(const Param& param, const Variables& var, array_t& force)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma wait // here is not ACC parallelized

    // Apply general stress (Neumann) boundary conditions
    for (int i = 0; i < 6; ++i) {
    // Check if there is a stress BC applied on this boundary
    if (var.stress_bc_types[i] == 0) continue;

    const auto& bdry = *(var.bfacets[i]);
    const int bound = static_cast<int>(bdry.size());

    // Loop over all boundary facets
    for (int n = 0; n < bound; ++n) {
        int e = bdry[n].first;      // This facet belongs to element e
        int f = bdry[n].second;     // This facet is the f-th facet of element e
        const int *conn = (*var.connectivity)[e];  // Element connectivity

        // Outward-normal vector and z-center for the facet
        double normal[NDIMS] = {0.0}; 
        double zcenter = 0.0;
        normal_vector_of_facet(f, conn, *var.coord, normal, zcenter);

        double traction[NDIMS] = {0.0};

#ifdef THREED
        // 3D case: Apply stress based on boundary type (x0, x1, y0, y1, z0, z1)
        switch (var.stress_bc_types[i]) {
            case 1: // Traction in the x-direction
                traction[0] = var.stress_bc_values[i]; // Set x-direction traction value
                break;
            case 2: // Traction in the y-direction
                traction[1] = var.stress_bc_values[i]; // Set y-direction traction value
                break;
            case 3: // Traction in the z-direction
                traction[2] = var.stress_bc_values[i]; // Set z-direction traction value
                break;
            default:
                continue; // Skip if no valid stress BC type
        }
#else
        // 2D case: Apply stress based on boundary type (x0, x1, z0, z1)
        switch (var.stress_bc_types[i]) {
            case 1: // Traction in the x-direction
                traction[0] = var.stress_bc_values[i]; // Set x-direction traction value
                break;
            case 3: // Traction in the z-direction (mapped to y in 2D)
                traction[1] = var.stress_bc_values[i]; // Set y-direction traction value in 2D
                break;
            default:
                continue; // Skip if no valid stress BC type
        }
#endif

        // Apply the traction to the boundary facet nodes
        for (int j = 0; j < NODES_PER_FACET; ++j) {
            int nn = conn[NODE_OF_FACET[f][j]];  // Get global node index for this facet
            double *force_node = force[nn];      // Pointer to the force array for this node

            // Distribute the traction (stress) to the nodes on the facet
            for (int d = 0; d < NDIMS; ++d) {
                // if(d == 1) std::cout<< force_node[d]<<","<< d <<","<< traction[d]<<std::endl;
                force_node[d] += traction[d] * normal[d] / NODES_PER_FACET; 
            }
        }
    }
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

namespace {

    void simple_diffusion(const Variables& var)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
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
            #pragma acc parallel loop async
            for (int i=0; i<var.nnode; i++) {
                total_dx[i] = 0.;
                total_slope[i] = 0.;
            }

#ifndef ACC
            #pragma omp for
#endif
            #pragma acc parallel loop async
            for (int i=0; i<var.surfinfo.etop; ++i) {
#ifdef THREED
                // this facet belongs to element e
                int e = top[i].first;
                // this facet is the f-th facet of e
                int f = top[i].second;

                const int *conn = (*var.connectivity)[e];
                int n0 = (*var.connectivity)[e][NODE_OF_FACET[f][0]];
                int n1 = (*var.connectivity)[e][NODE_OF_FACET[f][1]];
                int n2 = (*var.connectivity)[e][NODE_OF_FACET[f][2]];

                double projected_area;
                {
                    // normal vector of this facet
                    double normal[NDIMS];

                    // two vectors n0-n1 and n0-n2
                    // n is the cross product of these two vectors
                    // the length of n is 2 * triangle area
                    double x01, y01, z01, x02, y02, z02;
                    x01 = coord[n1][0] - coord[n0][0];
                    y01 = coord[n1][1] - coord[n0][1];
                    z01 = coord[n1][2] - coord[n0][2];
                    x02 = coord[n2][0] - coord[n0][0];
                    y02 = coord[n2][1] - coord[n0][1];
                    z02 = coord[n2][2] - coord[n0][2];

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
                shp2dx[0] = iv * (coord[n1][1] - coord[n2][1]);
                shp2dx[1] = iv * (coord[n2][1] - coord[n0][1]);
                shp2dx[2] = iv * (coord[n0][1] - coord[n1][1]);
                shp2dy[0] = iv * (coord[n2][0] - coord[n1][0]);
                shp2dy[1] = iv * (coord[n0][0] - coord[n2][0]);
                shp2dy[2] = iv * (coord[n1][0] - coord[n0][0]);

                double D[NODES_PER_FACET][NODES_PER_FACET];
                for (int j=0; j<NODES_PER_FACET; j++) {
                    for (int k=0; k<NODES_PER_FACET; k++) {
                        D[j][k] = (shp2dx[j] * shp2dx[k] +
                                shp2dy[j] * shp2dy[k]);
                    }
                }

                const int n[NODES_PER_FACET] = {n0, n1, n2};
                for (int j=0; j<NODES_PER_FACET; j++) {
                    double slope = 0;
                    for (int k=0; k<NODES_PER_FACET; k++)
                        slope += D[j][k] * coord[n[k]][2];
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
            #pragma acc parallel loop async
            for (int i=0; i<ntop; ++i) {
                int n = top_nodes[i];
#ifdef THREED
                int_vec &surf_sup = (*var.surfinfo.node_and_elems)[i];

                for (int j=0; j<surf_sup.size(); ++j) {
                    int k = surf_sup[j];
                    total_dx[n] += (*var.etmp)[k];

                    // this facet belongs to element e
                    int e = top[k].first;
                    // this facet is the f-th facet of e
                    int f = top[k].second;

                    const int *conn = (*var.connectivity)[e];
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
            #pragma acc parallel loop async
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
#ifdef USE_NPROF
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

    void get_surface_info(const Variables& var, \
        double_vec& top_base, double_vec& top_depth) {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
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
#ifdef USE_NPROF
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
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
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

#ifdef USE_NPROF
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
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
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

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void terrigenous_deposition(const Param& param,const Variables& var) {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
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
#ifdef USE_NPROF
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

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }
}


void surface_plstrain_diffusion(const Param &param, \
    const Variables& var, double_vec& plstrain)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double half_life = 1.e2 * YEAR2SEC;
    double lambha = 0.69314718056 / half_life; // ln2
    // #pragma omp parallel for default(none)      \
    //     shared(param, var, plstrain, lambha)
    for (auto e=(*var.top_elems).begin();e<(*var.top_elems).end();e++) {
        // Find the most abundant marker mattype in this element
        int_vec &a = (*var.elemmarkers)[*e];
        int mat = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
        if (mat != param.mat.mattype_oceanic_crust)
            plstrain[*e] -= plstrain[*e] * lambha * var.dt;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void correct_surface_element(const Variables& var, double_vec& volume, double_vec& volume_n, \
    tensor_t& stress, tensor_t& strain, tensor_t& strain_rate, double_vec& plstrain)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
#ifndef ACC
    #pragma omp parallel default(none) \
        shared(var, volume, volume_n, stress, strain, strain_rate, plstrain)
#endif
    {
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop async
        for (int i=0; i<var.ntop_elems; i++) {
            const double *coord1[NODES_PER_ELEM];

            const int e = (*var.top_elems)[i];

            for (int j=0; j<NODES_PER_ELEM;j++)
                coord1[j] = (*var.coord)[(*var.connectivity)[e][j]];
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
        #pragma acc parallel loop async
        for (int n=0;n<var.surfinfo.ntop;n++) {
            int nt = (*var.surfinfo.top_nodes)[n];
            int_vec &sup = (*var.support)[nt];
            volume_n[nt] = 0.;
            for (size_t i=0;i<sup.size();i++)
                volume_n[nt] += volume[sup[i]];
        }
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void surface_processes(const Param& param, const Variables& var, array_t& coord, tensor_t& stress, tensor_t& strain, \
                       tensor_t& strain_rate, double_vec& plstrain, double_vec& volume, double_vec& volume_n, \
                       SurfaceInfo& surfinfo, std::vector<MarkerSet*> &markersets, \
                       int_vec2D& elemmarkers, int_vec2D& markers_in_elem)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var)
#endif
    #pragma acc parallel loop async
    for (int i=0;i<var.surfinfo.ntop;i++)
        (*var.surfinfo.dh)[i] = 0.;

    switch (param.control.surface_process_option) {
    case 0:
        // no surface process
        break;
    case 1:
        simple_diffusion(var);
        break;
    case 101:
        custom_surface_processes(var, coord);
        break;
    case 102:
#ifdef THREED
        std::cout << "3D deposition of sediment processes is not ready yet.";
        exit(168);
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
        std::exit(1);
    }

#ifndef ACC
    #pragma omp parallel default(none) shared(param, var, coord)
#endif
    {
        // go through all surface nodes and abject all surface node by dh
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop async
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
        #pragma acc parallel loop async
        for (int i=0;i<var.surfinfo.etop;i++) {
            int_vec n(NDIMS);
            double dh_e = 0.;

            for (int j=0;j<NDIMS;j++) {
                int k = (*var.surfinfo.elem_and_nodes)[i][j];
                n[j] = (*var.surfinfo.top_nodes)[k];
                dh_e += (*var.surfinfo.dh)[k];
            }

#ifdef THREED
            double base = (((*var.coord)[n[1]][0] - (*var.coord)[n[0]][0])*((*var.coord)[n[2]][1] - (*var.coord)[n[0]][1]) \
                - ((*var.coord)[n[2]][0] - (*var.coord)[n[0]][0])*((*var.coord)[n[1]][1] - (*var.coord)[n[0]][1]))/2.0;
#else
            double base = (*var.coord)[n[0]][0] - (*var.coord)[n[1]][0];
#endif
            (*var.surfinfo.edvacc_surf)[i] += dh_e * base / NDIMS;
        }
    }

//     if (var.steps%1000 == 0) {
//         double sum = 0.;
// #ifndef ACC
//         #pragma omp parallel default(none) shared(var) reduction(+:sum)
// #endif
//         #pragma acc parallel loop async reduction(+:sum)
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
    #pragma acc parallel loop reduction(max:maxdh) async
    for (int i=0; i<var.surfinfo.ntop; ++i) {
        double tmp = fabs((*var.surfinfo.dh)[i]);

        if (maxdh < tmp)
            maxdh = tmp;
    }

    #pragma acc wait

    surfinfo.max_surf_vel = maxdh / var.dt;

    correct_surface_element(var, volume, volume_n, stress, strain, strain_rate, plstrain);

    if (var.steps != 0) {
        if ( var.steps % param.mesh.quality_check_step_interval == 0) {
            markersets[0]->correct_surface_marker(param, var, *var.surfinfo.dhacc, elemmarkers, markers_in_elem);
            std::fill(var.surfinfo.dhacc->begin(), var.surfinfo.dhacc->end(), 0.);

            // set marker of sediment.
            markersets[0]->set_surface_marker(param, var, param.mesh.smallest_size, param.mat.mattype_sed, *var.surfinfo.edvacc_surf, elemmarkers, markers_in_elem);
        }
    }

#ifdef THREED
    // todo
#else
    if (!(var.steps % param.mesh.quality_check_step_interval &&  var.steps != 0))
        // correct the plastic strain of urface element for preventing surface landslide.
        surface_plstrain_diffusion(param,var,plstrain);

    if ( param.control.is_reporting_terrigenous_info && var.steps%10000 == 0 &&  var.steps != 0  ) {
        double max_dh = 0., min_dh = 0., max_dh_oc = 0.;
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

#ifdef USE_NPROF
    nvtxRangePop();
#endif


}
