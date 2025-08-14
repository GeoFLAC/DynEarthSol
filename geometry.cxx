#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include <cmath>
#include <limits>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "geometry.hpp"
#include "bc.hpp"
#include "mesh.hpp"


/* Given two points, returns the distance^2 */
double dist2(const double* a, const double* b)
{
    double sum = 0;
    for (int i=0; i<NDIMS; ++i) {
        double d = b[i] - a[i];
        sum += d * d;
    }
    return sum;
}


/* Given four 3D points, returns the (signed) volume of the enclosed
   tetrahedron */
#pragma acc routine seq
static double tetrahedron_volume(const double *d0,
                                 const double *d1,
                                 const double *d2,
                                 const double *d3)
{
    double x01 = d0[0] - d1[0];
    double x12 = d1[0] - d2[0];
    double x23 = d2[0] - d3[0];

    double y01 = d0[1] - d1[1];
    double y12 = d1[1] - d2[1];
    double y23 = d2[1] - d3[1];

    double z01 = d0[2] - d1[2];
    double z12 = d1[2] - d2[2];
    double z23 = d2[2] - d3[2];

    return (x01*(y23*z12 - y12*z23) +
            x12*(y01*z23 - y23*z01) +
            x23*(y12*z01 - y01*z12)) / 6;
}


/* Given two points, returns the area of the enclosed triangle */
#pragma acc routine seq
static double triangle_area(const double *a,
                            const double *b,
                            const double *c)
{
    double ab0, ab1, ac0, ac1;

    // ab: vector from a to b
    ab0 = b[0] - a[0];
    ab1 = b[1] - a[1];
    // ac: vector from a to c
    ac0 = c[0] - a[0];
    ac1 = c[1] - a[1];

#ifndef THREED
    // area = norm(cross product of ab and ac) / 2
    return std::fabs(ab0*ac1 - ab1*ac0) / 2;
#else
    double ab2, ac2;
    ab2 = b[2] - a[2];
    ac2 = c[2] - a[2];

    // vector components of ab x ac
    double d0, d1, d2;
    d0 = ab1*ac2 - ab2*ac1;
    d1 = ab2*ac0 - ab0*ac2;
    d2 = ab0*ac1 - ab1*ac0;

    // area = norm(cross product of ab and ac) / 2
    return std::sqrt(d0*d0 + d1*d1 + d2*d2) / 2;
#endif
}

double compute_volume(const double **coord)
{
    const double *a = coord[0];
    const double *b = coord[1];
    const double *c = coord[2];
#ifdef THREED
    const double *d = coord[3];
    return tetrahedron_volume(a, b, c, d);
#else
    return triangle_area(a, b, c);
#endif
}

void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(coord, connectivity, volume)
#endif
    #pragma acc parallel loop 
    for (int e=0; e<volume.size(); ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        const double *d = coord[n3];
        volume[e] = tetrahedron_volume(a, b, c, d);
#else
        volume[e] = triangle_area(a, b, c);
#endif
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void compute_volume(const Variables &var,
                    double_vec &volume)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, volume)
#endif
    #pragma acc parallel loop 
    for (int e=0; e<var.nelem; ++e) {
        int n0 = (*var.connectivity)[e][0];
        int n1 = (*var.connectivity)[e][1];
        int n2 = (*var.connectivity)[e][2];

        const double *a = (*var.coord)[n0];
        const double *b = (*var.coord)[n1];
        const double *c = (*var.coord)[n2];

#ifdef THREED
        int n3 = (*var.connectivity)[e][3];
        const double *d = (*var.coord)[n3];
        volume[e] = tetrahedron_volume(a, b, c, d);
#else
        volume[e] = triangle_area(a, b, c);
#endif
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void compute_dvoldt(const Variables &var, double_vec &dvoldt, double_vec &etmp)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    /* dvoldt is the volumetric strain rate, weighted by the element volume,
     * lumped onto the nodes.
     */
//    std::fill_n(dvoldt.begin(), var.nnode, 0);

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(var, etmp)
#endif
    #pragma acc parallel loop
    for (int e=0;e<var.nelem;e++) {
        const int *conn = (*var.connectivity)[e];
        const double *strain_rate= (*var.strain_rate)[e];
        // TODO: try another definition:
        // dj = (volume[e] - volume_old[e]) / volume_old[e] / dt
        double dj = trace(strain_rate);
        etmp[e] = dj * (*var.volume)[e];
    }

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(var,dvoldt,etmp)
#endif
    #pragma acc parallel loop
    for (int n=0;n<var.nnode;n++) {
        dvoldt[n] = 0.;
        for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e)
	        dvoldt[n] += etmp[*e];
        dvoldt[n] /= (*var.volume_n)[n];
    }

    // std::cout << "dvoldt:\n";
    // print(std::cout, dvoldt);
    // std::cout << "\n";
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void compute_edvoldt(const Variables &var, double_vec &dvoldt,
                     double_vec &edvoldt)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    /* edvoldt is the averaged (i.e. smoothed) dvoldt on the element.
     * It is used in update_stress() to prevent mesh locking.
     */

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, edvoldt)
#endif
    #pragma acc parallel loop
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double dj = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dj += dvoldt[n];
        }
        edvoldt[e] = dj / NODES_PER_ELEM;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
    // std::cout << "edvoldt:\n";
    // print(std::cout, edvoldt);
    // std::cout << "\n";
}


void NMD_stress(const Param& param, const Variables &var,
    double_vec &dp_nd, tensor_t& stress, double_vec &etmp)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // dp_nd is the pressure change, weighted by the element volume,
    // lumped onto the nodes.

//    double **centroid = elem_center(*var.coord, *var.connectivity); // centroid of elements
/*
    // weight with inverse distance
    if(false) {
        #pragma omp parallel for default(none) shared(var,centroid,tmp_result)
        for (int e=0;e<var.nelem;e++) {
            const int *conn = (*var.connectivity)[e];
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                const double *d = (*var.coord)[conn[i]];
                tmp_result[i][e] = 1. / sqrt( dist2(d, centroid[e])  );
                tmp_result[i + NODES_PER_ELEM][e] = tmp_result[i][e] * (*var.dpressure)[e];
            }
        }

        #pragma omp parallel for default(none) shared(var,dp_nd,tmp_result)
        for (int n=0;n<var.nnode;n++) {
            double dist_inv_sum = 0.;
            for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e) {
                const int *conn = (*var.connectivity)[*e];
                for (int i=0;i<NODES_PER_ELEM;i++) {
                    if (n == conn[i]) {
                        dist_inv_sum += tmp_result[ i ][*e];
                        dp_nd[n] += tmp_result[i + NODES_PER_ELEM][*e];
                        break;
                    }
                }
            }
            dp_nd[n] /= dist_inv_sum;
        }

    // weight with volumn
    } else {
        */

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,etmp)
#endif
    #pragma acc parallel loop
    for (int e=0;e<var.nelem;e++) {
        const int *conn = (*var.connectivity)[e];
        double dp = (*var.dpressure)[e];
        etmp[e] = dp * (*var.volume)[e];
    }

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,dp_nd,etmp)
#endif
    #pragma acc parallel loop
    for (int n=0;n<var.nnode;n++) {
        dp_nd[n] = 0;
        for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e)
            dp_nd[n] += etmp[*e];
        dp_nd[n] /= (*var.volume_n)[n];
    }
//    }

    /* dp_el is the averaged (i.e. smoothed) dp_nd on the element.
     */
#ifndef ACC
    #pragma omp parallel for default(none) shared(param, var, dp_nd, stress)
#endif
    #pragma acc parallel loop
    for (int e=0; e<var.nelem; ++e) {

        double factor;
        switch (param.mat.rheol_type) {
        case MatProps::rh_viscous:
        case MatProps::rh_maxwell:
        case MatProps::rh_evp:
            if ((*var.viscosity)[e] < param.control.mixed_stress_reference_viscosity)
                factor = 0.;
            else
                factor = std::min((*var.viscosity)[e] / (param.control.mixed_stress_reference_viscosity * 10.), 1.);
            break;
        default:
            factor = 1;
        }

        const int *conn = (*var.connectivity)[e];
        double dp = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dp += dp_nd[n];
        }
        double dp_el = dp / NODES_PER_ELEM;

    	double* s = stress[e];
	    double dp_orig = (*var.dpressure)[e];
        double ddp = ( - dp_orig + dp_el ) / NDIMS * factor;
	    for (int i=0; i<NDIMS; ++i)
            s[i] += ddp;
    }

//    delete [] centroid[0];
//    delete [] centroid;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

double compute_dt(const Param& param, Variables& var)

{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // constant dt
    if (param.control.fixed_dt != 0) return param.control.fixed_dt;

    // dynamic dt
    double dt_maxwell = std::numeric_limits<double>::max();
    double dt_diffusion = std::numeric_limits<double>::max();
    double dt_hydro_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();
    int mattype_ref = param.mat.mattype_mantle;

    // Define element velocity arrays
    double_vec velocity_x_element(var.nelem, 0.0);
    double_vec velocity_y_element(var.nelem, 0.0);
    double_vec velocity_z_element(var.nelem, 0.0); // Used only for 3D
    double vx_element = 0.0, vy_element = 0.0, vz_element = 0.0;
    // Global max velocity for elements
    // double global_max_vem = std::max(std::max(var.max_vbc_val, 0.0), 1E-20);
    double global_max_vem = 0.0;  // 
    double global_dt_min = std::numeric_limits<double>::max(); // based on length and S wave velocity

#ifndef ACC
    #pragma omp parallel for reduction(min:minl, dt_maxwell, dt_diffusion, dt_hydro_diffusion, global_dt_min) \
        reduction(max: global_max_vem) \
        default(none) shared(param, var, velocity_x_element, velocity_y_element, velocity_z_element) \
        private(vx_element, vy_element, vz_element)
#endif
    #pragma acc parallel loop reduction(min:minl, dt_maxwell, dt_diffusion, dt_hydro_diffusion, global_dt_min) \
        reduction(max: global_max_vem)

    for (int e=0; e<var.nelem; ++e) {

        vx_element = 0.0, vy_element = 0.0, vz_element = 0.0;
        int n0 = (*var.connectivity)[e][0];
        int n1 = (*var.connectivity)[e][1];
        int n2 = (*var.connectivity)[e][2];

        // calculate maxium velocity in the element
        const array_t& vel = *var.vel;
        const int *conn = (*var.connectivity)[e];
        double weight = 1.0 / NODES_PER_ELEM; 

        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            vx_element += vel[conn[j]][0] * weight;
            vy_element += vel[conn[j]][1] * weight;
            #ifdef THREED
            vz_element += vel[conn[j]][2] * weight;
            #endif
        }

        velocity_x_element[e] = vx_element;
        velocity_y_element[e] = vy_element;
        #ifdef THREED
        velocity_z_element[e] = vz_element;
        #endif

        // min height of this element
        #ifdef THREED  
        double max_vem = std::sqrt(vx_element*vx_element + vy_element*vy_element + vz_element*vz_element);
        #else
        double max_vem = std::sqrt(vx_element*vx_element + vy_element*vy_element);
        #endif

        // Find global max velocity
        global_max_vem = std::max(global_max_vem, max_vem);

        // std::cout<< "Element: " << e << " max_vem: " << global_max_vem << std::scientific << std::setprecision(5) << std::endl;

        const double *a = (*var.coord)[n0];
        const double *b = (*var.coord)[n1];
        const double *c = (*var.coord)[n2];

        // min height of this element
        double minh;
#ifdef THREED
        {
            int n3 = (*var.connectivity)[e][3];
            const double *d = (*var.coord)[n3];

            // max facet area of this tet
            double maxa = std::max(std::max(triangle_area(a, b, c),
                                            triangle_area(a, b, d)),
                                   std::max(triangle_area(c, d, a),
                                            triangle_area(c, d, b)));
            minh = 3 * (*var.volume)[e] / maxa;
        }
#else
        {
            // max edge length of this triangle
            double maxl = std::sqrt(std::max(std::max(dist2(a, b),
                                                      dist2(b, c)),
                                             dist2(a, c)));
            minh = 2 * (*var.volume)[e] / maxl;
        }
#endif
        dt_maxwell = std::min(dt_maxwell,
                              0.5 * var.mat->visc_min / (1e-40 + var.mat->shearm(e)));
        if (param.control.has_thermal_diffusion)
            dt_diffusion = std::min(dt_diffusion,
                                    0.5 * minh * minh / var.mat->therm_diff_max);
        
        // Compute dt_hydro_diffusion (hydraulic)
        if (param.control.has_hydraulic_diffusion)
            if (var.mat->hydro_diff_max > 0) {
                dt_hydro_diffusion = std::min(dt_hydro_diffusion,
                                            0.5 * minh * minh / var.mat->hydro_diff_max);
            }
        minl = std::min(minl, minh);

        // Find global min delta t to meet CFL condition
        global_dt_min = std::min(global_dt_min, minl/std::sqrt(var.mat->shearm(e)/var.mat->rho(e)) /5.0);
    }

    double max_vbc_val;
    if (param.control.characteristic_speed == 0) {
        max_vbc_val = var.max_vbc_val;
        if (param.control.surface_process_option > 0)
            max_vbc_val = std::max(max_vbc_val, var.surfinfo.max_surf_vel*5e-1);
    }
    else
        max_vbc_val = param.control.characteristic_speed;

    double dt_advection = 0;
    double dt_elastic = 0;
    global_max_vem = std::max(global_max_vem, var.max_vbc_val);
    var.max_global_vel_mag = global_max_vem;
    var.global_dt_min = global_dt_min;
    
    {
        double max_vel = 0;
        #pragma omp parallel for reduction(max:max_vel) default(none) shared(var)
        for (auto n=var.surfinfo.top_nodes->begin(); n<var.surfinfo.top_nodes->end(); ++n) {
            double vel = 0;
            for (int i=0; i<NDIMS; ++i)
                vel += (*var.vel)[*n][i] * (*var.vel)[*n][i];
            max_vel = std::max(max_vel, vel);
        }
        max_vbc_val = std::max(max_vbc_val, std::sqrt(max_vel));
    }

    // Calculate dt_advection and dt_elastic 
    if(param.control.has_ATS)
    {
        dt_advection = 0.5 * minl / global_max_vem;
        dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (global_max_vem * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(param.mat.bulk_modulus[0] / param.mat.rho0[0]);
        dt_elastic = std::max(dt_elastic, global_dt_min); // Ensure dt_elastic is not smaller than global_dt_min
    }
    else
    {
        dt_advection = 0.5 * minl / max_vbc_val;
        dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (max_vbc_val * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(param.mat.bulk_modulus[mattype_ref] / param.mat.rho0[mattype_ref]);
    }

    // Combine dt calculations and incorporate dt_hydro_diffusion
    double dt = std::min({dt_elastic, dt_maxwell, dt_advection, dt_diffusion, dt_hydro_diffusion}) * param.control.dt_fraction;
    // double dt = std::min({dt_elastic, dt_maxwell, dt_advection, dt_diffusion}) * param.control.dt_fraction;
    if (param.debug.dt) {
        std::cout << "step #" << var.steps << "  dt: " << dt_maxwell << " " << dt_diffusion << " " 
                  << dt_hydro_diffusion << " " << dt_advection << " " << dt_elastic << " sec\n";
    }
    if (dt <= 0) {
        std::cerr << "Error: dt <= 0!  " << dt_maxwell << " " << dt_diffusion
                  << " " << dt_hydro_diffusion << " " << dt_advection << " " << dt_elastic << "\n";
        std::exit(11);
    }
    
#ifdef USE_NPROF
    nvtxRangePop();
#endif
    return dt;
}

double compute_dt_PT(const Param& param, const Variables& var)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // constant dt
    if (param.control.fixed_dt != 0) return param.control.fixed_dt;

    // dynamic dt
    double dt_maxwell = std::numeric_limits<double>::max();
    double dt_diffusion = std::numeric_limits<double>::max();
    double dt_hydro_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();
    int mattype_ref = param.mat.mattype_mantle;

    #pragma omp parallel for reduction(min:minl,dt_maxwell,dt_diffusion,dt_hydro_diffusion)    \
        default(none) shared(param,var)
    // #pragma acc parallel loop reduction(min:minl, dt_maxwell, dt_diffusion,dt_hydro_diffusion)
    for (int e=0; e<var.nelem; ++e) {
        int n0 = (*var.connectivity)[e][0];
        int n1 = (*var.connectivity)[e][1];
        int n2 = (*var.connectivity)[e][2];

        const double *a = (*var.coord)[n0];
        const double *b = (*var.coord)[n1];
        const double *c = (*var.coord)[n2];

        // min height of this element
        double minh;
#ifdef THREED
        {
            int n3 = (*var.connectivity)[e][3];
            const double *d = (*var.coord)[n3];

            // max facet area of this tet
            double maxa = std::max(std::max(triangle_area(a, b, c),
                                            triangle_area(a, b, d)),
                                   std::max(triangle_area(c, d, a),
                                            triangle_area(c, d, b)));
            minh = 3 * (*var.volume)[e] / maxa;
        }
#else
        {
            // max edge length of this triangle
            double maxl = std::sqrt(std::max(std::max(dist2(a, b),
                                                      dist2(b, c)),
                                             dist2(a, c)));
            minh = 2 * (*var.volume)[e] / maxl;
        }
#endif
        dt_maxwell = std::min(dt_maxwell,
                              0.5 * var.mat->visc_min / (1e-40 + var.mat->shearm(e)));
        // if (param.control.has_thermal_diffusion)
        //     dt_diffusion = std::min(dt_diffusion,
        //                             0.5 * minh * minh / var.mat->therm_diff_max);
        
        // // Compute dt_hydro_diffusion (hydraulic)
        // if (var.mat->hydro_diff_max > 0) {
        //     dt_hydro_diffusion = std::min(dt_hydro_diffusion,
        //                                   0.5 * minh * minh / var.mat->hydro_diff_max);
        // }
        minl = std::min(minl, minh);
    }


    // max_vbc_val is maximum boundary velocity
    double max_vbc_val;
    if (param.control.characteristic_speed == 0) {
        max_vbc_val = var.max_vbc_val; 

        if (param.control.surface_process_option > 0)
            max_vbc_val = std::max(max_vbc_val, var.surfinfo.max_surf_vel*5e-1);
    }
    else
        max_vbc_val = param.control.characteristic_speed;

    double dt_advection = 0.5 * minl / max_vbc_val;
    double dt_elastic = (param.control.is_quasi_static) ?
        0.5 * minl / (max_vbc_val * param.control.inertial_scaling) :
        0.5 * minl / std::sqrt(param.mat.bulk_modulus[mattype_ref] / param.mat.rho0[mattype_ref]);

    double dt = std::min({dt_elastic, dt_maxwell, dt_advection}) * param.control.dt_fraction;
    if (param.debug.dt) {
        std::cout << "step #" << var.steps << "  dt: " << dt_maxwell << " " << dt_advection << " " << dt_elastic << " sec\n";
    }
    if (dt <= 0) {
        std::cerr << "Error: dt <= 0!  " << dt_maxwell << " "  << dt_advection << " " << dt_elastic << "\n";
        std::exit(11);
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
    return dt;
}

void compute_mass(const Param &param, const Variables &var,
                  double max_vbc_val, double_vec &volume_n,
                  double_vec &mass, double_vec &tmass, double_vec &hmass, double_vec &ymass, elem_cache &tmp_result)

{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // volume_n is (node-averaged volume * NODES_PER_ELEM)
    // volume_n.assign(volume_n.size(), 0);
    // mass.assign(mass.size(), 0);
    // tmass.assign(tmass.size(), 0);

    const double pseudo_speed = max_vbc_val * param.control.inertial_scaling; // for non-ATP using max velocity on boundary
    const double pseudo_speed_ATP = var.max_global_vel_mag * param.control.inertial_scaling; // for ATP using global max velocity

    double diff_e;

    #pragma acc serial
    {
        // Retrieve hydraulic properties for the element
        double perm_e = var.mat->perm(0);                // Intrinsic permeability 
        double mu_e = var.mat->mu_fluid(0);              // Fluid dynamic viscosity
        double alpha_b = var.mat->alpha_biot(0);         // Biot coefficient
        double rho_f = var.mat->rho_fluid(0);            // Fluid density
        double phi_e = var.mat->phi(0);        // Element porosity
        double comp_fluid = var.mat->beta_fluid(0);        // fluid comporessibility
        double bulkm = var.mat->bulkm(0);
        double shearm = var.mat->shearm(0);
        double matrix_comp = 1.0 / (bulkm +4.0*shearm/3.0);

        rho_f = 1000.0; 
        double gamma_w = rho_f * param.control.gravity; // specific weight
        
        // Hydraulic conductivity using permeability and viscosity
        double hydraulic_conductivity = perm_e * gamma_w / mu_e;
        
        // Compute element diffusivity and update max using reduction
        diff_e = hydraulic_conductivity / (phi_e * comp_fluid + alpha_b * matrix_comp) / gamma_w;
    }

    if (pseudo_speed < diff_e && param.control.has_hydraulic_diffusion)
    {
        std::cout << "pseudo speed is too slow, increase mass scaling" << std::endl;
        std::exit(11);
    }

#ifndef ACC
#ifdef GPP1X
    #pragma omp parallel default(none) shared(var, param, volume_n, mass, tmass, hmass, ymass, pseudo_speed, pseudo_speed_ATP, tmp_result)
#else
    #pragma omp parallel default(none) shared(var, param, volume_n, mass, tmass, hmass, ymass, tmp_result)
#endif
#endif
    {
#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop
        for (int e=0;e<var.nelem;e++) {
            double *tr = tmp_result[e];
            double rho;

            if(param.control.has_ATS)
            {
                double apprent_speed = std::min(pseudo_speed_ATP, std::sqrt(var.mat->shearm(e)/var.mat->rho(e))); // minimum speed

                rho = (param.control.is_quasi_static) ?
                (*var.mat).bulkm(e) / (apprent_speed * apprent_speed) :  // pseudo density for quasi-static sim
                (*var.mat).rho(e);  // true density for dynamic sim

                if (param.control.has_hydraulic_diffusion && (param.control.is_quasi_static == false)) {
                    // Modified density considering porosity for hydraulic diffusion
                        rho = (*var.mat).rho(e) * (1 - (*var.mat).phi(e)) + 1000.0 * (*var.mat).phi(e);
                    }
            }
            else
            {    
                rho = (param.control.is_quasi_static) ?
                (*var.mat).bulkm(e) / (pseudo_speed * pseudo_speed) :  // pseudo density for quasi-static sim
                (*var.mat).rho(e);  // true density for dynamic sim

                if (param.control.has_hydraulic_diffusion && (param.control.is_quasi_static == false)) {
                    // Modified density considering porosity for hydraulic diffusion
                        rho = (*var.mat).rho(e) * (1 - (*var.mat).phi(e)) + 1000.0 * (*var.mat).phi(e);
                    }

            }

            double bulk_comp = 1.0/(*var.mat).bulkm(e); // lambda + 2G/3
            if(NDIMS == 2) bulk_comp = 1.0/((*var.mat).bulkm(e) + (*var.mat).shearm(e)/3.0); // lambda + G 

            double hm_coeff = (*var.mat).alpha_biot(e) + (*var.mat).phi(e) - (*var.mat).alpha_biot(e) * (*var.mat).phi(e); 
            double m = rho * (*var.volume)[e] / NODES_PER_ELEM;
            double tm = (*var.mat).rho(e) * (*var.mat).cp(e) * (*var.volume)[e] / NODES_PER_ELEM;
            double hm = (hm_coeff * bulk_comp + (*var.mat).phi(e) * (*var.mat).beta_fluid(e)) * (*var.volume)[e] / NODES_PER_ELEM;
            double ym = 9 * (*var.mat).bulkm(e) * (*var.mat).shearm(e) / (3 * (*var.mat).bulkm(e) + (*var.mat).shearm(e)) / NODES_PER_ELEM; // Young's modulus

            tr[0] = (*var.volume)[e];
            tr[1] = m;
            if (param.control.has_thermal_diffusion)
                tr[2] = tm;
            tr[3] = hm; // check
            tr[4] = ym; 
        }

#ifndef ACC
        #pragma omp for
#endif
        #pragma acc parallel loop
        for (int n=0;n<var.nnode;n++) {
            volume_n[n]=0;
            mass[n]=0;
            tmass[n]=0;
            hmass[n]=0;
            ymass[n]=0;
        
            for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e) {
                double *tr = tmp_result[*e];
                volume_n[n] += tr[0];
                mass[n] += tr[1];
                if (param.control.has_thermal_diffusion)
                    tmass[n] += tr[2];
                hmass[n] += tr[3];
                ymass[n] += tr[4];
            }
        }
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void compute_shape_fn(const Variables &var, shapefn &shpdx, shapefn &shpdy, shapefn &shpdz)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none)      \
        shared(var, shpdx, shpdy, shpdz)
#endif
    #pragma acc parallel loop
    for (int e=0;e<var.nelem;e++) {

        int n0 = (*var.connectivity)[e][0];
        int n1 = (*var.connectivity)[e][1];
        int n2 = (*var.connectivity)[e][2];

        const double *d0 = (*var.coord)[n0];
        const double *d1 = (*var.coord)[n1];
        const double *d2 = (*var.coord)[n2];

#ifdef THREED
        {
            int n3 = (*var.connectivity)[e][3];
            const double *d3 = (*var.coord)[n3];

            double iv = 1 / (6 * (*var.volume)[e]);

            double x01 = d0[0] - d1[0];
            double x02 = d0[0] - d2[0];
            double x03 = d0[0] - d3[0];
            double x12 = d1[0] - d2[0];
            double x13 = d1[0] - d3[0];
            double x23 = d2[0] - d3[0];

            double y01 = d0[1] - d1[1];
            double y02 = d0[1] - d2[1];
            double y03 = d0[1] - d3[1];
            double y12 = d1[1] - d2[1];
            double y13 = d1[1] - d3[1];
            double y23 = d2[1] - d3[1];

            double z01 = d0[2] - d1[2];
            double z02 = d0[2] - d2[2];
            double z03 = d0[2] - d3[2];
            double z12 = d1[2] - d2[2];
            double z13 = d1[2] - d3[2];
            double z23 = d2[2] - d3[2];

            shpdx[e][0] = iv * (y13*z12 - y12*z13);
            shpdx[e][1] = iv * (y02*z23 - y23*z02);
            shpdx[e][2] = iv * (y13*z03 - y03*z13);
            shpdx[e][3] = iv * (y01*z02 - y02*z01);

            shpdy[e][0] = iv * (z13*x12 - z12*x13);
            shpdy[e][1] = iv * (z02*x23 - z23*x02);
            shpdy[e][2] = iv * (z13*x03 - z03*x13);
            shpdy[e][3] = iv * (z01*x02 - z02*x01);

            shpdz[e][0] = iv * (x13*y12 - x12*y13);
            shpdz[e][1] = iv * (x02*y23 - x23*y02);
            shpdz[e][2] = iv * (x13*y03 - x03*y13);
            shpdz[e][3] = iv * (x01*y02 - x02*y01);
        }
#else
        {
            double iv = 1 / (2 * (*var.volume)[e]);

            shpdx[e][0] = iv * (d1[1] - d2[1]);
            shpdx[e][1] = iv * (d2[1] - d0[1]);
            shpdx[e][2] = iv * (d0[1] - d1[1]);

            shpdz[e][0] = iv * (d2[0] - d1[0]);
            shpdz[e][1] = iv * (d0[0] - d2[0]);
            shpdz[e][2] = iv * (d1[0] - d0[0]);
        }
#endif
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

double elem_quality(const array_t &coord, const conn_t &connectivity,
                    const double_vec &volume, int e)
{
    /* This function returns the quality (0~1) of the element.
     * The quality of an equidistant (i.e. best quality) tetrahedron/triangle is 1.
     */
    double quality;
    double vol = volume[e];
    int n0 = connectivity[e][0];
    int n1 = connectivity[e][1];
    int n2 = connectivity[e][2];

    const double *a = coord[n0];
    const double *b = coord[n1];
    const double *c = coord[n2];

#ifdef THREED
    {
        int n3 = connectivity[e][3];
        const double *d = coord[n3];
        double normalization_factor = 216 * std::sqrt(3);

        double area_sum = (triangle_area(a, b, c) +
                           triangle_area(a, b, d) +
                           triangle_area(c, d, a) +
                           triangle_area(c, d, b));
        quality = normalization_factor * vol * vol / (area_sum * area_sum * area_sum);
    }
#else
    {
        double normalization_factor = 4 * std::sqrt(3);

        double dist2_sum = dist2(a, b) + dist2(b, c) + dist2(a, c);
        quality = normalization_factor * vol / dist2_sum;
    }
#endif

    return quality;
}


double worst_elem_quality(const array_t &coord, const conn_t &connectivity,
                          const double_vec &volume, int &worst_elem)
{
    double q = 1;
    worst_elem = 0;
    for (std::size_t e=0; e<volume.size(); e++) {
        double quality = elem_quality(coord, connectivity, volume, e);
        if (quality < q) {
            q = quality;
            worst_elem = e;
        }
    }
    return q;
}


