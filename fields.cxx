#include <iostream>
#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "matprops.hpp"
#include "fields.hpp"
#include "utils.hpp"

void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);
    var.hmass = new double_vec(n);
    var.ymass = new double_vec(n);
    var.edvoldt = new double_vec(e);

//    var.marker_in_elem = new int_vec2D(e);

    {
        // these fields are reallocated during remeshing interpolation
        var.temperature = new double_vec(n);
        var.ppressure = new double_vec(n);
        var.dppressure = new double_vec(n);
        var.coord0 = new array_t(n);
        var.plstrain = new double_vec(e);
        var.delta_plstrain = new double_vec(e);
        var.vel = new array_t(n, 0);
        var.strain = new tensor_t(e, 0);
        var.stress = new tensor_t(e, 0);
        var.stressyy = new double_vec(e, 0);
        var.old_mean_stress = new double_vec(e, 0);
        // var.stress_old = new tensor_t(e, 0);
        var.radiogenic_source = new double_vec(e, 0);

//         // these fields are from Denis' code
//         var.MAX_shear = new double_vec(e);
//         var.MAX_shear_0 = new double_vec(e);
//         var.CL = new double_vec(e);
//         var.dl_min = new double_vec(e);
//         var.maxv = new double_vec(e);
//         var.strain_energy = new double_vec(e);
//         var.kinetic_energy = new double_vec(e);
// #ifdef RS
//         var.state1 = new double_vec(e);
//         var.slip_velocity = new double_vec(e);
//         var.friction_coefficient = new double_vec(e);
//         var.RS_shear = new double_vec(e);
//         var.Failure_mode = new double_vec(e);
// #endif
    }

    var.ntmp = new double_vec(n);
    var.dpressure = new double_vec(e, 0);

    var.viscosity = new double_vec(e,param.mat.visc_max);

    var.force = new array_t(n, 0);
    var.force_residual = new array_t(n, 0);

    var.strain_rate = new tensor_t(e, 0);

    var.shpdx = new shapefn(e);
    if (NDIMS == 3) var.shpdy = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);

    var.tmp_result = new elem_cache(e);
    var.tmp_result_sg = new double_vec(e);
}


void reallocate_variables(const Param& param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    delete var.volume;
    delete var.volume_old;
    delete var.volume_n;
    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    delete var.mass;
    delete var.tmass;
    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    delete var.hmass;
    var.hmass = new double_vec(n);
    delete var.ymass;
    var.ymass = new double_vec(n);
  
    delete var.edvoldt;
    var.edvoldt = new double_vec(e);

//    delete var.marker_in_elem;
//    var.marker_in_elem = new int_vec2D(e);

    delete var.ntmp;
    var.ntmp = new double_vec(n);
    delete var.dpressure;
    var.dpressure = new double_vec(e, 0);
    delete var.viscosity;
    var.viscosity = new double_vec(e,param.mat.visc_max);
    delete var.force;
    var.force = new array_t(n, 0);

    delete var.force_residual;
    var.force_residual = new array_t(n, 0);

    delete var.strain_rate;
    var.strain_rate = new tensor_t(e, 0);

    delete var.shpdx;
    delete var.shpdz;
    var.shpdx = new shapefn(e);
    if (NDIMS == 3) {
        delete var.shpdy;
        var.shpdy = new shapefn(e);
    }
    var.shpdz = new shapefn(e);

    delete var.mat;
    var.mat = new MatProps(param, var);

    delete var.tmp_result;
    var.tmp_result = new elem_cache(e);
    delete var.tmp_result_sg;
    var.tmp_result_sg = new double_vec(e);

}


void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot, elem_cache &tmp_result)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma omp parallel for default(none) shared(var,temperature,tmp_result)
    #pragma acc parallel loop
    for (int e=0;e<var.nelem;e++) {
        // diffusion matrix

        const int *conn = (*var.connectivity)[e];
        double *tr = tmp_result[e];
        double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volume
        double rh = (*var.radiogenic_source)[e] * (*var.volume)[e] * var.mat->rho(e) / NODES_PER_ELEM;
        const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
#endif
        const double *shpdz = (*var.shpdz)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double diffusion = 0.;
            for (int j=0; j<NODES_PER_ELEM; ++j) {
#ifdef THREED
                diffusion += (shpdx[i] * shpdx[j] + \
                            shpdy[i] * shpdy[j] + \
                            shpdz[i] * shpdz[j]) * temperature[conn[j]];
#else
                diffusion += (shpdx[i] * shpdx[j] + \
                            shpdz[i] * shpdz[j]) * temperature[conn[j]];
#endif
            }
            tr[i] = diffusion * kv - rh;
        }
    }

    #pragma omp parallel for default(none) \
        shared(param,var,tdot,temperature,tmp_result)
    #pragma acc parallel loop
    for (int n=0;n<var.nnode;n++) {
        tdot[n]=0;
        for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e) {
            const int *conn = (*var.connectivity)[*e];
            const double *tr = tmp_result[*e];
            for (int i=0;i<NODES_PER_ELEM;i++) {
                if (n == conn[i]) {
                    tdot[n] += tr[i];
                    break;
                }
            }
        }
    // Combining temperature update and bc in the same loop for efficiency,
    // since only the top boundary has Dirichlet bc, and all the other boundaries
    // have no heat flux bc.
        if ((*var.bcflag)[n] & BOUNDZ1)
            temperature[n] = param.bc.surface_temperature;
        else
            temperature[n] -= tdot[n] * var.dt / (*var.tmass)[n];
    }

    // Combining temperature update and bc in the same loop for efficiency,
    // since only the top boundary has Dirichlet bc, and all the other boundaries
    // have no heat flux bc.
//    #pragma omp parallel for default(none) shared(var, param, tdot, temperature)
//    for (int n=0; n<var.nnode; ++n) {
//        if ((*var.bcflag)[n] & BOUNDZ1)
//            temperature[n] = param.bc.surface_temperature;
//        else
//            temperature[n] -= tdot[n] * var.dt / (*var.tmass)[n];
//    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

// Function to check if a node is a boundary node
inline bool is_boundary_node_for_pp(const int n, const Variables &var) {
    int boundary_types[6] = {BOUNDX0, BOUNDX1, BOUNDY0, BOUNDY1, BOUNDZ0, BOUNDZ1};
    for (int i = 0; i < 6; ++i) {
        if (((*var.bcflag)[n] & boundary_types[i]) && (var.hbc_types[i] == 1)) {
            return true;
        }
    }
    return false;
}


void update_pore_pressure(const Param &param, const Variables &var,
                          double_vec &ppressure, double_vec &dppressure, double_vec &tdot, elem_cache &tmp_result, tensor_t& stress, double_vec& old_mean_stress)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    // Initialize diff_max_local for reduction
    double diff_max_local = 1.0e-38;

    #pragma omp parallel for default(none) shared(var, ppressure, tmp_result, stress, old_mean_stress, param, diff_max_local)
    #pragma acc parallel loop
    for (int e = 0; e < var.nelem; e++) {
        const int *conn = (*var.connectivity)[e];
        double *tr = tmp_result[e];
        double current_mean_stress = trace(stress[e])/NDIMS;
        double mean_stress_change = current_mean_stress - old_mean_stress[e];

        // Retrieve hydraulic properties for the element
        double perm_e = var.mat->perm(e);                // Intrinsic permeability 
        double mu_e = var.mat->mu_fluid(e);              // Fluid dynamic viscosity
        double alpha_b = var.mat->alpha_biot(e);         // Biot coefficient
        double rho_f = var.mat->rho_fluid(e);            // Fluid density
        double phi_e = var.mat->phi(e);        // Element porosity
        double comp_fluid = var.mat->beta_fluid(e);        // fluid comporessibility
        double bulkm = var.mat->bulkm(e);
        double shearm = var.mat->shearm(e);
        double matrix_comp = 1.0 / (bulkm +4.0*shearm/3.0);

        double bulk_comp = 1.0/(*var.mat).bulkm(e); // lambda + 2G/3
        if(NDIMS == 2) bulk_comp = 1.0/((*var.mat).bulkm(e) + (*var.mat).shearm(e)/3.0); // lambda + G 

        rho_f = 1000.0; 
        double gamma_w = rho_f * param.control.gravity; // specific weight
        
        // Hydraulic conductivity using permeability and viscosity
        double hydraulic_conductivity = perm_e * gamma_w / mu_e;
        double kv = hydraulic_conductivity * (*var.volume)[e];

        // Compute element diffusivity and update max using reduction
        double diff_e = hydraulic_conductivity / (phi_e * comp_fluid + alpha_b * matrix_comp) / gamma_w;
        diff_max_local = std::max(diff_max_local, diff_e);

        // volume term (poroelastic effect)
        double pe = alpha_b * mean_stress_change * bulk_comp * (*var.volume)[e] / NODES_PER_ELEM  / var.dt;
       
        const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
#endif
        const double *shpdz = (*var.shpdz)[e];

        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            double diffusion = 0.0;
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
#ifdef THREED
                diffusion += (shpdx[i] * shpdx[j] +
                             shpdy[i] * shpdy[j] +
                             shpdz[i] * shpdz[j]) * (ppressure[conn[j]]/gamma_w + (*var.coord)[conn[j]][NDIMS-1]); 
#else
                diffusion += (shpdx[i] * shpdx[j] +
                             shpdz[i] * shpdz[j]) * (ppressure[conn[j]]/gamma_w + (*var.coord)[conn[j]][NDIMS-1]);
#endif
                }

            // Add diffusion, compressibility, poroelastic effects, and source term
            tr[i] = kv * diffusion + pe;
            // tr[i] = kv * diffusion;
        }
    }

    // Update global hydro_diff_max after the loop
    var.mat->hydro_diff_max = diff_max_local;

    #pragma omp parallel for default(none) shared(param, var, tdot, ppressure, dppressure, tmp_result)
    #pragma acc parallel loop
    for (int n = 0; n < var.nnode; n++) {
        tdot[n] = 0.0;
        for (auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e) {
            const int *conn = (*var.connectivity)[*e];
            const double *tr = tmp_result[*e];
            for (int i = 0; i < NODES_PER_ELEM; i++) {
                if (n == conn[i]) {
                    tdot[n] += tr[i];
                    break;
                }
            }
        }

        // Update pore pressure for non-boundary nodes
        if (!is_boundary_node_for_pp(n, var)) {
            // Ensure mass is non-zero before division
            if ((*var.hmass)[n] > 0.0) {
                ppressure[n] -= tdot[n] * var.dt / (*var.hmass)[n]; // Update pore pressure
                dppressure[n] = tdot[n] * var.dt / (*var.hmass)[n]; // Update pressure change
            }
        }

        // if ((*var.bcflag)[n] & BOUNDZ1) {
        //     ppressure[n] = 0.0; 
        //     dppressure[n] = 0.0;
        // }
            
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void update_strain_rate(const Variables& var, tensor_t& strain_rate)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double *v[NODES_PER_ELEM];

    #pragma omp parallel for default(none) shared(var, strain_rate) private(v)
    #pragma acc parallel loop private(v)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
        const double *shpdz = (*var.shpdz)[e];
        double *s = strain_rate[e];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = (*var.vel)[conn[i]];

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
        // YY component
        n = 1;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][1] * shpdy[i];
#endif

        // ZZ component
#ifdef THREED
        n = 2;
#else
        n = 1;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][NDIMS-1] * shpdz[i];

#ifdef THREED
        // XY component
        n = 3;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdy[i] + v[i][1] * shpdx[i]);
#endif

        // XZ component
#ifdef THREED
        n = 4;
#else
        n = 2;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdz[i] + v[i][NDIMS-1] * shpdx[i]);

#ifdef THREED
        // YZ component
        n = 5;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][1] * shpdz[i] + v[i][2] * shpdy[i]);
#endif
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


static void apply_damping(const Param& param, const Variables& var, array_t& force)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // flatten 2d arrays to simplify indexing
    const double small_vel = 1e-13;

    switch (param.control.damping_option) {
    case 0:
        // no damping, stress field can become very noisy
        break;
    case 1:
        // damping when force and velocity are parallel
        // acclerating when force and velocity are anti-parallel
#ifdef GPP1X
        #pragma omp parallel for default(none) \
            shared(var, param, force, small_vel)
#else
        #pragma omp parallel for default(none) \
            shared(var, param, force)  firstprivate(small_vel)
#endif
        #pragma acc parallel loop
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0;j<NDIMS;j++)
                if (std::fabs((*var.vel)[i][j]) > small_vel) {
                    force[i][j] -= param.control.damping_factor * std::copysign(force[i][j], (*var.vel)[i][j]);
                }
        }
        break;
    case 2:
        // damping prop. to force
        #pragma omp parallel for default(none) \
            shared(var, param, force)
        #pragma acc parallel loop
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0;j<NDIMS;j++)
                force[i][j] -= param.control.damping_factor * force[i][j];
        }
        break;
    case 3:
        // damping when force and velocity are parallel
        // weakly acclerating when force and velocity are anti-parallel
        #pragma omp parallel for default(none) shared(var, param, force)
        #pragma acc parallel loop
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0;j<NDIMS;j++) {
                if ((force[i][j]<0) == ((*var.vel)[i][j]<0)) {
                    // strong damping
                    force[i][j] -= param.control.damping_factor * force[i][j], (*var.vel)[i][j];
                }
                else {
                    // weak acceleration
                    force[i][j] += (1 - param.control.damping_factor) * force[i][j];
                }
            }
        }
        break;
    case 4:
        // rayleigh damping
#ifdef GPP1X
        #pragma omp parallel for default(none) \
            shared(var, param, force, small_vel)
#else
        #pragma omp parallel for default(none) \
            shared(var, param, force)  firstprivate(small_vel)
#endif
        #pragma acc parallel loop
        for (int i=0; i<var.nnode; ++i) {
            double mass = (*var.mass)[i];
            double young = (*var.ymass)[i];
            double critical_coeff = 2.0 * std::sqrt(mass * young);
            for (int j=0;j<NDIMS;j++)
                if (std::fabs((*var.vel)[i][j]) > small_vel) {
                    double f_C = param.control.damping_factor * std::copysign(force[i][j], (*var.vel)[i][j]);
                    double f_V = critical_coeff * (*var.vel)[i][j];
                    double f_damping = (std::fabs(f_C) < std::fabs(f_V)) ? f_V : f_C;
                    force[i][j] -= f_damping;
                }
        }
        break;
    default:
        std::cerr << "Error: unknown damping_option: " << param.control.damping_option << '\n';
        std::exit(1);
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

/*
#pragma acc routine seq
static double rho(const conn_t &var_connectivity, \
    const double_vec &var_temperature, const int_vec2D &var_elemmarkers, \
    const double_vec &rho0, const double_vec &alpha, int nmat,  int e)
{
    const double celsius0 = 273;

    // average temperature of this element
    double T = 0;
    const int *conn = var_connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += var_temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    double TinCelsius = T - celsius0;
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        // TODO: compressibility
        result += rho0[m] * (1 - alpha[m] * TinCelsius) * var_elemmarkers[e][m];
        n += var_elemmarkers[e][m];
    }
    return result / n;
}
*/

void update_force(const Param& param, const Variables& var, array_t& force, array_t& force_residual, elem_cache& tmp_result)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma omp parallel for default(none) shared(var,param,tmp_result)
    #pragma acc parallel loop
    for (int e=0;e<var.nelem;e++) {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
#endif
        const double *shpdz = (*var.shpdz)[e];
        double *s = (*var.stress)[e];
        double vol = (*var.volume)[e];
        double *tr = tmp_result[e];

        double buoy = 0;
        if (param.control.gravity != 0) {
            // Calculate buoyancy based on the gravity and element properties
            // buoy = var.mat->rho(e) * param.control.gravity / NODES_PER_ELEM;
            buoy = (var.mat->rho(e) * (1 - var.mat->phi(e)) + 1000.0 * var.mat->phi(e)) * param.control.gravity / NODES_PER_ELEM;
        }


        for (int i=0; i<NODES_PER_ELEM; ++i) {
#ifdef THREED
            tr[i] = (s[0]*shpdx[i] + s[3]*shpdy[i] + s[4]*shpdz[i]) * vol;
            tr[i+NODES_PER_ELEM] = (s[3]*shpdx[i] + s[1]*shpdy[i] + s[5]*shpdz[i]) * vol;
            tr[i+NODES_PER_ELEM*2] = (s[4]*shpdx[i] + s[5]*shpdy[i] + s[2]*shpdz[i] + buoy) * vol;
#else
            tr[i] = (s[0]*shpdx[i] + s[2]*shpdz[i]) * vol;
            tr[i+NODES_PER_ELEM] = (s[2]*shpdx[i] + s[1]*shpdz[i] + buoy) * vol;
#endif
        }
    }

    #pragma omp parallel for default(none) shared(var,force,force_residual,tmp_result)
    #pragma acc parallel loop
    for (int n=0;n<var.nnode;n++) {
        std::fill_n(force[n],NDIMS,0); 
        double *f = force[n];
        std::fill_n(force_residual[n],NDIMS,0); 
        double *f_residual = force_residual[n];
        for( auto e = (*var.support)[n].begin(); e < (*var.support)[n].end(); ++e) {
            const int *conn = (*var.connectivity)[*e];
            const double *tr = tmp_result[*e];
            for (int i=0;i<NODES_PER_ELEM;i++) {
                if (n == conn[i]) {
                    for (int j=0;j<NDIMS;j++)
                    {
                        f[j] -= tr[i+NODES_PER_ELEM*j];
                        f_residual[j] = tr[i+NODES_PER_ELEM*j];
                    }
                    break;
                }
            }
        }
    }

    apply_stress_bcs(param, var, force);

    // if(var.time <= 1.0)
    // {
    //     apply_stress_bcs_neumann(param, var, force);
    // }

    // if (param.control.is_quasi_static) {
    //     apply_damping(param, var, force);
    // }

    if (!param.ic.has_body_force_adjustment) apply_stress_bcs_neumann(param, var, force);
    apply_damping(param, var, force);
    
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

double calculate_residual_force(const Variables& var, array_t& force_residual)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double l2 = 0.0;
    double num = var.nnode * NDIMS;

    #pragma omp parallel for default(none) shared(var, force_residual, num) reduction(+:l2)
    #pragma acc parallel loop
    for (int i = 0; i < var.nnode; ++i)
        for (int j = 0; j < NDIMS; ++j)
            l2 += std::pow(force_residual[i][j], 2) / num;

#ifdef USE_NPROF
    nvtxRangePop();
#endif

    return std::sqrt(l2);
}


void update_velocity(const Variables& var, array_t& vel)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma omp parallel for default(none) shared(var, vel)
    #pragma acc parallel loop
    for (int i=0; i<var.nnode; ++i)
        for (int j=0;j<NDIMS;j++)
            vel[i][j] += var.dt * (*var.force)[i][j] / (*var.mass)[i];

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void update_velocity_PT(const Variables& var, array_t& vel)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    #pragma omp parallel for default(none) shared(var, vel)
    #pragma acc parallel loop
    for (int i=0; i<var.nnode; ++i)
        for (int j=0;j<NDIMS;j++)
            vel[i][j] += var.dt_PT * (*var.force)[i][j] / (*var.mass)[i];

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void update_coordinate(const Variables& var, array_t& coord)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // double* x = var.coord->data();
    // const double* v = var.vel->data();

    // for gpu parallelization dt and bound need to be sent to 

    #pragma omp parallel for default(none) shared(var)
    #pragma acc parallel loop collapse(2) 
    for (int i=0; i<var.nnode; ++i) {
        for (int j=0 ; j<NDIMS; ++j){
            (*var.coord)[i][j] += (*var.vel)[i][j] * var.dt;
        }
        
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


namespace {

#ifdef THREED
    #pragma acc routine seq 
    void jaumann_rate_3d(double *s, double dt, double w3, double w4, double w5)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[3] * w3 - 2.0 * s[4] * w4;
        s_inc[1] = 2.0 * s[3] * w3 - 2.0 * s[5] * w5;
        s_inc[2] = 2.0 * s[4] * w4 + 2.0 * s[5] * w5;
        s_inc[3] = s[0] * w3 - s[1] * w3 - s[4] * w5 - s[5] * w4;
        s_inc[4] = s[0] * w4 - s[2] * w4 + s[3] * w5 - s[5] * w3;
        s_inc[5] = s[1] * w5 - s[2] * w5 + s[3] * w4 + s[4] * w3;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }

#else
    #pragma acc routine seq
    void jaumann_rate_2d(double *s, double dt, double w2)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[2] * w2;
        s_inc[1] = 2.0 * s[2] * w2;
        s_inc[2] = s[0] * w2 - s[1] * w2;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }

#endif

}


void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // The spin rate tensor, W, and the Cauchy stress tensor, S, are
    //     [  0  w3  w4]     [s0 s3 s4]
    // W = [-w3   0  w5],  S=[s3 s1 s5].
    //     [-w4 -w5   0]     [s4 s5 s2]
    //
    // Stress (and strain) increment based on the Jaumann rate is
    // dt*(S*W-W*S). S*W-W*S is also symmetric.
    // So, following the indexing of stress tensor,
    // we get
    // sj[0] = dt * ( -2 * s3 * w3 - 2 * s4 * w4)
    // sj[1] = dt * (  2 * s3 * w3 - 2 * s5 * w5)
    // sj[2] = dt * (  2 * s4 * w4 + 2 * s5 * w5)
    // sj[3] = dt * ( s0 * w3 - s1 * w3 - s4 * w5 - s5 * w4)
    // sj[4] = dt * ( s0 * w4 - s2 * w4 + s3 * w5 - s5 * w3)
    // sj[5] = dt * ( s1 * w5 - s2 * w5 + s3 * w4 + s4 * w3)

    #pragma omp parallel for default(none) shared(var, stress, strain)
    #pragma acc parallel loop
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];

#ifdef THREED

        double w3, w4, w5;
        {
            const double *shpdx = (*var.shpdx)[e];
            const double *shpdy = (*var.shpdy)[e];
            const double *shpdz = (*var.shpdz)[e];

            double *v[NODES_PER_ELEM];
            for (int i=0; i<NODES_PER_ELEM; ++i)
                v[i] = (*var.vel)[conn[i]];

            w3 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w3 += 0.5 * (v[i][0] * shpdy[i] - v[i][1] * shpdx[i]);

            w4 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w4 += 0.5 * (v[i][0] * shpdz[i] - v[i][NDIMS-1] * shpdx[i]);

            w5 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w5 += 0.5 * (v[i][1] * shpdz[i] - v[i][NDIMS-1] * shpdy[i]);
        }

        jaumann_rate_3d(stress[e], var.dt, w3, w4, w5);
        jaumann_rate_3d(strain[e], var.dt, w3, w4, w5);

#else

        double w2;
        {
            const double *shpdx = (*var.shpdx)[e];
            const double *shpdz = (*var.shpdz)[e];

            double *v[NODES_PER_ELEM];
            for (int i=0; i<NODES_PER_ELEM; ++i)
                v[i] = (*var.vel)[conn[i]];

            w2 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w2 += 0.5 * (v[i][NDIMS-1] * shpdx[i] - v[i][0] * shpdz[i]);
        }

        jaumann_rate_2d(stress[e], var.dt, w2);
        jaumann_rate_2d(strain[e], var.dt, w2);

#endif
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

