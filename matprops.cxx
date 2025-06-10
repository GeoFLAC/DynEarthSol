
#include <algorithm>
#include <cmath>

#include "constants.hpp"
#include "utils.hpp"
#include "matprops.hpp"

#ifdef ACC
// ACC version

namespace {

    #pragma acc routine seq
    double get_prem_pressure(double depth)
    {
        // reference pressure profile from isotropic PREM model
        const int nlayers = 46;
        const double
            ref_depth[] = { 0e3,    3e3,    15e3,   24.4e3, 40e3,
                            60e3,   80e3,   115e3,  150e3,  185e3,
                            220e3,  265e3,  310e3,  355e3,  400e3,
                            450e3,  500e3,  550e3,  600e3,  635e3,
                            670e3,  721e3,  771e3,  871e3,  971e3,
                            1071e3, 1171e3, 1271e3, 1371e3, 1471e3,
                            1571e3, 1671e3, 1771e3, 1871e3, 1971e3,
                            2071e3, 2171e3, 2271e3, 2371e3, 2471e3,
                            2571e3, 2671e3, 2741e3, 2771e3, 2871e3,
                            2891e3 };

        // pressure in PREM table is given in kilobar, converted to 10^8 Pa
        const double
            ref_pressure[] = { 0e8,      0.3e8,    3.3e8,    6.0e8,    11.2e8,
                               17.8e8,   24.5e8,   36.1e8,   47.8e8,   59.4e8,
                               71.1e8,   86.4e8,   102.0e8,  117.7e8,  133.5e8,
                               152.2e8,  171.3e8,  190.7e8,  210.4e8,  224.3e8,
                               238.3e8,  260.7e8,  282.9e8,  327.6e8,  372.8e8,
                               418.6e8,  464.8e8,  511.6e8,  558.9e8,  606.8e8,
                               655.2e8,  704.1e8,  753.5e8,  803.6e8,  854.3e8,
                               905.6e8,  957.6e8,  1010.3e8, 1063.8e8, 1118.2e8,
                               1173.4e8, 1229.7e8, 1269.7e8, 1287.0e8, 1345.6e8,
                               1357.5e8 };

        // PREM model doesn't work if depth is above sea level, always returns 0 pressure
        if (depth <= 0) return 0;

        int n;
        for (n=1; n<nlayers; n++) {
            if (depth <= ref_depth[n]) break;
        }

        // linear interpolation
        double pressure = ref_pressure[n-1] + (ref_pressure[n] - ref_pressure[n-1]) *
            (depth - ref_depth[n-1]) / (ref_depth[n] - ref_depth[n-1]);

        return pressure;
    }


    #pragma acc routine seq
    double get_prem_pressure_modified(double depth)
    {
        // reference pressure profile from isotropic PREM model, modified for
        // average continental crust (density 2800 kg/m^3, thickness 24.4 km)
        const int nlayers = 46;
        const double
            ref_depth[] = { 0e3,    3e3,    15e3,   24.4e3, 40e3,
                            60e3,   80e3,   115e3,  150e3,  185e3,
                            220e3,  265e3,  310e3,  355e3,  400e3,
                            450e3,  500e3,  550e3,  600e3,  635e3,
                            670e3,  721e3,  771e3,  871e3,  971e3,
                            1071e3, 1171e3, 1271e3, 1371e3, 1471e3,
                            1571e3, 1671e3, 1771e3, 1871e3, 1971e3,
                            2071e3, 2171e3, 2271e3, 2371e3, 2471e3,
                            2571e3, 2671e3, 2741e3, 2771e3, 2871e3,
                            2891e3 };

        // pressure in PREM table is given in kilobar, converted to 10^8 Pa
        const double
            ref_pressure[] = { 0e8,      0.82e8,    4.1e8,    6.7e8,    11.2e8,
                               17.8e8,   24.5e8,   36.1e8,   47.8e8,   59.4e8,
                               71.1e8,   86.4e8,   102.0e8,  117.7e8,  133.5e8,
                               152.2e8,  171.3e8,  190.7e8,  210.4e8,  224.3e8,
                               238.3e8,  260.7e8,  282.9e8,  327.6e8,  372.8e8,
                               418.6e8,  464.8e8,  511.6e8,  558.9e8,  606.8e8,
                               655.2e8,  704.1e8,  753.5e8,  803.6e8,  854.3e8,
                               905.6e8,  957.6e8,  1010.3e8, 1063.8e8, 1118.2e8,
                               1173.4e8, 1229.7e8, 1269.7e8, 1287.0e8, 1345.6e8,
                               1357.5e8 };

        // PREM model doesn't work if depth is above sea level, always returns 0 pressure
        if (depth <= 0) return 0;

        int n;
        for (n=1; n<nlayers; n++) {
            if (depth <= ref_depth[n]) break;
        }

        // linear interpolation
        double pressure = ref_pressure[n-1] + (ref_pressure[n] - ref_pressure[n-1]) *
            (depth - ref_depth[n-1]) / (ref_depth[n] - ref_depth[n-1]);

        return pressure;
    }

    using VectorBase = double_vec;

    double arithmetic_mean(const VectorBase &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        for (std::size_t i=0; i<s.size(); i++) {
            result += n[i] * s[i];
            m += n[i];
        }
        return result / m;
    }


    double harmonic_mean(const VectorBase &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        for (std::size_t i=0; i<s.size(); i++) {
            result += n[i] / s[i];
            m += n[i];
        }
        return m / result;
    }
}


double ref_pressure(const Param& param, double z)
{
    // Get pressure at this depth
    double depth = -z;
    double p;
    int mattype_ref = param.mat.mattype_mantle;
    if (param.control.ref_pressure_option == 0)
        if (param.control.has_hydraulic_diffusion) {
        // Modified density considering porosity for hydraulic diffusion
            p = (param.mat.rho0[mattype_ref] * (1 - param.mat.porosity[mattype_ref]) + 1000.0 * param.mat.porosity[mattype_ref]) * param.control.gravity * depth;
        } else {
            // Standard reference pressure without hydraulic diffusion
            p = param.mat.rho0[mattype_ref] * param.control.gravity * depth;
        }

    else if (param.control.ref_pressure_option == 1)
        p = get_prem_pressure(depth);
    else if (param.control.ref_pressure_option == 2)
        p = get_prem_pressure_modified(depth);
    return p;
}


MatProps::MatProps(const Param& p, const Variables& var) :
  rheol_type(p.mat.rheol_type),
  nmat(p.mat.nmat),
  is_plane_strain(p.mat.is_plane_strain),
  visc_min(p.mat.visc_min),
  visc_max(p.mat.visc_max),
  tension_max(p.mat.tension_max),
  therm_diff_max(p.mat.therm_diff_max),
  hydro_diff_max(1e-1),
  rho_magma(p.mat.rho_magma),
  has_magmatism(p.control.has_magmatism),
  coord(*var.coord),
  connectivity(*var.connectivity),
  temperature(*var.temperature),
  stress(*var.stress),
//   stress_old(*var.stress_old),
  strain_rate(*var.strain_rate),
  elemmarkers(*var.elemmarkers),
  ppressure(*var.ppressure),
  dppressure(*var.dppressure),
  eff_fmelt(*var.eff_fmelt),
  magma_fraction(*var.magma_fraction),
  hydro_diff(*var.hydro_diff)
{
    rho0 = p.mat.rho0;
    alpha = p.mat.alpha;
    bulk_modulus = p.mat.bulk_modulus;
    shear_modulus = p.mat.shear_modulus;
    visc_exponent = p.mat.visc_exponent;
    visc_coefficient = p.mat.visc_coefficient;
    visc_activation_energy = p.mat.visc_activation_energy;
    visc_activation_volume = p.mat.visc_activation_volume;
    heat_capacity = p.mat.heat_capacity;
    therm_cond = p.mat.therm_cond;
    pls0 = p.mat.pls0;
    pls1 = p.mat.pls1;
    cohesion0 = p.mat.cohesion0;
    cohesion1 = p.mat.cohesion1;
    friction_angle0 = p.mat.friction_angle0;
    friction_angle1 = p.mat.friction_angle1;
    dilation_angle0 = p.mat.dilation_angle0;
    dilation_angle1 = p.mat.dilation_angle1;

    // Hydraulic parameters
    porosity = p.mat.porosity;
    hydraulic_perm = p.mat.hydraulic_perm;
    fluid_rho0 = p.mat.fluid_rho0;
    fluid_alpha = p.mat.fluid_alpha;
    fluid_bulk_modulus = p.mat.fluid_bulk_modulus;
    fluid_visc = p.mat.fluid_visc;
    biot_coeff = p.mat.biot_coeff;
    bulk_modulus_s = p.mat.bulk_modulus_s;
    // Rate-and-state friction parameters
    direct_a = p.mat.direct_a;
    evolution_b = p.mat.evolution_b;
    characteristic_velocity = p.mat.characteristic_velocity;
    // static_friction_coefficient = p.mat.static_friction_coefficient;

}


MatProps::~MatProps()
{
    #pragma acc exit data delete(rho0,alpha,bulk_modulus,shear_modulus,visc_exponent)
    #pragma acc exit data delete(visc_coefficient,visc_activation_energy,visc_activation_volume,heat_capacity,therm_cond,pls0)
    #pragma acc exit data delete(pls1,cohesion0,friction_angle0,friction_angle1,dilation_angle0,dilation_angle1)
    #pragma acc exit data delete(porosity,hydraulic_perm,fluid_rho0,fluid_alpha,fluid_bulk_modulus,fluid_visc,biot_coeff,bulk_modulus_s)
    #pragma acc exit data delete(direct_a,evolution_b,characteristic_velocity)
}


double MatProps::bulkm(int e) const
{
    return harmonic_mean(bulk_modulus, elemmarkers[e]);
}


double MatProps::shearm(int e) const
{
    return harmonic_mean(shear_modulus, elemmarkers[e]);
}


double MatProps::visc(int e) const
{
    const double gas_constant = 8.3144;
    const double min_strain_rate = 1e-30;

    // average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    // stress
    const double *s = stress[e];
    double s0 = trace(s) / NDIMS;

    // strain-rate
    double edot = second_invariant(strain_rate[e]);
    // min strain rate to prevent viscosity -> inf
    edot = std::max(edot, min_strain_rate);

    // viscosity law from Chen and Morgan, JGR, 1990
    double result = 0;
    int n = 0;
    
        if (has_magmatism) {
        double fmagma = magma_fraction[e];
        double T_lava = 1473.; // 1200.+273.
        T = (1. - fmagma) * T + fmagma * T_lava;
    }

    for (int m=0; m<nmat; m++) {
        double pow = 1 / visc_exponent[m] - 1;
        double pow1 = -1 / visc_exponent[m];
        double visc0 = 0.25 * std::pow(edot, pow) * std::pow(0.75 * visc_coefficient[m], pow1)
            * std::exp((visc_activation_energy[m] + visc_activation_volume[m] * s0)
            / (visc_exponent[m] * gas_constant * T)) * 1e6;
        result += elemmarkers[e][m] / visc0;
        n += elemmarkers[e][m];
    }

    double visc = n / result;

    // applying min & max limits
    visc = std::min(std::max(visc, visc_min), visc_max);

    return visc;
}


void MatProps::plastic_weakening(int e, double pls,
                                 double &cohesion, double &friction_angle,
                                 double &dilation_angle, double &hardening) const
{
    double c, f, d, h;
    c = f = d = h = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        if (pls < pls0[m]) {
            // no weakening yet
            c += cohesion0[m] * k;
            f += friction_angle0[m] * k;
            d += dilation_angle0[m] * k;
            h += 0;
        }
        else if (pls < pls1[m]) {
            // linear weakening
            double p = (pls - pls0[m]) / (pls1[m] - pls0[m]);
            c += (cohesion0[m] + p * (cohesion1[m] - cohesion0[m])) * k;
            f += (friction_angle0[m] + p * (friction_angle1[m] - friction_angle0[m])) * k;
            d += (dilation_angle0[m] + p * (dilation_angle1[m] - dilation_angle0[m])) * k;
            h += (cohesion1[m] - cohesion0[m]) / (pls1[m] - pls0[m]) * k;
        }
        else {
            // saturated weakening
            c += cohesion1[m] * k;
            f += friction_angle1[m] * k;
            d += dilation_angle1[m] * k;
            h += 0;
        }
    }
    cohesion = c / n;
    friction_angle = f / n;
    dilation_angle = d / n;
    hardening = h / n;
}

void MatProps::plastic_weakening_rsf(int e, double pls,
                                 double &cohesion, double &friction_angle,
                                 double &dilation_angle, double &hardening, double &slip_rate) const
{
    double c, f, d, h;
    c = f = d = h = 0;

    double d_a, e_b, c_v, mu_0, mu_d;
    d_a = e_b = c_v = mu_0 = mu_d = 0;

    double static_friction_angle = 0;

    int n = 0;
    for (int m=0; m<nmat; m++) {
        int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        if (pls < pls0[m]) {
            // no weakening yet
            c += cohesion0[m] * k;
            f += friction_angle0[m] * k;
            d += dilation_angle0[m] * k;
            h += 0;
        }
        else if (pls < pls1[m]) {
            // linear weakening
            double p = (pls - pls0[m]) / (pls1[m] - pls0[m]);
            c += (cohesion0[m] + p * (cohesion1[m] - cohesion0[m])) * k;
            f += (friction_angle0[m] + p * (friction_angle1[m] - friction_angle0[m])) * k;
            d += (dilation_angle0[m] + p * (dilation_angle1[m] - dilation_angle0[m])) * k;
            h += (cohesion1[m] - cohesion0[m]) / (pls1[m] - pls0[m]) * k;
        }
        else {
            // saturated weakening
            c += cohesion1[m] * k;
            f += friction_angle1[m] * k;
            d += dilation_angle1[m] * k;
            h += 0;
        }

        // Rate-and-state friction parameters
        d_a += direct_a[m] * k;
        e_b += evolution_b[m] * k;
        c_v += characteristic_velocity[m] * k;
    }

    d_a /= n; // direct effect parameter
    e_b /= n; // evolution effect parameter   
    c_v /= n; // characteristic velocity
    static_friction_angle = f / n; // friction angle
    mu_0 = std::tan(DEG2RAD * static_friction_angle); // static friction coefficient
    mu_d = mu_0 + (d_a - e_b) * log(slip_rate / c_v); // dynamic friction angle
    friction_angle = std::atan(mu_d) / DEG2RAD;
    cohesion = c / n;
    dilation_angle = d / n;
    hardening = h / n;
}

void MatProps::plastic_props(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max) const
{
    double cohesion, phi, psi;

    plastic_weakening(e, pls, cohesion, phi, psi, hardn);

    // derived variables
    double sphi = std::sin(phi * DEG2RAD);
    double spsi = std::sin(psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? tension_max : std::min(tension_max, cohesion/std::tan(phi*DEG2RAD));
}

void MatProps::plastic_props_rsf(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max, double& slip_rate) const
{
    double cohesion, phi, psi;

    plastic_weakening_rsf(e, pls, cohesion, phi, psi, hardn, slip_rate);

    // derived variables
    double sphi = std::sin(phi * DEG2RAD);
    double spsi = std::sin(psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? tension_max : std::min(tension_max, cohesion/std::tan(phi*DEG2RAD));
}


double MatProps::rho(int e) const
{
    const double celsius0 = 273;

    // average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    double TinCelsius = T - celsius0;
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        // TODO: compressibility
        result += rho0[m] * (1 - alpha[m] * TinCelsius) * elemmarkers[e][m];
        n += elemmarkers[e][m];
    }
    double rho = result / n;
    if (has_magmatism) {
        double tt_melt = std::min(1., magma_fraction[e] + eff_fmelt[e]);
        return (1 - tt_melt) * rho + tt_melt * rho_magma;
    } else {
        return rho;
    }
}


void  MatProps::update_hydro_diff(const Param &param, const Variables &var,
                        double_vec &var_hydro_diff, double_vec &tdot, elem_cache &tmp_result) const
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double separate_rate = 10;
    double half_life = 1e6 * YEAR2SEC;
    double decay_rate = std::pow(0.5, var.dt/half_life);
    double_vec &etmp = *var.etmp;

    double diff_multi,depth_max,pls_min,t_max,xh1t,yh1t;
    pls_min = 1.;
    diff_multi = 6.;
    depth_max = -20.e3;
    t_max = 600. + 273.;
    xh1t = 30.e3;
    yh1t = 0.e3;

    #pragma omp parallel default(none) shared(param,var,tdot,var_hydro_diff, \
        tmp_result,separate_rate,decay_rate,coord,connectivity,temperature, \
        pls_min,diff_multi,depth_max,t_max,xh1t,yh1t,etmp)
    {
        #pragma omp for
        // #pragma acc parallel loop
        for (int e=0;e<var.nelem;e++) {
            // diffusion matrix
            const int *conn = (*var.connectivity)[e];
            double *tr = tmp_result[e];
            double kv = separate_rate * (*var.volume)[e];
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
                                shpdz[i] * shpdz[j]) * var_hydro_diff[conn[j]];
    #else
                    diffusion += (shpdx[i] * shpdx[j] + \
                                shpdz[i] * shpdz[j]) * var_hydro_diff[conn[j]];
    #endif
                }
                tr[i] = diffusion * kv;
            }
        }

        #pragma omp for
        // #pragma acc parallel loop
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
            var_hydro_diff[n] -= tdot[n] * var.dt / (*var.tmass)[n];
            var_hydro_diff[n] *= decay_rate;
        }

        #pragma omp for
        for (int e=0;e<var.nelem;e++) {
            if ((*var.plstrain)[e] <= pls_min) continue;

            const int *conn = connectivity[e];
            double T = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                T += temperature[conn[i]];
            T /= NODES_PER_ELEM;
            if (T > t_max) continue;

            double pe[NDIMS];
            for (int i=0; i<NDIMS; ++i) {
                for (int j=0; j<NODES_PER_ELEM; ++j)
                    pe[i] += coord[conn[j]][i];
                pe[i] /= NODES_PER_ELEM;
            }
            if (pe[0] > xh1t && pe[1] > yh1t && 
                pe[0] < (param.mesh.xlength - xh1t) && pe[1] < (param.mesh.ylength - yh1t) &&
                pe[2] >= depth_max) {
                etmp[i] = diff_multi;
            } else {
                etmp[i] = 0;
            }
        }

        #pragma omp for
        for (int n=0;n<var.nnode;n++) {
            int_vec &sup = (*var.supper)[n]
            for (int i=0; i<sup.size();i++) {
                if (etmp[sup[i]]>0.) {
                    var_hydro_diff[n] = etmp[sup[i]]
                    updated = true;
                    break;
                }
            }
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

double MatProps::cp(int e) const
{
    return arithmetic_mean(heat_capacity, elemmarkers[e]);
}


double MatProps::k(int e) const
{
    double ek = arithmetic_mean(therm_cond, elemmarkers[e]);

    if (has_magmatism) {
        const int* conn = connectivity[e];
        double ehydro = 0;
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            ehydro += hydro_diff[conn[i]];
        ehydro /= NODES_PER_ELEM;
        return ek * (1. + ehydro);
    } else {
        return ek;
    }
}

// hydraulic parameters
double MatProps::phi(int e) const
{
    return arithmetic_mean(porosity, elemmarkers[e]);
}

double MatProps::perm(int e) const
{
    return harmonic_mean(hydraulic_perm, elemmarkers[e]);
}

double MatProps::alpha_fluid(int e) const
{
    return arithmetic_mean(fluid_alpha, elemmarkers[e]);
}

double MatProps::beta_fluid(int e) const
{
    // Return the inverse of harmonic mean (compressibility)
    return 1.0 / harmonic_mean(fluid_bulk_modulus, elemmarkers[e]);
}

double MatProps::rho_fluid(int e) const
{
    const double p0 = 1.0e5; // Reference pressure in Pascals (1 atm)
    const double T0 = 273; // Reference temperature in Kelvin (0°C)

    // Average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    // Average pore pressure of this element
    double p = 0;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        p += ppressure[conn[i]];
    }
    p /= NODES_PER_ELEM;

    // Use the functions to get thermal expansivity and compressibility
    double alpha_f = alpha_fluid(e);  // Thermal expansivity for the fluid
    double beta_f = beta_fluid(e);    // Compressibility for the fluid

    // Calculate fluid density based on thermal expansivity and compressibility
    double result = 0;
    int n = 0;
    for (int m = 0; m < nmat; m++) {
        double rho_f = (fluid_rho0)[m];  // Reference fluid density
        result += rho_f * (1 - alpha_f * (T - T0) + beta_f * (p - p0)) * elemmarkers[e][m];
        n += elemmarkers[e][m];
    }

    // Return the averaged fluid density
    return result / n;
}

double MatProps::mu_fluid(int e) const
{
    return arithmetic_mean(fluid_visc, elemmarkers[e]);
}

double MatProps::alpha_biot(int e) const
{
    return arithmetic_mean(biot_coeff, elemmarkers[e]);
}

double MatProps::beta_mineral(int e) const
{
    // Return the inverse of harmonic mean (compressibility)
    return 1.0 / harmonic_mean(bulk_modulus_s, elemmarkers[e]);
}

// Rate-and-state friction parameters
double MatProps::d_a(int e) const
{
    return arithmetic_mean(direct_a, elemmarkers[e]);
}

double MatProps::e_b(int e) const
{
    return arithmetic_mean(evolution_b, elemmarkers[e]);
}

double MatProps::c_v(int e) const
{
    return arithmetic_mean(characteristic_velocity, elemmarkers[e]);
}


#else

namespace {

    double get_prem_pressure(double depth)
    {
        // reference pressure profile from isotropic PREM model
        const int nlayers = 46;
        const static double
            ref_depth[] = { 0e3,    3e3,    15e3,   24.4e3, 40e3,
                            60e3,   80e3,   115e3,  150e3,  185e3,
                            220e3,  265e3,  310e3,  355e3,  400e3,
                            450e3,  500e3,  550e3,  600e3,  635e3,
                            670e3,  721e3,  771e3,  871e3,  971e3,
                            1071e3, 1171e3, 1271e3, 1371e3, 1471e3,
                            1571e3, 1671e3, 1771e3, 1871e3, 1971e3,
                            2071e3, 2171e3, 2271e3, 2371e3, 2471e3,
                            2571e3, 2671e3, 2741e3, 2771e3, 2871e3,
                            2891e3 };

        // pressure in PREM table is given in kilobar, converted to 10^8 Pa
        const static double
            ref_pressure[] = { 0e8,      0.3e8,    3.3e8,    6.0e8,    11.2e8,
                               17.8e8,   24.5e8,   36.1e8,   47.8e8,   59.4e8,
                               71.1e8,   86.4e8,   102.0e8,  117.7e8,  133.5e8,
                               152.2e8,  171.3e8,  190.7e8,  210.4e8,  224.3e8,
                               238.3e8,  260.7e8,  282.9e8,  327.6e8,  372.8e8,
                               418.6e8,  464.8e8,  511.6e8,  558.9e8,  606.8e8,
                               655.2e8,  704.1e8,  753.5e8,  803.6e8,  854.3e8,
                               905.6e8,  957.6e8,  1010.3e8, 1063.8e8, 1118.2e8,
                               1173.4e8, 1229.7e8, 1269.7e8, 1287.0e8, 1345.6e8,
                               1357.5e8 };

        // PREM model doesn't work if depth is above sea level, always returns 0 pressure
        if (depth <= 0) return 0;

        int n;
        for (n=1; n<nlayers; n++) {
            if (depth <= ref_depth[n]) break;
        }

        // linear interpolation
        double pressure = ref_pressure[n-1] + (ref_pressure[n] - ref_pressure[n-1]) *
            (depth - ref_depth[n-1]) / (ref_depth[n] - ref_depth[n-1]);

        return pressure;
    }


    double get_prem_pressure_modified(double depth)
    {
        // reference pressure profile from isotropic PREM model, modified for
        // average continental crust (density 2800 kg/m^3, thickness 24.4 km)
        const int nlayers = 46;
        const static double
            ref_depth[] = { 0e3,    3e3,    15e3,   24.4e3, 40e3,
                            60e3,   80e3,   115e3,  150e3,  185e3,
                            220e3,  265e3,  310e3,  355e3,  400e3,
                            450e3,  500e3,  550e3,  600e3,  635e3,
                            670e3,  721e3,  771e3,  871e3,  971e3,
                            1071e3, 1171e3, 1271e3, 1371e3, 1471e3,
                            1571e3, 1671e3, 1771e3, 1871e3, 1971e3,
                            2071e3, 2171e3, 2271e3, 2371e3, 2471e3,
                            2571e3, 2671e3, 2741e3, 2771e3, 2871e3,
                            2891e3 };

        // pressure in PREM table is given in kilobar, converted to 10^8 Pa
        const static double
            ref_pressure[] = { 0e8,      0.82e8,    4.1e8,    6.7e8,    11.2e8,
                               17.8e8,   24.5e8,   36.1e8,   47.8e8,   59.4e8,
                               71.1e8,   86.4e8,   102.0e8,  117.7e8,  133.5e8,
                               152.2e8,  171.3e8,  190.7e8,  210.4e8,  224.3e8,
                               238.3e8,  260.7e8,  282.9e8,  327.6e8,  372.8e8,
                               418.6e8,  464.8e8,  511.6e8,  558.9e8,  606.8e8,
                               655.2e8,  704.1e8,  753.5e8,  803.6e8,  854.3e8,
                               905.6e8,  957.6e8,  1010.3e8, 1063.8e8, 1118.2e8,
                               1173.4e8, 1229.7e8, 1269.7e8, 1287.0e8, 1345.6e8,
                               1357.5e8 };

        // PREM model doesn't work if depth is above sea level, always returns 0 pressure
        if (depth <= 0) return 0;

        int n;
        for (n=1; n<nlayers; n++) {
            if (depth <= ref_depth[n]) break;
        }

        // linear interpolation
        double pressure = ref_pressure[n-1] + (ref_pressure[n] - ref_pressure[n-1]) *
            (depth - ref_depth[n-1]) / (ref_depth[n] - ref_depth[n-1]);

        return pressure;
    }


    double arithmetic_mean(const VectorBase &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        for (std::size_t i=0; i<s.size(); i++) {
            result += n[i] * s[i];
            m += n[i];
        }
        return result / m;
    }


    double harmonic_mean(const VectorBase &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        for (std::size_t i=0; i<s.size(); i++) {
            result += n[i] / s[i];
            m += n[i];
        }
        return m / result;
    }


    class Vector : public VectorBase {
    private:
        const double_vec &a;
        const int len;
    public:
        Vector(const double_vec &a_, int len_) :
            a(a_), len(len_)
        {}

        double operator[](std::size_t i) const
        {
            return a[i];
        }

        std::size_t size() const
        {
            return a.size();
        }
    };


    class Vector1 : public VectorBase {
    private:
        const double d;
        const int len;
    public:
        Vector1(const double_vec &a_, int len_) :
            d(a_[0]), len(len_)
        {}

        double operator[](std::size_t i) const
        {
            return d; // always return the same element
        }

        std::size_t size() const
        {
            return 1;
        }
    };


}


VectorBase* VectorBase::create(const double_vec &a, int len)
{
    if (static_cast<int>(a.size()) == len)
        return new Vector(a, len);
    if (a.size() == 1 && len != 1)
        return new Vector1(a, len);

    std::cerr << "Error: incorrect parameters received in VectorBase::create() at "
              << __FILE__ << ':' << __LINE__ << '\n';
    std::exit(12);
    return NULL;
}


double ref_pressure(const Param& param, double z)
{
    // Get pressure at this depth
    int mattype_ref = param.mat.mattype_mantle;
    double depth = -z;
    double p;
    if (param.control.ref_pressure_option == 0)
        if (param.control.has_hydraulic_diffusion) {
        // Modified density considering porosity for hydraulic diffusion
            p = (param.mat.rho0[mattype_ref] * (1 - param.mat.porosity[mattype_ref]) + 1000.0 * param.mat.porosity[mattype_ref]) * param.control.gravity * depth;
        } else {
            // Standard reference pressure without hydraulic diffusion
            p = param.mat.rho0[mattype_ref] * param.control.gravity * depth;
        }

    else if (param.control.ref_pressure_option == 1)
        p = get_prem_pressure(depth);
    else if (param.control.ref_pressure_option == 2)
        p = get_prem_pressure_modified(depth);
    return p;
}


MatProps::MatProps(const Param& p, const Variables& var) :
  rheol_type(p.mat.rheol_type),
  nmat(p.mat.nmat),
  is_plane_strain(p.mat.is_plane_strain),
  visc_min(p.mat.visc_min),
  visc_max(p.mat.visc_max),
  tension_max(p.mat.tension_max),
  therm_diff_max(p.mat.therm_diff_max),
  hydro_diff_max(1e-1),
  coord(*var.coord),
  connectivity(*var.connectivity),
  temperature(*var.temperature),
  stress(*var.stress),
//   stress_old(*var.stress_old),
  strain_rate(*var.strain_rate),
  elemmarkers(*var.elemmarkers),
  ppressure(*var.ppressure),
  dppressure(*var.dppressure)
{
    rho0 = VectorBase::create(p.mat.rho0, nmat);
    alpha = VectorBase::create(p.mat.alpha, nmat);
    bulk_modulus = VectorBase::create(p.mat.bulk_modulus, nmat);
    shear_modulus = VectorBase::create(p.mat.shear_modulus, nmat);
    visc_exponent = VectorBase::create(p.mat.visc_exponent, nmat);
    visc_coefficient = VectorBase::create(p.mat.visc_coefficient, nmat);
    visc_activation_energy = VectorBase::create(p.mat.visc_activation_energy, nmat);
    visc_activation_volume = VectorBase::create(p.mat.visc_activation_volume, nmat);
    heat_capacity = VectorBase::create(p.mat.heat_capacity, nmat);
    therm_cond = VectorBase::create(p.mat.therm_cond, nmat);
    pls0 = VectorBase::create(p.mat.pls0, nmat);
    pls1 = VectorBase::create(p.mat.pls1, nmat);
    cohesion0 = VectorBase::create(p.mat.cohesion0, nmat);
    cohesion1 = VectorBase::create(p.mat.cohesion1, nmat);
    friction_angle0 = VectorBase::create(p.mat.friction_angle0, nmat);
    friction_angle1 = VectorBase::create(p.mat.friction_angle1, nmat);
    dilation_angle0 = VectorBase::create(p.mat.dilation_angle0, nmat);
    dilation_angle1 = VectorBase::create(p.mat.dilation_angle1, nmat);

    // Hydraulic parameters
    porosity = VectorBase::create(p.mat.porosity, nmat);
    hydraulic_perm = VectorBase::create(p.mat.hydraulic_perm, nmat);
    fluid_rho0 = VectorBase::create(p.mat.fluid_rho0, nmat);
    fluid_alpha = VectorBase::create(p.mat.fluid_alpha, nmat);
    fluid_bulk_modulus = VectorBase::create(p.mat.fluid_bulk_modulus, nmat);
    fluid_visc = VectorBase::create(p.mat.fluid_visc, nmat);
    biot_coeff = VectorBase::create(p.mat.biot_coeff, nmat);
    bulk_modulus_s = VectorBase::create(p.mat.bulk_modulus_s, nmat);
    // Rate-and-state friction parameters
    direct_a = VectorBase::create(p.mat.direct_a, nmat);
    evolution_b = VectorBase::create(p.mat.evolution_b, nmat);
    characteristic_velocity = VectorBase::create(p.mat.characteristic_velocity, nmat);
    // static_friction_coefficient = VectorBase::create(p.mat.static_friction_coefficient, nmat);
}


MatProps::~MatProps()
{
    delete rho0;
    delete alpha;
    delete bulk_modulus;
    delete shear_modulus;
    delete visc_exponent;
    delete visc_coefficient;
    delete visc_activation_energy;
    delete visc_activation_volume;
    delete heat_capacity;
    delete therm_cond;
    delete pls0;
    delete pls1;
    delete cohesion0;
    delete cohesion1;
    delete friction_angle0;
    delete friction_angle1;
    delete dilation_angle0;
    delete dilation_angle1;

    // Deleting hydraulic properties
    delete porosity;
    delete hydraulic_perm;
    delete fluid_rho0;
    delete fluid_alpha;
    delete fluid_bulk_modulus;
    delete fluid_visc;
    delete biot_coeff;
    delete bulk_modulus_s;
    // Deleting rate-and-state friction properties
    delete direct_a;
    delete evolution_b;
    delete characteristic_velocity;
    // delete static_friction_coefficient;
}


double MatProps::bulkm(int e) const
{
    return harmonic_mean(*bulk_modulus, elemmarkers[e]);
}


double MatProps::shearm(int e) const
{
    return harmonic_mean(*shear_modulus, elemmarkers[e]);
}


double MatProps::visc(int e) const
{
    const double gas_constant = 8.3144;
    const double min_strain_rate = 1e-30;

    // average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    // stress
    const double *s = stress[e];
    double s0 = trace(s) / NDIMS;

    // strain-rate
    double edot = second_invariant(strain_rate[e]);
    // min strain rate to prevent viscosity -> inf
    edot = std::max(edot, min_strain_rate);

    // viscosity law from Chen and Morgan, JGR, 1990
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        double pow = 1 / (*visc_exponent)[m] - 1;
        double pow1 = -1 / (*visc_exponent)[m];
        double visc0 = 0.25 * std::pow(edot, pow) * std::pow(0.75 * (*visc_coefficient)[m], pow1)
            * std::exp(((*visc_activation_energy)[m] + (*visc_activation_volume)[m] * s0)
            / ((*visc_exponent)[m] * gas_constant * T)) * 1e6;
        result += elemmarkers[e][m] / visc0;
        n += elemmarkers[e][m];
    }

    double visc = n / result;

    // applying min & max limits
    visc = std::min(std::max(visc, visc_min), visc_max);

    return visc;
}


void MatProps::plastic_weakening(int e, double pls,
                                 double &cohesion, double &friction_angle,
                                 double &dilation_angle, double &hardening) const
{
    double c, f, d, h;
    c = f = d = h = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        if (pls < (*pls0)[m]) {
            // no weakening yet
            c += (*cohesion0)[m] * k;
            f += (*friction_angle0)[m] * k;
            d += (*dilation_angle0)[m] * k;
            h += 0;
        }
        else if (pls < (*pls1)[m]) {
            // linear weakening
            double p = (pls - (*pls0)[m]) / ((*pls1)[m] - (*pls0)[m]);
            c += ((*cohesion0)[m] + p * ((*cohesion1)[m] - (*cohesion0)[m])) * k;
            f += ((*friction_angle0)[m] + p * ((*friction_angle1)[m] - (*friction_angle0)[m])) * k;
            d += ((*dilation_angle0)[m] + p * ((*dilation_angle1)[m] -(* dilation_angle0)[m])) * k;
            h += ((*cohesion1)[m] - (*cohesion0)[m]) / ((*pls1)[m] - (*pls0)[m]) * k;
        }
        else {
            // saturated weakening
            c += (*cohesion1)[m] * k;
            f += (*friction_angle1)[m] * k;
            d += (*dilation_angle1)[m] * k;
            h += 0;
        }
    }
    cohesion = c / n;
    friction_angle = f / n;
    dilation_angle = d / n;
    hardening = h / n;
}

void MatProps::plastic_weakening_rsf(int e, double pls,
                                 double &cohesion, double &friction_angle,
                                 double &dilation_angle, double &hardening, double &slip_rate) const
{
    double c, f, d, h;
    c = f = d = h = 0;

    double d_a, e_b, c_v, mu_0, mu_d;
    d_a = e_b = c_v = mu_0 = mu_d = 0;

    double static_friction_angle = 0;

    int n = 0;
    for (int m=0; m<nmat; m++) {
        int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        if (pls < (*pls0)[m]) {
            // no weakening yet
            c += (*cohesion0)[m] * k;
            f += (*friction_angle0)[m] * k;
            d += (*dilation_angle0)[m] * k;
            h += 0;
        }
        else if (pls < (*pls1)[m]) {
            // linear weakening
            double p = (pls - (*pls0)[m]) / ((*pls1)[m] - (*pls0)[m]);
            c += ((*cohesion0)[m] + p * ((*cohesion1)[m] - (*cohesion0)[m])) * k;
            f += ((*friction_angle0)[m] + p * ((*friction_angle1)[m] - (*friction_angle0)[m])) * k;
            d += ((*dilation_angle0)[m] + p * ((*dilation_angle1)[m] -(* dilation_angle0)[m])) * k;
            h += ((*cohesion1)[m] - (*cohesion0)[m]) / ((*pls1)[m] - (*pls0)[m]) * k;
        }
        else {
            // saturated weakening
            c += (*cohesion1)[m] * k;
            f += (*friction_angle1)[m] * k;
            d += (*dilation_angle1)[m] * k;
            h += 0;
        }

        // Rate-and-state friction parameters
        d_a += (*direct_a)[m] * k;
        e_b += (*evolution_b)[m] * k;
        c_v += (*characteristic_velocity)[m] * k;
    }

    d_a /= n; // direct effect parameter
    e_b /= n; // evolution effect parameter   
    c_v /= n; // characteristic velocity
    static_friction_angle = f / n; // friction angle
    mu_0 = std::tan(DEG2RAD * static_friction_angle); // static friction coefficient
    mu_d = mu_0 + (d_a - e_b) * log(slip_rate / c_v); // dynamic friction angle
    friction_angle = std::atan(mu_d) / DEG2RAD;
    cohesion = c / n;
    dilation_angle = d / n;
    hardening = h / n;

    // std::cout << "d_a: " << d_a << " e_b: " << e_b << " c_v: " << c_v << " mu_0: " << mu_0 << " mu_d: " << mu_d << " friction_angle: " << friction_angle << std::endl;

    // std::cout << "d_a: " << d_a << std::endl;
    // std::cout << "e_b: " << e_b << std::endl;  
    // std::cout << "c_v: " << c_v << std::endl;
    // std::cout << "mu_0: " << mu_0 << std::endl;
    // std::cout << "mu_d: " << mu_d << std::endl;
    // std::cout << "friction_angle: " << friction_angle << std::endl;
}

void MatProps::plastic_props(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max) const
{
    double cohesion, phi, psi;

    plastic_weakening(e, pls, cohesion, phi, psi, hardn);

    // derived variables
    double sphi = std::sin(phi * DEG2RAD);
    double spsi = std::sin(psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? tension_max : std::min(tension_max, cohesion/std::tan(phi*DEG2RAD));
}

void MatProps::plastic_props_rsf(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max, double& slip_rate) const
{
    double cohesion, phi, psi;

    plastic_weakening_rsf(e, pls, cohesion, phi, psi, hardn, slip_rate);

    // derived variables
    double sphi = std::sin(phi * DEG2RAD);
    double spsi = std::sin(psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? tension_max : std::min(tension_max, cohesion/std::tan(phi*DEG2RAD));
}

double MatProps::rho(int e) const
{
    const double celsius0 = 273;

    // average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    double TinCelsius = T - celsius0;
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        // TODO: compressibility
        result += (*rho0)[m] * (1 - (*alpha)[m] * TinCelsius) * elemmarkers[e][m];
        n += elemmarkers[e][m];
    }
    return result / n;
}


double MatProps::cp(int e) const
{
    return arithmetic_mean(*heat_capacity, elemmarkers[e]);
}


double MatProps::k(int e) const
{
    return arithmetic_mean(*therm_cond, elemmarkers[e]);
}

// hydraulic parameters
double MatProps::phi(int e) const
{
    return arithmetic_mean(*porosity, elemmarkers[e]);
}

double MatProps::perm(int e) const
{
    return harmonic_mean(*hydraulic_perm, elemmarkers[e]);
}

double MatProps::alpha_fluid(int e) const
{
    return arithmetic_mean(*fluid_alpha, elemmarkers[e]);
}

double MatProps::beta_fluid(int e) const
{
    // Return the inverse of harmonic mean (compressibility)
    return 1.0 / harmonic_mean(*fluid_bulk_modulus, elemmarkers[e]);
}

double MatProps::rho_fluid(int e) const
{
    const double p0 = 1.0e5; // Reference pressure in Pascals (1 atm)
    const double T0 = 273; // Reference temperature in Kelvin (0°C)

    // Average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    // Average pore pressure of this element
    double p = 0;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        p += ppressure[conn[i]];
    }
    p /= NODES_PER_ELEM;

    // Use the functions to get thermal expansivity and compressibility
    double alpha_f = alpha_fluid(e);  // Thermal expansivity for the fluid
    double beta_f = beta_fluid(e);    // Compressibility for the fluid

    // Calculate fluid density based on thermal expansivity and compressibility
    double result = 0;
    int n = 0;
    for (int m = 0; m < nmat; m++) {
        double rho_f = (*fluid_rho0)[m];  // Reference fluid density
        result += rho_f * (1 - alpha_f * (T - T0) + beta_f * (p - p0)) * elemmarkers[e][m];
        n += elemmarkers[e][m];
    }

    // Return the averaged fluid density
    return result / n;
}

double MatProps::mu_fluid(int e) const
{
    return arithmetic_mean(*fluid_visc, elemmarkers[e]);
}

double MatProps::alpha_biot(int e) const
{
    return arithmetic_mean(*biot_coeff, elemmarkers[e]);
}

double MatProps::beta_mineral(int e) const
{
    // Return the inverse of harmonic mean (compressibility)
    return 1.0 / harmonic_mean(*bulk_modulus_s, elemmarkers[e]);
}

// Rate-and-state friction parameters
double MatProps::d_a(int e) const
{
    return arithmetic_mean(*direct_a, elemmarkers[e]);
}

double MatProps::e_b(int e) const
{
    return arithmetic_mean(*evolution_b, elemmarkers[e]);
}

double MatProps::c_v(int e) const
{
    return arithmetic_mean(*characteristic_velocity, elemmarkers[e]);
}

#endif