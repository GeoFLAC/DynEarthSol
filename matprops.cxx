#include <algorithm>
#include <cmath>

#include "constants.hpp"
#include "utils.hpp"
#include "matprops.hpp"


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

    #pragma acc routine seq
    static inline double param_at(const double_vec &v, int m)
    {
        const std::size_t n = v.size();
        if (n == 0) return 0.0;
        if (n == 1) return v[0];
        std::size_t idx = static_cast<std::size_t>(m);
        if (idx >= n) idx = n - 1;
        return v[idx];
    }


    double arithmetic_mean(const double_vec &s, const int_vec &n)
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


    double harmonic_mean(const double_vec &s, const int_vec &n)
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

    if (param.control.ref_pressure_option == 0)
        if (param.control.has_hydraulic_diffusion) {
        // Modified density considering porosity for hydraulic diffusion
            p = (param.mat.rho0[param.mat.mattype_ref] * (1 - param.mat.porosity[param.mat.mattype_ref]) + \
                1000.0 * param.mat.porosity[param.mat.mattype_ref]) * param.control.gravity * depth;
        } else {
            // Standard reference pressure without hydraulic diffusion
            p = param.mat.rho0[param.mat.mattype_ref] * param.control.gravity * depth;
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
  dppressure(*var.dppressure),
  log_table(*var.log_table),
  tan_table(*var.tan_table),
  sin_table(*var.sin_table)
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
    characteristic_distance = p.mat.characteristic_distance;
    // static_friction_coefficient = p.mat.static_friction_coefficient;
}


MatProps::~MatProps()
{
    #pragma acc exit data delete(rho0,alpha,bulk_modulus,shear_modulus,visc_exponent)
    #pragma acc exit data delete(visc_coefficient,visc_activation_energy,visc_activation_volume)
    #pragma acc exit data delete(heat_capacity,therm_cond,pls0,pls1,cohesion0)
    #pragma acc exit data delete(friction_angle0,friction_angle1,dilation_angle0,dilation_angle1)

    // Deleting hydraulic properties
    #pragma acc exit data delete(porosity,hydraulic_perm,fluid_rho0,fluid_alpha,fluid_bulk_modulus,fluid_visc,biot_coeff,bulk_modulus_s)
    // Deleting rate-and-state friction properties
    #pragma acc exit data delete(direct_a,evolution_b,characteristic_velocity,characteristic_distance)
    // #pragma acc exit data delete(static_friction_coefficient)
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

    for (int m=0; m<nmat; m++) {
        double pow = 1 / visc_exponent[m] - 1;
        double pow1 = -1 / visc_exponent[m];
        double visc0 = 0.25 * pow_safe(log_table,edot, pow) * pow_safe(log_table, 0.75 * visc_coefficient[m], pow1)
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
                                 double &cohesion, double &dynamic_friction_angle,
                                 double &dilation_angle, double &hardening, double &slip_rate,
                                 double& dyn_fric_coeff, double& state_variable,
                                 int state_model, double dt) const
{
    (void) dt;

    double c = 0.0;
    double f = 0.0;
    double d = 0.0;
    double h = 0.0;
    double d_a_avg = 0.0;
    double e_b_avg = 0.0;
    double c_v_avg = 0.0;
    double d_c_avg = 0.0;
    int n = 0;

    for (int m = 0; m < nmat; ++m) {
        const int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;

        const double pls0_m = param_at(pls0, m);
        const double pls1_m = param_at(pls1, m);
        const double cohesion0_m = param_at(cohesion0, m);
        const double cohesion1_m = param_at(cohesion1, m);
        const double friction0_m = param_at(friction_angle0, m);
        const double friction1_m = param_at(friction_angle1, m);
        const double dilation0_m = param_at(dilation_angle0, m);
        const double dilation1_m = param_at(dilation_angle1, m);

        if (pls < pls0_m) {
            c += cohesion0_m * k;
            f += friction0_m * k;
            d += dilation0_m * k;
        }
        else if (pls < pls1_m) {
            const double denom = std::max(pls1_m - pls0_m, 1e-30);
            const double p = (pls - pls0_m) / denom;
            c += (cohesion0_m + p * (cohesion1_m - cohesion0_m)) * k;
            f += (friction0_m + p * (friction1_m - friction0_m)) * k;
            d += (dilation0_m + p * (dilation1_m - dilation0_m)) * k;
            h += (cohesion1_m - cohesion0_m) / denom * k;
        }
        else {
            c += cohesion1_m * k;
            f += friction1_m * k;
            d += dilation1_m * k;
        }

        d_a_avg += param_at(direct_a, m) * k;
        e_b_avg += param_at(evolution_b, m) * k;
        c_v_avg += param_at(characteristic_velocity, m) * k;
        d_c_avg += param_at(characteristic_distance, m) * k;
    }

    if (n == 0) {
        cohesion = 0.0;
        dynamic_friction_angle = 0.0;
        dilation_angle = 0.0;
        hardening = 0.0;
        dyn_fric_coeff = 0.0;
        return;
    }

    d_a_avg /= n;
    e_b_avg /= n;
    c_v_avg /= n;
    d_c_avg /= n;

    const double static_friction_angle = f / n;
    const double mu_0 = std::tan(DEG2RAD * static_friction_angle);

    const double v_eff = std::max(slip_rate, 1e-30);
    const double cv_eff = std::max(c_v_avg, 1e-30);
    const double dc_eff = std::max(d_c_avg, 1e-30);
    const double theta_eff = std::max(state_variable, 1e-30);

    double mu_d;
    if (state_model == 0) {
        mu_d = mu_0 + (d_a_avg - e_b_avg) * std::log(v_eff / cv_eff);
    } else {
        const double log_v = std::log(v_eff / cv_eff);
        const double log_theta = std::log((cv_eff * theta_eff) / dc_eff);
        mu_d = mu_0 + d_a_avg * log_v + e_b_avg * log_theta;
    }
    mu_d = std::max(mu_d, 1e-6);

    dynamic_friction_angle = std::atan(mu_d) / DEG2RAD;
    cohesion = c / n;
    dilation_angle = d / n;
    hardening = h / n;
    dyn_fric_coeff = mu_d;
}

void MatProps::update_state_variable(int e, double slip_rate, double& state_variable,
                                     double dt, int state_model) const
{
    const double theta_min = 1e-12;
    const double theta_max = 1e12;
    const double ratio_min = 1e-10;

    switch (state_model) {
    case 0:
        // Steady-state model: keep theta fixed.
        break;
    case 1:
        {
            double d = 0.0;
            int n = 0;
            for (int m = 0; m < nmat; ++m) {
                const int k = elemmarkers[e][m];
                if (k == 0) continue;
                d += param_at(characteristic_distance, m) * k;
                n += k;
            }
            d = (n > 0) ? (d / n) : 0.0;
            if (d < 1e-12) return;

            const double dtheta = (1.0 - (slip_rate * state_variable / d)) * dt;
            if (!std::isfinite(dtheta)) return;

            double new_theta = state_variable + dtheta;
            if (new_theta < theta_min) new_theta = theta_min;
            if (new_theta > theta_max) new_theta = theta_max;
            state_variable = new_theta;
        }
        break;
    case 2:
        {
            double d = 0.0;
            int n = 0;
            for (int m = 0; m < nmat; ++m) {
                const int k = elemmarkers[e][m];
                if (k == 0) continue;
                d += param_at(characteristic_distance, m) * k;
                n += k;
            }
            d = (n > 0) ? (d / n) : 0.0;
            if (d < 1e-12) return;

            double theta = state_variable;
            if (theta < theta_min) theta = theta_min;
            if (theta > theta_max) theta = theta_max;

            double ratio = slip_rate * theta / d;
            if (ratio < ratio_min) ratio = ratio_min;

            const double dtheta = -ratio * std::log(ratio) * dt;
            if (std::isfinite(dtheta)) {
                double new_theta = theta + dtheta;
                if (!std::isfinite(new_theta) || new_theta <= 0.0) {
                    state_variable = d / std::max(slip_rate, 1e-30);
                } else {
                    if (new_theta < theta_min) new_theta = theta_min;
                    if (new_theta > theta_max) new_theta = theta_max;
                    state_variable = new_theta;
                }
            } else {
                state_variable = d / std::max(slip_rate, 1e-30);
            }
        }
        break;
    default:
        break;
    }
}

void MatProps::plastic_props(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max) const
{
    double cohesion, phi, psi;

    plastic_weakening(e, pls, cohesion, phi, psi, hardn);

    // derived variables
    double sphi = sin_safe(sin_table,phi * DEG2RAD);
    double spsi = sin_safe(sin_table,psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? tension_max : std::min(tension_max, cohesion/tan_safe(tan_table,phi*DEG2RAD));
}

void MatProps::plastic_props_rsf(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max, double& slip_rate,
                             double& dyn_fric_coeff, double& state_variable,
                             double dt, int state_model) const
{
    double cohesion, phi, psi;

    update_state_variable(e, slip_rate, state_variable, dt, state_model);
    plastic_weakening_rsf(e, pls, cohesion, phi, psi, hardn, slip_rate,
                          dyn_fric_coeff, state_variable, state_model, dt);

    // derived variables
    double sphi = std::sin(phi * DEG2RAD);
    double spsi = std::sin(psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? tension_max : std::min(tension_max, cohesion/std::tan(phi*DEG2RAD));
}

void MatProps::rsf_friction_from_state(int e, double pls, double slip_rate,
                                       double state_variable, double& dyn_fric_coeff,
                                       int state_model) const
{
    double cohesion, phi, psi, hardn;
    double theta = state_variable;
    plastic_weakening_rsf(e, pls, cohesion, phi, psi, hardn, slip_rate,
                          dyn_fric_coeff, theta, state_model, 0.0);
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
    return rho;
}


double MatProps::cp(int e) const
{
    return arithmetic_mean(heat_capacity, elemmarkers[e]);
}


double MatProps::k(int e) const
{
    return arithmetic_mean(therm_cond, elemmarkers[e]);
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

double MatProps::d_c(int e) const
{
    return arithmetic_mean(characteristic_distance, elemmarkers[e]);
}
