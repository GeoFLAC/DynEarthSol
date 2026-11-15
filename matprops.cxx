#include <algorithm>
#include <cmath>
#include <limits>

#include "constants.hpp"
#include "utils.hpp"
#include "matprops.hpp"


namespace {

    #pragma acc routine seq
    inline double rsf_slip_rate_eff(double slip_rate)
    {
        const double cutoff = 1e-16;
        return (slip_rate > cutoff) ? slip_rate : cutoff;
    }

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
                            prem_max_depth };

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

        // Deeper than prem_max_depth (or for a NaN depth) no row matches, n stays nlayers and the
        // interpolation reads one past both tables. mesh.zlength is checked against the limit at
        // input, but a run whose nodes later move deeper than that, such as a diverged one,
        // overruns here.
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
                            prem_max_depth };

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

        // Deeper than prem_max_depth (or for a NaN depth) no row matches, n stays nlayers and the
        // interpolation reads one past both tables. mesh.zlength is checked against the limit at
        // input, but a run whose nodes later move deeper than that, such as a diverged one,
        // overruns here.
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


    #pragma acc routine seq
    double arithmetic_mean(const double_vec &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        #pragma acc loop seq
        for (std::size_t i=0; i<s.size(); i++) {
            // Absent materials add 0*s[i] == 0 and 0 to m, so skipping is exact.
            if (n[i] == 0) continue;
            result += n[i] * s[i];
            m += n[i];
        }
        return result / m;
    }


    #pragma acc routine seq
    double harmonic_mean(const double_vec &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        #pragma acc loop seq
        for (std::size_t i=0; i<s.size(); i++) {
            // Absent materials add 0/s[i] == 0 and 0 to m, so skipping is exact --
            // and it avoids the 0/0 = NaN an absent material with s[i] == 0 gave.
            if (n[i] == 0) continue;
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
    double p = 0;

    if (param.control.ref_pressure_option == 0) {
        if (param.control.has_hydraulic_diffusion) {
        // Modified density considering porosity for hydraulic diffusion
            p = (param.mat.rho0[param.mat.mattype_ref] * (1 - param.mat.porosity[param.mat.mattype_ref]) + \
                1000.0 * param.mat.porosity[param.mat.mattype_ref]) * param.control.gravity * depth;
        } else {
            // Standard reference pressure without hydraulic diffusion
            p = param.mat.rho0[param.mat.mattype_ref] * param.control.gravity * depth;
        }
    }
    else if (param.control.ref_pressure_option == 1)
        p = get_prem_pressure(depth);
    else if (param.control.ref_pressure_option == 2)
        p = get_prem_pressure_modified(depth);
    return p;
}


// Only an NPROF build creates the lookup tables, and only there do the *_safe functions read them;
// elsewhere they are null, and binding a reference to null is undefined even if it is never read.
static const double_vec& table_or_empty(const double_vec *table)
{
    static const double_vec empty;
    return table ? *table : empty;
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
  log_table(table_or_empty(var.log_table)),
  tan_table(table_or_empty(var.tan_table)),
  sin_table(table_or_empty(var.sin_table))
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
    cohesion_slope.assign(nmat, 0);
    friction_slope.assign(nmat, 0);
    dilation_slope.assign(nmat, 0);
    // param_at: a model without plasticity may leave these vectors empty.
    for (int m = 0; m < nmat; ++m) {
        const double span = param_at(pls1, m) - param_at(pls0, m);
        if (!(span > 0)) continue;   // no weakening segment
        cohesion_slope[m] = (param_at(cohesion1, m) - param_at(cohesion0, m)) / span;
        friction_slope[m] = (param_at(friction_angle1, m) - param_at(friction_angle0, m)) / span;
        dilation_slope[m] = (param_at(dilation_angle1, m) - param_at(dilation_angle0, m)) / span;
    }

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

    // Material-only terms of the creep law, in the association visc() evaluates them in
    // so the per-element result is bit-identical at opt <= 2. Round-off-equivalent only
    // under -ffast-math and on openacc=1 (device -> host libm).
    {
        // J mol^-1 K^-1, as the Chen & Morgan creep law uses it.
        const double gas_constant = 8.3144;
        const int nm = nmat;   // == visc_exponent.size(); nmat is the bound visc() loops to
        visc_pow_edot.resize(nm);
        visc_coef_term.resize(nm);
        visc_nR.resize(nm);
        for (int m = 0; m < nm; ++m) {
            visc_pow_edot[m] = 1 / visc_exponent[m] - 1;
            const double pow1 = -1 / visc_exponent[m];
            visc_coef_term[m] = pow_safe(log_table, 0.75 * visc_coefficient[m], pow1);
            visc_nR[m] = visc_exponent[m] * gas_constant;
        }
    }

    refresh_elem_cache();
}


/* Rebuild the per-element property means only if elemmarkers moved; they are a pure
 * function of it. Every call site sits right after something that moves elemmarkers;
 * a missed one is silent staleness. */
void MatProps::refresh_elem_cache()
{
    const int nelem = elemmarkers.size();

    // A size change means a new mesh: nothing about the old contents survives it.
    const bool resized = int(cache_props.size()) != nelem * MP_STRIDE;
    if (resized) {
        cache_props.resize(std::size_t(nelem) * MP_STRIDE);
        cache_markers.assign(std::size_t(nelem) * nmat, -1);
    }

    // Compare while copying. `changed` tracks whether the CONTENT moved: on most steps
    // no marker moves and the whole rebuild below is skipped.
    int changed = resized ? 1 : 0;
#ifndef ACC
    #pragma omp parallel for default(none) firstprivate(nelem) reduction(|:changed)
#endif
    // NO `async`: the host reads `changed` immediately below, and appending async would
    // read the reduction before it lands and silently freeze the means.
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector reduction(|:changed)
    for (int e = 0; e < nelem; ++e) {
        for (int m = 0; m < nmat; ++m) {
            const std::size_t i = std::size_t(e) * nmat + m;
            const int v = elemmarkers[e][m];
            if (cache_markers[i] != v) { cache_markers[i] = v; changed = 1; }
        }
    }
    if (!changed) return;

#ifndef ACC
    #pragma omp parallel for default(none) firstprivate(nelem)
#endif
    #pragma acc parallel loop if(des_acc_offload_enabled()) gang vector
    for (int e = 0; e < nelem; ++e) {
        double* c = &cache_props[std::size_t(e) * MP_STRIDE];
        c[MPC_BULKM]      = harmonic_mean(bulk_modulus, elemmarkers[e]);
        c[MPC_SHEARM]     = harmonic_mean(shear_modulus, elemmarkers[e]);
        c[MPC_PHI]        = arithmetic_mean(porosity, elemmarkers[e]);
        c[MPC_ALPHA_BIOT] = arithmetic_mean(biot_coeff, elemmarkers[e]);
        c[MPC_CP]         = arithmetic_mean(heat_capacity, elemmarkers[e]);
        c[MPC_K]          = arithmetic_mean(therm_cond, elemmarkers[e]);
        // Stored as the compressibility, so beta_fluid() is a plain load like the others.
        c[MPC_BETA_FLUID] = 1.0 / harmonic_mean(fluid_bulk_modulus, elemmarkers[e]);
    }
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
    return cache_props[std::size_t(e) * MP_STRIDE + MPC_BULKM];
}


double MatProps::shearm(int e) const
{
    return cache_props[std::size_t(e) * MP_STRIDE + MPC_SHEARM];
}


double MatProps::visc(int e) const
{
    const double min_strain_rate = 1e-30;

    // average temperature of this element
    double T = 0;
    ConstConnAccessor conn = connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    // stress
    ConstTensorAccessor s = stress[e];
    double s0 = trace(s) / NDIMS;

    // strain-rate
    double edot = second_invariant(strain_rate[e]);
    // min strain rate to prevent viscosity -> inf
    edot = std::max(edot, min_strain_rate);

    // viscosity law from Chen and Morgan, JGR, 1990
    double result = 0;
    int n = 0;

    for (int m=0; m<nmat; m++) {
        // Absent materials add 0/visc0 == 0, so the skip is exact, and it keeps
        // visc0 == 0 from contributing a NaN.
        const int marker_count = elemmarkers[e][m];
        if (marker_count == 0) continue;

        double visc0 = 0.25 * pow_safe(log_table,edot, visc_pow_edot[m]) * visc_coef_term[m]
            * std::exp((visc_activation_energy[m] + visc_activation_volume[m] * s0)
            / (visc_nR[m] * T)) * 1e6;
        result += marker_count / visc0;
        n += marker_count;
    }

    double visc = n / result;

    // applying min & max limits
    visc = std::min(std::max(visc, visc_min), visc_max);

    return visc;
}


// Dynamic friction coefficient from the static friction angle and the rate
// and state terms.
static double rsf_coefficient(double static_friction_angle, const double terms[2])
{
    const double mu_0 = std::tan(DEG2RAD * static_friction_angle);
    return std::max(mu_0 + terms[0] + terms[1], 1e-6);
}

// The material vectors hold nmat entries (get_numbers expands one value
// for all), so they are indexed directly here, on every element's path.
void MatProps::weakening_means(int e, double pls,
                               double &cohesion, double &friction_angle,
                               double &dilation_angle, double &hardening, int &n) const
{
    double c, f, d, h;
    c = f = d = h = 0;
    n = 0;
    for (int m=0; m<nmat; m++) {
        const int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        if (pls < pls0[m]) {
            // no weakening yet
            c += cohesion0[m] * k;
            f += friction_angle0[m] * k;
            d += dilation_angle0[m] * k;
        }
        else if (pls < pls1[m]) {
            // linear weakening
            const double p = (pls - pls0[m]) / (pls1[m] - pls0[m]);
            c += (cohesion0[m] + p * (cohesion1[m] - cohesion0[m])) * k;
            f += (friction_angle0[m] + p * (friction_angle1[m] - friction_angle0[m])) * k;
            d += (dilation_angle0[m] + p * (dilation_angle1[m] - dilation_angle0[m])) * k;
            h += cohesion_slope[m] * k;
        }
        else {
            // saturated weakening
            c += cohesion1[m] * k;
            f += friction_angle1[m] * k;
            d += dilation_angle1[m] * k;
        }
    }
    if (n == 0) {
        cohesion = friction_angle = dilation_angle = hardening = 0;
        return;
    }
    cohesion = c / n;
    friction_angle = f / n;
    dilation_angle = d / n;
    hardening = h / n;
}

void MatProps::plastic_weakening(int e, double pls,
                                 double &cohesion, double &friction_angle,
                                 double &dilation_angle, double &hardening) const
{
    int n;
    weakening_means(e, pls, cohesion, friction_angle, dilation_angle, hardening, n);
}

// The piece's angle slopes and end, by the same branches as weakening_means.
void MatProps::complete_piece(int e, StrengthPiece& piece) const
{
    if (piece.has_slopes) return;
    const double pls = piece.pls;
    double df = 0, dd = 0, end = std::numeric_limits<double>::max();
    int n = 0;
    for (int m=0; m<nmat; m++) {
        const int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        if (pls < pls0[m])
            end = std::min(end, pls0[m]);
        else if (pls < pls1[m]) {
            df += friction_slope[m] * k;
            dd += dilation_slope[m] * k;
            end = std::min(end, pls1[m]);
        }
    }
    piece.end = end;
    piece.dfriction = n ? df / n : 0;
    piece.ddilation = n ? dd / n : 0;
    piece.has_slopes = true;
}

void MatProps::strength_here(int e, double pls, StrengthPiece& piece) const
{
    int n;
    weakening_means(e, pls, piece.cohesion, piece.friction, piece.dilation,
                    piece.dcohesion, n);
    piece.pls = pls;
    piece.rsf = 0;
    piece.rsf_terms[0] = piece.rsf_terms[1] = 0;
    // complete_piece sets these once the element yields; until then they hold
    // values, so that copying the piece copies no indeterminate ones.
    piece.end = std::numeric_limits<double>::max();
    piece.dfriction = piece.ddilation = 0;
    piece.has_slopes = false;
    finish_piece(piece);
}

void MatProps::strength_piece(int e, double pls, StrengthPiece& piece) const
{
    strength_here(e, pls, piece);
    complete_piece(e, piece);
}

// The angle factors and the strength at the start of the piece.
void MatProps::finish_piece(StrengthPiece& piece) const
{
    double phi = piece.friction;
    if (piece.rsf == 1) phi = std::atan(rsf_coefficient(phi, piece.rsf_terms)) / DEG2RAD;
    angle_factors(phi, piece.dilation, piece.rsf != 0, piece.factors);
    strength_from(piece.cohesion, piece.factors,
                  piece.amc, piece.anphi, piece.anpsi, piece.ten_max);
}

void MatProps::rsf_terms(int e, double slip_rate, double state_variable, int state_model,
                         double terms[2]) const
{
    double d_a_avg = 0.0;
    double e_b_avg = 0.0;
    double c_v_avg = 0.0;
    double d_c_avg = 0.0;
    int n = 0;
    for (int m = 0; m < nmat; ++m) {
        const int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        d_a_avg += param_at(direct_a, m) * k;
        e_b_avg += param_at(evolution_b, m) * k;
        c_v_avg += param_at(characteristic_velocity, m) * k;
        d_c_avg += param_at(characteristic_distance, m) * k;
    }
    d_a_avg /= n;
    e_b_avg /= n;
    c_v_avg /= n;
    d_c_avg /= n;

    const double v_eff = rsf_slip_rate_eff(slip_rate);
    const double cv_eff = std::max(c_v_avg, 1e-30);
    const double dc_eff = std::max(d_c_avg, 1e-30);
    const double theta_eff = std::max(state_variable, 1e-30);

    if (state_model == 0) {
        terms[0] = (d_a_avg - e_b_avg) * std::log(v_eff / cv_eff);
        terms[1] = 0;
    } else {
        const double log_v = std::log(v_eff / cv_eff);
        const double log_theta = std::log((cv_eff * theta_eff) / dc_eff);
        terms[0] = d_a_avg * log_v;
        terms[1] = e_b_avg * log_theta;
    }
}

void MatProps::plastic_weakening_rsf(int e, double pls,
                                 double &cohesion, double &dynamic_friction_angle,
                                 double &dilation_angle, double &hardening, double &slip_rate,
                                 double& dyn_fric_coeff, double& state_variable,
                                 int state_model, double dt, bool frozen_history) const
{
    (void) dt;

    double static_friction_angle;
    int n;
    weakening_means(e, pls, cohesion, static_friction_angle, dilation_angle, hardening, n);
    if (n == 0) {
        dynamic_friction_angle = 0.0;
        dyn_fric_coeff = 0.0;
        return;
    }

    if (frozen_history) {
        // Initialization solves an instantaneous skeleton correction. Its
        // numerical correction velocity is not a physical RSF slip rate.
        dynamic_friction_angle = std::atan(dyn_fric_coeff) / DEG2RAD;
        return;
    }

    double terms[2];
    rsf_terms(e, slip_rate, state_variable, state_model, terms);
    const double mu_d = rsf_coefficient(static_friction_angle, terms);
    dynamic_friction_angle = std::atan(mu_d) / DEG2RAD;
    dyn_fric_coeff = mu_d;
}

void MatProps::rsf_strength_piece(int e, double pls, double slip_rate, double state_variable,
                                  double dyn_fric_coeff, int state_model, bool frozen_history,
                                  StrengthPiece& piece) const
{
    int n;
    weakening_means(e, pls, piece.cohesion, piece.friction, piece.dilation,
                    piece.dcohesion, n);
    piece.pls = pls;
    piece.has_slopes = false;
    complete_piece(e, piece);
    if (frozen_history) {
        // The stored coefficient, as in plastic_weakening_rsf.
        piece.rsf = 2;
        piece.friction = std::atan(dyn_fric_coeff) / DEG2RAD;
        piece.dfriction = 0;
    }
    else {
        piece.rsf = 1;
        rsf_terms(e, slip_rate, state_variable, state_model, piece.rsf_terms);
    }
    finish_piece(piece);
}

void MatProps::piece_point(const StrengthPiece& piece, double offset,
                           double& cohesion, AngleFactors& factors) const
{
    cohesion = piece.cohesion + piece.dcohesion*offset;
    if (!piece.angles_change()) {
        factors = piece.factors;
        return;
    }
    double phi = piece.friction + piece.dfriction*offset;
    if (piece.rsf == 1) phi = std::atan(rsf_coefficient(phi, piece.rsf_terms)) / DEG2RAD;
    angle_factors(phi, piece.dilation + piece.ddilation*offset, piece.rsf != 0, factors);
}

void MatProps::piece_strength(const StrengthPiece& piece, double offset,
                              double& amc, double& anphi, double& anpsi, double& ten_max) const
{
    double cohesion;
    AngleFactors factors;
    piece_point(piece, offset, cohesion, factors);
    strength_from(cohesion, factors, amc, anphi, anpsi, ten_max);
}

// Cohesion changes linearly along a piece and so does the friction angle, which the
// rate-and-state terms raise monotonically; anphi, sqrt_anphi and tan_phi rise with it. Each
// therefore takes its extremes at x and y, and the strength those of their products.
//
// Where the friction angle phi itself is linear, at dphi per unit plastic strain, u =
// sqrt_anphi = tan(pi/4 + phi/2) and w = 1/tan_phi change by u' = dphi (1 + u^2)/2 and w' =
// -dphi (1 + w^2), so that
//   anphi'' = dphi^2 (1 + u^2)(1 + 3u^2)/2,
//   amc''   = (1 + u^2)(2 dcohesion dphi + cohesion dphi^2 u),
//   (cohesion w)'' = (1 + w^2)(2 cohesion dphi^2 w - 2 dcohesion dphi),
// each largest in size at the largest u, w and cohesion. The lookup tables of the GPU build
// interpolate sin and tan linearly, within h^2/8 of sin and h^2 tan (1 + tan^2)/4 of tan for
// a table step h, and anphi, sqrt_anphi and cohesion/tan_phi follow them through their
// derivatives, bounded here twice over.
void MatProps::piece_strength_range(const StrengthPiece& piece, double x, double y,
                                    StrengthRange& range) const
{
    double c[2];
    AngleFactors f[2];
    piece_point(piece, x, c[0], f[0]);
    piece_point(piece, y, c[1], f[1]);
    const double c_lo = std::min(c[0], c[1]), c_hi = std::max(c[0], c[1]);
    const double t_lo = std::min(f[0].tan_phi, f[1].tan_phi);
    const double t_hi = std::max(f[0].tan_phi, f[1].tan_phi);
    const double u_lo = std::min(f[0].sqrt_anphi, f[1].sqrt_anphi);
    const double u_hi = std::max(f[0].sqrt_anphi, f[1].sqrt_anphi);
    range.anphi[0] = std::min(f[0].anphi, f[1].anphi);
    range.anphi[1] = std::max(f[0].anphi, f[1].anphi);
    range.amc[0] = 2 * c_lo * u_lo;
    range.amc[1] = 2 * c_hi * u_hi;
    // As strength_from: cohesion/tan_phi up to tension_max, which zero friction (tan_phi 0) gives.
    range.ten_max[0] = t_hi > 0 ? std::min(tension_max, c_lo / t_hi) : tension_max;
    range.ten_max[1] = t_lo > 0 ? std::min(tension_max, c_hi / t_lo) : tension_max;

    range.curvature_amc = range.curvature_anphi = range.curvature_ten = -1;
    range.error_amc = range.error_anphi = range.error_ten = 0;
    if (piece.rsf == 1 && piece.dfriction != 0) return;
    const double dphi = (piece.rsf == 2 ? 0 : piece.dfriction) * DEG2RAD;
    const double dc = std::fabs(piece.dcohesion * dphi);
    range.curvature_anphi = 0.5 * dphi * dphi * (1 + u_hi * u_hi) * (1 + 3 * u_hi * u_hi);
    range.curvature_amc = (1 + u_hi * u_hi) * (2 * dc + c_hi * dphi * dphi * u_hi);
    if (t_lo > 0) {
        const double w = 1 / t_lo;
        range.curvature_ten = (1 + w * w) * 2 * (dc + c_hi * dphi * dphi * w);
    }
#ifdef NPROF
    if (piece.rsf != 0) return;
    const double sin_error = SIN_XDELTA * SIN_XDELTA / 8;
    const double tan_error = TAN_XDELTA * TAN_XDELTA * t_hi * (1 + t_hi * t_hi) / 4;
    // d anphi/d sin = (anphi + 1)^2/2, d sqrt_anphi = d anphi/(2 sqrt_anphi).
    range.error_anphi = (range.anphi[1] + 1) * (range.anphi[1] + 1) * sin_error;
    range.error_amc = 2 * c_hi * range.error_anphi / u_lo;
    if (t_lo > 0) range.error_ten = 2 * c_hi * tan_error / (t_lo * t_lo);
#endif
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

            const double slip_rate_eff = rsf_slip_rate_eff(slip_rate);
            const double dtheta = (1.0 - (slip_rate_eff * state_variable / d)) * dt;
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

            const double slip_rate_eff = rsf_slip_rate_eff(slip_rate);
            double ratio = slip_rate_eff * theta / d;
            if (ratio < ratio_min) ratio = ratio_min;

            const double dtheta = -ratio * std::log(ratio) * dt;
            if (std::isfinite(dtheta)) {
                double new_theta = theta + dtheta;
                if (!std::isfinite(new_theta) || new_theta <= 0.0) {
                    state_variable = d / slip_rate_eff;
                } else {
                    if (new_theta < theta_min) new_theta = theta_min;
                    if (new_theta > theta_max) new_theta = theta_max;
                    state_variable = new_theta;
                }
            } else {
                state_variable = d / slip_rate_eff;
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
    ep_derived(cohesion, phi, psi, amc, anphi, anpsi, ten_max);
}

void MatProps::ep_derived(double cohesion, double phi, double psi,
                          double& amc, double& anphi, double& anpsi, double& ten_max) const
{
    AngleFactors factors;
    angle_factors(phi, psi, false, factors);
    strength_from(cohesion, factors, amc, anphi, anpsi, ten_max);
}

void MatProps::angle_factors(double phi, double psi, bool rsf, AngleFactors& factors) const
{
    const double sphi = rsf ? std::sin(phi * DEG2RAD) : sin_safe(sin_table,phi * DEG2RAD);
    const double spsi = rsf ? std::sin(psi * DEG2RAD) : sin_safe(sin_table,psi * DEG2RAD);

    factors.anphi = (1 + sphi) / (1 - sphi);
    factors.anpsi = (1 + spsi) / (1 - spsi);
    factors.sqrt_anphi = std::sqrt(factors.anphi);
    factors.zero_friction = phi == 0;
    factors.tan_phi = factors.zero_friction ? 0 :
        rsf ? std::tan(phi*DEG2RAD) : tan_safe(tan_table,phi*DEG2RAD);
}

void MatProps::strength_from(double cohesion, const AngleFactors& factors,
                             double& amc, double& anphi, double& anpsi, double& ten_max) const
{
    anphi = factors.anphi;
    anpsi = factors.anpsi;
    amc = 2 * cohesion * factors.sqrt_anphi;
    ten_max = factors.zero_friction ? tension_max : std::min(tension_max, cohesion/factors.tan_phi);
}

void MatProps::plastic_props_rsf(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max, double& slip_rate,
                             double& dyn_fric_coeff, double& state_variable,
                             double dt, int state_model, bool frozen_history) const
{
    double cohesion, phi, psi;

    if (!frozen_history) update_state_variable(e, slip_rate, state_variable, dt, state_model);
    plastic_weakening_rsf(e, pls, cohesion, phi, psi, hardn, slip_rate,
                          dyn_fric_coeff, state_variable, state_model, dt, frozen_history);
    rsf_derived(cohesion, phi, psi, amc, anphi, anpsi, ten_max);
}

void MatProps::rsf_derived(double cohesion, double phi, double psi,
                           double& amc, double& anphi, double& anpsi, double& ten_max) const
{
    AngleFactors factors;
    angle_factors(phi, psi, true, factors);
    strength_from(cohesion, factors, amc, anphi, anpsi, ten_max);
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
    ConstConnAccessor conn = connectivity[e];
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
    return cache_props[std::size_t(e) * MP_STRIDE + MPC_CP];
}


double MatProps::k(int e) const
{
    return cache_props[std::size_t(e) * MP_STRIDE + MPC_K];
}

// hydraulic parameters
double MatProps::phi(int e) const
{
    return cache_props[std::size_t(e) * MP_STRIDE + MPC_PHI];
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
    // Already stored as the compressibility; see refresh_elem_cache().
    return cache_props[std::size_t(e) * MP_STRIDE + MPC_BETA_FLUID];
}

double MatProps::rho_fluid(int e) const
{
    const double p0 = 1.0e5; // Reference pressure in Pascals (1 atm)
    const double T0 = 273; // Reference temperature in Kelvin (0°C)

    // Average temperature of this element
    double T = 0;
    ConstConnAccessor conn = connectivity[e];
    #pragma acc loop seq
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    // Average pore pressure of this element
    double p = 0;
    #pragma acc loop seq
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
    #pragma acc loop seq
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
    return cache_props[std::size_t(e) * MP_STRIDE + MPC_ALPHA_BIOT];
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
