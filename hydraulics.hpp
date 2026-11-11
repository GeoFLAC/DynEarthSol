#ifndef DYNEARTHSOL3D_HYDRAULICS_HPP
#define DYNEARTHSOL3D_HYDRAULICS_HPP

#include <cmath>

#include "parameters.hpp"
#include "matprops.hpp"

// Constant-compressibility fluid density; rho0 is referenced to p=0 and T0.
#pragma acc routine seq
inline double fluid_density(const MatProps &mat, int e, double p, double t, double t0)
{
    return mat.rho_fluid_reference(e)*std::exp(mat.beta_fluid(e)*p-mat.alpha_fluid(e)*(t-t0));
}

struct FluidState {
    double porosity, density;
    double pressure_slope, temperature_slope; // derivatives of mass per bulk volume
};

// The saved element porosity anchors pressure/temperature at their element
// means. Nodal storage samples this same law at each quadrature pressure/T.
#pragma acc routine seq
inline FluidState fluid_state(const MatProps &mat, int e, double p, double t, double t0,
                              double phi, double reference_p, double reference_t,
                              double volume_strain = 0, double recovery = 0)
{
    const double alpha = mat.alpha_biot(e);
    const double beta_s = mat.beta_mineral(e), beta_f = mat.beta_fluid(e);
    const double alpha_s = mat.alpha_solid(e), alpha_f = mat.alpha_fluid(e);
    const double increment = -volume_strain-beta_s*(p-reference_p)+
                             alpha_s*(t-reference_t)+recovery;
    const double next = phi-(alpha-phi)*std::expm1(increment);
    const double density = fluid_density(mat,e,p,t,t0);
    return {next,density,density*(next*beta_f+(alpha-next)*beta_s),
                        -density*(next*alpha_f+(alpha-next)*alpha_s)};
}

// Compliance for the recorded mean total stress: in-plane mean in 2D,
// three-dimensional mean in 3D. Both relate that mean to volumetric strain.
#pragma acc routine seq
inline double hydraulic_mean_stress_compliance(const MatProps &mat, int e)
{
    double bulk_comp = 1.0 / mat.bulkm(e);
    if (NDIMS == 2) bulk_comp = 1.0 / (mat.bulkm(e) + mat.shearm(e) / 3.0);
    return bulk_comp;
}

// Biot fluid content is zeta = alpha*eps_v + p/M, with tension-positive strain.
#pragma acc routine seq
inline double hydraulic_inverse_biot_modulus(const MatProps &mat, int e)
{
    const double phi = mat.phi(e);
    return phi * mat.beta_fluid(e) + (mat.alpha_biot(e) - phi) * mat.beta_mineral(e);
}

// Temperature baseline belongs to the physical step, not to constitutive trials.
#pragma acc routine seq
inline double hydraulic_temperature_change(const Variables &var, int e)
{
    if (!var.old_temperature && !var.pending_temperature) return 0;
    double change = 0;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const int n = (*var.connectivity)[e][i];
        change += (var.pending_temperature ? (*var.pending_temperature)[n] :
                   (*var.temperature)[n]-(*var.old_temperature)[n])/NODES_PER_ELEM;
    }
    return change;
}

// Temperature changes the fluid weight when heat is transported.
#pragma acc routine seq
inline bool hydraulic_buoyancy(const Param &param)
{
    return param.control.has_thermal_diffusion || param.control.has_fluid_thermal_advection;
}

// Darcy drive with Boussinesq density rho_f0*(1 - alpha_f*(T - T0)). The
// Dirichlet conduction profile is a nodal potential; element means use only
// its temperature departure, keeping conduction at rest. T0 is surface_temperature.
#pragma acc routine seq
inline void hydraulic_drive(const Param &param, const Variables &var, const double_vec &pressure,
                            int e, const double *const grad[NDIMS], double drive[NDIMS])
{
    const double weight = var.mat->rho_fluid_reference(e) * param.control.gravity;
    const bool buoyant = hydraulic_buoyancy(param) && weight != 0 && var.mat->alpha_fluid(e) != 0;
    if (!buoyant) {
        // Without thermal buoyancy: grad p + rho_f*g*e_z.
        for (int d = 0; d < NDIMS; ++d) {
            drive[d] = 0;
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
                const int n = (*var.connectivity)[e][j];
                drive[d] += grad[d][j] * (pressure[n] + weight*(*var.coord)[n][NDIMS-1]);
            }
        }
        return;
    }
    const double buoy = var.mat->alpha_fluid(e) * weight;
    double slope = 0, excess = 0;
    if (param.bc.tbc[iboundz0] == 1 && param.bc.tbc[iboundz1] == 1) {
        const double bottom = std::isfinite(param.bc.tbc_value[iboundz0]) ?
            param.bc.tbc_value[iboundz0] : param.bc.mantle_temperature;
        const double top = std::isfinite(param.bc.tbc_value[iboundz1]) ?
            param.bc.tbc_value[iboundz1] : param.bc.surface_temperature;
        slope = (bottom - top) / param.mesh.zlength;
    }
    for (int d = 0; d < NDIMS; ++d) drive[d] = 0;
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        const int n = (*var.connectivity)[e][j];
        const double z = (*var.coord)[n][NDIMS-1];
        const double potential = pressure[n] + weight*z + 0.5*buoy*slope*z*z;
        for (int d = 0; d < NDIMS; ++d) drive[d] += grad[d][j] * potential;
        excess += (*var.temperature)[n] - (param.bc.surface_temperature - slope*z);
    }
    drive[NDIMS-1] -= buoy * excess / NODES_PER_ELEM;
}

// Fixed-stress pressure storage, in 1/Pa. Geometry supplies volume/shape weight.
#pragma acc routine seq
inline double hydraulic_storage(const MatProps &mat, int e)
{
    const double alpha = mat.alpha_biot(e);
    return hydraulic_inverse_biot_modulus(mat, e) + alpha * alpha * hydraulic_mean_stress_compliance(mat, e);
}

// Coupled explicit steps use the current strain increment before evaluating
// stress, including when the mechanical mass is scaled. PT and the reduced
// quasistatic flow model retain the fixed-stress split and pending increment.
#pragma acc routine seq
inline bool hydraulic_explicit_coupling(const Param &param)
{
    return param.control.has_hydraulic_diffusion && !param.control.has_PT &&
           (!param.control.is_quasi_static ||
            param.control.has_poroelastic_pressure_feedback);
}

#pragma acc routine seq
inline bool hydraulic_strain_form(const Param &param)
{
    return param.control.hydraulic_content_model == 0 && hydraulic_explicit_coupling(param);
}

#pragma acc routine seq
inline double hydraulic_storage(const Param &param, const MatProps &mat, int e)
{
    return hydraulic_explicit_coupling(param) ? hydraulic_inverse_biot_modulus(mat, e) :
                                       hydraulic_storage(mat, e);
}

// Minimum nodal tangent bounds the local undrained modulus. Density cancels
// from the mass equation; the wave storage has units of inverse pressure.
#pragma acc routine seq
inline double hydraulic_wave_storage(const Variables &var, int e)
{
    if (!var.fluid_split || !(var.fluid_reference_temperature > 0))
        return hydraulic_inverse_biot_modulus(*var.mat,e);
    const auto conn = (*var.connectivity)[e];
    double p = 0, t = 0;
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        p += (*var.ppressure)[conn[j]]/NODES_PER_ELEM;
        t += (*var.temperature)[conn[j]]/NODES_PER_ELEM;
    }
    double storage = hydraulic_inverse_biot_modulus(*var.mat,e);
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        const int n = conn[j];
        const auto state = fluid_state(*var.mat,e,(*var.ppressure)[n],(*var.temperature)[n],
            var.fluid_reference_temperature,var.mat->phi(e),p,t);
        storage = std::min(storage,state.pressure_slope/state.density);
    }
    return storage;
}

#pragma acc routine seq
inline double hydraulic_diffusivity(const Param &param, const MatProps &mat, int e)
{
    return mat.perm_bound * mat.perm(e) / mat.mu_fluid(e) / hydraulic_storage(param, mat, e);
}

// Bulk density with the pores filled to saturation sat. Gravity passes the
// element saturation; mass and wave-speed estimates keep the saturated default.
#pragma acc routine seq
inline double hydraulic_bulk_density(const MatProps &mat, int e, double sat = 1)
{
    const double phi = mat.phi(e);
    return mat.rho(e) * (1 - phi) + mat.rho_fluid_reference(e) * phi * sat;
}

// Owned nodal solid mass and actual fluid content determine physical weight.
// volume_n is the sum of patch volumes (without simplex lumping).
#pragma acc routine seq
inline double hydraulic_bulk_density(const Param &param, const Variables &var, int e, double sat = 1)
{
    if (!var.fluid_split) return hydraulic_bulk_density(*var.mat,e,sat);
    double p = 0, t = 0, strain = 0, density = 0;
    for (int d = 0; d < NDIMS; ++d) strain += (*var.strain)[e][d];
    const double ratio = param.control.has_moving_mesh ? 1 : std::exp(strain);
    const auto conn = (*var.connectivity)[e];
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        p += (*var.ppressure)[conn[j]]/NODES_PER_ELEM;
        t += (*var.temperature)[conn[j]]/NODES_PER_ELEM;
    }
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        const int n = conn[j];
        const auto state = fluid_state(*var.mat,e,(*var.ppressure)[n],(*var.temperature)[n],
                                      var.fluid_reference_temperature,var.mat->phi(e),p,t);
        density += (*var.solid_mass)[n]/(*var.volume_n)[n]+
                   ratio*state.density*state.porosity/NODES_PER_ELEM;
    }
    return density;
}

// Marker-weighted relative permeability of an element's materials at pressure p.
#pragma acc routine seq
inline double element_kr(const MatProps &mat, const int_vec &markers, double p)
{
    if (mat.nmat == 1) return retention_kr(mat.retention(0), p);
    double kr = 0;
    int total = 0;
    for (int m = 0; m < mat.nmat; ++m) {
        if (markers[m] == 0) continue;
        kr += markers[m] * retention_kr(mat.retention(m), p);
        total += markers[m];
    }
    return kr / total;
}

// Marker-weighted effective saturation Se of an element's materials at
// pressure p; sat, if given, receives their saturation S = Sr + (1-Sr)*Se.
// Bishop's chi = Se, so the element's Bishop pressure is Se*p. Both are
// exactly 1 where every present law is saturated.
#pragma acc routine seq
inline double element_se(const MatProps &mat, const int_vec &markers, double p, double *sat = nullptr)
{
    double se = 0, s = 0;
    int total = 0;
    for (int m = 0; m < mat.nmat; ++m) {
        if (markers[m] == 0) continue;
        const Retention &r = mat.retention(m);
        const double se_m = retention_se(r, p);
        se += markers[m] * se_m;
        s += markers[m] * (r.sr + (1 - r.sr) * se_m);
        total += markers[m];
    }
    if (sat) *sat = s / total;
    return se / total;
}

struct ElementPressure {
    double raw, bishop, strength;
};

// The hydraulic residual and constitutive increment use opposite signs.
// Keep each subtraction order, including the saturated raw-pressure path.
#pragma acc routine seq
template <bool residual>
inline ElementPressure element_pressure(const Variables &var, int e,
                                        const double_vec &p, const double_vec &dp,
                                        bool increment, bool strength)
{
    ConstConnAccessor conn = (*var.connectivity)[e];
    const int_vec &markers = (*var.elemmarkers)[e];
    const bool retention = increment && var.mat->has_retention;
    ElementPressure value = {0, 0, 0};
    bool saturated = true;
    #pragma acc loop seq
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        const int n = conn[j];
        if (increment) {
            if constexpr (residual) value.raw -= dp[n] / NODES_PER_ELEM;
            else value.raw += dp[n] / double(NODES_PER_ELEM);
        }
        if (!retention && !strength) continue;
        const double se = element_se(*var.mat, markers, p[n]);
        if (strength) value.strength += se * p[n] / NODES_PER_ELEM;
        if (!retention) continue;
        const double old = p[n] + dp[n];
        const double se_old = element_se(*var.mat, markers, old);
        if (se != 1 || se_old != 1) saturated = false;
        if constexpr (residual) value.bishop += (se * p[n] - se_old * old) / NODES_PER_ELEM;
        else value.bishop += (se_old * old - se * p[n]) / double(NODES_PER_ELEM);
    }
    if (saturated) value.bishop = value.raw;
    return value;
}

// Pressure fields remain allocated for advection/output when transport is off.
void allocate_hydraulic_state(Variables &var);
void validate_hydraulic_boundaries(const Param &param, const Variables &var);
// After nodal remapping, reset derived assembly data without replacing the
// transferred pressure or pending increment.
void reset_hydraulic_assembly(const Param &param, Variables &var);
// Seed the transport timestep before the first step and after material/remesh
// changes. The pressure update also refreshes this bound from current properties.
void refresh_hydraulic_diffusivity(const Param &param, const Variables &var);
// Physical Darcy volume flux (m/s); tensor, retention and gravity conventions
// match pressure assembly. Reuse caller-owned gradients when supplied.
#pragma acc routine seq
void hydraulic_element_flux(const Param &param, const Variables &var, int e, double *q,
                            const double *const gradient[NDIMS] = nullptr);
// Nodal well extraction at schedule ratio 1, var.fluid_source: each screen's
// rate integrated exactly against the linear shape functions. Mesh-derived,
// so it is rebuilt with the mass and needs no checkpoint.
void assemble_well_source(const Param &param, const Variables &var);
// Advance only accepted physical strain and pressure increments; etmp is scratch.
void update_porosity(const Variables &var);
void prepare_pore_recovery(const Variables &var);
double nodal_fluid_mass(const Param &param, const Variables &var, int n);
void initialize_fluid_content(const Param &param, Variables &var);
double_vec hydraulic_reference(const Param &param);

// Reused physical-step workspace. The generating dt is captured before mesh
// motion can select the following step's dt; PT trials never advance it.
class HydraulicStep {
public:
    void begin(const Param &param, const Variables &var);
    void advance(const Param &param, Variables &var);
    void correct_mechanics(const Param &param, Variables &var) const;
private:
    double interval = 0, start_time = 0;
    elem_cache reference, transport, heat_transport;
    double_vec pressure, temperature, initial_mass, initial_split;
    double_vec mass, target, lag, rate, row, outgoing, boundary_rate, source_rate;
    double_vec candidate, candidate_mass;
    double_vec energy, initial_energy, energy_rate, heat_row, heat_source, source_heat, boundary_heat;
    double_vec candidate_temperature, candidate_energy, pending_energy;
    array_t initial_velocity;
    double_vec initial_inertia;
    double content(const Param &param, const Variables &var, int n, double p, double t,
                   double fraction, double &slope, double *thermal_slope = nullptr) const;
    double invert(const Param &param, const Variables &var, int n, double fraction,
                  double value, double guess, double t) const;
    bool equilibrate(const Param &param, const Variables &var, int n, double fraction,
                     double mass_target, double energy_target, double &p, double &t) const;
};

// One exact-overlap transfer for fluid inventories and pending shear energy.
// The linear model retains its storage-weighted correction; finite-density
// mass and sensible energy use global constraints and explicit boundary flows.
class ContentRemap {
public:
    ContentRemap(const Param &param, const Variables &var);
    void transfer(const Variables &var, const array_t &old_coordinate,
                  const conn_t &old_connectivity);
    void restore(const Param &param, Variables &var) const;
private:
    bool linear, nonlinear, active;
    elem_cache old_content;
    shapefn outside;
    double total[3] = {};
    double magnitude[3] = {};
    double mapped[3] = {};
    Array2D<double,6*NODES_PER_ELEM> inventory;
    Array2D<double,6> projected;
    double physical_total[6] = {}, physical_mapped[6] = {};
    double_vec heat_density, mapped_heat;
    double heat_total = 0;
    void restore_inventory(const Param &param, Variables &var) const;
};

// Hydraulic state changes belong to initialization or accepted physical steps;
// PT constitutive trials must not advance transport or its stress reference.
void initial_hydrostatic_state(const Param &param, const Variables &var,
                              double_vec &ppressure, double_vec &dppressure);
void update_old_mean_stress(const Param &param, const Variables &var,
                           tensor_t &stress, double_vec &old_mean_stress);

// dppressure is the negative pressure change across one physical interval.
// Physical-inertia steps consume it in the current stress update; quasistatic/
// PT steps consume it in the next mechanical update. tdot and tmp_result are
// caller-owned assembly scratch and retain no hydraulic history.
void update_pore_pressure(const Param &param, const Variables &var,
                         double_vec &ppressure, double_vec &dppressure,
                         double_vec &tdot, elem_cache &tmp_result,
                         tensor_t &stress, double_vec &old_mean_stress);

#endif
