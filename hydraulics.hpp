#ifndef DYNEARTHSOL3D_HYDRAULICS_HPP
#define DYNEARTHSOL3D_HYDRAULICS_HPP

#include "parameters.hpp"
#include "matprops.hpp"
#include "fields.hpp"

// phi*fluid + (alpha - phi)*solid: a fluid and a grain property weighted by
// the pore fraction and the rest of the Biot coefficient.
#pragma acc routine seq
inline double pore_weighted(double phi, double alpha, double fluid, double solid)
{
    return phi*fluid+(alpha-phi)*solid;
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
    return pore_weighted(mat.phi(e), mat.alpha_biot(e), mat.beta_fluid(e), mat.beta_mineral(e));
}

// Temperature baseline belongs to the physical step, not to constitutive trials.
#pragma acc routine seq
inline double hydraulic_temperature_change(const Variables &var, int e)
{
    if (!var.old_temperature) return 0;
    double change = 0;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const int n = (*var.connectivity)[e][i];
        change += ((*var.temperature)[n]-(*var.old_temperature)[n])/NODES_PER_ELEM;
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
                            int e, const double grad[NDIMS][NODES_PER_ELEM], double drive[NDIMS])
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
    if (param.bc.tbc[iboundz0] == 1 && param.bc.tbc[iboundz1] == 1)
        slope = (thermal_boundary_value(param, iboundz0) - thermal_boundary_value(param, iboundz1)) /
                param.mesh.zlength;
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
// quasistatic flow model use the fixed-stress split and pending increment.
#pragma acc routine seq
inline bool hydraulic_strain_form(const Param &param)
{
    return param.control.has_hydraulic_diffusion && !param.control.has_PT &&
           (!param.control.is_quasi_static ||
            param.control.has_poroelastic_pressure_feedback);
}

#pragma acc routine seq
inline double hydraulic_storage(const Param &param, const MatProps &mat, int e)
{
    return hydraulic_strain_form(param) ? hydraulic_inverse_biot_modulus(mat, e) :
                                       hydraulic_storage(mat, e);
}

#pragma acc routine seq
inline double hydraulic_diffusivity(const Param &param, const MatProps &mat, int e)
{
    return mat.perm(e) / mat.mu_fluid(e) / hydraulic_storage(param, mat, e);
}

// Bulk density with the pores filled to saturation sat. Gravity passes the
// element saturation; mass and wave-speed estimates keep the saturated default.
#pragma acc routine seq
inline double hydraulic_bulk_density(const MatProps &mat, int e, double sat = 1)
{
    const double phi = mat.phi(e);
    return mat.rho(e) * (1 - phi) + mat.rho_fluid_reference(e) * phi * sat;
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
        s += markers[m] * retention_sat_of_se(r, se_m);
        total += markers[m];
    }
    if (sat) *sat = s / total;
    return se / total;
}

// Element means of the nodal saturation S and effective saturation Se
// (Bishop's chi), and of the relative permeability kr: one averaging rule for
// the body force, the Darcy flux, the pressure step and the output.
#pragma acc routine seq
inline double element_mean_saturation(const MatProps &mat, const int_vec &markers,
                                      ConstConnAccessor conn, const double_vec &p,
                                      double *se_mean = nullptr)
{
    double s = 0, se = 0;
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        double sat;
        se += element_se(mat, markers, p[conn[j]], &sat);
        s += sat;
    }
    if (se_mean) *se_mean = se / NODES_PER_ELEM;
    return s / NODES_PER_ELEM;
}

#pragma acc routine seq
inline double element_mean_kr(const MatProps &mat, const int_vec &markers,
                              ConstConnAccessor conn, const double_vec &p)
{
    double kr = 0;
    for (int j = 0; j < NODES_PER_ELEM; ++j) kr += element_kr(mat, markers, p[conn[j]]);
    return kr / NODES_PER_ELEM;
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
// Physical Darcy volume flux (m/s); retention and gravity conventions
// match pressure assembly. Reuse caller-owned gradients when supplied.
#pragma acc routine seq
void hydraulic_element_flux(const Param &param, const Variables &var, int e, double *q,
                            const double (*gradient)[NODES_PER_ELEM] = nullptr);
// Fix the well screens in the material of the whole mesh (var.well_pieces): the
// configured vertical screens, or, after a remesh, the previous pieces from the
// positions of their ends (2*NDIMS per piece, see well_screen_positions).
void place_well_screens(const Param &param, Variables &var, const double_vec *ends = nullptr);
void well_screen_positions(const Variables &var, double_vec &ends);
// Nodal well extraction at schedule ratio 1, var.fluid_source: each screen's
// rate integrated exactly against the linear shape functions along its pieces.
// The pieces move with the mesh, so it changes only when they are placed.
void assemble_well_source(const Param &param, const Variables &var);

// One exact-overlap transfer for fluid inventories and pending shear energy.
// The linear model applies a storage-weighted correction.
class ContentRemap {
public:
    ContentRemap(const Param &param, const Variables &var);
    void transfer(const Variables &var, const array_t &old_coordinate,
                  const conn_t &old_connectivity);
    void restore(const Param &param, Variables &var) const;
private:
    bool linear, active;
    elem_cache old_content;
    shapefn outside;
    double total[2] = {};
    double magnitude[2] = {};
    double mapped[2] = {};
    double_vec heat_density, mapped_heat;
    double heat_total = 0;
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
