#ifndef DYNEARTHSOL3D_HYDRAULICS_HPP
#define DYNEARTHSOL3D_HYDRAULICS_HPP

#include "parameters.hpp"
#include "matprops.hpp"

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

// Darcy drive grad p + rho_f*g*e_z.
#pragma acc routine seq
inline void hydraulic_drive(const Param &param, const Variables &var, const double_vec &pressure,
                            int e, const double grad[NDIMS][NODES_PER_ELEM], double drive[NDIMS])
{
    const double weight = var.mat->rho_fluid_reference(e) * param.control.gravity;
    for (int d = 0; d < NDIMS; ++d) {
        drive[d] = 0;
        for (int j = 0; j < NODES_PER_ELEM; ++j) {
            const int n = (*var.connectivity)[e][j];
            drive[d] += grad[d][j] * (pressure[n] + weight*(*var.coord)[n][NDIMS-1]);
        }
    }
}

// Fixed-stress pressure storage, in 1/Pa. Geometry supplies volume/shape weight.
#pragma acc routine seq
inline double hydraulic_storage(const MatProps &mat, int e)
{
    const double alpha = mat.alpha_biot(e);
    return hydraulic_inverse_biot_modulus(mat, e) + alpha * alpha * hydraulic_mean_stress_compliance(mat, e);
}

// Advance pressure before stress in physical-inertia or coupled explicit steps.
// PT and the reduced quasistatic flow model advance pressure after mechanics.
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
    // Without mechanical feedback, use the reduced fixed-stress capacity
    // at either inertia setting.
    return hydraulic_strain_form(param) && param.control.has_poroelastic_pressure_feedback
        ? hydraulic_inverse_biot_modulus(mat, e) : hydraulic_storage(mat, e);
}

#pragma acc routine seq
inline double hydraulic_mobility(const MatProps &mat, int e)
{
    return mat.perm(e) / mat.mu_fluid(e);
}

#pragma acc routine seq
inline double hydraulic_diffusivity(const Param &param, const MatProps &mat, int e)
{
    return hydraulic_mobility(mat, e) / hydraulic_storage(param, mat, e);
}

// Bulk density of the saturated solid and fluid mixture.
#pragma acc routine seq
inline double hydraulic_bulk_density(const MatProps &mat, int e)
{
    const double phi = mat.phi(e);
    return mat.rho(e) * (1 - phi) + mat.rho_fluid_reference(e) * phi;
}

struct ElementPressure {
    double increment, mean;
};

// The hydraulic residual and constitutive increment use opposite signs.
#pragma acc routine seq
template <bool residual>
inline ElementPressure element_pressure(const Variables &var, int e,
                                        const double_vec &p, const double_vec &dp,
                                        bool increment, bool mean)
{
    ConstConnAccessor conn = (*var.connectivity)[e];
    ElementPressure value = {0, 0};
    #pragma acc loop seq
    for (int j = 0; j < NODES_PER_ELEM; ++j) {
        const int n = conn[j];
        if (increment) {
            if constexpr (residual) value.increment -= dp[n] / NODES_PER_ELEM;
            else value.increment += dp[n] / double(NODES_PER_ELEM);
        }
        if (mean) value.mean += p[n] / NODES_PER_ELEM;
    }
    return value;
}

// Pressure fields remain allocated for advection/output when transport is off.
void allocate_hydraulic_state(Variables &var);
void validate_hydraulic_boundaries(const Param &param, const Variables &var);
// After nodal remapping, reset derived assembly data without replacing the
// transferred pressure or pending increment.
void reset_hydraulic_assembly(const Param &param, Variables &var);
// Fix the well screens in the material of the whole mesh (var.well_pieces): the
// configured vertical screens, or, after a remesh, the previous pieces from the
// positions of their ends (2*NDIMS per piece, see well_screen_positions).
void place_well_screens(const Param &param, Variables &var, const double_vec *ends = nullptr);
void well_screen_positions(const Variables &var, double_vec &ends);
// Nodal well supply at schedule ratio 1, var.fluid_source: each screen's
// rate integrated exactly against the linear shape functions along its pieces.
// The pieces move with the mesh, so it changes only when they are placed.
void assemble_well_source(const Param &param, const Variables &var);

// One exact-overlap transfer for the fluid inventories. The linear model
// applies a storage-weighted correction.
class ContentRemap {
public:
    ContentRemap(const Param &param, const Variables &var);
    void transfer(const Variables &var, const array_t &old_coordinate,
                  const conn_t &old_connectivity);
    void restore(const Param &param, Variables &var) const;
private:
    bool linear;
    elem_cache old_content;
    shapefn outside;
    double total[2] = {};
    double magnitude[2] = {};
    double mapped[2] = {};
};

// Hydraulic state changes belong to initialization or accepted physical steps;
// PT constitutive trials must not advance transport or its stress reference.
void initial_hydrostatic_state(const Param &param, const Variables &var,
                              double_vec &ppressure, double_vec &dppressure);
void update_old_mean_stress(const Param &param, const Variables &var,
                           tensor_t &stress, double_vec &old_mean_stress);

// dppressure is the negative pressure change across one physical interval.
// hydraulic_strain_form selects the pressure order. tdot and tmp_result
// are assembly scratch.
void update_pore_pressure(const Param &param, const Variables &var,
                         double_vec &ppressure, double_vec &dppressure,
                         double_vec &tdot, elem_cache &tmp_result,
                         tensor_t &stress, double_vec &old_mean_stress);

#endif
