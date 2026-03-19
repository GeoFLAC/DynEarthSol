#include <algorithm>
#include <cmath>
#include "pseudo_transient.hpp"

#include "bc.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "rheology.hpp"

namespace {

int pt_loop_rheology_type(int rheol_type)
{
    // PT is a mechanical re-equilibration loop, not a physical-time
    // constitutive update. Keep only the recoverable rheology during the inner
    // iterations so plastic softening / RSF evolution does not accumulate just
    // because more PT iterations were needed.
    switch (rheol_type) {
    case MatProps::rh_ep:
    case MatProps::rh_ep_rsf:
        return MatProps::rh_elastic;
    case MatProps::rh_evp:
    case MatProps::rh_evp_rsf:
        return MatProps::rh_maxwell;
    default:
        return rheol_type;
    }
}

void advance_relaxation_step(const Param& param, Variables& var)
{
    update_strain_rate(var, *var.strain_rate);
    compute_dvoldt(var, *var.ntmp, *var.etmp);
    compute_edvoldt(var, *var.ntmp, *var.edvoldt);
    update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
        *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
        *var.strain_rate,
        *var.ppressure, *var.dppressure, *var.vel,
        *var.dyn_fric_coeff, *var.state_variable);
    if (param.control.is_using_mixed_stress)
        NMD_stress(param, var, *var.ntmp, *var.stress, *var.etmp);
    update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
    update_velocity(var, *var.vel);
}

void run_relaxation_iterations(Param& param, Variables& var, bool zero_velocity_bc)
{
    const int restore_rheology_type = param.mat.rheol_type;
    param.mat.rheol_type = pt_loop_rheology_type(restore_rheology_type);

    const bool restore_hydraulic_diffusion = param.control.has_hydraulic_diffusion;
    if (restore_hydraulic_diffusion)
        param.control.has_hydraulic_diffusion = false;

    var.l2_residual = calculate_residual_force(var, *var.force_residual);
    double residual_old = var.l2_residual;

    for (int pt_step = 0; pt_step < param.control.PT_max_iter; ++pt_step) {
        apply_vbcs(param, var, *var.vel);

        advance_relaxation_step(param, var);

        var.l2_residual = calculate_residual_force(var, *var.force_residual);
        const double relative_change =
            std::fabs((var.l2_residual - residual_old) / std::max(std::fabs(residual_old), 1e-30));
        if (relative_change < param.control.PT_relative_tolerance) {
            break;
        }
        residual_old = var.l2_residual;
    }

    param.mat.rheol_type = restore_rheology_type;

    if (restore_hydraulic_diffusion)
        param.control.has_hydraulic_diffusion = true;
}

void run_time_step_pseudo_transient_iterations(Param& param, Variables& var)
{
    // Keep the relaxed PT velocity as the baseline state. 
    run_relaxation_iterations(param, var, false);
    apply_vbcs(param, var, *var.vel);
}

} // namespace

void run_pseudo_transient_loop(Param& param, Variables& var)
{
    if (!param.control.has_PT)
        return;

    run_time_step_pseudo_transient_iterations(param, var);
}

void run_initial_mechanical_equilibrium(Param& param, Variables& var)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // This stage only re-equilibrates the initial mechanical imbalance. 
    run_relaxation_iterations(param, var, true);
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}
