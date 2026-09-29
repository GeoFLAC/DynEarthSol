#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>

#include "pseudo_transient.hpp"
#include "bc.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "rheology.hpp"
#include "utils.hpp"

namespace {
struct PTWorkspace {
    MechanicalState start;
    VelocityConstraints constraints;
    tensor_t lagged_stress;
    array_t old_velocity, undamped;
    double_vec stress_fraction, mobility;
    PTWorkspace(int ne, int nn) : start(ne), constraints(nn), lagged_stress(ne),
        old_velocity(nn), undamped(nn), stress_fraction(ne), mobility(nn) {}
};

double residual_rms(const Variables& var, const VelocityConstraints& constraints,
                    const array_t& force)
{
    double sum = 0;
    int rank = 0, bad = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, constraints, force) reduction(+:sum,rank) reduction(|:bad)
#endif
    #pragma acc parallel loop gang vector reduction(+:sum,rank) reduction(|:bad)
    for (int n=0; n<var.nnode; ++n) {
        rank += constraints.free_rank[n];
        for (int d=0; d<NDIMS; ++d) {
            const double f = force[n][d];
            bad |= !std::isfinite(f) || !std::isfinite((*var.vel)[n][d]);
            sum += f*f;
        }
    }
    if (bad || !std::isfinite(sum)) return std::numeric_limits<double>::quiet_NaN();
    return rank ? std::sqrt(sum/rank) : 0.0;
}

bool finite_candidate(const Variables& var)
{
    int bad = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var) reduction(|:bad)
#endif
    #pragma acc parallel loop gang vector reduction(|:bad)
    for (int e=0; e<var.nelem; ++e) {
        for (int d=0; d<NSTR; ++d)
            bad |= !std::isfinite((*var.stress)[e][d]) || !std::isfinite((*var.strain)[e][d]);
        bad |= !std::isfinite((*var.stressyy)[e]) || !std::isfinite((*var.plstrain)[e]) ||
               !std::isfinite((*var.delta_plstrain)[e]) || !std::isfinite((*var.viscosity)[e]) ||
               !std::isfinite((*var.dpressure)[e]) ||
               !std::isfinite((*var.state_variable)[e]) || !std::isfinite((*var.dyn_fric_coeff)[e]);
    }
    return !bad;
}

PTResult solve(const Param& param, Variables& var, bool initial)
{
    #pragma acc wait
    // Heap ownership also keeps the workspace's nested arrays available to the
    // managed-memory accelerator path. No snapshot survives mesh motion/remesh.
    std::unique_ptr<PTWorkspace> owner(new PTWorkspace(var.nelem, var.nnode));
    PTWorkspace& w = *owner;
    w.start.capture(var);
    for (int n=0; n<var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) w.old_velocity[n][d] = (*var.vel)[n][d];
    for (int e=0; e<var.nelem; ++e)
        for (int d=0; d<NSTR; ++d) w.lagged_stress[e][d] = (*var.stress)[e][d];
    build_velocity_constraints(param, var, w.constraints, initial);
    if (initial)
        for (int n=0; n<var.nnode; ++n)
            for (int d=0; d<NDIMS; ++d) (*var.vel)[n][d] = 0;

    PTResult result;
    double best = std::numeric_limits<double>::infinity();
    int last_improvement = 0;
    for (int iteration=0; iteration<param.control.PT_max_iter; ++iteration) {
        project_free_vectors(var, w.constraints, *var.vel, true);
        evaluate_mechanical_trial(param, var, w.start);
        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result, &w.undamped);
        project_free_vectors(var, w.constraints, w.undamped);
        result.iterations = iteration+1;
        result.residual = residual_rms(var, w.constraints, w.undamped);
        var.l2_residual = result.residual;
        if (iteration == 0) result.initial_residual = result.residual;
        if (!std::isfinite(result.residual) || !finite_candidate(var)) {
            result.status = PTStatus::nonfinite;
            break;
        }
        const double threshold = param.control.PT_absolute_tolerance +
            param.control.PT_relative_tolerance*result.initial_residual;
        if (result.residual <= threshold) {
            result.status = PTStatus::converged;
            // The accepted stress is the constitutive/NMD candidate, never the
            // solver's lagged stress. No second constitutive update follows.
            for (int n=0; n<var.nnode; ++n)
                for (int d=0; d<NDIMS; ++d) (*var.force_residual)[n][d] = w.undamped[n][d];
            if (initial)
                for (int n=0; n<var.nnode; ++n)
                    for (int d=0; d<NDIMS; ++d) (*var.vel)[n][d] = w.old_velocity[n][d];
            return result;
        }
        if (result.residual < best*(1-1e-12)) {
            best = result.residual;
            last_improvement = iteration;
        } else if (param.control.PT_stagnation_window > 0 &&
                   iteration-last_improvement >= param.control.PT_stagnation_window) {
            result.status = PTStatus::stagnated;
            break;
        }
        if (iteration+1 == param.control.PT_max_iter) break;
        compute_pt_factors(param, var, w.stress_fraction, w.mobility);
#ifndef ACC
        #pragma omp parallel for default(none) shared(var, w)
#endif
        #pragma acc parallel loop gang vector async
        for (int e=0; e<var.nelem; ++e)
            for (int d=0; d<NSTR; ++d) {
                w.lagged_stress[e][d] += w.stress_fraction[e]*((*var.stress)[e][d]-w.lagged_stress[e][d]);
                (*var.stress)[e][d] = w.lagged_stress[e][d];
            }
        // Existing assembly supplies exactly the same physical loads to the
        // lagged solver stress. Artificial damping is excluded from the update.
        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result, &w.undamped);
        project_free_vectors(var, w.constraints, w.undamped);
#ifndef ACC
        #pragma omp parallel for default(none) shared(var, w)
#endif
        #pragma acc parallel loop gang vector async
        for (int n=0; n<var.nnode; ++n)
            for (int d=0; d<NDIMS; ++d)
                (*var.vel)[n][d] += w.mobility[n]*w.undamped[n][d];
        #pragma acc wait
    }
    // Failed candidates cannot become physical history. The caller stops before
    // transport, mesh motion, output or a new physical time step can consume them.
    w.start.restore(var);
    for (int n=0; n<var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) (*var.vel)[n][d] = w.old_velocity[n][d];
    update_strain_rate(var, *var.strain_rate);
    compute_dvoldt(var, *var.ntmp, *var.etmp);
    compute_edvoldt(var, *var.ntmp, *var.edvoldt);
    update_force(param, var, *var.force, *var.force_residual, *var.tmp_result);
    return result;
}
}

void evaluate_mechanical_trial(const Param& param, Variables& var,
                               const MechanicalState& physical_start)
{
    physical_start.restore(var);
    update_strain_rate(var, *var.strain_rate);
    compute_dvoldt(var, *var.ntmp, *var.etmp);
    compute_edvoldt(var, *var.ntmp, *var.edvoldt);
    update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
                  *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
                  *var.strain_rate, *var.ppressure, *var.dppressure, *var.vel,
                  *var.dyn_fric_coeff, *var.state_variable, true);
    if (param.control.is_using_mixed_stress)
        NMD_stress(var, *var.stress, *var.ntmp, *var.etmp);
    #pragma acc wait
}

PTResult run_physical_step_pt(const Param& param, Variables& var)
{
    return solve(param, var, false);
}

PTResult run_initial_equilibrium_pt(const Param& param, Variables& var)
{
    // Initial prestress correction is elastic, not a physical creep/aging step.
    // Do not silently use that contract for an inelastic initial state.
    if (param.mat.rheol_type != MatProps::rh_elastic)
        die(EXIT_CONFIG_VALUE, "Initial PT equilibrium currently requires elastic rheology; duration-based isostasy is separate.");
    Param initial = param;
    initial.control.has_hydraulic_diffusion = false;
    return solve(initial, var, true);
}

const char* pt_status_name(PTStatus status)
{
    switch (status) {
    case PTStatus::converged: return "converged";
    case PTStatus::max_iterations: return "max_iterations";
    case PTStatus::stagnated: return "stagnated";
    case PTStatus::nonfinite: return "nonfinite";
    }
    return "unknown";
}

void require_pt_convergence(const Variables& var, const PTResult& result)
{
    std::printf("PT step=%d status=%s iterations=%d residual=%.17g initial=%.17g\n",
                var.steps, pt_status_name(result.status), result.iterations, result.residual, result.initial_residual);
    if (result.status != PTStatus::converged)
        die(result.status == PTStatus::nonfinite ? EXIT_RUNTIME_NAN : EXIT_RUNTIME_NONCONVERGENCE,
            "PT mechanical equilibrium failed; physical-start mechanical history restored.");
}
