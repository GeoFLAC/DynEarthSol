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
}

double pt_residual_rms(const Variables& var, const VelocityConstraints& constraints,
                    const array_t& force)
{
    double scale = 0;
    long long rank = 0;
    int bad = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, constraints, force) reduction(max:scale) reduction(+:rank) reduction(|:bad)
#endif
    #pragma acc parallel loop gang vector reduction(max:scale) reduction(+:rank) reduction(|:bad)
    // Flatten components so every reduction update belongs to the parallel loop.
    for (long long k=0; k<static_cast<long long>(var.nnode)*NDIMS; ++k) {
        const int n = k/NDIMS, d = k%NDIMS;
        if (d == 0) rank += constraints.free_rank[n];
        const double f = force[n][d];
        bad |= !std::isfinite(f) || !std::isfinite((*var.vel)[n][d]);
        scale = std::max(scale, std::abs(f));
    }
    if (bad) return std::numeric_limits<double>::quiet_NaN();
    if (!rank || scale == 0) return 0;
    // Scale before squaring: large forces must not overflow and tiny imbalances
    // must not underflow to an artificial zero residual.
    double sum = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, force) firstprivate(scale) reduction(+:sum)
#endif
    #pragma acc parallel loop gang vector reduction(+:sum)
    for (int n=0; n<var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) {
            const double f = force[n][d]/scale;
            sum += f*f;
        }
    return scale*std::sqrt(sum/rank);
}

namespace {
// Upper magnitude estimate in the free subspace, from the already assembled
// element cache. Absolute projector coefficients avoid cancellation in oblique
// constraints. This is a roundoff scale, not a physical error bound.
double absolute_force_rms(const Variables& var, const VelocityConstraints& constraints,
                          array_t& magnitude)
{
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, constraints, magnitude)
#endif
    #pragma acc parallel loop gang vector
    for (int n=0; n<var.nnode; ++n) {
        double assembled[NDIMS] = {};
        const int* patch = var.support.patch(n);
        const int* local = var.support.local(n);
        for (int k=0; k<var.support.size(n); ++k) {
            ConstElemCacheAccessor tr = (*var.tmp_result)[patch[k]];
            for (int d=0; d<NDIMS; ++d)
                assembled[d] += std::abs(tr[local[k]+NODES_PER_ELEM*d]);
        }
        for (int i=0; i<NDIMS; ++i) {
            double value = 0;
            for (int j=0; j<NDIMS; ++j)
                value += std::abs(constraints.projector[n][i*NDIMS+j])*assembled[j];
            magnitude[n][i] = value;
        }
    }
    return pt_residual_rms(var, constraints, magnitude);
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

    // Physical transport uses this same dt once after convergence. Initial
    // equilibrium and deliberately fixed-grid runs retain reference loads.
    const double displacement_dt = !initial && param.control.has_moving_mesh ? var.dt : 0;
    const bool running_scale = !initial && param.control.PT_use_running_scale;
    double trial_reference = running_scale ? var.PT_initial_residual_max : 0;
    std::unique_ptr<array_t> magnitude;
    if (running_scale) magnitude.reset(new array_t(var.nnode));
    double roundoff_floor = 0;
    PTResult result;
    double best = std::numeric_limits<double>::infinity();
    int last_improvement = 0;
    for (int iteration=0; iteration<param.control.PT_max_iter; ++iteration) {
        project_free_vectors(var, w.constraints, *var.vel, true);
        evaluate_mechanical_trial(param, var, w.start, initial);
        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result, &w.undamped, displacement_dt);
        project_free_vectors(var, w.constraints, w.undamped);
        result.iterations = iteration+1;
        result.residual = pt_residual_rms(var, w.constraints, w.undamped);
        var.l2_residual = result.residual;
        if (iteration == 0) {
            result.initial_residual = result.residual;
            result.threshold = param.control.PT_absolute_tolerance +
                param.control.PT_relative_tolerance*result.initial_residual;
            if (running_scale) {
                if (!var.PT_skip_scale_update)
                    trial_reference = std::max(trial_reference, result.initial_residual);
                roundoff_floor = 1000*std::numeric_limits<double>::epsilon()*
                    absolute_force_rms(var, w.constraints, *magnitude);
                result.threshold = std::max(1e-6*trial_reference, roundoff_floor);
            }
        }
        if (!std::isfinite(result.residual) || !std::isfinite(roundoff_floor) || !finite_candidate(var)) {
            result.status = PTStatus::nonfinite;
            break;
        }
        if (result.residual <= result.threshold) {
            result.status = PTStatus::converged;
            if (running_scale) {
                std::printf("PT scale step=%d initial=%.17g reference=%.17g floor=%.17g threshold=%.17g skipped=%d\n",
                            var.steps, result.initial_residual, trial_reference, roundoff_floor,
                            result.threshold, static_cast<int>(var.PT_skip_scale_update));
                var.PT_initial_residual_max = trial_reference;
                var.PT_skip_scale_update = false;
            }
            // The accepted stress is the constitutive/NMD candidate, never the
            // solver's lagged stress. No second constitutive update follows.
            for (int n=0; n<var.nnode; ++n)
                for (int d=0; d<NDIMS; ++d) (*var.force_residual)[n][d] = w.undamped[n][d];
            if (initial) {
                for (int n=0; n<var.nnode; ++n)
                    for (int d=0; d<NDIMS; ++d) (*var.vel)[n][d] = w.old_velocity[n][d];
                // The correction velocity is numerical. Initial output and
                // physical transport must see rates from the restored velocity.
                update_strain_rate(var, *var.strain_rate);
                compute_dvoldt(var, *var.ntmp, *var.etmp);
                compute_edvoldt(var, *var.ntmp, *var.edvoldt);
                #pragma acc wait
            }
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
        // Geometry, physical dt and elastic moduli are fixed during this solve.
        // Only a viscous rheology adds candidate-dependent viscosity to factors.
        if (iteration == 0 || (param.mat.rheol_type & MatProps::rh_viscous))
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
        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result, &w.undamped, displacement_dt);
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
                               const MechanicalState& physical_start, bool initial_equilibrium)
{
    physical_start.restore(var);
    update_strain_rate(var, *var.strain_rate);
    compute_dvoldt(var, *var.ntmp, *var.etmp);
    compute_edvoldt(var, *var.ntmp, *var.edvoldt);
    update_stress(param, var, *var.stress, *var.stressyy, *var.dpressure,
                  *var.viscosity, *var.strain, *var.plstrain, *var.delta_plstrain,
                  *var.strain_rate, *var.ppressure, *var.dppressure, *var.vel,
                  *var.dyn_fric_coeff, *var.state_variable, true, initial_equilibrium);
    if (param.control.is_using_mixed_stress)
        NMD_stress(var, *var.stress, *var.ntmp, *var.etmp,
                   param.mat.is_plane_strain ? var.stressyy : nullptr);
    #pragma acc wait
}

PTResult run_physical_step_pt(const Param& param, Variables& var)
{
    return solve(param, var, false);
}

PTResult run_initial_equilibrium_pt(const Param& param, Variables& var)
{
    // Instantaneous limit: no viscous relaxation and no RSF aging. Plastic
    // correction remains part of the skeleton's initial equilibrium state.
    if (!(param.mat.rheol_type & MatProps::rh_elastic))
        die(EXIT_CONFIG_VALUE, "Initial instantaneous PT equilibrium requires an elastic skeleton; use a physical-time process for pure viscous material.");
    std::unique_ptr<Param> initial_owner(new Param(param));
    Param& initial = *initial_owner;
    // Keep the physical load definition: the hydraulic option also controls
    // sidewall support pressure. This solve never advances fluid transport;
    // update_stress uses initial_equilibrium to leave pending pressure alone.
    initial.mat.rheol_type &= ~MatProps::rh_viscous;
    if (initial.mat.rheol_type & MatProps::rh_rsf) {
        update_strain_rate(var, *var.strain_rate);
        #pragma acc wait
        refresh_rsf_friction(initial, var, *var.dyn_fric_coeff, *var.state_variable);
        #pragma acc wait
    }
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
    if (result.status != PTStatus::converged) {
        std::fprintf(stderr, "PT stopping threshold=%.17g; no candidate was accepted.\n",
                     result.threshold);
        die(result.status == PTStatus::nonfinite ? EXIT_RUNTIME_NAN : EXIT_RUNTIME_NONCONVERGENCE,
            "PT mechanical equilibrium failed; physical-start mechanical history restored.");
    }
}
