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

struct PTWorkspace {
    MechanicalState start;
    VelocityConstraints constraints;
    tensor_t lagged_stress;
    array_t old_velocity, undamped;
    double_vec stress_fraction, mobility, height, effective_viscosity;
    // Adaptive dynamic relaxation only (PT_option=1); empty otherwise.
    double_vec mass_rows;
    array_t mass, momentum, previous_force;

    PTWorkspace(int nelem, int nnode, bool adaptive)
        : start(nelem), constraints(nnode), lagged_stress(adaptive ? 0 : nelem),
          old_velocity(nnode), undamped(nnode), stress_fraction(adaptive ? 0 : nelem),
          mobility(adaptive ? 0 : nnode), height(adaptive ? 0 : nelem), effective_viscosity(nelem),
          mass_rows(adaptive ? static_cast<std::size_t>(nelem)*NODES_PER_ELEM*NDIMS : 0),
          mass(adaptive ? nnode : 0), momentum(adaptive ? nnode : 0),
          previous_force(adaptive ? nnode : 0) {}
};

// One adaptive dynamic relaxation update (Underwood 1983; Papadrakakis 1981).
// The Gershgorin mass bounds the linear reference operator, not the full
// nonlinear response. Estimate a low mode with a Rayleigh quotient along
// the latest increment, F_prev - F = K dv. Momentum is dropped when it opposes
// the current force (kinetic damping; gradient restart, O'Donoghue & Candes 2015).
// Only the route to equilibrium changes; acceptance uses the physical residual.
void update_adaptive_dr(const Variables& var, PTWorkspace& w, int iteration, double& lowest)
{
    double num = 0, den = 0, drive = 0;
    if (iteration > 0) {
#ifndef ACC
        #pragma omp parallel for default(none) shared(var, w) reduction(+:num,den,drive)
#endif
        #pragma acc parallel loop gang vector reduction(+:num,den,drive)
        for (int n=0; n<var.nnode; ++n)
            for (int d=0; d<NDIMS; ++d) {
                const double dv = w.momentum[n][d];
                num += dv*(w.previous_force[n][d] - w.undamped[n][d]);
                den += w.mass[n][d]*dv*dv;
                drive += dv*w.undamped[n][d];
            }
        if (num > 0 && den > 0)
            lowest = std::max(1e-12, std::min(1.0, num/den));
    }
    const double s = std::sqrt(lowest);
    const double alpha = 4/((1+s)*(1+s));
    const double beta = (iteration > 0 && drive >= 0) ? ((1-s)/(1+s))*((1-s)/(1+s)) : 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, w) firstprivate(alpha, beta)
#endif
    #pragma acc parallel loop gang vector
    for (int n=0; n<var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) {
            w.momentum[n][d] = (beta > 0 ? beta*w.momentum[n][d] : 0.0) +
                alpha*w.undamped[n][d]/w.mass[n][d];
            w.previous_force[n][d] = w.undamped[n][d];
        }
    project_free_vectors(var, w.constraints, w.momentum);
#ifndef ACC
    #pragma omp parallel for default(none) shared(var, w)
#endif
    #pragma acc parallel loop gang vector
    for (int n=0; n<var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) (*var.vel)[n][d] += w.momentum[n][d];
}

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

bool is_finite_candidate(const Variables& var)
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

const char* pt_status_name(PTStatus status)
{
    switch (status) {
    case PTStatus::converged: return "converged";
    case PTStatus::max_iterations: return "max_iterations";
    case PTStatus::nonfinite: return "nonfinite";
    }
    return "unknown";
}

PTResult solve_pt(const Param& param, Variables& var, bool initial)
{
    #pragma acc wait
    // Heap ownership also keeps the workspace's nested arrays available to the
    // managed-memory accelerator path. No snapshot survives mesh motion/remesh.
    const bool adaptive = param.control.PT_option == 1;
    std::unique_ptr<PTWorkspace> owner(new PTWorkspace(var.nelem, var.nnode, adaptive));
    PTWorkspace& w = *owner;
    w.start.capture(var);
    for (int n=0; n<var.nnode; ++n)
        for (int d=0; d<NDIMS; ++d) w.old_velocity[n][d] = (*var.vel)[n][d];
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
    bool stagnation_reported = false;
    double lowest_mode = 1;
    for (int iteration=0; iteration<param.control.PT_max_iter; ++iteration) {
        project_free_vectors(var, w.constraints, *var.vel, true);
        evaluate_mechanical_trial(param, var, w.start, initial);
        update_force(param, var, *var.force, *var.force_residual, *var.tmp_result, &w.undamped, displacement_dt);
        project_free_vectors(var, w.constraints, w.undamped);
        result.iterations = iteration+1;
        result.residual = pt_residual_rms(var, w.constraints, w.undamped);
        var.l2_residual = result.residual;
        if (iteration == 0) {
            // Match the auxiliary solver stress to the first constitutive trial.
            if (!adaptive)
                for (int e=0; e<var.nelem; ++e)
                    for (int d=0; d<NSTR; ++d) w.lagged_stress[e][d] = (*var.stress)[e][d];
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
        if (!std::isfinite(result.residual) || !std::isfinite(roundoff_floor) || !is_finite_candidate(var)) {
            result.status = PTStatus::nonfinite;
            break;
        }
        if (result.residual <= result.threshold) {
            result.status = PTStatus::converged;
            if (running_scale) {
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
        } else if (!stagnation_reported && param.control.PT_stagnation_window > 0 &&
                   iteration-last_improvement >= param.control.PT_stagnation_window) {
            std::fprintf(stderr, "Warning: PT step=%d has not improved its best residual for %d iterations; "
                         "continuing (residual=%.17g, best=%.17g, threshold=%.17g).\n",
                         var.steps, param.control.PT_stagnation_window, result.residual, best, result.threshold);
            stagnation_reported = true;
        }
        if (iteration+1 == param.control.PT_max_iter) break;
        // Geometry, physical dt and elastic moduli are fixed during this solve.
        // Only a viscous rheology adds candidate-dependent viscosity to factors.
        if (iteration == 0 || (param.mat.rheol_type & MatProps::rh_viscous)) {
            compute_pt_factors(param, var, w.stress_fraction, w.mobility,
                               w.height, w.effective_viscosity);
            if (adaptive)
                compute_pt_mass(var, w.effective_viscosity, w.mass_rows, w.mass);
        }
        if (adaptive) {
            update_adaptive_dr(var, w, iteration, lowest_mode);
            continue;
        }
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

} // anonymous namespace


PTResult run_physical_step_pt(const Param& param, Variables& var)
{
    return solve_pt(param, var, false);
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
    return solve_pt(initial, var, true);
}


void require_pt_convergence(const Param& param, const Variables& var, const PTResult& result)
{
    // Keep failures visible even between ordinary progress reports.
    if (result.status != PTStatus::converged || var.steps == 0 ||
        (var.steps >= var.info_display_next_step &&
         var.steps % param.mesh.quality_check_step_interval == 0))
        std::fprintf(result.status == PTStatus::converged ? stdout : stderr,
                     "PT step=%d status=%s iterations=%d residual=%.17g initial=%.17g threshold=%.17g\n",
                     var.steps, pt_status_name(result.status), result.iterations,
                     result.residual, result.initial_residual, result.threshold);
    if (result.status != PTStatus::converged) {
        die(result.status == PTStatus::nonfinite ? EXIT_RUNTIME_NAN : EXIT_RUNTIME_NONCONVERGENCE,
            "PT mechanical equilibrium failed; physical-start mechanical history restored.");
    }
}
