#ifndef PSEUDO_TRANSIENT_HPP
#define PSEUDO_TRANSIENT_HPP

#include "parameters.hpp"

struct MechanicalState;
// Reusable trial evaluation for future outer pressure/mechanics iterations.
// Velocity and pressure inputs live in var; physical_start is never modified.
// Geometry and physical dt must remain fixed for the lifetime of the baseline.
void evaluate_mechanical_trial(const Param& param, Variables& var,
                               const MechanicalState& physical_start, bool initial_equilibrium = false);

struct VelocityConstraints;
// Force must already be projected into the admissible velocity subspace.
double pt_residual_rms(const Variables& var, const VelocityConstraints& constraints,
                       const array_t& force);

enum class PTStatus { converged, max_iterations, stagnated, nonfinite };
struct PTResult {
    PTStatus status = PTStatus::max_iterations;
    int iterations = 0;
    double residual = 0, initial_residual = 0;
};
const char* pt_status_name(PTStatus status);
// Called before any physical constitutive update. One accepted candidate is
// left in Variables; failed solves restore physical-start mechanical history.
PTResult run_physical_step_pt(const Param& param, Variables& var);
// Instantaneous skeleton correction on fixed geometry, with homogeneous supports and
// no physical aging, transport or pressure increment consumption.
PTResult run_initial_equilibrium_pt(const Param& param, Variables& var);
void require_pt_convergence(const Variables& var, const PTResult& result);

#endif
