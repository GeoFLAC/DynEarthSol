#ifndef DYNEARTHSOL3D_PSEUDO_TRANSIENT_HPP
#define DYNEARTHSOL3D_PSEUDO_TRANSIENT_HPP

#include "parameters.hpp"

struct VelocityConstraints;

// Force must already be projected into the admissible velocity subspace.
double pt_residual_rms(const Variables& var, const VelocityConstraints& constraints,
                       const array_t& force);

enum class PTStatus { converged, max_iterations, nonfinite };

struct PTResult {
    PTStatus status = PTStatus::max_iterations;
    int iterations = 0;
    double residual = 0, initial_residual = 0;
    double threshold = 0;
};

// Called before any physical constitutive update. One accepted candidate is
// left in Variables; failed solves restore physical-start mechanical history.
PTResult run_physical_step_pt(const Param& param, Variables& var);

// Instantaneous skeleton correction on fixed geometry, with homogeneous supports and
// no physical aging, transport or pressure increment consumption.
PTResult run_initial_equilibrium_pt(const Param& param, Variables& var);

void require_pt_convergence(const Param& param, const Variables& var, const PTResult& result);

#endif
