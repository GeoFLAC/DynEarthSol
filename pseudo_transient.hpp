#ifndef PSEUDO_TRANSIENT_HPP
#define PSEUDO_TRANSIENT_HPP

#include "parameters.hpp"

// Existing definition lives in dynearthsol.cxx. Declared here so the
// pseudo-transient orchestration below (a separate translation unit) can reuse it.
void update_mesh(const Param &param, Variables &var);

// Pseudo-transient relaxation of the current physical step, extracted verbatim
// from main() in dynearthsol.cxx. hydraulic_diffusion_switch is owned by the
// caller and is set by this function when it disables
// param.control.has_hydraulic_diffusion for the PT loop; this function restores
// the flag when the switch is set before returning.
void run_physical_step_pt(Param &param, Variables &var, bool &hydraulic_diffusion_switch);

#endif
