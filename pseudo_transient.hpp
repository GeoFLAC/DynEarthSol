#ifndef DYNEARTHSOL_PSEUDO_TRANSIENT_HPP
#define DYNEARTHSOL_PSEUDO_TRANSIENT_HPP

#include "parameters.hpp"

void run_pseudo_transient_loop(Param& param, Variables& var);
void run_initial_mechanical_equilibrium(Param& param, Variables& var);

#endif
