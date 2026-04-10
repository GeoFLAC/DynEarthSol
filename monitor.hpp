#ifndef DYNEARTHSOL3D_MONITOR_HPP
#define DYNEARTHSOL3D_MONITOR_HPP

#include "parameters.hpp"

void monitor_initialize(const Param& param, Variables& var);
void monitor_write_if_due(const Param& param, Variables& var);
void monitor_before_remesh(const Param& param, Variables& var);
void monitor_remesh_update(const Param& param, Variables& var);
void monitor_finalize(Variables& var);

#endif
