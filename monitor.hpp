#ifndef DYNEARTHSOL3D_MONITOR_HPP
#define DYNEARTHSOL3D_MONITOR_HPP

#include "parameters.hpp"

void monitor_initialize(const Param& param, Variables& var);
void monitor_write_if_due(const Param& param, Variables& var);
void monitor_before_remesh(const Param& param, Variables& var);
void monitor_remesh_update(const Param& param, Variables& var);
void monitor_finalize(Variables& var);
// Save and restore the tracer's exact reference in the checkpoint mesh.
void monitor_checkpoint(const Param& param, const Variables& var,
                        double_vec& positions, int_vec& elements, double_vec& weights);
void monitor_restore(const Variables& var, const double_vec& positions,
                     const int_vec& elements, const double_vec& weights);

#endif
