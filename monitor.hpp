#ifndef DYNEARTHSOL3D_MONITOR_HPP
#define DYNEARTHSOL3D_MONITOR_HPP

#include "parameters.hpp"

// The tracers' exact reference in the checkpoint mesh; empty when there is none.
struct MonitorCheckpoint {
    double_vec positions;   // NDIMS per point
    int_vec elements;       // containing element, -1 outside the mesh
    double_vec weights;     // NODES_PER_ELEM per point
};

// A restart passes the checkpoint's state; an empty one locates the configured points.
void monitor_initialize(const Param& param, Variables& var, const MonitorCheckpoint& saved);
void monitor_write_if_due(const Param& param, Variables& var);
void monitor_before_remesh(const Param& param, Variables& var);
void monitor_remesh_update(const Param& param, Variables& var);
void monitor_finalize(Variables& var);
MonitorCheckpoint monitor_checkpoint(const Param& param, const Variables& var);
#endif
