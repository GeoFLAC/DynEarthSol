#ifndef DYNEARTHSOL3D_POINT_SOURCES_HPP
#define DYNEARTHSOL3D_POINT_SOURCES_HPP

#include "parameters.hpp"

// Visit each current-mesh contribution once. The rate has volume/time units
// for fluid sources and energy/time units for heat sources; the receiver owns
// conversion to physical inventories. No point locator is duplicated there.
using PointSourceSink = void (*)(int point, int node, double rate, void *context);
void deposit_point_sources(const Injection &injection, const Variables &var,
                           double start, double stop, PointSourceSink sink, void *context,
                           const double_vec *scale = nullptr, bool peak_rate = false);

// Mean source rate over an accepted (possibly subcycled) interval. Positive
// rates supply fluid volume or heat; shape weights conserve the requested sum.
// scale, if given, multiplies each point's rate (one entry per point).
// An optional second field reuses the same mapping and weights each nodal
// contribution by point_weight; both output totals are checked.
// peak_rate bounds the active rate over the interval instead of averaging it;
// use it only for a timestep bound, never for an integrated source amount.
void assemble_point_sources(const Injection &injection, const Variables &var,
                           double start, double stop, double_vec &source,
                            const double_vec *scale = nullptr, bool peak_rate = false,
                            double_vec *weighted = nullptr, const double_vec *point_weight = nullptr);

#endif
