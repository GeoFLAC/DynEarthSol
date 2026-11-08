#ifndef DYNEARTHSOL3D_GWB_HPP
#define DYNEARTHSOL3D_GWB_HPP

#include "parameters.hpp"

// Initial conditions from a Geodynamic World Builder (GWB) world, ic.gwb_filename:
// ic.temperature_option = 936 takes the nodal temperature from it, ic.mattype_option = 936
// the initial material of each marker, GWB composition i being mattype i (the world numbers
// its compositions in the cfg's material order, and must define one at every marker: a
// marker where none reaches 0.5 stops the run). Needs a build with usegwb=1;
// validate_parameters rejects option 936 otherwise. ic.weakzone_option = 936 sets the weak zone,
// weakzone_plstrain times gwb_fault_weight(), in the elements whose centroid lies in a GWB
// "fault" feature.
//
// GWB is queried at a point x of the initial mesh (z = 0 at the surface): depth -z, and in
// 2D the point lies on the world's cross section, x along it.
//
// init() calls gwb_build() before it meshes (MMG refines around the weak zone), which reads
// ic.gwb_filename and builds the world in memory if any of the three options is 936, and
// gwb_release() once the initial fields are set. A query outside that span is a bug.

void gwb_build(const Param &param);
double gwb_temperature(const Param &param, const double *x);
int gwb_mattype(const Param &param, const double *x);
// 0 outside the world's faults; inside one, 1, or with ic.gwb_fault_is_gaussian,
// exp(-d^2 / (2 sigma^2)) for the distance d from the fault's plane (the largest over the
// faults), sigma being ic.weakzone_standard_deviation.
double gwb_fault_weight(const Param &param, const double *x);
void gwb_release();

#endif
