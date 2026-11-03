#ifndef DYNEARTHSOL3D_GWB_HPP
#define DYNEARTHSOL3D_GWB_HPP

#include "parameters.hpp"

// Initial conditions from a Geodynamic World Builder (GWB) world, ic.gwb_filename:
// ic.temperature_option = 936 takes the nodal temperature from it, ic.mattype_option = 936
// the initial material of each marker, GWB composition i being mattype i (the world numbers
// its compositions in the cfg's material order, and must define one at every marker: a
// marker where none reaches 0.5 stops the run). Needs a build with usegwb=1;
// validate_parameters rejects option 936 otherwise. ic.weakzone_option = 936 sets the weak zone,
// weakzone_plstrain, in the elements whose centroid lies in a GWB "fault" feature.
//
// GWB is queried at a point x of the initial mesh (z = 0 at the surface): depth -z, and in
// 2D the point lies on the world's cross section, x along it. The world is opened at the
// first query and stays open until gwb_close(), which init() calls once the initial
// fields are set.

double gwb_temperature(const Param &param, const double *x);
int gwb_mattype(const Param &param, const double *x);
bool gwb_in_fault(const Param &param, const double *x);
void gwb_close();

#endif
