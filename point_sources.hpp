#ifndef DYNEARTHSOL3D_POINT_SOURCES_HPP
#define DYNEARTHSOL3D_POINT_SOURCES_HPP

#include "parameters.hpp"

void assemble_fluid_point_sources(const Param& param, const Variables& var,
                                  double_vec& fluid_source);

#endif
