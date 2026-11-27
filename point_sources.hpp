#ifndef DYNEARTHSOL3D_POINT_SOURCES_HPP
#define DYNEARTHSOL3D_POINT_SOURCES_HPP

#include "parameters.hpp"

// Where the point sources and well-screen piece ends sit in the material, taken
// on the whole mesh before a remesh: NDIMS per point, 2*NDIMS per piece.
struct SourcePositions {
    double_vec fluid, wells;
};
SourcePositions source_positions(const Variables &var);
// Fix the point sources and well screens in the material of the whole mesh: at
// their configured places, or at the positions taken before a remesh. Then
// assemble the well source on it.
void place_sources(const Param &param, Variables &var, const SourcePositions *at = nullptr);
// The element of the whole mesh holding each point of one set and its barycentric
// weights there: at the configured coordinates, or at positions (NDIMS per point).
void place_point_sources(const Injection &injection, const Variables &var,
                         PointSourceOwners &owners, const double_vec *positions = nullptr);

// Mean source rate over an accepted (possibly subcycled) interval. Positive
// rates supply fluid volume; shape weights conserve the requested sum.
void assemble_point_sources(const Injection &injection, const PointSourceOwners &owners,
                            const Variables &var, double start, double stop, double_vec &source);

#endif
