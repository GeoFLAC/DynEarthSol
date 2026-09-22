#ifndef DYNEARTHSOL3D_PHASECHANGES_HPP
#define DYNEARTHSOL3D_PHASECHANGES_HPP

void phase_changes(const Param&, Variables&);

// True when the run builds and reads the serpentinization grid (var.vtopo):
// mat.phase_change_option 3 with mat.has_serpentinization on. The one test
// for building, restoring and validating it.
bool uses_serpentinization_grid(const Param&);

// Dies (EXIT_CONFIG_VALUE) when a mat.phase_change_option 3 rule would misfire:
// two active rules from one material, or any rule from mattype_mantle.
// Called at input.
void check_phase_transition_rules(const Param&);

// The mat.mattype_* keys (without "mat.") that mat.phase_change_option 3 reads:
// the source and target of each active rule, mattype_mantle when any rule is
// active, and the grid's when uses_serpentinization_grid(). Called at input.
std::vector<const char*> phase_change_keys(const Param&);

// (Re)build the background grid for geometric serpentinization; no-op
// unless uses_serpentinization_grid().
// Called eagerly after init() and at the end of remeshing. A restart
// restores the grid from the checkpoint instead (restart() in
// dynearthsol.cxx); subduction_phase_changes() calls this only when that
// left none it can use: no grid, or no crust extents while the margin
// feature is on (a checkpoint written with the margin off).
void build_serpentinization_grid(const Param&, Variables&);

#endif
