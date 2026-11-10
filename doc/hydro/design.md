# Poromechanics implementation and state ownership

The [overview](index.md) links the physical formulations and reproducible
benchmarks. This note describes the current code boundaries and time levels.
Pressure is compression-positive; stress and volumetric strain are
tension/expansion-positive.

## Source owners

| Source | Responsibility |
| --- | --- |
| `hydraulics.cxx/.hpp` | Pressure initialization, storage and Darcy assembly, transport substeps, accepted hydraulic history and inventory remapping. |
| `point_sources.cxx/.hpp` | Current-mesh point mapping and schedule integration shared by fluid and heat sources. |
| `unsaturated.hpp` | Retention laws and local pressure/content conversion. |
| `matprops.cxx/.hpp` | Material means, evolving porosity and permeability, and constitutive accessors. |
| `rheology.cxx` | Effective-stress increments, constitutive trials and accepted nonelastic work. |
| `fields.cxx/.hpp` | Shared shape gradients, force/velocity updates and temperature transport. |
| `pseudo_transient.cxx/.hpp` | Mechanical equilibrium trials and acceptance against the caller-owned physical baseline. |
| `output.cxx`, `remeshing.cxx` | Checkpoint and mesh-transfer integration with those state owners. |

There is one Darcy driving gradient for pressure and advected heat, one
point-source mapper, and one material-property path. Scratch assembly arrays
are not accepted physical history. Both dimensions use the same sources,
DES dimensional constants, and paired OpenMP/OpenACC loops.

## Physical-step contract

Transport advances over an accepted physical interval. Its stable substeps
advance pressure, optional matrix pressure and advected heat; they do not
repeat mechanical, porosity, plastic-strain or RSF history updates.
Scheduled source amounts are integrated over the actual subintervals.

For the fixed-stress PT path, mechanics reads the preceding pressure
increment, trials the physical constitutive interval, and commits once.
Pressure then advances from the accepted mechanical change. For the
strain-form explicit path, pressure is advanced before the force used by
the accepted momentum update; see [transport](transport.md) and
[physical inertia](dynamics.md). Coupling time accuracy must be measured
separately from transport substep stability.

The optional finite-density model preserves the current interval, old
velocity and old inertial mass during its moving-geometry predictor. The
accepted force supplies one momentum update; a later GVS mass/timestep
change belongs to the next interval. A checkpoint records whether the
pressure increment has already been consumed.

## Pressure and effective stress

`dppressure` stores the raw increment `p_old - p_new`, with the physical
interval already included. It is not a pressure rate. Each element applies
its own Biot coefficient and retention law. The Bishop increment is
`Se(p_old)*p_old - Se(p_new)*p_new`; an element saturated at both states
keeps the existing saturated arithmetic. Matching fixed-stress terms use
the same element-aware retention law.

Pressure corrections that change total stress while leaving effective
stress unchanged must not invoke a second plastic return. An increment
consumed before momentum must not be consumed again by the following
constitutive step. The checkpoint and remapping paths preserve that
ownership; see [unsaturated coupling](unsaturated.md) and the
[finite-density model](mass.md).

## PT trials and energy across timestep changes

PT iterations use the physical constitutive interval while searching for
mechanical equilibrium. A rejected trial restores the baseline. Only the
accepted caller transfers nonelastic work to the pending heat reservoir.
Hydraulic subcycling never creates another mechanical trial or deposits the
same work twice.

Pending heat is extensive energy. Changing the next physical timestep
changes its conversion to a rate, not the stored energy. Checkpoints and
remapping preserve that reservoir independently of temperature. Physical
water mass supplies the fluid heat capacity; the numerical fixed-stress
buffer contributes neither density nor sensible heat. The reference scope
is the declared sensible-energy model, not full pressure-work accounting.
See [plastic heating](plastic-heat.md), [fault heating](fault-heat.md) and
[mass and energy](mass.md).

## Mesh and restart boundaries

Remapping rebuilds current-mesh assembly and source maps after transferring
accepted state. The linear model conserves its declared Biot excess content;
the finite-density model accounts for physical water, numerical buffer,
solid heat capacity and sensible energy with explicit domain and boundary
exchanges. These are model-specific global constraints, not exact local
inventory matching. Retentive remapping is outside that qualification.

A restart must satisfy the documented dimensional, constitutive and
checkpoint-schema guards. Mechanical-method transitions have their own
restrictions and reference checks; passing an ordinary restart does not
establish an arbitrary PT/explicit conversion or a complete earthquake
cycle. Reproduction commands and the limits of each comparison are in the
[benchmark guide](../../benchmarks/poromechanics/README.md).
