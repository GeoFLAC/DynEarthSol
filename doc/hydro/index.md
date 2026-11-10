# Poromechanics in DES

DES couples deformation, Darcy flow and heat transport through shared
material properties and accepted physical-step state. Transport remains
explicit and uses stable substeps within the mechanical interval. The
[implementation guide](design.md) describes source ownership and time levels;
[transport](transport.md), [unsaturated flow](unsaturated.md) and
[finite-density mass and energy](mass.md) give the equations and restrictions.

The [benchmark guide](../../benchmarks/poromechanics/README.md) provides
commands, inputs and independent references. A generated review gives one
selected result per problem, with a setup diagram, reference derivation,
measured error, input deck, portable figures and source hashes. The
[DES manual](https://geoflac.github.io/des3d/) is the user-facing documentation.

## Formulations and comparisons

| Subject | Notes and references |
| --- | --- |
| Groundwater flow | [Leaky aquifer](hantush.md), [several wells](theis-wells.md), [open drainage](open-flow.md), [anisotropy](anisotropic-flow.md), [scheduled sources](sources.md). |
| Unsaturated coupling | [Retention and Bishop stress](unsaturated.md), [coupled drainage](drainage.md). |
| Consolidation and inventory | [Water balance](content.md), [Terzaghi load changes](terzaghi-pulse.md), [averaged aquifer](aquifer.md), [regional subsidence](regional.md). |
| Mechanical methods | [Explicit/GVS/PT comparisons](mechanics.md), [physical inertia and method transitions](dynamics.md). |
| Strength and openings | [Effective-stress material paths](strength.md), [pressure-driven yield](pressure-yield.md), [RSF velocity steps](rsf.md), [Kirsch](kirsch.md), [Salençon](salencon.md), [borehole injection](borehole.md). |
| Heat and buoyancy | [Heat transport](heat.md), [thermal pressurization](thermal-cell.md), [surface heating](mctigue.md), [injection temperature](injection-temperature.md), [Nusselt number](nusselt.md), [linear 3D convection](convection-3d.md). |
| Plastic work and properties | [Plastic heating](plastic-heat.md), [fault-zone heating](fault-heat.md), [plastic porosity and recovery](plastic-content.md). |
| Slope applications | [Embankment strength reduction](embankment.md), [pressure-triggered sliding](slope.md), [rainfall](rainfall.md), [rapid drawdown](drawdown.md). |

## Interpreting the evidence

A parent regression measures preservation of numerical behavior. An
independent reference measures physical accuracy. Each comparison states
its normalization, spatial/time resolution, source version and validity
limits; a completed run or a closed discrete budget is not sufficient by
itself to establish an accurate solution.

The default linear storage model and optional finite-density model have
different inventory and energy meanings. The latter tracks physical water
and retained sensible heat separately from a numerical fixed-stress buffer.
Scaled fluid storage in convection studies must be retained in the matching
reference. Application comparisons retain early-time and resolution limits;
these examples do not establish complete earthquake cycles, general
retentive remapping or fully qualified thermal plumes.

Build and package a review using the commands in the benchmark guide.
`provenance.json` identifies the data and source versions used by its figures;
`EXPORT.json` records the website chapter destinations and exported hashes.
Generated fields, experimental histories and local coordination notes are
not source-repository inputs.
