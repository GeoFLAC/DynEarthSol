# Consolidation refinement and two-well qualification

The explicit core integration is followed by two independent executable
checks. These add no production solver changes. Raw outputs, configurations,
logs and executable hashes are retained under
`benchmarks/poromechanics/results/` and are ignored by Git.

## Terzaghi spatial and temporal refinement

`benchmarks/poromechanics/consolidation_refinement.py` uses the existing
independent finite-grain Biot series in `terzaghi_reference.py`. It starts
from the traction-column configurations, with a regular fixed mesh, zero
gravity, disabled thermal diffusion and PT option 1. The PT absolute and
relative tolerances are `1e-8` and `1e-12`. Material properties and the
100 kPa load come from the input configuration, not duplicated constants
in the analysis. The column is 10 m high, 1 m wide and, in 3D, 1 m thick.
Lateral rollers confine it; its base is fixed vertically and sealed, and
its top is drained and traction-loaded.

The spatial sequence uses spacings 1, 0.5 and 0.25 m with a common physical
step of 0.125 s. The temporal sequence holds the 0.5 m mesh fixed and halves
the step from 1 to 0.125 s. All runs sample exactly 100, 200, 300 and 400 s.
Frame zero is pre-load and excluded. Step-aligned output avoids comparing
slightly different physical times. Duplicate mesh/dt combinations run once.

An independent finite-element assembly reconstructs the lumped storage and
absolute element row sums. Every physical timestep is below that hydraulic
bound, so no sequence changes its temporal resolution through subcycling.
The finest 3D bound is 0.2537355 s, above the common 0.125 s step.
Input and save-file checksums are checked before replay; `--plot-only`
reads the recorded configurations and retains the original run provenance.

Results in `results/terzaghi-refinement-20261107/`:

| Dimension | Spacing (m) | Elements | Maximum pressure error / load | Maximum settlement error / drained settlement |
| --- | ---: | ---: | ---: | ---: |
| 2D | 1 | 20 | 1.77455% | 0.688285% |
| 2D | 0.5 | 80 | 0.419347% | 0.140578% |
| 2D | 0.25 | 320 | 0.0853354% | 0.0113060% |
| 3D | 1 | 50 | 2.70225% | 1.00745% |
| 3D | 0.5 | 400 | 0.683366% | 0.208055% |
| 3D | 0.25 | 3,200 | 0.124053% | 0.0207812% |

Both error measures decrease at every spatial refinement. The pressure
orders are 2.08/2.30 in 2D and 1.98/2.46 in 3D. Differences between
successive temporal solutions, normalized by the same load/drained scales,
give orders 1.00–1.03. Time plots use these differences because spatial and
temporal errors can cancel in the error against the exact solution.
The observed orders belong to this sample window and these meshes; they
are not universal convergence guarantees. Settlement orders above two can
also reflect cancellation with the remaining temporal error.

The runner requires decreasing spatial errors and decreasing successive
temporal differences, plus finest-grid pressure/load below 1% and settlement/
drained below 0.5%. Both dimensions pass. Separate PNG/SVG/PDF figures, CSV
histories, full metrics and an HTML gallery are saved with the runs.

The additional [Terzaghi loading/unloading study](terzaghi-pulse.md) checks
three material regimes through $T_v=2.5$, with its own signed load-history
reference and separate spatial/temporal sequences.

## Pumping/injection superposition

`tests/functional/hydromechanics/check_multiwell.py` runs one pumping well,
one injection well and both together on identical fixed, sealed meshes.
The shared schedule pumps/injects, shuts in, then reverses direction; its
events fall between 37 s physical steps. Each run has 600 steps and seven
saved frames. Comparisons cover pressure, stress, strain and the discrete
fluid-balance diagnostic. Independently, the integrated nodal storage and
mean-stress contribution must close against the exact scheduled well volume.
The shared balance helper uses `abs(Q)` for its normalization, so an
injection-only case cannot pass through a negative error scale.

DES still applies objective stress/strain rotation on a fixed mesh. Exact
finite-amplitude superposition therefore is not the model. An initial 2D
case at the original well-test rates gave a relative pressure discrepancy
of `1.39e-7` and stress discrepancy of `9.79e-7`. Reducing rates by 1,000
reduced those discrepancies by about 1,000, consistent with the small-strain
limit of the objective update. That initial case is retained in
`results/multiwell-2d-20261107/`; it is not counted as a passing test.

The accepted fixture uses the smaller rates and requires maximum absolute
strain below `1e-5`, relative superposition discrepancies below `1e-8` and
relative content-balance errors below `1e-9`. Results in
`results/multiwell-small-strain-{2,3}d-20261107/` are:

| Dimension | Maximum pressure discrepancy | Maximum discrepancy across checked fields | Maximum content-balance error |
| --- | ---: | ---: | ---: |
| 2D | 1.39e-10 | 9.84e-10 | 7.81e-16 |
| 3D | 6.64e-11 | 3.17e-9 | 9.39e-16 |

The largest strains are `1.51e-6` in 2D and `3.77e-7` in 3D. Both
dimensions pass. These checks do not establish Theis radial-flow accuracy,
finite-domain error, nonlinear superposition, or evolving-property/remap
conservation. The [overview](index.md) links the separate flow and
application comparisons.
