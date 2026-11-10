# Hydromechanics figures and analytical references

Open the generated `results/review/index.html` after building the review:
one best-verified
snapshot per benchmark, with a setup schematic, analytical equations,
reference curves, measured errors and scope. Its offline analytical appendix
includes the derivations and reference limits for every selected case.
The [review document](../../doc/hydro/index.md)
links the formulations, source ownership and reference comparisons.

`build_review.py` redraws the selected snapshots from their frozen evidence.
`review.py` owns the pages and analytical derivations, and
`setup_drawings.py` supplies the vector setup drawings;
the builder rejects a new selected case without that material. Each problem has a separate HTML/Markdown page, with short names such as
`heat.html` and `wave.html`. Setup/result figures export to PNG/SVG/PDF;
equations use individual SVGs and editable LaTeX. The heat case is described in
[groundwater heat transport](../../doc/hydro/heat.md).

`benchmarks/poromechanics/results` holds generated results and is not
part of the repository; it can be a link to a separate results directory.

## Export the checked review for the manual

After a review passes its data, figure and link checks, export it to a new
upload folder:

```sh
python3 benchmarks/poromechanics/export_manual.py REVIEW NEW_UPLOAD
```

The exporter follows the current DES3D structure: poromechanics benchmarks
by subject, groundwater/slope applications under Tutorials, and RSF/fault
heating under fault benchmarks. It keeps supporting assets in one static
folder and rewrites links between chapters. `EXPORT.json` lists every
chapter destination and file checksum. Runtime identities are removed from
portable metadata; physical checkpoint payloads remain unchanged. The
source review stays intact. Test the upload in a separate manual checkout
with the site's production build before installing it; exporting does not
publish the site.
The analysis and replay commands read results through that path.

Build both CPU executables from the repository root, then run:

```sh
make -j2 ndims=2
make -j2 ndims=3
python3 tests/functional/hydromechanics/check_theis_reference.py -v
python3 tests/functional/hydromechanics/check_terzaghi_reference.py -v
python3 tests/functional/hydromechanics/check_averaged_aquifer_reference.py -v
python3 tests/functional/hydromechanics/check_unsaturated_reference.py -v
python3 tests/functional/hydromechanics/check_wells.py ./dynearthsol2d benchmarks/poromechanics/results/wells-2d
python3 tests/functional/hydromechanics/check_wells.py ./dynearthsol3d benchmarks/poromechanics/results/wells-3d
python3 tests/functional/hydromechanics/check_poroperm.py ./dynearthsol2d benchmarks/poromechanics/results/poroperm-2d
python3 tests/functional/hydromechanics/check_poroperm.py ./dynearthsol3d benchmarks/poromechanics/results/poroperm-3d
python3 benchmarks/poromechanics/prepare_figures.py --output benchmarks/poromechanics/results/baseline
```

NumPy, SciPy and Matplotlib are required for this benchmark runner. It uses
one OpenMP thread and a 180-second timeout per DES run. An existing output
directory is refused; use a new destination for another run. To redraw the
figures without repeating the simulations:

```sh
python3 benchmarks/poromechanics/prepare_figures.py --output benchmarks/poromechanics/results/baseline --plot-only
```

Open `results/baseline/index.html` for the gallery. Every figure is saved as
PNG, SVG and PDF. The same directory includes CSV histories, `summary.json`,
the exact run configurations, complete logs and raw DES outputs. Generated
results are ignored by Git.

The [implementation guide](../../doc/hydro/design.md) describes physical-step
ownership and restart boundaries. `check_poroperm.py`
checks evolving properties and exact restart on fixed and remeshed domains.
It does not qualify conservative fluid-content transfer across remeshing.
`check_boundaries.py` verifies scheduled top/side flux and pressure, balance
and checkpoint continuation. The [explicit formulation](../../doc/hydro/transport.md)
describes the component tests and diagnostic fields.

The [refinement and multi-well qualification](../../doc/hydro/refinement.md)
records separate mesh/time trends, the small-strain superposition limits,
acceptance budgets and measured results. Run the cases sequentially:

```sh
(
  ulimit -v 6291456
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  python3 benchmarks/poromechanics/consolidation_refinement.py --output /new/refinement-results
  python3 tests/functional/hydromechanics/check_multiwell.py ./dynearthsol2d /new/multiwell-2d
  python3 tests/functional/hydromechanics/check_multiwell.py ./dynearthsol3d /new/multiwell-3d
)
```

The refinement runner saves an `index.html` gallery, PNG/SVG/PDF plots,
CSV data and checksummed provenance. Its `--plot-only` option rebuilds the
report from recorded inputs and raw saves without rerunning DES. The finest
3D case took about 324 s on the recorded host; larger studies need a separate
resource/runtime budget. All output directories must be new.

`pumping_slab.py` compares a uniform thin slab against the averaged-aquifer
reference with pressure/strain error limits and executable/config hashes.
For example, run `python3 benchmarks/poromechanics/pumping_slab.py /new/slab-results
--resolution 20` under a memory ceiling. Refinement studies remain separate
from passing this coarse comparison.

For a behavior-preserving production refactor, keep the parent output and
compare every save/checkpoint physical payload, for each dimension:

```sh
python3 benchmarks/poromechanics/compare_runs.py benchmarks/poromechanics/results/baseline/terzaghi-2d benchmarks/poromechanics/results/hydraulic-module/terzaghi-2d --report benchmarks/poromechanics/results/hydraulic-module/comparison-2d.json
```

The comparison requires identical output-file sets and exact payload equality;
only runtime/provenance records are excluded. It is separate from analytical
accuracy and conservation checks.

For the intentional Biot storage correction, preserve parent executables in a
separate directory **before rebuilding**, then run the same finite-grain column
with each pair of binaries:

```sh
python3 benchmarks/poromechanics/compare_biot_storage.py --parent-bindir /path/to/parent-binaries --output benchmarks/poromechanics/results/biot-storage-comparison
```

The default is `alpha=0.6`, with `K_s=K/(1-alpha)`. Each dimension uses the
same input before/after and is compared at its own exact saved times. Open
`results/biot-storage-comparison/index.html`; the figures, logs, hashes and
CSV histories are kept together. `--plot-only` preserves original provenance.
Use `--biot 1 --fixed-dt-seconds 0.5` in a new output directory to isolate
the small finite-grain correction from changes in automatic timestep size.
The alpha=1 run retains the original `K_s=37e9 Pa`.

Executable input validation can be checked independently in both dimensions:

```sh
python3 tests/functional/hydromechanics/check_hydraulic_input.py ./dynearthsol2d
python3 tests/functional/hydromechanics/check_hydraulic_input.py ./dynearthsol3d
```

Non-water reference densities have a separate full executable check, including
hydrostatic pressure/stress stationarity and exact continuation from checkpoint
frame 1. It retains raw data and logs in a new directory:

```sh
python3 tests/functional/hydromechanics/check_reference_density.py ./dynearthsol2d benchmarks/poromechanics/results/reference-density-2d
python3 tests/functional/hydromechanics/check_reference_density.py ./dynearthsol3d benchmarks/poromechanics/results/reference-density-3d
```

- `terzaghi_comparison`: actual DES 2D/3D output against classical Biot
  consolidation pressure and settlement, with normalized errors.
- `terzaghi_pressure_fields`: actual 2D pressure maps at three saved times.
- `theis_reference`: analytical pumping, injection and recovery curves. DES
  wells have separate executable balance and Theis comparisons in
  `tests/functional/hydromechanics/check_theis.py`; reference curves alone
  do not qualify the DES solver.
- `averaged_aquifer_reference.py`: a thickness-averaged Biot pumping reference
  for an infinite aquifer or finite annulus, with numerical Laplace inversion
  for the finite case. Bear-Corapcioglu is a limiting case. See the
  [derivation and storage convention](../../doc/hydro/aquifer.md).
  `pumping_slab.py` supplies a DES pressure/strain comparison within its
  small-strain slab and finite-domain assumptions.
- `unsaturated_reference.py`: retention laws, conservative explicit column
  updates and transient/steady infiltration references. The
  [unsaturated formulation](../../doc/hydro/unsaturated.md) describes
  retention, element-aware Bishop coupling and their separate reference
  checks, including the U0-U9 levels.
  `check_unsaturated_column.py EXECUTABLE OUTPUT_DIRECTORY` runs transient,
  steady, mechanics-mode and fixed/remeshed restart checks. The fixtures
  explicitly disable heat for these isothermal comparisons.

For a review of actual unsaturated DES fields against the references, export
the completed two-dimension verification results:

```sh
python3 benchmarks/poromechanics/unsaturated_figures.py benchmarks/poromechanics/results/unsaturated-integration-20261107 benchmarks/poromechanics/results/unsaturated-integration-20261107/figures
```

This writes an HTML gallery, PNG/SVG/PDF profiles, compact CSV data and hashes
of the consumed inputs/saves. The profile errors must reproduce the test
summaries. The export requires completed `unsaturated-2d` and `unsaturated-3d`
subdirectories and refuses an existing output directory. These plots qualify
transport comparisons; production Bishop mechanics remains unfinished.

The existing Terzaghi traction configurations are used with `PT_option=1`,
`PT_absolute_tolerance=1e-6` and explicit `gravity=10`, matching the completed
PT regression audit. Other inputs are preserved. Pressure references take
properties directly from the run configuration rather than the old notebooks.
The analytical model uses `1/M = phi/K_f + (alpha-phi)/K_s`, constrained modulus
`K+4G/3`, and diffusivity `(k/mu)/(1/M+alpha^2/(K+4G/3))`.

DES frame zero is pre-load. The reference at zero time is the undrained
post-load state, so only positive-time frames enter error metrics. Pressure
errors are normalized by the applied load; settlement errors are normalized
by final drained settlement. These coarse baseline runs precede separate
mesh/time convergence and fluid-balance qualification.

Reference: [PyLith Terzaghi solution](https://github.com/geodynamics/pylith/blob/a44a527dea7b4fb96a27e54791283c53a4a41340/tests/fullscale/poroelasticity/terzaghi/terzaghi_soln.py).

## Ordinary explicit, GVS and PT comparisons

`mechanical_methods.py` compares the same sudden-traction Terzaghi column
with PT0/PT1, ordinary explicit mechanics and GVS shear/bulk. Defaults sweep
S=1e4/1e5/1e6 and initial velocity scales 1e-12/1e-6 m/s. The column has
K/G=2 to distinguish the GVS ceilings, zero gravity, a fixed mesh and no
heat diffusion. Comparisons use actual saved physical times. This is an
exploration, not an acceptance test: a `complete` execution may have large
analytical or force errors. Failures and incomplete time coverage remain
in the report. See [method qualification and limits](../../doc/hydro/mechanics.md).

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 benchmarks/poromechanics/mechanical_methods.py benchmarks/poromechanics/results/methods-2d --dimension 2
python3 benchmarks/poromechanics/mechanical_methods.py benchmarks/poromechanics/results/methods-2d --analyze-only
```

Run dimensions sequentially. `--spacing`, `--dt-fraction`, `--duration`,
`--initial-speeds`, `--scalings` and `--methods` support controlled followups;
`--max-steps` and `--timeout` bound each run. GVS uses adaptive `fixed_dt=0`.
The output contains `index.html`, PNG/SVG/PDF error histories, CSV histories,
JSON metrics and provenance, input/executable/output hashes and raw saves.
Replay refuses changed inputs or raw saves. Reported physical-reference
kinetic and skeleton strain energies are diagnostics, not a coupled energy
balance or a substitute for force and analytical-error checks.

## Rebuild the selective review

```sh
python3 benchmarks/poromechanics/build_review.py benchmarks/poromechanics/results benchmarks/poromechanics/results/review --refresh
```

The review records its selected cases in `provenance.json`, with a plain-language
explanation and a black reference curve overlaid on DES markers. Selection
uses the predefined finest grid and target time, without parameter fitting.
It includes portable Markdown, PNG/SVG/PDF figures, compact profiles, input
decks and provenance. Unfinished applications and exploratory mechanics
results are not presented as completed application validations. Raw runs
remain available through the supporting evidence rather than the main page.

The [undrained-wave qualification](../../doc/hydro/dynamics.md)
explains the physical-inertia formulation and reproducible wave runner.
Its selected 3D result is one figure with pressure and velocity; separate
mesh/time studies remain in the evidence appendix.
The portable review includes that physical-inertia derivation and checkpoint
compatibility note under `appendix/`.
