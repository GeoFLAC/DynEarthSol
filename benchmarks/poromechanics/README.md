# Poromechanics benchmarks and analytical references

Inputs, analytical references and study runners live here. The complete theory,
benchmark comparisons, application figures and figure scripts live in the
[DES manual](https://geoflac.github.io/des3d/). Figure sources are in des3d's
`docs/benchmarks/poromechanics/img/_src/` and
`docs/tutorial/groundwater-slopes/img/_src/`; set `DES_DIR` to this checkout.
The [implementation notes](https://geoflac.github.io/des3d/docs/theory/poromechanics) describe the solver.
Generated study results are ignored by Git. Run and check them on a compute host.

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
python3 benchmarks/poromechanics/consolidation.py --output benchmarks/poromechanics/results/baseline
```

NumPy and SciPy are required for these study runners. It uses
one OpenMP thread and a 180-second timeout per DES run. An existing output
directory is refused; use a new destination for another run. To check saved runs without repeating the simulations:

```sh
python3 benchmarks/poromechanics/consolidation.py --output benchmarks/poromechanics/results/baseline --analyze-only
```

Each study saves its configurations, logs, raw fields, metrics and provenance.

The [implementation guide](https://geoflac.github.io/des3d/docs/developer/source-layout) describes physical-step
ownership and restart boundaries. `check_poroperm.py`
checks evolving properties and exact restart on fixed and remeshed domains.
It does not qualify conservative fluid-content transfer across remeshing.
`check_boundaries.py` verifies scheduled top/side flux and pressure, balance
and checkpoint continuation. The [explicit formulation](https://geoflac.github.io/des3d/docs/theory/poromechanics)
describes the component tests and diagnostic fields.

The [refinement and multi-well qualification](https://geoflac.github.io/des3d/docs/benchmarks/poromechanics/poroelastic-consolidation/terzaghi)
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

The refinement runner saves CSV metrics and checksummed JSON provenance. Its
`--analyze-only` option checks recorded inputs and raw saves.

`pumping_slab.py` compares a uniform thin slab against the averaged-aquifer
reference with pressure/strain error limits and executable/config hashes.
For example, run `python3 benchmarks/poromechanics/pumping_slab.py /new/slab-results
--resolution 20` under a memory ceiling. Refinement studies remain separate
from passing this coarse comparison.

For a behavior-preserving production refactor, keep the parent output and
compare every save/checkpoint physical payload, for each dimension:

```sh
python3 tests/functional/compare_runs.py benchmarks/poromechanics/results/baseline/terzaghi-2d benchmarks/poromechanics/results/hydraulic-module/terzaghi-2d --report benchmarks/poromechanics/results/hydraulic-module/comparison-2d.json
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
same input before/after and is compared at its own exact saved times.
JSON metrics, logs, hashes and CSV histories are kept together.
`--analyze-only` preserves original provenance.
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
  [derivation and storage convention](https://geoflac.github.io/des3d/docs/benchmarks/poromechanics/poroelastic-consolidation/aquifer).
  `pumping_slab.py` supplies a DES pressure/strain comparison within its
  small-strain slab and finite-domain assumptions.
- `unsaturated_reference.py`: retention laws, conservative explicit column
  updates and transient/steady infiltration references. The
  [unsaturated formulation](https://geoflac.github.io/des3d/docs/theory/poromechanics#partly-saturated-flow) describes
  retention, element-aware Bishop coupling and their separate reference
  checks, including the U0-U9 levels.
  `check_unsaturated_column.py EXECUTABLE OUTPUT_DIRECTORY` runs transient,
  steady, mechanics-mode and fixed/remeshed restart checks. The fixtures
  explicitly disable heat for these isothermal comparisons.

The manual figure scripts read saved unsaturated 2D/3D profiles.

## Ordinary explicit, GVS and PT comparisons

`mechanical_methods.py` compares the same sudden-traction Terzaghi column
with PT0/PT1, ordinary explicit mechanics and GVS shear/bulk. Defaults sweep
S=1e4/1e5/1e6 and initial velocity scales 1e-12/1e-6 m/s. The column has
K/G=2 to distinguish the GVS ceilings, zero gravity, a fixed mesh and no
heat diffusion. Comparisons use actual saved physical times. This is an
exploration, not an acceptance test: a `complete` execution may have large
analytical or force errors. Failures and incomplete time coverage remain
in the report. See [method qualification and limits](https://geoflac.github.io/des3d/docs/benchmarks/poromechanics/poroelastic-consolidation/terzaghi#coupled-gvs-and-pt-in-3d).

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 benchmarks/poromechanics/mechanical_methods.py benchmarks/poromechanics/results/methods-2d --dimension 2
python3 benchmarks/poromechanics/mechanical_methods.py benchmarks/poromechanics/results/methods-2d --analyze-only
```

Run dimensions sequentially. `--spacing`, `--dt-fraction`, `--duration`,
`--initial-speeds`, `--scalings` and `--methods` support controlled followups;
`--max-steps` and `--timeout` bound each run. GVS uses adaptive `fixed_dt=0`.
The output contains CSV histories, JSON metrics and provenance,
input/executable/output hashes and raw saves. Figure scripts live in the manual.
Replay refuses changed inputs or raw saves. Reported physical-reference
kinetic and skeleton strain energies are diagnostics, not a coupled energy
balance or a substitute for force and analytical-error checks.
