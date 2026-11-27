# Poromechanics benchmarks and analytical references

Inputs, analytical references and study runners live here. The complete theory,
benchmark comparisons, application figures and figure scripts live in the
[DES manual](https://geoflac.github.io/des3d/). Set `DES_DIR` to this checkout when using its
`docs/benchmarks/poromechanics/img/_src/` figure scripts.
The [implementation notes](https://geoflac.github.io/des3d/docs/theory/poromechanics) describe the solver.
Generated study results are ignored by Git. Run and check them on a compute host.

Build both CPU executables from the repository root, then run:

```sh
make -j2 ndims=2
make -j2 ndims=3
python3 tests/functional/hydromechanics/check_theis_reference.py -v
python3 tests/functional/hydromechanics/check_terzaghi_reference.py -v
python3 tests/functional/hydromechanics/check_averaged_aquifer_reference.py -v
python3 tests/functional/hydromechanics/check_wells.py ./dynearthsol2d benchmarks/poromechanics/results/wells-2d
python3 tests/functional/hydromechanics/check_wells.py ./dynearthsol3d benchmarks/poromechanics/results/wells-3d --ndims 3
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
ownership and restart boundaries. `check_boundaries.py` verifies scheduled
top/side flux and pressure, balance and checkpoint continuation. The [explicit formulation](https://geoflac.github.io/des3d/docs/theory/poromechanics)
describes the component tests and diagnostic fields.

The [multi-well qualification](https://geoflac.github.io/des3d/docs/benchmarks/poromechanics/poroelastic-consolidation/terzaghi)
records the small-strain superposition limits, acceptance budgets and
measured results. Run the cases sequentially:

```sh
(
  ulimit -v 6291456
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  python3 tests/functional/hydromechanics/check_multiwell.py ./dynearthsol2d /new/multiwell-2d
  python3 tests/functional/hydromechanics/check_multiwell.py ./dynearthsol3d /new/multiwell-3d --ndims 3
)
```

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
