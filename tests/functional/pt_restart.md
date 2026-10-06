# PT checkpoint and remeshing regression

From the repository root, use a CPU executable built with either `hdf5=0` or
`hdf5=1` (HDF5 comparisons require Python h5py and NumPy):

```sh
OMP_NUM_THREADS=4 python3 tests/functional/pt_initial_checkpoint.py ./dynearthsol2d tests/functional/pt_remesh_restart.cfg /tmp/des-pt-restart-2d --restart-frame 1 --require-remesh
OMP_NUM_THREADS=4 python3 tests/functional/pt_initial_checkpoint.py ./dynearthsol3d tests/functional/pt_remesh_restart.cfg /tmp/des-pt-restart-3d --restart-frame 1 --require-remesh
```

Use new output directories for each run. Omit `--restart-frame` to check the
initial equilibrium checkpoint; use `--restart-frame 2` to restart after a prior
mesh change. Logs, inputs and compared outputs remain in the output directory.
Use `--timeout SECONDS` to extend the per-process wall-clock limit on slower
backends; the default is 240 seconds. This does not change solver tolerances or
the iteration cap.

To exercise adaptive dynamic relaxation, add `PT_option = 1` to the fixture's
`[control]` section. Use `OMP_NUM_THREADS=1` for its byte-exact restart check:
the adaptive damping uses floating-point reductions whose OpenMP summation
order is not deterministic. In the 2026-10-06 four-thread remesh/hydrostatic
fixtures, exact comparisons failed while the largest field-relative difference
was 1.36e-14. Keep the exact checker strict; do not interpret that measurement as
a guarantee for other models or thread counts. See the finalization record in
`../../doc/pt-validation-20261006.md` for the tested scope.

The test checks the restored velocity before the next solve and compares every
subsequent save/checkpoint payload exactly, excluding wall
clock duration. HDF5 datasets must also have matching shapes and types and finite
numeric values. `--require-remesh` requires completed remeshing in both runs and
changed connectivity. The fixture deliberately triggers existing mesh-quality
loop-limit warnings; a pass does not establish acceptable mesh quality. The
one-year output interval avoids unrelated output-clock rounding in this short,
step-driven test. The harness suppresses extra remesh output frames so frame
numbers match physical step numbers.

## Hydrostatic initial equilibrium

```sh
OMP_NUM_THREADS=2 python3 tests/functional/pt_initial_checkpoint.py ./dynearthsol2d tests/functional/pt_hydrostatic_equilibrium.cfg /tmp/des-pt-hydrostatic-2d --require-stationary
OMP_NUM_THREADS=2 python3 tests/functional/pt_initial_checkpoint.py ./dynearthsol3d tests/functional/pt_hydrostatic_equilibrium.cfg /tmp/des-pt-hydrostatic-3d --require-stationary
```

Repeat with `--restart-frame 1` and new output directories. This no-load fixture
requires stress and pore pressure to remain within 1e-8 of their initial maximum
absolute magnitude (scale floor 1 Pa), with exactly fixed coordinates. Initial
PT equilibration retains physical sidewall support and body forces; it suppresses
transport and pressure-increment consumption, not the hydraulic load definition.
It does not establish general poroelastic analytical accuracy. Binary and HDF5
executables use the same checker; the stationarity check requires NumPy.

## Dynamic time interval and oblique boundaries

```sh
OMP_NUM_THREADS=2 python3 tests/functional/pt_initial_checkpoint.py ./dynearthsol2d tests/functional/pt_dynamic_restart.cfg /tmp/des-pt-dynamic-2d --restart-frame 1
OMP_NUM_THREADS=2 python3 tests/functional/pt_initial_checkpoint.py ./dynearthsol3d tests/functional/pt_dynamic_restart.cfg /tmp/des-pt-dynamic-3d --restart-frame 1
OMP_NUM_THREADS=2 python3 tests/functional/pt_initial_checkpoint.py ./dynearthsol2d tests/functional/pt_oblique_restart.cfg /tmp/des-pt-oblique-2d --restart-frame 1
```

The stationary elastic fixtures use automatic global velocity scaling. The axis
fixture has a nonzero mechanical response; the 2D oblique fixture uses rigid
translation to isolate type-11 boundary velocities. A repeated legacy boundary
operation must not alter accepted PT velocity, and restart must not select a new
time interval ahead of the normal
physical update cadence. These cases complement the fixed-dt schedule tests.
Relative polygon paths are resolved against the source configuration directory.

## Regular output schedule

```sh
OMP_NUM_THREADS=1 python3 tests/functional/pt_output_schedule.py ./dynearthsol2d tests/functional/pt_output_schedule.cfg /tmp/des-pt-schedule --schedule mixed
```

Repeat with `--schedule step`, `time`, and `catchup`, new output directories,
and the 3D executable. The same quiescent cell tests binary and HDF5 checkpoints.
The comparison includes frame numbers, physical steps, output times and all
subsequent save/checkpoint payloads. Cadence settings remain unchanged on restart.
PT checkpoints retain the original schedule anchors and the next regular output
index. Older checkpoints without these fields remain readable with a warning
and legacy scheduling; incomplete new schedule records are rejected. This test
does not qualify earthquake-event history or averaged-field accumulator restart.

## Nonconvergence diagnostics

A failed PT solve reports its requested stopping threshold separately from the
existing status line. The threshold is the configured absolute tolerance plus
the relative tolerance times the initial residual of that physical step. Near
equilibrium, this relative-only target may fall below attainable arithmetic
accuracy. A small reported residual alone does not authorize acceptance: failed
candidates are still rolled back, and the solver does not silently relax the
threshold. Diagnose the force scale and conditioning before changing controls.
