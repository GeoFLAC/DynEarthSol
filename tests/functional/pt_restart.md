# PT checkpoint and remeshing regression

## Changed timestep controls

```sh
OMP_NUM_THREADS=1 python3 tests/functional/check_restart_dt.py EXE OUT --legacy-executable OLD_EXE
```

Use a fresh output directory and executables of the same dimension and
container. Twelve tiny-cell cases check exact unchanged PT/explicit
continuation and the first interval after changing method, physical
inertia, scaling, characteristic speed, step fraction or fixed timestep.
The explicit velocity increment also checks the mass actually used.

New checkpoints retain a versioned record of the timestep controls.
Unchanged controls preserve the saved next interval; changed controls
select the current interval and matching mass before the first new step.
Older files detect PT/explicit conversion through the existing output
clock schema and honor an explicitly changed positive `fixed_dt`.
An older ordinary-explicit file does not record its previous scaling
settings: resume it with its original settings to write a current-format
checkpoint before changing those settings. GVS already reselects its
adaptive timestep/mass pair at restart.

The test supports native binary and HDF5. The recorded 2D/3D cells select
0.05/0.0333 s after a PT checkpoint carrying 1 s; the parent incorrectly
keeps 1 s. This check does not qualify an earthquake-cycle energy budget.

## Unchanged physical continuation

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
subsequent save/checkpoint physical payload exactly, excluding wall-clock duration
and per-process provenance. HDF5 datasets must also have matching shapes and types and finite
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
and the 3D executable. Use `--schedule step-default` to omit the time interval
and exercise the disabled-time default; the checkpoint clock must remain
finite. The runner selects elastic rheology so `--method pt`, `explicit` and
`gvs` exercise their requested scheduler. Add `--rsf` with `pt` or `gvs` to
retain the original RSF fixture, which forces global velocity scaling on.
Repeat PT with both options. The same cell tests binary and HDF5 checkpoints.
The comparison includes frame numbers, physical steps, output times and all
subsequent save/checkpoint payloads. The fixed physical timestep isolates the
output clock; it does not qualify adaptive GVS or changing physical timesteps.

All methods retain the original schedule anchors, next regular output index
and configured step/time intervals. Unchanged cadence preserves the original
frame sequence. An explicitly changed interval starts a new cadence at the
restart step/time. Check decreasing and increasing step intervals with
`--schedule step --step-after-restart 2` and `6`; use
`--schedule time --time-after-restart 3.7` to change the time interval to
3.7 fixed timesteps. These checks compare the physical state against a separate
run saving every step; output-clock metadata is allowed to differ.

PT's existing schedule field names remain readable. New non-PT checkpoints use
generic names; the reader detects the saved schema even for a permitted
step-zero method conversion. Older checkpoints without schedule fields remain
readable with a warning and legacy scheduling. Without saved interval fields,
an interval change cannot be detected reliably. Incomplete new records are
rejected. These checks do not qualify earthquake-event history, averaged-field
accumulator restart or changing `fixed_dt` during continuation.

## Nonconvergence diagnostics

A failed PT solve reports its requested stopping threshold separately from the
existing status line. The threshold is the configured absolute tolerance plus
the relative tolerance times the first warm-trial residual (`warm_initial`).
The reported `initial` is the selected starting candidate, which can be smaller.
`evaluations` counts full trials; `iterations` counts loop passes. Both residuals
are reported before physical transport, mesh motion and stress rotation. Near
equilibrium, this relative-only target may fall below attainable arithmetic
accuracy. A small reported residual alone does not authorize acceptance: failed
candidates are still rolled back, and the solver does not silently relax the
threshold. Diagnose the force scale and conditioning before changing controls.
