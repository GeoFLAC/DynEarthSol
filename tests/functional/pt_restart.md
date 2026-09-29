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

The test compares every final save/checkpoint payload exactly, excluding wall
clock duration. HDF5 datasets must also have matching shapes and types and finite
numeric values. `--require-remesh` requires completed remeshing in both runs and
changed connectivity. The fixture deliberately triggers existing mesh-quality
loop-limit warnings; a pass does not establish acceptable mesh quality. The
one-year output interval avoids unrelated output-clock rounding in this short,
step-driven test. The harness suppresses extra remesh output frames so frame
numbers match physical step numbers.

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
