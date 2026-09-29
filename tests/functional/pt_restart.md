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
