# Experimental CPU MPI validation

Build separate MPI executables in both dimensions:

```sh
make ndims=2 usempi=1 usemmg=1 hdf5=1 BOOST_ROOT_DIR=/usr
make ndims=3 usempi=1 usemmg=1 hdf5=1 BOOST_ROOT_DIR=/usr
OMP_NUM_THREADS=1 mpiexec -n 4 ./dynearthsol2d.mpi examples/core-complex-mmg.cfg
```

The backend advances compact owned/ghost meshes, exchanges neighbor halos and
repartitions balanced element ownership after every MMG remesh. Rank zero
retains mesh/marker management, serial MMG, interpolation and output.
This is an experimental CPU implementation; GPU timesteps are not integrated.

Use parent HDF5/MMG executables as references:

```sh
python tests/functional/mpi/run_mpi.py \
  --serial /path/to/parent/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --serial3d /path/to/parent/dynearthsol3d --mpi3d ./dynearthsol3d.mpi \
  --launcher 'mpiexec --bind-to none' --ranks 1 2 3 4 6 8
```

The driver checks complete physical datasets, checkpoint/restart, remeshing,
hydraulic diffusion, hybrid threads, averaged output and empty partitions.
Wall-clock metadata is excluded. Configurations, logs and comparisons are
retained. Parent `make set` and changed `make cmp` in `benchmarks-cores/` remain
required; repeat with valid 2D and 3D cases.
