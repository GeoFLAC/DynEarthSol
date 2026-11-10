# Independent DES and GoSPL MPI coupling

GoSPL already supports MPI; its [HPC guide](https://gospl.readthedocs.io/en/v2026.7.14/getting_started/installHPC.html)
shows MPI execution. This bridge connects that existing capability to DES.
`usempi=1 usegospl=1` now supports a DES MPI world of N ranks coupled to a
persistent GoSPL world of M ranks. N and M are independent. GoSPL's installed
erosion/diffusion and native advection algorithms are unchanged.

## Execution

Build the usual GoSPL prerequisites, then both DES dimensions with an MPI wrapper:

```sh
make -j8 ndims=2 usempi=1 usegospl=1 usemmg=1 hdf5=1 MPI_CXX=/des/mpi/bin/mpicxx
make -j8 ndims=3 usempi=1 usegospl=1 usemmg=1 hdf5=1 MPI_CXX=/des/mpi/bin/mpicxx
```

GPU builds additionally use `openacc=1 openmp=0 nofma=1 GPU_CC=90` and an
NVHPC-backed wrapper for H100. Set the documented GoSPL, Python and HDF5 paths
for the host. GPU MPI uses one physical device per DES rank. The tested PETSc
environment uses CPUs; GPU acceleration here applies to DES.

Run the launcher from the model directory, using absolute paths to executables
and the real `gospl_extensions` Python interface:

```sh
/gospl/env/bin/python /src/DynEarthSol/gospl_driver/run_coupled.py \
  --gospl-python /gospl/env/bin/python \
  --interface /opt/gospl_extensions/cpp_interface/gospl_python_interface.py \
  --gospl-launcher '/gospl/env/bin/mpiexec --bind-to none -n 2' \
  --workdir /mounted/results/model \
  -- /des/mpi/bin/mpiexec --bind-to none -n 3 \
       /src/DynEarthSol/dynearthsol3d.mpi model.cfg
```

The launcher keeps the current model working directory; `--workdir` specifies
where launcher/worker logs go. YAML/mesh paths retain their usual model-relative
meaning. The DES root and GoSPL root must run on the launcher host because their
control socket is local. Explicit host/binding options belong in the respective
launcher prefixes. Reserve slots for **both** worlds, including OpenMP threads;
the launcher does not obtain extra scheduler allocations or partition CPU affinity.
Use scheduler allocations and the site's designated model storage. Both MPI
worlds must fit within the allocated host/CPU/GPU budget.

When DES and GoSPL use different MPI/SDK installations, select the worker's
libraries separately, for example:

```sh
--worker-env LD_LIBRARY_PATH=/gospl/env/lib \
--worker-env OPAL_PREFIX=/gospl/env
```

These flags are launcher arguments before `--`. DES inherits the launcher's
environment. Every GoSPL rank must see its Python environment and model files.
MPI-enabled DES requires this launcher/remote-client setup when using GoSPL;
directly embedding GoSPL in the DES MPI world is rejected at initialization.
No internet access on the compute hosts is required.

When the Python environment also contains MPI, keep the DES wrapper's library
directory ahead of that environment **at link time**, as well as at runtime:
`MPI_CXX='/des/mpi/bin/mpicxx -L/des/mpi/lib -Wl,-rpath-link,/des/mpi/lib'`
on Linux. Inspect the resulting dependencies. Mixing incompatible MPI libraries
at link time can cause initialization failures; setting only runtime
`LD_LIBRARY_PATH` does not repair the link.

## Data movement and model lifetime

The C++ bridge registers `des_gospl_client` before importing its Python interface.
That client imports no mpi4py, PETSc or GoSPL. One authenticated, private local
socket carries array messages to the GoSPL root; the worker broadcasts operations
within its own MPI world. Model processes and the connection persist for the
whole run. Neither model files nor a fresh process are created per exchange.

DES sends coordinates, velocities and elapsed time through the existing API.
The root-side surface command gathers only surface-related fields on ordinary
coupling events. Idle steps retain the original mechanical surface-volume/stress
correction on compact partitions, with neighbor halos, and do not gather state
to the DES root or call the GoSPL worker. Mesh-check marker maintenance, output
and remeshing retain their existing global management paths.

The installed `gospl_extensions` methods assume that `hGlobal.getArray()` is
a complete, sequential mesh array. In PETSc MPI it contains only owned rows.
`parallel_surface.py` adapts the coupling boundary:

- Seed elevation is interpolated to mesh coordinates, then assigned through
  `glbIDs` to owned rows and `locIDs` to local/ghost rows.
- Each rank advances its existing GoSPL partition. Owned erosion increments are
  gathered with a cached unique global-ID map; the root interpolates to the DES
  query points using a cached tree. Tectonic uplift is removed exactly once.
- Mapping buffers persist. DES surface queries are refreshed as coordinates
  move; GoSPL's surface model is not reset when DES remeshes.

This adapter also applies at GoSPL np=1, where PETSc ordering can differ from
input mesh ordering. Consequently the corrected worker route is not expected
to reproduce every result of the legacy direct embedded bridge. In the retained
fixture its first erosion array differs by up to `2.0885e-5 m` from that older
bridge. The new rank-local seed-placement checks verify the intended mapping;
the old direct embedded route retains the external extension's behavior.

Calls are synchronous at coupling events. Returned erosion shape and finite
values are checked before delivery. Coordinated API failures return an error
to DES and terminate the worker world with status 53 after collective cleanup.
The launcher reaps both child process groups. `--timeout` bounds startup and
RPC waits; `--run-timeout` optionally limits total DES execution.

## Coupled continuation

The remote worker route saves GoSPL state whenever DES writes a continuation
checkpoint. Keep the DES `.save`/`.chkpt` files and the adjacent
`<model>.chkpt.<frame>.vtkhdf.gospl/` directory together. Its immutable generation
contains rank-local model state, PETSc vector/matrix data, solver settings and
the original GoSPL mesh; a generation token in the DES checkpoint selects it.
Checksums are validated before restoring state. Binary output uses the analogous
sidecar directory without the `.vtkhdf` suffix.

Use the ordinary `sim.is_restarting`, `restarting_from_modelname` and
`restarting_from_frame` settings with the same coupled launcher. The checkpoint
also retains the pending coupling clock, step counter, velocity-reference
surface and its remesh generation. Restoration does not reseed the surface or
replay earlier GoSPL steps. A nonzero checkpoint without coupled state is rejected
instead of silently starting a new surface model.

GoSPL's rank count, partition, physical configuration and PETSc version must
match the saved state; DES's rank count can change independently. Output settings
may change. The verified configuration is the bounded tutorial-derived continental
erosion/diffusion case. Other GoSPL process combinations need their own continuation
check; this work does not establish every external model's restart behavior.
As in ordinary DES, choose checkpoints between mesh checks: output precedes a
same-step remesh, and such a checkpoint does not record the subsequent remesh.

## Verification scope

Checks cover a bounded 3D continental erosion/diffusion configuration with marine
deposition disabled. They exercise independent DES/GoSPL rank counts, CPU/GPU
execution, MMG, time/step triggering and collective error cleanup. Within the
same GoSPL rank count, 36 complete output/checkpoint comparisons were bit-exact.
The boundary audits verify owned/ghost seed placement, units, pending clocks and
exactly-once erosion application. These checks do not establish equivalence
between different GoSPL rank counts or unrelated process configurations.

Nonzero CPU/GPU continuation, including changing DES rank count across MMG,
passed 178 complete file comparisons at the existing `1e-8` relative-field
criterion; integer fields/topology were exact. A two-node coupled check passed
nine complete comparisons bit-exact. Retain the GoSPL rank count and partition
at restart. Mesh-independent sampling supports one retained CPU/GPU remeshing
case; it is not a universal remeshing tolerance.

In one fixed 204,800-tetrahedron model, four DES OpenMP threads plus one GoSPL
rank were faster end to end than one or three A6000 GPUs. Choose resources using
full launcher wall time, including initialization and output. The tested PETSc
build computes on CPUs; GPUs accelerate DES mechanics. Reserve a separate CPU
for the worker and set `--worker-env OMP_NUM_THREADS=1` and single-threaded BLAS.

Reproduce with the [quickstart](mpi-quickstart.md) and
[functional drivers](../tests/functional/gospl/README.md). Keep detailed run logs,
source/executable hashes and machine-specific records in ignored `*.local/`
directories.
