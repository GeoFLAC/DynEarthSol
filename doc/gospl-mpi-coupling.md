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
On godo/godo-s, launch inside `qsub` allocations. Earthsim scheduling is excluded
from this task. Keep modeling files on mounted `/data/slee91` on the other hosts.

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
on Linux. Inspect the resulting dependencies. Mixing Open MPI 1.8's C++ library
with conda's MPI 4 library caused an initialization crash in the retained Godo
probe; setting only runtime `LD_LIBRARY_PATH` did not fix that link.

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

## Verification on 10 November 2026

The real installed GoSPL/PETSc and C++ bridge ran on nemex-gpu's mounted storage.
Both 2D/3D CPU and GPU MPI GoSPL executables compiled. Coupled model checks are
3D and use eight DES steps/four exchanges per case:

| DES execution | DES ranks | GoSPL ranks | Check |
|---|---:|---:|---|
| CPU | 1 | 1 | reference boundary audit |
| CPU | 3 | 1 | DES rank-count equivalence |
| CPU | 3 | 2 | independent MPI sizes |
| H100 GPU | 2 | 1 | CPU/GPU coupling equivalence |
| H100 GPU | 2 | 2 | CPU/GPU coupling equivalence |
| CPU, four MMG remeshes | 3 | 2 | surface correspondence refresh |
| H100 GPU, four MMG remeshes | 2 | 2 | GPU synchronization and refreshed correspondence |
| CPU, time-triggered coupling | 3 | 2 | clock gating equivalent to the step fixture |

All eight cases passed. Across the same GoSPL rank count, all **36 complete
saved/checkpoint comparisons** passed bit for bit: DES CPU1/CPU3, CPU/GPU for
GoSPL1, CPU/GPU for GoSPL2 and time/step triggers. GoSPL np1/np2 scientific
equivalence is not the target of these tests.

Each case verified one model initialization, one elevation seed, four real
exchanges, consistent times/velocity units, exactly-once erosion application
and destruction. Seed-placement errors across owned and ghost arrays were at
most `1.39e-17 m`; returned height application errors were zero. The two GoSPL
ranks owned 58 and 72 of the 130 mesh vertices, establishing an actual distributed
model rather than duplicate full models. All 99 DES surface query points were
audited at each event. Separate per-rank audit files avoid duplicate writers.

For the fixed-mesh distributed cases, profiles record four `Surface` commands,
four local `SurfaceIdle` commands and one initial `State` distribution. There
is no repeated full-state distribution for ordinary coupling. Four real GoSPL
step calls together took about 0.13–0.37 s across these instrumented fixtures;
the tiny two-rank GoSPL case was slower than the one-rank case. This is functional
and overhead evidence, not a scaling benchmark or a reason to allocate more
GoSPL ranks to small models.

The error-path check sends an invalid operation to a two-rank worker, verifies
an error reaches the client, and requires worker status 53 with both groups
reaped. The original simultaneous-abort attempt and its MPI launcher failure
remain in the artifacts; the corrected collective-cleanup repeat is the
acceptance check.

Non-coupled regressions also passed: parent `make set` / changed four-rank
`make cmp` in 2D/3D through godo-s job **883**, plus four complete 2D/3D H100
fixed-mesh CPU/GPU checks at one/two GPU ranks. No full hardware/scaling matrix
was repeated. The subsequent [completion report](validation-2026-11-10.md)
establishes nonzero continuation, two-node coupled execution, mesh-independent
physical comparisons and a controlled coupled cost measurement. Existing standalone
DES multi-node results remain in the [earlier report](validation-2026-11-09.md).

Evidence is under `benchmarks-cores/gospl-functional.local/`: full short-run
arrays and audit records in `results/mpi-coupling-pdjlv12j/`, build logs in
`results/mpi-build-Y2j9q8/`, default checks in `results/default-gpu-regression/`,
and job 883 accounting/regression logs. Exact executable and source hashes are
retained separately from the historical November 9 GPU reports.

For the measured 204,800-tetrahedron, 100-step case on Kersti, one DES MPI rank
with four OpenMP threads and one GoSPL rank was fastest end to end: median
25.995 s, versus 51.062 s for one A6000 and 53.411 s for three. Three GPUs
slightly reduced the compute phase but added initialization/output cost. Start
small coupling jobs with that CPU configuration, reserve an additional physical
CPU for the single-threaded worker (`--worker-env OMP_NUM_THREADS=1`), and keep
both MPI worlds alive. More nodes/GPUs require a workload-specific benefit;
three-dimensional DES alone is not evidence of a speedup. The tested GoSPL/PETSc
build computes on CPUs.
