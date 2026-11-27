# MPI validation

`run_mpi.py` runs short correctness checks against a parent serial executable and
requires NumPy and h5py. Use CPU HDF5/MMG executables built from the parent and
changed sources.

## PT mechanics

```sh
python tests/functional/mpi/run_pt_mpi.py \
  --serial /path/to/current-master/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --serial3d /path/to/current-master/dynearthsol3d --mpi3d ./dynearthsol3d.mpi \
  --ranks 1 3 4 8 --workdir benchmarks-cores/mpi-pt.local
```

Build all executables with HDF5 and MMG. The driver runs both upstream PT methods
in both dimensions. It reuses the upstream force/RMS and held-load/reload oracle,
100-step analytical Maxwell creep, hydrostatic initialization, running-scale
history across MMG remeshing, pending-remesh
restart, dynamic timestep and 2D oblique-boundary checks. Forced iteration-cap
failures must terminate the MPI job with status 54.

Every saved/checkpoint dataset is compared with the current-master serial run.
Integer topology/state and layouts require exact equality; floating fields require
`abs(error) <= 1e-12 + 1e-8 * max(abs(reference field))`. The projected PT residual
uses the loaded nodal-force scale because it measures cancellation; independent
upstream checks also require its actual RMS to satisfy the configured convergence
target. Pending pore-pressure increments use the pore-pressure level as their
scale: subtracting pressure levels near hydrostatic equilibrium can lose relative
accuracy in an almost zero increment. Byte equality is recorded separately. Each run's own fresh/restart
continuation must pass the upstream byte-exact checker; this is stricter than
the comparison between serial and MPI partitions.

The forced-remesh fixtures explicitly increase `PT_max_iter` from 5,000 to 50,000
for MMG meshes, preserving their stopping tolerances. Configurations, process
commands and logs are retained. The two remesh cases exercise boundary options
0 and 11; with MMG enabled both use MMG. Mesh-quality loop-limit warnings in this
fixture do not establish acceptable production mesh quality. `--threads 2` checks the hybrid execution path;
ADR restart can lose byte equality with multiple OpenMP threads, as documented
by upstream, so that failure must not be hidden by relaxing the restart checker.

## Full NVIDIA GPU MPI checks

`usempi=1 openacc=1 openmp=0` builds separate `.mpi.gpu` executables.
Use an NVHPC-backed MPI wrapper and the architecture of the target device.
Initialization, MMG and field/marker remapping remain root CPU operations;
compact timestep kernels and halo pack/unpack run on GPUs with pinned host MPI
staging. Local ranks must select distinct physical GPU UUIDs. Stage dependencies
before using offline hosts, and use an allocation on scheduler-managed hosts.
Keep model workspaces on the site's designated storage.

On one host where every rank sees the same GPUs, at least one per rank (local
rank N runs on device N):

```sh
python3 tests/functional/mpi/run_gpu_mpi.py \
  --cpu <parent build>/dynearthsol2d --gpu ./dynearthsol2d.mpi.gpu \
  --cpu3d <parent build>/dynearthsol3d --gpu3d ./dynearthsol3d.mpi.gpu \
  --ranks 1 2 --restart --workdir <workdir> \
  --launcher "$NVHPC_DIR/comm_libs/mpi/bin/mpiexec --bind-to none"
```

The GPU driver checks fixed, repeated-MMG, hydraulic and averaged-output cases in
both dimensions. It compares complete finite floating datasets within `1e-8`,
integer/topology datasets exactly, every saved frame/checkpoint, and multi-GPU
outputs against the one-GPU control. Restart checks resume non-averaged MMG runs
from their own checkpoints and compare equivalent CPU continuations. Averaged
parent restarts reset their averaging window and are not included as trajectory
restart controls. The reports retain physical UUIDs, and every GPU run logs its
OpenACC kernel launches: each rank's GPU must run the force update, or the check
fails. That log names no host, so the driver requires all GPU ranks on one host.

## Partition and overlap checks

Workers keep compact owned/ghost meshes and exchange neighbor halos; root retains
global management state. Ownership is rebuilt after every MMG remesh. See the
[architecture guide](../../../doc/mpi.md) for central operations and scaling limits.

`run_mpi.py` defaults to ranks 1/2/3/4/6/8 and also checks hybrid OpenMP,
averaging across remeshing, empty partitions, and the per-side BC table (letters,
tractions, an oblique boundary and a restart that changes a condition).
No extra simulation cfg parameter is needed. The `[mpi][partition]` log rows
report each rank's ownership ([reading the MPI log](../../../doc/mpi.md#reading-the-mpi-log)).

## Accelerator halo checks

The standalone harness exercises the production SoA pack/unpack indexing through
OpenACC allocations, then the shared MPI transport. CPU-only validation is explicit:

```sh
python tests/functional/mpi/run_accelerator_checks.py --mode host \
  --workdir /tmp/des-accelerator-host
```

This builds both dimensions as C++17 and checks 1/2/3/4/6/8 ranks with GCC's
`-fopenacc`. Use `--std c++20` for the C++20 compatibility run.
It records `gpu_verified: false`. On a separate NVIDIA system with NVHPC and
CUDA-aware Open MPI, expose all local GPUs and use one rank per GPU:

```sh
OMPI_CXX=nvc++ python tests/functional/mpi/run_accelerator_checks.py --mode gpu \
  --ranks 1 2 --workdir /tmp/des-accelerator-gpu
```

GPU compilation defaults to `-acc -cuda`; host compilation uses `-fopenacc`.
`--acc-flags` overrides these for a particular compiler/architecture. Choose rank
counts available on that system. GPU mode requires real devices and MPI's CUDA
capability query and fails if either is unavailable. It never silently switches
to the host. Every local rank selects its GPU before MPI initialization, and the
test verifies placement using a shared-memory communicator and physical CUDA
UUIDs. It exercises small and large irregular buffers, asymmetric empty messages,
and scalar/vector/tensor widths in both dimensions.

Reports retain compiler/run commands, executable and local source hashes, topology,
and every rank's host, local rank, UUID and PCI address. `--min-nodes 2` additionally
requires two distinct processor names in every invocation. `--prebuilt-dir DIR`
uses existing `halo-2d`/`halo-3d` binaries on a host without a compiler; retain the
producer's build report and verify that executable hashes match. Compatible MPI,
CUDA and NVHPC runtime libraries must already be available on the target.

These checks verify halo infrastructure. Use the full GPU MPI driver above for
timestep and CPU MMG/device-state transitions. Staged transfers do not establish
GPU IPC, GPUDirect or NVLink routing, and passing checks do not imply speedup.
