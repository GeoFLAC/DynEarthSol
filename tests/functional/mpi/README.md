# MPI validation and example reproduction

`run_mpi.py` runs short correctness checks.
`reproduce_examples.py` runs the existing scientific examples and retains their
configurations, checks every saved physical dataset/checkpoint against a parent
serial executable, compares complete monitor CSV records, and optionally plots
results. Both drivers require NumPy and h5py; plotting also needs Matplotlib.
Use CPU HDF5/MMG executables built from the parent and changed sources.

## PT mechanics

```sh
python tests/functional/mpi/run_pt_mpi.py \
  --serial /path/to/current-master/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --serial3d /path/to/current-master/dynearthsol3d --mpi3d ./dynearthsol3d.mpi \
  --ranks 1 3 4 8 --plot --workdir benchmarks-cores/mpi-pt.local
```

Build all executables with HDF5 and MMG. The driver runs both upstream PT methods
in both dimensions. It reuses the upstream force/RMS and held-load/reload oracle,
100-step analytical Maxwell creep, hydrostatic initialization, running-scale
history across MMG remeshing, pending-remesh
restart, dynamic timestep and 2D oblique-boundary checks. Forced iteration-cap
failures must terminate the MPI job with status 54. Run with
`DES_MPI_PARTITIONER=metis` to repeat with METIS-enabled executables.

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
fixture do not establish acceptable production mesh quality. `pt-creep.png`
overlays the checked serial/MPI strain histories; it is a correctness figure,
not a speedup measurement. `--threads 2` checks the hybrid execution path;
ADR restart can lose byte equality with multiple OpenMP threads, as documented
by upstream, so that failure must not be hidden by relaxing the restart checker.

## Full NVIDIA GPU MPI checks

`usempi=1 openacc=1 openmp=0 nofma=1` builds separate `.mpi.gpu` executables.
Use an NVHPC-backed MPI wrapper and the architecture of the target device.
Initialization, MMG and field/marker remapping remain root CPU operations;
compact timestep kernels and halo pack/unpack run on GPUs with pinned host MPI
staging. Local ranks must select distinct physical GPU UUIDs. Stage dependencies
before using offline hosts, and use an allocation on scheduler-managed hosts.
Keep model workspaces on the site's designated storage.

Inside a single-host Grid Engine allocation, with both A100s visible:

```sh
python3 tests/functional/mpi/run_gpu_mpi.py \
  --cpu /parent/dynearthsol2d --gpu ./dynearthsol2d.mpi.gpu \
  --cpu3d /parent/dynearthsol3d --gpu3d ./dynearthsol3d.mpi.gpu \
  --ranks 1 2 --restart --workdir /shared/results \
  --launcher 'python3 /shared/source/tests/functional/mpi/launch_sge.py --local-only'
```

`--local-only` requires this host to be the sole scheduler-allocated host; it
supports MPI builds lacking Grid Engine integration while retaining slot, host,
job-ID and executable-hash checks. Native Grid Engine CPU MPI launchers use the
normal launcher without this option. No scheduler bypass is permitted.

The GPU driver checks fixed, repeated-MMG, hydraulic and averaged-output cases in
both dimensions. It compares complete finite floating datasets within `1e-8`,
integer/topology datasets exactly, every saved frame/checkpoint, and multi-GPU
outputs against the one-GPU control. Restart checks resume non-averaged MMG runs
from their own checkpoints and compare equivalent CPU continuations. Averaged
parent restarts reset their averaging window and are not included as trajectory
restart controls. The reports retain physical UUIDs and actual kernel traces.

`reproduce_examples.py --gpu` separately checks longer scientific examples
against the parent CPU. `--gpu --gpu-reference` instead checks a one-GPU control
against GPU MPI and records that reference backend explicitly. Keep both reports:
host compiler differences can change the initial Triangle mesh, and adaptive
CPU/GPU trajectories can diverge. Never relabel a failed CPU comparison as passed
because the GPU MPI comparison succeeds.

## Scientific reproductions

```sh
python tests/functional/mpi/reproduce_examples.py \
  --serial /path/to/parent/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --serial3d /path/to/parent/dynearthsol3d --mpi3d ./dynearthsol3d.mpi \
  --cases core-mmg rsf-eq maxwell-2d maxwell-3d \
  --launcher 'mpiexec --bind-to none' --plot
```

Default mode is `full`: the original simulation duration, mesh, materials,
loading, remeshing and output settings are retained. Only input asset paths are
resolved against the source directory. Runs are isolated in a fresh temporary
directory, even when `--workdir` supplies a parent directory. All overrides and
source hashes are written to `provenance.json`; executable hashes are recorded
with the commands. Full core-complex and earthquake cases can take tens of minutes
or longer. `--timeout` defaults to 7200 seconds per executable; a timeout fails the
run and terminates its process group, including MPI workers.

| Case | Source | Full duration / purpose |
|---|---|---|
| `core-mmg` | `examples/core-complex-mmg.cfg` | 500,000 yr; plastic localization and repeated MMG adaptation |
| `rsf-eq` | `examples/rate_and_state_friction/ep_rsf_eq.cfg` | 1,000 yr; velocity-weakening earthquake example with adaptive stepping |
| `rsf-creep` | `examples/rate_and_state_friction/ep_rsf_creep.cfg` | 1,000 yr; velocity-strengthening control |
| `rsf-large-dc` | `examples/rate_and_state_friction/ep_rsf_creep_largeDc.cfg` | 1,000 yr; large characteristic-distance control |
| `evp-rifting` | `examples/rifting-2d.cfg` | 2,000,000 yr plus the original 50,000 yr isostasy adjustment |
| `maxwell-2d`, `maxwell-3d` | `benchmarks/maxwell-current.cfg` | Five Maxwell relaxation times; analytical constitutive check |

For a bounded regression first:

```sh
python tests/functional/mpi/reproduce_examples.py \
  --serial /path/to/parent/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --cases core-mmg rsf-eq evp-rifting --mode regression \
  --launcher 'mpiexec --bind-to none' --plot
```

Regression caps are 20,000 core-complex steps, 2,000 RSF steps, and 1,000 rifting
steps. It changes step/output limits, preserving mesh, material and loading
settings, and keeps the original isostasy adjustment. `--max-steps` supplies a
custom cap. Such runs verify agreement only over the completed interval; an
RSF startup check is not a reproduction of earthquake cycles. `sampled_dynamic_frames`
in the report counts saved frames with model-wide speed above 1e-6 m/s and
nonzero plastic activity; it is not a complete event catalog or a recurrence test.
The monitor plot is needed to inspect the evolution between full-field outputs.

Results include `comparison.json`, a top-level `results.json`, commands, logs,
all inputs/outputs, and (with `--plot`) `reproduction.png`/`.pdf`. Core-complex and
rifting also produce `plastic-strain.png`/`.pdf`. RSF cases produce a four-panel
`monitor-history.png`/`.pdf` showing speed, stress, friction and displacement;
a dynamic burst also produces `earthquake-zoom.png`/`.pdf`. Monitor speed plots
use the existing RSF plotting script's display floor of 1e-12 m/s on logarithmic axes. Every physical VTKHDF dataset must
be finite and byte-identical, excluding wall-clock measurements. Saved frame
sets, checkpoint frame sets, monitor records, completed steps and remesh counts
must match. The normal `benchmarks-cores/compare.py` also checks the final saved
frame. Failure leaves artifacts and `completed: false`; it never marks an
unfinished run reproduced. Timings from concurrently launched manual runs are
not controlled performance benchmarks.

## Maxwell adaptation

`benchmarks/maxwell.cfg` is a historical input requiring BC type `100`; its
companion `maxwell.diff` supplied special source changes. Current DES rejects
that boundary type. The historical input and patch are retained unchanged.

`maxwell-current.cfg` is a separate, documented constitutive reproduction. It
keeps the original unit box, viscosity 1e22 Pa s, shear modulus 1e10 Pa, bulk
modulus 1e12 Pa, and imposed strain-rate components. A fixed cell and ordinary
velocity BCs impose exx=-1e-14/s and ezz=+1e-14/s (eyy=0 in 3D). Its coarse
boundary-only mesh removes interior velocity transients. With tau=eta/mu=1e12 s,

    sigma_zz(t) = 2 * eta * 1e-14 * (1 - exp(-t/tau)).

It runs 5,000 fixed 1e9-second steps through five relaxation times, saving every
100 steps. The checker requires every element's vertical stress to agree with
the analytical curve within 1e-6 relative error at each nonzero saved time,
in addition to exact serial/MPI agreement. This verifies the constitutive
response; it does not reproduce the historical deforming mesh experiment.

The standard reference workflow works in a parent checkout's `benchmarks-cores`
with `maxwell-current.cfg` copied there under a local case name: run `make set`
with the parent binary, then `make cmp FRAME=50` with the changed MPI binary.
Perform both 2D and 3D comparisons; see `AGENTS.md` and `CONTRIBUTING.md`.

## Partition and overlap checks

Workers keep compact owned/ghost meshes and exchange neighbor halos; root retains
global management state. Ownership is rebuilt after every MMG remesh. See the
[architecture guide](../../../doc/mpi.md) for central operations and scaling limits.

`run_mpi.py` defaults to ranks 1/2/3/4/6/8 and also checks hybrid OpenMP,
averaging across remeshing and empty partitions. Repeat the same physical-field
checks with `DES_MPI_PARTITIONER=metis` using a METIS-enabled executable.
`DES_MPI_OVERLAP=0` and `1` select blocking and overlapped halo assembly; use the
same reference and inputs for both. No extra simulation cfg parameter is needed.

`plot_partitions.py RUN [RUN ...] --output FIGURES` reads actual partition traces
and matching `run.log` files, checks unique complete ownership against per-rank
counts, and draws 2D ownership maps and load/ghost bars. See
[trace generation](../../../doc/mpi.md#inspecting-actual-partitions).

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

## CPU correctness through Grid Engine

On Grid Engine clusters, submit multi-CPU work through `qsub`, including MPI
and multi-threaded OpenMP on one node. `launch_sge.py` requires `JOB_ID`,
`NSLOTS` and `PE_HOSTFILE`, verifies the slot total, and lets Open MPI's Grid Engine
support select hosts and start tasks. It refuses rank/thread counts exceeding the
allocation and audits each worker's actual host, local rank count, job ID and
executable SHA-256. `sge-launch.json` retains the allocation, commands and workers.
The source, input assets, working directory and dependencies must be visible at
the same absolute paths on the allocated hosts.

For the full CPU matrix, put this command in a batch script after selecting the
matching MPI runtime and setting `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE` and
`OPENBLAS_NUM_THREADS=1`:

```sh
python3 /path/to/DynEarthSol/tests/functional/mpi/run_mpi.py \
  --serial /path/to/parent/dynearthsol2d --mpi /path/to/changed/dynearthsol2d.mpi \
  --serial3d /path/to/parent/dynearthsol3d --mpi3d /path/to/changed/dynearthsol3d.mpi \
  --ranks 1 2 3 4 6 8 \
  --launcher 'python3 /path/to/DynEarthSol/tests/functional/mpi/launch_sge.py' \
  --workdir /shared-path/validation
```

Submit that script with eight slots using the site's queue and PE, for example
`qsub -S /bin/bash -q all.q -pe smp 8 -l h_rt=00:20:00 job.sh`.
The eight-slot allocation accommodates the driver's six-rank empty-partition
case and two-rank/two-thread cases. Do not reduce slots without accounting for
those cases. Python 3.6 is supported; errors do not require `shlex.join`.

For a separate multi-node matrix, put `run_sge_matrix.py` in a batch script:

```sh
python3 /path/to/DynEarthSol/tests/functional/mpi/run_sge_matrix.py \
  --serial /path/to/parent/dynearthsol2d --mpi /path/to/changed/dynearthsol2d.mpi \
  --serial3d /path/to/parent/dynearthsol3d --mpi3d /path/to/changed/dynearthsol3d.mpi \
  --ranks 2 4 --workdir /shared-path/multi-node-validation
```

A multi-node PE, for example
`qsub -S /bin/bash -q all.q -pe ompi_pf 4 -l h_rt=00:15:00 job.sh`,
requests slots across compute nodes. Replace queue and PE names for the site.
The driver requires at least two actual allocated nodes. Its 24 comparisons cover 2D/3D
fixed meshes, forced MMG remeshing and restart, with halo overlap off/on. Fresh
frames/checkpoints 0–4 must match the parent exactly. Restart frames 3/4 must
reproduce uninterrupted physics, and every checkpoint field must match an
equivalent parent serial restart. The parent also reinitializes output scheduling
on restart, so comparing restart checkpoints to an uninterrupted checkpoint is
not a valid equality criterion; no checkpoint field is omitted to hide this.

Retain the job ID, `PE_HOSTFILE`, job log, `qacct -j JOB_ID` and comparison JSON.
The host/slot budget checks can be exercised locally without MPI or a scheduler:
`python3 tests/functional/mpi/test_launch_sge.py`.

On offline compute nodes, stage interpreters and compatible dependencies before
submission. Match dependency MPI libraries to the solver's MPI implementation.
Compile and run inside the scheduler allocation when required by the site.

## CPU correctness across SSH hosts

`launch_cluster.py` plugs into the existing correctness driver when nodes have
separate working directories and do not require scheduler submission:

Use `launch_sge.py` above on Grid Engine clusters. A static SSH host list is not
a resource allocation; follow the site's scheduler requirements.

```sh
python tests/functional/mpi/run_mpi.py \
  --serial /shared-path/parent/dynearthsol2d --mpi /shared-path/changed/dynearthsol2d.mpi \
  --serial3d /shared-path/parent/dynearthsol3d --mpi3d /shared-path/changed/dynearthsol3d.mpi \
  --ranks 2 3 6 \
  --launcher "python3 $PWD/tests/functional/mpi/launch_cluster.py --hosts host-a,host-b" \
  --workdir /tmp/des-cross-node
```

Install the solver, this launcher, input assets and a compatible MPI/runtime at
identical absolute paths on every node first; this does not require a shared
filesystem. SSH must authenticate noninteractively and verify pinned host keys.
For offline hosts, prepare sources, submodules and binary-compatible dependencies
on a connected workstation, then transfer the complete bundle before building
or running checks.
`--known-hosts FILE` selects an existing task-specific key file without changing
the user's SSH configuration. A forwarded temporary agent can supply SSH
authentication without copying private keys to workers.

The launcher checks executable SHA-256 hashes, stages the cfg and restart inputs,
and creates the working directory on each host. It collects rank-zero outputs
back to the driver host when needed. `cluster-launch.json` records commands,
input hashes, actual hosts/ranks, OpenMP thread settings, CPU masks and collection
status. Missing/duplicate ranks or fewer than two actual nodes fail the check.
`--min-nodes` changes that requirement. `--map-by seq` assigns contiguous groups
in the supplied host order; the default maps ranks by node. CPU masks are recorded
without pinning, so these runs are correctness checks, not scaling benchmarks.
