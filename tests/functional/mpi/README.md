# MPI validation and example reproduction

`run_mpi.py` runs short correctness checks and optional repeatable timings.
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

For DES/GoSPL MPI coupling with independent rank counts, see the
[coupling tests](../gospl/README.md#independent-mpi-coupling).

`usempi=1 openacc=1 openmp=0 nofma=1` builds separate `.mpi.gpu` executables.
Use an NVHPC-backed MPI wrapper and the architecture of the target device.
Initialization, MMG and field/marker remapping remain root CPU operations;
compact timestep kernels and halo pack/unpack run on GPUs with pinned host MPI
staging. Local ranks must select distinct physical GPU UUIDs. Host intranet-only
dependencies must be staged offline, and multi-CPU builds/runs on
godo/godo-s/earthsim-01 must be submitted through `qsub`. Model workspaces on
earthsim-01/nemex-gpu/kersti belong under mounted `/data/slee91/` storage.

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

`benchmark_devices.py --variants variants.json --workdir results` compares
OpenMP, CPU MPI and GPU MPI on an allocated single host. Variant records specify
`name`, `kind` (`openmp`, `mpi`, `gpu`), `workers`, `executable`, `launcher` and
`env`, keeping different MPI ABIs and library paths separate. All variants start
from the same CPU frame-zero checkpoint, giving identical geometry and fields.
One warmup and at least three randomized measured rounds are checked in full.
Actual CPU affinities, GPU UUIDs, executable/checkpoint hashes and phase/MPI
timings are retained. Total simulation time includes restart initialization and
output. Overlapping operation/communication timers are not additive phases.

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

## Spatial partitions and CPU scaling

The MPI backend now stores compact owned/ghost meshes on workers, exchanges
neighbor halos and rebuilds balanced element ownership after every MMG remesh.
Rank zero still owns global management state. See [the architecture guide](../../../doc/mpi.md)
for central operations, load metrics, trace figures and the unimplemented GPU path.

`run_mpi.py` defaults to 1, 2, 3, 4, 6 and 8 ranks and also exercises hybrid
OpenMP, averaged fields across remeshing and empty partitions. It retains every
physical field comparison. To compare equal CPU budgets with verified Linux
CPU affinity, warmups and three interleaved measured rounds:

```sh
python tests/functional/mpi/benchmark_scaling.py \
  --serial /path/to/parent/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --workdir /tmp/des-scaling
python tests/functional/mpi/plot_scaling.py /tmp/des-scaling/scaling-XXXX/results.json \
  --output /tmp/des-scaling-figures
```

`pin_cpu.py` is the benchmark's process wrapper. Each MPI rank gets one allowed
CPU ID and one OpenMP thread; the OpenMP process gets the same N allowed IDs.
CPU masks, thread counts, commands, sample times, executable hashes and physical
comparison reports are retained. Read the hardware/clock limitations in the guide;
short fixed-mesh timings do not establish full adaptive simulation performance.

`plot_partitions.py RUN [RUN ...] --output FIGURES` reads actual partition traces
and matching `run.log` files. It checks ownership against per-rank counts before
drawing 2D maps and load/ghost bars. `transport_smoke.cxx` exercises the shared
contiguous-buffer transport with unequal and zero-length messages. Its optional
CUDA mode requires a separate CUDA-aware MPI system; CPU success does not validate
GPU execution or NVLink routing.

## MPI overhead measurements

### RCB versus METIS

Use the same METIS-enabled executable for both methods, with parent executables
for physical reference outputs:

```sh
python tests/functional/mpi/benchmark_partitioners.py \
  --serial /path/to/parent/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --serial3d /path/to/parent/dynearthsol3d --mpi3d ./dynearthsol3d.mpi \
  --ranks 3 4 --repeats 5 --workdir /tmp/des-partitioners
python tests/functional/mpi/plot_partitioners.py \
  /tmp/des-partitioners/partitioners-XXXX/results.json --output /tmp/des-partitioner-figures
```

The default six cases cover small/large fixed meshes and repeated MMG remeshing
in both dimensions. `--cases` can select a subset; omit the 3D executables only
when selecting exclusively 2D cases. `--ranks` also accepts other local process
counts. At least three repeats are required. Use Linux/Open MPI with enough
allowed CPU IDs; run without competing builds or simulations.

Each method is warmed up. Within each round, process-count order and method
order within each pair are randomized. Both methods use identical pinned CPU
IDs, one OpenMP thread per process and enabled halo overlap. Actual CPU masks,
thread count and requested/actual partitioner names are verified. Every saved
frame and checkpoint in every run must be finite and exactly match the parent;
remesh and partition-generation counts must match too. A failed/incomplete run
retains `completed: false` and cannot be plotted as a successful benchmark.

Reports retain total simulation time, complete launcher walltime, all partition
rebuild times, communication timers, logical halo/root payloads, per-rank loads,
paired sample IDs, commands and source/executable hashes. Partition time includes
owner/halo-map construction and is already included in total simulation time.
The plotter also extracts initialization, timestep-computation and output times
from the retained solver logs, then writes a standalone `timings.json` with those
phases for later plotting without the original run directories. Initialization
can dominate short runs, so inspect timestep time as well as the complete run.
The communication timer on rank zero includes packing, MPI calls and waiting;
it is not an isolated network-latency measurement. Counts/bytes can be compared
directly, while elapsed time also depends on scheduling, imbalance and cache
behavior. The plotter reports medians, full sample ranges and how many paired
rounds each method won. Local shared-memory timings do not predict multi-server
network performance or establish that a custom RCB refinement is worthwhile.

### Before/after halo optimizations

Compare the preceding MPI executable with the updated executable in blocking and
overlap modes, using the same pinned CPU IDs and three interleaved measured rounds:

```sh
python tests/functional/mpi/benchmark_overhead.py \
  --before /path/to/previous/dynearthsol2d.mpi --after ./dynearthsol2d.mpi \
  --serial /path/to/parent/dynearthsol2d --workers 3 4 --repeats 3 \
  --workdir /tmp/des-overhead
python tests/functional/mpi/plot_overhead.py /tmp/des-overhead/overhead-XXXX/results.json \
  --output /tmp/des-overhead-figures
```

Do not run other simulations/builds during timing. The fixed-mesh cases use
4,001 and 64,177 elements for 1,000 steps; the small adaptive case uses 196 steps.
Every measured run checks all saved physical fields against the parent. Reports
include executable hashes, validated affinity, sample ranges, logical transfer
payloads, and per-rank operation/communication profiles. Seconds outside timed
communication include scheduling and other overhead, so they are not pure compute
times. A short test with balanced element counts does not prove balanced runtime.

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

On 2026-11-09, the 2D/3D checks passed on godo-s's two A100s (1/2 ranks),
nemex-gpu's H100 PCIe and H100 NVL (1/2 ranks), and kersti's three A6000s
(1/2/3 ranks). UUID checks confirmed distinct physical devices. The tested H100
CUDA IPC configuration returned incorrect values for larger messages; NVIDIA's
HPC-X runtime passed with `UCX_TLS=tcp,sm,self,cuda_copy`, which stages transfers
through host memory. Kersti used the same executable hashes, compiled with
`-gpu=cc86,cc90,nordc`. The A100 checks used NVHPC 21.3 and its Open MPI 3.1.5;
H100/A6000 checks used NVHPC 25.7 and HPC-X Open MPI 4.1.7rc1.

These verify GPU halo infrastructure. Full GPU MPI timesteps and CPU MMG/device
state transitions still require integration. GPU IPC, GPUDirect and NVLink routing
are not established by the passing staged transfers. No GPU speedup is claimed.

## CPU correctness through Grid Engine

On godo, godo-s and earthsim-01, submit multi-CPU work through `qsub`. This includes
MPI and multi-threaded OpenMP on one node. `launch_sge.py` requires `JOB_ID`,
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

On godo-s, submit that script with eight slots, for example
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

On godo, `qsub -S /bin/bash -q all.q -pe ompi_pf 4 -l h_rt=00:15:00 job.sh`
uses the site's round-robin PE to request slots across compute nodes. The driver
requires at least two actual allocated nodes. Its 24 comparisons cover 2D/3D
fixed meshes, forced MMG remeshing and restart, with halo overlap off/on. Fresh
frames/checkpoints 0–4 must match the parent exactly. Restart frames 3/4 must
reproduce uninterrupted physics, and every checkpoint field must match an
equivalent parent serial restart. The parent also reinitializes output scheduling
on restart, so comparing restart checkpoints to an uninterrupted checkpoint is
not a valid equality criterion; no checkpoint field is omitted to hide this.
In the tested parent, the only difference was one ULP in `reference_frame_time`
at continuation frame 3; frame 4 checkpoints were identical again. Saved physical
frames were exact throughout the checked continuation.

Retain the job ID, `PE_HOSTFILE`, job log, `qacct -j JOB_ID` and comparison JSON.
The host/slot budget checks can be exercised locally without MPI or a scheduler:
`python3 tests/functional/mpi/test_launch_sge.py`.

The compute nodes are offline. Stage interpreters and binary-compatible wheels
from the workstation when system Python is missing. Stage a serial HDF5 build
for an Open MPI solver when the site's installed HDF5 links Intel MPI. Compile
dependencies and solver variants inside a scheduler allocation as well.
Earthsim-01's scheduler was unavailable on 2026-11-09: qmaster port 6444 refused
connections after a logged controlled shutdown. Its parallel verification stays
pending until the scheduler is restored.

### Scheduler verification, 2026-11-09

Godo-s job 834 built serial HDF5 1.12.2 and parent/changed serial and MPI solvers
in both dimensions using eight allocated slots, with no compute-host downloads.
Job 835 passed all 78 CPU comparisons at 1/2/3/4/6/8 ranks, including hybrid
OpenMP, and its 79 retained launch audits all matched the eight-slot allocation.
Job 833 passed the four A100 halo checks at 1/2 ranks in 2D/3D using two slots.

Godo job 8614 passed all 24 multi-node comparisons and both standard parent
`make set` / MPI `make cmp` checks. Its four slots were allocated on
`godo_node20`, `godo_node21`, `godo_node34` and `godo_node44`; the two-/four-rank
tests used two/four actual nodes, verified by the worker audits. All successful
validation jobs have retained `qacct` records with zero failure and exit status.
Compute nodes without Python used a staged Python 3.6 interpreter, NumPy 1.19.5,
h5py 3.1.0 and `cached_property` 1.5.2 wheels prepared offline.

Reports, batch scripts, dependency hashes, accounting records and the parent
restart control are retained in `benchmarks-cores/mpi-sge.local/`. Full inputs and
snapshots remain under `/home/slee91/des-mpi-sge-20261109/` on the shared cluster
filesystem. These are correctness checks; no scaling or GPU speedup is inferred.

## CPU correctness across SSH hosts

`launch_cluster.py` plugs into the existing correctness driver when nodes have
separate working directories and do not require scheduler submission:

Use `launch_sge.py` above on godo, godo-s and earthsim-01. A static SSH host list
is not a resource allocation. If the scheduler is unavailable, retain prepared
jobs for later submission instead of running them directly.

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
Godo, godo-s, nemex-gpu, earthsim-01 and kersti are intranet machines without
internet access. Download sources, submodules and any missing dependency wheels
on the connected workstation, then transfer the complete bundle. Builds and
checks on these machines must use installed or staged dependencies without
network downloads or package-index access.
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

The 2026-11-09 nemex-gpu/earthsim-01 matrix passed all 54 physical comparisons
at 2/3/6 ranks, including two-thread hybrid runs, 2D/3D remeshing, restart,
hydraulics, averaged fields and empty partitions. All checked physical datasets
were bit-exact to parent `cce59ab`. Error termination and solver help also passed.
Parent `make set` / changed six-rank `make cmp` passed in both dimensions.
The working configuration launched from earthsim-01, used sequential placement
with nemex-gpu first, selected `ens6f0`, and used Open MPI's `ob1`/TCP transport.
Other launch directions and the kersti cross-subnet pair had TCP connectivity
failures; SSH connectivity alone does not establish MPI connectivity. Raw reports
and logs are retained in `benchmarks-cores/mpi-cluster.local/` and the isolated
remote `/tmp/des-mpi-20261109/` directories.

Those historical checks were launched directly before the user specified the
`qsub` requirement. They establish the recorded numerical agreement but are not
scheduler-submitted validation; future parallel runs must follow the requirement
above.

Godo/godo-s also passed four-rank 2D/3D parent `make set` / MPI `make cmp`
comparisons using Open MPI 1.8.7 with `ob1`/TCP over `ib0` (IP over InfiniBand).
Those builds used native DES binary output (`hdf5=0`): the installed HDF5 was
linked to Intel MPI, so it was not mixed with the Open MPI solver. Ordinary
Ethernet/default-interface attempts stalled; selecting `ib0` completed the run.

## Parallel-remeshing compatibility prototypes

These tools run outside the DES timestep backend. Production MMG, Chase's metric
and freezing rules, quality gates, interpolation and marker remapping remain in
place. METIS partitioning is independently usable with `make usempi=1 usemetis=1`
and `DES_MPI_PARTITIONER=metis`; see [the guide](../../../doc/mpi.md).

`build_remeshers.py` builds pinned ParMmg2D/ParMmg sources and their selected MMG
revisions under a supplied work directory. It never installs globally. Requirements:
C/C++ compilers, Git, CMake, MPI and METIS 5. Example with a local METIS install:

```sh
python tests/functional/mpi/build_remeshers.py \
  --workdir /tmp/des-remeshers --metis-prefix /path/to/metis
```

The pinned revisions are in the script and `build-manifest.json` records revisions,
commands, executable hashes and source-patch hashes. Build logs are retained.
The tested 2D revision uses MMG 5.7 and the 3D revision uses MMG 5.8. A separate
`--patch-2d-merge` option applies `parmmg2d-merge.patch`: a prototype fix for
one-based merge indexing and copying scalar/tensor metrics into centralized output.
It does not repair boundary labels or required-entity semantics. Use a separate
workdir when comparing unmodified and patched libraries. No library is bundled,
and no upstream issue/PR is created by these tools.

For a first external compatibility matrix using existing DES HDF5 snapshots:

```sh
python tests/functional/mpi/remesh_prototype.py \
  --input2d /path/to/2d/result.save.000000.vtkhdf \
  --parmmg2d /tmp/des-remeshers/build-2d/bin/parmmg2d_O3 \
  --mmg2d /path/to/mmg2d_O3 \
  --input3d /path/to/3d/result.save.000000.vtkhdf \
  --parmmg3d /tmp/des-remeshers/build-3d/bin/parmmg_O3 \
  --mmg3d /path/to/mmg3d_O3 --ranks 2 3 4 --cycles 2 \
  --workdir /tmp/des-remeshing-checks
```

This matrix constructs diagnostic scalar/tensor metrics and an optional frozen
cell. It checks signed element validity, area/volume, boundary references, exported
metric finiteness/positive definiteness, and required vertices/cells. Diagnostic
transfer verifies affine nodal interpolation and relocates markers while retaining
identity/material/history. Cell fields are sampled at centroids: this is **not
conservative DES field remapping**. Integral drift is reported, not treated as
proof of conservation. `--restore-boundaries2d` explicitly enables a prototype
geometric label adapter within a recorded distance tolerance; it must not be
confused with native label preservation. Results identify patched executable hashes.

For the stronger test, capture the **actual** prepared MMG inputs. Build serial
2D/3D DES objects with MMG/HDF5 in an existing build checkout, then relink them with
the test-only capture shim (Linux/GNU linker):

```sh
python tests/functional/mpi/build_capture.py \
  --build-dir /path/to/built/DynEarthSol --output /tmp/des-capture
mkdir -p /tmp/des-captured-inputs
DES_MMG_CAPTURE=/tmp/des-captured-inputs OMP_NUM_THREADS=1 \
  /tmp/des-capture/dynearthsol2d.capture model.cfg
```

`capture_mmg.c` wraps the MMG entry point and calls the original function. It
writes each prepared input mesh/metric, equivalent CLI controls and original MMG
output to a fresh `call-XXXXXX` directory. `--wrap` also works with static MMG;
LD_PRELOAD only works with dynamically linked MMG. Capture relinking reads the
real Makefile object/library settings and optional MMG VTK dependencies. Compiler
and dependency settings must match the prebuilt objects. Capture runs were checked
against the parent: all physical frames were exact in 2D/3D. The helper reproduced
all captured mesh/metric files exactly.

Replay selected calls without constructing another metric or freezing policy:

```sh
python tests/functional/mpi/replay_mmg.py \
  --captures /tmp/des-captured-inputs/call-XXXXXX \
  --mmg2d /path/to/mmg2d_O3 \
  --parmmg2d /tmp/des-remeshers/build-2d/bin/parmmg2d_O3 \
  --workdir /tmp/des-mmg-replay
```

Supply the 3D pair too when replaying 3D captures. The replay retains commands,
logs, output meshes, missing required-entity counts and `results.json`. A failed
compatibility check exits nonzero, even when the remesher returned success.
`completed: true` means the matrix finished; `all_passed` reports its outcome.
The captured inputs include the existing SDK state; different remeshers can apply
different internal policies and are not expected to produce identical meshes.

The tested parallel CLIs changed required entities on actual DES inputs. ParMmg2D
also changed boundary labels and had invalid native metric exports; serial MMG
controls preserved the tested constraints. Therefore these prototypes do not enable
parallel remeshing in production. `plot_prototypes.py --help` generates compatibility,
boundary-label, ownership and halo figures from retained results; these are not
large-machine or GPU performance benchmarks.
