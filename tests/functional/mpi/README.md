# MPI validation and example reproduction

`run_mpi.py` runs short correctness checks and optional repeatable timings.
`reproduce_examples.py` runs the existing scientific examples and retains their
configurations, checks every saved physical dataset/checkpoint against a parent
serial executable, compares complete monitor CSV records, and optionally plots
results. Both drivers require NumPy and h5py; plotting also needs Matplotlib.
Use CPU HDF5/MMG executables built from the parent and changed sources.

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
  --acc-flags=-acc --ranks 1 2 3 4 --workdir /tmp/des-accelerator-gpu
```

Choose rank counts available on that system. GPU mode requires real devices and
MPI's CUDA capability query and fails if either is unavailable. It never silently
switches to the host. The runner retains topology, compiler/device availability,
commands and logs. It verifies halo infrastructure only: a full GPU DES timestep,
CPU MMG/device-state transitions, distinct physical devices and actual NVLink
routing still require separate integration and hardware validation. The local
workspace has no exposed Intel/NVIDIA device or GPU compiler, and its MPI lacks
CUDA support; no GPU speedup is claimed.

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
