# Experimental MPI backends

`usempi=1` builds separate `dynearthsol2d.mpi` / `dynearthsol3d.mpi` executables.
Regular timestep kernels advance compact spatial partitions and exchange neighbor
halos. The project defaults to C++17 and also supports C++20, matching upstream.
Temperature, force and mass assembly can overlap element halo transfers with
computation on interior nodes. `usempi=1 openacc=1` additionally builds
`dynearthsol2d.mpi.gpu` / `dynearthsol3d.mpi.gpu` for NVIDIA devices.

## Build and run

```sh
make -j6 ndims=2 usempi=1 usemmg=1 hdf5=1
make -j6 ndims=3 usempi=1 usemmg=1 hdf5=1
cd benchmarks-cores
OMP_NUM_THREADS=1 mpiexec -n 4 ../dynearthsol2d.mpi test-tiny.cfg
```

`MPI_CXX=/path/to/mpicxx` selects the wrapper. Configuration parsing uses
the vendored CLI11 header; Boost is no longer required. Existing cfg files
apply without new simulation parameters. One CPU rank uses the ordinary
full-range kernels. A one-rank GPU run uses the same compact
backend as multi-GPU runs. NVHPC profiling cannot currently be combined with
`usempi=1`. `usempi=1 usegospl=1` couples to an independent persistent GoSPL
MPI worker world; see [MPI coupling](gospl-mpi-coupling.md) for the launcher and
separate DES/GoSPL rank counts. Serial/OpenMP and non-MPI OpenACC executables
remain separate.

## NVIDIA GPU MPI

The completed intranet validation, measured overheads and remaining limitations
are recorded in [the November 2026 report](validation-2026-11-09.md).

Build with an NVHPC-backed MPI wrapper; submit builds and runs through `qsub` on
godo/godo-s/earthsim-01. The compute hosts have no internet access: stage complete
dependencies and runtimes first. On earthsim-01, nemex-gpu and kersti, model files,
checkpoints and results belong on mounted `/data/slee91/` storage.

```sh
export OMPI_CXX=/path/to/nvc++
make -j8 ndims=2 usempi=1 openacc=1 openmp=0 nofma=1 GPU_CC=80 \
  MPI_CXX=/path/to/mpicxx usemmg=1 hdf5=1
make -j8 ndims=3 usempi=1 openacc=1 openmp=0 nofma=1 GPU_CC=80 \
  MPI_CXX=/path/to/mpicxx usemmg=1 hdf5=1
```

Use one rank per visible physical GPU. Device selection precedes managed
allocations; local ranks must have distinct physical UUIDs. The implementation
packs compact SoA halo rows on devices, copies contiguous buffers through pinned
host memory, and supplies **host pointers** to MPI. GPU producers finish before
packing; receives finish before upload/unpack and buffer reuse. This path needs
ordinary MPI transport, and makes no claim about NVLink, CUDA IPC or GPUDirect.

Compact timestep fields use managed allocation and are prefetched after
distribution. Initialization, MMG, double-precision nearest-neighbor remapping,
global marker/surface management and output run on the root CPU. Conditional
OpenACC regions enforce this boundary; command boundaries drain pending kernels
before host state changes. MMG is followed by collection, repartitioning and
managed-state refresh. Managed allocation also covers the root state container.

GPU MPI builds use IEEE division; `nofma=1` disables host and device FMA. NVHPC
21.3 was exercised on A100/A6000 in the earlier source matrix; H100 was built natively with NVHPC 25.7 targeting
`GPU_CC=90` after the staged 21.3 runtime could not enumerate it. Keep each
MPI/SDK runtime separate when staging libraries. In
particular, an HDF5 library directory containing an older SDK can accidentally
make a newer compiler link the old OpenACC runtime.
The current C++17 integration and reviewer GPU examples were rerun with NVHPC
25.7 on A6000; use the project's supported SDK/toolchain for reproduction.

The validation drivers retain full-field and checkpoint comparisons at `1e-8`,
integer/topology equality, GPU identities and actual kernel traces. Short 2D/3D
fixed, MMG, hydraulic and averaged-output cases exercise GPU MPI. Non-averaged MMG
checkpoints also exercise each run's own restart. The parent resets averaging
windows at restart, so averaged-output trajectory restart is outside this matrix.
See the test README for commands and retained reports. Longer CPU/GPU comparisons
must be assessed separately: host compilers can generate different Triangle
meshes, and adaptive trajectories can diverge. Once meshes differ, compare
physical observables or interpolate fields to common spatial locations rather
than require equal node/element arrays. Array-comparator rejection after remeshing
alone does not identify a solver defect. A passing short matrix does not
establish complete scientific equivalence or a speedup.

## Physical-step PT mechanics

Both upstream PT methods run on the compact compute meshes: localized relaxation
(`PT_option=0`) and adaptive dynamic relaxation (`PT_option=1`). Initial body-force
equilibration uses the same distributed path. Each rank snapshots and restores its
owned material history; no mechanical baseline crosses a remesh. Trial velocities
and element assembly/factor data use neighbor halos. Physical mesh motion, transport
and remeshing still happen only after convergence.

Residual norms count owned free degrees of freedom once, using a global force
scale before squaring. The squared sum and free-rank count share one reduction;
ADR's three secant sums share another. Candidate selection, nonfinite rejection
and convergence therefore make the same decision on every rank, including empty
partitions. The accepted diagnostic, running-scale history and post-remesh exclusion
are retained for output/restart. Iterations do not gather the global mesh.

ADR changes floating-point summation order across partitions. Parent comparisons
must measure field error with an explicit tolerance; byte-exact agreement is
reported only when observed. One thread per rank avoids an additional source of
reduction-order changes. This does not guarantee that every nonlinear PT problem
converges; the physical and numerical limits in [pt-solver.md](pt-solver.md) apply.

The MPI PT regression reuses upstream's independent force, Maxwell creep, initial
equilibrium and restart checks. See [the test guide](../tests/functional/mpi/README.md#pt-mechanics).

## Ownership, halos and remeshing

### Optional METIS partitioning

RCB is still the default. Build METIS support explicitly, then select it through
the MPI environment (no new simulation cfg parameter):

```sh
make ndims=2 usempi=1 usemetis=1 usemmg=1 hdf5=1
DES_MPI_PARTITIONER=metis OMP_NUM_THREADS=1 \
  mpiexec -n 4 ./dynearthsol2d.mpi model.cfg
```

METIS 5 headers/libraries are required. `METIS_INCLUDE_DIR` and `METIS_LIB_DIR`
override compiler/linker default paths; `make config` prints them. Repeat the
build with `ndims=3` for 3D. Requesting METIS from a build without it fails clearly;
unknown partitioner names fail too. One-rank execution uses the original kernels
without partitioning. Tiny meshes with no more elements than ranks use an
explicitly logged RCB fallback, preserving empty-partition support.

The METIS graph connects elements sharing any vertex, matching DES node-star
dependencies more closely than a face-only graph. K-way partitioning uses the
communication-volume objective, seed 0 and a 3% element-balance tolerance. Node
ownership and halo construction then use the same code as RCB. The graph objective
is a proxy: actual halo bytes also depend on node ownership and exchanged fields.
METIS does not currently use measured constitutive/marker weights or migrate
incrementally. It is rebuilt after each MMG remesh, as RCB is.

Both dimensions compiled with/without METIS and in serial. All 78 parent comparisons
passed with METIS at 1/2/3/4/6/8 ranks, including averaged remeshing. Parent `make set`
and changed `make cmp` passed in 2D/3D. A failed 3D averaged-remeshing case exposed
an old-size inactive hydraulic scratch array; resizing it on state distribution
fixed the bounds error and the case passed under AddressSanitizer.

On the small 40-step comparison, final total ghost-element counts for RCB/METIS
were 44/35 (2D, 3 ranks), 50/54 (2D, 4 ranks), 257/264 (3D, 3 ranks), and 449/401
(3D, 4 ranks). Every physical output matched. Thus METIS did not uniformly reduce
the halo. These are structural comparisons, not performance claims. Figures and
raw reports are retained in `benchmarks-cores/mpi-prototypes.local/` (ignored).

For a runtime comparison, use `benchmark_partitioners.py` and
`plot_partitioners.py` as described in the
[test guide](../tests/functional/mpi/README.md#rcb-versus-metis). Both methods use
the same MPI executable, one thread per rank, identical verified CPU masks,
warmups and randomized paired rounds. Every saved field/checkpoint is checked
against the parent. Reports distinguish total simulation time, complete launcher
walltime, partition rebuild cost and logical halo payload. A smaller halo alone
does not establish a faster run, and local shared-memory results do not predict
multi-server network performance.

#### Controlled laptop comparison (2026-11-06)

Six inputs cover small/large fixed meshes and repeated remeshing in 2D/3D, using
three/four ranks and five paired measured rounds per method. All 120 measured
runs and 24 warmups matched all saved physical datasets and checkpoints exactly.
Median total simulation times for the larger meshes were:

| Case | Ranks | RCB (s) | METIS (s) | METIS halo-byte reduction |
|---|---:|---:|---:|---:|
| 2D, 64,177 elements, 1,000 steps | 3 | 4.5796 | 4.7473 | -1.2% |
| 2D, 64,177 elements, 1,000 steps | 4 | 4.1035 | 4.1359 | 9.6% |
| 3D, 36,549 elements, 1,000 steps | 3 | 12.1349 | 12.1964 | 19.2% |
| 3D, 36,549 elements, 1,000 steps | 4 | 12.2053 | 12.2353 | 15.2% |

Reduced halo traffic did not produce a consistent runtime advantage. Three-rank
large-3D paired results favored METIS in four of five rounds, even though the
ratio of separate method medians slightly favored RCB; both statistics and the
full sample ranges are retained. Close results therefore do not establish a
general winner. Large-3D setup accounts for roughly 70–75% of these short runs;
the separate timestep-computation medians also favored RCB. RCB partitioning was
cheaper in every case. Repeated-remeshing total-runtime medians favored RCB in
both dimensions, with the 2D case faster in all five pairs at both rank counts.

Keep RCB as the default and METIS optional. These measurements do not justify
adding a custom connectivity-refinement step to RCB. Larger adaptive workloads
and actual multi-server networks may have different tradeoffs. Figures, paired
samples, phase timings, exact comparisons and compact audit logs are retained in
`benchmarks-cores/mpi-partitioners.local/` (ignored). Solver C++ and binaries were
unchanged for this comparison.

### Default geometric partitioning

Recursive coordinate bisection splits element centroids along the longest spatial
axis. Element counts differ by at most one, including non-power-of-two process
counts and empty partitions. Interior nodes belong to the incident element rank;
shared nodes are assigned to an incident rank with the smallest current node load.
This balances element counts first. It does not guarantee equal node counts,
marker counts, constitutive work or execution time.

Every compute rank stores owned nodes/elements first, followed by ghost rows.
The local mesh includes owned elements, the complete element stars of owned nodes,
and their vertices. Local indices map to stable global IDs until remeshing.
Node patch sums retain global element ordering to preserve the serial arithmetic.
Workers do not receive the full mesh. Rank zero also retains a global management
copy, so root memory and centralized operations remain scaling limits.

Ownership and halo maps are rebuilt after **every** MMG remesh, even if node and
element counts are unchanged. Owned and ghost marker subsets are reconstructed,
including marker history and barycentric coordinates, with consistency checks.
MMG and interpolation run centrally, preserving existing 2D/3D behavior without
requiring ParMMG. Repartitioning occurs after remeshing, not periodically between
remeshes. Cost-weighted partitioning remains a possible improvement if measured
load warrants it; optional METIS already supplies a communication-volume objective.

## Distributed and central work

Distributed kernels include temperature/pore-pressure updates, strain rate,
constitutive stress/strain/plasticity and friction, force, velocity, coordinates,
stress rotation, volumes, volume-change rates, mixed discretization, masses and
timestep bounds. Boundary conditions use global wall extents; node owners supply
authoritative values. Residual norms use a global sum, whose rounding can vary
with process count.

Rank zero drives initialization/isostasy, scheduling, phase changes, hydrous-marker
advection, marker creation/correction, MMG/interpolation, surface management and
output. Workers wait during this work. Thus **not every operation is parallel**.
A global state is collected only when needed by these operations, checkpoints,
output or optional earthquake moment accumulation. Ordinary surface options 0/1
collect/scatter a top/boundary strip. Quality checks collect coordinates and
element volumes; advanced surface work can require full state. The optional
serpentinization grid collects state only when that model feature is active.
Monitor output collects queried rows and their dependencies.
Output averaging keeps local accumulators and gathers them when writing.
Earthquake speed/plastic-rate diagnostics use small global maximum reductions.

Regular dependency exchanges pack only needed ghost fields into reusable
contiguous peer buffers. Receives/sends use `MPI_Irecv`/`MPI_Isend`, followed by
`MPI_Waitall` and unpacking. Related fields share transfers; buffers and maps
persist until repartitioning. With `DES_MPI_OVERLAP=1`, temperature, force and mass
assembly first posts element transfers, computes owned nodes whose incident
elements are all local, waits/unpacks, then computes interface nodes. Per-node
element summation order is unchanged. Overlap is enabled by default.
`DES_MPI_OVERLAP=0` uses a blocking exchange
and one owned-node pass. Other exchanges still complete before dependent work.
MPI progress during computation depends on the MPI implementation; posting a
transfer does not by itself prove simultaneous data movement. Full-field
`MPI_Allgatherv` is no longer used on each regular timestep operation.

## Inspecting actual partitions

Set an output prefix before launching; its parent directory must exist:

```sh
mkdir -p /tmp/des-partition-run
DES_MPI_PARTITION_TRACE=/tmp/des-partition-run/partition OMP_NUM_THREADS=1 \
  mpiexec -n 4 ./dynearthsol2d.mpi model.cfg > /tmp/des-partition-run/run.log 2>&1
python tests/functional/mpi/plot_partitions.py /tmp/des-partition-run \
  --output /tmp/des-partition-figures
```

Trace files `partition.0001.csv`, `.0002.csv`, etc. contain actual global IDs,
owners, node positions, centroids and connectivity at each partition generation.
`[mpi][partition]` log rows report step, owned nodes/elements/markers, ghosts and
neighbors for every rank. The plotter validates complete, unique ownership and
matches trace counts against log rows before making figures and metrics JSON.
The plotter currently draws 2D triangular meshes; 3D traces need a separate viewer.
Plots show count balance and ghost storage; they do not establish equal runtime.

`[mpi][kernel]` reports operation counts and root elapsed seconds including
synchronization. Management operations have counters too. `[mpi][communication]`
reports timed communication, `halo_sent_bytes` summed over ranks,
`root_collection_bytes` and `state_bytes` sent to workers. These logical payload
counts exclude headers, reductions, MPI protocol traffic and local root copies;
they are not a complete network profile.

`[mpi][rank-profile]` reports operation time and instrumented communication time
on every rank. Their difference includes local work, scheduling and uninstrumented
overhead; it is not a hardware compute counter. `[mpi][overlap]` reports aggregate
wait time and the interval between posting and completion, including blocking
exchanges' small intervals. These counters help diagnose imbalance but do not
establish asynchronous MPI progress. Reductions are now included in communication
timing, so seconds cannot be compared directly with earlier instrumentation.

## Correctness checks

For godo/godo-s/earthsim-01, use the
[Grid Engine workflow](../tests/functional/mpi/README.md#cpu-correctness-through-grid-engine).
The scheduler-submitted 2026-11-09 checks passed 78 CPU comparisons on godo-s and
24 two-/four-node comparisons on godo, plus the standard 2D/3D comparisons. Rank
audits verified scheduler job IDs, allocated hosts/slots and executable hashes.
The two A100s also passed their halo checks through `qsub`. Earthsim-01's qmaster
was unavailable, so scheduler-submitted checks there remain pending.

Multi-node CPU and physical GPU halo checks are now available through the
[cluster validation commands](../tests/functional/mpi/README.md#cpu-correctness-across-ssh-hosts).
The 2026-11-09 two-node CPU matrix passed 54 exact parent comparisons and the
2D/3D `make set`/`make cmp` workflow. A100, H100 and A6000 halo checks passed in
both dimensions with distinct GPU UUIDs. H100 larger-message checks required
disabling CUDA IPC and using host staging. These earlier halo checks are
correctness results, not full-solver scaling or NVLink measurements; the NVIDIA
GPU MPI section describes the subsequent solver integration.

Use HDF5/MMG executables built from the parent commit and changed sources:

```sh
python tests/functional/mpi/run_mpi.py \
  --serial /path/to/parent/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --serial3d /path/to/parent/dynearthsol3d --mpi3d ./dynearthsol3d.mpi \
  --launcher 'mpiexec --bind-to none' --ranks 1 2 3 4 6 8
```

The 3D pair is optional. Results, logs, configurations and snapshots are retained
in a fresh directory. The driver uses the defaults, anisotropy, core-complex,
core-complex/MMG and subduction examples, capping duration and adjusting output.
Additional cases cover hybrid MPI/OpenMP, hydraulic diffusion, pseudo-transient
iterations, thermal-off, hydration, repeated remeshing, restart, isostasy, averaged
outputs, empty partitions and error termination. At least eight local MPI slots
are needed for the default rank matrix.

Saved physical VTKHDF datasets, including complete vectors/tensors, connectivity
and marker records, must be finite and byte-identical to the parent; walltime is
excluded. The standard comparator also runs. The repository reference workflow
remains required:

```sh
# Parent checkout, before changing sources:
cd benchmarks-cores
make set NDIMS=2 CASE=test-tiny.cfg
# Changed executable, reusing that reference:
make cmp NDIMS=2 CASE=test-tiny.cfg \
  EXE='mpiexec -n 3 /path/to/changed/dynearthsol2d.mpi'
```

Repeat with a valid 3D case. Output precedes same-step remeshing in DES; restarting
such a checkpoint can omit that remesh, also without MPI. Restart tests therefore
place checkpoints between mesh checks.

The [scientific reproduction guide](../tests/functional/mpi/README.md) covers
full-duration core-complex/MMG, earthquake cycles, bounded rifting, and a separate
current-format Maxwell benchmark with an analytical constitutive comparison.
Full examples and bounded regression runs are explicitly distinguished.

## OpenMP versus MPI on one laptop

Use equal CPU budgets: N OpenMP threads versus N MPI ranks with one thread each.
The following Linux/Open MPI driver pins the OpenMP process to the first N allowed
CPU IDs and each MPI rank to one of those same IDs. It records and validates actual
CPU masks and runtime thread counts, warms every variant, interleaves randomized
measured rounds, and checks physical outputs against the one-thread parent.

```sh
python tests/functional/mpi/benchmark_scaling.py \
  --serial /path/to/parent/dynearthsol2d --mpi ./dynearthsol2d.mpi \
  --workers 1 2 3 4 6 8 --repeats 3 --workdir /tmp/des-scaling
```

No competing builds/simulations should run during timing. The fixed-mesh inputs
derive from core-complex/MMG, with 1,000 steps at 500 m / 125 m resolution.
Reports retain simulation time and complete walltime including MPI launch.
These short cases isolate timestep scaling; they do not predict complete adaptive
run performance. On hybrid CPUs, Linux/WSL CPU IDs do not necessarily expose the
physical P/E/LP core type. Affinity gives repeatable CPU sets, not a guarantee of
identical physical cores or stable laptop clocks.

### Laptop measurement (2026-11-05)

These timings are historical measurements on the pre-CLI11/C++17 source base
`cce59ab`; they are not a new timing study of the current upstream integration.

Intel Core Ultra 7 356H under WSL, GCC 13.3, Open MPI 4.1.6, MMG 5.8,
system Boost/HDF5, parent `cce59ab`. Median simulation seconds from three
interleaved runs; launcher startup excluded. Both modes use the same first N
allowed CPU IDs. Full sample ranges and launcher walltimes are in the local JSON.

| CPU budget | OpenMP: 4,001 elements | MPI: 4,001 elements | OpenMP: 64,177 elements | MPI: 64,177 elements |
|---|---:|---:|---:|---:|
| 1 | 0.495 | 0.499 | 8.877 | 8.526 |
| 2 | 0.303 | 0.331 | 4.961 | 5.938 |
| 3 | 0.244 | 0.269 | 3.996 | 4.717 |
| 4 | 0.208 | 0.242 | 3.366 | 4.228 |
| 6 | 0.218 | 0.252 | 3.173 | 4.659 |
| 8 | 0.226 | 0.249 | 3.073 | 3.824 |

OpenMP wins at every tested parallel CPU budget. The fastest small case uses four
OpenMP threads; the fastest large case uses eight. Laptop clock variation is
visible in the sample ranges; these are model-specific measurements.

In a separate, unpinned three-repeat comparison, the preceding replicated MPI
backend took 13.445 s on two ranks for 64,177 elements, versus 5.676 s with neighbor
halos. Logical payload fell from approximately 23.45 GB of field/state transfers
to 0.89 GB of halo/root/state transfers. This is an approximate aggregate payload
comparison, not a network byte measurement; the two implementations count different
transfer categories. Do not compare those unpinned times directly with the table.

The final short-case matrix passed 78 configuration comparisons against the parent,
including averaged remeshing in 2D/3D. Every checked physical dataset was exact.
Separate parent `make set` / changed `make cmp` checks passed in 2D and 3D, as did
20 forced MMG remeshes at three/four ranks. A further 50,000-step continuation of
the mature core-complex model passed every saved/checkpoint comparison at three
and four ranks, including two natural remeshes. Final element ownership was
1,165/1,165/1,165 and 873/874/874/874, respectively. These count balances do not
prove equal per-rank execution time. The complete 1,000-year earthquake example
also passed on four ranks (98,378 steps, 40 saved frames, 12 checkpoints and
complete monitor records), as did the 2D/3D Maxwell constitutive checks and the
bounded rifting case on three ranks.

Local figures/results are in `benchmarks-cores/mpi-halo.local/` (Git-ignored).
Use the measured best configuration for each model. More workers can be slower
because of communication, memory traffic and the laptop's mixed core speeds.

### Overhead measurement (2026-11-06)

`benchmark_overhead.py` compares the preceding neighbor-halo executable with the
updated executable in blocking and overlap modes. The same CPU IDs, one thread
per rank, warmups and three randomized interleaved rounds apply. All 54 measured
runs matched every saved physical field exactly. Median simulation seconds:

| Case | Ranks | Previous MPI | Updated blocking | Updated overlap |
|---|---:|---:|---:|---:|
| 4,001 elements, fixed, 1,000 steps | 3 | 0.2644 | 0.2528 | 0.2568 |
| 4,001 elements, fixed, 1,000 steps | 4 | 0.2334 | 0.2251 | 0.2267 |
| 64,177 elements, fixed, 1,000 steps | 3 | 4.9237 | 4.4548 | 4.2908 |
| 64,177 elements, fixed, 1,000 steps | 4 | 4.5710 | 4.2345 | 4.1711 |
| Small adaptive case, 196 steps | 3 | 0.3274 | 0.3291 | 0.3289 |
| Small adaptive case, 196 steps | 4 | 0.3404 | 0.3314 | 0.3320 |

For the larger fixed mesh, runtime fell by 12.9% at three ranks and 8.7% at four.
Logical root collection payload fell from 1,123.7 to 22.0 MB and from 1,263.9
to 24.8 MB, respectively (98.0%); halo bytes were unchanged. The fixed mesh needs
no geometry quality check, so skipping the inactive serpentinization collection
accounts for that payload reduction. These figures measure the combined changes,
not the isolated effect of any one optimization. Overlap improves the larger case
over the updated blocking variant, but costs a little on the small case. The
adaptive case shows little change; it does not establish a general adaptive-run
speedup. Four-rank sample ranges overlap substantially on this laptop.

For the larger case, the maximum per-rank median time outside instrumented
communication was within 1.0% of the rank mean for stress, temperature and force
at both rank counts. The short adaptive case was within 2.9% for the plotted
kernels, over only milliseconds of work. These measurements do not justify adding
constitutive weights for these inputs. Other materials, GPUs and larger adaptive
meshes need their own profiles. The new figures, sample ranges and raw JSON are
in `benchmarks-cores/mpi-overhead.local/` (ignored); see the
[timing commands](../tests/functional/mpi/README.md#mpi-overhead-measurements).

The updated executable also matched the complete earthquake example (98,378
steps, 40 saved frames, 12 checkpoints, 9,838 monitor records) and 50,000-step
mature core-complex continuations at three/four ranks, each with two natural MMG
remeshes. Across these longer runs, the sum of stress/temperature/force/mass time
outside instrumented communication had maximum/mean rank ratios of 1.005, 1.009
and 1.008, respectively. Individual kernels can differ more (core-complex mass
recalculation reached 1.076 on four ranks); the aggregate measurements still do
not justify adding material weights for these examples.

## Accelerator preparation and hardware choice

`mpi_transport.hpp` accepts contiguous buffers and never dereferences their
pointers. The CPU solver uses host buffers. The standalone transport test exercises
unequal and asymmetric zero-length messages with the same helper:

```sh
mpicxx -std=c++17 -O2 -DOMPI_SKIP_MPICXX=1 -I. -Inanoflann/include \
  tests/functional/mpi/transport_smoke.cxx -o /tmp/des-transport
mpiexec -n 4 /tmp/des-transport
```

A `-DDES_TEST_CUDA` build additionally needs CUDA headers/runtime libraries and
CUDA-aware Open MPI. It selects one GPU per local rank, rejects duplicate device
UUIDs, checks `MPIX_Query_cuda_support`, and compares received device-buffer data.
This mode has **not** been compiled or run on the CPU-only development workspace.
It is a transport harness, not a working GPU solver, NVLink benchmark or proof of
GPUDirect. The query checks MPI capability, not the actual path taken; see the
[Open MPI capability API](https://docs.open-mpi.org/en/main/man-openmpi/man3/MPIX_Query_cuda_support.3.html).

The NVIDIA solver backend packs/unpacks on GPUs, synchronizes producers before
MPI and consumers after MPI, and initially uses one rank per GPU with pinned host
staging. A separate CUDA-aware transport may use NVLink/NVSwitch or
PCIe depending on hardware and transport; verify the path with profiling. Noncontiguous
GPU MPI datatypes can be expensive; see [Open MPI CUDA support](https://docs.open-mpi.org/en/main/tuning-apps/accelerators/cuda.html).
CPU MMG would require a collect/remesh/repartition/upload stage.

The CPU solver and `accelerator_halo.cxx` now share the exact SoA pack/unpack
indexing in `mpi_halo_kernels.hpp`. The harness uses OpenACC allocations and
synchronous device kernels, then passes contiguous device addresses to the same
MPI post/wait helpers. It covers irregular rows, scalar/vector/tensor widths,
asymmetric empty messages and reusable buffers in both dimensions. Explicit host
mode passed at 1/2/3/4/6/8 ranks. Host success verifies indexing and transport only.
The local workspace exposes no GPU. Remote 2D/3D halo checks passed on A100,
H100 and A6000 devices on 2026-11-09, verifying distinct physical UUIDs and larger
buffers. The tested H100 CUDA IPC configuration failed for larger messages;
CUDA-aware MPI with host staging passed. Full solver verification is now separate
from that harness; passing transport checks alone does not establish solver
correctness, performance or NVLink routing. See the
[accelerator test commands](../tests/functional/mpi/README.md#accelerator-halo-checks).

The existing OpenACC build targets NVIDIA GPUs. An Intel integrated GPU would
require a supported compiler/backend and kernel port, such as SYCL or Intel
OpenMP target. First query the actual device's double-precision support and
measure representative DES kernels; SYCL exposes `aspect::fp64` as shown in
[Intel's device discovery guide](https://www.intel.com/content/www/us/en/developer/articles/technical/device-discovery-with-sycl.html).
No Intel GPU runtime/device is exposed in this workspace. NPU AI TOPS describes
neural-network throughput and is not a measurement of DES double-precision mesh
operations; it does not justify an NPU port.

Full timestep and MMG-transition checks accompany the GPU backend. Compare
complete scientific runs, including MMG, initialization and output, and retain
CPU comparisons separately from one-GPU/multi-GPU comparisons. Passing halo or
short solver checks does not predict multi-GPU solver speedup.

## Partitioning and remeshing choices: research snapshot, 2026-11-05

Partitioning assigns ownership; remeshing changes vertices/connectivity. A graph
partitioner cannot replace MMG, and adding a parallel remesher does not by itself
parallelize DES physics, interpolation, marker management or output. The current
implementation still uses the in-tree RCB partitioner and centralized MMG.
METIS is now optional; ParMETIS, ParMmg2D and Omega_h are not production backends.

### Partitioners

- **METIS:** a practical first graph-partitioning candidate while the full mesh
  already resides on root. It supports edge-cut and communication-volume objectives
  and load constraints. DES communicates through complete node stars, so a graph
  containing only face-neighbor elements would omit some dependencies. Model the
  actual exchanges and compare resulting halo bytes and maximum rank time, rather
  than relying on the library's edge-cut score alone. See the
  [METIS manual](https://karypis.github.io/glaros/files/sw/metis/manual.pdf).
- **ParMETIS:** appropriate for already distributed graphs, with incremental
  `AdaptiveRepart` and multiple balance constraints. Initial partitioning and
  adaptive repartitioning solve different cost problems; migration must be priced
  against communication saved before the next remesh. See the
  [ParMETIS manual](https://karypis.github.io/glaros/files/sw/parmetis/manual.pdf).
- **Scotch/PT-Scotch:** credible graph alternatives; PT-Scotch operates on
  distributed graphs and Scotch also supports mapping onto weighted hardware
  topologies. This is relevant to unequal compute units or communication links,
  but does not establish DES performance. See the
  [project description](https://www.labri.fr/perso/pelegrin/scotch/).
- **Zoltan RCB/PHG:** RCB provides geometric, multiweight balancing; PHG supports
  hypergraph partitioning and migration-aware repartitioning. Shared-node groups
  may be represented more naturally as hyperedges than pairwise edges. That is a
  modeling hypothesis for DES, not a demonstrated speedup. PHG generally spends
  more partitioning time for potentially better communication quality. See
  [RCB](https://sandialabs.github.io/Zoltan/ug_html/ug_alg_rcb.html),
  [PHG](https://sandialabs.github.io/Zoltan/ug_html/ug_alg_phg.html), and the
  [graph/hypergraph tradeoff](https://sandialabs.github.io/Zoltan/ug_html/ug_graph_vs_hg.html).

MFEM's documented nonconforming `Rebalance()` uses a space-filling curve over
refinement trees; initial conforming partitioning uses METIS. It therefore does
not establish that METIS is always the best dynamic rebalancer. DES's general
MMG remeshing does not retain that refinement-tree structure. See
[MFEM's load-balancing explanation](https://mfem.org/howto/ncmesh/#load-balancing).

### Remeshers and frameworks

- **MMG2D/MMG3D:** the current DES wrapper calls these sequential libraries on
  root. Keeping them initially preserves the validated workflow, but the root
  memory requirement and remeshing pause remain.
- **ParMmg:** distributed CPU adaptation for 3D volume meshes, using iterative
  local MMG remeshing and repartitioning. The latest tagged release checked was
  v1.5.0 (2024-11-01); repository activity also exists after that release. An old
  release date alone does not prove abandonment. The release notes emphasize
  distributed memory capacity and explicitly do not promise time scalability.
  See [the project](https://github.com/MmgTools/ParMmg) and
  [v1.5.0 notes](https://github.com/MmgTools/ParMmg/releases/tag/v1.5.0).
- **ParMmg2D:** a separate newer project supplies distributed 2D adaptation,
  also built around MMG. It is community-supported and had no published releases
  when checked. Its build pins an MMG commit and accesses MMG internals, so verify
  version compatibility and library interfaces before adopting it. See the
  [README](https://github.com/MmgTools/parMmg2D/blob/main/README.md) and
  [build definition](https://github.com/MmgTools/parMmg2D/blob/main/CMakeLists.txt).
- **Omega_h:** triangle/tetrahedral adaptation with MPI and accelerator support.
  The actively maintained SCOREC fork has Kokkos backend integration, including
  SYCL configuration, making it a useful longer-term portability candidate.
  Porting DES mesh/state handling and validating field transfer would be substantial
  work; this is not a drop-in MMG replacement or proof of support/performance on
  the laptop iGPU. See the [SCOREC repository](https://github.com/SCOREC/omega_h)
  and [backend configuration](https://github.com/SCOREC/omega_h/blob/master/CMakeLists.txt).
- **PRAgMaTIc:** another 2D/3D anisotropic MPI/OpenMP mesh-adaptation library.
  **PETSc DMPlex** is a mesh-management framework exposing adaptation backends
  including MMG and ParMmg; adopting it is a larger architectural change, not an
  independent cure for sequential remeshing. See
  [PRAgMaTIc](https://meshadaptation.github.io/) and
  [DMPlex](https://petsc.org/main/manual/dmplex/).

In the measured 50,000-step mature core-complex parent continuation, MMG itself
took 0.102 s of 25.118 s (two remeshes, about 0.41%). That small CPU case does not
predict a large or GPU-accelerated case, but it does not justify a costly remesher
replacement for present laptop performance.

The isolated MMG-family prototypes are described in the
[test guide](../tests/functional/mpi/README.md#parallel-remeshing-compatibility-prototypes).
Capture/replay uses the exact inputs prepared by the existing DES MMG pipeline,
including Chase's metric/freezing/interface policies; none of these policies or
the production field/marker remapping were replaced. Serial MMG preserves the
captured required entities, while the tested parallel CLIs fail that contract
on multiple rank counts. ParMmg2D also loses boundary references in centralized
output and its native metric export can be invalid. Keep the serial-MMG backend;
parallel integration requires these compatibility issues to be resolved first.
The prototype retains an experimental metric/indexing patch and an explicit,
geometry-bounded 2D label reconstruction mode, neither of which resolves all
required-entity failures. No upstream issue/PR was opened. Omega_h remains a
separate research candidate if GPU-resident remeshing becomes a requirement.

Before changing the remesher, check required/frozen entities, free surfaces and
material interfaces, mesh validity/quality, stress/plastic-strain and marker
transfer, conservation and bounded scientific observables across repeated remeshes.
Parallel remeshing can produce a different valid mesh, so byte-identical output
against serial MMG is not a suitable universal acceptance criterion. Retain exact
checks for changes that preserve the same mesh/arithmetic, as the current CPU
halo implementation does.
