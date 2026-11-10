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

Build with an NVHPC-backed MPI wrapper and its compatible runtime. Submit
parallel work through the site's scheduler and keep model data on the storage
allocated for the job. Stage dependencies first on offline compute hosts.

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

GPU MPI builds use IEEE division; `nofma=1` disables host and device FMA.
Use the project's supported NVHPC SDK and target the installed GPU architecture.
Keep MPI/SDK runtimes separate when staging libraries: an HDF5 library directory
containing an older SDK can make a newer compiler link the old OpenACC runtime.

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

For a runtime comparison, use `benchmark_partitioners.py` and
`plot_partitioners.py` as described in the
[test guide](../tests/functional/mpi/README.md#rcb-versus-metis). Both methods use
the same MPI executable, one thread per rank, identical verified CPU masks,
warmups and randomized paired rounds. Every saved field/checkpoint is checked
against the parent. Reports distinguish total simulation time, complete launcher
walltime, partition rebuild cost and logical halo payload. A smaller halo alone
does not establish a faster run, and local shared-memory results do not predict
multi-server network performance.

RCB remains the default. Existing paired checks found no consistent end-to-end
advantage for METIS. A smaller halo alone does not justify switching partitioners;
measure partition rebuilds, initialization and timestep costs on the target model.

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

Use the [Grid Engine workflow](../tests/functional/mpi/README.md#cpu-correctness-through-grid-engine)
for scheduler-managed resources, or the
[SSH-host workflow](../tests/functional/mpi/README.md#cpu-correctness-across-ssh-hosts)
for explicitly allocated hosts with separate working directories.

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

### Verification scope

Native Apple Silicon checks with Apple clang, Open MPI, OpenMP, HDF5 and MMG
built both dimensions as serial and MPI executables. Against master `6195cc2`,
four parent/changed regressions and 34 MPI matrix comparisons were bit-exact.
All 20 PT checks passed without changing the existing tolerances; six ADR field
comparisons differed within `rtol=1e-8`, `atol=1e-12`. This verifies short CPU
cases on one Mac; GPU and GoSPL coupling were not exercised there.

GPU checks cover bounded fixed-mesh and remeshing cases with physical device
identity and kernel execution checks. Adaptive CPU/GPU trajectories can produce
different meshes; compare physical observables when topology or coordinates
change. These checks do not establish universal scientific equivalence or speedup.
See the [coupling scope](gospl-mpi-coupling.md#verification-scope) for GoSPL.

Keep per-run logs, hashes, reports and hardware notes under an ignored `*.local/`
directory. Version reusable commands, acceptance criteria and scoped conclusions.

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

## Communication overhead

`benchmark_overhead.py` compares two MPI executables with blocking and overlap
modes using the same CPU masks, one thread per rank, warmups and randomized
paired rounds. See the [timing commands](../tests/functional/mpi/README.md#mpi-overhead-measurements).
Overlap can help larger fixed meshes and cost time on small cases. Include root
collections and adaptive mesh management in end-to-end measurements. Rank timings
and halo payloads diagnose bottlenecks; neither alone establishes a speedup.

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
This transport harness does not establish full-solver correctness, NVLink routing
or GPUDirect. The query checks MPI capability, not the actual path taken; see the
[Open MPI capability API](https://docs.open-mpi.org/en/main/man-openmpi/man3/MPIX_Query_cuda_support.3.html).

The NVIDIA solver backend packs/unpacks on GPUs, synchronizes producers before
MPI and consumers after MPI, and initially uses one rank per GPU with pinned host
staging. A separate CUDA-aware transport may use NVLink/NVSwitch or
PCIe depending on hardware and transport; verify the path with profiling. Noncontiguous
GPU MPI datatypes can be expensive; see [Open MPI CUDA support](https://docs.open-mpi.org/en/main/tuning-apps/accelerators/cuda.html).
CPU MMG uses a collect/remesh/repartition/upload stage.

The CPU solver and `accelerator_halo.cxx` now share the exact SoA pack/unpack
indexing in `mpi_halo_kernels.hpp`. The harness uses OpenACC allocations and
synchronous device kernels, then passes contiguous device addresses to the same
MPI post/wait helpers. It covers irregular rows, scalar/vector/tensor widths,
asymmetric empty messages and reusable buffers in both dimensions. Explicit host
mode passed at 1/2/3/4/6/8 ranks. Host success verifies indexing and transport only.
Full solver verification is separate from halo checks; passing transport checks
alone does not establish solver correctness, performance or NVLink routing. See
the [accelerator test commands](../tests/functional/mpi/README.md#accelerator-halo-checks).

The existing OpenACC build targets NVIDIA GPUs. Other accelerators need a
supported backend and separate kernel/solver validation.

Full timestep and MMG-transition checks accompany the GPU backend. Compare
complete scientific runs, including MMG, initialization and output, and retain
CPU comparisons separately from one-GPU/multi-GPU comparisons. Passing halo or
short solver checks does not predict multi-GPU solver speedup.

## Remeshing limits

Partitioning assigns ownership; remeshing changes vertices/connectivity. The
production backend uses optional METIS or RCB partitioning and sequential
MMG2D/MMG3D on root. Parallel remeshers are not integrated. Root memory and the
remeshing pause therefore remain scaling limits.

A replacement remesher must preserve required/frozen entities, free surfaces and
material interfaces, and validate mesh quality, field/marker transfer,
conservation and physical observables across repeated remeshes. Different valid
meshes require physical comparisons rather than universal byte-exact equality.
