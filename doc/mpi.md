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
the vendored CLI11 header, so Boost is not needed. Existing cfg files
apply without new simulation parameters. One CPU rank uses the ordinary
full-range kernels. A one-rank GPU run uses the same compact
backend as multi-GPU runs. NVHPC profiling cannot currently be combined with
`usempi=1`. GoSPL coupling is not yet available with `usempi=1`; the Makefile
rejects `usempi=1 usegospl=1`. Serial/OpenMP and non-MPI OpenACC executables
remain separate.

## NVIDIA GPU MPI

Build with an NVHPC-backed MPI wrapper and its compatible runtime. Submit
parallel work through the site's scheduler and keep model data on the storage
allocated for the job. Stage dependencies first on offline compute hosts.

```sh
export OMPI_CXX=/path/to/nvc++
make -j8 ndims=2 usempi=1 openacc=1 openmp=0 MPI_CXX=/path/to/mpicxx usemmg=1 hdf5=1
make -j8 ndims=3 usempi=1 openacc=1 openmp=0 MPI_CXX=/path/to/mpicxx usemmg=1 hdf5=1
```

The Makefile detects the GPU's compute capability; set `GPU_CC` to build for
another device.

Use one rank per physical GPU. Local rank N on a host runs on that host's visible
device N, so every rank on a host must see the same GPUs, at least one per rank.
A launcher that gives each rank a different single GPU stops local rank 1 with
exit code 52. Device selection precedes managed allocations; local ranks must
have distinct physical UUIDs. The implementation
packs compact SoA halo rows on devices, copies contiguous buffers through pinned
host memory, and supplies **host pointers** to MPI. GPU producers finish before
packing; receives finish before upload/unpack and buffer reuse. This path needs
ordinary MPI transport, and makes no claim about NVLink, CUDA IPC or GPUDirect.
Mesh states sent after remeshing, gathered and scattered field rows and array
reductions also go through pinned host buffers; scalars go from stack variables.
No MPI data buffer is managed memory: a CUDA-aware transport (UCX) corrupted
managed buffers on two GPUs.

Compact timestep fields use managed allocation and are prefetched after
distribution. Initialization, MMG, double-precision nearest-neighbor remapping,
global marker/surface management and output run on the root CPU. Conditional
OpenACC regions enforce this boundary; command boundaries drain pending kernels
before host state changes. MMG is followed by collection, repartitioning and
managed-state refresh. Managed allocation also covers the root state container.

GPU MPI builds use IEEE division. As in other builds, `nofma=1` (opt-in) disables
host and device FMA, which brings GPU results closer to CPU ones.
Use the project's supported NVHPC SDK and target the installed GPU architecture.
Keep MPI/SDK runtimes separate when staging libraries: an HDF5 library directory
containing an older SDK can make a newer compiler link the old OpenACC runtime.

The validation drivers retain full-field and checkpoint comparisons at `1e-8`,
integer/topology equality and GPU identities, and require each rank's GPU to run
the force update in its OpenACC launch log. Short 2D/3D
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

### METIS partitioning

Every MPI build partitions the mesh with METIS, so `usempi=1` needs METIS 5:

```sh
make ndims=2 usempi=1 usemmg=1 hdf5=1
OMP_NUM_THREADS=1 mpiexec -n 4 ./dynearthsol2d.mpi model.cfg
```

`METIS_INCLUDE_DIR` and `METIS_LIB_DIR` override compiler/linker default paths;
`make config` prints them. Repeat the build with `ndims=3` for 3D. One-rank CPU
execution uses the original kernels without partitioning. A one-rank GPU run uses
the compact backend with every element on its rank, without calling METIS. Tiny
meshes with no more elements than ranks use an explicitly logged recursive
coordinate bisection, preserving empty-partition support.

The METIS graph connects elements sharing any vertex, matching DES node-star
dependencies more closely than a face-only graph. K-way partitioning uses the
communication-volume objective, seed 0 and a 3% element-balance tolerance. Node
ownership and halo construction then use the same code as for tiny meshes. The graph objective
is a proxy: actual halo bytes also depend on node ownership and exchanged fields.
METIS does not currently use measured constitutive/marker weights or migrate
incrementally. It is rebuilt after each MMG remesh.

A smaller halo alone does not establish a faster run; include partition rebuilds,
initialization, output and timestep costs when assessing the target workload.

### Geometric partitioning of tiny meshes

For a mesh with no more elements than ranks, recursive coordinate bisection splits
element centroids along the longest spatial axis. Element counts differ by at most one, including non-power-of-two process
counts and empty partitions. Interior nodes belong to the incident element rank;
shared nodes are assigned to an incident rank with the smallest current node load.
This balances element counts first. It does not guarantee equal node counts,
marker counts, constitutive work or execution time.

Every compute rank stores owned nodes/elements first, followed by ghost rows.
The local mesh includes owned elements, the complete element stars of owned nodes,
and their vertices. Local indices map to stable global IDs until remeshing.
Node patch sums retain global element ordering to preserve the serial arithmetic.
Workers do not receive the full mesh. Rank zero also retains a global management
copy, so root memory and centralized operations remain scaling limits. Root sends
each rank's compact state as one MPI message, which must stay under 2 GiB
(`INT_MAX` bytes); a larger partition stops the run with exit code 52 before the
transfer. Use more ranks to make each partition smaller.

Ownership and halo maps are rebuilt after **every** MMG remesh, even if node and
element counts are unchanged. Compute ranks retain per-element material counts;
individual marker histories, barycentric coordinates and index lists stay on
the management mesh, where marker consistency is checked after rebalancing.
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
Rank zero runs this work with its OpenMP threads, so models with phase
changes, surface option 1 or frequent remeshing usually run faster with fewer
ranks per node and several threads each than with one single-threaded rank per core.
A global state is collected only when needed by these operations, checkpoints,
output or optional earthquake moment accumulation. Surface option 1
collects/scatters a top/boundary strip. Option 0 uses compact mechanical surface
correction outside mesh checks and requested surface reports, with neighbor
halos. Quality checks collect coordinates and
element volumes; advanced surface work can require full state. The optional
serpentinization grid collects state only when that model feature is active.
Monitor output collects queried rows and their dependencies.
Output averaging keeps local accumulators and gathers them when writing.
Earthquake speed/plastic-rate diagnostics use small global maximum reductions.

Regular dependency exchanges pack only needed ghost fields into reusable
contiguous peer buffers. Receives/sends use `MPI_Irecv`/`MPI_Isend`, followed by
`MPI_Waitall` and unpacking. Related fields share transfers; buffers and maps
persist until repartitioning. Temperature, force and mass assembly first posts
element transfers, computes owned nodes whose incident elements are all local,
waits/unpacks, then computes interface nodes. Per-node element summation order is
unchanged. Other exchanges still complete before dependent work.
MPI progress during computation depends on the MPI implementation; posting a
transfer does not by itself prove simultaneous data movement. Regular
timesteps do not gather full fields with `MPI_Allgatherv`.

## Adding a distributed kernel

Rank 0 runs the ordinary time loop on its global mesh. When the loop calls a
kernel that has an MPI operation, rank 0 broadcasts a `Command` naming the
operation, and every rank, rank 0 included, runs that kernel on its compact
mesh. The workers wait for commands in `MpiSimulation::worker_loop()`; both
sides run `Impl::execute()` in `mpi_simulation.cxx`.

**1. The dispatch line.** The kernel starts with

```cpp
if (var.mpi && var.mpi->dispatch(MpiOp::X, var)) return;
```

On rank 0's global variables, `dispatch()` broadcasts the command, runs the
operation and returns `true`, so the global body is skipped. Inside
`execute()` the same kernel is called again on the compact variables. That call
falls through to the body, because `dispatch()` returns `false` while an
operation is running, for variables other than rank 0's global ones, and in runs
that are not distributed (CPU runs with one rank).

**2. The operation.** Add `X` to `enum class MpiOp` and its name to `op_name()`,
and add a `case MpiOp::X:` to `execute()`. The case calls the kernel on the
compact variables and then exchanges every field the kernel writes that a
neighbor reads, for example

```cpp
case MpiOp::Temperature:
    update_temperature(param, var, *var.temperature, *var.tmp_result);
    exchange(MpiEntity::Node, {{var.temperature->data(), 1}});
    break;
```

Element loops run over owned elements (`mpi_element_end(var)`), node loops over
owned nodes (`mpi_node_end(var)`); the exchange fills the ghosts. Node assembly
that can overlap communication uses `mpi_gather_passes()` and
`mpi_gather_nodes()`, as `update_temperature()` does.

**3. Scalars.** `Command` carries the step, time, `dt`, the velocity and
timestep bounds, the surface velocity, the operation's argument and two
switches; `Impl::command()` copies them into each rank's variables. A new scalar
that the kernel reads goes into both. A global scalar result comes from a
reduction (`minimum()`, `sum()`, `maximum()`). If rank 0's time loop reads it
afterwards, copy it to `master` in the case (`if (rank == 0) master->... = ...`),
as `TimeStep` and `Residual` do.

**4. Work that only rank 0 can do.** Operations that need the global mesh,
markers or MMG are listed in `management()`. Their case first brings rank 0's
copy up to date (`full()` for every field, `geometry()` for coordinates and
volumes, `boundary()` for the surface strip), runs the function on `*master`
under `if (rank == 0)`, and then sends the result back with `distribute(false)`,
`distribute(true)` when the mesh changed (this also repartitions), or
`scatter_surface()`. OpenACC offload is off during these operations.

**5. The version rule.** After every operation, `execute()` advances `version`
unless the operation is on the list at its end (`Collect`, `Geometry`,
`Boundary`, `Sample`, `Average`, `Metrics`, `TimeStep`, `Residual`, `Report`).
`full()`, `geometry()` and `boundary()` gather again only when the version has
moved since their last gather. An operation that changes any field must
therefore not be on the list; one that only reads fields or computes scalars
may be added to it.

**6. Fields and offload.** A new per-node or per-element field of `Variables`
goes into `mpi_local_fields.hpp`, or compact meshes never get it.
`tests/functional/mpi/check_mpi_fields.py` fails for a field allocated with
`new` that is in neither that list nor its own list of fields that never leave
rank 0. An OpenACC loop that can also run on rank 0's management mesh needs
`if(des_acc_offload_enabled())`; in a GPU MPI build that mesh stays on the host.

**Options.** An MPI build accepts only the checked values of the options that
choose code inside the distributed kernels (`check_mpi_support()` in
`input.cxx`). Add a new value there once the steps above and a test cover it.
Kernel arguments other than `var` and the one scalar are not forwarded to the
compute ranks: a kernel that takes more must carry them in `Command`.

**7. Tests.** `tests/functional/mpi/run_mpi.py` compares serial and MPI runs
field by field; add a case to its list that uses the kernel. PT kernels are
checked by `run_pt_mpi.py`, halo packing by `run_accelerator_checks.py`.

## Reading the MPI log

`[mpi][partition]` log rows report step, owned nodes/elements/markers, ghosts and
neighbors for every rank at each partition generation. Counts show balance and
ghost storage; they do not establish equal runtime.

`[mpi][kernel]` reports operation counts and root elapsed seconds including
synchronization. Management operations have counters too. `[mpi][communication]`
reports timed communication, `halo_sent_bytes` summed over ranks,
`root_collection_bytes` and `state_bytes` sent to workers. These logical payload
counts exclude headers, reductions, MPI protocol traffic and local root copies;
they are not a complete network profile.

`[mpi][rank-profile]` reports operation time and instrumented communication time
on every rank. Their difference includes local work, scheduling and uninstrumented
overhead; it is not a hardware compute counter. These counters help diagnose
imbalance but do not establish asynchronous MPI progress.

## Correctness checks

On scheduler-managed resources, run the checks inside an allocation and pass the
site's MPI launcher with `--launcher`. Use HDF5/MMG executables built from the
parent commit and changed sources:

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

Repeat with a valid 3D case.

### Verification scope

The CPU checks above also build and pass natively on macOS (Apple clang, Open
MPI, OpenMP, HDF5 and MMG); they cover short CPU cases, not GPU runs.

GPU checks cover bounded fixed-mesh and remeshing cases with physical device
identity and kernel execution checks. Adaptive CPU/GPU trajectories can produce
different meshes; compare physical observables when topology or coordinates
change. These checks do not establish universal scientific equivalence or speedup.

Keep per-run logs, hashes, reports and hardware notes under an ignored `*.local/`
directory. Version reusable commands, acceptance criteria and scoped conclusions.

## Performance assessment

Compare equal CPU budgets and retain the same geometry, initial state and physical
outputs. Measure complete wall time including initialization, remeshing, output
and MPI startup. Use warmups and repeated runs without competing jobs; record CPU
placement and GPU identities. Each run records its ranks, nodes, ranks per node,
threads per rank, MPI thread level and MPI library in the
`[runtime.mpi]` section of `<modelname>.manifest` and prints them at start as
`[runtime][mpi]`. Local shared-memory measurements do not predict
multi-node or adaptive-model speedup. Keep workload-specific measurement
scripts with their results under ignored `*.local/` directories.

The [halo checks](../tests/functional/mpi/README.md#accelerator-halo-checks) exercise
the shared SoA pack/unpack indexing and MPI post/wait helpers, including irregular
rows and empty messages. They do not establish full-solver correctness or NVLink
routing. Use the full solver CPU/GPU checks for timestep and MMG transitions.

## Remeshing limits

Partitioning assigns ownership; remeshing changes vertices/connectivity. The
production backend uses METIS partitioning and sequential
MMG2D/MMG3D on root. Parallel remeshers are not integrated. Root memory and the
remeshing pause therefore remain scaling limits.

A replacement remesher must preserve required/frozen entities, free surfaces and
material interfaces, and validate mesh quality, field/marker transfer,
conservation and physical observables across repeated remeshes. Different valid
meshes require physical comparisons rather than universal byte-exact equality.
