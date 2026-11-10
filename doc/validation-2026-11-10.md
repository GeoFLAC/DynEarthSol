# DES–GoSPL coupling completion, 10 November 2026

The four follow-up items are complete for the bounded tutorial-derived coupling
case: nonzero continuation, actual two-node coupling, mesh-independent physical
comparisons and one controlled resource/cost comparison. GoSPL's existing MPI
solvers were used unchanged. This is coupling verification, not independent
validation of GoSPL or a general scaling survey.

## Continuation

DES checkpoints now reference an immutable GoSPL state generation. Rank-local
arrays, PETSc vectors/matrices and solver settings, forcing/counter data and the
original surface mesh are retained. DES also saves the pending coupling clock,
velocity-reference coordinates and remesh generation. Restart restores these
without reseeding or replaying earlier surface steps. Checksums, configuration,
PETSc version, mesh and partition compatibility are checked before restoration.

Every case ran 12 DES steps with an exchange every three steps, or a 50-year
time threshold. Checkpoints at steps 1, 5 and 7 are between mesh checks; they
exercise the initial pending interval and later partial intervals.

| Execution | DES ranks, continuous → restart | GoSPL ranks | Restart steps | Complete save/checkpoint comparisons |
|---|---:|---:|---|---:|
| CPU, fixed mesh | 3 → 3 | 2 | 1, 5, 7 | 46 |
| A6000, fixed mesh | 3 → 3 | 2 | 1, 5 | 36 |
| CPU, six MMG remeshes | 3 → 2 | 2 | 1, 5, 7 | 46 |
| CPU, time triggering | 3 → 3 | 2 | 1, 5 | 36 |
| A6000, six MMG remeshes | 3 → 2 | 2 | 5 | 14 |

All **178 file comparisons** passed the normal finite, `1e-8` relative-field
criterion; integer fields/topology were exact. The checked final owned/local
GoSPL physical vectors were bit-exact in every continuation, and GoSPL time plus
the pending interval agreed with DES time. Some DES floating fields differed
below tolerance, so continuation is not claimed to be universally bit-exact.

Private copies with a damaged rank-state checksum and a changed GoSPL rank
count were both rejected for the expected reason. Each worker world exited 53
after cleanup; original fixtures were unchanged. The first development attempt
failed because the C++ input vector was not allocated before reading; that
failed artifact remains separate from the accepted runs.

GoSPL rank count/partition stays fixed across continuation; DES rank count may
change. The tested processes are continental erosion/diffusion with marine
deposition disabled. Other external GoSPL configurations need a targeted
continuation check. Ordinary DES output precedes same-step remeshing, so select
restart checkpoints between mesh checks.

## Actual two-node coupling

Godo `qsub` job **8655**, `ompi_pf` four-slot allocation, used two slots each on
`godo_node01` and `godo_node02`. GoSPL used one real distributed rank on each
node. DES ran first with one rank, then one rank per node. Per-rank host records,
allocation files and GoSPL collective host records establish placement.

Both runs completed eight steps/four real exchanges. All **nine** complete
save/checkpoint comparisons between DES1 and DES2 were bit-exact. Audits checked
99 surface points at each exchange, one initialization/seed, correct ownership,
clock/velocity units and exactly-once height application. Seed error was at most
`1.39e-17 m`; height application error was zero. Verification job **8657** reused
the finished models after a missing postprocessing module was staged; modeling
was not repeated for that packaging error.

The old cluster required three environment adaptations: rebuilding the same
GoSPL v2026.7.14 Fortran extension against glibc 2.17 (job **892**), a job-local
hostname adapter for Open MPI 4's handling of underscores, and linking DES to
its own MPI libraries before the conda environment. System hostnames, SSH
credentials and firewall settings were unchanged. Independent MPI worlds retain
separate library environments. The earlier Nemex–Kersti bootstrap reached the
remote daemon but its non-SSH TCP path was blocked; native scheduler networking
provided the accepted route.

These tiny multi-node cases establish function, not speedup. In particular the
DES2 run's four worker erosion calls took about 19.25 s versus 0.30 s in DES1;
this allocation/transport result argues against spreading a small coupled model
across nodes for performance.

## Physical comparison after remeshing

The existing 400-step 3D CPU/GPU outputs from job **867**, after two MMG
remeshes, were reused. Each had 3,978 nodes/12,500 tetrahedra at the same physical
time, 945.93327 years. Different mesh coordinates do not constitute a coupling
failure. The new comparison does not require matching node or element IDs.

Both meshes had finite physical fields and consistently oriented, nondegenerate
tetrahedra. All 2,048 common interior points and 529 common top-surface points
were covered. Nodal quantities were interpolated barycentrically; cell fields
were sampled in their containing tetrahedra. Volume integrals use tetrahedron
weights, with exact linear-field temperature integrals.

| Observable | CPU/GPU difference |
|---|---:|
| Total volume | 0 at stored precision |
| Volume-weighted mean stress, relative maximum | `9.01e-15` |
| Volume-weighted mean strain, relative maximum | `2.55e-14` |
| Volume-weighted mean strain rate, relative maximum | `3.00e-11` |
| Common-point velocity, relative RMS | `3.20e-11` |
| Common-point stress, relative RMS | `1.26e-13` |
| Common-point strain rate, relative RMS | `3.90e-10` |
| Top elevation, maximum absolute | `7.06e-10 m` |

These observations support agreement for this retained case. They do not impose
a universal remeshing tolerance or prove behavior for unrelated geometries. The
sampling utility assumes a box interior and a single-valued top surface.

## Representative coupled cost

Kersti ran one 3D tutorial-derived problem: 47,385 DES nodes, **204,800
tetrahedra**, 500 GoSPL vertices, 100 steps/five exchanges, fixed mesh. All variants
restarted the same initial checkpoint and used GoSPL1. DES and the worker were
bound to distinct physical CPUs; BLAS used one thread. Each variant had one
warmup, followed by three interleaved measured runs. The workload was run once
as this controlled comparison, without a hardware matrix or full tutorial sweep.

| DES execution | Full launcher wall, median (s) | Initialize (s) | Compute (s) | Output (s) |
|---|---:|---:|---:|---:|
| CPU, 1 thread | 54.676 | 32.161 | 14.815 | 4.516 |
| CPU, 4 threads | **25.995** | 13.841 | 4.538 | 4.486 |
| CPU, 3 MPI ranks | 48.599 | 34.383 | 6.259 | 4.890 |
| A6000, 1 GPU | 51.062 | 37.351 | 4.902 | 4.623 |
| A6000, 3 GPUs | 53.411 | 39.907 | **4.304** | 5.348 |

Columns are separate medians and need not sum exactly. DES's `Execute` timer
includes initialization/compute/output; the table's wall timer also includes
launch and worker shutdown. Overlapping RPC, kernel and MPI timers are not added.
Wall-time ranges across repeats were 25.994–26.043 s for CPU4 and
53.371–53.467 s for GPU3. The OMP4 variant's worker inherited `OMP_NUM_THREADS=4`
while bound to one physical CPU; the other variants inherited 1. All used
`OPENBLAS_NUM_THREADS=1`; worker erosion-call times remained 0.27–0.29 s.

All 15 measured runs passed **45** complete save/checkpoint comparisons;
the four other warmup variants passed another 12. All were bit-exact relative
to warmup CPU1. The CPU1 warmup's three self-comparisons are excluded from these
counts. Only five `run_and_get_erosion` calls occurred per 100-step run. The worker
persisted, and ordinary idle steps did not reseed or invoke the surface solver.
Five real worker erosion calls totaled median 0.27–0.29 s, about 1% of CPU4 wall
time. MPI/GPU `Surface` timers, which also include root work and transfers, were
about 0.45–0.53 s per run. These are inclusive observations, not a pure transport
overhead estimate.

For this workload the recommended default is **one DES rank with four OpenMP
threads plus one GoSPL rank**, reserving an additional CPU and explicitly setting
the worker's OpenMP/BLAS threads to one. It was 2.10× faster end to end than CPU1.
Three GPUs slightly improved compute time but took 2.05× CPU4's full wall time.
Use more GPUs/nodes only after demonstrating an end-to-end benefit on a larger
or different workload. The tested GoSPL/PETSc build uses CPUs; these GPUs
accelerate DES mechanics.

## Regression and evidence

Final CPU/GPU MPI GoSPL builds compiled both dimensions on Nemex (build
`mpi-build-PBxBJn`) and Godo-s (CPU job **897**, GPU job **900**). Final coupled
GPU continuation used the correctly isolated SDK MPI libraries. Parent `make
set` / changed four-rank `make cmp` passed in both dimensions in Godo-s job
**896** after checkpoint I/O changes. Non-HDF5 binary I/O compiled in both
dimensions. Python modules compiled, and the mesh/timing utilities processed
the retained real results. Earlier default GPU and independent-rank checks
remain in the [coupling guide](gospl-mpi-coupling.md).

Artifacts are under `benchmarks-cores/gospl-functional.local/`, including
`completion-results.json` and `completion-manifest.json`. Full short-run arrays
are retained in `results/restart-cpu-02`, `restart-gpu-01`, `continuation-remesh`,
`continuation-time`, `continuation-gpu-remesh`, `checkpoint-rejection` and
`nodes-8655`. Mesh-independent inputs/report are in `mesh-observables/`.
Cost logs/reports are local in `results/coupled-cost`; its full model arrays remain
on Kersti at `/data/slee91/des-gospl-functional-20261109/results/coupled-cost-pa6tj027`.

Exact build/executable hashes are retained per generation. The cost/fixed-mesh
continuation binaries precede the final surface-reference length-guard change;
the later remesh/time continuation and final builds include it. That guard only
changes validation of malformed checkpoint lengths. Historical November 9 GPU
results and their source hashes were not relabeled as current results.

All `/data` work used the mounted storage. Godo/Godo-s parallel work used `qsub`
and allocated slots. Earthsim scheduling was excluded as requested. Our jobs,
workers and bootstrap relay have finished; unrelated user jobs were left running.

## PR integration and runnable reviewer examples

Before publication, the GPU/coupling extension was integrated with the existing
CPU MPI history and current `master` (`68cc2a0`), preserving upstream C++17/CLI11
and physical-step PT support. The earlier matrices above retain their original
source/executable provenance; the following are the refreshed current-source checks:

- Godo-s job **901**: current-master parent `make set` / changed four-rank
  `make cmp` passed in both dimensions; serial and MPI builds passed.
- Godo-s job **905**: the quickstart CPU command passed all **34** cases
  bit-exact, including PT/hydraulic, MMG, restart and six-rank empty partitions.
  Four C++20 host halo checks passed in 2D/3D at one/two ranks. The driver needs
  six allocated CPU slots even when its main matrix selects ranks 1 and 2.
- Nemex: both CPU/GPU dimensions built with the current sources. CPU remeshing
  continuation (46 file comparisons) and time-triggered continuation (36) passed.
  The new checked-in `smoke.cfg`/`smoke.yml` example started from fresh meshes,
  completed DES3 + GoSPL2, and supplied a successful 46-comparison remeshing restart.
- Kersti, NVHPC **25.7**, SM86: all four current-source fixed-mesh GPU cases
  (2D/3D, GPU1/GPU2) passed **40** recorded complete file comparisons. Physical
  UUIDs and device kernel execution were verified. A further GPU3 → GPU2 coupled
  continuation with six MMG remeshes passed 14 complete file comparisons.
- Godo job **8664**: the refreshed two-node DES1/DES2 + GoSPL2 models passed
  all nine complete file comparisons bit-exact and the coupling-boundary audits.

The supported-SDK GPU repeat is the acceptance result for the C++17 integration.
A legacy NVHPC 21.3/GNU11 configuration compiled but failed CUDA symbol
initialization on Kersti; it is retained separately, not counted as a pass.
Performance timings above were not rerun or relabeled after the upstream merge.

The [reviewer quickstart](mpi-quickstart.md) includes prerequisites, copyable build
and launch commands, slot/GPU counts, expected outputs and continuation steps.
Refreshed artifacts are under `results/final-regression-901`, `reviewer-cpu-905`,
`reviewer-gpu25`, `reviewer-continuation`, `review-gpu-continuation` and
`nodes-8664` in the local evidence directory. Large GPU arrays remain on mounted
Kersti storage in `results/gpu-mpi-k73skkzj`; the PR contains source, fixtures,
reproduction tools and summaries.

## Automated review follow-up

The MPI/C++20 CI job exposed a startup crash with `usemmg=0`: the root's two
MMG-only mesh-size hint vectors are empty, but compact selection indexed them
with node IDs. Those fields now participate in MPI state only under `USEMMG`.
The failing three-rank case was reproduced before the fix and completes after it;
the MMG-enabled field list is unchanged.

Non-MMG 2D/3D builds passed parent (`c81fb3b`) `make set` / changed three-rank
`make cmp`, both bit-exact, with remeshing disabled in the retained tiny configs.
All 40 existing PT checks passed across dimensions 2/3, methods 0/1 and ranks
3/4, covering elasticity, Maxwell relaxation, hydrostatic/dynamic steps and
collective failure. Existing tolerances were retained. C++17/C++20 host halo
checks passed at ranks 1/2 in both dimensions. The Python review cleanup also
passed two-rank protocol failure, rank-local KeyboardInterrupt/SystemExit cleanup
and 14 coupled remeshing-continuation comparisons. Review artifacts are retained
under `results/pr-review-*`, `results/pr128-*` and the mounted Nemex run directories.
