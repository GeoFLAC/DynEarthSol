# Intranet MPI, GPU and DES–GoSPL validation

Validation completed on 9 November 2026 against parent
`cce59ab5b2cfcf2520519ef2b818361eb0ae7e7c` and the working-tree MPI/GPU changes.
CPU MPI, GPU MPI functional checks and bounded DES–GoSPL coupling are verified
within the scopes below. GPU MPI remains experimental. CPU/GPU post-remeshing
arrays differ in the cases recorded below; pointwise array equality does not
establish or disprove physical equivalence once the meshes can differ.

Interpretation updated on 10 November 2026: these array-comparator rejections
are reproducibility observations, not demonstrated solver defects. Physical
equivalence across different meshes still requires suitable observables.
Earthsim scheduler availability is outside the remaining work at the user's request.

The subsequent [independent MPI coupling implementation and checks](gospl-mpi-coupling.md)
supersede the initial single-process GoSPL support scope below.

## Resources and execution

| Resource | Work completed | Execution |
|---|---|---|
| godo | CPU MPI across 4 and 8 allocated nodes | Grid Engine `qsub` |
| godo-s | CPU regression, 2 A100 GPU MPI, controlled scaling | Grid Engine `qsub` |
| nemex-gpu | 2 H100 GPU MPI; real CPU/GPU GoSPL coupling | Detached jobs on mounted `/data/slee91` |
| kersti | 3 A6000 GPU MPI; kernel traces | Detached jobs on mounted `/data/slee91` |
| earthsim-01 | Mount and scheduler inspection | qmaster port 6444 refused connections; no computation bypass |

The hosts have no internet access. Complete dependencies and a packed GoSPL
environment were staged from the connected workstation. Jobs use persistent
working directories and scheduler/detached execution independent of the client
SSH/VPN connection; result collection still requires connectivity. No other
users' or pre-existing jobs were stopped.
One H100 was shared with an unrelated simulation, so H100 timings are excluded
from performance claims. All two-H100 correctness checks still identified
distinct physical GPUs.

## CPU MPI

The godo-s regression job **854** passed **78 comparisons** across both
dimensions, 1/2/3/4/6/8 MPI ranks, hybrid threading, MMG, restart, optional physics
and averaged output. Parent `make set` and changed `make cmp` passed.

Godo jobs **8633/8634** passed 24 comparisons on four allocated nodes; jobs
**8636/8637** passed 12 on eight allocated nodes. These **36 multi-node
comparisons** cover both dimensions, overlap on/off, remeshing and each run's
own restart. Recorded allocation files establish actual execution hosts;
concurrent jobs can share nodes, so these are not 16 unique nodes.

After the last initialization and coupling fixes, job **877** rebuilt serial
and MPI executables in both dimensions and repeated independent parent
`make set` / four-rank `make cmp` checks. Both passed. Job **881** repeats these
checks after the final GPU boundary-reduction synchronization fix; both passed.
Eight scheduler-launcher unit tests and Python/shell syntax checks also passed.

## GPU MPI

The implementation uses one physical device per local MPI rank, compact GPU
timestep kernels, persistent device buffers and pinned host halo staging. MPI
receives host pointers. Initialization, MMG, remapping, markers and output
remain CPU management work. This verifies neither NVLink nor GPUDirect.

| Hardware | Primary short matrix | Result |
|---|---|---|
| 2 A100 | godo-s job 868 | 20 checks passed |
| 2 H100 | `hardware-TNUhnE` | 20 checks passed |
| 3 A6000 | `hardware-3163cI` | 30 checks passed |

The **70 checks** span 2D/3D fixed mesh, MMG, hydraulic and averaged-output
cases. They compare complete saved fields and checkpoints with finite floating
values at `1e-8`, exact integer/topology data and distinct GPU UUIDs.
Non-averaged MMG cases also resume their own checkpoints. Averaged trajectory
restart is excluded because the parent resets averaging windows on restart.
A6000 runtime traces record actual temperature, strain, stress, force and
velocity GPU kernels in both dimensions; passing output checks alone were not
used as evidence of offload.

NVHPC 21.3 builds target SM80 on A100/A6000; H100 requires the native NVHPC 25.7
SM90 build. IEEE division and disabled FMA are part of the tested configuration.
The old staged runtime could not enumerate H100 and was replaced for that host.

Final review also added a wait before CPU MPI consumes the asynchronous 2D
boundary-height reduction. After rebuilding both dimensions, eight-step fixed
mesh checks compare CPU/GPU and GPU rank counts again: A100 job **879**,
H100 `final-sync-sSFgtT` and A6000 `final-sync-FU3VEF`. All 12 checks passed,
covering both dimensions and one versus the host's maximum tested GPU rank
count. These targeted repeats are kept separate from the 70 primary matrix checks.

Longer checks separate scientific behavior from GPU rank-count consistency:

- A100 CPU versus two GPUs: 2D/3D Maxwell relaxation, 5,000 steps, passed exactly
  and matched the analytical tolerance `1e-6`. The 2,000-step RSF case passed
  all physical fields/checkpoints at `1e-8` and its complete monitor CSV exactly.
- A6000 one versus three GPUs: Maxwell 5,000 steps, core-complex MMG 20,000
  steps/two remeshes, and RSF 2,000 steps all passed bit for bit.
- H100 one versus two GPUs: core-complex MMG and RSF passed bit for bit.
  The original Maxwell run failed an initial-checkpoint diagnostic comparison
  because fixed-timestep startup left `max_global_vel_mag` indeterminate.
  Explicit diagnostic initialization fixed it. Bounded four-step Maxwell
  repeats in both dimensions passed every saved/checkpoint field bit for bit
  on H100 1/2 GPUs and A6000 1/3 GPUs. The original failing report remains;
  the full 5,000-step H100 Maxwell run was not repeated after this scalar fix.
- Standard repository tests on A100 one versus two GPUs passed bit for bit in
  2D/3D (job 871), including full 3D saved/checkpoint comparisons. The latter
  records one legitimate matching positive-infinity output-cadence metadata
  value; it does not exclude physical datasets or tolerate other nonfinite data.

**Post-remeshing array-comparison rejections:** A100 core-complex MMG first
exceeds the pointwise threshold after remeshing at frame 2, mesh-quality relative difference
`1.2818153157184855e-4`. The standard 3D CPU/GPU comparison (job 867, 400 steps,
two remeshes) also exceeds `1e-8` in coordinate components, up to `3.61e-6`.
Other fields in that standard comparison meet the tolerance. The raw comparator
reports are retained unchanged. Adaptive meshing can change coordinates,
connectivity and node/element identities after small arithmetic differences;
requiring equal indexed arrays across those trajectories is not a general
scientific acceptance criterion. These observations alone do not demonstrate
incorrect CPU or GPU physics, and do not automatically establish agreement.

Fixed-mesh and same-mesh arithmetic checks retain their existing exact/`1e-8`
requirements. Across different meshes, assess physical fields at common spatial
locations, volume-weighted stress/strain observables, topographic profiles and
valid mesh quality. Account for the model's deformation, thermal expansion and
erosion/sedimentation when assessing conservation. A remapping error would need
evidence such as an unexplained jump in those observables, invalid elements or
an error that does not decrease with refinement. That mesh-independent assessment
has not yet been performed for the two array-comparator rejections above.

## Measured overhead and scaling

Godo-s job **873** used eight physical CPU cores, a shared CPU initial checkpoint,
one warmup and three randomized measured rounds per variant. All 80 runs had
physical outputs checked. The table shows median simulation seconds for 1,000
steps, including restart initialization and output, excluding launcher startup.

| Execution | ~4,001 elements | ~64,000 elements |
|---|---:|---:|
| OpenMP 1 | 1.135 | 17.356 |
| OpenMP 2 | 0.689 | 10.660 |
| OpenMP 4 | 0.479 | 7.369 |
| OpenMP 8 | 0.371 | 5.252 |
| CPU MPI 1 | 1.102 | 17.993 |
| CPU MPI 2 | 0.790 | 11.853 |
| CPU MPI 4 | 0.655 | 8.532 |
| CPU MPI 8 | 0.633 | 6.618 |
| A100 MPI 1 | 2.140 | 7.208 |
| A100 MPI 2 | 3.789 | 8.999 |

That allocation bound CPUs to NUMA node 0, remote from the GPUs. Rather than
repeat all variants, job **874** tested only CPU1/GPU1/GPU2 with CPUs local to
the two GPUs (NUMA nodes 1/2). All 24 warmup/measured runs passed physical checks.

| Execution | Small simulation / full wall seconds | Large simulation / full wall seconds |
|---|---:|---:|
| OpenMP 1 | 1.110 / 1.228 | 17.189 / 17.290 |
| A100 MPI 1 | 1.975 / 8.603 | 6.834 / 13.726 |
| A100 MPI 2 | 3.761 / 10.863 | 8.482 / 15.682 |

The large one-GPU case is **2.515× CPU1 during simulation**, but only **1.260×
including launcher/device startup**. GPU startup contributes about seven
seconds in these invocations. Two GPUs are slower than one at both measured
sizes; their reported communication interval grows to about 2.36 seconds in
the large local case versus 0.85 for one GPU. These timers include synchronization
and overlap operation timers, so they must not be summed as independent costs.
The tested eight-thread OpenMP execution is faster than the tested one-GPU
execution. Different affinities and these short 2D workloads limit extrapolation.

For these sizes, use OpenMP or one GPU according to measured total runtime.
There is no demonstrated multi-GPU speedup or H100/A6000 scaling result.
Timings retain exact executable hashes. They precede the final startup-diagnostic
initialization fix; their adaptive-timestep cases initialize those diagnostics
and already passed full physical comparisons. The hot kernels were unchanged.

## DES–GoSPL coupling boundary

The scope is coupling only. The installed GoSPL 2026.7.14/PETSc 3.21.2 CPU
environment and real C++/Python bridge were reused; no GoSPL algorithm validation
or standalone scaling is claimed. PETSc reports CUDA support disabled.

`coupling-4y_zaw90` ran eight DES steps and four actual exchanges on a coarse
3D tutorial-derived uniform box: CPU with two OpenMP threads, then one H100
using the CPU initial checkpoint and the same GoSPL mesh. Marine deposition
was disabled to keep the check focused on the boundary. Each exchange covered
99 surface points. Audit records verified:

- one-time elevation seeding and one embedded GoSPL MPI rank;
- finite surface arrays, first instantaneous and subsequent averaged velocities;
- matching DES/GoSPL clocks at 38.464494, 76.928989, 115.393483 and 153.857978 yr;
- nonzero returned erosion, applied exactly once with zero measured height error;
- successful destruction with no remaining bridge model handles;
- all five saved frames and four continuation checkpoints bit for bit CPU/GPU.

Two coupling fixes were verified: use DES `YEAR2SEC` (365.2422 days) for elapsed
time and velocity conversion, and finish pending GPU work before host exchange.
The earlier hard-coded 365-day conversion introduced about 0.0664% disagreement.

The four real GoSPL step calls together took **0.125 s CPU / 0.154 s GPU-run**;
velocity transfer calls together took about **0.00128 s** in each run. Audit
writes took **0.00176 / 0.00275 s**. Full process wall times were **4.622 / 4.722 s**;
these tiny instrumented runs establish functionality, not coupling performance
or GPU acceleration. Actual model clocks are audited directly; GoSPL file-output
cadence is not used as the clock reference.

**Initial November 9 scope:** one DES process with OpenMP or one GPU, one GoSPL
rank. At that point `usempi=1 usegospl=1` was rejected by the build. This was an
implementation gap. It was subsequently removed by the independent-worker
implementation described below.

### MPI extension follow-up, 10 November

GoSPL's existing MPI capability is now connected through an independent
persistent worker world. DES and GoSPL rank counts are separate, including
DES 3 / GoSPL 2 and GPU DES 2 / GoSPL 2. Coupling transfers use owned PETSc IDs
and retain model state through DES remeshing. Eight bounded CPU/GPU coupling
cases, 36 full physical comparisons, distributed seed placement, time/step
gating and controlled failure cleanup are recorded in the
[implementation and verification guide](gospl-mpi-coupling.md).

## Evidence and reproducibility

Repository instructions and drivers are in [MPI tests](../tests/functional/mpi/README.md)
and [coupling tests](../tests/functional/gospl/README.md). Machine-specific
artifacts are intentionally ignored under:

- `benchmarks-cores/mpi-gpu.local/`: source hashes, executable hashes, allocation
  files, scheduler accounting, JSON comparisons, kernel logs and scaling samples.
  `progress.json` records final disposition; `validation-summary.json` indexes
  the primary reports. `source-manifest.json` identifies the final local sources.
- `benchmarks-cores/gospl-functional.local/results/coupling-4y_zaw90/`: complete
  short coupling fixture, saved fields/checkpoints, JSONL call audit and NPZ arrays.
- Remote GPU results: `/home/slee91/des-mpi-gpu-20261109` on godo-s;
  `/data/slee91/des-mpi-gpu-20261109` on nemex-gpu/kersti.
- Remote coupling environment/results:
  `/data/slee91/des-gospl-functional-20261109` on nemex-gpu.

Failed exploratory/build runs remain available but are not counted as successful
verification. Standalone GoSPL probes are excluded from coupling counts.
No commit, pull request, release or deployment was created.
