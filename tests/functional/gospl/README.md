# Bounded DES–GoSPL coupling check

This check exercises the real embedded bridge and installed GoSPL. It verifies
the DES coupling boundary, not GoSPL's surface algorithms or parallel scaling.
Run it with a working Python/GoSPL environment and separately built CPU/OpenMP
and single-GPU DES executables (`usegospl=1 usempi=0`).

```sh
python tests/functional/gospl/run_coupling_smoke.py \
  --cpu /path/to/dynearthsol3d.cpu.gospl \
  --gpu /path/to/dynearthsol3d.gpu.gospl \
  --extension /path/to/gospl_extensions --workdir /mounted/results
```

The driver coarsens the tutorial to a uniform 3D box, runs eight DES steps and
four coupling events, and disables marine deposition. The CPU run uses two
OpenMP threads. The GPU run reuses its initial checkpoint and GoSPL mesh to
avoid a duplicate mesh generation and control initialization differences.
Each embedded GoSPL model has one MPI rank. The audit module delegates to the
real bridge and records arrays, clocks, operation durations and cleanup.

Checks cover one-time elevation seeding, finite surface arrays, instantaneous
and interval-averaged velocity transfer, consistent year units, exactly one
application of each returned erosion increment, and model destruction. All
saved physical fields and continuation checkpoints are compared at `1e-8` with
exact topology. JSON reports retain executable hashes, timings and overrides.
Timings include instrumentation and startup; this tiny check is not a GPU
performance benchmark. GoSPL output-file cadence is independent of the actual
coupling clock, which is audited directly.

## Independent MPI coupling

`run_mpi_coupling.py` uses the persistent worker launcher and real GoSPL at
independent DES/GoSPL rank counts. It requires CPU/GPU `usempi=1 usegospl=1`
executables and a retained initial CPU fixture from the check above. Run inside
the allocated resources, with models on mounted storage where required:

```sh
python tests/functional/gospl/run_mpi_coupling.py \
  --cpu /bin/dynearthsol3d.cpu.mpi.gospl \
  --gpu /bin/dynearthsol3d.gpu.mpi.gospl \
  --seed /mounted/coupling-reference/cpu --extension /opt/gospl_extensions \
  --gospl-python /gospl/env/bin/python \
  --cpu-launcher '/des/mpi/bin/mpiexec --bind-to none' \
  --gpu-launcher '/gpu/mpi/bin/mpiexec --bind-to none' \
  --gospl-launcher '/gospl/env/bin/mpiexec --bind-to none' \
  --workdir /mounted/results
```

Use `--worker-env NAME=VALUE` and `--gpu-env NAME=VALUE` to keep MPI/SDK runtimes
separate. The worker interface audit delegates to the real bridge and records
per-rank seed/erosion arrays, clocks, durations and cleanup. The matrix covers
DES CPU 1/3, GPU 2, GoSPL 1/2, time/step triggering and CPU/GPU remeshing. It
compares all physical fields at fixed GoSPL rank count, not GoSPL's own algorithms
across rank counts. Owned and ghost seed placement is checked independently.
An injected protocol error must terminate and reap both process groups, with
worker status 53. `--resume /existing/artifact --cases NAME` appends selected
checks without rerunning completed models; executable hashes must match.

The [MPI coupling guide](../../../doc/gospl-mpi-coupling.md) describes allocation,
root placement and measured scope. The tested PETSc build uses CPUs, so the
GPU case accelerates DES mechanics only.

## Continuation and physical observations

`run_continuation.py` runs 12 steps continuously and resumes checkpoints before
the first exchange and between later exchanges. Supply `--exe`, `--seed`,
`--extension`, `--gospl-python`, `--workdir`, `--launcher` and `--gospl-launcher`
as in the MPI check. Its launcher prefixes omit `-n`; the driver adds `--ranks`
and `--gospl-ranks` (defaults 3 and 2). `--restart-ranks 2` changes only DES's
continuation rank count; `--remesh` and `--time-mode` select the targeted cases.
`--worker-env` and `--des-env` separate runtimes. It compares every subsequent
DES field/checkpoint plus final owned/local GoSPL physical vectors and clocks.
`--verify-existing PATH` rechecks retained results without launching a model.

`reject_checkpoint.py` takes an existing `--continuous` result, the executable,
environment and launcher arguments above, and tests private copies for checksum
damage and a changed GoSPL rank count. Both must fail for the expected reason
and return worker status 53 after cleanup; the source checkpoint is unchanged.

`mesh_observables.py REFERENCE.vtkhdf CANDIDATE.vtkhdf --output report.json`
compares 3D tetrahedral outputs at the same physical time without requiring node
or element identity. It checks finite fields, nondegenerate consistently oriented
tetrahedra, common interior/top coverage, volume-weighted cell statistics and
interpolated fields/elevation. Sampling uses an interior box and a single-valued
top surface, appropriate to the retained fixture. Reported differences need
case-specific scientific interpretation; they are not a universal tolerance.

`benchmark_coupled.py` runs one shared initial checkpoint with CPU OMP1/OMP4,
CPU MPI3 and GPU1/GPU3, all coupled to GoSPL1. It reserves distinct physical CPUs
for DES and the worker, warms each variant once, then interleaves three measured
rounds. Every output is compared. `--cpu`, `--gpu`, `--seed-config`, `--extension`,
`--gospl-python`, `--workdir` and the three launcher prefixes are required.
Use `--worker-env`/`--gpu-env` to isolate MPI libraries. Initialization, compute,
output, full launcher wall time and worker-call timings are separate; overlapping
timers must not be added. `--summarize-existing PATH` parses retained timing logs
without rerunning models.

Verification scope and limitations are in the
[coupling guide](../../../doc/gospl-mpi-coupling.md#verification-scope).
