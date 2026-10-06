# PT refactor and adaptive relaxation submission record

Date: 2026-10-06 (Asia/Seoul).

## Scope and source identity

`refactor/pt` continues the refactor based on `308d8c3`, through reference
`796fad2` and adaptive relaxation `69c8712`. It is not based on the old public
`feature/pseudo-transient` branch (`911a1b5`). That branch is preserved as
`archive/pseudo-transient-20261006` in the public repository.

The finalization adds ADR rejection coverage and a selectable EP/RSF benchmark
mode, repairs two fixture-owned aliases exposed by sanitizers, and clarifies
algorithm limitations. It does not change the relaxation arithmetic or add
another acceleration scheme. The default is `PT_option=0`; ADR is opt-in with
`PT_option=1`. The existing constitutive update remains authoritative, physical
time is fixed within a solve, and only an accepted trial commits material
history. Failure restores the physical-start state. Reference geometry and
physical mesh/transport stages retain the ownership described in [pt-solver.md](pt-solver.md).

The public master was checked at `7cd88a1`; its newer provenance changes have
not been replayed into this branch. The three-way textual merge check found no
conflict. Integration with those changes remains a PR/CI review item.

## Environment and acceptance criteria

- WSL2 Linux 6.6.87.2, GCC 11.4.0, Python 3.11.14.
- CPU serial debug: `opt=0 openmp=0 hdf5=0`; OpenMP: `opt=2 openmp=1 hdf5=0`.
- GPU: NVIDIA RTX 4070 SUPER, NVHPC 25.9, `openacc=1 GPU_CC=89 hdf5=0`.
- Existing focused assertions retain their analytical tolerances and exact
  rollback checks. ASan/UBSan runs halt on undefined behavior.
- Restart harnesses retain byte-exact comparisons (excluding wall-clock data).
  A failed exact comparison is reported as such below, not silently relaxed.
- Existing `benchmarks-cores/compare.py` uses its established 1e-8 criterion.
  EP/RSF uses the unchanged analytical limits, including 2e-5 stress error.
- These checks establish scoped correctness, not a new performance claim.

## Results reproduced during finalization

| Check | Scope | Result |
| --- | --- | --- |
| Focused force/constraints/material/rollback | 2D/3D ASan/UBSan and OpenMP, both PT options and both tolerance modes | 5,740 / 8,547 assertions pass |
| Serial restart | Both options, 2D/3D dynamic dt, hydraulic remesh and hydrostatic stationarity; 2D oblique constraints | Exact continuation; remesh connectivity changes required |
| Initial checkpoint | Both options, 2D/3D serial hydrostatic fixture | Exact continuation and stationary physical fields |
| Output scheduling | Both options, 2D/3D serial and OpenMP; step/time/mixed/catch-up | Exact continuation |
| OpenMP restart | Four threads, original relaxation; ADR dynamic/oblique cases | Exact continuation |
| OpenMP ADR remesh/hydrostatic | Four threads, 2D/3D | Four byte-exact checks fail; maximum field-relative difference 1.36e-14 |
| OpenMP ADR remesh/hydrostatic | One thread, 2D/3D | Exact continuation; hydrostatic stationarity and remeshing checks pass |
| Default PT versus reference `796fad2` | 2D/3D dynamic and hydraulic-remesh cases, all saved/checkpoint payloads, one thread | Byte-identical apart from wall-clock metadata |
| Dry non-PT regression | Parent/base `make set`, final `make cmp`; 2D `test-rect-tiny.cfg`, 3D `test-3d-equ-tiny.cfg` | BIT-EXACT against both `308d8c3` and `796fad2` |
| EP/RSF analytical suite | Both options, nine 200,000-step shear cases plus healing, CPU one thread | Pass; maximum stress error 1.156e-5 as a fraction, rate error 5e-8, healing error zero |
| OpenACC focused tests | 2D/3D, both options and tolerance modes | 5,740 / 8,547 assertions pass |
| OpenACC dry restart | Both options, 2D/3D dynamic case; 2D oblique case | Exact continuation |
| OpenACC ADR dry remesh | 2D/3D, three physical steps and restart from frame 1 | Exact continuation with changed connectivity |
| OpenACC original-relaxation dry remesh | 2D/3D | Incomplete: original 240-second fresh-run limits reached |

The longer original-relaxation 2D GPU attempt used `--timeout 900` and the same
three-step input, iteration cap and tolerance. It produced frame 1 but was
stopped before finishing; it is not a passing restart result. The full legacy
GPU remesh workload remains unqualified here. Some GPU correctness runs
shared the device, so none of their wall times is a performance measurement.

The four-thread ADR differences are consistent with reduction-order sensitivity
in the adaptive damping. For the
reported runs, `max(abs(a-b))/max(max(abs(a)),max(abs(b)))` over each changed
physical field is at most 1.355654502864176e-14. Integer/topology fields and the
restored velocity before continuation agree. Serial and one-thread continuation
are exact. This is evidence for these fixtures, not a determinism guarantee for
other workloads or thread counts.

The EP/RSF cells are prescribed-motion constitutive checks and accept their
first candidate; they do not establish multi-iteration convergence. The
traction-driven dynamic/remesh checks supply that separate coverage.

## Reproduction

Run dimensions sequentially in one worktree. See the existing fixture
[README](../tests/residual-assembly-test/README.md) and
[restart instructions](../tests/functional/pt_restart.md) for all commands.

```sh
make -C tests/residual-assembly-test check ndims=2 backend=asan
make -C tests/residual-assembly-test check ndims=3 backend=asan
OMP_NUM_THREADS=4 make -C tests/residual-assembly-test check ndims=2 backend=omp
OMP_NUM_THREADS=4 make -C tests/residual-assembly-test check ndims=3 backend=omp
make -C tests/residual-assembly-test check ndims=2 backend=acc
make -C tests/residual-assembly-test check ndims=3 backend=acc

OMP_NUM_THREADS=1 python3 benchmarks/simple_shear_rsf/check_simple_shear_benchmark.py --exe ./dynearthsol2d --all --pt --pt-option 0 --keep-output
OMP_NUM_THREADS=1 python3 benchmarks/simple_shear_rsf/check_simple_shear_benchmark.py --exe ./dynearthsol2d --all --pt --pt-option 1 --keep-output
```

For functional cases, copy the supplied cfg to a new directory, set
`control.PT_option` to 0 or 1, and use the existing checkpoint and schedule
runners. Resolve the oblique fixture's polygon path before moving its cfg.
GPU remesh cases set `has_hydraulic_diffusion=false` and only claim dry coverage.
The byte-exact ADR restart checks use `OMP_NUM_THREADS=1`; the four-thread
results above are kept separately. All runs use new output directories.

For dry baseline comparisons, use separate parent/final executables and the
same `benchmarks-cores` inputs. `make set` runs the base executable, then
`make cmp` runs the final executable with `NDIMS=2 CASE=test-rect-tiny.cfg` or
`NDIMS=3 CASE=test-3d-equ-tiny.cfg`, `OMP=1`, and the corresponding `EXE` paths.
The baseline executables used `opt=2 openmp=0 hdf5=0`; final production
executables used `opt=2 openmp=1 hdf5=0`, one thread.

Local full logs, generated configs, executables, comparison JSON and orchestration
scripts are retained under
`DynEarthSol-worktrees/test-runs/pt-final-20261006/`. Initial failing attempts
are retained alongside the corrected fixture and isolated comparison reruns.
The failed initial fixture used unbound scratch/strain-rate aliases; the final
fixture fixes those aliases. Only the isolated `make set/cmp` reruns are cited
above. No raw simulation payload is included in the source branch.

## CPU/GPU field comparison

The dynamic fixture applies a nonzero top traction on a fixed common mesh.
Compare both relaxation options across all saved physical frames. For each
field, use the maximum absolute difference divided by the maximum magnitude
over the compared history, with the existing 1e-8 numerical criterion. This
avoids using only a near-zero terminal velocity as the normalization scale.
Topology must agree exactly; report absolute differences as well.

Both options pass in 2D/3D. The largest history-normalized difference is
6.548e-16; coordinates and connectivity agree exactly. Absolute differences,
scales and per-field results are retained in `gpu-history-comparison.json` in
the local evidence directory.

## Earlier measurements and limits

The performance numbers in the original reference/ADR commits are historical
implementation measurements. Their H100 and serial CPU timings were not
reproduced here, and should not be interpreted as RTX 4070 SUPER measurements.
This submission does not qualify large-mesh performance, general nonlinear
convergence, full hydraulic GPU coupling, or long-time core-complex morphology.
HDF5, MMG, Exodus and GoSPL were not rerun in this finalization; earlier coverage
is described separately in the solver documentation. Earthquake-event history
and averaged-output accumulator restart remain outside the checkpoint fixtures.

The reference stiffness row-sum bound does not prove stability of the full
nonlinear constrained operator. Softening may remove a stable quasi-static
equilibrium. Running-scale convergence remains opt-in and is not a general
qualification of history-dependent materials. The separate pressure-unit Darcy
flux correction in this series changes hydraulic behavior; the dry non-PT
BIT-EXACT results do not claim unchanged hydraulic results.
