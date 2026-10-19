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

The initial submission did not include master's newer provenance changes.
The follow-up described below merges public master `7cd88a1` without rewriting
the shared PT history; the merge had no textual conflicts.

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

## Follow-up cleanup and master integration

The first cleanup removes ADR's unused `lagged_stress` allocation and initialization,
along with obsolete commented-out force and main-loop code. It leaves the shared
factor calculation and its material/mesh guards unchanged, as well as the
relaxation arithmetic, acceptance criteria and rollback state.

The merged master adds per-process provenance to binary output/checkpoints and
HDF5 attributes. PT restart and schedule checkers exclude the binary provenance
record, alongside wall time, from physical-state equality. They continue to
compare all physical fields exactly. This matches their existing HDF5 behavior,
where provenance attributes are outside the physical datasets.

The pre-cleanup submission binaries (`b6d233f`, with the same production code
as the recorded finalization builds) provide the comparison baseline. Both PT
options, 2D/3D serial dynamic and hydraulic-remesh cases, retain byte-identical
physical payloads across all 64 saved/checkpoint files. The dry non-PT
`make set/cmp` checks also remain BIT-EXACT for `test-rect-tiny.cfg` and
`test-3d-equ-tiny.cfg`, using OpenMP executables with one thread.

The integrated sources build in 2D/3D serial debug, OpenMP and OpenACC modes.
Focused ASan/UBSan, four-thread OpenMP and OpenACC checks again pass 5,740 / 8,547
assertions. Both options pass serial initial hydrostatic equilibrium, dynamic
restart, hydraulic remesh, all four output schedules and 2D oblique constraints;
one-thread OpenMP dynamic restart and hydraulic remesh also pass in 2D/3D.
Additional 2D/3D HDF5 OpenMP builds pass both options' dynamic restart and mixed
output-schedule checks. Physical continuation remains byte-exact throughout.

OpenACC passes both options' 2D/3D dynamic restart, both 2D oblique cases and
ADR dry remesh in both dimensions. The corresponding 64 pre/post-cleanup GPU
save/checkpoint files also have byte-identical physical payloads. A separate
audit checks the fresh/restarted source identity in 1,238 binary/HDF5 provenance
records across the 54 functional runs; all agree with their run configuration.

Follow-up commands, binaries, source hashes and comparison records are retained
locally under `DynEarthSol-worktrees/test-runs/pt-cleanup-20261006/`.

The subsequent cleanup, based on `c66634c`, also leaves ADR's `stress_fraction`,
`mobility` and `height` arrays empty and skips their unused coefficient outputs.
Element-height/theta and nodal-support calculations needed for the existing
validation remain; the rejection conditions and effective-viscosity arithmetic
are unchanged. One shared refresh condition now covers factors and ADR mass.
No new solver controls or helper functions are introduced.

The factor cleanup repeats 2D/3D serial, OpenMP and OpenACC builds and focused
checks (5,740 / 8,547 assertions for each backend, including ASan/UBSan).
Dry `make set/cmp` remains BIT-EXACT. Both options' small serial Maxwell cases
retain identical physical payloads, iteration counts, residuals and thresholds
in both dimensions. This is a short regression, not a Maxwell/Rayleigh–Taylor
qualification. The 25 completed CPU/GPU comparison cases retain byte-identical
physical state across 200 save/checkpoint files. Binary CPU restart/remesh and
GPU fresh dynamic/ADR-remesh comparisons are retained under
`test-runs/pt-factors-cleanup-20261006/`.

One initial 3D serial ADR restart aborted with a stack-smashing diagnostic in
pending remeshing, before the next PT solve. Its eight fresh physical payloads
match the parent exactly. The same checkpoint then passed ten restart attempts
with each of the parent and changed binaries, all with exact continuation.
A full AddressSanitizer build also passes fresh/remesh/restart and matches the
parent's physical output. The original abort has not been reproduced or assigned
a cause; its failed log/result remains in the evidence and is not counted as a
passing restart. A separate full ASan/UBSan attempt stops earlier in unchanged
TetGen (`tetgen.h:3225`, left shift of negative value during initial meshing).
That full UBSan check is incomplete; the focused PT ASan/UBSan suite passes.

## Earlier measurements and limits

The performance numbers in the original reference/ADR commits are historical
implementation measurements. Their H100 and serial CPU timings were not
reproduced here, and should not be interpreted as RTX 4070 SUPER measurements.
This submission does not qualify large-mesh performance, general nonlinear
convergence, full hydraulic GPU coupling, or long-time core-complex morphology.
MMG, Exodus and GoSPL were not rerun in this finalization or follow-up; earlier
coverage is described separately in the solver documentation. HDF5 follow-up
coverage is limited to the output/restart cases above. Earthquake-event history
and averaged-output accumulator restart remain outside the checkpoint fixtures.

The reference stiffness row-sum bound does not prove stability of the full
nonlinear constrained operator. Softening may remove a stable quasi-static
equilibrium. Running-scale convergence remains opt-in and is not a general
qualification of history-dependent materials. The separate pressure-unit Darcy
flux correction in this series changes hydraulic behavior; the dry non-PT
BIT-EXACT results do not claim unchanged hydraulic results.


## Review follow-up: force output, state ownership and example (2026-10-19)

This follow-up addresses unused PT damping, the indexed mechanical snapshot,
and the missing runnable example. Trial force assembly skips artificial damping;
the legacy damped `force` is computed once on acceptance. The accepted projected
undamped vector is saved as optional `PT residual force`. Diagnostic validity is
cleared on failure and buffer reallocation; save files may restore it without
adding a checkpoint requirement. The existing hydraulic restart policy is not
changed by this follow-up. `MechanicalState` now lives in the PT module and uses
seven named scalar arrays rather than numbered columns.

The new `examples/pt-elastic-load.cfg` is a dry, fixed-mesh elastic compression
with a nonzero initial imbalance, three 0.01 s physical steps, and both PT options.
For 2D, the computed equivalent shortening is approximately 4.286 mm and the
vertical stress is approximately -100 kPa. See the solver notes for the setup,
reference modulus and fixed-mesh displacement interpretation.

Validation used GNU 11.4, NVHPC 25.9 and an RTX 4070 SUPER (cc89):

- Both dimensions compile in serial debug (binary output), OpenMP O2
  (binary/HDF5), and OpenACC (HDF5) configurations.
- The focused assembly/material/rollback suite passes 6,832 assertions in 2D
  and 9,639 in 3D for each of OpenMP (four threads), serial ASan/UBSan, and
  OpenACC. These sanitizer results cover the focused suite, not a new full
  application sanitizer qualification.
- The new force-output regression passes in both dimensions and both PT options
  for serial debug, OpenMP binary, OpenMP HDF5 and GPU HDF5. It independently
  removes damping, applies the example's axis constraints, checks free-rank RMS
  against the reported residual (`rtol=1e-12`, `atol=1e-14`), and checks exact
  diagnostic persistence on restart.
- Compared with the parent executable, all existing physical save/checkpoint
  payloads, legacy force fields, PT iteration counts and residual reports are
  bit-exact for the four single-thread elastic examples. Only runtime and
  provenance records are excluded.
- Parent `make set` then changed `make cmp` passes BIT-EXACT for
  `NDIMS=2 CASE=test-rect-tiny.cfg` and
  `NDIMS=3 CASE=test-3d-equ-tiny.cfg`. Additional full-payload comparisons of
  their six and ten save/checkpoint files also find no physical differences.
- Twelve single-thread binary lifecycle cases pass: both dimensions and both
  options, each with dynamic dt, remesh/restart, and stationary hydrostatic
  initialization. Separate checks confirm diagnostic presence before remeshing
  and omission immediately afterward. Parent dry saves without the new field
  restart successfully for both dimensions/options with bit-exact continuation
  of existing fields. This does not change or qualify older hydraulic restart
  compatibility.
- CPU/GPU examples use identical reference coordinates and connectivity.
  Across all four examples and saved frames, the maximum normalized difference
  in stress, strain or velocity is 2.91e-16. Normalizers are 1e5 Pa, 0.01, and
  1 m/s respectively; the declared acceptance limit was 1e-5. Both backends
  independently satisfy their residual criteria.
- Binary-to-VTK conversion handles the optional vector, including its absence
  from the unsolved initial frame.

Reproduce the new output check from the repository root, using fresh output
paths and the appropriate executable for each build:

```sh
OMP_NUM_THREADS=1 python3 tests/functional/pt_residual_output.py ./dynearthsol2d examples/pt-elastic-load.cfg /tmp/pt-review-2d-0 --option 0
OMP_NUM_THREADS=1 python3 tests/functional/pt_residual_output.py ./dynearthsol2d examples/pt-elastic-load.cfg /tmp/pt-review-2d-1 --option 1
```

Repeat with `dynearthsol3d`. GPU checks use `dynearthsol2d.gpu` or
`dynearthsol3d.gpu` and `--timeout 180`; the initial 60 s 2D GPU check exceeded
its wall-clock limit, while all final checks passed. This changes only the test
process timeout, not the model's iteration cap or convergence settings.
Focused-suite commands and lifecycle commands are in
`tests/residual-assembly-test/README.md` and `tests/functional/pt_restart.md`.

No long runs or new acceleration experiments are part of this follow-up.
The short elastic results do not qualify general nonlinear convergence,
long-time plastic/RSF behavior or hydraulic accuracy. MMG, Exodus and GoSPL
configurations were not rerun for this follow-up.
