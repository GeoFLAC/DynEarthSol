# Force assembly, PT constraints and material trials

Failed-trial checks cover both tolerance modes. Running-scale PT must preserve
the committed residual maximum, pending remesh exclusion, velocities and material
history when an iteration cap rejects the candidate.

Run from this directory:

```sh
make check ndims=2 backend=asan
make check ndims=3 backend=asan
OMP_NUM_THREADS=4 make check ndims=2 backend=omp
OMP_NUM_THREADS=4 make check ndims=3 backend=omp
make check ndims=2 backend=acc
make check ndims=3 backend=acc
```

OpenACC requires NVHPC, the production build dependencies and a GPU. `GPU_CC`
defaults to the detected device and can be overridden. Accelerator checks use the
root build graph; do not run them concurrently with another root build. The test
entry point replaces only application `main`; all physical routines are real.

Two positively oriented simplices share one vertex. Checks cover cancellation,
Neumann traction, gravity with porous bulk density, and artificial damping.
The optional output is the unprojected assembled physical force before damping.
Supplying it must preserve both existing outputs bit for bit. Analytic assertions
use `abs(x-y) <= 1e-9 + 1e-14*abs(x) + 1e-14*abs(y)` in the fixture's force units.

`mode=BASE` calls the old five-argument API and measures the legacy residual
overwrite defect. The default `mode=B0` exercises the optional undamped output.
Both modes now require the current source tree for the added PT tests. The fixture
also checks affine constraint fixed points, projector symmetry/idempotence/rank,
oblique horizontal normals and imposed edge motion; pressure-only responses for
elastic/Maxwell/viscous/EP/EVP at Biot coefficients 0, 0.6 and 1; plane-strain yy;
repeated shear and yielded plastic trials; and frozen initial RSF properties.
It does not qualify full-solve convergence, remeshing or physical RSF aging.

`../functional/pt_initial_checkpoint.py EXECUTABLE CONFIG OUTPUT_DIRECTORY`
uses a supplied short binary-output PT case to check that frame-zero initialization
is complete, then compares all physical save/checkpoint fields after a frame-zero
restart. Only wall-clock metadata is excluded from byte comparisons.

The initial fixture was authored by Qwen3.8-27B-FP8 under supervision. Codex added
the real material/body-force fixture, null-output equivalence checks and build
drivers, and independently validated the baseline and changed implementations.

The PT constraint/material and initial-checkpoint tests were added by Codex.

The PT residual reduction is tested directly with first/last-node loads, unit,
1e200 and 1e-200 force scales, zero free rank and nonfinite inputs. Integration
comparisons must also use a common mesh and a nonzero initial imbalance: a solver
exit code alone cannot rule out falsely reported convergence on a backend.

The material fixture also calls the actual PT solver with iteration caps 1 and 3,
zero stopping tolerance, nonzero shear velocity and a pending pressure increment.
For elastic/Maxwell/viscous/EP/EVP it requires a positive residual and iteration-limit
failure, then checks exact restoration of stress, strain, all seven captured scalar
histories and velocity. Pressure inputs, physical time/dt and array identities must
remain unchanged. The three-iteration case exercises numerical updates before
rollback. These are rejection/lifecycle checks, not convergence benchmarks.

A planar Winkler fixture checks the trial-height load against the independent
spring resultant `Delta Fz = -rho*g*area*dt*vz`, for dry and hydraulic boundary
densities in 2D/3D. Repeated dt values and a return to dt=0 must recover the same
forces without changing coordinates. The failed-solve lifecycle fixture enables
surface diffusion and requires both the per-step `dh` and accumulated `dhacc`
sentinels to remain untouched by PT iterations. Accepted physical surface updates
are checked separately by the integration lifecycle, not by these force tests.

The Darcy fixture checks the zero-gravity limit of the actual hydraulic update.
Constant hydraulic potential must remain stationary; a nonconstant pressure field
must diffuse with zero net fluid change on the closed two-element mesh. Adding a
hydrostatic background at nonzero gravity must leave the pressure increments and
diffusivity unchanged. These checks run in both dimensions and all test backends;
they do not establish the accuracy of the coupled consolidation time integration.
