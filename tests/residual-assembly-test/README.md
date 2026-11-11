# Force assembly, PT constraints and material trials

Failed-trial checks cover both relaxation options and both tolerance modes. Running-scale PT must preserve
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
Supplying it must preserve both existing outputs bit for bit. The PT damping
bypass must return the raw force; damping it once must recover the legacy result. Analytic assertions
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

The constant-load integration regression in `../functional/pt_residual_output.py`
checks both held equilibrium and a load doubled on restart. With the supplied
running-scale elastic example, steps 2/3 must be accepted without relaxation
updates. `../functional/pt_maxwell_creep.py` checks continued traction-driven
Maxwell creep for 100 steps against a homogeneous discrete reference. Commands
and the meaning of candidate evaluations are in `../../doc/pt-solver.md`.

Mohr-Coulomb checks prescribe plastic increments and reconstruct elastic trial
stresses from final stresses on a single face or either edge. Positive and
negative cohesion slopes, dilation angles 0/30 degrees, and equal/unequal edge
multipliers must return to those targets. Rebuilding strength from the recorded
plastic history and applying a second zero-strain step must preserve yield and
history. Plane-strain cases cover tensile clipping and states at and on either
side of the tensile cap. These are local 3D and plane-strain tests within a linear cohesion segment,
with fixed friction/dilation. The pure 2D, two-principal-stress formulation is
separate; general apex and fully coupled shear/tension evolution are not qualified.

`check.sh` also runs a strongly softening case in a separate process. An
inadmissible return must report its diagnostic and exit with status 54 on every
test backend; an out-of-order candidate must not continue silently. The material
kernel completes before the host checks this return status.

For a fixed friction angle, the return uses `H = 2*sqrt(Nphi)*dc/d(plstrain)`.
Its strength increment is `H*depls`, using exactly the increment stored in history.
On an edge, isotropic elasticity fixes the difference `d` of the two face
multipliers. Their sum `s` satisfies `E*s + H*depls(s,d) = R`. The recorded
plane-strain measure is piecewise linear in `s,d`; the 3D measure has
`depls^2 = A*s^2 + B*d^2`. The return solves that scalar consistency equation,
then checks multiplier signs, stress ordering, tensile bounds and yield.
