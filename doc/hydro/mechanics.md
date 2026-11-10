# Mechanical method qualification for coupled models

Ordinary explicit mechanics and global velocity scaling (GVS) are primary
comparison methods alongside PT. An instantaneous load does not by itself
require PT, and successful earthquake-cycle integration does not by itself
establish quasi-static Mandel/Terzaghi accuracy. Compare the equations,
initial state, mechanical relaxation time, drainage time and physical
observables. Keep inaccurate, unstable and timed-out cases in the report.

## Selected 3D comparison

The coupled 3D check uses the same Terzaghi column and independent series
reference as the 2D study: a 1 m by 1 m by 10 m confined column, sealed sides
and base, drained top and a sudden 100 kPa top traction. K=20 MPa,
G=10 MPa, zero gravity, fixed mesh and no heat. The
[column schematic](../../benchmarks/poromechanics/results/review/fig/terzaghi/setup.svg)
shows its vertical section; both horizontal directions are confined here.
The [reference and convergence definitions](refinement.md) give pressure
and settlement equations. Reference values use each frame's actual time.

The selected comparison uses PT1 with physical dt=1 s and shear-ceiling
GVS with S=1e5, initial velocity scale 1e-5 m/s, `fixed_dt=0` and
`dt_fraction=0.5`. The independent force diagnostic assembles internal
forces from saved stresses and removes constrained components.

| Method | Mesh spacing | Late pressure/load error | Late settlement/drained error | Final free-node force/total load |
| --- | --- | --- | --- | --- |
| PT1 | 1 m | 0.4915% | 0.1679% | 9.62e-9 |
| PT1 | 0.5 m | 0.1706% | 0.1011% | 5.39e-9 |
| GVS shear, S=1e5 | 1 m | 0.8916% | 0.1551% | 5.00e-7 |
| GVS shear, S=1e5 | 0.5 m | 0.2196% | 0.03967% | 2.51e-7 |

Late means saved frames from 200 s through the final frame near 400 s.
The drainage time is 2,933.74 s. Maximum pressure/load error across the
entire saved 40–400 s window on the fine mesh is 0.8685% for GVS and
0.8417% for PT1. Halving the spacing also changes the automatic GVS
timestep: this is combined mesh/step evidence, not independent temporal
convergence. The finer mesh has 189 nodes and 400 tetrahedra.

![Selected 3D pressure and settlement errors](../../benchmarks/poromechanics/results/runs/gvs-qualification-20261108/three-d-h05/mechanical-methods.svg)

The retained coarse matrix contains all eleven PT0/PT1, ordinary explicit
and GVS shear/bulk cases at S=1e4–1e6; all finish. Ordinary explicit at
S=1e4 has 2.086% late pressure error and force/load 0.01245. That larger
inertial error remains part of the result. The full matrix is linked from
its [study index](../../benchmarks/poromechanics/results/runs/gvs-qualification-20261108/three-d-h1/index.html).
This check does not qualify cold starts, the full drainage history,
Mandel/Cryer in 3D or earthquake-cycle transitions.

The study uses a frozen executable identified by its qualification manifest,
separately from the checkout used to run the analysis. Reproduce with new
output directories:

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
V=benchmarks/poromechanics/results/runs/gvs-qualification-20261108
python3 benchmarks/poromechanics/mechanical_methods.py /tmp/des-gvs-3d-h1 \
  --dimension 3 --executable "$V/final-regression/binaries/dynearthsol3d" \
  --binary-manifest "$V/qualification.json" --spacing 1 --duration 400 \
  --dt-fraction .5 --pt-dt 1 --initial-speeds 1e-5 \
  --scalings 1e4 1e5 1e6 --methods pt0 pt1 explicit gvs-shear gvs-bulk \
  --timeout 300 --max-steps 10000000
python3 benchmarks/poromechanics/mechanical_methods.py /tmp/des-gvs-3d-h05 \
  --dimension 3 --executable "$V/final-regression/binaries/dynearthsol3d" \
  --binary-manifest "$V/qualification.json" --spacing .5 --duration 400 \
  --dt-fraction .5 --pt-dt 1 --initial-speeds 1e-5 \
  --scalings 1e5 --methods pt1 gvs-shear --timeout 600 --max-steps 10000000
```

## Pressure and strength sampling in PT output

In the fixed-stress PT ordering, mechanics evaluates strength with $p_n$,
then flow advances to $p_{n+1}$. The saved pair therefore combines the
new mechanical stress with the new pressure, although the return used
the preceding pressure. This is the one-step coupling lag; hydraulic
subcycling does not remove it. The explicit/GVS strain-form ordering is
different, as described below.

For saturated Mohr–Coulomb material with tension-positive total stress,
let $N=(1+\sin\phi)/(1-\sin\phi)$ and $q_c=2c\sqrt N$. With the
pressure mechanics used, the yield slack is

$$
f_n=\sigma_{\min}+\alpha p_n
       -N(\sigma_{\max}+\alpha p_n)+q_c.
$$

Using the subsequent pressure instead changes that diagnostic by

$$
f_{\mathrm{saved}}=f_n-(N-1)\alpha(p_{n+1}-p_n).
$$

Thus a rising pressure produces a negative saved-frame slack even when
the return is admissible at its own sampled pressure. Both the lag error
and the constitutive return need their own checks.

A read-only screen of the recorded [slope application](slope.md), with
$h=0.25$ m, $\phi=35^\circ$ and $\alpha=1$, gives the following minimum
**in-plane** slack in the middle third, at the last accepted state:

| Physical steps | Using saved pressure | Using preceding pressure |
| --- | --- | --- |
| 404 | -205.617 Pa | -0.000744 Pa |
| 808 | -102.843 Pa | -0.000226 Pa |

The saved-pressure deficit halves with the step. The $h=0.5$ and
$0.125$ m states show the same timing effect; the largest remaining
preceding-pressure discrepancy is 0.003824 Pa. This is bounded evidence
of first-order sampling lag, not a full coupling-convergence study.
No late checkpoint was written, so $\sigma_{yy}$ is unavailable and
full sorted-principal admissibility cannot be established from those
files. The independent affine [edge-return failure](salencon.md) is a
separate defect, reproduced without smoothing or a free force residual.

The exact inputs, saved frames, reader and analysis hashes are in
`results/runs/slope-application-20261108/yield-screen.json`. Replay without
running DES:

```sh
python3 benchmarks/poromechanics/results/runs/slope-application-20261108/audit-yield.py
```

## Restart output clock

The 3D adaptive GVS restart check exposed an extra frame immediately after
restart. At common physical steps, 311 saved/checkpoint payloads were exact;
the output cursor had not been persisted for non-PT runs. The corrected
output state retains the epoch, next regular frame and configured intervals
for every method. An unchanged cadence preserves the frame sequence; a
changed cadence starts at the restart time/step. PT keeps its existing
checkpoint field names. This change does not change the mechanical or
hydraulic update, fictitious mass, timestep or energy equation.

The corrected adaptive 3D replay restarts at step 251,498, about 200 s.
All six continuation frames have the original frame number, time and step;
305 payloads in ten subsequent save/checkpoint files match exactly. All
396 physical save payloads across eleven fresh frames also match the
previous solver. Fixed-dt schedule checks cover PT, ordinary explicit, GVS
and RSF, both containers and dimensions: 160 cases compare 2,128 files,
with eight additional PT1 cases. Both dry core comparisons remain exact.
The combined injection-temperature build also passes twelve schedule checks.
See the
[restart commands and compatibility limits](../../tests/functional/pt_restart.md#regular-output-schedule).
Evidence is retained under `results/runs/output-clock-20261108/`.

Legacy checkpoints lacking schedule/cadence fields retain their documented
fallback. Earthquake-event history and averaged-field accumulators remain
separate checks. No changing-mass energy balance or general
variable-timestep qualification is claimed here.

### Changing mechanical method or timestep controls

A PT gravity stage can have a much larger physical step than the following
explicit stage permits. The checkpoint now records the timestep controls;
changing method, physical inertia, scaling, characteristic speed, step
fraction or fixed step reselects the current timestep and mass before the
first continued step. Unchanged controls retain exact physical continuation.

Forty-eight native/HDF5 cases cover tiny 2D/3D cells, including older PT
checkpoints. Independent mesh-height wave limits and the explicit
velocity/mass identity check the selected interval and mass. A saved 1 s
PT interval becomes 0.05 s in the 2D explicit cell and 0.0333 s in 3D;
the parent wrongly keeps 1 s. Physical-inertia conversion uses the real
density and its wave limit. Both parent dry core comparisons remain exact.
See [commands and legacy limits](../../tests/functional/pt_restart.md#changed-timestep-controls).
Evidence: `results/runs/restart-timestep-20261108/qualification.json`.

## What the scaling parameter controls

In `geometry.cxx`, let `S=control.inertial_scaling` and `v` be the velocity
scale. Quasi-static mass uses `rho_fict=K/(S*v)^2`, subject to the GVS ceiling.
Thus increasing S reduces fictitious mass at fixed v. S is not the density
multiplier commonly called mass scaling in other explicit codes.

Without GVS, v is the boundary/characteristic velocity. With GVS it is the
maximum norm of the **element-mean** velocity, bounded below by that initial
scale. A traction-only model with no `control.characteristic_speed` starts
from the `1e-12 m/s` boundary-velocity floor. Record this value: S alone
does not establish the mechanical timescale. A starting scale must be
qualified against the load response, not selected only to obtain a pass.

| GVS setting | Speed ceiling | Fictitious density floor |
| --- | --- | --- |
| `control.mass_scaling_reference_speed=shear` | `sqrt(G/rho)` | `rho*K/G` |
| `control.mass_scaling_reference_speed=bulk` | `sqrt(K/rho)` | `rho` |

Here rho is `MatProps::rho`, the solid reference/thermal density. With
hydraulics enabled, the physical dynamic branch instead uses the saturated
mixture density. Consequently the bulk ceiling alone does not establish
correct saturated coseismic inertia. It is also a bulk-modulus speed,
not the full compressional-wave speed including shear or undrained coupling.
Check those wave speeds and the coupled timestep independently.

The physical-inertia branch now uses a strain-driven fluid equation and an
automatic undrained-wave timestep limit. The
[closed-column qualification](dynamics.md) establishes its
small-strain saturated wave response in both dimensions. Coupled explicit/GVS
now uses the same strain-driven pressure equation, with its existing mass
ceiling and an undrained timestep based on fictitious density. The physical
wave check does not qualify a GVS earthquake-cycle transition.

`control.fixed_dt` bypasses `compute_dt`'s adaptive velocity-scale update.
Use `fixed_dt=0` for this GVS comparison. Heat and hydraulic substeps have
their own bounds; those do not establish mechanical stability or accuracy.

## Required comparisons

| Case | Methods / settings | Evidence needed |
| --- | --- | --- |
| Dry confined sudden traction | Ordinary explicit, GVS shear/bulk, PT 0/1; S=1e4/1e5/1e6 | Velocity/mass identity, force relaxation, stress/displacement, step refinement and restart |
| Terzaghi | Same method/scaling matrix; separate cold and load-informed starting scales | Pressure and settlement at actual physical times; force imbalance; initial loading delay; mechanical/drainage time separation; independent mesh/time refinement |
| Mandel / Cryer | Same comparisons with the correct plate/resultant or spherical boundary conditions | Overshoot amplitude/time, displacement and resultant; repeat with resolved mechanical transients |
| Tectonic loading | S=1e4/1e5/1e6, ordinary explicit/GVS and a converged equilibrium reference | Insensitivity of slow observables to scaling, timestep and damping; force and kinetic/strain-energy histories |
| Earthquake cycle | Adaptive GVS through interseismic/coseismic transitions and a resolved physical-inertia reference | Rupture speed, slip, stress drop, recurrence, mass/dt transitions, damping sensitivity and energy budget |

For a quasi-static analytical reference, mechanical settling must be fast
relative to drainage/loading at the times used for comparison. For a
coseismic reference, inertia is part of the target response; forcing static
force balance would remove that response. Internal stress release and
external traction changes can both generate waves. Their amplitude,
timescale, initial stress balance and boundary conditions determine what
must be resolved. Existing PT fixtures using S=1e11 do not establish that
ordinary explicit/GVS applications need that value.

Kinetic/strain-energy ratios are useful alongside force and observable
errors, but a small ratio alone is not sufficient, especially near zero
initial strain energy. The [Abaqus mass-scaling reference](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEANLRefMap/simaanl-c-massscaling.htm)
likewise distinguishes quasi-static and dynamic uses and calls for checking
the inertial effect of modified mass. Its mass multiplier is not DES's S.
For changing GVS masses, distinguish physical and fictitious kinetic energy
and account for mass-change work before claiming a discrete energy balance.

## Fixed-mesh GVS implementation check

The old fixed-mesh path could refresh dt at a slow update without rebuilding
the velocity-dependent mass. A new dry sudden-traction check independently
assembles nodal mass and verifies `delta_v=dt*force/mass` on free components
at every step. It also checks that the stored next-step velocity scale
matches the accepted velocity. The old code fails at step 2 with an 84%
velocity-increment discrepancy. The fixed path refreshes mass and timestep
together after each accepted physical step, as the moving-mesh path does.

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 tests/functional/check_gvs_mass.py ./dynearthsol2d /tmp/gvs-mass-2d
python3 tests/functional/check_gvs_mass.py ./dynearthsol3d /tmp/gvs-mass-3d
```

Both ceilings pass 64 steps in both dimensions, with maximum relative
velocity-identity error `7.19e-15`. Velocity scales span more than seven
orders of magnitude. Every saved continuation frame matches the fresh
coordinate, velocity, stress, strain, temperature and physical clock
exactly after a step-16 restart. Both dimensions build and parent dry core
`make set` / changed `make cmp` remain bit-exact. Evidence:
`results/gvs-mass-validation-20261107/final-verification.json`.

These checks qualify mass/dt synchronization, not full coupled transient
accuracy. The method-comparison runner records those analytical errors
separately.

## Sudden-traction Terzaghi pilot, 2026-11-07

The runner `benchmarks/poromechanics/mechanical_methods.py` retains 43
exploratory 2D runs. All use the same 1 m by 10 m confined column, 100 kPa
sudden traction, K=20 MPa, G=10 MPa, zero gravity, fixed mesh and no heat.
The reference drainage time H^2/c is 2,933.74 s. Inputs, executable/output
hashes, error histories, force assembly and energy diagnostics are retained.

- The 20-case 400 s matrix includes both PT options and ordinary explicit /
  GVS shear/bulk at all three scalings and both initial velocity scales.
  Cold starts at 1e-12 m/s can finish with almost no load response. Other
  cases grow large errors, fail or stop before the requested output window.
  None of those outcomes establishes that GVS intrinsically needs PT.
- At h=1 m, reducing `dt_fraction` from 0.1 to 0.01 allows all six seeded
  GVS cases to finish the 20 s window. At S=1e6, pressure/settlement errors
  approach the coarse PT errors.
- At h=0.5 m, fraction 0.01 is inaccurate. Fraction 0.005 gives the results
  below at S=1e6 and a 1e-6 m/s initial scale. The first two-million-step
  trials stop before 20 s; fresh runs with a four-million-step cap complete
  it. A step-cap exit is recorded separately from nonfinite failure.

| Method, h=0.5 m | Max pressure/load error, saved times 2–20 s | Max pressure/load error, 10–20 s | Max settlement/drained error, 10–20 s |
| --- | --- | --- | --- |
| PT0 / PT1, physical dt=0.125 s | 9.057% | 2.000% | 0.286% |
| GVS shear, S=1e6, fraction=0.005 | 10.027% | 2.233% | 0.379% |
| GVS bulk, S=1e6, fraction=0.005 | 10.026% | 2.234% | 0.378% |

The final GVS maximum free-node force / total applied load is 1.57e-4
(shear) or 3.67e-4 (bulk); PT is about 4.77e-8. These are saved-stress
diagnostics after objective rotation, not the PT solver's acceptance norm.
The physical-reference kinetic/strain-energy ratios are small in these
final states, but that test alone would be misleading: an inaccurate
ordinary explicit case at h=0.5 m and fraction 0.01 has a ratio below 0.5%
while its pressure/load error exceeds 11.

This is positive evidence that non-PT GVS can carry this sudden traction
over the tested early window. It does not qualify the full consolidation
history, initial sub-second transient, separate spatial/time convergence,
3D coupled GVS, Mandel/Cryer, tectonic applications or earthquake cycles.
The successful fraction is not a universal setting. In particular, the
dry elastic/Gershgorin limits and stable hydraulic substeps do not by
themselves establish stability of the lagged coupled mechanical-pressure
update. Derive/check the undrained wave and coupled timestep limits next.

Results under `benchmarks/poromechanics/results/`:

- `mechanical-methods-2d-20261107/`: full matrix, including failed runs.
- `mechanical-methods-short-dt-2d-20261107/`: all scalings, smaller step.
- `mechanical-methods-h05-2d-20261107/`: finer mesh, larger fraction.
- `mechanical-methods-h05-halfdt-2d-20261107/`: initial step-cap trials.
- `mechanical-methods-h05-complete-2d-20261107/`: completed fine GVS cases.

Each directory has an HTML gallery and PNG/SVG/PDF figures. Reanalysis
reproduces every numerical metric and rejects altered inputs or raw saves.


## Coupled explicit correction, 2026-11-08

The extended 400 s test at h=0.5 m and an initial velocity scale of
1e-5 m/s finds nonfinite failure in all nine legacy coupled explicit/GVS
cases at S=1e4, 1e5 and 1e6, despite stable hydraulic substeps. Dry and
pressure-feedback-disabled controls remain finite. At the same fixed
mechanical dt=0.001 s over 40 s, replacing lagged fixed-stress pressure by
the strain-driven equation reduces maximum pressure/load error from 613%
to 9.29%, and final saved-stress force/load imbalance from 0.364 to
0.000752. The old kinetic/strain-energy ratio is only 5.08e-6: a small
kinetic ratio alone would incorrectly suggest equilibrium.

The candidate uses the existing strain-driven equation before stress and
limits the undrained wave using the actual scaled mass. All nine runs
finish at automatic fraction 0.5, including the GVS mass ceiling. PT's
algorithm, the mass formula and decoupled flow are unchanged.

| Method | S | Late pressure/load error | Late settlement/drained error | Final saved-stress force/total load |
| --- | --- | --- | --- | --- |
| PT0 / PT1, dt=1 s | 1e5 | 0.153% | 0.0881% | 7.65e-9 |
| Ordinary explicit | 1e4 | 4.78% | 1.21% | 8.86e-3 |
| Ordinary explicit | 1e5 | 0.234% | 0.0588% | 1.15e-5 |
| Ordinary explicit | 1e6 | 0.208% | 0.0339% | 1.00e-8 |
| GVS shear | 1e4 | 0.209% | 0.0854% | 3.39e-4 |
| GVS shear | 1e5 | 0.207% | 0.0331% | 4.17e-7 |
| GVS shear | 1e6 | 0.206% | 0.0329% | 4.17e-9 |

Late means saved times at or after 0.05 H^2/c (about 147 s); references
use each actual saved time. For GVS at S>=1e5 the maximum pressure/load
error over the entire saved 40–400 s window is about 0.85%. Bulk-ceiling
results approach the same errors at these scalings. The weaker S=1e4
ordinary result remains in the report rather than being discarded.

This is a finite-window comparison, not full mesh/time convergence or a
universal choice of initial speed, scaling or damping. Separate cold-start,
mesh/time and Mandel/Cryer checks remain necessary. The proposed correction
also changes coupled non-PT Richards steps from fixed-stress to strain
storage; fully wet laws must remain identical to the saturated model, and
unsaturated fluid content and restart are checked independently. Older
nonzero-step coupled quasistatic checkpoints require their archived solver.

Evidence: `results/runs/gvs-qualification-20261108/`, including the legacy
failures, diagnostic controls, matched-step comparison, source patch,
prototype matrix and final qualification. These exploratory results stay
outside the selected benchmark figures until the refinement work is done.

The final implementation passes both builds and C++20 syntax checks;
serial/OpenMP assembly passes 7,355 assertions in 2D and 10,871 in 3D.
The 84 unchanged-path checks preserve 260 files / 7,930 physical payloads
exactly, including PT wet/unsaturated and fixed/remeshed restart. Thirty
sealed explicit/GVS cases check exact wet-law payloads, physical restart
and nodal fluid content (maximum residual 2.86e-16). Twelve executable
probes independently verify the undrained timestep below and above both
GVS speed ceilings. The final 400 s shear-GVS replay at S=1e5 preserves
all 616 prototype payloads exactly. Nonzero-step legacy ordering is
rejected rather than silently converted. This is CPU evidence; no new GPU
qualification is claimed.
