# Physical-inertia hydraulic steps

The closed-column wave verifies the saturated elastic response when
`control.is_quasi_static=false` and `control.has_PT=false`. Pressure and
velocity agree with the independent standing-wave solution within 0.3%
of their amplitudes over one period on the finest 3D grid. The selected
figure belongs in the [guided review](index.md);
this page records the formulation and its limits.

## Why the pressure update depends on the mechanical formulation

For tension-positive strain and compression-positive pressure, the small-strain
Biot fluid content is

$$
\zeta=\alpha\epsilon_v+S_0p-\beta\Delta T,\qquad
S_0=M^{-1}=\frac{\phi}{K_f}+\frac{\alpha-\phi}{K_s},\qquad
\beta=\phi\alpha_f+(\alpha-\phi)\alpha_s.
$$

When the accepted strain rate is available, the fluid equation uses pressure
capacity $S_0$ and source $-\alpha\dot\epsilon_v+\beta\dot T$ directly.
The volumetric rate is `edvoldt`, the projected rate also used by the existing
anti-locking constitutive update. Using the unprojected velocity divergence
would give pressure and recorded strain different fluid contents.

PT and the reduced quasistatic flow model retain the fixed-stress split. Eliminating
volumetric strain with the recorded mean total stress gives capacity
$S_B=S_0+\alpha^2/K_b$, where $K_b=K+G/3$ in 2D and $K_b=K$ in 3D.
Its mean-stress source and thermal correction remain unchanged. That split
is tested separately against consolidation; stable Darcy substeps alone do
not establish its coupled wave accuracy.

The coupled explicit branch advances temperature, computes the projected
strain rate, advances pressure, then evaluates stress, force and velocity.
Its new pressure increment contributes to the current stress update.
PT and reduced-flow steps keep their ordering and consume the pressure
increment in the next mechanical update. Both store raw
`dppressure = -Delta p`; hydraulic substeps integrate the source once over
the full physical interval.

| Mechanical formulation | Pressure capacity | Pressure source and ordering |
| --- | --- | --- |
| Physical inertia, no PT | $S_0$ | Projected strain source when enabled; pressure before stress |
| Coupled explicit, scaled mass or GVS | $S_0$ | Projected strain; pressure before stress |
| PT | $S_B$ | Mean total stress; pending increment for the next mechanical update |
| Quasistatic, pressure feedback disabled | $S_B$ | Reduced flow model; pending increment |

The capacity selection is independent of retention law. In the saturated
range, laws 1–3 must reproduce model 0 exactly. In an unsaturated material
the conservative content is $\phi S(p)+s p$, with $s$ selected from this
table. Changing the numerical ordering also changes the meaning of a saved
pending pressure increment; see the checkpoint restrictions below.

## Wave speed and timestep

In a homogeneous sealed column with zero permeability and no thermal change,
$\zeta=0$ implies $p=-\alpha\epsilon_v/S_0$. The constrained longitudinal
modulus and physical mixture density are

$$
H_u=K+\frac43G+\frac{\alpha^2}{S_0},\qquad
\rho_b=(1-\phi)\rho_s+\phi\rho_f,\qquad c_u=\sqrt{H_u/\rho_b}.
$$

Automatic physical-inertia timesteps include the geometric limit
$h_{\min}/(5c_u)$, multiplied by `control.dt_fraction`, alongside the other
limits. This limit is exercised through the actual executable in both
dimensions. It is not a general spectral proof for every coupled mesh,
constitutive law or boundary condition. `control.fixed_dt` retains its
explicit override; hydraulic subcycling cannot stabilize an oversized
mechanical step.

For scaled explicit mass, the same undrained modulus is combined with the
actual fictitious density $K/(S v)^2$, including the existing GVS speed
ceiling. The automatic limit uses $h\sqrt{\rho_{\mathrm{fict}}/H_u}/5$.
The heterogeneous implementation conservatively combines the minimum
unclipped bound with the minimum mass-floor bound. It introduces no new
scaling parameter. See [the sudden-load comparison](mechanics.md) for the
measured effect and remaining accuracy limits.

For a column of length $L$, the independently prescribed initial velocity
$v(z,0)=V\sin(\pi(z+L)/L)$, zero initial pressure/strain and fixed ends give

$$
v(z,t)=V\sin(\pi(z+L)/L)\cos(2\pi t/T),\qquad
p(z,t)=-\frac{\alpha V}{S_0c_u}\cos(\pi(z+L)/L)\sin(2\pi t/T),\qquad
T=2L/c_u.
$$

The runner initializes this velocity through a zero-step binary checkpoint.
It fits no phase, amplitude or material parameter. All-node errors use the
prescribed amplitudes; conservation compares $S_0p$ with the volume-weighted
nodal projection of the recorded volumetric strain.

## Qualification record

Evidence is in `results/runs/physical-pressure-fix-20261107/`.
`wave-qualification.json` records raw-data hashes and replayed metrics.

- Mesh spacings 0.05, 0.025 and 0.0125 m at 4096 steps per period give
  approximately second-order pressure convergence in both dimensions.
- On the fixed 0.0125 m grid, 1024/2048/4096/8192 steps give velocity
  self-convergence orders 1.001–1.003 and pressure orders 2.000 for this
  standing mode. The complete method has first-order time accuracy;
  the pressure result is specific to this problem and sampling.
- The selected 3D, 729-node, 8192-step result has maximum velocity error
  0.171% and pressure error 0.262%, across all nodes and 65 saved times.
  Maximum local fluid-content defect is 1.06e-14 of the strain amplitude.
- Physical fields and clocks continue exactly from a quarter-period
  checkpoint in all eight fine-grid runs. The test does not claim identical
  legacy output-scheduler metadata. Incompatible pressure ordering is rejected.
- Both dimensions build. Parent dry `make set` / changed `make cmp` remain
  bit-exact. The regression pass has 94 successful checks and zero changes
  in 360 PT/quasistatic/GVS save/checkpoint files (10,714 physical payloads),
  including active unsaturated flow and fixed/remeshed restart.
- Direct production assembly tests pass 7,217 assertions in 2D and 10,649
  in 3D, including physical strain/heating sources and forced subcycling.

The earlier lagged fixed-stress physical-inertia implementation fails the
closed-column content and phase checks on the tested 0.1 m grid. Refinement
to 4096 steps per period still leaves large errors. Historical failures are
preserved as research evidence, outside the selected-result gallery.

## Checkpoint scope and limits

Strain-form explicit checkpoints carry `physical hydraulic step=1`, distinguishing
an already consumed increment from the fixed-stress pending increment.
A new solver rejects mismatched ordering with exit code 22. Conversion is
allowed only at step zero with an exactly zero pressure increment.
The allowed conversion reselects the timestep and rebuilds masses for the
chosen formulation, instead of retaining an obsolete saved wave limit.
Continue nonzero-step checkpoints from older physical-inertia or coupled
quasistatic solvers with their
matching archived solver; they are not silently converted.

This qualification covers a small-strain, saturated, fixed-mesh elastic wave
over one period. GVS still uses its existing fictitious mass and solid-density
speed ceiling; reaching that ceiling does not establish correct saturated
coseismic inertia. Coupled GVS scaling/damping, complete energy accounting,
earthquake-cycle transitions, full consolidation histories and nonlinear/
unsaturated dynamics lie outside this wave comparison.

Reproduce a fine-grid study, using a new output directory:

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 benchmarks/poromechanics/undrained_wave.py ./dynearthsol3d /new/wave-3d \
  --spacing .0125 --steps 1024 2048 4096 8192 --samples 64 --restart-check
```

Run 2D and 3D sequentially. The output includes the input decks, initialized
checkpoint, runner snapshot, raw fields, restart checks and automatic-step
test. `qualify_wave.py --help` describes replay with separate
coarse/middle mesh studies at 4096 steps.

## Manual transition to physical inertia

A separate restart control now switches ordinary explicit or GVS mechanics
at a nonzero coupled state to physical inertia. Both phases already use
the strain-form hydraulic ordering; this does not convert a nonzero PT
pending-pressure history or implement an automatic earthquake trigger.

![Coupled inertia-switch control](../../benchmarks/poromechanics/results/runs/inertia-switch-20261108/fig/setup.svg)

The fixed-mesh unit column has a fixed base, side-normal rollers and a
constant −1 kPa top traction. Water boundaries are sealed; internal Darcy
redistribution remains active ($k=10^{-11}$ m², $\mu=10^{-3}$ Pa s).
Gravity, damping and thermal transport are zero. Material properties are
$K=20$ MPa, $G=10$ MPa, $\alpha=0.8$, $\phi=0.25$, $K_f=100$ MPa and
$K_s=1$ GPa. The warm phase takes 32 accepted steps at scaling
$10^4$, $10^5$ or $10^6$, with ordinary explicit, shear-ceiling GVS or
bulk-ceiling GVS. Its pressure and strain are already nonzero at the switch.

The physical phase sets `is_quasi_static=false`,
`use_global_velocity_scaling=false` and `mass_scaling_reference_speed=shear`.
The last setting restores the valid neutral value of the GVS-only option.
Automatic timestep selection and nodal mass must use the physical mixture:

$$
\rho_b=(1-\phi)\rho_s+\phi\rho_f=2237.5\ {\rm kg/m^3},
$$

$$
K_u=K+\frac43G+\frac{\alpha^2}{S},\qquad c_u=\sqrt{K_u/\rho_b}.
$$

The independent minimum-simplex-altitude calculation gives the new
hydraulic-wave limit $\Delta t=f h_{\min}/(5c_u)$, with $f=0.2$.
It controls these cases: **42.8985 µs in 2D** and **35.0265 µs in 3D**.
At the first accepted physical interval, every unconstrained velocity
component must satisfy

$$
v_{i,n+1}-v_{i,n}=\frac{\Delta t}{\rho_b V_i}F_{i,n+1}.
$$

This checks the mass actually used by mechanics, rather than reading a
configured mass or only inspecting the saved timestep. Across 36 cases in
2D/3D native/HDF5 output, its largest relative error is
$4.23\times10^{-16}$. All 3,618 checked saved-state payloads are unchanged
at the switch, excluding the intentionally reselected timestep and runtime
metadata. The closed global Biot-content budget stays within
$2.27\times10^{-16}$ of its stored/deformation scale. A further physical
restart preserves 24,840 continuation payloads exactly.

Changing artificial mass does not release physical heat. At fixed saved
velocity, the physical kinetic energy is continuous. The change in the
numerical kinetic-energy expression is reported separately:

$$
\Delta K_{\rm num}=\frac12\sum_i
(m_{i,\rm physical}-m_{i,\rm artificial})|v_i|^2.
$$

For example, at scaling $10^6$ the 2D shear-ceiling GVS mass uses an
effective density of 5300 kg/m³; bulk-ceiling GVS reaches 2650 kg/m³.
Neither equals the saturated mixture density. The switch preserves the
stress, strain, pressure and velocity fields and reselects the mass and
step; it does not rescale velocity or deposit the numerical mass-energy
change as heat. Full mechanical/thermal work, event-trigger hysteresis,
RSF coseismic paths and earthquake-cycle restart/remapping remain open.

Reproduce the matrix using the corresponding frozen candidate executable:

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
solver=benchmarks/poromechanics/results/runs/plastic-source-20261108/candidate/dynearthsol2d
python3 tests/functional/check_inertia_switch.py "$solver" /new/inertia-switch
python3 tests/functional/check_inertia_switch.py "$solver" /new/inertia-switch --replay
```

Repeat sequentially with the native/HDF5 2D/3D binaries. The inputs,
every-frame output, independent references and hashes are retained in
[`inertia-switch-20261108`](../../benchmarks/poromechanics/results/runs/inertia-switch-20261108/).
