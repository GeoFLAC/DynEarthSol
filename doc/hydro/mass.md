# Finite-density water and sensible heat

`control.hydraulic_content_model=1` selects the saturated finite-density
inventory model. The default `0` keeps the linear Biot/retention model and
its existing arithmetic. This is an optional constitutive extension, with
separate qualification; it is not a change to the default model's EOS.
It requires `control.is_quasi_static=true`: scaled-mass ordinary explicit
mechanics, GVS or PT. Unscaled physical inertia is excluded from this mode.

The explicit/GVS pressure-ordering correction passes independent free-mechanics
checks in 2D and 3D on CPU and H100. Native and HDF5 CPU cell, transport,
restart and remap checks pass, as do clean parent-core comparisons and
C++17/C++20 CPU builds. This evidence covers the declared saturated model
and loading paths; it does not qualify arbitrary finite-strain mechanics.

The selected controls are a sealed compression cell, a drained heated cell
and a warm-water mixing cell. Their equations and results are in the worked
examples. Run their independent checks with:

```bash
python3 tests/functional/hydromechanics/check_mass.py EXE OUT
```

Use a new output folder. The executable determines the dimension and output
container. The suite includes physical-step and transport-substep refinement,
restart with unchanged and changed timesteps, remeshing, and accepted plastic
heat. Spatial transport is checked on the accepted mechanical strain path
against an independent Radau reference; that check does not independently
solve mechanical equilibrium.

The additional free-mechanics check uses a small-amplitude, undamped
undrained standing wave with declared fictitious density:

```bash
python3 benchmarks/poromechanics/undrained_wave.py EXE OUT \
  --content-model 1 --scaling 1e5 --amplitude 7.938948919340288e-9 --spacing .025 \
  --steps 1024 2048 4096 --samples 16 --cycles 4 --restart-check
```

Use one initial displacement amplitude, $u_0=1$ nm, across inertial scalings.
The sinusoidal initial velocity is $A=\pi c u_0$, where $c$ is the declared
undrained wave speed. The tested pairs are:

| Inertial scaling | Initial velocity amplitude (m/s) |
| --- | --- |
| `1e4` | `7.938948919340288e-10` |
| `1e5` | `7.938948919340288e-9` |
| `1e6` | `7.938948919340288e-8` |

Each pair passes ordinary explicit and constant-floor GVS in both dimensions
on CPU and H100; the middle pair also passes with `--moving-mesh`. These
28 studies retain 1024/2048/4096-step refinement, four periods, unchanged
error criteria and exact restart checks. A fixed velocity at all scalings
instead shrinks the strain signal by 100 between the endpoints; the original
`1e-8` velocity test at scaling `1e6` reaches the inventory roundoff floor
and misses its temporal-convergence criterion. Its failed record is retained.

`--moving-mesh` checks the geometry predictor. `--gvs` keeps the wave below
a declared constant velocity floor to retain an analytical reference.
Adaptive GVS has a separate momentum
identity check, `check_gvs_mass.py EXE OUT --content-model 1`, under
`tests/functional/`; add `--moving-mesh` to check original-mass ownership
after geometry prediction. Both geometries pass in 2D/3D on CPU and H100;
the measured momentum identity closes within $6.6\times10^{-15}$ and their
restarts preserve the physical fields exactly.

## Constitutive equations

For constant fluid bulk modulus and volumetric thermal expansivity,

$$
d\rho_f/\rho_f=dp/K_f-\alpha_f\,dT.
$$

The reference density is at zero pressure and the saved temperature origin:

$$
\rho_f(p,T)=\rho_{f0}\exp[p/K_f-\alpha_f(T-T_0)].
$$

This is a constant-compressibility EOS, not a tabulated water-property model.
The Eulerian Biot porosity law, with constant material coefficients, is

$$
d\phi=(\alpha-\phi)(d\varepsilon_v+dp/K_s-\alpha_s\,dT-dh).
$$

Here strain is tension positive; positive recovery $h$ closes accepted pore
dilation. Its integrated form is anchored at the physical step's initial
pressure, temperature and porosity. With a homogeneous reference state,

$$
\phi=\alpha-(\alpha-\phi_0)
\exp[-\varepsilon_v-p/K_s+\alpha_s(T-T_0)+h].
$$

The fluid inventory is $M=V\rho_f\phi$. A moving mesh supplies its actual
volume. A fixed mesh uses the virtual volume $V_0\exp(\mathrm{tr}\,\varepsilon)$;
this convention does not establish finite-strain mechanical accuracy.
The state must retain positive density and $0<\phi\leq\alpha\leq1$.

## Explicit time integration and local state conversion

Mechanics keeps its existing ordinary explicit, GVS or PT option. Transport
uses forward updates of extensive inventories. At each transport substep,
Darcy mass flux, upwind sensible-heat flux and conduction use the **current
substep's starting state**. Each interior face has one signed flux, added to
one node and subtracted from the other. Scheduled sources are integrated over
the actual substep interval.

For an interior node, the update has the form

$$
M^{j+1}+L_n(p^{j+1}-p_n)=A^j+\delta t_j F_M(p^j,T^j),
$$

$$
E^{j+1}=E^j+\delta t_j F_E(p^j,T^j).
$$

$A^j$ includes the initial physical mass, the released previous numerical
buffer and cumulative real mass flux. $L_n$ stays frozen during the physical
coupling step. Local scalar pressure/temperature roots convert these updated
inventories to a thermodynamic state. They introduce no solve for neighboring
nodes' future pressure or temperature. Local nonlinear state conversion is
separate from implicit time integration of spatial transport.

Pressure-controlled nodes instead impose their boundary pressure. Their
water reaction and matching $c_fT$ sensible-heat exchange are measured
explicitly. Fixed-temperature nodes have a separate energy reaction.

The substep is restricted by the current hydraulic row bound, outgoing mass,
thermal capacity/transport bound and deformation increment. Invalid local
states cause a shorter substep before any state is committed. The controls
check temporal convergence and positive admissible states; they do not prove
unconditional stability, arbitrary nonlinear material laws or exact buoyant
hydrostatic balance.

## Physical mass and the fixed-stress lag

PT retains an extensive numerical pressure stabilizer:

$$
L_n=\sum_{e\ni n}\frac{\rho_eV_e}{N_e}\frac{\alpha_e^2}{B_e},
\qquad B_{2D}=K+G/3,\quad B_{3D}=K.
$$

At a free-pressure node it leaves $R=L_n(p_{n+1}-p_n)$. The physical-step
ledger closes **$M+R$**, not physical $M$ alone. In a sealed cell the physical
mass defect is $-R$; the pressure and mass defects decrease at first order
with the physical timestep in the reported compression study. $R$ contributes
neither weight nor heat capacity. Resetting it at a pressure boundary is a
numerical reset, not a water or heat flux.

Ordinary explicit and GVS use $L=0$, so their ledger closes physical $M$.
They apply the current transport pressure increment to total stress before
the accepted velocity update. This isotropic shift leaves effective stress
unchanged; it does not repeat the plastic return. The undrained wave tangent
$K+4G/3+\alpha^2/S$ restricts the mechanical timestep, using the current
minimum nodal storage $S$ in each element.

A fixed mesh needs one force evaluation. A moving mesh predicts its geometry
with the existing mechanical step, accepts transport on that geometry, and
recomputes force from the current pressure. Its accepted velocity starts
from the original velocity, mass and physical interval; the predictor is
not a second momentum increment. Nodal mixed stress is applied once after
the pressure shift. Remapping and checkpoints distinguish a raw pressure
increment already consumed by stress from one still pending. The new
checkpoint version rejects earlier finite-density snapshots that lack this
ownership information; default linear-model checkpoints keep their schema.

The physical interval and accepted shear energy are captured before GVS
selects the next timestep. Transport substeps and PT trials never redeposit
that accepted energy. The integrated energy target is carried through
substeps without accumulating the local root solver's roundoff residuals.

## Sensible heat and remapping

Solid mass and solid heat capacity $C_s$ are material inventories. Draining
water does not recreate either from the current porosity. Retained energy is

$$
E=(C_s+c_fM)T.
$$

Inlet water uses its injection temperature and EOS density; extraction uses
local temperature and density. The heat ledger distinguishes supplied heat,
well enthalpy, carried boundary heat, fixed-temperature reaction and retained
energy change. Only actual fluid mass contributes to heat capacity.

A shared geometric-overlap remap carries $M$, $R$, solid mass, $C_s$, $E$ and
the mechanical temperature reference. Accepted pending shear energy uses the
same overlap geometry. Removed and imported material have separate inventory
budgets. Pressure/temperature constraints account for their reservoir
reactions. Checkpoints preserve these inventories, the EOS reference and
pending coupling increments; changing the timestep does not reinterpret them.

This is a **sensible-heat model**. Full pressure work, thermoelastic caloric
energy, phase conversion and recovery energy are not included. The selected
zero-pressure rigid open-cell controls have no mechanical pressure work;
they do not by themselves qualify a general thermal plume or earthquake
energy balance. Unsaturated retention, matrix exchange, surface processes
and phase conversion are excluded from this optional mode.
