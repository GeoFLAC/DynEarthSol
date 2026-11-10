# Injection temperature of fluid point sources

A fluid point source (`injection.*`) adds volume to the pressure equation.
Before this option it carried no heat: hot-water injection needed a separate
heat source, which neither follows the injected rate nor depends on the
temperature it meets. `injection.temperature` gives the temperature of the
injected fluid of each point.

## Equation

The default linear-content model (`control.hydraulic_content_model=0`)
integrates heat in the advective form

$$
C\,T_t+\rho_fc_f\,\mathbf q\cdot\nabla T=\nabla\cdot(\lambda\nabla T)+\rho_fc_f\,Q\,(T_{in}-T),
$$

with $C$ the bulk heat capacity and $Q$ the injected volume rate per unit
volume. The conservative form, $\partial_t(CT)+\nabla\cdot(\rho_fc_f\mathbf
qT)=\dots+\rho_fc_fQT_{in}$, gives this once fluid mass conservation,
$\nabla\cdot\mathbf q=Q$, removes $\rho_fc_fT\,\nabla\cdot\mathbf q$: the
source term is the injected enthalpy relative to the local temperature.
Extraction ($Q<0$) removes fluid at the local temperature, which the
advective term already carries, so extracting points add nothing.

Per node, $\rho_fc_f(\sum_pQ_{pn}T_{in,p}-T_n\sum_pQ_{pn})$ over injecting
points, with the point-source shape weights and $\rho_fc_f$ of the reference
material. The term relaxes the node toward $T_{in}$ at rate
$\rho_fc_fQ_n/C_n$, which enters the explicit heat-step bound.
`assemble_point_sources` takes an optional per-point scale, so the same
mapping gives $\sum Q_{pn}$ (scale 1 on injecting points) and
$\sum Q_{pn}T_{in,p}$ (scale $T_{in,p}$). Unset, nothing changes: core and
hydraulic regressions are bit-exact (2D/3D `make cmp`; 40 comparisons,
8378 payloads) and the point-source checks pass (23 cases in 2D, 22 in 3D).

## Assumptions and limits

- Heat capacity is DES's fixed $C$: the heat carried by fluid that goes into
  storage is not tracked. For the fluid stored by pressurization this is of
  relative order $\rho_fc_f\,\delta V/(C V)$, the injected fraction of the
  stored heat capacity, small within linear storage. A sealed cell that
  stores all of its injection, with no flow, mixes exponentially in DES and
  rationally in the physical ledger $(C V+\rho_fc_fQt)\,\theta$; this sealed
  mixing case is outside the scope.
- No pressure work, viscous dissipation or thermal expansion of the
  injected water; $c_f$ constant.
- Sequential splitting: temperature is updated before pressure, so thermal
  advection uses the Darcy flux of the pressure at the start of the step.
  A source switched on from rest has no advective flux in its first step.
  When that step needs several thermal substeps (the inflow term shortens
  the bound), the injection node saturates and part of the first step's
  heat is lost. The deficit is a fixed fraction of one step's heat, so its
  share of the total falls with the step and with time.

## Verification

**Radial hot-water injection (Avdonin’s problem).** A line source at $T_{in}$
in a layer at $T_0$ with no heat loss has the similarity solution
$\theta=\Gamma(\nu,u)/\Gamma(\nu)$, $u=r^2C/(4\lambda t)$,
$\nu=\rho_fc_fq/(4\pi\lambda)$ (`avdonin_reference.py`, which checks that it
satisfies the equation and stores all injected heat). DES runs a quarter
plane with the source at its corner, drained far sides and a stiff skeleton
($\nu=3.3$; the front reaches 2.3 m of an 8 m domain).

Largest nodal error at radius above 0.25 m, relative to $T_{in}-T_0$, at
1.4, 2.9, 4.3 and 5.8 days, and stored against injected heat:

| Mesh, steps | Largest error | rms (5.8 d) | Stored/injected − 1 |
| --- | --- | --- | --- |
| 0.4 m, 500 | 11.7 / 8.5 / 7.2 / 5.8% | 1.15% | 6e-6 |
| 0.2 m, 500 | 5.5 / 4.0 / 3.3 / 2.8% | 0.55% | 6e-6 |
| 0.1 m, 500 | 3.0 / 2.1 / 1.7 / 1.5% | 0.28% | 0.43 → 0.11% |
| 0.2 m, 200 | 5.8 / 4.1 / 3.4 / 2.9% | 0.55% | 0.69 → 0.17% |
| 0.2 m, 1000 | 5.7 / 4.1 / 3.4 / 2.9% | 0.56% | 7e-6 |

Temperature converges at first order with the mesh, as expected of the
clipped (upwind-like) advection; the largest error sits at the front. The
step barely matters. The heat balance is exact except where the first step
needs several thermal substeps (0.1 m at 500 steps, 0.2 m at 200 steps):
then the start deficit of about half a step's heat appears, and its share
falls as $1/t$. With one substep, forward Euler from $T_0$ puts in the full
heat and the balance holds to $10^{-5}$.

**Closed-box ledger, hot doublet.** A sealed, insulated box with an
injector at $T_{in}$ and a producer at the same rate. No heat or fluid
crosses the boundary, so by energy conservation alone
$E(t)-E(0)=\rho_fc_fq\int_0^t(T_{in}-T_{prod})\,ds$, with $T_{prod}$
interpolated at the producer. The hot front reaches the producer, so the
ledger includes outflow at a changing temperature.

Over 400 steps of 1500 s ($6\times10^5$ s; the producer reaches 77-78% of
$T_{in}-T_0$), stored heat against the ledger, relative to the injected
heat $\rho_fc_fq(T_{in}-T_0)t$:

| Mesh | Start deficit | Ledger after the first step, largest / final | Total, final |
| --- | --- | --- | --- |
| 0.2 m | 0.75 step | 2.4e-4 / 1.1e-4 | 1.8e-3 |
| 0.1 m | 0.92 step | 1.9e-4 / 2.7e-5 | 2.3e-3 |

After the first step the ledger closes to a few $10^{-4}$, including the
outflow at a changing producer temperature. The start deficit grows on
the finer mesh, whose first step needs more thermal substeps; it is a
fixed amount, so its share of the total falls as $1/t$. At 0.4 m it is
0.75 step with 200 steps and 0.53 with 400.

The active-rate substep bound for short pulses preserves these constant-rate
cases exactly. Separate reruns check Avdonin at 0.4 m and the doublet at
0.4 m with 200 steps; their recorded inputs and fields are unchanged.

## Reproduce

```
python3 benchmarks/poromechanics/avdonin_reference.py
python3 tests/functional/hydromechanics/check_injection_temperature.py EXE OUT --spacings 0.2 0.1
python3 tests/functional/hydromechanics/check_injection_temperature.py EXE OUT --spacings 0.4
python3 tests/functional/hydromechanics/check_injection_temperature.py EXE OUT --spacings 0.2 --steps 200   # and 1000
python3 tests/functional/hydromechanics/check_injection_energy.py EXE OUT --spacings 0.2 0.1 --steps 400
```

Recorded studies: `results/runs/injection-temperature-20261108/`.

## Completed 3D slab

The 0.2 m slab completes at source `4ec082f`. Maximum temperature errors
at 125,000, 250,000, 375,000 and 500,000 s are 5.851, 4.064, 3.240 and
2.767% of the imposed temperature rise. Stored heat is lower than the
reference by a fixed 71.259 kJ startup amount: 0.2728% of supplied heat
at the first time and 0.06819% at the last. This is the default fixed-capacity
model and the startup limitation described above.

The optional [finite-density inventory model](mass.md) instead carries
physical water mass and retained sensible heat. Its separate drained-heating
and mixing controls do not retrospectively qualify the default model or a
general spatial thermal plume.
