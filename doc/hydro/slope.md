# Application: pore pressure triggers a translational slide

The first slope application (status item 6): a soil layer on inclined
bedrock fails as pore pressure rises from the bedrock. There is no exact
solution, so DES is compared with three engineering estimates, and the
gaps between them are part of the result.

## Model

A layer of vertical thickness $D=2$ m on bedrock inclined at
$\beta=20^\circ$, 40 m long with vertical ends, in plane strain. The base
is fixed (no slip) and supplies water; the surface is free and drained
($p=0$); the ends are rollers and sealed. The 4 m at each end are the same
rock with cohesion 10 MPa, so the corners stay elastic. Without them the
corners yield first under gravity and PT stalls at 0.25 m; the
comparison uses only the middle third. The rest is Mohr–Coulomb with
$c=2$ kPa, $\phi=35^\circ$, no dilation, $K=50$ MPa, $G=5.17$ MPa
($\nu=0.45$), saturated, $\rho_s=2650$ kg/m$^3$, $\phi_n=0.3$.

Gravity is applied as a load: a zero-step run without gravity or stress
writes a checkpoint, and the slope restarts from it with gravity on, so
its first step finds the gravitational equilibrium of the real geometry.
The legacy initial stress is isotropic and measured from $z=0$, which on a
slope is far from equilibrium. A $\nu=0.45$ material is needed because an
elastic slope at rest has slope-parallel stress $K_0\sigma_n$,
$K_0=\nu/(1-\nu)$: with $\nu=0.25$ at $30^\circ$ that state already
violates Mohr–Coulomb throughout. The base pressure then rises linearly
over 2000 s; permeability $10^{-11}$ m$^2$ keeps the seepage quasi-steady.

## Estimates

Statics of an infinite slope give, on the bedrock plane,
$\sigma_n=\gamma D\cos^2\beta=37.3$ kPa and
$\tau=\gamma D\sin\beta\cos\beta=13.6$ kPa ($\gamma$ saturated unit weight).

- **Onset of yield** (pointwise Mohr–Coulomb on the poroelastic path): a
  pore-pressure rise $u$ lowers the slope-normal effective stress by $u$ and
  the slope-parallel one by $K_0u$, while $\tau$ stays; yield at
  $u=14.1$ kPa.
- **Non-associated slip plane** (Davis's reduction, as in Drescher and Detournay 1993): a zero-dilation shear band
  fails as if $\tan\phi^*=\sin\phi$, $c^*=c\cos\phi$; $u=16.5$ kPa.
- **Classical infinite slope**, $FS=[c+(\sigma_n-u)\tan\phi]/\tau=1$:
  $u=20.8$ kPa.

## Results

"Basal plastic strain" is the mean plastic strain of the elements in the
basal 0.6 m of the middle third.

| Mesh | Steps | Onset of basal yield | Plastic strain 0.001 | Plastic strain 0.002 | Last PT equilibrium |
| --- | --- | --- | --- | --- | --- |
| 0.5 m | 404 | 12.94 kPa | 16.47 kPa | 17.22 kPa | 18.08 kPa |
| 0.25 m | 404 | 13.56 kPa | 16.54 kPa | 17.19 kPa | 18.00 kPa |
| 0.125 m | 404 | 13.87 kPa | 16.42 kPa | 16.95 kPa | 17.30 kPa |
| 0.25 m | 808 | 13.52 kPa | 16.51 kPa | 17.16 kPa | 18.55 kPa |

Yield starts at the base and converges at first order toward the
pointwise estimate (14.1 kPa). Above it the basal slip grows along one
curve on every mesh and both step counts: the base pressure at a basal
plastic strain of 0.001 varies by 0.7%, and it lies at Davis's
non-associated estimate. The slip then accelerates, and PT stops finding
an equilibrium at 17.3–18.6 kPa, still 11–17% below the classical factor
of safety. That end point is not converged: it rises with the step count
and falls on the finest mesh. It records where those PT runs cease to
satisfy their equilibrium criterion; it does not identify a physical
collapse pressure. The comparison therefore uses a specified basal slip
for this layer with strong ends. Zero-dilation plasticity makes the slope
weaker than limit equilibrium assumes, consistent with Davis's correction.

Limits: no dynamic run-out; an associated-flow (dilating) variant became
non-finite in an earlier design and was not repeated; 2D only.

A weaker soil ($c=1$ kPa, $\phi=30^\circ$, gravity applied elastically
first), used by the rainfall application, reaches a basal plastic strain of
0.001 at 11.38, 11.49 and 11.41 kPa on the 0.5, 0.25 and 0.125 m meshes, against
estimates of 8.8 (pointwise), 11.9 (non-associated) and 15.5 kPa
(classical). Its basal yield starts at 7.1, 7.9 and 8.3 kPa, approaching
the pointwise estimate. Rainfall onto the slope is [rainfall.md](rainfall.md).

## Fluid content during plastic yielding

While the base yields, the water stored in the slope must still match the
water that crossed its boundaries. An independent volume integral of the
saturated Biot content, plus the pressure increment PT has not yet passed
to the mechanics, closes that budget to about $1.7\times10^{-17}$ m³ per
metre on both soils (0.25 m mesh, through 1090 s at $\phi=35^\circ$ and
835 s at $\phi=30^\circ$). See the [content derivation](plastic-content.md)
for the pending term. On these runs the base pressure at a basal plastic
strain of 0.001 is 16.539 kPa ($c=2$ kPa, $\phi=35^\circ$) and 11.486 kPa
($c=1$ kPa, $\phi=30^\circ$), as in the table above.

PT stops after the fixed-slip criterion, at 16.680 and 12.705 kPa on these
runs. Like the last PT equilibrium in the table, this is where the solver
stops, not a collapse load.

The independent checker reads a run's saved frames and checkpoints, for
example the 35° slope through 1090 s:

```sh
python3 tests/functional/hydromechanics/check_content.py OUT/slope.cfg --until 1090
```

The checker requires every accepted step, constant properties, saturated
flow and fixed geometry. It includes the initial pending increment at an
elastic-gravity restart and does not claim an unsaturated or remapped
content law.

## Reproduce

```
python3 benchmarks/poromechanics/slope_application.py EXE OUT --spacings 0.5 0.25 0.125
python3 benchmarks/poromechanics/slope_application.py EXE OUT --spacings 0.25 --steps 808
```

The review’s slope page ships the 0.125 m input, its seed checkpoint and
the study's numbers (`record.json`).
Reference: Drescher and Detournay (1993), Limit load in translational
failure mechanisms for associative and non-associative materials,
Géotechnique 43(3), 443–456, which uses Davis's reduced strength
parameters and shows that the non-associative (coaxial) limit load lies
below the associative one.
