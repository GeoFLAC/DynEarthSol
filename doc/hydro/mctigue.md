# Thermal pressurization: McTigue surface heating

The drained, traction-free surface of a fluid-saturated poroelastic
half-space is raised by $\Delta T$ at $t=0^+$. Pore water expands more than
the grains, so pressure rises below the surface while draining upward, and
the solid heaves. This checks the thermal pressurization source, the thermal
stress and their two-way coupling with drainage and deformation.

![Setup](../../benchmarks/poromechanics/results/review/fig/mctigue/setup.svg)

DES models a laterally confined column 1 m deep. Its sides and base are
sealed, insulated and fixed normally; its top is drained ($p=0$),
traction-free and held at $T_0+\Delta T$ with $\Delta T=10$ K. Heat
conducts without Darcy advection; gravity is off; PT keeps equilibrium.
The material has $K=G=1$ GPa, $\alpha=0.8$, $\phi=0.2$, $K_f=2.2$ GPa,
$K_s=5$ GPa, $k=10^{-19}\ \mathrm{m^2}$, $\mu=10^{-3}$ Pa s, volumetric
expansions $\alpha_s=3\times10^{-5}$ and $\alpha_f=2.07\times10^{-4}$ K$^{-1}$,
and thermal diffusivity $\kappa=10^{-6}\ \mathrm{m^2/s}$.

## Analytical solution

Uniaxial strain with a traction-free surface keeps the vertical total
stress zero (tension positive, $H=K+4G/3$):

$$
\varepsilon=\frac{\alpha p+K\alpha_sT}{H}.
$$

With fluid content $\zeta=\alpha\varepsilon+S_0p-\beta T$,
$S_0=\phi/K_f+(\alpha-\phi)/K_s$ and $\beta=\phi\alpha_f+(\alpha-\phi)\alpha_s$
(the DES law in [transport.md](transport.md)), Darcy flow gives

$$
Sp_t-\frac{k}{\mu}p_{zz}=S\gamma T_t,\qquad
S=S_0+\frac{\alpha^2}{H},\quad
\gamma=\frac{\beta-\alpha K\alpha_s/H}{S},\quad
c=\frac{k}{\mu S}.
$$

For $T=\Delta T\,\mathrm{erfc}(z/2\sqrt{\kappa t})$ (McTigue, 1986),

$$
p=\frac{\gamma\Delta T}{1-c/\kappa}\left[\mathrm{erfc}\frac{z}{2\sqrt{\kappa t}}-\mathrm{erfc}\frac{z}{2\sqrt{ct}}\right],
$$

and integrating the strain gives the surface heave. Here $\gamma=101$ kPa/K
and $c/\kappa=0.21$. `benchmarks/poromechanics/mctigue_reference.py` checks
the equation, the boundary and initial conditions, the heave integral and an
independent Crank–Nicolson solution (agreement within $2\times10^{-4}$).

## The step must be a load

DES applies a fixed boundary temperature at initialization. A surface that
is hot in frame 0 is therefore part of the initial state: its nodes carry no
thermal stress and no thermal pressurization. Run that way, the top element
layer misses its load and the errors stop converging (0.025 m mesh:
pressure 10%, strain 24%, heave 10%). The test instead writes an initial
checkpoint at $T_0$ and restarts with the hot surface, so the step is applied
at $t=0^+$ as in the reference (same mesh: 3.9%, 5.6%, 2.3%). Applications
that heat a boundary suddenly need the same care.

## Results

Largest errors over eight saved times to 40 ks; pressure relative to its
peak, temperature to $\Delta T$, strain to its largest value:

| Run | Pressure | Temperature | Strain | Heave |
| --- | --- | --- | --- | --- |
| 2D 0.05 m, 200 steps | 15.8% | 1.04% | 19.4% | 5.32% |
| 2D 0.025 m, 200 steps | 3.86% | 0.27% | 5.57% | 2.35% |
| 2D 0.0125 m, 200 steps | 0.92% | 0.067% | 2.10% | 1.90% |
| 2D 0.00625 m, 200 steps | 0.99% | 0.041% | 2.19% | 1.81% |
| 2D 0.00625 m, 400 steps | 0.46% | 0.041% | 1.03% | 0.91% |
| 3D 0.025 m, 200 steps | 5.82% | 0.22% | 6.48% | 2.92% |
| 3D 0.0125 m, 200 steps | 1.68% | 0.054% | 2.54% | 2.33% |

The errors fall at about second order until the mesh reaches the time-step
floor near 0.0125 m. Doubling the steps then halves the pressure, strain and
heave errors: the one-step lag of the sequential pressure-mechanics
coupling. Temperature, which no lag affects, reaches 0.04%.

## Heat ledger

The heated surface puts in exactly the heat stored in the column:
$\rho c_p\int_0^L(T-T_0)\,dz$ from the image series of the insulated
column (indistinguishable from the half-space value
$2\rho c_p\Delta T\sqrt{\kappa t/\pi}$ here, $3\times10^{-9}$). DES's
stored heat, from its lumped nodal volumes, exceeds it at 5 ks by 1.505,
0.366, 0.091 and 0.022% on the 0.05, 0.025, 0.0125 and 0.00625 m meshes:
second order, unchanged by 400 steps, and smaller at later times because
the surplus is a fixed amount near the heated node. This closes the heat
part of the energy accounting; pressure work and mechanical energy are
not included.

## Reproduce

```
python3 benchmarks/poromechanics/mctigue_reference.py
python3 tests/functional/hydromechanics/check_mctigue.py ./dynearthsol2d OUT --spacings 0.05 0.025 0.0125 0.00625
python3 tests/functional/hydromechanics/check_mctigue.py ./dynearthsol2d OUT --spacings 0.00625 --steps 400
python3 tests/functional/hydromechanics/check_mctigue.py ./dynearthsol3d OUT --ndims 3 --spacings 0.025 0.0125
```

Recorded studies with input and executable hashes:
`results/runs/mctigue-20261108/{mesh-2d,steps400-2d,mesh-3d}/`. The finest
2D runs take 2 (200 steps) and 3 minutes (400 steps) on one thread.
