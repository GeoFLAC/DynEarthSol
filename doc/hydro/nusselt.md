# Steady hydrothermal convection: Nusselt number

The onset check (heat.md) verifies when convection starts. This check
verifies how much heat a steady roll carries once it runs: the Nusselt
number of the Horton–Rogers–Lapwood square cell at Ra = 100 and 200,
about 2.6 and 3.8 times the conductive heat flow.

## Problem

A unit square porous cell has impermeable walls, a base at $T_0+\Delta T$,
a top at $T_0$ and insulated sides. With $\theta=(T-T_0)/\Delta T$, lengths
scaled by $H$ and a streamfunction $\psi$, Boussinesq Darcy flow and heat
transport read

$$
\nabla^2\psi=-\mathrm{Ra}\,\theta_x,\qquad
u\theta_x+w\theta_z=\nabla^2\theta,\qquad u=\psi_z,\ w=-\psi_x,
$$

with $\mathrm{Ra}=\rho_f^2c_f g\alpha_f\Delta T kH/(\mu\lambda)$. The
steady heat balance at every height gives

$$
\mathrm{Nu}=1+\frac{\rho_fc_fH}{\lambda\Delta T}\langle q_zT\rangle ,
$$

the volume form used on DES fields, with the element Darcy flux rebuilt from
the saved pressure and temperature. A top-gradient estimate is reported
alongside.

## Reference

There is no closed form above onset.
`benchmarks/poromechanics/convection_steady.py` solves the equations above
independently: second-order central differences, relaxed Picard iteration,
three grids ($32^2$, $64^2$, $128^2$, observed order 1.99) and Richardson
extrapolation. It returns Nu = 1 below onset and balances top and base heat
flows. Nu = 2.6459 at Ra = 100 (published values are near 2.65) and 3.8098
at Ra = 200.

## DES runs

The cell of `check_convection.py` is seeded with a roll and run for 1.6 Ms
until Nu changes by less than $10^{-8}$. Storage does not enter the steady
state, so the pore fluid is made compressible ($K_f$ = 10 Pa): pressure then
diffuses about 20 times faster than heat instead of 2000, and the explicit
runs are 100 times cheaper. The transient is not compared.

| Mesh $h/H$ | Ra = 100 | wall | Ra = 200 | wall |
| --- | --- | --- | --- | --- |
| 0.05 | +4.30% | −7.82% | +17.89% | −13.37% |
| 0.025 | +0.26% | −1.87% | +2.43% | −3.40% |
| 0.0125 | +0.18% | −0.67% | +0.83% | −1.41% |

The volume form converges fast and the wall estimate from below, so the two
bracket the reference. At Ra = 100 the volume form flattens near +0.17%;
the Darcy flux rebuilt from saved fields is not DES's own discrete flux,
which is the likely floor. The 0.0125 runs take 12 and 21 minutes on one
thread.

## Reproduce

```
python3 benchmarks/poromechanics/convection_steady.py
python3 tests/functional/hydromechanics/check_convection_nusselt.py ./dynearthsol2d OUT --spacings 0.05 0.025 0.0125
```

Recorded study: `results/runs/convection-nusselt-20261108/study/` (see its
NOTE.txt for the executable).
