# Pressure-driven yielding of a confined layer

Rising pore pressure lowers the effective stress until the ground yields:
the mechanism behind injection-induced and rainfall-induced failure. This
checks DES's effective-stress Mohr–Coulomb plasticity in a boundary-value
problem with an exact solution: a laterally confined saturated layer whose
pore pressure rises.

## Problem

A layer 2 m deep is saturated, with its water table at the top, the DES
legacy isotropic lithostatic initial stress and hydrostatic pressure, and a
drained surcharge $q=20$ kPa on its top. Sides and base are fixed normally
and sealed. The top pressure rises linearly; with a high permeability the
whole layer follows. The material is elastic–perfectly plastic,
$K=100$ MPa, $G=60$ MPa, $c=10$ kPa, $\phi=20^\circ$, no dilation,
$\alpha=1$, $\phi_n=0.3$.

## Exact solution

Under uniaxial strain both horizontal stresses stay equal, and the total
vertical stress is fixed by the weight, so a uniform rise $\Delta p$ gives
(tension positive, $d$ depth, $\kappa=\nu/(1-\nu)$,
$\gamma'=(\rho_b-\alpha\rho_f)g$)

$$
\sigma'_v=-(\gamma'd+q)+\alpha\Delta p,\qquad
\sigma'_h=-(\gamma'd+\kappa q)+\kappa\alpha\Delta p .
$$

Yielding, $(\sigma_1-\sigma_3)/2+(\sigma_1+\sigma_3)/2\,\sin\phi=c\cos\phi$,
starts at the first root $\Delta p^*(d)$ of the elastic trial state. It
increases with depth, so a yield front moves down. Each depth is
independent under uniaxial strain, so the front stays exact after the
layer above has yielded. The yielded top reaches the apex of the
Mohr–Coulomb cone, mean effective tension $c\cot\phi$, at
$\Delta p=c\cot\phi+q$: beyond it no equilibrium exists and the layer would
lift off, which bounds the test at 95% of that pressure (45.1 kPa). The
surcharge keeps the top in compression; without it the apex follows the
onset within a few kPa. `pressure_yield_reference.py` checks the onset on
the yield surface, the elastic state just below it, the closed form of the
surcharge-free case and the ordering of onset and apex.

The infinite-slope version of this test needs a traction-free initial
stress on an inclined surface. The DES legacy initial stress is isotropic
and measured from $z=0$, so it does not provide one; that is a missing
initial-stress option, not a plasticity limit.

## Results

Each element's yield pressure is its own pressure rise at the first saved
step with plastic strain, which removes any diffusion lag; it is compared
with $\Delta p^*$ at its centroid depth. Onset is 36.7 kPa at the top and
the front reaches 1.2 m.

| Mesh | Steps | Yielded | Max error | Mean error | One increment |
| --- | --- | --- | --- | --- | --- |
| 0.2 m | 400 | 11/20 | 1.98% | 1.25% | 0.30% |
| 0.1 m | 400 | 23/40 | 1.33% | 0.90% | 0.31% |
| 0.05 m | 400 | 47/80 | 0.93% | 0.70% | 0.31% |
| 0.05 m | 1600 | 48/80 | 0.56% | 0.35% | 0.08% |
| 0.05 m | 6400 | 48/80 | 0.43% | 0.28% | 0.02% |
| 0.025 m | 1600 | 96/160 | 0.34% | 0.23% | 0.08% |
| 0.025 m | 6400 | 96/160 | 0.23% | 0.15% | 0.02% |

The 0.025 m runs use PT option 1: the default option went non-finite in the
first (surcharge) step on that mesh. On 0.05 m with 1600 steps the two
options agree (mean error 0.351% against 0.354%).

DES yields slightly late everywhere. Before yield its effective stresses
follow the exact uniaxial path within 0.05% (19 Pa of 40 kPa), with no
difference between the two triangle orientations, so the error is the
onset itself: detection on the first saved step past onset and the
one-step lag of the fixed-stress coupling, both of which shrink with the
step, plus a spatial part that shrinks with the mesh. Plotted as raw yield
pressures on a narrow axis the late bias looks large; the review figure
therefore shows the yield-front depth, the DES stress paths and the signed
error by step and mesh.

The complete 3D sanity check at `3774df0` (0.2 m, 400 steps) yields
28 of 50 elements. Maximum onset error is 1.9504%, mean error 1.1285%, and
the largest saved pressure increment is 0.3040% of the analytical onset.
The first onset is 37.207 kPa against 37.085 kPa at its element centroid.
This completes the loading ramp after the Mohr–Coulomb edge-return correction;
the earlier stalled run is excluded from the qualification. This is one
3D mesh, with the spatial and timestep convergence studies supplied by 2D.

Ordinary explicit and GVS runs on the same 3D mesh, at inertial scaling
`1e5`, also complete the ramp and yield 28 of 50 elements. Their maximum
onset errors are 1.65982% and 1.65981%, respectively; mean errors are
0.86363%. The first onset is 37.094 kPa in both, against 37.085 kPa in the
reference. This agreement is specific to this load, scaling and mesh.

## Reproduce

```
python3 benchmarks/poromechanics/pressure_yield_reference.py
python3 tests/functional/hydromechanics/check_pressure_yield.py EXE3D OUT --ndims 3 --spacings 0.2 --steps 400
python3 tests/functional/hydromechanics/check_pressure_yield.py EXE OUT --spacings 0.2 0.1 0.05 --steps 400
python3 tests/functional/hydromechanics/check_pressure_yield.py EXE OUT --spacings 0.05 --steps 1600
python3 tests/functional/hydromechanics/check_pressure_yield.py EXE OUT --spacings 0.05 0.025 --steps 1600 --pt-option 1
python3 tests/functional/hydromechanics/check_pressure_yield.py EXE OUT --spacings 0.025 --steps 6400 --pt-option 1
```

Recorded studies (frozen `e77de6e` binary):
`results/runs/pressure-yield-20261108/{mesh-2d,steps1600-2d,steps6400-2d,pt1-steps1600-2d,pt1-steps6400-2d,sanity-3d}/`.
