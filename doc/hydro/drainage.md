# Coupled drainage of a partly saturating column (U9)

A saturated sand column drains through its base: the top desaturates, the
pore pressure falls below zero and the column settles under the Bishop
effective stress and its lighter weight. This is the Liakopoulos (1965)
experiment in shape. It is the first check of the Bishop coupling with flow
during drainage; the drained column (U8) and the sealed cell test the end
states and the transient coupling without flow.

## Problem

A laterally confined column 1 m high is saturated and hydrostatic with its
water table at the top. At $t=0^+$ its base is opened to $p=0$ (a 1 ms ramp
of the scheduled boundary pressure); the top is sealed and traction-free,
the base fixed. Material after Liakopoulos: $E=1.3$ MPa, $\nu=0.4$,
$\alpha=1$, $\phi=0.2975$, $\kappa=4.5\times10^{-13}\ \mathrm{m^2}$, and a
van Genuchten law ($\alpha_{vG}=1/3000$ Pa$^{-1}$, $n=1.8$, $S_r=0.1$).
Gravity acts; pressure feedback is on.

## Reference

There is no closed form. `drainage_reference.py` solves the same continuum
independently: Richards flow with $\chi=S_e$, uniaxial strain with Bishop
stress and the weight change, which makes the storage nonlocal (drying
above a point unloads it),

$$
c(p)\,p_t+\alpha S_e\varepsilon_t=-q_z,\qquad
\varepsilon_t=\frac{\sigma_t+\alpha\Pi'(p)\,p_t}{K+\tfrac43G},\qquad
\sigma_t(z)=-g\rho_f\phi(1-S_r)\int_z^H\frac{dS_e}{dp}p_t\,dz',
$$

with vertex-centred finite volumes, the dense coupled storage matrix and an
implicit BDF integrator. Grid-to-grid pressure differences fall at second
order (5.2 → 1.3 Pa at 300 s, of a 9.8 kPa scale); the comparison uses 200
cells. Over two hours the column settles 2.07 mm and the top desaturates
to $S=0.79$.

## Results

Largest nodal pressure error relative to $\rho_f gH$, and settlement error,
at 300, 1200, 3600 and 7200 s (2D, PT, 480 steps unless stated):

| Run | Pressure | Settlement |
| --- | --- | --- |
| 0.05 m | 0.58 / 0.34 / 0.22 / 0.16% | 0.24 / 0.06 / 0.03 / 0.01% |
| 0.025 m | 0.17 / 0.09 / 0.06 / 0.04% | 0.24 / 0.05 / 0.02 / 0.01% |
| 0.0125 m | 0.14 / 0.02 / 0.02 / 0.01% | 0.24 / 0.05 / 0.02 / 0.01% |
| 0.0125 m, 960 steps | 0.07 / 0.02 / 0.02 / 0.01% | 0.12 / 0.03 / 0.01 / 0.01% |

Late-time pressure converges at second order with the mesh. The first
frame and the settlement are limited by the time step: the base opens
instantly and the sequential coupling lags one step, and doubling the
steps halves both (first order).

The explicit strain form (dynamic relaxation, inertial scaling $10^5$,
characteristic speed $10^{-5}$ m/s) on the 0.05 m mesh gives pressure
errors 0.55 / 0.33 / 0.21 / 0.16%, comparable to PT, and settlement within
0.0071% over the four compared times. Its small stable step greatly
reduces the coupling-lag error. An independent replay of both recorded
PT/explicit cases reproduces all eight pressure/settlement comparisons;
the 54 input/output hashes and recomputed reference are recorded in
`drainage-20261108/explicit-2d/independent-replay.json`. A first explicit
run left `control.characteristic_speed` at 0
with no boundary velocity, so the skeleton never moved; it is kept with a
note (`mesh-2d/NOTE-explicit.txt`) and excluded.

The 3D sanity check at `4ec082f` uses 0.05 m spacing and 480 PT steps.
At the same four times its pressure errors are 0.527 / 0.300 / 0.188 / 0.135%
of $\rho_f gH$, and settlement errors are 0.190 / 0.059 / 0.029 / 0.014%.
These agree with the 2D sequence at the same spacing; this single 3D mesh
does not replace the 2D convergence study.

## Reproduce

```
python3 benchmarks/poromechanics/drainage_reference.py
python3 tests/functional/hydromechanics/check_unsaturated_drainage.py EXE OUT --spacings 0.05 0.025 0.0125
python3 tests/functional/hydromechanics/check_unsaturated_drainage.py EXE OUT --spacings 0.05 --explicit
python3 tests/functional/hydromechanics/check_unsaturated_drainage.py EXE OUT --spacings 0.0125 --steps 960
python3 tests/functional/hydromechanics/check_unsaturated_drainage.py EXE3D OUT --ndims 3 --spacings 0.05
```

Recorded studies (frozen `e77de6e` binary):
`results/runs/drainage-20261108/{mesh-2d,steps960-2d,explicit-2d}/`.
