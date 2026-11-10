# Consolidation and recovery after partial unloading

A confined elastic column first carries a compressive load, then loses
half of it. The material settles while it drains and recovers after the
load reduction. Signed superposition of the Terzaghi step solution gives
both pressure and settlement without fitting any parameters.

![Loaded and drained column](../../benchmarks/poromechanics/results/runs/terzaghi-regimes-20261109/fig/setup.svg)

The height is 1 m. The base is sealed and fixed vertically; sides are
sealed with normal-motion rollers. The top holds zero pore pressure and
carries uniform total traction $-F(t)$. Gravity is zero and the mesh stays
fixed. All cases use an elastic skeleton and the saturated linear Biot law.
The 2D column has width $2h$, keeping two cells across the confined width;
the one-dimensional continuum solution is independent of this width.

## Load history and analytical solution

Let depth $d$ increase from the drained top. Define

$$
H_m=K+\frac43G,\qquad
S=\frac{\phi}{K_f}+\frac{\alpha-\phi}{K_s},
$$

$$
C=S+\frac{\alpha^2}{H_m},\qquad
D=\frac{k}{\mu C},\qquad T_v=\frac{Dt}{H^2}.
$$

The pressure response per unit compressive traction is

$$
\mathcal P(d,t)=\frac{\alpha}{H_m C}
\sum_{n=0}^{\infty}\frac{4}{(2n+1)\pi}
\sin(\lambda_n d)e^{-D\lambda_n^2t},
$$

$$
\lambda_n=\frac{(2n+1)\pi}{2H}.
$$

The sine basis enforces $p(0,t)=0$ and $p_d(H,t)=0$. Its initial interior
pressure is the undrained Biot response. Integrating pressure gives the
mean response and settlement per unit traction:

$$
\overline{\mathcal P}(t)=\frac{\alpha}{H_m C}
\sum_{n=0}^{\infty}\frac{8e^{-D\lambda_n^2t}}{(2n+1)^2\pi^2},
$$

$$
\mathcal U(t)=\frac{H}{H_m}
[1-\alpha\overline{\mathcal P}(t)].
$$

For $t<0$ both response functions are zero. At $t=0^+$, pressure and
settlement take their undrained values. For the load reduction at
$t_j=0.5H^2/D$,

$$
p(d,t)=F_0[\mathcal P(d,t)-\tfrac12\mathcal P(d,t-t_j)],
$$

$$
\delta(t)=F_0[\mathcal U(t)-\tfrac12\mathcal U(t-t_j)].
$$

The final load is $F_0/2$; pressure tends to zero and settlement tends to
$F_0H/(2H_m)$. The single-step solution follows from the diffusion equation
$Cp_t=(k/\mu)p_{dd}$. Constant coefficients make its pressure and displacement
operators linear, which establishes the signed superposition above.
`TerzaghiReference.load_history` evaluates arbitrary signed load increments.
Nine reference checks include the unloading jump, independent pressure
quadrature for settlement, the diffusion equation and both drainage limits.

The initial DES frame is pre-load. The checkpoint at $T_v=0.5$ is also
before unloading; its comparison uses the earlier load history. The
restarted phase changes traction and compares subsequent accepted states
with the complete two-event solution. This distinguishes a saved initial
state from the analytical right-hand limit after a load change.

## Prescribed material regimes

All regimes have $K_s=1$ GPa and $\mu=10^{-3}$ Pa s.
The Biot coefficient obeys $\alpha=1-K/K_s$ in each case.

| Regime | $K,G$ (MPa) | $K_f$ (MPa) | $\phi$ | $\alpha$ | $F_0$ (Pa) | $k$ (m²) |
| --- | --- | --- | --- | --- | --- | --- |
| Compressible fluid | 20, 10 | 1 | 0.25 | 0.98 | 10,000 | $10^{-11}$ |
| Water | 20, 10 | 2200 | 0.25 | 0.98 | 10,000 | $10^{-11}$ |
| Soft skeleton | 0.2, 0.1 | 2200 | 0.40 | 0.9998 | 100 | $10^{-14}$ |

The characteristic drainage times $H^2/D$ are 27.9542 s, 2.96556 s and
299,958 s. Every case runs through $T_v=2.5$ with the same dimensionless
load history. Loads keep strains small; these are prescribed physical
regimes, not changes made to fit a curve.

## Selected comparison and refinement

![Pressure and recovery in the soft-skeleton case](../../benchmarks/poromechanics/results/runs/terzaghi-regimes-20261109/fig/result.svg)

The selected soft-skeleton result uses $h=0.05$ m and 4000 steps per unit
$T_v$. Lines are the analytical series; circles are DES. Negative excess
pressure after unloading and the recovery toward half the original drained
settlement are part of the linear reference.

| Regime | Maximum pressure error / $F_0$ | Maximum settlement error / $F_0H/H_m$ | Temporal order: pressure, settlement |
| --- | --- | --- | --- |
| Compressible fluid | 0.01344% | 0.00312% | 1.006, 1.005 |
| Water | 0.16838% | 0.05713% | 1.026, 1.015 |
| Soft skeleton | 0.17138% | 0.05945% | 1.027, 1.015 |

The spatial sequence uses 0.2, 0.1 and 0.05 m at 4000 steps per unit
$T_v$. Pressure and settlement errors decrease on each refinement.
Time comparisons hold the 0.1 m mesh fixed and use differences between
successive full pressure fields and settlement histories. Absolute errors
can cancel between spatial and temporal contributions.

The compressible-fluid run at 1000 steps per unit $T_v$ takes two hydraulic
substeps. The 2000/4000-step runs take one; the mixed sequence has pressure
order 0.62. An independently reconstructed FE row bound is 0.0243242 s.
The 2000/4000/8000 sequence stays within that bound with one substep and
recovers first-order time convergence. Its earlier coarse results remain
in the evidence. The other regimes use the 1000/2000/4000 sequence.

Twenty retained cases provide 1020 selected observations, with all raw
hashes and analytical comparisons replayed. The matching native/HDF5
control preserves 416 field arrays exactly. Tiny 3D columns (45 nodes,
80 tetrahedra) execute all three regimes, but their 0.25 m mesh has up to
4.34% pressure/load error and 2.56% settlement error: these are dimensional
controls, not refined 3D accuracy qualification. That remains open.

These tests use PT option 1, fixed material properties and fixed geometry.
They extend the load/material coverage of consolidation; they do not
qualify evolving storage, conservative remapping, plasticity or the full
explicit/GVS comparison matrix. Each saved time is compared on its actual
mesh; no phase or parameter is fitted.

## Reproduce

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
solver=benchmarks/poromechanics/results/runs/plastic-source-20261108/candidate/dynearthsol2d
python3 benchmarks/poromechanics/terzaghi_pulse.py "$solver" /new/terzaghi-pulse
python3 benchmarks/poromechanics/terzaghi_pulse.py "$solver" /new/terzaghi-pulse --replay
python3 benchmarks/poromechanics/terzaghi_pulse.py "$solver" /new/terzaghi-time \
  --regimes soft-fluid --spacings .1 --counts 8000 --time-spacing .1
```

The [retained evidence](../../benchmarks/poromechanics/results/runs/terzaghi-regimes-20261109/)
includes the independent convergence analysis, exact inputs, source/binary
hashes, every comparison table, schematics and the selected vector figures.
