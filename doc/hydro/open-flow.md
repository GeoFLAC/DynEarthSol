# Open boundaries: draining a smooth pressure pulse

A smooth pressure pulse drains through the two zero-pressure ends of a
10 m column. The lateral faces are sealed and the cross-section is 1 m².
Constant storage and disabled mechanical pressure feedback isolate the
transport ledger from constitutive coupling.

![Setup](../../benchmarks/poromechanics/results/review/fig/open-flow/setup.svg)

## Analytical pressure, stored fluid and outflow

For $p_t=Dp_{xx}$, $p(0,t)=p(L,t)=0$ and
$p(x,0)=P\sin(\pi x/L)$, separation of variables gives

$$
p(x,t)=P\sin(\pi x/L)e^{-\lambda t},\qquad
\lambda=D\pi^2/L^2.
$$

The test uses $P=1000$ Pa, $D=1$ m²/s, $L=10$ m and the reduced model's
constant pressure capacity $S_B$. Integration gives stored fluid

$$
M(t)=S_BA\int_0^Lp\,dx=M_0e^{-\lambda t},\qquad
M_0=2S_BPAL/\pi.
$$

Darcy flux at each end is outward. Its integrated volume is

$$
V_{\rm out}(t)=M_0[1-e^{-\lambda t}].
$$

Thus stored fluid plus cumulative outflow is constant. The reference
tests both quantities, rather than judging conservation from the saved
residual alone.

## Results and discrete accounting

![Result](../../benchmarks/poromechanics/results/review/fig/open-flow/result.svg)

All 200 physical steps are saved through 10 s on each 3D mesh:

| Spacing, m | Max pressure error / $P$ | Max stored-volume error / $M_0$ |
| --- | --- | --- |
| 1 | 0.6433% | 0.8162% |
| 0.5 | 0.1679% | 0.2038% |
| 0.25 | 0.04240% | 0.05093% |

Both errors approach second-order spatial convergence. The outer step is
0.05 s; stable hydraulic substeps are automatic. At every frame, a separate
tetrahedral-volume integration computes the change in lumped fluid content.
It is compared with the sum of all accepted boundary-volume increments.
The largest error is $1.20\times10^{-21}$ m³, or
$1.64\times10^{-15}M_0$.

The numerical initial integral differs from the continuum $M_0$ by its
spatial quadrature error. The exact discrete conservation check therefore
uses the measured initial mesh integral, while the accuracy table uses
the continuum reference. These are distinct checks.

The initial pulse is placed in a step-zero checkpoint with zero pending
pressure increment; the modified seed and checkpoint are hashed. This
does not establish conservation under remapping, varying storage or
unsaturated mechanical coupling. Those remain separate required checks.

## Reproduce

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
transport_solver=benchmarks/poromechanics/results/runs/gvs-qualification-20261108/final-regression/binaries/dynearthsol3d
transport_manifest=benchmarks/poromechanics/results/runs/gvs-qualification-20261108/qualification.json
python3 tests/functional/hydromechanics/check_open_flow.py \
  "$transport_solver" benchmarks/poromechanics/results/runs/transport-qualification-20261108/replay-open-flow \
  --binary-manifest "$transport_manifest"
```

The recorded study, input/seed hashes and all-frame budgets are in
[summary.json](../../benchmarks/poromechanics/results/runs/transport-qualification-20261108/open-flow/summary.json).
The executable is bound to the frozen solver qualification for `3774df0`.
