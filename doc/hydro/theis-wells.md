# Several wells: analytical superposition and recovery

Two wells extract water and one injects it into a confined, homogeneous
slab. All stop at 40 s; pressure then recovers until 80 s. This is a
constant-property transport check with mechanical pressure feedback
disabled. It complements the existing coupled small-strain superposition
test, `check_multiwell.py`, which remains unchanged.

![Setup](../../benchmarks/poromechanics/results/review/fig/theis-wells/setup.svg)

The full slab is 160 by 160 by 5 m, with all faces sealed. Wells penetrate
the full thickness. With extraction positive and $Q=10^{-4}$ m³/s:

| Well | $(x,y)$, m | Rate until 40 s |
| --- | --- | --- |
| 1 | $(60,80)$ | $Q$ |
| 2 | $(100,80)$ | $0.6Q$ |
| 3 | $(80,100)$ | $-0.25Q$ |

The selected reduced flow model has linear storage $S_B$ and diffusivity
$D=1$ m²/s. Its pressure transmissivity and storativity are

$$
\mathcal T_p=\kappa B/\mu,\qquad S_p=B S_B,\qquad D=\mathcal T_p/S_p.
$$

## Analytical solution and balance

Linearity of the pressure diffusion equation permits spatial and temporal
superposition of the [Theis solution](https://doi.org/10.1029/TR016i002p00519).
Let $r_i$ be distance from well $i$, and define $W(r,t)=0$ for $t\le0$ and

$$
W(r,t)=E_1\!\left(\frac{r^2}{4Dt}\right),\qquad t>0.
$$

Pressure drawdown after pumping and shut-in is therefore

$$
s_p(\mathbf x,t)=\sum_i\frac{Q_i}{4\pi\mathcal T_p}
 [W(r_i,t)-W(r_i,t-40)].
$$

The reference is for an infinite aquifer. At the finite slab's far sides
its value remains below $10^{-5}$ of the stated pressure scale during the
comparison. Observations exclude every well axis: distance to the nearest
well is 10–30 m. Error is normalized by
$C_\Sigma=\sum_i|Q_i|/(4\pi\mathcal T_p)$, without fitting properties.

Integration over the sealed mesh independently gives

$$
S_B\int_V[p(t)-p(0)]\,dV=-\sum_i Q_i\min(t,40).
$$

## Results

![Result](../../benchmarks/poromechanics/results/review/fig/theis-wells/result.svg)

The 5 and 2.5 m meshes at a 1 s outer step give maximum pressure errors
of 2.547% and 0.639% of $C_\Sigma$. The 2.5 m mesh with a 0.125 s outer
step gives 0.688%. The slight increase exposes partial mesh/time error
cancellation; the smaller error is not selected as a fitted result.
Hydraulic subcycling has its own step bound.

Separate single-well runs reproduce the combined pressure to
$1.14\times10^{-15}C_\Sigma$. The largest independently integrated
volume error among all six runs is $1.13\times10^{-17}$ m³. These checks
establish constant-property transport and its signed source budget, not
nonlinear or evolving-property superposition.

## Reproduce

Run from the repository, using a fresh output directory:

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
transport_solver=benchmarks/poromechanics/results/runs/gvs-qualification-20261108/final-regression/binaries/dynearthsol3d
transport_manifest=benchmarks/poromechanics/results/runs/gvs-qualification-20261108/qualification.json
python3 tests/functional/hydromechanics/check_theis_wells.py \
  "$transport_solver" benchmarks/poromechanics/results/runs/transport-qualification-20261108/replay-multiwell \
  --binary-manifest "$transport_manifest"
```

The recorded study and output hashes are in
[summary.json](../../benchmarks/poromechanics/results/runs/transport-qualification-20261108/multiwell/summary.json).
It binds the executable to the frozen solver qualification for `c5fc173`.
The existing coupled small-strain test remains under
`check_multiwell.py`; this analytical check is `check_theis_wells.py`.
