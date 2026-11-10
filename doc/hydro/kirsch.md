# Elastic circular opening: the Kirsch limit

A dry elastic solid with a circular opening receives isotropic external
compression. This is the first tunnel component check: stress, elastic
strain, symmetry constraints and normal traction on curved facets. It
does not yet qualify anisotropic far loading, plasticity or poroelastic
tunnel injection.

![Setup](../../benchmarks/poromechanics/results/review/fig/kirsch/setup.svg)

DES meshes a quarter annulus with opening radius $a=1$ m and outer radius
$R=20a$. Both straight cuts have zero normal displacement; the inner arc
is free and the outer arc carries compressive normal stress $P=100$ kPa.
Gravity and transport are off. The material has $K=G=1$ GPa in plane
strain. PT supplies static equilibrium; it is not a claim that this load
requires PT.

## Analytical solution

Radial equilibrium and isotropic Hooke's law give displacement
$u_r=C r+D/r$. With tension-positive stress, define

$$
A_0=-\frac{PR^2}{R^2-a^2},\qquad
B_0=-\frac{Pa^2R^2}{R^2-a^2}.
$$

The exact finite-annulus solution is

$$
\sigma_{rr}=A_0-B_0/r^2,\qquad
\sigma_{\theta\theta}=A_0+B_0/r^2,\qquad
\sigma_{r\theta}=0.
$$

It satisfies $d\sigma_{rr}/dr+(\sigma_{rr}-\sigma_{\theta\theta})/r=0$,
$\sigma_{rr}(a)=0$ and $\sigma_{rr}(R)=-P$. For
$\lambda=K-2G/3$, plane-strain displacement and strains are

$$
u_r=\frac{A_0r}{2(\lambda+G)}+\frac{B_0}{2Gr},\qquad
\epsilon_{rr}=\frac{du_r}{dr},\quad
\epsilon_{\theta\theta}=\frac{u_r}{r}.
$$

The out-of-plane stress is
$\sigma_{yy}=\nu(\sigma_{rr}+\sigma_{\theta\theta})$.
These strains are measured from an unloaded annulus. An excavation
displacement relative to an already compressed intact solid needs that
background deformation subtracted; the excavation history is not tested here.

As $R\to\infty$, the solution approaches the isotropic
[Kirsch solution](https://docs.itascacg.com/itasca950/3dec/docproject/source/examples/CylindricalHoleInAnIEM.html):

$$
\sigma_{rr}=-P(1-a^2/r^2),\qquad
\sigma_{\theta\theta}=-P(1+a^2/r^2).
$$

The finite model's wall hoop concentration is $2R^2/(R^2-a^2)$ rather
than exactly 2. Its 0.251% truncation difference is kept separate from
mesh error by comparing DES with the finite-annulus reference.

## Results

![Result](../../benchmarks/poromechanics/results/review/fig/kirsch/result.svg)

Errors use all element centroids in $1.1a\le r\le5a$. Curved boundaries
are faceted; local maximum errors need not decrease monotonically on
independently generated meshes.

| Inner spacing, m | Radial stress RMS / $P$ | Hoop stress RMS / $P$ | Largest local stress error / $P$ |
| --- | --- | --- | --- |
| 0.2 | 0.911% | 0.835% | 3.12% |
| 0.1 | 0.561% | 0.511% | 2.36% |
| 0.05 | 0.299% | 0.342% | 2.59% |
| 0.025 | 0.148% | 0.169% | 2.15% |

The finest radial/hoop strain RMS errors are 0.0719% / 0.0802% of $P/G$.
Out-of-plane stress RMS error is 0.0168% of $P$. A separate NumPy force
assembly from stress, simplex gradients and outer-edge tractions gives
maximum free-node imbalance $1.57\times10^{-8}PR$. The outer resultant
also matches $(-PR,+PR)$ per unit out-of-plane thickness.

## Reproduce

From the repository, using a fresh output directory:

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 benchmarks/poromechanics/kirsch_reference.py
application_solver=benchmarks/poromechanics/results/runs/gvs-qualification-20261108/final-regression/binaries/dynearthsol2d
application_manifest=benchmarks/poromechanics/results/runs/gvs-qualification-20261108/qualification.json
application_run=benchmarks/poromechanics/results/runs/application-qualification-20261108/replay-kirsch
python3 tests/functional/hydromechanics/check_kirsch.py \
  "$application_solver" "$application_run" \
  --spacings .2 .1 .05 .025 --binary-manifest "$application_manifest"
python3 benchmarks/poromechanics/qualify_kirsch.py \
  "$application_run" "$application_run/review"
```

The [qualification record](../../benchmarks/poromechanics/results/runs/application-qualification-20261108/kirsch/qualification.json)
binds the source, immutable solver, polygon inputs and field hashes.
The finest refinement lives below `kirsch/finer/`; all earlier meshes and
their original runners remain available. Salençon and poroelastic injection
are the remaining tunnel stages.
