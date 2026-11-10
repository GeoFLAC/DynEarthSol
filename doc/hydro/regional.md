# Application: seasonal pumping and land subsidence

A well field pumps a soft confined aquifer for six months, rests for six
months and repeats once. The ground subsides around the wells, partly
rebounds in each rest season and keeps a slowly growing subsidence farther
out. This applies the qualified pumping and thin-slab (averaged-aquifer)
physics at field scale, with a realistic schedule and rebound.

![Setup](../../benchmarks/poromechanics/results/review/fig/regional/setup.svg)

This is a preliminary application study. The spatial/time limits below
remain part of its acceptance record.

## Model

DES models a quarter of a 6 km square slab, $B=50$ m thick, with the well
field on its corner edge pumping $Q/4$, $Q=0.02\ \mathrm{m^3/s}$. Symmetry
rollers hold $x=0$ and $y=0$, the base is a roller and the top is
traction-free; the outer sides hold $p=0$ and zero normal displacement.
Gravity is off, so all fields are changes from the initial state. The
sediment has $K=5$ MPa, $G=3$ MPa, $\alpha=0.9$, $\phi=0.3$ and
$\kappa=10^{-12}\ \mathrm{m^2}$; the consolidation coefficient is
$7.4\times10^{-3}\ \mathrm{m^2/s}$, so the pressure front reaches 1.4 km in
two years. PT uses adaptive dynamic relaxation (`PT_option = 1`), 104 weekly
steps; the default accelerated scheme stalls near $10^{-2}$ of the initial
residual on a domain this size. Structured meshes are used: the unstructured
mesher fails on a slab 60 times wider than thick.

## Reference: a closed basin

The averaged (generalised Bear–Corapcioglu) model gives the areal strain as
the local pressure plus a uniform function of time,

$$
e_h=\frac{\alpha p}{2(\lambda+G)}+h(t),\qquad
\varepsilon_{zz}=\frac{\alpha p-\lambda e_h}{\lambda+2G},\qquad
\delta=-B\varepsilon_{zz}.
$$

In a closed basin the fixed outer walls make $h(t)$ nonzero, and through the
poroelastic coupling it deepens the drawdown everywhere. The infinite-domain
line well sets $h=0$ and departs from the DES basin by up to about 20% in the
rest seasons; the finite-annulus solution of
`averaged_aquifer_reference.py` (fixed outer pressure and displacement,
Talbot inversion) is used instead, on the circle of the same area
($R=2L/\sqrt\pi=3.39$ km), superposed over the rate changes.

## Results

Over $250\le r\le1500$ m, 17 frames to two years, against the closed-basin
reference (largest error over the band, relative to that frame's largest
value, after the first frame; worst absolute error at the observation wells
over all frames; peak subsidence 24.1 cm):

| Run | Pressure | Subsidence | Well 250 m | 500 m | 1000 m |
| --- | --- | --- | --- | --- | --- |
| 100 m, 104 steps | 6.3% | 20.5% | 1.57 cm | 0.52 cm | 0.17 cm |
| 50 m, 104 steps | 3.8% | 18.3% | 1.46 cm | 0.50 cm | 0.08 cm |
| 25 m, 104 steps | 3.7% | 13.0% | 1.41 cm | 0.42 cm | 0.09 cm |
| 50 m, 416 steps | 0.9% | 12.1% | 0.41 cm | 0.23 cm | 0.01 cm |
| 25 m, 416 steps | 0.4% | 7.0% | 0.32 cm | 0.10 cm | 0.01 cm |

Away from rate changes the pressure agrees within 0.2-1%. The largest
errors sit in the frames just after each pump-on and shut-in, where weekly
steps are too coarse for the sequential coupling: quartering the step on
the 50 m mesh cuts the worst well error from 1.46 to 0.41 cm. The first
frame, at six weeks, is also under-resolved in space (31.6% of its small
subsidence on the 25 m mesh at 104 steps, 12.8% at 416). On the 25 m mesh
with 416 steps the wells agree within 0.32 cm of a 24 cm peak, and that run
supplies the review page. Its runner's final binary-hash check failed only
because the bare binary was rebuilt while it ran; the process kept its
original image, the same binary as the other two studies, and identical
copies are archived (see `final-h25-s416/NOTE.txt`).

## Reproduce

```
python3 benchmarks/poromechanics/regional_subsidence.py ./dynearthsol3d OUT --spacings 100 50 25
python3 benchmarks/poromechanics/regional_subsidence.py ./dynearthsol3d OUT --spacings 50 --steps 416
python3 benchmarks/poromechanics/regional_subsidence.py ./dynearthsol3d OUT --spacings 25 --steps 416
```

Recorded studies: `results/runs/regional-subsidence-20261108/{study,steps416,final-h25-s416}/`.
