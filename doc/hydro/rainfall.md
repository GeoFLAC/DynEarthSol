# Application: rainfall triggers a shallow slide

The second slope application: rain falls on a partly
saturated soil layer over impermeable bedrock, a perched water table rises
on the bedrock, and the slope fails. It uses the verified unsaturated flow
(Srivastava–Yeh), the Bishop coupling (drained column, sealed cell,
draining column) and effective-stress plasticity (confined layer) together.
There is no exact solution.

## Model

The slope of [slope.md](slope.md): 2 m of soil on 20° bedrock, 40 m long,
strong 4 m ends, gravity applied as a load by restart from an unstressed,
gravity-free checkpoint. The soil is weaker than on that page, $c=1$ kPa and
$\phi=30^\circ$ (no dilation), so that it fails before the layer saturates.
Gardner retention ($\alpha=2.4\times10^{-4}$ Pa$^{-1}$, $S_r=0.05$), Bishop
$\chi=S_e$, initial uniform suction 5 kPa ($S_e=0.30$, Bishop stress
$-1.5$ kPa). Rain $q=10^{-5}$ m/s (36 mm/h) enters through the surface; the
bedrock is sealed; the toe drains ($p=0$); the crest end is sealed.

DES has no runoff or ponding: rain that meets a saturated surface would be
forced in and build an unphysical overpressure. A first design with the
stronger soil of [slope.md](slope.md) did exactly that before failing
(surface pressure +11 kPa at 37,500 s), so the soil was weakened instead of
the rain being changed; the runner records the surface pressure at every
step, and the page refuses a run in which it turns positive.

## Comparison

Failure is read where the mean plastic strain of the basal 0.6 m (middle
third) reaches 0.001, as on [slope.md](slope.md): that curve converges with
mesh and step, while the last PT equilibrium does not. The same soil loaded
by bedrock pressure (gravity applied elastically first) reaches that slip at
11.38, 11.49 and 11.41 kPa on the 0.5, 0.25 and 0.125 m meshes. Infinite-slope estimates
for it: pointwise onset 8.8 kPa, non-associated slip plane 11.9 kPa,
classical factor of safety 15.5 kPa.

## Results

| Mesh | Onset of basal yield | Plastic strain 0.001 | Last PT equilibrium |
| --- | --- | --- | --- |
| 0.5 m | 6.7 h, 9.05 kPa | 8.72 h, 12.32 kPa | 9.86 h, 15.77 kPa |
| 0.25 m | 6.5 h, 8.94 kPa | 8.71 h, 12.39 kPa | 9.31 h, 13.52 kPa |

The surface never saturates (largest surface pressure 0). The time to
failure agrees to 0.1% between the meshes; the base pressure there is
0.9 kPa (8%) above the bedrock-loaded value.

The rain does not arrive as a front: with $k=10^{-11}$ m$^2$ and
$S_e\approx0.3$ it spreads through the 2 m column within hours, and the
layer fills from the bedrock up at a nearly uniform rate along the slope.
Basal yield starts close to the pointwise estimate, but the slope needs a
higher bedrock pressure to slip: the perched table leaves suction, and
Bishop strength, above it and drives seepage parallel to the slope, while
pressure from the bedrock drives an upward seepage force through the
whole layer.

## Current-source coupling audit

On a rerun of the 0.5 m input with DES's fluid source for plastic yielding,
basal plastic strain 0.001 occurs at 31,407.6 s (8.72 h), with mean basal
pressure 12.323 kPa, matching the table to round-off. Saved surface
pressures stay at or below zero throughout. The bedrock-pressure loading
on the same 0.5 m mesh reaches the criterion at 11.379 kPa, so rainfall
needs 0.944 kPa more mean basal pressure. A 0.25 m rerun has not been
completed. The saturated water budget of the pressure-driven slides closes
to round-off ([slope.md](slope.md)); a general unsaturated mass budget is
not checked here.

## Reproduce

```
python3 benchmarks/poromechanics/slope_application.py EXE OUT --spacings 0.5 0.25 0.125 --cohesion 1000 --friction 30 --elastic-gravity
python3 benchmarks/poromechanics/rainfall_slope.py EXE OUT --spacings 0.5 0.25
```

The review’s rainfall page ships the 0.25 m input, its seed checkpoint and
the study's numbers (`record.json`).
