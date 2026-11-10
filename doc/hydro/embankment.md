# Factor of safety of an embankment by strength reduction

The standard slope-stability problem: a homogeneous embankment whose
factor of safety DES finds by strength reduction and Bishop's simplified
method of slices finds on a circular slip surface. It is the mechanical
base for the hydromechanical slope cases (reservoir, drawdown). The
translational slide and rainfall pages ([slope.md](slope.md),
[rainfall.md](rainfall.md)) are a different, shallow landslide on bedrock.

## Model

A slope $H=10$ m high with a 2:1 face on a foundation $H$ deep; toe
platform 0–20 m at $z=-10$ m, face 20–40 m, crest 40–60 m at $z=0$, base
$z=-20$ m. Plane strain; base fixed, sides on rollers. Mohr–Coulomb
$c=10$ kPa, $\phi=20^\circ$, $\psi=0$; $E=100$ MPa, $\nu=0.3$;
$\gamma=20$ kN/m$^3$; dry. The geometry and strength are those of the
classic homogeneous slope analysed by strength reduction in Griffiths and
Lane (1999); their paper could not be retrieved here, so their numbers are
not used.

## Method

1. Seed: no gravity, zero stress, zero steps (checkpoint).
2. Gravity: two elastic PT steps (checkpoint). The second leaves the mesh
   at rest. The step is 0.02 s per metre of mesh spacing, below the
   explicit stable step, so the explicit trials can restart from it
   directly.
3. Trials: for each $F$, restart with $c/F$ and $\tan\phi/F$ and run
   explicit DES (no PT), quasi-static mass scaling (inertial scaling
   $10^5$, characteristic speed $10^{-5}$ m/s), DES's default damping 0.8,
   for 2000 s on the fixed original mesh. A trial is coming to rest when
   its peak speed over the last five frames is below that three quarters
   through the run, and sliding otherwise.

The factor of safety lies between the largest $F$ at rest and the smallest
that slides. Bishop's FS is the minimum over circles through any entry and
exit point (`bishop_reference.py`); it reproduces $\tan\phi/\tan\beta$ for
a shallow circle in sand and is converged in the number of slices.

## Results

| Mesh | Largest $F$ at rest | Smallest $F$ sliding | FS | vs Bishop 1.369 |
| --- | --- | --- | --- | --- |
| 1 m | 1.38 (20.6 mm) | 1.39 | 1.385 | +1.2% |
| 0.5 m | 1.36 (21.6 mm) | 1.37 | 1.365 | −0.3% |
| 0.25 m | 1.35 (50.9 mm) | 1.36 | 1.355 | −1.0% |

The midpoints fall by 0.020 and 0.010 as the mesh halves, toward about
1.345, 1.7% below Bishop: strength reduction finds the mechanism, which need
not be a circle, and a finer mesh resolves the shear band better. The
mechanism passes through the toe and leaves the crest where Bishop's
critical circle does (centre (23.4, 12.6) m, radius 22.9 m). Below the
factor of safety the displacement grows smoothly (0.25 m: 4.5 mm at $F=1$,
12.6 mm at 1.34) and jumps by an order of magnitude at it; sliding trials
are still accelerating after 2000 s.

Under elastic gravity ($K_0=\nu/(1-\nu)=0.43$) the deep foundation already
yields slightly at $F=1$; it does not take part in the mechanism.

## Practical notes

- Run the trials explicitly. With zero dilation, PT at DES's default
  tolerance can stop converging well below the factor of safety (on the
  0.5 m mesh already at $F=1.1$, where the slope clearly stands), so a PT
  stall does not mark failure.
- Keep the mesh fixed (`has_moving_mesh = false`). On a moving mesh a
  failing slope slides, flattens and stops (1 m mesh: 44 mm at $F=1.39$ up
  to 1.07 m at $F=1.5$), so every trial ends at rest. That is run-out, not
  a factor of safety, which refers to the original geometry.

## Reproduce

```
python3 benchmarks/poromechanics/bishop_reference.py
python3 benchmarks/poromechanics/strength_reduction_slope.py EXE OUT --spacings 1 0.5 0.25
```

The review's embankment page ships one trial with its restart files and
the study's numbers (`record.json`).

References: Bishop (1955), Géotechnique 5(1), 7–17; Griffiths and Lane
(1999), Géotechnique 49(3), 387–403.
