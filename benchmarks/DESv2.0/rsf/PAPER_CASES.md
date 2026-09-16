# Rate-and-state friction: the calculations reported in the paper

These are the configurations for Sect. 5.3, Sect. 5.4 and Appendices C and E of the DES v2.0
paper. For a quick introduction, start with the smaller cases in
`examples/rate_and_state_friction/` in the source repository. The parameter study here contains
74 calculations.

| Directory | Cases | Paper |
|---|---|---|
| `shear_box/` | 9 | Sect. 5.3, Fig. 10 — simple shear against analytical solutions |
| `decollement/grid/` | 60 | Sect. 5.4.2 — the `a/b`–`D_c` parameter grid |
| `decollement/contour/` | 14 | Appendix E.2 — constant `L_inf` at 325 and 450 m |
| `decollement/resolution/` | 2 | Appendix E.3 — coarse and fine meshes |
| `decollement/density_floor/` | 2 | Appendix E.4 — the paired mass-scaling test |
| `decollement/wavefield/` | 2 | Appendix C — see the note at the end |
| `decollement/other/` | 7 | classified in Appendix E.2, outside the 74 |
| `strike_slip/` | 9 | Sect. 5.4.3 — comparison with Herrendörfer et al. (2018) |
| `mesh/` | 3 | the meshes the cases name |

Case names carry their parameters: `ab0p5_dc0p005` is `a/b = 0.5` with `D_c = 5` mm,
`ct_L650_ab0p45` sits on the `L_inf = 325` m contour at `a/b = 0.45` (the name carries twice
`L_inf`), and `tr_*` and `dc25_*` are further grid points.

## Running a case

The application configurations name their meshes as `../../mesh/<file>` (`../mesh/` for
`strike_slip/`); the shear-box benchmark generates its unit-square mesh internally. The
application cases write into `output/`, which has to exist first. Run cases in separate copies
of their directory because many configurations share output paths. For example:

```bash
cd decollement/grid
mkdir -p output
OMP_NUM_THREADS=8 /path/to/dynearthsol2d ab0p5_dc0p005.cfg
```

Set `OMP_NUM_THREADS` yourself. Left unset, OpenMP takes every core it can see, which on a large
node is far slower for these problems than a modest thread count.

The output base path is set by `sim.modelname`, and the monitor path by
`monitor.output_prefix`. The décollement cases monitor one node at
`x = 50` km on the top of the fault zone, `y = -4.95` km. The three strike-slip `*prof*`
configurations monitor 150 nodes, giving the 75 station pairs used for the slip profiles; the
other six use eight pairs at selected locations.

These are long. A décollement case integrates 2000 to 3000 model years and a strike-slip case
up to 300 yr. The time step is set by the mass-scaled elastic limit and the state-evolution limit
together, so wall-clock time runs to hours or days depending on the case and the machine.

## What the configurations set

| Entry | Meaning |
|---|---|
| `state_var_model = 1` | the aging law |
| `direct_a`, `evolution_b` | `a` and `b`, one value per material |
| `characteristic_velocity` | the reference velocity `V_0` |
| `characteristic_distance` | `D_c` |
| `control.rsf_slip_rate_projection_option = 1` | the continuum rate measure of Eq. (5.6), `V = 2 w eps_II` from the total deviatoric strain rate |
| `control.rsf_dtheta_max = 0.2` | for adaptive aging-law cases, bounds the fractional state change over a step, Eq. (5.9) |
| `control.inertial_scaling` | the mass-scaling coefficient `c` |
| `control.mass_scaling_reference_speed` | the modulus of the apparent-speed ceiling; `shear` floors the fictitious density at `rho K / G`, `bulk` lets it return to `rho` |
| `control.damping_option = 1` | FLAC damping |

Across the décollement study `a` is held at 0.003 and `a/b` is set through `evolution_b`, so
`ab0p5_*` has `b = 0.006` and `ab0p2_*` has `b = 0.015`.

## Changing a case

The four knobs the paper varies:

- **Response class.** `evolution_b` sets `a/b` and `characteristic_distance` sets `D_c`. Raising
  either moves a case from rapid episodes toward slow transients and then toward motion near the
  loading rate. Appendix E.2 gives the classification and the peak speeds that separate the
  classes.
- **Resolution.** The element size inside the fault zone comes from the regional size entries in
  the `.poly`, not from the configuration; `decollement/resolution/` differs from the grid cases
  in `mesh` section entries and in the `.poly` regions. The rate measure `V` scales with element
  size, so a resolution change moves it directly — that is the point of Appendix E.3, not a
  side effect.
- **Mass scaling.** `control.inertial_scaling` and `control.mass_scaling_reference_speed`. The
  pair in `density_floor/` differs only in the second, and Appendix E.4 reports what that does to
  peak speed, episode duration, stress drop and recurrence variability.
- **How long.** `max_time_in_yr`, with `output_time_interval_in_yr` and the `[monitor]`
  `step_interval` controlling how much is written. The décollement cases classify on the record
  after the first 1000 yr, so a shorter run will not reproduce the classification.

Two things to leave alone unless they are the subject of the test. In the adaptive aging-law
cases, `control.rsf_dtheta_max` prevents a step from spanning the state transient; it defaults
to off. The fixed-time-step shear-box cases do not use this adaptive limiter. The fault-zone
geometry in the `.poly` controls its cross-band discretization, so mesh changes also change the
element-scale rate measure.

## What these need

Use DES commit `b94e994d694e980e6a2a638fda00f6100acbd150` or a release that contains it. The
adaptive aging-law cases require `control.rsf_dtheta_max`, and the `density_floor/` pair requires
`control.mass_scaling_reference_speed`.

## The two Appendix C cases

`wavefield/wf1_baseline.cfg` and `wf2_baseline.cfg` do not start from time zero. They continue
saved states to 950.2 and 951.05 yr to write dense output across one rapid episode. Appendix C
uses divergence and the out-of-plane curl of the velocity field as compressional- and
shear-sensitive diagnostics. The required checkpoints are not included in this subset, and the
absolute `restarting_from_modelname` values record their original locations. These two inputs
also predate `control.rsf_dtheta_max` and therefore omit it.
