# RSF configurations

Choose a configuration by the problem you want to run. There are 108 cfg files
and three external meshes. The [main guide](README.md) gives the compatible
solver and representative commands for the four main problem groups.

| Directory | Configurations | Purpose |
|---|---:|---|
| [shear_box/](shear_box/) | 9 | Verify EP, steady-state RSF, and aging-law RSF against simple-shear reference solutions |
| [localization/](localization/README.md) | 3 | Compare shear-band development in EP and strengthening/weakening EP-RSF |
| [decollement/grid/](decollement/grid/) | 60 | Vary the friction-parameter ratio `a/b` and characteristic distance `D_c` |
| [decollement/contour/](decollement/contour/) | 14 | Vary `a/b` at constant nucleation half-length `L_inf = 325` or `450 m` |
| [decollement/resolution/](decollement/resolution/) | 2 | Compare coarse and fine fault-zone discretizations |
| [decollement/density_floor/](decollement/density_floor/) | 2 | Compare shear- and bulk-based mass-scaling reference speeds |
| [decollement/wavefield/](decollement/wavefield/) | 2 | Sample velocity disturbances across an episode; requires unavailable restart checkpoints |
| [decollement/other/](decollement/other/) | 7 | Additional baseline, friction-parameter, and resolution configurations |
| [strike_slip/](strike_slip/) | 9 | Compare band discretization, characteristic distance, and slip-monitoring setups |

## Meshes and output paths

Three external mesh files are supplied in [mesh/](mesh/). The decollement cfgs
reference them through `../../mesh/`, and the strike-slip cfgs through `../mesh/`.
Keep that directory structure when copying a case. Shear-box and localization
meshes are generated internally.

Output paths are set by `sim.modelname` and `monitor.output_prefix`. Many cases
share output names, so use a separate working copy for each run as described in
the main guide. The decollement cases run for 2000-3000 model years, strike-slip
cases for up to 300 years, and localization cases for 12000 years.

## Parameter conventions

- `ab0p5_dc0p005` denotes `a/b = 0.5` and `D_c = 0.005 m`. The decollement grid
  fixes `a = 0.003` and varies `b` through `evolution_b`.
- `ct_L650_ab0p45` denotes the `L_inf = 325 m` contour at `a/b = 0.45`; the filename
  carries twice `L_inf`. Names beginning with `tr_` and `dc25_` are further grid points.
- All RSF configurations explicitly select `control.rsf_slip_rate_projection_option = 1`,
  the total-deviatoric-strain invariant rate. The compatible source retains
  option 0 as its default, so keep the explicit option in the cfgs.
- Adaptive aging-law cases use `control.rsf_dtheta_max = 0.2`. Fixed-step shear-box
  cases do not use this adaptive limiter; the two historical wavefield restarts
  also omit it.
- The density-floor pair varies `control.mass_scaling_reference_speed` between
  `shear` and `bulk`. No source patch needs to be applied to the compatible solver.

## Wavefield restart requirements

`decollement/wavefield/wf1_baseline.cfg` and `wf2_baseline.cfg` continue saved
states to 950.2 and 951.05 years to sample a rapid episode. Their required
checkpoints are absent, and their absolute `restarting_from_modelname` values
record the original locations. To run these inputs, obtain the matching
checkpoints and update those paths. They cannot run from this collection alone.
These configurations predate the adaptive state-step limiter and omit it.
