# RSF paper-case inventory

The 108 configurations support Sections 5.3-5.4 and Appendices D-G of the revised
DES v2.0 paper. The parameter study consists of the 60 grid cases and 14 contour
cases; the remaining inputs are verification, localization, comparison, or
sensitivity cases. See the [main guide](README.md) for the compatible solver and
representative run commands.

| Directory | Configurations | Paper |
|---|---:|---|
| [shear_box/](shear_box/) | 9 | Section 5.3.1, Figure 10, Appendix E: simple-shear verification |
| [localization/](localization/README.md) | 3 | Section 5.3.2, Figure 11, Appendix F: EP and strengthening/weakening EP-RSF |
| [decollement/grid/](decollement/grid/) | 60 | Section 5.4.2: `a/b`-`D_c` grid |
| [decollement/contour/](decollement/contour/) | 14 | Appendix G.2: constant `L_inf` at 325 and 450 m |
| [decollement/resolution/](decollement/resolution/) | 2 | Appendix G.3: mesh-resolution sensitivity |
| [decollement/density_floor/](decollement/density_floor/) | 2 | Appendix G.4: mass-scaling density-floor comparison |
| [decollement/wavefield/](decollement/wavefield/) | 2 | Appendix D: restart inputs; checkpoints unavailable |
| [decollement/other/](decollement/other/) | 7 | Appendix G.2: additional cases outside the 74-case study |
| [strike_slip/](strike_slip/) | 9 | Section 5.4.3: Herrendoerfer et al. (2018) comparison |

Three external mesh files are supplied in [mesh/](mesh/). The decollement cfgs
reference them through `../../mesh/`, and the strike-slip cfgs through `../mesh/`.
Keep that directory structure when copying a case. Shear-box and localization
meshes are generated internally.

## Input conventions

- `ab0p5_dc0p005` denotes `a/b = 0.5` and `D_c = 0.005 m`. The decollement study
  fixes `a = 0.003` and varies `b` through `evolution_b`.
- `ct_L650_ab0p45` denotes the `L_inf = 325 m` contour at `a/b = 0.45`; the filename
  carries twice `L_inf`. Names beginning with `tr_` and `dc25_` are further grid points.
- All RSF configurations explicitly select `control.rsf_slip_rate_projection_option = 1`,
  the total-deviatoric-strain invariant rate. The compatible PR source retains
  option 0 as its default, so keep the explicit option in the cfgs.
- Adaptive aging-law cases use `control.rsf_dtheta_max = 0.2`. Fixed-step shear-box
  cases do not use this adaptive limiter; the two historical wavefield restarts
  also omit it.
- The density-floor pair varies `control.mass_scaling_reference_speed` between
  `shear` and `bulk`. No source patch needs to be applied to the compatible solver.

Output paths are set by `sim.modelname` and `monitor.output_prefix`. Run cases in
separate copies as described in the main guide. The decollement cases run for
2000-3000 model years, strike-slip cases for up to 300 years, and localization
cases for 12000 years. Short runs do not establish the reported response classes.

## Appendix D: missing restart checkpoints

`decollement/wavefield/wf1_baseline.cfg` and `wf2_baseline.cfg` continue saved
states to 950.2 and 951.05 years to sample a rapid episode. They cannot run from
this collection alone: the required checkpoints are absent, and their absolute
`restarting_from_modelname` values record the original locations. These inputs
predate the adaptive state-step limiter. Appendix D uses velocity divergence
and out-of-plane curl as compressional- and shear-sensitive diagnostics.

## Missing strike-slip sensitivity input

The final-condition three-element comparison with per-element `D_c = 0.01 m`
is not included. The nine archived strike-slip configurations should not be
substituted for that missing sensitivity run.

## Validation scope

The compatible PR source passed the nine-case simple-shear and zero-rate healing
checks. The three localization inputs passed three-step startup checks on their
original meshes. The full long-running campaigns have not been repeated on that
source revision; these checks do not establish complete paper reproduction.
