# Shear-band localization: EP and EP-RSF

These three extensional models compare shear-band development in elastoplastic
(EP) and rate-and-state friction (EP-RSF) materials, using an initial friction
angle of 30 degrees and zero dilation.

| Configuration | Rheology | a | b | a/b |
| --- | --- | --- | --- | --- |
| `ep_phi30_psi0.cfg` | EP | - | - | - |
| `ep_rsf_strengthening_phi30_psi0.cfg` | EP-RSF | 0.018 | 0.015 | 1.2 |
| `ep_rsf_weakening_phi30_psi0.cfg` | EP-RSF | 0.003 | 0.015 | 0.2 |

The 20 by 10 km domain uses an internally generated uniform triangular grid with
100 m spacing. No external mesh file is needed. A circular seed of radius 500 m,
centered at (10, -9.5) km, starts with accumulated plastic strain 1; the surroundings
start at zero. Cohesion decreases from 30 to 0.1 MPa over plastic strain 0 to 0.5.
The side boundaries move outward at 0.5 cm/yr each, the top is traction-free, and the
base uses Winkler support. The elastic moduli are K = 50 GPa and G = 30 GPa, density
is 2700 kg/m3, and gravity is 10 m/s2. Initial stress is lithostatic.

All three models use global-velocity mass scaling with c = 2e5, the shear reference
speed, and FLAC damping factor 0.8. Both RSF cases use the aging law, Dc = 0.005 m,
V0 = 1e-9 m/s, and the total deviatoric strain-rate invariant
(`control.rsf_slip_rate_projection_option = 1`).

## Run

Use the compatible 2D executable specified in the [main guide](../README.md);
the solver retained on this benchmark branch predates the required RSF controls.
Run from this directory with that executable. Each case has a distinct output prefix.

```bash
mkdir -p output
OMP_NUM_THREADS=8 /path/to/dynearthsol2d ep_phi30_psi0.cfg
OMP_NUM_THREADS=8 /path/to/dynearthsol2d ep_rsf_strengthening_phi30_psi0.cfg
OMP_NUM_THREADS=8 /path/to/dynearthsol2d ep_rsf_weakening_phi30_psi0.cfg
```

Each configuration integrates 12000 years and writes fields every 500 years.
Monitor records are written every 1000 steps. The middle monitoring point starts
at the seed center; the other two points lie near the side boundaries.

## Inspecting the results

Compare the accumulated plastic strain with its initial field to identify the
shear bands developing from the weak seed. Use the same color scale for all
three models. If estimating band dips, use the same depth interval in each model.
The Arthur angle for the initial friction and dilation angles is 52.5 degrees.

The monitor CSVs provide coordinates, velocity, stress, friction, and state
outputs. Seed displacement is the magnitude of the monitored seed node's change
in position from its initial location. The second invariant of deviatoric stress
in the monitored seed element is
`sqrt(0.25 * (stress_0 - stress_1)^2 + stress_2^2)`; divide by `1e6` for MPa.
Use these histories to compare smooth deformation with episodic displacement
and stress drops.

The inputs passed three-step startup checks on their original meshes with the
source revision linked in the main guide. Those checks do not establish the
final shear-band geometry or the full 12000-year response.
