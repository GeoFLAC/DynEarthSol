# GoSPL Coupling Example: Gaussian Weak Zone 3D Rift

This directory contains example input files for running a 3D Gaussian weak zone
rift model with DynEarthSol coupled to GoSPL for surface processes.

## Input Files

### `gaussian-weakzone-3d-with-gospl.cfg`

DynEarthSol configuration file for the Gaussian weak zone 3D rift model.

Key settings:
- **Domain**: 100 km (x) × 80 km (y) × 10 km (z), 1 km base resolution
- **Run time**: 1 Myr with 20 kyr output interval
- **Boundary conditions**: ±1.5 cm/yr extension in x-direction
- **GoSPL coupling**: enabled via `surface_process_option = 11`
  - Config file: `gospl_config_gaussian_weakzone_3D.yml`
  - Coupling frequency: every 200 steps
  - GoSPL mesh resolution: 500 m
  - Velocity coupling: all 3 DES velocity components passed to GoSPL

### `gospl_config_gaussian_weakzone_3D.yml`

GoSPL configuration file paired with the DynEarthSol `.cfg` above. Parameters
are chosen to match the ASPECT + FastScape reference model
(`gaussian_weakzone_3D.prm`).

Key settings:
- **Flow direction**: multi-flow (`flowdir: 6`, equivalent to FastScape `p=-1`)
- **Boundary conditions**: east and west boundaries open, north and south closed (`bc: '1010'`)
- **Stream power law**: K = 1×10⁻⁵ m^(1-2m)/yr, m = 0.4, n = 1
- **Hillslope diffusion**: Ka = 1×10⁻² m²/yr
- **Sea level**: −2000 m (prevents marine flooding of the initial surface)
- **Output**: `output_gaussian_weakzone_3D/` every 20 kyr

## Usage

1. Build `dynearthsol-gospl` first in the root directory (see `gospl_driver/README.md`).
2. Run from this directory so the relative path to the GoSPL config resolves correctly:

```bash
conda activate gospl
../../dynearthsol-gospl gaussian-weakzone-3d-with-gospl.cfg
```

GoSPL is driven automatically by the coupling layer; no separate GoSPL invocation
is needed.

## Known Issues

### Intermittent PETSc error during coupling

When `seadepo: true` is set in the GoSPL config, you may occasionally see the following during
a coupling step:

```
Error in run_and_get_erosion: error code 77
[0] TSSolve() ...
[0] TSAdaptChoose() ...
[0] Unexpected state: bad hmax in TSAdaptChoose()
Error: GoSPL run_and_get_erosion failed
```

This is a floating-point issue in PETSc's Rosenbrock W-scheme adaptive time stepper inside
GoSPL's nonlinear marine deposition solver. The stepper occasionally lands just past the
integration end time (`self.dt`), making the remaining interval `hmax = max_time -
current_time` slightly negative. The coupled run recovers automatically: DynEarthSol skips
applying erosion for that one step and proceeds normally.

If the sea level is well below the model surface (e.g., `sea: position: -2000.`), marine
deposition is physically inactive and `seadepo: false` can be set to avoid this code path
entirely.
