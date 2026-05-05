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

1. `dynearthsol3d-gospl` should be built first in the root direcotry.
2. Run from the directory containing both input files: e.g., in the current directory,

```bash
../../dynearthsol3d-gospl gaussian-weakzone-3d-with-gospl.cfg
```

GoSPL is driven automatically by the coupling layer; no separate GoSPL invocation
is needed.

## Scripts

The `scripts/` subdirectory contains helper utilities (`umeshFcts.py`) used by
the GoSPL driver to build and manipulate the unstructured mesh. They are imported
at runtime and do not need to be invoked manually.
