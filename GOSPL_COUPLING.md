# GoSPL Coupling in DynEarthSol

This document describes how to enable and configure the coupling between DynEarthSol (DES) and GoSPL (Global Surface Processes Landscape).

## Prerequisites

- **GoSPL**: Must be installed and accessible.
- **Python**: DynEarthSol must be built with Python support (`HAS_GOSPL_CPP_INTERFACE` defined).
- **Libraries**: `numpy`, `scipy` must be available in the Python environment.

## Configuration

To enable GoSPL coupling, set the following in your DynEarthSol input parameter file (`.cfg`):

```cfg
[Control]
surface_process_option = 11
surface_process_gospl_config_file = gospl_config.yml
gospl_coupling_interval_in_yr = 1000.0  # Coupling time step
```

## Mesh Generation

DynEarthSol automatically generates a compatible mesh for GoSPL at the start of the simulation.

- **Process**:
    - The current surface nodes of DynEarthSol are extracted.
    - A regular generic mesh (strip) is created covering the bounding box of the surface nodes.
    - **Dimensions**: `N x N` grid where `N = sqrt(num_surface_nodes)`.
    - **Coordinates**: `z` is set to 0. `x` and `y` cover the extent.
    - **Output**: `gospl_mesh.npz` is saved in the simulation directory.

## GoSPL Configuration (`gospl_config.yml`)

Your GoSPL configuration file must point to the generated mesh:

```yaml
...
mesh:
  filename: gospl_mesh.npz
...
```

## Notes

- The coupling assumes a 2D plane-strain DES model coupled to a pseudo-3D strip in GoSPL.
- Elevation changes in GoSPL are applied back to DES surface nodes via interpolation.
