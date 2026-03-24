# GoSPL Coupling in DynEarthSol

This document describes how to enable and configure the coupling between DynEarthSol (DES) and GoSPL (Global Surface Processes Landscape) for landscape evolution modeling.

## Prerequisites

- **GoSPL**: Must be installed and accessible via the Python environment.
- **Python**: DynEarthSol must be built with Python support (`HAS_GOSPL_CPP_INTERFACE` defined).
- **Libraries**: `numpy`, `scipy` must be available in the Python environment.

## Configuration

To enable GoSPL coupling, set the following in your DynEarthSol input parameter file (`.cfg`):

```cfg
[Control]
surface_process_option = 11
surface_process_gospl_config_file = gospl_config.yml
```

### Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `surface_process_option` | 0 | Set to `11` for GoSPL coupling |
| `surface_process_gospl_config_file` | "" | Path to GoSPL YAML configuration file |
| `gospl_coupling_frequency` | 1 | Run GoSPL every N DES steps (higher values reduce overhead) |
| `gospl_mesh_resolution` | -1 | Mesh node spacing in meters (-1 = auto-sized based on DES surface nodes) |
| `gospl_mesh_perturbation` | 0.3 | Fraction of grid spacing to randomly perturb node positions (0-1) |

### Example Configuration

```cfg
[control]
surface_process_option = 11
surface_process_gospl_config_file = ./gospl_config.yml
gospl_coupling_frequency = 200
gospl_mesh_resolution = 500
gospl_mesh_perturbation = 0.3
```

## Mesh Generation

DynEarthSol automatically generates a compatible mesh for GoSPL at the start of the simulation.

### Process

1. Surface nodes of DynEarthSol are extracted.
2. A regular grid mesh is created covering the bounding box of the surface nodes.
3. Grid dimensions are determined by:
   - **Auto-sizing** (when `gospl_mesh_resolution = -1`): `N x N` grid where `N = sqrt(num_surface_nodes)`
   - **User-specified resolution**: Grid nodes calculated from domain extent / resolution
4. If `gospl_mesh_perturbation > 0`, interior node positions are randomly offset to reduce grid artifacts.
5. A Delaunay triangulation is performed and saved as `gospl_mesh.npz`.

### Output

- **File**: `gospl_mesh.npz` (saved in the simulation directory)
- **Contents**: Vertices (x, y, z coordinates) and triangular cells

## GoSPL Configuration (`gospl_config.yml`)

Your GoSPL configuration file must point to the generated mesh:

```yaml
mesh:
  filename: gospl_mesh.npz
```

Refer to GoSPL documentation for other configuration options (erosion rates, precipitation, etc.).

## Coupling Mechanism

### Two-Way Coupling Algorithm

The coupling uses a "pure erosion" approach to avoid interpolation smoothing errors:

1. Before running GoSPL: Interpolate GoSPL elevations to DES surface points → `z_before`
2. Run GoSPL erosion/deposition processes for the accumulated time
3. After running GoSPL: Interpolate updated GoSPL elevations to DES points → `z_after`
4. Compute erosion: `dz = z_after - z_before`
5. Apply `dz` to DES surface nodes

This differencing cancels out interpolation smoothing, ensuring only true erosion/deposition changes are applied.

### Coupling Frequency

When `gospl_coupling_frequency = N`:
- GoSPL runs every N DES time steps
- GoSPL receives the accumulated time (N × dt)
- The full elevation change is applied in one step

Use higher coupling frequencies for erosion-dominated systems to reduce computational overhead.

### Boundary Buffer

A 5% buffer on each edge of the GoSPL mesh prevents extrapolation errors for nodes that move outside the original mesh domain during simulation.

## Notes

- The coupling assumes a 2D plane-strain DES model coupled to a pseudo-3D strip in GoSPL.
- Elevation changes in GoSPL are applied back to DES surface nodes via IDW interpolation.
- The mesh is generated once at initialization (static mesh approach).
- For large deformation simulations, consider the limitations of the static mesh approach (see `gospl_driver/DESIGN_DECISIONS.md`).

## See Also

- [`gospl_driver/DESIGN_DECISIONS.md`](gospl_driver/DESIGN_DECISIONS.md) - Design decisions and deferred improvements
- [`gospl_driver/examples/core-complex-3d-with-gospl.cfg`](gospl_driver/examples/core-complex-3d-with-gospl.cfg) - Example configuration
