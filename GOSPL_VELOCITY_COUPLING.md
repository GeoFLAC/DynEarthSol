# GoSPL Velocity Coupling

This document describes the **velocity coupling** extension to the DynEarthSol (DES) – GoSPL interface. With velocity coupling enabled, DES surface node vertical velocities are passed to GoSPL at each coupling step as its internal uplift/subsidence field (`upsub`). This replaces the need for GoSPL to read its own tectonic file and ensures that the landscape evolution model sees the same tectonic forcing as the geodynamic model.

## Background

### Basic elevation-only coupling (default)

Without velocity coupling, the exchange at each GoSPL call is:

```
DES surface topography  →  GoSPL (sets elevation field)
                           GoSPL runs erosion/deposition
GoSPL erosion change    →  DES (updates surface node positions)
```

GoSPL's internal tectonic uplift/subsidence (`upsub`) is either:
- read from a separate tectonic file specified in the GoSPL YAML, or
- absent (zero), treating topography as purely erosion-modified.

In either case, the uplift/subsidence that DES resolves through its FEM solve is **not** communicated to GoSPL. GoSPL therefore computes river incision, hillslope diffusion, and drainage reorganisation on a surface that does not rise or sink in response to the actual geodynamic deformation. This can cause drainage networks to form incorrectly or not respond dynamically to active deformation.

### Velocity coupling

With `gospl_velocity_coupling = 1`, the exchange becomes:

```
DES surface topography   →  GoSPL (sets elevation field)
DES surface vertical Vz  →  GoSPL upsub field (m/yr)
                            GoSPL runs erosion WITH correct uplift rate
Pure erosion only        →  DES (uplift contribution removed to avoid double-counting)
```

GoSPL now sees the correct tectonic context, so slope-driven processes (stream power erosion, hillslope diffusion) respond dynamically to the geodynamic velocity field.

---

## The Double-Counting Problem and Its Solution

When DES passes vertical velocities to GoSPL, both models apply the same uplift:

| Model | What it applies |
|-------|----------------|
| DES FEM solve | Moves surface nodes by `vz × dt` (tectonic uplift) |
| GoSPL (via `upsub`) | Internally raises its elevation field by `upsub × dt` |

If we naively take `gospl_after − gospl_before` as the elevation change to apply back to DES, we would apply the uplift a second time, creating runaway surface growth.

**Fix:** The pure erosion signal is isolated by subtracting the uplift contribution:

```
Δh_erosion = (gospl_after − gospl_before) − vz_yr × Δt
```

where `vz_yr` is the surface node vertical velocity in m yr⁻¹ and `Δt` is the accumulated coupling interval in years. Only `Δh_erosion` is applied back to DES, exactly as in the elevation-only coupling scheme.

---

## Algorithm

At each **coupling step** the sequence is:

1. **Push DES topography → GoSPL**
   `apply_elevation_data()` interpolates current DES surface elevations onto the GoSPL mesh.

2. **Push DES vertical velocities → GoSPL `upsub`**
   `apply_velocity_data()` sends surface node vertical velocities (converted to m yr⁻¹) to GoSPL's uplift/subsidence field. Only executed when `gospl_velocity_coupling = 1`.

3. **Query elevation before erosion**
   `interpolate_elevation_to_points()` snapshots the GoSPL elevation at DES node positions (`z_before`). This baseline captures any interpolation smoothing so it cancels out in the differencing step.

4. **Run GoSPL**
   `run_processes_for_dt()` advances GoSPL by the accumulated coupling time `Δt`.
   When velocity coupling is on, `skip_tectonics = true` is passed, which prevents GoSPL's internal `getTectonics()` from overwriting the `upsub` values just set in step 2.

5. **Query elevation after erosion**
   `interpolate_elevation_to_points()` snapshots the updated GoSPL elevation at DES positions (`z_after`).

6. **Compute pure erosion**
   ```
   Δh_total  = z_after − z_before                 (erosion + uplift)
   Δh_erosion = Δh_total − vz_yr × Δt             (subtract uplift; velocity coupling only)
   erosion_rate = Δh_erosion / Δt_seconds          (m s⁻¹)
   ```

7. **Apply erosion to DES**
   The first step's worth of erosion (`erosion_rate × dt`) is applied immediately. On subsequent non-coupling steps the rate is spread gradually with linear extrapolation, exactly as in the base coupling scheme (see [GOSPL_COUPLING.md](GOSPL_COUPLING.md)).

---

## Configuration

Add the following to your DynEarthSol input parameter file (`.cfg`):

```cfg
[control]
surface_process_option = 11
surface_process_gospl_config_file = ./gospl_config.yml
gospl_velocity_coupling = 1
```

### Full parameter reference (GoSPL-related)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `surface_process_option` | 0 | Set to `11` for GoSPL coupling |
| `surface_process_gospl_config_file` | `""` | Path to GoSPL YAML configuration file |
| `gospl_coupling_frequency` | 1 | Run GoSPL every N DES steps |
| `gospl_rate_change_tolerance` | 0.3 | Adaptive coupling threshold (0 = fixed N) |
| `gospl_mesh_resolution` | -1 | Mesh spacing in metres (-1 = auto) |
| `gospl_initial_topo_amplitude` | 0.0 | Random initial topography amplitude (m) |
| `gospl_mesh_perturbation` | 0.3 | Node jitter fraction (0 = regular grid) |
| `gospl_velocity_coupling` | 0 | **Set to `1` to enable velocity coupling** |

### Interaction with other parameters

| `gospl_coupling_frequency` | `gospl_velocity_coupling` | Behaviour |
|:-:|:-:|---|
| 1 | 0 | Elevation-only coupling every step (original scheme) |
| N > 1 | 0 | Gradual erosion application with linear extrapolation, no velocity exchange |
| 1 | 1 | Velocity coupling every step |
| N > 1 | 1 | Velocity coupling at coupling steps; velocity at the coupling step is held constant across the N-step interval for the extrapolation |

Adaptive coupling frequency (`gospl_rate_change_tolerance`) works unchanged alongside velocity coupling.

---

## GoSPL Configuration Notes

When `gospl_velocity_coupling = 1`:

- GoSPL's **`tecdata`** key in the YAML config can be left empty or omitted. The `skip_tectonics` flag prevents `getTectonics()` from being called, so any pre-existing tectonic file is silently bypassed for each run.
- GoSPL's `upsub` (vertical) and `hdisp` (horizontal) fields are both updated by `apply_velocity_data()`. Currently only the vertical component is populated from DES data; horizontal components are set to zero.
- If you previously had a tectonic file that you want GoSPL to retain (e.g. for a background plate motion), do **not** enable velocity coupling — the two approaches are mutually exclusive.

---

## Unit Conventions

| Quantity | DES internal unit | GoSPL unit | Conversion |
|----------|------------------|------------|------------|
| Time step `dt` | seconds (s) | years (yr) | `÷ 3.1536 × 10⁷` |
| Vertical velocity `vz` | m s⁻¹ | m yr⁻¹ | `× 3.1536 × 10⁷` |
| Elevation | metres (m) | metres (m) | — |
| Erosion rate (internal) | m s⁻¹ | — | applied as `rate × dt_s` |

---

## Implementation Details

### Files modified

| File | Change |
|------|--------|
| `gospl_driver/gospl-driver.hpp` | Added `bool velocity_coupling` member; added `skip_tectonics` parameter to `run_processes_for_dt()`; declared `apply_velocity_data()` |
| `gospl_driver/gospl-driver.cxx` | Constructor initialises `velocity_coupling = false`; `run_processes_for_dt()` forwards `skip_tectonics` to C API; new `apply_velocity_data()` wraps C interface |
| `parameters.hpp` | Added `bool gospl_velocity_coupling` to `Control` struct |
| `input.cxx` | Registered `control.gospl_velocity_coupling` option (default `false`) |
| `dynearthsol.cxx` | Passes `gospl_velocity_coupling` to driver at startup; prints status message |
| `bc.cxx` | **Step 1b** (new): builds velocity array in m yr⁻¹ and calls `apply_velocity_data()` before querying before-elevations; **Step 3**: passes `skip_tectonics = velocity_coupling`; **erosion loop**: subtracts `vz_yr × Δt` from raw GoSPL diff when coupling is on |

### Key functions

```cpp
// gospl-driver.hpp / gospl-driver.cxx
int  GoSPLDriver::apply_velocity_data(const double* coords,
                                      const double* velocities,  // m/yr, [vx, vy, vz] per point
                                      int num_points, double timer,
                                      int k = 3, double power = 1.0);

double GoSPLDriver::run_processes_for_dt(double dt,
                                         bool verbose = false,
                                         bool skip_tectonics = false);
```

The `skip_tectonics` parameter maps directly to the `skip_tectonics` argument of the underlying `gospl_extensions` C API function `run_processes_for_dt()`.

---

## Limitations and Future Work

- **Single-step velocity snapshot**: The surface velocity is sampled once per coupling step and held constant for the duration of the coupling interval. This is consistent with the piecewise-constant rate approximation used for erosion in the gradual-application scheme.
- **Vertical only**: Horizontal velocities (`hdisp`) are zeroed. Horizontal advection of the GoSPL mesh by DES horizontal surface velocities is not yet implemented.
- **Static GoSPL mesh**: The GoSPL mesh is generated once at initialisation. Large lateral displacements of the DES surface may eventually carry surface nodes outside the GoSPL domain; the 5% boundary buffer mitigates edge-interpolation errors.
- **Velocity used is from the current coupling step**: If the simulation undergoes remeshing between coupling steps, the velocity at the coupling step is used. Post-remeshing changes to node numbering are handled because the velocity array is rebuilt from scratch using current `top_nodes` at each coupling step.

---

## See Also

- [GOSPL_COUPLING.md](GOSPL_COUPLING.md) — Overview of the DES–GoSPL coupling including mesh generation and the gradual erosion application scheme
- [`gospl_driver/DESIGN_DECISIONS.md`](gospl_driver/DESIGN_DECISIONS.md) — Design rationale and deferred improvements
- [`gospl_driver/coupling_improvements_explained.md`](gospl_driver/coupling_improvements_explained.md) — ASCII diagrams explaining the gradual application and adaptive frequency schemes
- [GoSPL documentation](https://gospl.readthedocs.io) — GoSPL configuration reference
