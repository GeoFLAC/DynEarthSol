# GoSPL Velocity Coupling

This document describes the **velocity coupling** extension to the DynEarthSol (DES) – GoSPL interface. With velocity coupling enabled, the expected tectonic displacement of DES surface nodes is pre-added to the elevation sent to GoSPL at each coupling step. GoSPL then runs its landscape evolution on the correctly-uplifted surface, so erosion rates respond dynamically to the geodynamic velocity field without any risk of double-counting.

## Background

### Basic elevation-only coupling (default)

Without velocity coupling, the exchange at each GoSPL call is:

```
DES surface topography (current)  →  GoSPL (sets elevation field)
                                      GoSPL runs erosion/deposition
GoSPL erosion change              →  DES (updates surface node positions)
```

GoSPL computes river incision, hillslope diffusion, and drainage reorganisation on the **current** DES topography without knowing that the crust is actively rising or sinking. This means:

- In a fast-uplifting region, GoSPL's slopes are underestimated → erosion is too slow.
- Drainage networks do not respond dynamically to deformation.

### Velocity coupling

With `gospl_velocity_coupling = 1`, the exchange becomes:

```
DES topography + expected uplift  →  GoSPL (pre-uplifted elevation)
                                      GoSPL runs erosion on uplifted surface
GoSPL erosion change              →  DES (pure erosion, no double-counting)
```

GoSPL sees the topography as it will appear after the current coupling interval's tectonic motion, so stream power incision and hillslope diffusion are driven by the correct slopes.

---

## Why Pre-Uplift Rather Than Passing `upsub`?

An alternative approach would be to use GoSPL's `apply_velocity_data()` API to set its internal `upsub` (uplift/subsidence rate) field and then subtract the resulting displacement from the before/after difference. This was the initial implementation but it was abandoned because:

1. **`apply_velocity_data` uses a `timer` argument** that GoSPL interprets as the time point of the velocity data. Passing the current simulation time (e.g. 500 000 yr) caused GoSPL to compute a cumulative displacement from t = 0, producing a displacement orders of magnitude larger than the one-interval increment — and making the subtraction correction wildly incorrect.

2. **GoSPL applies `upsub` progressively** across its own internal sub-steps, so `gospl_after − gospl_before` does not simply equal `upsub × Δt + erosion`; the interaction between uplift and erosion within sub-steps makes the correction inexact.

The pre-uplift approach sidesteps both issues: GoSPL's internal state is never touched for velocity, and no correction is needed.

---

## Algorithm

At each **coupling step** the sequence is:

1. **Build pre-uplifted elevation**
   For each surface node `i`:
   ```
   z_preuplifted[i] = z_DES[i] + vz[i] × Δt_seconds
   ```
   where `vz[i]` is the DES vertical velocity (m s⁻¹) and `Δt_seconds` is the accumulated coupling interval in seconds.
   When `gospl_velocity_coupling = 0` this is just `z_DES[i]`.

2. **Push pre-uplifted elevation → GoSPL**
   `apply_elevation_data()` interpolates the pre-uplifted elevations onto the GoSPL mesh.

3. **Query elevation before erosion**
   `interpolate_elevation_to_points()` snapshots GoSPL's elevation at DES node positions (`z_before`). Any interpolation smoothing is recorded in this baseline so it cancels in the differencing step.

4. **Run GoSPL**
   `run_processes_for_dt()` advances GoSPL by the accumulated coupling time `Δt` (years). GoSPL performs erosion/deposition only — no `upsub` is set, so no correction is needed.

5. **Query elevation after erosion**
   `interpolate_elevation_to_points()` gives `z_after`.

6. **Compute pure erosion**
   ```
   Δh_erosion = z_after − z_before        (pure erosion; GoSPL added no uplift)
   erosion_rate = Δh_erosion / Δt_seconds  (m s⁻¹)
   ```

7. **Apply erosion to DES**
   The first step's worth of erosion (`erosion_rate × dt`) is applied immediately. On subsequent non-coupling steps the rate is distributed gradually with linear extrapolation, exactly as in the base coupling scheme (see [GOSPL_COUPLING.md](GOSPL_COUPLING.md)).

   The DES FEM solve simultaneously applies the tectonic deformation (`vz × dt` per step) — this is the DES side of the uplift. The two contributions are additive and independent.

### Why there is no double-counting

| What applies the uplift | Source |
|-------------------------|--------|
| DES FEM solve | Moves surface nodes by `vz × dt` each step |
| GoSPL (pre-uplift approach) | **None** — GoSPL only adds erosion to the pre-uplifted surface |

Because GoSPL never moves its elevation field due to an uplift rate, `gospl_after − gospl_before` is strictly erosion/deposition.

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
| N > 1 | 0 | Gradual erosion application with linear extrapolation |
| 1 | 1 | Velocity coupling every step; `z_preuplifted = z_DES + vz × dt` |
| N > 1 | 1 | Velocity coupling at coupling steps; uplift uses `vz` at that step times the full accumulated `Δt` |

Adaptive coupling frequency (`gospl_rate_change_tolerance`) works unchanged alongside velocity coupling.

---

## GoSPL Configuration Notes

When `gospl_velocity_coupling = 1`, GoSPL's `tecdata` key in the YAML config may be left empty or omitted. The pre-uplift approach does not interact with GoSPL's internal tectonic machinery at all (`getTectonics()` is called normally), so a pre-existing tectonic file will still run if present. If you have a tectonic file that provides background plate motion or other forcing, it will add on top of the DES-derived pre-uplift.

---

## Unit Conventions

| Quantity | DES internal unit | Used as | Conversion |
|----------|------------------|---------|------------|
| Time step `dt` | seconds (s) | `total_dt_seconds` | `accumulated_dt_yr × 3.1536 × 10⁷` |
| Vertical velocity `vz` | m s⁻¹ | m s⁻¹ (multiplied by `total_dt_seconds`) | — |
| Pre-uplift displacement | metres (m) | added to z before `apply_elevation_data` | — |
| Elevation | metres (m) | metres (m) | — |
| Erosion rate (internal) | m s⁻¹ | applied as `rate × dt_s` per step | — |

---

## Implementation Details

### Files modified

| File | Change |
|------|--------|
| `gospl_driver/gospl-driver.hpp` | Added `bool velocity_coupling` member; added `skip_tectonics` parameter to `run_processes_for_dt()` (kept for future use); declared `apply_velocity_data()` (kept for future use) |
| `gospl_driver/gospl-driver.cxx` | Constructor initialises `velocity_coupling = false`; `run_processes_for_dt()` forwards `skip_tectonics` to C API; `apply_velocity_data()` wraps C interface |
| `parameters.hpp` | Added `bool gospl_velocity_coupling` to `Control` struct |
| `input.cxx` | Registered `control.gospl_velocity_coupling` option (default `false`) |
| `dynearthsol.cxx` | Passes `gospl_velocity_coupling` to driver at startup; prints status message |
| `bc.cxx` | **Step 1** (modified): adds `vz × Δt_seconds` to `des_elevations[i]` when velocity coupling is on before calling `apply_elevation_data()`; no other changes to the coupling loop |

### Key code change in `bc.cxx`

```cpp
// Step 1: build elevation to hand to GoSPL
for (std::size_t i = 0; i < ntop; ++i) {
    int n = top_nodes[i];
    double z = (*var.coord)[n][NDIMS-1];
    if (var.gospl_driver->velocity_coupling) {
        z += (*var.vel)[n][NDIMS-1] * total_dt_seconds;  // pre-add tectonic uplift
    }
    des_elevations[i] = z;
}
apply_elevation_data(coords, des_elevations, ntop, 3, 1.0);

// Steps 2–6 are unchanged; no correction in the erosion loop.
erosion = gospl_elev_after[i] - gospl_elev_before[i];  // pure erosion
```

---

## Limitations and Future Work

- **Single-step velocity snapshot**: `vz` is sampled once per coupling step and assumed constant for the full interval `Δt`. This is the same piecewise-constant approximation used for the erosion rate in the gradual-application scheme.
- **Vertical only**: The pre-uplift applies only to the vertical coordinate. Horizontal advection of the GoSPL mesh by DES horizontal surface velocities is not implemented.
- **Static GoSPL mesh**: The GoSPL mesh is generated once at initialisation. Large lateral displacements of DES surface nodes may eventually carry them outside the GoSPL domain; the 5% boundary buffer mitigates edge-interpolation errors.
- **`apply_velocity_data` / `skip_tectonics` retained**: These driver functions are preserved in the API for potential future use (e.g. a hybrid scheme), but are not called in the current coupling loop.

---

## See Also

- [GOSPL_COUPLING.md](GOSPL_COUPLING.md) — Overview of the DES–GoSPL coupling including mesh generation and the gradual erosion application scheme
- [`gospl_driver/DESIGN_DECISIONS.md`](gospl_driver/DESIGN_DECISIONS.md) — Design rationale and deferred improvements
- [`gospl_driver/coupling_improvements_explained.md`](gospl_driver/coupling_improvements_explained.md) — ASCII diagrams explaining the gradual application and adaptive frequency schemes
- [GoSPL documentation](https://gospl.readthedocs.io) — GoSPL configuration reference
