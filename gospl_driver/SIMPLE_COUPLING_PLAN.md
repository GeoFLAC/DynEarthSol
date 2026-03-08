# DES–GoSPL Simple Coupling (ASPECT–FastScape Scheme)

> **Status: Implemented.** All steps below have been completed. See the
> [Implementation Notes](#implementation-notes) section at the bottom for key
> divergences from the original plan.

## Motivation

The current `gospl_coupling` scheme resets GoSPL's internal topography (`hGlobal`) from DES
at every coupling step via `apply_elevation_data`. This prevents GoSPL from accumulating
hydrological state (drainage network, flow directions) across steps.

Neuharth et al. (2022, *Tectonics*, doi:10.1029/2021TC007166) couple ASPECT and FastScape
every 5 ASPECT timesteps using a simpler, more principled scheme:

- **FastScape is initialized once** from ASPECT's topography.
- At each coupling step, ASPECT sends all three surface velocity components (vx, vy, vz)
  to FastScape.
- FastScape applies vz as uplift and vx/vy for lateral advection of the topography, then
  runs its own erosion/deposition.
- The topographic change Δh = h_after − h_before is applied to ASPECT as a mesh velocity.
- FastScape's topography is **never reset** from ASPECT (except after remeshing).

This eliminates drift corrections, the double-ownership problem, and the need for extrapolation
between coupling steps. GoSPL retains its drainage network across coupling steps.

---

## Current State (gospl_coupling branch)

### `use_gospl()` in `bc.cxx` (lines 1590–1719):
```
Every N DES steps:
  Step 1: apply_elevation_data  ← resets GoSPL hGlobal from DES every time
  Step 2: interpolate h_before
  Step 3: run_processes_for_dt
  Step 4: interpolate h_after
  Step 5: Δh = h_after − h_before → apply to DES nodes
```

### `GoSPLDriver` (`gospl-driver.hpp`) — fields that exist:
- `coupling_frequency`, `step_counter`, `accumulated_dt`
- `mesh_bounds[4]`, `mesh_bounds_valid`

### `GoSPLDriver` — methods that exist:
- `apply_elevation_data`, `interpolate_elevation_to_points`
- `run_processes_for_dt` (no `skip_tectonics` param yet)

### `gospl_extensions` stack — what already exists:

| Layer | File | Function |
|---|---|---|
| Python model | `enhanced_model.py` | `set_uplift_rate(src_pts, vz_yr)` |
| Python model | `enhanced_model.py` | `run_and_get_erosion(dt, query_pts)` |
| Python interface | `gospl_python_interface.py` | `set_uplift_rate(handle, coords, vz_yr, ...)` |
| C bridge | `gospl_extensions.cpp` | `set_uplift_rate(...)` |
| C header | `gospl_extensions.h` | `set_uplift_rate(...)` |

### What is missing for full 3-component velocity coupling:
- No `set_surface_velocity(vx, vy, vz)` anywhere in the stack
- No horizontal advection of `hGlobal` in `run_and_get_erosion`
- `run_processes_for_dt` in `GoSPLDriver` does not forward `skip_tectonics`
- No `needs_elevation_reset` flag on `GoSPLDriver`

---

## Target Scheme

```
Initialization (once):
  apply_elevation_data  →  GoSPL owns topography from here on

Every N DES steps (coupling step):
  if needs_elevation_reset (post-remesh only):
      apply_elevation_data  →  re-sync after node repositioning
      needs_elevation_reset = false
  set_surface_velocity(vx_yr, vy_yr, vz_yr from DES)
  h_before = interpolate_elevation_to_points
  run_processes_for_dt(total_dt, skip_tectonics=true)
      [inside GoSPL: apply horiz advection Δx=vx*dt, Δy=vy*dt to hGlobal
                     apply vz as upsub, then run erosion/deposition]
  h_after  = interpolate_elevation_to_points
  Δh = h_after − h_before  →  apply to DES nodes

Non-coupling steps:
  DES moves normally; GoSPL completely idle.
  No erosion applied to DES between coupling steps.
```

---

## Implementation Steps

### Step 1: Add `set_surface_velocity` to `EnhancedModel` (Python)

**File:** `/home/echoi2/opt/gospl_extensions/gospl_model_ext/enhanced_model.py`

Add a new method after `set_uplift_rate`:

```python
def set_surface_velocity(self, src_pts, vx_yr, vy_yr, vz_yr, k=3, power=1.0):
    """
    Interpolate all three DES surface velocity components onto the GoSPL mesh
    and store for the next run_and_get_erosion() call.

    - vz_yr is stored as _upsub_override (vertical uplift/subsidence, m/yr)
    - vx_yr, vy_yr are stored as _vx_override, _vy_override (horizontal, m/yr)
      and applied as semi-Lagrangian advection of hGlobal before running GoSPL.

    :param src_pts: (N, 3) DES surface node coordinates
    :param vx_yr:   (N,)  x-velocity at each DES node in m/yr
    :param vy_yr:   (N,)  y-velocity at each DES node in m/yr
    :param vz_yr:   (N,)  z-velocity (uplift) at each DES node in m/yr
    :param k:       number of IDW neighbours (default 3)
    :param power:   IDW power exponent (default 1.0)
    """
    from scipy.spatial import cKDTree
    src_pts = np.asarray(src_pts, dtype=np.float64)
    vx_yr   = np.asarray(vx_yr,   dtype=np.float64)
    vy_yr   = np.asarray(vy_yr,   dtype=np.float64)
    vz_yr   = np.asarray(vz_yr,   dtype=np.float64)
    k = max(1, min(int(k), src_pts.shape[0]))
    src_tree = cKDTree(src_pts, leafsize=10)
    dists, idxs = src_tree.query(self.mCoords, k=k)
    if k == 1:
        dists = dists[:, None]
        idxs  = idxs[:, None]
    eps = 1.0e-20
    weights = 1.0 / np.maximum(dists, eps) ** power
    onIDs = np.where(dists[:, 0] < eps)[0]
    if onIDs.size > 0:
        weights[onIDs] = 0.0
        weights[onIDs, 0] = 1.0
    weights /= weights.sum(axis=1, keepdims=True)
    self._upsub_override = (weights * vz_yr[idxs]).sum(axis=1)
    self._vx_override    = (weights * vx_yr[idxs]).sum(axis=1)
    self._vy_override    = (weights * vy_yr[idxs]).sum(axis=1)
```

### Step 2: Apply horizontal advection in `run_and_get_erosion` (Python)

**File:** `/home/echoi2/opt/gospl_extensions/gospl_model_ext/enhanced_model.py`

Modify `run_and_get_erosion` to apply semi-Lagrangian horizontal advection when
`_vx_override`/`_vy_override` are set, before taking the `h_before` snapshot:

```python
def run_and_get_erosion(self, dt, query_pts, k=3, power=1.0):
    from scipy.spatial import cKDTree
    skip = hasattr(self, '_upsub_override') and self._upsub_override is not None
    has_horiz = (hasattr(self, '_vx_override') and self._vx_override is not None)

    # --- Horizontal advection (semi-Lagrangian) ---
    # h_new(x) = h_old(x - v_horiz * dt)
    # Trace each GoSPL node back along its horizontal velocity to find
    # where the material came from, then interpolate elevation from there.
    if has_horiz:
        displaced = self.mCoords.copy()
        displaced[:, 0] -= self._vx_override * dt   # vx in m/yr, dt in yr → metres
        displaced[:, 1] -= self._vy_override * dt
        mesh_tree = self._get_mesh_tree()
        k_adv = max(1, min(int(k), self.mCoords.shape[0]))
        adv_dists, adv_idxs = mesh_tree.query(displaced, k=k_adv)
        if k_adv == 1:
            adv_dists = adv_dists[:, None]
            adv_idxs  = adv_idxs[:, None]
        eps = 1.0e-20
        adv_w = 1.0 / np.maximum(adv_dists, eps) ** power
        on = np.where(adv_dists[:, 0] < eps)[0]
        if on.size > 0:
            adv_w[on] = 0.0; adv_w[on, 0] = 1.0
        adv_w /= adv_w.sum(axis=1, keepdims=True)
        h = self.hGlobal.getArray()
        h_advected = (adv_w * h[adv_idxs]).sum(axis=1)
        h[:] = h_advected
        if hasattr(self, 'hLocal') and hasattr(self, 'locIDs'):
            self.hLocal.getArray()[:] = h[self.locIDs]
        self._vx_override = None
        self._vy_override = None

    # --- Vertical uplift (upsub) ---
    if skip:
        old_upsub = self.upsub.copy() if hasattr(self, 'upsub') and self.upsub is not None else None
        self.upsub = self._upsub_override[self.locIDs]

    # Snapshot, run, diff on native mesh
    h_before = self.hGlobal.getArray().copy()
    self.runProcessesForDt(dt, verbose=False, skip_tectonics=skip)
    h_after  = self.hGlobal.getArray()
    delta_h  = h_after - h_before

    if skip:
        self._upsub_override = None
        if old_upsub is not None:
            self.upsub = old_upsub

    # IDW: native-mesh delta_h → query_pts
    query_pts = np.asarray(query_pts, dtype=np.float64)
    tree = self._get_mesh_tree()
    k_q = max(1, min(int(k), self.mCoords.shape[0]))
    dists, idxs = tree.query(query_pts, k=k_q)
    if k_q == 1:
        dists = dists[:, None]; idxs = idxs[:, None]
    eps = 1.0e-20
    weights = 1.0 / np.maximum(dists, eps) ** power
    on = np.where(dists[:, 0] < eps)[0]
    if on.size > 0:
        weights[on] = 0.0; weights[on, 0] = 1.0
    weights /= weights.sum(axis=1, keepdims=True)
    return (weights * delta_h[idxs]).sum(axis=1)
```

### Step 3: Add Python interface wrapper

**File:** `/home/echoi2/opt/gospl_extensions/cpp_interface/gospl_python_interface.py`

Add after `set_uplift_rate`:

```python
def set_surface_velocity(handle: int, coords, vx_yr, vy_yr, vz_yr,
                         num_points: int, k: int = 3, power: float = 1.0) -> int:
    """Pass all three DES velocity components to GoSPL."""
    model = _models.get(handle)
    if model is None:
        return -1
    try:
        coords_array = np.asarray(coords).reshape(num_points, 3)
        vx_array = np.asarray(vx_yr).reshape(num_points)
        vy_array = np.asarray(vy_yr).reshape(num_points)
        vz_array = np.asarray(vz_yr).reshape(num_points)
        model.set_surface_velocity(coords_array, vx_array, vy_array, vz_array,
                                   k=k, power=power)
        return 0
    except Exception as e:
        print(f"Error in set_surface_velocity: {e}")
        return -1
```

### Step 4: Add C bridge function

**File:** `/home/echoi2/opt/gospl_extensions/include/gospl_extensions.h`

Add declaration:

```c
/**
 * Interpolate all three DES surface velocity components (m/yr) onto the GoSPL
 * mesh and store internally. Consumed by the next run_and_get_erosion() call.
 * vz is applied as uplift/subsidence; vx and vy drive semi-Lagrangian horizontal
 * advection of the elevation field before erosion runs.
 *
 * @param handle     Model handle
 * @param coords     DES surface node coordinates (num_points * 3)
 * @param vx_yr      X-velocity at each node in m/yr (num_points)
 * @param vy_yr      Y-velocity at each node in m/yr (num_points)
 * @param vz_yr      Z-velocity (uplift) at each node in m/yr (num_points)
 * @param num_points Number of DES surface nodes
 * @param k          IDW nearest-neighbour count
 * @param power      IDW power exponent
 * @return 0 on success, -1 on error
 */
int set_surface_velocity(ModelHandle handle,
                         const double* coords,
                         const double* vx_yr,
                         const double* vy_yr,
                         const double* vz_yr,
                         int num_points, int k, double power);
```

**File:** `/home/echoi2/opt/gospl_extensions/cpp_interface/gospl_extensions.cpp`

Add static function pointer, load in `initialize_gospl_extensions`, release in
`finalize_gospl_extensions`, and implement the C bridge function following the
same pattern as `set_uplift_rate` (lines 388–415), passing four numpy arrays
(coords, vx, vy, vz).

### Step 5: Add `set_surface_velocity` to `GoSPLDriver`

**File:** `gospl_driver/gospl-driver.hpp`

Add field and method:
```cpp
bool needs_elevation_reset;  // true at init and after remeshing
bool velocity_coupling;      // if true, send all 3 DES velocity components

int set_surface_velocity(const double* coords,
                         const double* vx_yr,
                         const double* vy_yr,
                         const double* vz_yr,
                         int num_points, int k = 3, double power = 1.0);
```

Update `run_processes_for_dt` signature:
```cpp
double run_processes_for_dt(double dt, bool verbose = false,
                             bool skip_tectonics = false);
```

**File:** `gospl_driver/gospl-driver.cxx`

Implement `set_surface_velocity` delegating to `::set_surface_velocity(...)`.
Update `run_processes_for_dt` to pass `skip_tectonics` as int to the C API.
Initialize new fields in the constructor:
```cpp
needs_elevation_reset(true), velocity_coupling(false)
```

### Step 6: Add `gospl_velocity_coupling` parameter

**File:** `parameters.hpp`
```cpp
bool gospl_velocity_coupling;  // Send all 3 DES velocity components to GoSPL (default: false)
```

**File:** `input.cxx`
```cpp
("control.gospl_velocity_coupling",
 po::value<bool>(&p.control.gospl_velocity_coupling)->default_value(false),
 "Send DES surface velocities (vx, vy, vz) to GoSPL each coupling step. "
 "vz drives uplift; vx/vy drive lateral advection of GoSPL topography.")
```

**File:** `dynearthsol.cxx`
```cpp
var.gospl_driver->velocity_coupling = param.control.gospl_velocity_coupling;
```

### Step 7: Rewrite `use_gospl()` in `bc.cxx`

Replace the current body with the new scheme. Key changes from current code:

- **Remove** `apply_elevation_data` from every coupling step (only on reset events).
- **Add** `set_surface_velocity(vx, vy, vz)` before running GoSPL (when `velocity_coupling`).
- **Replace** `run_processes_for_dt` + before/after interpolation with
  `run_and_get_erosion` (single call that snapshots, runs, diffs on native mesh).
- **Keep** boundary node skipping logic.

```cpp
void use_gospl(const Param& param, const Variables& var)
{
    if (!var.gospl_driver || !var.gospl_driver->is_initialized()) return;

    const double SEC_PER_YR = 3.1536e7;
    double dt_years = var.dt / SEC_PER_YR;
    var.gospl_driver->accumulated_dt += dt_years;
    var.gospl_driver->step_counter++;

    if (var.gospl_driver->step_counter < var.gospl_driver->coupling_frequency)
        return;  // GoSPL idle this step

    double total_dt = var.gospl_driver->accumulated_dt;
    var.gospl_driver->accumulated_dt = 0.0;
    var.gospl_driver->step_counter   = 0;

    const int_vec& top_nodes = *var.bnodes[iboundz1];
    const std::size_t ntop = top_nodes.size();
    std::vector<double> coords(ntop * 3);
    for (std::size_t i = 0; i < ntop; ++i) {
        int n = top_nodes[i];
        for (int d = 0; d < NDIMS; d++)
            coords[i * 3 + d] = (*var.coord)[n][d];
    }

    // Re-initialize GoSPL topography only on first step or after remeshing
    if (var.gospl_driver->needs_elevation_reset) {
        std::vector<double> des_elev(ntop);
        for (std::size_t i = 0; i < ntop; ++i)
            des_elev[i] = (*var.coord)[top_nodes[i]][NDIMS-1];
        var.gospl_driver->apply_elevation_data(coords.data(), des_elev.data(), ntop, 3, 1.0);
        var.gospl_driver->needs_elevation_reset = false;
    }

    // Pass DES surface velocities to GoSPL
    if (var.gospl_driver->velocity_coupling) {
        std::vector<double> vx_yr(ntop), vy_yr(ntop), vz_yr(ntop);
        for (std::size_t i = 0; i < ntop; ++i) {
            int n = top_nodes[i];
            vx_yr[i] = (*var.vel)[n][0]        * SEC_PER_YR;
            vy_yr[i] = (*var.vel)[n][1]        * SEC_PER_YR;  // 0 in 2D
            vz_yr[i] = (*var.vel)[n][NDIMS-1]  * SEC_PER_YR;
        }
        var.gospl_driver->set_surface_velocity(
            coords.data(), vx_yr.data(), vy_yr.data(), vz_yr.data(), ntop, 3, 1.0);
    }

    // Run GoSPL and get Δh in one call (h_before snapshot + run + h_after on native mesh)
    std::vector<double> net_erosion(ntop, 0.0);
    int run_result = gd->run_and_get_erosion(
        total_dt, coords.data(), ntop, net_erosion.data(), 3, 1.0);
    if (run_result != 0) {
        std::cerr << "Error: GoSPL run_and_get_erosion failed" << std::endl;
        return;
    }

    // Apply Δh to DES surface nodes (skip boundary fringe)
    bool check_bounds = var.gospl_driver->mesh_bounds_valid;
    double x_min = var.gospl_driver->mesh_bounds[0], x_max = var.gospl_driver->mesh_bounds[1];
    double y_min = var.gospl_driver->mesh_bounds[2], y_max = var.gospl_driver->mesh_bounds[3];
    double x_mg = 0.05*(x_max-x_min), y_mg = 0.05*(y_max-y_min);

    double min_dh = 0, max_dh = 0, sum_dh = 0; int skipped = 0;
    for (std::size_t i = 0; i < ntop; ++i) {
        int n = top_nodes[i];
        double x = (*var.coord)[n][0], y = (*var.coord)[n][1];
        if (check_bounds && (x < x_min+x_mg || x > x_max-x_mg ||
                             y < y_min+y_mg || y > y_max-y_mg)) { skipped++; continue; }
        double dh = net_erosion[i];
        min_dh = std::min(min_dh, dh); max_dh = std::max(max_dh, dh); sum_dh += dh;
        (*var.coord)[n][NDIMS-1] += dh;
    }
    std::cout << "GoSPL: dh [" << min_dh << ", " << max_dh << "]"
              << " mean=" << sum_dh/ntop;
    if (skipped) std::cout << " | " << skipped << " boundary nodes skipped";
    std::cout << std::endl;
}
```

Note: `run_and_get_erosion` already exists in `gospl_extensions` (both C and Python layers);
only `set_surface_velocity` (replacing `set_uplift_rate`) is new. The C wrapper for
`run_and_get_erosion` also needs to be added to `gospl-driver.hpp`/`.cxx`.

### Step 8: Set `needs_elevation_reset = true` after remeshing

Find where DES remeshing occurs (likely in `remesh.cxx` or `dynearthsol.cxx`) and add:
```cpp
if (var.gospl_driver) var.gospl_driver->needs_elevation_reset = true;
```

---

## Files to Modify

| Repo | File | Change |
|---|---|---|
| `gospl_extensions` | `gospl_model_ext/enhanced_model.py` | Add `set_surface_velocity`; modify `run_and_get_erosion` for horiz advection |
| `gospl_extensions` | `cpp_interface/gospl_python_interface.py` | Add `set_surface_velocity` wrapper |
| `gospl_extensions` | `cpp_interface/gospl_extensions.cpp` | Add C bridge for `set_surface_velocity` |
| `gospl_extensions` | `include/gospl_extensions.h` | Add C declaration |
| `DynEarthSol-velcoupling` | `gospl_driver/gospl-driver.hpp` | Add fields + method declarations |
| `DynEarthSol-velcoupling` | `gospl_driver/gospl-driver.cxx` | Implement new methods |
| `DynEarthSol-velcoupling` | `bc.cxx` | Rewrite `use_gospl()` |
| `DynEarthSol-velcoupling` | `parameters.hpp` | Add `gospl_velocity_coupling` |
| `DynEarthSol-velcoupling` | `input.cxx` | Register parameter |
| `DynEarthSol-velcoupling` | `dynearthsol.cxx` | Set flag; set `needs_elevation_reset` after remesh |

---

## What Is Removed vs. Current `gospl_coupling` Branch

| Removed | Reason |
|---|---|
| `apply_elevation_data` at every coupling step | GoSPL owns topography after init |
| `run_processes_for_dt` + before/after interpolation | Replaced by `run_and_get_erosion` |

## What Is Not Ported from `feature/gospl-velocity-coupling`

| Not ported | Reason |
|---|---|
| `pending_erosion_rate` / `prev_erosion_rate` | No interpolation between coupling steps |
| Linear extrapolation on non-coupling steps | GoSPL idle; DES moves freely |
| `elevation_sync_interval` / `elevation_drift_alpha` | No drift problem to correct |
| `adaptive_coupling_frequency` / `rate_change_tolerance` | Simplification; can add later |

---

## Open Questions

1. **Remesh site**: Where exactly in `dynearthsol.cxx` or `remesh.cxx` should
   `needs_elevation_reset = true` be set? Needs investigation before Step 8.
2. **2D models**: `vy_yr` is always 0 in 2D DES models (no Y-extent). The semi-Lagrangian
   advection will still work correctly — zero vy means no Y displacement.
3. **Non-coupling step behavior**: With no extrapolation, DES topography is frozen by
   GoSPL's contribution between coupling steps. Acceptable for typical N=200?

---

## Implementation Notes

The following deviations from the original plan occurred during implementation.

### 1. Boundary artifact fix: mesh padding instead of buffer zone

The plan showed a 5 % boundary margin where DES nodes near the GoSPL mesh edge would
be skipped (the `x_mg`/`y_mg` logic in Step 7). This was tried but produced large
artificial relief at the buffer boundary from the very first coupling step, because
GoSPL truncates drainage networks at its mesh boundary.

**Actual fix:** `generate_mesh()` now extends the GoSPL mesh by a configurable fraction
(`gospl_mesh_padding`, default 0.1 = 10 %) beyond the DES domain on every side.
All DES surface nodes therefore sit in the interior of the GoSPL mesh where drainage
is never truncated. The boundary-skip logic in `use_gospl()` was removed entirely;
`dh` is applied uniformly to all DES surface nodes.

New DES parameter (registered in `input.cxx`, stored in `parameters.hpp`):
```
gospl_mesh_padding = 0.1   # extend GoSPL mesh 10% beyond DES domain on each side
```
Pass-through in `dynearthsol.cxx` → `generate_mesh(padding)`.

### 2. Horizontal advection: GoSPL native `_varAdvector()` instead of custom IDW

The plan described a custom semi-Lagrangian IDW remap of `hGlobal` before running
GoSPL. In practice this was replaced with GoSPL's built-in face-velocity advection
for `advscheme > 0` (the default):

```python
# enhanced_model.py — run_and_get_erosion(), when has_horiz and advscheme > 0
from gospl._fortran import getfacevelocity
nodeVel = np.zeros((self.lpoints, 3), dtype=np.float64)
nodeVel[:, 0] = self._vx_override[self.locIDs]
nodeVel[:, 1] = self._vy_override[self.locIDs]
self.hdisp = nodeVel
getfacevelocity(self.lpoints, nodeVel)
old_dt = self.dt
self.dt = dt
self._varAdvector()
self.dt = old_dt
self.hdisp = None
```

The custom semi-Lagrangian path is kept as a fallback when `advscheme == 0`.
This reduces IDW calls per coupling step from 3 to 2.

### 3. Vertical uplift bug fix: `skip_tectonics=False` + `tecdata=None`

The plan used `skip_tectonics=True` in `runProcessesForDt`. This silently disabled
`applyTectonics()`, so `self.upsub` was stored but never applied.

**Fix:** Always call `runProcessesForDt(skip_tectonics=False)`. To prevent GoSPL from
re-loading tectonic data from its config file, temporarily set `self.tecdata = None`
before the call. `getTectonics()` returns early when `tecdata is None`, while
`applyTectonics()` still runs and correctly applies `self.upsub`.

```python
old_tecdata = self.tecdata
self.upsub = self._upsub_override[self.locIDs]
self.tecdata = None          # prevents config-file reload
self.runProcessesForDt(dt, verbose=False, skip_tectonics=False)
self.tecdata = old_tecdata
```

### 4. GoSPL SPL parameter tuning

Default parameters in `gospl_config_slower.yml` produced essentially no incision on
100–220 kyr timescales. Three compounding issues were identified and corrected in a
new file `gospl_config_K3e-5_Ka0.02.yml`:

| Parameter | Old | New | Reason |
|---|---|---|---|
| `spl.K` | 3e-8 | 3e-5 | ~1000× too small for regional 100-kyr scale |
| `diffusion.hillslopeKa` | 0.1 | 0.02 | Hillslope diffusion was dominating over fluvial incision |
| `sea.position` | -10. | -2000. | Rift basin nodes were going "marine", suppressing incision |
| `domain.seadepo` | True | False | Belt-and-suspenders: disable marine deposition |

DES config updated: `surface_process_gospl_config_file = ./gospl_config_K3e-5_Ka0.02.yml`

### 5. Monitoring output additions

`use_gospl()` in `bc.cxx` now prints:
- Whether velocity coupling is on/off
- Min/max of vx, vy, vz at each coupling step (when velocity_coupling=true)
- Min/max/mean of `dh` applied to DES nodes

### 6. Resolved open questions

- **Q1 (remesh site):** `needs_elevation_reset = true` is set in `remesh.cxx` at the
  point where DES node coordinates are repositioned.
- **Q2 (2D):** Confirmed — zero `vy` works correctly.
- **Q3 (N=200):** Accepted; no interpolation between coupling steps.
