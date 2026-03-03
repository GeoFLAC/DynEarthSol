# DES–GoSPL Coupling: Design Review and Proposed Redesign

This document reviews the current coupling scheme between DynEarthSol (DES) and GoSPL
from three perspectives — **straightforwardness**, **accuracy**, and **efficiency** —
and proposes a cleaner design assuming the freedom to modify both `gospl_extensions`
and the DES `gospl_driver`.

---

## 1. What the Current Scheme Does

At each coupling step the sequence in `bc.cxx / use_gospl()` is:

```
(A) apply_elevation_data(DES_coords, DES_z)     DES nodes → GoSPL mesh  (IDW, builds KDTree from DES nodes)
(B) interpolate_elevation_to_points(DES_coords)  GoSPL mesh → DES nodes  (IDW, builds KDTree from GoSPL mesh) → z_before
(C) run_processes_for_dt(total_dt)
(D) interpolate_elevation_to_points(DES_coords)  GoSPL mesh → DES nodes  (IDW, builds KDTree from GoSPL mesh) → z_after
(E) erosion[i] = z_after[i] - z_before[i]
```

With velocity coupling enabled, step (A) pre-adds `vz × Δt_seconds` to the elevations
before giving them to GoSPL, so GoSPL sees the "future" topography.

Between coupling steps, the stored `erosion_rate` is applied at each DES step with
linear extrapolation.

---

## 2. Problems Identified

### 2.1 Straightforwardness

**The before/after trick is non-obvious.**
Steps (B) and (D) exist only to cancel the IDW smoothing error introduced by step (A).
The cancellation reasoning ("same smoothing in both directions") is subtle and breaks
down if anything changes the GoSPL mesh between the two calls. It also obscures the
real intent: *get the change GoSPL made to the surface*.

**Velocity coupling is a workaround for a missing API.**
The `apply_velocity_data` function comes from `DataDrivenTectonics`, which is designed
for time-series velocity archives identified by a `timer` key. This is completely the
wrong abstraction for "set a constant uplift rate for the next run". The fragility of
this mismatch caused the flat-surface bug (GoSPL interpreted `timer = var.time` as a
cumulative-displacement look-up, not a per-interval rate). The pre-uplift workaround
fixes the symptom but is physically approximate (see §2.2).

**`runProcessesForDt` manipulates `tEnd` to get one step.**
Internally it sets `self.tEnd = self.tNow + dt` to make GoSPL's main loop run for
exactly one iteration. Any GoSPL logic that checks `tNow ≥ tEnd` for file output,
stratigraphy saves, or other side effects is silently triggered.

---

### 2.2 Accuracy

**GoSPL's drainage state is reset every coupling step.**
Step (A) overwrites all of GoSPL's `hGlobal` with IDW-smoothed DES elevation. This
forces GoSPL to recompute flow accumulation, chi values, and drainage area from scratch
at the next run. Real river networks have geological memory — once a major trunk stream
captures a watershed, it does not reset every N DES steps. Discarding drainage history
every coupling step systematically underestimates the efficiency of the drainage
network and misrepresents how rivers reorganise in response to deformation.

**IDW smoothing of the input elevation distorts GoSPL's drainage.**
Step (A) is a low-pass filter (IDW with k = 3). Local slope features that drive stream
power erosion are attenuated *before* GoSPL runs. The before/after cancellation in
step (E) removes the smoothing from the *output* signal, but it cannot undo the effect
that smoother input slopes had on GoSPL's *internal* flow routing and erosion during
step (C). The more often the elevation is reset, the more this attenuation accumulates.

**Pre-uplift applies tectonic displacement as a step function.**
When velocity coupling is on, all of `vz × Δt` is added to the elevation *before*
GoSPL starts. GoSPL then computes erosion on that elevated surface for the full
coupling interval. In reality, the surface rises *continuously* during the interval;
the slope and drainage area at each GoSPL sub-step should reflect the incrementally
rising surface. For large intervals (N = 200) or high uplift rates, the step-function
approximation systematically overestimates slopes early and underestimates them late.

---

### 2.3 Efficiency

**Two KDTree builds per coupling step are redundant.**
`interpolate_elevation_to_points` builds a `cKDTree(self.mCoords)` on every call.
Steps (B) and (D) call it with identical `DES_coords` and an identical GoSPL mesh. The
tree is built and discarded twice. Since GoSPL's mesh (`mCoords`) never changes, this
tree should be built once and cached.

**Three IDW passes per coupling step where one would suffice.**
Step (A): IDW from DES nodes to GoSPL mesh (M GoSPL nodes × k DES neighbours).
Step (B): IDW from GoSPL mesh to DES nodes (N DES nodes × k GoSPL neighbours).
Step (D): same as (B).

The purpose of (B) and (D) is to produce a single quantity: the change GoSPL made to
its elevation field, mapped to DES nodes. This requires only ONE interpolation of the
elevation difference, computed directly in GoSPL's native mesh space. If the elevation
reset (step A) is done less frequently, the cost drops further (see §3.4).

---

## 3. Proposed Redesign

The redesign has two independent but complementary parts:

- **§3.1–3.3**: replace the before/after trick with `run_and_get_erosion`, and replace
  the pre-uplift hack with a proper `set_uplift_rate` API.
- **§3.4**: preserve GoSPL's drainage state between coupling steps by calling
  `apply_elevation_data` only when necessary (initialization, remeshing) and using a
  gentle drift correction otherwise.

Together they eliminate the three IDW passes, the before/after trick, and the drainage
state resets.

---

### 3.1 New `gospl_extensions` Python API (`EnhancedModel`)

```python
def set_uplift_rate(self, src_pts, vz_yr, k=3, power=1.0):
    """
    Interpolate DES vertical velocities onto GoSPL mesh nodes and store as upsub.
    Applied (and cleared) during the very next run_and_get_erosion() call.

    :param src_pts: (N, 3) DES surface node coordinates
    :param vz_yr:   (N,) vertical velocity in m/yr at those nodes
    """
    tree = spatial.cKDTree(src_pts)
    distances, idx = tree.query(self.mCoords, k=k)
    weights = 1.0 / np.maximum(distances, 1e-20) ** power
    wsum = weights.sum(axis=1)
    self._upsub_override = (weights * vz_yr[idx]).sum(axis=1) / wsum
    on = np.where(distances[:, 0] < 1e-20)[0]
    if on.size:
        self._upsub_override[on] = vz_yr[idx[on, 0]]


def run_and_get_erosion(self, dt, query_pts, k=3, power=1.0):
    """
    Run GoSPL for dt years and return net erosion/deposition at query_pts.

    The difference hGlobal_after - hGlobal_before is computed directly on the
    GoSPL mesh (no IDW error in the difference itself), then interpolated once
    to query_pts using the cached KDTree.

    If set_uplift_rate() was called beforehand, that upsub field replaces
    getTectonics(); it is consumed and cleared after the run.

    :param dt:        coupling interval in years
    :param query_pts: (N, 3) DES surface node coordinates
    :return:          (N,) net elevation change in metres (negative = erosion)
    """
    # 1. Snapshot hGlobal in GoSPL's native space — no IDW
    h_before = self.hGlobal.getArray().copy()

    # 2. Apply upsub override if set_uplift_rate() was called
    skip = hasattr(self, '_upsub_override')
    if skip:
        self.upsub.setArray(self._upsub_override)
        del self._upsub_override

    # 3. Run GoSPL; upsub (tectonic) + erosion both update hGlobal
    self.runProcessesForDt(dt, skip_tectonics=skip)

    # 4. Elevation change in GoSPL's native space (exact — IDW error is zero here)
    delta_h = self.hGlobal.getArray() - h_before

    # 5. Single IDW interpolation of delta_h from GoSPL mesh to DES nodes
    tree = self._get_mesh_tree()               # cached; never rebuilt
    distances, idx = tree.query(query_pts, k=k)
    weights = 1.0 / np.maximum(distances, 1e-20) ** power
    wsum = weights.sum(axis=1)
    erosion = (weights * delta_h[idx]).sum(axis=1) / wsum
    on = np.where(distances[:, 0] < 1e-20)[0]
    if on.size:
        erosion[on] = delta_h[idx[on, 0]]
    return erosion


def apply_drift_correction(self, src_pts, des_elevation, alpha=0.1, k=3, power=1.0):
    """
    Gently nudge GoSPL's hGlobal toward DES's current elevation to prevent
    slow drift accumulation without fully resetting the drainage state.

    h_gospl_new = h_gospl + alpha * (h_des_on_mesh - h_gospl)

    alpha = 1.0 is a full reset; alpha = 0.1 is a 10% nudge.
    Recommended: alpha in [0.1, 0.3] applied every 10–20 coupling steps.

    :param src_pts:      (N, 3) DES surface node coordinates
    :param des_elevation (N,) DES surface node elevations
    :param alpha:        correction strength in [0, 1]
    """
    tree = spatial.cKDTree(src_pts)
    distances, idx = tree.query(self.mCoords, k=k)
    weights = 1.0 / np.maximum(distances, 1e-20) ** power
    wsum = weights.sum(axis=1)
    h_des_on_mesh = (weights * des_elevation[idx]).sum(axis=1) / wsum
    on = np.where(distances[:, 0] < 1e-20)[0]
    if on.size:
        h_des_on_mesh[on] = des_elevation[idx[on, 0]]

    h = self.hGlobal.getArray()
    h += alpha * (h_des_on_mesh - h)   # in-place blend

    if hasattr(self, 'hLocal') and hasattr(self, 'locIDs'):
        self.hLocal.getArray()[:] = h[self.locIDs]


def _get_mesh_tree(self):
    """Return a cached KDTree of GoSPL mesh node coordinates (built once)."""
    if not hasattr(self, '_mesh_kdtree'):
        from scipy.spatial import cKDTree
        self._mesh_kdtree = cKDTree(self.mCoords, leafsize=10)
    return self._mesh_kdtree
```

---

### 3.2 New C API additions in `gospl_extensions.h`

```c
/**
 * Interpolate DES vertical velocities onto GoSPL's upsub field.
 * Applied during the next run_and_get_erosion() call; cleared after.
 */
int set_uplift_rate(ModelHandle handle,
                    const double* coords, const double* vz_yr,
                    int num_points, int k, double power);

/**
 * Run GoSPL for dt years and return net erosion at query coordinates.
 * The elevation difference is computed on GoSPL's mesh (no IDW error),
 * then interpolated once to the query points.
 * If set_uplift_rate() was called, that upsub is used and then consumed.
 */
int run_and_get_erosion(ModelHandle handle,
                        double dt,
                        const double* coords, int num_points,
                        double* erosion,
                        int k, double power);

/**
 * Gently blend GoSPL's elevation toward DES's elevation to limit drift.
 * alpha = 1.0 is a full reset; alpha ≈ 0.1–0.2 is a gentle nudge.
 */
int apply_drift_correction(ModelHandle handle,
                           const double* coords, const double* des_elevation,
                           int num_points, double alpha, int k, double power);
```

---

### 3.3 New `GoSPLDriver` state and wrappers

New members in `gospl-driver.hpp`:

```cpp
// --- Drainage state continuity ---
bool   needs_elevation_reset;       // true on first step and after remeshing
int    elevation_sync_counter;      // coupling steps since last full reset
int    elevation_sync_interval;     // do drift correction every N steps (0 = disabled)
double elevation_drift_alpha;       // drift correction strength [0, 1]; default 0.2
```

New wrapper methods:

```cpp
int  set_uplift_rate(const double* coords, const double* vz_yr,
                     int num_points, int k = 3, double power = 1.0);

int  run_and_get_erosion(double dt,
                         const double* coords, int num_points,
                         double* erosion_m,
                         int k = 3, double power = 1.0);

int  apply_drift_correction(const double* coords, const double* des_elevation,
                             int num_points, double alpha,
                             int k = 3, double power = 1.0);
```

New input parameters (in `parameters.hpp` / `input.cxx`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gospl_elevation_sync_interval` | 10 | Apply gentle drift correction every N coupling steps (0 = only on reset events) |
| `gospl_drift_correction_alpha` | 0.2 | Fraction of the DES–GoSPL elevation gap corrected per sync event |

---

### 3.4 Drainage State Continuity: When to Reset Elevation

The fundamental insight enabling drainage continuity is that with correct velocity
coupling, GoSPL's elevation naturally tracks DES's elevation without needing a full
reset every step:

| Source of elevation change | DES | GoSPL (proposed) |
|----------------------------|-----|-----------------|
| Tectonic uplift/subsidence | FEM solve: `Δz = vz × dt` | `upsub = vz_yr`; GoSPL applies same `Δz = vz_yr × dt` |
| Erosion/deposition | GoSPL feedback applied to nodes | Computed by GoSPL directly |

Both models advance by the same tectonic increment each coupling step.
The remaining sources of drift are small and slow:

- **IDW interpolation error at initialization** — one-time, decays quickly.
- **Horizontal material advection** — as DES surface nodes move laterally, the material
  at a fixed GoSPL mesh location changes. The drift rate is approximately
  `v_horiz × ∇z`, typically < 0.01 m per coupling step for geological extension
  rates and moderate slopes.

This means a full elevation reset is needed only in two situations:

1. **Initialization** — GoSPL has no elevation information yet.
2. **After DES remeshing** — surface node positions change discontinuously;
   GoSPL's spatial state no longer corresponds to the DES surface.

Between these events, a **gentle drift correction** (`apply_drift_correction` with
`alpha ≈ 0.2`) applied every `elevation_sync_interval` coupling steps is sufficient
to keep the models aligned without disrupting the drainage network.

**`use_gospl()` coupling step with drainage continuity:**

```cpp
// --- Detect reset conditions ---
// Post-remeshing: pending_erosion_rate size changes when ntop changes.
if (var.gospl_driver->pending_erosion_rate.size() != ntop)
    var.gospl_driver->needs_elevation_reset = true;

// --- Elevation management ---
if (var.gospl_driver->needs_elevation_reset) {
    // Full reset: DES is authoritative (first step or post-remesh)
    gospl_driver->apply_elevation_data(coords.data(), des_elevations.data(),
                                        ntop, 3, 1.0);
    gospl_driver->needs_elevation_reset = false;
    gospl_driver->elevation_sync_counter = 0;

} else if (gospl_driver->elevation_sync_interval > 0 &&
           gospl_driver->elevation_sync_counter % gospl_driver->elevation_sync_interval == 0) {
    // Gentle drift correction — does NOT reset drainage state
    gospl_driver->apply_drift_correction(coords.data(), des_elevations.data(),
                                          ntop, gospl_driver->elevation_drift_alpha,
                                          3, 1.0);
}
// (otherwise: GoSPL evolves freely — drainage state fully preserved)
gospl_driver->elevation_sync_counter++;

// --- Velocity coupling ---
if (gospl_driver->velocity_coupling) {
    std::vector<double> vz_yr(ntop);
    for (size_t i = 0; i < ntop; i++)
        vz_yr[i] = (*var.vel)[top_nodes[i]][NDIMS-1] * SEC_PER_YR;
    gospl_driver->set_uplift_rate(coords.data(), vz_yr.data(), ntop, 3, 1.0);
}

// --- Run GoSPL and get erosion in one call ---
std::vector<double> net_erosion(ntop);
gospl_driver->run_and_get_erosion(total_dt, coords.data(), ntop,
                                   net_erosion.data(), 3, 1.0);

// --- Apply to DES (gradual application and adaptive frequency unchanged) ---
for (size_t i = 0; i < ntop; i++) {
    if (/* inside boundary buffer */) {
        erosion_rate[i] = net_erosion[i] / total_dt_seconds;
        (*var.coord)[top_nodes[i]][NDIMS-1] += erosion_rate[i] * dt_seconds;
    }
    gospl_driver->pending_erosion_rate[i] = erosion_rate[i];
}
```

---

## 4. Comparison

| Aspect | Current scheme | Proposed scheme |
|--------|---------------|-----------------|
| **Coupling step logic** | Before-snapshot → run → after-snapshot → subtract | Push elev (when needed) → set uplift → run-and-get-erosion |
| **IDW calls per normal coupling step** | 3 | 1 (only GoSPL→DES of Δh) |
| **IDW calls when drift correction fires** | 3 | 2 (DES→GoSPL blend + GoSPL→DES Δh) |
| **KDTree builds per coupling step** | 3 (no caching) | 0–1 (cached GoSPL tree, DES tree only at reset) |
| **Differencing location** | After IDW (IDW error added to signal) | On GoSPL's native mesh (no IDW error in Δh) |
| **Velocity coupling** | Pre-uplift step function (approximate) | Continuous `upsub` throughout GoSPL run (exact) |
| **GoSPL drainage state** | Fully reset every coupling step | Preserved between coupling steps; only nudged gently |
| **Drainage network memory** | None — discarded every N DES steps | Maintained across coupling steps |
| **IDW smoothing of erosion physics** | Repeated every coupling step | Only at initialization and remeshing events |

### What remains unchanged

- Gradual application of erosion rate with linear extrapolation between coupling steps.
- Adaptive coupling frequency driven by the L2 rate-change metric.
- Boundary buffer excluding DES nodes near the GoSPL mesh edge.

---

## 5. Remaining Longer-Term Improvements

**5.1 Deforming GoSPL mesh**
The GoSPL mesh is static. As DES deforms, surface nodes drift laterally and eventually
reach the 5% boundary buffer. For simulations > ~1 Myr with high horizontal velocities,
the effective domain shrinks noticeably. Periodically regenerating the GoSPL mesh from
the current DES surface node positions — re-interpolating GoSPL's state including
`hGlobal`, drainage areas, and sediment thicknesses onto the new mesh — would
eliminate the boundary problem entirely.

**5.2 Replace `tEnd` manipulation in `runProcessesForDt`**
`EnhancedModel.run_one_step(dt)` already calls individual process methods directly
without the `runProcesses()` main-loop overhead. Using `run_one_step` inside
`runProcessesForDt` in place of the `tEnd` trick would avoid silently triggering
GoSPL's output/stratigraphy/flexure side effects.
