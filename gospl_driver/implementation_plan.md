# Improved GoSPL Coupling: Linear Rate Extrapolation + Adaptive Frequency

## Problem

The current coupling scheme applies the full elevation change `dh` as a lump sum every N DES steps. This creates impulsive perturbations that grow with N, generating spurious transient responses. We propose distributing the erosion continuously using linear rate extrapolation and adapting N based on how rapidly erosion rates change.

## Proposed Changes

### GoSPL Driver State

#### [MODIFY] [gospl-driver.hpp](file:///home/echoi2/opt/DynEarthSol/gospl_driver/gospl-driver.hpp)

Add new state members to `GoSPLDriver` for rate history and adaptive frequency:

```cpp
// --- Gradual application state ---
std::vector<double> pending_erosion_rate;    // dh/dt rate to apply per DES step (per node, in m/s)
std::vector<double> prev_erosion_rate;       // Previous cycle's rate (for extrapolation)
bool has_prev_rate;                          // Whether prev_erosion_rate is valid (false on first cycle)
int remaining_steps;                         // Steps remaining to apply pending rate

// --- Adaptive coupling frequency ---
int base_coupling_frequency;                 // User-specified baseline N
int adaptive_coupling_frequency;             // Currently active N (may differ from base)
double rate_change_tolerance;                // Threshold for adapting N (default: 0.3 = 30% change)
double rate_change_metric;                   // Last computed rate-change metric (for diagnostics)
```

Update constructor initialization list in `gospl-driver.cxx` to initialize these.

---

### Configuration Parameters

#### [MODIFY] [parameters.hpp](file:///home/echoi2/opt/DynEarthSol/parameters.hpp)

Add to `struct Control`:

```cpp
double gospl_rate_change_tolerance;          // Rate-change threshold for adaptive N (default: 0.3)
```

#### [MODIFY] [input.cxx](file:///home/echoi2/opt/DynEarthSol/input.cxx)

Add parameter declaration:

```cpp
("control.gospl_rate_change_tolerance", po::value<double>(&p.control.gospl_rate_change_tolerance)->default_value(0.3),
 "Rate-change tolerance for adaptive GoSPL coupling frequency. "
 "When relative change exceeds this, N is halved; when below tolerance/4, N is doubled. "
 "Set to 0 to disable adaptive frequency.")
```

#### [MODIFY] [dynearthsol.cxx](file:///home/echoi2/opt/DynEarthSol/dynearthsol.cxx)

Pass new parameter to driver after initialization (near line 606):

```cpp
var.gospl_driver->rate_change_tolerance = param.control.gospl_rate_change_tolerance;
var.gospl_driver->base_coupling_frequency = param.control.gospl_coupling_frequency;
var.gospl_driver->adaptive_coupling_frequency = param.control.gospl_coupling_frequency;
```

---

### Core Algorithm Changes

#### [MODIFY] [bc.cxx](file:///home/echoi2/opt/DynEarthSol/bc.cxx)

Rewrite `use_gospl()` (lines 1589–1719) with two distinct code paths within the function:

**Non-coupling steps** (step_counter < adaptive_coupling_frequency): Instead of returning immediately, apply the linearly-extrapolated erosion rate:

```cpp
// Apply gradual erosion from previous coupling cycle
if (!var.gospl_driver->pending_erosion_rate.empty()) {
    double dt_seconds = var.dt;
    int n = top_nodes[i];
    
    // Linear extrapolation: rate + (rate - prev_rate) * progress
    double progress = (double)var.gospl_driver->step_counter / 
                      var.gospl_driver->adaptive_coupling_frequency;
    double rate_i = var.gospl_driver->pending_erosion_rate[i];
    if (var.gospl_driver->has_prev_rate) {
        double prev_rate_i = var.gospl_driver->prev_erosion_rate[i];
        rate_i += (rate_i - prev_rate_i) * progress;
    }
    
    // Limiter: clamp extrapolated rate to 2x magnitude of base rate
    double base_rate = var.gospl_driver->pending_erosion_rate[i];
    double limit = 2.0 * std::abs(base_rate) + 1e-20;  // small epsilon for zero rates
    rate_i = std::max(-limit, std::min(limit, rate_i));
    
    coord[n][NDIMS-1] += rate_i * dt_seconds;
}
return;
```

**Coupling steps** (step_counter == adaptive_coupling_frequency): Execute the full goSPL call (same as current code), then:

1. Compute the erosion rate: `rate[i] = erosion[i] / total_dt_seconds` (in m/s)
2. Store `prev_erosion_rate = pending_erosion_rate` (save previous)
3. Set `pending_erosion_rate = rate` (new rate for next cycle)
4. Set `has_prev_rate = true`
5. Apply the first step's worth of erosion: `rate[i] * var.dt`
6. Compute adaptive coupling frequency

**Adaptive frequency logic** (after computing new rates):

```cpp
// Compute rate-change metric (L2 relative change)
if (var.gospl_driver->has_prev_rate && var.gospl_driver->rate_change_tolerance > 0) {
    double sum_diff_sq = 0, sum_rate_sq = 0;
    for (size_t i = 0; i < ntop; ++i) {
        double diff = var.gospl_driver->pending_erosion_rate[i] 
                    - var.gospl_driver->prev_erosion_rate[i];
        sum_diff_sq += diff * diff;
        sum_rate_sq += var.gospl_driver->pending_erosion_rate[i] 
                     * var.gospl_driver->pending_erosion_rate[i];
    }
    double metric = (sum_rate_sq > 1e-30) ? std::sqrt(sum_diff_sq / sum_rate_sq) : 0.0;
    var.gospl_driver->rate_change_metric = metric;
    
    int N = var.gospl_driver->adaptive_coupling_frequency;
    int N_base = var.gospl_driver->base_coupling_frequency;
    double tol = var.gospl_driver->rate_change_tolerance;
    
    if (metric > tol && N > 1) {
        N = std::max(1, N / 2);       // Halve: couple more often
    } else if (metric < tol / 4 && N < 4 * N_base) {
        N = std::min(4 * N_base, N * 2);  // Double: couple less often
    }
    var.gospl_driver->adaptive_coupling_frequency = N;
    
    std::cout << "GoSPL: rate_change=" << metric 
              << " | adaptive N=" << N << std::endl;
}
```

---

### Documentation Update

#### [MODIFY] [GOSPL_COUPLING.md](file:///home/echoi2/opt/DynEarthSol/GOSPL_COUPLING.md)

Update "Coupling Frequency" section (lines 89–96) to describe the new behavior:

- Gradual application with linear extrapolation
- Adaptive frequency with configurable tolerance
- New parameter `gospl_rate_change_tolerance`

---

## User Review Required

> [!IMPORTANT]
> **Backward Compatibility**: The default behavior changes. With `gospl_coupling_frequency = 1`, the behavior is identical to current (every step, no extrapolation needed). For N > 1, the **default** now applies erosion gradually instead of as a lump sum. Users who relied on the lump-sum behavior for specific reasons would need to be aware. This should be strictly better for all practical purposes.

> [!IMPORTANT]
> **Adaptive N bounds**: The adaptive scheme bounds N between 1 and 4× the user-specified value. The 4× upper bound is somewhat arbitrary. Should we make this configurable, or is 4× a reasonable cap?

> [!IMPORTANT]
> **Unit system**: The erosion rates are stored in m/s (SI), matching DES internal units. The conversion from goSPL's output (which works in years) happens during the coupling step. Please confirm this is consistent with how `var.dt` is defined (I see it's in seconds based on the `dt_years = var.dt / 3.1536e7` conversion at line 1598).

## Verification Plan

### Manual Verification

Since the coupling requires both DynEarthSol and goSPL to be properly set up with mesh files and configuration:

1. **Build test**: Run `make` to verify the code compiles without errors
2. **Regression test**: Run an existing goSPL-coupled simulation (e.g., using `test_gospl_coupling.cfg` or `gospl_driver/examples/core-complex-3d-with-gospl.cfg`) with `gospl_coupling_frequency = 1` to confirm behavior is unchanged (N=1 bypasses both extrapolation and adaptive logic)
3. **Gradual application test**: Run the same simulation with `gospl_coupling_frequency = 100` and compare output against the lump-sum behavior. The gradual scheme should produce smoother topography evolution
4. **Adaptive behavior test**: Run with `gospl_coupling_frequency = 50` and `gospl_rate_change_tolerance = 0.3`. Monitor the console output for `adaptive N=` messages to verify N adjusts up/down based on rate changes
5. **Diagnostic comparison**: Could you suggest an existing test case or setup that exercises the coupling, so we can define concrete pass/fail criteria?
