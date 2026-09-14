# PT Solver Improvements — Duretz et al. (2026)

Reference: Duretz, T., de Montserrat, A., Sevilla, R., Räss, L., Utkin, I., Spang, A.
"Automatic tuning of iterative pseudo-transient solvers for modeling the deformation
of heterogeneous media." *Geosci. Model Dev.* 19, 5343–5362 (2026).
https://doi.org/10.5194/gmd-19-5343-2026

## Motivation

The Räss et al. (2022) PT parameters — Re, CFL — are derived analytically for a
homogeneous uniform medium.  In DES they are fixed scalars set once before the PT
loop.  Duretz et al. (2026) show that for models with large viscosity contrasts or
evolving plasticity, *automatically adapting* the parameters during the iteration
gives 2–12× speed-ups over fixed-parameter PT.  The paper's two main tools are:

1. **Gershgorin circle theorem** — cheap upper bound on λ_max → sets Δτ (CFL).
2. **Rayleigh quotient** from successive iterate differences → lower bound on λ_min
   → sets damping (Re).

Both require re-evaluating parameters periodically during the iteration, not just once.

## Implementation plan

### Phase 1 — Periodic re-evaluation (DONE, this branch)

**What it does.**  `update_pt_params()` is called every `PT_retune_interval` PT
iterations (default 100) *inside* the PT loop, immediately after `update_stress_PT()`.
The existing function already contains a yield-aware μ_ve softening path that adjusts
the per-element damping when elements approach yield; re-evaluating it periodically
keeps the PT stepping factors in sync with the evolving plastic state.

The mesh is frozen during PT (no advection), so element heights and topology are
stable between retune calls; only stress and the derived μ_ve softening can change.

**New config parameter.**

```
control.PT_retune_interval = 100   # re-evaluate PT params every N iterations
                                   # 0 = disable (legacy behaviour)
```

**Cost.** `update_pt_params()` is O(nelem + nnode), roughly 1–2% of one iteration
cost at interval 100.

**Files changed.**
- `parameters.hpp` — added `int PT_retune_interval`
- `input.cxx` — registered param with default 100
- `examples/defaults.cfg` — documented default
- `dynearthsol.cxx` — retune call in both PT loops (init + main), after
  `update_stress_PT()` and before `update_force()`

### Phase 2 — Rayleigh quotient for adaptive Re (TODO)

**What it does.**  Every `PT_retune_interval` iterations, estimate the actual λ_min
of the assembled system from the Rayleigh quotient

$$
\lambda_\min \approx \frac{\sum_{i,j} \Delta v_{ij}\,\Delta f_{ij}}
                          {\sum_i \frac{1}{\texttt{PT\_dtau\_rho}[i]}\sum_j(\Delta v_{ij})^2}
$$

where Δv = v^k − v^{k−1} and Δf = f^k − f^{k−1}.  Use λ_min to update Re so the
damping tracks the actual slowest mode, not the analytical estimate for uniform material.

**Derivation needed before coding.**  The mapping Re → η (damping rate) in the
telegraph equation is:
  η_e = G̃Δτ_e / (G_e Δt) = Re·CFL·h_e·μ_ve_e / ((r+2)·L·G_e·Δt)

The Rayleigh quotient gives λ_min in units of [s⁻²] (eigenvalue of the assembled
stiffness / inertia).  Converting λ_min to the optimal Re requires matching these
units through:

  η* = 2√(λ_min · λ_max)   (geometric-mean critical damping)
  Re_opt = η* · (r+2) · L · G_e · Δt / (CFL · h_e · μ_ve_e)

The remaining issue is that Re is a scalar while the formula involves element-local
quantities; a global weighted average (h_mean, μ_ve_mean) is probably sufficient.

**New data structures needed.**
- `array_t* PT_vel_prev`  — velocity snapshot at last retune step (nnode × NDIMS)
- `array_t* PT_force_prev` — force_residual snapshot at last retune step (nnode × NDIMS)

**Files to change.**  `parameters.hpp`, `dynearthsol.cxx`, `geometry.cxx`
(`update_pt_params` or a thin `update_pt_params_adaptive` wrapper).

### Phase 3 — Gershgorin λ_max refinement (TODO)

**What it does.**  Replace or supplement the element-local Δτ estimate with a
Gershgorin-based λ_max that accounts for off-diagonal stiffness coupling across
viscosity-contrast boundaries.

**When to implement.**  Only if Phases 1–2 are insufficient for problems with >3–4
orders-of-magnitude viscosity contrast across adjacent elements.  The existing
per-element CFL is already a reasonable Gershgorin approximation for smooth fields.

**Files to change.**  `geometry.cxx` (new `compute_gershgorin_lambda_max()`),
`dynearthsol.cxx`.

## Testing strategy

1. **Regression (Phase 1):** `tests/functional/3d-evp-regular.cfg` —
   iteration counts and final stress must match the baseline within tolerance.
2. **Convergence improvement (Phase 1+2):** Gaussian weak-zone setup
   (`weakzone_option = 4`, high viscosity contrast) — plot residual vs. iteration
   with retune on vs. off; expect fewer iterations to tolerance.
3. **GoSPL coupling (Phase 1):** Short coupled run (`examples/coupled_N0200/`)
   with retune enabled — confirm no interaction with `copy_stress_PT` bookkeeping.
