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

## Implementation

### Phase 1 — Periodic re-evaluation (DONE, commit 3c6f15d)

**What it does.**  `update_pt_params()` is called every `PT_retune_interval` PT
iterations (default 100) *inside* the PT loop, immediately after `update_stress_PT()`.
The existing function already contains a yield-aware μ_ve softening path that adjusts
the per-element damping when elements approach yield; re-evaluating it periodically
keeps the PT stepping factors in sync with the evolving plastic state.

The mesh is frozen during PT (no advection), so element heights and topology are
stable between retune calls; only stress and the derived μ_ve softening can change.

**Convergence guard.**  The retune is skipped when
`l2_residual / force_scale < 100 × PT_relative_tolerance` (= 1e-4 with the default
tolerance of 1e-6).  This prevents disrupting a nearly-converged descent.

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

### Phase 2 — Rayleigh quotient for adaptive Re (DONE, commit a80bec5)

**What it does.**  Every `PT_retune_interval` iterations, `rayleigh_update_Re()`
estimates the actual λ_min of the assembled system from the Rayleigh quotient

$$
\lambda_\min \approx \frac{-\sum_{i,j} \Delta v_{ij}\,\Delta f_{ij}}
                          {\sum_i \frac{1}{\texttt{PT\_dtau\_rho}[i]}\sum_j(\Delta v_{ij})^2}
$$

where Δv = v_current − v_snapshot and Δf = f_current − f_snapshot over the last
`PT_retune_interval` iterations.  The critical-damping Reynolds number is then

$$
\text{Re}_\text{opt} = 2\sqrt{\lambda_\min \cdot \lambda_\max} \cdot
    \frac{(r+2)\,L\,\bar{G}\,\Delta t}{\text{CFL}\,\bar{h}\,\bar{\mu}_{ve}}
$$

with λ_max = CFL²/(r+2) from the stability limit, and with G̅, h̅, μ̅_ve
volume-weighted global means.  Re is clamped to [0.5, 2.0] × `PT_Re` to limit
overshoot on the first few calls before the iterate differences are representative.

After computing the quotient, the snapshots are updated to the current state for
the next retune call.

**Observed λ_min behavior (2D EP benchmark).**

- Early iterations (fast-mode transient): λ_min ≈ 0.1–0.2, Re increases to ~20–28.
- Steady convergence regime: λ_min ≈ 1e-3, Re_new ≈ 2 → clamped to 0.5 × Re_ref ≈ 7.5.
- Near the end of a long solve (dominant modes nearly exhausted): λ_min can jump
  above λ_max; Re is clamped at 2 × Re_ref ≈ 30.

The lower-bound clamping to 0.5 × Re_ref is the operative outcome in most cases:
Re is reduced from the default ~14.9, which gives LARGER dtau_rho (more aggressive
velocity stepping) at the cost of smaller Gdtau.  For problems where the default
Re is too conservative this adaptation enables convergence.

**Benchmark — 2D EP, pt-retune-demo, 10000 iter budget:**

| PT loop | Init residual | retune=0 iters | retune=100 iters |
|---------|--------------|----------------|------------------|
| init    | 6.4e9        | **10000 (FAIL)** | **5979** |
| step 1  | 8.3e10       | **10000 (FAIL)** | **5266** |
| step 2  | 8.0e10       | 2430           | 3378 (+39%) |
| step 3  | 8.5e10       | 908            | 810 (−11%) |
| step 4  | 8.5e10       | 1013           | 1055 (+4%) |
| **TOTAL** |             | **24351**       | **16488 (−32%)** |

The adaptive Re enables convergence of the init loop and step 1 that completely
fail (hit max_iter) without retuning.  Step 2 is 39% slower because its first
100-iteration window captures a fast-mode transient (200× residual drop), giving
λ_min ≈ 10 and boosting Re to the upper bound; this cost is outweighed by the
savings on the previously-failing loops.

**New data structures.**
- `array_t* PT_vel_prev` — velocity snapshot at last retune call (nnode × NDIMS)
- `array_t* PT_force_prev` — force snapshot at last retune call (nnode × NDIMS)
- `double PT_Re_adaptive` — current adaptive Re (reset to `PT_Re` at each PT loop start)
- `double PT_mu_ve_mean`, `PT_G_mean` — volume-weighted means for the Re formula

**Files changed.**
- `parameters.hpp` — new fields listed above
- `fields.cxx` — allocate snapshot arrays and `stress_old`; restore var.support accessor
- `geometry.cxx` — `rayleigh_update_Re()` function; `update_pt_params()` reads
  `PT_Re_adaptive` and accumulates volume-weighted means
- `geometry.hpp` — declaration of `rayleigh_update_Re()`
- `dynearthsol.cxx` — both PT loops: initialize `PT_Re_adaptive`, take initial
  snapshot, call `rayleigh_update_Re` then `update_pt_params` at retune points
- `input.cxx` — restored PT param registrations after merge

### Phase 3 — Gershgorin λ_max refinement (TODO)

**What it does.**  Replace or supplement the analytical λ_max = CFL²/(r+2) estimate
with a Gershgorin-based per-element bound that accounts for off-diagonal stiffness
coupling across viscosity-contrast boundaries.

**When to implement.**  Only if Phases 1–2 are insufficient for problems with >3–4
orders-of-magnitude viscosity contrast across adjacent elements.  The existing
analytical λ_max is already a reasonable approximation for smooth fields; and the
current [0.5, 2] × Re bounds limit how much a wrong λ_max can hurt.

**Files to change.**  `geometry.cxx` (new `compute_gershgorin_lambda_max()`),
`dynearthsol.cxx`.

## Testing strategy

1. **Regression (Phase 1+2):** `tests/functional/3d-evp-regular.cfg` —
   iteration counts and final stress must match the baseline within tolerance.
2. **Convergence improvement (Phase 1+2):** Gaussian weak-zone setup
   (`weakzone_option = 4`, high viscosity contrast) — plot residual vs. iteration
   with retune on vs. off; expect fewer iterations to tolerance.
3. **GoSPL coupling (Phase 1):** Short coupled run (`examples/coupled_N0200/`)
   with retune enabled — confirm no interaction with `copy_stress_PT` bookkeeping.
