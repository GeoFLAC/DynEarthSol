# PT boundary-pluck: symptom, analysis, and remedies

## Symptom

Running `core-complex-mmg-pt.cfg` (`has_PT = yes`) with thermal diffusion
disabled, the model is supposed to open a plastic shear zone in the
prescribed initial "weak zone" and accommodate the extension imposed via
`vbc_val_x0`/`vbc_val_x1` by transmitting that kinematics through the
domain. Instead, the boundary nodes carrying the prescribed velocity appear
to decouple from the interior — the imposed velocity shows up right at the
boundary but does not ramp smoothly across the domain the way a converged
elastostatic solution requires. Informally: the boundary gets "plucked
off" rather than driving deformation inward.

A quick check of `ccmmg.save.*.vtkhdf` (via `h5py`) supported this: the
x-velocity profile across the domain showed a real but heavily
under-shot gradient (edge values roughly an order of magnitude below the
nominal ±1e-9 m/s boundary condition), rather than the smooth,
domain-spanning linear ramp expected at equilibrium.

## Analysis

Both DR (dynamic relaxation) and PT reach the quasi-static balance the
same way: by letting a signal propagate from the boundary into the
interior over iterations, then damping it once it arrives. Neither method
solves the balance directly.

- In DR, that signal is a genuine elastic wave advancing in real,
  mass-scaled time.
- In PT, per the telegraph-equation derivation in `accelerated-PT.md`, the
  pseudo-time relaxation is *also* a damped wave equation,
  `∂²σ/∂τ² + η ∂σ/∂τ = V̂² ∂²σ/∂x²`, with the numerical modulus `G̃`
  standing in for real inertia. That second-order term is the
  "pseudo-inertia": a boundary change cannot be felt at the interior
  instantaneously — it must cross the mesh at the pseudo-wave speed `V̂`,
  one iteration at a time, just as a real elastic disturbance crosses
  physical space at the real wave speed, one CFL-limited step at a time.

Consequently, the number of PT iterations needed for a boundary condition
to reach the far side of the domain is bounded below by the same `L/h`
scaling derived in `accelerated-PT.md`
(`N ≈ ln(1/ε)·L/(π·CFL·h)`), regardless of whether the iteration variable
is called a "time step" (DR) or a "pseudo-time step" (PT). If the PT inner
loop exits — by hitting `PT_max_iter`, or, more relevant here, by
tripping the stagnation-detection early exit — before that many
iterations have elapsed, the resulting state is exactly what a
partially-propagated wave looks like: boundary nodes already relaxed onto
the imposed velocity, interior not yet reached. That is indistinguishable
from "the boundary got plucked off."

DR has no analogous early exit — it just keeps taking many small physical
steps, so the wave gets as many iterations as real time allows, amortized
over the whole run. PT buys its speed by taking one much larger physical
step and asking a single inner pseudo-time solve to do all the
propagating and damping that DR would spread across hundreds of small
steps — and that inner solve has a finite, apparently insufficient,
iteration budget here.

Three concrete, code-level gaps compound the effect, all inspected in the
current source:

1. **Stagnation-window floor uses the wrong length scale.**
   `pt_stagnation_window()` (`dynearthsol.cxx:541-548`) sizes its default
   window as `2·PT_L/(CFL·h_mean)`, using the *domain-average* element
   size (`var.PT_h_mean`). But `update_pt_params()` (`geometry.cxx:1086`)
   already computes local pseudo-timestep factors precisely because a
   graded mesh needs the disturbance to reach its *finest* elements,
   governed by `h_min`, not `h_mean` (see the comment at
   `geometry.cxx:1088-1091`). The weak zone is exactly the region MMG
   refines most aggressively (`remeshing_option = 11`), i.e. exactly
   where `h ≪ h_mean` — so the stagnation window is shortest precisely
   where the true domain-crossing requirement is longest.

2. **Stagnation test only watches a scalar residual, not whether the
   plastic front is still moving.** Two consecutive windows of <0.1%
   residual improvement trigger acceptance (`stag_strikes >= 2`,
   `dynearthsol.cxx:910-923`). The code's own comment
   (`dynearthsol.cxx:874-877`) already acknowledges that the nonsmooth
   plastic return-mapping projection can floor the residual above
   tolerance — a plateaued residual does not necessarily mean nothing is
   still propagating.

3. **η/Δτ tuning uses elastic stiffness, not the actual (yield-softened)
   tangent stiffness.** In `update_pt_params()`, `mu_ve_e` is built from
   `G_e = var.mat->shearm(e)` (the elastic shear modulus, unaffected by
   yielding) and the current viscosity — never from a
   plasticity-consistent tangent modulus. Once an element yields, its
   effective stiffness for propagating a disturbance is lower than what
   its damping was tuned against, biasing exactly the weak-zone
   iterations toward appearing "resolved" before they physically are.

## Follow-up: solutions 1/4 tested, found insufficient

Setting `PT_stagnation_window = 3000` explicitly (the config-only version
of fix 1) was tried. Result: the "stagnated at iter N" message moved to a
much later `N` (e.g. iter 9000 instead of ~2000-4500), confirming the
setting took effect — but the residual did **not** trend toward
tolerance. It oscillated in a bounded band. Concretely, from a step-7
log with the widened window:

```
iter 3400: 3.119e12   iter 3700: 3.125e12   iter 4000: 3.124e12   iter 4300: 3.129e12   iter 4600: 3.137e12
iter 3500: 2.769e12   iter 3800: 2.742e12   iter 4100: 2.729e12   iter 4400: 2.730e12   iter 4700: 2.721e12
iter 3600: 2.191e12   iter 3900: 2.135e12   iter 4200: 2.133e12   iter 4500: 2.136e12   iter 4800: 2.138e12
```

This is a stable **limit cycle with a ~300-iteration period**, not slow
convergence — it held through 9000 iterations with no downward trend.
Widening the stagnation window only changes how long the loop waits
before giving up; it cannot fix a state that will never converge at the
current damping. **Solutions 1 and 4 are therefore ruled out as fixes**
(they remain useful only as a compute-saving/diagnostic tool, since they
reveal whether a step is oscillating vs. still genuinely converging).

A residual limit cycle like this is the signature of a fixed-point
iteration bouncing off a nonsmooth yield-surface switch: the state
overshoots into (or out of) yielding, the pseudo-relaxation corrects,
overshoots back, repeat. Two mechanisms in the source can produce
exactly this:

- Gap 3 (elastic-only `mu_ve_e`, above): a newly-yielded element is
  damped as if still fully elastic — too stiff a correction for the
  actually-softened material.
- **The plastic-weakening step-size cap lags by one step.**
  `compute_dt_PT`'s `dt_weakening` limit (`geometry.cxx:1028-1034`)
  estimates the allowed step from *last step's* `delta_plstrain/dt`. An
  element yielding for the first time has `delta_plstrain = 0` from the
  previous step, so it contributes nothing to the cap — the very step
  that first activates the weak zone can take an oversized, uncapped
  leap through a large, brand-new plastic disturbance in one shot.

### Diagnostic test: does a smaller step avoid the cycle?

To isolate step-size as the trigger (independent of any code change),
`core-complex-mmg-pt-smalldt.cfg` was run with `dt_fraction = 0.1`
(all else identical to `core-complex-mmg-pt.cfg`, `PT_stagnation_window =
3000`), producing steps roughly 10x smaller (`dt ≈ 385.6 yr` vs. the
original `dt ≈ 3855.8 yr`).

Result: every step through the equivalent activation point (step 6→7,
`time ≈ 2313-2699 yr`) converged with a **smooth, strictly monotonic**
residual decay toward a plateau — no oscillation:

```
iter   300: residual/residual_0 = 0.001463
iter  1000: residual/residual_0 = 0.001365
iter  5000: residual/residual_0 = 0.001312
iter 10000: residual/residual_0 = 0.001306
iter 18000: residual/residual_0 = 0.001304   (stagnated, still monotonic)
```

Still above tolerance (so `PT_max_iter`/stagnation logic is still doing
real work), but qualitatively a different, healthier regime than the
limit cycle at full step size. This is a fairly direct confirmation that
the oscillation is driven by the size of the per-step plastic increment,
not an inherent property of the PT scheme or an artifact of insufficient
iteration budget in general.

(This run later hit an unrelated MMG3D remeshing failure —
`BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH` — a separate meshing-library
robustness issue, out of scope for this note.)

## Potential solutions

Given the above, the priority is revised: 1 and 4 are demoted to
diagnostic tools; the fixes that actually address a non-convergent limit
cycle are 2 (partial: saves compute, doesn't fix physics), 3, and a new
item 5. None requires reverting to DR-scale small time-stepping across
the whole run — the goal is to spend extra iterations/smaller steps only
on the few steps that newly activate plastic yielding, leaving the many
already-near-equilibrium steps to convergence quickly via warm-starting.

1. ~~Fix the stagnation-window length scale~~ (`var.PT_h_min` instead of
   `var.PT_h_mean` in `pt_stagnation_window()`). **Tested via config
   (solution 4) and found insufficient** — see follow-up above. Keep only
   as a way to let a genuinely-converging step run long enough to prove
   it either converges or cycles.

2. **Make the stagnation break physics-aware.** Detect the limit cycle
   directly (e.g. track residual variance/periodicity over a window, or
   whether the actively-yielding element set is still changing) instead
   of just "no improvement." Would let the loop bail out in ~1000
   iterations instead of burning thousands sampling the same cycle — a
   compute-saving fix, not a convergence fix.

3. **Feed plastic softening into `mu_ve_e`.** Couple a yield-aware
   tangent modulus into `update_pt_params()`'s `mu_ve_e` formula so local
   damping matches the material's actual (softened) post-yield stiffness.
   Directly targets the overshoot mechanism behind the cycle.

5. ~~Fix the one-step lag in `dt_weakening`~~ **— tried, found
   insufficient; see negative result below.** Do not re-attempt the
   "persist the last known rate across a gap" version of this fix without
   also addressing the true-first-activation case described there.

## Negative result: solution 5 (persisted-rate fallback) does not fix it

**What was implemented.** A new `Variables` field,
`PT_steps_min_prev` (`parameters.hpp`), persisted the most recent finite
`steps_min` computed in `compute_dt_PT()` (`geometry.cxx`) across calls.
When a call found no element with a fresh `delta_plstrain > 0` (i.e.
`steps_min` stayed at its "unconstrained" sentinel), the persisted value
was substituted instead of allowing an uncapped `dt_weakening`:

```cpp
double steps_min_effective = steps_min;
if (check_weakening) {
    if (var.PT_steps_min_prev > 0)
        steps_min_effective = std::min(steps_min_effective, var.PT_steps_min_prev);
    if (steps_min < 1e300)
        var.PT_steps_min_prev = steps_min;
}
double dt_weakening = (check_weakening && steps_min_effective < 1e300) ?
    steps_min_effective * var.dt : std::numeric_limits<double>::max();
```

**Test.** Rebuilt and reran the original full-step config
(`core-complex-mmg-pt.cfg`, `dt_fraction = 1.0`, `PT_stagnation_window =
3000`) — the exact case that previously produced the step-7 limit cycle.

**Result: no effect whatsoever.** Step 7's output was bit-for-bit
identical before and after the change:

```
Before: PT stagnated at iter 9000: residual = 3523569472439.132812 (residual/residual_0 = 0.003629)
After:  PT stagnated at iter 9000: residual = 3523569472439.132812 (residual/residual_0 = 0.003629)
```

Same iteration count, same residual to the last digit — `dt` for step 7
was computed identically with and without the fallback in place.

**Why it did nothing.** The fallback only helps a *reactivation after a
gap*: some element must have reported a finite rate at least once before,
so there's something to persist through a later gap. Steps 1-6 in this
model never had any element with `delta_plstrain > 0` — step 7 is a
**genuine first-ever activation across the whole domain**, with zero rate
history anywhere to persist. `var.PT_steps_min_prev` never became finite
before the step that needed it, so the substitution was a no-op exactly
in the one case this note cares about.

**Takeaway for future attempts.** Any fix framed as "reuse/extrapolate a
previously observed rate" cannot help a true first-ever activation — there
is nothing to extrapolate from. A fix here needs a mechanism that doesn't
depend on *any* prior observed rate, e.g.:
- comparing current stress against the yield surface directly (a
  distance-to-yield estimate rather than a rate estimate), or
- a generic conservative fallback cap applied whenever `check_weakening`
  is on but no rate has *ever* been recorded domain-wide (not just this
  step) — arbitrary/untuned, but at least closes the true-cold-start gap
  that the persisted-rate approach cannot.

Neither was implemented as part of this investigation. The code changes
for this attempt were fully reverted (`git diff`/`git status` clean on
`parameters.hpp` and `geometry.cxx`); `compute_dt_PT()` and `Variables`
are back to their pre-solution-5 state.

## Negative result: softening-curve shape (`pls0` plateau) is not the cause

**Idea.** `mat.pls0` (currently `0` in `core-complex-mmg-pt.cfg`) sets the
width of a flat, unweakened plateau before the cohesion/friction ramp
begins — schematically `‾\_` (flat at `cohesion0`/`friction_angle0` for
`plstrain < pls0`, linear ramp down to `pls1`, flat at `cohesion1`
beyond) instead of `\_` (ramp starting immediately at `plstrain = 0`).
The hypothesis: with `pls0 = 0`, the very first plastic increment already
starts shrinking the yield surface, so the first (large, PT-sized)
activation step has to solve a return-mapping problem against a target
that moves *within that same step* — a harder, less-standard nonlinearity
than fixed-yield-surface (non-softening) plasticity. A nonzero `pls0`
would let the first activation step(s) be pure elastic-perfectly-plastic
(fixed yield surface), which is generally well-posed, and would also
seed real `delta_plstrain` rate history before softening ever starts —
closing the exact gap that made the solution-5 negative result above
unfixable (no rate history at the moment it's needed).

**Tests run**, all at the original full step size (`dt_fraction = 1.0`,
`PT_stagnation_window = 3000`), varying only `mat.pls0` (config-only,
no rebuild):

| `pls0` | step 6→7 stagnation | character |
|---|---|---|
| `0` (baseline) | iter 9000, residual/residual_0 = 0.003629 | ~300-iter limit cycle |
| `0.05` | iter 15000, residual/residual_0 = 0.003719 | still oscillating, e.g. 0.00356→0.00527→0.00318 across iters 1000-11000 |
| `0.2` | iter 9000, residual/residual_0 = 0.002755 | still oscillating, e.g. 0.00266→0.00304→0.00276 across iters 1000-9000 |

Raising `pls0` 4x (`0.05 → 0.2`, i.e. up to 40% of `pls1 = 0.5`) barely
moved the stagnation magnitude and did not change the oscillating
character at all. If "the step blows through the whole plateau into the
ramp in one shot" were the operative mechanism, a much wider plateau
should have shown at least some improvement — it didn't.

**Decisive follow-up: disable softening entirely.** Set
`cohesion1 = cohesion0 = 4.4e7` (`friction_angle1` was already equal to
`friction_angle0 = 30` in this config), removing the weakening curve
altogether — the material is now purely elastic-perfectly-plastic, yield
surface never shrinks, at any `plstrain`. Same full step size, same
weak-zone IC (`weakzone_plstrain = 0.5`, now without any accompanying
strength drop).

**Result: the oscillation is still there, same magnitude.** Step 6→7
stagnated at iter 9000, residual/residual_0 = 0.002660 — statistically
indistinguishable from both `pls0` tests above — and the iteration trace
is the same bounded, non-converging fluctuation (~0.0026-0.0036), not
monotonic convergence:

```
iter 1000: 0.002931     iter 4000: 0.002875     iter 7000: 0.002925
iter 2000: 0.002606     iter 5000: 0.002903     iter 8000: 0.002956
iter 3000: 0.002625     iter 6000: 0.002912     iter 9000: 0.002660
```

**Conclusion.** The softening-curve shape (`pls0`/`pls1` ramp) is **ruled
out** as the cause of the limit cycle. The instability appears even when
the yield surface never shrinks at all, so it cannot be about a moving
target during the return-mapping solve. It is much more likely driven by
the nonsmooth elastic→plastic *onset* itself — the switching from
unyielded to yielded response — combined with gap 3 above (`mu_ve_e` in
`update_pt_params()` damped from the *elastic* modulus `G_e`, regardless
of whether the element has ever yielded). That mismatch would bite on
first yield whether or not any softening ever follows.

**Takeaway for future attempts.** Do not re-attempt any variant of
"delay/reshape the weakening curve" (wider `pls0` plateau, smoother ramp,
etc.) expecting it to fix this oscillation — it was tested directly,
including the limiting case of no softening at all, and made no
difference. The remaining, still-untested candidate is solution 3
(couple a yield-aware, not purely elastic, tangent modulus into
`mu_ve_e`) — that is the next thing to try.

All config changes made for this test (`pls0`, `cohesion1`) were
diagnostic only, made directly in `examples/core-complex-mmg-pt.cfg`;
restore `pls0 = 0` and `cohesion1 = 4e6` before using that file for
anything other than reproducing this note.

## Positive (partial) result: solution 3, then damping the plastic return-map too

### Solution 3: yield-aware `mu_ve_e`

**What was implemented** (`geometry.cxx`, `update_pt_params()`). Before
computing `mu_ve_e`, each element's elastic shear modulus `G_e` is
softened when its current stress is close to or past yield:

```cpp
if (has_plastic) {
    ConstTensorAccessor s = (*var.stress)[e];
    double mean = 0;
    for (int d = 0; d < NDIMS; ++d) mean += s[d];
    mean /= NDIMS;
    double J2 = 0;
    for (int d = 0; d < NDIMS; ++d) { double dev = s[d] - mean; J2 += 0.5 * dev * dev; }
    for (int d = NDIMS; d < NSTR; ++d) J2 += s[d] * s[d];
    double tau_eq = std::sqrt(J2);

    double amc, anphi, anpsi, hardn, ten_max;
    var.mat->plastic_props(e, (*var.plstrain)[e], amc, anphi, anpsi, hardn, ten_max);
    double yield_ratio = (amc > 0) ? tau_eq / amc : 0.0;

    if (yield_ratio > yield_onset) {   // yield_onset = 0.7, yield_floor = 0.05
        double t = std::min((yield_ratio - yield_onset) / (1.0 - yield_onset), 1.0);
        G_e *= 1.0 - t * (1.0 - yield_floor);
    }
}
```

`yield_ratio` is a deviatoric-stress-only (friction/pressure term
ignored) proxy for "how close to Mohr-Coulomb failure" — evaluated fresh
from the current stress every call, so (unlike solution 5) it has no
cold-start blind spot; it doesn't need any history of prior yielding.

**Result (`pls0 = 0`, `cohesion1 = 4e6` restored, `dt_fraction = 1.0`):**
step 6→7 showed **real, partial convergence** instead of pure oscillation
— residual dropped from 0.00473 (iter 900) to a minimum of ~0.00284
(iter 15000), then drifted back up to 0.00319 by the iter-24000
stagnation. Genuine ~40% improvement over the initial value, but not
monotonic all the way and still short of tolerance. Earlier steps were
also visibly better than baseline (e.g. step 4: 0.001534 vs. baseline's
0.001895).

**Why only partial.** Rereading `update_stress_PT()` while implementing
this revealed the more precise mechanism: each iteration relaxes stress
toward the elastic target at a *damped* rate `theta` (built from
`mu_ve_e`), but then — regardless of `theta` — applies a **full,
undamped** Mohr-Coulomb return-map correction every single iteration if
that relaxed stress violates yield (`elasto_plastic()`/
`elasto_plastic2d()`, called with zero strain increment). Softening
`mu_ve_e` only slows the *creep toward* yield; it does nothing to the
*snap-back itself*, which stays instantaneous.

### Damping the return-map correction itself

**What was implemented** (`rheology.cxx`, `update_stress_PT()`). The
fully-projected stress from `elasto_plastic()`/`elasto_plastic2d()` is
now blended back toward its pre-projection value using the *same*
`w0`/`w1` weights already used for the elastic relaxation immediately
above it in the same loop iteration, instead of being applied in full:

```cpp
double s_pre[NSTR];
for (int i = 0; i < NSTR; ++i) s_pre[i] = s[i];

if (var.mat->is_plane_strain) {
    double& syy = (*var.stressyy)[e];
    double syy_pre = syy;
    elasto_plastic2d(..., s, syy, ...);
    syy = w0 * syy_pre + w1 * syy;
}
else {
    elasto_plastic(..., s, ...);
}
for (int i = 0; i < NSTR; ++i)
    s[i] = w0 * s_pre[i] + w1 * s[i];
```

This was tested **combined with** solution 3 (both left active together
— they address complementary halves of the same mechanism).

**Result:** a genuine qualitative change. Step 6→7:

```
iter     0: 0.12830   (initial)
iter  2000: 0.01274
iter  6000: 0.01077
iter 12000: 0.01056
iter 18000: 0.01049
iter 24000: 0.01047
iter 26000: 0.01047   ← stagnated (2-strike "no improvement")
```

**The oscillation is gone.** This is smooth, strictly monotonic
convergence toward a tight plateau — the same healthy character the
small-`dt` diagnostic showed, not the ~30-70%-swing limit cycle every
prior full-step test produced. But the plateau itself
(`residual/residual_0 ≈ 0.0105`) is *higher* than the oscillating
baseline's floor (~0.0036) or solution-3-alone's floor (~0.0032).
Earlier steps (2-4) improved further under the combination (e.g. step 2:
0.000085 vs. baseline's 0.000627) — the combination is not uniformly
worse, just worse at this specific step's plateau value.

The tail of the trace above is still (very slowly) decreasing when the
2-consecutive-window stagnation check fires (7122935661769 →
7120316510430 → 7122298924266 — sub-0.01% wiggles, not yet flat) — so it
is not confirmed whether this trajectory is truly stuck at ~0.0105 or
would keep closing the gap toward tolerance given a wider stagnation
window. That is the immediate next thing to check, before concluding
anything about the final floor value.

**Status:** both changes are left in place in the working tree
(`geometry.cxx`, `rheology.cxx`) pending that follow-up test. This is the
first result in this investigation that changes the qualitative
character of the failure (oscillation → convergence) rather than just
shifting a stagnation-iteration count or a residual magnitude within the
same oscillating regime.

(Run hit a different MMG3D crash after this step —
`MMG3D_consistency_error_message` — same unrelated meshing-library
territory as the earlier crashes in this document.)

## Major finding: `dt` was never actually being recomputed per step

While trying to widen `PT_stagnation_window` further and testing a new
`dt_yield` limiter (below), a much more fundamental issue surfaced:
**every diagnostic run in this entire investigation up to this point had
`dt` frozen at its initial, pre-deformation value.**

`dynearthsol.cxx:1006-1019` batches `dt` recomputation together with
genuinely expensive, slowly-varying physics:

```cpp
const int slow_updates_interval = 10;
if (var.steps % slow_updates_interval == 0) {
    phase_changes(param, var);
    if (param.control.has_hydration_processes) advect_hydrous_markers(...);
    var.dt = (param.control.has_PT) ? compute_dt_PT(param, var) : compute_dt(param, var);
}
```

`var.steps % 10 == 0` is only true at step 0. Every test run this session
crashed (via the recurring, unrelated MMG3D issue) well before reaching
`var.steps = 10`, so `dt` was computed exactly once, from the pristine,
undeformed initial state, and held fixed for the entire step 1-9 range
every test in this document was conducted over. Confirmed directly: an
unconditional `fprintf` at `compute_dt_PT()`'s entry/return showed it was
only ever called twice total (both at `var.steps=0`) in a run that
reached "Step = 6" — and the constant `dt = 3855.79` yr seen identically
across *every* test in this document (baseline, `pls0` sweep,
no-softening, solution 3, damped correction, wide stagnation window) is
the direct fingerprint of this.

**Consequence:** `dt_weakening` (pre-existing) and `dt_yield` (below,
new) both compute fresh, history-free per-step limits from the *current*
stress/strain-rate/plastic-strain state — but that state was never
actually current. Every conclusion in this document about "the step 6→7
limit cycle" is still valid (it's about the PT *inner iteration*'s
behavior for a given, fixed `dt`), but none of the dt-*selection*-side
experiments (existing `dt_weakening`, and the `dt_yield` attempt below)
had a fair test until this was fixed.

**Fix applied:** moved `var.dt = ...` out of the `slow_updates_interval`
gate so it recomputes every step, leaving `phase_changes()` /
`advect_hydrous_markers()` on their existing 10-step cadence:

```cpp
if (var.steps % slow_updates_interval == 0) {
    phase_changes(param, var);
    if (param.control.has_hydration_processes) advect_hydrous_markers(...);
}
var.dt = (param.control.has_PT) ? compute_dt_PT(param, var) : compute_dt(param, var);
```

This fix is kept (validated below as stable on its own, independent of
`dt_yield`).

### Related, independent bug found and fixed: `var.max_global_vel_mag` stuck at 0 in PT mode

While investigating the above, the user noticed the printed `vmax` was
*always* exactly `0.00000e+00` in every PT run. Root cause:
`compute_dt()` (DR mode) sets `var.max_global_vel_mag = global_max_vem;`
(`geometry.cxx:885`); `compute_dt_PT()` computed the equivalent local
`global_max_vem` for its own `dt_advection` but never wrote it back to
`var.max_global_vel_mag` — so the field stayed at its zero-initialized
value for the entire run, in every PT test in this document.

This isn't just a display bug: `compute_mass()`'s `pseudo_speed_ATP`
(`geometry.cxx:1307`) is built directly from `var.max_global_vel_mag`,
and is used (in place of the boundary-velocity-based `pseudo_speed`)
whenever `control.use_global_velocity_scaling` is true. Had that flag
been enabled here, `apprent_speed = min(pseudo_speed_ATP, ...) = min(0,
...) = 0` would have produced `rho = bulkm/0 = inf` in the mass/inertia
scaling. It wasn't enabled in this config (default `false`; only
auto-forced true for rate-and-state-friction rheology, which this model
doesn't use), so it didn't corrupt these particular runs — but it's a
live hazard for any PT config that does use
`use_global_velocity_scaling`.

**Fix applied** (`geometry.cxx`, `compute_dt_PT()`): added
`var.max_global_vel_mag = global_max_vem;` right after `global_max_vem`
is finalized, mirroring `compute_dt()`. Kept permanently — this is an
unambiguous, independent correctness fix, not an experimental change.

## Negative result: `dt_yield`, a time-to-yield dt limiter

With per-step `dt` recomputation now real, a further idea (motivated by
the earlier "why not scale PT's `dt` off DR's `dt_elastic`" discussion)
was tried: cap `dt` directly from how close each element's current
elastic stress state is to yielding, using only the *current*
stress/strain-rate (no history, unlike `dt_weakening`):

```cpp
dt_yield_e = yield_f * headroom / (2 * G_e * eps_eq)
```

where `headroom = amc - tau_eq` (remaining distance to the deviatoric,
friction-ignored Mohr-Coulomb strength) and `eps_eq` is the equivalent
strain rate. New config knob `control.PT_yield_approach_fraction`
(default `0.1`), implemented in `geometry.cxx`'s `compute_dt_PT()` plus a
small shared `deviatoric_equiv()` helper.

**Result: a genuine, self-inflicted runaway**, not the intended
stabilization. With `dt_yield` active (and the per-step `dt` fix and
`max_global_vel_mag` fix both in place so this was a fair test), `vmax`
exploded:

```
step 1: vmax = 7.50000e-10 m/s   (sensible, matches boundary condition scale)
step 2: vmax = 2.60068e-07 m/s   (~350x jump)
step 3: vmax = 8.99454e-04 m/s   (~3500x jump again)
```

Isolated by re-running with `PT_yield_approach_fraction = 0` (same
per-step `dt` fix, same `max_global_vel_mag` fix): `vmax` stayed sensible
and `dt_advection`/`dt_weakening` varied smoothly, confirming `dt_yield`
specifically — not the per-step `dt` fix, not the `vmax` fix — was the
trigger.

**Mechanism (not a formula bug — a real structural interaction).**
`update_pt_params()`'s `mu_ve_e = 1/(1/(G·dt) + 1/μ_s)` (Räss et al. 2022,
Eq. 35) is *supposed* to depend on `dt`: in the elastic-dominated regime
here, `mu_ve_e ≈ G·dt`, which makes `θ_e = G̃Δτ_e/(G·dt)` come out
`dt`-independent (by design — PT's relative convergence rate shouldn't
depend on what physical `dt` was chosen). But `PT_dtau_rho` (the
per-node pseudo-timestep driving `update_velocity_PT()`,
`vel[i] += dtau_rho[i]·force[i]`) scales as `~1/mu_ve_e`, i.e. `~1/dt` —
also by design, to keep the *velocity* side of the iteration in step with
an elastic target that itself shrinks with `dt`. That `1/dt` scaling is
correct and shouldn't be "fixed."

The likely actual mechanism: this problem's PT solves chronically don't
reach true convergence (documented throughout this file — every step is
accepted via the stagnation heuristic well above the `1e-6` tolerance,
even after solution 3 + the damped return-map). At the original, larger
`dt`, that leftover fractional error produces a modest absolute velocity
error. Once `dt_yield` shrinks `dt`, `dtau_rho` grows as `1/dt` to
compensate (correctly, per the theory) — but the *same* fractional
non-convergence, multiplied by a larger `dtau_rho`, produces a larger
absolute velocity error at the same "stagnated" exit point. That inflated,
still-not-equilibrated velocity feeds the *next* step's `dt_yield`
strain-rate estimate, shrinking `dt` further, growing `dtau_rho` further
— a compounding loop that only requires the pre-existing convergence gap
to exist, not any flaw in `mu_ve_e`/`dtau_rho`.

**Takeaway for future attempts.** Don't try to patch the `1/dt` scaling
in `mu_ve_e`/`PT_dtau_rho` — it's theoretically load-bearing (Räss et
al.'s design for `dt`-independent relative convergence), not a defect.
Any new `dt`-selection criterion that reads back the PT solve's own
output (stress, strain rate) to decide how much to shrink `dt` is only
safe if that output is trustworthy (near-equilibrium) to begin with — it
isn't, yet, in this model. Fix PT's convergence-quality problem first
(continue from solution 3 / the damped return-map); revisit a
time-to-yield-style `dt` limiter only after that's solid.

**Status:** `dt_yield` and `control.PT_yield_approach_fraction` were
fully reverted (`geometry.cxx`, `parameters.hpp`, `input.cxx`,
`examples/core-complex-mmg-pt.cfg`). The per-step `dt` recomputation fix
and the `max_global_vel_mag` fix were both kept — both are validated,
independent corrections, confirmed stable on their own with `dt_yield`
disabled.

## Staggered elastic-solve / return-map scheme

Motivated by the question: if the elastic solve is allowed to fully
equilibrate *before* any plastic return-map correction is applied (rather
than interleaving a damped correction into every PT iteration, as the
existing code does), does the boundary-pluck disappear? Implemented as
`control.PT_staggered` (default off) + `control.PT_max_outer_iter`:

- Outer loop, only active when `rheol_type & rh_plastic` and not RSF:
  1. Run the normal PT inner loop with `update_stress_PT(..., apply_correction=false)`
     — pure elastic relaxation toward the Maxwell target, no return-map at
     all — until it converges or stagnates.
  2. Apply one full, undamped Mohr-Coulomb correction to every element via
     a new `apply_return_map_PT()` (no `s_pre`/blend — a real, not
     damped, projection).
  3. If nothing yielded, stop. Otherwise re-solve the elastic problem
     (step 1 again) from the corrected stress state, repeat up to
     `PT_max_outer_iter` rounds.
- Each outer round gets a *fresh* `PT_max_iter` budget — it's a distinct
  elastic re-solve from a newly-corrected state, not a continuation.

**Bug found and fixed: `stress_old` re-anchoring.** Initially every outer
round after the first stagnated at the same ~174-element/round residual
plateau. Root cause: `update_stress_PT`'s elastic target is built from
`var.stress_old` (the stress at the *start of the physical step*), which
was never updated between outer rounds — so round 2's "elastic solve" kept
relaxing toward a target anchored to the pre-correction state, undoing
round 1's correction. Fix: `copy_stress_PT(*var.stress, *var.stress_old)`
after each correction round, so each new elastic re-solve is anchored to
the just-corrected state. (Checked with the user first whether `target`
should depend on pseudo-time `dτ` rather than physical `dt` — confirmed
`target`'s `dt`-dependence is intentional, representing the elastic
predictor for the full physical step; `dτ` only controls the `w0`/`w1`
blend rate. Not the bug — a sanity check before the real fix.)

With this fix, the staggered scheme achieves genuine (non-stagnation)
convergence every outer round tested — a first in this investigation.
**But the underlying boundary-pluck symptom was unchanged**: even at
genuine convergence, `ccmmg-pt.save.000001.vtkhdf` (steps=3) still showed
only 44 of 3018 nodes with any meaningful velocity, all exactly at the two
boundary faces.

## Force-scale / gravity discovery (and a fix that didn't resolve the pluck)

Suspected the residual-convergence check itself was masking the problem:
`calculate_characteristic_force()` normalizes the residual against the sum
of `|tr|` (raw nodal force contributions, non-cancelling), which includes
the lithostatic background stress and gravity — both enormous and
present even in perfect equilibrium. A `force_scale` dominated by this
background could let the relative-residual test pass while the *actual*
unbalanced force (the part that should drive deformation) is still huge
in absolute terms but tiny by comparison.

First attempt (subtracting only the explicit `buoy` term) changed
`force_scale` by <1% — showing gravity's explicit addition isn't the
dominant term; the *locally-balancing* background lithostatic stress is.

Second, more general fix: `calculate_characteristic_force()` now takes a
lazy, one-time snapshot of `tmp_result` into a new `var.tmp_result_ref`
(`Variables` member, `fields.hpp`/`fields.cxx`) on its first call ever, and
computes the RMS of `|tr − tr_ref|` (element-wise delta from that
reference) instead of raw `|tr|`. This removes the constant background
contribution from the scale. (Known limitation, accepted: `tmp_result_ref`
isn't reallocated/interpolated across remeshing; not triggered in these
test runs.)

**Result: changed the convergence bookkeeping (more outer rounds needed,
different residual magnitudes) but did not change the converged end
state.** `ccmmg-pt.save.000001.vtkhdf`, re-checked after this fix, still
showed the identical 44/3018-node boundary-only pattern. **Conclusion:
the pluck is not a residual-normalization artifact** — some other
mechanism is preventing interior velocity from developing even when the
elastic sub-solve reports genuine convergence.

## Decisive control experiment: does a genuinely yield-free elasto-plastic run still pluck?

User's proposal: keep the full elasto-plastic code structure (staggered
wrapper, yield-aware `mu_ve_e`, return-map machinery all active/gated on
`rh_plastic`) but set material strength so high that yielding can never
actually trigger. If the pluck persists, it's caused by something in the
staggered/PT machinery itself, independent of plasticity. If it
disappears, actual yielding is the cause.

New config `examples/core-complex-mmg-pt-noyield.cfg`: identical to
`core-complex-mmg-pt.cfg` (`rheology_type = elasto-plastic`,
`PT_staggered = yes`) but `cohesion0 = cohesion1 = 1e32`.

**First attempt was contaminated.** Raising only `cohesion0`/`cohesion1`
left `mat.max_tension` (`input.cxx:821`) at its default `1e9` Pa. Cohesion
failure was suppressed, but the kinematically-imposed boundary velocity
creates a local stress concentration that easily exceeds a 1 GPa *tensile*
cutoff — so real return-map correction was still happening via the
tension-failure mode, concentrated at/near both boundaries (confirmed:
~171 elements with plastic strain *not* equal to the exact weak-zone IC
value of `0.5`, spanning from each boundary inward). This produced a
velocity field that looked mostly smooth in the interior but still showed
an abrupt ~50% drop right at the boundary elements — caught by visual
inspection in ParaView (plastic strain clearly nonzero at the boundary),
not by the coarse node-count metric used earlier, which only checked
`|vx|` thresholds and missed the small-scale boundary drop.

**Corrected run:** also set `max_tension = 1e32`. Result: **zero elements
deviate from the initial condition** (every nonzero-`pls` element is
exactly the untouched `0.5` weak-zone IC value — genuinely no yielding
anywhere). Velocity field at `ccmmg-pt-noyield.save.000001.vtkhdf`
(steps=3):

```
x = -365 to 9708:     vx_mean = -9.05e-10
x = 9708 to 19781:    vx_mean = -7.00e-10
...
x = 90292 to 100365:  vx_mean = +8.76e-10

near left boundary (within 2 km):   vx ≈ -1.0e-9   (= exact BC value)
just inside (2–6 km from boundary): vx ≈ -8.5e-10 to -9.3e-10  (smooth, no cliff)
```

A clean, monotonic, domain-spanning ramp with no abrupt boundary drop —
matching the pure-elastic (`rheology_type = elastic`) reference behavior.
2385/3018 nodes have `|vx| > 1e-10` (vs. 44/3018 in every plastic run with
real softening enabled).

**Conclusion: the boundary-pluck is caused specifically by actual plastic
yielding/return-mapping, not by the staggered wrapper, the yield-aware
`mu_ve_e` code, `dt`, or `force_scale` normalization** — all of that
machinery runs identically in this test and produces a normal result when
yielding doesn't actually trigger. This confirms (on much firmer footing
than the earlier "solution 3" / damped-correction evidence) that the
mechanism lives specifically in how the return-map correction interacts
with the (elastic or staggered) velocity/force update — most likely the
same limit-cycle-style interaction documented above (damped elastic
relaxation vs. an undamped correction), now shown to hold even in the
staggered architecture where the correction is no longer *interleaved*
per-iteration but applied once per fully-converged outer round.

**Open question, next to investigate:** why does the same elastic-solve
machinery that spreads velocity smoothly across the whole domain (this
control, and the standalone pure-elastic `.cfg`) fail to do so once
real return-map corrections start accumulating at yielded elements —
i.e., what specifically about a corrected stress state at a subset of
elements suppresses net velocity development in the *rest* of the domain,
rather than just locally at the yielded elements themselves.

**Takeaway for future attempts.** When constructing a "disable yielding"
control for this material model, remember there are *two* independent
failure surfaces gated by *two* independent parameters — shear/cohesion
failure (`cohesion0`/`cohesion1`/`friction_angle*`) and tensile failure
(`mat.max_tension`, default `1e9` Pa, easy to miss since it's not part of
the `[mat]` cohesion/friction block most edits focus on). Raising only
the cohesion terms is not sufficient to fully suppress yielding in a
boundary-velocity-driven setup — check plastic-strain-vs-IC directly
(not just velocity thresholds) to confirm a "no-yield" control is
actually yield-free before trusting its result.

## `dt_bc_yield`: a boundary-yield-aware `dt` limiter (per-face)

Motivated by the question the yield-free control raised: *in the staggered
scheme, why do the boundary elements yield first at all?* Mechanism: every
PT iteration `apply_vbcs()` pins the boundary node velocities to their full
imposed value. Before the deformation has genuinely propagated across the
domain, a boundary-adjacent element sees an artificially concentrated local
strain rate (`~v_bc/h_min` rather than the eventual equilibrium `~v_bc/L`),
which drives it past yield transiently — a numerical, not physical,
first-yield at the boundary. That premature return-mapping is exactly what
suppresses inward transmission.

Fix idea (distinct from the reverted `dt_yield`, which read *back* the PT
solve's own output): bound `dt` using only **a-priori** quantities — the
imposed `vbc` velocity, the boundary mesh size, and the current yield
headroom — so that even the worst case (all of `v_bc` absorbed by the
boundary-adjacent element) cannot drive that element past yield within one
step. Implemented in `compute_dt_PT()` (`geometry.cxx`), gated on
`rh_plastic`.

**Bug found and fixed: global min-h / max-vbc mismatch.** The first version
paired the *globally* smallest boundary element size with the *globally*
largest `|vbc_val|`. In this model those come from different faces: the
finest boundary mesh is on the y-faces (`h ≈ 243 m`), where `vbc_val_y0 =
vbc_val_y1 = 0` (no imposed velocity), while the largest `|vbc_val|` is on
the coarsely-meshed x-faces. Combining them produced an unphysical, ~10x
over-conservative worst case that needlessly throttled `dt` domain-wide —
which in turn starved the *weak zone's* own elastic stress build-up, the
opposite of the intended behavior (the predefined weak zones are supposed to
yield almost immediately when the deformation wave reaches them, precisely
to prevent pervasive yielding elsewhere).

**Fix:** compute `h_min` and `eps_worst = |vbc_face| / h_min_face` **per
face**, pairing each face's own imposed velocity with its own local mesh
size, and take the max `eps_worst` only over the faces an element actually
touches. Result: `dt` jumped from 2.66 yr (buggy global pairing) to 84.08
yr; physical progress reached `t ≈ 606 yr` by step 10 (vs. `~26 yr`), and
the first genuine yielding events resolved quickly instead of cascading.

**Status:** the per-face `dt_bc_yield` is a real improvement over the
original boundary-first-yielding bug for the early elasto-plastic steps, but
it is *not* a complete fix for the whole simulation — later steps still
misbehaved, which ultimately traced to a separate, deeper problem (the
Winkler instability below), not to the boundary-yield mechanism.
(`dt_bc_yield` uses the same deviatoric, friction-ignored `tau_eq` proxy as
solution 3 for the headroom estimate — see the caveat about that proxy under
"solution 3" above; under gravity/confining pressure the *true*
Mohr-Coulomb threshold includes a friction×pressure term the proxy omits, so
the proxy is conservative, which is acceptable for a `dt` cap.)

## Negative result: disabling solution 3 under the staggered scheme

Hypothesis: solution 3's yield-aware `mu_ve_e` softening — designed to tame
the *interleaved* damped-relaxation-vs-undamped-correction limit cycle of
the pre-staggered code — is redundant and possibly harmful under
`PT_staggered`, where the elastic sub-solve is correction-free
(`apply_correction=false`) and there is no interleaved correction to protect
against. Gated it off under staggered mode (`&& !param.control.PT_staggered`
in `update_pt_params()`).

**Result: negative.** The step-11 residual convergence floor persisted at
essentially the same magnitude (~0.0042 vs. 0.004223 before). So solution 3
is *not* the cause of that floor. The gate was left in place (its rationale
is genuinely correct under staggered mode — it's just not the fix for this
floor).

## Residual metric A/B: raw-sum vs delta-from-reference

The user recalled that an *earlier* implementation converged to `1e-6` every
step and suspected the delta-from-reference `calculate_characteristic_force`
change (see "Force-scale / gravity discovery" above) was responsible for the
later degradation. Tested directly by temporarily reverting
`calculate_characteristic_force()` to the original raw-`|tr|`-sum formula.

**Result: nuanced.** With the raw-sum metric, steps 1–4 *did* converge
cleanly (partially confirming the memory), but the same progressive
degradation still appeared from step 5 onward at a similar per-step growth
rate. So the metric choice affects *when* the degradation becomes visible,
but is **not** the root cause — that was the Winkler vertical instability
(below).

**Decision: keep raw-sum; delta removed.** The delta variant was built to
solve a problem it turned out not to cause, and it carried two hazards raw-sum
doesn't:

1. **Broken across remeshing.** Its `tmp_result_ref` snapshot is not
   reallocated/interpolated on remesh, so `force_scale` becomes garbage after
   the first remesh — and the production config is `remeshing_option = 11`.
2. **Denominator collapse near the reference state.** `|tr − tr_ref| → 0`
   whenever the stress state is near the reference (early steps, or any step
   returning near it), making the relative residual blow up spuriously. This
   is what made the delta version *show* degradation earlier in the A/B — a
   less stable denominator, not worse physics.

Raw-sum's only real downside — the locally-balancing lithostatic/gravity
background inflates the scale, so the effective tolerance on the *driving*
problem is looser than the nominal `1e-6` — is a tunable looseness
(`PT_relative_tolerance`), not a correctness bug. If a genuinely
background-free scale is ever wanted, the right construction is a
remesh-invariant, stateless scale from the boundary/imposed-BC tractions
(what actually drives the problem), *not* a one-time global snapshot.

The delta machinery was removed: `calculate_characteristic_force()`
(`fields.cxx`) is the raw-`|tr|`-sum, and the now-unused `tmp_result_ref`
member was deleted from `parameters.hpp`.

## Major finding: Winkler-driven vertical instability (root cause) and its fix

The above threads all circled a degradation that also appears in the
**pure-elastic** case — the simplest possible configuration, with none of
the plastic/staggered/solution-3 machinery active. Reducing to that case
made the true root cause visible.

### Diagnosis: an exponential *vertical* blowup

Reconstructing the force residual per component from output frames (raw
`force` PointData + `bcflag`, zeroing x at x0/x1 and y at y0/y1, never z —
since `vbc_z0 = vbc_z1 = 0` leaves z unconstrained everywhere) on
`core-complex-mmg-pt-elastic-noremesh.cfg` (`remeshing_option = 0`, to run
past the recurring MMG3D crash):

- z-residual is **exactly zero at t=0** (confirming
  `initial_body_force_adjustment()` establishes lithostatic balance), then
  grows roughly geometrically, dominating x/y by 2–6 orders of magnitude,
  and is roughly *uniform* across the whole domain.
- It is **purely vertical**: at every unstable frame `vz_max == vmax` to
  full precision. `vmax` grows `4e-3 → 4.3 → 7.7e3 → … → 5.7e24 m/s`,
  **sign-flipping every step** (mesh coords fly to ±100s of km, alternating).
- Physical time **freezes at 16003.83 yr** — the advection CFL `dt ~ h/vmax
  → 0` as `vmax → ∞`, so real time stalls while the instability races ahead
  in step count.

Sign-oscillating exponential growth is the signature of an *under-damped /
unstable mode* in a pseudo-time integrator — the numerical opposite of the
DR-style settling we were hoping for. The "geometrically growing z-residual"
was the *onset* of this blowup, not slow settling.

### Control experiment: a hard vertical pin removes it entirely

The vertical direction has no Dirichlet constraint (`vbc_z0 = vbc_z1 = 0`);
the only thing resisting vertical motion is gravity + the soft Winkler
foundation. Test: replace the Winkler support with a hard pin — but note
**`input.cxx:1276` makes Winkler take priority over `vbc_z0`**: with Winkler
on, setting `vbc_z0 = 1` is *silently forced back to 0*. You must set
`has_winkler_foundation = no` **and** `vbc_z0 = 1` together (config
`core-complex-mmg-pt-elastic-zpin.cfg`). (A first attempt that set only
`vbc_z0 = 1` was a no-op — bit-identical to the unpinned run — which is how
this priority rule was found.)

**Result:** with the bottom hard-pinned and Winkler off, the pure-elastic
case is completely stable — it ran the **full 500 kyr (126 frames)**,
`vmax` pinned at `1.001e-9 m/s` throughout, `vz_max` **monotonically
decaying** (`9.0e-11 → 4.69e-11`), PT converging to the `1e-6` tolerance
every step (~750–840 iters). So the vertical rigid-body mode *can* settle
DR-style — the Winkler foundation specifically is what destabilizes it.

### Mechanism (code trace)

`apply_stress_bcs()` (`bc.cxx`) computes the bottom Winkler pressure `p =
compensation_pressure − (ρ+Δρ)·g·(zcenter + zlength)`, where `zcenter` is
the **current** facet-center z-position — this is the spring. It's called
from `update_force()` (`fields.cxx:705`), every PT iteration.

The decisive fact: **the mesh is frozen during PT iterations** (explicit
comment in the `dynearthsol.cxx` PT loop; no `update_coordinate` call
inside — the mesh advects once, after the inner loop). So `zcenter`, and
therefore `p`, is **constant across all PT iterations of a step** — the
Winkler spring contributes *zero stiffness* during the solve; it is a
constant load, not a spring. It only behaves as a spring *across physical
steps*, once the outer loop finally moves the mesh.

With the mesh frozen, a uniform vertical translation produces zero internal
elastic strain → zero internal restoring force → it is a rigid-body,
zero-energy mode. With no Dirichlet z-constraint and no Winkler stiffness
during PT, **that vertical rigid-body mode is completely unconstrained in
the frozen-mesh PT problem.** `update_velocity_PT()` (`fields.cxx:930`)
integrates the *raw* force (Winkler load + gravity included) into velocity
for every DOF including z; `apply_vbcs()` re-pins x/y but nothing resets z.
The net vertical force (Winkler load + gravity + internal stress not exactly
balancing once the domain deforms away from the frozen position)
continuously kicks `vel_z`; an unconstrained mode with a persistent net
force cannot be zeroed → PT stalls (the ~3% floor), leaving a net vertical
drift that advects the mesh, shifts `zcenter`, and compounds into the 1e24
blowup. A hard z-pin cures it precisely because Dirichlet-constraining the
bottom z-DOFs *removes the rigid-body mode*.

So the Winkler treatment is **worse than explicit**: during the PT iteration
it provides literally zero stiffness for the very mode it is supposed to
support.

### Fix: semi-implicit (predicted-position) Winkler

Evaluate the Winkler pressure at the **predicted end-of-step position** `z +
vz·dt` — exactly where `update_coordinate` will move the boundary — instead
of the frozen current position. This is backward-Euler in the vertical: it
makes `p` depend on the current PT iterate's velocity, restoring the
physical spring stiffness *during* the frozen-mesh iteration. Change is
localized to the Winkler branch of `apply_stress_bcs()` (`bc.cxx`), gated on
`param.control.has_PT` (the DR path moves the mesh every step, so it keeps
the explicit position):

```cpp
double zc = zcenter;
if (param.control.has_PT) {
    double vz = 0;
    for (int j=0; j<NODES_PER_FACET; ++j)
        vz += (*var.vel)[idx[j]][NDIMS-1];
    zc += (vz / NODES_PER_FACET) * var.dt;
}
p = var.compensation_pressure
  - (rho_effective + param.bc.winkler_delta_rho) * param.control.gravity
    * (zc + param.mesh.zlength);
```

The new term linearizes to `δp = −(ρ+Δρ)·g·vz·dt` — a velocity-proportional
restoring (damping) force on the vertical mode, converting the zero-energy
mode into a damped, convergent one. Stability margin is **`dt`-independent**:
the effective damping `c ∝ dt` is multiplied in `update_velocity_PT` by
`dtau_rho ∝ 1/(G·dt)`, so `dtau_rho·c` is `dt`-independent (`~10⁻³` for this
model — gentle and stable, yet nonzero enough to remove the singular mode).

**Result (Winkler ON, no z-pin — the exact config that previously blew up):**
PT now converges to the `1e-6` tolerance every step, `vmax` stays physical
(~`9.9e-10 m/s`), no blowup. Iteration count shows a one-time startup
transient (step 1: 11255, step 2: 12735 — the vertical mode relaxing from
the t=0 balanced state into its steady deformation-driven drift) then
**settles to ~700 iters/step by step ~10** (`3544 → 1472 → 1062 → 812 → 751
→ … → 693`) — comparable to, even slightly better than, the hard-pin
control. So the anticipated slow-convergence concern is only a transient,
and the pseudo-mass escalation fallback (folding `c` into the boundary-node
pseudo-density in `geometry.cxx`) is **not needed**.

This retains the physical Winkler foundation (no artificial pin) while
curing the root-cause instability. It is the first fix in this investigation
that resolves the actual failure mechanism for the pure-elastic case rather
than shifting a residual magnitude within the same regime.

## Plastic-case validation: the boundary-pluck is resolved

Ran the real elasto-plastic target `core-complex-mmg-pt.cfg` (staggered
scheme, Winkler on, `remeshing_option = 11`) with the semi-implicit Winkler
fix + the per-face `dt_bc_yield` limiter both active. It reached **396 output
frames** — far past where the recurring MMG3D crash used to cut runs off.

**The two things this investigation targeted are fixed:**

1. **No Winkler vertical blowup.** `vmax` stayed physical (~`1e-9 m/s`)
   through the elastic phase; PT converged to the `1e-6` tolerance **every
   step** — both the staggered elastic sub-solve and the return-map
   correction rounds (**3566 "converged", zero "stagnated"** across the whole
   run, including large yielding cascades that corrected 2400–2560 elements
   and still hit `residual = 0.000001`).

2. **The boundary-pluck is resolved.** In the elastic phase (frame 8,
   t=505 yr), **2374/3018 nodes carry velocity** in a smooth, domain-spanning
   ramp (`−8.7e-10` at the left face → ~0 at center → `+8.7e-10` at the
   right) — versus the old broken pattern of **44/3018 nodes, boundary-only**.
   Then, mapping actively-yielding element locations over time (via the
   `plastic strain-rate` field), the **weak zone activates and becomes the
   dominant, localizing deformation zone**: boundary yielding stays low and
   roughly constant (~30 → ~170 elements), while at t≈1280 yr the weak zone
   ignites and takes over — the median yielding location locks onto **x ≈
   52–53 km (weak-zone center at 50 km)** with 1000–1700 elements yielding in
   the 35–65 km band vs. ~170 near the boundaries (~10:1). That is the
   intended core-complex shear-zone behavior, the opposite of the pluck.

### Remaining issue (separate, benign): `dt_weakening` throughput under localization

Once the weak-zone band localizes into sustained vigorous plastic flow
(~2500 elements re-yielding per step), physical time advances very slowly —
`dt` drops from ~84 yr to a plateau of ~1.5–2.6e-3 yr and time creeps at
~1978 yr. **This is not a correctness or stability failure, and not the
`dt_yield` compounding loop** (that hinged on PT *not* converging; here PT
converges cleanly to `1e-6` every step, 0 stagnations, even in the slow
regime). It is `dt_weakening` (`geometry.cxx`) operating **as designed**:
`dt_weakening = allowance/rate` (with `allowance` scaled by
`PT_dpls_fraction = 0.1`), so an actively-weakening band pins `dt` small to
resolve its strain-softening curve. `dt` plateaus (does not → 0, in fact
rises slightly late — consistent with the band beginning to saturate), and
`vmax` fluctuates bounded at ~100–600× the BC velocity — a real, bounded
shear-band slip, not a numerical runaway.

The limitation is **throughput**: during active localization PT collapses to
DR-like small steps, partially defeating its large-step advantage — even
though the solution stays converged, stable, and physically correct. Levers:
- Raise `PT_dpls_fraction` to allow larger plastic increments per step
  (trades weakening-curve resolution for larger `dt`).
- Confirm the "relaxes after it saturates" design intent by running long
  enough for the band's `plstrain` to exceed `pls1 = 0.5` everywhere (no
  further weakening → the cap should release and `dt` recover). The late
  slight `dt` rise is suggestive but unconfirmed.

**Bottom line:** the semi-implicit Winkler fix + per-face `dt_bc_yield`
together get the elasto-plastic model to do the physically right thing —
transmit deformation inward and localize in the weak zone, with PT converging
cleanly throughout. The boundary-pluck this whole document is about is
resolved. What remains is a distinct `dt_weakening` *throughput* question
under active localization, not a stability or convergence defect.

> **CORRECTION (see "Root cause" below).** The "localize in the weak zone"
> claim above was WRONG — it came from reading the `plastic strain-rate`
> output field with a loose threshold, which registered the weak zone merely
> *sitting at* its IC yield surface (pls=0.5) as if it were accumulating.
> Direct ParaView inspection of the *accumulated* `plastic strain` field
> (user-verified) showed the opposite: the weak zone stays inert and the real
> plastic strain accumulates at the **top corners of the velocity boundaries**.
> The boundary-pluck was NOT resolved by the Winkler + `dt_bc_yield` fixes.
> The real root cause was found afterward and is documented below.

## Root cause: the staggered scheme itself causes the pluck

After the correction above, a clean isolation finally pinned the mechanism.
The `plastic strain-rate` field is misleading for this question (it is nonzero
wherever material sits on the yield surface, including the pre-weakened
`pls=0.5` weak zone, whether or not strain is accumulating). Judge yielding by
the *accumulated* `plastic strain` field instead.

**Why the corners yield first (a real effect, but a red herring for the
pluck).** Computing the Mohr-Coulomb `fs = p[0] − p[NDIMS-1]·anphi + amc`
per element at the perfect-ramp frame showed the elements closest to yield are
the velocity-BC **top corners** (near-surface, near-zero confining pressure, so
strength ≈ cohesion only), not the weak zone (at depth, where the
`friction×pressure` term makes it *strong* despite its lower cohesion). This
motivated a series of config experiments — moving the weak zone shallow
(`weakzone_depth 0→0.5`), lowering weak-zone friction (`friction_angle1 30→10`),
and raising `cohesion0` uniformly (`4.4e7→4.4e8`, "strong crust, weak
detachment"). Each shifted *where* strength was lowest, but **none stopped the
pluck**: the velocity ramp still collapsed to boundary-only once anything
yielded.

**The `force_scale` red herring.** One run showed the velocity going
boundary-only (44 nodes at ±vbc, interior at ~1e-15 = the `1e-6` tolerance ×
boundary) with *nothing yielding much*, which looked like a convergence-
criterion artifact: the raw-sum `force_scale ≈ 7.7e14` is lithostatic-
dominated, so `1e-6·force_scale` is a loose absolute tolerance. **This was
tested and refuted:** the pure-elastic case (`rheology_type = elastic`), same
geometry / `dt` / `force_scale` / tolerance, produces a **perfect ramp**
(2915 interior nodes) at `1e-6`. So the raw-sum metric is fine (the earlier
"keep raw-sum" decision stands) and the pluck is not a `force_scale` problem.

**The decisive A/B.** Same config, only the correction scheme differs:

| configuration | velocity field | weak zone | `dt` / time |
|---|---|---|---|
| `rheology = elastic` | ramp (2915 nodes) | — | stable |
| `elasto-plastic`, **`PT_staggered = no`** | **ramp (2930 nodes)** | **localizes, pls 0.5→0.54** | **stable, advances** |
| `elasto-plastic`, `PT_staggered = yes` | **pluck (44 boundary-only)** | inert / never loads | collapses, freezes |

With `PT_staggered = no` — the older **damped per-iteration** return map (the
`w0`/`w1` blend inside `update_stress_PT`) — the ramp survives yielding and the
weak zone localizes cleanly: at step 50, 209 elements with `pls>0.51`, **all in
the weak-zone band (x=46–52 km), zero elsewhere**, max pls 0.536, `dt` steady
(~76 yr), time advancing, `vmax` physical. That is the intended core-complex
behavior.

**So the boundary-pluck is an artifact of the STAGGERED scheme's full,
undamped return-map projection + `stress_old` re-anchoring** — the abrupt
projection breaks the discretized force balance the elastic solve satisfies,
and re-equilibrating from it collapses the interior velocity to zero. The
irony: the staggered scheme (above) was the *experimental attempt to fix* the
pluck, but it is what *causes* it. The pre-existing non-staggered damped
correction handles yielding gracefully. Everything chased earlier —
`force_scale`, `dt`/`dt_weakening`, `dt_bc_yield`, corner vs weak-zone
placement, friction, cohesion — was downstream symptom or red herring; the
pluck rode on the staggered return map the whole time.

**Fix:** run with `PT_staggered = no`. The staggered scheme should be
considered broken for this problem (or fixed so its projection preserves force
balance) rather than used.

**Open follow-up (not yet done):** isolate what is *essential*. `PT_staggered
= no` fixes the *pluck*; the shallow weak zone / low `friction_angle1` / high
`cohesion0` changes only steer *where* localization nucleates (weak zone vs
BC corner). Worth testing whether the plain original config (deep weak zone,
`cohesion0 = 4.4e7`, friction 30) also localizes correctly with just
`PT_staggered = no`, to separate the pluck fix from the localization-steering
knobs.

## Resolution: isolation confirmed, and the staggered scheme removed from the code

The isolation follow-up above was run: the plain **original** config (deep weak
zone, `cohesion0 = 4.4e7`, friction 30, `dt_fraction = 1.0`) with *only*
`PT_staggered = no` keeps the ramp intact (2900+ interior nodes moving) **even
while the BC corner actively yields** (33 → 79 elements over the first steps).
So the two effects are cleanly separable and independent:

- **`PT_staggered = no` fixes the pluck** — unconditionally, regardless of where
  yielding happens. The ramp survives yielding.
- **Strength/geometry (weak-zone depth, `friction_angle1`, `cohesion0`) steers
  *where* deformation localizes** — with the original deep weak zone it stays
  inert (strong at depth) and the corner yields; the shallow, low-friction weak
  zone + strong crust move localization onto the weak zone (the intended
  core-complex shear band). This is Mohr-Coulomb physics, not a solver artifact.

**The staggered scheme has since been deleted from the code.** Because it only
ever *caused* the pluck and the non-staggered damped correction is strictly
better here, keeping it as a toggle was a footgun. Removed: `control.PT_staggered`
and `control.PT_max_outer_iter` (input parsing + `Param` struct), the staggered
outer loop in `dynearthsol.cxx`, `apply_return_map_PT()` in `rheology.cxx`, and
the now-redundant `apply_correction` argument of `update_stress_PT()` (the
damped return map is now always applied every PT iteration for plastic, non-RSF
rheologies). The `!PT_staggered` gate on the yield-aware `mu_ve_e` softening
("solution 3") in `update_pt_params()` was likewise dropped — that softening is
active for all plastic rheologies now, its original design. So the mode this
whole document recommends (`PT_staggered = no`) is no longer a setting; it is
simply how the PT solver works.

**Net outcome of the investigation.** The boundary-pluck is resolved. The PT
solver transmits boundary deformation inward and localizes it in the predefined
weak zone (given weak-zone parameters that make it the genuine weakest point at
its depth), converging cleanly to tolerance. Two separate, non-blocking items
were characterized along the way and remain open: the `dt_weakening` throughput
cost during vigorous localization (PT falls back to DR-scale small steps), and
the broader PT-vs-DR cost picture (PT is currently ~2× slower than DR in
wall-clock on this CPU build, dominated by per-step iteration count that grows
with localization — the main lever there is running the GPU/`ACC` build, which
is where accelerated PT's advantage actually lives).
