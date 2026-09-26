# Smooth weakening law + implicit self-consistent return map for `rh_ep`

## Context

`dt_weakening` (`geometry.cxx`'s `compute_dt_PT()`) limits PT's time step
because cohesion/friction/dilation (a piecewise-linear function of
accumulated plastic strain `plstrain`) are evaluated once, at start-of-step
`plstrain`, and held frozen through the whole step's Mohr-Coulomb return-map
solve. A step whose plastic-strain increment sweeps past a large fraction of
the ramp makes those frozen parameters stale, and the ramp's kinks at
`pls0`/`pls1` have no well-defined derivative, ruling out a local Newton
correction. Duvaut-Lions regularization (`doc/pt-duvaut-lions-viscoplasticity.md`)
attacked the symptom (how much `plstrain` a step is allowed to accumulate);
this change attacks the root cause directly: evaluate cohesion/friction
**self-consistently at the converged end-of-step `plstrain`**, via a smooth
weakening law (so the parameter is always differentiable) and a fixed-point
iteration on the return map itself (so "converged end-of-step" is actually
solved for, not just approximated once).

**Goal:** remove the staleness that motivates `dt_weakening`, gated behind an
off-by-default flag with a structural guarantee that the default path is
bit-identical to today's code, then validate on the same PT weak-zone
benchmark used for the Duvaut-Lions work.

## Design

### Smooth weakening law (`control.has_smooth_weakening`)

`MatProps::plastic_weakening()` (`matprops.cxx`) gained an exponential-decay
branch alongside the existing piecewise-linear ramp:

```
X(pls) = X1 + (X0 - X1)*exp(-pls/scale)      dX/dpls = -(X0-X1)/scale*exp(-pls/scale)
```

for cohesion, friction angle, and dilation angle, using the existing
`cohesion0/1`, `friction_angle0/1`, `dilation_angle0/1` virgin/residual
values unchanged, plus a new per-material parameter `mat.pls_scale` (the
e-folding strain) in place of `pls1` as the ramp-width control. `pls0` is
unused by this law (no onset plateau) but left in place for the
piecewise-linear path. `pls_scale` is a deliberately new parameter name, not
a reinterpretation of `pls1` — silently repurposing `pls1`'s meaning would
change what existing `.cfg` files mean with no syntax change.

### Implicit self-consistent return map (`elasto_plastic{,2d}_implicit`, `rheology.cxx`)

Rather than hand-deriving the analytic Newton Jacobian for the yield
consistency condition (worked out partially already in
`doc/deferred-improvements.md`, for the cohesion+friction terms only), the
implementation uses a **fixed-point (Picard) iteration**: re-solve the
existing closed-form return map (`elasto_plastic`/`elasto_plastic2d`,
unchanged) with `hardn=0` at successively refined estimates of end-of-step
`plstrain`, until self-consistent:

```
kappa = pls_old
repeat up to 20 times:
    amc, anphi, ... = plastic_props(e, kappa)      # smooth law, evaluated at the current guess
    (re-run the return map from the saved pre-step stress, hardn=0)
    kappa_new = pls_old + depls
    if |kappa_new - kappa| < tol: converged
    kappa = kappa_new
```

`hardn=0` is intentional: with `amc`/`anphi` held fixed for one inner solve,
the closed-form return map is then the *exact* projection onto that (locally
non-hardening) surface — the outer iteration, not the closed form's own
linear-hardening correction term, is what supplies the end-of-step
softening. This sidesteps the risk of an algebra error in a hand-derived
Jacobian for safety-critical constitutive code, at the cost of linear
(rather than quadratic) convergence — acceptable given the warm start from
`pls_old` and the small per-step `plstrain` increments typical of this
solver.

Per `doc/deferred-improvements.md`'s own finding ("dilation-angle softening
affects the flow direction but not the yield consistency condition"),
`anpsi` (flow direction) stays frozen at `pls_old` throughout the iteration —
not part of the self-consistent solve. This also keeps the iteration a clean
1D root-find in `plstrain`, not a coupled system.

Wired into `update_stress()`'s `rh_ep` case for both the plain and
Duvaut-Lions-blended paths, gated on `control.has_smooth_weakening`; the
`!implicit_on` branch is exactly today's code, byte-for-byte.

### `dt_weakening` bypass

`compute_dt_PT()`'s `check_weakening` gate now also requires
`!control.has_smooth_weakening` — with self-consistent end-of-step
parameters there is no staleness left for this limiter to protect against.
`pls_weakening_allowance()` itself is untouched (still used by the
piecewise-linear/default path).

### Scope

`rh_ep` only, 3D (`elasto_plastic`) and 2D plane-strain (`elasto_plastic2d`)
shear and tensile branches. Out of scope: `rh_evp`/RSF variants,
`update_stress_PT()`'s zero-strain re-projection, 2D's hard tension-cutoff
clamps, `dt_bc_yield` and `update_pt_params()`'s yield-ratio softening
(both still read the ordinary, now-smooth, start-of-step `plastic_props()`
unchanged — see Findings below for why extending the fix to `dt_bc_yield`
would not help).

## Validation

- Build: clean, no new warnings.
- Gate-off: bit-identical by construction (`if (implicit_on) ... else <unchanged code>`,
  default `false`) — confirmed by direct diff review of the `else` branches.
- Gate-on: ran `tests/functional/2d-ep-irregular.cfg` (both `is_plane_strain=yes`,
  exercising `elasto_plastic2d_implicit`, and `=no`, exercising
  `elasto_plastic_implicit`) full-length (40000 steps). No NaNs, no crashes,
  both wrapper functions exercised. Plastic strain/stress/velocity fields
  differ from the piecewise-linear baseline by a small, physically plausible
  amount (e.g. max plstrain 0.53826 → 0.53961; max stress differs by ~1%),
  confirming the new code path is genuinely exercised and stable, not just
  falling through to a no-op.

## Findings: `dt_weakening` was never the bottleneck on the PT weak-zone benchmark

Ran the full `examples/gaussian-weakzone-3d-PT` benchmark (same one used for
the Duvaut-Lions `tau` sweep) with `has_smooth_weakening` off vs on
(`pls_scale=0.15`, comparable ramp width to the baseline's `pls0=0,
pls1=0.5`):

| | steps | sim. time reached | wall-clock | PT stagnation events |
|---|---|---|---|---|
| baseline (piecewise-linear, `dt_weakening` active) | 117 | 8025 yr | 30:10 | 12 |
| `has_smooth_weakening=yes` (`dt_weakening` bypassed) | 118 | 8033 yr | 32:41 | 12 |

No improvement — if anything, ~8% slower from the Picard iteration's
per-element overhead, with an identical stagnation count. Enabling
`debug.dt=yes` on the **baseline** (`has_smooth_weakening` off — this is
today's unmodified code path, not a Stage 1+2 artifact) and reading the raw
`dt_maxwell`/`dt_advection`/`dt_diffusion`/`dt_hydro_diffusion`/
`dt_weakening`/`dt_bc_yield` breakdown per step showed `dt_bc_yield` is the
binding term at *every single step* of the run, matching the actual `dt`
used almost exactly (including its decline from ~46 yr to ~38 yr and
recovery to ~40 yr around steps 97-106). `dt_weakening` sat at 10,000-21,000
years throughout — roughly 300-500x looser than the actual `dt` — never
remotely close to binding, even before this change existed.

So on this benchmark, `dt_weakening`'s staleness was never actually
constraining anything, in the original code or with Stage 1+2 enabled;
`dt_bc_yield` (`geometry.cxx:1801-1936`, a worst-case elastic-overshoot bound
protecting PT's inner-loop stability at velocity-Dirichlet boundaries, see
its own in-code comment for the feedback-loop history) was and remains the
actual ceiling. Stage 1+2 is a correct, validated fix for a real staleness
problem in the return map and is being kept for that reason (better-posed
than the piecewise-linear/frozen-parameter scheme on its own merits), but it
does not — and was never going to — move this benchmark's wall-clock time.

**Why the same treatment does not transfer to `dt_bc_yield`:** `dt_bc_yield`
also reads `amc` at start-of-step `plstrain`, but the staleness there works
in the *opposite* direction from `dt_weakening`'s problem. Since cohesion
only decreases with weakening, the start-of-step `amc` is always \>= the true
(already-softened) end-of-step value, so `headroom = amc - tau_eq` is
currently *overestimated* — `dt_bc_yield` is already more permissive than a
self-consistent version would be. Making its `amc` self-consistent would
shrink `headroom` and make `dt_bc_yield` *more* restrictive, not less. A
genuine improvement there would have to relax the worst-case assumption
itself (100% of a step's imposed boundary velocity absorbed by one
boundary-adjacent element, zero elastic propagation into the interior) or
address the underlying PT inner-loop stability concern it exists to avoid —
neither is a weakening-law problem.

## Follow-up: the implicit update now reaches the PT solve (commit 418eda8)

The implicit return map was applied only by the post-PT corrector
(`update_stress()`); `update_stress_PT()` built PT's target with the ordinary
return map at start-of-step parameters, so PT converged against a stress
different from the one stored. The target now uses
`elasto_plastic{,2d}_implicit()` when `has_smooth_weakening` is on. Even so,
it does not fix the reduced localization of large-`dt` PT steps -- the
plastic-strain increment per step is too small for the parameters to change
much within a step. See "Step size vs. strain localization" at the end of
`doc/pt-boundary-pluck.md`.
