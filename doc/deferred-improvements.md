# Deferred Improvements

This document records improvements that were identified but deferred.
Each entry states the issue, the correct fix, and why it was deferred.

---

## PT solver

### Phase 3 — Gershgorin λ_max refinement

**Issue.** The analytical stability ceiling $\lambda_\max = \mathrm{CFL}^2/(r+2)$ is
derived for a uniform medium.  In models with large viscosity contrasts across
adjacent elements, the true λ_max can exceed this bound, making the PT step size
too aggressive and slowing convergence.

**Fix.** Compute a per-element Gershgorin upper bound on λ_max from the assembled
stiffness matrix entries, and use the global max as the stability limit.  New
function `compute_gershgorin_lambda_max()` in `geometry.cxx`, called from
`rayleigh_update_Re()`.  Also requires updating `dynearthsol.cxx`.

**Deferral reason.** The [0.5, 2] × Re clamping in the Rayleigh-quotient update
(Phase 2) limits the damage from a wrong λ_max.  The Gershgorin bound is only
necessary if Phase 1+2 remain insufficient for models with >3–4 orders-of-magnitude
viscosity contrast.  Left as Phase 3 in the implementation notes (pt-duretz2026.md).

---

### Rayleigh retune: step-2 regression

**Issue.** When the first PT retune window (iterations 1–100 of a new time step)
captures a fast-mode transient (residual drops ~200×), the Rayleigh quotient yields
$\lambda_\min \approx 10$ (fast-mode eigenvalue, not the slow mode).  This drives Re
to the upper clamp (2 × Re_ref ≈ 30), slowing slow-mode convergence for the remainder
of that PT loop.  Observed as a 39% iteration-count increase in step 2 of the 2D EP
benchmark.

**Fix options.**
1. Delay the first retune to iteration 500 (or a configurable `PT_first_retune`
   parameter), so the fast transient has passed before the first Rayleigh estimate.
2. Multi-window exponential averaging of successive λ_min estimates to dampen
   transient spikes.

**Deferral reason.** Step 2 regression is outweighed by savings in the init loop and
step 1 (total: −32% iterations); net effect is positive.  The clamping to
[0.5, 2] × Re_ref bounds the regression.  Revisit if a problem class consistently
triggers this pattern.

---

## Return mapping / plasticity

### Correct hardening modulus for Mohr-Coulomb with softening friction angle

**Issue.** `hardn` in `plastic_weakening()` (`matprops.cxx:404`) stores only the
cohesion slope $dc/d\kappa$.  The return-map consistency condition
(`elasto_plastic()`, `rheology.cxx:392`) uses `2√(anphi)·hardn` in the `alam`
denominator, which accounts for cohesion softening but ignores friction-angle
softening.

The full consistency condition for Mohr-Coulomb yield
$F = \sigma_1 - \alpha_\phi\,\sigma_3 + 2c\sqrt{\alpha_\phi}$
(with $\alpha_\phi = (1+\sin\phi)/(1-\sin\phi)$) gives

$$\frac{dF}{d\kappa} = 2\sqrt{\alpha_\phi}\,\frac{dc}{d\kappa}
  + \left(\frac{c}{\sqrt{\alpha_\phi}} - \sigma_3\right)
    \frac{2\cos\phi}{(1-\sin\phi)^2}\,\frac{d\phi}{d\kappa}$$

The second term is missing.  In models where $d\phi/d\kappa$ is large relative to
$dc/d\kappa$, the increment size $\lambda$ (alam) is underestimated, slightly
over-returning.

**Fix.**
1. Add `dphi_dpls` output to `plastic_weakening()` (`matprops.cxx`), computed as
   `(friction_angle1[m] - friction_angle0[m]) / (pls1[m] - pls0[m])`.
2. Pass the current minimum principal stress into `elasto_plastic()` (already
   available as `p[NDIMS-1]`).
3. Replace `2*std::sqrt(anphi)*hardn` in the `alam` denominator with the full
   two-term expression above, using `dphi_dpls * (2*cos(phi)/(1-sin(phi))^2)`.

Dilation-angle softening ($d\psi/d\kappa$) affects the flow direction but not the
consistency condition $F$, so it does not change `alam`.

**Deferral reason.** The effect is second-order for typical geodynamic softening
laws where cohesion dominates; friction-angle softening is usually slower and smaller
in absolute stress change per unit plastic strain.  The existing approximation is
consistent with what most published MC implementations do.  Revisit if a benchmark
comparison shows systematic over-softening in friction-dominated regimes.

---
