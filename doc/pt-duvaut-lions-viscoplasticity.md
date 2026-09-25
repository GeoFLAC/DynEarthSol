# Duvaut-Lions viscoplastic regularization for `rh_ep`

## Context

The accelerated Pseudo-Transient (PT) solver is currently ~2x slower than the
classic Dynamic Relaxation (DR) solver on the 3D weak-zone benchmark
(`examples/gaussian-weakzone-3d-PT*.cfg`) used throughout this session's
DR-vs-PT investigation. The bottleneck is `dt_weakening`
(`geometry.cxx`'s `compute_dt_PT()`), which collapses PT's time step once
plastic localization becomes vigorous. `dt_weakening` exists because
cohesion/friction (a piecewise-linear function of accumulated plastic strain
`plstrain`) are frozen for a whole step and only updated once, post-
convergence — a step that lets `plstrain` jump too far makes those frozen
parameters stale. Under today's rate-independent (instantaneous Mohr-Coulomb
return-map) plasticity, a large step's elastic trial stress can overshoot the
yield surface arbitrarily far, producing an unbounded per-step plastic-strain
increment — forcing `dt_weakening` to shrink `dt` hard, precisely defeating
PT's large-step advantage.

Duvaut-Lions viscoplastic regularization blends the elastic trial stress with
the rate-independent ("inviscid") plastically-corrected stress, weighted by
`dt/(tau+dt)` for a new material relaxation time `tau`. This makes the
plastic-strain increment intrinsically proportional to `dt` (bounded by
`tau`), which should relax `dt_weakening` substantially for PT's large steps
while barely perturbing DR (whose steps are already `dt >> tau`, close to the
rate-independent limit). It also algorithmically reuses a blend structure
(`s_pre`/`w0`/`w1`) already present in this codebase's PT return-map damping
(`update_stress_PT()`), just with a physical relaxation time instead of a
numerical PT pseudo-damping ratio.

**Goal of this change:** implement the blend where physical plastic strain is
actually accumulated, gated behind an off-by-default flag, with a structural
guarantee that the default (and `tau=0`) path is bit-identical to today's
code — then validate that it relaxes `dt_weakening` and narrows the PT/DR
wall-time gap without corrupting the physical answer.

## Scope

**In scope:** `MatProps::rh_ep` only — the rheology the benchmark
(`rheology_type = elasto-plastic`) actually uses.

**Explicitly out of scope** (documented, not attempted):
- `rh_evp` (`rheology.cxx:872-920`) — already has its own rate-dependence via
  a Maxwell-viscoelastic candidate competing on J2; layering DL on top raises
  an unresolved question (does DL's `tau` apply on top of or instead of the
  Maxwell time?) best left to a follow-up.
- `rh_ep_rsf` / `rh_evp_rsf` — already rate-dependent via state-variable
  evolution.
- `update_stress_PT()` (`rheology.cxx:1032-1148`) — its `w0`/`w1` blend is an
  unrelated *numerical* PT pseudo-damping device (built from
  `PT_Gdtau_e`/`Gdt`), and the `depls` it computes is already discarded
  (zero strain increment, pure re-projection). Physical plastic strain is
  accumulated only in `update_stress()`'s corrector call — that's the only
  function this change touches.
- Full OpenACC data-lifetime wiring for `MatProps` — the existing
  `#pragma acc exit data delete(...)` lists in `matprops.cxx` (lines
  308-317) are already incomplete (e.g. `cohesion1` is missing from the list
  that includes `cohesion0`) and have no matching `enter data create` calls
  anywhere in the file. This is a pre-existing gap, and this session's build
  target is exclusively the non-ACC CPU path (`openacc=0`). Add the new
  vector to the existing (already-incomplete) delete list for consistency,
  but do not attempt to fix ACC support as part of this change.
- 2D plane-strain tension-cutoff proportionality (see validation below) —
  documented as an accepted, bounded approximation.

## Naming (matches existing conventions)

- Gate: `control.has_duvaut_lions` (bool, default `false` — matches
  `control.has_hydraulic_diffusion`/`has_PT` naming, `parameters.hpp:251`).
- Per-material array: `mat.relaxation_time` (matches `mat.cohesion0` style).
- `MatProps` private storage: `double_vec relaxation_time;`
- `MatProps` public accessor: `double tau_dl(int e) const;` — a short,
  distinct name from the storage member (matching the existing
  `direct_a`/`d_a()`, `evolution_b`/`e_b()` abbreviation pattern), avoiding
  any collision with `rh_evp`'s pre-existing "elasto-**visco**-plastic"
  naming (a user must not read `mat.relaxation_time` and assume it affects
  `rh_evp`).

## Concrete changes

### `parameters.hpp`
- `Control` struct: `bool has_duvaut_lions;` next to `has_hydraulic_diffusion`
  (line 251).
- `Mat` struct: `double_vec relaxation_time;` next to `cohesion0, cohesion1;`
  (line 472).

### `input.cxx`
- Register the gate alongside the other `control.has_*` options:
  ```cpp
  ("control.has_duvaut_lions", po::value<bool>(&p.control.has_duvaut_lions)->default_value(false),
   "Enable Duvaut-Lions viscoplastic regularization of the rh_ep return map")
  ```
- Register the array right after `mat.cohesion1` (~line 863-864), default
  `"[0]"` so an enabled-but-unconfigured gate degrades safely to the
  rate-independent limit rather than picking an arbitrary nonzero `tau`:
  ```cpp
  ("mat.relaxation_time", po::value<std::string>()->default_value("[0]"),
   "Duvaut-Lions viscoplastic relaxation time of the materials '[d0, d1, ...]' "
   "(seconds; only used when control.has_duvaut_lions is true; 0 recovers "
   "rate-independent plasticity)")
  ```
- Parse call right after `get_numbers(vm, "mat.cohesion1", ...)` (~line 1476):
  ```cpp
  get_numbers(vm, "mat.relaxation_time", p.mat.relaxation_time, p.mat.nmat, -1);
  ```

### `matprops.hpp`
- Public: `#pragma acc routine seq` + `double tau_dl(int e) const;` next to
  `pls_weakening_allowance` (~line 67-68).
- Private: `double_vec relaxation_time;` next to `cohesion0, cohesion1;`
  (line 133).

### `matprops.cxx`
- Constructor: `relaxation_time = p.mat.relaxation_time;` right after
  `cohesion1 = p.mat.cohesion1;` (~line 211).
- New accessor (on-demand marker-weighted mean — same idiom as the `tmax`
  computation inside `plastic_props()`, ~lines 597-606; *not* the cached
  `bulkm`/`shearm` idiom, since this is read once per element per
  `update_stress()` call, not a hot per-PT-iteration path):
  ```cpp
  double MatProps::tau_dl(int e) const
  {
      double tau = 0;
      int n = 0;
      for (int m = 0; m < nmat; m++) {
          int k = elemmarkers[e][m];
          if (k == 0) continue;
          n += k;
          tau += relaxation_time[m] * k;
      }
      return (n > 0) ? tau / n : 0.0;
  }
  ```
  (`get_numbers(..., p.mat.nmat, -1)` already expands `relaxation_time` to
  length `nmat` at parse time, so direct `relaxation_time[m]` indexing is
  safe — same as the existing `tmax` loop.)
- Add `relaxation_time` to the existing `#pragma acc exit data delete(...)`
  list at line 310 (alongside `cohesion0`), for consistency with (already
  incomplete) existing practice — not a fix to that pre-existing gap.

### `rheology.cxx` — two changes

**(a) Prerequisite, isolated, behavior-identical refactor:** extract
`elasto_plastic2d`'s inline elastic-trial block (lines 518-532) into a
reusable helper so the DL blend can compute a plane-strain trial without
duplicating the formula:
```cpp
#pragma acc routine seq
static void elastic_trial2d(double bulkm, double shearm, const double* de,
                             double* s, double& syy,
                             bool has_hydraulic_diffusion, double dpp)
{
    double a1 = bulkm + 4./3*shearm;
    double a2 = bulkm - 2./3*shearm;
    double sxx = s[0] + de[1]*a2 + de[0]*a1;
    double szz = s[1] + de[0]*a2 + de[1]*a1;
    double sxz = s[2] + de[2]*2*shearm;
    syy += (de[0] + de[1]) * a2;
    if (has_hydraulic_diffusion) { sxx += dpp; syy += dpp; szz += dpp; }
    s[0] = sxx; s[1] = szz; s[2] = sxz;
}
```
Replace lines 518-532 with a call to it. Ship and verify this alone first
(Stage 0 below) before adding any DL logic — isolates "refactor broke
something" from "DL logic broke something" if a bisect is ever needed.

**(b) Modify the `rh_ep` branch** (lines 849-871):
```cpp
case MatProps::rh_ep:
    {
        double depls = 0;
        double bulkm = var.mat->bulkm(e);
        double shearm = var.mat->shearm(e);
        double amc, anphi, anpsi, hardn, ten_max;
        var.mat->plastic_props(e, plstrain[e], amc, anphi, anpsi, hardn, ten_max);
        int failure_mode;

        const bool dl_on = param.control.has_duvaut_lions;
        const double tau = dl_on ? var.mat->tau_dl(e) : 0.0;

        if (!dl_on || tau <= 0.0 || var.dt <= 0.0) {
            // Exactly today's path -- zero extra cost, zero behavior change.
            if (var.mat->is_plane_strain) {
                elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                 de, depls, s, syy, failure_mode,
                                 has_hydraulic_diffusion, dpp);
            } else {
                elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                               de, depls, s, failure_mode,
                               has_hydraulic_diffusion, dpp);
            }
        } else {
            // Duvaut-Lions blend: weight -> 1 as dt >> tau (rate-independent
            // limit), -> 0 as dt << tau (strong regularization toward the
            // elastic trial). depls scales by the same weight as the stress
            // blend -- see "Blend-consistency validation" below.
            const double weight = var.dt / (tau + var.dt);

            if (var.mat->is_plane_strain) {
                double s_trial[NSTR]; for (int i=0;i<NSTR;++i) s_trial[i]=s[i];
                double syy_trial = syy;
                elastic_trial2d(bulkm, shearm, de, s_trial, syy_trial,
                                 has_hydraulic_diffusion, dpp);

                elasto_plastic2d(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                 de, depls, s, syy, failure_mode,   // s,syy -> inviscid
                                 has_hydraulic_diffusion, dpp);

                for (int i=0;i<NSTR;++i) s[i] = (1-weight)*s_trial[i] + weight*s[i];
                syy = (1-weight)*syy_trial + weight*syy;
            } else {
                double s_trial[NSTR]; for (int i=0;i<NSTR;++i) s_trial[i]=s[i];
                if (has_hydraulic_diffusion) elastic_effective(bulkm, shearm, de, s_trial, dpp);
                else                          elastic(bulkm, shearm, de, s_trial);

                elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                               de, depls, s, failure_mode,          // s -> inviscid
                               has_hydraulic_diffusion, dpp);

                for (int i=0;i<NSTR;++i) s[i] = (1-weight)*s_trial[i] + weight*s[i];
            }
            depls *= weight;
        }
        plstrain[e] += depls;
        delta_plstrain[e] = depls;
    }
    break;
```
`dpp` is already computed once per element earlier in the loop (~lines
731-738) and is safe to reuse for both the trial and inviscid calls (it
depends only on nodal pore pressure, not on `s`/`de`). `elastic()` /
`elastic_effective()` are templated on the accessor type, so a raw
`double[NSTR]` satisfies them without signature changes.

**`update_stress_PT()` and the `rh_evp` branch: no changes**, per Scope above.

## Blend-consistency validation (already checked against the actual return-map code)

The claim `depls_new = weight * depls_inviscid` being consistent with the
stress blend requires the return map's stress correction to be exactly
proportional to `depls` (a "radial" projection). Verified directly:

- **3D `elasto_plastic()`** (both shear failure, lines 388-423, and tensile
  failure, lines 425-457): a single scalar (`alam`) linearly scales the
  principal-stress correction *and* `depls` is computed as
  `|alam| * sqrt(const)` — exactly proportional, in both branches. The
  no-failure early-return trivially satisfies the blend too
  (`s_inviscid==s_trial`, `depls==0`).
- **Plane-strain `elasto_plastic2d()`, pure-shear branch** (lines 651-658):
  same structure, exact.
- **Plane-strain `elasto_plastic2d()`, tension-cutoff branches** (lines
  614-619, 621-635, 664-688): these are **hard, absolute clamps** (e.g.
  `s[0]=s[1]=syy=ten_max`), not proportional to any `alam`-like scalar, and
  in the pure-clamp cases `depls` is left at 0 — a decoupling that already
  exists in today's rate-independent code (not introduced by this change).
  **Verdict:** accept as a known, bounded approximation — it doesn't affect
  the `tau=0`/disabled bit-identity guarantee (every branch's `weight=1`
  case reduces to exactly today's values), the benchmark this feature
  targets is a 3D problem (exact case) or a compressive shear-localization
  regime where multiaxial-tension elements are a small non-localizing
  subset. Do not attempt to fix `elasto_plastic2d`'s legacy tension-cutoff
  proportionality as part of this change.

## Staged validation

1. **Build sanity.** Land the `elastic_trial2d` refactor (1a) alone first;
   confirm 2D and 3D builds produce unchanged output on an existing plastic
   example (pure no-op refactor, isolated commit).
2. **Regression safety, gate off (default).** Re-run every config already
   used this session (`examples/gaussian-weakzone-3d.cfg`,
   `-PT.cfg`, `-PT-1e4/1e5/1e6.cfg`, `-PT-loose.cfg`) with
   `has_duvaut_lions` unset. Diff `run.log`'s `debug.dt=1` line
   (`geometry.cxx:1936-1940`, the `dt_maxwell dt_advection dt_diffusion
   dt_hydro_diffusion dt_weakening dt_bc_yield` breakdown) and the
   `.vtkhdf` output (reuse this session's h5py node-level velocity/stress
   read pattern) against pre-change baselines — expect bit-identical
   results.
3. **Gate-mechanism cross-check.** Re-run with `has_duvaut_lions=true`,
   `mat.relaxation_time=[0]` (exercises the arithmetic `tau<=0` short-circuit
   instead of the structural off-gate) — must also match Stage 2's baseline
   bit-for-bit, confirming both "disabled" paths agree.
4. **Intended-effect validation.** Enable with a swept `tau` (pick a
   starting order of magnitude from the unmodified run's `dt_weakening`
   column once it starts collapsing) on the PT configs and check: (a)
   `dt_weakening` stops collapsing as hard during active localization, (b)
   PT inner-iteration counts are essentially unaffected (this change only
   touches the once-per-step corrector, not `update_stress_PT()`'s
   convergence loop), (c) wall-clock time narrows the PT/DR gap, (d) the
   `.vtkhdf`-based physical-sanity checks (plastic-strain/shear-band
   pattern, velocity field, `delta_plstrain` smoothness per step) stay
   qualitatively consistent with the DR reference — differences from DR are
   *expected* (DL is a genuine physical change at large `dt`, not a purely
   numerical one) but should shrink as `tau→0`, and shear-band location
   should not qualitatively change; if it does, `tau` is too large for the
   chosen `dt` regime.
5. **(Future work, not this change) `tau` calibration** — systematically
   tie `tau` to a physically- or numerically-motivated scale (e.g. a Maxwell
   time, or the DR step size) once the mechanism is validated qualitatively.
