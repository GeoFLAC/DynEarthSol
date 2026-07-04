# The Accelerated Pseudo-Transient (PT) Solver in DynEarthSol

This document describes the accelerated pseudo-transient solver available on
the `feature/accelerated-pt` branch: the theory it implements, how the theory
was adapted to DynEarthSol's explicit FEM architecture, and how to use and
tune it in practice.

The method follows

> Räss, L., Utkin, I., Duretz, T., Omlin, S., and Podladchikov, Y. Y. (2022),
> *Assessing the robustness and scalability of the accelerated pseudo-transient
> method*, Geosci. Model Dev. 15, 5757–5786 (Sec. 2.4),

generalized here to compressible elasto-plasticity on unstructured, graded
tetrahedral/triangular meshes with local pseudo-time stepping.

---

## 1. The problem being solved

At every physical time step of size `Δt`, DynEarthSol needs the quasi-static
momentum balance

```
∇·σ + ρg = 0
```

with the constitutive update

```
σ = τ_old + C : ε̇(v) Δt      (elastic increment)
σ ← return_map(σ)             (Mohr–Coulomb yield + tension cutoff)
```

where `τ_old` is the stress at the beginning of the step, `C` the elastic
stiffness (`C:ε̇Δt = λΔt tr(ε̇) I + 2GΔt ε̇`), and `v` the unknown velocity of
this step. Boundary conditions are the usual velocity BCs (`bc.vbc_*`),
Winkler foundation, water loading, and traction BCs.

Classic DynEarthSol reaches quasi-static balance with FLAC-style dynamic
relaxation: scaled inertia (`control.inertial_scaling`) plus local damping,
advancing thousands of pseudo-dynamic steps. The PT solver replaces this with
an iteration *within* each time step that converges to the static solution and
only then advances the physical state.

## 2. Theory: the accelerated PT method

### 2.1 First-order vs. second-order pseudo-transient iteration

The naive ("first-order") PT method augments the static equation with a
pseudo-time derivative, `ρ̃ ∂v/∂τ = ∇·σ + ρg`, and integrates in pseudo-time
`τ` until the right-hand side vanishes. This is a diffusion-like process: an
error of wavelength `L` decays on the pseudo-time scale `L²`, so the iteration
count scales as `O((L/h)²)` — prohibitive for fine meshes.

The *accelerated* (second-order) method instead makes the system a **damped
wave equation**: errors propagate through the domain at a pseudo-wave speed
(crossing it in `O(L/h)` iterations) while being damped at a rate tuned so the
slowest, longest-wavelength mode decays optimally. The damping is not applied
to the momentum equation; it enters through a Maxwell-type relaxation in the
pseudo-time evolution of the stress. The iteration count then scales as
`O(L/h)` — the same scaling that makes multigrid-like methods fast, achieved
with purely local, explicitly parallel updates (ideal for GPUs).

### 2.2 Optimal numerical parameters

The scheme has three numerical parameters (Räss et al. 2022, Sec. 2.4):

| symbol | meaning | control name | default |
|---|---|---|---|
| `Re`  | numerical "Reynolds number": ratio of pseudo-inertia to viscous relaxation; sets the damping of the slowest mode | `control.PT_Re` | `3√10/2·π ≈ 14.93` (optimal for the Stokes-like problem) |
| `CFL` | pseudo-wave Courant number | `control.PT_CFL` | `0.9/√NDIMS` |
| `r`   | ratio of numerical bulk to shear modulus `K̃/G̃` | `control.PT_r` | `0.5` |

From these, with `L` a characteristic domain length (`control.PT_char_length`,
auto = largest mesh dimension) and `h` a local element size, the method defines

```
ρ̃    = Re · μ_ve / (V̂ L)         pseudo-density
Δτ    = CFL · h / V̂               pseudo-time step
G̃     = ρ̃ V̂² / (r + 2)           numerical shear modulus
```

`V̂` is a free velocity scale that **cancels** in the two combinations the
algorithm actually needs:

```
Δτ/ρ̃  = CFL · h · L / (Re · μ_ve)                    (velocity update)
G̃Δτ   = Re · CFL · h · μ_ve / ((r+2) · L)            (stress update)
```

`μ_ve` is the effective visco-elastic viscosity of the *physical* problem over
one time step (Räss et al. Eq. 35):

```
μ_ve = 1 / ( 1/(G·Δt) + 1/μ_s )
```

For DynEarthSol's elastic and elasto-plastic rheologies (`μ_s` effectively
infinite) this reduces to `μ_ve = G·Δt`.

### 2.3 The two updates

Each PT iteration `k` performs, in order:

**Stress relaxation** (`update_stress_PT()`, rheology.cxx). The stress is
relaxed toward the physical elastic target of the current velocity iterate:

```
τ*      = τ_old + C : ε̇(v^k) Δt
τ^{k+1} = ( τ^k + θ τ* ) / ( 1 + θ ),     θ_e = G̃Δτ / (G_e Δt)
```

This is Räss et al. Eq. 36 written for a total-stress, compressible
formulation. The essential point — and the reason a naive implementation
diverges — is that the relaxation rate must be built from the **numerical**
modulus `G̃`, not the physical one:

```
θ_e = Re·CFL·h_e / ((r+2)·L) · (μ_ve,e / (G_e Δt))  ≈  Re·CFL·h_e/((r+2)L)  ~  h/L  ≪ 1
```

The stress must build up *slowly* (a fraction `~h/L` of the physical increment
per iteration). That slow build-up is exactly what turns the coupled iteration
into a damped wave and allows the large velocity step below. Applying the full
physical increment `C:ε̇Δt` every iteration is equivalent to `θ = 1`, which
violates the pseudo-wave CFL condition by a factor `~L/h` and blows up within
tens of iterations.

**Plastic projection** (same kernel). After the relaxation, the stress is
projected back onto the yield surface (Mohr–Coulomb + tension cutoff) by
calling the standard return map with a *zero* strain increment. The projection
is non-expansive, so it does not destabilize the iteration, and its fixed
point satisfies both equilibrium and the yield condition. Plastic strain is
**not** accumulated during PT iterations (see §2.6).

**Velocity update** (`update_velocity_PT()`, fields.cxx):

```
v^{k+1}_i = v^k_i + dτ_ρ[i] · f_i ,
dτ_ρ[i]   = CFL · h_i · L · NODES_PER_ELEM / ( Re · μ_ve,i · volume_n[i] )
```

where `f_i` is the assembled nodal out-of-balance force (internal tractions +
gravity + boundary forces) and `volume_n[i]/NODES_PER_ELEM` is the nodal
volume. There is *no* velocity damping and *no* FLAC local damping in this
update — all dissipation comes from the stress relaxation, whose rate the
optimal `Re` controls.

### 2.4 Local pseudo-time stepping

Räss et al. derive the parameters for a uniform grid spacing `h`. On graded
unstructured meshes a single global `h = h_min` throttles the entire domain to
the pace of the worst element: a model with `h_min/L ~ 10⁻⁴` (one sliver in a
7 km domain) needed `>5·10⁴` iterations. DynEarthSol therefore uses **local**
factors, precomputed once per PT loop by `update_pt_params()` (geometry.cxx):

- per element: `G̃Δτ|_e = Re·CFL·h_e·μ_ve,e / ((r+2)·L)` with `h_e` the
  element's minimum height (stored in `var.PT_Gdtau_e`);
- per node: `dτ_ρ[i]` built from `h_i = min(h_e)` and `μ_ve,i = max(μ_ve,e)`
  over the node's adjacent elements (stored in `var.PT_dtau_rho`).

The min/max choices at nodes are conservative: at any element–node pair the
product of stress step and velocity step stays below the local stability
limit, including across element-size and material jumps. Each region then
relaxes at its own optimal rate, and the iteration count is governed by the
*typical* element count across the domain, not by the smallest element. In
practice this reduced a >50 000-iteration solve to a few hundred iterations.

### 2.5 What is (deliberately) frozen during PT iterations

- **The mesh.** The PT velocity is a relaxation iterate, not a physical
  velocity; advecting nodes with it (`coord += v Δt`) progressively distorts
  elements and destabilizes the iteration. The Lagrangian mesh update happens
  once, after convergence, with the converged velocity.
- **FLAC local damping** (`control.damping_option`). It is skipped while
  `PT_jump` is set; superposing it on the tuned wave dynamics degrades or
  destroys convergence.
- **Surface processes, thermal/hydraulic diffusion, plastic-strain
  bookkeeping** — all applied once per physical step, outside the PT loop.

### 2.6 Predictor–corrector structure of a time step

For `control.has_PT = yes` the main loop executes each step as:

1. Save `τ_old` (`copy_stress_PT` → `var.stress_old`).
2. Assemble forces of the `τ_old` state; record the initial residual and the
   characteristic force scale (§2.7).
3. `update_pt_params()` — recompute the local PT factors.
4. **PT loop** (up to `PT_max_iter`): `apply_vbcs → update_strain_rate →
   update_stress_PT (relax + project) → update_force → update_velocity_PT →
   residual check`.
5. **Corrector**: restore `τ_old`, then run the full physical
   `update_stress()` once with the converged velocity (followed by NMD if
   enabled). This is the only place where total strain, `plstrain`,
   `delta_plstrain`, viscosity, etc. are updated — so the constitutive
   bookkeeping is done exactly once per step and is consistent with the
   equilibrated velocity field. The pre-PT "predictor" stress update that the
   non-PT path performs is skipped entirely in PT mode: its plastic-strain
   increments would come from the stale previous-step velocity and be
   double-counted by the corrector.

The corrector's stress is (up to the relaxation tolerance) the same operation
the PT loop converged on, so the post-correction imbalance remains small.

### 2.7 Convergence criterion: unbalanced-force ratio

The residual is the L2 (RMS-per-DOF) norm of the **true out-of-balance
force**: everything assembled into `force` — internal tractions, gravity,
Winkler foundation, water loading, Neumann/stress BCs — taken *before* the
artificial damping, and with velocity-constrained DOFs masked out (their
reactions are not imbalance; the masking maps `bc.vbc_*` types to fixed
components per boundary).

Convergence is declared when

```
residual / force_scale < control.PT_relative_tolerance
```

where `force_scale` is the **characteristic gross force**: the RMS over DOFs
of `Σ_e |element nodal contribution|`, computed once per PT loop
(`calculate_characteristic_force()`). This is the FLAC-style *unbalanced force
ratio*: a value of `10⁻⁶` means the large opposing forces at every node cancel
to six significant digits — the same physical accuracy statement at every
step, for every model.

Normalizing by the *initial* residual instead (the previous behavior) fails in
both directions: a step that starts at equilibrium (e.g. any boundary-driven
model after its first step) has `residual_0` at round-off level, so the ratio
can never decrease by the tolerance factor and the loop spins for
`PT_max_iter`; a badly unbalanced start makes the test pass while large forces
remain. The attainable floor of the new ratio is machine precision
(`~10⁻¹⁶`), so tolerances down to `~10⁻¹²` are meaningful.

## 3. Practical usage

### 3.1 Minimal configuration

```ini
[control]
has_PT = yes
PT_max_iter = 50000          # safety cap per time step
PT_relative_tolerance = 1e-6 # unbalanced-force ratio (see below)
PT_info_interval = 100       # print residual every N PT iterations (0 = quiet)
```

Everything else has sensible defaults. The optional tuning knobs:

```ini
PT_Re = 14.93          # numerical Reynolds number (default 3*sqrt(10)/2*pi)
PT_CFL = 0.5196        # pseudo-wave CFL (default 0.9/sqrt(NDIMS))
PT_r = 0.5             # K~/G~ ratio
PT_char_length = 0     # characteristic length L; 0 = max(xlength,ylength,zlength)
```

You should rarely need to change these. If a model shows late-stage residual
oscillation or a very slow tail, the first thing to try is lowering `PT_CFL`
(e.g. 0.3); the second is increasing `PT_Re` slightly.

### 3.2 Reading the log

With `PT_info_interval > 0` each step prints:

```
PT start: initial residual = 1.13969e+11, force scale = 7.72484e+14
PT params: h_min=243.354 m, mu_ve_max=3.65030e+17 Pa*s, Gdtau=2.75137e+15 Pa*s
PT iter 0:   residual = ..., residual/residual_0 = ...
...
PT converged at iter 193: residual = 7.53511e+08 (residual/residual_0 = 0.000001)
```

- `force scale` is the gross-force denominator. `initial residual / force
  scale` tells you how far from equilibrium the step starts; values near
  `10⁻¹⁶` mean "already at equilibrium" and the loop exits immediately.
- `h_min`, `mu_ve_max`, `Gdtau` are global diagnostics (the solver itself uses
  the local per-element/per-node values).
- A healthy run converges in `O(100)` iterations per step (a few hundred for
  strongly graded meshes or after large load changes). If steps routinely hit
  `PT_max_iter`, see §3.4.

### 3.3 Choosing the tolerance

`PT_relative_tolerance` is a fraction of the *gross* force, not of the
per-step forcing. In a boundary-driven model where the driving stress
increment per step is a small fraction `s` of the lithostatic state (e.g.
`s ~ 5·10⁻⁵` for cm/yr extension of a 10-km crust with `Δt ~ 0.4 yr`), a
tolerance of `10⁻⁶` resolves each step's increment to `10⁻⁶/s ≈ 2%`. Tighten
to `10⁻⁷`–`10⁻⁸` if you need each increment resolved more finely — each decade
costs only a roughly fixed number of extra iterations (the convergence is
geometric), so this is cheap.

### 3.4 Performance characteristics and troubleshooting

- **Iteration count scales with `L/h_typical`**, thanks to local pseudo-time
  stepping — a few slivers in the mesh do not slow the solve down (they only
  relax their own neighborhood more carefully). Still, near-degenerate
  elements harm stress accuracy locally; if the log reports an `h_min` orders
  of magnitude below your nominal resolution, inspect the mesh.
- **Residual stalls at a plateau well above tolerance.** Likely causes, in
  order: (a) large regions at plastic yield — the coupled
  equilibrium-with-yield problem is genuinely harder; consider more steps
  (smaller `Δt`) or slightly relaxed tolerance; (b) an unbalanced boundary
  condition the residual legitimately measures (e.g. a traction BC
  inconsistent with the interior state).
- **Residual grows.** Should not happen with default parameters; if it does,
  reduce `PT_CFL`. Persistent growth indicates a configuration the stability
  analysis does not cover — please report it.
- **Steps that hit `PT_max_iter` leave the last (unconverged) velocity
  iterate in `vel`**, which then advects the mesh. A run where many steps do
  this will produce distorted meshes and remeshing failures — treat
  `PT_max_iter` as a safety net, not an operating mode.

### 3.5 Interaction with other features

| feature | behavior with `has_PT = yes` |
|---|---|
| Winkler foundation / water loading / traction BCs | included in both the force and the residual; converges normally |
| velocity BCs (`vbc_*` types 0–3) | enforced every PT iteration; constrained DOFs are excluded from the residual (types 4/5 are masked conservatively as fully constrained) |
| plasticity (`elasto-plastic`, `elasto-visco-plastic`) | yield enforced every PT iteration by projection; plastic strain accumulated once per step in the corrector |
| Maxwell viscosity | **not** relaxed during PT iterations (elastic target only); applied once per step in the corrector with viscosity lagged by one step |
| rate-and-state friction | not applied during PT iterations; PT with RSF rheologies is untested |
| remeshing (in-PT quality check) | supported: PT factors and `τ_old` are refreshed after the remesh; the force scale stays that of the pre-remesh mesh for the remainder of that step (a few-percent effect) |
| thermal/hydraulic diffusion, surface processes | suspended during PT iterations, applied once per step |
| `inertial_scaling`, `damping_option` | irrelevant / disabled inside the PT loop (used only by the non-PT dynamic-relaxation path) |

### 3.6 Implementation map

| function | file | role |
|---|---|---|
| `update_pt_params()` | geometry.cxx | per-element `G̃Δτ`, per-node `dτ_ρ`, diagnostics |
| `update_stress_PT()` | rheology.cxx | Eq. 36 relaxation + plastic projection |
| `update_velocity_PT()` | fields.cxx | accelerated velocity update |
| `copy_stress_PT()` | fields.cxx | save/restore `τ_old` |
| `calculate_characteristic_force()` | fields.cxx | gross-force scale for convergence test |
| `calculate_residual_force()` | fields.cxx | RMS out-of-balance force |
| `update_force()` (residual section) | fields.cxx | BC-inclusive, DOF-masked residual assembly |
| PT loops + predictor–corrector | dynearthsol.cxx | step orchestration |

## 4. References

- Räss, L., Utkin, I., Duretz, T., Omlin, S., and Podladchikov, Y. Y. (2022).
  Assessing the robustness and scalability of the accelerated pseudo-transient
  method. *Geoscientific Model Development*, 15, 5757–5786.
  https://doi.org/10.5194/gmd-15-5757-2022
- Cundall, P. A. (1987). Distinct element models of rock and soil structure.
  In *Analytical and Computational Methods in Engineering Rock Mechanics* —
  origin of the dynamic-relaxation / unbalanced-force-ratio practice used by
  FLAC and inherited by DynEarthSol's convergence metric.
- Duretz, T., Räss, L., Podladchikov, Y. Y., and coworkers' PT viscoplasticity
  papers for the relax-then-project treatment of plastic yield inside
  pseudo-transient iterations.
