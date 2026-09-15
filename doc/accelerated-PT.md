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

$$\nabla \cdot \sigma + \rho g = 0$$

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

**Where the "order" hides.**  The classification cannot be read off the
momentum equation alone.  Both variants can be written in the same form,

$$\tilde{\rho}\partial v/\partial\tau = \nabla \cdot \sigma + f,$$

which is first order in time and second order in space *as written* — yet
one variant is parabolic and the other hyperbolic.  The "order" refers to
the highest pseudo-time derivative acting on the primary field in the
equivalent single-field equation, i.e. **after the stress has been
eliminated** — and what elimination yields depends entirely on what σ
depends on:

- If σ is **memoryless** — recomputed fresh each iteration as an algebraic
  function of the *current* field — elimination introduces no new time
  derivative.  One time derivative survives: parabolic.
- If σ has **memory** — updated from its *own previous value* by its own
  evolution equation — eliminating it (differentiate the momentum equation
  in pseudo-time, substitute the stress update) produces a second time
  derivative: hyperbolic.

Equivalently, count the independent state variables marched in pseudo-time.
One (the field itself, stress slaved to it) → first order.  Two (field
*and* stress, each carrying history) → second order.  Two coupled
first-order ODEs are one second-order system.

**First-order PT: the constitutive law is enforced exactly, every iteration.**
The stress is algebraically slaved to the current velocity — recomputed in
full as `σ = τ_old + C:ε̇(v^k)Δt` each iteration (`τ_old` is the converged
stress of the *previous physical step*, a constant during the iteration —
it contributes no pseudo-time dynamics).  Substituting this into the
momentum update leaves a single pseudo-time derivative,

$$\tilde{\rho}\partial v/\partial\tau = \nabla \cdot (C\Delta t : \dot{\varepsilon}(v)) + (\text{known terms}),$$

a **parabolic** (diffusion) equation in pseudo-time with diffusivity
`D ~ GΔt/ρ̃`.  There is no wave: the stress has no independent dynamics, so
nothing stores and returns kinetic information.  This is the natural form
for Stokes flow (the setting of Räss et al.), where velocity *is* the
equilibrium unknown and `τ = 2μ ε̇(v)` is instantaneous.  Explicit stability caps the step at `Δτ ≤ h²/(2d·D)` (here *d* is
the number of spatial dimensions).  To see why the slowest mode decays so
slowly, consider the 1D case.  For a mode of wavenumber *k*, explicit Euler
multiplies its amplitude each iteration by the factor

$$g(k) = 1 - D k^2 \Delta\tau,$$

so the fractional decay per iteration is $1 - g(k) = D k^2 \Delta\tau$.
Stability requires $|g| \le 1$ for all *k*; the binding constraint comes from
the fastest (highest-*k*) mode, $k_{\max} = \pi/h$ (Nyquist), giving

$$\Delta\tau \lesssim \frac{1}{D k_{\max}^2} \sim \frac{h^2}{D}.$$

At this ceiling, the amplification factor for the *slowest* mode $k_{\min} = \pi/L$ is

$$g(k_{\min}) = 1 - D\left(\frac{\pi}{L}\right)^{2}\frac{h^2}{D}
            = 1 - \left(\frac{\pi h}{L}\right)^2
            \approx 1 - \left(\frac{h}{L}\right)^2,$$

so each iteration removes only a fraction $\sim(h/L)^2$ of the slow-mode error,
and the iteration count scales as `O((L/h)² · ln(1/tol))` — prohibitive for
fine meshes.  Information spreads the way diffusion spreads: incoherently,
over a radius `~h√k` after `k` iterations.  To see why, note that after $k$
iterations of size $\Delta\tau$, the elapsed pseudo-time is $k\Delta\tau$, so
a diffusive disturbance spreads over

$$r \sim \sqrt{D \cdot k\Delta\tau} \sim \sqrt{D \cdot k \cdot \frac{h^2}{D}} = h\sqrt{k}.$$

This is the random-walk result: each degree of freedom nudges its neighbours
by a small amount each step, with no organised wave front, so the spread
grows as $\sqrt{k}$ rather than linearly in $k$.  To reach across the domain
($r \sim L$) requires $k \sim (L/h)^2$, confirming the scaling above.  The
second-order scheme escapes this because its stress dynamics carry a genuine
wave front that advances $\sim h$ per iteration, crossing the domain
coherently in $O(L/h)$ steps.

**Second-order PT: the stress gets its own pseudo-time evolution.**  Instead
of enforcing the constitutive law, each iteration *relaxes toward* it (the
θ-update of §2.4 — a Maxwell element in pseudo-time).  The stress now
depends on its own previous value — it has memory — so eliminating it
produces a **second** pseudo-time derivative of the primary field: the
damped wave equation analyzed in §2.2.  Information travels
coherently — one element per iteration, radius `~h·k·CFL` — so a domain
crossing costs `O(L/h)` iterations, and a correctly tuned damping makes the
total count `O((L/h) · ln(1/tol))`, achieved with purely local, explicitly
parallel updates (ideal for GPUs).

**Why FLAC-style dynamic relaxation is second-order despite its first-order
momentum equation.**  DynEarthSol's native mode marches
`ρ ∂v/∂t = ∇·σ + f` — first order in v as written.  But its hypoelastic
update `∂σ/∂t = C:ε̇(v)` makes σ a second history-carrying state variable,
not a function of the current field.  Differentiating the momentum equation
in time and substituting gives

$$\rho\,\partial^2 v/\partial t^2 = \nabla \cdot (C : \nabla^s v),$$

a genuine wave equation (equivalently `ρ ü = ∇·(C:∇ˢu)` in displacement):
the second time derivative was hidden in the stress memory.  So dynamic
relaxation is *also* a second-order method in this taxonomy, and what
distinguishes *accelerated* PT is not the wave character but the parameter
choice.  Dynamic relaxation runs with the physical moduli and heuristic
damping (`inertial_scaling`, local damping), generally far from optimal;
the accelerated method derives `ρ̃`, `G̃` and the damping rate analytically
from the dispersion relation so that the slowest mode is critically damped
(§2.2).  In that sense Räss et al.'s "acceleration" of the Stokes solver
consists of promoting τ from a slaved quantity to a second state variable —
recovering the structure FLAC-style codes had all along — and then tuning
it.

### 2.2 Dispersion analysis: how critical damping is enforced

Strip the scheme to its 1D scalar essence (elastic limit, `μ_ve = GΔt`).
Per unit pseudo-time, the two updates read

$$\begin{aligned}
\tilde{\rho}\partial v/\partial\tau &= \partial\sigma/\partial x && \text{(momentum — no damping)}\\
\partial\sigma/\partial\tau &= \tilde{G}\partial v/\partial x - \eta(\sigma - \tau^*) && \text{(stress relaxation)}
\end{aligned}$$

The second line is the discrete θ-update of §2.4 in continuous form, with
relaxation rate `η = θ/Δτ = G̃/(GΔt)`.  The numerical modulus `G̃` plays two
roles at once: it is the spring of the pseudo-wave *and*, through η, the
dial on the Maxwell-type dissipation.  That is the structural trick of the
method — all dissipation lives in the constitutive update, none in the
momentum equation.

**Damped wave equation.**  Differentiate the stress equation in τ, substitute
`∂v/∂τ` from momentum (the target τ* is constant for the error), and the
perturbation about the converged state obeys

$$\partial^2\sigma/\partial\tau^2 + \eta\partial\sigma/\partial\tau = \hat{V}^2\partial^2\sigma/\partial x^2, \qquad \hat{V}^2 \propto \tilde{G}/\tilde{\rho}.$$

**Dispersion relation.**  A Fourier mode `σ ∝ e^{λτ} e^{ikx}` gives

$$\lambda^2 + \eta\lambda + \hat{V}^2 k^2 = 0 \quad\Rightarrow\quad \lambda = -\eta/2 \pm \sqrt{\eta^2/4 - \hat{V}^2 k^2}.$$

Each error component, labeled by its wavenumber `k`, decays as:

- **underdamped** (`η < 2V̂k`): complex roots, `Re(λ) = −η/2` — the mode
  oscillates (a wave sloshing through the domain) while decaying at rate
  `η/2`, *independent of k*;
- **overdamped** (`η > 2V̂k`): both roots are real; expanding
  $\sqrt{\eta^2/4 - \hat{V}^2k^2} \approx \eta/2 - \hat{V}^2k^2/\eta$
  gives the slow root $\lambda \approx -\hat{V}^2k^2/\eta$ — a
  diffusive rate `−Dk²` with `D = V̂²/η`.  Overdamping does not merely slow
  a mode down; it reverts that mode to first-order (parabolic) behavior.
  In the limit `η → ∞` (stress slaved to velocity) the whole spectrum is
  diffusive and the `O((L/h)²)` method of §2.1 is recovered;
- **critical** (`η = 2V̂k`): double root `λ = −V̂k` — the fastest decay that
  wavenumber can achieve.

**Critical damping as a min–max choice.**  There is one knob, η, and a
spectrum of modes from `k_min = π/L` (the largest error structure the domain
can hold) up to `k_max ~ π/h`; one η can critically damp only one k.  The
iteration converges at the pace of its *slowest* mode, so η is chosen to
maximize the minimum decay rate.  Underdamped modes (`k > η/2V̂`) decay at
`η/2` — raising η helps them; overdamped modes decay at `V̂²k²/η` — raising η
hurts them.  The optimum places the boundary between the regimes exactly at
the smallest wavenumber:

$$\eta^* = 2\hat{V} k_{\min} = \frac{2\pi\hat{V}}{L}.$$

Then the `k_min` mode is exactly critical (rate `V̂ k_min`) and every other
mode is underdamped with the *same* rate `η*/2 = V̂ k_min`: the entire
spectrum decays uniformly, the best achievable bound.  Nudge η either way
and the minimum drops.

**Iteration count.**  The discrete scheme must also resolve the fastest wave
stably — that is where CFL enters: `Δτ = CFL·h/V̂`.  The decay per iteration
of every mode is then

$$|\lambda|\,\Delta\tau = \hat{V} k_{\min} \cdot \frac{\mathrm{CFL}\cdot h}{\hat{V}} = \frac{\pi\cdot\mathrm{CFL}\cdot h}{L}$$

so tolerance ε is reached in `N ≈ ln(1/ε) · L/(π·CFL·h)` iterations.  Note that V̂ cancels here exactly
as it does in the implementation: only the combinations `Δτ/ρ̃` and `G̃Δτ`
are ever needed, which is why the code never chooses a wave speed.

**Where Re comes from.**  From the definitions in §2.3,

$$\eta = \frac{\tilde{G}}{G\Delta t} = \frac{\tilde{\rho}\hat{V}^2/(r+2)}{\mu_{ve}} = \frac{\mathrm{Re}\cdot\hat{V}}{(r+2)\cdot L}.$$

Setting this equal to `η* = 2πV̂/L` gives

$$\mathrm{Re}^* = 2\pi(r+2) \approx 15.7 \qquad (r = 0.5),$$

within 5% of the default `Re = 3√10·π/2 ≈ 14.93`, which comes from the same
analysis carried out for the full coupled vector system (velocity +
pressure, with r-dependent P- and S-wave branches) instead of the scalar
scalar damped wave equation.  The identity is the point: **`Re` is the
damping-to-wave-crossing ratio `η·L/V̂` (up to the `(r+2)` factor), and its
optimal value is the one that places the slowest mode at critical damping.**

### 2.3 Optimal numerical parameters

The scheme has three numerical parameters (Räss et al. 2022, Sec. 2.4):

| symbol | meaning | control name | default |
|---|---|---|---|
| `Re`  | numerical "Reynolds number": ratio of pseudo-inertia to viscous relaxation; sets the damping of the slowest mode at its critical value (derivation in §2.2) | `control.PT_Re` | `3√10/2·π ≈ 14.93` (optimal for the Stokes-like problem) |
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

$$\mu_{ve} = \frac{1}{1/(G\Delta t) + 1/\mu_s}$$

For DynEarthSol's elastic and elasto-plastic rheologies (`μ_s` effectively
infinite) this reduces to `μ_ve = G·Δt`.

### 2.4 The two updates

Each PT iteration `k` performs, in order:

**Stress relaxation** (`update_stress_PT()`, rheology.cxx). The stress is
relaxed toward the physical elastic target of the current velocity iterate:

```
τ*      = τ_old + C : ε̇(v^k) Δt
τ^{k+1} = ( τ^k + θ τ* ) / ( 1 + θ ),     θ_e = G̃Δτ / (G_e Δt)
```

**Derivation from Räss et al.**  Eq. 34 of Räss et al. (2022) is the
dual-time PT augmentation of the visco-elastic Stokes stress equation
(Eq. 33): it adds a pseudo-time derivative (1/2G̃) ∂τ/∂τ alongside the
physical Maxwell relaxation:

$$\frac{1}{2\tilde{G}}\frac{\partial\tau_{ij}}{\partial\tau} + \frac{\tau_{ij} - \hat{\tau}_{ij}}{2G\Delta t} + \frac{\tau_{ij}}{2\mu_s} = \dot{\varepsilon}_{ij}(v)$$

where $\hat{\tau}_{ij} = \tau_\text{old}$.  Discretising $\partial\tau/\partial\tau$ forward-in-pseudo-time and
collecting $\tau^{k+1}$:

$$\tau^{k+1}\!\left[\frac{1}{2\tilde{G}\Delta\tau} + \frac{1}{2\mu_{ve}}\right] = \dot{\varepsilon}_{ij} + \frac{\tau^k}{2\tilde{G}\Delta\tau} + \frac{\tau_\text{old}}{2G\Delta t}$$

using $1/\mu_{ve} = 1/(G\Delta t) + 1/\mu_s$ (Eq. 35).  Multiplying through by $2\tilde{G}\Delta\tau$
and letting $\theta = \tilde{G}\Delta\tau/\mu_{ve}$:

$$\tau^{k+1}(1+\theta) = \tau^k + 2\tilde{G}\Delta\tau \dot{\varepsilon}_{ij} + \theta\cdot\frac{\mu_{ve}}{G\Delta t}\cdot\tau_\text{old}$$

In the **elastic limit** ($\mu_s \to \infty$, $\mu_{ve} = G\Delta t$, $\theta = \tilde{G}\Delta\tau/(G\Delta t)$):

$$\begin{aligned}
\tau^{k+1}(1+\theta) &= \tau^k + \theta\cdot 2G\Delta t \dot{\varepsilon}_{ij} + \theta\cdot\tau_\text{old}\\
                     &= \tau^k + \theta\cdot(\tau_\text{old} + 2G\Delta t \dot{\varepsilon}_{ij})\\
                     &= \tau^k + \theta\cdot\tau^*
\end{aligned}$$

which gives the update formula above exactly.  DynEarthSol's elasto-plastic
rheology satisfies this limit ($\mu_s \to \infty$).  For a true Maxwell viscoelastic
run (finite $\mu_s$) the term $\theta\cdot(\mu_{ve}/G\Delta t)\cdot\tau_\text{old}$ no longer equals $\theta\cdot\tau_\text{old}$,
and a modified $\theta = \tilde{G}\Delta\tau/\mu_{ve}$ (using the effective $\mu_{ve}$) would keep the same
$(\tau^k + \theta\tau^*)/(1+\theta)$ form but with a different $\tau^*$ definition.

The derivation above is for the deviatoric Räss et al. formulation
(incompressible, deviatoric τ).  DynEarthSol uses total stress with
compressibility: `C:ε̇Δt` includes the volumetric term `λΔt tr(ε̇) I`.
The algebra is identical per component; only the physical elastic target
τ* is broader.

The essential point is that the relaxation rate must be built from the **numerical**
modulus `G̃` (see §2.3), not the physical one:

```
θ_e = Re·CFL·h_e / ((r+2)·L) · (μ_ve,e / (G_e Δt)),
```
where h_e and G_e are element size and element stiffness. 
In the elastic limit, μ_ve,e = G_e Δt. So,

```
θ_e ≈  Re·CFL·h_e/((r+2)L)  ~  h/L  ≪ 1.
```
The stress must build up *slowly* (a fraction `~h/L` of the physical increment
per iteration). That slow build-up is exactly what turns the coupled iteration
into a damped wave and allows the large velocity step in the accelerated PT. 
Applying the full physical increment `C:ε̇Δt` every iteration is equivalent to 
`θ = 1`, which violates the pseudo-wave CFL condition by a factor `~L/h`. 
Such a test case blew up within tens of iterations.

**Plastic projection** (`update_stress_PT()`, rheology.cxx but using one the existing functions, `elasto_plastic` or `elastic_plastic2D`). After the relaxation, 
only the stress is projected back onto the yield surface 
(Mohr–Coulomb + tension cutoff) by calling the standard return mapping. Plastic strain is **not** accumulated during PT iterations (see §2.7).

The term *projection* is precise: the yield surface (a Mohr–Coulomb cone in
stress space) bounds a convex feasible region, and the return map finds the
closest stress on that surface: i.e., the orthogonal projection $P(\sigma^{trial})$
onto the convex set. Thus, the projection is *non-expansive*, which means that for 
any two stress states, the projection never increases the distance between them,

$$\|P(\sigma^a) - P(\sigma^b)\| \leq \|\sigma^a - \sigma^b\|,$$

a standard result for projections onto convex sets.  This is why the
projection cannot amplify the error $\|\sigma^k - \sigma^*\|$ and therefore
cannot destabilize the convergence driven by the relaxation step.

The projection operation also has a *fixed point*. The fixed point of the combined 
operator (elastic or viscoelastic guess then project) is the state
where neither step moves the stress: the guessed stress is stationary when momentum
is balanced and the constitutive law is satisfied, and the projection is
stationary when the stress is already on the yield surface.  Both hold
simultaneously only at the physical solution, so the fixed point satisfies both
equilibrium and the yield condition.  The argument relies on the yield surface
being convex, which holds for Mohr–Coulomb and Drucker–Prager.

**Velocity update** (`update_velocity_PT()`, fields.cxx):

$$v_i^{k+1} = v_i^k + \Delta\tau_{\rho,i}\, f_i,$$
where $\Delta\tau_{\rho,i}$ denotes $\Delta\tau/\tilde{\rho}$ in §2.3.

$$\Delta\tau_{\rho,i} = \frac{\mathrm{CFL}\cdot h_i \cdot L \cdot N_\text{npe}}{\mathrm{Re}\cdot\mu_{ve,i}\cdot V_i}$$

where $f_i$ is the assembled nodal out-of-balance force (internal tractions +
gravity + boundary forces), $V_i = \texttt{volume\_n}[i]$ is the summed nodal
volume, and $N_\text{npe} = \texttt{NODES\_PER\_ELEM}$ so that
$V_i/N_\text{npe}$ is the nodal volume. There is *no* velocity damping and
*no* FLAC local damping in this update — all dissipation comes from the stress
relaxation, whose rate the optimal `Re` controls.

### 2.5 Local pseudo-time stepping

Räss et al. derive the parameters for a uniform grid spacing `h`. On graded
unstructured meshes a single global `h = h_min` throttles the entire domain to
the pace of the worst element: a model with `h_min/L ~ 10⁻⁴` (e.g., one sliver in a
7 km domain) needed `>5·10⁴` iterations. DynEarthSol therefore uses **local**
factors, precomputed once per PT loop by `update_pt_params()` (geometry.cxx):

- per element: `G̃Δτ|_e = Re·CFL·h_e·μ_ve,e / ((r+2)·L)` with `h_e` the
  element's minimum height (stored in `var.PT_Gdtau_e`);
- per node: `dτ_ρ[i]` (=Δτ/ρ̃ ) built from `h_i = min(h_e)` and `μ_ve,i = max(μ_ve,e)`
  over the node's adjacent elements (stored in `var.PT_dtau_rho`).

The min/max choices at nodes are conservative: at any element–node pair the
product of stress step and velocity step stays below the local stability
limit, including across element-size and material jumps. Each region then
relaxes at its own optimal rate, and the iteration count is governed by the
*typical* element count across the domain, not by the smallest element. In
practice this reduced a >50 000-iteration solve to a few hundred iterations.

### 2.6 What is (deliberately) frozen during PT iterations

- **The mesh.** The PT velocity is a relaxation iterate, not a physical
  velocity; advecting nodes with it (`coord += v Δt`) progressively distorts
  elements and destabilizes the iteration. The Lagrangian mesh update happens
  once, after convergence, with the converged velocity.
- **FLAC local damping** (`control.damping_option`). It is skipped while
  `PT_jump` is set; superposing it on the tuned wave dynamics degrades or
  destroys convergence.
- **Surface processes, thermal/hydraulic diffusion, plastic-strain
  bookkeeping** — all applied once per physical step, outside the PT loop.

### 2.7 Predictor–corrector structure of a time step

For `control.has_PT = yes` the main loop executes each step as:

1. Save `τ_old` (`copy_stress_PT` → `var.stress_old`).
2. Assemble forces of the `τ_old` state; record the initial residual and the
   characteristic force scale (§2.8).
3. `update_pt_params()` — recompute the local PT factors.
4. **Predictor — PT loop** (up to `PT_max_iter`): `apply_vbcs →
   update_strain_rate → update_stress_PT (relax + project) → update_force →
   update_velocity_PT → residual check`.  The stress update here uses only
   the elastic target `τ* = τ_old + C:ε̇(v^k)Δt`; viscous evolution, plastic
   strain accumulation, and other constitutive bookkeeping are deliberately
   omitted so that the PT stress can be relaxed freely across many iterations.
   The loop converges to the velocity field that balances the forces implied
   by the full elastic constitutive law.
5. **Corrector**: restore `τ_old`, then run the full physical
   `update_stress()` once with the converged velocity (followed by NMD if
   enabled). This is the only place where total strain, `plstrain`,
   `delta_plstrain`, viscosity, etc. are updated — so the constitutive
   bookkeeping is done exactly once per step and is consistent with the
   equilibrated velocity field.

The corrector's stress is (up to the relaxation tolerance) the same operation
the PT loop converged on, so the post-correction imbalance remains small.

### 2.8 Convergence criterion: unbalanced-force ratio

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

### 2.9 The physical time step in PT mode

In dynamic relaxation the physical step `Δt` does double duty: it is the
*stability* step of the pseudo-time march (the mass-scaled wave CFL,
`dt_elastic = 0.5·h_min/(v_bc·inertial_scaling)`, usually the binding limit
in quasi-static runs) **and** the *physical* step for constitutive updates
and mesh advection.  PT separates those roles: the iteration's stability
comes entirely from its own pseudo-time parameters, so `dt_elastic` no
longer has a job.

In PT mode, `compute_dt_PT()` therefore caps `Δt` only by the physics the
outer loop integrates explicitly:

- **mesh advection** — `0.5·h_min / max(v_actual, v_bc)`, so no node moves
  more than half the smallest element per step;
- **Maxwell relaxation** — `0.5·visc_min/G`;
- **thermal / hydraulic diffusion** when those (explicit) updates are on;
- **plastic weakening** — see below;
- times `dt_fraction`; `fixed_dt` overrides everything.  If no finite limit
  applies, it falls back to the DR step.

The weakening limit exists because piecewise-linear strain softening is
integrated forward-Euler in the weakening variable: the yield parameters
are evaluated at the *start-of-step* plastic strain, the return map's
`hardn` term linearizes only the cohesion slope *within the current
segment*, friction-angle weakening inside the increment is ignored
entirely, and crossing a segment boundary (`pls0` or `pls1`) invalidates
even that linearization.  Under DR the tiny `dt_elastic` masked all of
this; with the PT step, an advection-bound `Δt` can push a nucleating
shear band across its *entire* ramp in one step using intact-strength
parameters.  The cap is adaptive and evaluated **per element**:

```
dt ≤ Δt_prev · min_e( A_e / Δε_pl,e )
```

where `Δε_pl,e` is the previous step's `delta_plstrain` (elements with
none impose no limit) and `A_e` is the element's *allowance* — the
admissible increment given where it sits on the weakening curves of its
constituent materials (`pls_weakening_allowance()`, evaluated per material
present in the element via `elemmarkers`, so mixed-material elements whose
effective curve has breakpoints from several ramps are handled):

- **on a ramp** (`pls0 ≤ pls < pls1`): `A = f·(pls1−pls0)` — at most a
  fraction `f` of the ramp per step;
- **below a ramp** (`pls < pls0`): `A = (pls0−pls) + f·(pls1−pls0)` — the
  flat run-up is free, plus the same entry allowance, so distant elements
  never throttle `dt` and the cap can never collapse to zero at a
  breakpoint;
- **saturated** (`pls ≥ pls1`) or **zero-width ramp** (`pls0 = pls1`, a
  parameter jump no step size can resolve): no limit from that material.

Since `Δε_pl ≈ rate·Δt_prev`, the cap settles directly at
`Δt ≈ A/rate` for the tightest element — it tightens while a band
actively weakens and relaxes once the band saturates.  The one-step
reaction lag at band *initiation* is absorbed by the safety fraction
(`f = PT_dpls_fraction`, default 0.1; 0 disables).

The resulting `Δt` is typically orders of magnitude larger than the DR step
(×`inertial_scaling` when advection binds).  This is safe for the *solver*
because `Δt` cancels out of the PT convergence rate by construction: both
pseudo-time factors contain `μ_ve` in opposite positions
(`Δτ/ρ̃ ∝ 1/μ_ve`, `G̃Δτ ∝ μ_ve`), so the pseudo-wave CFL and the
per-iteration decay `~π·CFL·h/L` (§2.2) are `Δt`-free, as is the
relaxation weight `θ_e = Re·CFL·h_e/((r+2)L)`.  A larger `Δt` only enlarges
the load increment per step, which costs iterations *logarithmically* (a
larger starting residual relative to the gross force) — the rate is
untouched.  What a large `Δt` does affect is **accuracy**: the plastic
path is resolved in coarser stress increments (see §2.10 for the
consequence), and viscoelastic elements shift along the `μ_ve` spectrum
between `GΔt` and `μ`.  Throttle with `dt_fraction` when that matters.

Dynamic-relaxation mode is untouched: all of this is gated on `has_PT`,
and `isostasy_adjustment` (which marches DR internally) keeps the DR step
even in PT runs.

### 2.10 Stagnation exit: the nonlinear residual floor

The convergence theory of §2.2 is linear.  With plasticity, each iteration
ends in a non-expansive yield projection (§2.4), and when the load
increment per step is large — precisely the regime the PT time step of
§2.9 enables — a sizeable population of elements sits *at* yield.  The
iteration can then enter a **limit cycle**: elements alternate between the
elastic relaxation pulling stress outside the yield surface and the
projection snapping it back, and the residual oscillates on a plateau
instead of decaying.  The plateau level (typically `10⁻⁶`–`10⁻⁴` of gross
force, rising with increment size) is the accuracy the nonsmooth problem
admits at that increment; iterating past it accomplishes nothing, and
without an exit every such step burns all of `PT_max_iter`.

The loop therefore detects stagnation: it tracks the best residual seen
and, at the end of each *check window*, tests whether it improved by at
least 0.1%; two consecutive failed windows end the loop with a warning
that reports the achieved level.  Two design points matter:

- **The window is measured in domain crossings, not iterations.**  A
  near-critically damped iteration *rings*: the residual oscillates while
  its envelope decays, so the best-so-far improves in bursts, roughly once
  per ring period — which is of the order of a pseudo-wave domain
  crossing, `L/(CFL·h)` iterations.  A short fixed window (say 200
  iterations) reads the flat stretch between bursts as stagnation and
  exits while the envelope is still falling (observed: exits during 4%-
  per-200-iterations decay).  The auto window is two crossings,
  `2L/(CFL·h_mean)`, with `h_mean` the mean element height — not `h_min`,
  because with local pseudo-time stepping (§2.5) the crossing time is set
  by typical elements, not the worst sliver.
- **The trackers reset after an in-PT remesh**, where the residual
  legitimately jumps.

The exit is controlled by `control.PT_stagnation_window`
(`-1` = auto as above, `0` = disabled, `N` = fixed window of `N`
iterations).  Steps that can reach the tolerance still do — the exit needs
at least two windows of no progress to fire.  If the reported stagnation
levels are higher than you can accept, reduce the increment size
(`dt_fraction`): smaller steps mean smaller excursions past yield and a
lower floor.

### 2.11 Periodic retuning and Rayleigh-quotient adaptive Re

The analytical parameters Re, CFL, r are derived for a homogeneous medium
(§2.3).  For problems with spatially varying viscosity (plasticity, large
thermal contrasts) the effective λ_min of the assembled system can differ
substantially from the value the uniform-medium formula assumed, causing either
overdamping (wasted iterations) or convergence failure.  Two coordinated
mechanisms, both gated on `PT_retune_interval > 0`, adapt the parameters
during the iteration.

**Periodic `update_pt_params()`.** Every `PT_retune_interval` PT iterations
the per-element `G̃Δτ|_e` and per-node `dτ_ρ[i]` factors are recomputed.
During elastic phases these are nearly constant;
during progressive yielding μ_ve,e softens and the local stepping factors
follow.  Cost per call is O(nelem + nnode), ~1–2% of one iteration at interval 100.

**Rayleigh-quotient Re estimate.** Every `PT_retune_interval` iterations
`rayleigh_update_Re()` also estimates λ_min from the ratio of iterate
differences over the preceding window (Duretz et al., 2026, Sect. 3):

$$\lambda_\min \approx \frac{-\sum_{i,j} \Delta v_{ij}\,\Delta f_{ij}}{\sum_i \frac{1}{\Delta\tau_{\rho,i}}\sum_j (\Delta v_{ij})^2}$$

where $\Delta\tau_{\rho,i}$ is the per-node pseudo-time step (`var.PT_dtau_rho[i]`),
Δv = v_current − v_snapshot and Δf = f_current − f_snapshot captured at
the last retune call.  The optimal Re is then

$$\mathrm{Re}_\text{new} = 2\sqrt{\lambda_\min \cdot \lambda_\max} \cdot \frac{(r+2)\,L\,\bar{G}\,\Delta t}{\mathrm{CFL}\,\bar{h}\,\bar{\mu}_{ve}}$$

with $\lambda_\max = \mathrm{CFL}^2/(r+2)$ from the stability limit.  To limit overshoot during
the first few windows (when the iterate differences may not yet represent the
slow-mode eigenvalue), Re_new is clamped to [0.5, 2.0] × `PT_Re`.

**Retune guard.** The retune is skipped when `l2_residual / force_scale <
100 × PT_relative_tolerance` (= 1e-4 with the default tolerance).  Near
convergence the quotient would estimate numerical noise rather than physics.

**Resetting between PT loops.** `PT_Re_adaptive` is reset to `param.PT_Re` at
the start of each PT loop (init and main), so adaptation does not carry across
physical time steps.

**Observed behavior and benchmark.** In a 2D elasto-plastic benchmark
(`examples/pt-retune-demo.cfg`, 5 time steps, `PT_max_iter = 10000`):

- Without retuning: the init loop and step 1 both fail to converge (hit
  max_iter = 10000).  Total iterations: 24 351.
- With `PT_retune_interval = 100`: the init loop and step 1 converge at 5979
  and 5266 iterations respectively.  Total iterations: 16 488 (**−32%**).

The adaptive mechanism enables convergence by reducing Re from the default
~14.9 to ~7.5 (lower clamping bound) once slow modes dominate, which gives
larger `dτ_ρ` (more aggressive velocity stepping).  Step 2 is ~39% slower
because its first 100-iteration window captures a fast-mode transient (200×
residual drop), biasing λ_min high; this is outweighed by the savings on the
previously-failing loops.

## 3. Practical usage

### 3.1 Minimal configuration

```ini
[control]
has_PT = yes
PT_max_iter = 50000          # safety cap per time step
PT_relative_tolerance = 1e-6 # unbalanced-force ratio (see below)
PT_info_interval = 100       # print residual every N PT iterations (0 = quiet)
PT_retune_interval = 100     # re-evaluate PT params every N iterations (0 = off)
```

Everything else has sensible defaults. The optional tuning knobs:

```ini
PT_Re = 14.93              # numerical Reynolds number (default 3*sqrt(10)/2*pi)
PT_CFL = 0.5196            # pseudo-wave CFL (default 0.9/sqrt(NDIMS))
PT_r = 0.5                 # K~/G~ ratio
PT_char_length = 0         # characteristic length L; 0 = max(xlength,ylength,zlength)
PT_stagnation_window = -1  # early exit on residual plateau (§2.10);
                           # -1 = auto (2L/(CFL*h_mean)), 0 = off, N = fixed
PT_dpls_fraction = 0.1     # max fraction of the weakening ramp (pls1-pls0)
                           # traversed per step (§2.9); 0 = off
```

Note that in PT mode the physical time step is no longer limited by the
dynamic-relaxation stability step (§2.9) — expect `Δt` (and hence the load
increment per step) to be orders of magnitude larger than in a DR run of
the same model.  `dt_fraction` and `fixed_dt` still apply and are the way
to trade step size against plastic-path accuracy.

You should rarely need to change these.  The one parameter worth a conscious
choice is `PT_char_length`: it stands for the longest error wavelength the
damping is tuned to (`k_min = π/L`, §2.2), so it should approximate the
largest linear dimension of the *deforming* region.  The default (largest
mesh dimension) is right unless deformation is confined to a small part of a
large domain.  `PT_CFL` is a stability margin — lower it (e.g. 0.3) only if
the residual grows.  For the under-/over-damping symptoms and which way to
move `PT_char_length` or `PT_Re`, see §3.4.

### 3.2 Reading the log

With `PT_info_interval > 0` each step prints:

```
PT start: initial residual = 1.13969e+11, force scale = 7.72484e+14
PT params: h_min=243.354 m, mu_ve_max=3.65030e+17 Pa*s, Gdtau=2.75137e+15 Pa*s
PT iter 0:   residual = ..., residual/residual_0 = ...
...
PT converged at iter 193: residual = 7.53511e+08 (residual/residual_0 = 0.000001)
```

With `PT_retune_interval > 0`, retune calls appear on stderr between the
iteration lines:

```
PT iter 100: residual = 6.33204e+07, residual/force_scale = 8.19e-08
[retune] λ_min=1.24e-03, Re=2.21→clamped to 7.45, PT_Gdtau updated
PT iter 200: residual = 1.41283e+07, residual/force_scale = 1.83e-08
```

The `→clamped` token appears when `Re_new` falls outside the [0.5, 2] ×
`PT_Re` guard range (§2.11).

- `force scale` is the gross-force denominator. `initial residual / force
  scale` tells you how far from equilibrium the step starts; values near
  `10⁻¹⁶` mean "already at equilibrium" and the loop exits immediately.
- `h_min`, `mu_ve_max`, `Gdtau` are global diagnostics (the solver itself uses
  the local per-element/per-node values).
- A healthy run converges in `O(100)` iterations per step (a few hundred for
  strongly graded meshes or after large load changes). If steps routinely hit
  `PT_max_iter`, see §3.4.
- A step may instead end with

  ```
  PT stagnated at iter 8000: residual = 2.83557e+09 (residual/residual_0 = 0.000004 > tolerance 1e-06); accepting current state.
  ```

  — the residual plateaued (§2.10) and the loop accepted the achieved
  level, here `4·10⁻⁶` of gross force.  Occasional stagnation slightly
  above the tolerance is normal for elasto-plastic runs with large steps;
  if the reported levels are too high for your purposes, reduce
  `dt_fraction`.

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
- **Residual oscillates ("rings") while decaying slowly** → the slow modes
  are underdamped: `PT_char_length` is larger than the actual deforming
  region (or `PT_Re` is too small).  **Residual drops fast at first, then
  crawls through a long smooth tail** → the slowest mode is overdamped:
  `PT_char_length` is too small (or `PT_Re` too large).  These are the two
  sides of the critical-damping optimum of §2.2; adjust `PT_char_length`
  first, `PT_Re` second.
- **Residual stalls at a plateau well above tolerance.** Likely causes, in
  order: (a) large regions at plastic yield — the coupled
  equilibrium-with-yield problem is genuinely harder; consider more steps
  (smaller `Δt` via `dt_fraction`) or slightly relaxed tolerance; (b) an
  unbalanced boundary condition the residual legitimately measures (e.g. a
  traction BC inconsistent with the interior state).  The stagnation exit
  (§2.10) detects case (a) and ends the step at the plateau instead of
  burning `PT_max_iter`; watch the reported levels across steps — a slow
  upward drift tracks the growth of the localized plastic zone.
- **Residual grows.** Should not happen with default parameters; if it does,
  reduce `PT_CFL`. Persistent growth indicates a configuration the stability
  analysis does not cover — please report it.
- **Adaptive Re (retune) behaves unexpectedly.**  The first retune call at
  iter `PT_retune_interval` can observe a fast-mode transient — a large initial
  residual drop — which gives a high λ_min and pushes Re toward its upper
  bound (2 × `PT_Re`), *slowing* that particular PT loop.  This is normal: the
  upper bound prevents instability, and later retune calls (once slow modes
  dominate) will return Re to the lower range.  If the first retune is
  consistently mis-estimating in a problematic way, try doubling
  `PT_retune_interval` to skip the transient window.  To disable adaptive Re
  entirely while keeping periodic `update_pt_params`, set
  `PT_retune_interval = 0` (disables both) or — if you want periodic
  param updates only — you can set `PT_retune_interval` to a value but note
  that both features share the same interval guard.
- **Steps that hit `PT_max_iter` leave the last (unconverged) velocity
  iterate in `vel`**, which then advects the mesh. A run where many steps do
  this will produce distorted meshes and remeshing failures — treat
  `PT_max_iter` as a safety net, not an operating mode.  With the
  stagnation exit enabled (the default) this should essentially never
  happen: plateaued steps end at the plateau with a near-equilibrated
  velocity.  If you disabled the exit (`PT_stagnation_window = 0`) and see
  `PT_max_iter` hits, re-enable it or lower `dt_fraction`.

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
| `update_pt_params()` | geometry.cxx | per-element `G̃Δτ`, per-node `dτ_ρ`, volume-weighted μ_ve / G means; reads `PT_Re_adaptive` when set |
| `rayleigh_update_Re()` | geometry.cxx | Rayleigh-quotient λ_min estimate → adaptive Re (§2.11) |
| `update_stress_PT()` | rheology.cxx | Sect. 2.4 (Eq. 34) relaxation + plastic projection |
| `update_velocity_PT()` | fields.cxx | accelerated velocity update |
| `copy_stress_PT()` | fields.cxx | save/restore `τ_old` |
| `calculate_characteristic_force()` | fields.cxx | gross-force scale for convergence test |
| `calculate_residual_force()` | fields.cxx | RMS out-of-balance force |
| `update_force()` (residual section) | fields.cxx | BC-inclusive, DOF-masked residual assembly |
| `compute_dt_PT()` | geometry.cxx | physical `Δt` without the DR stability limit (§2.9) |
| `pt_stagnation_window()` | dynearthsol.cxx | auto check window, `2L/(CFL·h_mean)` (§2.10) |
| PT loops + predictor–corrector | dynearthsol.cxx | step orchestration, retune guard, snapshot management |

## 4. References

- Räss, L., Utkin, I., Duretz, T., Omlin, S., and Podladchikov, Y. Y. (2022).
  Assessing the robustness and scalability of the accelerated pseudo-transient
  method. *Geoscientific Model Development*, 15, 5757–5786.
  https://doi.org/10.5194/gmd-15-5757-2022
- Cundall, P. A. (1987). Distinct element models of rock and soil structure.
  In *Analytical and Computational Methods in Engineering Rock Mechanics* —
  origin of the dynamic-relaxation / unbalanced-force-ratio practice used by
  FLAC and inherited by DynEarthSol's convergence metric.
- Duretz, T., de Montserrat, A., Sevilla, R., Räss, L., Utkin, I., and Spang, A.
  (2026). Automatic tuning of iterative pseudo-transient solvers for modeling
  the deformation of heterogeneous media. *Geoscientific Model Development*, 19,
  5343–5362. https://doi.org/10.5194/gmd-19-5343-2026
