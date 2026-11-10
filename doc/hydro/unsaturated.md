# Unsaturated flow: explicit formulation and verification ladder

DES implements explicit unsaturated transport, element-aware Bishop stress
and saturation-dependent gravity. The Python reference and prototype are checked
(`benchmarks/poromechanics/unsaturated_reference.py`,
`tests/functional/hydromechanics/check_unsaturated_reference.py`, 19 checks
including the U8 Bishop column, about 5 s). The retention header
`unsaturated.hpp` is checked against that reference
(`check_unsaturated_header.py`, nine checks, g++ with warnings as errors). Its
host and parallel evaluation are compared in `tests/unsaturated-test/`
(serial, OpenMP and ASan pass; OpenACC runs on a separate GPU machine).
The production solver includes retention parameters and laws, unsaturated
transport and water-table initialization (U5-U7 below). Mechanics uses an element-aware Bishop increment in the
effective stress, the same chi_e in the fixed-stress transport term and
gravity on the partly saturated bulk density (U8 below). Retention
laws are not yet combined with evolving porosity/permeability, thermal
pressurization or matrix exchange; those inputs are rejected.

The DES hydraulic update stays explicit. This formulation adds no implicit
time integration and no global Picard or Newton solve. The only nonlinear
operation is a scalar root at each node, which inverts a monotone
pointwise function.

## Physics and conventions

Richards flow with the gas phase at constant (zero gauge) pressure; no
hysteresis. DES conventions are kept: z points up, the pore pressure p is
positive in compression, suction is -p, and

    q = -(k*kr(p)/mu) * (grad p + rho_f*g*e_z).

Saturation laws, with s = max(-p, 0), Se the effective saturation and
S = Sr + (1 - Sr)*Se:

| Model | Se(p) | kr(Se) | Saturated for |
| --- | --- | --- | --- |
| van Genuchten-Mualem | (1 + (alpha*s)^n)^-m, m = 1 - 1/n | Se^0.5 (1 - (1 - Se^(1/m))^m)^2 | p >= 0 |
| Brooks-Corey-Burdine | (alpha*s)^-n for alpha*s > 1 | Se^(3 + 2/n) | p >= -1/alpha |
| Gardner exponential | exp(-alpha*s) | Se | p >= 0 |

`alpha` is in 1/Pa (1/alpha is the air-entry pressure for Brooks-Corey).
For Brooks-Corey, n is the pore-size index lambda. The Gardner law is kept
because it has an exact transient solution (Srivastava and Yeh, 1991), which
gives the transport operator a true analytical test.

Bishop effective stress (tension positive) uses chi = Se (Lu et al., 2010):

    sigma' = sigma + alpha_B * pi,   pi = Se(p) * p   (pi = p when saturated).

## Conserved quantity

The nodal fluid volume per bulk volume is

    F(p) = phi*S(p) + s*p,

where `s` is the linear storage selected by `hydraulic_storage(param, mat, e)`.
PT and reduced quasistatic flow use fixed-stress storage
`1/M + alpha_B^2 * compliance`; coupled explicit/GVS uses `1/M` and the
projected strain source before stress. This selection does not depend on
retention law. F is strictly increasing (`dF/dp >= s > 0`), so it has a
unique inverse at every node. When saturated, `F = phi + s*p`, exactly the
linear update for the same mechanical formulation. The
[mechanical formulation table](dynamics.md) records the ordering and
checkpoint contract. The one-dimensional transport references prescribe
their chosen linear storage explicitly.

The lumped nodal content follows the existing `hmass` assembly. For a node i
with patch elements e,

    M_i(p) = hmass_i * p + sum_e (V_e/N) * phi_e * S_e(p),

each element uses its own retention law at the nodal pressure. This keeps
material interfaces exact without nodal material properties.

## Discrete update (mapping onto `update_pore_pressure`)

The transport and coupling update the element loop as follows:

1. Darcy mobility: multiply `kv_pressure`/`kv_gravity` (now the tensor flux)
   by `kr_e` = the arithmetic mean of the nodal `kr_e(p_j)`.
2. Bishop coupling: multiply the poroelastic term `pe` by
   `chi_e` = the mean of the nodal `Se_e(p_j)`. The same chi enters the
   effective stress (`rheology.cxx`), so the coupling stays symmetric.
   The mixture weight in `update_force` uses the element saturation
   `S_e`, the mean of the nodal `S_e(p_j)`.

At a free node, the existing loop computes `increment = tdot*dt/hmass`
(minus the change in p). The unsaturated update replaces this with a content
increment `dm = -tdot*dt` and solves

    sum_e (V_e/N)*phi_e*(S_e(p') - S_e(p)) + hmass*(p' - p) = dm

for p' in increment form:

- Saturated fast path. If p and the trial `p + dm/hmass` are both at or above
  every patch element's `p_sat`, then p' is that trial. This is bit-identical
  to the current code. `retention_model = 0` keeps the current code path
  untouched.
- Otherwise, a safeguarded Newton iteration. The function is monotone, and
  the trial `p + dm/hmass` and the old p bracket the root. The first guess is
  the closed-form inverse of Se with storage neglected. Stop when the bracket
  collapses to a few ulps; the prototype's residual is <= 1e-14 in content
  units. Writing the residual as an increment, with no absolute content,
  avoids the cancellation of `phi + s_b*p` when s_b*p is small. The C++
  version stops when |g| <= 4*eps*(|dm| + sum w): the residual cannot be
  resolved below that, and an ulp-width bracket test never terminates there.
  Measured over 1e5 random updates per law, the mean is about 3 iterations;
  the maximum is 3 for van Genuchten and Brooks-Corey and 10 for Gardner.

`dppressure` stores raw pressure increments. Each element forms its own
Bishop increment in `rheology.cxx`, as described in
[element-aware coupling](#element-aware-bishop-coupling).

The fluid balance (`var.fluid_balance`) counts `stored` as the change in
M_i(p) rather than `hmass*Delta p`. Boundary reactions keep their current
definition.

No new history is introduced: p is the only state. S, kr and chi are
derived, so the existing checkpoints, nodal remapping and restart carry over
unchanged. Remapping p does not conserve fluid content exactly across
saturation fronts. That belongs to the separate conservative-remap job
(queue item 11).

## Explicit diffusion bound and nonlinear limits

At a free node, write the update as
`p' = p + (dt/(V_i*c))*sum_j T_ij*lambda*(p_j - p_i)`
for source-free diffusion with nonnegative, fixed transmissibilities.
Here c is the chord slope
(F(p') - F(p))/(p' - p). F is monotone with slope >= s_b, so c >= s_b, and
kr <= 1 gives lambda <= k/mu. The step is therefore a convex combination
(monotone, with no overshoot) whenever

    dt <= V_i*s_b / sum_j T_ij*(k/mu)

This conditional estimate explains retaining the saturated diffusion bound.
Production uses the assembled tensor absolute-row bound in
`hydraulics.cxx::pressure_step_limit`, not only the scalar geometry estimate
`0.5*h^2/D_sat`. The one-dimensional prototype's Courant check supports its
own column calculation; it is not a proof for arbitrary DES patches.
Pressure-dependent mobility, gravity, sources, interfaces and obtuse elements
need separate nonlinear stability/monotonicity qualification. Stable separate
transport and mechanical steps also do not establish coupled stability.

A sharper local bound, based on the current kr and dF/dp with an a-posteriori
chord check and step halving, is implemented in the prototype (`mode='local'`).
It only pays off when no node in the mesh is near saturation, because van
Genuchten and Brooks-Corey give dF/dp -> s_b at the water table, and that
node then sets dt anyway. DES uses the assembled row bound for subcycling;
the prototype's local optimization is not used by the production selector.

P1 elements with obtuse angles do not give an M-matrix. That caveat is the
same as in the saturated update and is not new here.

## Gravity and density

When gravity is on, the bulk density becomes
`rho_s*(1 - phi) + rho_f*phi*S` (the element mean of nodal S) in the
body-force assembly (`fields.cxx`, via `hydraulic_bulk_density`). The
inertial/scaled mass in `compute_mass` may keep the saturated value. In
quasi-static runs it only sets the relaxation speed, and rebuilding it every
step is not justified.

## Parameters

| Name | Units | Meaning |
| --- | --- | --- |
| `mat.retention_model` | - | 0 saturated (default, current code), 1 van Genuchten-Mualem, 2 Brooks-Corey-Burdine, 3 Gardner |
| `mat.retention_alpha` | 1/Pa | per material |
| `mat.retention_n` | - | VG n > 1, or the BC pore-size index > 0; unused for Gardner |
| `mat.residual_saturation` | - | per material, 0 <= Sr < 1 |
| `ic.water_table_z` | m | optional hydrostatic water table: p = -rho_f*g*(z - z_wt), negative above it |

## Reading saturation and mobility

With hydraulic diffusion and an active retention law, each saved frame
includes four element fields:

| Field | Units | Meaning |
| --- | --- | --- |
| `saturation` | 1 | Water-filled pore fraction S, including residual saturation |
| `effective saturation` | 1 | Se, the Bishop coefficient used by the mechanical coupling |
| `relative permeability` | 1 | kr multiplying the intrinsic hydraulic permeability |
| `bulk density` | kg/m³ | Partly saturated mixture density used for gravity |

Each field averages the material laws using the element's marker fractions,
then averages over the element's nodes. In particular, averaging Se(p) or
kr(p) differs from evaluating the law at the mean pressure. Bulk density
uses the same element-mean porosity and material densities as the body-force
assembly. The existing `density` field is the solid material density.
Inertial mass and wave-speed estimates still use the saturated mixture
density; this output addition does not change that choice.

These four fields describe the saved pressure state, including when stress
averaging is enabled. They are derived outputs, with no additional
checkpoint history. They are absent when all retention models are zero.
Both the Python reader and `2vtk.py` support them; native VTK-HDF writes
them as cell data even when node and element counts coincide. The global
five-entry `fluid balance` remains separate from cell and point data.

For example, after a run named `column`:

```python
from Dynearthsol import Dynearthsol

model = Dynearthsol('column')
frame = model.frames[-1]
saturation = model.read_field(frame, 'saturation')
relative_permeability = model.read_field(frame, 'relative permeability')
```

The interface check uses six-node, six-element meshes in both dimensions,
with van Genuchten, Brooks-Corey and Gardner materials. Wet, partly
saturated, dry and nonuniform-pressure cases compare against the independent
Python laws and saved marker fractions. A disabled-law case checks that the
fields are omitted. Run it for either a binary or HDF5 executable:

```bash
OUT=benchmarks/poromechanics/results/runs/retention-check
OMP_NUM_THREADS=1 python3 tests/functional/hydromechanics/check_retention_output.py ./dynearthsol2d "$OUT/2d" --ndims 2
OMP_NUM_THREADS=1 python3 tests/functional/hydromechanics/check_retention_output.py ./dynearthsol3d "$OUT/3d" --ndims 3
```

The qualification record is in
`benchmarks/poromechanics/results/runs/retention-output-20261108/`.
These output checks do not close the U9 drainage or application requirements.

## Implementation contracts

The five parameters above follow the input registry and
`examples/defaults.cfg`. `retention_alpha` is in inverse pressure.
Material values use the existing
`nmat` length/broadcast validation, finite values, and model-specific ranges.
The model is a per-material integer list (default `[0]`), supporting different
laws at interfaces. Omitting the water-table option preserves the existing
initialization; prescribed type-2 pressure takes precedence at boundary nodes.

The transport and constitutive updates keep these contracts:

1. `update_pore_pressure` resets the pending increment once, before all
   hydraulic substeps. The next physical mechanical update consumes that
   increment; PT constitutive trials must only read it. Preserve this
   checkpoint/remapping contract. Never reset it inside `pressure_step`.
2. In the model-0/all-saturated fast path, preserve
   `increment = tdot*dt/hmass`, `p -= increment`, `dp += increment`, including
   the existing type-2 boundary arithmetic. Recomputing `dp` as `p_old-p_new`
   changes rounding: at `p_old=1e8`, an increment of `1e-8` gives a rounded
   pressure difference of about `1.49e-8`. Algebraic equality alone does not
   establish U5 bit-exactness for the pending increment or mechanics.
3. A nodal pressure is shared, but `Se_e(p)` can differ between the patch's
   materials. A single nodal Bishop increment cannot represent every
   element's law. Specify an element-aware increment or an explicit,
   validated common-law restriction; do not silently average away this
   distinction. Test a shared node with different retention laws.
4. `update_porosity` currently uses the negative raw pressure increment for
   the grain-compressibility term. Changing its input to Bishop pressure
   changes that constitutive law. Derive and verify that change separately,
   or preserve access to the raw pressure increment. The arithmetic
   telescoping argument does not establish constitutive consistency.
5. Retain the current assembled tensor row-sum stability bound, including
   matrix exchange. The old scalar `0.5*h*h/D` estimate is not the complete
   subcycling criterion. Qualify changing mobility under gravity and
   heterogeneous/obtuse patches before claiming nonlinear monotonicity;
   the constant-mobility M-matrix argument alone is insufficient.

Scalar reference tests check the content law separately from connected DES
runs. The U5-U8 checks below cover the saturated path, transport and Bishop
coupling; the [coupled drainage example](drainage.md) supplies the bounded
U9 mechanical comparison.

Infiltration uses boundary type `hbc = 3` (outward Darcy flux in m/s), with
a negative value for rainfall. A prescribed flux continues through
saturation; no surface-ponding reservoir or automatic seepage-face switch
is implemented. Rainfall applications must qualify those surface conditions
separately. The proposed ponding-stop diagnostic is not present in this stage.

## Element-aware Bishop coupling

`dppressure` keeps its meaning,
the accumulated raw -Delta p over the hydraulic substeps. The transport fast
path, the type-2 boundary arithmetic, checkpointing, remapping and the
porosity grain term (requirement 4) are therefore unchanged. Bishop pressure
is not stored at the nodes.

Each element computes its own Bishop increment in `rheology.cxx`, where
`dpp` is formed today (requirement 3). With nodal p_j and the reconstructed
old pressure p_j + dppressure_j,

    dpp_e = -alpha_e * mean_j [pi_e(p_j) - pi_e(p_j + dppressure_j)],

where pi_e is that element's law. When the element's law is saturated at
both pressures of every node, the current expression
`alpha_e * mean_j dppressure_j` is used unchanged, so model 0 and
all-saturated meshes stay bit-exact (requirement 2). The ppressure gather it
needs runs only for elements with `retention_model > 0`.

Reconstruction error: p_j + dppressure_j
differs from the stored old pressure by the rounding accumulated over the
step's substeps, a few ulps of the pressure level |p_j|. That error is not
small relative to the increment: when the Bishop increment is only a few
ulps of |p| itself, its relative error can be O(1). The absolute error stays
at the ulp level of p and does not accumulate across mechanical steps,
because each step restarts from the stored p. No extra history is added only
to remove this ordinary rounding. Coupling checks cover several substeps,
mixed-law shared nodes, saturation
crossings in both directions; and repeated PT trials, which must leave p and
dppressure untouched. Those tolerance-based unsaturated checks stay separate
from the bit-exact saturated check.

The fixed-stress transport term uses chi_e, the mean of the element's own
nodal Se_e(p_j). The content update sums each patch element's own law
(`pore_content`). No quantity is averaged across different laws.

Stability (requirement 5): the assembled tensor and matrix row-sum bound is
kept as it is. The chord-slope argument shows that unsaturated rows add no
restriction beyond the saturated row sums (chord >= s_b, kr <= 1). This
holds for constant mobility within a step. The monotonicity claim with
gravity, mobility that changes during the step, and heterogeneous or obtuse
patches is left to U6 and U7, which are measured, not assumed.

## Approximations, stated

- Elastic storage acts as s_b*p in the unsaturated zone instead of S*s_b*p.
  The content difference is at most s_b*|p|*(1 - S): 1e-4 for s_b = 1e-9 /Pa
  at 100 kPa suction, but 1e-2 for a soft soil with s_b = 1e-7 /Pa. Compare
  this with phi*Delta S before using strongly compressible, very dry
  materials. The same choice is what makes F globally monotone.
- Bishop coupling uses chi = Se in both effective stress and the
  poroelastic transport term. With Sr = 0 this equals chi = S.
- With pressure feedback, rheologies that include a plastic branch use
  the accepted projected strain directly: alpha*chi*deps, with the
  stabilizing raw-pressure term -L*D. They do not apply the elastic
  mean-stress/Bishop correction a second time. The
  [yielding-cell reference](plastic-content.md) checks a homogeneous
  Gardner wetting path against independent scalar content inversion.
- Elastic fixed-stress correction: for the implemented
  fluid-rate term alpha*chi*deps_v, the source reads a mechanical stress change
  K_b*deps - alpha*E that already contains the Bishop increment E of the
  last physical step, while the stabilising storage stays S_B = S0 + L,
  L = alpha^2*compliance. Without a correction PT converges to a different
  transient, capacity S0 + L + phi*dS/dp - L*chi*Pi'. At the start of each
  hydraulic update, before the pending increment is cleared,
  `pressure_increments` freezes per element D = mean(p - p_old) and
  E = mean(Pi(p) - Pi(p_old)), Pi = Se*p, with the element's own law, the
  increments that mechanics has just consumed. The source gains
  L*(chi*E - D)/dt with the substep's chi, which leaves
  S0*dp + phi*dS = -alpha*chi*deps + L*(D - dp): a lag that vanishes as the
  step shrinks, with no change of the continuum or of the step bound.
  If both pressures are saturated, E is set equal to raw D exactly; the
  correction is skipped when chi*E - D is exactly zero. D must remain
  available because a later hydraulic substep can enter suction. The
  sealed wet-transition regression (`check_bishop_history.py`) starts
  wet with nonzero pending pressure history and crosses zero within a
  physical step. Refining hydraulic substeps from 160 to 1600 (2D) and
  240 to 2400 (3D) reduces the largest interval error against an independent
  scalar ODE from 7.02e-6 to 7.02e-7 Pa and 6.88e-6 to 6.87e-7 Pa.
  Restart immediately before that transition preserves 4,962 physical
  save/checkpoint payloads exactly across the four runs. Zeroing D/E at
  the start instead leaves errors of 0.00170/0.00646 Pa after refinement.
  A clipped storage,
  max(chi*Pi', 0), was tried and rejected: where Pi' < 0 (Gardner beyond
  alpha*suction = 1) it converges to a different continuum.
  The sealed Bishop cell (`check_bishop_cell.py`) tests the implemented rate
  equation. Uncorrected, PT reached 7.8 Pa against an exact 37.0 Pa at every
  step. Corrected, the moderate state (chi*Pi' > 0) falls 39.4 -> 22.5 ->
  11.4 -> 5.73 -> 2.86 -> 1.43% at 10..320 steps, and the strong-suction
  state (Pi' < 0) 20.1 -> 10.3 -> 5.13 -> 2.57 -> 1.28 -> 0.64%, both first
  order and extrapolating to the exact 37.037 and 3.7037 Pa. The explicit
  strain form also converges to that rate equation. A restart with the same
  settings reproduces the fresh run exactly. The versioned timestep-control
  contract reselects changed controls before the first resumed interval;
  broader unsaturated variable-step accuracy remains a separate check.
- Legacy checkpoints do not identify whether their mechanical history used
  chi=1 or Bishop stress. Changing that constitutive law can change both
  the interpretation of a pending increment and the historical stress/strain
  state; a drained or equilibrated state alone does not establish compatibility.
  A constitutive-version guard or an explicit migration remains open.
  Continue such unfinished runs with their recorded source and executable.
  The current-model restart checks above do not prove legacy-model migration.
- Gravity uses the partly saturated weight. Lumped mass and wave-speed
  estimates keep the saturated density, a slightly larger mass. The
  lithostatic initial stress still uses the saturated weight, so above the
  water table the initial PT equilibrium (`has_body_force_adjustment`)
  removes a real imbalance of phi*(1 - S)*rho_f*g per volume and can need
  more iterations than a saturated start: the 2.5 cm 3D U6 column needs
  more than 10,000.
- A mixed-material element weights each material's saturation by its marker
  fraction times the element's marker-mean porosity, the same element-mean
  convention as `hmass`, not by each material's own porosity.
- Element kr is the arithmetic mean of nodal kr. Infiltration into very dry
  soil may need upstream weighting. That would be a separate, tested
  numerical change.
- Explicit mobility evaluates kr at the old pressure and never evaluates
  its derivative. This keeps the flux update explicit; it does not by
  itself establish nonlinear stability near saturation.

## Verification ladder

| Level | Check | Where | Status |
| --- | --- | --- | --- |
| U0 | Constitutive laws: saturated limits, monotonicity, dSe/dp against FD, closed forms, inverses | Python | Done |
| U1 | Content inversion: saturated branch identical to p + dm/s, residual <= 1e-14, wetting through saturation and back | Python | Done |
| U2 | Srivastava-Yeh series against an independent matrix-exponential solution (< 3e-5) | Python | Done |
| U3 | Explicit column against Srivastava-Yeh: second-order convergence, balance <= 1e-15, Courant <= 1 | Python | Done (errors 4.2e-3, 1.1e-3, 2.9e-4 at 21/41/81 nodes) |
| U4 | Explicit column steady VG/BC profiles against the exact ODE profile, under the existing saturated dt | Python | Done (VG second order; BC about first order at the air-entry kink) |
| U5 | DES `retention_model = 0` dry cores and hydraulic tests bit-exact with the parent; models 1-3 with an all-saturated mesh bit-exact with model 0 | DES 2D/3D | Done: make cmp and all 28 hydraulic frames BIT-EXACT; Brooks-Corey (1 MPa air entry) BIT-EXACT with model 0 |
| U6 | DES thin column against Srivastava-Yeh (Gardner, hbc = 3 top, hbc = 2 bottom) | DES 2D/3D (`check_unsaturated_column.py`) | Done: see "Transport verification" |
| U7 | DES steady VG/BC infiltration profiles; exact restart and remeshed restart | DES 2D/3D (`check_unsaturated_column.py`) | Done: see "Transport verification" |
| U8 | Bishop column: lower the water table from 9 m to 5 m in an elastic, laterally confined 10 m column; settlement u(z) from the Bishop term and the drained overburden | `check_bishop_column.py`, DES 2D/3D, PT; van Genuchten, Brooks-Corey, Gardner and a two-law column (shared interface nodes) | Done for drained states: 2D errors at 1/0.5/0.25 m are 0.255/0.055/0.012% (VG), 0.584/0.127/0.024% (BC), 0.716/0.155/0.040% (Gardner), order 2.1-2.3; two-law 1.02/0.38/0.16%, order 1.3 (saturation jumps at the interface). The pre-Bishop build errs by 26-99%. Transient Bishop coupling is checked separately by the sealed cell and U9 |
| U9 | Coupled drainage (Liakopoulos-type) and the rainfall/slope applications in benchmarks.md | DES | Bounded 2D PT and ordinary-explicit drainage is compared with an independent Richards–Bishop reference; see [drainage](drainage.md). The sealed Bishop cell has an exact ODE reference. Broader GVS/3D regimes and reproduction of the Liakopoulos experiment remain open |

## Transport verification

`tests/functional/hydromechanics/check_unsaturated_column.py EXE OUT` runs a
laterally confined column, z in [-1, 0], with `ic.water_table_z = -1`, top
`hbc = 3` and base `hbc = 2` (p = 0). The soft elastic skeleton (K = G = 10
MPa) gives s_b ~ 7.5e-8 /Pa. That keeps the explicit substep count near 1e5
and shifts the solution from storage-free Srivastava-Yeh by only
s_b*|p|/dtheta ~ 0.3%.

Transient (Gardner, 0.1/cm, Ks = 1 cm/h, top flux 0.9 Ks, T = 0.5 at 7200 s);
maximum error of K/Ks:

| Dimension | Resolution | T = 0.13 | T = 0.25 | T = 0.38 | T = 0.50 | Balance residual / step inflow |
| --- | --- | --- | --- | --- | --- | --- |
| 2D | 5 cm (26 nodes) | 0.109 | 0.081 | 0.057 | 0.044 | 5.4e-13 |
| 2D | 2.5 cm (50 nodes) | 0.037 | 0.024 | 0.020 | 0.016 | 9.3e-12 |
| 2D | 1.25 cm (98 nodes, exploratory run, not in the test) | 0.012 | 0.0078 | 0.0056 | 0.0047 | not measured |
| 3D | 5 cm (166 nodes) | 0.136 | 0.090 | 0.054 | 0.036 | 8.2e-11 |
| 3D | 2.5 cm (309 nodes) | 0.025 | 0.015 | 0.0099 | 0.0071 | 1.1e-9 |

That is a ratio of 2.7-3.4 per halving in 2D (approaching second order like
the 1D prototype, U3) and 3.6-6.0 in 3D. Pressure falls to -5.5 kPa, so the
unsaturated path runs. The steps are 60 s in 2D and 30 s in 3D, which the
elastic PT mechanics requires (see the end of this section).

Steady (flux Ks/2, 2e6 s; a compressible pore fluid, K_f = 40 kPa, shortens
the run because storage does not enter the steady state):

| Law | 5 cm | 2.5 cm | Ratio | Profile range |
| --- | --- | --- | --- | --- |
| van Genuchten (5e-4 /Pa, n = 2, Sr = 0.1), 2D | 25.2 Pa | 5.7 Pa | 4.4 | 595 Pa |
| Brooks-Corey (1e-3 /Pa, n = 1.5, Sr = 0.05), 2D | 33.5 Pa | 18.8 Pa | 1.8 | 1113 Pa |
| van Genuchten, 3D (120 s steps) | 9.7 Pa | - | - | 595 Pa |
| Brooks-Corey, 3D (120 s steps) | 38.7 Pa | - | - | 1113 Pa |

On the same node levels, the 1D prototype has an error of 8.2 Pa (van
Genuchten, 5 cm). The excess in DES comes from P1 triangles whose element kr
averages nodes on two levels, and it converges.

Restart: `pt_initial_checkpoint.py` continues a column with its water table
at mid-depth (van Genuchten) exactly, both on a fixed mesh and across a real
remesh.

Fluid balance: each node's stored term is its actual content change. The
remaining residual is round-off in that content. Each update rounds at about
eps times the node's pore volume, with random sign (signed/absolute sum 0.03
over 1e6 updates), because the Newton-corrected point is returned at the
floor. Relative to one step's inflow, this floor grows with nodes and
substeps: 9e-12 (2D) and 1.1e-9 (3D) at 2.5 cm, against the test's 1e-8.
The saturated path keeps its existing ~1e-16 balance.

The reported dry-column failure was reproduced and traced to explicit heat
diffusion exceeding its stability limit, with temperature diverging before
stress and velocity. It is not an established PT timestep restriction.
The integrated core bounds heat with explicit substeps; the dry 600 s column
then retains exactly constant temperature/stress and zero velocity for 100
steps with both PT options in both dimensions. These isothermal infiltration
fixtures now explicitly disable unused heat diffusion. Their original time
resolutions and single-resolution 3D steady scope are retained so this
integration does not silently change their convergence study.

## Mechanics modes every coupled test must cover

How DES reaches mechanical equilibrium differs between applications, and the
column tests above inherit PT with `inertial_scaling = 1e11` from
`pt_hydrostatic_equilibrium.cfg`. That is far from production settings:

| Mode | Used for | Settings |
| --- | --- | --- |
| PT | Quasi-static equilibrium comparison for sudden/slow loading, including Terzaghi, Mandel and U8 | `has_PT = yes`, both PT options, controlled physical dt |
| Ordinary explicit / dynamic relaxation | Slow and sudden loading comparisons, including tectonic cases | `has_PT = no`, `inertial_scaling` 1e4-1e6, dt from DES |
| Explicit with GVS | Slow and sudden loading comparisons; earthquake-cycle transitions | `use_global_velocity_scaling = yes`, `inertial_scaling` 1e4-1e6, adaptive dt and mass |

Sudden traction does not itself require PT. The integrated early Terzaghi
pilot includes non-PT GVS; its timestep, initial velocity scale and coupled
refinement limits are recorded in
[mechanical-method qualification](mechanics.md).
Earthquake cycles additionally need an accurate inertial response. Internal
stress release alone does not establish equivalence to a quasi-static
analytical load step. Measured transport-only comparisons so far:

- In the checked transport-only cases, the 5 cm Srivastava-Yeh
  column gives an error of 0.0431 without PT (inertial_scaling 1e5,
  GVS on or off) and 0.0438 with PT in 2D, and 0.0353 without PT against
  0.0358 with PT in 3D near the target T = 0.5. Each error uses its run's
  actual saved time. The adaptive run steps are governed by hydraulic
  diffusion, and their counts differ between dimensions.
  `check_unsaturated_column.py` keeps the no-PT/GVS case as a guard.
- The earlier divergence at 600 s across scaling 1e4-1e11 was caused by the
  heat step described above; it does not select the mechanics method.

Checks for broader coupled applications:
- Compare instantaneous- and slow-load benchmarks with PT and ordinary
  explicit/GVS at inertial_scaling 1e4, 1e5 and 1e6. Qualify mechanical
  settling versus drainage/loading and refine the coupled timestep.
- Validate cycle cases without PT against a resolved physical-inertia
  reference, including scaling ceilings, damping and mass-change effects.
- With GVS, check that hydraulic subcycling and the unsaturated content
  update stay stable and conservative while dt jumps by orders of magnitude.
- Report which mode each result comes from.

## References

- van Genuchten (1980), Soil Sci. Soc. Am. J. 44, 892-898.
- Mualem (1976), Water Resour. Res. 12(3), 513-522.
- Brooks and Corey (1964), Hydrology Paper 3, Colorado State University.
- Gardner (1958), Soil Sci. 85, 228-232.
- Srivastava and Yeh (1991), Water Resour. Res. 27(5), 753-762.
- Lu, Godt and Wu (2010), Water Resour. Res. 46, W05515.
- Ippisch, Vogel and Bastian (2006), Adv. Water Resour. 29, 1780-1789.
- Celia, Bouloutas and Zarba (1990), Water Resour. Res. 26(7), 1483-1496
  (mixed-form mass conservation).
