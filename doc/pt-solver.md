# Physical-step pseudo-transient mechanics

Implementation draft: final numerical and backend qualification is pending.

With `control.has_PT=yes`, mechanics starts from the physical-step baseline,
solves on fixed reference geometry, and leaves one accepted constitutive state.
The no-PT path keeps its existing constitutive and velocity-update ordering.

`MechanicalState` captures stress, strain, out-of-plane stress, plastic history,
per-step plastic increment, pressure correction, viscosity and RSF history.
`evaluate_mechanical_trial` restores that baseline in place before computing the
candidate; material-property aliases remain valid. The physical dt and current
pressure inputs are fixed during each mechanical solve. NMD applies to every
candidate before force evaluation. Trial Maxwell volume strain uses the corrected
strain-rate trace times physical dt on the fixed reference mesh.

The ordinary `update_stress` return maps remain the constitutive authority. For PT,
Maxwell normal stresses receive the pending isotropic pressure increment once.
Plane strain includes yy in mean/deviatoric stress, the EVP branch comparison and
NMD; the selected EVP candidate supplies yy as well as the in-plane components.
Pure viscous replacement stress uses the current pore-pressure level rather than
repeated increments. Non-PT callers retain the existing defaults. These paths still
require the final material and coupling tests before numerical qualification.

The residual is the RMS assembled undamped force projected onto admissible velocity
perturbations, normalized by the free rank. Gravity, Neumann traction and foundation
loads remain in the force. Only constrained reactions and artificial damping are
excluded. All-constrained models have zero free residual. Convergence requires

```
R <= PT_absolute_tolerance + PT_relative_tolerance * R_initial
```

The relative scale is the first candidate's imbalance, not the gross prestress.
Absolute tolerance has nodal-force units: N/m in 2D, N in 3D. Stagnation is failure,
not convergence. Results distinguish `converged`, `max_iterations`, `stagnated` and
`nonfinite`; a failed solve restores mechanical history and stops before output,
transport or mesh movement can consume the candidate.

Boundary constraints use the same ordered dispatch as `apply_vbcs`, evaluated in
affine and homogeneous modes. Corrected PT normalizes oblique type-11 and edge
operations, then derives an orthogonal free-space projector from the fixed-point
constraints. The default legacy boundary calls retain their previous operations.
The final PT velocity is already constrained and is not passed through the legacy
map a second time. Incompatible affine fixed-point constraints are rejected.

Localized dual-time stress/velocity relaxation uses mesh heights, effective
viscosity, `PT_CFL` and `PT_Re`. Lagged stress is solver scratch. Acceptance always
uses the full constitutive candidate, so a small residual of lagged stress cannot
commit an inconsistent material state. No plastic-yield heuristic is introduced.
The factor construction follows the effective-viscosity/dual-time approach in
[Räss et al. (2022)](https://gmd.copernicus.org/articles/15/5757/2022/); transferring
it to this discretization still requires the final convergence tests.

Physical mesh motion, stress rotation, transport and remeshing remain outside the
trial loop. Equilibrium is established on the reference geometry; it does not
claim equilibrium after moving that geometry. No live baseline crosses a remesh.

Initial body-force equilibration uses the instantaneous elastic/plastic skeleton
with homogeneous supports and no transport, viscous relaxation or physical aging.
RSF friction and state are frozen at the initial physical values, independent of
numerical correction velocities. Plastic history is committed only for the accepted
initial candidate. Maxwell and EVP initialization use their instantaneous skeletons.
Pure viscous material has no such skeleton and is explicitly rejected for this
initialization contract. Duration-based isostasy retains its physical-time process.
Initial equilibrium precedes frame-zero output and monitor initialization.
Completion state is checkpointed so a restart does not repeat initial adjustment.

Hydraulic checkpoints now store the pending signed pressure increment consumed
by the next mechanical step. Old hydraulic checkpoints lacking it are rejected:
the missing increment cannot be inferred exactly from a pressure level alone.
Dry restart behavior is unchanged. Existing barycentric remapping transfers the
pending increment with the other nodal hydraulic fields.

The hydraulic transport equation and its existing in-plane mean-stress coupling
remain unchanged; this refactor does not implement a new poroelastic formulation.
At corners the ordered boundary dispatch retains its precedence rules; the affine
fixed-point construction is not a simultaneous enforcement of every conflicting
boundary prescription. These cases require explicit final validation.
