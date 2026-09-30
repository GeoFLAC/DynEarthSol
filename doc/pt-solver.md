# Physical-step pseudo-transient mechanics

This refactor defines the mechanical trial/commit boundary for future coupling.
Validation coverage and convergence limits are stated below; it does not replace
the existing constitutive integration algorithms.

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
repeated increments. Non-PT callers retain the existing defaults. Material trial
tests cover these pressure and plane-strain paths. This does not qualify a new
coupled hydraulic formulation.

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
[Räss et al. (2022)](https://gmd.copernicus.org/articles/15/5757/2022/); convergence
remains dependent on the material response, mesh, loading and requested tolerance.

Physical mesh motion, stress rotation, transport and remeshing remain outside the
trial loop. Equilibrium is established on the reference geometry; it does not
claim equilibrium after moving that geometry. No live baseline crosses a remesh.

For a moving-mesh physical PT step, height-dependent pressure boundary laws
(Winkler foundation, water and side support) use the candidate facet height
`z_reference + physical_dt * mean(candidate_vertical_velocity)`. Facet normals,
areas, internal-force gradients and body-force volumes remain at the reference
configuration. This is an incremental pressure-load correction, not a complete
finite-deformation follower-load formulation. Initial equilibration, fixed-grid
runs and non-PT calls retain reference-height loads. The candidate height is
computed afresh; no trial writes coordinates or advances surface diffusion.

Surface processes run in `update_mesh` after accepted mechanics. Their `dh`
buffer is cleared at each physical call; `dhacc` intentionally records physical
surface changes for marker correction. Neither is advanced by mechanical trials.
This boundary-load correction removes the frozen-foundation force-compatibility
obstruction, but does not guarantee convergence with the default PT iteration
limit. Core-complex long-time morphology remains unqualified until the remaining
slow-convergence/residual-floor behavior is resolved.

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
Dry checkpoints do not require this hydraulic field. Existing barycentric remapping
transfers the pending increment with the other nodal hydraulic fields.

PT restart preserves the saved velocity, next-step dt and velocity scale. It does
not reapply the legacy velocity boundary map or select a new dt merely because
execution resumed. The next PT solve imposes its constraints at the next physical
time; mesh and slow-update stages retain ownership of subsequent dt selection.
Exact continuation assumes unchanged controls. Changed timestep controls take
effect at those existing selection stages, rather than in the restart loader.
A pending remesh can therefore still select dt before the next step, just as in
uninterrupted execution. Non-PT restart retains its existing boundary/dt behavior.

The hydraulic transport equation and its existing in-plane mean-stress coupling
remain unchanged; this refactor does not implement a new poroelastic formulation.
At corners the ordered boundary dispatch retains its precedence rules; the affine
fixed-point construction is not a simultaneous enforcement of every conflicting
boundary prescription. The focused constraint tests cover the stated projection and precedence semantics;
arbitrary contradictory boundary specifications are not supported.

## Validation and limits

The [focused tests](../tests/residual-assembly-test/README.md) exercise actual force,
constraint, material, residual and failed-solve paths in 2D/3D on ASan/UBSan,
OpenMP and OpenACC. Separate [functional tests](../tests/functional/pt_restart.md)
cover initial checkpoints, pending remesh stages and regular output scheduling.
The existing [EP/RSF benchmark](../benchmarks/simple_shear_rsf/README.md) has an
opt-in PT mode. Maxwell relaxation has also been checked against its discrete
update and continuum time-convergence reference. Those checks do not establish
convergence for every material/loading combination.

- A solve may reach the iteration limit or stagnate. Retain its failed status;
  it must not commit a discarded candidate or silently relax its tolerance.
- A relative-only target can fall below floating-point cancellation accuracy near
  equilibrium. No automatic residual floor or time-step retry is introduced.
- Weak volumetric modes in some 2D viscous cases require more than the default
  iteration cap. Difficult finite-increment 3D plastic cases can also fail.
  The existing incremental Mohr–Coulomb return is preserved; nonsmooth behavior
  is not a reason to change that algorithm within this refactor.
- Checkpoint tests cover the regular output schedule with unchanged output controls.
  Earthquake-event history and averaged-output accumulators are outside that scope.
- Dry GPU mechanical correctness and dry remesh/restart have scoped coverage;
  full hydraulic GPU qualification and realistic-size GPU performance are not
  claimed. Small-mesh launch overhead remains a known cost.

No experimental secant extrapolation, adaptive coefficient tuning or new plastic
corner-return algorithm is part of this implementation. Such changes require
separate numerical contracts and validation.
