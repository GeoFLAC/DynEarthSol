---
title: Plastic shear heating
---

# Plastic work, temperature and pore pressure

`control.has_shear_heating=true` converts accepted nonelastic work to heat.
The default is false. This integrates Eh Tan's
[existing PR #88](https://github.com/GeoFLAC/DynEarthSol/pull/88) and its
corrected follow-ups into the current poromechanics owners. The priority
is irreversible plastic deformation during earthquakes, including EP-RSF.
It does not introduce a second constitutive or PT solver.

The [worked example](../../benchmarks/poromechanics/results/review/plastic-heat.html)
contains the setup drawing, input, aging-law solution, plastic-work integral,
physical temperature and pressure, and portable tables. It prescribes
spin-free pure shear in a sealed cell; it is a material-level validation.

## What becomes heat

With tension-positive stress and compression-positive pore pressure,
the effective stress is $\sigma'=\sigma+\alpha p_B I$, where $p_B=p$ for
saturated elements and the element's Bishop pressure otherwise.
The accepted plastic work density is

$$
\Delta W^p=\bar\sigma':\Delta\epsilon^p,\qquad
\Delta\epsilon^p=\Delta\epsilon-\mathsf{C}^{-1}:\Delta\sigma'_e.
$$

`rheology.cxx` uses the Hooke inverse of the actual branch, including
out-of-plane stress for plane strain. The reversible thermal correction
precedes the baseline snapshot; the Biot/Bishop stress increment is
removed from the inferred elastic increment. Pressure-only and
thermal-only elastic changes generate no heat. Normal contributions are
summed before clamping the total to zero. Existing viscous and Maxwell
branches retain their own deviatoric dissipation; the selected example
uses plasticity.

All computed nonelastic work is converted: the heat fraction is one.
Stored hardening, fracture energy and a variable heat fraction are not
separate reservoirs in this implementation.

## Physical steps, PT and changed timesteps

Constitutive trials overwrite candidate power in W/m³. PT captures and
restores that field with its existing mechanical state, including failed
trials. Initial PT equilibration generates no heat. After the mechanical
step is accepted, the existing physical caller adds

$$
E_{e,\mathrm{pending}}\mathrel{+}=V_e\,\Delta t_n\,q_{e,n}.
$$

The next thermal update deposits this stored energy, split across its
thermal substeps, then clears it. It never multiplies old power by a new
timestep. A physical temperature frame therefore has a one-step heat
lag. For a closed cell of fixed capacity, the separate ledger is

$$
E_T+E_{\mathrm{pending}}=\sum_{n,e}V_e\Delta W^p_{e,n}.
$$

The fluid-content ledger separately includes PT's pending pressure.
Artificial kinetic energy from inertial scaling or GVS is not a heat
source; changing to physical mass does not deposit that energy.
This uses the existing PT acceptance boundary, without further PT redesign.

## Restart and remapping

Native and HDF5 checkpoints save version 1, accepted power and extensive
pending energy. Missing members, unknown versions and nonfinite/negative
values are rejected. Disabling heating with undeposited energy is rejected.
Older checkpoints without heating start with zero energy. A legacy
nonzero power-only checkpoint cannot identify the generating interval;
continue it using its matching solver rather than inventing pending heat.

The existing remap stencils transfer pending energy with weights normalized
per old element; omitted old elements transfer to their nearest new owner.
This conserves the global pending-energy total. The output field
`pending shear heat` is extensive: J per element in 3D and J/m in 2D.
It is not an energy density. `shear heating` is power density in W/m³.

## Measured evidence and limits

The selected aging EP-RSF run reaches a reference temperature rise of
1.380832 K and pressure of 0.569950 MPa at 0.2 s. At $\Delta t=0.0002$ s,
maximum temperature and pressure errors, normalized by these final rises,
are 0.124% and 0.161%. Halving the timestep halves both errors. The heat
and split fluid budgets close to round-off.

Additional checks cover elastic zero heat, ordinary explicit/GVS/PT,
pressure and thermal exclusions, failed PT rollback, changed timesteps,
exact restarts, and actual topology changes in both dimensions and output
containers. Both parent tiny cores remain bit-exact with heating disabled;
the assembly suite passes 19,790 assertions. The qualification record binds
the builds and sources and distinguishes the initial matrix from the
final legacy-checkpoint guard.

The prescribed cell holds effective confinement constant: total normal
stress reacts to pressure. A separate [fixed-normal-traction patch](fault-heat.md)
demonstrates thermal weakening with physical mixture inertia and zero
damping, including a manual slow-GVS-to-physical transition. Neither cell
contains a rupture front or a complete earthquake cycle.
General enthalpy conservation with changing fluid storage and remapped
temperature is also outside this check. The role of heating and thermal
pressurization in earthquake slip motivates further fault-scale models;
see [Rice (2006)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JB004006).

Return to the [overview](index.md).
