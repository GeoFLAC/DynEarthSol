# Groundwater heat transport

`control.has_fluid_thermal_advection=true` activates the integrated private
DES heat-transport model. It requires hydraulic transport and saturated
materials. The default is false. With hydraulic and thermal transport active,
fluid thermal expansivity also contributes Boussinesq buoyancy to Darcy flow.
Storage, thermal capacity and inertia retain the reference fluid density;
tensor permeability, boundary conventions and hydraulic subcycling are shared.

## Equation and scope

With temperature-dependent work terms omitted, the implemented transport is

$$
C\dot T+\rho_{f0}c_f\mathbf q\cdot\nabla T
=\nabla\cdot(\lambda\nabla T)+r,
\qquad C=(1-\phi)\rho_s c_s+\phi\rho_{f0}c_f.
$$

`mat.fluid_heat_capacity` supplies $c_f$ in J/(kg K). The Darcy flux is

$$
\mathbf q=-\frac{\boldsymbol\kappa k_r}{\mu}
(\nabla p+\rho_{f0}[1-\alpha_f(T-T_0)]g\mathbf e_z),
\qquad T_0=T_{\mathrm{surface}}.
$$

This is Darcy flux, not pore velocity. For constant properties in a column,
the thermal-front speed is $v_T=\rho_{f0}c_f q_x/C$. The flux helper uses
the same pressure and gravity operator as hydraulics. Buoyancy uses
`mat.fluid_alpha` when hydraulic transport and either thermal diffusion or
fluid thermal advection are enabled. It introduces no additional switch.
Zero expansivity or zero gravity uses the existing pressure/gravity arithmetic.
This Boussinesq approximation is distinct from a general fluid equation of state.

The primitive-temperature advection equation is equivalent to conservative
sensible-heat transport for constant capacity and divergence-free flow.
General enthalpy conservation with changing porosity, compressible storage,
fluid sources and mechanical work requires additional terms and separate
qualification. Unsaturated advection is rejected pending a consistent
saturation-dependent heat capacity.

## Boundaries, initialization and stable steps

The imported `bc.tbc_x0` through `bc.tbc_z1` select insulated (0) or
prescribed temperature (1) faces. Their `bc.tbc_val_*` values are in kelvin;
omitted values use the existing surface/mantle defaults. Inconsistent
temperatures on intersecting prescribed faces are rejected. Prescribed
temperatures are restored after remeshing.

The imported `ic.temperature_option=4` supports the registered Gaussian
temperature parameters. `ic.pore_pressure_option=1` initializes an affine
pressure field from the opposite prescribed x-face pressures, while option
2 initializes a uniform `ic.pore_pressure_value`.

`ic.temperature_option=5` adds the same Gaussian anomaly to the linear
surface-to-mantle conduction profile. For flow, the background gradient comes
from the actual active top/bottom prescribed thermal values, including their
default fallbacks. An insulated end supplies no prescribed background slope.
The background hydrostatic potential and temperature departure are assembled
separately so that a homogeneous conductive box can remain at rest.

The advection/conduction element operator adds symmetric graph diffusion,
limited to the advective coupling. At zero Darcy flux it recovers ordinary
Galerkin conduction, including on obtuse elements. Consequently its
off-diagonal entries need not be nonpositive on every mesh: the method does
not guarantee a maximum principle for arbitrary conduction meshes. A
lumped-capacity absolute-row bound selects thermal substeps, including when
the outer physical step uses `control.fixed_dt`. Changing Darcy flux refreshes
that bound. Mechanical and hydraulic histories still advance on their own
intervals.

The port also restores `control.has_poroelastic_pressure_feedback` (default
true), `mat.derive_biot_coeff_from_bulk_moduli` (default false), and private
`bc.hbc_val_*` pressure aliases. The feedback switch disables only the
mechanical pressure source, not the selected storage. Derived Biot
coefficient is $1-K_d/K_s$; invalid moduli and porosity/Biot combinations
are rejected. Modern hydraulic type 2 remains prescribed pressure; the
legacy type-1/value pair is translated explicitly rather than changing
modern type-1 drainage.

## Qualification

The [guided review](../../benchmarks/poromechanics/results/review/index.html#ogata_banks)
includes a setup schematic, the Ogata–Banks reference, its derivation and
one selected 3D result. The reference solves a constant-coefficient,
semi-infinite inlet-step problem. The finite DES column is compared before
the front reaches its downstream boundary. Parameters are deliberately
scaled; a very small fluid bulk modulus enlarges hydraulic storage without
changing the steady Darcy flux. This isolates transport rather than
representing an aquifer application.

Recorded evidence under `results/runs/private-des-integration-20261107/`:

- Both dimensions build; assembly checks pass 7,249/10,693 assertions.
- Parent `make set` and candidate `make cmp` pass in 2D/3D. The 94-check
  regression preserves 360 save/checkpoint files and 10,714 physical
  payloads exactly with the new controls disabled.
- Separate spatial sequences reach maximum all-node/saved-time errors of
  1.304% (2D) and 1.306% (3D) of a 1 K inlet jump at spacing 0.003125 m.
  Final spatial orders are 0.997 and 1.000. Independent time halving gives
  orders 1.020 and 1.010.
- Steady pressure remains exact. Physical restart is exact for every
  recorded transport run. Fluid-only capacity limits pass in both
  dimensions. A large fixed-step case remains bounded and accurate.
- Twenty-six input/compatibility checks cover invalid parameters,
  conflicting temperature corners, legacy pressure aliases and the
  independently derived Biot coefficient.

Recompute these metrics from the frozen, checksummed outputs with:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 \
  benchmarks/poromechanics/qualify_heat.py \
  benchmarks/poromechanics/results/runs/private-des-integration-20261107
```

`tests/functional/hydromechanics/check_heat_transport.py` produces new
solver runs and records executable, input, source-diff and output hashes.
Scheduled [fluid/heat point sources](sources.md) have separate transport
and input-temperature checks. These heat-transport comparisons do not
establish a complete earthquake-cycle energy model.
