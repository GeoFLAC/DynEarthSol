# Explicit hydraulic extensions

The transport operator, its material properties and its accepted-step history
stay in the existing DES modules. CPU loops retain OpenMP `default(none)`;
there is no implicit integration, global nonlinear solve or separate solver
framework. Component checks link the production routines.

## Darcy transport and stability

`mat.permeability_model=2` selects the power law
`k=k0*(phi/phi0)^n`, with positive dimensionless material values in
`mat.permeability_exponent` (default `[3]`). It uses the existing `MatProps`
owner and current porosity. For mixed elements, `k0` retains the harmonic
material mean, while `n` uses the marker-count arithmetic mean; `phi0` is
the existing marker-weighted reference porosity. Model 0 remains fixed
permeability and model 1 remains Kozeny–Carman. A zero reference permeability
stays zero. This power-law choice is also used in
[Yang et al. (2021)](https://doi.org/10.1029/2021JB021787).
The [plastic property check](plastic-content.md#plastic-pore-recovery-with-finite-grain-stiffness)
qualifies both signs of dilation with finite grain stiffness and adds recovery
of accepted plastic pore opening. It retains the current plastic Biot
coefficient and does not establish finite-density mass or energy conservation.

`mat.permeability_tensor = [xx,yy,zz,xy,xz,yz]` is a dimensionless symmetric
positive-definite tensor multiplying the material's scalar `hydraulic_perm`
and optional Kozeny-Carman factor. Its orientation is common to all materials;
material magnitudes remain distinct. In 2D its x-z submatrix is used.

Transport computes one gradient of `p + rho_f*g*z` per element, then applies
the permeability tensor and distributes its residual through the existing
nodal support. This replaces the repeated pressure/gravity node-pair sums.
A row-sum bound on the assembled, lumped operator limits explicit substeps,
including when the operator sets `fixed_dt`. The bound is cached until mass,
material or mesh assembly invalidates it. Only pressure and optional matrix
pressure advance during these substeps; mechanics and porosity advance once
per accepted physical step. Scheduled sources use exact subinterval integrals.

Coupled explicit mechanics (`has_PT=false`, pressure feedback enabled) uses
Biot pressure capacity `1/M` and the accepted projected volumetric strain
rate, with pressure advanced before stress and force. This includes scaled
mass and GVS. Physical inertia uses this ordering also when pressure feedback
is disabled. PT and quasistatic flow with feedback disabled retain the
fixed-stress capacity and ordering. The
[closed-column wave qualification](dynamics.md) derives the
two forms, automatic undrained-wave timestep and checkpoint compatibility.

## Boundary data and balance

Boundary types are 0 sealed, 1 initial pressure held, 2 prescribed absolute
pressure, and 3 prescribed outward Darcy flux. `bc.hbc_pressure` (Pa) and
`bc.hbc_flux` (m/s) each contain six values ordered x0 x1 y0 y1 z0 z1.
`bc.hbc_schedule_time_in_yr` / `bc.hbc_schedule_ratio` define one shared
piecewise-linear multiplier, with `bc.num_hbc_schedule` entries; repeated times
are steps. Values are held outside the table. Type 2 applies at initialization
and at substep endpoints; type 3 uses exact interval averages and true facet
measure, including side faces. Type-2 values must agree at shared tagged
nodes. Imported mesh labels do not imply an intersection: disjoint tagged
boundaries may prescribe different pressures, including the two arcs of an
annulus. The check runs after boundary flags are rebuilt at initialization,
restart and remeshing; conflicting values report the node and tag indices
with configuration error 11. Built-in boxes also retain the early face check.
Type 2 takes precedence where it meets a legacy held-pressure face.

The save/checkpoint field `fluid balance` contains five last-step increments:

1. Lumped pressure storage change, including any matrix reservoir.
2. Accepted mechanical/thermal source contribution.
3. Well extraction volume (positive outward).
4. Net imposed/inferred boundary volume (positive inward).
5. Residual: storage + source + wells - boundary.

This diagnoses conservation of the actual discrete update with its frozen
coefficients. It is not a claim of exact nonlinear fluid-content conservation
when porosity/storage evolve, or across remeshing. The five values are restored
from checkpoints to preserve exact diagnostic continuation.

### Remapping linear Biot content

The fixed-property saturated remap constrains the existing linear model's
excess fluid content,

$$
\zeta=S_0p+\alpha\epsilon_v-\beta(T-T_0),
\qquad S_0=\frac{\phi}{K_f}+\frac{\alpha-\phi}{K_s}.
$$

With mechanical feedback disabled, it uses the reduced flow model's
effective capacity and thermal coefficient.
Here `beta = phi*alpha_f + (alpha-phi)*alpha_s` when thermal pressurization
is enabled, and `T0` is `bc.surface_temperature`. The recorded mean-stress
modulus `Kb` is `K+G/3` in 2D and `K` in 3D.

This is an excess-content measure; it does not supply the total retained
mass needed by the evolving-storage energy model. The numerical fixed-stress
buffer and any immobile reservoir are transferred as separate constraints:

$$
R=-\sum_e V_e L_e\overline{\mathtt{dppressure}}_e,
\qquad L_e=\alpha_e^2/K_{b,e}.
$$

`FluidRemap` snapshots these quantities before interpolation. Exact
old/new simplex intersections determine content retained in the new domain.
Uncovered new regions use the existing remapper's provisional state;
removed regions keep their own content outflow. The intersection integrates
products of linear shape functions. Tetrahedra retain face topology during
clipping, with extended-precision intermediates for nearly coincident faces.

After rebuilding markers, properties and masses, the correction minimizes
the storage-weighted squared pressure change subject to the global content
constraint. It therefore adds a common pressure offset to free-pressure
nodes. Prescribed nodes remain fixed; an all-Dirichlet mismatch is reported
as a boundary reaction. The pending buffer and matrix content have their
own corrections. Total stress follows the reconstructed mechanical pressure
reference, preserving effective stress at that reconstruction stage.

This is global conservation of the linear content measure. It does not
impose local conservation or a pressure maximum principle. Retention and
evolving porosity retain their existing remap until their inventory laws
are extended; thermal energy has its own unfinished qualification. The
remap log records geometric removal, addition and Dirichlet reaction
separately from the five physical-step balance entries above.

Forty focused native/HDF5 cases cover both dimensions, explicit/GVS/PT,
matrix storage, thermal sources, reduced flow and swept domains. The worst
normalized content defect is `1.01e-12`; 35,590 continuation payloads match
exactly. The unchanged-domain parent control changes integrated content by
`3.585e-6` m² in 2D. Independent overlap tests include the two nearly
coincident-face cases found during actual remeshing. The
[qualification record](../../benchmarks/poromechanics/results/runs/finish-20261109/fluid-remap/qualified2/qualification.json)
retains source/binary hashes and the separate monitor issue identified on
the parent commit.

The [deforming-column water budget](content.md) additionally integrates
physical Biot content from saved pressure and strain. It checks basal inflow,
an open drain and changing GVS timesteps, and identifies the finite pending
increment that distinguishes the PT split's conserved quantity from physical
content. Its independent mixed-boundary analytical solution and mesh/time
controls retain that first-order coupling error explicitly.

## Groundwater heat transport

The integrated saturated heat model uses physical Darcy flux from
`hydraulics.cxx::hydraulic_element_flux`, mixed solid/fluid capacity and a
stabilized conduction/advection operator in `fields.cxx`. See
[heat transport](heat.md) for the registered controls, equation and measured
Ogata–Banks mesh/time/restart qualification. [Point sources](sources.md)
use one shared source mapper and exact interval integration for fluid and heat.

For constant capacity and divergence-free flow, the primitive-temperature
advection term is consistent with conservative sensible-heat transport.
A general fluid enthalpy flux has the form

$$
\mathbf j_E=-\boldsymbol\lambda_{\mathrm{eff}}\nabla T
             +\rho_f h_f\mathbf q,
\qquad
C_{\mathrm{bulk}}=(1-\phi)\rho_s c_s+\phi\rho_f c_f.
$$

General enthalpy conservation with changing porosity/compressible storage,
pressure work, injection temperature and outflow energy accounting still
needs a consistent extension and separate qualification. The physical flux
helper includes both pressure and gravity; the volume-weighted hydraulic
residual cannot be used directly as an advection velocity. Unsaturated
heat capacity, thermal dispersion and nonlinear thermal applications remain
outside the completed constant-capacity benchmark scope.

## Thermal pressurization

Heat diffusion also uses explicit substeps bounded by the assembled lumped
heat-capacity/operator row ratio. This applies to `control.fixed_dt` as well
as an automatically selected physical step. The conservative bound is half
the forward-Euler spectral limit; it is not a monotonicity guarantee on an
arbitrary simplex mesh. The cached bound is invalidated when masses are
assembled, including remeshing and material changes. Temperature advances
through the substeps, while the saved old temperature remains the start of
the physical step, so mechanical and fluid thermal increments are consumed
once. Isothermal benchmarks should set `control.has_thermal_diffusion=false`.

`control.has_thermal_pressurization=true` requires thermal and hydraulic
transport and an elastic bulk response. The existing material `alpha` is
volumetric, consistent with the density law. The thermal stress increment is
`-K*alpha_s*Delta T` on each normal component, including out-of-plane stress
in plane strain. A nodal temperature baseline is taken before the physical
heat update; repeated mechanical trials read the same baseline.

For `beta = phi*alpha_f + (alpha_b-phi)*alpha_s`, the adopted fluid content is
`zeta = alpha_b*eps_v + p/M - beta*T`. For an elastic fixed-stress source,
eliminating strain from the recorded mean total stress gives the hydraulic thermal source coefficient
`beta - alpha_b*K*alpha_s/K_b`, where `K_b` is the existing mean-stress modulus.
The physical strain-driven branch uses `beta` directly.
The fixed-stress plastic path also takes accepted strain directly, uses
`beta`, and retains its raw pending-pressure buffer; see the
[yielding-cell balance and thermal reference](plastic-content.md).
Porosity history includes the grain-expansion increment `-alpha_s*Delta T`.
The direct constitutive checks cover sealed fixed-stress heating, drained
outflow, fixed-strain thermal stress and repeated trials. Executable lifecycle
checks combine heat, porosity/permeability and matrix exchange. Nonlinear
plastic/RSF thermal applications need their own qualification.

The [uniform heated cell](thermal-cell.md) adds analytical pressure, strain
constraints and a supplied-heat ledger in both dimensions and output formats.
It measures the PT pressure lag separately, and gives a drained-cell
counterexample to treating the fixed-capacity heat model as a complete
retained-fluid enthalpy balance.

See [MOOSE's governing equations](https://mooseframework.inl.gov/falcon/modules/porous_flow/governing_equations.html)
for the tensor-Darcy and thermo-poroelastic conventions. The fixed-stress
elimination above follows DES's recorded stress convention.

## Local matrix exchange

`mat.matrix_storage` (1/Pa per bulk volume) and `mat.matrix_exchange`
(1/(Pa s)) enable an immobile matrix reservoir in elements with nonzero
exchange. Its pressure starts at the element's mean transport pressure.
For transfer `q = exchange*(p - p_matrix)`, the matrix gains `q*dt`, and the
transport residual loses exactly the same volume. Substeps obey both the
transport row bound and the reservoir-capacity bound.

There is one optional element pressure history, written to output/checkpoints
as `matrix pressure`. Restart requires that field when exchange is enabled.
It follows the existing element remapper. This is local dual-storage exchange,
not a separately meshed fracture/matrix mechanical model or a proof of
conservative reservoir remapping. Checks compare homogeneous two-reservoir
relaxation with its exponential solution, verify first-order time convergence,
conservation and stability under a large physical step.

## Verification commands

Build dimensions sequentially with `make -j2 ndims=2` and `make -j2 ndims=3`.
For each dimension:

```sh
OMP_NUM_THREADS=1 make -C tests/residual-assembly-test check ndims=2 backend=omp
OMP_NUM_THREADS=1 python3 tests/functional/hydromechanics/check_boundaries.py ./dynearthsol2d /new/boundary-results
OMP_NUM_THREADS=1 python3 tests/functional/hydromechanics/check_poroperm.py ./dynearthsol2d /new/property-results
python3 tests/functional/hydromechanics/check_hydraulic_input.py ./dynearthsol2d
```

Replace 2d/ndims=2 with 3d/ndims=3 for the other dimension. Results directories
must be new. See the [overview](index.md) for the reference comparisons
and their limits.
