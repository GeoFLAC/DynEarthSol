# Scheduled fluid and heat point sources

The recovered private point-source module is integrated with the current
hydraulic and thermal subcycling. `injection` supplies fluid volume and
`heat_injection` supplies heat. Both are optional, default off, and share
one schedule/mapping implementation. Fluid rates are positive for supply;
the existing well convention remains positive for extraction. Their
contributions coexist in the hydraulic balance diagnostic as net extraction.

Each section registers `enabled`, `num_points`, `points_x/y/z`,
`points_unit`, `rate_model`, `rate`, `total_amount`, `start_time_in_yr`,
`end_time_in_yr` and `duration_in_yr`. Coordinates use the DES coordinate
system (z negative below the top); units are mm, cm, m or km. In 2D,
`points_z` is the vertical coordinate and `points_y` is a mutually exclusive
legacy alias. Rates and amounts may broadcast; coordinate arrays must have
exactly `num_points` entries. Source times use years at input and seconds
internally. End time and duration are alternatives.

`constant_rate` supplies the signed rate over its active interval.
`total_amount` divides its signed amount by a finite positive duration.
Fluid units are m³/s and m³ in 3D; in 2D they are per unit out-of-plane
thickness. Heat uses W and J with the same 2D thickness convention.
Fluid sources require hydraulic transport. Heat sources also work with
thermal diffusion disabled. Their capacity follows the selected thermal
model: existing solid capacity by default, mixed solid/fluid capacity when
fluid thermal advection is enabled. A supplied fluid volume carries no
enthalpy unless `injection.temperature` gives its temperature $T_{in}$ per
point (K; fluid thermal advection required). An injecting node then gains
$\rho_f c_f\,(\sum_p Q_{pn}T_{in,p}-T_n\sum_p Q_{pn})$, with $\rho_f c_f$
of the reference material; extracting points add nothing, since the
extracted fluid leaves at the local temperature, which the advective
transport already carries. The explicit heat-step bound uses each source's
active rate; its contribution uses the exact overlap with each substep.
This prevents a short pulse from bypassing the bound through a small
whole-step average. With fixed heat capacity, this mixing term does not
account for heat capacity carried into retained-fluid storage. Changing-storage
energy and deformation work remain unqualified. The current injection/outflow
tests use constant capacity, a warm region inside the domain and no deformation
work; the startup flux lag is reported separately. See the
[advective heat equation and scope](heat.md). Without the option, a heat
source remains an independent input.

## Reference and discrete balance

For a signed rate $R_i$ active over $[a_i,b_i]$, the exact supplied amount is

$$
A_i(t)=R_i\max(0,\min(t,b_i)-a_i).
$$

On an integration interval $[t_0,t_1]$, use
$(A_i(t_1)-A_i(t_0))/(t_1-t_0)$ as the mean rate. Hydraulic and thermal
substeps each use their actual interval. Subcycling does not repeat the
whole physical-step source amount.

At source position $\mathbf x_i$ in a simplex, distribute the rate as
$f_n=R_iN_n(\mathbf x_i)$. The shape functions sum to one. Shared faces
therefore map to the first containing element once, with deterministic
ownership. Small barycentric roundoff is clamped and the weights normalized;
the assembled signed total is independently checked. The current mesh is
used after motion/remeshing. Unmappable points and sources that reach
prescribed-pressure/temperature nodes are rejected. A fixed node must not
silently discard part of a prescribed source.

For the isolated closed-box verification, mechanical pressure feedback,
thermal expansion and transport fluxes are disabled. Fixed capacity gives

$$
\sum_n m_{h,n}(p_n-p_{n,0})=\sum_i[A_i(t)-A_i(0)],\qquad
\sum_n m_{T,n}(T_n-T_{n,0})=\sum_i[E_i(t)-E_i(0)].
$$

The reference masses are assembled independently from simplex volumes,
$m_{h,n}=\sum_{e\ni n}V_eS_B/N_e$ and
$m_{T,n}=\sum_{e\ni n}V_eC/N_e$. The quasistatic pressure split uses
$S_B=S_0+\alpha^2/K_*$, where $K_*=K$ in 3D and $K_*=K+G/3$ for the
existing 2D in-plane mean-stress convention. This checks the source and
storage bookkeeping; it is not a spatial singular-source solution or a
validation of general enthalpy conservation with changing storage.

The [selected review case](../../benchmarks/poromechanics/results/review/index.html#point_sources)
includes its setup schematic, exact schedule integrals and one comparison
of stored volume/heat. Separate remesh tests establish deterministic source
relocation and continuation, rather than conservative fluid/heat remapping.

## Inlet pulses and the thermal bound

The inlet-temperature term is a relaxation coefficient. Averaging a short
pulse over the whole physical step can underestimate that coefficient in
a thermal substep. The bound therefore uses the active rate of every inlet
whose window intersects the physical step. The actual source still uses
the exact substep/window overlap. This bound is conservative when different
inlets have disjoint active windows; resolving a short pulse with the
physical timestep also improves temporal accuracy.

A clamped, insulated cell with zero permeability isolates this coefficient.
With injection on a vertex, its nodal temperature obeys the implemented
fixed-capacity mixing equation. Let $A$ be the total injected amount and
$C_n$ the nodal heat capacity. The independent final-temperature reference is

$$
T_f=T_{in}-(T_{in}-T_0)\exp[-\rho_f c_f A/C_n].
$$

For $T_0=283$ K, $T_{in}=333$ K and $\rho_f c_f A/C_n=2.25$, the
reference is 327.73004 K. The old whole-step bound gives 395.5 K for a
pulse lasting 1% of one physical step. The active-rate bound keeps early,
late and total-amount pulses within 283–333 K in both dimensions. Reducing
the physical step from 0.001 s to 0.0001 s reduces the final ODE error from
2.7233% to 0.2674% of the temperature jump. Coarse boundedness alone does
not establish temperature accuracy.

This is a numerical regression for the implemented mixing equation. Its
large injected amount deliberately exposes the coefficient bound; it is
outside the small-storage physical energy qualification. It does not close
the changing-storage energy requirement or remove the separate startup
Darcy-flux lag.

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 tests/functional/hydromechanics/check_injection_energy.py \
  ./dynearthsol2d /new/inlet-pulse-2d --pulse --ndims 2
python3 tests/functional/hydromechanics/check_injection_energy.py \
  ./dynearthsol3d /new/inlet-pulse-3d --pulse --ndims 3
```

Inputs, raw frames, parent failures and source/binary hashes are retained
under `results/runs/inlet-pulse-20261108/`.

## Recorded qualification

Evidence is in `results/runs/private-des-integration-20261107/`:

- Both production builds and parent `make set` / candidate `make cmp` pass.
  Assembly checks pass 7,270/10,722 assertions in 2D/3D.
- The disabled-feature 94-check regression preserves 360 save/checkpoint
  files and 10,714 physical payloads exactly against the heat-port parent.
- The source runner passes 23 cases in 2D and 22 in 3D: signed fluid/heat
  balances, partial intervals, rate/amount and end/duration equivalence,
  coordinate units/2D alias, subcycling, mixed capacity, source/well
  coexistence, invalid inputs and fixed-node overlap rejection.
- Physical state and actual saved times remain exact on restart. Forced
  remeshing also preserves every continuation field exactly in both
  dimensions. The half-step non-PT run exposes a pre-existing few-ulp
  difference in `reference_frame_time`; this output-clock control value is
  recorded separately in that historical study. The subsequent
  [output-clock correction](mechanics.md#restart-output-clock) has exact
  continuation metadata in the new source checks.

```sh
python3 tests/functional/hydromechanics/check_point_sources.py \
  ./dynearthsol2d /new/points-2d --dimension 2
python3 tests/functional/hydromechanics/check_point_sources.py \
  ./dynearthsol3d /new/points-3d --dimension 3
```

Use one OpenMP/BLAS thread and the 6 GiB process limit for these runs.
