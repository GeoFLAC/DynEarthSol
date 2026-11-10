# Thermal pressurization of a heated material cell

A sealed, rigid porous cell warms by 10 K. Its fluid cannot escape and
its bulk volume cannot change, so the fluid's thermal expansion raises
pore pressure. The selected explicit result matches the linear reference
and its supplied-heat ledger. Separate controls measure the PT coupling
lag and the missing open-system energy extension.

![Sealed heated cell](../../benchmarks/poromechanics/results/runs/thermal-cell-20261108/native2d/fig/setup.svg)

The cell is 1 m on each side; the 2D calculation uses unit thickness.
Gravity, permeability, conduction and solid thermal expansion are zero.
The uniform heat source acts for 0.25 s; the run then holds through 1 s.
Temperature starts at 273 K. The material has $K=20$ MPa, $G=10$ MPa,
$K_s=1000$ MPa, $K_f=100$ MPa and $\phi=0.25$.
Here $\alpha=1-K/K_s=0.98$ and $\alpha_f=2\times10^{-4}$ K$^{-1}$.

## Temperature, pressure and stress references

Define the temperature rise $\theta=T-T_0$, inverse Biot modulus $S$,
fluid expansion coefficient $\beta$, and reference heat capacity $C_0$:

$$
S=\frac{\phi}{K_f}+\frac{\alpha-\phi}{K_s},\qquad
\beta=\phi\alpha_f.
$$

$$
C_0=(1-\phi)\rho_s c_s+\phi\rho_f c_f
=3.0325\times10^6\ {\rm J/(m^3 K)}.
$$

The adopted linear fluid content is

$$
\zeta=S p+\alpha\varepsilon_v-\beta\theta.
$$

For a sealed, fixed-strain cell, $\zeta=\varepsilon_v=0$. Therefore

$$
p=\frac{\beta}{S}\theta,\qquad
\Delta\sigma_{ii}=-\alpha p.
$$

The normal total stress is negative in the tension-positive convention.
No mechanical strain is needed to build this pressure. The source supplies
$\dot Q/V=C_0(40\ {\rm K/s})$ during the pulse, giving

$$
\theta(t)=40\min(t,0.25)\ {\rm K},\qquad
E-E_0=Q_{\rm supplied}=VC_0\theta.
$$

Nodal source strengths equal the uniform power density times independently
integrated nodal volumes. Thus uniform heating comes from actual input
sources, rather than modifying saved temperatures.

For the additional cell with zero recorded mean total stress, elasticity
gives $B\varepsilon_v-\alpha p=0$, hence

$$
p=\frac{\beta\theta}{S+\alpha^2/B}.
$$

The recorded in-plane mean uses $B=K+G/3$ in 2D; the 3D mean uses $B=K$.
The 2D free faces do not prescribe the out-of-plane stress to zero.

## Selected sealed-cell result

![Explicit pressure and closed heat balance](../../benchmarks/poromechanics/results/runs/thermal-cell-20261108/native2d/fig/result.svg)

At the end of the pulse, the fixed-strain reference is **154.799 kPa**.
The explicit pressure, temperature and fixed-capacity heat ledger agree
to round-off. Native and HDF5 output in tiny 2D/3D cells give the same
physical results. The free-stress PT controls also match their analytical
pressure and zero recorded normal stresses.

A separate explicit restart changes the step from 0.01 to 0.025 s during
heating. Temperature, pressure and supplied heat remain within
$3\times10^{-14}$ of their respective reference scales. This checks
accepted thermal increments across a step change; it does not qualify
pending shear-heating energy, which is a separate required feature.

## PT timestep control

PT consumes the previous pressure increment. The fixed-strain cell
therefore retains the split's numerical capacity $L=\alpha^2/B$:

$$
S p_n+L(p_n-p_{n-1})=\beta\theta_n.
$$

Consequently its exact discrete reference is

$$
p_n=\frac{\beta\theta_n+L p_{n-1}}{S+L},\qquad p_0=0.
$$

The solver matches this recurrence. It approaches the physical pressure
as the physical step shrinks, and keeps relaxing the pending increment
after heating stops. During the heating ramp, at 0.2 s:

| PT step, s | Pressure error / physical reference |
| --- | --- |
| 0.010 | 49.65% |
| 0.005 | 30.31% |
| 0.00125 | 7.964% |
| 0.000625 | 3.982% |

The small-step error falls at first order. The 0.01 s PT case is not an
accurate representation of this rapid sealed heating pulse, even though
its mechanical force converges. The selected explicit cell produces the
physical pressure directly. This is a coupling-time error, distinct from
a failed local material return or a lack of mechanical equilibrium.

![PT and open-energy controls](../../benchmarks/poromechanics/results/runs/thermal-cell-20261108/native2d/fig/controls.svg)

## Drained-cell energy limit and its analytical deficit

In the drained control every material-cell vertex holds $p=0$, while
the bulk strain remains zero. The fluid content and cumulative outward
water volume are then

$$
\zeta=-\beta\theta,\qquad V_o=V\beta\theta.
$$

The boundary reaction measures this volume even though a single
uniform-pressure material cell has no resolved Darcy boundary layer.
This is a drained constitutive limit, not a spatial drainage benchmark.

An open-system energy ledger must use the mass still in the cell and
the enthalpy carried away. These conventions are described in
[MOOSE's heat-flow equations](https://mooseframework.inl.gov/modules/porous_flow/governing_equations.html#heat-flow).
The following deficit is derived here for DES's linear fluid-content law.

With $a=\rho_f c_f\beta$ and constant specific heats, the retained
capacity would be $C(\theta)=C_0-a\theta$. Relative to the initial state,

$$
\Delta E_{\rm retained}
=V[(C_0-a\theta)(T_0+\theta)-C_0T_0].
$$

Integrating the outflow enthalpy using its changing temperature gives

$$
H_o=Va\left(T_0\theta+\frac12\theta^2\right).
$$

The current fixed-capacity temperature equation supplies $VC_0\theta$.
Combining these independent quantities leaves

$$
\Delta E_{\rm retained}+H_o-Q_{\rm supplied}
=-\frac12Va\theta^2.
$$

For 10 K this is **−10,450 J**, or **0.03446%** of supplied heat.
The measured deficit is the same at steps 0.02, 0.01 and 0.005 s and in
both dimensions/output formats. Time refinement cannot remove it.
The implemented fixed-capacity ledger closes; the retained-fluid ledger
does not. This isolates the previously unqualified storage extension.

For this rigid, zero-pressure control there is no bulk mechanical work or
pressure-boundary work to explain the difference. A consistent retained
capacity would instead satisfy

$$
C_0\theta-\frac12a\theta^2=Q_{\rm supplied}/V.
$$

This control does not justify replacing a capacity in isolation for
general deforming flow. The extension still needs shared mass/energy
accounting, inlet and outflow enthalpy, volume changes, remapping,
restart and PT's pending coupling contribution. Thermal-plume
applications remain unqualified until those checks pass. Nonzero solid
thermal expansion and shear heating retain their separate references.

## Reproduce

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
solver=benchmarks/poromechanics/results/runs/restart-timestep-20261108/candidate/dynearthsol2d
python3 benchmarks/poromechanics/thermal_cell.py "$solver" /new/thermal-cell \
  --dt .01 --figures
python3 benchmarks/poromechanics/thermal_cell.py "$solver" /new/thermal-fine \
  --methods pt --boundaries fixed --dt .000625
python3 benchmarks/poromechanics/thermal_cell.py "$solver" /new/thermal-cell \
  --replay --figures
```

Use the corresponding frozen HDF5/3D executable and `--ndims 3` for the
tiny dimensional/backend controls. Every input, source recipe and saved
frame is retained in
[`thermal-cell-20261108`](../../benchmarks/poromechanics/results/runs/thermal-cell-20261108/).
The solver qualification is the same frozen `50c4001` source/binary
manifest used by the [coupled water-budget chapter](content.md).
