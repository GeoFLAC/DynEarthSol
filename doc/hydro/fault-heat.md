---
title: Fault-zone thermal weakening
---

# Fault-zone heating at fixed normal traction

A sealed, prestressed porous patch slides at 1 m/s under 50 MPa total
normal compression. Plastic work heats the mixture; expanding water raises
pore pressure and lowers shear resistance. Ordinary explicit mechanics uses
physical mixture mass and zero damping. PT is a separate comparison.

The worked example supplies the setup, selected figure, input and tables.
This is a homogeneous local fault-zone model with one free normal platen
mode. It has no spatial rupture front or complete earthquake cycle.

## Setup and input

The reference square is 10 mm high and wide, with two triangles and four
nodes. Its base is fixed. Both side faces prescribe
$v_x=V(z+H)/H$, with $z\in[-H,0]$ and $V=1$ m/s. The top has a free,
uniform normal velocity (`vbc_z1=8`) under prescribed total traction
$\sigma_{zz}=-N$. Every hydraulic and thermal face is sealed and insulated;
gravity, permeability, conductivity and solid thermal expansion are zero.
The mesh stays fixed, so this checks the integrated constitutive path in
the reference geometry, rather than deformation of a finite fault layer.

| Setting | Value |
| --- | --- |
| Drained bulk and shear moduli, $K,G$ | 20 and 10 GPa |
| Cohesion, friction, dilation | 1 MPa, 30°, 0° |
| Biot coefficient, porosity | 0.8, 0.2 |
| Grain and water bulk moduli | 100 and 2.2 GPa |
| Grain and water densities | 2700 and 1000 kg/m³ |
| Grain and water heat capacities | 1000 and 4180 J/(kg K) |
| Water volumetric thermal expansion | $2\times10^{-4}$ K⁻¹ |
| Initial $T,p$ | 300 K, 0 Pa |

Initially $\sigma_{xx}=\sigma_{yy}=\sigma_{zz}=-N$ and
$\tau_0=N\sin\phi+c\cos\phi=25.866$ MPa, on the smooth Mohr–Coulomb
face. The mixture density is 2360 kg/m³ and its fixed volumetric heat
capacity is $C=2.996\times10^6$ J/(m³ K). The thermal-pressure coefficient
is $\beta=\varphi\alpha_f=4\times10^{-5}$ K⁻¹. The elastic storage is
$S=\varphi/K_f+(\alpha-\varphi)/K_s$.

From the DES repository, generate and check the selected case:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 benchmarks/poromechanics/fault_heat.py \
  --exe ./dynearthsol2d --output OUT --mode physical --strain 10 --samples 600
```

`OUT/cell.cfg` is the complete input; `cell.info`, saves and checkpoints
contain the numerical fields. `history.csv` pairs every saved state with
the independent reference; `summary.json` records errors and binary/source
identity. The column `sampled_work_J_m3` is a sparse power integral in this
case. It is not an exact accepted-step heat ledger.

## Reference and derivation

The reference integrates the continuous Mohr–Coulomb/Jaumann equations,
independently of the DES return map. Prime denotes differentiation with
respect to engineering shear strain $\gamma=Vt/H$; set
$\ell=K-2G/3$, $s=\sin\phi$ and

$$
d=\frac{\sigma_{xx}+N}{2},\quad
m=\frac{\sigma_{xx}-N}{2},\quad
q=\sqrt{d^2+\tau^2},\quad a=d/q,\quad b=\tau/q.
$$

The smooth-face yield condition is
$q+s(m+\alpha p)-c\cos\phi=0$. With zero dilation, write the plastic
increments as $\epsilon^p_{xx}{}'=ra/2$,
$\epsilon^p_{zz}{}'=-ra/2$ and $\gamma^p{}'=rb$; their work is $rq$.
Let $e=\epsilon_v'$ and let the heat fraction $\chi=1$ (zero for the
no-heating control). Fluid content and energy give

$$
S p'+\alpha e=\beta T',\qquad C T'=\chi r q.
$$

Normal equilibrium and yield consistency determine $e,r$ at each state:

$$
\left(\ell+2G+\frac{\alpha^2}{S}\right)e
+\left(Ga-\frac{\alpha\beta\chi q}{SC}\right)r=\tau,
$$

$$
\big[s(\ell+G)-Ga\big]e-Gr=-Gb.
$$

Then integrate

$$
\sigma_{xx}'=\ell e-Gr a+\tau-\alpha p',\qquad
\tau'=G-Gr b-d,\qquad W^{p}{}'=rq.
$$

The reference also checks that
$\sigma_{yy}=-N+\ell\epsilon_v-\alpha p$ remains the intermediate
principal stress. The solution stays on that face through $\gamma=10$.
Tightening the integration tolerance from $10^{-10}$ to $10^{-11}$
changes normalized stress by less than $1.6\times10^{-10}$.

An instructive closed-form limit follows when the elastic stress
differences are negligible: $q\simeq\tau$, $r\simeq1$ and
$\epsilon_v\simeq\alpha p/(\ell+G)$. Define

$$
\Lambda=\frac{\beta}{S+\alpha^2/(\ell+G)}.
$$

The leading thermal-weakening law is then

$$
\tau\simeq\tau_0 e^{-\alpha s\Lambda\gamma/C},\qquad
p\simeq\frac{\tau_0-\tau}{\alpha s},\qquad
T-T_0\simeq p/\Lambda.
$$

This limiting law explains the trend. The plotted reference is the finite
modulus differential solution above, including rotation and normal
compliance. Thermal pressurization as a mechanism for rapid fault weakening
is discussed by [Rice (2006)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2005JB004006);
the Mohr–Coulomb patch reference here is derived for these DES boundaries.

## Results and method checks

At 100 mm slip (0.1 s), temperature rises 70.147 K, pressure reaches
22.559 MPa, and shear resistance falls to 16.829 MPa: 34.9% weakening.
Maximum saved-frame shear and normal-stress errors are 0.0011% of
$\tau_0$ and 0.0022% of $N$, respectively. The normal reference imposes
equilibrium; physical inertia gives a small measured departure. These
metrics sample the saved frames, rather than bounding every wave peak.

At 30 mm slip, separate PT runs with $\Delta t=20,10,5$ µs have temperature
errors 0.0711%, 0.0356%, 0.0178%; halving the step halves the error. Their
accepted heat ledger includes undeposited plastic heat and closes within
$4.2\times10^{-13}$ of $\tau_0$ per unit volume. No-heating and
zero-water-expansion controls separate plastic heating from thermal
pressurization. Native and HDF5 restarts reproduce 18 endpoint field
payloads exactly; changed-step PT restarts retain the heat ledger.

Fast-slip GVS runs at scalings $10^4,10^5,10^6$ coincide because all reach
the configured shear-wave speed cap. Their fictitious density is 5400
kg/m³; this cap does not select the physical mixture density. A separate
slow phase at $10^{-5}$ m/s is restarted at 1 m/s with physical mass and
without GVS or PT. The new step is $6.88\times10^{-8}$ s, versus old GVS
steps $2.00\times10^{-3},2.00\times10^{-4},2.00\times10^{-5}$ s. Stress,
strain, pressure, temperature and pending heat are retained exactly.
The first 20 physical intervals close the accepted heat ledger within
$1.4\times10^{-14}$ of $\tau_0$. Prescribed acceleration changes kinetic
energy; artificial kinetic energy is not deposited as heat.

For example, reproduce the slow-to-fast check with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 benchmarks/poromechanics/fault_heat.py \
  --exe ./dynearthsol2d --output OUT-switch --mode gvs --scaling 1e5 \
  --velocity 1e-5 --strain .02 --physical-transition
```

## Limits

The material uses constant friction and zero dilation, with all accepted
plastic work converted to heat. There is no weakening from RSF, flash
heating, conduction, drainage or a propagating rupture in this case.
The fixed reference capacity excludes the general evolving-storage
enthalpy problem. This manual switch at a specified restart does not
qualify automatic earthquake events or long tectonic-time clock handling.

Return to [plastic shear heating](plastic-heat.md) or the [overview](index.md).
