# Water entering, stored and drained in a deforming column

Water enters the base of a laterally confined elastic column and leaves
through its drained top. The solid expands as pressure rises. This checks
the water budget with mechanical feedback, including a changing GVS
timestep, rather than checking only the pressure-equation residual.

![Column and boundary conditions](../../benchmarks/poromechanics/results/runs/fluid-content-20261108/spatial-pt/fig/setup.svg)

The height is 1 m, gravity is zero, the base is fixed vertically and the
top carries zero total traction. Basal inflow is $q_i=10^{-6}$ m/s;
$p=0$ at the top. The mesh stays fixed, while the material strain evolves.
The cross-section is $2h$ in the 2D refinement study; all volumes below
are normalized by area, so the physical one-dimensional problem is the same.

| Property | Value |
| --- | --- |
| Drained bulk and shear moduli $K,G$ | 20, 10 MPa |
| Biot coefficient $\alpha$; porosity $\phi$ | 0.8; 0.25 |
| Fluid and grain moduli $K_f,K_s$ | 100, 1000 MPa |
| Permeability $k$; viscosity $\mu$ | $10^{-11}$ m²; $10^{-3}$ Pa s |

## Analytical solution

Let $d$ measure depth from the top, with tension-positive strain and
compression-positive pore pressure. The pressure capacity and constrained
modulus are

$$
S=\frac{\phi}{K_f}+\frac{\alpha-\phi}{K_s},\qquad
H_m=K+\frac43G.
$$

Vertical equilibrium gives $\sigma_{zz}=H_m\varepsilon_v-\alpha p=0$.
The physical increment in Biot fluid content is therefore

$$
\zeta=S p+\alpha\varepsilon_v=C_o p,\qquad
C_o=S+\frac{\alpha^2}{H_m}.
$$

With mobility $m=k/\mu$, mass balance reduces to

$$
C_o p_t=m p_{dd},\qquad
p(0,t)=0,\qquad m p_d(H,t)=q_i.
$$

Initially pressure and strain are zero. Subtract the steady profile
$q_i d/m$ and expand the remainder in eigenfunctions that vanish at the
top and have zero derivative at the base. Defining

$$
k_n=\frac{(2n+1)\pi}{2H},\qquad D=\frac{m}{C_o},
$$

the sine coefficients follow by integrating $-q_i d/m$ over the column:

$$
p(d,t)=\frac{q_i}{m}\left[
d-\frac2H\sum_{n=0}^{\infty}
\frac{(-1)^n\sin(k_n d)}{k_n^2}e^{-D k_n^2t}\right].
$$

Integrating this pressure gives the independently predicted stored volume.
Differentiation at the drain gives the outward Darcy flux:

$$
\bar p=\frac{q_i}{m}\left[
\frac H2-\frac{2}{H^2}\sum_{n=0}^{\infty}
\frac{(-1)^n}{k_n^3}e^{-D k_n^2t}\right],
$$

$$
V_s=A H C_o\bar p,\qquad
V_o=Aq_i t-V_s.
$$

$$
q_o=q_i\left[1-\frac2H\sum_{n=0}^{\infty}
\frac{(-1)^n}{k_n}e^{-D k_n^2t}\right].
$$

The runner checks the series with finite differences of the diffusion
equation, the boundary derivatives and independent quadrature of pressure
and outflow. It uses the exact zero state at $t=0$.

## Physical content and the PT coupling lag

For every saved step, a separate mesh-volume integration computes
$Z=\sum_e V_e(S\bar p_e+\alpha\varepsilon_{v,e})$. Each pressure and
mechanical contribution is also reconstructed from consecutive saved
states and compared with the reported last-step balance.

In explicit mechanics and GVS, strain and pressure advance in the same
physical step. Thus $Z$ equals accumulated net boundary inflow to round-off.
An inertial transient may still differ from the quasistatic analytical
solution; conservation alone does not establish equilibrium accuracy.

PT consumes the previous step's pressure increment. Its saved mean total
stress consequently uses $p_{n-1}$, while the saved pressure is $p_n$.
For fixed material coefficients, the quantity conserved by that split is

$$
Z_{\rm split}=Z+\sum_e V_e L\overline{(p_n-p_{n-1})}_e,
\qquad L=\frac{\alpha^2}{B}.
$$

Here $B=K+G/3$ for the recorded in-plane mean stress in 2D and $B=K$ in
3D. The extra term is a coupling lag, not additional physical water.
The checkpoint stores $p_{n-1}-p_n$ as the raw pending increment.
Its contribution is shown separately; reporting only the split budget
would hide the finite physical-content error.

## Selected result and numerical controls

![Pressure and water budget on the selected mesh](../../benchmarks/poromechanics/results/runs/fluid-content-20261108/spatial-pt/fig/result.svg)

The main figure uses $h=0.05$ m and a 0.005 s PT step through 4 s.
The maximum pressure error during 1–4 s is **0.0286%** of the 100 Pa
steady basal pressure. Integrated Darcy flux is measured on actual top
facets, separately from the boundary reaction used by the solver.

| Spacing, m | Pressure error in 1–4 s / 100 Pa | Final integrated Darcy outflow error / injected volume |
| --- | --- | --- |
| 0.20 | 0.649% | 0.165% |
| 0.10 | 0.134% | 0.0296% |
| 0.05 | 0.0286% | 0.00975% |

At 1 s on the 0.05 m mesh, reducing the timestep separates the PT lag
from a water-budget defect:

| PT step, s | Physical content minus net inflow / injected volume | Budget including pending increment |
| --- | --- | --- |
| 0.020 | −1.055% | round-off |
| 0.010 | −0.523% | round-off |
| 0.005 | −0.260% | round-off |

The lag falls at first order. At 4 s the finest case's pending term is
0.00232% of injected water. Physical storage and reference drainage can
show error cancellation; the table therefore retains the independent
strong-flux measurement rather than selecting the smallest single error.

Short ordinary explicit controls cover inertial scalings $10^4$, $10^5$
and $10^6$. An additional 2D GVS case uses $q_i=10^{-3}$ m/s to exercise
adaptive mass and timestep changes: its step ranges from 0.000610 s down
to 0.0000187 s. Its physical-content budget closes within
$2.4\times10^{-15}$ of injected volume. Its early pressure discrepancy
is retained as an inertial/discretization limit, not presented as a
quasistatic accuracy result. Tiny native/HDF5 cells additionally check the
2D/3D stress conventions, and eight coupled restart controls retain all
physical save/checkpoint payloads exactly.

This qualification uses constant porosity and compressibility, saturated
flow, an elastic skeleton and fixed geometry. Evolving storage,
unsaturated content, local conservative remapping and the corresponding
enthalpy budget remain separate requirements.

The [yielding-column budget](plastic-content.md) extends the constant-property
check to plastic and EP-RSF material paths and open drainage. It identifies
the correction needed when an elastic mean-stress conversion is used after
a plastic return, while retaining the same pending-pressure stabilization.

## Reproduce

The frozen native/HDF5 solvers and their source/binary hashes are bound by
`restart-timestep-20261108/qualification.json` (implementation `50c4001`).
Use a new output directory for every simulation:

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
solver=benchmarks/poromechanics/results/runs/restart-timestep-20261108/candidate/dynearthsol2d
python3 benchmarks/poromechanics/fluid_content.py "$solver" /new/content-spatial \
  --spacings .2 .1 .05 --dt .005 --duration 4 --figures
python3 benchmarks/poromechanics/fluid_content.py "$solver" /new/content-time \
  --spacings .05 --dt .01 --duration 1
python3 benchmarks/poromechanics/fluid_content.py "$solver" /new/content-gvs \
  --method gvs-bulk --spacings .2 --duration .1 --inflow .001
python3 benchmarks/poromechanics/fluid_content.py "$solver" /new/content-spatial \
  --replay --figures
```

The retained studies, input files, every-step CSV tables, all-frame hashes
and numerical limits are in
[`fluid-content-20261108`](../../benchmarks/poromechanics/results/runs/fluid-content-20261108/).
