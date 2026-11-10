# Water balance while a confined material yields

A saturated column is compressed at a prescribed top velocity. The
fluid stays inside the sealed column while the skeleton yields. Its mean
pressure follows from water conservation independently of the plastic
stress return, providing a direct check of their coupling.

![Sealed compression cell](../../benchmarks/poromechanics/results/runs/plastic-source-20261108/fig/setup.svg)

The height and width are 1 m. The base and side-normal motion are fixed;
the top moves down at $v=10^{-3}$ m/s for 1 s. Gravity is zero and the
mesh stays fixed. The 2D cell is plane strain with unit thickness.
The baseline has $K=20$ MPa, $G=10$ MPa, $\alpha=1$, $\phi_n=0.25$,
$K_f=1$ MPa and effectively incompressible grains ($K_s=10^{20}$ Pa).
Mohr–Coulomb cohesion is 100 Pa, friction is 20° and dilation is varied
between −10°, 0° and 10°. Porosity and pressure capacity stay fixed.

## Independent pressure and water-volume reference

With tension-positive strain, the imposed total volume change is
$\int\varepsilon_v\,dV=-Avt$. Linear Biot content is

$$
\zeta=\alpha\varepsilon_v+S p,\qquad
S=\frac{\phi_n}{K_f}+\frac{\alpha-\phi_n}{K_s}.
$$

For the sealed cell, integrating mass balance therefore gives

$$
\bar p(t)=\frac{\alpha vt}{HS}.
$$

This relation holds through yielding: it uses total volume change,
not an elastic stress-to-strain conversion. Ordinary explicit mechanics
uses the current projected volumetric strain before updating stress and
matches the global content balance directly.

The one-step PT split retains the pending pressure contribution described
in the [deforming-column budget](content.md). Its conserved quantity is

$$
Z_{\rm split}=\int\zeta\,dV+VL(\bar p_n-\bar p_{n-1}),
\qquad L=\frac{\alpha^2}{B}.
$$

Here $B=K+G/3$ for the recorded 2D mean stress and $B=K$ in 3D.
The independent discrete pressure reference follows immediately:

$$
\bar p_n=\frac{\alpha vt_n/H+L\bar p_{n-1}}{S+L}.
$$

The physical continuous reference reaches 4 kPa at 1 s. With a 0.01 s
PT step, the remaining pending-pressure lag is about 0.17% of that value.
The budget including the pending contribution closes to round-off.

## Correct source after a plastic return

An elastic mean-stress increment can be converted to a volume increment.
A plastic return adds stress relaxation, and in plane strain includes an
out-of-plane plastic contribution. Treating that stress change as elastic
volume change creates or removes fluid in the pressure update.

The corrected fixed-stress plastic path uses the accepted projected
volumetric strain rate. Its stabilizing pressure capacity stays $S+L$;
the source retains the previous raw pressure increment $D$:

$$
(S+L)\Delta p+\alpha\Delta\varepsilon_v-LD
=\text{net flow per bulk volume}.
$$

With heating, the left side also contains $-\beta\Delta T$, where
$\beta=\phi_n\alpha_f+(\alpha-\phi_n)\alpha_s$.
Pressure increments are frozen before the pending field is cleared and
remain fixed across hydraulic substeps. They are increments, not rates
rescaled from the previous timestep. The accepted mechanical strain
comes from the existing PT trial/commit solver; no new PT owner,
constitutive state or parameter is introduced.

Purely elastic rheology keeps its existing mean-stress path. Saturated
and partly saturated elastic regression runs preserve 10,120 physical
payloads exactly in native/HDF5 output and both dimensions.

## Selected result and controls

![Pressure and volume balance after yielding](../../benchmarks/poromechanics/results/runs/plastic-source-20261108/fig/result.svg)

The main figure shows the corrected 2D PT case with 10° dilation. The
second panel includes the pending increment explicitly; it is not counted
as physical water. Every saved pressure and strain is integrated on the
actual mesh. Both the global boundary-volume change and the pressure
recurrence are checked independently of the reported pressure residual.

The parent defect is reproduced in mechanically converged cells:

| Cell | Parent budget drift / imposed volume | Corrected budget |
| --- | --- | --- |
| 2D, 10° dilation | +1.238% | round-off |
| 2D, zero dilation | −3.105% | round-off |
| 3D, 10° dilation | +4.489% | round-off |

Other controls cover negative dilation, ordinary explicit mechanics,
internal Darcy flow and open drainage. The finite-permeability controls
take 17–39 hydraulic substeps per physical step. A 0.01-to-0.02 s PT
restart retains the pressure-increment balance; unchanged restart
histories remain exact. In the 2D drained PT control, 0.992116 L per metre
thickness leaves by 1 s; pressure storage plus pending contribution
accounts for the rest of the imposed 1 L pore-volume loss.
This checks the boundary reaction and material content, rather than an
independent spatial Darcy-flux reference.

Older plastic checkpoints load with their saved physical fields unchanged.
The corrected continuation preserves any fluid-content offset accumulated
by the older source; it does not reconstruct fluid lost before the restart.
Final validation therefore uses fresh runs with the corrected source.

A rerun of the separate pressure-driven yielding benchmark on the
0.05 m mesh with 1600 steps gives 0.559% maximum and 0.347% mean onset
error, consistent with its recorded convergence trend. Both parent/candidate
core comparisons are bit-exact; the 2D/3D component suites pass 7,579 and
11,171 checks, respectively.

## Coupled retention, RSF and thermal controls

The homogeneous Gardner control starts at −5 kPa suction and yields
while wetting under compression. With $S_e=e^{ap}$, its independent scalar
content inversion solves, for each step,

$$
(S+L)\Delta p+\phi_n\Delta S_e
+\alpha S_{e,\rm old}\Delta\varepsilon_v-LD=0.
$$

This verifies the implemented nonlinear retention-rate coupling and raw
pressure buffer. It does not establish a finite porosity/deformation
mass law or conservative remapping.

The EP-RSF control has evolving pressure and friction state, rather than
a prescribed constant pore pressure. For its constant affine slip rate,
the aging reference is

$$
\theta_{n+1}=\theta_n+(1-V\theta_n/D_c)\Delta t,
$$

$$
\mu=\mu_0+a\log(V/V_0)+b\log(V_0\theta/D_c).
$$

The coefficients are $a=0.01$, $b=0.015$, $D_c=0.01$ m and
$V_0=3\times10^{-4}$ m/s. Projection option 1 uses
$V=2w\dot\varepsilon_{\rm II}$; the reference computes the minimum
simplex altitude $w$ independently from the mesh.

An independent two-face Mohr–Coulomb reference checks the effective
stress, including the hidden 2D out-of-plane component. At the compression
edge, let $\sigma_A$ be axial stress and $\sigma_L$ the two equal lateral
stresses. For zero dilation, the sum of the two flow multipliers is

$$
N_\varphi=\frac{1+\sin\varphi}{1-\sin\varphi},
\qquad \varphi=\arctan\mu.
$$

$$
\Delta\gamma=\frac{[N_\varphi\sigma_L^t-\sigma_A^t
-2c\sqrt{N_\varphi}]_+}{G(2+N_\varphi)}.
$$

$$
\sigma_A=\sigma_A^t+2G\Delta\gamma,\qquad
\sigma_L=\sigma_L^t-G\Delta\gamma.
$$

Trial stresses come from the imposed elastic strain increment. Total
stress then subtracts the pressure already consumed by mechanics.
Friction state and coefficient match their discrete reference exactly;
the stress error is about $1.1\times10^{-15}$ of $K v t/H$.

The thermal control combines imposed compression and supplied heating,
nonzero solid expansion and a plastic return. Its pressure, stress and
fluid-content references agree within $2.2\times10^{-14}$ of their scales.
This checks the full thermal fluid source on the plastic path. It adds
no shear-heating law or open-system enthalpy qualification.

Finite-grain plastic porosity, permeability enhancement/recovery,
conservative remapping and general open-system energy remain separate
requirements. A correct last-step pressure residual alone does not
establish any of them.

## Application replay

The same correction is checked in the saturated pressure-driven slope,
with gravity, open pressure boundaries and localized plastic strain.
Both 0.25 m soils retain their pressure at basal strain 0.001: shifts are
−1.60 Pa around 16.54 kPa (35° friction) and +0.014 Pa around 11.49 kPa
(30°). Matched-window physical-content-plus-pending defects fall from
14.3 and 40.4 mL per metre thickness to about $1.7\times10^{-17}$ m³.
The [application audit](slope.md#current-source-coupling-audit) states
the windows, boundary-volume normalization and numerical stopping limits.

The 0.5 m rainfall replay retains its 8.72433 h, 12.32327 kPa fixed-slip
point to round-off. Its accepted surface pressures remain nonpositive.
This confirms the selected application trend; it does not establish a
general finite-deformation unsaturated content law.

## Plastic pore recovery with finite grain stiffness

The additional property check uses $\alpha=0.8$, $K_s=100$ MPa and
$K_f=100$ MPa. It retains DES's plastic effective stress
$\boldsymbol\sigma+\alpha p\mathbf I$ and its matching total-strain fluid
source. A separate plastic Biot coefficient would define a different
constitutive model; it cannot be introduced by changing porosity alone.
The distinction between elastic and plastic porosity coefficients is
described in [Hu (2009), sections 1.3.2.2 and 5.1.3](https://pepite-depot.univ-lille.fr/LIBRE/Autres/2009/50376-2009-Hu-2.pdf).

![Drained opening and holding cell](../../benchmarks/poromechanics/results/runs/finish-20261109/healing/qualified/fig/setup.svg)

The selected cell opens at $10^{-3}$ m/s for 1 s, then holds its strain.
All vertices have $p=0$; temperature is 293 K. The 2D cell is plane strain.
The hold raises cohesion to prevent further yielding and changes the
physical step from 0.01 to 0.025 s. This isolates recovery of the accepted
plastic opening. Compression, negative/zero dilation and a purely elastic
opening are separate controls.

Recovery uses signed plastic volume, rather than the scalar plastic-strain
invariant. In 3D and plane strain it is

$$
\Delta\varepsilon_v^p=\Delta\varepsilon_v
-\frac{\operatorname{tr}_3\Delta\boldsymbol\sigma_{\rm elastic}}{3K}.
$$

The stress change excludes reversible pressure and thermal increments;
plane strain includes the hidden out-of-plane stress. Only the accepted
constitutive branch supplies this increment. A remaining recoverable
dilation $\eta\ge0$ has a recovery time $\tau(T)$. At each physical step,

$$
h=\eta_n(1-e^{-\Delta t/\tau}),\qquad
\eta_{n+1}=\max(0,\eta_n-h+\Delta\varepsilon_v^p).
$$

Recovery acts before adding the new plastic increment. Compaction reduces
the remaining opening and is never itself healed. Elastic volume changes
do not enter $\eta$. The existing porosity history includes the closure:

$$
\Delta\xi=\Delta\varepsilon_v+\Delta p/K_s-\alpha_s\Delta T-h,
\qquad \phi=\alpha-(\alpha-\phi_0)e^{-\xi}.
$$

At frozen pressure, strain and temperature, closure removes pore volume
$\delta\phi_c=(\alpha-\phi_n)(e^h-1)$ per bulk volume. The hydraulic
source includes its negative increment. Thus closure expels fluid through
a drain or raises pressure in a sealed cell; changing permeability alone
would omit this coupling. This source is qualified in the existing
tangent Biot update, not as a finite-density mass law.

Set `mat.pore_healing_time` to a positive reference time in seconds;
zero, the default, disables recovery. The material temperature law is

$$
\tau(T)=\tau_0\frac{T}{T_0}
\exp\!\left[\frac{E_a}{R}(T^{-1}-T_0^{-1})\right].
$$

`mat.pore_healing_temperature` gives $T_0$ in K and
`mat.pore_healing_energy` gives $E_a$ in J/mol, with $R=8.3144$ J/(mol K).
Mixed elements average rates by marker count. The temperature dependence
and porosity power law follow the modeling choices discussed by
[Yang et al. (2021)](https://doi.org/10.1029/2021JB021787); the DES recovery
state and coupling are specified above and do not reproduce their complete
fault model. No porosity clamp or unstable Euler decay is applied.

For the hold, let $s$ be time since loading stopped. The analytical solution
uses the accepted loaded state as its initial condition, without fitting:

$$
\eta(s)=\eta_*e^{-s/\tau},\qquad
\phi(s)=\alpha-(\alpha-\phi_*)e^{\eta_*[1-e^{-s/\tau}]}.
$$

Together with $k=k_0(\phi/\phi_0)^3$, this supplies the plotted reference.

![Plastic opening recovers during the hold](../../benchmarks/poromechanics/results/runs/finish-20261109/healing/qualified/fig/result.svg)

The selected hold has porosity error $1.12\times10^{-16}$ and permeability
error $1.35\times10^{-15}$ relative to $k_0$. The compression controls also
check plastic volume against an independent two-face Mohr–Coulomb flow
rule, including both signs of dilation. The elastic hold preserves its
porosity exactly. Input, checkpoint, timestep-change and actual-remesh
checks are recorded with the dimensional/container matrix in
[`healing/qualified`](../../benchmarks/poromechanics/results/runs/finish-20261109/healing/qualified/).

PT trials freeze recovery and restore the candidate plastic increment;
`update_porosity` commits $\eta$ once after acceptance. Checkpoints save
the recoverable opening and the signed increment. Enabling recovery from
an older checkpoint without this history is rejected. Permeability remains
derived by `MatProps`; no second porosity or permeability array is added.

This validates the chosen property law and tangent hydraulic coupling.
Finite-density evolving-content remapping, recovery energy and the general
retained-fluid heat balance remain open. It is a material check rather
than a calibrated earthquake-cycle application.

## Reproduce

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
solver=benchmarks/poromechanics/results/runs/plastic-source-20261108/candidate/dynearthsol2d
python3 tests/functional/hydromechanics/check_plastic_content.py "$solver" /new/plastic-content
python3 tests/functional/hydromechanics/check_plastic_content.py --replay /new/plastic-content
python3 benchmarks/poromechanics/results/runs/plastic-source-20261108/render.py \
  .
```

Repeat with the corresponding native/HDF5 frozen 3D solver and
`--ndims 3` for the tiny dimensional controls. Inputs, every-frame tables,
source and binary hashes are retained in
[`plastic-source-20261108`](../../benchmarks/poromechanics/results/runs/plastic-source-20261108/).

For pore recovery, build the chosen dimension/container, then run:

```sh
python3 tests/functional/hydromechanics/check_recovery.py ./dynearthsol2d /new/recovery
# Plotting instructions and scripts are in the DES manual.
```

The checker writes exact inputs and the summary; the plotter exports the
hold table, schematic and PNG/SVG/PDF figures with source/file hashes.
