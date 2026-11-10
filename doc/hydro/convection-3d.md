# Three-dimensional porous convection cells

A unit cube heated from below tests a Darcy–Boussinesq temperature mode
that varies in all three directions. All faces are impermeable. The top
and base hold 273 and 283 K; vertical walls are insulated. Mechanics is
stiff, the mesh is fixed and pressure feedback is off. Darcy and thermal
transport keep their explicit substeps.

The initial conductive profile has a 0.1 mK Gaussian disturbance. Its
volume projection onto

$$
\delta T=A(\tau)\cos(\pi x/H)\cos(\pi y/H)\sin[\pi(z/H+1)]
$$

measures the $(1,1,1)$ mode. A separate unseeded cube checks conduction
at rest. The selected result uses 35,937 nodes and 163,840 tetrahedra.

## Analytical reference

Let $\kappa_T=\lambda/C$, $\tau=\kappa_Tt/H^2$ and
$D=D_p/\kappa_T$, where $D_p=k/(\mu S_B)$. The test uses the declared
fixed-stress storage

$$
S_B=\phi/K_f+(\alpha-\phi)/K_s+\alpha^2/K.
$$

The fluid bulk modulus is deliberately scaled to 10 Pa to make this
verification economical. The reference retains that storage; replacing
it with an incompressible solution would test a different transient.
It is not a realistic water-compressibility application.

Scale the pressure-mode amplitude by $\rho_fg\alpha_fH$, giving $P$ in
temperature units. Linear pressure diffusion and heat advection give

$$
\dot P=-DK_m^2P-D\pi A,
$$

$$
\dot A=Ra\pi P+(Ra-K_m^2)A,
$$

where dots mean derivatives with respect to $\tau$, $a^2=2\pi^2$ and
$K_m^2=3\pi^2$. The first equation follows from the vertical derivative
of the temperature perturbation in the divergence of Darcy flux. The
second combines thermal diffusion with advection of the conductive
temperature gradient. The characteristic polynomial is

$$
\lambda_m^2+b_1\lambda_m+b_0=0,
$$

$$
b_1=(D+1)K_m^2-Ra,\qquad b_0=D(K_m^4-Ra\,a^2).
$$

Its larger real root gives the late growth rate
$\sigma=\kappa_T\lambda_+/H^2$. Zero growth gives
$Ra_c=K_m^4/a^2=9\pi^2/2$. As $D\to\infty$, elimination of $P$
recovers $\lambda_+=Ra\,a^2/K_m^2-K_m^2$. These limits, neutral
stability and zero-Rayleigh diffusion are checked independently.

## Results

All three meshes use the same properties, physical timestep of 500 s and
200,000 s duration. Late rates are measured over dimensionless times
0.1–0.2; the reference slopes have no fitted material parameter.

| Mesh spacing | Rate error at $Ra/Ra_c=0.7$ | At 1.2 | At 1.4 |
| --- | --- | --- | --- |
| 1/8 m | 9.19% | 36.79% | 23.33% |
| 1/16 m | 2.43% | 9.65% | 6.11% |
| 1/32 m | 0.62% | 2.44% | 1.55% |

The spatial sequence approaches second order. Halving the physical
timestep to 250 s changes the selected $Ra/Ra_c=1.4$ rate by 0.0067%.
The unseeded finest cube retains the conductive temperature exactly.
All eleven histories were independently replayed from raw fields and
their file hashes checked. Earlier 1,000 Pa runs remain separate evidence;
they are not mixed into this same-property convergence table.

## Reproduce

```bash
python3 tests/functional/hydromechanics/check_convection_3d.py \
  dynearthsol3d OUT --spacing .03125 --ratio 1.4 \
  --fluid-bulk 10 --steps 400 --dt 500
```

Use a new output folder for each mesh and Rayleigh ratio. Add `--seed 0`
for the rest control. The temporal control uses `--steps 800 --dt 250`.
The review page provides the cube schematic, selected input, figures,
rate table and source/binary/reference hashes. It qualifies linear 3D
growth and decay; nonlinear steady cells and deformable convection are
outside this comparison.
