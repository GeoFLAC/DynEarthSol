# Averaged aquifer reference and DES storage

`benchmarks/poromechanics/averaged_aquifer_reference.py` implements a
thickness-averaged Biot reference for infinite and finite-annulus aquifers.
Its independent checks are in
`tests/functional/hydromechanics/check_averaged_aquifer_reference.py`.
The assumptions and paper conventions below define what can be compared
with the slab and regional DES studies; they do not imply that a general
three-dimensional or phreatic aquifer has been qualified.

## Conventions

Stress and strain are tension-positive. The excess pore pressure `p` is
positive in compression, so extraction makes it negative. With Biot
coefficient `alpha`, total stress is

    sigma = lambda eps_v I + 2 G eps - alpha p I

Fluid content is `zeta = alpha eps_v + p/M`, with
`1/M = phi/K_f + (alpha - phi)/K_s`. Mass balance is
`d zeta/dt = (k/mu) lap p` minus the extraction. The well rate `Q` (m^3/s,
total) is positive for extraction. Subsidence is positive downward. `U_r` is
positive outward.

## Thickness-averaged model

The model is that of Bear and Corapcioglu (1981b), with `alpha`, `K_s` and
fluid compressibility kept. It is a confined, homogeneous, isotropic layer
of constant thickness `B`. The well is fully penetrating, and the top and
bottom are sealed, so `p = p(r, t)`. The model makes two mechanical
assumptions:

- **A1, plane incremental total stress:** `sigma_zz`, `sigma_xz` and
  `sigma_yz` do not change.
- **A2:** horizontal displacement is uniform over the thickness.

From these assumptions:

1. `sigma_zz = 0` gives `eps_zz = (alpha p - lambda e_h)/(lambda + 2G)`. Here
   `e_h = (1/r) d(r U_r)/dr` is the areal strain.
2. The horizontal stresses are plane stress, with
   `lambda* = 2G lambda/(lambda+2G)` and `alpha* = 2G alpha/(lambda+2G)`.
3. Radial equilibrium of an irrotational field gives
   `(lambda*+2G) e_h - alpha* p = (lambda*+2G) h(t)`, uniform in `r`:

       e_h = alpha p / (2(lambda+G)) + h(t)
       eps_v = alpha p / (lambda+G) + 2G h / (lambda+2G)
       r U_r = c_in + alpha/(2(lambda+G)) * int_rw^r rho p drho + h (r^2 - rw^2)/2

4. Flow:

       S_s dp/dt + c dh/dt = (k/mu) lap p
       S_s = 1/M + alpha^2/(lambda+G)
       c   = 2G alpha/(lambda+2G)

In an infinite domain with `U_r -> 0`, `h = 0`. The pressure is then Theis
with storativity `rho g B S_s`, and the vertical strain equals the areal
strain. This is the paper's `alpha'' = 1/(lambda+G)` (eq. 94) and
`Delta z/B = eps/2` (eq. 103).

`mechanics='uniaxial'` is the vertical-compressibility model of Jacob and
Bear and Corapcioglu (1981a): `e_h = 0`, `eps_zz = alpha p/(lambda+2G)` and
`S_s = 1/M + alpha^2/(lambda+2G)`. For the same pressure change, plane-stress
subsidence is `(1 - nu)` times the uniaxial value. If both models are
calibrated to the same storativity, it is one half, which is the paper's
"approximately half".

## Finite annulus solution

`AveragedAquifer.finite_well` solves the same equations exactly on
`rw <= r <= R`.

- **Hydraulic conditions:**
  - At the well face, the total rate is `2 pi rw B (k/mu) dp/dr = Q`.
  - At the outer radius, the pressure is fixed (`p = 0`) or there is no flow.
- **Mechanical conditions:** at each radius, either `U_r = 0` or
  `sigma_rr = (lambda*+2G) h - 2G U_r/r = 0`.
- **Solution method:** in the Laplace domain,
  `p^ = A I0(qr) + B K0(qr) - (c/S_s) h^`, with `q^2 = s S_s mu/k`. `A`, `B`
  and `h^` come from a 3x3 linear system. Exponentially scaled Bessel bases
  keep it well conditioned. The result is inverted with the fixed Talbot
  method (Abate and Valko 2004), using 32 nodes.

`h(t)` is the integration function `g(t)` that the paper sets to zero for an
infinite domain (eqs. 89-92). With a fixed outer boundary, it is nonzero from
the first instant and grows with the pumped volume. While the pressure front
is far from `R`, it is approximately

    h ~ alpha Q t / (2 pi (lambda+G) B R^2 (1/M + alpha^2/(lambda+2G)))

The approximation is accurate to 1-5% for `sqrt(4Dt)/R < 0.2`. `h` feeds back
into the pressure through `c dh/dt`, so a finite-domain pressure is not Theis.
Size a DES domain from `finite_well(...)['h']`, and do not assume that a
distant boundary is harmless.

The checks compare the solution with each of the following:

- the printed closed forms, in the limit `alpha=1`, `K_f=K_s=inf`,
  identical to rounding
- the paper's example, with Fig. 3's peak `U_r` reproduced
- an independent radial finite-volume solver, which matches to 1e-4
- the sealed fluid balance `int zeta dV = -Q t`, which holds to 1e-5
- the boundary conditions and the plane-stress constraint
- equilibrium and kinematics, by finite differences
- the large-`R` limit of the line well
- the steady Thiem profile

Run them with:

```sh
python3 tests/functional/hydromechanics/check_averaged_aquifer_reference.py -v
```

## What the reference does not remove

The rederivation removes the approximations the paper makes in solving its
own equations:

| Approximation | Paper location | Removed? |
|---|---|---|
| line sink, `x K1(x) -> 1` | Part 2 (98)-(99) | yes, finite `rw` |
| fluid compressibility dropped in the strain amplitude, `alpha'' >> n beta` | (101)-(104) | yes |
| `u_w W(u_w) ~ 0`, `exp(-u_w) ~ 1` in `U_r` | (107) | yes, `U_r(rw)=0` exactly |
| `g(t) = 0`, infinite domain | (89)-(92) | yes, `h(t)` solved |
| `alpha = 1`, incompressible grains | throughout | yes |

The following are part of the model definition, so a rederivation cannot
remove them:

- **A1 and A2 are exact only in the thin-layer limit.** In a 3D slab with a
  traction-free top and a roller base, `eps_zz(r)` varies, which creates
  `sigma_rz ~ G z d(eps_zz)/dr`. The plane-stress result is then accurate to
  `O((B/L)^2)`, where `L` is the radial length scale of `p`. Compare
  thickness-averaged DES quantities at `r >> B` and `sqrt(4Dt) >> B`, with a
  `B/r` sweep that shows the gap shrinking. Overburden stiffness, or a base
  fixed in all directions, violates A1.
- **Near-well, early-time coupling.** Full 3D Biot shows a transient pressure
  rise near the well (the Noordbergum effect, Verruijt 1969, Hsieh 1996). The
  averaged model cannot show it. This is a modelling difference, not a DES
  error.
- **Layering and aquitard delay.** Part 1 lumps aquifers and aquitards with a
  common head. A layered DES model with resolved aquitards is a different
  problem.
- **2D slices.** A 2D vertical DES slice cannot represent radial plane flow
  together with subsidence.

## Paper reproduction and errata

Notes from a full reading of all three PDFs, with equations checked against
page images:

- **Part 1 (WRR 17(4) 937):**
  - Eq. (83) reproduces Fig. 2 to within 0.5 cm, using the stated
    `C_v = 0.6 m^2/s`.
  - Table 2 with `rho g = 9810 N/m^3` gives `T = 9.50e-3 m^2/s`,
    `sum alpha' B = 1.725e-4 m^3/N` and `C_v = 56 cm^2/s`, not `6e3`. The
    stated value requires `rho g ~ 100 N/m^3`, the figure printed next to
    eq. (44). Reproduce the figure with the stated `C_v`, but do not use the
    example as physical data.
  - The signs in (42), (63) against (62), and (71) against (80)-(81) are
    inconsistent. A self-consistent form is
    `d delta/dt = div(C_v grad delta) + Q` (with `S_v = 1`) and
    `delta = rho g alpha' B s`.
- **Part 2 (WRR 17(4) 947):**
  - Fig. 3 is reproduced: `U_r,max = 17.06 cm` at `8.545 km`, against about
    17.1 cm at 8.564 km in the paper.
  - The peak satisfies `W(u) = (1 - e^-u)/u`, so `u = 0.3214`. The paper
    gives 0.323.
  - A quoted ratio, 51093, is about 2% off its own eq. (107), which gives
    50076. The figure follows the formula.
- **Part 3 (WRR 19(4) 895, phreatic; not implemented here):**
  - Eqs. (80), (85) and (89) reproduce the numbers in Figs. 4-8.
  - The printed `u = 0.232` should be 0.3214. The printed 1.1367 corresponds
    to 0.323.
  - Eq. (74) adds the two consolidation coefficients. Dividing (72) by `T`
    combines them harmonically instead: `C_v* = C_v/(1 + S_y/2) = 5000 cm^2/s`,
    not 36000. As a result, `C_v**` changes from 235.9 to 226.7 cm^2/s.
  - The overburden decrease `S_y gamma_w s` is entered with a
    compression-positive sign inside the tension-positive
    `sigma = sigma* - p I`. Physically, a water-table drop raises the
    effective vertical stress by `(1 - S_y) gamma_w s`, but the paper uses
    `(1 + S_y)`. This affects every `S_y`-dependent mechanical term in
    (68)-(89). Rederive it independently before any phreatic benchmark uses
    it.

## DES pressure storage

The linear saturated Biot capacity is

$$
S_0=\frac{\phi}{K_f}+\frac{\alpha-\phi}{K_s}.
$$

The fixed-stress split adds its mechanical stabilization capacity
$L=\alpha^2/K_b$. The strain-form explicit update and the fixed-stress PT
update use their own consistent source terms and timestep limits; see
[transport](transport.md). These coefficients must not be interchanged
when matching a pumping reference. The optional finite-density inventory
model has separate [constitutive equations and limits](mass.md).

## References

- Bear, J. and Corapcioglu, M. Y. (1981a, 1981b), WRR 17(4), 937-946 and
  947-958, doi:10.1029/WR017i004p00937, doi:10.1029/WR017i004p00947.
- Corapcioglu, M. Y. and Bear, J. (1983), WRR 19(4), 895-908,
  doi:10.1029/WR019i004p00895.
- Abate, J. and Valko, P. P. (2004), Int. J. Numer. Meth. Eng. 60, 979-993.
- Verruijt, A. (1969), in Flow Through Porous Media, ed. R. J. M. De Wiest;
  Hsieh, P. A. (1996), Ground Water 34(6), 1082-1089.
