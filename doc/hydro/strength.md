---
title: Pore pressure and material strength
sidebar_label: Effective-stress strength
---

# Pore pressure and material strength

This material test checks whether pore pressure changes Mohr–Coulomb
strength through effective stress. It isolates the constitutive update
before a slope or injection calculation introduces flow, geometry and
equilibrium errors. Every nodal velocity component is prescribed.

![Unit-cell setup](../../benchmarks/poromechanics/results/review/fig/strength/setup.svg)

The unit cell starts at total stress $\boldsymbol\sigma=-C\mathbf I$.
Normal strains are $\epsilon_x=e$ and $\epsilon_z=-e$, with zero shear
strain and, in 3D, zero $\epsilon_y$. Thus total and plastic volume changes
are zero. The 2D calculation uses plane strain. Pressure is uniform and
constant: permeability is zero and pressure feedback is disabled for
this material test. Gravity, mesh motion and thermal diffusion are off.

The material has $K=G=1$ GPa, cohesion $c=100$ kPa, friction angle
$\phi=30^\circ$ and dilation angle zero. There is no strain weakening.
The rate-and-state cases use $a=b=0$: they check the constant-friction
limit while retaining state-variable checkpoint continuation.

## Analytical response

For tension-positive total stress and compression-positive pore pressure,

$$
\boldsymbol\sigma'=\boldsymbol\sigma+q\mathbf I,
\qquad q=\alpha\langle S_e(p)p\rangle_e.
$$

In a saturated material $q=\alpha p$. The extreme effective principal
stresses are $-(C-q)-t$ and $-(C-q)+t$. Substitution into the
Mohr–Coulomb condition gives

$$
t_y=c\cos\phi+(C-q)\sin\phi.
$$

With zero dilation the mean stress stays constant. The elastic and
perfectly plastic portions of this imposed path are therefore

$$
t(e)=\min(2Ge,t_y),\qquad
\epsilon_p=\max(0,e-t_y/(2G)).
$$

All reference paths stay on the compressive shear branch. This test
does not qualify tensile or mixed tensile/shear returns.

DES stores total stress. Its existing return can evaluate this same
skeleton criterion by translating the shear intercept $A=2c\sqrt{N_\phi}$
and tensile cap $T$, with $N_\phi=(1+\sin\phi)/(1-\sin\phi)$:

$$
A_{\rm total}=A-(N_\phi-1)q,\qquad T_{\rm total}=T-q.
$$

This leaves the return algorithm and stored tensor representation intact.
Absolute pressure is needed even when its pending increment is zero.
Initial equilibrium also evaluates the effective-stress strength; it
continues to avoid consuming a pending physical pressure increment.

## Selected result and supporting checks

![Pressurized material path](../../benchmarks/poromechanics/results/review/fig/strength/result.svg)

The selected 3D EP case has $C=1$ MPa and $p=0.8$ MPa. DES follows the
elastic path and yields at **186.603 kPa**, matching the reference. The
previous implementation returned 586.603 kPa because it omitted the
absolute pressure from strength. All saved stress and plastic-strain
values are checked, including the elastic portion and the yield plateau.

| State | Total confinement (kPa) | Pressure (kPa) | Biot coefficient | Reference yield load (kPa) |
| --- | ---: | ---: | ---: | ---: |
| Dry equivalent effective stress | 200 | 0 | 1 | 186.603 |
| Pressurized | 1000 | 800 | 1 | 186.603 |
| Dry equivalent total stress | 1000 | 0 | 1 | 586.603 |
| Smaller Biot coefficient | 1000 | 800 | 0.6 | 346.603 |
| Gardner suction | 200 | -200 | 1 | 200.136 |

The Gardner case uses retention coefficient $10^{-5}$ Pa$^{-1}$ and
residual saturation 0.1, so $S_e=\exp(-2)$ and $q=-27.067$ kPa.
Both EP and constant-friction EP-RSF pass in 2D and 3D, in binary and
native HDF5 builds. The largest normalized stress/plastic-strain error
is below $10^{-13}$. Zero-pressure cases retain 8,424 physical payloads
in 272 files exactly relative to the parent. Current-model restart
before yielding also reproduces the continuation exactly.

## Reproduction and scope

From the repository root, using executables built from the recorded commit:

```bash
OUT=benchmarks/poromechanics/results/runs/strength-check
OMP_NUM_THREADS=1 python3 tests/functional/hydromechanics/check_pore_strength.py ./dynearthsol2d "$OUT/2d" --ndims 2
OMP_NUM_THREADS=1 python3 tests/functional/hydromechanics/check_pore_strength.py ./dynearthsol3d "$OUT/3d" --ndims 3
```

`runs/ep-pressure-20261108/qualification.json` records the source/executable
hashes and measured acceptance. Each study retains inputs, frames and
its reference comparisons. The setup and result also have PNG/PDF exports.

This is a fixed-pressure, zero-dilation constitutive check. Pressure-driven
failure with flow, evolving RSF friction/state accuracy, plastic dilation,
mixed-material suction, tensile returns and the full slope/injection
applications remain separate requirements. Changed strength also changes
the interpretation of old plastic histories: continue unfinished legacy
plastic runs with their recorded source and executable. General
constitutive-checkpoint migration is still open.
