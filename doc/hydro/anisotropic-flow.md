# Anisotropic flow with a prescribed permeability change

The [selected problem page](../../benchmarks/poromechanics/results/review/anisotropic-time.html)
contains the setup, analytical solution and the finest pressure comparison.
The test uses a fixed, drained 1 m square. It changes scalar permeability
and directional anisotropy at 0.5 s, then checks a second restart within
the new phase. Material parameters are prescribed; none are fitted.

The pressure equation uses physical strain-form storage, with strain
feedback into pressure disabled. Pressure still loads the solid, so the
physical timestep also respects the elastic limit. This isolates transport
from deformation feedback without changing the production solver.

## Analytical reference

The two initial modes have relative amplitudes 1 and 0.35. Let
$\varphi_n=\sin(n\pi x)\sin[\pi(z+1)]$ for $n=1,2$, with lengths divided
by the 1 m side. All four faces have zero pressure. The equation and storage are

$$
S_0\partial_t p=\mu^{-1}\nabla\cdot(\mathbf{k}\nabla p),
\qquad S_0=\phi/K_f+(\alpha-\phi)/K_s.
$$

Substitution gives independent decay rates

$$
\lambda_{n,j}=\pi^2(n^2k_{x,j}+k_{z,j})/(\mu S_0).
$$

Write $\tau_1=\min(t,t_s)$ and $\tau_2=\max(0,t-t_s)$, with $t_s=0.5$ s.
The continuous solution through the jump is

$$
p=P_0\sum_{n=1}^{2}a_n(0)\exp(-\lambda_{n,1}\tau_1-\lambda_{n,2}\tau_2)\varphi_n.
$$

Here $P_0=1000$ Pa, $k_0=10^{-15}$ m²,
$\mathbf{k}_1=k_0\operatorname{diag}(4,1)$ and
$\mathbf{k}_2=k_0\operatorname{diag}(2,18)$. The implementation sets the
second scalar permeability to $2k_0$ and its dimensionless tensor to
$\operatorname{diag}(1,9)$. The two horizontal wave numbers distinguish
the directional coefficients; their exchange cannot preserve both decays.
The material uses $\phi=0.3$, $\alpha=1$, $K_f=2.17$ GPa, $K_s=37$ GPa
and $\mu=1.002\times10^{-3}$ Pa s.

## Measured accuracy and scope

At a fixed 0.00005 s step, maximum pressure error over the saved history
falls from 1.4644% to 0.36409% to 0.090124% of $P_0$ on 0.1, 0.05 and
0.025 m meshes. The independent time study uses the finest mesh with
20,000, 40,000 and 80,000 steps. Successive final-field differences are
0.0059876 and 0.0029937 Pa: their ratio is 2.00004, consistent with first
order in time. The selected 80,000-step result has 0.089961% pressure error.

All physical state and output-clock payloads match exactly after the
second-phase restart. These checks establish directional diffusion and
operator reconstruction after a prescribed coefficient change. Continuous
state-dependent permeability, rotated tensors changing in time, deformation
feedback and general fluid-content budgets remain separate checks. No 3D
qualification is inferred from this 2D study.

## Reproduction

Use a fresh output directory for each command. The manifest identifies the
frozen solver separately from the analysis checkout.

```sh
ulimit -v 6291456
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
B=benchmarks/poromechanics/results/runs/inlet-pulse-20261108
python3 tests/functional/hydromechanics/check_anisotropic_flow.py \
  "$B/candidate/dynearthsol2d" /tmp/des-anisotropic-mesh \
  --binary-manifest "$B/qualification.json"
python3 tests/functional/hydromechanics/check_anisotropic_flow.py \
  "$B/candidate/dynearthsol2d" /tmp/des-anisotropic-steps40k \
  --spacings .025 --steps 40000 --binary-manifest "$B/qualification.json"
python3 tests/functional/hydromechanics/check_anisotropic_flow.py \
  "$B/candidate/dynearthsol2d" /tmp/des-anisotropic-steps80k \
  --spacings .025 --steps 80000 --binary-manifest "$B/qualification.json"
```

The portable problem input uses an included analytical initial-state
checkpoint. In its data directory, run `input.cfg` and then
`continuation.cfg` with the recorded solver to reproduce both phases.
Evidence is under `results/runs/transport-qualification-20261108/anisotropic-time*`.
