# Rate-and-state friction under a velocity step

The [material-strength page](strength.md) checks EP-RSF only at constant
friction ($a=b=0$). This checks the evolving law: a pressurized cell sheared
at yield sees its loading rate rise tenfold.

## Problem

The strength unit cell (every node prescribed; normal extension rates
$\pm\dot e$; plane strain; initial stress $-C\mathbf I$, $C=1$ MPa; fixed
pore pressure 0.8 MPa, $\alpha=1$; cohesion 100 kPa; no dilation), with
EP-RSF and the aging law: $a=0.01$, $b=0.015$, $D_c=1$ mm,
$\mu_0=\tan30^\circ$, $V_0=V_1$, and slip rate $V=2w\dot\epsilon_{II}$
(`control.rsf_slip_rate_projection_option = 1`; $w$ is the element's
minimum altitude, $1/\sqrt2$ m here). The cell loads for 64 s at
$\dot e_1=7.8\times10^{-6}$ s$^{-1}$, starting and staying at steady state,
then restarts from its last checkpoint at $10\dot e_1$ for 64 s.

## Exact responses

$$
\mu=\mu_0+a\ln\frac{V}{V_0}+b\ln\frac{V_0\theta}{D_c},\qquad
\theta(t)=\frac{D_c}{V_2}+\left(\frac{D_c}{V_1}-\frac{D_c}{V_2}\right)e^{-V_2(t-t_1)/D_c}.
$$

Friction jumps by $a\ln10=0.0230$ and relaxes over $D_c/V_2\approx9$ s to
a net change of $(a-b)\ln10=-0.0115$. With every node prescribed the
stress is $t=\min(2Ge,\ c\cos\phi+(C-\alpha p)\sin\phi)$, $\tan\phi=\mu$.
The discrete reference repeats DES's update: forward Euler for $\theta$
once per step, then $\mu$, then the stress.

## Results

| Step | DES vs discrete ($\theta$, $\mu$) | Stress vs discrete | $\mu$ vs continuum |
| --- | --- | --- | --- |
| 1 s | 0 | 9e-16 | 9.56e-4 |
| 0.5 s | 0 | 9e-16 | 4.67e-4 |
| 0.25 s | 0 | 9e-16 | 2.31e-4 |
| 0.125 s | 0 | 9e-16 | 1.15e-4 |

DES reproduces its discrete recursion exactly, so under PT the state is
advanced once per physical step, not per iteration, and it survives the
restart. The continuum error is first order in the step, as forward Euler
implies, and 4% of the direct jump at 1 s.

Not covered: the slip law, the maximum-shear projection (option 0),
stick–slip instability (no elastic interaction here) and coupling of
friction to flow or dilation.

## Reproduce

```
python3 tests/functional/hydromechanics/check_rsf_step.py EXE OUT --refinements 1 2 4 8
```

Recorded study: `results/runs/rsf-step-20261108/steps-2d/`. References:
Dieterich (1979, JGR 84, 2161–2168); Ruina (1983, JGR 88, 10359–10370).
