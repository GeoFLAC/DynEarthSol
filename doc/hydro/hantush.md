# Leaky aquifer: Hantush–Jacob pumping

A well pumps a confined aquifer covered by a leaky aquitard whose top is
held at the initial pressure. Water drawn down through the aquitard limits
the drawdown, which tends to a steady cone instead of the ever-growing
Theis drawdown. This checks a layered aquifer system: a second material,
a fixed-pressure top, the well, storage and radial diffusion. Mechanical
pressure feedback is off, as in the Theis check.

![Setup](../../benchmarks/poromechanics/results/review/fig/hantush/setup.svg)

DES meshes a quarter box, $0\le x,y\le75$ m. The aquifer is $B=2.5$ m
thick; above it the aquitard is one element thick ($b'=h$), and its top
holds $p=0$. The base and the far sides are sealed. A well on the corner
edge, screened through the aquifer only, pumps $Q/4$ with
$Q=10^{-4}\ \mathrm{m^3/s}$ for 200 s and then shuts in until 320 s, in
160 steps of 2 s. The aquifer has diffusivity 1 m²/s; the aquitard stores
a thousandth of the aquifer's water per volume.

## Analytical solution

With aquifer transmissivity and storativity in pressure units,
$\mathcal{T}_p=\kappa B/\mu$ and $S_p=BS_B$, an aquitard without storage
carries a linear vertical pressure profile and leaks
$(\kappa_z'/\mu b')\,s_p$ into the aquifer, so the drawdown
$s_p=p_{\rm initial}-p$ obeys

$$
S_p\,\partial_t s_p=\mathcal{T}_p\left(\partial_r^2 s_p+\frac{1}{r}\partial_r s_p\right)-\frac{\mathcal{T}_p}{\lambda^2}s_p,
\qquad \lambda^2=\frac{\kappa B\,b'}{\kappa_z'} .
$$

For a constant rate from $t=0$ (Hantush and Jacob, 1955),

$$
s_p=\frac{Q}{4\pi\mathcal{T}_p}\,W\!\left(\frac{r^2S_p}{4\mathcal{T}_pt},\frac{r}{\lambda}\right),\qquad
W(u,\beta)=\int_u^\infty\frac{1}{y}\exp\!\left(-y-\frac{\beta^2}{4y}\right)dy .
$$

$W(u,0)=E_1(u)$ is Theis and $W(0,\beta)=2K_0(\beta)$ the steady cone.
Shut-in is superposed. `benchmarks/poromechanics/hantush_reference.py`
evaluates $W$ by quadrature and checks both limits, the identity
$W(u,\beta)+W(\beta^2/4u,\beta)=2K_0(\beta)$ and the equation itself.

## Matching the assumptions

The theory neglects aquitard storage, horizontal aquitard flow and vertical
resistance in the aquifer. The aquitard leakance is held fixed on every
mesh, so $\lambda=15$ m throughout. Vertical permeability is ten times
horizontal in both layers. The remaining departures are aquitard storage
$S'_B b'/S_p\le0.12\%$, where $S'_B$ is aquitard storage per unit
volume and $S_p=BS_B$ is aquifer storage per unit area, horizontal
aquitard transmissivity
$(b'/\lambda)^2/10\le0.28\%$ (both falling with $h$) and the aquifer's
vertical pressure drop $B^2/(20\lambda^2)=0.14\%$.

## Results

Drawdown at the aquifer nodes with $10\le r\le30$ m, at eight times through
pumping and recovery, relative to $C_Q=Q/(4\pi\mathcal{T}_p)$:

| Mesh $h$ | Nodes | Max error | RMS error |
| --- | --- | --- | --- |
| 2.5 m | 2,883 | 6.53% | 1.38% |
| 1.25 m | 14,884 | 1.49% | 0.30% |
| 0.625 m | 87,846 | 0.41% | 0.076% |

The error falls with order 2.00. It is largest near the well during
pumping and falls below 0.1% of $C_Q$ after shut-in on the two finer meshes. Leakage changes the drawdown by
up to 53% of $C_Q$ from Theis, so the comparison tests the leakage rather
than the Theis part alone. The skeleton strain stays below $10^{-14}$.

## Reproduce

```
python3 benchmarks/poromechanics/hantush_reference.py
python3 tests/functional/hydromechanics/check_hantush.py ./dynearthsol3d OUTPUT --spacings 2.5 1.25 0.625
```

The recorded study, with input and executable hashes, is
`results/runs/bishop-20261108/hantush/`. The 0.625 m mesh takes about
45 minutes on one thread.
