# Fluid injection into a poroelastic borehole

The poroelastic step of the tunnel ladder (Kirsch, Salençon, then this):
fluid at pressure $p_w$ enters a borehole wall, loads it and diffuses into
the rock, which changes the stresses around the hole as it goes. This is
the axisymmetric part of Detournay and Cheng (1988), solved here on the
finite annulus that DES models.

## Problem

A quarter annulus, $a=1$ m and $R=20a$, in plane strain, with symmetry
rollers on both cuts and meshes refined out to $5a$ (the Kirsch geometry).
The rock starts unstressed at $p=0$; gravity and heat are off; pressure
feedback is on. At $t=0^+$ the wall takes pore pressure $p_w=1$ MPa (a 1 ms
ramp of the scheduled boundary pressure) and normal traction $-p_w$. The
outer arc is sealed and traction-free. $K=200$ MPa, $G=100$ MPa,
$\alpha=1$, $1/M=\phi/K_f=10^{-9}$ Pa$^{-1}$, $k/\mu=10^{-9}$
m$^2$/(Pa s), so $c=0.25$ m$^2$/s; the front reaches about $3a$ at 36 s.

The recorded study uses a sealed outer arc. The reference also handles a
drained outer arc. Prescribed-pressure validation now checks actual shared
tagged nodes on imported meshes, so the disjoint wall and outer arc may
take different pressures despite their $x_1$ and $z_1$ labels. The results
below still refer to the sealed study; a drained DES comparison is separate.

## Reference

With radial symmetry $u=mJ/r+C_1r+C_2/r$, $J=\int_a^rp\,r'dr'$,
$m=\alpha/(\lambda+2G)$, so the volumetric strain is $mp+2C_1$ and

$$
\sigma_{rr}=-2Gm\frac{J}{r^2}+2(\lambda+G)C_1-2G\frac{C_2}{r^2},\qquad
\sigma_{\theta\theta}=2Gm\left(\frac{J}{r^2}-p\right)+2(\lambda+G)C_1+2G\frac{C_2}{r^2},
$$

with $C_1$, $C_2$ from $\sigma_{rr}(a)=-p_w$, $\sigma_{rr}(R)=0$. Fluid
mass gives $S\,p_t+2\alpha\dot C_1=(k/\mu)\nabla^2p$,
$S=1/M+\alpha m$. In an infinite medium $C_1=0$ and this is plain
diffusion; on the annulus $C_1$ depends on $J(R)$ and couples the pressure
to its own integral, moving it by about 1% of $p_w$. In Laplace space the
pressure is a Bessel solution plus a uniform part $-2\alpha\bar C_1/S$,
with $\bar C_1$ from one scalar equation, inverted by Stehfest.
`borehole_reference.py` checks the wall and outer conditions, radial
equilibrium, the steady limits, and the inversion against an independent
Crank–Nicolson solution of the coupled equation (agreement $5\times10^{-6}$
of $p_w$, drained and sealed).

## Results

Largest nodal pressure error and element total-stress errors for
$1.05a\le r\le4a$, relative to $p_w$ (144 steps of 0.25 s unless stated):

| Mesh | Pressure, 9 / 18 / 27 / 36 s | Pressure rms, 36 s | Hoop stress max / rms | Radial stress max |
| --- | --- | --- | --- | --- |
| 0.2 a | 0.111 / 0.117 / 0.125 / 0.130% | 0.090% | 4.57 / 0.416% | 3.14% |
| 0.1 a | 0.067 / 0.037 / 0.035 / 0.034% | 0.021% | 1.85 / 0.165% | 1.59% |
| 0.05 a | 0.058 / 0.017 / 0.009 / 0.007% | 0.002% | 1.06 / 0.084% | 1.03% |
| 0.1 a, 576 steps | 0.039 / 0.036 / 0.040 / 0.039% | 0.027% | 1.88 / 0.165% | 1.59% |

Late-time pressure converges at about second order with the mesh; the
first frame is limited by the step (0.067% falls to 0.039% with four times
the steps on the 0.1 a mesh). Stress converges at first order, as expected
of constant-stress triangles; the largest element errors sit at the wall,
where the hoop stress is steepest.

## Reproduce

```
python3 benchmarks/poromechanics/borehole_reference.py
python3 tests/functional/hydromechanics/check_borehole.py EXE OUT --spacings 0.2 0.1 0.05
python3 tests/functional/hydromechanics/check_borehole.py EXE OUT --spacings 0.1 --steps 576
```

Recorded studies: `results/runs/borehole-20261108/{mesh-2d,steps576-2d}/`.
The independent `replay.py` in their parent directory recomputes all
pressure and stress metrics for all four cases and 16 saved comparisons.
`independent-replay.json` records the reader/reference hashes, the original
solver hash, the study artifact hashes and the independent reference-check
log. This is a replay of those recorded
studies, not a new solver run or a plastic-tunnel qualification.
The field comparison is 2D; it does not establish 3D borehole accuracy.
