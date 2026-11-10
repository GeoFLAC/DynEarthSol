# Rapid drawdown of a reservoir against an embankment

An application: the [embankment](embankment.md) holds back a reservoir at
its crest level, and the reservoir is drawn down quickly. The water no
longer supports the face, while the pore pressure inside falls only as the
slope drains. This is the critical case for the upstream slope of a dam
(Morgenstern 1963). There is no analytical solution, so the reference is
Bishop's simplified method applied at each time to the pore pressure DES
has computed.

## Model

- **Geometry and strength:** those of the embankment page: a 2:1 face with
  $H=10$ m, a foundation $H$ deep, $c=10$ kPa, $\phi=20^\circ$, $\psi=0$,
  $E=100$ MPa and $\nu=0.3$.
- **Saturated state:** porosity 0.3 and saturated unit weight 20 kN/m³
  (grain density 2484 kg/m³).
- **Retention and permeability:** Gardner retention with
  $\alpha=2.4\times10^{-4}$ Pa$^{-1}$, $S_r=0.05$, Bishop $\chi=S_e$, and a
  permeability of $10^{-13}$ m².
- **Selected mesh:** 0.5 m, 5,822 nodes and 11,337 triangles. The completed
  1 m study remains a separate comparison. Computation used the remote
  cluster, with 12 independent strength trials running concurrently.

## Method

1. **Seed:** no gravity and zero stress (checkpoint).
2. **Reservoir:** 25000 s, explicit throughout.
   - DES's water loading (`has_water_loading`, `surf_base_level=0`) loads
     the face and the toe platform.
   - The crest is held at $p=0$, the water level, and the base at its
     hydrostatic pressure. The slope fills to hydrostatic; the misfit is
     0.248% of the base pressure.
3. **Drawdown:** 20000 s.
   - The water load is removed and the face and toe platform drain
     ($p=0$). The crest and base are sealed.
   - The mechanics stays elastic, as limit equilibrium assumes.
4. **Trials:** at 250, 1000, 3000, 6000, 12000 and 20000 s, the
   strength-reduction trials of the embankment page restart from the
   drainage checkpoint, with the same settings and the same at-rest
   criterion.
   - The pore pressure is held: the permeability is frozen ($10^{-25}$ m²)
     and the pore fluid is made compressible (1 kPa bulk modulus).
   - A stiff fluid would make the trial undrained and add excess pore
     pressures of order 10 kPa that limit equilibrium does not have.

Every stage after the seed is explicit, with inertial scaling $10^5$ and
DES's default damping. DES refuses to change the hydraulic formulation
between PT and explicit at a restart, so the reservoir and drawdown stages
cannot be PT while the trials are explicit.

At each time, `bishop()` interpolates DES's $\chi p$ and saturation onto
the slices. The slice weight uses the unit weight
$g(\rho_s(1-n)+\rho_f n S)$, and the critical circle is searched again.

## Results

| Time after drawdown | Bishop on DES's pore pressure | DES bracket |
| --- | --- | --- |
| 250 s | 0.857 | 0.80–0.85 |
| 1000 s | 0.877 | 0.85–0.90 |
| 3000 s | 0.890 | 0.85–0.90 |
| 6000 s | 0.902 | 0.85–0.90 |
| 12000 s | 0.913 | 0.85–0.90 |
| 20000 s | 0.918 | 0.90–0.95 |

**Three of six brackets contain Bishop.** At 250, 6000 and 12000 s,
Bishop exceeds the bracket's upper edge by 0.00685, 0.00202 and 0.01341,
respectively (0.2–1.5% of Bishop). DES is slightly conservative at all
three times. The largest midpoint difference is 4.20%; this is a different
measure from the upper-edge miss. The trial spacing is 0.05. The original
criterion of containing Bishop at every time is therefore not met on this
mesh; it is retained rather than relaxed. All 72 trials and six independent
Bishop searches were replayed from their raw fields and reproduced the
recorded metrics exactly.

**The slope fails on drawdown.** It stands dry at 1.369 but drops below 1
when drawn down, and drainage restores only a little of its strength in
20000 s.

**The mechanism matches Bishop's circle.** It passes through the toe and
leaves the crest about 2 m behind its edge. The circle is centred at
(25.0, 8.0) m with a radius of 18.7 m at 250 s, and it widens slightly
as the slope drains (19.8 m at 20000 s).

The older 1 m drawdown study contains Bishop at every time, with a 2.5%
maximum midpoint difference. Refinement preserves the failure-and-recovery
trend but does not improve every bracket. This application is presented
with its small conservative difference; no further mesh or factor tuning
is used to match Bishop.

## Reproduce

```
python3 benchmarks/poromechanics/drawdown_slope.py EXE OUT \
  --spacings .5 --workers 12 --flow-threads 12
```

The review's drawdown page ships the reservoir and drawdown inputs with
the restart files, the steps to set up a trial, and the study's numbers
(`record.json`).

## References

- Bishop (1955), Géotechnique 5(1), 7–17.
- Morgenstern (1963), Géotechnique 13(2), 121–131.
