# Rapid drawdown, selected 0.5 m mesh

Use a 2D DynEarthSol executable from source revision 3774df00ba1708a4f70a836d2073f845653aac64.
The supplied reservoir frame is the hydrostatic starting state at 25000 s.
Preserve its save, checkpoint and info files together with slope.poly.
Run from this directory:

```bash
export OMP_NUM_THREADS=1
dynearthsol2d input.cfg
```

This continues drainage for 20000 s. The reference compares Bishop's method
on the DES pore pressure with strength reduction; it is not an exact analytical
solution. Three of six strength brackets contain the reference. The other
three are slightly conservative; all-times containment is not achieved.
The 72 trials and six Bishop searches were replayed remotely from saved fields.
record.json contains the selected results and solver hash.

To regenerate the whole study, run from the matching DynEarthSol source:

```bash
python3 benchmarks/poromechanics/drawdown_slope.py /absolute/path/to/dynearthsol2d OUT --spacings 0.5
```

This bundle contains only the selected restart input. The study and reference
code already live in the matching source checkout.
