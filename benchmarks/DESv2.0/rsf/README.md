# RSF models for the revised DES v2.0 paper

The revised paper uses four main model groups: simple shear, extensional
localization, a decollement parameter study, and a strike-slip comparison.
This directory contains 108 configurations and three meshes. The full inventory
and appendix mapping are in [PAPER_CASES.md](PAPER_CASES.md).

## Setup

Use a 2D executable built from the RSF source in
[PR #91](https://github.com/GeoFLAC/DynEarthSol/pull/91), for example public commit
[`1b4e94f`](https://github.com/GeoFLAC/DynEarthSol/tree/1b4e94f1c333fdade9dc512cf9cb20e10e8e9aad).
Follow that source revision's build instructions in a separate source checkout.
The solver and simple-shear scripts retained elsewhere on this benchmark branch
predate the revised RSF controls. The linked source is a PR revision, not a
published v2.0.1 release.

From the root of this benchmark checkout, prepare a working copy that preserves
the relative mesh paths:

```bash
export SOURCE_DIR=/absolute/path/to/compatible-source
export DYNEXE="$SOURCE_DIR/dynearthsol2d"
export OMP_NUM_THREADS=8
RUN_DIR=$(mktemp -d)
cp -R benchmarks/DESv2.0/rsf "$RUN_DIR/rsf"
export RSF_RUN="$RUN_DIR/rsf"
```

Choose the thread count for your machine. Use a fresh copy for another run;
many configurations in the same directory share monitor output names.

## 1. Simple shear — Figure 10

The nine cases verify EP, steady-state RSF, and aging-law RSF against reference
solutions. Use the checker from the compatible source checkout; it also runs
the separate zero-velocity healing test:

```bash
cd "$SOURCE_DIR/benchmarks/simple_shear_rsf"
python3 check_simple_shear_benchmark.py --exe "$DYNEXE" --all
```

The paper's individual inputs are in [shear_box/](shear_box/). For figure
production, follow the source checkout's
[benchmark guide](https://github.com/GeoFLAC/DynEarthSol/blob/1b4e94f1c333fdade9dc512cf9cb20e10e8e9aad/benchmarks/simple_shear_rsf/README.md).
The checker uses the Python standard library; plotting requires NumPy and
Matplotlib.

## 2. Extensional localization — Figure 11

Three models compare EP, velocity-strengthening EP-RSF, and velocity-weakening
EP-RSF. They generate their meshes internally and integrate 12000 model years.
Each configuration has a distinct output prefix.

```bash
cd "$RSF_RUN/localization"
mkdir -p output
"$DYNEXE" ep_phi30_psi0.cfg
"$DYNEXE" ep_rsf_strengthening_phi30_psi0.cfg
"$DYNEXE" ep_rsf_weakening_phi30_psi0.cfg
```

See [localization/README.md](localization/README.md) for the model setup and
field definitions.

## 3. Decollement parameter study — Section 5.4.2

The parameter study varies `a/b` and `D_c` in a 100 m fault zone. It contains
60 grid cases and 14 constant-nucleation-length contour cases. The example below
runs one representative grid case, `a/b = 0.5` and `D_c = 0.005 m`:

```bash
cd "$RSF_RUN/decollement/grid"
mkdir -p output
"$DYNEXE" ab0p5_dc0p005.cfg
```

Other study cases use the same procedure with their own cfg and a fresh working
copy. This one run does not reproduce the entire parameter map. Classifications
use the record after the first 1000 model years. The required mesh is included.

## 4. Strike-slip comparison — Section 5.4.3

These models compare a 150 km domain containing a 250 m band with Herrendoerfer
et al. (2018). For example, run the three-element-band case with per-element
`D_c` approximately 0.00333 m:

```bash
cd "$RSF_RUN/strike_slip"
mkdir -p output
"$DYNEXE" herrendorfer2018_A_c3dc3.cfg
```

The configurations include the paired monitoring points used to calculate
fault-relative slip. The required meshes are included under [mesh/](mesh/).

## Reproduction limits

The two Appendix D inputs in `decollement/wavefield/` require checkpoints that
are not included. One final-condition strike-slip sensitivity input is also
unavailable; see [PAPER_CASES.md](PAPER_CASES.md). The full application campaigns
have not been repeated on the linked PR revision. The localization inputs have
passed short startup checks, which do not reproduce their final 12 kyr fields.
