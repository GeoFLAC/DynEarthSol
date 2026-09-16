# Rate-and-State Friction Benchmarks

This directory provides rate-and-state friction (RSF) benchmarks and example
problems for DynEarthSol: simple shear, extensional localization, a decollement
parameter study, and a strike-slip comparison. It contains 108 configurations
and three meshes. See [CASES.md](CASES.md) to choose a configuration and check
its input requirements.

## Setup

Use a 2D executable built from a compatible RSF source revision, such as
[`1b4e94f`](https://github.com/GeoFLAC/DynEarthSol/tree/1b4e94f1c333fdade9dc512cf9cb20e10e8e9aad).
Follow that source revision's build instructions in a separate source checkout.
The solver and simple-shear scripts retained elsewhere on this benchmark branch
predate the revised RSF controls. The link pins the source revision used for
these instructions.

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

## 1. Simple-shear verification

The nine cases verify elastoplasticity (EP), steady-state RSF, and aging-law RSF
against reference solutions. Use the checker from the compatible source
checkout; it also runs the separate zero-velocity healing test:

```bash
cd "$SOURCE_DIR/benchmarks/simple_shear_rsf"
python3 check_simple_shear_benchmark.py --exe "$DYNEXE" --all
```

Individual simple-shear inputs are in [shear_box/](shear_box/). For figure
production, follow the source checkout's
[benchmark guide](https://github.com/GeoFLAC/DynEarthSol/blob/1b4e94f1c333fdade9dc512cf9cb20e10e8e9aad/benchmarks/simple_shear_rsf/README.md).
The checker uses the Python standard library; plotting requires NumPy and
Matplotlib.

## 2. Extensional localization

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

## 3. Decollement parameter study

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

## 4. Strike-slip comparison

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

## Restart inputs and validation

The two inputs in `decollement/wavefield/` require checkpoints that are not
included; see [CASES.md](CASES.md). They cannot start from the supplied cfgs alone.

The simple-shear and healing checks passed on the linked source revision.
The localization inputs passed short startup checks; their full 12000-year runs
and the other long application runs have not been repeated on that revision.
