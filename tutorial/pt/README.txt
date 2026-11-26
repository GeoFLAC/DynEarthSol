DES3D pseudo-transient mechanics tutorials
========================================

Manual pages:
https://geoflac.github.io/DynEarthSol/docs/tutorial/pt-elastic
https://geoflac.github.io/DynEarthSol/docs/tutorial/pt-maxwell

Reference source: GeoFLAC/DynEarthSol
Version for reproducing the figures: 707e481fea957971d5cb6ca67d8d24bd689a3467
Reference figures: GCC 13.3, Linux, OpenMP with one thread, binary output.
Both dimensions and both PT options were checked. Published figures use 2D.

Create pt-tutorial/ inside a checkout of that revision and extract this
bundle there. From the checkout root, build:

    make -j2 ndims=2 hdf5=0
    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install numpy matplotlib
    export OMP_NUM_THREADS=1
    export OMP_DYNAMIC=FALSE
    cd pt-tutorial

Run and analyze:

    ../dynearthsol2d elastic-0.cfg > elastic-0.log 2>&1
    ../dynearthsol2d elastic-1.cfg > elastic-1.log 2>&1
    python analyze_pt.py elastic elastic-0 elastic-1 --source ..
    ../dynearthsol2d maxwell-0.cfg > maxwell-0.log 2>&1
    ../dynearthsol2d maxwell-1.cfg > maxwell-1.log 2>&1
    python analyze_pt.py maxwell maxwell-0 maxwell-1 --source ..

The same inputs work with a 3D build. Use a fresh output directory and
dynearthsol3d. The analysis script detects the dimension. hdf5=1 output
also needs h5py in the Python environment; this bundle was checked with
hdf5=0 and needs no HDF5 installation.

Inputs derive from examples/pt-elastic-load.cfg in the checked revision.
The elastic configs differ in modelname and PT_option only. Maxwell adds
100 steps, maxwell rheology, plane strain in 2D and fixed viscosity 1e6 Pa s.

All examples explicitly use PT_use_running_scale=yes. This is an opt-in
convergence criterion; the solver default is no. It has limited validation
for history-dependent materials. See the manual for the limitations.

The fixed mesh does not move. Equivalent top displacement is reconstructed
by integrating velocity. Negative stress, strain and displacement denote
compression. The analysis checks the supplied 0.01 s timestep and 3/100-step
histories, all components of stress/strain/coordinates/velocity for
finiteness, axial stress/strain/vertical-velocity reference values, and PT reports.
Reference CSVs in reference/ contain the checked 2D results. Figures are
generated from solver output by analyze_pt.py; no output data are fabricated.

To regenerate the manual figures, copy elastic-response.svg,
elastic-convergence.svg and maxwell-response.svg from a checked 2D run into
static/tutorial/pt/ in the documentation repository.

Expected responses:
elastic: strain_zz about -0.004285714, stress_zz about -100000 Pa;
         subsequent held steps retain strain without further compression.
Maxwell: strain_zz about -0.00991968 at 1 s, stress_zz about -100000 Pa;
         creep continues after the first loading step.

These are homogeneous dry fixed-grid checks. They do not qualify arbitrary
hydraulic coupling, moving meshes or strain-softening localization.
