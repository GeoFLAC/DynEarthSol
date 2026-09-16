<!-- Citation & License -->
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20293557-blue.svg)](https://doi.org/10.5281/zenodo.20293557)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<br> <!-- Core Platforms & Compilers -->
[![Basic build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/basic-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/basic-build.yml)
[![macOS build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/macos-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/macos-build.yml)
[![nvc build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/nvc-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/nvc-build.yml)
[![OMP g++ matrix](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/omp-gcc-matrix.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/omp-gcc-matrix.yml)
[![Functional tests](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/functional-tests.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/functional-tests.yml)
<br> <!-- External Modules & Integrations -->
[![Exodus build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/exodus-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/exodus-build.yml)
[![MMG build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/mmg-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/mmg-build.yml)
[![GoSPL build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/gospl-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/gospl-build.yml)

# Overview

DynEarthSol, DES in short, is a finite element code that solves the momentum
balance and the heat transfer in Lagrangian form using unstructured meshes, in
two or three dimensions. It can be used to study the long-term deformation of
Earth's lithosphere and problems alike.

# Building DES

## Getting the Source Code
This repository uses Git submodules for three libraries: `nanoflann` (always),
`mmg` (with `usemmg=1`) and `knn-bvh` (with `openacc=1`). `make` initializes the
ones a build needs, fetching them over the internet before compiling.

For environments without internet access (such as certain HPC compute nodes), pre-downloading all dependencies is highly recommended. To ensure all necessary source files are downloaded upfront, clone the repository with its submodules using the `--recurse-submodules` flag:

```bash
git clone --recurse-submodules https://github.com/GeoFLAC/DynEarthSol.git
```

If you have already cloned the repository without the submodules, you can initialize and update them by running the following command inside the DynEarthSol directory:

```bash
git submodule update --init --recursive
```

## Requirements

Install the dependencies with your package manager; the build finds them on its
own, with no paths to edit and none to pass on the command line.

```bash
# macOS (Homebrew) -- libomp because Apple clang ships no OpenMP runtime
brew install boost libomp

# Debian / Ubuntu
sudo apt install g++ make libboost-program-options-dev

# Fedora / RHEL
sudo dnf install gcc-c++ make boost-devel
```

The requirements in detail:
* You will need a C++11 compiler. CI builds with g++ 8 through 15, Apple clang
  and nvc++.
* You will need `Boost::Program_options` (1.42 or newer). Packaged versions are
  found automatically. On macOS the search order is the Homebrew `boost` keg,
  then MacPorts, then the Homebrew prefix itself, then an activated conda
  environment; elsewhere it is the conda environment, then the system Boost the
  compiler finds on its own. A prefix is accepted only if it holds both the
  headers and the library. Package managers come before conda deliberately: a
  shell that auto-activates `base` has conda on for everything, so it is the
  weaker signal of intent, but a conda Boost is still used when it is the only
  one installed.
  * To build Boost yourself: download the source from www.boost.org, run
    `./bootstrap.sh` then `./b2 --with-program_options -q` in the untarred
    directory, and point the build at it with `make BOOST_ROOT_DIR=/path/to/boost`.
    A build directory is recognised by its `stage/` subdirectory, so name the
    untarred directory, not `stage/lib`.
* The Python tools (`2vtk.py`, `Dynearthsol.py`) need Python 3 with NumPy and
  SciPy, plus h5py to read vtkhdf output.
* **macOS users**: Apple clang has no built-in OpenMP, so the runtime has to
  come from somewhere; `brew install libomp` is the whole answer for most
  people. See [OpenMP on macOS](#llvm) for the search order, the overrides, and
  how to build LLVM OpenMP from source when a package manager is not an option.

## Building

From the repository root:

```bash
make ndims=2    # 2D executable, dynearthsol2d
make            # 3D executable, dynearthsol3d
```

That is the whole procedure. `make config` prints the compiler and every
dependency path the build resolved to; it is the quickest way to see that
something was found somewhere you did not expect, and the most useful thing to
paste into a bug report. The subsections below cover the cases where the
defaults are not what you want.

### Build options

* Options go on the command line (`make ndims=2 hdf5=1`) or into the
  corresponding variable at the top of the `Makefile`; the command line wins.
  Changing a flag rebuilds what it reaches, so no `make clean` is needed.
* `make check-deps` verifies the external libraries, including their
  architecture, without building.
* `BOOST_ROOT_DIR`, `HDF5_INCLUDE_DIR`, `HDF5_LIB_DIR` and `NVHPC_DIR` sit
  together in an **Optional paths** block near the top of the `Makefile`. All
  may stay blank; set one only to force a specific install, which skips
  detection for that dependency. Every other dependency's path is declared
  beside the feature that uses it (`OPENMP_*`, `EXO_*`, `MMG_*`, the GoSPL group).
  An *exported* value is honoured too, which is how a module file or a conda
  activation hook names a prefix once; the catch is that a stale export in a
  shell profile silently defeats a good package-manager install. `make config`
  shows which prefix was used, and an empty `BOOST_ROOT_DIR=` on the command
  line forces detection.
* `opt=2` is the default optimized build; `opt=3` adds `-march=native -O3`;
  `opt=0` is for debugging; `opt=-1` adds `-fsanitize=address`, which reports
  where a memory error occurs and where the memory was allocated without gdb or
  valgrind, but cannot be combined with valgrind.
* `openmp=0` builds a single-threaded executable, which is what valgrind needs.

### Submodules
* [nanoflann](https://github.com/jlblancoc/nanoflann): header-only KD-tree for
  nearest-neighbour searches on the CPU. Always required.
* [mmg](https://github.com/MmgTools/mmg): mesh optimization during remeshing.
  Required with `usemmg=1`; `make` configures and builds it under `mmg/build`.
* [knn-bvh](https://github.com/GeoFLAC/knn-bvh): K-nearest-neighbour searches
  on Bounding Volume Hierarchies, the GPU counterpart of nanoflann. Required with
  `openacc=1`.

### Optional packages
* [Exodus](https://github.com/sandialabs/seacas/), `useexo=1`, for importing a
  mesh in the ExodusII format. 3D only.
  * Suggested building procedure, run in the root directory of DES; it
    downloads and builds NetCDF and HDF5, then Exodus, with headers and libraries
    landing in `./seacas/include` and `./seacas/lib`:
    ```BASH
    git clone https://github.com/sandialabs/seacas.git
    cd seacas && export ACCESS=`pwd`
    COMPILER=gnu MATIO=NO GNU_PARALLEL=NO CGNS=NO FMT=NO ./install-tpl.sh
    mkdir build; cd build
    ../cmake-exodus
    make; make install
    ```
  * `EXO_INCLUDE` and `EXO_LIB_DIR` override those defaults.
* [MMG](https://www.mmgtools.org/), `usemmg=1`, for mesh optimization during
  remeshing, in 2D and 3D. It is the `mmg` submodule: `make usemmg=1` initializes
  it and builds `libmmg2d` or `libmmg3d` under `mmg/build`, so nothing is
  installed by hand. The configure step ignores an exported `CFLAGS`; a tree
  configured earlier keeps its cached flags until `mmg/build` is removed.
  `MMG_INCLUDE` and `MMG_LIB_DIR` point at a different build.
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/), `hdf5=1`, for writing
  results in the HDF5-based [vtkhdf](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format)
  format, which is compressed (up to 50 % smaller) and opens directly in ParaView.
  * Install it with `brew install hdf5` (macOS), `apt install libhdf5-dev`
    (Debian/Ubuntu) or `dnf install hdf5-devel` (Fedora/RHEL). The build locates
    it via `pkg-config` and then the usual install layouts, picking the prefix
    that matches the host architecture: on a Mac that has used both Homebrew
    prefixes, `h5cc` and `pkg-config` on `PATH` are often the x86_64 ones.
  * `HDF5_INCLUDE_DIR` and `HDF5_LIB_DIR` force one specific HDF5, an HPC module,
    a conda env or a hand-built copy:
    `make hdf5=1 HDF5_INCLUDE_DIR=/prefix/include HDF5_LIB_DIR=/prefix/lib`.
    Either one alone is enough; the other is still detected.
* [GoSPL](https://gospl.readthedocs.io/) (Global Scalable Paleo Landscape
  Evolution), `use_gospl=1`, for two-way coupling with a surface process model.
  GoSPL handles erosion, deposition and hillslope diffusion; DES handles
  tectonics. At each coupling event DES surface velocities are passed to GoSPL,
  which returns the erosion/deposition increment applied to DES surface nodes.
  * Requires a GoSPL conda environment, following
    [the GoSPL installation procedure](https://gospl.readthedocs.io/en/latest/getting_started/installConda.html),
    and the [gospl_extensions](https://github.com/GeoFLAC/gospl_extensions) C++
    interface library:
    ```bash
    git clone https://github.com/GeoFLAC/gospl_extensions.git ~/opt/gospl_extensions
    cd ~/opt/gospl_extensions/cpp_interface
    conda activate gospl
    make install-local
    ```
  * `GOSPL_EXT_DIR` and `CONDA_ENV_PATH` override the defaults
    (`~/opt/gospl_extensions` and `~/miniconda3/envs/gospl`). `make use_gospl=1`
    also generates the `dynearthsol-gospl` wrapper script.
  * See `gospl_driver/README.md` for build, runtime and coupling details, and
    `gospl_driver/examples/` for example configs. The Docker image below ships
    the whole stack ready to run.

<div id="llvm"></div>

### OpenMP on macOS

Apple clang implements the OpenMP *pragmas* but ships no OpenMP *runtime*, so
one has to be installed. Almost always this is enough:

```bash
brew install libomp
```

Nothing else to configure: `make` locates the header and the library itself.
`make config` reports which one it picked.

* **Search order**, highest priority first:
  1. `external/openmp-install` in the source tree — an LLVM OpenMP built here by
     hand (see below) keeps taking priority over anything installed system-wide
  2. Homebrew `libomp`
  3. Homebrew `llvm`
  4. MacPorts `libomp`
  5. the activated conda environment, `$CONDA_PREFIX`

  `omp.h` and `libomp.dylib` are searched separately over that same list,
  because Homebrew's `llvm` keeps its `omp.h` under `lib/clang/<version>/include`
  rather than beside the library. Whichever provider is installed therefore
  supplies both — but a provider holding only one half is skipped for that half
  alone, so `make config` is worth a glance if you have several installed.

  Homebrew is looked for at the prefix matching the machine's architecture --
  `/opt/homebrew` on Apple Silicon, `/usr/local` on Intel -- and never via
  `PATH`. On a Mac that has run both, `PATH` regularly offers the other one's
  tools, and building against those produces libraries of the wrong
  architecture. Pass `BREW_PREFIX=/your/prefix` if yours is somewhere else.

* **Overrides**, when the runtime is somewhere the search does not look, or when
  you want a specific one. Set these on the command line, or edit them in the
  clang++ branch of the `Makefile`, where the OpenMP search lives:

  ```bash
  make OPENMP_ROOT_DIR=/prefix                        # /prefix/include + /prefix/lib
  make OPENMP_INCLUDE_DIR=... OPENMP_LIB_DIR=...      # when they are not siblings
  ```

  The second form is what Homebrew's `llvm` needs if named explicitly, for the
  reason given above. The exported-variable behaviour is the same as for the
  Optional paths block.

* **Thread-wait performance on Apple Silicon.** LLVM `libomp` hardcodes
  hybrid-CPU detection for all Apple Silicon, setting thread *blocktime* to
  **0 µs**: threads yield immediately after each parallel region instead of
  spin-waiting, which cuts CPU utilisation to ~250 % on a 6-core Mac against
  ~600 % on Linux. DES therefore sets `OMP_WAIT_POLICY=active` at startup on
  macOS unless `OMP_WAIT_POLICY` or `KMP_BLOCKTIME` is already set, and prints a
  notice saying so. To override:
  ```bash
  export OMP_WAIT_POLICY=passive   # yield immediately, lower power
  export KMP_BLOCKTIME=20          # fine-grained control (ms, libomp only)
  ```

* **Building LLVM OpenMP from source**, if a package manager is not an option.
  The result lands in `external/openmp-install`, which the search prefers, so no
  build variable has to be set afterwards. (LLVM OpenMP 19.1.7 or newer will
  suffice; replace `arm64` with `x86_64` on an Intel Mac.)

  ```BASH
  mkdir -p external && cd external
  curl -L https://github.com/llvm/llvm-project/releases/download/llvmorg-19.1.7/cmake-19.1.7.src.tar.xz -o cmake-19.1.7.src.tar.xz
  tar xf cmake-19.1.7.src.tar.xz
  curl -L https://github.com/llvm/llvm-project/releases/download/llvmorg-19.1.7/openmp-19.1.7.src.tar.xz -o openmp-19.1.7.src.tar.xz
  tar xf openmp-19.1.7.src.tar.xz

  mkdir -p openmp-19.1.7.src/build && cd openmp-19.1.7.src/build
  cmake -DCMAKE_INSTALL_PREFIX=$(pwd)/../../openmp-install \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_MODULE_PATH=$(pwd)/../../cmake-19.1.7.src/Modules \
        -DCMAKE_OSX_ARCHITECTURES=arm64 \
        -DLIBOMP_INSTALL_ALIASES=OFF \
        ..
  make -j4 && make install && cd ../../
  ```

### Docker

* Build the image. `build.sh` sets the dimension (`NDIMS=2`) and the compiler
  (`gcc-11`; also `gcc-8`, `clang-14`) at its top; edit them there or pass them
  in the environment.
  ```bash
  ./build.sh
  ```
* Run it
  ```bash
  docker run --rm -it dynearthsol/gcc-11
  ```
* `GOSPL=1 ./build.sh` builds `dynearthsol/gcc-11-gospl` instead: a 3D
  executable with GoSPL coupling, the `gospl` conda environment and
  `gospl_extensions` included, so no conda setup is needed on the host. Mount a
  directory holding the cfg and the GoSPL YAML and run `dynearthsol-gospl`
  inside it; see `gospl_driver/README.md`.
* `docker/Dockerfile.cuda` is the NVHPC image the GPU CI builds in.

### Common make invocations

Here are a few practical examples for common build configurations (run these from the project root):

```bash
# default optimized 3D build
make

# show the compiler, flags and dependency paths the build resolved to
make config

# check the external libraries are present, without building
make check-deps

# debugging build (no optimizations, no OpenMP)
make opt=0 openmp=0

# build 2D version
make ndims=2

# enable MMG mesh optimization (builds the mmg submodule)
make usemmg=1

# enable Exodus mesh import, 3D only (requires seacas/exodus libs)
make useexo=1

# enable HDF5-based vtkhdf output support (requires HDF5)
make hdf5=1

# enable GoSPL surface process coupling (requires gospl_extensions and gospl conda env)
make use_gospl=1

# NVHPC/profiler build (uses nvc++ when set)
make nprof=1

# embed the uncommitted code changes in the executable (off by default)
make snapshot_diff=1

# OpenACC build (NVHPC compiler)
make openacc=1

# OpenACC targeting a specific GPU (else from nvidia-smi, default 80; see make config)
make openacc=1 GPU_CC=90
```

# Running DES
* Execute `dynearthsol2d input.cfg` (or `dynearthsol3d`). The input file is
  required; `-h` or `--help` lists every parameter with its description.
* The input format is documented in `examples/defaults.cfg`, and working cases
  are under `examples/`.
* Process exit codes are two digits whose first digit says who fixes it: 1x the
  input, 2x the environment, 3x-6x the code. The table is in
  [CONTRIBUTING.md](CONTRIBUTING.md).
* The [input file generator](https://geoflac.github.io/des-inputgen) builds a
  cfg from a web form.
* Benchmark cases with analytical solutions are under `benchmarks/`; the
  regression cases the developers compare against are under `benchmarks-cores/`.
* **Build provenance**: every executable embeds a `build.snapshot` block -- revision,
  make options, dependency versions and providers, toolchain, compile and link flags:
  ```bash
  strings <exe> | grep '^build\.snapshot\.'
  ```
  A run prints the same block at start as `[build][...]` lines, ahead of the `[Runtime]`
  host and device lines; `has_runtime_info_display = no` in the cfg silences both.
  A toolchain that also emits the literal into `.debug_*` may print it more than
  once; the copies are identical.
* **Uncommitted changes**: `make snapshot_diff=1` (off by default) also embeds the
  working-tree diff:
  ```bash
  strings <exe> | sed -n '/^build\.code-changes\.begin :$/,/^build\.code-changes\.end   :$/p'
  ```
  Scope is `*.c *.h *.cxx *.hpp *.cpp *.cu`, `Makefile` and `3x3-C/Makefile`, tracked in
  `HEAD` only; anything excluded is counted in the payload. `rev`'s `-dirty` suffix and
  the `dirty=` counts use the same scope.
* **Running with GoSPL**: set `surface_process_option = 11` in the cfg and use
  the generated wrapper; in the Docker image the environment is already active.
  ```bash
  conda activate gospl
  ./dynearthsol-gospl your_input.cfg
  ```
  See `gospl_driver/README.md` and `gospl_driver/examples/` for the details.

# Visualizing DES outputs
* Run `2vtk.py modelname` to convert the binary output to VTK files; `2vtk.py -h`
  lists the options, among them markers (`-m`), principal stresses (`-p`) and
  full tensors (`-t`).
* What a frame contains is set in the input file (`has_marker_output`,
  `is_outputting_averaged_fields`, ... in the `[sim]` section of
  `examples/defaults.cfg`), not by editing sources.
* With `hdf5=1` the output is already `.vtkhdf`, which
  [ParaView](https://www.paraview.org/download/) opens directly; `2vtk.py -u` adds
  the derived fields to it in place. Plain VTK files open in ParaView or
  [VisIt](https://visit-dav.github.io/visit-website/).

# Contributing and bug reports

Bug reports, comments and suggestions are welcome on the
[issue tracker](https://github.com/GeoFLAC/DynEarthSol/issues); a template
asks for what a fix needs. Development conventions -- building, regression
tests, commit messages, pull requests, releases -- are in
[CONTRIBUTING.md](CONTRIBUTING.md), and user-visible changes per release in
[CHANGELOG.md](CHANGELOG.md).

# Citing

Cite the [DynEarthSol v2.0 paper](https://doi.org/10.5194/egusphere-2026-2922)
and the version you used. [CITATION.cff](CITATION.cff) carries the version DOI,
and GitHub renders it under *Cite this repository*; the badge above is the
concept DOI that always resolves to the latest release.

# License

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT / X Windows System license. See
[LICENSE](https://github.com/GeoFLAC/DynEarthSol/blob/master/LICENSE) for the full text.

The files under `3x3-C/`, `tetgen/` and `triangle/` are distributed under their
own licenses. The submodules `knn-bvh`, `mmg` and `nanoflann` are separate
projects with their own.

