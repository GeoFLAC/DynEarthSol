<!-- Citation & License -->
[![DOI](https://zenodo.org/badge/838469917.svg)](https://doi.org/10.5281/zenodo.20293557)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<br> <!-- Core Platforms & Compilers -->
[![Basic build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/basic-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/basic-build.yml)
[![macOS build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/macos-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/macos-build.yml)
[![nvc build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/nvc-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/nvc-build.yml)
<br> <!-- External Modules & Integrations -->
[![Exodus build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/exodus-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/exodus-build.yml)
[![MMG build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/mmg-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/mmg-build.yml)
[![GoSPL build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/gospl-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/gospl-build.yml)

# Overview

DynEarthSol3D, DES3D in short, is a finite element code that solves the momentum balance and 
the heat transfer in Lagrangian form using unstructured meshes. It can be
used to study the long-term deformation of Earth's lithosphere and problems
alike.

# Building DES3D

## Getting the Source Code
This repository uses Git submodules (e.g., `knn-bvh` for GPU version) to manage certain internal libraries. By default, the Makefile automatically prepares these submodules, fetching them via the internet before the build process begins. 

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

Then build:

```bash
make ndims=2    # 2D executable, dynearthsol2d
make            # 3D executable, dynearthsol3d
```

`make config` prints the compiler and every dependency path the build resolved
to. It is the quickest way to see that something was found somewhere you did not
expect, and the most useful thing to paste into a bug report.

In more detail:
* You will need a recent C++ compiler that supports C++11 standard. (GNU g++
  4.4 or newer version will suffice.)
* You will need a recent version of `Boost::Program_options` library (1.42 or
  newer version). Packaged versions are found automatically. On macOS the search
  order is the Homebrew `boost` keg, then MacPorts, then the Homebrew prefix
  itself, then an activated conda environment; elsewhere it is the conda
  environment, then the system Boost the compiler finds on its own. A prefix is
  accepted only if it holds both the headers and the library, so a headers-only
  tree is skipped rather than selected and then failed at the link.
  * Package managers come before conda deliberately. A shell that
    auto-activates `base` has conda on for everything, not just for the work
    that wants it, so it is the weaker signal of intent — but a conda Boost is
    still used when it is the only one installed.
  * To build Boost yourself instead:
    * Download the source code from www.boost.org
    * In the untarred source directory, run `./bootstrap.sh`
    * In the same directory, run `./b2 --with-program_options -q` to build
       the library.
    * Then point the build at it: `make BOOST_ROOT_DIR=/path/to/boost`. A build
      directory is recognised by its `stage/` subdirectory, so the untarred
      source directory is the right thing to name, not `stage/lib`.
* You will need Python 2.6+ or 3.2+ and the Numpy package.
* **macOS users**: Apple clang has no built-in OpenMP, so the runtime has to
  come from somewhere; `brew install libomp` is the whole answer for most
  people. See [OpenMP on macOS](#llvm) for the search order, the overrides, and
  how to build LLVM OpenMP from source when a package manager is not an option.
### Submodules
* [knn-bvh](https://github.com/GeoFLAC/knn-bvh): A library for K-Nearest Neighbors (KNN) searches utilizing Bounding Volume Hierarchies (BVH). This submodule provides GPU acceleration for intensive spatial queries (e.g., particle locating, mesh interpolation, or contact detection)
  * Required when compiling with `openacc=1` (GPU acceleration enabled).

### Optional packages
* [Exodus](https://github.com/gsjaardema/seacas/) for importing a mesh in the ExodusII format
  * Suggested building procedure
    * Run the following in the root directory of DES3D:
      ```BASH
      git clone https://github.com/sandialabs/seacas.git
      cd seacas && export ACCESS=`pwd`
      COMPILER=gnu MATIO=NO GNU_PARALLEL=NO CGNS=NO FMT=NO ./install-tpl.sh
      mkdir build; cd build
      ../cmake-exodus
      make; make install
      ```
  * The above procedure will download and build NetCDF and HDF5; and then build EXODUS.
  * The header files and built shared library will be in `./seacas/include` and `./seacas/lib`. 
* [MMG3D](https://www.mmgtools.org/mmg-remesher-downloads) for mesh optimization during remeshing in three-dimensional models
  * Suggested building procedure
    * Run the following in the root directory of DES3D:
      ```BASH
      git clone https://github.com/MmgTools/mmg.git
      cd mmg; mkdir build; cd build
      cmake ..
      make
      ```
    * The header files and built shared library will be in `mmg/build/include` and `mmg/build/lib`. 
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/) for outputting model results in HDF5-based vtkhdf format, which is compressed (reducing size by up to 50%) and can be visualized directly in Paraview.
  * The HDF5 Library is generally pre-installed on modern computer operating systems; otherwise `brew install hdf5` (macOS), `apt install libhdf5-dev` (Debian/Ubuntu) or `dnf install hdf5-devel` (Fedora/RHEL). `make hdf5=1` locates it on its own, via `pkg-config` and then the usual install layouts.
    * Prefer that search over `which h5cc`: on a Mac that has ever used both Homebrew prefixes, `h5cc` and `pkg-config` on `PATH` are often the `/usr/local` (x86_64) ones while the build is arm64, and their paths link an HDF5 of the wrong architecture. The Makefile picks the prefix matching the host architecture instead.
  * The HDF5-based vtkhdf format follows the data structure of VTK, which can be visualized directly in Paraview. Please refer to the official [VTKHDF File Format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format) documentation for more information.

* [GoSPL](https://gospl.readthedocs.io/) (Global Scalable Paleo Landscape Evolution) for two-way coupling with a surface process model
  * GoSPL handles erosion, deposition, and hillslope diffusion; DES handles tectonics. At each coupling event DES surface velocities are passed to GoSPL, and GoSPL returns the erosion/deposition increment, which is applied to DES surface nodes.
  * Requires the [gospl_extensions](https://github.com/GeoFLAC/gospl_extensions) C++ interface library:
    ```bash
    git clone https://github.com/GeoFLAC/gospl_extensions.git ~/opt/gospl_extensions
    cd ~/opt/gospl_extensions/cpp_interface
    conda activate gospl
    make install-local
    ```
  * Also requires a GoSPL conda environment. Follow [the GoSPL installation procedure](https://gospl.readthedocs.io/en/latest/getting_started/installConda.html).
  * Set `use_gospl = 1` in the Makefile and run `make`. A `dynearthsol-gospl` wrapper script is generated automatically.
  * See `gospl_driver/README.md` for full build, runtime, and coupling details, and `gospl_driver/examples/` for example configs.

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
  `make check-deps` checks the architecture of what it found and says so
  outright, rather than leaving it to the linker.

* **Overrides**, when the runtime is somewhere the search does not look, or when
  you want a specific one. Set these on the command line, or edit them in the
  clang++ branch of the `Makefile`, where the OpenMP search lives:

  ```bash
  make OPENMP_ROOT_DIR=/prefix                        # /prefix/include + /prefix/lib
  make OPENMP_INCLUDE_DIR=... OPENMP_LIB_DIR=...      # when they are not siblings
  ```

  An *exported* `OPENMP_ROOT_DIR` or `BOOST_ROOT_DIR` is deliberately ignored.
  These are DES's own variable names, not ones any package manager or module
  file sets, so an exported value is nearly always a forgotten line in a shell
  profile — and since naming a prefix skips detection, one stale export would
  silently override a perfectly good install and be reported as a missing
  dependency. Pass an empty `BOOST_ROOT_DIR=` on the command line to force
  detection even when the `Makefile` has a prefix written into it.

  The second form is what Homebrew's `llvm` needs if named explicitly, for the
  reason given above.

* **No OpenMP at all** is a supported configuration -- `make openmp=0` builds a
  single-threaded executable, which is also what you want for valgrind. The
  OpenMP runtime calls in the source are guarded, so nothing else is needed.

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

  * macOS / Apple Silicon — OpenMP thread-wait performance

    LLVM `libomp` hardcodes hybrid-CPU detection for all Apple Silicon, setting
    thread *blocktime* to **0 µs**. Threads yield immediately after each parallel
    region instead of spin-waiting, reducing CPU utilisation to ~250% on a 6-core
    Mac vs. ~600% on Linux. DES3D automatically sets `OMP_WAIT_POLICY=active` at startup on macOS (unless
    `OMP_WAIT_POLICY` or `KMP_BLOCKTIME` is already set), restoring near-Linux
    throughput. A notice is printed at startup confirming this.

    To override:
    ```bash
    export OMP_WAIT_POLICY=passive   # yield immediately, lower power
    export KMP_BLOCKTIME=20          # fine-grained control (ms, libomp only)
    ```

## Or, using docker

* Build docker image
  ```bash
  ./build.sh
  ```
* Run docker
  ```bash
  docker run --rm -it dynearthsol/gcc-11 # default compiler
  ```
## Building procedure

Boost and, on macOS, the OpenMP runtime are located automatically, so a plain
`make` is normally the whole procedure. `make config` shows what was found;
everything below is for the cases where the defaults are not what you want.

* Options can be set either on the command line (`make ndims=2 hdf5=1`) or by
  editing the corresponding variable at the top of the `Makefile`. A value given
  on the command line wins.
* `BOOST_ROOT_DIR`, `HDF5_INCLUDE_DIR`, `HDF5_LIB_DIR` and `NVHPC_DIR` are
  declared together in an **Optional paths** block near the top of the `Makefile`;
  all may stay blank, and are only there to force one specific install. Every
  other dependency's path is declared beside the code that uses it — `OPENMP_*` in
  the clang++ branch, and `EXO_*`, `MMG_*` and the GoSPL group with their own
  features — so it is read together with the flags it feeds.
* `make config` prints the resolved compiler, flags and dependency paths.
  `make check-deps` verifies the external libraries without building anything.
* Edit `Makefile`
  * `BOOST_ROOT_DIR` selects a specific Boost, overriding the search. Set it to
  the untarred boost directory if you built `Boost::Program_options` yourself
  following the instructions above; a build directory with a `stage/`
  subdirectory is recognised as such.
  * If importing an exodus mesh:
    * Set `useexo = 1` and `ndims = 3`. Only 3D exodus mesh can be imported.
    * Set `EXO_INCLUDE` and `EXO_LIB_DIR` paths if it differs from the default values.
  * If mesh optimization with mmg is desired for remeshing:
    * Set `usemmg = 1`.
    * Set `MMG_INCLUDE` and `MMG_LIB_DIR` paths if it differs from the default values.
  * If coupling with GoSPL surface processes:
    * Set `use_gospl = 1`.
    * Set `GOSPL_EXT_DIR` and `CONDA_ENV_PATH` if they differ from the defaults (`~/opt/gospl_extensions` and `~/miniconda3/envs/gospl`).
  * If outputing in HDF5-based vtkhdf format:
    * set `hdf5 = 1`. The paths are found automatically on Linux and macOS, so
      normally there is nothing else to set.
    * `HDF5_INCLUDE_DIR` and `HDF5_LIB_DIR` override that search. Set them
      (in the Makefile or on the command line) to build against one specific
      HDF5 -- an HPC module, a conda env, a hand-built copy:
      `make hdf5=1 HDF5_INCLUDE_DIR=/prefix/include HDF5_LIB_DIR=/prefix/lib`.
      Either one alone is enough; the other is still detected.
    * Install python HDF5 lib by `pip install h5py` for further analyzed vtk visualization.
  * OpenMP on macOS needs no setup beyond `brew install libomp`; see
    [OpenMP on macOS](#llvm) for the search order and for `OPENMP_ROOT_DIR` /
    `OPENMP_INCLUDE_DIR` / `OPENMP_LIB_DIR` if you need to name one explicitly.
* Run `make` to build optimized executable.
* Or run `make opt=0` to build a debugging executable.
* Or run `make openmp=0` to build the executable without OpenMP. This is
  necessary to debug the code under valgrind.
* Or run `make opt=-1` to build a memory-specific debugging executable using `-fsanitize=address`, a compiler flag for detacting memory address issues. It can show where the issue occurs and where variables are allocated during execution, without needing additional tools such as GDB or Valgrind. However, valgrind cannot easily coexist with -fsanitize=address. as using both together may cause library-related errors.

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

# enable MMG mesh optimization (requires MMG headers/libs)
make usemmg=1

# enable Exodus input support (requires seacas/exodus libs)
make useexo=1

# enable HDF5-based vtkhdf output support (requires HDF5)
make hdf5=1

# enable GoSPL surface process coupling (requires gospl_extensions and gospl conda env)
make use_gospl=1

# NVHPC/profiler build (uses nvc++ when set)
make nprof=1

# OpenACC build (NVHPC compiler)
make openacc=1

# OpenACC targeting a specific GPU (else from nvidia-smi, default 80; see make config)
make openacc=1 GPU_CC=90
```

# Running DES3D
* Execute `dynearthsol2d [inputfile: examples/defaults.cfg by default]`.
* Pay attention to any warnings. For instance, if a warning about potential 
  race condition is printed on screen, do follow the given suggestions.
* Several example input files are provided under `examples/` directory. The
  format of the input file is described in `examples/defaults.cfg`.
* Use the [simple input file generator](https://geoflac.github.io/des-inputgen) to create input files for your simulations. This tool provides an easy-to-use interface for generating configuration files tailored to your specific needs.
* Benchmark cases with analytical solution can be found under `benchmarks/`
  directory.
* Execute the executable with `-h` flag to see the available input parameters
  and their descriptions.
* **Running with GoSPL**: Use the auto-generated wrapper script and set `surface_process_option = 11` in your config file:
  ```bash
  conda activate gospl
  ./dynearthsol-gospl your_input.cfg
  ```
  See `gospl_driver/README.md` and `gospl_driver/examples/` for configuration details.

# Visualizing DES3D outputs
* Run `2vtk.py [modelname: 'results' by default]` to convert the binary output to VTK files.
* Execute `2vtk.py -h` to see more usage information.
* Some of the simulation outputs might be disabled. Edit `2vtk.py` and
  `output.cxx` to disable/enable them.
* Plot the VTK files with [Paraview](https://www.paraview.org/download/) or [Visit](https://visit-dav.github.io/visit-website/).

# Bug reports
      
Bug reports, comments, and suggestions are always welcome. The best 
channel is to create an issue on the Issue Tracker here:
  https://github.com/GeoFLAC/DynEarthSol/issues

# License

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT / X Windows System license. See
[LICENSE](https://github.com/GeoFLAC/DynEarthSol/blob/master/LICENSE) for the full text.

The files under the subdirectories `3x3-C/`, `nanoflann/`, `tetgen/`
and `triangles/` are distributed by their own license(s).

