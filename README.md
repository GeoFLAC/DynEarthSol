[![Basic build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/basic-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/basic-build.yml)
[![Exodus build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/exodus-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/exodus-build.yml)
[![MMG build](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/mmg-build.yml/badge.svg)](https://github.com/GeoFLAC/DynEarthSol/actions/workflows/mmg-build.yml)

# Overview

DynEarthSol3D, DES3D in short, is a finite element code that solves the momentum balance and 
the heat transfer in Lagrangian form using unstructured meshes. It can be
used to study the long-term deformation of Earth's lithosphere and problems
alike.

# Building DES3D
## Requirements
* You will need a recent C++ compiler that supports C++11 standard. (GNU g++
  4.4 or newer version will suffice.)
* You will need a recent version of `Boost::Program_options` library (1.42 or
  newer version). Instructions for building the library:
  * Download the source code from www.boost.org
  * In the untarred source directory, run `./bootstrap.sh`
  * In the same directory, run `./b2 --with-program_options -q` to build
     the library.
* You will need Python 2.6+ or 3.2+ and the Numpy package.
* **macOS users**: For OpenMP support on macOS, see the [LLVM OpenMP library build instructions](#llvm) below.
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
  * The HDF5 Library is generally pre-installed on modern computer operating systems. User can use `which h5cc` to find the path to the HDF5 Library.
  * The HDF5-based vtkhdf format follows the data structure of VTK, which can be visualized directly in Paraview. Please refer to the official [VTKHDF File Format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format) documentation for more information.

<div id="llvm"></div>

* [LLVM](https://github.com/llvm/llvm-project) OpenMP library for macOS requires special setup due to Apple Clang lacking built-in OpenMP. 
  * Suggested building procedure
    * Download LLVM CMake Modules and OpenMP. (LLVM OpenMP 15.0.7 or newer version will suffice.)
      ```BASH
      mkdir -p external && cd external
      # Download CMake modules
      curl -L https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.7/cmake-15.0.7.src.tar.xz -o cmake-15.0.7.src.tar.xz
      tar xf cmake-15.0.7.src.tar.xz

      # Download OpenMP source
      curl -L https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.7/openmp-15.0.7.src.tar.xz -o openmp-15.0.7.src.tar.xz
      tar xf openmp-15.0.7.src.tar.xz
      ```
    * Build OpenMP
      ```BASH
      # Configure and build for ARM64 (Apple Silicon) or x86_64 (Intel Mac)
      mkdir -p openmp-15.0.7.src/build && cd openmp-15.0.7.src/build
      cmake -DCMAKE_INSTALL_PREFIX=$(pwd)/../../openmp-install \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_MODULE_PATH=$(pwd)/../../cmake-15.0.7.src/Modules \
            -DCMAKE_OSX_ARCHITECTURES=arm64 \
            -DLIBOMP_INSTALL_ALIASES=OFF \
            ..
      make -j4 && make install && cd ../../
      ```
    * Installed LLVM OpenMP will be in `external/openmp-install`.
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
* Edit `Makefile` 
  * Modify `BOOST_ROOT_DIR` if you manually built or installed 
  boost library.
    * If you followed the instructions above to build 
  `Boost::Program_options` library, set `BOOST_ROOT_DIR` to the untarred boost
  directory.
  * If importing an exodus mesh:
    * Set `useexo = 1` and `ndims = 3`. Only 3D exodus mesh can be imported.
    * Set `EXO_INCLUDE` and `EXO_LIB_DIR` paths if different from the default values.
  * If mesh optimization with mmg is desired for remeshing:
    * Set `usemmg = 1`.
    * Set `MMG_INCLUDE` and `MMG_LIB_DIR` paths if different from the default values.
  * If outputing in HDF5-based vtkhdf format:
    * set `hdf5 = 1`.
    * set `HDF5_INCLUDE_DIR` to the HDF5 header file directory.
    * set `HDF5_LIB_DIR` to the HDF5 library directory.
    * Install python HDF5 lib by `pip install h5py` for further analyzed vtk visualization.
  * if enabling openMP on macOS:
    *  set `OPENMP_ROOT_DIR` path if different from the default values.
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

# NVHPC/profiler build (uses nvc++ when set)
make nprof=1

# OpenACC build (NVHPC compiler)
make openacc=1
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

