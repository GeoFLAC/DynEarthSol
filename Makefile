# -*- Makefile -*-
#
# Makefile for DynEarthSol3D
#
# Author: Eh Tan <tan2@earth.sinica.edu.tw>
#
## Build notes
## - Run simply `make` to build the optimized production executable.
## - For a debugging build, run for example: `make opt=0 openmp=0`.
## Common configuration variables (set on the make command line or edit below):
##  - ndims = 3 : build 3D code; set to 2 for the 2D code.
##  - opt = 0..3 : optimization level. 0 = debugging (no optimizations),
##       1 = low optimization, 2 = default optimized build, 3 = aggressive
##       optimizations (-march=native, -O3, etc.).
##  - openacc = 1 : enable OpenACC compilation (NVHPC).
##  - openmp = 1 : enable OpenMP parallelization.
##  - nprof = 1 : enable NVHPC nprof profiling build (uses nvc++ when set),
##       1 = main dynearthsol loop profiling, 2 = detailed profiling.
##  - gprof = 1 : enable GNU gprof instrumentation (-pg).
##  - usemmg = 1 : enable MMG mesh optimization support (requires MMG headers/libs).
##  - hdf5 = 1 : enable HDF5-based vtkhdf output support (requires hdf5).
##  - adaptive_time_step = 1 : enable adaptive time stepping.
##  - use_R_S = 1 : enable Rate-and-State friction (requires adaptive_time_step).
##  - useexo = 1 : enable ExodusII import support (3D only; requires seacas/exodus libs).

ndims = 3
opt = 2
openacc = 0
openmp = 1
nprof = 0
gprof = 0
usemmg = 0
adaptive_time_step = 0
use_R_S = 0
useexo = 0
use_gospl = 0
hdf5 = 0
nofma = 0   # disable FMA instructions when using nvc++, may help if using mixed precision

ifeq ($(ndims), 2)
	useexo = 0    # for now, can import only 3d exo mesh
endif

ifneq ($(adaptive_time_step), 1)
	use_R_S = 0   # Rate - State friction law only works with adaptive time stepping technique
endif

OSNAME := $(shell uname -s)

## Prefixes every dependency search below starts from. Chosen by host
## architecture, never from PATH: a Mac that has run both Homebrews keeps the
## other one's pkg-config and headers on PATH, and building against those links
## the wrong architecture. `brew --prefix <formula>` is not used either -- it
## prints the would-be path and exits 0 for a formula that is not installed.
ifeq ($(OSNAME), Darwin)
	ifeq ($(shell uname -m), arm64)
		BREW_PREFIX ?= /opt/homebrew
	else
		BREW_PREFIX ?= /usr/local
	endif
	MACPORTS_PREFIX ?= /opt/local
endif

## $(call firstdir,<file>,<candidate dirs>) -- first directory in the list that
## really holds <file>, else empty. Globs allowed in both arguments
## (`lib/clang/*/include`, `libhdf5.*`); the candidate list is a priority list.
##
## $(wildcard) expands the globs, `test -f` decides. Not `ls` and not $(wildcard)
## alone: both accept a DANGLING SYMLINK, and $(BREW_PREFIX)/opt/<formula> is a
## symlink into the Cellar, so a half-removed keg would stop the search on a dud
## instead of falling through to a provider that works.
firstdir = $(patsubst %/,%,$(dir $(call firstfile,$(1),$(2))))
firstfile = $(firstword $(shell for f in $(wildcard $(addsuffix /$(1),$(2))); do \
	test -f "$$f" && { echo $$f; break; }; done))

## Drop a candidate whose prefix is empty, so an unset CONDA_PREFIX contributes
## nothing rather than a bare, plausible-looking /include or /lib.
under = $(if $(strip $(1)),$(addprefix $(strip $(1)),$(2)))

## Select C++ compiler and set paths to necessary libraries
ifeq ($(openacc), 1)
	CXX = nvc++
	suffix = .gpu
else
	ifneq ($(nprof), 0)
		CXX = nvc++
	else
		# Select compiler based on platform
		ifeq ($(OSNAME), Darwin)
			# clang++ is the default optimized compiler for macOS
			CXX = clang++
		else
			CXX = g++
		endif
	endif
endif
CXX_BACKEND = ${CXX}

## path to HDF5's base directory. Leave both blank to auto-detect; set them
## (here or on the command line) to force one specific install -- an HPC
## module, a conda env, a hand-built HDF5 -- which skips detection entirely.
## The layouts the search below knows about, if you would rather set them:
##   Debian/Ubuntu  /usr/include/hdf5/serial  /usr/lib/<arch>-linux-gnu/hdf5/serial
##   Fedora/RHEL    /usr/include              /usr/lib64
##   Homebrew       $(brew --prefix hdf5)/include  $(brew --prefix hdf5)/lib
HDF5_INCLUDE_DIR = #/usr/include/hdf5/serial
HDF5_LIB_DIR = #/usr/lib/x86_64-linux-gnu/hdf5/serial

ifeq ($(hdf5), 1)
## Every stage below fills only what is still blank, so a path set by hand is
## never overwritten, and setting just one of the two still leaves the other
## detected instead of empty (a blank -L would eat the -lhdf5 that follows it).
HDF5_DETECT :=
ifeq ($(strip $(HDF5_INCLUDE_DIR)),)
	HDF5_DETECT := yes
endif
ifeq ($(strip $(HDF5_LIB_DIR)),)
	HDF5_DETECT := yes
endif
ifeq ($(HDF5_DETECT), yes)
	ifeq ($(OSNAME), Darwin)
		## Use the architecture-matched Homebrew prefix (see BREW_PREFIX above)
		## and that prefix's own pkg-config, not whichever one PATH offers.
		HDF5_BREW = $(BREW_PREFIX)
		PKGCONFIG := $(BREW_PREFIX)/bin/pkg-config
	else
		HDF5_BREW =
		PKGCONFIG := pkg-config
	endif

	## Ask pkg-config for directories rather than parsing --cflags/--libs:
	## h5cc -show and --libs mix in inline .a paths and -lsz/-lz/-ldl, which
	## are painful to filter and differ per build. Debian and Ubuntu name the
	## serial module hdf5-serial; everyone else ships plain hdf5.
	HDF5_PC := $(shell for m in hdf5 hdf5-serial; do \
		$(PKGCONFIG) --exists $$m 2>/dev/null && { echo $$m; break; }; done)
	ifneq ($(HDF5_PC),)
		ifeq ($(strip $(HDF5_INCLUDE_DIR)),)
			HDF5_INCLUDE_DIR := $(shell $(PKGCONFIG) --variable=includedir $(HDF5_PC))
		endif
		ifeq ($(strip $(HDF5_LIB_DIR)),)
			HDF5_LIB_DIR := $(shell $(PKGCONFIG) --variable=libdir $(HDF5_PC))
		endif
	endif

	## Without a .pc file, probe the layouts real installs use -- Fedora and
	## RHEL need this branch unconditionally, as their hdf5-devel ships no
	## pkg-config file at all. Debian splits its serial build across two trees
	## (/usr/include/hdf5/serial and .../hdf5/serial under the multiarch libdir),
	## so the header and the library are searched separately, not under one
	## prefix. /usr/local comes before the system directories so a hand-built
	## HDF5 installed there wins over the distro package.
	HDF5_INC_DIRS = $(call under,$(HDF5_BREW),/opt/hdf5/include) \
		/usr/include/hdf5/serial \
		/usr/local/include /usr/include
	HDF5_LIB_DIRS = $(call under,$(HDF5_BREW),/opt/hdf5/lib) \
		/usr/lib/*-linux-gnu/hdf5/serial \
		/usr/local/lib /usr/local/lib64 \
		/usr/lib/*-linux-gnu \
		/usr/lib64 /usr/lib
	ifeq ($(strip $(HDF5_INCLUDE_DIR)),)
		HDF5_INCLUDE_DIR := $(call firstdir,hdf5.h,$(HDF5_INC_DIRS))
	endif
	ifeq ($(strip $(HDF5_LIB_DIR)),)
		HDF5_LIB_DIR := $(call firstdir,libhdf5.*,$(HDF5_LIB_DIRS))
	endif

	## Checked separately: finding one without the other is still unbuildable,
	## and failing here beats a confusing link error hundreds of lines later.
	## Indented with SPACES, not a tab: a tab-indented $(error) this early
	## reads as a recipe line and make reports "commands commence before first
	## target" instead of the message below.
	ifeq ($(strip $(HDF5_INCLUDE_DIR)),)
        $(error hdf5=1 but no HDF5 headers found. Install HDF5 (macOS: brew install hdf5; \
Debian/Ubuntu: apt install libhdf5-dev; Fedora/RHEL: dnf install hdf5-devel), \
or set the paths by hand: \
make hdf5=1 HDF5_INCLUDE_DIR=/prefix/include HDF5_LIB_DIR=/prefix/lib)
	endif
	ifeq ($(strip $(HDF5_LIB_DIR)),)
        $(error hdf5=1: found HDF5 headers in $(HDF5_INCLUDE_DIR) but no libhdf5. \
Set HDF5_LIB_DIR to the directory holding libhdf5, \
e.g. make hdf5=1 HDF5_LIB_DIR=/prefix/lib)
	endif
endif
endif

## path to cuda's base directory
NVHPC_DIR = # /cluster/nvidia/hpc_sdk/Linux_x86_64/21.2

## path to Boost's base directory. Blank = auto-detect; set it here or on the
## command line to force one install (an HPC module, a conda env, a hand-built
## Boost), which skips detection.
##
## Plain `=`, not `?=`: an EXPORTED BOOST_ROOT_DIR is ignored on purpose. It is
## DES's own variable name, so an exported value is almost always a stale shell
## profile line -- and since setting it skips detection, one forgotten export
## defeats a good install and reports it as "Boost missing".
BOOST_ROOT_DIR = # /path/to/boost

## Candidate prefixes, in priority order. macOS has no system Boost, so /usr is
## never the answer there and these tiers are what remove the old requirement to
## pass BOOST_ROOT_DIR=$(brew --prefix boost) by hand.
##
## Package managers come BEFORE conda: a shell that auto-activates `base` has
## conda on for everything, not just the work that wants it, so it is the weaker
## signal of intent. A conda Boost is still used when it is the only one.
## The `opt/boost` keg beats the bare $(BREW_PREFIX) farm -- it pins the version.
## Linux keeps the tiers it had: conda, else the /usr default below.
ifeq ($(OSNAME), Darwin)
	BOOST_PREFIXES = $(BREW_PREFIX)/opt/boost $(MACPORTS_PREFIX) $(BREW_PREFIX) \
		$(CONDA_PREFIX)
else
	BOOST_PREFIXES = $(CONDA_PREFIX)
endif

## Both halves required: a headers-only tree (extracted tarball, half-removed
## keg, conda's libboost-headers without libboost) would otherwise be selected
## and then fail at the link. The `libboost_program_options.*` glob matches
## exactly what -lboost_program_options can link, and not MacPorts' `-mt` variant.
## `override` so a command-line `BOOST_ROOT_DIR=` still falls through to
## detection -- a command-line value otherwise outranks even the /usr fallback.
ifndef BOOST_ROOT_DIR
	override BOOST_ROOT_DIR := $(shell for d in $(BOOST_PREFIXES); do \
		test -f $$d/include/boost/program_options.hpp || continue; \
		set -- $$d/lib/libboost_program_options.*; \
		test -f "$$1" || continue; \
		echo $$d; break; done)
endif

## Nothing found: keep the historical /usr default so that `make clean`, `make
## config` and a hand-set BOOST_LIB_DIR all still work, and record that this is
## a failed search rather than a deliberate choice. On macOS the difference is
## decisive -- there is no system Boost under /usr (this macOS has no
## /usr/include at all), so the fallback there always means "not found", and
## check-deps says so instead of letting the link produce a wall of undefined
## symbols.
ifndef BOOST_ROOT_DIR
	override BOOST_ROOT_DIR = /usr
	## Only on macOS does the fallback mean "the search failed" -- on Linux /usr
	## is the normal, correct answer and the compiler finds Boost there on its
	## own. Keeping that distinction is what lets check-deps be decisive on
	## macOS without misreporting a healthy Linux build.
	ifeq ($(OSNAME), Darwin)
		BOOST_NOT_FOUND := yes
	endif
endif

########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

BOOST_LDFLAGS = -lboost_program_options
# Detect multiarch tuple (e.g. x86_64-linux-gnu on Debian/Ubuntu); empty on macOS/other.
# Skipped on macOS, where there is no multiarch layout and `gcc` is Apple clang,
# which has no -print-multiarch: two shell-outs per make invocation for nothing.
ifneq ($(OSNAME), Darwin)
	MULTIARCH := $(shell gcc -print-multiarch 2>/dev/null)
	ifeq ($(MULTIARCH),)
		MULTIARCH := $(shell dpkg-architecture -qDEB_HOST_MULTIARCH 2>/dev/null)
	endif
endif
ifdef BOOST_ROOT_DIR
	# check existence of stage/ directory
	has_stage_dir = $(wildcard $(BOOST_ROOT_DIR)/stage)
	ifeq (, $(has_stage_dir))
		# no stage dir, BOOST_ROOT_DIR is the installation directory
		ifeq ($(BOOST_ROOT_DIR), /usr)
			# Conda environments can disrupt gcc's default search order, so explicitly
			# add the multiarch path before /usr/include (keeps bits/stdint-least.h
			# resolvable) and use the multiarch lib dir where boost is installed.
			ifneq ($(MULTIARCH),)
				BOOST_CXXFLAGS = -I/usr/include/$(MULTIARCH) -I/usr/include
				BOOST_LIB_DIR = /usr/lib/$(MULTIARCH)
			else
				BOOST_LIB_DIR = /usr/lib
			endif
		else
			BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)/include
			BOOST_LIB_DIR = $(BOOST_ROOT_DIR)/lib
		endif
	else
		# with stage dir, BOOST_ROOT_DIR is the build directory
		BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)
		BOOST_LIB_DIR = $(BOOST_ROOT_DIR)/stage/lib
	endif
	## Full path on macOS, for the same reason as libomp above: a -L dir also
	## answers the implicit -lc++, and a conda prefix ships its own libc++.dylib
	## (left on -L, the executable loads two libc++ images and says nothing).
	## Only when the .dylib exists -- a static-only Boost still needs -l -- and
	## via firstfile so a dangling Cellar symlink falls back to -l too.
	## Linux keeps -L/-l, right for .so versioning and the stage layout.
	BOOST_DYLIB = $(call firstfile,libboost_program_options.dylib,$(BOOST_LIB_DIR))
	ifeq ($(OSNAME)$(if $(BOOST_DYLIB),yes,), Darwinyes)
		BOOST_LDFLAGS = $(BOOST_DYLIB)
	else
		BOOST_LDFLAGS += -L$(BOOST_LIB_DIR)
	endif
	ifneq ($(OSNAME), Darwin)
		# Apple's ld supports -rpath only in the comma form, so Darwin gets its
		# rpath from the install_name_tool step after the link instead.
		BOOST_LDFLAGS += -Wl,-rpath=$(BOOST_LIB_DIR)
	endif
endif

ifeq ($(useexo), 1)
	# path to exodus header files
	EXO_INCLUDE = ./seacas/include

	# path of exodus library files, if not in standard system location
	EXO_LIB_DIR = ./seacas/lib

	EXO_CXXFLAGS = -I$(EXO_INCLUDE) -DUSEEXODUS
	EXO_LDFLAGS = -L$(EXO_LIB_DIR) -lexodus
	ifneq ($(OSNAME), Darwin)  # Apple's ld doesn't support -rpath
		EXO_LDFLAGS += -Wl,-rpath=$(EXO_LIB_DIR)
	endif
endif

ifeq ($(usemmg), 1)
	# path to MMG3D header files
	MMG_INCLUDE = ./mmg/build/include

	# path of MMG3D library files, if not in standard system location
	MMG_LIB_DIR = ./mmg/build/lib

	MMG_CXXFLAGS = -I$(MMG_INCLUDE) -DUSEMMG
	ifeq ($(ndims), 3)	
		MMG_LDFLAGS = -L$(MMG_LIB_DIR) -lmmg3d
	else
		MMG_LDFLAGS = -L$(MMG_LIB_DIR) -lmmg2d
	endif
	ifneq ($(OSNAME), Darwin)  # Apple's ld doesn't support -rpath
		MMG_LDFLAGS += -Wl,-rpath=$(MMG_LIB_DIR)
	endif

	MMG_LIB = $(MMG_LIB_DIR)/libmmg$(ndims)d.a
endif

ifeq ($(use_gospl), 1)
	# GoSPL extensions library configuration
	# Path to locally built gospl_extensions (clone from GitHub)
	GOSPL_EXT_DIR = $(HOME)/opt/gospl_extensions
	GOSPL_INCLUDE = $(GOSPL_EXT_DIR)/include
	GOSPL_LIB_DIR = $(GOSPL_EXT_DIR)/lib
	GOSPL_PYTHONPATH = $(GOSPL_EXT_DIR)/cpp_interface
	
	# Export PYTHONPATH for GoSPL extensions
    	export PYTHONPATH := $(GOSPL_PYTHONPATH):$(PYTHONPATH)
  
	GOSPL_CXXFLAGS = -I$(GOSPL_INCLUDE) -DHAS_GOSPL_CPP_INTERFACE -Igospl_driver
	GOSPL_LDFLAGS = -L$(GOSPL_LIB_DIR) -lgospl_extensions
	ifneq ($(OSNAME), Darwin)  # Apple's ld doesn't support -rpath
		GOSPL_LDFLAGS += -Wl,-rpath=$(GOSPL_LIB_DIR)
	endif
	
	# Python include/lib dirs — override on the command line for non-conda builds:
	#   make use_gospl=1 PYTHON_VERSION=3.12 PYTHON_INCLUDE_DIR=/usr/include/python3.12 PYTHON_LIB_DIR=/usr/lib/x86_64-linux-gnu
	CONDA_ENV_PATH = $(HOME)/miniconda3/envs/gospl
	PYTHON_VERSION     ?= 3.11
	PYTHON_INCLUDE_DIR ?= $(CONDA_ENV_PATH)/include/python$(PYTHON_VERSION)
	PYTHON_LIB_DIR     ?= $(CONDA_ENV_PATH)/lib

	GOSPL_CXXFLAGS += -I$(PYTHON_INCLUDE_DIR)
	GOSPL_LDFLAGS += -lpython$(PYTHON_VERSION) -L$(PYTHON_LIB_DIR)
	ifneq ($(OSNAME), Darwin)
		GOSPL_LDFLAGS += -Wl,-rpath=$(PYTHON_LIB_DIR)
	endif
endif


ifneq (, $(findstring clang++, $(CXX)))
	CXX_SUPPORTED = yes
	CXXFLAGS = -g -std=c++0x -DGPP1X
	LDFLAGS = -lm
	TETGENFLAG = 
	
	# macOS needs extra headerpad for install_name_tool
	ifeq ($(OSNAME), Darwin)
		LDFLAGS += -Wl,-headerpad_max_install_names
	endif 

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	else ifeq ($(opt), 3)
		CXXFLAGS += -march=native -O3 -ffast-math -funroll-loops
	else # debugging
		CXXFLAGS += -O0 -Wall -Wno-unused-variable -Wno-unused-function -Wno-unknown-pragmas
		ifeq ($(opt), -1)
			CXXFLAGS += -fsanitize=address
			LDFLAGS += -fsanitize=address
		endif
	endif
 
	ifeq ($(openmp), 1)
	ifeq ($(OSNAME), Darwin)
		## Apple clang ships no OpenMP runtime, so one has to be found. Set
		## OPENMP_ROOT_DIR to force a prefix (its include/ and lib/), or
		## OPENMP_INCLUDE_DIR / OPENMP_LIB_DIR when they are not siblings.
		## Plain `=`: an exported value is ignored, as for BOOST_ROOT_DIR.
		OPENMP_ROOT_DIR = # /path/to/openmp
		ifneq ($(strip $(OPENMP_ROOT_DIR)),)
			ifeq ($(strip $(OPENMP_INCLUDE_DIR)),)
				OPENMP_INCLUDE_DIR := $(OPENMP_ROOT_DIR)/include
			endif
			ifeq ($(strip $(OPENMP_LIB_DIR)),)
				OPENMP_LIB_DIR := $(OPENMP_ROOT_DIR)/lib
			endif
		endif

		## Candidate lists in priority order, named once and reused by the
		## search, the check-deps `Searched:` line and `make config`, so the
		## three cannot disagree. external/openmp-install is first: it is where
		## this project's build-from-source instructions put LLVM OpenMP.
		## Header and library are searched separately because Homebrew's llvm
		## keeps omp.h under lib/clang/<ver>/include, not beside the library.
		OPENMP_INC_DIRS = external/openmp-install/include \
			$(BREW_PREFIX)/opt/libomp/include \
			$(BREW_PREFIX)/opt/llvm/lib/clang/*/include \
			$(MACPORTS_PREFIX)/include/libomp \
			$(call under,$(CONDA_PREFIX),/include)
		OPENMP_LIB_DIRS = external/openmp-install/lib \
			$(BREW_PREFIX)/opt/libomp/lib \
			$(BREW_PREFIX)/opt/llvm/lib \
			$(MACPORTS_PREFIX)/lib/libomp \
			$(call under,$(CONDA_PREFIX),/lib)

		ifeq ($(strip $(OPENMP_INCLUDE_DIR)),)
			OPENMP_INCLUDE_DIR := $(call firstdir,omp.h,$(OPENMP_INC_DIRS))
		endif
		ifeq ($(strip $(OPENMP_LIB_DIR)),)
			OPENMP_LIB_DIR := $(call firstdir,libomp.dylib,$(OPENMP_LIB_DIRS))
		endif

		## Flags only when both halves were found, so a half-detected runtime
		## cannot produce a truncated command line; check-deps does the
		## reporting, which keeps `make clean` and `make config` working.
		ifneq ($(strip $(OPENMP_INCLUDE_DIR)),)
		ifneq ($(strip $(OPENMP_LIB_DIR)),)
			## rpath must be absolute, so a relative OPENMP_LIB_DIR
			## (external/... above) is anchored at the build directory. It is
			## needed when the dylib's install name is @rpath-relative --
			## conda's libomp is, Homebrew's is not.
			ifeq ($(filter /%,$(OPENMP_LIB_DIR)),$(OPENMP_LIB_DIR))
				OPENMP_RPATH := $(OPENMP_LIB_DIR)
			else
				OPENMP_RPATH := $(CURDIR)/$(OPENMP_LIB_DIR)
			endif
			## -idirafter, not -I: placed after the system directories it can
			## only supply a header the toolchain lacks, i.e. omp.h. With -I,
			## Homebrew llvm's omp.h dir -- a compiler *resource* dir holding
			## its own stdint.h -- breaks the SDK ("unknown type name
			## 'uint8_t'" from sys/resource.h). Conda shadows the same way.
			##
			## Full path, not -L/-lomp: every -L dir also answers the implicit
			## -lc++/-lz/-lm, and a conda prefix ships its own libc++.dylib.
			CXXFLAGS += -Xpreprocessor -fopenmp -idirafter $(OPENMP_INCLUDE_DIR)
			LDFLAGS += $(OPENMP_LIB_DIR)/libomp.dylib
		endif
		endif
	else
		## Non-Apple clang drives OpenMP itself: -fopenmp both enables the
		## pragmas and links the runtime, with no paths to find.
		CXXFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif
	endif

else ifneq (, $(findstring g++, $(CXX_BACKEND))) # if using any version of g++
	CXX_SUPPORTED = yes
	CXXFLAGS = -g -std=c++0x
	LDFLAGS = -lm
	TETGENFLAG = -Wno-unused-but-set-variable -Wno-int-to-pointer-cast

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	else ifeq ($(opt), 3) # experimental, use at your own risk :)
		CXXFLAGS += -march=native -O3 -ffast-math -funroll-loops
	else # debugging flags
		CXXFLAGS += -O0 -Wall -Wno-unused-variable -Wno-unused-function -Wno-unknown-pragmas -fbounds-check -ftrapv
		ifeq ($(opt), -1)
			CXXFLAGS += -fsanitize=address
			LDFLAGS += -fsanitize=address
		endif
	endif

	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif

	ifeq ($(gprof), 1)
		CXXFLAGS += -pg
		LDFLAGS += -pg
	endif

	GCCVERSION = $(shell $(CXX) --version | grep g++ | sed 's/^.* //g' | cut -d. -f1)

	ifeq ($(shell expr $(GCCVERSION) \> 10), 1)
		CXXFLAGS += -DGPP1X
	endif

else ifneq (, $(findstring icpc, $(CXX_BACKEND))) # if using intel compiler, tested with v14
	CXX_SUPPORTED = yes
	CXXFLAGS = -g -std=c++0x
	LDFLAGS = -lm

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	else ifeq ($(opt), 3) # experimental, use at your own risk :)
		CXXFLAGS += -fast -fast-transcendentals -fp-model fast=2
	else # debugging flags
		CXXFLAGS += -O0 -check=uninit -check-pointers=rw -check-pointers-dangling=all -fp-trap-all=all
	endif

	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif

else ifneq (, $(findstring nvc++, $(CXX)))
	CXX_SUPPORTED = yes
	CXXFLAGS = -g -Minfo=mp,accel
	LDFLAGS =
	TETGENFLAGS = 

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	endif

	ifeq ($(nofma), 1)
		CXXFLAGS += -mno-fma
	endif

	ifeq ($(openacc), 1)
		CXXFLAGS += -acc=gpu -cuda -DACC
		LDFLAGS += -acc=gpu -gpu=mem:managed -cuda
		# CXXFLAGS += -acc=gpu -Mcuda -DACC
		# LDFLAGS += -acc=gpu -gpu=managed -Mcuda
		ifeq ($(nofma), 1)
			CXXFLAGS += -gpu=mem:managed,nofma
			# CXXFLAGS += -gpu=managed,nofma
		else
			CXXFLAGS += -gpu=mem:managed
			# CXXFLAGS += -gpu=managed
		endif
	endif

	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif

	ifneq ($(nprof), 0)
		CXXFLAGS += -I$(NVHPC_DIR)/cuda/include -DNPROF
		ifeq ($(nprof), 2)
			CXXFLAGS += -DNPROF_DETAIL
		endif
		LDFLAGS += -g
	endif
else ifneq (, $(findstring pgc++, $(CXX)))
	CXX_SUPPORTED = yes
	CXXFLAGS = -march=core2
	LDFLAGS = 
	TETGENFLAGS = 

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2 -silent
	else ifeq ($(opt), 3)
		CXXFLAGS += -O3 -fast -silent
	endif
 
	ifeq ($(openmp), 1)
		CXXFLAGS += -mp
		LDFLAGS += -mp
	endif

	ifneq ($(nprof), 0)
			CXXFLAGS += -Minfo=mp -I$(NVHPC_DIR)/cuda/include -DNPROF
			ifeq ($(nprof), 2)
				CXXFLAGS += -DNPROF_DETAIL
			endif
			LDFLAGS += -L$(NVHPC_DIR)/cuda/lib64 -Wl,-rpath,$(NVHPC_DIR)/cuda/lib64 -lnvToolsExt
	endif
else
## Unsupported compiler: leave CXX_SUPPORTED unset. `make config` says so, and
## check-deps refuses the build -- which used to be done with a second `all:`
## target here, whose recipe was silently overridden by the real `all:` further
## down (make warned "overriding commands for target 'all'" and used the other).
endif

ifeq ($(hdf5), 1)
	CXXFLAGS += -DHDF5 -I$(HDF5_INCLUDE_DIR)
	LDFLAGS += -L$(HDF5_LIB_DIR) -lhdf5
	ifneq ($(OSNAME), Darwin)
		LDFLAGS += -Wl,-rpath,$(HDF5_LIB_DIR)
	endif
endif

## Is git in the path?
HAS_GIT := $(shell git --version 2>/dev/null)
ifneq ($(HAS_GIT),)
        ## Is this a git repository?
        IS_REPO := $(shell git rev-parse --is-inside-work-tree 2>/dev/null)
endif

SRCS =	\
	barycentric-fn.cxx \
	ats_output_scheduler.cxx \
	brc-interpolation.cxx \
	bc.cxx \
	binaryio.cxx \
	dynearthsol.cxx \
	earthquake_state.cxx \
	fields.cxx \
	geometry.cxx \
	ic.cxx \
	ic-read-temp.cxx \
	input.cxx \
	matprops.cxx \
	mesh.cxx \
	monitor.cxx \
	nn-interpolation.cxx \
	output.cxx \
	phasechanges.cxx \
	remeshing.cxx \
	rheology.cxx \
	runtime_info.cxx \
	markerset.cxx \
	knn.cxx

ifeq ($(use_gospl), 1)
	SRCS += gospl_driver/gospl-driver.cxx
endif

INCS =	\
	array2d.hpp \
	ats_output_scheduler.hpp \
	barycentric-fn.hpp \
	binaryio.hpp \
	constants.hpp \
	earthquake_state.hpp \
	monitor.hpp \
	parameters.hpp \
	matprops.hpp \
	sortindex.hpp \
	utils.hpp \
	mesh.hpp \
	markerset.hpp \
	output.hpp \
	knn.hpp

OBJS = $(SRCS:.cxx=.$(ndims)d$(suffix).o)

EXE = dynearthsol$(ndims)d$(suffix)


## Libraries

TET_SRCS = tetgen/predicates.cxx tetgen/tetgen.cxx
TET_INCS = tetgen/tetgen.h
TET_OBJS = $(TET_SRCS:.cxx=$(suffix).o)

TRI_SRCS = triangle/triangle.c
TRI_INCS = triangle/triangle.h
TRI_OBJS = $(TRI_SRCS:.c=$(suffix).o)

M_SRCS = $(TRI_SRCS)
M_INCS = $(TRI_INCS)
M_OBJS = $(TRI_OBJS)

ifeq ($(ndims), 3)
	M_SRCS += $(TET_SRCS)
	M_INCS += $(TET_INCS)
	M_OBJS += $(TET_OBJS)
	CXXFLAGS += -DTHREED
endif

ifeq ($(adaptive_time_step), 1)
	CXXFLAGS += -DATS
ifeq ($(use_R_S), 1)
	CXXFLAGS += -DRS
endif
endif

ifeq ($(useexo), 1)
	CXXFLAGS += $(EXO_CXXFLAGS)
	LDFLAGS += $(EXO_LDFLAGS)
endif

ifeq ($(usemmg), 1)
	CXXFLAGS += $(MMG_CXXFLAGS)
	LDFLAGS += $(MMG_LDFLAGS)
endif

ifeq ($(use_gospl), 1)
	CXXFLAGS += $(GOSPL_CXXFLAGS)
	LDFLAGS += $(GOSPL_LDFLAGS)
endif

C3X3_DIR = 3x3-C
C3X3_LIBNAME = 3x3$(suffix)

ANN_DIR = nanoflann
CXXFLAGS += -I$(ANN_DIR)/include

GOSPL_DIR = gospl_driver
CXXFLAGS += -I$(GOSPL_DIR)

KNN_BVH_DIR = knn-bvh
ifeq ($(openacc), 1)
	CXXFLAGS += -I$(KNN_BVH_DIR)/include
	KNN_BVH_LIB = $(KNN_BVH_DIR)/lib/libknn_bvh.$(ndims)d.a
	LDFLAGS += $(KNN_BVH_LIB)
endif

# Enable Array2D structure of Array
CXXFLAGS += -DSOA

## Action

.PHONY: all clean take-snapshot prepare build check-deps config

all: prepare
	$(MAKE) build

## Name the missing dependency and the command that installs it, rather than
## letting it surface as a linker error hundreds of lines down. A recipe, not
## $(error), so `make clean` and `make config` still work on a machine that
## cannot build yet.
##
## Fatal on macOS only: there this detection is the sole source of -I/-L and
## there is no system Boost nor Apple omp.h, so a miss is certain. Elsewhere the
## toolchain has paths no filesystem probe can see (CPATH, LIBRARY_PATH, a
## module sysroot), so a failed guess must not block a working build.
ifeq ($(OSNAME), Darwin)
DEP_FATAL = exit 1
BOOST_INSTALL_HINT = brew install boost
else
DEP_FATAL = echo "         (continuing: your compiler may find it via CPATH or a module)"
BOOST_INSTALL_HINT = apt install libboost-program-options-dev / dnf install boost-devel
endif

check-deps:
	@test -n "$(CXX_SUPPORTED)" || { \
		echo "Unsupported compiler '$(CXX)': no flags would be set. Use clang++,"; \
		echo "g++, icpc, nvc++ or pgc++, e.g. make CXX=g++."; exit 1; }
	@test -f "$(BOOST_ROOT_DIR)/include/boost/program_options.hpp" \
	  || test -f "$(BOOST_ROOT_DIR)/boost/program_options.hpp" \
	  || test -f /usr/include/boost/program_options.hpp || { \
		echo "Missing: Boost.ProgramOptions.   Install: $(BOOST_INSTALL_HINT)"; \
		echo "         Searched: $(strip $(BOOST_PREFIXES)) /usr"; \
		echo "         Or name one: make BOOST_ROOT_DIR=/prefix   (now $(BOOST_ROOT_DIR))"; \
		$(DEP_FATAL); }
ifeq ($(OSNAME), Darwin)
ifeq ($(openmp), 1)
ifneq (, $(findstring clang++, $(CXX)))
	@test -f "$(OPENMP_INCLUDE_DIR)/omp.h" \
	  && test -f "$(OPENMP_LIB_DIR)/libomp.dylib" || { \
		echo "Missing: an OpenMP runtime -- Apple clang has the pragmas, ships none."; \
		echo "         Install: brew install libomp        (or https://brew.sh)"; \
		echo "         Searched: $(strip $(OPENMP_LIB_DIRS))"; \
		echo "         Or: make openmp=0        make OPENMP_ROOT_DIR=/prefix"; \
		echo "             make BREW_PREFIX=/your/prefix   (now $(BREW_PREFIX))"; \
		exit 1; }
endif
endif
endif

## Print what the build resolved to -- the first thing to ask for when a build
## misbehaves elsewhere. Each library's architecture is shown because a path
## alone does not tell you whether it matches the host.
config:
	@echo "target             $(EXE)   (ndims=$(ndims), opt=$(opt))"
	@echo "platform           $(OSNAME) $(shell uname -m)"
	@echo "compiler           $(CXX)$(if $(CXX_SUPPORTED),,   ** UNSUPPORTED -- no flags will be set **)"
	@echo "openmp             $(strip $(openmp))"
ifeq ($(OSNAME), Darwin)
	@echo "brew prefix        $(BREW_PREFIX)"
ifeq ($(openmp), 1)
	@echo "  omp.h from       $(if $(strip $(OPENMP_INCLUDE_DIR)),$(OPENMP_INCLUDE_DIR),** NOT FOUND **)"
	@echo "  libomp from      $(if $(strip $(OPENMP_LIB_DIR)),$(OPENMP_LIB_DIR) [`lipo -archs $(OPENMP_LIB_DIR)/libomp.dylib 2>/dev/null || echo missing`],** NOT FOUND **)"
endif
	@echo "boost prefix       $(BOOST_ROOT_DIR)$(if $(BOOST_NOT_FOUND),   ** NOT FOUND -- fallback **)"
else
	@echo "boost prefix       $(BOOST_ROOT_DIR)"
endif
	@echo "  boost lib dir    $(BOOST_LIB_DIR)$(if $(filter Darwin,$(OSNAME)), [`lipo -archs $(BOOST_LIB_DIR)/libboost_program_options.dylib 2>/dev/null || echo 'no .dylib'`])"
	@echo "  boost -I         $(if $(strip $(BOOST_CXXFLAGS)),$(BOOST_CXXFLAGS),none -- relying on the compiler default search path)"
	@echo "hdf5               $(strip $(hdf5))"
ifeq ($(hdf5), 1)
	@echo "  hdf5 include     $(HDF5_INCLUDE_DIR)"
	@echo "  hdf5 lib dir     $(HDF5_LIB_DIR)"
endif
	@echo "usemmg             $(strip $(usemmg))"
	@echo "useexo             $(strip $(useexo))"
	@echo "use_gospl          $(strip $(use_gospl))"
	@echo "openacc            $(strip $(openacc))"
	@echo
	@echo "CXXFLAGS           $(CXXFLAGS) $(BOOST_CXXFLAGS)"
	@echo "LDFLAGS            $(LDFLAGS) $(BOOST_LDFLAGS)"

ifeq ($(use_gospl), 1)
.PHONY: install-gospl-wrapper
install-gospl-wrapper: 
	@echo '#!/bin/bash' > dynearthsol-gospl
	@echo '# Auto-generated wrapper for DynEarthSol with GoSPL support' >> dynearthsol-gospl
	@echo 'export PYTHONPATH="$(GOSPL_PYTHONPATH):$$PYTHONPATH"' >> dynearthsol-gospl
	@echo 'exec "$$(dirname "$$0")/$(EXE)" "$$@"' >> dynearthsol-gospl
	@chmod +x dynearthsol-gospl
endif

## check-deps first, as a prerequisite rather than another line in this recipe,
## so that the ordering holds under `make -j` and `make check-deps` still works
## on its own.
prepare: check-deps
	@if command -v git >/dev/null 2>&1 && git rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
		if git submodule status $(ANN_DIR) | grep -q '^[-+]'; then \
			echo "   Status mismatch. Updating submodule $(ANN_DIR)..."; \
			git submodule update --init --recursive $(ANN_DIR); \
		fi; \
	elif [ -f "$(ANN_DIR)/include/nanoflann.hpp" ]; then \
		:; \
	else \
		echo "Error: DES build requires $(ANN_DIR), but git is unavailable or this source tree is not a git checkout."; \
		echo "       Please initialize/provide $(ANN_DIR) before building with 'git submodule update --init --recursive $(ANN_DIR)'."; \
		exit 1; \
	fi
ifeq ($(openacc), 1)
	@if command -v git >/dev/null 2>&1 && git rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
		if git submodule status $(KNN_BVH_DIR) | grep -q '^[-+]'; then \
			echo "   Status mismatch. Updating submodule $(KNN_BVH_DIR)..."; \
			git submodule update --init --recursive $(KNN_BVH_DIR); \
		fi; \
	elif [ -f "$(KNN_BVH_DIR)/Makefile" ]; then \
		:; \
	else \
		echo "Error: OpenACC build requires $(KNN_BVH_DIR), but git is unavailable or this source tree is not a git checkout."; \
		echo "       Please initialize/provide $(KNN_BVH_DIR) before building with openacc=1."; \
		exit 1; \
	fi
endif
ifeq ($(usemmg), 1)
	@if command -v git >/dev/null 2>&1 && git rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
		if git submodule status mmg | grep -q '^[-+]'; then \
			echo "   Status mismatch. Updating submodule mmg..."; \
			git submodule update --init --recursive mmg; \
		fi; \
	elif [ -f "mmg/CMakeLists.txt" ]; then \
		:; \
	else \
		echo "Error: usemmg requires the mmg submodule, but git is unavailable or this source tree is not a git checkout."; \
		echo "       Please initialize/provide the mmg submodule before building with usemmg=1."; \
		exit 1; \
	fi

	@mkdir -p mmg/build
	@if [ ! -f "mmg/build/Makefile" ]; then \
		echo "   Configuring MMG..."; \
		cd mmg/build && LDFLAGS="" CXXFLAGS="" cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..; \
	fi
	@if [ ! -f "$(MMG_LIB)" ]; then \
		$(MAKE) -C mmg/build; \
	fi
endif

build: $(EXE) tetgen/tetgen triangle/triangle take-snapshot

$(EXE): $(M_OBJS) $(OBJS) $(C3X3_DIR)/lib$(C3X3_LIBNAME).a $(KNN_BVH_LIB) $(MMG_LIB)
		$(CXX) $(M_OBJS) $(OBJS) $(LDFLAGS) $(BOOST_LDFLAGS) \
			-L$(C3X3_DIR) -l$(C3X3_LIBNAME) \
			-o $@

ifeq ($(use_gospl), 1)
	@+$(MAKE) install-gospl-wrapper
	@echo ""
	@echo "=============================================="
	@echo "✅ DynEarthSol built with GoSPL support!"
	@echo "=============================================="
	@echo "🚀 To run with GoSPL support:"
	@echo ""
	@echo "Use the wrapper script (PYTHONPATH is set automatically):"
	@echo "  ./dynearthsol-gospl your_input.cfg"
	@echo ""
	@echo "Or set PYTHONPATH manually and use the regular executable:"
	@echo "  PYTHONPATH=$(GOSPL_PYTHONPATH):\$$PYTHONPATH ./$(EXE) your_input.cfg"
	@echo "=============================================="
	@echo ""
endif
ifeq ($(OSNAME), Darwin)  # fix for dynamic library problem on Mac
		@# A Boost built by hand with b2 records a bare install name, which the
		@# loader cannot resolve; rewrite it to the full path. A no-op for a
		@# packaged Boost, whose install name is already absolute.
		install_name_tool -change libboost_program_options.dylib $(BOOST_LIB_DIR)/libboost_program_options.dylib $@
		@# Record each dependency's directory as an rpath, skipping any that is
		@# already there: install_name_tool fails outright on a duplicate, and
		@# these two are the same directory whenever Boost and libomp come from
		@# one prefix (an activated conda env, say). An empty OPENMP_RPATH --
		@# openmp=0, or a runtime that was not found -- drops out of the loop.
		@for d in $(BOOST_LIB_DIR) $(OPENMP_RPATH); do \
			otool -l $@ | grep -qF " path $$d (" || install_name_tool -add_rpath $$d $@; \
		done
ifeq ($(useexo), 1)  # fix for dynamic library problem on Mac
		install_name_tool -change libexodus.dylib $(EXO_LIB_DIR)/libexodus.dylib $@
endif
ifeq ($(usemmg), 1)  # fix for dynamic library problem on Mac
ifeq ($(ndims), 3)
		install_name_tool -change libmmg3d.dylib $(MMG_LIB_DIR)/libmmg3d.dylib $@
else
		install_name_tool -change libmmg2d.dylib $(MMG_LIB_DIR)/libmmg2d.dylib $@
endif
endif # end of usemmg
endif # end of Darwin

take-snapshot:
	@# snapshot of the code for building the executable
	@echo Flags used to compile the code: > snapshot.diff
	@echo '  '  CXX=$(CXX) opt=$(opt) openmp=$(openmp) >> snapshot.diff
	@echo '  '  CXXFLAGS=$(CXXFLAGS) >> snapshot.diff
	@echo '  '  LDFLAGS=$(LDFLAGS) >> snapshot.diff
	@echo '  '  PATH="$(PATH)" >> snapshot.diff
	@echo '  '  LD_LIBRARY_PATH="$(LD_LIBRARY_PATH)" >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
ifneq ($(HAS_GIT),)
ifneq ($(IS_REPO),)
	@echo '==== Summary of the code ====' >> snapshot.diff
	@git show -s >> snapshot.diff
	@echo >> snapshot.diff
	@echo >> snapshot.diff
	@git status >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (not checked-in) ==' >> snapshot.diff
	@echo >> snapshot.diff
	@git diff >> snapshot.diff
	@echo >> snapshot.diff
	@echo '== Code modification (checked-in but not in "origin") ==' >> snapshot.diff
	@echo >> snapshot.diff
	@git log --patch -- origin..HEAD >> snapshot.diff
else
	@echo "Warning: Not a git repository. Cannot take code snapshot." | tee -a snapshot.diff
	@echo "Warning: Use 'git clone' to copy the code!" | tee -a snapshot.diff
endif
else
	@echo "'git' is not in path, cannot take code snapshot." >> snapshot.diff
endif

$(OBJS): %.$(ndims)d$(suffix).o : %.cxx $(INCS)
	$(CXX) $(CXXFLAGS) $(BOOST_CXXFLAGS) -c $< -o $@

$(KNN_BVH_LIB):
	$(MAKE) -C $(KNN_BVH_DIR) NDIM=$(ndims)

$(TRI_OBJS): %$(suffix).o : %.c $(TRI_INCS)
	@# Triangle cannot be compiled with -O2
	$(CXX) $(CXXFLAGS) -O1 -DTRILIBRARY -DREDUCED -DANSI_DECLARATORS -c $< -o $@

triangle/triangle: triangle/triangle.c
	$(CXX) $(CXXFLAGS) -O1 -DREDUCED -DANSI_DECLARATORS triangle/triangle.c -o $@

tetgen/predicates$(suffix).o: tetgen/predicates.cxx $(TET_INCS)
	@# Compiling J. Shewchuk predicates, should always be
	@# equal to -O0 (no optimization). Otherwise, TetGen may not
	@# work properly.
	$(CXX) $(CXXFLAGS) -DTETLIBRARY -O0 -c $< -o $@

tetgen/tetgen$(suffix).o: tetgen/tetgen.cxx $(TET_INCS)
	$(CXX) $(CXXFLAGS) -DNDEBUG -DTETLIBRARY $(TETGENFLAG) -c $< -o $@

tetgen/tetgen: tetgen/predicates.cxx tetgen/tetgen.cxx
	$(CXX) $(CXXFLAGS) -O0 -DNDEBUG $(TETGENFLAG) tetgen/predicates.cxx tetgen/tetgen.cxx -o $@

$(C3X3_DIR)/lib$(C3X3_LIBNAME).a:
	@# CUDA_DIR only when there is an NVHPC_DIR to build it from; otherwise this
	@# handed the sub-make a literal "/cuda" on every single build.
	@+$(MAKE) -C $(C3X3_DIR) openacc=$(openacc) nofma=$(nofma) \
		$(if $(strip $(NVHPC_DIR)),CUDA_DIR=$(NVHPC_DIR)/cuda)

clean-submodules:
	@echo "   Cleaning external submodules..."
	@if [ -d "mmg/build" ]; then \
		echo "   Removing MMG build directory..."; \
		rm -rf mmg/build; \
	fi
	@if [ -f "$(KNN_BVH_DIR)/Makefile" ]; then \
		$(MAKE) -C $(KNN_BVH_DIR) clean NDIM=$(ndims) >/dev/null 2>&1; \
	fi

deepclean: 
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
ifeq ($(use_gospl), 1)
	@rm -f dynearthsol-gospl
endif
	@+$(MAKE) -C $(C3X3_DIR) clean openacc=$(openacc)
	
cleanall: clean
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
	@+$(MAKE) -C $(C3X3_DIR) clean openacc=$(openacc)
ifeq ($(use_gospl), 1)
	@rm -f dynearthsol-gospl
endif
	@+$(MAKE) -C $(C3X3_DIR) clean

clean:
	@rm -f $(OBJS) $(EXE)
