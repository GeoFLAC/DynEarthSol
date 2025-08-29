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
##  - nprof = 1 : enable NVHPC nprof profiling build (uses nvc++ when set).
##  - gprof = 1 : enable GNU gprof instrumentation (-pg).
##  - usemmg = 1 : enable MMG mesh optimization support (requires MMG headers/libs).
##  - netcdf = 1 : enable NetCDF output support (requires netCDF and netcdf-cxx4).
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
netcdf = 0

ifeq ($(ndims), 2)
	useexo = 0    # for now, can import only 3d exo mesh
endif

ifneq ($(adaptive_time_step), 1)
	use_R_S = 0   # Rate - State friction law only works with adaptive time stepping technique
endif

## Select C++ compiler and set paths to necessary libraries
ifeq ($(openacc), 1)
	CXX = nvc++
	suffix = .gpu
else
	ifeq ($(nprof), 1)
		CXX = nvc++
	else
		CXX = g++
	endif
endif
CXX_BACKEND = ${CXX}

## path to netCDF's base directory, if not in standard system location
NETCDF_DIR = # /path/to/netcdf-c
NETCDFCXX_DIR = # /path/to/netcdf-cxx4

## path to cuda's base directory
NVHPC_DIR = # /cluster/nvidia/hpc_sdk/Linux_x86_64/21.2

## path to Boost's base directory, if not in standard system location
BOOST_ROOT_DIR =

########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

OSNAME := $(shell uname -s)

BOOST_LDFLAGS = -lboost_program_options
ifdef BOOST_ROOT_DIR
	# check existence of stage/ directory
	has_stage_dir = $(wildcard $(BOOST_ROOT_DIR)/stage)
	ifeq (, $(has_stage_dir))
		# no stage dir, BOOST_ROOT_DIR is the installation directory
		BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)/include
		BOOST_LIB_DIR = $(BOOST_ROOT_DIR)/lib
	else
		# with stage dir, BOOST_ROOT_DIR is the build directory
		BOOST_CXXFLAGS = -I$(BOOST_ROOT_DIR)
		BOOST_LIB_DIR = $(BOOST_ROOT_DIR)/stage/lib
	endif
	BOOST_LDFLAGS += -L$(BOOST_LIB_DIR)
	ifneq ($(OSNAME), Darwin)  # Apple's ld doesn't support -rpath
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
endif


ifneq (, $(findstring clang++, $(CXX)))
	CXXFLAGS = -v -DGPP1X
	LDFLAGS = -v

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	endif
 
	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif

else ifneq (, $(findstring g++, $(CXX_BACKEND))) # if using any version of g++
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
	CXXFLAGS = -g -Minfo=mp,accel
	LDFLAGS =
	TETGENFLAGS = 

	ifeq ($(opt), 1)
		CXXFLAGS += -O1
	else ifeq ($(opt), 2)
		CXXFLAGS += -O2
	endif

	ifeq ($(openacc), 1)
		CXXFLAGS += -acc=gpu -cuda -DACC
		LDFLAGS += -acc=gpu -gpu=mem:managed -cuda
		# CXXFLAGS += -acc=gpu -Mcuda -DACC
		# LDFLAGS += -acc=gpu -gpu=managed -Mcuda
		ifeq ($(nprof), 1)
			CXXFLAGS += -gpu=mem:managed,nofma -mno-fma
			# CXXFLAGS += -gpu=managed,nofma -mno-fma
		else
			CXXFLAGS += -gpu=mem:managed
			# CXXFLAGS += -gpu=managed
		endif
	endif

	ifeq ($(openmp), 1)
		CXXFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif

	ifeq ($(nprof), 1)
		CXXFLAGS += -I$(NVHPC_DIR)/cuda/include -DUSE_NPROF
		LDFLAGS += -L$(NVHPC_DIR)/cuda/lib64 -Wl,-rpath,$(NVHPC_DIR)/cuda/lib64 -lnvToolsExt -g
	endif
else ifneq (, $(findstring pgc++, $(CXX)))
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

	ifeq ($(nprof), 1)
			CXXFLAGS += -Minfo=mp -I$(NVHPC_DIR)/cuda/include -DUSE_NPROF
			LDFLAGS += -L$(NVHPC_DIR)/cuda/lib64 -Wl,-rpath,$(NVHPC_DIR)/cuda/lib64 -lnvToolsExt
	endif
else
# the only way to display the error message in Makefile ...
all:
	@echo "Unknown compiler, check the definition of 'CXX' in the Makefile."
	@false
endif

ifeq ($(netcdf), 1)
	CXXFLAGS += -DNETCDF  -I$(NETCDFCXX_DIR)/build/include
	LDFLAGS += -L$(NETCDF_DIR)/lib -lnetcdf -L$(NETCDFCXX_DIR)/build/lib64 -lnetcdf-cxx4 
endif

## Is git in the path?
HAS_GIT := $(shell git --version 2>/dev/null)
ifneq ($(HAS_GIT),)
        ## Is this a git repository?
        IS_REPO := $(shell git rev-parse --s-inside-work-tree 2>/dev/null)
endif

SRCS =	\
	barycentric-fn.cxx \
	brc-interpolation.cxx \
	bc.cxx \
	binaryio.cxx \
	dynearthsol.cxx \
	fields.cxx \
	geometry.cxx \
	ic.cxx \
	ic-read-temp.cxx \
	input.cxx \
	matprops.cxx \
	mesh.cxx \
	nn-interpolation.cxx \
	output.cxx \
	phasechanges.cxx \
	remeshing.cxx \
	rheology.cxx \
	markerset.cxx \
	knn.cxx

INCS =	\
	array2d.hpp \
	barycentric-fn.hpp \
	binaryio.hpp \
	constants.hpp \
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

C3X3_DIR = 3x3-C
C3X3_LIBNAME = 3x3$(suffix)

ANN_DIR = nanoflann
CXXFLAGS += -I$(ANN_DIR)/include

## Action

.PHONY: all clean take-snapshot

all: $(EXE) tetgen/tetgen triangle/triangle take-snapshot

$(EXE): $(M_OBJS) $(OBJS) $(C3X3_DIR)/lib$(C3X3_LIBNAME).a
		$(CXX) $(M_OBJS) $(OBJS) $(LDFLAGS) $(BOOST_LDFLAGS) \
			-L$(C3X3_DIR) -l$(C3X3_LIBNAME) \
			-o $@
ifeq ($(OSNAME), Darwin)  # fix for dynamic library problem on Mac
		install_name_tool -change libboost_program_options.dylib $(BOOST_LIB_DIR)/libboost_program_options.dylib $@
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
	@+$(MAKE) -C $(C3X3_DIR) openacc=$(openacc) nprof=$(nprof) CUDA_DIR=$(NVHPC_DIR)/cuda

deepclean: 
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
	@+$(MAKE) -C $(C3X3_DIR) clean openacc=$(openacc)
	
cleanall: clean
	@rm -f $(TET_OBJS) $(TRI_OBJS) $(OBJS) $(EXE)
	@+$(MAKE) -C $(C3X3_DIR) clean openacc=$(openacc)

clean:
	@rm -f $(OBJS) $(EXE)