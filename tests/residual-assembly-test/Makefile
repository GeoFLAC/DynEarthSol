# Directly link the real assembly/BC/material routines. Unused solver sections
# are discarded; no production routine is replaced by a test stub.
CXX = g++
ndims ?= 2
mode ?= B0
backend ?= serial
ROOT := ../..
TARGET := residual-$(ndims)d-$(backend)-$(mode)
CPPFLAGS += -I$(ROOT) -I$(ROOT)/nanoflann/include -DGPP1X -DSOA -D$(mode)
ifeq ($(ndims),3)
CPPFLAGS += -DTHREED
endif
CXXFLAGS += -std=c++17 -O0 -g -ffunction-sections -fdata-sections
LDFLAGS += -Wl,--gc-sections
ifeq ($(backend),omp)
CXXFLAGS += -fopenmp
endif
ifeq ($(backend),asan)
CXXFLAGS += -fsanitize=address,undefined
endif
ifeq ($(backend),acc)
CXX = nvc++
GPU_CC ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.')
CPPFLAGS += -DACC -I$(ROOT)/knn-bvh/include
CXXFLAGS += -acc=gpu -gpu=cc$(GPU_CC),mem:managed -cuda
endif
SOURCES := main.cxx $(addprefix $(ROOT)/,fields.cxx bc.cxx matprops.cxx utils.cxx rheology.cxx pseudo_transient.cxx)
HEADERS := $(wildcard $(ROOT)/*.hpp)
.PHONY: check clean
ifeq ($(backend),acc)
# NVHPC retains device registrations even for otherwise unused functions.
# Reuse the production dependency graph instead of substituting test stubs.
check:
	env -u MAKEFLAGS -u MFLAGS mode=$(mode) $(MAKE) -C $(ROOT) -f Makefile -f tests/residual-assembly-test/acc.mk ndims=$(ndims) openacc=1 hdf5=0 GPU_CC=$(GPU_CC) residual-assembly-check
else
.DEFAULT_GOAL := $(TARGET)
$(ROOT)/3x3-C/lib3x3.a: $(wildcard $(ROOT)/3x3-C/*.c $(ROOT)/3x3-C/*.h)
	$(MAKE) -C $(ROOT)/3x3-C openacc=0
$(TARGET): $(SOURCES) $(HEADERS) $(ROOT)/3x3-C/lib3x3.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SOURCES) $(LDFLAGS) $(ROOT)/3x3-C/lib3x3.a -o $@
check: $(TARGET)
	./$(TARGET)
endif
clean:
	rm -f residual-2d-* residual-3d-* driver.*.o des-main.*.o
