# Backend check for the header-only unsaturated.hpp. No production object is
# linked; the header is the code under test.
#   make check backend=serial|omp|asan|acc
CXX = g++
backend ?= serial
ROOT := ../..
TARGET := unsaturated-$(backend)
CPPFLAGS += -I$(ROOT)
CXXFLAGS += -std=c++17 -O2 -Wall -Wextra -Wno-unknown-pragmas -Werror
ifeq ($(backend),omp)
CXXFLAGS += -fopenmp
endif
ifeq ($(backend),asan)
CXXFLAGS += -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all
endif
ifeq ($(backend),acc)
CXX = nvc++
GPU_CC ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.')
ifeq ($(strip $(GPU_CC)),)
GPU_CC := 80
endif
nofma ?= 0
CPPFLAGS += -DACC
# The production nvc++ accelerator flags (root Makefile, openacc=1). -Minfo
# reports which loops and routines were offloaded, so a silent host fallback
# shows in the build log.
CXXFLAGS = -std=c++17 -O2 -Minfo=accel --diag_suppress=declared_but_not_referenced -acc=gpu -cuda
ifeq ($(nofma),1)
CXXFLAGS += -mno-fma -gpu=cc$(GPU_CC),mem:managed,nofma
else
CXXFLAGS += -gpu=cc$(GPU_CC),mem:managed
endif
endif

.PHONY: check clean
.DEFAULT_GOAL := check
$(TARGET): main.cxx $(ROOT)/unsaturated.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.cxx -o $@
check: $(TARGET)
	./$(TARGET)
clean:
	rm -f unsaturated-serial unsaturated-omp unsaturated-asan unsaturated-acc
