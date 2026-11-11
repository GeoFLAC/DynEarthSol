CXX = g++
ndims ?= 2
ROOT := ../..
TARGET := moment-$(ndims)d
CPPFLAGS += -I$(ROOT) -I$(ROOT)/nanoflann/include -DGPP1X -DSOA
CXXFLAGS += -std=c++17 -O0 -g -fopenmp -ffunction-sections -fdata-sections -Wno-unknown-pragmas
LDFLAGS += -fopenmp -Wl,--gc-sections
ifeq ($(ndims),3)
CPPFLAGS += -DTHREED
endif
.PHONY: check clean
$(TARGET): main.cxx $(ROOT)/earthquake_state.cxx $(ROOT)/matprops.cxx $(ROOT)/utils.cxx $(wildcard $(ROOT)/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.cxx $(ROOT)/matprops.cxx $(ROOT)/utils.cxx $(LDFLAGS) -o $@
check: $(TARGET)
	@task_dir=$$(mktemp -d) || exit 1; \
	trap 'rm -rf "$$task_dir"' EXIT; \
	cd "$$task_dir" && OMP_NUM_THREADS=1 $(CURDIR)/$(TARGET)
clean:
	rm -f moment-2d moment-3d
