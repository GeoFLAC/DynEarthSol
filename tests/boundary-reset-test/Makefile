CXX = g++
ROOT := ../..
CXXFLAGS += -std=c++17 -O0 -g -ffunction-sections -fdata-sections -fopenmp -Wno-unknown-pragmas
CPPFLAGS += -I$(ROOT) -I$(ROOT)/nanoflann/include -DGPP1X -DSOA
LDFLAGS += -Wl,--gc-sections
.PHONY: check clean
check: reset-test
	OMP_NUM_THREADS=1 ./reset-test
reset-test: main.cxx $(ROOT)/remeshing.cxx $(wildcard $(ROOT)/*.hpp)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.cxx $(LDFLAGS) -o $@
clean:
	rm -f reset-test
