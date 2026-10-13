# Include after the repository Makefile, running from the repository root.
TEST_DIR := tests/residual-assembly-test
TEST_EXE := $(TEST_DIR)/residual-$(ndims)d-acc-$(mode)
TEST_DRIVER := $(TEST_DIR)/driver.$(ndims)d$(suffix).$(mode).o
TEST_MAIN := $(TEST_DIR)/des-main.$(ndims)d$(suffix).o
TEST_OBJS := $(filter-out dynearthsol.$(ndims)d$(suffix).o,$(OBJS))

$(TEST_DRIVER): $(TEST_DIR)/main.cxx $(INCS) $(BUILD_STAMP)
	$(CXX) $(CXXFLAGS) $(BOOST_CXXFLAGS) -I. -D$(mode) -c $< -o $@
$(TEST_MAIN): dynearthsol.cxx $(INCS) $(BUILD_STAMP)
	$(CXX) $(CXXFLAGS) $(BOOST_CXXFLAGS) -Dmain=des_application_main -c $< -o $@
$(TEST_EXE): $(TEST_DRIVER) $(TEST_MAIN) $(TEST_OBJS) $(M_OBJS) $(C3X3_DIR)/lib$(C3X3_LIBNAME).a $(KNN_BVH_LIB) $(MMG_LIB) $(LINK_STAMP)
	$(CXX) $(TEST_DRIVER) $(TEST_MAIN) $(TEST_OBJS) $(M_OBJS) $(LDFLAGS) $(BOOST_LDFLAGS) -L$(C3X3_DIR) -l$(C3X3_LIBNAME) -o $@
.PHONY: residual-assembly-check
residual-assembly-check: $(TEST_EXE)
	./$(TEST_EXE)
