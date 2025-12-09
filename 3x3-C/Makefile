# C compiler
# ----------
CC            = g++
CFLAGS        = -Wall -O3 -ffast-math #-funsafe-loop-optimizations

ifeq ($(nprof), 1)
	CC = nvc++
endif
ifeq ($(openacc), 1)
	CC = nvc++
	suffix = .gpu
endif

ifneq (, $(findstring nvc++, $(CC)))
	CFLAGS = -Wall -O3
	ifeq ($(nofma), 1)
		CFLAGS += -mno-fma
	endif

	ifeq ($(nprof), 1)
		CFLAGS += -Minfo=mp,accel -I$(CUDA_DIR)/include
	endif
	ifeq ($(openacc), 1)
		ifeq ($(nofma), 1)
			CFLAGS += -acc=gpu -gpu=mem:managed,nofma -cuda
		else
			CFLAGS += -acc=gpu -gpu=mem:managed -cuda
		endif
	endif

endif

# Archiver
AR            = ar
ARFLAGS       = -crs

# Compiler and linker options
# ---------------------------
INCFLAGS      =
LOCAL_LDFLAGS =

# Project files
# -------------
BIN    = lib3x3$(suffix).a
OBJ    = dsytrd3$(suffix).o dsyevc3$(suffix).o dsyevq3$(suffix).o dsyevh3$(suffix).o
HEADER = dsytrd3.h dsyevc3.h dsyevq3.h dsyevh3.h


# Rules
# -----
$(BIN): Makefile $(OBJ)
	$(AR) $(ARFLAGS) $(BIN) $(OBJ)

%$(suffix).o : %.c Makefile $(HEADER)
	$(CC) $(CFLAGS) -c $(INCFLAGS) $< -o $@

.PHONY: clean
clean:
	@rm -f $(BIN) $(OBJ)


