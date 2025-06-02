# C compiler
# ----------
CC            = g++
CFLAGS        = -Wall -O3 -ffast-math #-funsafe-loop-optimizations

ifeq ($(nprof), 1)
	CC = nvc++
endif
ifeq ($(openacc), 1)
	CC = nvc++
endif

ifneq (, $(findstring nvc++, $(CC)))
	CFLAGS = Wall -O3 -mno-fma

	ifeq ($(nprof), 1)
		CFLAGS += -Minfo=mp,accel -I$(CUDA_DIR)/include
	endif
	ifeq ($(openacc), 1)
		CFLAGS += -acc=gpu -gpu=mem:managed,nofma -cuda
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
BIN    = lib3x3.a
OBJ    = dsytrd3.o dsyevc3.o dsyevq3.o dsyevh3.o
HEADER = dsytrd3.h dsyevc3.h dsyevq3.h dsyevh3.h


# Rules
# -----
$(BIN): Makefile $(OBJ)
	$(AR) $(ARFLAGS) $(BIN) $(OBJ)

%.o : %.c Makefile $(HEADER)
	$(CC) $(CFLAGS) -c $(INCFLAGS) $<

.PHONY: clean
clean:
	@rm -f $(BIN) $(OBJ)


