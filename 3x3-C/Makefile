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
		## GPU_CC = compute capability, from the parent make. Blank lets nvc++ infer it
		## from the build host, which is what makes a tree shared between hosts link
		## objects built for the wrong arch. Assembled into ONE -gpu=: repeating the
		## option leaves it to nvc++ whether occurrences merge or the last one wins.
		GPU_OPTS = mem:managed
		ifneq ($(strip $(GPU_CC)),)
			GPU_OPTS := cc$(strip $(GPU_CC)),$(GPU_OPTS)
		endif
		ifeq ($(nofma), 1)
			GPU_OPTS := $(GPU_OPTS),nofma
		endif
		CFLAGS += -acc=gpu -gpu=$(GPU_OPTS) -cuda
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
# The stamp DELETES the objects when flags change: make 3.81 compares whole
# seconds, so mtime is unreliable here, and the delete must run in an earlier
# process than the build -- hence the re-entry in `all`.
STAMP = .build-flags-3x3$(suffix)

.PHONY: all clean FORCE
all: $(STAMP)
	@$(MAKE) -q --no-print-directory $(BIN) || $(MAKE) --no-print-directory $(BIN)

$(STAMP): FORCE
	@echo '$(CC)|$(CFLAGS)|$(INCFLAGS)' > $@.tmp
	@if cmp -s $@.tmp $@; then rm -f $@.tmp; else \
		test -f $@ && echo "   3x3-C flags changed since the last build -- rebuilding"; \
		mv $@.tmp $@; rm -f $(OBJ) $(BIN); \
	fi

FORCE:

$(BIN): Makefile $(OBJ)
	$(AR) $(ARFLAGS) $(BIN) $(OBJ)

%$(suffix).o : %.c Makefile $(HEADER)
	$(CC) $(CFLAGS) -c $(INCFLAGS) $< -o $@

clean:
	@rm -f $(BIN) $(OBJ) $(STAMP)


