
CASE := test-3d-equ-tiny.cfg
#CASE := test-tiny-long.cfg
#CASE := test-big.cfg
#CASE := test-topo.cfg
FRAME := 4
OMP := 1
ACC := 0
TFLAG = # time -v
DEBUG := 0
NSYS := 0
VTK := 0
NDIMS := 3
EXE := ../dynearthsol$(NDIMS)d
FUNC := compute_mass
NUM := 0
DATE := $(shell date +%Y-%m-%d)
HOST := $(shell hostname -s)
KERNEL := CPU
ifeq ($(ACC), 1)
	KERNEL := GPU
	EXE := ../dynearthsol$(NDIMS)d.gpu
endif

MODELNAME := $(shell awk -F '=' '/^modelname/ {gsub(/^[ \t]+|[ \t]+/, "", $$2); print $$2}' $(CASE))
JOBNAME = ${CASE}.${OMP}-${DATE}-${HOST}-${KERNEL}-${NUM}

GITID = $(shell git rev-parse HEAD | cut -c1-7)

ifeq ($(DEBUG), 1)
	ifeq ($(ACC), 1)
		DEBUGFLAG = CUDA_LAUNCH_BLOCKING=1 CUDA_MEMCHECK=1 # NV_ACC_NOTIFY=3 NV_ACC_DEBUG=0x800  # compute-sanitizer #  # 
	else
		DEBUGFLAG = valgrind --tool=helgrind --log-file=valgrind-out.txt -v
	endif
endif
ifeq ($(NSYS), 1)
	DEBUGFLAG = nsys profile --trace=cuda,nvtx,osrt,openmp,openacc --output=${MODELNAME}.${GITID}.${CASE}.${OMP}
	# NSYSSUM = nsys stats -o ${MODELNAME}.${GITID}.${CASE}.${OMP}.csv ${MODELNAME}.${GITID}.${CASE}.${OMP}.nsys-rep
endif
ifeq ($(VTK), 1)
	VTKFLAG = ${MAKE} vtk
endif


# Use single cpu for omp to avoid round-off errors
RUN_TEST = rm -f ${MODELNAME}.*; OMP_NUM_THREADS=${OMP} ${TFLAG} ${DEBUGFLAG} ${EXE} ${CASE}
RESTART_TEST = OMP_NUM_THREADS=${OMP} ${TFLAG} ${DEBUGFLAG} ${EXE} ${CASE}

ORIG = orig-${CASE}
CURR = ~/data/jobs/${CASE}

all: cmp

test:
	${RUN_TEST}

job:
	# run jobs on server
	${RUN_TEST} > screen-${JOBNAME}.log &
	cp ${CASE} ${JOBNAME}.cfg

auto:
	${MAKE} CASE=test-big.cfg OMP=4
	${MAKE} CASE=test-topo.cfg OMP=4
	${MAKE} CASE=test-model.cfg OMP=4

set:
	# run test, store test output to old/
	${RUN_TEST}
	${MAKE} store

detect-orig:
	$(if $(wildcard ${ORIG}),,$(error Error: Original results directory '${ORIG}' doesn't exist. Please run 'make set' first.))

cmp: detect-orig
	${RUN_TEST}
	${NSYSSUM} > /dev/null
	${VTKFLAG}
	python ./compare.py ${ORIG} ${FRAME} ${MODELNAME}
	cp ${MODELNAME}.info ${MODELNAME}.${GITID}.${CASE}.${OMP}.info

restart: detect-orig
	${RESTART_TEST}
	${NSYSSUM} > /dev/null
	${VTKFLAG}
	python ./compare.py ${ORIG} ${FRAME} ${MODELNAME}
	cp ${MODELNAME}.info ${MODELNAME}.${GITID}.${CASE}.${OMP}.info

# restart reproducibility test; config file with restart information is needed
IS_RESTART_EXISTS := $(shell grep -i "^[[:space:]\#]*is_restarting[[:space:]]*=" $(CASE))
RESTART_FRAME := $(shell awk -F '=' '/^[[:space:]]*restarting_from_frame/ {gsub(/^[ \t]+|[ \t]+/, "", $$2); print $$2}' $(CASE))
MODELNAME_EXISTS := $(shell grep -i "^[[:space:]\#]*restarting_from_modelname[[:space:]]*=" $(CASE))
FRC_CASE = $(CASE:%.cfg=%-frc.cfg)
fresh-restart-cmp:
	$(if $(strip $(IS_RESTART_EXISTS)),,$(error Error: 'sim.is_restarting' does not exist in ${CASE}...))
	$(if $(strip $(RESTART_FRAME)),,$(error Error: 'sim.restarting_from_frame' is not set in ${CASE}...))
	$(if $(strip $(MODELNAME_EXISTS)),,$(error Error: 'sim.restarting_from_modelname' does not exist in ${CASE}...))
	@sed 's|^[[:space:]\#]*is_restarting[[:space:]]*=.*|is_restarting = no|' ${CASE} > ${FRC_CASE}
	${MAKE} set CASE=${FRC_CASE}
	@sed 's|^[[:space:]\#]*is_restarting[[:space:]]*=.*|is_restarting = yes|' ${CASE} > ${FRC_CASE}
	@sed 's|^[[:space:]\#]*restarting_from_modelname[[:space:]]*=.*|restarting_from_modelname = orig-${FRC_CASE}/${MODELNAME}|' ${FRC_CASE} > ${FRC_CASE}.tmp
	@mv ${FRC_CASE}.tmp ${FRC_CASE}
	${MAKE} CASE=${FRC_CASE}
	rm -f ${FRC_CASE}

store:
	rm -rf ${ORIG}
	mkdir -p ${ORIG}
	mv ${MODELNAME}.* ${ORIG}
	cp ../snapshot.diff ${ORIG}

setjob:
	# st job to compare
	rm -rf orig-job
	ln -s ~/data/jobs/${CASE} orig-job

cmpjob:
	python ./compare.py orig-job ${CURR} ${FRAME} ${MODELNAME}

vtk:
	python ../2vtk.py -m  ${MODELNAME}

.PHONY: all set cmp store test job fresh-restart-cmp
