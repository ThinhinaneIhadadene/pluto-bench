
# common Makefile
BASEDIR=$(dir $(lastword $(MAKEFILE_LIST)))

CC=gcc

MPI=openmpi
NPROCS=4
NTHREADS=4
NTHREADS_WITH_MPI=2
NTHREADS_LIST = 1 2 4 8 16 32
NPROCS_LIST = 1 2 4 8 16 32
TILE_SIZES = 16 32 64 128 256
#paths
POLYRTINCDIR=$(BASEDIR)../../../../../pluto/polyrt
POLYBENCHINCDIR=$(BASEDIR)polybench/utilities
POLYBENCHSRC=$(BASEDIR)polybench/utilities/polybench.c
HOSTS_FILE=$(BASEDIR)hosts
PLC=$(BASEDIR)../../../../../pluto/polycc

# Intel MKL and AMD ACML library paths
MKL=/opt/intel/mkl
ACML=/usr/local/acml

#Metis libray path
METIS=/usr/local/lib/libmetis.a

#Intel scalapack library paths and compiler flags
S_CFLAGS +=  -DMKL_ILP64 -I$(MKLROOT)/include
S_LFLAGS +=  $(MKLROOT)/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_cdft_core.a $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lpthread -lm


ifeq ($(CC), icc)
	CXX           := icpc
	OPT_FLAGS     := -O3 -xHost -ansi-alias -ipo -fp-model precise
	PAR_FLAGS     := -parallel
	OMP_FLAGS     := -openmp
	ifeq ($(MPI), mvapich)
		MPICC        := mpicc -cc=icc -D__MPI
		MPICXX       := mpicc -cc=icpc -D__MPI
	else ifeq ($(MPI), intelmpi)
		MPICC        := mpiicc -D__MPI
		MPICXX       := mpiicpc -D__MPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
	else #openmpi
		MPICC        := OMPI_CC=icc mpicc -D__MPI
		MPICXX       := OMPI_CC=icpc mpicc -D__MPI
	endif
else
	# for gcc
	CXX           := g++
	OPT_FLAGS     := -O3 -march=native -mtune=native -ftree-vectorize
	PAR_FLAGS     := -ftree-parallelize-loops=$(NTHREADS)
	OMP_FLAGS     := -fopenmp
	ifeq ($(MPI), mvapich)
		MPICC        := mpicc -cc=gcc -D__MPI
		MPICXX       := mpicc -cc=g++ -D__MPI
	else ifeq ($(MPI), intelmpi)
		MPICC        := mpicc -D__MPI
		MPICXX       := mpicxx -D__MPI
	else #openmpi
		MPICC        := OMPI_CC=gcc mpicc -D__MPI
		MPICXX       := OMPI_CC=g++ mpic++ -D__MPI
	endif
endif

ifeq ($(NON_CLUSTER), true)
	MPI = openmpi
endif

ifeq ($(MPI), mvapich)
	MPIARGS       = -np $(NPROCS) -hostfile $(HOSTS_FILE)
	MPIENV        = MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS_WITH_MPI)
	MPIRUN        = mpirun_rsh $(MPIARGS) $(MPIENV)
else ifeq ($(MPI), intelmpi)
	MPIARGS       = -ppn=1 -np $(NPROCS) -hostfile $(HOSTS_FILE)
	MPIENV        = -env I_MPI_FABRICS ofa -env OMP_NUM_THREADS $(NTHREADS_WITH_MPI)
	MPIRUN        = mpirun $(MPIARGS) $(MPIENV)
else #openmpi
	MPIARGS       = -np $(NPROCS)
	MPIENV        = OMP_NUM_THREADS=$(NTHREADS_WITH_MPI)
	MPIRUN        = $(MPIENV) mpirun $(MPIARGS)
endif

CFLAGS += -DTIME
LDFLAGS += -lm
PLCFLAGS += --isldep --lastwriter --indent
DISTOPT_FLAGS += --cloogsh --timereport
TILEFLAGS +=

ifdef PERFCTR
	CFLAGS += -DPERFCTR -L/usr/local/lib64 -lpapi
endif

ifdef POLYBENCH
	CFLAGS += -DPOLYBENCH_USE_SCALAR_LB -DPOLYBENCH_TIME -I $(POLYBENCHINCDIR) $(POLYBENCHSRC)
	DISTOPT_FLAGS += --variables_not_global
endif

all: orig par dist

$(SRC).par.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel $(TILEFLAGS) $(PLCFLAGS)  -o $@

#for distributed version with no communication optimisation
$(SRC).dist.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --nocommopt --tile $(TILEFLAGS) $(PLCFLAGS)  -o $@
#for distributed version with comm opt
$(SRC).distopt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@
#for distributed version with comm opt using foifi
$(SRC).distopt_foifi.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --commopt_foifi --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@
#for distributed version with comm opt using fop
$(SRC).distopt_fop.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --commopt_fop --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

orig: $(SRC).c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).c -o $@ $(LDFLAGS)

par: $(SRC).par.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par.c -o $@  $(LDFLAGS)

dist: $(SRC).dist.c pi_$(SRC).dist.c sigma_$(SRC).dist.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).dist.c pi_$(SRC).dist.c sigma_$(SRC).dist.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt: $(SRC).distopt.c sigma_$(SRC).distopt.c pi_$(SRC).distopt.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt.c sigma_$(SRC).distopt.c pi_$(SRC).distopt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt_foifi: $(SRC).distopt_foifi.c sigma_$(SRC).distopt_foifi.c pi_$(SRC).distopt_foifi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_foifi.c sigma_$(SRC).distopt_foifi.c pi_$(SRC).distopt_foifi.c\
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt_fop: $(SRC).distopt_fop.c sigma_$(SRC).distopt_fop.c pi_$(SRC).distopt_fop.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_fop.c sigma_$(SRC).distopt_fop.c pi_$(SRC).distopt_fop.c\
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)


clean:
	rm -f out_*  *.par.* \
		*.dist.c  *.distopt.c  *.distopt_foifi.c *.distopt_fop.c \
	  orig opt  par dist distopt distopt_fop distopt_foifi \
		hopt hopt *.par2d.c *.out.* *.dist*.h \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
		.distmem .srcfilename .outfilename .vectorize par2d parsetab.py *.body.c *.pluto.c \
		*.par.cloog *.tiled.cloog *.pluto.cloog sigma.cloog is_receiver.cloog sigma_fop.cloog sigma_check_fop.cloog is_receiver_fop.cloog packunpack.cloog \
		count_remote_dep_tasks.cloog remote_update_dep_tasks.cloog remote_count_dep_tasks.cloog \
		count_local_dep_tasks.cloog local_update_dep_tasks.cloog add_outgoing_edges.cloog local_init_remote_dep_tasks.cloog \
		count_remote_src_tasks.cloog count_local_src_tasks.cloog count_sending_tasks.cloog \
		init_tasks.cloog write_out.cloog compute_task.cloog \
		pi_*.c sigma.c sigma_*.c packunpack.c pi1.c tau.c pi_defs.h *.append.c .appendfilename *.cloog debug_print_node*
