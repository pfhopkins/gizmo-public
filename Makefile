
#-----------------------------------------------------------------
#
# You might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh 
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the 
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the 
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=GIZMO
# 
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype 
#     This file has priority over your shell variable.:
#
#     Uncomment your system in  "Makefile.systype".
#
# If you add an ifeq for a new system below, also add that systype to
# Template-Makefile.systype
#
###########
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel. The code has been modified
#   substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO
#   (dealing with new files and filename conventions, libraries, parser, logic)
#
#############

CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat GIZMO_config.h)

HG_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null)
HG_REPO := $(shell git config --get remote.origin.url)
HG_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
BUILDINFO = "Build on $(HOSTNAME) by $(USER) from $(HG_BRANCH):$(HG_COMMIT) at $(HG_REPO)"
OPT += -DBUILDINFO='$(BUILDINFO)'
OPT += -DGIZMO_SOURCE_DIR='"$(CURDIR)/"'


# initialize some default flags -- these will all get re-written below
CC	= mpicc		# sets the C-compiler (default, will be set for machine below)
CXX	= mpiCC		# sets the C++-compiler (default, will be set for machine below)
FC	= mpif90	# sets the fortran compiler (default, will be set for machine below)
OPTIMIZE = -Wall  -g   # optimization and warning flags (default)
MPICHLIB = -lmpich	# mpi library (arbitrary default, set for machine below)
CHIMESINCL = # default to empty, will only be used below if called
CHIMESLIBS = # default to empty, will only be used below if called
HYPRE_INCL = # hypre library for AMG-preconditioned solver in MG gradient correction
HYPRE_LIBS = # hypre library for AMG-preconditioned solver in MG gradient correction



## read the systype information to use the blocks below for different machines
## precedence: environment SYSTYPE overrides; otherwise ~/.gizmo is checked before Makefile.systype
HOME_GIZMO := $(wildcard $(HOME)/.gizmo)
ifdef SYSTYPE
SYSTYPE := "$(patsubst "%",%,$(SYSTYPE))"
else ifneq ($(HOME_GIZMO),)
include $(HOME_GIZMO)
SYSTYPE := "$(patsubst "%",%,$(SYSTYPE))"
else
include Makefile.systype
endif

ifneq ($(HOME_GIZMO),)
INCL = $(HOME_GIZMO)
else ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif
FINCL =



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Orion-OpenMPI")
CC       =  mpicc
CXX      =  mpicxx -std=c++11
FC       =  mpifort
OPTIMIZE = -g -O3 -march=native -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB = -lstdc++
OPT     += -DHDF5_DISABLE_VERSION_CHECK
## modules: module load gcc openmpi gsl hwloc fftw-mpi hdf5-mpi
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Frontera")
CC       =  mpicc
CXX      =  mpicxx -std=c++17
FC       =  mpif90 -nofor_main
OPTIMIZE = -ggdb -O2 -xCORE-AVX2 -Wno-unknown-pragmas -Wall -Wno-format-security -qopenmp
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I$(TACC_SUNDIALS_INC)
CHIMESLIBS = -L$(TACC_SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial
endif
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
HYPRE_VERSION := $(shell ls /opt/homebrew/Cellar/hypre/ 2>/dev/null | sort -V | tail -n 1)
HYPRE_INCL = -I/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/include/
HYPRE_LIBS = -L/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/lib/ -lHYPRE
endif
MKL_INCL = -I$(TACC_MKL_INC)
MKL_LIBS = -L$(TACC_MKL_LIB) -mkl=sequential
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB) -lgsl -lgslcblas
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
## compiles with module set: intel/19 impi hdf5 fftw3 gsl valgrind python3
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"CaltechHPC")
CC       =  mpicc
CXX      =  mpic++
FC       =  mpif90 -nofor_main
OPTIMIZE = -O2 -xCORE-AVX2
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
endif
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
HYPRE_VERSION := $(shell ls /opt/homebrew/Cellar/hypre/ 2>/dev/null | sort -V | tail -n 1)
HYPRE_INCL = -I/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/include/
HYPRE_LIBS = -L/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/lib/ -lHYPRE
endif
MKL_INCL = -I$(CPATH)
MKL_LIBS = -L$(LIBRARY_PATH) -mkl=sequential
GSL_INCL = -I$(CPATH)
GSL_LIBS = -L$(LIBRARY_PATH)
FFTW_INCL= -I$(CPATH)
FFTW_LIBS= -L$(LIBRARY_PATH)
HDF5INCL = -I$(CPATH) -DH5_USE_16_API
HDF5LIB  = -L$(LIBRARY_PATH) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
# Compiles with following modules:   1) intel/20.1    2) hdf5/1.10.1   3) gsl/2.4       4) fftw/3.3.7
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"BigRed200")
CC       =  cc # For Cray use this instead of mpicc
CXX      =  CC # mpic++
FC       =  $(CC)
OPTIMIZE =  -O2

# Extra compile time warning flags
# OPTIMIZE += -Wall -Wextra -Wuninitialized -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -Wno-unused-variable -Wno-unused-but-set-variable

ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
endif
# All the library paths are already setup
MKL_INCL =
MKL_LIBS = -mkl=sequential
GSL_INCL = -I/N/soft/sles15sp6/gsl/gnu/2.8/include
GSL_LIBS = -L/N/soft/sles15sp6/gsl/gnu/2.8/lib
FFTW_INCL=
FFTW_LIBS=
HDF5INCL =
HDF5LIB  = -lhdf5 -lz
MPICHLIB =

# For mpi intel compilers and HDF5
OPT += -DUSE_MPI_IN_PLACE -DH5_USE_16_API -DNO_ISEND_IRECV_IN_DOMAIN -DHDF5_DISABLE_VERSION_CHECK
## modules to load: (August 2025)
## module swap PrgEnv-gnu PrgEnv-intel
## module swap cray-mpich-ucx/9.0.0 cray-mpich-ucx/8.1.32  # currently defaults to a pre-release mpi version
## module load gsl cray-fftw cray-hdf5
## srun --cpus-per-task=$SLURM_CPUS_PER_TASK --ntasks-per-node=$SLURM_NTASKS_PER_NODE --nodes=$SLURM_JOB_NUM_NODES ./GIZMO ./gizmo_parameters.txt
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Expanse")
CC       = mpicc
CXX      = mpicxx
FC       = $(CC)
OPTIMIZE = -Ofast
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
endif
GSL_INCL = -I$(GSLHOME)/include
GSL_LIBS = -L$(GSLHOME)/lib
FFTW_INCL= -I$(FFTWHOME)/include
FFTW_LIBS= -L$(FFTWHOME)/lib
HDF5INCL = -I$(HDF5HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += #
## modules to load
## module load slurm intel openmpi_ib fftw/2.1.5 gsl hdf5
## run job with
## mpirun -v -x LD_LIBRARY_PATH ./GIZMO params.txt
endif

#----------------------------------------------------------------------------------------------
# Environment for building GIZMO on a macbook with libraries installed via homebrew. But note
# that the specific GSL and HDF5 versions are hardcoded here...
ifeq ($(SYSTYPE),"MacBookCellar")
CC       =  mpicc
CXX      =  mpicxx -std=c++17
FC       =  $(CC) #mpifort  ## change this to "mpifort" for packages requiring linking secondary fortran code, currently -only- the helmholtz eos modules do this, so I leave it un-linked for now to save people the compiler headaches
OPTIMIZE = -O3 -funroll-loops -ffast-math -march=native -flto 
OPTIMIZE += -Wno-unused-command-line-argument ## -g -Wall # compiler warnings
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
#CXX     = mpic++
CHIMESINCL = -I/usr/local/include/sundials
CHIMESLIBS = -L/usr/local/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
OPTIMIZE += -L/opt/homebrew/opt/libomp/lib -lomp
endif
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
ifneq (MHD_MODIFIED_GRADIENT_CG_ONLY,$(findstring MHD_MODIFIED_GRADIENT_CG_ONLY,$(CONFIGVARS)))
HYPRE_VERSION := $(shell ls /opt/homebrew/Cellar/hypre/ 2>/dev/null | sort -V | tail -n 1)
HYPRE_INCL = -I/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/include/
HYPRE_LIBS = -L/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/lib/ -lHYPRE
endif
endif
ifeq (MHD_MODIFIED_GRADIENT_USE_PARDISO,$(findstring MHD_MODIFIED_GRADIENT_USE_PARDISO,$(CONFIGVARS)))
MKL_INCL = -I$(MKLROOT)/include
MKL_LIBS = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else
MKL_INCL = #
MKL_LIBS = #
endif
GSL_INCL = -I/opt/homebrew/Cellar/gsl/2.8/include #-I$(PORTINCLUDE)
GSL_LIBS = -L/opt/homebrew/Cellar/gsl/2.8/lib #-L$(PORTLIB)
FFTW_INCL= -I/opt/homebrew/Cellar/fftw/3.3.10_3/include
FFTW_LIBS= -L/opt/homebrew/Cellar/fftw/3.3.10_3/lib
HDF5_VERSION := $(shell ls /opt/homebrew/Cellar/hdf5/ 2>/dev/null | sort -V | tail -n 1)
HDF5INCL = -I/opt/homebrew/Cellar/hdf5/$(HDF5_VERSION)/include -DH5_USE_16_API  #-I$(PORTINCLUDE) -DH5_USE_16_API
HDF5LIB  = -L/opt/homebrew/Cellar/hdf5/$(HDF5_VERSION)/lib -lhdf5 -lz  #-L$(PORTLIB)
MPICHLIB = #
OPT     += -DDISABLE_ALIGNED_ALLOC -DCHIMES_USE_DOUBLE_PRECISION #
endif

#----------------------------
ifeq ($(SYSTYPE),"github-ubuntu")
CC       =  mpicc
CXX      =  mpicxx
FC       =  $(CC)
OPTIMIZE = -g -fcommon -O1 -funroll-loops -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -fipa-cp-clone  ## optimizations for gcc compilers (1/2)
OPTIMIZE += -ftree-loop-distribute-patterns -fvect-cost-model -ftree-partial-pre   ## optimizations for gcc compilers (2/2)
OPTIMIZE += -g -Wall # compiler warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CXX     = mpic++
CHIMESINCL = -I/usr/include/sundials
CHIMESLIBS = -L/usr/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp # openmp required compiler flags
FC       = $(CC)
endif
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib 
FFTW_INCL= -I/usr/include
FFTW_LIBS= -L/usr/lib
HDF5INCL = -I/usr/include/hdf5/openmpi -DH5_USE_16_API
HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ -lhdf5 -lz
MPICHLIB = #
OPT     += -DDISABLE_ALIGNED_ALLOC -DCHIMES_USE_DOUBLE_PRECISION
## to get required packages: sudo apt install libhdf5-openmpi-dev libgsl-dev libopenmpi-dev
endif

#----------------------------
# Should work on any Flatiron institute linux cluster environment: rusty, popeye and linux workstations
ifeq ($(SYSTYPE),"RUSTY")
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpicxx
endif
FC      = mpifort
OPTIMIZE =  -O3 -ffast-math -funroll-loops -march=native -g -Wall
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
GSL_INCL = -I$(GSL_BASE)/include
GSL_LIBS = -L$(GSL_BASE)/lib -Xlinker -R -Xlinker $(GSL_BASE) -lgsl -lgslcblas
FFTW3_BASE= /mnt/sw/nix/store/bjzkf3pwcw0gy54db19kd4rl0xdiq98s-fftw-3.3.10/.
FFTW_INCL= -I$(FFTW3_BASE)/include
FFTW_LIBS= -L$(FFTW3_BASE)/lib -Xlinker -R -Xlinker $(FFTW3_BASE)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_BASE)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_BASE)/lib -Xlinker -R -Xlinker $(HDF5_BASE)/lib -lhdf5 -lz
endif


#================================================================================
# Rusty (Flatiron) -- per-architecture systypes.
#
# The generic "RUSTY" systype above uses -march=native, which requires building on
# the same node type you run on. The three blocks below name the target explicitly
# so you can build anywhere (e.g. a workstation or a login node) and still get a
# correctly-tuned binary. Pick the one matching your --constraint.
#
#   systype         --constraint     cores/node   RAM/node   ISA target
#   RUSTY_ROME      rome                128        1.0 TB    znver2  (Zen 2)
#   RUSTY_GENOA     genoa                96        1.5 TB    znver4  (Zen 4)
#   RUSTY_ICELAKE   icelake&cpu          64        1.0 TB    icelake-server
#
# ---- HOW TO RUN (measured 2026-08-02, 2.5e6-cell STARFORGE full-physics box) ----
#
# 1. USE PURE MPI. One rank per core, OPENMP off. Throughput in a fixed wall-clock
#    budget, best config per node type:
#        Genoa  96x1 -> 24 steps      (best absolute AND best per core)
#        Rome  128x1 -> 16 steps
#        Ice    64x1 -> 13 steps
#    Per-core cost (core-seconds/step, lower better): Genoa ~500, Ice ~700, Rome ~1150.
#    Genoa wins outright despite having fewer cores than Rome -- prefer it when free.
#
# 2. HYBRID MPI+OpenMP COSTS YOU THROUGHPUT. Relative to pure MPI on the same node:
#    2 threads/rank is roughly free (Rome ties), then it degrades -- ~1.5-2x slower
#    by 8 threads/rank, ~2.5x by 16. Only go hybrid if you need fewer ranks for
#    memory reasons; 2 threads/rank is the sweet spot if you must.
#
# 3. IF YOU BUILD WITH OPENMP, YOU MUST PASS --cpu-bind=cores TO srun:
#        srun --cpus-per-task=$SLURM_CPUS_PER_TASK --cpu-bind=cores ./GIZMO ...
#    OMP_PROC_BIND/OMP_PLACES are honored only by OpenMP-enabled binaries. Without
#    a per-task cpuset each rank's libgomp independently sees the whole node and
#    pins its threads to the same first N cores, so only N cores of the node do any
#    work while they are R-way oversubscribed. This cost us a factor of 10-49x and
#    is completely silent -- pure-MPI builds are immune, so it masquerades as
#    "OpenMP is slow". Always dump the placement and check it:
#        srun --cpu-bind=cores --label grep Cpus_allowed_list /proc/self/status
#    Every rank must report a DISTINCT cpu list.
#
# 4. SCALE MaxMemSize WITH RANK COUNT. It is per-MPI-process, so a value tuned for
#    128 ranks starves an 8-rank run of the same problem. Budget ~75% of node RAM
#    divided by ranks-per-node.
#
# 5. COMPILERS. Spread is small (~8% end to end): intel-oneapi 2025 was fastest on
#    all three architectures, aocc/5.0 a close second, gcc/13-14 ~5-8% behind.
#    Load a compiler and point the MPI wrappers at it:
#        module load modules gcc/14.2.0 openmpi gsl hdf5/mpi-1.12.3
#        export OMPI_CC=gcc   OMPI_CXX=g++         # gcc
#        export OMPI_CC=clang OMPI_CXX=clang++     # aocc/5.0.0
#        export OMPI_CC=icx   OMPI_CXX=icpx        # intel-oneapi-compilers/2025.0.4
#    NOTE gcc/13.2.0 only resolves under modules/2.3-20240529; load that explicitly
#    or the build and the run will link against different GSL and fail at startup.
#    NOTE icpx at -O2 and above miscompiled the H2 rotational partition sum into an
#    infinite loop (fixed in eos/hydrogen_molecule.cc). If you use the Intel
#    compiler, sanity-check that a run actually advances past the first few steps.
#
# 6. EXTRA FLAGS. Append with OPT_EXTRA, e.g.
#        make SYSTYPE='"RUSTY_GENOA"' OPT_EXTRA="-flto"
#    Later flags win over earlier ones, so this also un-sets defaults:
#        OPT_EXTRA="-O2"                  # drop to -O2
#        OPT_EXTRA="-fno-fast-math"       # disable fast-math entirely
#        OPT_EXTRA="-fno-finite-math-only" # keep fast-math but restore NaN/Inf semantics
#        OPT_EXTRA="-mprefer-vector-width=256"  # avoid AVX-512 downclocking
#================================================================================

#----------------------------
# Rusty AMD Rome (EPYC, Zen 2) -- 128 cores/node, 1 TB. Run: --constraint=rome, 128x1.
ifeq ($(SYSTYPE),"RUSTY_ROME")
CC       =   mpicc
CXX      =   mpicxx
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpicxx
endif
FC       =   mpifort
OPTIMIZE =  -O3 -ffast-math -funroll-loops -march=znver2 -mtune=znver2 -g -Wall $(OPT_EXTRA)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp   # remember --cpu-bind=cores on srun (see note 3 above)
endif
GSL_INCL = -I$(GSL_BASE)/include
GSL_LIBS = -L$(GSL_BASE)/lib -Xlinker -R -Xlinker $(GSL_BASE) -lgsl -lgslcblas
FFTW3_BASE= /mnt/sw/nix/store/bjzkf3pwcw0gy54db19kd4rl0xdiq98s-fftw-3.3.10/.
FFTW_INCL= -I$(FFTW3_BASE)/include
FFTW_LIBS= -L$(FFTW3_BASE)/lib -Xlinker -R -Xlinker $(FFTW3_BASE)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_BASE)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_BASE)/lib -Xlinker -R -Xlinker $(HDF5_BASE)/lib -lhdf5 -lz
endif


#----------------------------
# Rusty AMD Genoa (EPYC, Zen 4) -- 96 cores/node, 1.5 TB. Run: --constraint=genoa, 96x1.
# FASTEST OF THE THREE, both in absolute throughput and per core. Prefer this.
# Genoa has AVX-512; if a kernel turns out to be downclock-limited, compare
# OPT_EXTRA="-mprefer-vector-width=256".
ifeq ($(SYSTYPE),"RUSTY_GENOA")
CC       =   mpicc
CXX      =   mpicxx
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpicxx
endif
FC       =   mpifort
OPTIMIZE =  -O3 -ffast-math -funroll-loops -march=znver4 -mtune=znver4 -g -Wall $(OPT_EXTRA)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp   # remember --cpu-bind=cores on srun (see note 3 above)
endif
GSL_INCL = -I$(GSL_BASE)/include
GSL_LIBS = -L$(GSL_BASE)/lib -Xlinker -R -Xlinker $(GSL_BASE) -lgsl -lgslcblas
FFTW3_BASE= /mnt/sw/nix/store/bjzkf3pwcw0gy54db19kd4rl0xdiq98s-fftw-3.3.10/.
FFTW_INCL= -I$(FFTW3_BASE)/include
FFTW_LIBS= -L$(FFTW3_BASE)/lib -Xlinker -R -Xlinker $(FFTW3_BASE)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_BASE)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_BASE)/lib -Xlinker -R -Xlinker $(HDF5_BASE)/lib -lhdf5 -lz
endif


#----------------------------
# Rusty Intel Ice Lake -- 64 cores/node, 1 TB. Run: --constraint='icelake&cpu', 64x1.
# NB use 'icelake&cpu': the bare 'icelake' feature also matches the A100/H100 GPU
# nodes, and a CPU-only job should not tie those up.
ifeq ($(SYSTYPE),"RUSTY_ICELAKE")
CC       =   mpicc
CXX      =   mpicxx
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpicxx
endif
FC       =   mpifort
OPTIMIZE =  -O3 -ffast-math -funroll-loops -march=icelake-server -mtune=icelake-server -g -Wall $(OPT_EXTRA)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp   # remember --cpu-bind=cores on srun (see note 3 above)
endif
GSL_INCL = -I$(GSL_BASE)/include
GSL_LIBS = -L$(GSL_BASE)/lib -Xlinker -R -Xlinker $(GSL_BASE) -lgsl -lgslcblas
FFTW3_BASE= /mnt/sw/nix/store/bjzkf3pwcw0gy54db19kd4rl0xdiq98s-fftw-3.3.10/.
FFTW_INCL= -I$(FFTW3_BASE)/include
FFTW_LIBS= -L$(FFTW3_BASE)/lib -Xlinker -R -Xlinker $(FFTW3_BASE)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_BASE)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_BASE)/lib -Xlinker -R -Xlinker $(HDF5_BASE)/lib -lhdf5 -lz
endif


#----------------------------
# PSMN (ENS Lyon) - Debian 13, e5-2667v4 nodes
# Compiler: system OpenMPI 5.0.7 (Debian); libs: EasyBuild gompi-2025a toolchain
ifeq ($(SYSTYPE),"PSMN")
CC       =  mpicc
CXX      =  mpicxx
FC       =  mpifort
OPTIMIZE = -O2 -g -Wall -funroll-loops -finline-functions
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas
FFTW_BASE= /applis/PSMN/debian13/E5/software/FFTW.MPI/3.3.10-gompi-2025a
FFTW_INCL= -I$(FFTW_BASE)/include
FFTW_LIBS= -L$(FFTW_BASE)/lib -Wl,-rpath,$(FFTW_BASE)/lib
HDF5_BASE= /applis/PSMN/debian13/E5/software/HDF5/1.14.6-gompi-2025a
HDF5INCL = -I$(HDF5_BASE)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_BASE)/lib -Wl,-rpath,$(HDF5_BASE)/lib -lhdf5 -lz
MPICHLIB =
OPT     += -DDISABLE_ALIGNED_ALLOC
endif


#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------


#
# different code groups that need to be compiled. the groupings below are
# arbitrary (they will all be added to OBJS and compiled, and if they are
# un-used it should be through use of macro flags in the source code). But
# they are grouped below for the sake of clarity when adding/removing code
# blocks in the future
#
CORE_OBJS =	core/main.o core/accel.o core/timestep.o core/init.o file_io/restart.o file_io/io.o \
			core/predict.o declarations/global.o core/begrun.o core/run.o declarations/allvars.o \
			file_io/read_ic.o domain/domain.o core/driftfac.o core/kicks.o mesh/ngb.o \
			compile_time_info.o mesh/merge_split.o \
			core/transport_subcycle.o

SYSTEM_OBJS =   system/system.o \
				system/allocate.o \
				system/mymalloc.o \
				system/parallel_sort.o \
                system/peano.o \
                system/parallel_sort_special.o \
                system/mpi_util.o \
                system/pinning.o

GRAVITY_OBJS  = gravity/forcetree.o \
                gravity/forcetree_update.o \
                gravity/gravtree.o \
				gravity/cosmology.o \
				gravity/potential.o \
				gravity/pm_periodic.o \
                gravity/pm_nonperiodic.o \
                gravity/longrange.o \
                gravity/ags_rkern.o \
                gravity/binary.o \
                gravity/star_direct_gravity.o

HYDRO_OBJS = 	hydro/hydro_toplevel.o \
				hydro/density.o \
				hydro/gradients.o \
				hydro/mg_gradient_correction.o \
				turb/dynamic_diffusion.o \
				turb/dynamic_diffusion_velocities.o \
				turb/turb_driving.o \
				turb/turb_powerspectra.o

EOSCOOL_OBJS =  cooling/cooling.o \
				cooling/grackle.o \
				cooling/simple_chemistry.o \
				eos/eos.o \
				eos/hydrogen_molecule.o \
				eos/cosmic_ray_fluid/cosmic_ray_alfven.o \
				eos/cosmic_ray_fluid/cosmic_ray_utilities.o \
				solids/elastic_physics.o \
				solids/grain_physics.o \
				solids/ism_dust_chemistry.o

STARFORM_OBJS = galaxy_sf/sfr_eff.o \
                galaxy_sf/stellar_evolution.o \
                galaxy_sf/mechanical_fb.o \
                galaxy_sf/thermal_fb.o \
                galaxy_sf/radfb_local.o \
                galaxy_sf/dm_dispersion_rkern.o

SINK_OBJS = sinks/sink.o \
            sinks/sink_util.o \
            sinks/sink_environment.o \
            sinks/sink_feed.o \
            sinks/sink_swallow_and_kick.o

RHD_OBJS =  radiation/rt_utilities.o \
			radiation/rt_CGmethod.o \
			radiation/rt_source_injection.o \
			radiation/rt_chem.o \
			radiation/rt_dust_opacity.o 

FOF_OBJS =	structure/fof.o \
			structure/subfind/subfind.o \
			structure/subfind/subfind_vars.o \
			structure/subfind/subfind_collective.o \
			structure/subfind/subfind_serial.o \
			structure/subfind/subfind_so.o \
			structure/subfind/subfind_cont.o \
			structure/subfind/subfind_distribute.o \
			structure/subfind/subfind_findlinkngb.o \
			structure/subfind/subfind_nearesttwo.o \
			structure/subfind/subfind_loctree.o \
			structure/subfind/subfind_potential.o \
			structure/subfind/subfind_density.o \
			structure/twopoint.o \
			structure/lineofsight.o

MISC_OBJS = sidm/cbe_integrator.o \
			sidm/dm_fuzzy.o \
			sidm/sidm_core.o

## hypre, for MHD_MODIFIED_GRADIENT's AMG-preconditioned solver. mg_gradient_correction.cc
## includes HYPRE.h unless MHD_MODIFIED_GRADIENT_CG_ONLY or _USE_PARDISO is set, and there is no
## hypre module on Rusty (nor one bundled with petsc), so fall back to a local build --
## scripts/build_hypre_rusty.sh makes one. Only fills in when no systype block above already set a
## path, so the homebrew prefixes stay intact. Harmless when the prefix does not exist: only a
## config that actually includes HYPRE.h cares.
ifeq ($(strip $(HYPRE_INCL)),)
ifeq ($(strip $(HYPRE_PATH)),)
HYPRE_PATH := $(HOME)/opt/hypre
endif
HYPRE_INCL = -I$(HYPRE_PATH)/include
HYPRE_LIBS = -L$(HYPRE_PATH)/lib -L$(HYPRE_PATH)/lib64 -Xlinker -R -Xlinker $(HYPRE_PATH)/lib -Xlinker -R -Xlinker $(HYPRE_PATH)/lib64 -lHYPRE
endif

## name of executable and optimizations
EXEC   = GIZMO
OPTIONS = $(OPTIMIZE) $(OPT)

.DEFAULT_GOAL := $(EXEC)

ifeq (JACO,$(findstring JACO,$(CONFIGVARS)))
JACO_MODEL := $(strip $(shell grep '^\#define JACO' GIZMO_config.h | awk '{print $$3}'))
ifeq ($(JACO_MODEL),)
JACO_MODEL := wind_comparison
endif
JACO_GENERATED = cooling/jaco_eos.cc cooling/microphysics_func_jac.cc cooling/microphysics_func_jac.h cooling/jaco_interp.h
EOSCOOL_OBJS += cooling/jaco.o cooling/jaco_eos.o cooling/microphysics_func_jac.o
endif

## combine all the objects above
OBJS  = $(CORE_OBJS) $(SYSTEM_OBJS) $(GRAVITY_OBJS) $(HYDRO_OBJS) \
		$(EOSCOOL_OBJS) $(STARFORM_OBJS) $(SINK_OBJS) $(RHD_OBJS) \
		$(FOF_OBJS) $(MISC_OBJS)

## fortran recompiler block
FOPTIONS = $(OPTIMIZE) $(FOPT)
FOBJS =

## include files needed at compile time for the above objects
INCL    += 	declarations/allvars.h \
			core/proto.h \
			gravity/forcetree.h \
			gravity/myfftw3.h \
			domain/domain.h \
			system/myqsort.h \
			mesh/kernel.h \
			eos/eos.h \
			sinks/sink.h \
			structure/fof.h \
			structure/subfind/subfind.h \
			cooling/cooling.h \
			Makefile

ifeq (JACO,$(findstring JACO,$(CONFIGVARS)))
INCL += cooling/microphysics_func_jac.h cooling/jaco_interp.h
JACO_PYTHON ?= python3
$(JACO_GENERATED): $(CONFIG)
	$(JACO_PYTHON) -m jaco.codegen.gizmo.gizmo $(JACO_MODEL) --language c --ext .cc --output-dir cooling
endif


## now we add special cases dependent on compiler flags. normally we would
##  include the files always, and simply use the in-file compiler variables
##  to determine whether certain code is compiled [this allows us to take
##  advantage of compiler logic, and makes it easier for the user to
##  always specify what they want]. However special cases can arise, if e.g.
##  there are certain special libraries needed, or external compilers, for
##  certain features

# helmholtz eos routines need special treatment here because they are written
#  in fortran and call the additional fortran compilers and linkers. these could
#  be written to always compile and just be ignored, but then the large majority
#  of cases that -don't- need the fortran linker would always have to go
#  through these additional compilation options and steps (and this 
#  can cause additional problems on some machines). so we sandbox it here.
ifeq (EOS_HELMHOLTZ,$(findstring EOS_HELMHOLTZ,$(CONFIGVARS)))
OBJS    += eos/eos_interface.o
INCL    += eos/helmholtz/helm_wrap.h
FOBJS   += eos/helmholtz/helm_impl.o eos/helmholtz/helm_wrap.o
FINCL   += eos/helmholtz/helm_const.dek eos/helmholtz/helm_implno.dek eos/helmholtz/helm_table_storage.dek eos/helmholtz/helm_vector_eos.dek
endif

ifeq (EOS_ANEOS,$(findstring EOS_ANEOS,$(CONFIGVARS)))
OBJS    += eos/aneos.o
INCL    += eos/aneos.h
endif

# chimes files are treated as special for now because they require special external libraries (e.g. sundials) that are otherwise not
#   used anywhere else in the code, and have not had their macro logic cleaned up to allow appropriate compilation without chimes flags enabled
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
OBJS    += cooling/chimes/chimes.o cooling/chimes/chimes_cooling.o cooling/chimes/init_chimes.o cooling/chimes/rate_equations.o cooling/chimes/update_rates.o 
INCL    += cooling/chimes/chimes_interpol.h cooling/chimes/chimes_proto.h cooling/chimes/chimes_vars.h 
endif

# if grackle libraries are installed they must be a shared library as defined here
ifeq (COOL_GRACKLE,$(findstring COOL_GRACKLE,$(CONFIGVARS)))
OPTIONS += -DCONFIG_BFLOAT_8
GRACKLEINCL =
GRACKLELIBS = -lgrackle
else
GRACKLEINCL =
GRACKLELIBS =
endif

# KETJU regularized integrator for BH and stellar dynamics
ifeq (KETJU_REGULARIZATION,$(findstring KETJU_REGULARIZATION,$(CONFIGVARS)))
OBJS    += galaxy_sf/ketju_coupling.o
KETJU_INCL = -Iketju-integrator/include
KETJU_LIBS = -Lketju-integrator/lib -lketju-integrator -lstdc++
else
KETJU_INCL =
KETJU_LIBS =
endif

# linking libraries (includes machine-dependent options above)
CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) \
         $(GRACKLEINCL) $(CHIMESINCL) $(HYPRE_INCL) $(MKL_INCL) $(KETJU_INCL)



# one annoying thing here is the FFTW libraries, since they are named differently depending on
#  whether they are compiled in different precision levels, or with different parallelization options, so we
#  have to have a big block here 'sorting them out'.
#
fftw_on_key = # default to 'off'
ifeq (PMGRID,$(findstring PMGRID, $(CONFIGVARS)))
  fftw_on_key = yes # needed for this module
endif
ifeq (TURB_DRIVING_SPECTRUMGRID,$(findstring TURB_DRIVING_SPECTRUMGRID, $(CONFIGVARS)))
  fftw_on_key = yes  # needed for this module
endif
FFTW_LIBNAMES = # default to 'off'
ifdef fftw_on_key
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -lfftw3_mpi -lfftw3 # double-precision libraries
else
  FFTW_LIBNAMES = -lfftw3f_mpi -lfftw3f # single-precision libraries
endif
endif


LIBS = $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas \
	   $(FFTW_LIBS) $(FFTW_LIBNAMES) -lm $(GRACKLELIBS) $(CHIMESLIBS) $(HYPRE_LIBS) $(MKL_LIBS) $(KETJU_LIBS)


$(EXEC): $(OBJS)
	$(CXX) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): %.o: %.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(CXX) $(CFLAGS) -c $< -o $@

#$(FOBJS): %.o: %.f90
#	$(FC) $(OPTIMIZE) -c $< -o $@

compile_time_info.cc: $(CONFIG)
	$(PERL) file_io/prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) *.oo *.c~ compile_time_info.cc GIZMO_config.h
	rm -f cooling/jaco_eos.cc cooling/microphysics_func_jac.cc cooling/microphysics_func_jac.h cooling/jaco_interp.h


