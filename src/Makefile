
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
CONFIGVARS += OPENMP

HG_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null)
HG_REPO := $(shell git config --get remote.origin.url)
HG_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
BUILDINFO = "Build on $(HOSTNAME) by $(USER) from $(HG_BRANCH):$(HG_COMMIT) at $(HG_REPO)"
OPT += -DBUILDINFO='$(BUILDINFO)'
OPT += -DGIZMO_SOURCE_DIR='"$(CURDIR)/"'


# initialize some default flags -- these will all get re-written below
CC	= mpicc		# sets the C-compiler (default, will be set for machine below)
CXX	= mpiCC		# sets the C++-compiler (default, will be set for machine below)
OPTIMIZE = -Wall  -g   # optimization and warning flags (default)
MPICHLIB = -lmpich	# mpi library (arbitrary default, set for machine below)
GPU_CFLAGS = # GPU offload compiler flags (set for GPU systypes below)
GPU_LDFLAGS = # GPU offload linker flags (set for GPU systypes below)
KOKKOS_PATH = # path to Kokkos installation (set for GPU systypes below)
KOKKOS_CPPFLAGS = # populated by Kokkos Makefile.kokkos include below
KOKKOS_CXXFLAGS = #
KOKKOS_LDFLAGS  = #
KOKKOS_LIBS     = #
# CHIMES SUNDIALS/CVODE requirement (READ BEFORE ENABLING CHIMES):
#   The bundled CHIMES network uses the pre-6.0 SUNDIALS API -- context-free
#   CVodeCreate()/N_VNew_Serial()/SUNDenseMatrix()/SUNLinSol_Dense() and the
#   'realtype' type. It therefore requires SUNDIALS in the 2.7 - 5.x range.
#   SUNDIALS 6.x made a SUNContext argument mandatory on those calls, and 7.x
#   removed 'realtype' (now 'sunrealtype'), so 6.0+ will NOT compile against the
#   current CHIMES sources without source-level API updates. Point CHIMESINCL/
#   CHIMESLIBS below at a SUNDIALS install in the supported range.
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



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Frontera")
CC       =  mpicc
CXX      =  mpicxx -std=c++17
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
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
## compiles with module set: intel/19 impi hdf5 fftw3 valgrind python3
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Vista")
CC       =  mpicc
CXX      =  mpicxx -std=c++17
FC       =  mpif90
OPTIMIZE = -O2 -Wall
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
## Kokkos GPU offload for cooling.cc (and future GPU-ported files).
## cooling.cc and eos/eos.cc are compiled via nvcc_wrapper → nvcc for device code.
## TACC's kokkos/4.5.01-cuda module sets TACC_KOKKOS_DIR/INC/LIB/BIN.
## We set flags manually (no Makefile.kokkos — TACC's Kokkos is CMake-installed).
## Load modules: nvidia/25.9 cuda kokkos/4.5.01-cuda openmpi hdf5/2.0.0 fftw3
KOKKOS_PATH    = $(TACC_KOKKOS_DIR)
## Compile flags: Kokkos includes + CUDA relaxed-constexpr/lambda extensions + sm_90 arch
KOKKOS_CPPFLAGS = -I$(TACC_KOKKOS_INC)
KOKKOS_CXXFLAGS = --expt-relaxed-constexpr --expt-extended-lambda -arch=sm_90
## Link flags: Kokkos libs (core + containers) + CUDA runtime
KOKKOS_LDFLAGS  = -L$(TACC_KOKKOS_LIB) -Wl,-rpath,$(TACC_KOKKOS_LIB) -L$(TACC_CUDA_LIB) -Wl,-rpath,$(TACC_CUDA_LIB)
KOKKOS_LIBS     = -lkokkoscore -lkokkoscontainers -lcudart -lcuda
## nvcc_wrapper reads NVCC_WRAPPER_DEFAULT_COMPILER as its host compiler.
## Set it to mpicxx so nvcc_wrapper → nvcc -ccbin mpicxx, bringing in MPI headers/libs.
export NVCC_WRAPPER_DEFAULT_COMPILER = mpicxx
GPU_CXX    = $(TACC_KOKKOS_BIN)/nvcc_wrapper --std=c++17
GPU_CFLAGS = $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS)
GPU_LDFLAGS = $(KOKKOS_LDFLAGS)
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I$(TACC_SUNDIALS_INC)
CHIMESLIBS = -L$(TACC_SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial
endif
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
## Vista: `module load hypre` provides TACC_HYPRE_INC + TACC_HYPRE_LIB.
## TACC_HYPRE_GPU_DIR honored if set (future CUDA-Hypre build).
ifdef TACC_HYPRE_GPU_DIR
HYPRE_INCL = -I$(TACC_HYPRE_GPU_DIR)/include
HYPRE_LIBS = -L$(TACC_HYPRE_GPU_DIR)/lib -lHYPRE
else ifdef TACC_HYPRE_INC
HYPRE_INCL = -I$(TACC_HYPRE_INC)
HYPRE_LIBS = -L$(TACC_HYPRE_LIB) -lHYPRE
else ifdef TACC_HYPRE_DIR
HYPRE_INCL = -I$(TACC_HYPRE_DIR)/include
HYPRE_LIBS = -L$(TACC_HYPRE_DIR)/lib -lHYPRE
endif
endif
MKL_INCL = #-I$(TACC_MKL_INC)
MKL_LIBS = #-L$(TACC_MKL_LIB) -lmkl_rt
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
## submit to 'gh' or 'gh-dev' queue on Vista (NVIDIA Grace Hopper H200, sm_90)
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Frontier")
## Frontier (OLCF), AMD MI250X (gfx90a) via the Kokkos HIP backend.
## Prebuilt module kokkos/4.6.02-gpu (Kokkos_ENABLE_HIP, ARCH_VEGA90A+ZEN3, built with
## rocm-6.4.2 hipcc, shared libs, RDC off). Device TUs (GPU_OBJS) are compiled by hipcc
## (defines __HIPCC__ => GIZMO_GPU_COMPILER, like nvcc_wrapper on Vista); host TUs + the
## final link use the Cray CC wrapper (amd PE), which auto-adds GPU-aware MPI (GTL,
## -lmpi_gtl_hsa), -lamdhip64, and cray-hdf5/fftw/gsl paths for the host objects.
## Load: amd/6.4.2 rocm/6.4.2 cray-mpich/9.0.1 craype-accel-amd-gfx90a kokkos/4.6.02-gpu
##       cray-hdf5 cray-fftw gsl hypre/3.0.0-mpi
## (hypre only needed for MHD_MODIFIED_GRADIENT; see the HYPRE block at the end of this
## section. cray-mpich/9.0.1 is what the hypre module is built against.)
CC       =  cc
CXX      =  CC -std=c++17
FC       =  ftn
OPTIMIZE = -O2 -g   # -g does not change codegen at -O2; without it a Frontier core file has no line numbers
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
KOKKOS_PATH     = $(OLCF_KOKKOS_ROOT)
KOKKOS_CPPFLAGS = -I$(OLCF_KOKKOS_ROOT)/include
## hipcc device flags: target gfx90a + cray-mpich headers (device TUs are -c only; CC links).
## <mpi.h> is included directly by device TUs but hipcc lacks the Cray wrapper's auto include.
## gfx90a:xnack+ (not bare gfx90a): the code uses HIP managed/unified memory
## (static __managed__ mirrors + Kokkos SharedSpace=HIPManagedSpace), which is
## only host-accessible on MI250X with XNACK. Requires HSA_XNACK=1 at runtime too.
KOKKOS_CXXFLAGS = --offload-arch=gfx90a:xnack+ -march=znver3 -I$(CRAY_MPICH_DIR)/include
KOKKOS_LDFLAGS  = -L$(OLCF_KOKKOS_ROOT)/lib64 -Wl,-rpath,$(OLCF_KOKKOS_ROOT)/lib64
KOKKOS_LIBS     = -lkokkoscore -lkokkoscontainers
GPU_CXX     = hipcc --std=c++17
GPU_CFLAGS  = $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS)
GPU_LDFLAGS = $(KOKKOS_LDFLAGS)
MKL_INCL =
MKL_LIBS =
## cray-hdf5 include set explicitly so hipcc finds <hdf5.h> (allvars.h includes it
## unconditionally); the Cray CC wrapper auto-adds it for host TUs, but hipcc does not.
HDF5INCL = -I$(HDF5_ROOT)/include -DH5_USE_16_API
HDF5LIB  = -lhdf5 -lz
FFTW_INCL=
FFTW_LIBS=
GSL_INCL = -I$(OLCF_GSL_ROOT)/include
GSL_LIBS = -L$(OLCF_GSL_ROOT)/lib
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
## Frontier: `module load hypre/3.0.0-mpi` sets OLCF_HYPRE_ROOT (the OLCF module
## convention, as for gsl/kokkos above); its libraries are under lib64. That version
## builds against amd/6.4.2 + cray-mpich/9.0.1, which is the compiler the
## kokkos/4.6.02-gpu module also requires; the newer hypre modules need amd/7.x and so
## cannot be combined with that Kokkos build. HYPRE_ROOT/HYPRE_DIR are honored as
## fallbacks for a hand-built install; both lib and lib64 are searched.
ifdef OLCF_HYPRE_ROOT
HYPRE_INCL = -I$(OLCF_HYPRE_ROOT)/include
HYPRE_LIBS = -L$(OLCF_HYPRE_ROOT)/lib -L$(OLCF_HYPRE_ROOT)/lib64 -lHYPRE
else ifdef HYPRE_ROOT
HYPRE_INCL = -I$(HYPRE_ROOT)/include
HYPRE_LIBS = -L$(HYPRE_ROOT)/lib -L$(HYPRE_ROOT)/lib64 -lHYPRE
else ifdef HYPRE_DIR
HYPRE_INCL = -I$(HYPRE_DIR)/include
HYPRE_LIBS = -L$(HYPRE_DIR)/lib -L$(HYPRE_DIR)/lib64 -lHYPRE
endif
endif
## submit to Frontier batch (AMD MI250X, 8x GCD per node, gfx90a)
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Frontier_CPU")
## Frontier (OLCF) with the Kokkos OpenMP backend: the same Kokkos source path as
## "Frontier", dispatched to CPU threads instead of the MI250X. Host compiler, host
## libraries and host flags are IDENTICAL to the "Frontier" block, so a Frontier vs
## Frontier_CPU pair differs only in the Kokkos backend and in which compiler builds
## GPU_OBJS -- that is the point of this systype, and it is why OPTIMIZE deliberately
## does NOT add -march=znver3 even though this is a host-only build.
## hipcc is not used, so __HIPCC__ (and therefore GIZMO_GPU_COMPILER) is undefined:
## the per-TU AllDeviceMirror blocks are skipped and All stays the normal extern,
## exactly as in Vista_CPU and Frontera_Kokkos.
## Modules:
##   module reset
##   module load amd/6.4.2
##   module load cray-mpich/9.0.1
##   module load kokkos/4.6.01-omp
##   module load cray-hdf5 cray-fftw
##   module load hypre/3.0.0-mpi
## i.e. the "Frontier" module set MINUS rocm and craype-accel-amd-gfx90a, with the
## OpenMP-backend Kokkos module in place of kokkos/4.6.02-gpu. (hypre is only needed
## for MHD_MODIFIED_GRADIENT; no gsl module -- nothing here consumes GSL_INCL/GSL_LIBS.)
## Do NOT export HSA_XNACK=1 for this build; it has no device allocations.
## To compare against the "Frontier" GPU build, run at the same rank and thread counts.
CC       =  cc
CXX      =  CC -std=c++17
FC       =  ftn
OPTIMIZE = -O2 -g   # matches the "Frontier" block exactly; see the note above
## OpenMP is required by the Kokkos OpenMP backend itself, so it is unconditional here
## rather than gated on the OPENMP config flag (same reasoning as Frontera_Kokkos).
OPTIMIZE += -fopenmp
KOKKOS_PATH     = $(OLCF_KOKKOS_ROOT)
KOKKOS_CPPFLAGS = -I$(OLCF_KOKKOS_ROOT)/include
KOKKOS_LDFLAGS  = -L$(OLCF_KOKKOS_ROOT)/lib64 -Wl,-rpath,$(OLCF_KOKKOS_ROOT)/lib64
KOKKOS_LIBS     = -lkokkoscore -lkokkoscontainers
## GPU-ported TUs: the host compiler, plus the Kokkos includes. No offload arch.
GPU_CXX     = CC -std=c++17
GPU_CFLAGS  = $(KOKKOS_CPPFLAGS)
GPU_LDFLAGS = $(KOKKOS_LDFLAGS)
MKL_INCL =
MKL_LIBS =
HDF5INCL = -I$(HDF5_ROOT)/include -DH5_USE_16_API
HDF5LIB  = -lhdf5 -lz
FFTW_INCL=
FFTW_LIBS=
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
ifdef OLCF_HYPRE_ROOT
HYPRE_INCL = -I$(OLCF_HYPRE_ROOT)/include
HYPRE_LIBS = -L$(OLCF_HYPRE_ROOT)/lib -L$(OLCF_HYPRE_ROOT)/lib64 -lHYPRE
else ifdef HYPRE_ROOT
HYPRE_INCL = -I$(HYPRE_ROOT)/include
HYPRE_LIBS = -L$(HYPRE_ROOT)/lib -L$(HYPRE_ROOT)/lib64 -lHYPRE
else ifdef HYPRE_DIR
HYPRE_INCL = -I$(HYPRE_DIR)/include
HYPRE_LIBS = -L$(HYPRE_DIR)/lib -L$(HYPRE_DIR)/lib64 -lHYPRE
endif
endif
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Vista_CPU")
## Vista Grace ARM with Kokkos OpenMP backend (debug oracle for the GPU build).
## Same source code path as Vista (Kokkos kernels), just dispatched to OpenMP
## CPU threads instead of CUDA. GIZMO_GPU_COMPILER is NOT defined (no nvcc),
## so the per-TU AllDeviceMirror blocks are skipped and All remains the normal extern.
## Load modules: nvidia openmpi hdf5/2.0.0 fftw3 gsl kokkos/4.5.01-omp
## (i.e. kokkos/4.5.01-omp NOT kokkos/4.5.01-cuda; the rest matches the CUDA build).
## mpicxx on Vista wraps nvc++ (NVIDIA HPC SDK), not gcc — no -foffload needed.
CC       =  mpicc
CXX      =  mpicxx -std=c++17
FC       =  mpif90
OPTIMIZE = -O2 -Wall
OPTIMIZE += -fopenmp
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I$(TACC_SUNDIALS_INC)
CHIMESLIBS = -L$(TACC_SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial
endif
## Kokkos paths (TACC module sets TACC_KOKKOS_INC/LIB)
KOKKOS_INCL     = -I$(TACC_KOKKOS_INC)
KOKKOS_CPPFLAGS = -I$(TACC_KOKKOS_INC)
KOKKOS_LIBS_PATH= -L$(TACC_KOKKOS_LIB) -Wl,-rpath,$(TACC_KOKKOS_LIB)
KOKKOS_LIBS     = -lkokkoscore -lkokkoscontainers
## GPU TU files: same compiler (mpicxx → gcc), just add Kokkos includes
GPU_CXX    = mpicxx -std=c++17
GPU_CFLAGS = $(KOKKOS_INCL)
GPU_LDFLAGS= $(KOKKOS_LIBS_PATH)
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Lonestar6")
## Lonestar6 (TACC), NVIDIA A100 GPU nodes (sm_80; 3 per node on 2x AMD EPYC 7763 "Milan"), Kokkos
## CUDA backend. Device TUs (GPU_OBJS) are compiled by nvcc_wrapper with gcc as host compiler,
## reached through mpicxx so the MPI headers come with it; host TUs and the link use mpicxx.
## Load modules: gcc/13.2.0 impi/21.12 cuda/12.8 hdf5/1.14.6 fftw3/3.3.10
##               (+ hypre/2.31.0 for MHD_MODIFIED_GRADIENT, sundials/7.3.0 for CHIMES)
## TACC's only CUDA Kokkos module (kokkos/4.1.00-cuda) predates sort_by_key with a comparator
## (Kokkos 4.3), so each user builds Kokkos 4.6.01 for this machine once, into $WORK (login node,
## ~10 minutes); LS6_KOKKOS_CUDA_DIR below is its install prefix:
##   module purge; module load TACC gcc/13.2.0 cuda/12.8 cmake
##   KV=4.6.01; mkdir -p $WORK/src && cd $WORK/src
##   git clone --branch $KV --depth 1 https://github.com/kokkos/kokkos.git kokkos-$KV
##   cmake -S kokkos-$KV -B kokkos-$KV/build-cuda-a100 -DCMAKE_BUILD_TYPE=Release \
##     -DCMAKE_INSTALL_PREFIX=$WORK/apps/kokkos/$KV-cuda-a100 \
##     -DCMAKE_CXX_COMPILER=$WORK/src/kokkos-$KV/bin/nvcc_wrapper -DCMAKE_CXX_STANDARD=17 \
##     -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON
##   cmake --build kokkos-$KV/build-cuda-a100 -j 8 && cmake --install kokkos-$KV/build-cuda-a100
## Login nodes have no GPU driver, so libcuda is linked from the toolkit's stub directory, which is
## deliberately not on the run-time path: the compute nodes supply the real driver.
## The H100 nodes (sm_90) need a Kokkos built for HOPPER90 and are not covered by this block.
## submit to gpu-a100, gpu-a100-dev or gpu-a100-small
CC       =  mpicc
CXX      =  mpicxx -std=c++17
FC       =  mpif90
OPTIMIZE = -O2 -g -march=znver3 -Wall
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
LS6_KOKKOS_CUDA_DIR ?= $(WORK)/apps/kokkos/4.6.01-cuda-a100
KOKKOS_PATH     = $(LS6_KOKKOS_CUDA_DIR)
KOKKOS_CPPFLAGS = -I$(LS6_KOKKOS_CUDA_DIR)/include
KOKKOS_CXXFLAGS = --expt-relaxed-constexpr --expt-extended-lambda -arch=sm_80
KOKKOS_LDFLAGS  = -L$(LS6_KOKKOS_CUDA_DIR)/lib64 -L$(TACC_CUDA_LIB) -Wl,-rpath,$(TACC_CUDA_LIB) -L$(TACC_CUDA_LIB)/stubs
KOKKOS_LIBS     = -lkokkoscore -lkokkoscontainers -lcudart -lcuda
export NVCC_WRAPPER_DEFAULT_COMPILER = mpicxx
GPU_CXX     = $(LS6_KOKKOS_CUDA_DIR)/bin/nvcc_wrapper --std=c++17
GPU_CFLAGS  = $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS)
GPU_LDFLAGS = $(KOKKOS_LDFLAGS)
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I$(TACC_SUNDIALS_INC)
CHIMESLIBS = -L$(TACC_SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial -lsundials_core
endif
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
HYPRE_INCL = -I$(TACC_HYPRE_INC)
HYPRE_LIBS = -L$(TACC_HYPRE_LIB) -lHYPRE
endif
MKL_INCL =
MKL_LIBS =
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Lonestar6_CPU")
## Lonestar6 with the Kokkos OpenMP backend: the same Kokkos source path as "Lonestar6", dispatched
## to CPU threads. Host compiler, MPI, libraries and flags are identical to "Lonestar6" and the Kokkos
## version is the same (4.6.01), so the pair differs only in the Kokkos backend (as Frontier and
## Frontier_CPU do). Runs on the CPU nodes (normal, development) or on the GPU nodes as a CPU control.
## No nvcc, so GIZMO_GPU_COMPILER is undefined and the per-TU AllDeviceMirror blocks are skipped.
## Load modules: gcc/13.2.0 impi/21.12 kokkos/4.6.01-omp hdf5/1.14.6 fftw3/3.3.10
##               (+ hypre/2.31.0 for MHD_MODIFIED_GRADIENT, sundials/7.3.0 for CHIMES)
## i.e. the "Lonestar6" set with kokkos/4.6.01-omp in place of cuda and the self-built CUDA Kokkos.
CC       =  mpicc
CXX      =  mpicxx -std=c++17
FC       =  mpif90
OPTIMIZE = -O2 -g -march=znver3 -Wall
## the Kokkos OpenMP backend needs OpenMP whatever the Config says
OPTIMIZE += -fopenmp
KOKKOS_PATH     = $(TACC_KOKKOS_DIR)
KOKKOS_CPPFLAGS = -I$(TACC_KOKKOS_INC)
KOKKOS_LDFLAGS  = -L$(TACC_KOKKOS_LIB) -Wl,-rpath,$(TACC_KOKKOS_LIB)
KOKKOS_LIBS     = -lkokkoscore -lkokkoscontainers
GPU_CXX     = mpicxx -std=c++17
GPU_CFLAGS  = $(KOKKOS_CPPFLAGS)
GPU_LDFLAGS = $(KOKKOS_LDFLAGS)
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I$(TACC_SUNDIALS_INC)
CHIMESLIBS = -L$(TACC_SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial -lsundials_core
endif
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
HYPRE_INCL = -I$(TACC_HYPRE_INC)
HYPRE_LIBS = -L$(TACC_HYPRE_LIB) -lHYPRE
endif
MKL_INCL =
MKL_LIBS =
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Frontera_Kokkos")
## Frontera (CPU-only, Cascade Lake) with Kokkos OpenMP backend.
## Same Kokkos source path as Vista_CPU (CPU dispatch), Intel compiler + TACC libs.
## GIZMO_GPU_COMPILER is NOT defined (no nvcc) -> device-mirror blocks skipped,
## GPU-ported TUs build with mpicxx + Kokkos includes, dispatch to OpenMP threads.
## Load modules: intel/23.1.0 impi/21.9.0 kokkos/4.6.01-omp fftw3/3.3.10 hdf5/1.14.6 hypre/2.31.0
## (kokkos/4.6.01-omp sets TACC_KOKKOS_DIR/INC/LIB/BIN)
CC       =  mpicc
CXX      =  mpicxx -std=c++17
OPTIMIZE = -ggdb -O2 -xCORE-AVX2 -Wno-unknown-pragmas -Wall -Wno-format-security
## The Intel compilers default to a FAST floating-point model, unlike gcc and clang.
## That default alone breaks the position->key contract in declarations/typedefs.h, which is
## why no -ffast-math appears above and the build still refused to compile: the contract needs
## reciprocal substitution AND fused multiply-add both held off, and the two are separate
## switches.  Correctness before speed here -- a key that disagrees between the domain
## exchange and the tree build detaches a particle from every rank's tree, silently.
OPTIMIZE += -fp-model=precise -no-fma
## OpenMP required for the Kokkos OpenMP backend (unconditional, not gated on OPENMP config)
OPTIMIZE += -qopenmp
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I$(TACC_SUNDIALS_INC)
CHIMESLIBS = -L$(TACC_SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial
endif
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
## hypre module sets TACC_HYPRE_INC + TACC_HYPRE_LIB
ifdef TACC_HYPRE_INC
HYPRE_INCL = -I$(TACC_HYPRE_INC)
HYPRE_LIBS = -L$(TACC_HYPRE_LIB) -lHYPRE
else ifdef TACC_HYPRE_DIR
HYPRE_INCL = -I$(TACC_HYPRE_DIR)/include
HYPRE_LIBS = -L$(TACC_HYPRE_DIR)/lib -lHYPRE
endif
endif
## Kokkos OpenMP backend (TACC module sets TACC_KOKKOS_INC/LIB)
KOKKOS_INCL     = -I$(TACC_KOKKOS_INC)
KOKKOS_CPPFLAGS = -I$(TACC_KOKKOS_INC)
KOKKOS_LIBS_PATH= -L$(TACC_KOKKOS_LIB) -Wl,-rpath,$(TACC_KOKKOS_LIB)
KOKKOS_LIBS     = -lkokkoscore -lkokkoscontainers
## GPU-ported TUs: same compiler (mpicxx), just add Kokkos includes
GPU_CXX    = mpicxx -std=c++17
GPU_CFLAGS = $(KOKKOS_INCL)
GPU_LDFLAGS= $(KOKKOS_LIBS_PATH)
MKL_INCL = #
MKL_LIBS = #
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"CaltechHPC")
CC       =  mpicc
CXX      =  mpic++
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
FFTW_INCL= -I$(CPATH)
FFTW_LIBS= -L$(LIBRARY_PATH)
HDF5INCL = -I$(CPATH) -DH5_USE_16_API
HDF5LIB  = -L$(LIBRARY_PATH) -lhdf5 -lz
MPICHLIB = #
OPT     += -DHDF5_DISABLE_VERSION_CHECK
# Compiles with following modules:   1) intel/20.1    2) hdf5/1.10.1   3) fftw/3.3.7
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"BigRed200")
CC       =  cc # For Cray use this instead of mpicc
CXX      =  CC # mpic++
OPTIMIZE =  -O2

# Extra compile time warning flags
# OPTIMIZE += -Wall -Wextra -Wuninitialized -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -Wno-unused-variable -Wno-unused-but-set-variable

ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
endif
# All the library paths are already setup
MKL_INCL =
MKL_LIBS = -mkl=sequential
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
## module load cray-fftw cray-hdf5
## srun --cpus-per-task=$SLURM_CPUS_PER_TASK --ntasks-per-node=$SLURM_NTASKS_PER_NODE --nodes=$SLURM_JOB_NUM_NODES ./GIZMO ./gizmo_parameters.txt
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Expanse")
CC       = mpicc
CXX      = mpicxx
OPTIMIZE = -O3   # not -Ofast: that implies -ffast-math, see the position->key contract in declarations/typedefs.h
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
endif
FFTW_INCL= -I$(FFTWHOME)/include
FFTW_LIBS= -L$(FFTWHOME)/lib
HDF5INCL = -I$(HDF5HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += #
## modules to load
## module load slurm intel openmpi_ib fftw/2.1.5 hdf5
## run job with
## mpirun -v -x LD_LIBRARY_PATH ./GIZMO params.txt
endif

## MacBookCellar with Kokkos (OpenMP backend) — for testing GPU code paths without a GPU.
## Kokkos::parallel_for dispatches to OpenMP threads, SharedSpace = HostSpace.
## GPU TU files (cooling.cc, density_gpu.cc, etc.) are compiled with the same mpicxx
## compiler but with Kokkos include flags. GIZMO_GPU_COMPILER is NOT defined (no nvcc),
## so the per-TU AllDeviceMirror blocks are skipped and All remains the normal extern global.
## Install: brew install kokkos (needs libomp: brew install libomp)
ifeq ($(SYSTYPE),"MacBookCellar_Kokkos")
CC       =  mpicc
CXX      =  mpicxx -std=c++20
OPTIMIZE = -O3 -funroll-loops -march=native   # no -ffast-math: see the position->key contract in declarations/typedefs.h
OPTIMIZE += -Wno-unused-command-line-argument
## OpenMP is required for Kokkos OpenMP backend
OPTIMIZE += -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
OPTIMIZE += -L/opt/homebrew/opt/libomp/lib -lomp
## Kokkos paths (homebrew)
KOKKOS_INCL = -I/opt/homebrew/opt/kokkos/include
KOKKOS_CPPFLAGS = -I/opt/homebrew/opt/kokkos/include
KOKKOS_LIBS_PATH = -L/opt/homebrew/opt/kokkos/lib
KOKKOS_LIBS = -lkokkoscore -lkokkoscontainers
## GPU TU files: same compiler, just add Kokkos includes
GPU_CXX    = mpicxx -std=c++20
GPU_CFLAGS = $(KOKKOS_INCL)
GPU_LDFLAGS = $(KOKKOS_LIBS_PATH)
ifeq (MHD_MODIFIED_GRADIENT,$(findstring MHD_MODIFIED_GRADIENT,$(CONFIGVARS)))
ifneq (MHD_MODIFIED_GRADIENT_CG_ONLY,$(findstring MHD_MODIFIED_GRADIENT_CG_ONLY,$(CONFIGVARS)))
HYPRE_VERSION := $(shell ls /opt/homebrew/Cellar/hypre/ 2>/dev/null | sort -V | tail -n 1)
HYPRE_INCL = -I/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/include/
HYPRE_LIBS = -L/opt/homebrew/Cellar/hypre/$(HYPRE_VERSION)/lib/ -lHYPRE
endif
endif
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
# NOTE: set these to a SUNDIALS install in the CHIMES-supported 2.7-5.x range (see the
# note at the CHIMESINCL/CHIMESLIBS defaults above). Homebrew currently ships SUNDIALS
# 7.x (too new: no 'realtype', context-mandatory CVode API), so `brew`'s default path
# will not compile CHIMES; the paths below are the legacy Intel-Homebrew location and
# must be pointed at a supported SUNDIALS to build CHIMES on this platform.
CHIMESINCL = -I/usr/local/include/sundials
CHIMESLIBS = -L/usr/local/lib -lsundials_cvode -lsundials_nvecserial
endif
MKL_INCL = #
MKL_LIBS = #
FFTW_INCL= -I/opt/homebrew/Cellar/fftw/3.3.10_3/include
FFTW_LIBS= -L/opt/homebrew/Cellar/fftw/3.3.10_3/lib
HDF5_VERSION := $(shell ls /opt/homebrew/Cellar/hdf5/ 2>/dev/null | sort -V | tail -n 1)
HDF5INCL = -I/opt/homebrew/Cellar/hdf5/$(HDF5_VERSION)/include -DH5_USE_16_API
HDF5LIB  = -L/opt/homebrew/Cellar/hdf5/$(HDF5_VERSION)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += -DCHIMES_USE_DOUBLE_PRECISION -DGIZMO_GPU_ARENA_DEBUG
endif

#----------------------------
ifeq ($(SYSTYPE),"github-ubuntu")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE = -g -fcommon -O1 -funroll-loops -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -fipa-cp-clone  ## optimizations for gcc compilers (1/2)
OPTIMIZE += -ftree-loop-distribute-patterns -fvect-cost-model -ftree-partial-pre   ## optimizations for gcc compilers (2/2)
OPTIMIZE += -g -Wall # compiler warnings
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CXX     = mpic++
CHIMESINCL = -I/usr/include/sundials
CHIMESLIBS = -L/usr/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp # openmp required compiler flags
endif
MKL_INCL = #
MKL_LIBS = #
FFTW_INCL= -I/usr/include
FFTW_LIBS= -L/usr/lib
HDF5INCL = -I/usr/include/hdf5/openmpi -DH5_USE_16_API
HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ -lhdf5 -lz
MPICHLIB = #
OPT     += -DCHIMES_USE_DOUBLE_PRECISION
## to get required packages: sudo apt install libhdf5-openmpi-dev libopenmpi-dev
endif

#----------------------------
# Should work on any Flatiron institute linux cluster environment: rusty, popeye and linux workstations
ifeq ($(SYSTYPE),"RUSTY")
CC       =   mpicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC       =   mpicxx
endif
OPTIMIZE =  -O3 -funroll-loops -march=native -g -Wall   # no -ffast-math: see the position->key contract in declarations/typedefs.h
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
FFTW3_BASE= /mnt/sw/nix/store/bjzkf3pwcw0gy54db19kd4rl0xdiq98s-fftw-3.3.10/.
FFTW_INCL= -I$(FFTW3_BASE)/include
FFTW_LIBS= -L$(FFTW3_BASE)/lib -Xlinker -R -Xlinker $(FFTW3_BASE)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_BASE)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_BASE)/lib -Xlinker -R -Xlinker $(HDF5_BASE)/lib -lhdf5 -lz
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
			core/predict.o declarations/global.o core/begrun.o core/run.o core/gizmo_fatal.o declarations/allvars.o \
			declarations/lifecycle_counters.o \
			file_io/read_ic.o file_io/hdf5_deflate_filter.o domain/domain.o core/driftfac.o core/kicks.o \
			mesh/ghost_exchange.o mesh/ghost_writeback.o mesh/neighbor_list.o mesh/sfc_tiles.o mesh/gpu_dirty_tracker.o mesh/mode_b_local_walker.o mesh/mode_b_p2p_transport.o compile_time_info.o mesh/merge_split.o \
			core/transport_subcycle.o

SYSTEM_OBJS =   system/system.o \
				system/allocate.o \
				system/mymalloc.o \
				system/parallel_sort.o \
                system/peano.o \
                system/parallel_sort_special.o \
                system/mpi_util.o \
                system/memory_ledger.o \
                system/pinning.o

GRAVITY_OBJS  = gravity/forcetree.o \
                gravity/forcetree_update.o \
                gravity/gravtree.o \
				gravity/cosmology.o \
				gravity/pm_periodic.o \
                gravity/pm_nonperiodic.o \
                gravity/longrange.o \
                gravity/ags_rkern.o \
                gravity/star_direct_gravity.o \
                gravity/let_pack.o

HYDRO_OBJS = 	hydro/hydro_toplevel.o \
				hydro/density.o \
				hydro/gradients.o \
				hydro/hydro_corridor.o \
				hydro/mg_gradient_correction.o \
				turb/dynamic_diffusion.o \
				turb/dynamic_diffusion_velocities.o \
				turb/turb_powerspectra.o
## turb/turb_driving.o is in GPU_OBJS (compiled by nvcc_wrapper for Kokkos offload)

## GPU_OBJS: files compiled by nvcc_wrapper (Kokkos CUDA backend) with GPU_CFLAGS.
## Must NOT also appear in OBJS/EOSCOOL_OBJS or the pattern rule will create duplicate symbols.
## eos/eos.o is here because it contains yhelium/Get_Gas_Mean_Molecular_Weight_mu/
## Get_Gas_Molecular_Mass_Fraction which are called from device cooling functions.
GPU_OBJS = core/gpu_drift.o cooling/cooling.o eos/eos.o hydro/density_gpu.o hydro/density_loop.o hydro/cellcorrections_loop.o hydro/gradients_loop.o hydro/hydro_force_loop.o mesh/gpu_neighbor_list.o mesh/neighbor_loop_runner.o radiation/rt_chem.o turb/turb_driving.o turb/difffilter_loop.o solids/grain_drag.o galaxy_sf/dm_dispersion.o galaxy_sf/dm_dispersion_loop.o gravity/ags_density_loop.o gravity/ags_force_loop.o sidm/dm_fuzzy_loop.o sidm/cbe_integrator_gpu.o sidm/cbe_integrator_gradients.o radiation/rt_source_injection_loop.o galaxy_sf/thermal_fb.o galaxy_sf/thermal_fb_loop.o sinks/sink_feed.o sinks/sink_feed_loop.o sinks/sink_environment.o sinks/sink_env1_loop.o sinks/sink_env2_loop.o galaxy_sf/mechfb_loop.o solids/grain_physics_loop.o sinks/sink_swallow_and_kick.o sinks/sink_swk_loop.o galaxy_sf/radfb_rp_loop.o system/gpu_particles_arena.o system/kokkos_mem_telemetry.o gravity/gpu_gravity_tree.o gravity/gpu_gravtree.o gravity/gpu_moment_refresh.o gravity/gpu_nextnode_thread.o gravity/gpu_morton.o gravity/gpu_peano_walk.o gravity/gpu_topology_build.o gravity/gpu_topology_finalize.o gravity/gpu_pseudo_update.o gravity/gpu_force_drift.o gravity/gpu_force_update.o
## Nuclear network files are added to GPU_OBJS below (conditional on NUCLEAR_NETWORK)
EOSCOOL_OBJS =  \
				cooling/grackle.o \
				cooling/simple_chemistry.o \
				eos/hydrogen_molecule.o \
				eos/cosmic_ray_fluid/cosmic_ray_utilities.o \
				solids/elastic_physics.o \
				solids/grain_physics.o \
				solids/grain_promotion.o \
				solids/ism_dust_chemistry.o

STARFORM_OBJS = galaxy_sf/sfr_eff.o \
                galaxy_sf/stellar_evolution.o \
                galaxy_sf/mechanical_fb.o \
                galaxy_sf/radfb_local.o \
                galaxy_sf/dm_heating.o
# galaxy_sf/thermal_fb.o is in GPU_OBJS (3e.2): the runner-template caller
# transitively includes thermal_fb_loop.h → Kokkos_Core.hpp, so the TU must
# build with $(GPU_CC) (nvcc_wrapper on Vista) like sinks/sink_feed.o.
# galaxy_sf/radfb_local.o stays in STARFORM_OBJS (host-only): radfb_local.cc
# contains only HII_heating_singledomain + do_the_local_ionization (legacy
# host code with many bare All.* reads). The RP runner-template caller
# `radiation_pressure_winds_consolidated` lives in radfb_rp_loop.cc (a GPU
# TU). The macro-poisoning rationale for the split is gone under the
# per-TU AllDeviceMirror scheme (device-pass-only macro), but keeping
# radfb_local.cc out of GPU_OBJS still saves compile time + avoids
# pulling Kokkos headers into the unrelated host code.
# galaxy_sf/dm_dispersion.o and galaxy_sf/dm_dispersion_loop.o are in GPU_OBJS
# (3d.D): dm_dispersion.cc includes dm_dispersion_loop.h → Kokkos_Core.hpp.

SINK_OBJS = sinks/sink.o \
            sinks/sink_util.o
## sinks/sink_swallow_and_kick.o moved to GPU_OBJS (3d.3) — the .cc now
## includes sink_swk_loop.h which pulls in Kokkos_Core.hpp via the inline
## pair body's atomics, so CUDA backend requires nvcc_wrapper compilation.
## sinks/sink_swk_loop.o is in GPU_OBJS for the same reason.
## sinks/sink_env2_loop.o is in GPU_OBJS (3d.2): the .cc uses
## Kokkos::kokkos_malloc / kokkos_free in populate/cleanup_device_context
## (per-active Jgas/Jstar UVM staging), so the CUDA backend requires
## nvcc_wrapper compilation. The header is Kokkos-free (pure i-side
## pair body) so consumers (sink_environment.cc) stay host-pattern.
## sinks/sink_feed.o and sinks/sink_feed_loop.o are in GPU_OBJS (3d.1) —
## both pull in Kokkos_Core.hpp via sink_feed_loop.h's atomics, so the
## CUDA backend requires nvcc_wrapper compilation. On Mac OpenMP backend
## GPU_CC falls back to CXX so the same rule still works.
## sinks/sink_environment.o and sinks/sink_env1_loop.o moved to GPU_OBJS
## (descope-fix audit, 2026-05-21) — sink_env1_loop.h now includes
## Kokkos_Core.hpp unconditionally (Plan #1's Kokkos::atomic_min on
## P[j].SwallowTime under SINGLE_STAR_SINK_DYNAMICS+SINK_GRAVCAPTURE_GAS).
## sink_environment.cc transitively pulls in the header, so it also requires
## nvcc_wrapper to avoid the catastrophic KOKKOS_ENABLE_CUDA / __CUDACC__
## diagnostic that mpicxx triggers on Kokkos_Setup_Cuda.hpp.

RHD_OBJS =  radiation/rt_utilities.o \
			radiation/rt_source_injection.o \
			radiation/rt_dust_opacity.o
## radiation/rt_chem.o is in GPU_OBJS (compiled by nvcc_wrapper for Kokkos offload)

FOF_OBJS =	structure/fof.o \
			structure/group_search.o \
			structure/subfind/subfind.o \
			structure/subfind/subfind_vars.o \
			structure/subfind/subfind_collective.o \
			structure/subfind/subfind_serial.o \
			structure/subfind/subfind_so.o \
			structure/subfind/subfind_cont.o \
			structure/subfind/subfind_distribute.o \
			structure/subfind/subfind_findlinkngb.o \
			structure/subfind/subfind_nearesttwo.o \
			structure/subfind/subfind_modern_search.o \
			structure/subfind/subfind_loctree.o \
			structure/subfind/subfind_potential.o \
			structure/subfind/subfind_density.o \
			structure/twopoint.o

MISC_OBJS = sidm/cbe_integrator.o \
			sidm/cbe_collisions.o \
			sidm/dm_fuzzy.o \
			sidm/sidm_core.o

## name of executable and optimizations
EXEC   = GIZMO
OPTIONS = $(OPTIMIZE) $(OPT)

## combine all the objects above
OBJS  = $(CORE_OBJS) $(SYSTEM_OBJS) $(GRAVITY_OBJS) $(HYDRO_OBJS) \
		$(EOSCOOL_OBJS) $(STARFORM_OBJS) $(SINK_OBJS) $(RHD_OBJS) \
		$(FOF_OBJS) $(MISC_OBJS)
## GPU_OBJS are kept separate so the pattern rule does not override their specific compile rules


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


## now we add special cases dependent on compiler flags. normally we would
##  include the files always, and simply use the in-file compiler variables
##  to determine whether certain code is compiled [this allows us to take
##  advantage of compiler logic, and makes it easier for the user to
##  always specify what they want]. However special cases can arise, if e.g.
##  there are certain special libraries needed, or external compilers, for
##  certain features

# helmholtz eos — now pure C++ (ported from Fortran), no Fortran compiler needed
ifeq (EOS_HELMHOLTZ,$(findstring EOS_HELMHOLTZ,$(CONFIGVARS)))
OBJS    += eos/eos_interface.o eos/helmholtz/helmholtz.o
INCL    += eos/helmholtz/helmholtz.h
endif

ifeq (EOS_ANEOS,$(findstring EOS_ANEOS,$(CONFIGVARS)))
OBJS    += eos/aneos.o
INCL    += eos/aneos.h
endif

# nuclear reaction network — nuclear.o is GPU-compiled; physics functions are header-only
# (nuclear_physics_functions.h, included by nuclear.cc) so no separate physics .o is needed.
ifeq (NUCLEAR_NETWORK,$(findstring NUCLEAR_NETWORK,$(CONFIGVARS)))
GPU_OBJS += nuclear/nuclear.o
INCL     += nuclear/nuclear.h nuclear/nuclear_physics_functions.h
endif
ifeq (NUCLEAR_NETWORK_SOLVER=1,$(findstring NUCLEAR_NETWORK_SOLVER=1,$(CONFIGVARS)))
OBJS    += nuclear/nuclear_skynet.o
SKYNETLIBS = -lskynet
else
SKYNETLIBS =
endif
ifeq (NUCLEAR_NETWORK_SOLVER=2,$(findstring NUCLEAR_NETWORK_SOLVER=2,$(CONFIGVARS)))
OBJS    += nuclear/nuclear_torch.o
TORCHLIBS = -ltorch_nuclear
else
TORCHLIBS =
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


# linking libraries (includes machine-dependent options above)
CFLAGS = $(OPTIONS) $(FFTW_INCL) $(HDF5INCL) \
         $(GRACKLEINCL) $(CHIMESINCL) $(HYPRE_INCL) $(MKL_INCL) $(KOKKOS_CPPFLAGS)



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


LIBS = $(HDF5LIB) -g $(MPICHLIB) \
	   $(FFTW_LIBS) $(FFTW_LIBNAMES) -lm $(GRACKLELIBS) $(CHIMESLIBS) $(SKYNETLIBS) $(TORCHLIBS) $(HYPRE_LIBS) $(MKL_LIBS)


$(EXEC): $(OBJS) $(GPU_OBJS)
	$(CXX) $(OPTIMIZE) $(GPU_LDFLAGS) $(OBJS) $(GPU_OBJS) $(KOKKOS_LIBS) $(LIBS) -o $(EXEC)

## GPU-offloaded files: compiled via nvcc_wrapper so nvcc handles device code.
## GPU_CXX = $(KOKKOS_PATH)/bin/nvcc_wrapper on GPU systypes; falls back to $(CXX) otherwise.
## GPU_CFLAGS carries the Kokkos include/feature flags populated by Makefile.kokkos.
GPU_CC = $(if $(GPU_CXX),$(GPU_CXX),$(CXX))
cooling/cooling.o: cooling/cooling.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
eos/eos.o: eos/eos.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
nuclear/nuclear.o: nuclear/nuclear.cc nuclear/nuclear_physics_functions.h $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
hydro/density_gpu.o: hydro/density_gpu.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
mesh/gpu_neighbor_list.o: mesh/gpu_neighbor_list.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
mesh/neighbor_loop_runner.o: mesh/neighbor_loop_runner.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
radiation/rt_chem.o: radiation/rt_chem.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
turb/turb_driving.o: turb/turb_driving.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
turb/difffilter_loop.o: turb/difffilter_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
solids/grain_drag.o: solids/grain_drag.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## galaxy_sf/dm_dispersion_gpu.o retired by DMDispersionSpec.
galaxy_sf/dm_dispersion.o: galaxy_sf/dm_dispersion.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
galaxy_sf/dm_dispersion_loop.o: galaxy_sf/dm_dispersion_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/ags_density_loop.o: gravity/ags_density_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
hydro/density_loop.o: hydro/density_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
hydro/cellcorrections_loop.o: hydro/cellcorrections_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
hydro/gradients_loop.o: hydro/gradients_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
hydro/hydro_force_loop.o: hydro/hydro_force_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/ags_force_loop.o: gravity/ags_force_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## sidm/dm_fuzzy_gpu.o retired by DMGradSpec (dm_fuzzy port). dm_fuzzy.o stays
## host-only (do_dm_fuzzy_drift_kick / _initialization); DMGrad loop lives in
## dm_fuzzy_loop.o.
sidm/dm_fuzzy_loop.o: sidm/dm_fuzzy_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## sinks/sink_environment_gpu.o retired; legacy
## sink_environment_second_evaluate_gpu absorbed by SinkEnv2Spec runner port.
sidm/cbe_integrator_gpu.o: sidm/cbe_integrator_gpu.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
sidm/cbe_integrator_gradients.o: sidm/cbe_integrator_gradients.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
radiation/rt_source_injection_loop.o: radiation/rt_source_injection_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
galaxy_sf/thermal_fb.o: galaxy_sf/thermal_fb.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
galaxy_sf/thermal_fb_loop.o: galaxy_sf/thermal_fb_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## sinks/sink_feed_gpu.o retired; legacy `sink_feed_gpu.cc` out of build,
## file deletion deferred to cleanup commit.
sinks/sink_feed.o: sinks/sink_feed.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
sinks/sink_feed_loop.o: sinks/sink_feed_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
sinks/sink_env2_loop.o: sinks/sink_env2_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## sinks/sink_environment.o + sink_env1_loop.o promoted to nvcc_wrapper after
## the descope-fix audit: sink_env1_loop.h now includes
## Kokkos_Core.hpp unconditionally for Kokkos::atomic_min on P[j].SwallowTime
## under SINGLE_STAR_SINK_DYNAMICS + SINK_GRAVCAPTURE_GAS. mpicxx compilation
## would trip the catastrophic KOKKOS_ENABLE_CUDA / __CUDACC__ diagnostic.
sinks/sink_environment.o: sinks/sink_environment.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
sinks/sink_env1_loop.o: sinks/sink_env1_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
galaxy_sf/mechfb_loop.o: galaxy_sf/mechfb_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## grain_physics_loop.o: runner-template Spec host hooks + toplevels for
## GrainBackrxSpec + GrainRTGasSpec + GrainRTGrainSpec. Replaced the bespoke
## GPU evaluators formerly in grain_physics_gpu.cc (retired).
solids/grain_physics_loop.o: solids/grain_physics_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## sinks/sink_swallow_and_kick_gpu.cc + sinks/sink_swallow_and_kick_functions.h
## deleted (port + cleanup folded into one commit per SSOT).
## New build rules for 3d.3:
sinks/sink_swk_loop.o: sinks/sink_swk_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
sinks/sink_swallow_and_kick.o: sinks/sink_swallow_and_kick.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## radfb_local_gpu.o retired in radfb_local cleanup commit
## (legacy radiation_pressure_winds_gpu replaced by runner-template port).
galaxy_sf/radfb_rp_loop.o: galaxy_sf/radfb_rp_loop.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
system/gpu_particles_arena.o: system/gpu_particles_arena.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
system/kokkos_mem_telemetry.o: system/kokkos_mem_telemetry.cc $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_gravity_tree.o: gravity/gpu_gravity_tree.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_gravtree.o: gravity/gpu_gravtree.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_moment_refresh.o: gravity/gpu_moment_refresh.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_nextnode_thread.o: gravity/gpu_nextnode_thread.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_morton.o: gravity/gpu_morton.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_peano_walk.o: gravity/gpu_peano_walk.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_topology_build.o: gravity/gpu_topology_build.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_topology_finalize.o: gravity/gpu_topology_finalize.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_pseudo_update.o: gravity/gpu_pseudo_update.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_force_drift.o: gravity/gpu_force_drift.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
core/gpu_drift.o: core/gpu_drift.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
gravity/gpu_force_update.o: gravity/gpu_force_update.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@
## gravity/let_pack.o is CPU-only (no Kokkos): the $(OBJS) static-pattern rule
## builds it with $(CXX). It must NOT have its own $(GPU_CC) rule (a duplicate
## explicit recipe make would warn about and override anyway).
declarations/allvars_gpu.o: declarations/allvars_gpu.cu declarations/global_data_all_struct.h $(CONFIG) compile_time_info.cc
	$(GPU_CC) $(CFLAGS) $(GPU_CFLAGS) -c $< -o $@

$(OBJS): %.o: %.cc $(INCL) $(CONFIG) compile_time_info.cc
	$(CXX) $(CFLAGS) -c $< -o $@

# Vista NVC++ 25.5 ICE workaround: gen_llvm_expr() unknown opcode at -O2 inside
# addFB_evaluate in mechanical_fb.cc. Compile this one CPU-only object at -O1 on
# Vista via a target-specific CFLAGS override (NOT a second explicit recipe, which
# collides with the $(OBJS) static-pattern rule and warns). The static-pattern
# recipe then builds it with the -O1 CFLAGS. GPU path lives in mechfb_loop.cc;
# negligible performance impact.
ifeq ($(SYSTYPE),"Vista")
galaxy_sf/mechanical_fb.o: CFLAGS := $(subst -O2,-O1,$(CFLAGS))
endif

compile_time_info.cc: $(CONFIG)
	$(PERL) file_io/prepare-config.perl $(CONFIG)

clean:
	find . -name '*.o' -type f -delete
	rm -f $(EXEC) *.oo *.c~ compile_time_info.cc GIZMO_config.h
## NOTE: nuke ALL object files generically (not just $(OBJS)/$(GPU_OBJS)/$(FOBJS)
## for the active Config). The OBJS lists depend on the current Config.sh, so a
## Config-scoped clean leaves stale .o behind when switching between Configs
## (the Makefile has no header-dependency tracking) -- those then silently link
## into a mixed/broken binary. A full .o sweep is always the intended policy.
