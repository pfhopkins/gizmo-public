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
#   slighty by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly 
#   dealing with new files and filename conventions)
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


# initialize some default flags -- these will all get re-written below
CC	= mpicc		# sets the C-compiler (default, will be set for machine below)
CXX	= mpiCC		# sets the C++-compiler (default, will be set for machine below)
FC	= mpif90	# sets the fortran compiler (default, will be set for machine below)
OPTIMIZE = -Wall  -g   # optimization and warning flags (default)
MPICHLIB = -lmpich	# mpi library (arbitrary default, set for machine below)
CHIMESINCL = # default to empty, will only be used below if called
CHIMESLIBS = # default to empty, will only be used below if called

# one annoying thing here is the FFTW libraries, since they are named differently depending on
#  whether they are compiled in different precision levels, or with different parallelization options, so we
#  have to have a big block here 'sorting them out'.
#
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
    FFTW_LIBNAMES =  #-lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
    FFTW_LIBNAMES =  #-ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
    FFTW_LIBNAMES =  #-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
# we only need fftw if PMGRID is turned on
ifneq (USE_FFTW3, $(findstring USE_FFTW3, $(CONFIGVARS)))
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
# or if TURB_DRIVING_SPECTRUMGRID is activated
ifeq (TURB_DRIVING_SPECTRUMGRID, $(findstring TURB_DRIVING_SPECTRUMGRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
  FFTW_LIBNAMES = #
endif
endif
else # use FFTW3 instead of FFTW2.?
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -lfftw3_mpi -lfftw3
else #single precision 
  FFTW_LIBNAMES = -lfftw3f_mpi -lfftw3f
endif
else 
# or if TURB_DRIVING_SPECTRUMGRID is activated
ifeq (TURB_DRIVING_SPECTRUMGRID, $(findstring TURB_DRIVING_SPECTRUMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -lfftw3_mpi -lfftw3
else #single precision  
  FFTW_LIBNAMES = -lfftw3f_mpi -lfftw3f
endif
else 
  FFTW_LIBNAMES = #
endif
endif
endif


## read the systype information to use the blocks below for different machines
ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif
FINCL =




#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Stampede")
CC       =  mpicc
CXX      =  mpic++
FC       =  mpif90 -nofor_main
OPTIMIZE = -O3 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -g -Wall # compiler warnings
#OPTIMIZE += -parallel -openmp  # openmp (comment out this line if OPENMP not used)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp  # openmp required compiler flags
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(TACC_MKL_INC)
MKL_LIBS = -L$(TACC_MKL_LIB) -mkl=sequential
##MKL_LIBS = -L$(TACC_MKL_LIB) -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
#MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
## modules to load: 
## module load intel mvapich2 gsl hdf5 fftw2
##  -- performance is very similar with impi (intel-mpi) instead of mpavich2, 
##   if preferred use that with MPICHLIB line uncommented
## newest version of code needed for compatibility with calls in MPI-2 libraries
##
endif



ifeq ($(SYSTYPE),"Stampede2")
CC       =  mpicc
CXX      =  mpic++
FC       =  mpif90 -nofor_main
OPTIMIZE = -O3 $(TACC_VEC_FLAGS) -ipo -funroll-loops -no-prec-div -fp-model fast=2
## above is preferred, $(TACC_VEC_FLAGS) automatically incorporates the TACC preferred flags for both KNL or SKX nodes
#OPTIMIZE = -O3 -xMIC-AVX512 -ipo -funroll-loops -no-prec-div -fp-model fast=2  # (deprecated, -xMIC-AVX512 is specific to the KNL nodes)
OPTIMIZE += -g -Wall # compiler warnings
#OPTIMIZE += -parallel -openmp  # openmp (comment out this line if OPENMP not used)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -qopenmp  # openmp required compiler flags
endif
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I${MY_SUNDIALS_PARENT_DIR}/sundials/instdir/include
CHIMESLIBS = -L${MY_SUNDIALS_PARENT_DIR}/sundials/instdir/lib64 -lsundials_cvode -lsundials_kinsol -lsundials_nvecserial
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(TACC_MKL_INC)
MKL_LIBS = -L$(TACC_MKL_LIB) -mkl=sequential
##MKL_LIBS = -L$(TACC_MKL_LIB) -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
ifeq (USE_FFTW3, $(findstring USE_FFTW3, $(CONFIGVARS)))
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
endif
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
##
## module load TACC intel impi hdf5 gsl fftw2
##  - note you can choose to use FFTW3 now instead of FFTW2, but you will need to load that module and change the compiler link appropriately
## note is you are using the KNL system it has a large number of slow cores, so some changes to 'usual' compilation parameters are advised:
##  - recommend running with ~16 mpi tasks/node. higher [32 or 64] usually involves a performance hit unless the problem is more scale-able;
##     use the remaining nodes in OPENMP. Do not use >64 MPI tasks/node [need ~4 cores free for management] and do not use >2 threads/core
##     [should never have >128 threads/node] -- the claimed 4 hardware threads/core includes non-FP threads which will severely slow performance.
##     so 'default' would be ~16 tasks/node, OMP_NUM_THREADS=8.
##  - because of the large core/thread count, MULTIPLEDOMAINS should be set low, MULTIPLEDOMAINS=1 ideally [already problem is heavily-divided].
##     - likewise be careful with domain decomposition, TreeDomainUpdateFrequency param [so don't spend very long running domain decompositions]
##  - memory is large per node: for 16 tasks/node, large MaxMemSize=5450 is reasonable, with BufferSize=450, and large PartAllocFactor=40 can be used
##  - run job with "tacc_affinity" on.
##
## for SUNDIALS above [CHIMES], you'll need to install your own version, but its easy to do: instructions from Alex Gurvitch which worked well for me:
##  - create parent folder: mkdir MY_SUNDIALS_PARENT_DIR/sundials
##  - copy sundials-5.1.0.tgz to folder MY_SUNDIALS_PARENT_DIR/sundials
##  - cd MY_SUNDIALS_PARENT_DIR/sundials
##  - tar -xvzf sundials-5.1.0.tgz
##  - mkdir instdir
##  - mkdir builddir
##  - cmake -DCMAKE_INSTALL_PREFIX=MY_SUNDIALS_PARENT_DIR/sundials/instdir sundials-5.1.0
##  - make
##  - make install
##  - finally make sure to define the directory MY_SUNDIALS_PARENT_DIR and export it in your .bashrc file, so its visible to the compilers and linkers
##
endif




ifeq ($(SYSTYPE),"Frontera")
CC       =  mpicc
CXX      =  mpic++
FC       =  mpif90 -nofor_main
OPTIMIZE = -O2 -xCORE-AVX2
#OPTIMIZE = -O3 $(TACC_VEC_FLAGS) -ipo -funroll-loops -no-prec-div -fp-model fast=2
#OPTIMIZE = -O3 -xCORE-AVX512 -ipo -funroll-loops -no-prec-div -fp-model fast=2
## above is preferred, $(TACC_VEC_FLAGS) automatically incorporates the TACC preferred flags for both KNL or SKX nodes, but gives tiny performance hit
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -qopenmp
endif
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I$(TACC_SUNDIALS_INC)
CHIMESLIBS = -L$(TACC_SUNDIALS_LIB) -lsundials_cvode -lsundials_nvecserial
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(TACC_MKL_INC)
MKL_LIBS = -L$(TACC_MKL_LIB) -mkl=sequential
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
ifeq (USE_FFTW3, $(findstring USE_FFTW3, $(CONFIGVARS)))
FFTW_INCL= -I$(TACC_FFTW3_INC)
FFTW_LIBS= -L$(TACC_FFTW3_LIB)
endif
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE -DNO_ISEND_IRECV_IN_DOMAIN -DHDF5_DISABLE_VERSION_CHECK
##
# UPDATE (9/19): Intel/19.0.5 is now working, and Intel/18 is actually sometimes running slower now because of some of the changes made to the impi installation.
#          Depending on when your code was compiled and exactly which flags you used, you may notice a performance drop with intel/18, and should switch to 19.
#          For intel/19: module load intel/19 impi hdf5 fftw3 gsl
#
# Previous: presently must use intel/18.x versions. 19.x versions compile and work, but lots of problems (+slower), esp. for high Ntasks or OpenMP
#  e.g.: module load intel/18.0.5 impi hdf5 fftw3 gsl
#  until recently, GSL module did -not- support intel/18.x, so needed to build it yourself (see update below). example instructions below:
#    -- 1. get newest GSL: ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
#       2. unpack, 3. then in folder run: "./configure --prefix=$HOME/gsl-2.5 CC=icc" followed by 4. "make" and 5. "make all"
#           (here I'm setting "$HOME/gsl-2.5" as the local install directory, you set yours appropriately)
#       6. in your .bashrc file, add "export HOME_GSL_DIR=$HOME/gsl-2.5" and
#           "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME_GSL_DIR:$HOME_GSL_DIR/.libs:$HOME_GSL_DIR/cblas/.libs"
#           (obviously if you use a different parent install directory, change the directory name here accordingly).
#       7. when you submit jobs, make sure you include a "source $HOME/.bashrc" in your run script or the export flags above, to link the libraries. I was using
#			GSL_INCL = -I$(HOME_GSL_DIR)
#			GSL_LIBS = -L$(HOME_GSL_DIR)/.libs -L$(HOME_GSL_DIR)/cblas/.libs
# [update: GSL module is now installed for intel/18.0.5, so you can simply load the module. but I'll keep the install instructions above, they can be useful]
#
# As usual include "umask 022" and "ulimit -s unlimited" in your .bashrc file to save headaches later
# fftw2/3 work equally well. usual intuition re: multipledomains, pmgrid, treedomainfreq, etc, apply.
# The different code optimizations above make very tiny differences. for stability I am for now using -O2 -xCORE-AVX2, nothing 'fancy' but this doesn't cost us
# Run scripts are simple SBATCH, like on Stampede and many other machines. Examples of several appear in this file. Example run script:
#                    #!/bin/bash
#                    #SBATCH -J (NAME) -p normal -N (NUMBER_OF_NODES) --ntasks-per-node (56/OPENMP_NUMBER) -t (RUNTIME_REQUEST) -A (ACCOUNT_NAME_TO_CHARGE)
#                    export OMP_NUM_THREADS=(OPENMP_NUMBER)
#                    source $HOME/.bashrc
#                    ibrun ./GIZMO ./params.txt (GIZMO_STARTUP_FLAG) 1>gizmo.out 2>gizmo.err
#     where quantities in (X) are the things you want to set.
# With these options, hybrid MPI+OpenMP works well. Because of the node configuration, optimal hybrid performance will typically use either
#   OPENMP=4 (ntasks-per-node=14) or OPENMP=7 (ntasks-per-node=8). Small jobs (<200 cores) might be better with smaller/no OPENMP, very large jobs higher,
#   (OPENMP can be any integer, ntasks-per-node must be even or severe performance hits apply).
#   Intel/19 now functional seems to favor slightly lower OPENMP number, shifting to perhaps OPENMP=2 (ntasks-per-node=28) for small jobs, =4 for medium, =7 for very large
#
# Note that the Frontera setup is NOT built for hyperthreading, even though the CLX nodes are supposed to support it. If you ask for 112 threads/node (insteady of 56),
#   the code will actually work, but very slowly. Stick to 56 for now.
#
# [old: There are still odd memory issues. The machine should have 3.3gb/core available after OS, etc, but in practice we need to allocate less than this. MPI errors
#   have also been appearing in large runs (for almost all users) related to memory. Be careful for now, and communicate to TACC support staff re: memory issues.]
#   I am using ~3gb/core for low task numbers, lower still for higher task numbers. 
##
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Bridges2")
CC       = mpicc
CXX      = mpic++
FC       = mpif90  # gcc/clang
#FC       = mpif90 -nofor_main  # intel
OPTIMIZE = -O3
OPTIMIZE += -march=znver1 -mfma -fvectorize -mfma -mavx2 -m3dnow -floop-unswitch-aggressive -fcommon  # aocc/clang
#OPTIMIZE += -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer -fcommon  # gcc
#OPTIMIZE += -march=core-avx2 -fma -ftz -fomit-frame-pointer -ipo -funroll-loops -no-prec-div -fp-model fast=2  # intel
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp  # gcc/clang
#OPTIMIZE += -qopenmp  # intel
endif
MKL_INCL = -I$(INCLUDE)
MKL_LIBS = -L$(LIBRARY_PATH) -mkl=sequential
GSL_INCL = -I$(INCLUDE)
GSL_LIBS = -L$(LIBRARY_PATH)
FFTW_INCL= -I$(INCLUDE)
FFTW_LIBS= -L$(LIBRARY_PATH)
HDF5INCL = -I$(INCLUDE) -DH5_USE_16_API
HDF5LIB  = -L$(LIBRARY_PATH) -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
## modules to load for aocc/clang:
## module load aocc
## module load openmpi/4.0.2-clang2.1
## modules to load for gcc:
## module load gcc
## module load openmpi/3.1.6-gcc8.3.1
## modules to load for intel:
## module load intel
## module load openmpi/4.0.2-intel20.4
## additional modules to load for any compiler:
## module load mkl
## module load hdf5
## module load fftw
# NOTE: this machine does not appear to always parallelize correctly in a hybrid MPI/OpenMP (i.e. when OPENMP is enabled) setup - it can potentially assign multiple threads to the same physical core and get lousy performance. A surefire way to get the correct thread affinity is to generate a rankfile and include as an argument to mpirun, e.g. https://github.com/mikegrudic/make_rankfile
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Anvil")
CC       = mpicc
CXX      = mpic++
FC       = mpif90  # gcc/clang
#FC       = mpif90 -nofor_main  # intel
OPTIMIZE = -O3
OPTIMIZE += -march=znver1 -mfma -fvectorize -mfma -mavx2 -m3dnow -floop-unswitch-aggressive -fcommon  # aocc/clang
#OPTIMIZE += -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer -fcommon  # gcc
#OPTIMIZE += -march=core-avx2 -fma -ftz -fomit-frame-pointer -ipo -funroll-loops -no-prec-div -fp-model fast=2  # intel
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp  # gcc/clang
#OPTIMIZE += -qopenmp  # intel
endif
MKL_INCL = -I$(C_INCLUDE_PATH)
MKL_LIBS = -L$(LIBRARY_PATH) -mkl=sequential
GSL_INCL = -I$(C_INCLUDE_PATH)
GSL_LIBS = -L$(LIBRARY_PATH)
FFTW_INCL= -I$(C_INCLUDE_PATH)
FFTW_LIBS= -L$(LIBRARY_PATH)
HDF5INCL = -I$(C_INCLUDE_PATH) -DH5_USE_16_API
HDF5LIB  = -L$(LIBRARY_PATH) -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
# modules to load: aocc openmpi hdf5 gsl fftw
# NOTE: this machine does not appear to always parallelize correctly in a hybrid MPI/OpenMP (i.e. when OPENMP is enabled) setup - it can potentially assign multiple threads to the same physical core and get lousy performance. A surefire way to get the correct thread affinity is to generate a rankfile and include as an argument to mpirun, e.g. https://github.com/mikegrudic/make_rankfile
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Pleiades")
CC       =  icc -lmpi
CXX      =  icc -lmpi -lmpi++
FC       =  ifort -nofor_main -lmpi
OPTIMIZE = -O3 -axCORE-AVX512,CORE-AVX2 -xAVX -ip -funroll-loops -no-prec-div -fp-model fast=2 -diag-disable 3180  # all core types
OPTIMIZE += -Wall  # compiler warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -qopenmp
endif
GMP_INCL =
GMP_LIBS =
GSL_INCL = -I$(PKGSRC_BASE)/include/gsl
GSL_LIBS = -L$(PKGSRC_BASE)/lib
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
ifeq (USE_FFTW3, $(findstring USE_FFTW3, $(CONFIGVARS)))
FFTW_INCL= -I$(FFTW3_HOME)/include
FFTW_LIBS= -L$(FFTW3_HOME)/lib
endif
HDF5INCL = -I$(PKGSRC_BASE)/include -DH5_USE_16_API
HDF5LIB  = -L$(PKGSRC_BASE)/lib -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
endif
## modules to load:
##     module load comp-intel mpi-hpe/mpt pkgconfig
## also add this to your .bashrc (after loading the modules above):
##     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PKGSRC_BASE/lib
## FFTW: Pleiades does not have an MPI-enabled FFTW, so install your own version in your home directory
##     for FFTW2 - to build:
##         ./configure --prefix=$HOME/local/fftw2 --enable-mpi --enable-type-prefix --enable-float
##     for FFTW3 - to build:
##         ./configure --prefix=$HOME/local/fftw3 --enable-mpi --enable-float
##     then define the following variable(s) in your .bashrc:
##         export FFTW2_HOME=$HOME/local/fftw2
##         and/or
##         export FFTW3_HOME=$HOME/local/fftw3
##     so Gizmo can find your FFTW installation
## in your job submission script, include the following before submitting your job:
##     export MPI_DSM_DISTRIBUTE=0
##     export KMP_AFFINITY=disabled


#----------------------------
ifeq ($(SYSTYPE),"MacBookPro")
CC       =  mpicc
CXX      =  mpiccxx
FC       =  $(CC) #mpifort  ## change this to "mpifort" for packages requiring linking secondary fortran code, currently -only- the helmholtz eos modules do this, so I leave it un-linked for now to save people the compiler headaches
OPTIMIZE = -O1 -funroll-loops
OPTIMIZE += -g -Wall # compiler warnings
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CXX     = mpic++
CHIMESINCL = -I/usr/local/include/sundials
CHIMESLIBS = -L/usr/local/lib -lsundials_cvode -lsundials_nvecserial
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I/usr/local/include #-I$(PORTINCLUDE)
GSL_LIBS = -L/usr/local/lib #-L$(PORTLIB)
FFTW_INCL= -I/usr/local/include
FFTW_LIBS= -L/usr/local/lib
HDF5INCL = -I/usr/local/include -DH5_USE_16_API #-I$(PORTINCLUDE) -DH5_USE_16_API
HDF5LIB  = -L/usr/local/lib -lhdf5 -lz #-L$(PORTLIB)
MPICHLIB = #
OPT     += -DDISABLE_ALIGNED_ALLOC -DCHIMES_USE_DOUBLE_PRECISION #
##
## update 2020: on more recent macs, MacPorts is not as useful a library installer.
##  I [PFH] switched over to homebrew. First you still need to install the extended XCode developer
##  tools and make sure all the appropriate extended tools, permissions, etc, are installed on your mac.
##  You can find tutorials online with simple searches like "how to install gcc and mpicc on osx"
##  rather than looking for something GIZMO-specific. Then look for how to install homebrew or
##  another package manager.
## Most of the compilation tools will then be available. You may need to install your own HDF5 libraries,
##  but this can be done with homebrew. Similarly some special code sub-modules require their own
##  packages. For example, chimes modules require sundials, which can be easily installed with
##  "brew install sundials". In the above, /usr/local/lib etc reflect the default homebrew install locations.
## FFTW is trickier, see the instructions below [special treatment is still needed]
##
## PFH: this is my own laptop installation (2013 MacBook Pro running Yosemite)
## --
## I have installed GSL and HDF5 through MacPorts (once you have it installed, just use:
## sudo port install gsl
## sudo port install hdf5
## then the shortcut PORTINCLUDE/PORTLIB are just my own links to the macports installation
##  directories. in my case they are the default:
## PORTLIB=/opt/local/lib
## PORTINCLUDE=/opt/local/include
## --
## Unfortunately, FFTW is more complicated, since macports, fink, and other repository systems
## do not support direct installation of the MPI version of FFTW2, which is what GIZMO needs
## if you want to run with PMGRID or POWERSPEC enabled (if not, it should just compile without
## FFTW just fine). Be sure to install FFTW 2.1.5: get it from http://www.fftw.org/
## then unpack it, go into the unpacked directory, and configure it with:
## ./configure --enable-mpi --enable-type-prefix --enable-float
## (this set of commands is important to install the correct version)
## then "make" and finally "sudo make install"
## that should install it to its default location, /usr/local/, which is where FFTW_INCL/FFW_LIBS
## are set to point (to the respective include and lib sub-directories). check to make sure you
## have the fftw libraries correctly installed.
## --
## With this done, and the code successfully compiled, you should be able to run it with
## mpirun -np X ./GIZMO 1>gizmo.out 2>gizmo.err &
## (here "X" is the number of processes you want to use, I'm assuming youre running from the
##  same directory with the code so ./GIZMO is just in the local directory, and GIZMO is the
##  compiled file, and the 1> and 2> commands route stdin and stderr to the desired files)
##--
## If you're having trouble, I recommend the excellent guides to installing GADGET-2 at:
## http://astrobites.org/2011/04/02/installing-and-running-gadget-2/
## and
## https://gauge.wordpress.com/2009/06/16/pitp-2009-installing-gadget2/
## (by Nathan Goldbaum and Javiera Guedes, respectively) -- the installation should be
## nearly identical here
##
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Wheeler")
CC       = mpicc ## gcc compilers, for intel replace this with mpiicc
CXX      = mpicpc ## gcc compilers, for intel replace this with mpiicpc
FC       = $(CC)
#OPTIMIZE = -Wall -g -O3 -xHOST -ipo -no-prec-div -fp-model fast=2 -fast-transcendentals -funroll-loops ## optimizations for intel compilers
##OPTIMIZE += -pg ## profiling for intel compilers
OPTIMIZE = -g -O1 -ffast-math -funroll-loops -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -fipa-cp-clone  ## optimizations for gcc compilers (1/2)
OPTIMIZE += -ftree-loop-distribute-patterns -fvect-cost-model -ftree-partial-pre   ## optimizations for gcc compilers (2/2)
#OPTIMIZE += -ftree-loop-distribute-patterns -ftree-slp-vectorize -fvect-cost-model -ftree-partial-pre   ## optimizations for gcc compilers (2/2)
#OPTIMIZE += -pg -fprofile -fprofile-arcs -ftest-coverage -fprofile-generate ## full profiling, for gcc compilers
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp # openmp required compiler flags
FC       = $(CC)
endif
GMP_INCL =
GMP_LIBS =
MKL_INCL =
MKL_LIBS =
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += -DUSE_MPI_IN_PLACE
## modules to load (intel compilers):
## module load intel/17 gsl/2.1 hdf5/1.8.17
## or for gcc compilers:
## module load gcc/5.3.0 openmpi/2.0.1 gsl/2.1 hdf5/1.8.17
## -- currently fftw2 is running from a custom install, but it should soon be fully module-supported (current module doesnt have mpi)
##     it is built in my directory with the config flags:
##      ./configure --prefix=$HOME/fftw_intel --enable-mpi --enable-type-prefix --enable-float CC=mpiicc CFLAGS='-O3 -fstrict-aliasing -malign-double -fomit-frame-pointer'
##      linked via the above FFTW2_HOME=$HOME/fftw_intel (where the libraries are installed)
##      (for a gcc compiler version, just omit the "CC" and "CFLAGS" flags above)
## in your job submission script, be sure to run gizmo with the following (if using intel compilers, otherwise this is irrelevant):
##   export I_MPI_DAPL_TRANSLATION_CACHE=0
##   before your "mpirun", (or include it in your .bashrc and source that before running): this is necessary or else the communication over DAPL will generate MPI memory errors
##
## note that with the newer builds of HDF5 on the machine, you may need to add the line
## export HDF5_DISABLE_VERSION_CHECK=1
##  to your .bashrc file, or it will think the wrong HDF5 file is linked and crash (even though it is fine)
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Zwicky")
CC       =  mpicc
CXX      =  mpicpc
FC       =  $(CC) ##mpiifort -nofor_main
OPTIMIZE = -O3 -funroll-loops
OPTIMIZE += -g -Wall # compiler warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp # openmp required compiler flags
FC       = $(CC)
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKL_HOME)/include
MKL_LIBS = -L$(MKL_HOME)/lib/em64t -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(GSL_HOME)/include
GSL_LIBS = -L$(GSL_HOME)/lib
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB = #
OPT     += # -DUSE_MPI_IN_PLACE
## modules to load: 
## module load intel/2011.4.191 impi/4.0.2.003 gsl/1.15-gcc HDF5 
##  -- the machine is quite picky, impi seems to be the only working mpi option right now
##  --  currently fftw2 isnt pre-installed, built library in my directory, with config flags:
##       ./configure --prefix=/home/phopkins/fftw --enable-mpi --enable-type-prefix --enable-float --with-gcc
##      linked via the above FFTW2_HOME=/home/phopkins/fftw (where the libraries are installed)
endif


#------------------------------------------------------------------------------
ifeq ($(SYSTYPE), "Edison")
CC       =  cc #instead ofmpicc
CXX      =  CC #instead of mpipc
FC       =  ftn #instead of $(CC)
OPTIMIZE =  -O3 -funroll-loops -ffast-math -finline-functions -funswitch-loops
OPTIMIZE += -g -Wall -fpredictive-commoning -fgcse-after-reload -fvect-cost-model
ifeq (OPENMP, $(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
endif
GMP_INCL =
GMP_LIBS =
MKL_INCL = -I$(INCLUDE)
MKL_LIBS = -L$(LIBRARY_PATH) -mkl=sequential
GSL_INCL = -I$(GSL_DIR)/include
GSL_LIBS = -L$(GSL_DIR)/lib
FFTW_INCL= -I$(FFTW_INC)
FFTW_LIBS= -L$(FFTW_DIR)
HDF5INCL = -I$(HDF5_INCLUDE_OPTS) -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_DIR)/lib -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
##
## modules to load: intel, impi, gsl, fftw/2.1.5.9, cray-hdf5
## note: there is a module called "hdf5" which will not work. Use cray-hdf5.
endif
#-----------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"SciNet")
CC       =  mpicc     # sets the C-compiler
OPTIMIZE =  -O1 -xHost -funroll-loops -no-prec-div -fast-transcendentals -fp-model fast=2 -ipo  # speed
#OPTIMIZE += -openmp -parallel -par-num-threads=4  # for openmp mode
OPTIMIZE += -g -debug parallel -Wall  # warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE +=-openmp -parallel  # openmp required compiler flags
endif
FC       =  $(CC)
GSL_INCL =  -I${SCINET_GSL_INC}
GSL_LIBS =  -L${SCINET_GSL_LIB} #-limf
FFTW_INCL=  -I${SCINET_FFTW_INC}
FFTW_LIBS=  -L${SCINET_FFTW_LIB}
MPICHLIB =
HDF5INCL =  -I${SCINET_HDF5_INC} -DH5_USE_16_API
HDF5LIB  =  -L${SCINET_HDF5_LIB} -lhdf5 -lz
MPICHLIB =
##
## Notes:
## 
### benchmarking suggests these optimizations, 256 cores with omp=4 or 2, no DOUBLE, multidomain=16 or 32, altogether gives best performance in
###   simple galaxy merger experiment (6x speedup over 16-core run with old-but-highly-optimized code).
##
## module load intel use.experimental openmpi/intel/1.6.0 gsl fftw/2.1.5-intel-openmpi hdf5/intel-openmpi/1.8.9
## NOTYPEPREFIX_FFTW should not be set on this machine
## 
## flags: 
## OPT      += -DNOCALLSOFSYSTEM -DMPICH_IGNORE_CXX_SEEK -DNO_ISEND_IRECV_IN_DOMAIN
##   -- these used to be recommended, with new compiler settings they don't seem necessary, but may help
## If memory problems crash the code, recommend small-scale chunking: MPISENDRECV_SIZELIMIT=10-100 (but this costs speed!)
##
## old options with intelmpi (not as good as openmpi):
##    module load intel intelmpi gsl fftw/2.1.5-intel-intelmpi4 use.experimental hdf5/intelmpi/1.8.9
##    OPTIMIZE =  -O2 -m64 -mt_mpi -openmp -xhost -g -debug parallel -mcmodel=medium -funroll-loops -Wall
##
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Iron")
CC       =   mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
#OPTIMIZE =   -std=c99 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -march=native
#OPTIMIZE =   -std=c99 -O3 -fno-tree-vectorize -march=native -Wno-implicit-function-declaration -Wno-cpp
OPTIMIZE =   -std=c99 -O2 -fno-tree-vectorize -march=native -Wno-implicit-function-declaration -Wno-cpp
OPTIMIZE += -g   #-Wall # compiler warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
endif
GSL_INCL =  -I$(GSL_BASE)/include
GSL_LIBS =  -L$(GSL_BASE)/lib
FFTW_INCL=  -I$(FFTW2_BASE)/include
FFTW_LIBS=  -L$(FFTW2_BASE)/lib
MPICHLIB =
HDF5INCL =  -I$(HDF5_BASE)/include -DH5_USE_16_API
HDF5LIB  =  -L$(HDF5_BASE)/lib -lhdf5 -lz
GMP_INCL =  #-I$(GMPDIR)/include
GMP_LIBs =  #-L$(GMPDIR)/lib
#module load slurm
#module add gcc
#module load openmpi2/2.0.2-hfi
#module add lib/hdf5/1.8.21-openmpi2
#module add lib/fftw2/2.1.5-openmpi2
#module add lib/gsl
endif
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Gordon")
CC       =  mpicc
CXX      =  mpicxx
FC       =  $(CC)  #mpif90 -nofor-main
OPTIMIZE =  -O3 -no-prec-div -xHOST
OPTIMIZE += -g -Wall # compiler warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -openmp # openmp required compiler flags
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I/opt/intel/composer_xe_2013_sp1.2.144/mkl/include
MKL_LIBS = -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib -mkl=sequential
GSL_INCL = -I/opt/gsl/2.1/intel/include
GSL_LIBS = -L/opt/gsl/2.1/intel/lib
FFTW_INCL= -I/opt/fftw/2.1.5/intel/mvapich2_ib/include
FFTW_LIBS= -L/opt/fftw/2.1.5/intel/mvapich2_ib/lib
HDF5INCL = -I/opt/hdf5/intel/mvapich2_ib/include -DH5_USE_16_API
HDF5LIB  = -L/opt/hdf5/intel/mvapich2_ib/lib -lhdf5 -lz
MPICHLIB = -L/opt/mvapich2/intel/ib/lib
OPT     += -DUSE_MPI_IN_PLACE
## modules to load:
## module load intel mvapich2_ib
## module load hdf5
## module load fftw/2.1.5
## module load gsl
endif
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Comet")
CC       =  mpicc
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -O3 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -g -Wall # compiler warnings
#OPTIMIZE += -parallel -openmp  # openmp (comment out this line if OPENMP not used)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp  # openmp required compiler flags
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I/opt/intel/composer_xe_2013_sp1.2.144/mkl/include
MKL_LIBS = -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib -mkl=sequential
##MKL_LIBS = -L/opt/mvapich2/intel/ib/lib -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I/opt/gsl/2.1/intel/include
GSL_LIBS = -L/opt/gsl/2.1/intel/lib
FFTW_INCL= -I/opt/fftw/2.1.5/intel/mvapich2_ib/include
FFTW_LIBS= -L/opt/fftw/2.1.5/intel/mvapich2_ib/lib
HDF5INCL = -I/opt/hdf5/intel/mvapich2_ib/include -DH5_USE_16_API
HDF5LIB  = -L/opt/hdf5/intel/mvapich2_ib/lib -lhdf5 -lz
MPICHLIB = -L/opt/mvapich2/intel/ib/lib
#MPICHLIB = -L/opt/openmpi/intel/ib/lib
OPT     += -DUSE_MPI_IN_PLACE
## modules to load:
## module load gsl intel hdf5 mvapich2_ib fftw/2.1.5
##  -- performance is very similar with impi (intel-mpi) instead of mpavich2,
##   if preferred use that with MPICHLIB line uncommented
## newest version of code needed for compatibility with calls in MPI-2 libraries
endif
#----------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Darter")
CC       =  cc
CXX      =  CC
FC       =  ftn -nofor_main
OPTIMIZE = -O3 -ipo -no-prec-div -static -xHost  # speed
OPTIMIZE += -g  # compiler warnings
#OPTIMIZE += -parallel -openmp  # openmp (comment out this line if OPENMP not used)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp  # openmp required compiler flags
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = #
GSL_LIBS = #
FFTW_INCL= #
FFTW_LIBS= #
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = #
MPICHLIB = #
OPT     += -DUSE_MPI_IN_PLACE
## modules to load:
## module swap PrgEnv-cray PrgEnv-intel
## module load intel gsl cray-hdf5-parallel fftw/2.1.5.9
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Quest")
CC       =  mpicc
CXX      =  mpic++
FC       =  $(CC)
OPTIMIZE = -O2 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp # openmp required compiler flags 
endif
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I/home/ajr882/sundials/include  
CHIMESLIBS = -L/home/ajr882/sundials/lib -lsundials_cvode -lsundials_nvecserial 
endif 
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKLROOT)/include
MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I/projects/b1026/pascal/software/gsl/1.16/include
GSL_LIBS = -L/projects/b1026/pascal/software/gsl/1.16/lib -lgsl -lgslcblas -lm
FFTW_INCL= -I/projects/b1026/pascal/software/fftw/2.1.5-mvp/include
FFTW_LIBS= -L/projects/b1026/pascal/software/fftw/2.1.5-mvp/lib
HDF5INCL = -I/projects/b1026/pascal/software/hdf5/1.8.12/include -DH5_USE_16_API
HDF5LIB  = -L/projects/b1026/pascal/software/hdf5/1.8.12/lib -lhdf5 -lz
MPICHLIB = -lmpich 
OPT     += -DUSE_MPI_IN_PLACE
#### modules to load:
#module load mpi/mvapich2-intel2013.2
#module use /projects/b1026/pascal/software/modules
#module load hdf5/1.8.12.1 gsl/1.16 fftw/2.1.5-mvp
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Quest-intel")
CC       =  mpicc
CXX      =  mpic++
FC       =  $(CC)
OPTIMIZE = -O2 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp -mt_mpi 
endif
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
CHIMESINCL = -I/home/sundials/include  
CHIMESLIBS = -L/home/sundials/lib -lsundials_cvode -lsundials_nvecserial 
endif 
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKLROOT)/include
MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= -I/home/libraries/fftw-2.1.5_install/include 
FFTW_LIBS= -L/home/libraries/fftw-2.1.5_install/lib 
HDF5INCL = -DH5_USE_16_API 
HDF5LIB  = -lhdf5 -lz
MPICHLIB = 
OPT     += -DUSE_MPI_IN_PLACE
#### modules to load:
#module load intel/2013.2
#module load mpi/intel-mpi-4.1.0
#module load hdf5/1.8.12-serial
#module load gsl/1.16-intel
#module load fftw/2.1.5-intel
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Titan")
CC       =  cc
CXX      =  CC
FC       =  $(CC) #ftn
OPTIMIZE = -O3 -ipo -funroll-loops -no-prec-div -fp-model fast=2 -static
OPTIMIZE += -g
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -openmp # (intel) openmp required compiler flags
FC       = $(CC)
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I$(GSL_DIR)/include
GSL_LIBS = -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm
FFTW_INCL= -I/opt/cray/fftw/2.1.5.8/include
FFTW_LIBS= -L/opt/cray/fftw/2.1.5.8/lib
HDF5INCL = -I$(HDF5_DIR)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_DIR)/lib -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
endif
## in your .bashrc file, include
## module swap PrgEnv-pgi PrgEnv-intel
## module load cray-hdf5-parallel fftw/2.1.5.8 gsl mercurial



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"BlueWaters")
CC       =  cc
CXX      =  CC
FC       =  $(CC) #ftn
#OPTIMIZE = -O3 -ipo -funroll-loops -no-prec-div -fp-model fast=2 -static
OPTIMIZE = -fast -no-ipo
#OPTIMIZE += -g
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -qopenmp -fopenmp # (intel) openmp required compiler flags
FC       = $(CC)
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I$(GSL_DIR)/include
GSL_LIBS = -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm
FFTW_INCL= -I$(FFTW_DIR)/include
FFTW_LIBS= -L$(FFTW_DIR)/lib
HDF5INCL = -I$(HDF5_DIR)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_DIR)/lib -lhdf5 -lz
MPICHLIB = -lmpich_intel
OPT     += -DUSE_MPI_IN_PLACE
endif
## in your .bashrc file, include
##    module swap PrgEnv-cray PrgEnv-intel
##    module load cray-hdf5-parallel fftw/2.1.5.9 gsl bwpy szip
##    module swap intel intel/17.0.4.196
##    umask 022
##    ulimit -s unlimited
##
## NOTE: everything above is for intel compilers. PGI performs almost as well (within 5%) and has OpenACC, much better for GPU/xk nodes.
##  to use PGI, replace PrgEnv-intel with PrgEnv-pgi everywhere here, and replace the optimize line with:
##  OPTIMIZE = -m64 -mcmodel=medium -Mdalign -Mllalign -Munroll -O3 -fastsse -Mipa=fast -fast -Msafeptr
##
## example run submission script:
##
##  #!/bin/bash
##  #PBS -N NAME                 # name as desired
##  #PBS -q debug                # debug or normal or whatever queue
##  #PBS -l nodes=NNODES:ppn=32:xe   # PPN=32, essentially always; nodes=NNODES (set to number)
##  ##PBS -l flags=preemptee     # makes job pre-emptable: receives discount for time used
##  ##PBS -l minwclimit=23:10:00  -l walltime=48:00:00  # flexible clock time (set min/max) - often gets on machine faster
##  #PBS -l walltime=00:15:00    # fixed run-time, set appropriately
##  #PBS -j oe                   # outputs
##  #PBS -o run.out              # outputs
##  cd $PBS_O_WORKDIR            # run in submission directory
##  . /opt/modules/default/init/bash # NEEDED to add module commands to shell
##  source $HOME/.bashrc         # source your bash shortcuts, etc
##  module swap PrgEnv-cray PrgEnv-intel    # make sure correct compilers [intel here] loaded
##  module load intel cray-hdf5-parallel fftw/2.1.5.9 gsl szip   # make sure correct modules loaded
##  module swap intel intel/17.0.4.196
##  export KMP_AFFINITY=disabled  # include if and only if you are using intel17+ compilers, needed to leverage newer cache optimizations
##  export OMP_NUM_THREADS=NOMP        # turn on for OpenMP runs, NOMP=number of threads (=1 if not using OpenMP, in everything below)
##  aprun -S SNODENUM -cc numa_node -N NMPINODE -d NOMP -n NTOTAL ./GIZMO ./params.txt 1>gizmo.out 2>gizmo.err
##  ## (in the above:
##  ##     -S SNODENUM : set SNODENUM = NMPINODE/4 = 8/NOMP. This is how many MPI tasks per die (4 dies per node) -- prevents memory issues
##  ##        which can dramatically slow down code by threads spawning without local access to memory.
##  ##        If using >8 OpenMP threads (so this would be <1), this command should be removed.
##  ##     -N NMPINODE : set NMPINODE = 32/NOMP = number of MPI tasks per node. BW has 32 integer-cores per node, with this structure all can be used; even though every pair shares 1 FP core.
##  ##     -d NOMP = specifies spacing/skipping of CPUs for multi-threaded alignment: set equal to NOMP
##  ##     NTOTAL = NNODES*NMPINODE = total number of MPI tasks (not number of processors)
##  ##     (in previous notes, had "-ss" do memory-binding, and "-j 1" prevents integer-cores sharing the same FPU, but this is better handled by the numa_node memory strategy)
##  ##
##
##  For FIRE runs, generally reasonable with OPENMP off or set =2,4,8. Any higher over-runs the die, not good.
##  May be some (small) performance improvement from using newest Intel MPI: in .bashrc and run-file add line "module swap intel intel/17.0.4.196" after the lines
##   that load intel compilers, and in run-file, before OPENMP definition, include line "export KMP_AFFINITY=disabled" (this is necessary with OpenMP in intel17+ or performance will suffer)
##


#------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Mira")
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
CC       = mpixlc_r # xl compilers appear to give significant speedup vs gcc
CXX      = mpixlcxx_r # _r for thread-safe versions, desired with openmp
OPTIMIZE = -openmp -qsmp=omp:noauto # -fopenmp for gcc or bgclang
else
CC       = mpixlc
CXX      = mpixlcxx
OPTIMIZE =
endif
FC       =  $(CC)
OPTIMIZE += -O3 -static -qnostrict -lm -ldl #-lpthread
OPTIMIZE += -g #-qlist -qsource -qreport -qlistopt # compiler warnings: qlist, etc produce list of opts
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I$(MIRA_GSL_INC)
GSL_LIBS = -lm -ldl -lpthread -L$(MIRA_GSL_LIB) -lgsl
FFTW_INCL= -I$(MIRA_FFTW2_INC)
FFTW_LIBS= -L$(MIRA_FFTW2_LIB)
HDF5INCL = -I$(MIRA_HDF5_INC) -DH5_USE_16_API -I$(MIRA_SZIP_INC) -I$(MIRA_LZIP_INC)
HDF5LIB  = -L$(MIRA_SZIP_LIB) -lszip -L$(MIRA_LZIP_LIB) -lz -L$(MIRA_HDF5_LIB) -lhdf5 -lz -lszip
MPICHLIB = #
OPT     += -DUSE_MPI_IN_PLACE -DREDUCE_TREEWALK_BRANCHING
##
## in .bashrc, need to define environmental variables:
##   export MIRA_HDF5_INC=/soft/libraries/hdf5/current/cnk-xl/current/include
##   export MIRA_HDF5_LIB=/soft/libraries/hdf5/current/cnk-xl/current/lib
##   export MIRA_GSL_INC=/soft/libraries/3rdparty/gsl/1.9/xl/include
##   export MIRA_GSL_LIB=/soft/libraries/3rdparty/gsl/1.9/xl/lib
##   export MIRA_SZIP_INC=/soft/libraries/alcf/current/xl/SZIP/include
##   export MIRA_SZIP_LIB=/soft/libraries/alcf/current/xl/SZIP/lib
##   export MIRA_LZIP_INC=/soft/libraries/alcf/current/xl/ZLIB/include
##   export MIRA_LZIP_LIB=/soft/libraries/alcf/current/xl/ZLIB/lib
##   export MIRA_FFTW2_INC=/home/phopkins/fftw/include
##   export MIRA_FFTW2_LIB=/home/phopkins/fftw/lib
##   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MIRA_LZIP_LIB:$MIRA_SZIP_LIB:$MIRA_FFTW2_LIB:$MIRA_GSL_LIB:$MIRA_HDF5_LIB:/bgsys/drivers/ppcfloor/comm/lib
##
## for HDF5,GSL,LZIP,SZIP these link to the current general-use versions of these libraries. the last command (adding these to the LD_LIBRARY_PATH)
##   is also critical, since they are not in the paths by default and it will be unable to find them even with the links above.
## for FFTW2, the pre-compiled libraries do not contain the correct mpi libraries, so you will have to compile your own. FFTW2 is installed and
##   compiled in the directory shown for me: you have to install it and link it yourself since the directory cannot be shared. for the gcc
##   compilers above, FFTW should be compiled with the following settings:
##     for xl compilers:
##       ./configure --prefix=$HOME/fftw --enable-mpi --enable-type-prefix --enable-float LDFLAGS=-L$HOME/lib CFLAGS=-I$HOME/include CC=mpixlc
##
## also in your .soft file, you want to enable:
##   for XL compilers:
##     +mpiwrapper-xl
##     +python
##     @default
## to load the mpi compilers and mpi wrappers, and MPICH libraries (python there is optional)
## xl appears to provide some improvement over gcc; xl-ndebug provides no noticeable further improvement, despite being more unsafe
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ranger_intel")
CC       =  mpicc
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -O3 -xO -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -parallel -openmp  # openmp
OPTIMIZE += -g -Wall -debug parallel # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB =      # must be empty if using openmpi
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
##
## Notes:
## 
## include the following in your .bashrc file (there is no default fftw2 module):
## module load intel/10.1 openmpi/1.2.4 gmp gsl hdf5 #now have to add fftw2 manually
## export TACC_FFTW2_INC=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/include
## export TACC_FFTW2_LIB=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/lib
## export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
##
## Options
## OPT += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
##   are not necessary, but may improve stability in some cases
##
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ranger_pgi")
CC       =  mpicc 
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -tp barcelona-64 -fast -Mipa=fast,inline -Munroll -Mvect -O4
OPTIMIZE += -mp -Mconcur  # openmp
OPTIMIZE += -Wall  # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPT     += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
## 
## Notes:
##
## include the following in your .bashrc file:
##   module load pgi mvapich gmp gsl fftw2 hdf5
##   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
## 
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK 
FC       =  $(CC)
OPTIMIZE = -g -O2 -Wall -Wno-unused-but-set-variable
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"antares")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
FC       =  $(CC)
OPTIMIZE = -g -O2 -Wall -Wno-unused-but-set-variable
GSL_INCL = -I/home/ptorrey/local/gsl-2.1/include
GSL_LIBS = -L/home/ptorrey/local/gsl-2.1/lib -lgsl -lm
FFTW_INCL= -I/home/ptorrey/local/fftw-2.1.5/include
FFTW_LIBS= -L/home/ptorrey/local/fftw-2.1.5/libGSL_INCL
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"CITA")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/include/gsl
GSL_LIBS =  -L/usr/lib/libgsl
FFTW_INCL=  -I/opt/fftw-2.1.5/include
FFTW_LIBS=  -L/opt/fftw-2.1.5/lib
MPICHLIB =  -L/usr/lib/libmpi
HDF5INCL =  -I/usr/include
HDF5LIB  =  -L/usr/lib/libhdf5 -static -lhdf5 -lz
endif 


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Sauron-gcc")
CC       =   mpicc.gcc   # sets the C-compiler
OPTIMIZE =   -O3 -funroll-loops -march=k8 -msse2 -static
GSL_INCL =   -I/usr/local/gsl.gcc/include
GSL_LIBS =   -L/usr/local/gsl.gcc/lib -static -lgsl -lgslcblas
FFTW_INCL=   -I/usr/local/fftw.gcc/include
FFTW_LIBS=   -L/usr/local/fftw.gcc/lib -static -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
MPICHLIB =
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Sauron")
CC       =  mpicc  -m64 # sets the C-compiler
CXX      =  mpiCC  -m64
OPTIMIZE =   -g
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   # sets the C-compiler
CXX      =  mpiCC
OPTIMIZE =   -g -Wall -fopenmp
# GSL_INCL =  -I/usr/common/pdsoft/include
# GSL_LIBS =  -L/usr/common/pdsoft/lib
GSL_INCL =  -I/afs/mpa/home/volker/Libs/include
GSL_LIBS =  -L/afs/mpa/home/volker/Libs/lib
FFTW_INCL=  -I/afs/mpa/home/volker/Libs/include
FFTW_LIBS=  -L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
MPICHLIB =
HDF5INCL =  -I/afs/mpa/home/volker/Libs/include
HDF5LIB  =  -L/afs/mpa/home/volker/Libs/lib -lhdf5 -lz 
OPT     +=  -DOLD_HDF5
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
CORE_OBJS =	main.o accel.o  timestep.o init.o restart.o io.o \
			predict.o global.o begrun.o run.o allvars.o read_ic.o \
			domain.o driftfac.o kicks.o ngb.o compile_time_info.o merge_split.o

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
                gravity/ags_hsml.o \
                gravity/binary.o

HYDRO_OBJS = 	hydro/hydro_toplevel.o \
				hydro/density.o \
				hydro/gradients.o \
				turb/dynamic_diffusion.o \
				turb/dynamic_diffusion_velocities.o \
				turb/turb_driving.o \
				turb/turb_powerspectra.o

EOSCOOL_OBJS =  cooling/cooling.o \
				cooling/grackle.o \
				eos/eos.o \
				eos/cosmic_ray_fluid/cosmic_ray_alfven.o \
				eos/cosmic_ray_fluid/cosmic_ray_utilities.o \
				solids/elastic_physics.o \
				solids/grain_physics.o \
				nuclear/nuclear_network_solver.o \
				nuclear/nuclear_network.o 

STARFORM_OBJS = galaxy_sf/sfr_eff.o \
                galaxy_sf/stellar_evolution.o \
                galaxy_sf/mechanical_fb.o \
                galaxy_sf/thermal_fb.o \
                galaxy_sf/radfb_local.o \
                galaxy_sf/dm_dispersion_hsml.o

SINK_OBJS = galaxy_sf/blackholes/blackhole.o \
            galaxy_sf/blackholes/blackhole_util.o \
            galaxy_sf/blackholes/blackhole_environment.o \
            galaxy_sf/blackholes/blackhole_feed.o \
            galaxy_sf/blackholes/blackhole_swallow_and_kick.o

RHD_OBJS =  radiation/rt_utilities.o \
			radiation/rt_CGmethod.o \
			radiation/rt_source_injection.o \
			radiation/rt_chem.o \
			radiation/rt_cooling.o

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

## name of executable and optimizations
EXEC   = GIZMO
OPTIONS = $(OPTIMIZE) $(OPT)

## combine all the objects above
OBJS  = $(CORE_OBJS) $(SYSTEM_OBJS) $(GRAVITY_OBJS) $(HYDRO_OBJS) \
		$(EOSCOOL_OBJS) $(STARFORM_OBJS) $(SINK_OBJS) $(RHD_OBJS) \
		$(FOF_OBJS) $(MISC_OBJS)

## fortran recompiler block
FOPTIONS = $(OPTIMIZE) $(FOPT)
FOBJS =

## include files needed at compile time for the above objects
INCL    += 	allvars.h \
			proto.h \
			gravity/forcetree.h \
			gravity/myfftw3.h \
			domain.h \
			system/myqsort.h \
			kernel.h \
			eos/eos.h \
			galaxy_sf/blackholes/blackhole.h \
			structure/fof.h \
			structure/subfind/subfind.h \
			cooling/cooling.h \
			nuclear/nuclear_network.h \
			Makefile


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

# chimes files are treated as special for now because they require special external libraries (e.g. sundials) that are otherwise not
#   used anywhere else in the code, and have not had their macro logic cleaned up to allow appropriate compilation without chimes flags enabled
ifeq (CHIMES,$(findstring CHIMES,$(CONFIGVARS)))
OBJS    += cooling/chimes/chimes.o cooling/chimes/chimes_cooling.o cooling/chimes/init_chimes.o cooling/chimes/rate_equations.o cooling/chimes/update_rates.o 
INCL    += cooling/chimes/chimes_interpol.h cooling/chimes/chimes_proto.h cooling/chimes/chimes_vars.h 
endif

# if HDF5 explicitly disabled, remove the linked libraries
ifeq (IO_DISABLE_HDF5,$(findstring IO_DISABLE_HDF5,$(CONFIGVARS)))
HDF5INCL =
HDF5LIB  =
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
CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(GMP_INCL) \
         $(GRACKLEINCL) $(CHIMESINCL)

LIBS = $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas \
	   $(FFTW_LIBS) $(FFTW_LIBNAMES) -lm $(GRACKLELIBS) $(CHIMESLIBS)

ifeq (PTHREADS_NUM_THREADS,$(findstring PTHREADS_NUM_THREADS,$(CONFIGVARS))) 
LIBS += -lpthread
endif

$(EXEC): $(OBJS) $(FOBJS)  
	$(FC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)  $(CONFIG)  compile_time_info.c

$(FOBJS): %.o: %.f90
	$(FC) $(OPTIMIZE) -c $< -o $@

compile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) *.oo *.c~ compile_time_info.c GIZMO_config.h


