#ifndef __MY_FFTW3_H__
#define __MY_FFTW3_H__

#include <fftw3-mpi.h>

#define cmplx_re(c) ((c)[0])
#define cmplx_im(c) ((c)[1])

#ifdef DOUBLEPRECISION_FFTW 
  #define fftw_real			    double 
  #define MPI_TYPE_FFTW			    MPI_DOUBLE 
#else 
  #define fftw_real			    float 
  #define MPI_TYPE_FFTW			    MPI_FLOAT 
  #define fftw_complex			    fftwf_complex 
  #define fftw_mpi_local_size_3d_transposed fftwf_mpi_local_size_3d_transposed
  #define fftw_mpi_init			    fftwf_mpi_init
  #define fftw_plan			    fftwf_plan 
  #define fftw_mpi_local_size_3d	    fftwf_mpi_local_size_3d 
  #define fftw_mpi_plan_dft_r2c_3d	    fftwf_mpi_plan_dft_r2c_3d 
  #define fftw_mpi_plan_dft_c2r_3d	    fftwf_mpi_plan_dft_c2r_3d 
  #define fftw_execute			    fftwf_execute 
  #define fftw_destroy_plan		    fftwf_destroy_plan
#endif

#endif
