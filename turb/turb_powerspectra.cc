#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include "../declarations/gpu_rng.h"

#include "../declarations/allvars.h"
#include "../core/proto.h"

/*
 *  This code was originally written for GADGET3 by Andreas Bauer; it has been
 *   modified slightly by Phil Hopkins for GIZMO, but is largely intact.
 */

#if defined(TURB_DRIVING_SPECTRUMGRID) && defined(BOX_PERIODIC) && (defined(TURB_DRIVING))
#define  TURB_DRIVING_SPECTRUMGRID2 (2*(TURB_DRIVING_SPECTRUMGRID/2 + 1))
#include "../gravity/myfftw3.h"

#include <vector>
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "../system/gpu_particles_arena.h"

#if (TURB_DRIVING_SPECTRUMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif
static fftw_plan fft_velx_plan, fft_vely_plan, fft_velz_plan;
static fftw_plan fft_svelx_plan, fft_svely_plan, fft_svelz_plan;
static fftw_plan fft_vrhox_plan, fft_vrhoy_plan, fft_vrhoz_plan;
static fftw_plan fft_vortx_plan, fft_vorty_plan, fft_vortz_plan;
static fftw_plan fft_dis1field_plan, fft_dis2field_plan; 
static fftw_plan fft_rand_plan; 
static fftw_plan fft_dens_plan; 
static ptrdiff_t slabstart_x, nslab_x, slabstart_y, nslab_y;
static ptrdiff_t fftsize, maxfftsize;
static MPI_Datatype MPI_TYPE_PTRDIFF; 
static fftw_real *velfield[3];
#ifdef TURB_DIFF_DYNAMIC
static fftw_real *velbarfield[3];
static fftw_real *velhatfield[3];
#endif
static fftw_real *smoothedvelfield[3];
static fftw_real *vorticityfield[3];
static fftw_real *velrhofield[3];
static fftw_real *dis1field;
static fftw_real *dis2field;
static fftw_real *densityfield;
static fftw_real *randomfield;
static fftw_real *workspace;
static float    *RandomValue;
static fftw_complex *fft_of_field;
static float *powerspec_turb_nearest_distance, *powerspec_turb_nearest_rkern;
void powerspec_turb_calc_and_bin_spectrum(fftw_plan plan, fftw_real *field, int flag);

#define BINS_PS  2000	                 	/* number of bins for power spectrum computation */

static long long CountModes[BINS_PS];
static double    SumPower[BINS_PS];
static double    Power[BINS_PS];
static double    Kbin[BINS_PS];
static double    K0, K1;
static double    binfac;
static double    vel_disp[3];
#ifdef TURB_DIFF_DYNAMIC
static double    velbar_disp[3];
static double    velhat_disp[3];
#endif
static double    velrho_disp[3];
static double    empty_disp[3] = {0, 0, 0};




void powerspec_turb(int filenr)
{
  int i;
  char fname[DEFAULT_PATH_BUFFERSIZE_TOUSE];

  if(ThisTask == 0)
    printf("Start turbulent powerspec computation\n");

  double tstart, tend;
  tstart = my_second();

  // define MPI_TYPE_PTRDIFF */
  if (sizeof(ptrdiff_t) == sizeof(long long)) {
    MPI_TYPE_PTRDIFF = MPI_LONG_LONG; 
  } else if (sizeof(ptrdiff_t) == sizeof(long)) {
    MPI_TYPE_PTRDIFF = MPI_LONG; 
  } else if (sizeof(ptrdiff_t) == sizeof(int)) {
    MPI_TYPE_PTRDIFF = MPI_INT; 
  }

  fftsize = fftw_mpi_local_size_3d_transposed(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID/2 + 1, 
	  MPI_COMM_WORLD, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y); 
  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_TYPE_PTRDIFF, MPI_MAX, MPI_COMM_WORLD); 

  /* allocate the memory to hold the FFT fields */

  velfield[0] = (fftw_real *) mymalloc("velfield[0]", maxfftsize * sizeof(fftw_real));
  velfield[1] = (fftw_real *) mymalloc("velfield[1]", maxfftsize * sizeof(fftw_real));
  velfield[2] = (fftw_real *) mymalloc("velfield[2]", maxfftsize * sizeof(fftw_real));

#ifdef TURB_DIFF_DYNAMIC
  velbarfield[0] = (fftw_real *) mymalloc("velbarfield[0]", maxfftsize * sizeof(fftw_real));
  velbarfield[1] = (fftw_real *) mymalloc("velbarfield[1]", maxfftsize * sizeof(fftw_real));
  velbarfield[2] = (fftw_real *) mymalloc("velbarfield[2]", maxfftsize * sizeof(fftw_real));

  velhatfield[0] = (fftw_real *) mymalloc("velhatfield[0]", maxfftsize * sizeof(fftw_real));
  velhatfield[1] = (fftw_real *) mymalloc("velhatfield[1]", maxfftsize * sizeof(fftw_real));
  velhatfield[2] = (fftw_real *) mymalloc("velhatfield[2]", maxfftsize * sizeof(fftw_real));
#endif

  smoothedvelfield[0] = (fftw_real *) mymalloc("smoothedvelfield[0]", maxfftsize * sizeof(fftw_real));
  smoothedvelfield[1] = (fftw_real *) mymalloc("smoothedvelfield[1]", maxfftsize * sizeof(fftw_real));
  smoothedvelfield[2] = (fftw_real *) mymalloc("smoothedvelfield[2]", maxfftsize * sizeof(fftw_real));

  velrhofield[0] = (fftw_real *) mymalloc("velrhofield[0]", maxfftsize * sizeof(fftw_real));
  velrhofield[1] = (fftw_real *) mymalloc("velrhofield[1]", maxfftsize * sizeof(fftw_real));
  velrhofield[2] = (fftw_real *) mymalloc("velrhofield[2]", maxfftsize * sizeof(fftw_real));

  vorticityfield[0] = (fftw_real *) mymalloc("vorticityfield[0]", maxfftsize * sizeof(fftw_real));
  vorticityfield[1] = (fftw_real *) mymalloc("vorticityfield[1]", maxfftsize * sizeof(fftw_real));
  vorticityfield[2] = (fftw_real *) mymalloc("vorticityfield[2]", maxfftsize * sizeof(fftw_real));

  dis1field = (fftw_real *) mymalloc("dis1field", maxfftsize * sizeof(fftw_real));
  dis2field = (fftw_real *) mymalloc("dis2field", maxfftsize * sizeof(fftw_real));
  randomfield = (fftw_real *) mymalloc("randomfield", maxfftsize * sizeof(fftw_real));

  densityfield = (fftw_real *) mymalloc("densityfield", maxfftsize * sizeof(fftw_real));

  fft_velx_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID,
	  velfield[0], (fftw_complex *) velfield[0], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_vely_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  velfield[1], (fftw_complex *) velfield[1], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_velz_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  velfield[2], (fftw_complex *) velfield[2], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_svelx_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  smoothedvelfield[0], (fftw_complex *) smoothedvelfield[0], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_svely_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  smoothedvelfield[1], (fftw_complex *) smoothedvelfield[1], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_svelz_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  smoothedvelfield[2], (fftw_complex *) smoothedvelfield[2], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_vrhox_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  velrhofield[0], (fftw_complex *) velrhofield[0], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_vrhoy_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  velrhofield[1], (fftw_complex *) velrhofield[1], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_vrhoz_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  velrhofield[2], (fftw_complex *) velrhofield[2], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_vortx_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  vorticityfield[0], (fftw_complex *) vorticityfield[0], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_vorty_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  vorticityfield[1], (fftw_complex *) vorticityfield[1], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_vortz_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  vorticityfield[2], (fftw_complex *) vorticityfield[2], 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_dis1field_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  dis1field, (fftw_complex *) dis1field, 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_dis2field_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  dis2field, (fftw_complex *) dis2field, 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_rand_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  randomfield, (fftw_complex *) randomfield, 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_dens_plan = fftw_mpi_plan_dft_r2c_3d(TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, TURB_DRIVING_SPECTRUMGRID, 
	  densityfield, (fftw_complex *) densityfield, 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  workspace = (fftw_real *) mymalloc("workspace", maxfftsize * sizeof(fftw_real));

  memset(velfield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(velfield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(velfield[2], 0, maxfftsize * sizeof(fftw_real));

#ifdef TURB_DIFF_DYNAMIC
  memset(velbarfield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(velbarfield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(velbarfield[2], 0, maxfftsize * sizeof(fftw_real));

  memset(velhatfield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(velhatfield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(velhatfield[2], 0, maxfftsize * sizeof(fftw_real));
#endif

  memset(smoothedvelfield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(smoothedvelfield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(smoothedvelfield[2], 0, maxfftsize * sizeof(fftw_real));

  memset(velrhofield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(velrhofield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(velrhofield[2], 0, maxfftsize * sizeof(fftw_real));

  memset(vorticityfield[0], 0, maxfftsize * sizeof(fftw_real));
  memset(vorticityfield[1], 0, maxfftsize * sizeof(fftw_real));
  memset(vorticityfield[2], 0, maxfftsize * sizeof(fftw_real));

  memset(dis1field, 0, maxfftsize * sizeof(fftw_real));
  memset(dis2field, 0, maxfftsize * sizeof(fftw_real));
  memset(randomfield, 0, maxfftsize * sizeof(fftw_real));

  memset(densityfield, 0, maxfftsize * sizeof(fftw_real));

  RandomValue = (float *) mymalloc("RndField", N_gas * sizeof(float));

  for(i=0; i < N_gas; i++)
    RandomValue[i] = (float)gizmo_gpu_rand_gaussian((uint64_t)(42 + ThisTask), (uint64_t)i);

  powerspec_turb_obtain_fields();
 
  powerspec_turb_calc_dispersion();



  /* Now compute the power spectrum of the velocities */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_velx_plan, velfield[0], 1);   /* only here the modes are counted */
  powerspec_turb_calc_and_bin_spectrum(fft_vely_plan, velfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(fft_velz_plan, velfield[2], 0);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_vel_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, vel_disp);


#ifdef TURB_DIFF_DYNAMIC
  /* Now compute the power spectrum of the velbar quantities */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(velbarfield[0], 1);   /* only here the modes are counted */
  powerspec_turb_calc_and_bin_spectrum(velbarfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(velbarfield[2], 0);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_velbar_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, velbar_disp);

  /* Now compute the power spectrum of the velhat quantities */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(velhatfield[0], 1);   /* only here the modes are counted */
  powerspec_turb_calc_and_bin_spectrum(velhatfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(velhatfield[2], 0);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_velhat_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, velhat_disp);
#endif

  /* Now compute the power spectrum of the smoothed velocities */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_svelx_plan, smoothedvelfield[0], 1);   /* only here the modes are counted */
  powerspec_turb_calc_and_bin_spectrum(fft_svely_plan, smoothedvelfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(fft_svelz_plan, smoothedvelfield[2], 0);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_smoothedvel_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, vel_disp);





  /* now compute the power spectrum of the sqrt(rho)-weighted veloicty */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_vrhox_plan, velrhofield[0], 1);
  powerspec_turb_calc_and_bin_spectrum(fft_vrhoy_plan, velrhofield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(fft_vrhoz_plan, velrhofield[2], 0);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_velrho_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, velrho_disp);



  /* now compute the power spectrum of the vorticity */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_vortx_plan, vorticityfield[0], 1);
  powerspec_turb_calc_and_bin_spectrum(fft_vorty_plan, vorticityfield[1], 0);
  powerspec_turb_calc_and_bin_spectrum(fft_vortz_plan, vorticityfield[2], 0);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_vorticity_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, velrho_disp);




  /* Now compute the power spectrum of the dissipation1 */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_dis1field_plan, dis1field, 1);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_dis1_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, empty_disp);




  /* Now compute the power spectrum of the dissipation2 */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_dis2field_plan, dis2field, 1);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_dis2_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, empty_disp);


  /* Now compute the power spectrum of the random field */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_rand_plan, randomfield, 1);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_random_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, empty_disp);

  /* Now compute the power spectrum of the density field */

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      CountModes[i] = 0;
    }

  powerspec_turb_calc_and_bin_spectrum(fft_dens_plan, densityfield, 1);

  powerspec_turb_collect();

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/powerspec_density_%03d.txt", All.OutputDir, filenr);
  powerspec_turb_save(fname, empty_disp);

  myfree(RandomValue);

  myfree(workspace);
  myfree(densityfield);
  myfree(randomfield);
  myfree(dis2field);
  myfree(dis1field);
  myfree(vorticityfield[2]);
  myfree(vorticityfield[1]);
  myfree(vorticityfield[0]);
  myfree(velrhofield[2]);
  myfree(velrhofield[1]);
  myfree(velrhofield[0]);
  myfree(smoothedvelfield[2]);
  myfree(smoothedvelfield[1]);
  myfree(smoothedvelfield[0]);
#ifdef TURB_DIFF_DYNAMIC
  myfree(velhatfield[2]);
  myfree(velhatfield[1]);
  myfree(velhatfield[0]);
  myfree(velbarfield[2]);
  myfree(velbarfield[1]);
  myfree(velbarfield[0]);
#endif
  myfree(velfield[2]);
  myfree(velfield[1]);
  myfree(velfield[0]);

  fftw_destroy_plan(fft_dens_plan);
  fftw_destroy_plan(fft_rand_plan); 

  fftw_destroy_plan(fft_dis2field_plan); 
  fftw_destroy_plan(fft_dis1field_plan); 

  fftw_destroy_plan(fft_vortz_plan); 
  fftw_destroy_plan(fft_vorty_plan); 
  fftw_destroy_plan(fft_vortx_plan); 

  fftw_destroy_plan(fft_vrhoz_plan); 
  fftw_destroy_plan(fft_vrhoy_plan); 
  fftw_destroy_plan(fft_vrhox_plan); 

  fftw_destroy_plan(fft_svelz_plan); 
  fftw_destroy_plan(fft_svely_plan); 
  fftw_destroy_plan(fft_svelx_plan); 

  fftw_destroy_plan(fft_velz_plan); 
  fftw_destroy_plan(fft_vely_plan); 
  fftw_destroy_plan(fft_velx_plan); 

  tend = my_second();
  
  PRINT_STATUS("end turbulent power spectra  took %g seconds", timediff(tstart, tend));
}


void powerspec_turb_calc_and_bin_spectrum(fftw_plan fplan, fftw_real *field, int flag)
{
  double k2, kx, ky, kz;
  int x, y, z, zz, ip;
  
  K0 = 2 * M_PI / All.BoxSize;	                        /* minimum k */
  K1 = K0 * TURB_DRIVING_SPECTRUMGRID / 2;	                                /* maximum k */
  binfac = BINS_PS / (log(K1) - log(K0));

  /* Do the FFT of the velocity_field */  /* rhogrid -> velfield */
  
  fftw_execute(fplan); 
  
  fft_of_field = (fftw_complex *) field;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) \
  reduction(+:SumPower[0:BINS_PS]) reduction(+:CountModes[0:BINS_PS]) \
  private(z,zz,kx,ky,kz,k2,ip)
#endif
  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < TURB_DRIVING_SPECTRUMGRID; x++)
      for(z = 0; z < TURB_DRIVING_SPECTRUMGRID; z++)
	{
	  zz = z;
	  if(z >= TURB_DRIVING_SPECTRUMGRID / 2 + 1)
	    zz = TURB_DRIVING_SPECTRUMGRID - z;

	  if(x > TURB_DRIVING_SPECTRUMGRID / 2)
	    kx = x - TURB_DRIVING_SPECTRUMGRID;
	  else
	    kx = x;
	  if(y > TURB_DRIVING_SPECTRUMGRID / 2)
	    ky = y - TURB_DRIVING_SPECTRUMGRID;
	  else
	    ky = y;
	  if(z > TURB_DRIVING_SPECTRUMGRID / 2)
	    kz = z - TURB_DRIVING_SPECTRUMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  ip = TURB_DRIVING_SPECTRUMGRID * (TURB_DRIVING_SPECTRUMGRID / 2 + 1) * (y - slabstart_y) + (TURB_DRIVING_SPECTRUMGRID / 2 + 1) * x + zz;

	  double po = (cmplx_re(fft_of_field[ip]) * cmplx_re(fft_of_field[ip])
		       + cmplx_im(fft_of_field[ip]) * cmplx_im(fft_of_field[ip])) / pow(TURB_DRIVING_SPECTRUMGRID, 6);

	  if(k2 > 0)
	    {
	      if(k2 < (TURB_DRIVING_SPECTRUMGRID / 2.0) * (TURB_DRIVING_SPECTRUMGRID / 2.0))
		{
		  double k = sqrt(k2) * 2 * M_PI / All.BoxSize;

		  if(k >= K0 && k < K1)
		    {
		      int bin = log(k / K0) * binfac;

		      SumPower[bin] += po;

		      if(flag)
			CountModes[bin] += 1;
		    }
		}
	    }
	}
}



void powerspec_turb_collect(void)
{
  int i, n;
  long long int *countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_PS * sizeof(long long));
  double *powerbuf = (double *) mymalloc("powerbuf", NTask * BINS_PS * sizeof(double));

  MPI_Allgather(CountModes, BINS_PS * sizeof(long long), MPI_BYTE,
		countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      CountModes[i] = 0;
      for(n = 0; n < NTask; n++)
	CountModes[i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPower, BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[i] = 0;
      for(n = 0; n < NTask; n++)
	SumPower[i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);

  for(i = 0; i < BINS_PS; i++)
    {
      Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[i] > 0)
	Power[i] = SumPower[i] / CountModes[i];
      else
	Power[i] = 0;
    }
}



void powerspec_turb_save(char *fname, double *disp)
{
  FILE *fd;
  char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
  int i;
  if(ThisTask == 0)
    {
      if(!(fd = fopen(fname, "w"))) {snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "can't open file `%s`\n", fname); terminate(buf);}
      fprintf(fd, "%.16g\n", All.Time);
      i = TURB_DRIVING_SPECTRUMGRID;
      fprintf(fd, "%d\n", i);
      i = BINS_PS;
      fprintf(fd, "%d\n", i);
      fprintf(fd, "%g\n", disp[0]); fprintf(fd, "%g\n", disp[1]); fprintf(fd, "%g\n", disp[2]);
      for(i = 0; i < BINS_PS; i++) {fprintf(fd, "%g %g %g %g\n", Kbin[i], Power[i], (double) CountModes[i], SumPower[i]);}
      fclose(fd);
    }
}



/* this function determines the velocity fields by using the nearest cell's values 
 */ 
double powerspec_turb_obtain_fields(void)
{
  int j, dummy;
  long long ntot, npleft;
  int ndone, ndone_flag, ngrp, sendTask, recvTask, place, nexport, nimport, iter;

  double tstart = my_second();

  PRINT_STATUS("Start finding nearest gas-particle for mesh-cell centers (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
  large_array_offset i, n, Ncount = ((large_array_offset)nslab_x) * (TURB_DRIVING_SPECTRUMGRID * TURB_DRIVING_SPECTRUMGRID);  /* number of grid points on the local slab */

  powerspec_turb_nearest_distance = (float *) mymalloc("powerspec_turb_nearest_distance", sizeof(float) * Ncount);
  powerspec_turb_nearest_rkern = (float *) mymalloc("powerspec_turb_nearest_rkern", sizeof(float) * Ncount);

  for(n = 0; n < Ncount; n++)
    {
      powerspec_turb_nearest_distance[n] = 1.0e30;
      powerspec_turb_nearest_rkern[n] = All.BoxSize / pow(All.TotN_gas, 1.0/3);
    }

    /* Modern path: replicates the legacy GLOBAL nearest-gas search semantic
     * (each FFT slab cell gets the globally-nearest gas particle's data, found
     * across ALL ranks — not just the slab-owning rank's local domain).
     *
     * Algorithm per iteration:
     *   1. Each rank packs its still-unresolved slab cells into a pending list.
     *   2. MPI_Allgatherv pending lists -> every rank holds the global pending
     *      list (with each cell tagged by its slab_owner_rank).
     *   3. Each rank builds a local NL with the global pending positions as
     *      arbitrary-source overrides (uses gpu_ngb_list_build's
     *      source_positions_host extension). Walks NL and considers only
     *      home gas (j < num_local) so each cell has at most one candidate
     *      per rank, with the candidate residing as HOME on that rank.
     *   4. MPI_Allreduce(MIN_LOC) on (dist, ThisTask) per pending cell ->
     *      every rank knows the winning rank for each cell.
     *   5. Winning ranks pack full field data (Vel, Vorticity, SmoothedVel,
     *      RandomValue, Density, Velocity_bar/Hat, DuDt_diss) into Alltoallv
     *      send buffers keyed by slab_owner.
     *   6. Slab owners receive data, write velfield[ip] / vorticityfield[ip]
     *      / etc., and mark cells resolved.
     *   7. Cells still unresolved on slab owners get rkern *= 2; iterate.
     *
     * Matches legacy semantic exactly: global nearest, iterative MAXITER
     * radius doubling, terminate on non-convergence. No descope: a slab cell
     * far from this rank's particle domain will be resolved by whichever rank
     * has nearby home gas via Allreduce.
     *
     * Memory peak: O(N_total_pending) cells held briefly during Allgatherv +
     * Alltoallv. For typical SPECTRUMGRID <= 128 this is <100 MB per rank;
     * for 256+ chunking would be needed. */
    {
        int max_iter = MAXITER;
        int iter_count = 0;
        int ghost_imported = 0;
        {
            int need_import_local = (ghost_get_num_ghosts() <= 0) ? 1 : 0;
            int need_import = 0;
            MPI_Allreduce(&need_import_local, &need_import, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if(need_import) {
                if(ghost_get_num_ghosts() > 0) ghost_exchange_cleanup();
                gizmo_density_prep_ghosts(gizmo_ghost_safety_factor());
                ghost_imported = 1;
            }
        }
        int num_local = ghost_get_num_local();
        if(num_local <= 0) num_local = NumPart;
        int num_all = num_local + ghost_get_num_ghosts();
        if(num_all <= 0) num_all = NumPart;

        struct WinnerData {
            double Vel[3];
            double SmoothedVel[3];
            double Vorticity[3];
            double Density;
            double DuDt_diss;
            double RandomValue;
#ifdef TURB_DIFF_DYNAMIC
            double Velocity_bar[3];
            double Velocity_hat[3];
#endif
        };
        struct AlltoallvSend {
            large_array_offset slab_cell_n;     /* index into slab owner's local nearest_distance[] */
            double dist;
            WinnerData data;
        };

        long long npleft;
        long long ntot;
        do {
            /* Stage 1: pack THIS rank's pending slab cells */
            int my_pending = 0;
            for(large_array_offset nn = 0; nn < Ncount; nn++) {
                if(powerspec_turb_nearest_distance[nn] > 1.0e29) my_pending++;
            }

            /* Stage 2: gather pending counts globally */
            std::vector<int> all_counts(NTask, 0);
            MPI_Allgather(&my_pending, 1, MPI_INT, all_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
            std::vector<int> all_offsets(NTask + 1, 0);
            for(int t = 0; t < NTask; t++) all_offsets[t+1] = all_offsets[t] + all_counts[t];
            int total_pending = all_offsets[NTask];

            if(total_pending == 0) { npleft = 0; ntot = 0; break; }

            /* Stage 3: build my pending entries (positions, radii, slab-cell-id) */
            std::vector<double> my_pos(my_pending * 3);
            std::vector<double> my_rad(my_pending);
            std::vector<int>    my_slabn(my_pending); /* index into THIS rank's nearest_distance[] */
            int p_at = 0;
            for(large_array_offset nn = 0; nn < Ncount; nn++) {
                if(powerspec_turb_nearest_distance[nn] > 1.0e29) {
                    int xx = nn / (TURB_DRIVING_SPECTRUMGRID * TURB_DRIVING_SPECTRUMGRID);
                    int yy = (nn - (large_array_offset)xx * TURB_DRIVING_SPECTRUMGRID * TURB_DRIVING_SPECTRUMGRID) / TURB_DRIVING_SPECTRUMGRID;
                    int zz = (nn - (large_array_offset)xx * TURB_DRIVING_SPECTRUMGRID * TURB_DRIVING_SPECTRUMGRID
                                 - (large_array_offset)yy * TURB_DRIVING_SPECTRUMGRID);
                    int xx_glob = xx + slabstart_x;
                    my_pos[3*p_at+0] = (xx_glob + 0.5) / TURB_DRIVING_SPECTRUMGRID * boxSize_X;
                    my_pos[3*p_at+1] = (yy      + 0.5) / TURB_DRIVING_SPECTRUMGRID * boxSize_Y;
                    my_pos[3*p_at+2] = (zz      + 0.5) / TURB_DRIVING_SPECTRUMGRID * boxSize_Z;
                    my_rad[p_at]     = powerspec_turb_nearest_rkern[nn];
                    my_slabn[p_at]   = (int)nn;
                    p_at++;
                }
            }

            /* Stage 4: Allgatherv the pending lists */
            std::vector<double> g_pos(total_pending * 3);
            std::vector<double> g_rad(total_pending);
            std::vector<int>    g_slabn(total_pending);
            std::vector<int>    g_owner(total_pending);
            std::vector<int> rdisp_pos(NTask), rcnt_pos(NTask);
            std::vector<int> rdisp_d  (NTask), rcnt_d  (NTask);
            for(int t = 0; t < NTask; t++) {
                rcnt_d  [t] = all_counts[t];
                rdisp_d [t] = all_offsets[t];
                rcnt_pos[t] = all_counts[t] * 3;
                rdisp_pos[t] = all_offsets[t] * 3;
            }
            MPI_Allgatherv(my_pos.data(), my_pending * 3, MPI_DOUBLE,
                           g_pos.data(), rcnt_pos.data(), rdisp_pos.data(), MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(my_rad.data(), my_pending, MPI_DOUBLE,
                           g_rad.data(), rcnt_d.data(), rdisp_d.data(), MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgatherv(my_slabn.data(), my_pending, MPI_INT,
                           g_slabn.data(), rcnt_d.data(), rdisp_d.data(), MPI_INT, MPI_COMM_WORLD);
            for(int t = 0; t < NTask; t++) {
                for(int k = 0; k < all_counts[t]; k++) g_owner[all_offsets[t] + k] = t;
            }

            /* Stage 5: build NL on this rank with global pending positions as override */
            gpu_neighbor_list_t gnl = {};
            std::vector<int> gnl_neighbors_host;
            std::vector<int> sentinel_active(total_pending > 0 ? total_pending : 1, 0);
            if(total_pending > 0) {
                gpu_particles_arena_acquire(num_all, P, CellP);
                struct particle_data *P_gpu = gpu_particles_arena_P();
                gpu_ngb_list_build(P_gpu, num_all,
                                   sentinel_active.data(), total_pending,
                                   NGB_SEARCH_ONEWAY, 1 /* gas only */,
                                   &gnl, NULL, 1.0,
                                   g_rad.data(),
                                   g_pos.data() /* arbitrary-source override */);
                /* gnl.neighbors is DEVICE_SPACE; host loop below indexes it. */
                if(gnl.total_pairs > 0) {
                    gnl_neighbors_host.resize((size_t)gnl.total_pairs);
                    gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
                }
            }
            const int *gnl_neighbors = gnl_neighbors_host.empty() ? NULL : gnl_neighbors_host.data();

            /* Stage 6: per global pending cell, compute THIS rank's home-gas best */
            std::vector<double> local_dist(total_pending, 1.0e30);
            std::vector<int>    local_idx (total_pending, -1);
            for(int aa = 0; aa < total_pending; aa++) {
                double sx = g_pos[3*aa+0], sy = g_pos[3*aa+1], sz = g_pos[3*aa+2];
                int64_t n_off = gnl.offsets[aa], n_off_end = gnl.offsets[aa+1];
                int best_index = -1;
                double best_r2 = MAX_REAL_NUMBER;
                MyDouble xtmp = 0;
                for(int64_t kk = n_off; kk < n_off_end; kk++) {
                    int p_idx = gnl_neighbors[kk];
                    if(p_idx >= num_local) continue; /* skip ghosts: each rank only contributes its HOME gas as candidate, so MIN_LOC tiebreaks cleanly */
                    if(p_idx >= N_gas) continue;
                    if(P[p_idx].Type != 0 || P[p_idx].Mass <= 0) continue;
                    double dx_raw = P[p_idx].Pos[0] - sx;
                    double dy_raw = P[p_idx].Pos[1] - sy;
                    double dz_raw = P[p_idx].Pos[2] - sz;
                    double dx = NGB_PERIODIC_BOX_LONG_X(dx_raw, dy_raw, dz_raw, 1);
                    double dy = NGB_PERIODIC_BOX_LONG_Y(dx_raw, dy_raw, dz_raw, 1);
                    double dz = NGB_PERIODIC_BOX_LONG_Z(dx_raw, dy_raw, dz_raw, 1);
                    double r2 = dx*dx + dy*dy + dz*dz;
                    if(r2 < best_r2) { best_r2 = r2; best_index = p_idx; }
                }
                if(best_index >= 0) { local_dist[aa] = sqrt(best_r2); local_idx[aa] = best_index; }
            }

            if(total_pending > 0) {
                gpu_ngb_list_free(&gnl, NULL);
                gpu_particles_arena_invalidate();
            }

            /* Stage 7: MPI_Allreduce(MIN_LOC) per cell -> (global_dist, winner_rank) */
            struct DRpair { double d; int r; };
            std::vector<DRpair> local_dr(total_pending), global_dr(total_pending);
            for(int aa = 0; aa < total_pending; aa++) {
                local_dr[aa].d = local_dist[aa];
                /* sentinel rank for "no contribution" must be larger than any real rank
                 * so MIN_LOC tiebreak picks any real contributor; use NTask. */
                local_dr[aa].r = (local_idx[aa] >= 0) ? ThisTask : NTask;
            }
            MPI_Allreduce(local_dr.data(), global_dr.data(), total_pending, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

            /* Stage 8: For cells where I won, pack data per slab_owner */
            std::vector<int> send_count(NTask, 0);
            for(int aa = 0; aa < total_pending; aa++) {
                if(global_dr[aa].r == ThisTask && local_idx[aa] >= 0) send_count[g_owner[aa]]++;
            }
            std::vector<int> send_disp(NTask, 0);
            int total_send = 0;
            for(int t = 0; t < NTask; t++) { send_disp[t] = total_send; total_send += send_count[t]; }

            std::vector<AlltoallvSend> send_buf(total_send);
            std::vector<int> cursor(send_disp);
            for(int aa = 0; aa < total_pending; aa++) {
                if(global_dr[aa].r != ThisTask || local_idx[aa] < 0) continue;
                int idx = local_idx[aa];
                int owner = g_owner[aa];
                AlltoallvSend &e = send_buf[cursor[owner]++];
                e.slab_cell_n = (large_array_offset)g_slabn[aa];
                e.dist = local_dist[aa];
                e.data.Vel[0] = P[idx].Vel[0]; e.data.Vel[1] = P[idx].Vel[1]; e.data.Vel[2] = P[idx].Vel[2];
                e.data.SmoothedVel[0] = CellP[idx].SmoothedVel[0];
                e.data.SmoothedVel[1] = CellP[idx].SmoothedVel[1];
                e.data.SmoothedVel[2] = CellP[idx].SmoothedVel[2];
                e.data.Vorticity[0] = CellP[idx].Vorticity[0];
                e.data.Vorticity[1] = CellP[idx].Vorticity[1];
                e.data.Vorticity[2] = CellP[idx].Vorticity[2];
                e.data.Density = CellP[idx].Density;
                e.data.DuDt_diss = CellP[idx].DuDt_diss;
                e.data.RandomValue = RandomValue[idx];
#ifdef TURB_DIFF_DYNAMIC
                e.data.Velocity_bar[0] = CellP[idx].Velocity_bar[0];
                e.data.Velocity_bar[1] = CellP[idx].Velocity_bar[1];
                e.data.Velocity_bar[2] = CellP[idx].Velocity_bar[2];
                e.data.Velocity_hat[0] = CellP[idx].Velocity_hat[0];
                e.data.Velocity_hat[1] = CellP[idx].Velocity_hat[1];
                e.data.Velocity_hat[2] = CellP[idx].Velocity_hat[2];
#endif
            }

            std::vector<int> recv_count(NTask, 0);
            MPI_Alltoall(send_count.data(), 1, MPI_INT, recv_count.data(), 1, MPI_INT, MPI_COMM_WORLD);
            std::vector<int> recv_disp(NTask, 0);
            int total_recv = 0;
            for(int t = 0; t < NTask; t++) { recv_disp[t] = total_recv; total_recv += recv_count[t]; }

            /* Convert counts/disps to bytes for MPI_BYTE Alltoallv */
            std::vector<int> sc_b(NTask), sd_b(NTask), rc_b(NTask), rd_b(NTask);
            for(int t = 0; t < NTask; t++) {
                sc_b[t] = send_count[t] * sizeof(AlltoallvSend);
                sd_b[t] = send_disp [t] * sizeof(AlltoallvSend);
                rc_b[t] = recv_count[t] * sizeof(AlltoallvSend);
                rd_b[t] = recv_disp [t] * sizeof(AlltoallvSend);
            }
            std::vector<AlltoallvSend> recv_buf(total_recv > 0 ? total_recv : 1);
            MPI_Alltoallv(send_buf.data(), sc_b.data(), sd_b.data(), MPI_BYTE,
                          recv_buf.data(), rc_b.data(), rd_b.data(), MPI_BYTE, MPI_COMM_WORLD);

            /* Stage 9: unpack received data into MY slab's velfield arrays */
            for(int k = 0; k < total_recv; k++) {
                AlltoallvSend &e = recv_buf[k];
                large_array_offset n_target = e.slab_cell_n;
                int ii = n_target / (TURB_DRIVING_SPECTRUMGRID * TURB_DRIVING_SPECTRUMGRID);
                int jj = (n_target - (large_array_offset)ii * TURB_DRIVING_SPECTRUMGRID * TURB_DRIVING_SPECTRUMGRID) / TURB_DRIVING_SPECTRUMGRID;
                int kk = (n_target - (large_array_offset)ii * TURB_DRIVING_SPECTRUMGRID * TURB_DRIVING_SPECTRUMGRID
                                   - (large_array_offset)jj * TURB_DRIVING_SPECTRUMGRID);
                int ip = TURB_DRIVING_SPECTRUMGRID2 * (TURB_DRIVING_SPECTRUMGRID * ii + jj) + kk;
                powerspec_turb_nearest_distance[n_target] = e.dist;
                velfield[0][ip] = e.data.Vel[0]; velfield[1][ip] = e.data.Vel[1]; velfield[2][ip] = e.data.Vel[2];
#ifdef TURB_DIFF_DYNAMIC
                velbarfield[0][ip] = e.data.Velocity_bar[0];
                velbarfield[1][ip] = e.data.Velocity_bar[1];
                velbarfield[2][ip] = e.data.Velocity_bar[2];
                velhatfield[0][ip] = e.data.Velocity_hat[0];
                velhatfield[1][ip] = e.data.Velocity_hat[1];
                velhatfield[2][ip] = e.data.Velocity_hat[2];
#endif
                smoothedvelfield[0][ip] = e.data.SmoothedVel[0];
                smoothedvelfield[1][ip] = e.data.SmoothedVel[1];
                smoothedvelfield[2][ip] = e.data.SmoothedVel[2];
                velrhofield[0][ip] = sqrt(e.data.Density) * e.data.Vel[0];
                velrhofield[1][ip] = sqrt(e.data.Density) * e.data.Vel[1];
                velrhofield[2][ip] = sqrt(e.data.Density) * e.data.Vel[2];
                vorticityfield[0][ip] = e.data.Vorticity[0];
                vorticityfield[1][ip] = e.data.Vorticity[1];
                vorticityfield[2][ip] = e.data.Vorticity[2];
                if(e.data.DuDt_diss >= 0) { dis1field[ip] = sqrt(e.data.DuDt_diss); dis2field[ip] = 0; }
                else                      { dis1field[ip] = 0; dis2field[ip] = sqrt(-e.data.DuDt_diss); }
                randomfield[ip] = e.data.RandomValue;
                densityfield[ip] = e.data.Density;
            }

            /* Stage 10: bookkeeping for next iter */
            npleft = 0;
            for(large_array_offset nnn = 0; nnn < Ncount; nnn++) {
                if(powerspec_turb_nearest_distance[nnn] > 1.0e29) {
                    npleft++;
                    powerspec_turb_nearest_rkern[nnn] *= 2.0;
                } else {
                    powerspec_turb_nearest_distance[nnn] = 0; /* skip on subsequent iters */
                }
            }
            sumup_longs(1, &npleft, &ntot);
            if(ntot > 0) {
                iter_count++;
                if(iter_count > 0 && ThisTask == 0) PRINT_STATUS("powespec_vel nearest iteration %d: need to repeat for %lld particles", iter_count, ntot);
                if(iter_count > max_iter) terminate("failed to converge");
            }
        } while(ntot > 0);

        if(ghost_imported) ghost_exchange_cleanup();
    }

  myfree(powerspec_turb_nearest_rkern);
  myfree(powerspec_turb_nearest_distance);

    if(ThisTask == 0) {printf("done finding velocity field\n");}

  double tend = my_second();
  return timediff(tstart, tend);
}


void powerspec_turb_calc_dispersion(void)
{
  int dim, i, j, k;

  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i=0; i < nslab_x;i++)
	for(j=0; j< TURB_DRIVING_SPECTRUMGRID; j++)
	  for(k=0; k< TURB_DRIVING_SPECTRUMGRID; k++)
	    {
	      int ip = TURB_DRIVING_SPECTRUMGRID2 * (TURB_DRIVING_SPECTRUMGRID * i + j) + k;

	      vsum += velfield[dim][ip];
	    }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(TURB_DRIVING_SPECTRUMGRID, 3);
      
      for(i=0; i < nslab_x;i++)
	for(j=0; j< TURB_DRIVING_SPECTRUMGRID; j++)
	  for(k=0; k< TURB_DRIVING_SPECTRUMGRID; k++)
	    {
	      int ip = TURB_DRIVING_SPECTRUMGRID2 * (TURB_DRIVING_SPECTRUMGRID * i + j) + k;
	      
	      velfield[dim][ip] -= vmean;
	    }

      for(i=0; i < nslab_x;i++)
	for(j=0; j< TURB_DRIVING_SPECTRUMGRID; j++)
	  for(k=0; k< TURB_DRIVING_SPECTRUMGRID; k++)
	    {
	      int ip = TURB_DRIVING_SPECTRUMGRID2 * (TURB_DRIVING_SPECTRUMGRID * i + j) + k;
	      
	      vdisp += velfield[dim][ip] * velfield[dim][ip];
	    }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      vel_disp[dim] = vdisp_all / pow(TURB_DRIVING_SPECTRUMGRID, 3);      
    }

#ifdef TURB_DIFF_DYNAMIC
  /* velbar */
  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i=0; i < nslab_x;i++)
  for(j=0; j< POWERSPEC_GRID; j++)
    for(k=0; k< POWERSPEC_GRID; k++)
      {
        int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;
        
        vsum += velbarfield[dim][ip];
      }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(POWERSPEC_GRID, 3);

      for(i=0; i < nslab_x;i++)
  for(j=0; j< POWERSPEC_GRID; j++)
    for(k=0; k< POWERSPEC_GRID; k++)
      {
        int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

        velbarfield[dim][ip] -= vmean;
      }

      for(i=0; i < nslab_x;i++)
  for(j=0; j< POWERSPEC_GRID; j++)
    for(k=0; k< POWERSPEC_GRID; k++)
      {
        int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

        vdisp += velbarfield[dim][ip] * velbarfield[dim][ip];
      }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      velbar_disp[dim] = vdisp_all / pow(POWERSPEC_GRID, 3);
    }

  /* velhat */
  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i=0; i < nslab_x;i++)
  for(j=0; j< POWERSPEC_GRID; j++)
    for(k=0; k< POWERSPEC_GRID; k++)
      {
        int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

        vsum += velhatfield[dim][ip];
      }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(POWERSPEC_GRID, 3);

      for(i=0; i < nslab_x;i++)
  for(j=0; j< POWERSPEC_GRID; j++)
    for(k=0; k< POWERSPEC_GRID; k++)
      {
        int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

        velhatfield[dim][ip] -= vmean;
      }

      for(i=0; i < nslab_x;i++)
  for(j=0; j< POWERSPEC_GRID; j++)
    for(k=0; k< POWERSPEC_GRID; k++)
      {
        int ip = POWERSPEC_GRID2 * (POWERSPEC_GRID * i + j) + k;

        vdisp += velhatfield[dim][ip] * velhatfield[dim][ip];
      }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      velhat_disp[dim] = vdisp_all / pow(POWERSPEC_GRID, 3);
    }
#endif
  for(dim = 0; dim < 3; dim++)
    {
      double vsum = 0, vsum_all, vmean, vdisp = 0, vdisp_all;

      for(i=0; i < nslab_x;i++)
	for(j=0; j< TURB_DRIVING_SPECTRUMGRID; j++)
	  for(k=0; k< TURB_DRIVING_SPECTRUMGRID; k++)
	    {
	      int ip = TURB_DRIVING_SPECTRUMGRID2 * (TURB_DRIVING_SPECTRUMGRID * i + j) + k;

	      vsum += velrhofield[dim][ip];
	    }

      MPI_Allreduce(&vsum, &vsum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      vmean = vsum_all / pow(TURB_DRIVING_SPECTRUMGRID, 3);
      
      for(i=0; i < nslab_x;i++)
	for(j=0; j< TURB_DRIVING_SPECTRUMGRID; j++)
	  for(k=0; k< TURB_DRIVING_SPECTRUMGRID; k++)
	    {
	      int ip = TURB_DRIVING_SPECTRUMGRID2 * (TURB_DRIVING_SPECTRUMGRID * i + j) + k;
	      
	      velrhofield[dim][ip] -= vmean;
	    }

      for(i=0; i < nslab_x;i++)
	for(j=0; j< TURB_DRIVING_SPECTRUMGRID; j++)
	  for(k=0; k< TURB_DRIVING_SPECTRUMGRID; k++)
	    {
	      int ip = TURB_DRIVING_SPECTRUMGRID2 * (TURB_DRIVING_SPECTRUMGRID * i + j) + k;
	      
	      vdisp += velrhofield[dim][ip] * velrhofield[dim][ip];
	    }

      MPI_Allreduce(&vdisp, &vdisp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      velrho_disp[dim] = vdisp_all / pow(TURB_DRIVING_SPECTRUMGRID, 3);      
    }
}







#endif



