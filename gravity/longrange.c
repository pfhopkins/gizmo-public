#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef USE_FFTW3
#include <fftw3-mpi.h>
#include "myfftw3.h"
#endif

/*! \file longrange.c
 *  \brief driver routines for computation of long-range gravitational PM force
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * significantly by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#ifdef PMGRID

/*! Driver routine to call initializiation of periodic or/and non-periodic FFT
 *  routines.
 */
void long_range_init(void)
{
#ifdef USE_FFTW3
  fftw_mpi_init(); 
#endif
#ifdef BOX_PERIODIC
  pm_init_periodic();
#ifdef PM_PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif
#else
  pm_init_nonperiodic();
#endif
}


void long_range_init_regionsize(void)
{
#ifdef BOX_PERIODIC
#ifdef PM_PLACEHIGHRESREGION
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
#else
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
}


/*! This function computes the long-range PM force for all particles.
 */
void long_range_force(void)
{
  int i;

#ifndef BOX_PERIODIC
  int j;
  double fac;
#endif


  for(i = 0; i < NumPart; i++)
    {
      P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;
#ifdef EVALPOTENTIAL
      P[i].PM_Potential = 0;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
      P[i].tidal_tensorpsPM[0][0] = P[i].tidal_tensorpsPM[0][1] = P[i].tidal_tensorpsPM[0][2] = 0;
      P[i].tidal_tensorpsPM[1][0] = P[i].tidal_tensorpsPM[1][1] = P[i].tidal_tensorpsPM[1][2] = 0;
      P[i].tidal_tensorpsPM[2][0] = P[i].tidal_tensorpsPM[2][1] = P[i].tidal_tensorpsPM[2][2] = 0;
#endif
    }

#ifdef SELFGRAVITY_OFF
  return;
#endif


#ifdef BOX_PERIODIC
  pmforce_periodic(0, NULL);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE   /* choose what kind of tidal field calculation you want (for this step use Fourier method; the direct-difference method is buggy still) */
    pmtidaltensor_periodic_fourier(0); pmtidaltensor_periodic_fourier(1); pmtidaltensor_periodic_fourier(2); pmtidaltensor_periodic_fourier(3); pmtidaltensor_periodic_fourier(4); pmtidaltensor_periodic_fourier(5); /* fourier */
    //pmtidaltensor_periodic_diff(); /* finite-difference */
#endif
#ifdef PM_PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE   /* choose what kind of tidal field calculation you want (fourier method is disfavored in current code) */
    //pmtidaltensor_nonperiodic_fourier(1, 0); pmtidaltensor_nonperiodic_fourier(1, 1); pmtidaltensor_nonperiodic_fourier(1, 2); pmtidaltensor_nonperiodic_fourier(1, 3); pmtidaltensor_nonperiodic_fourier(1, 4); pmtidaltensor_nonperiodic_fourier(1, 5); /* fourier */
    pmtidaltensor_nonperiodic_diff(1); /* finite-difference */
#endif
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1);	/* try again */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE   /* choose what kind of tidal field calculation you want (fourier method is disfavored in current code) */
        //pmtidaltensor_nonperiodic_fourier(1, 0); pmtidaltensor_nonperiodic_fourier(1, 1); pmtidaltensor_nonperiodic_fourier(1, 2); pmtidaltensor_nonperiodic_fourier(1, 3); pmtidaltensor_nonperiodic_fourier(1, 4); pmtidaltensor_nonperiodic_fourier(1, 5); /* fourier */
        pmtidaltensor_nonperiodic_diff(1); /* finite-difference */
#endif
    }
  if(i == 1)
    endrun(68686);
#endif
#else
  i = pmforce_nonperiodic(0);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
  //pmtidaltensor_nonperiodic_fourier(0, 0); pmtidaltensor_nonperiodic_fourier(0, 1); pmtidaltensor_nonperiodic_fourier(0, 2); pmtidaltensor_nonperiodic_fourier(0, 3); pmtidaltensor_nonperiodic_fourier(0, 4); pmtidaltensor_nonperiodic_fourier(0, 5);
  pmtidaltensor_nonperiodic_diff(0);
#endif

  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
        pm_init_regionsize();
        pm_setup_nonperiodic_kernel();
        i = pmforce_nonperiodic(0);    /* try again */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        //pmtidaltensor_nonperiodic_fourier(0, 0); pmtidaltensor_nonperiodic_fourier(0, 1); pmtidaltensor_nonperiodic_fourier(0, 2); pmtidaltensor_nonperiodic_fourier(0, 3); pmtidaltensor_nonperiodic_fourier(0, 4); pmtidaltensor_nonperiodic_fourier(0, 5);
        pmtidaltensor_nonperiodic_diff(0);
#endif
    }
  if(i == 1)
    endrun(68687);
#ifdef PM_PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
  //pmtidaltensor_nonperiodic_fourier(1, 0); pmtidaltensor_nonperiodic_fourier(1, 1); pmtidaltensor_nonperiodic_fourier(1, 2); pmtidaltensor_nonperiodic_fourier(1, 3); pmtidaltensor_nonperiodic_fourier(1, 4); pmtidaltensor_nonperiodic_fourier(1, 5);
  pmtidaltensor_nonperiodic_diff(1);
#endif
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

        /* try again */

        for(i = 0; i < NumPart; i++) {P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;}
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        P[i].tidal_tensorpsPM[0][0] = P[i].tidal_tensorpsPM[0][1] = P[i].tidal_tensorpsPM[0][2] = P[i].tidal_tensorpsPM[1][0] = P[i].tidal_tensorpsPM[1][1] = P[i].tidal_tensorpsPM[1][2] = P[i].tidal_tensorpsPM[2][0] = P[i].tidal_tensorpsPM[2][1] = P[i].tidal_tensorpsPM[2][2] = 0;
#endif
        i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        //pmtidaltensor_nonperiodic_fourier(0, 0); pmtidaltensor_nonperiodic_fourier(0, 1); pmtidaltensor_nonperiodic_fourier(0, 2); pmtidaltensor_nonperiodic_fourier(0, 3); pmtidaltensor_nonperiodic_fourier(0, 4); pmtidaltensor_nonperiodic_fourier(0, 5);
        //pmtidaltensor_nonperiodic_fourier(1, 0); pmtidaltensor_nonperiodic_fourier(1, 1); pmtidaltensor_nonperiodic_fourier(1, 2); pmtidaltensor_nonperiodic_fourier(1, 3); pmtidaltensor_nonperiodic_fourier(1, 4); pmtidaltensor_nonperiodic_fourier(1, 5);
        pmtidaltensor_nonperiodic_diff(0); pmtidaltensor_nonperiodic_diff(1); /* two-iterations here */
#endif
    }
  if(i != 0)
    endrun(68688);
#endif
#endif


#ifndef BOX_PERIODIC
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits * All.OmegaMatter;
      for(i = 0; i < NumPart; i++) {for(j = 0; j < 3; j++) {P[i].GravPM[j] += fac * P[i].Pos[j];}}
    }

  if(All.ComovingIntegrationOn == 0) /* special factor as in gravtree for cases where we want to run a non-cosmological simulation but with dark energy terms */
    {
      fac = All.OmegaLambda * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits;
      for(i = 0; i < NumPart; i++) {for(j = 0; j < 3; j++) {P[i].GravPM[j] += fac * P[i].Pos[j];}}
    }
#endif

}


#endif
