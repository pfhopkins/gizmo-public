#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"


/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO. The modifications
 * mostly center on added functionality for new modules, changed variables for 
 * cosmology, and consolidating the relevant functions into a single file.
 */



/* These are critical factors used throughout for co-moving integrations. Set them here and
   call THESE, rather than trying to come up with the factors throughout, since that makes debugging a nightmare */
void set_cosmo_factors_for_current_time(void)
{
    if(All.ComovingIntegrationOn)
    {
        /* All.cf_atime = a = 1/(1+z), the cosmological scale factor */
        All.cf_atime = All.Time;
        /* All.cf_a2inv is just handy */
        All.cf_a2inv = 1 / (All.Time * All.Time);
        /* All.cf_a3inv * Density_code = Density_physical */
        All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
        /* Pressure_code/Density_code = All.cf_afac1 * Pressure_physical/Density_physical */
         All.cf_afac1 = 1; // -- no longer explicitly needed -- retained only for historical compatibility reasons with users importing from gadget/arepo
        /* All.cf_afac2 * Pressure_code/Density_code * 1/r_code = Pressure_physical/Density_physical * 1/r_physical */
        All.cf_afac2 = 1 / (All.Time * All.cf_afac1);
        /* All.cf_afac3 * sqrt(Pressure_code/Density_code) = sqrt(Pressure_phys/Density_phys) = cs_physical */
        All.cf_afac3 = 1 / sqrt(All.cf_afac1); // -- no longer explicitly needed -- retained only for historical compatibility reasons with users importing from gadget/arepo
        /* time units: proper time dt_phys = 1/hubble_function(a) * dz/(1+z) = dlna / hubble_function(a)
            code time unit in comoving is dlna, so dt_phys = dt_code / All.cf_hubble_a   */
        All.cf_hubble_a = hubble_function(All.Time); /* hubble_function(a) = H(a) = H(z) */
        /* dt_code * v_code/r_code = All.cf_hubble_a2 * dt_phys * v_phys/r_phys */
        All.cf_hubble_a2 = All.Time * All.Time * hubble_function(All.Time);
        /* set custom tabulated values */
#ifdef GR_TABULATED_COSMOLOGY_G
        All.G = All.Gini * dGfak(All.Time);
#endif
    }
    else {All.cf_atime = All.cf_a2inv = All.cf_a3inv = All.cf_afac1 = All.cf_afac2 = All.cf_afac3 = All.cf_hubble_a = All.cf_hubble_a2 = 1;}
#ifdef CHIMES
    ChimesGlobalVars.cmb_temperature = (ChimesFloat) (2.725 / All.cf_atime);
#endif
}




/* this function gets called regardless of the cosmology choices: 
    anything which modifies the growth history should live here. 
    This is the usual hubble function H0 * E(z); so for example 
    the proper time:
        dt = 1/hubble_function(a) * dz/(1+z) = dlna / hubble_function(a)
*/
double INLINE_FUNC hubble_function(double a)
{
    double hubble_a;
    
#ifdef GR_TABULATED_COSMOLOGY_H
    hubble_a = All.Hubble_H0_CodeUnits * hubble_function_external(a);
#else
    hubble_a = All.OmegaRadiation / (a*a*a*a) + All.OmegaMatter / (a*a*a) + (1 - All.OmegaMatter - All.OmegaLambda - All.OmegaRadiation) / (a*a)
#ifdef GR_TABULATED_COSMOLOGY
    + DarkEnergy_a(a);
#else
    + All.OmegaLambda;
#endif
    hubble_a = All.Hubble_H0_CodeUnits * sqrt(hubble_a);
#endif
#ifdef GR_TABULATED_COSMOLOGY_G
    hubble_a *= dHfak(a);
#endif
    return (hubble_a);
}



#ifdef GR_TABULATED_COSMOLOGY
#if defined(GR_TABULATED_COSMOLOGY_W) || defined(GR_TABULATED_COSMOLOGY_G) || defined(GR_TABULATED_COSMOLOGY_H)

#define ANZ_W_A_IN 4000
#define ANZ_W_A 10000

static MyFloat atab[ANZ_W_A_IN];
static MyFloat wtab[ANZ_W_A_IN];
static MyFloat intwtab[ANZ_W_A + 1];
static MyFloat intwatab[ANZ_W_A + 1];

#ifdef GR_TABULATED_COSMOLOGY_G
static MyFloat dHtab[ANZ_W_A_IN];
static MyFloat dGtab[ANZ_W_A_IN];
static MyFloat intdHtab[ANZ_W_A + 1];
static MyFloat intdGtab[ANZ_W_A + 1];
#endif

#ifdef GR_TABULATED_COSMOLOGY_H
static MyFloat Htab[ANZ_W_A_IN];
static MyFloat intHtab[ANZ_W_A + 1];
#endif

/* set-up table with Exp(-Integral_a^1 (q+w(a'))/a' da'
 * needed in hubble function for time dependent w
 */

void fwa_init(void)
{
  int count = 0, i;
  char buf[200], buf1[200], buf2[200], buf3[200], buf4[200], buf5[200];
  FILE *fd;
  MyFloat a_first, w_first, a, w, sum;


  if((fd = fopen(All.TabulatedCosmologyFile, "r")))
    {
      if(ThisTask == 0)
	{
	  printf("\nreading w of a from file `%s'\n", All.TabulatedCosmologyFile);
	}
      atab[0] = -1.0;		/* we have to extrapolate wtab[0] later ! */
      count = 1;
      while(!feof(fd) && count < ANZ_W_A_IN)
	{
	  if(fgets(buf, 200, fd))
	    {
	      if(sscanf(buf, "%s%s%s%s%s", buf1, buf2, buf3, buf4, buf5) < 5)
		{
		  if(ThisTask == 0)
		    {
		      printf("Wrong syntax in file `%s', line %d\n", All.TabulatedCosmologyFile, count);
		      fflush(stdout);
		    }
		  endrun(0);
		}
	      a = atof(buf1); /* column is scale factor */
	      if(a == 0.0 && count == 1) count--; /* w(0) present in file, so fill the first element ! */
	      atab[count] = a;
	      wtab[count] = atof(buf2);
#ifdef GR_TABULATED_COSMOLOGY_H
          Htab[count] = atof(buf3);
#endif
#ifdef GR_TABULATED_COSMOLOGY_G
	      dGtab[count] = atof(buf4);
          dHtab[count] = atof(buf5);
#endif
	      count++;
	    }
	}
      fclose(fd);
      if(count >= ANZ_W_A_IN - 1)
	{
	  if(ThisTask == 0)
	    {
	      printf("File `%s' contains to many datapoints, increase ANZ_W_A_IN !\n", All.TabulatedCosmologyFile);
	      fflush(stdout);
	    }
	  endrun(0);
	}
      if(count <= 2)
	{
	  if(ThisTask == 0)
	    {
	      printf("File `%s' has to less Data Points (%d) !\n", All.TabulatedCosmologyFile, count);
	      fflush(stdout);
	    }
	  endrun(0);
	}

      if(atab[0] < 0.)		/* We still have to extrapolate w to a = 0 (w[0]) */
	{
	  atab[0] = 0.;
	  wtab[0] = wtab[1] - (wtab[2] - wtab[1]) / (atab[2] - atab[1]) * (atab[1] - atab[0]);
#ifdef GR_TABULATED_COSMOLOGY_G
	  dHtab[0] = dHtab[1] - (dHtab[2] - dHtab[1]) / (atab[2] - atab[1]) * (atab[1] - atab[0]);
	  dGtab[0] = dGtab[1] - (dGtab[2] - dGtab[1]) / (atab[2] - atab[1]) * (atab[1] - atab[0]);
#endif
#ifdef GR_TABULATED_COSMOLOGY_H
	  Htab[0] = Htab[1] - (Htab[2] - Htab[1]) / (atab[2] - atab[1]) * (atab[1] - atab[0]);
#endif
	}

/* Setp back if tables go bejond a=1 */
      while(atab[count - 1] > 1.0)
	count--;

/* Calculate w(1) if needed */
      if(atab[count - 1] < 1.)
	{
	  atab[count] = 1.0;
	  wtab[count] = wtab[count - 1] + (wtab[count - 1] - wtab[count - 2]) / (atab[count - 1] - atab[count - 2]) * (1. - atab[count - 1]);
#ifdef GR_TABULATED_COSMOLOGY_G
	  dHtab[count] = dHtab[count - 1] + (dHtab[count - 1] - dHtab[count - 2]) / (atab[count - 1] - atab[count - 2]) * (1. - atab[count - 1]);
	  dGtab[count] = dGtab[count - 1] + (dGtab[count - 1] - dGtab[count - 2]) / (atab[count - 1] - atab[count - 2]) * (1. - atab[count - 1]);
#endif
#ifdef GR_TABULATED_COSMOLOGY_H
	  Htab[count] = Htab[count - 1] + (Htab[count - 1] - Htab[count - 2]) / (atab[count - 1] - atab[count - 2]) * (1. - atab[count - 1]);
#endif
/*            if(ThisTask ==0) 
              {
                printf("%d %f %f %f\n",count,atab[count-2],atab[count-1],atab[count]);
                printf("%d %f %f %f\n",count,wtab[count-2],wtab[count-1],wtab[count]);
              }*/
	  count++;
	}

/* Now calculated the integral (starting from last to first to save Time !
 * Explicit asume that a[0]=0. and a[count-1]=1. , which is enshured by   
 * the loading precedure !                                                  */

/* Set todays values in the tables */
      intwtab[ANZ_W_A] = All.OmegaLambda;
      intwatab[ANZ_W_A] = wtab[count - 1];
#ifdef GR_TABULATED_COSMOLOGY_G
      intdHtab[ANZ_W_A] = dHtab[count - 1];
      intdGtab[ANZ_W_A] = dGtab[count - 1];
#endif
#ifdef GR_TABULATED_COSMOLOGY_H
      intHtab[ANZ_W_A] = Htab[count - 1];
#endif

/* Place count on last entry in table */
      count--;

      a_first = atab[count];	/* Startinv value should be 1.0 ! */
      w_first = wtab[count];	/* Starting value from table ! */
      sum = 0.0;		/* set int to 0.0 */
      for(i = ANZ_W_A - 1; i >= 1; i--)
	{
	  a = (MyFloat) i / (MyFloat) ANZ_W_A;
	  if(count > 1)		/* Still inside the table */
	    {
	      while(atab[count - 1] > a && count > 0)
		{
		  sum += 0.5 * ((1. + w_first) / a_first + (1. + wtab[count - 1]) / atab[count - 1])
		    * (a_first - atab[count - 1]);
		  count--;
		  a_first = atab[count];
		  w_first = wtab[count];
		}
	      w = w_first - (wtab[count] - wtab[count - 1]) / (atab[count] - atab[count - 1]) * (a_first - a);
	      sum += 0.5 * ((1. + w_first) / a_first + (1. + w) / a) * (a_first - a);
	      w_first = w;
	      a_first = a;
	    }
	  else
	    {
	      w = w_first - (wtab[count] - wtab[count - 1]) / (atab[count] - atab[count - 1]) * (a_first - a);
	      sum += 0.5 * ((1. + w_first) / a_first + (1. + w) / a) * (a_first - a);
	      w_first = w;
	      a_first = a;
	    }
	  intwtab[i] = All.OmegaLambda * exp(3. * sum);
	  intwatab[i] = wtab[count - 1] + (wtab[count] - wtab[count - 1]) /
	    (atab[count] - atab[count - 1]) * (a - atab[count - 1]);
#ifdef GR_TABULATED_COSMOLOGY_G
	  intdHtab[i] = dHtab[count - 1] + (dHtab[count] - dHtab[count - 1]) / (atab[count] - atab[count - 1]) * (a - atab[count - 1]);
	  intdGtab[i] = dGtab[count - 1] + (dGtab[count] - dGtab[count - 1]) / (atab[count] - atab[count - 1]) * (a - atab[count - 1]);
#endif
#ifdef GR_TABULATED_COSMOLOGY_H
	  intHtab[i] = Htab[count - 1] + (Htab[count] - Htab[count - 1]) / (atab[count] - atab[count - 1]) * (a - atab[count - 1]);
#endif
	}
      /* artificially define value for a=0 */
      intwtab[0] = intwtab[1];
      intwatab[0] = intwatab[1];
#ifdef GR_TABULATED_COSMOLOGY_G
      intdHtab[0] = intdHtab[1];
      intdGtab[0] = intdGtab[1];
#endif
#ifdef GR_TABULATED_COSMOLOGY_H
      intHtab[0] = intHtab[1];
#endif
    }
  else
    {
      if(ThisTask == 0)
	{
	  printf("\nFile `%s' not found !\n", All.TabulatedCosmologyFile);
	  fflush(stdout);
	}
      endrun(0);
    }

  if(ThisTask == 0)
    {
      printf("Integrating w(a) finisched.\n");
    }
}

/* This function the integral w(a) therm for the actual time
 * needed in the hubble function.
 */
double INLINE_FUNC fwa(double a)
{
  int ai;
  double fwa = 0.0;

  ai = a * ANZ_W_A;
  if(ai >= ANZ_W_A)
    {
      fwa = intwtab[ANZ_W_A];
    }
  else
    {
      fwa = intwtab[ai] + (intwtab[ai + 1] - intwtab[ai]) * (a * (double) ANZ_W_A - (double) ai);
    }
/*   if(ThisTask==0) printf("%f %f %f %f\n",a,intwtab[ai],fwa,intwtab[ai+1]);*/
  return (fwa);
}

#endif



double DarkEnergy_a(double a)	/* only needed for comoving integration */
{
#ifdef GR_TABULATED_COSMOLOGY_W
  return fwa(a);
#else
  return (All.OmegaLambda * pow(a, -3. * (1 + All.DarkEnergyConstantW)));
#endif
}


double DarkEnergy_t(double Time)	/* only needed for physical integration */
{
  return All.DarkEnergyConstantW;
}

#ifdef GR_TABULATED_COSMOLOGY_W

/* This function returns the interpolated equation of state parameter.
This is only used for information in the log files.
 */
double INLINE_FUNC get_wa(double a)
{
  int ai;
  double fw = 0.0;

  ai = a * ANZ_W_A;
  if(ai >= ANZ_W_A)
    {
      fw = intwatab[ANZ_W_A];
    }
  else
    {
      fw = intwatab[ai] + (intwatab[ai + 1] - intwatab[ai]) * (a * (double) ANZ_W_A - (double) ai);
    }
  return (fw);
}

#ifdef GR_TABULATED_COSMOLOGY_G

/* This function returns the interpolated correction for the Hubble function
 */
double INLINE_FUNC dHfak(double a)
{
  int ai;
  double fdH = 0.0;

  ai = a * ANZ_W_A;
  if(ai >= ANZ_W_A)
    {
      fdH = intdHtab[ANZ_W_A];
    }
  else
    {
      fdH = intdHtab[ai] + (intdHtab[ai + 1] - intdHtab[ai]) * (a * (double) ANZ_W_A - (double) ai);
    }
  return (fdH);
}

/* This function returns the interpolated correction for the Gravitational constant
 */
double INLINE_FUNC dGfak(double a)
{
  int ai;
  double fdG = 0.0;

  ai = a * ANZ_W_A;
  if(ai >= ANZ_W_A)
    {
      fdG = intdGtab[ANZ_W_A];
    }
  else
    {
      fdG = intdGtab[ai] + (intdGtab[ai + 1] - intdGtab[ai]) * (a * (double) ANZ_W_A - (double) ai);
    }
  return (fdG);
}

#endif

#ifdef GR_TABULATED_COSMOLOGY_H

/* This function returns the interpolated correction for the Hubble function
 */
double INLINE_FUNC hubble_function_external(double a)
{
  int ai;
  double H = 0.0;

  ai = a * ANZ_W_A;
  if(ai >= ANZ_W_A)
    {
      H = intHtab[ANZ_W_A];
    }
  else
    {
      H = intHtab[ai] + (intHtab[ai + 1] - intHtab[ai]) * (a * (double) ANZ_W_A - (double) ai);
    }
  return (H);
}

#endif


#endif
#endif

