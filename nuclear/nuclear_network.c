#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_const_cgs.h>
#include <gsl/gsl_const_num.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*
 *  This code is place-holder, inherited from GADGET3
 */

#ifdef NUCLEAR_NETWORK
//#include "./nuclear_network.h"

/* to use this network you need 5 files:
 * one file containing information about which species your network should use [species.txt]
 * one file containing the partition functions for the species [part.txt]
 *     you can download it e.g. from http://www.nscl.msu.edu/~nero/db/library.php?action=download | http://www.nscl.msu.edu/~nero/db/docs/part_frdm.asc
 * one file containing the reaction rates in the REACLIB 2 format (described in http://www.nscl.msu.edu/~nero/db/docs/reaclibFormat.pdf) [rates.txt]
 *     you can download it e.g. from http://www.nscl.msu.edu/~nero/db/library.php?action=download
 * one file containing the binding energies of the nuclei [masses.txt]
 *     you can download it e.g. from http://www.nndc.bnl.gov/masses/ | http://www.nndc.bnl.gov/masses/mass.mas03
 * one file containing weak reactions e.g. by langanke 2001
 */


static const double conv = 1.602177e-12 * 1.0e3 * 6.0221367e23;	/* eV2erg * 1.0e3 [keV] * avogadro */

void network_normalize(double *x, double *e, const struct network_data *nd, struct network_workspace *nw)
{
    double sum, xnew;
    int i;
    
    sum = 0;
    for(i = 0; i < nd->nuc_count; i++)
    {
        sum += x[i];
    }
    
    if(e)
    {
        for(i = 0; i < nd->nuc_count; i++)
        {
            xnew = x[i] / sum;
            *e -= (xnew - x[i]) * nd->nucdata[i].exm * conv;
            x[i] = xnew;
        }
    }
    else
    {
        for(i = 0; i < nd->nuc_count; i++)
            x[i] /= sum;
    }
}

int network_integrate(double temp, double rho, const double *x, double *dx, double dt, double *dedt,
                      double *drhodt, const struct network_data *nd, struct network_workspace *nw)
{
    double *y;
    double sum;
    int i;
    
    if(dt == 0 || temp < 1e7)
    {
        for(i = 0; i < nd->nuc_count; i++)
            dx[i] = 0;
        *dedt = 0;
        if(drhodt)
            *drhodt = 0;
        return 0;
    }
    
    /* calculate number densities */
    y = nw->y;
#ifdef NETWORK_SEP_YZ
    y[nd->iYz] = 0.0;
#endif
    for(i = 0; i < nd->nuc_count; i++)
    {
        y[i] = x[i] / nd->nucdata[i].na;
#ifdef NETWORK_SEP_YZ
        y[nd->iYz] += y[i] * gsl_pow_2(nd->nucdata[i].nz);
#endif
    }
    
#if NETWORK_VARIABLE
    y[nd->iTemp] = temp;
#endif
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
    y[nd->iRho] = rho;
#endif
    
    /* run network */
    network_solver_integrate(temp, rho, y, dt, nd, nw);
    
    /* normalise */
    sum = 0.0;
    for(i = 0; i < nd->nuc_count; i++)
    {
        if(y[i] > 1.0)
            y[i] = 1.0;
        if(y[i] < 1e-30)
            y[i] = 1e-30;
        sum += y[i] * nd->nucdata[i].na;
    }
    for(i = 0; i < nd->nuc_count; i++)
    {
        y[i] /= sum;
    }
    /* calculate change of mass fractions and energy release */
    *dedt = 0;
    for(i = 0; i < nd->nuc_count; i++)
    {
        dx[i] = (y[i] * nd->nucdata[i].na - x[i]) / dt;
        *dedt -= dx[i] / nd->nucdata[i].na * nd->nucdata[i].exm;
    }
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
    if(drhodt)
        *drhodt = (y[nd->iRho] - rho) / dt;
#else
    if(drhodt)
        *drhodt = 0.0;
#endif
    *dedt *= conv;
    
    return 0;
}


#define ELECTRON_CHARGE_ESU (4.80320427e-10)
/* for numerical derivatives */
#define NETWORK_DIFFVAR (1e-6)
/* threshold for contructing the sparse matrices */
#define NETWORK_SPARSE_THRESHOLD (0.0)


#define set_zero(x) memset(&x, 0, sizeof(x))
#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

#if !(_BSD_SOURCE || _SVID_SOURCE || _XOPEN_SOURCE >= 500 || __STDC_VERSION__ >= 199901L)
#define cbrt(x) pow(x, 1./3.)
#endif


#if NETWORK_SPARSE
static void store_csr(jacob_t * jacob, int row_ind, const double *row, size_t n_row);
#endif /* NETWORK_SPARSE */

#ifdef NETWORK_SCREENING
void screening_apply(double rho, double temp, double ye, double yz, const struct network_data *nd,
		     struct network_workspace *nw);
static void compute_constant_screening_factors(struct network_rate *rate, const struct network_nucdata *nuc1,
					       const struct network_nucdata *nuc2,
					       const struct network_data *nd);
static network_var compute_screening(const struct network_rate *rate, network_var itemp13,
				     network_var Gamma_prime, network_var lGamma_prime,
				     network_var Gamma_prime14, network_var H_w_common);
#endif /* NETWORK_SCREENING */

/* functions for computing derivatives */
static __inline network_var network_var_add(const network_var a, const network_var b);
static __inline __attribute__ ((unused))
network_var network_var_const_add(const network_var a, const double b);
static __inline network_var network_var_sub(const network_var a, const network_var b);
static __inline network_var network_var_mul(const network_var a, const network_var b);
static __inline network_var network_var_const_mul(const network_var a, const double b);
static __inline __attribute__ ((unused))
network_var network_var_div(const network_var a, const network_var b);
#if (NETWORK_VARIABLE == NETWORK_VAR_TEMP)
static __inline network_var network_var_inverse(const network_var a);
#endif
static __inline network_var network_var_exp(const network_var a);
#ifdef NETWORK_SCREENING
static __inline network_var network_var_sqrt(const network_var a);
//static __inline network_var network_var_cbrt(const network_var a);
#endif
static __inline __attribute__ ((unused))
network_var network_var_log(const network_var a);
static __inline __attribute__ ((unused))
network_var network_var_pow_int(const network_var a, int n);

/* abbreviations to make the code more readable */
#define n_a(a, b) network_var_add(a, b)
#define n_ca(a, b) network_var_const_add(a, b)
#define n_s(a, b) network_var_sub(a, b)
#define n_m(a, b) network_var_mul(a, b)
#define n_cm(a, b) network_var_const_mul(a, b)
#define n_d(a, b) network_var_div(a, b)
#define n_inv(a) network_var_inverse(a)
#define n_exp(a) network_var_exp(a)
#define n_sqrt(a) network_var_sqrt(a)
#define n_cbrt(a) network_var_cbrt(a)
#define n_log(a) network_var_log(a)
#define n_powi(a, n) network_var_pow_int(a, n)

static network_var network_var_add(const network_var a, const network_var b)
{
  network_var res;

  res.v = a.v + b.v;
  res.drho = a.drho + b.drho;
  res.dT = a.dT + b.dT;
  res.dYe = a.dYe + b.dYe;
  res.dYz = a.dYz + b.dYz;
  res.dabar = a.dabar + b.dabar;
  res.dzbar = a.dzbar + b.dzbar;

  return res;
}

static network_var network_var_const_add(const network_var a, const double b)
{
  network_var res = a;

  res.v = a.v + b;

  return res;
}

static network_var network_var_sub(const network_var a, const network_var b)
{
  network_var res;

  res.v = a.v - b.v;
  res.drho = a.drho - b.drho;
  res.dT = a.dT - b.dT;
  res.dYe = a.dYe - b.dYe;
  res.dYz = a.dYz - b.dYz;
  res.dabar = a.dabar - b.dabar;
  res.dzbar = a.dzbar - b.dzbar;

  return res;
}

static network_var network_var_mul(const network_var a, const network_var b)
{
  network_var res;

  res.v = a.v * b.v;
  res.drho = a.v * b.drho + a.drho * b.v;
  res.dT = a.v * b.dT + a.dT * b.v;
  res.dYe = a.v * b.dYe + a.dYe * b.v;
  res.dYz = a.v * b.dYz + a.dYz * b.v;
  res.dabar = a.v * b.dabar + a.dabar * b.v;
  res.dzbar = a.v * b.dzbar + a.dzbar * b.v;

  return res;
}

static network_var network_var_const_mul(const network_var a, const double b)
{
  network_var res;

  res.v = a.v * b;
  res.drho = a.drho * b;
  res.dT = a.dT * b;
  res.dYe = a.dYe * b;
  res.dYz = a.dYz * b;
  res.dabar = a.dabar * b;
  res.dzbar = a.dzbar * b;

  return res;
}

static network_var network_var_div(const network_var a, const network_var b)
{
  network_var res;
  double bv2i = 1.0 / (b.v * b.v);

  res.v = a.v / b.v;
  res.drho = (a.drho * b.v - a.v * b.drho) * bv2i;
  res.dT = (a.dT * b.v - a.v * b.dT) * bv2i;
  res.dYe = (a.dYe * b.v - a.v * b.dYe) * bv2i;
  res.dYz = (a.dYz * b.v - a.v * b.dYz) * bv2i;
  res.dabar = (a.dabar * b.v - a.v * b.dabar) * bv2i;
  res.dzbar = (a.dzbar * b.v - a.v * b.dzbar) * bv2i;

  return res;
}

#if (NETWORK_VARIABLE == NETWORK_VAR_TEMP)
static network_var network_var_inverse(const network_var a)
{
  network_var res;
  double nav2i;

  res.v = 1.0 / a.v;
  nav2i = -res.v * res.v;
  res.drho = nav2i * a.drho;
  res.dT = nav2i * a.dT;
  res.dYe = nav2i * a.dYe;
  res.dYz = nav2i * a.dYz;
  res.dabar = nav2i * a.dabar;
  res.dzbar = nav2i * a.dzbar;

  return res;
}
#endif

static network_var network_var_exp(const network_var a)
{
  network_var res;

  res.v = exp(a.v);
  res.drho = res.v * a.drho;
  res.dT = res.v * a.dT;
  res.dYe = res.v * a.dYe;
  res.dYz = res.v * a.dYz;
  res.dabar = res.v * a.dabar;
  res.dzbar = res.v * a.dzbar;

  return res;
}

#ifdef NETWORK_SCREENING
static network_var network_var_sqrt(const network_var a)
{
  network_var res;
  double fac;

  res.v = sqrt(a.v);
  fac = 0.5 / res.v;
  res.drho = fac * a.drho;
  res.dT = fac * a.dT;
  res.dYe = fac * a.dYe;
  res.dYz = fac * a.dYz;
  res.dabar = fac * a.dabar;
  res.dzbar = fac * a.dzbar;

  return res;
}
#endif

/*
static network_var network_var_cbrt(const network_var a)
{
  network_var res;
  double fac;

  res.v = cbrt(a.v);
  fac = (1.0 / 3.0) * res.v / a.v;
  res.drho = fac * a.drho;
  res.dT = fac * a.dT;
  res.dYe = fac * a.dYe;
  res.dYz = fac * a.dYz;
  res.dabar = fac * a.dabar;
  res.dzbar = fac * a.dzbar;

  return res;
}
*/

static network_var network_var_log(const network_var a)
{
  network_var res;
  double avi = 1.0 / a.v;

  res.v = log(a.v);
  res.drho = avi * a.drho;
  res.dT = avi * a.dT;
  res.dYe = avi * a.dYe;
  res.dYz = avi * a.dYz;
  res.dabar = avi * a.dabar;
  res.dzbar = avi * a.dzbar;

  return res;
}

static network_var network_var_pow_int(const network_var a, int n)
{
  network_var res;
  double f = gsl_pow_int(a.v, n - 1);

  res.v = f * a.v;
  f = n * f;
  res.drho = f * a.drho;
  res.dT = f * a.dT;
  res.dYe = f * a.dYe;
  res.dYz = f * a.dYz;
  res.dabar = f * a.dabar;
  res.dzbar = f * a.dzbar;

  return res;
}

static void azbar(const double y[], double *abar, double *zbar, const struct network_data *nd);

static int network_init_internal(char *speciesfile, char *ratesfile, char *partfile, char *massesfile,
				 char *weakratesfile, struct network_data *nd, int onlyweak);

int network_init(char *speciesfile, char *ratesfile, char *partfile, char *massesfile, char *weakratesfile,
		 struct network_data *nd)
{
  return network_init_internal(speciesfile, ratesfile, partfile, massesfile, weakratesfile, nd, 0);
}

int network_init_onlyweak(char *speciesfile, char *ratesfile, char *partfile, char *massesfile,
			  char *weakratesfile, struct network_data *nd)
{
  return network_init_internal(speciesfile, ratesfile, partfile, massesfile, weakratesfile, nd, 1);
}

static int network_init_internal(char *speciesfile, char *ratesfile, char *partfile, char *massesfile,
				 char *weakratesfile, struct network_data *nd, int onlyweak)
{
  FILE *fd;
  char cdummy[200], cdummy2[200];
  int i, j, k, found, needed, rc;
  int na, nz, nn, nminz;
  float spin, exm, q;
  char nucnames[6][6];
  int ratetype, nallocated;
  int nucids[6];
  int perm, mult, *count, *ratecount;
  char *masses, *spins, missing;
  char nucEduct[6], nucProduct[6];
  int nucEductId, nucProductId;
  int isFirst, isWeak, zsum;

  myprintf("Network init with onlyweak: %d\n", onlyweak);

  /*
   *
   * read species data
   *
   */

  if(!(fd = fopen(speciesfile, "r")))
    {
      myprintf("can't open file `%s' for reading species information.\n", speciesfile);
      endrun(1);
    }
  safe_fgets(cdummy, 200, fd);
  sscanf(cdummy, "%zu", &nd->nuc_count);

  {
#ifdef NETWORK_SEP_YZ
    const size_t begin = nd->nuc_count + 1;
    nd->iYz = nd->nuc_count;
#else
    const size_t begin = nd->nuc_count;
#endif /* NETWORK_SEP_YZ */
#if NETWORK_VARIABLE == NETWORK_VAR_TEMP
    nd->iTemp = begin;
    nd->n_matrix = begin + 1;
#elif NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
    nd->iTemp = begin;
    nd->iRho = nd->iTemp + 1;
    nd->n_matrix = begin + 2;
#else
    nd->n_matrix = begin;
#endif /* NETWORK_VARIABLE */
  }

  /* allocate memory for species data */
  nd->nucdata = malloc(nd->nuc_count * sizeof(struct network_nucdata));

  for(i = 0; i < nd->nuc_count; i++)
    {
      safe_fgets(cdummy, 200, fd);
      sscanf(cdummy, "%5c%d%d", nd->nucdata[i].name, &nd->nucdata[i].na, &nd->nucdata[i].nz);
      nd->nucdata[i].name[5] = 0;
      nd->nucdata[i].nn = nd->nucdata[i].na - nd->nucdata[i].nz;
      nd->nucdata[i].nrates = 0;
      nd->nucdata[i].nweakrates = 0;
      nd->nucdata[i].spin = 0.0;
    }
  fclose(fd);

  /*
   *
   * read mass excess
   *
   */

  masses = calloc(nd->nuc_count, sizeof(char));

  if(!(fd = fopen(massesfile, "r")))
    {
      myprintf("can't open file `%s' for reading mass excess.\n", massesfile);
      endrun(2);
    }

  /* skip 39 lines */
  for(i = 0; i < 39; i++)
    {
      safe_fgets(cdummy, 200, fd);
    }

  while(!feof(fd))
    {
      if(fgets(cdummy, 200, fd) == NULL)
	break;
      sscanf(&cdummy[1], "%d%d%d%d\n", &nminz, &nn, &nz, &na);

      for(i = 0; i < nd->nuc_count; i++)
	{
	  if(nd->nucdata[i].na == na && nd->nucdata[i].nz == nz)
	    {
	      sscanf(&cdummy[29], "%f", &exm);
	      nd->nucdata[i].exm = exm;
	      sscanf(&cdummy[54], "%f", &q);
	      nd->nucdata[i].q = q * 1.602177e-12 * 1.0e3 * nd->nucdata[i].na;	/* values are in KeV / nucleon, converting to erg */
	      masses[i] = 1;
	      break;
	    }
	}
    }

  fclose(fd);

  missing = 0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      if(!masses[i])
	{
	  myprintf("Nucleus %s missing in mass file.\n", nd->nucdata[i].name);
	  missing = 1;
	}

      nd->nucdata[i].m =
	nd->nucdata[i].na * GSL_CONST_CGS_UNIFIED_ATOMIC_MASS +
	nd->nucdata[i].exm * conv / GSL_CONST_NUM_AVOGADRO / (GSL_CONST_CGS_SPEED_OF_LIGHT *
							      GSL_CONST_CGS_SPEED_OF_LIGHT);
    }

  if(missing)
    {
      endrun(3);
    }

  free(masses);

  /*
   *
   * read partition functions
   *
   */

  spins = calloc(nd->nuc_count, sizeof(char));

  if(!(fd = fopen(partfile, "r")))
    {
      myprintf("can't open file `%s' for reading partition functions.\n", partfile);
      endrun(4);
    }

  /* default value for the partition function is 1.0 => log(part) = 0.0 */
  for(i = 0; i < nd->nuc_count; i++)
    {
      for(j = 0; j < 24; j++)
	{
	  nd->nucdata[i].part[j] = 0.0;
	}
    }

  /* skip 4 lines */
  for(i = 0; i < 4; i++)
    {
      safe_fgets(cdummy, 200, fd);
    }

  while(!feof(fd))
    {
      if(fgets(cdummy, 200, fd) == NULL)
	break;			/* skip name of the nucleus */
      if(fgets(cdummy, 200, fd) == NULL)
	break;
      sscanf(cdummy, "%d%d%f\n", &nz, &na, &spin);

      found = 0;
      for(i = 0; i < nd->nuc_count; i++)
	{
	  if(nd->nucdata[i].na == na && nd->nucdata[i].nz == nz)
	    {
	      found = 1;

	      /* read partition function data */
	      for(j = 0; j < 3; j++)
		{
		  char *tmp, *next_field;

		  safe_fgets(cdummy, 200, fd);
		  next_field = cdummy;
		  for(k = 0; k < 8; k++)
		    {
		      nd->nucdata[i].part[j * 8 + k] = strtod(next_field, &tmp);
		      if(tmp == next_field)
			{
			  myprintf("Error reading partition function data for element %s\n",
				   nd->nucdata[i].name);
			  endrun(5);
			}
		      next_field = tmp;
		      nd->nucdata[i].part[j * 8 + k] = log(nd->nucdata[i].part[j * 8 + k]);
		    }
		}

	      nd->nucdata[i].spin = spin;
	      spins[i] = 1;
	      break;
	    }
	}

      /* skip 3 lines containing partition function data */
      if(!found)
	{
	  for(i = 0; i < 3; i++)
	    {
	      safe_fgets(cdummy, 200, fd);
	    }
	}
    }

  fclose(fd);

  for(i = 0; i < nd->nuc_count; i++)
    {
      if(!spins[i])
	{
	  myprintf
	    ("There are no partition function and spin data for nucleus %s. Assuming spin 0 and constant partition function of 1.\n",
	     nd->nucdata[i].name);
	}
    }
  free(spins);

  /*
   *
   * read weak rate data
   *
   */

  nd->weakrate_count = 0;
  nallocated = 1024;
  nd->weakrates = malloc(nallocated * sizeof(struct network_weakrate));

  if(!(fd = fopen(weakratesfile, "r")))
    {
      myprintf("can't open file `%s' for reading weak rate data.\n", weakratesfile);
      endrun(8);
    }

  isFirst = 1;
  while(!feof(fd))
    {
      if(fgets(cdummy, 200, fd) == NULL)
	break;

      memcpy(nucEduct, cdummy, 5);
      nucEduct[5] = 0;
      memcpy(nucProduct, &cdummy[5], 5);
      nucProduct[5] = 0;

      nucEductId = -1;
      nucProductId = -1;

      for(i = 0; i < nd->nuc_count; i++)
	{
	  if(!strcmp(nucEduct, nd->nucdata[i].name))
	    {
	      nucEductId = i;
	      if(nucProductId >= 0)
		break;
	    }
	  if(!strcmp(nucProduct, nd->nucdata[i].name))
	    {
	      nucProductId = i;
	      if(nucEductId >= 0)
		break;
	    }
	}

      if(nucEductId >= 0 && nucProductId >= 0)
	{
	  /* we need this rate */
	  rc = nd->weakrate_count;

	  /* check if array to store rates is large enough */
	  if(rc == nallocated)
	    {
	      nallocated += 1024;
	      nd->weakrates = realloc(nd->weakrates, nallocated * sizeof(struct network_weakrate));
	    }

	  nd->weakrates[rc].input = nucEductId;
	  nd->weakrates[rc].output = nucProductId;

	  sscanf(&cdummy[12], "%lf", &nd->weakrates[rc].q1);
	  sscanf(&cdummy[24], "%lf", &nd->weakrates[rc].q2);
	  sscanf(&cdummy[36], "%d", &nd->weakrates[rc].isReverse);

	  for(i = 0; i < 143; i++)
	    {
	      safe_fgets(cdummy, 200, fd);
	      sscanf(&cdummy[13], "%lf", &nd->weakrates[rc].lambda1[i]);
	      sscanf(&cdummy[22], "%lf", &nd->weakrates[rc].lambda2[i]);

	      if(isFirst)
		{
		  if(i < 13)
		    sscanf(&cdummy[0], "%f", &nd->weakTemp[i]);
		  if(i % 13 == 0)
		    sscanf(&cdummy[6], "%f", &nd->weakRhoYe[i / 13]);
		}
	    }

	  nd->nucdata[nucEductId].nweakrates++;
	  nd->nucdata[nucProductId].nweakrates++;

	  nd->weakrate_count++;

	  isFirst = 0;
	}
      else
	{
	  /* skip rate data, 13*11 = 143 lines */
	  for(i = 0; i < 143; i++)
	    safe_fgets(cdummy, 200, fd);
	}
    }

  /* free memory that is not needed */
  nd->weakrates = realloc(nd->weakrates, nd->weakrate_count * sizeof(struct network_weakrate));

  fclose(fd);

  /*
   *
   * read rate data
   *
   */

  count = malloc(nd->nuc_count * sizeof(int));

  nd->rate_count = 0;
  nallocated = 1024;
  nd->rates = malloc(nallocated * sizeof(struct network_rate));

  if(!(fd = fopen(ratesfile, "r")))
    {
      myprintf("can't open file `%s' for reading rate data.\n", ratesfile);
      endrun(6);
    }

  /* first read performed outside of the loop to ensure we hit EOF after the last rate */
  safe_fgets(cdummy, 200, fd);

  while(!feof(fd))
    {
#ifdef REACLIB1
      if(isdigit(cdummy[0]))
	{
	  sscanf(cdummy, "%d", &ratetype);

	  /* skip 2 lines and get new line */
	  for(i = 0; i < 3 && !feof(fd); i++)
	    {
	      safe_fgets(cdummy, 200, fd);
	    }
	}
#else
      /* in the REACLIB2 format the ratetype is given for each rate individually */
      if(sscanf(cdummy, "%d", &ratetype) != 1)
	break;
      safe_fgets(cdummy, 200, fd);
#endif
      for(i = 0; i < 6; i++)
	{
	  sscanf(&cdummy[(i + 1) * 5], "%5c", nucnames[i]);
	  nucnames[i][5] = 0;
	}

      /* check if we need this rate. we need it if it only contains nuclei we consider */
      needed = 1;

      for(i = 0; i < 6; i++)
	{
	  nucids[i] = -1;
	  found = 0;
	  if(strncmp(nucnames[i], "     ", 5))
	    {
	      /* do not check one nuclei again */
	      for(j = 0; j < i; j++)
		{
		  if(!strcmp(nucnames[i], nucnames[j]))
		    {
		      nucids[i] = nucids[j];
		      found = 1;
		      break;
		    }
		}

	      /* if found, we do not have to check it again */
	      if(found)
		{
		  continue;
		}

	      /* if new nuclei, check table of all nuclei used */
	      for(j = 0; j < nd->nuc_count; j++)
		{
		  if(!strcmp(nd->nucdata[j].name, nucnames[i]))
		    {
		      nucids[i] = j;
		      found = 1;
		      break;
		    }
		}

	      /* if not found, we do not need this rate */
	      if(!found)
		{
		  needed = 0;
		  break;
		}
	    }
	}


      /* check whether this is a weak rate and we already have it */
      if(needed && ratetype == 1)
	{
	  for(i = 0; i < nd->weakrate_count; i++)
	    {
	      nucEductId = nd->weakrates[i].input;
	      nucProductId = nd->weakrates[i].output;

	      if(nucids[0] == nucEductId && nucids[1] == nucProductId)
		{
		  /* found it */
		  needed = 0;
		  break;
		}
	    }
	}


      if(needed)
	{
	  rc = nd->rate_count;

	  /* check if array to store rates is large enough */
	  if(rc == nallocated)
	    {
	      nallocated += 1024;
	      nd->rates = realloc(nd->rates, nallocated * sizeof(struct network_rate));
	    }

	  switch (ratetype)
	    {
	    case 1:
	      nd->rates[rc].ninput = 1;
	      nd->rates[rc].noutput = 1;
	      break;
	    case 2:
	      nd->rates[rc].ninput = 1;
	      nd->rates[rc].noutput = 2;
	      break;
	    case 3:
	      nd->rates[rc].ninput = 1;
	      nd->rates[rc].noutput = 3;
	      break;
	    case 4:
	      nd->rates[rc].ninput = 2;
	      nd->rates[rc].noutput = 1;
	      break;
	    case 5:
	      nd->rates[rc].ninput = 2;
	      nd->rates[rc].noutput = 2;
	      break;
	    case 6:
	      nd->rates[rc].ninput = 2;
	      nd->rates[rc].noutput = 3;
	      break;
	    case 7:
	      nd->rates[rc].ninput = 2;
	      nd->rates[rc].noutput = 4;
	      break;
	    case 8:
	      nd->rates[rc].ninput = 3;
	      if(nucids[4] == -1)
		{
		  nd->rates[rc].noutput = 1;
		}
	      else
		{
		  nd->rates[rc].noutput = 2;
		}
	      break;
	    case 9:
	      nd->rates[rc].ninput = 3;
	      nd->rates[rc].noutput = 2;
	      break;
	    case 10:
	      nd->rates[rc].ninput = 4;
	      nd->rates[rc].noutput = 2;
	      break;
	    case 11:
	      nd->rates[rc].ninput = 1;
	      nd->rates[rc].noutput = 4;
	      break;
	    default:
	      myprintf("Ratetype %d is not supported\n", ratetype);
	      needed = 0;
	      break;
	    }

	  if(needed)
	    {
	      for(i = 0; i < nd->rates[rc].ninput; i++)
		nd->rates[rc].input[i] = nucids[i];
	      for(i = 0; i < nd->rates[rc].noutput; i++)
		nd->rates[rc].output[i] = nucids[i + nd->rates[rc].ninput];
	    }
	}


      if(needed)
	{
	  /* check if this rate is a weak rate
	     it is a weak rate exactly if the reaction changes the total number of protons */
	  zsum = 0;
	  for(i = 0; i < nd->rates[rc].ninput; i++)
	    zsum += nd->nucdata[nd->rates[rc].input[i]].nz;
	  for(i = 0; i < nd->rates[rc].noutput; i++)
	    zsum -= nd->nucdata[nd->rates[rc].output[i]].nz;

	  if(zsum != 0)
	    {
	      isWeak = 1;
	    }
	  else
	    {
	      isWeak = 0;
	      /* if we want to load weak rates only, we do not need it */
	      if(onlyweak)
		needed = 0;
	    }
	}


      if(needed)
	{
	  rc = nd->rate_count;

	  nd->rates[rc].type = ratetype;
	  nd->rates[rc].index = rc;
	  nd->rates[rc].isWeak = isWeak;	/* (cdummy[47] == 'w') seems to be unused */
	  nd->rates[rc].isReverse = (cdummy[48] == 'v');
	  nd->rates[rc].isElectronCapture = (strcmp(&cdummy[43], "  ec") == 0) ? 1 : 0;
	  sscanf(&cdummy[52], "%lf\n", &nd->rates[rc].q);

	  safe_fgets(cdummy2, 200, fd);
	  for(i = 0; i < 4; i++)
	    {
	      sscanf(&cdummy2[i * 13], "%lf", &nd->rates[rc].data[i]);
	    }
	  safe_fgets(cdummy2, 200, fd);
	  for(i = 4; i < 7; i++)
	    {
	      sscanf(&cdummy2[(i - 4) * 13], "%lf", &nd->rates[rc].data[i]);
	    }

	  /* mark nuclei that are affected by this rate */
	  memset(count, 0, nd->nuc_count * sizeof(int));

	  for(i = 0; i < 6; i++)
	    {
	      if(nucids[i] != -1)
		{
		  if((ratetype < 4 && i < 1) || (ratetype >= 4 && ratetype < 8 && i < 2)
		     || (ratetype == 8 && i < 3))
		    {
		      count[nucids[i]]--;
		    }
		  else
		    {
		      count[nucids[i]]++;
		    }
		}
	    }

	  /* nuclei are only affected if their number is really changed
	   * this is i.e. not the case for catalysts like protons in the reaction
	   * p + he4 + ru88 => p + pd92 */
	  for(i = 0; i < nd->nuc_count; i++)
	    {
	      if(count[i] != 0)
		nd->nucdata[i].nrates++;
	    }

	  nd->rate_count++;
	}
      else
	{
	  /* skip 2 lines */
	  for(i = 0; i < 2; i++)
	    {
	      safe_fgets(cdummy, 200, fd);
	    }
	}

      /* This read is for the next step. If we hit EOF, the loop will terminate. */
      if(fgets(cdummy, 200, fd) == NULL)
	break;
    }

  /* free memory that is not needed */
  nd->rates = realloc(nd->rates, nd->rate_count * sizeof(struct network_rate));

  fclose(fd);

  for(i = 0; i < nd->nuc_count; i++)
    {
      nd->nucdata[i].rates = malloc(nd->nucdata[i].nrates * sizeof(struct network_rate));
      nd->nucdata[i].prates = malloc(nd->nucdata[i].nrates * sizeof(struct network_rate *));
      nd->nucdata[i].irates = malloc(nd->nucdata[i].nrates * sizeof(size_t));
      nd->nucdata[i].w = malloc(nd->nucdata[i].nrates * sizeof(double));

      nd->nucdata[i].weakrates = malloc(nd->nucdata[i].nweakrates * sizeof(struct network_weakrate));
      nd->nucdata[i].pweakrates = malloc(nd->nucdata[i].nweakrates * sizeof(struct network_weakrate *));
      nd->nucdata[i].iweakrates = malloc(nd->nucdata[i].nweakrates * sizeof(size_t));
      nd->nucdata[i].wweak = malloc(nd->nucdata[i].nweakrates * sizeof(double));
    }

  ratecount = calloc(nd->nuc_count, sizeof(int));

  for(i = 0; i < nd->rate_count; i++)
    {
      memset(count, 0, nd->nuc_count * sizeof(int));
      for(j = 0; j < nd->rates[i].ninput; j++)
	count[nd->rates[i].input[j]]--;
      for(j = 0; j < nd->rates[i].noutput; j++)
	count[nd->rates[i].output[j]]++;

      /* permutation factor that is important when more than one nucleus of a kind is involved */
      perm = 1;
      for(j = 0; j < nd->rates[i].ninput; j++)
	{
	  mult = 1;
	  for(k = 0; k < nd->rates[i].ninput; k++)
	    {
	      if(nd->rates[i].input[j] == nd->rates[i].input[k])
		{
		  if(k < j)
		    {
		      /* we already did this kind */
		      break;
		    }
		  else if(k > j)
		    {
		      /* another one of this kind */
		      mult++;
		    }
		}
	    }

	  for(k = 1; k <= mult; k++)
	    {
	      perm *= k;
	    }
	}

      for(j = 0; j < nd->nuc_count; j++)
	{
	  if(count[j] != 0)
	    {
	      rc = ratecount[j];
	      if(rc == nd->nucdata[j].nrates)
		{
		  /* this should never happen */
		  myprintf("rate %d, %s: ratecount too low (%d <> %zu), stopping.\n", i, nd->nucdata[j].name,
			   rc, nd->nucdata[j].nrates);
		  endrun(9);
		}
	      nd->nucdata[j].irates[rc] = i;
	      nd->nucdata[j].rates[rc] = nd->rates[i];
	      nd->nucdata[j].prates[rc] = &nd->rates[i];
	      nd->nucdata[j].w[rc] = (double) count[j] / (double) perm;
	      ratecount[j]++;
	    }
	}
    }

  /* do the same for weak rates */
  memset(ratecount, 0, nd->nuc_count * sizeof(int));
  for(i = 0; i < nd->weakrate_count; i++)
    {
      nd->nucdata[nd->weakrates[i].input].iweakrates[ratecount[nd->weakrates[i].input]] = i;
      nd->nucdata[nd->weakrates[i].input].weakrates[ratecount[nd->weakrates[i].input]] = nd->weakrates[i];
      nd->nucdata[nd->weakrates[i].input].pweakrates[ratecount[nd->weakrates[i].input]] = &nd->weakrates[i];
      nd->nucdata[nd->weakrates[i].input].wweak[ratecount[nd->weakrates[i].input]] = -1.0;
      ratecount[nd->weakrates[i].input]++;

      nd->nucdata[nd->weakrates[i].output].iweakrates[ratecount[nd->weakrates[i].output]] = i;
      nd->nucdata[nd->weakrates[i].output].weakrates[ratecount[nd->weakrates[i].output]] = nd->weakrates[i];
      nd->nucdata[nd->weakrates[i].output].pweakrates[ratecount[nd->weakrates[i].output]] = &nd->weakrates[i];
      nd->nucdata[nd->weakrates[i].output].wweak[ratecount[nd->weakrates[i].output]] = 1.0;
      ratecount[nd->weakrates[i].output]++;
    }

  free(ratecount);
  free(count);

  /* do electron capture rates externally */
  nd->electronCapture_count = 0;
  for(i = 0; i < nd->rate_count; i++)
    {
      if(nd->rates[i].isElectronCapture)
	nd->electronCapture_count++;
    }

  nd->electronCapture_rates = malloc(nd->electronCapture_count * sizeof(size_t));
  nd->electronCapture_count = 0;
  for(i = 0; i < nd->rate_count; i++)
    {
      if(nd->rates[i].isElectronCapture)
	{
	  nd->electronCapture_rates[nd->electronCapture_count] = i;
	  nd->electronCapture_count++;
	}
    }

  /* count the screened rates */
  ratecount = calloc(nd->rate_count, sizeof(int));
  nd->screened_count = 0;
  nd->tripleAlpha_rate = NULL;
  for(i = 0; i < nd->rate_count; i++)
    {
      if(nd->rates[i].ninput == 2)
	{
	  if(nd->nucdata[nd->rates[i].input[0]].nz != 0 && nd->nucdata[nd->rates[i].input[1]].nz != 0)
	    {
	      ratecount[i] = 1;
	      nd->screened_count++;
	    }
	}
      else if(nd->rates[i].ninput == 3)
	{
	  struct network_nucdata *in0, *in1, *in2;

	  in0 = &nd->nucdata[nd->rates[i].input[0]];
	  in1 = &nd->nucdata[nd->rates[i].input[1]];
	  in2 = &nd->nucdata[nd->rates[i].input[2]];

	  if(in0->nz == 2 && in0->nn == 2 && in1->nz == 2 && in1->nn == 2 && in2->nz == 2 && in2->nn == 2)
	    {
	      nd->tripleAlpha_rate = &nd->rates[i];
	    }
	}
    }

  if(nd->tripleAlpha_rate == NULL && !onlyweak)
    {
      myprintf("triple alpha rate not found\n");
    }

  /* build lists for the screened and tripleAlpha rates */
  nd->screened_rates = malloc(nd->screened_count * sizeof(struct network_rate *));
  nd->screened_count = 0;
  for(i = 0; i < nd->rate_count; i++)
    {
      if(ratecount[i])
	{
	  nd->screened_rates[nd->screened_count] = &nd->rates[i];
	  nd->screened_count++;
	}
    }

#ifdef NETWORK_SCREENING
  /* compute the constant factors for screening */
  for(i = 0; i < nd->screened_count; i++)
    {
      compute_constant_screening_factors(nd->screened_rates[i], NULL, NULL, nd);
    }
  {
    const struct network_nucdata he4 = {.nz = 2,.na = 4 };
    const struct network_nucdata be8 = {.nz = 4,.na = 8 };

    compute_constant_screening_factors(&nd->alpha_alpha_rate, &he4, &he4, nd);
    compute_constant_screening_factors(&nd->Be8_alpha_rate, &he4, &be8, nd);
  }
#endif /* NETWORK_SCREENING */

  free(ratecount);

  nd->initialized = 1;

  /* check if there are species without any rates */
  if(!onlyweak)
    {
      int any = 0;
      for(i = 0; i < nd->nuc_count; i++)
	{
	  if(nd->nucdata[i].nrates == 0 && nd->nucdata[i].nweakrates == 0)
	    {
	      myprintf("There are no rates connecting nucleus %d (%6s)\n", i, nd->nucdata[i].name);
	      any = 1;
	    }
	}
      if(any)
	endrun(11);
    }
  myprintf("Network init done, %zd species found, %zd rates loaded.\n", nd->nuc_count,
	   nd->rate_count + nd->weakrate_count);

  return 0;
}

int network_workspace_init(const struct network_data *nd, struct network_workspace *nw)
{
  int i, j;

  nw->nuc_count = nd->nuc_count;

  nw->y = malloc(nd->n_matrix * sizeof(double));
  nw->gg = malloc(nd->nuc_count * sizeof(network_var));
  nw->matrix = malloc(nd->n_matrix * nd->n_matrix * sizeof(double));

  nw->old_temp = -1.0;
  nw->old_rho = -1.0;
  nw->old_ye = -1.0;

  nw->baserate = malloc(nd->rate_count * sizeof(network_var));
  nw->unscreened_rate = malloc(nd->rate_count * sizeof(network_var));
  nw->rates = malloc(nd->rate_count * sizeof(network_var));
  nw->weakrates = malloc(nd->weakrate_count * sizeof(network_var));
  nw->yrates = malloc(nd->rate_count * sizeof(network_var));
  nw->yweakrates = malloc(nd->weakrate_count * sizeof(network_var));

  nw->prates = malloc(nd->nuc_count * sizeof(network_var *));
  nw->pweakrates = malloc(nd->nuc_count * sizeof(network_var **));

  for(i = 0; i < nd->nuc_count; i++)
    {
      nw->prates[i] = malloc(nd->nucdata[i].nrates * sizeof(network_var *));
      for(j = 0; j < nd->nucdata[i].nrates; j++)
	{
	  nw->prates[i][j] = &nw->yrates[nd->nucdata[i].irates[j]];
	}
      nw->pweakrates[i] = malloc(nd->nucdata[i].nweakrates * sizeof(network_var *));
      for(j = 0; j < nd->nucdata[i].nweakrates; j++)
	{
	  nw->pweakrates[i][j] = &nw->yweakrates[nd->nucdata[i].iweakrates[j]];
	}
    }

  nw->dTdYi = malloc(nd->nuc_count * sizeof(network_var));
  nw->drhodYi = malloc(nd->nuc_count * sizeof(network_var));

  nw->nsd = network_solver_init(1e-6, nd->n_matrix, nd->nuc_count, nd->nucdata);

  nw->initialized = 1;
  return 0;
}

void network_deinit(struct network_data *nd)
{
  int i;

  for(i = 0; i < nd->nuc_count; i++)
    {
      free(nd->nucdata[i].rates);
      free(nd->nucdata[i].prates);
      free(nd->nucdata[i].irates);
      free(nd->nucdata[i].w);

      free(nd->nucdata[i].weakrates);
      free(nd->nucdata[i].pweakrates);
      free(nd->nucdata[i].iweakrates);
      free(nd->nucdata[i].wweak);
    }

  free(nd->electronCapture_rates);
  free(nd->screened_rates);

  free(nd->nucdata);
  free(nd->rates);
  free(nd->weakrates);

  nd->initialized = 0;
}

void network_workspace_deinit(struct network_workspace *nw)
{
  int i;

  network_solver_deinit(nw->nsd);

  for(i = 0; i < nw->nuc_count; i++)
    {
      free(nw->prates[i]);
      free(nw->pweakrates[i]);
    }

  free(nw->y);
  free(nw->gg);
  free(nw->matrix);

  free(nw->baserate);
  free(nw->unscreened_rate);
  free(nw->rates);
  free(nw->weakrates);
  free(nw->yrates);
  free(nw->yweakrates);
  free(nw->prates);
  free(nw->pweakrates);
  free(nw->dTdYi);
  free(nw->drhodYi);

  nw->initialized = 0;
}

int network_getrhs(double rho, double temp, const double y[], int compute_derivs,
		   const struct network_data *nd, struct network_workspace *nw, double *rhs,
		   network_var * deriv)
{
  int j;
  double ne, nn, ye, yz;
  struct network_nucdata *nucdata, *nucend;
  struct network_rate *ratedata, *rateend;
  struct network_weakrate *weakratedata, *weakrateend;
  const double *yy;
  network_var *yrate, *rateval, *rhsiter, ***pratelist, ***pweakratelist;

  /* the following statements are important as the temp and rho arguments are not updated */
#ifdef NETWORK_SEP_YZ
  yz = y[nd->iYz];
#endif
#ifdef NETWORK_VARIABLE
  temp = y[nd->iTemp];
#endif
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  rho = y[nd->iRho];
#endif

  nn = 0.0;
  ne = 0.0;
  yy = y;
#ifndef NETWORK_SEP_YZ
  yz = 0.0;
#endif
  nucend = &nd->nucdata[nd->nuc_count];
  for(nucdata = nd->nucdata; nucdata != nucend; nucdata++)
    {
      nn += (*yy) * (*nucdata).na;
      ne += (*yy) * (*nucdata).nz;

#ifndef NETWORK_SEP_YZ
      yz += gsl_pow_2(nucdata->nz) * (*yy);
#endif
      yy++;
    }
  ye = ne / nn;

  getrates(rho, temp, ye, yz, compute_derivs, nd, nw);

  /* dy_i/dt */
  rateend = &nd->rates[nd->rate_count];
  for(yrate = nw->yrates, ratedata = nd->rates, rateval = nw->rates;
      ratedata != rateend; yrate++, ratedata++, rateval++)
    {
      /* The derivative w.r.t. y[i] is treated separately. */
      *yrate = n_cm(*rateval, y[ratedata->input[0]]);
      for(j = 1; j < ratedata->ninput; j++)
	*yrate = n_cm(*yrate, y[ratedata->input[j]]);
    }

  weakrateend = &nd->weakrates[nd->weakrate_count];
  for(yrate = nw->yweakrates, weakratedata = nd->weakrates, rateval = nw->weakrates;
      weakratedata != weakrateend; yrate++, weakratedata++, rateval++)
    {
      /* The derivative w.r.t. y[i] is treated separately. */
      *yrate = n_cm(*rateval, y[weakratedata->input]);
    }

  nucend = &nd->nucdata[nd->nuc_count];
  for(nucdata = nd->nucdata, rhsiter = deriv, pratelist = nw->prates, pweakratelist = nw->pweakrates;
      nucdata != nucend; nucdata++, rhsiter++, pratelist++, pweakratelist++)
    {
      network_var a = { 0 };
      network_var **prate;
      double *w, *wend;

      /* strong rates */
      wend = &nucdata->w[nucdata->nrates];
      for(w = nucdata->w, prate = *pratelist; w != wend; w++, prate++)
	{
	  a = n_a(a, n_cm(**prate, *w));
	}

      /* weak rates */
      wend = &nucdata->wweak[nucdata->nweakrates];
      for(w = nucdata->wweak, prate = *pweakratelist; w != wend; w++, prate++)
	{

	  a = n_a(a, n_cm(**prate, *w));
	}

      *rhsiter = a;
    }

#ifdef NETWORK_SEP_YZ
  {
    int i;
    network_var val = { 0 };

    for(i = 0; i < nd->nuc_count; i++)
      {
	val = n_a(val, n_cm(deriv[i], gsl_pow_2(nd->nucdata[i].nz)));
      }
    deriv[nd->iYz] = val;
  }
#endif /* NETWORK_SEP_YZ */
#if NETWORK_VARIABLE
  {
    int i;
    network_var dedT_var = { 0 }, dedtime =
    {
    0};
    network_var dp_drho = { 0 }, dpdT =
    {
    0}, dedrho =
    {
    0};
#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
    network_var dedA = { 0 }, dedZ =
    {
    0}, dpdA =
    {
    0}, dpdZ =
    {
    0};
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */
    struct eos_result res;

#endif /* NETWORK_VARIABLE */

#if NETWORK_VARIABLE == NETWORK_VAR_TEMP || NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP

    {
      double abar, zbar;
      double diff, old;
      struct helm_eos_cache cache;

      azbar(y, &abar, &zbar, nd);
      helm_eos_update_cache(rho, abar, zbar, &cache);

      eos_calc_tgiven_azbar(rho, &cache, temp, &res, 0);
      dedT_var.v = res.e.dtemp;
      dp_drho.v = res.p.drho;
      dpdT.v = res.p.dtemp;
      dedrho.v = res.e.drho;
#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
      dedA.v = res.e.dabar;
      dedZ.v = res.e.dzbar;
      dpdA.v = res.p.dabar;
      dpdZ.v = res.p.dzbar;
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */

      if(compute_derivs)
	{
	  /* numerical derivative of dedT */

	  diff = max(NETWORK_DIFFVAR, fabs(abar) * NETWORK_DIFFVAR);
	  old = abar;
	  abar += diff;
	  helm_eos_update_cache(rho, abar, zbar, &cache);
	  eos_calc_tgiven_azbar(rho, &cache, temp, &res, 0);
	  abar = old;
	  dedT_var.dabar = (res.e.dtemp - dedT_var.v) / diff;
	  dp_drho.dabar = (res.p.drho - dp_drho.v) / diff;
	  dpdT.dabar = (res.p.dtemp - dpdT.v) / diff;
	  dedrho.dabar = (res.e.drho - dedrho.v) / diff;
#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
	  dedA.dabar = (res.e.dabar - dedA.v) / diff;
	  dedZ.dabar = (res.e.dzbar - dedZ.v) / diff;
	  dpdA.dabar = (res.p.dabar - dpdA.v) / diff;
	  dpdZ.dabar = (res.p.dzbar - dpdZ.v) / diff;
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */

	  diff = max(NETWORK_DIFFVAR, fabs(zbar) * NETWORK_DIFFVAR);
	  old = zbar;
	  zbar += diff;
	  helm_eos_update_cache(rho, abar, zbar, &cache);
	  eos_calc_tgiven_azbar(rho, &cache, temp, &res, 0);
	  zbar = old;
	  dedT_var.dzbar = (res.e.dtemp - dedT_var.v) / diff;
	  dp_drho.dzbar = (res.p.drho - dp_drho.v) / diff;
	  dpdT.dzbar = (res.p.dtemp - dpdT.v) / diff;
	  dedrho.dzbar = (res.e.drho - dedrho.v) / diff;
#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
	  dedA.dzbar = (res.e.dabar - dedA.v) / diff;
	  dedZ.dzbar = (res.e.dzbar - dedZ.v) / diff;
	  dpdA.dzbar = (res.p.dabar - dpdA.v) / diff;
	  dpdZ.dzbar = (res.p.dzbar - dpdZ.v) / diff;
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */

	  diff = max(NETWORK_DIFFVAR, fabs(temp) * NETWORK_DIFFVAR);
	  old = temp;
	  temp += diff;
	  helm_eos_update_cache(rho, abar, zbar, &cache);
	  eos_calc_tgiven_azbar(rho, &cache, temp, &res, 0);
	  temp = old;
	  dedT_var.dT = (res.e.dtemp - dedT_var.v) / diff;
	  dp_drho.dT = (res.p.drho - dp_drho.v) / diff;
	  dpdT.dT = (res.p.dtemp - dpdT.v) / diff;
	  dedrho.dT = (res.e.drho - dedrho.v) / diff;
#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
	  dedA.dT = (res.e.dabar - dedA.v) / diff;
	  dedZ.dT = (res.e.dzbar - dedZ.v) / diff;
	  dpdA.dT = (res.p.dabar - dpdA.v) / diff;
	  dpdZ.dT = (res.p.dzbar - dpdZ.v) / diff;
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */

#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
	  diff = max(NETWORK_DIFFVAR, fabs(rho) * NETWORK_DIFFVAR);
	  old = rho;
	  rho += diff;
	  helm_eos_update_cache(rho, abar, zbar, &cache);
	  eos_calc_tgiven_azbar(rho, &cache, temp, &res, 0);
	  rho = old;
	  dedT_var.drho = (res.e.dtemp - dedT_var.v) / diff;
	  dp_drho.drho = (res.p.drho - dp_drho.v) / diff;
	  dpdT.drho = (res.p.dtemp - dpdT.v) / diff;
	  dedrho.drho = (res.e.drho - dedrho.v) / diff;
#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
	  dedA.drho = (res.e.dabar - dedA.v) / diff;
	  dedZ.drho = (res.e.dzbar - dedZ.v) / diff;
	  dpdA.drho = (res.p.dabar - dpdA.v) / diff;
	  dpdZ.drho = (res.p.dzbar - dpdZ.v) / diff;
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_T */
	}
    }

    /* dT/dt = dT/dE * dE/dt + sum_i ( dT/dy_i * dy_i/dt )
     * dE/dt = sum_i ( ebind_i * dy_i/dt ) */
    for(i = 0; i < nd->nuc_count; i++)
      {
	dedtime = n_a(dedtime, n_cm(deriv[i], nd->nucdata[i].exm));
      }

    dedtime = n_cm(dedtime, -conv);
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
    /* case of fixed pressure
     * dT/dE = dp_drho / (dp_drho * dEdT - dPdT * dEdrho)
     */
    nw->temp_deriv_fac = n_d(dp_drho, n_s(n_m(dp_drho, dedT_var), n_m(dpdT, dedrho)));
#else
    /* case of fixed density
     * dT/dE = 1 / (dEdT)
     */
    nw->temp_deriv_fac = network_var_inverse(dedT_var);
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP */

    deriv[nd->iTemp] = n_m(dedtime, nw->temp_deriv_fac);

#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
    if(compute_derivs)
      {
	const network_var abar = {.v = res.abar,.dabar = 1.0 };
	const network_var zbar = {.v = res.zbar,.dzbar = 1.0 };
	network_var dTdZ, dTdA, dAdyi, dTdA_dAdyi;

#if NETWORK_VARIABLE == NETWORK_VAR_TEMP
	/* fixed density */
	/* dT/dZ|e = - (de/dZ) / (de/dT) ... */
	dTdA = n_cm(n_d(dedA, dedT_var), -1.0);
	dTdZ = n_cm(n_d(dedZ, dedT_var), -1.0);
#elif NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
	/* fixed pressure */
	/* dT/dA|P = (de/dA * dp/drho - de/drho * dp/dA) / (de/drho * dp/dT - de/dT * dp/drho) ... */
	{
	  const network_var denominator = n_s(n_m(dedrho, dpdT), n_m(dedT_var, dp_drho));

	  dTdA = n_d(n_s(n_m(dedA, dp_drho), n_m(dedrho, dpdA)), denominator);
	  dTdZ = n_d(n_s(n_m(dedZ, dp_drho), n_m(dedrho, dpdZ)), denominator);
	}
#else
#error "unknown mode"
#endif /* NETWORK_VARIABLE == NETWORK_VAR_TEMP */
	dAdyi = n_cm(n_powi(abar, 2), -1.0);
	dTdA_dAdyi = n_m(dTdA, dAdyi);
	/* sum_i ( dT/dy_i * dy_i/dt )
	 * Z == zbar, A == abar
	 * dT/dy_i = dT/dZ * dZ/dy_i + dT/dA * dA/dy_i
	 * dA/dy_i = - A**2
	 * dZ/dy_i = - A * Z + A * Z_i
	 */
	for(i = 0; i < nd->nuc_count; i++)
	  {
	    network_var dZdy = {.v = abar.v * (nd->nucdata[i].nz - zbar.v),
	      .dabar = nd->nucdata[i].nz - zbar.v,
	      .dzbar = -abar.v
	    };

	    nw->dTdYi[i] = n_a(dTdA_dAdyi, n_m(dTdZ, dZdy));
	    deriv[nd->iTemp] = n_a(deriv[nd->iTemp], n_m(deriv[i], nw->dTdYi[i]));
	  }
      }
    else
      {
#if NETWORK_VARIABLE == NETWORK_VAR_TEMP
	const double dTdA = -dedA.v / dedT_var.v;
	const double dTdZ = -dedZ.v / dedT_var.v;
#elif NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
	const double idenominator = 1.0 / (dedrho.v * dpdT.v - dedT_var.v * dp_drho.v);
	const double dTdA = (dedA.v * dp_drho.v - dedrho.v * dpdA.v) * idenominator;
	const double dTdZ = (dedZ.v * dp_drho.v - dedrho.v * dpdZ.v) * idenominator;
#else
#error "unknown mode"
#endif /* NETWORK_VARIABLE == NETWORK_VAR_TEMP */
	const double scr1 = (dTdA * res.abar + dTdZ * res.zbar) * res.abar;
	const double scr2 = dTdZ * res.abar;
	double sum = 0.0;

	for(i = 0; i < nd->nuc_count; i++)
	  {
	    sum += deriv[i].v * (scr2 * nd->nucdata[i].nz - scr1);
	  }

	deriv[nd->iTemp].v += sum;
      }
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */

#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
    {
      /* so far only fixed pressure is implemented
       * Here rho is a function of E, P, Abar, Zbar
       * drho/dt|P = drho/dE|(P,A,Z) * dE/dt + drho/dA|(E,P,Z) * dA/dt
       *             + drho/dZ|(E,P,A) * dZ/dt
       * dE/dt see above (saved in dedtime)
       * drho/dE|(P,A,Z) = dP/dT / ( dP/dT * dE/rho - dE/dT * dP/drho )
       */

      nw->rho_deriv_fac = n_d(dpdT, n_s(n_m(dpdT, dedrho), n_m(dedT_var, dp_drho)));
      deriv[nd->iRho] = n_m(dedtime, nw->rho_deriv_fac);
#ifndef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
      {
	const network_var abar = {.v = res.abar,.dabar = 1.0 };
	const network_var zbar = {.v = res.zbar,.dzbar = 1.0 };
	network_var drhodZ, drhodA, dAdyi, drhodA_dAdyi;

	/* drho/dA|P = (de/dT * dp/dA - dp/dT * de/dA) / (de/drho * dp/dT - dp/drho * de/dT) */

	{
	  const network_var denominator = n_s(n_m(dedrho, dpdT), n_m(dedT_var, dp_drho));

	  drhodA = n_d(n_s(n_m(dedT_var, dpdA), n_m(dpdT, dedA)), denominator);
	  drhodZ = n_d(n_s(n_m(dedT_var, dpdZ), n_m(dpdT, dedZ)), denominator);
	}
	dAdyi = n_cm(n_m(abar, abar), -1.0);
	drhodA_dAdyi = n_m(drhodA, dAdyi);
	/* sum_i ( dT/dy_i * dy_i/dt )
	 * Z == zbar, A == abar
	 * dT/dy_i = dT/dZ * dZ/dy_i + dT/dA * dA/dy_i
	 * dA/dy_i = - A**2
	 * dZ/dy_i = - A * Z + A * Z_i
	 */
	for(i = 0; i < nd->nuc_count; i++)
	  {
	    network_var dZdy = {.v = abar.v * (nd->nucdata[i].nz - zbar.v),
	      .dabar = nd->nucdata[i].nz - zbar.v,
	      .dzbar = -abar.v
	    };

	    nw->drhodYi[i] = n_a(drhodA_dAdyi, n_m(drhodZ, dZdy));
	    deriv[nd->iRho] = n_a(deriv[nd->iRho], n_m(deriv[i], nw->drhodYi[i]));
	  }
      }
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */
    }
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP */
  }
#endif /* NETWORK_VARIABLE */
  {
    int i;

    for(i = 0; i < nd->n_matrix; i++)
      {
	rhs[i] = deriv[i].v;
      }
  }

  return 0;
}

int network_getjacob(const double y[], const network_var * rhs, const struct network_data *nd,
		     struct network_workspace *nw, jacob_t * jacob)
{
  double abar, zbar, ye;
  double dabardyi;
  size_t n_matrix2 = nd->n_matrix * nd->n_matrix;

#if NETWORK_SPARSE
  double *matrix = nw->matrix;	/* temporary space */

  jacob->n_elements = 0;
#else
  double *matrix = jacob->matrix;
#endif

  memset(matrix, 0, n_matrix2 * sizeof(double));

  azbar(y, &abar, &zbar, nd);
  ye = zbar / abar;
  dabardyi = -abar * abar;

  {
    int i;

    for(i = 0; i < nd->nuc_count; i++)
      {
	int j, k, l;
	/* row by row */
	double *row = &matrix[i * nd->n_matrix];

	/* drhs[y_i]/dy_j (j != the index in the next loop) */
	for(j = 0; j < nd->nucdata[i].nrates; j++)
	  {
	    struct network_rate *rate = &nd->nucdata[i].rates[j];
	    network_var rateval = n_cm(nw->rates[nd->nucdata[i].irates[j]], nd->nucdata[i].w[j]);

	    for(k = 0; k < rate->ninput; k++)
	      {
		double tempderiv = rateval.v;

		/* this ensures correct treatment of rates depending on Y**2, Y**3, ... */
		for(l = 0; l < rate->ninput; l++)
		  {
		    if(l != k)
		      tempderiv *= y[rate->input[l]];
		  }

		row[rate->input[k]] += tempderiv;
	      }

	  }

	/* weak rates */
	for(j = 0; j < nd->nucdata[i].nweakrates; j++)
	  {
	    struct network_weakrate *rate = &nd->nucdata[i].weakrates[j];
	    network_var rateval = n_cm(nw->weakrates[nd->nucdata[i].iweakrates[j]], nd->nucdata[i].wweak[j]);

	    row[rate->input] += rateval.v;
	  }

#if NETWORK_VARIABLE
	/* use the row to compute the iTemp and iRho rows
	 * before adding the other derivatives */
	{
	  double *temp_row = &matrix[nd->iTemp * nd->n_matrix];
#ifdef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
	  const double temp_fac = nd->nucdata[i].exm * (-conv) * nw->temp_deriv_fac.v;
#else
	  const double temp_fac = nd->nucdata[i].exm * (-conv) * nw->temp_deriv_fac.v + nw->dTdYi[i].v;
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
	  double *rho_row = &matrix[nd->iRho * nd->n_matrix];
#ifdef NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS
	  const double rho_fac = nd->nucdata[i].exm * (-conv) * nw->rho_deriv_fac.v;
#else
	  const double rho_fac = nd->nucdata[i].exm * (-conv) * nw->rho_deriv_fac.v + nw->drhodYi[i].v;
#endif /* NETWORK_NUCLEARNET_NEGLECT_DTDY_TERMS */
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP */

	  for(j = 0; j < nd->nuc_count; j++)
	    {
	      temp_row[j] += temp_fac * row[j];
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
	      rho_row[j] += rho_fac * row[j];
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP */
	    }
	}
#endif /* NETWORK_VARIABLE */

	/* now the derivatives of the rates w.r.t. the abundances
	 * (introduced by e.g. screening)
	 * this causes the Jacobian to be full instead of sparse
	 * This should be zero without screening.
	 * If NETWORK_SEP_YZ is defined, we neglect the Ye derivatives
	 * as Ye is only slowly changing.
	 * Ther abar and zbar derivatives are zero anyway.
	 */
#if defined(NETWORK_SCREENING) && !defined(NETWORK_SEP_YZ)
	{
	  network_var this_rhs = rhs[i];

	  for(k = 0; k < nd->nuc_count; k++)
	    {
	      const double dyedyi = (nd->nucdata[k].nz - nd->nucdata[k].na * ye);
#ifndef NETWORK_SEP_YZ
	      const double dyzdyi = gsl_pow_2(nd->nucdata[k].nz);
#endif
	      const double dzbardyi = (nd->nucdata[k].nz - zbar) * abar;

	      row[k] += this_rhs.dYe * dyedyi
#ifndef NETWORK_SEP_YZ
		+ this_rhs.dYz * dyzdyi
#endif
		+ this_rhs.dabar * dabardyi + this_rhs.dzbar * dzbardyi;
	    }
	}
#endif /* defined(NETWORK_SCREENING) && !defined(NETWORK_SEP_YZ) */

#ifdef NETWORK_SEP_YZ
	{
	  double *Yz_row = &matrix[nd->iYz * nd->n_matrix];
	  const double fac = gsl_pow_2(nd->nucdata[i].nz);

	  for(j = 0; j < nd->nuc_count; j++)
	    Yz_row[j] += fac * row[j];
	}

	/* drhs[y_i]/dYz */
	row[nd->iYz] = rhs[i].dYz;
#endif /* NETWORK_SEP_YZ */
#if NETWORK_VARIABLE
	/* drhs[y_i]/dT */
	row[nd->iTemp] = rhs[i].dT;
#endif /* NETWORK_VARIABLE */
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
	/* drhs[y_i]/drho */
	row[nd->iRho] = rhs[i].drho;
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP */

#if NETWORK_SPARSE
	store_csr(jacob, i, row, nd->n_matrix);
#endif /* NETWORK_SPARSE */
      }
  }

#ifdef NETWORK_SEP_YZ
  {
    /* drhs[Yz]/dyi have already been computed */
    double *row = &matrix[nd->iYz * nd->n_matrix];

    row[nd->iYz] = rhs[nd->iYz].dYz;
#if NETWORK_VARIABLE
    row[nd->iTemp] = rhs[nd->iYz].dT;
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
    row[nd->iRho] = rhs[nd->iYz].drho;
#endif
#endif

#if NETWORK_SPARSE
    store_csr(jacob, nd->iYz, row, nd->n_matrix);
#endif /* NETWORK_SPARSE */
  }
#endif

#if NETWORK_VARIABLE
  {
    /* drhs[T]/dyi */
    double *row = &matrix[nd->iTemp * nd->n_matrix];
    const network_var *rhsT = &rhs[nd->iTemp];
    int i;

    for(i = 0; i < nd->nuc_count; i++)
      {
	const double dyedyi = (nd->nucdata[i].nz - nd->nucdata[i].na * ye);
#ifndef NETWORK_SEP_YZ
	const double dyzdyi = gsl_pow_2(nd->nucdata[i].nz);
#endif
	const double dzbardyi = (nd->nucdata[i].nz - zbar) * abar;

	row[i] += rhsT->dYe * dyedyi
#ifndef NETWORK_SEP_YZ
	  + rhsT->dYz * dyzdyi
#endif
	  + rhsT->dabar * dabardyi + rhsT->dzbar * dzbardyi;
      }

#ifdef NETWORK_SEP_YZ
    /* drhs[T]/dYz */
    row[nd->iYz] = rhsT->dYz;
#endif
    /* drhs[T]/dT */
    row[nd->iTemp] = rhsT->dT;
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
    /* drhs[T]/drho */
    row[nd->iRho] = rhsT->drho;
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP */

#if NETWORK_SPARSE
    store_csr(jacob, nd->iTemp, row, nd->n_matrix);
#endif /* NETWORK_SPARSE */
  }
#endif /* NETWORK_VARIABLE */

#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  {
    /* drhs[rho]/dyi */
    double *row = &matrix[nd->iRho * nd->n_matrix];
    int i;

    for(i = 0; i < nd->nuc_count; i++)
      {
	const double dyedyi = (nd->nucdata[i].nz - nd->nucdata[i].na * ye);
#ifndef NETWORK_SEP_YZ
	const double dyzdyi = gsl_pow_2(nd->nucdata[i].nz);
#endif
	const double dzbardyi = (nd->nucdata[i].nz - zbar) * abar;
	const network_var *rhs_rho = &rhs[nd->iRho];

	row[i] += rhs_rho->dYe * dyedyi
#ifndef NETWORK_SEP_YZ
	  + rhs_rho->dYz * dyzdyi
#endif
	  + rhs_rho->dabar * dabardyi + rhs_rho->dzbar * dzbardyi;
      }

#ifdef NETWORK_SEP_YZ
    /* drhs[rho]/dYz */
    row[nd->iYz] = rhs[nd->iRho].dYz;
#endif
    /* drhs[rho]/dT */
    row[nd->iTemp] = rhs[nd->iRho].dT;
    /* drhs[rho]/drho */
    row[nd->iRho] = rhs[nd->iRho].drho;

#if NETWORK_SPARSE
    store_csr(jacob, nd->iRho, row, nd->n_matrix);
#endif /* NETWORK_SPARSE */
  }
#endif /* NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP */

#if NETWORK_SPARSE
  {
#if defined(NETWORK_SUPERLU)
    const int offset = 0;	/* C-style indices */
#elif defined(NETWORK_PARDISO)
    const int offset = 1;	/* Fortran-style indices */
#else
#error "unknown sparse matrix library"
#endif

    /* record the final size of the sparse matrix */
    jacob->rowstart[nd->n_matrix] = jacob->n_elements + offset;
#if DEBUG > 2
    myprintf("sparseness: %d / %zd (%g %%)\n", jacob->rowstart[nd->n_matrix], nd->n_matrix * nd->n_matrix,
	     100.0 * (double) jacob->rowstart[nd->n_matrix] / (double) (nd->n_matrix * nd->n_matrix));
#endif /* DEBUG */
  }
#endif /* NETWORK_SPARSE */

  return 0;
}

int getrates(double rho, double temp, double ye, double yz, int compute_derivs, const struct network_data *nd,
	     struct network_workspace *nw)
{
  double t9 = temp * 1e-9;
  network_var temp9[7] = { {0} };
  int i, j, changed;

  /* do not calculate rates again if nothing changed */
  if(nw->old_temp == temp && nw->old_rho == rho && nw->old_ye == ye)
    return 0;

  changed = 0;

  if(nw->old_temp != temp)
    {
      double t9l, tempi;

      network_part(temp, nd, nw);

      t9l = log(t9);
      tempi = 1.0 / temp;

      temp9[0].v = 1.0;
      temp9[0].dT = 0.0;
      temp9[1].v = 1.0 / t9;
      temp9[1].dT = -temp9[1].v * tempi;
      temp9[2].v = exp(-t9l / 3.0);
      temp9[2].dT = (-1.0 / 3.0) * temp9[2].v * tempi;
      temp9[3].v = 1.0 / temp9[2].v;
      temp9[3].dT = (1.0 / 3.0) * temp9[3].v * tempi;
      temp9[4].v = t9;
      temp9[4].dT = 1.0e-9;
      temp9[5].v = gsl_pow_5(temp9[3].v);	/* exp( t9l * 5.0 / 3.0 ) */
      temp9[5].dT = (5.0 / 3.0) * temp9[5].v * tempi;
      temp9[6].v = t9l;
      temp9[6].dT = tempi;

      for(i = 0; i < nd->rate_count; i++)
	{
	  struct network_rate *rate = &nd->rates[i];
	  double *data = rate->data;
	  network_var baserate = { 0 };

	  for(j = 0; j < 7; j++)
	    {
	      baserate.v += temp9[j].v * data[j];
	      baserate.dT += temp9[j].dT * data[j];
	    }
	  baserate = n_exp(baserate);

	  /* account for reverse rate */
	  if(rate->isReverse)
	    {
	      /* divide by input and multiply with output
	       * only normalised temperature dependent partition functions are used here
	       * as the groundstate is already part of the rate parametrization (included in a0) */
	      for(j = 0; j < rate->ninput; j++)
		{
		  baserate = n_d(baserate, nw->gg[nd->rates[i].input[j]]);
		}
	      for(j = 0; j < nd->rates[i].noutput; j++)
		{
		  baserate = n_m(baserate, nw->gg[nd->rates[i].output[j]]);
		}
	    }

	  nw->baserate[i] = baserate;
	}

      nw->old_temp = temp;
      changed = 1;
    }

  if(nw->old_rho != rho || changed)
    {
      const network_var rho_v = {.v = rho,.drho = 1.0 };
      const network_var *baserate = nw->baserate;

      for(i = 0; i < nd->rate_count; i++, baserate++)
	{
	  network_var *rate = &nw->unscreened_rate[i];
	  int n = nd->rates[i].ninput - 1;

	  *rate = *baserate;
	  /* multiply with rho (n-1)-times for n-body rate */
	  for(j = 0; j < n; j++)
	    *rate = n_m(*rate, rho_v);
	}

      nw->old_rho = rho;
      changed = 1;
    }

  if(nw->old_ye != ye || changed)
    {
      const network_var ye_v = {.v = ye,.dYe = 1.0 };

      for(i = 0; i < nd->electronCapture_count; i++)
	{
	  nw->unscreened_rate[nd->electronCapture_rates[i]] =
	    n_m(nw->unscreened_rate[nd->electronCapture_rates[i]], ye_v);
	}

      nw->old_ye = ye;
      changed = 1;
    }

  /* if something has changed, apply the screening and recalculate weak rates */
  if(changed)
    {
      memcpy(nw->rates, nw->unscreened_rate, nd->rate_count * sizeof(network_var));

#ifdef NETWORK_SCREENING
      screening_apply(rho, temp, ye, yz, nd, nw);
#endif /* NETWORK_SCREENING */

      /* compute weak rates
       * only when there actually are weak rates in the network */
      if(nd->weakrate_count > 0)
	{
	  network_var logrhoye = {.v = log10(rho * ye),.drho = 1.0 / (rho * log(10)),.dYe =
	      1.0 / (ye * log(10)) };
	  network_var dt = { 0 }
	  , drhoye =
	  {
	  0};
	  network_var b1, b2, rate1, rate2;
	  size_t iTempLow, iTempHigh, iRhoYeLow, iRhoYeHigh;
	  size_t idx1, idx2, idx3, idx4;
	  double at, arhoye;

	  iTempLow = 0;
	  while(t9 > nd->weakTemp[iTempLow])
	    iTempLow++;

	  if(iTempLow > 0)
	    iTempLow--;

	  iTempHigh = min(iTempLow + 1, 12);

	  iRhoYeLow = 0;
	  while(logrhoye.v > nd->weakRhoYe[iRhoYeLow])
	    iRhoYeLow++;

	  if(iRhoYeLow > 0)
	    iRhoYeLow--;

	  iRhoYeHigh = min(iRhoYeLow + 1, 10);

	  idx1 = iRhoYeLow * 13 + iTempLow;
	  idx2 = iRhoYeHigh * 13 + iTempLow;
	  idx3 = iRhoYeLow * 13 + iTempHigh;
	  idx4 = iRhoYeHigh * 13 + iTempHigh;

	  dt.v = t9 - nd->weakTemp[iTempLow];
	  dt.dT = 1e-9;
	  drhoye = logrhoye;
	  drhoye.v = logrhoye.v - nd->weakRhoYe[iRhoYeLow];
	  at = nd->weakTemp[iTempHigh] - nd->weakTemp[iTempLow];
	  arhoye = nd->weakRhoYe[iRhoYeHigh] - nd->weakRhoYe[iRhoYeLow];

	  for(i = 0; i < nd->weakrate_count; i++)
	    {
	      b1 =
		n_ca(n_cm(dt, (nd->weakrates[i].lambda1[idx3] - nd->weakrates[i].lambda1[idx1]) / at),
		     nd->weakrates[i].lambda1[idx1]);
	      b2 =
		n_ca(n_cm(dt, (nd->weakrates[i].lambda1[idx4] - nd->weakrates[i].lambda1[idx2]) / at),
		     nd->weakrates[i].lambda1[idx2]);
	      rate1 = n_a(b1, n_cm(n_m(n_s(b2, b1), drhoye), 1.0 / arhoye));

	      b1 =
		n_ca(n_cm(dt, (nd->weakrates[i].lambda2[idx3] - nd->weakrates[i].lambda2[idx1]) / at),
		     nd->weakrates[i].lambda2[idx1]);
	      b2 =
		n_ca(n_cm(dt, (nd->weakrates[i].lambda2[idx4] - nd->weakrates[i].lambda2[idx2]) / at),
		     nd->weakrates[i].lambda2[idx2]);
	      rate2 = n_a(b1, n_cm(n_m(n_s(b2, b1), drhoye), 1.0 / arhoye));

	      nw->weakrates[i] = n_a(n_exp(n_cm(rate1, log(10))), n_exp(n_cm(rate2, log(10))));
	    }
	}
    }

  return 0;
}

int network_part(double temp, const struct network_data *nd, struct network_workspace *nw)
{
  /* interpolates partition functions, given the temperature */
  /* to do: implement partition functions for T > 1e10 K (cf. Rauscher paper?) */
  int index, i;
  double tempLeft, tempRight;
  double dlgLeft, dlgRight;
  double grad;

  static const double network_parttemp[24] = { 1.0e8, 1.5e8, 2.0e8, 3.0e8, 4.0e8, 5.0e8, 6.0e8, 7.0e8,
    8.0e8, 9.0e8, 1.0e9, 1.5e9, 2.0e9, 2.5e9, 3.0e9, 3.5e9,
    4.0e9, 4.5e9, 5.0e9, 6.0e9, 7.0e9, 8.0e9, 9.0e9, 1.0e10
  };

  index = 0;
  temp = min(max(temp, network_parttemp[0]), network_parttemp[23]);

  while(temp > network_parttemp[index])
    {
      index++;
    }
  if(index > 0)
    index--;

  tempLeft = network_parttemp[index];
  tempRight = network_parttemp[index + 1];

  for(i = 0; i < nd->nuc_count; i++)
    {
      dlgLeft = nd->nucdata[i].part[index];
      dlgRight = nd->nucdata[i].part[index + 1];

      grad = (dlgRight - dlgLeft) / (tempRight - tempLeft);
      set_zero(nw->gg[i]);
      nw->gg[i].v = exp(dlgLeft + (temp - tempLeft) * grad);
      nw->gg[i].dT = nw->gg[i].v * grad;
    }

  return 0;
}

static void azbar(const double y[], double *abar, double *zbar, const struct network_data *nd)
{
  int i;
  /* compute abar, zbar */
  *abar = 0.0;
  *zbar = 0.0;
  for(i = 0; i < nd->nuc_count; i++)
    {
      *abar += y[i];
      *zbar += nd->nucdata[i].nz * y[i];
    }
  *abar = 1.0 / *abar;
  *zbar = *zbar * *abar;
}

#if NETWORK_SPARSE
static void store_csr(jacob_t * jacob, int row_ind, const double *row, size_t n_row)
{
  /* row_ind is always the C-style row index */
  int col_ind;
  const double *row_p;
  double *value_p;
  int *columns_p;
#if defined(NETWORK_SUPERLU)
  const int offset = 0;		/* C-style indices */
#elif defined(NETWORK_PARDISO)
  const int offset = 1;		/* Fortran-style indices */
#else
#error "unknown sparse matrix library"
#endif

  jacob->rowstart[row_ind] = jacob->n_elements + offset;

  value_p = &jacob->values[jacob->n_elements];
  columns_p = &jacob->columns[jacob->n_elements];
  for(col_ind = offset, row_p = row; col_ind < n_row + offset; col_ind++, row_p++)
    {
      if(fabs(*row_p) > NETWORK_SPARSE_THRESHOLD || row_ind + offset == col_ind)
	{
	  /* we include diagonal values even if they are zero
	   * because 1 is added to these values later in prepare_matrix */
	  *value_p = *row_p;
	  *columns_p = col_ind;
	  value_p++;
	  columns_p++;
	  jacob->n_elements++;
	}
    }
}
#endif /* NETWORK_SPARSE */

#ifdef NETWORK_SCREENING
static void compute_constant_screening_factors(struct network_rate *rate, const struct network_nucdata *nuc1,
					       const struct network_nucdata *nuc2,
					       const struct network_data *nd)
{
  /* If nuc1 or nuc2 != NULL, we use the this nucleus instead of the one specified in rate.
   * This is needed for triple alpha screening, as Be-8 is not an element included in the network */
  double z1, z2, a1, a2;

  if(nuc1 != NULL)
    {
      z1 = nuc1->nz;
      a1 = nuc1->na;
    }
  else
    {
      z1 = nd->nucdata[rate->input[0]].nz;
      a1 = nd->nucdata[rate->input[0]].na;
    }
  if(nuc2 != NULL)
    {
      z2 = nuc2->nz;
      a2 = nuc2->na;
    }
  else
    {
      z2 = nd->nucdata[rate->input[1]].nz;
      a2 = nd->nucdata[rate->input[1]].na;
    }

  rate->z_tilde1 = pow(z1 + z2, 5. / 3.) - pow(z1, 5. / 3.) - pow(z2, 5. / 3.);
  rate->z_tilde2 = pow(z1 + z2, 5. / 12.) - pow(z1, 5. / 12.) - pow(z2, 5. / 12.);

  rate->C_const = -0.5551 * 5. / 3. * log(z1 * z2 / (z1 + z2)) - 2.996;
  rate->tau_fac =
    pow(27. / 2. * gsl_pow_2(M_PI) * (a1 * a2) / (a1 + a2) * gsl_pow_2(z1 * z2) *
	gsl_pow_4(ELECTRON_CHARGE_ESU) / (GSL_CONST_NUM_AVOGADRO * GSL_CONST_CGS_BOLTZMANN *
					  gsl_pow_2(GSL_CONST_CGS_PLANCKS_CONSTANT_HBAR)), 1. / 3.);
  /* Note that hbar is correct here. The h in (A2) is a misprint. */
  /* This is taken from ApJ 181:439 DeWitt et al.
   * The formula from Wallace is missing a sqrt
   */
  rate->H_w_fac =
    gsl_pow_3(ELECTRON_CHARGE_ESU) * sqrt(4.0 * M_PI * GSL_CONST_NUM_AVOGADRO /
					  gsl_pow_3(GSL_CONST_CGS_BOLTZMANN)) * z1 * z2;
  rate->Gamma_fac = pow(2. / (z1 + z2), 1. / 3.) * z1 * z2;
}

void screening_apply(double rho, double temp, double ye, double yz, const struct network_data *nd,
		     struct network_workspace *nw)
{
  network_var H, Gamma_prime, lGamma_prime, Gamma_prime14, H_w_common;
  network_var itemp13 = { 0 }, ne;
  const network_var temp_var = {.v = temp,.dT = 1.0 }, rho_var =
  {
  .v = rho,.drho = 1.0}, ye_var =
  {
  .v = ye,.dYe = 1.0}, yz_var =
  {
  .v = yz,.dYz = 1.0};
  struct network_rate **rate, **rateend;

  /* This implements enhancement of the rates by screening as explained in
   * the appendix of Wallace, Woosley & Weaver 1982, ApJ 258
   * Note that there is an erratum ApJ 264:746 that corrects a coefficient
   */

  itemp13.v = 1.0 / cbrt(temp);
  itemp13.dT = (-1.0 / 3.0) * itemp13.v / temp;

  ne = n_cm(n_m(ye_var, rho_var), 1.0 / GSL_CONST_CGS_UNIFIED_ATOMIC_MASS);

  Gamma_prime =
    n_cm(n_d(n_cbrt(ne), temp_var),
	 gsl_pow_2(ELECTRON_CHARGE_ESU) / GSL_CONST_CGS_BOLTZMANN * cbrt(4. / 3. * M_PI));
  lGamma_prime = n_log(Gamma_prime);
  Gamma_prime14 = n_sqrt(n_sqrt(Gamma_prime));

  /* cf. ApJ 181:439 DeWitt et al. 1973
   * there is a square root missing in Wallace et al. 1982
   * theta_e is set to 1 following advice from Ivo Seitenzahl
   * it is inconsequential for large Z
   *
   * This common factor is only computed once for all rates.
   */
  H_w_common = n_sqrt(n_d(n_m(n_a(yz_var, ye_var), rho_var), n_powi(temp_var, 3)));

  rateend = &nd->screened_rates[nd->screened_count];
  for(rate = &nd->screened_rates[0]; rate != rateend; rate++)
    {
      int ind = (*rate)->index;

      H = compute_screening(*rate, itemp13, Gamma_prime, lGamma_prime, Gamma_prime14, H_w_common);

      nw->rates[ind] = n_m(nw->unscreened_rate[ind], n_exp(H));

    }
  if(nd->tripleAlpha_rate != NULL)
    {
      int ind = nd->tripleAlpha_rate->index;

      H =
	n_a(compute_screening
	    (&nd->alpha_alpha_rate, itemp13, Gamma_prime, lGamma_prime, Gamma_prime14, H_w_common),
	    compute_screening(&nd->Be8_alpha_rate, itemp13, Gamma_prime, lGamma_prime, Gamma_prime14,
			      H_w_common));

      nw->rates[ind] = n_m(nw->unscreened_rate[ind], n_exp(H));
    }
}

static network_var compute_screening(const struct network_rate *rate, network_var itemp13,
				     network_var Gamma_prime, network_var lGamma_prime,
				     network_var Gamma_prime14, network_var H_w_common)
{
  network_var Gamma, H, H_s, H_w;

  Gamma = n_cm(Gamma_prime, rate->Gamma_fac);

  if(Gamma.v >= 0.3)
    {
      /* not only weak */
      network_var Cvar, b, tau;

      Cvar = n_ca(n_a(n_a(n_cm(Gamma_prime, 0.896434 * rate->z_tilde1),
		       n_cm(Gamma_prime14, -3.44740 * rate->z_tilde2)),
		   n_cm(lGamma_prime, -0.5551)), rate->C_const);
      tau = n_cm(itemp13, rate->tau_fac);
      b = n_cm(n_d(Gamma, tau), 3.);

      /* note that there is an error in the coefficient of b^5 in the original paper */
      H_s = n_a(Cvar,
		n_a(n_a(n_a(n_m(n_powi(b, 3),
				n_cm(tau, -5. / (32. * 3.))),
			    n_m(n_powi(b, 4),
				n_s(n_cm(tau, 0.014 / 3.),
				    n_cm(Gamma, 0.0055)))),
			n_m(n_powi(b, 5),
			    n_a(n_cm(tau, 0.0128 / 3.),
				n_cm(Gamma, 0.0098)))), n_m(n_powi(b, 6), n_cm(Gamma, 0.0048))));
    }

  if(Gamma.v <= 0.8)
    {
      /* not only strong */
      H_w = n_cm(H_w_common, rate->H_w_fac);
    }

  if(Gamma.v < 0.3)
    H = H_w;
  else if(Gamma.v > 0.8)
    H = H_s;
  /* this linear interpolation between the two regimes ensures continuous behavior
   * the original suggestion by Wallace et al. 1982 was bad when the values differ by orders of magnitude
   * this was suggested in the erratum
   */
  else
    H = n_a(n_m(H_w, n_cm(n_ca(Gamma, -0.8), -2.0)), n_m(H_s, n_cm(n_ca(Gamma, -0.3), 2.0)));


  /* If H is implausibly large, this is probably due to the timestep being too large */
  if(H.v >= 30.0)
    {
      network_var zero_var = { 0 };
#if DEBUG > 3
      if(rate != &Be8_alpha_rate)
          PRINT_WARNING("H is %g for Z_1 = %d, Z_2 = %d\n", H.v, nd->nucdata[rate->input[0]].nz,nd->nucdata[rate->input[1]].nz);
      else
          PRINT_WARNING("H is %g for Be-8, He-4\n", H.v);
#endif /* DEBUG */
      H = zero_var;
    }
  assert(H.v >= 0.0);

  return H;
}

#endif /* NETWORK_SCREENING */




/* some really basic utility functions used in the network subroutines */


void myprintf(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}

char *util_fgets(char *str, int num, FILE * stream, char *file, int line)
{
    char *ret = fgets(str, num, stream);
    if(ret == NULL)
    {
        printf("error: fgets in file %s at line %d\n", file, line); endrun(200);
    }
    return ret;
}

size_t util_fread(void *ptr, size_t size, size_t count, FILE * stream, char *file, int line)
{
    size_t result = fread(ptr, size, count, stream);
    if(result != count)
    {
        printf("error: fread in file %s at line %d\n", file, line); endrun(201);
    }
    return result;
}

double SwapDouble(double Val)
{
    double nVal;
    int i;
    const char *readFrom = (const char *) &Val;
    char *writeTo = ((char *) &nVal) + sizeof(nVal);
    
    for(i = 0; i < sizeof(Val); ++i)
    {
        *(--writeTo) = *(readFrom++);
    }
    return nVal;
}

float SwapFloat(float Val)
{
    float nVal;
    int i;
    const char *readFrom = (const char *) &Val;
    char *writeTo = ((char *) &nVal) + sizeof(nVal);
    
    for(i = 0; i < sizeof(Val); ++i)
    {
        *(--writeTo) = *(readFrom++);
    }
    return nVal;
}

int SwapInt(int Val)
{
    int nVal;
    int i;
    const char *readFrom = (const char *) &Val;
    char *writeTo = ((char *) &nVal) + sizeof(nVal);
    
    for(i = 0; i < sizeof(Val); ++i)
    {
        *(--writeTo) = *(readFrom++);
    }
    return nVal;
}

int CheckSwap(char *fname, int *swap)
{
    FILE *fd;
    off_t fsize, fpos;
    int blocksize, blockend;
    
    if(!(fd = fopen(fname, "r")))
    {
        printf("can't open file `%s'.\n", fname);
        return -1;
    }
    
    fseeko(fd, 0, SEEK_END);
    fsize = ftello(fd);
    
    *swap = 0;
    fpos = 0;
    fseeko(fd, 0, SEEK_SET);
    safe_fread(&blocksize, sizeof(int), 1, fd);
    while(!feof(fd))
    {
        if(fpos + blocksize + 4 > fsize)
        {
            *swap += 1;
            break;
        }
        fpos += 4 + blocksize;
        fseeko(fd, fpos, SEEK_SET);
        safe_fread(&blockend, sizeof(int), 1, fd);
        if(blocksize != blockend)
        {
            *swap += 1;
            break;
        }
        fpos += 4;
        if(!fread(&blocksize, sizeof(int), 1, fd))
            break;
    }
    
    if(*swap == 0)
    {
        fclose(fd);
        return 0;
    }
    
    fpos = 0;
    fseeko(fd, 0, SEEK_SET);
    safe_fread(&blocksize, sizeof(int), 1, fd);
    while(!feof(fd))
    {
        blocksize = SwapInt(blocksize);
        if(fpos + blocksize + 4 > fsize)
        {
            *swap += 1;
            break;
        }
        fpos += 4 + blocksize;
        fseeko(fd, fpos, SEEK_SET);
        safe_fread(&blockend, sizeof(int), 1, fd);
        blockend = SwapInt(blockend);
        if(blocksize != blockend)
        {
            *swap += 1;
            break;
        }
        fpos += 4;
        if(!fread(&blocksize, sizeof(int), 1, fd))
            break;
    }
    
    fclose(fd);
    return 0;
}








#endif /* NUCLEAR_NETWORK */
