#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_const_cgs.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"


/*
 *  This code is place-holder, inherited from GADGET3
 */

#ifdef NUCLEAR_NETWORK
//#include "./nuclear_network.h"

#ifdef NETWORK_SUPERLU
#include "slu_ddefs.h"
#endif

#ifdef NETWORK_PARDISO
#include "mkl_pardiso.h"
#include "mkl_service.h"
#endif

#include "../eos/eos.h"

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

static void
  __attribute__ ((unused)) print_matrix(const int n, const double *matrix, const int *rowstart,
					const int *columns, const double *rhs, const int offset);

static void network_solver_integrate_internal(double temp, double rho, double *y, double maxtime,
					      const struct network_data *nd, struct network_workspace *nw,
					      struct network_solver_trajectory *traj);
static void network_solver_interpolate_traj_rhotemp(struct network_solver_trajectory *traj, double time,
						    double *y, double *rho, double *temp);
static void network_solver_calc_dy(double temp, double rho, double dt, int steps, double *y, double *yn,
				   const double *first_rhs, const jacob_t * old_jacob,
				   const struct network_data *nd, struct network_workspace *nw,
				   struct network_solver_trajectory *traj);
static void network_solver_extrapolate(int iter, int matrixsize, double *x, double *y, double *dy,
				       double *qcol);
static void normalize(double *y, const struct network_solver_data *nsd);
static void sanitize(double *y, const struct network_data *nd);
static void prepare_matrix(jacob_t * mod_jacob, const jacob_t * old_jacob, size_t n_matrix, double h);

struct network_solver_data *network_solver_init(double tolerance, int matrixsize, int nelements,
						const struct network_nucdata *nucdata)
{
  int i, k;
  struct network_solver_data *nsd;

  nsd = malloc(sizeof(struct network_solver_data));

  nsd->tolerance = tolerance;
  nsd->matrixsize = matrixsize;
  nsd->matrixsize2 = matrixsize * matrixsize;
  nsd->nelements = nelements;

  nsd->nsteps = 12;
  nsd->maxstep = nsd->nsteps - 1;
  nsd->steps = (int *) malloc(nsd->nsteps * sizeof(int));
  /* These steps fulfill the Toeplitz condition for alpha = 5/7
   * cf. Bader & Deuflhard 1983 eq. (3.8) */
  nsd->steps[0] = 2;
  nsd->steps[1] = 6;
  nsd->steps[2] = 10;
  nsd->steps[3] = 14;
  nsd->steps[4] = 22;
  nsd->steps[5] = 34;
  nsd->steps[6] = 50;
  nsd->steps[7] = 70;
  nsd->steps[8] = 98;
  nsd->steps[9] = 138;
  nsd->steps[10] = 194;
  nsd->steps[11] = 272;

  nsd->aion = malloc(nelements * sizeof(double));
  nsd->dy = (double *) malloc(nsd->matrixsize * sizeof(double));
  nsd->ynew = (double *) malloc(nsd->matrixsize * sizeof(double));
  nsd->yscale = (double *) malloc(nsd->matrixsize * sizeof(double));
  nsd->rhs_deriv = malloc(nsd->matrixsize * sizeof(network_var));
  nsd->first_rhs = malloc(nsd->matrixsize * sizeof(double));
  nsd->rhs = malloc(nsd->matrixsize * sizeof(double));
  nsd->x = (double *) malloc(nsd->maxstep * sizeof(double));
  nsd->err = (double *) malloc(nsd->maxstep * sizeof(double));
  nsd->qcol = (double *) malloc(nsd->matrixsize * nsd->nsteps * sizeof(double));
  nsd->a = (double *) malloc(nsd->nsteps * sizeof(double));
  nsd->alf = (double *) malloc((nsd->nsteps - 1) * (nsd->nsteps - 1) * sizeof(double));

  for(i = 0; i < nelements; i++)
    nsd->aion[i] = nucdata[i].na;

#if NETWORK_SPARSE
  nsd->jacob.values = malloc(nsd->matrixsize2 * sizeof(double));
  nsd->jacob.columns = malloc(nsd->matrixsize2 * sizeof(int));
  nsd->jacob.rowstart = malloc((nsd->matrixsize + 1) * sizeof(int));
  nsd->mod_jacob.values = malloc(nsd->matrixsize2 * sizeof(double));
  nsd->mod_jacob.columns = malloc(nsd->matrixsize2 * sizeof(int));
  nsd->mod_jacob.rowstart = malloc((nsd->matrixsize + 1) * sizeof(int));
#else
  nsd->jacob.matrix = malloc(nsd->matrixsize2 * sizeof(double));
  nsd->mod_jacob.matrix = malloc(nsd->matrixsize2 * sizeof(double));
#endif /* NETWORK_SPARSE */

  /* performance checker */
  nsd->a[0] = nsd->steps[0] + 1.0;
  /* Deuflhard 1983 eq. (2.11) */
  for(i = 0; i < nsd->maxstep; i++)
    nsd->a[i + 1] = nsd->a[i] + nsd->steps[i + 1];

  for(i = 1; i < nsd->maxstep; i++)
    {
      for(k = 0; k < i; k++)
	{
	  nsd->alf[i * nsd->maxstep + k] = pow(nsd->tolerance * 0.25, (nsd->a[k + 1] - nsd->a[i + 1]) /
					       ((nsd->a[i + 1] - nsd->a[0] + 1.0) * (2. * (k + 1.0) + 1.0)));
	}
    }

  nsd->a[0] += nsd->matrixsize;
  for(i = 0; i < nsd->maxstep; i++)
    nsd->a[i + 1] = nsd->a[i] + nsd->steps[i + 1];

  for(i = 1; i < nsd->maxstep - 1; i++)
    {
      if(nsd->a[i + 1] > nsd->a[i] * nsd->alf[i * nsd->maxstep + i - 1])
	break;
    }
  nsd->maxiter = i + 1;

  myprintf("Network solver initialization done\n");

  return nsd;
}

void network_solver_deinit(struct network_solver_data *nsd)
{
  free(nsd->steps);
  free(nsd->aion);
  free(nsd->dy);
  free(nsd->ynew);
  free(nsd->yscale);
  free(nsd->rhs_deriv);
  free(nsd->rhs);
  free(nsd->first_rhs);
  free(nsd->x);
  free(nsd->err);
  free(nsd->qcol);
  free(nsd->a);
  free(nsd->alf);

#if NETWORK_SPARSE
  free(nsd->jacob.values);
  free(nsd->jacob.columns);
  free(nsd->jacob.rowstart);
  free(nsd->mod_jacob.values);
  free(nsd->mod_jacob.columns);
  free(nsd->mod_jacob.rowstart);
#else
  free(nsd->jacob.matrix);
  free(nsd->mod_jacob.matrix);
#endif /* NETWORK_SPARSE */

  free(nsd);
}

void network_solver_interpolate_trajectory(struct network_solver_trajectory *traj, double time, double *rho,
					   double *energy)
{
  integertime timestep;
  double dt, tstep;

  timestep = traj->timestep;

  if(timestep > 0 && traj->timesteps[timestep - 1] > time)
    timestep = 0;

  while(time >= traj->timesteps[timestep] && timestep < traj->ntimesteps - 1)
    timestep++;

  traj->timestep = timestep;

  if(traj->timestep == 0)
    {
      *rho = traj->rho[0];
      *energy = traj->energy[0];
      return;
    }

  dt = time - traj->timesteps[timestep - 1];
  tstep = traj->timesteps[timestep] - traj->timesteps[timestep - 1];

  *rho = traj->rho[timestep - 1] + (traj->rho[timestep] - traj->rho[timestep - 1]) * dt / tstep;
  *energy = traj->energy[timestep - 1] + (traj->energy[timestep] - traj->energy[timestep - 1]) * dt / tstep;

  return;
}

static void network_solver_interpolate_traj_rhotemp(struct network_solver_trajectory *traj, double time,
						    double *y, double *rho, double *temp)
{
  double energy;
  network_solver_interpolate_trajectory(traj, time, rho, &energy);
  //eos_calc_egiven_y(*rho, y, energy, temp, NULL);
}

void network_solver_integrate_traj(const struct network_data *nd, struct network_workspace *nw,
				   struct network_solver_trajectory *traj)
{
  int i;

#ifdef NETWORK_VARIABLE
  myprintf("Integrate trajectory does not work with NETWORK_VARIABLE-switch.\n");
  endrun(70);
#endif

  for(i = 0; i < nd->nuc_count; i++)
    traj->x[i] /= nd->nucdata[i].na;
  network_solver_integrate_internal(0, 0, traj->x, 0, nd, nw, traj);
  for(i = 0; i < nd->nuc_count; i++)
    traj->x[i] *= nd->nucdata[i].na;
}

void network_solver_integrate(double temp, double rho, double *y, double dt, const struct network_data *nd,
			      struct network_workspace *nw)
{
  network_solver_integrate_internal(temp, rho, y, dt, nd, nw, NULL);
}

static void network_solver_integrate_internal(double temp, double rho, double *y, double maxtime,
					      const struct network_data *nd, struct network_workspace *nw,
					      struct network_solver_trajectory *traj)
{
  double *dy, *ynew, *yscale;
  double *x, *err, *qcol;
  double *a, *alf;
  double time, dttry, dtnext;
  double errmax, dum, red;
  double workmin, work, fac, scale;
  int i, m, iter, lastiter, iterm, iteropt;
  int first, reduce, step;
  double *first_rhs;
  network_var *rhs_deriv;
  jacob_t *jacob;
  struct network_solver_data *nsd = nw->nsd;
  double start;
#ifdef NUCLEARNET_OUTPUT_TIMEEVOLUTION
  FILE *fp;

  fp = fopen("network_output.dat", "w");
  if(fp == NULL)
    {
      perror("error opening file `network_output.dat'");
      endrun(1001);
    }
  fwrite(&nsd->nelements, sizeof(int), 1, fp);
  fwrite(&nsd->matrixsize, sizeof(int), 1, fp);
  for(i = 0; i < nsd->nelements; i++)
    fwrite(&nsd->aion[i], sizeof(double), 1, fp);
#endif

  /* to suppress compiler warnings that wrongly indicate uninitialized variables */
  red = 0.0;
  iterm = 0.0;

  dy = nsd->dy;
  ynew = nsd->ynew;
  yscale = nsd->yscale;
  x = nsd->x;
  err = nsd->err;
  qcol = nsd->qcol;
  a = nsd->a;
  alf = nsd->alf;

  first_rhs = nsd->first_rhs;
  rhs_deriv = nsd->rhs_deriv;
  jacob = &nsd->jacob;

  if(traj)
    {
      time = traj->time;
      maxtime = traj->maxtime;
    }
  else
    {
      time = 0;
    }
  m = 0;
  step = 0;
  first = 1;
  dttry = maxtime - time;
  iteropt = nsd->maxiter - 1;

  if(traj && dttry > 1e-3)
    dttry = 1e-3;

  while(time < maxtime)
    {
      start = (double) clock();

#if defined (DEBUG) && defined (NETWORK_SEP_YZ)
      {
	int i;
	double yz = 0.0;

	for(i = 0; i < nd->nuc_count; i++)
	  yz += y[i] * gsl_pow_2(nd->nucdata[i].nz);
	myprintf("relative deviation between yz and real yz: %g\n", (y[nd->iYz] - yz) / yz);
      }
#endif

      sanitize(y, nd);
      normalize(y, nsd);

#ifdef NUCLEARNET_OUTPUT_TIMEEVOLUTION
      printf("t(%03d): %11.5e, dt (%02d): %11.5e, x:", step, time, m, dttry);
      for(i = 0; i < nsd->nelements; i++)
	{
	  printf(" %8.1e", y[i] * nsd->aion[i]);
	}
      if(nsd->nelements < nsd->matrixsize)
	{
	  printf(" %8.1e", y[nsd->nelements]);
	}
      printf("\n");
      fwrite(&time, sizeof(double), 1, fp);
      fwrite(y, sizeof(double), nsd->matrixsize, fp);
#endif

      if(traj)
	{
	  network_solver_interpolate_traj_rhotemp(traj, time, y, &rho, &temp);

	  if((traj->mintemp && temp < traj->mintemp) || (traj->maxtemp && temp > traj->maxtemp))
	    {
	      break;
	    }
	}

      for(i = 0; i < nsd->matrixsize; i++)
	{
	  yscale[i] = max(fabs(y[i]), nsd->tolerance);
	}

      /* compute the Jacobian only once per successful step and use prepare_matrix
         for the different factors later */

      network_getrhs(rho, temp, y, 1, nd, nw, first_rhs, rhs_deriv);
      network_getjacob(y, rhs_deriv, nd, nw, jacob);

      reduce = 0;
      for(iter = 0; iter < nsd->maxiter; iter++)
	{
	  if(dttry == 0)
	    {
	      myprintf("dt is zero.\n");
	      return;
	    }

	  if(reduce == 1)
	    reduce = 2;

	  network_solver_calc_dy(temp, rho, dttry, nsd->steps[iter], y, ynew, first_rhs, jacob, nd, nw, traj);
	  x[iter] = gsl_pow_2(dttry / nsd->steps[iter]);
	  network_solver_extrapolate(iter, nsd->matrixsize, x, ynew, dy, qcol);

	  /* compute normalized error estimate */
	  if(iter > 0)
	    {
	      errmax = 1.0e-30;
	      for(i = 0; i < nsd->matrixsize; i++)
		{
		  dum = fabs(dy[i] / yscale[i]);
		  if(dum > errmax)
		    {
		      errmax = dum;
		      m = i;
		    }
		}

	      errmax /= nsd->tolerance;
	      iterm = iter - 1;
	      /* cf. Deuflhard 1983 eq. (2.12)
	       * 0.25 is a safety factor (cf. Timmes' torch Network) */
	      err[iterm] = pow(errmax / 0.25, 1.0 / (2. * (iterm + 1.) + 1.));
	    }

	  if(iter > 0 && (iter >= iteropt - 1 || first))
	    {
	      /* if converged, leave */
	      if(errmax < 1.0)
		break;

	      if(iter == nsd->maxiter - 1 || iter == iteropt + 1)
		{
		  red = 0.7 / err[iterm];
		  reduce = 1;
		}
	      else if(iter == iteropt)
		{
		  if(alf[iteropt * nsd->maxstep + iteropt - 1] < err[iterm])
		    {
		      red = 1.0 / err[iterm];
		      reduce = 1;
		    }
		}
	      else if(iteropt == nsd->maxiter - 1)
		{
		  if(alf[(nsd->maxiter - 2) * nsd->maxstep + iterm] < err[iterm])
		    {
		      red = alf[(nsd->maxiter - 2) * nsd->maxstep + iterm] * 0.7 / err[iterm];
		      reduce = 1;
		    }
		}
	      else if(alf[iteropt * nsd->maxstep + iterm] < err[iterm])
		{
		  red = alf[(iteropt - 1) * nsd->maxstep + iterm] / err[iterm];
		  reduce = 1;
		}

	      /* reduce step size */
	      if(reduce == 1)
		{
		  /* step not successful */
		  red = red > 1e-5 ? red : 1e-5;
		  red = red < 0.7 ? red : 0.7;
		  dttry = dttry * red;
		  iter = -1;
		}
	    }
	}

      lastiter = iter;
      /* step successful */
      for(i = 0; i < nsd->matrixsize; i++)
	{
	  y[i] = ynew[i];
	}

      first = 0;

      /* get a new estimate for the optimal number of iterations and corresponding step size */
      scale = 1.0;
      workmin = 1.0e35;
      /* Deuflhard 1983 eq. (2.13) */
      for(iter = 0; iter < iterm + 1; iter++)
	{
	  fac = err[iter] > 0.1 ? err[iter] : 0.1;
	  work = fac * a[iter + 1];
	  if(work < workmin)
	    {
	      scale = fac;
	      workmin = work;
	      iteropt = iter + 1;
	    }
	}

      dtnext = dttry / scale;
      if(iteropt >= lastiter && iteropt != nsd->maxiter - 1 && !reduce)
	{
	  fac = scale / alf[iteropt * nsd->maxstep + iteropt - 1];
	  fac = fac > 0.1 ? fac : 0.1;
	  if(a[iteropt + 1] * fac <= workmin)
	    {
	      dtnext = dttry / fac;
	      iteropt += 1;
	    }
	}

      time += dttry;

      if(time + dtnext > maxtime)
	{
	  first = 1;
	  iteropt = nsd->maxiter - 1;
	  dttry = maxtime - time;
	}
      else
	{
	  dttry = dtnext;
	}
      step++;

      if(traj)
	{
	  traj->time = time;
	}
    }

#ifdef NUCLEARNET_OUTPUT_TIMEEVOLUTION
  printf("t(%03d): %11.5e, dt (%02d): %11.5e, x:", step, time, m, dttry);
  for(i = 0; i < nsd->nelements; i++)
    {
      printf(" %8.1e", y[i] * nsd->aion[i]);
    }
  if(nsd->nelements < nsd->matrixsize)
    {
      printf(" %8.1e", y[nsd->nelements]);
    }
  printf("\n");
  fwrite(&time, sizeof(double), 1, fp);
  fwrite(y, sizeof(double), nsd->matrixsize, fp);
  fclose(fp);
#endif
}

#if defined(NETWORK_SUPERLU)
static void network_solver_calc_dy(double temp, double rho, double dt, int steps, double *y, double *yn,
				   const double *first_rhs, const jacob_t * old_jacob,
				   const struct network_data *nd, struct network_workspace *nw,
				   struct network_solver_trajectory *traj)
{
  double *dy, *rhs;
  jacob_t *jacob;
  network_var *deriv;
  int nMatrix;
  int i, j;
  double h;
  SuperMatrix A, B, L, U;
  superlu_options_t options;
  SuperLUStat_t stat;
  int *perm_c, *perm_r, info;
  struct network_solver_data *nsd = nw->nsd;

  /* take variables from global struct */
  dy = nsd->dy;
  rhs = nsd->rhs;
  deriv = nsd->rhs_deriv;
  jacob = &nsd->mod_jacob;
  nMatrix = nsd->matrixsize;

  h = dt / steps;		/* subtimestep */
  memcpy(rhs, first_rhs, sizeof(*first_rhs) * nMatrix);	/* because rhs is used for in-place output */
  prepare_matrix(jacob, old_jacob, nMatrix, h);

  dCreate_CompRow_Matrix(&A, nMatrix, nMatrix, jacob->n_elements, jacob->values, jacob->columns,
			 jacob->rowstart, SLU_NR, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&B, nMatrix, 1, rhs, nMatrix, SLU_DN, SLU_D, SLU_GE);

  perm_r = intMalloc(nMatrix);
  perm_c = intMalloc(nMatrix);

  set_default_options(&options);
  options.PrintStat = NO;
  StatInit(&stat);

  /* first step */
  for(i = 0; i < nMatrix; i++)
    {
      rhs[i] = h * rhs[i];
    }
  /* do decomposition only once, only here */
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
  options.Fact = FACTORED;	/* Indicate the factored form of A is supplied. */
  for(i = 0; i < nMatrix; i++)
    {
      dy[i] = rhs[i];
      yn[i] = y[i] + dy[i];
    }

  /* middle step */
  for(j = 1; j < steps; j++)
    {
      sanitize(yn, nd);
      if(traj)
	network_solver_interpolate_traj_rhotemp(traj, traj->time + j * h, yn, &rho, &temp);
      network_getrhs(rho, temp, yn, 0, nd, nw, rhs, deriv);
      for(i = 0; i < nMatrix; i++)
	{
	  rhs[i] = h * rhs[i] - dy[i];
	}
      dgstrs(TRANS, &L, &U, perm_c, perm_r, &B, &stat, &info);
      for(i = 0; i < nMatrix; i++)
	{
	  dy[i] = dy[i] + 2.0 * rhs[i];
	  yn[i] = yn[i] + dy[i];
	}
    }

  /* last step */
  sanitize(yn, nd);
  if(traj)
    network_solver_interpolate_traj_rhotemp(traj, traj->time + steps * h, yn, &rho, &temp);
  network_getrhs(rho, temp, yn, 0, nd, nw, rhs, deriv);
  for(i = 0; i < nMatrix; i++)
    {
      rhs[i] = h * rhs[i] - dy[i];
    }
  dgstrs(TRANS, &L, &U, perm_c, perm_r, &B, &stat, &info);
  for(i = 0; i < nMatrix; i++)
    {
      yn[i] = yn[i] + rhs[i];
    }

  StatFree(&stat);

  SUPERLU_FREE(perm_c);
  SUPERLU_FREE(perm_r);

  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);

  Destroy_SuperMatrix_Store(&A);
  Destroy_SuperMatrix_Store(&B);
}
#elif defined(NETWORK_PARDISO)
static void network_solver_calc_dy(double temp, double rho, double dt, int steps, double *y, double *yn,
				   const double *first_rhs, const jacob_t * old_jacob,
				   const struct network_data *nd, struct network_workspace *nw,
				   struct network_solver_trajectory *traj)
{
  double *dy, *rhs;
  network_var *deriv;
  jacob_t *jacob;
  double *values;
  int *rowstart, *columns;
  int matrixsize, matrixsize2;
  struct network_solver_data *nsd = nw->nsd;
  int i, j;
  double h;
  int mtype = 11;		/* real unsymmetric matrix */
  void *pt[64] = { 0 };		/* internal solver memory pointer (must be initialized to zero) */
  int iparm[64] = { 0 };
  int maxfct = 1;		/* maximum number of numerical factorization */
  int mnum = 1;			/* which factorization to use */
  int msglvl = 0;		/* print no statistical information */
  int nrhs = 1;			/* number of right hand sides to be solved for */
  int error = 0;		/* error flag */
  int idum;			/* dummy */
  double ddum;			/* dummy */
  double *ddumarray;		/*dummy */
  int phase;
  const int refinement = 2;

  /* Pardiso control parameters */
  iparm[0] = 1;			/* No solver default */
  iparm[1] = 0;			/* Minimum degree fill-in reordering */
  /*iparm[1] = 2; *//* Fill-in reordering from METIS */
  iparm[2] = mkl_get_max_threads();	/* number of processors */
  iparm[3] = 0;			/* No iterative-direct algorithm */
  iparm[4] = 0;			/* No user fill-in reducing permutation */
  iparm[5] = 1;			/* Write solution into rhs */
  iparm[6] = 0;			/* Not in use */
  iparm[7] = refinement;	/* Max numbers of iterative refinement steps */
  iparm[8] = 0;			/* Not in use */
  iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;		/* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;		/* Not in use */
  iparm[12] = 0;		/* Not in use */
  iparm[13] = 0;		/* Output: Number of perturbed pivots */
  iparm[14] = 0;		/* Not in use */
  iparm[15] = 0;		/* Not in use */
  iparm[16] = 0;		/* Not in use */
  iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1;		/* Output: Mflops for LU factorization */
  iparm[19] = 0;		/* Output: Numbers of CG Iterations */
#ifdef DEBUG
  iparm[26] = 1;		/* Check supplied matrix indices */
#endif /* DEBUG */

  /* take variables from global struct */
  dy = nsd->dy;
  rhs = nsd->rhs;
  deriv = nsd->rhs_deriv;
  jacob = &nsd->mod_jacob;
  values = jacob->values;
  columns = jacob->columns;
  rowstart = jacob->rowstart;
  matrixsize = nsd->matrixsize;
  matrixsize2 = nsd->matrixsize2;

  h = dt / steps;		/* subtimestep */
  memcpy(rhs, first_rhs, sizeof(*first_rhs) * matrixsize);	/* because rhs is used for in-place output */
  prepare_matrix(jacob, old_jacob, matrixsize, h);

  /* We don't use ddumarray for the result but PARDISO requests its existence. */
  if(!(ddumarray = malloc(matrixsize * sizeof(double))))
    {
      perror("error allocating temporary space");
      endrun(1002);
    }

  /* do decomposition only once */
  phase = 12;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &matrixsize, values, rowstart, columns, &idum, &nrhs, iparm,
	  &msglvl, &ddum, &ddum, &error);
  if(error != 0)
    {
      myprintf("error during LU decomposition\n");
      endrun(1003);
    }

  /* first step */
  for(i = 0; i < matrixsize; i++)
    rhs[i] = h * rhs[i];
  iparm[7] = refinement;	/* Max numbers of iterative refinement steps. */
  phase = 33;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &matrixsize, values, rowstart, columns, &idum, &nrhs, iparm,
	  &msglvl, rhs, ddumarray, &error);
  if(error != 0)
    {
      myprintf("error during first step\n");
      endrun(1004);
    }

  for(i = 0; i < matrixsize; i++)
    {
      dy[i] = rhs[i];
      yn[i] = y[i] + dy[i];
    }

  /* middle step */
  for(j = 1; j < steps; j++)
    {
      sanitize(yn, nd);
      if(traj)
	network_solver_interpolate_traj_rhotemp(traj, traj->time + j * h, yn, &rho, &temp);
      network_getrhs(rho, temp, yn, 0, nd, nw, rhs, deriv);
      for(i = 0; i < matrixsize; i++)
	rhs[i] = h * rhs[i] - dy[i];
      iparm[7] = refinement;	/* Max numbers of iterative refinement steps. */
      PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &matrixsize, values, rowstart, columns, &idum, &nrhs, iparm,
	      &msglvl, rhs, ddumarray, &error);
      if(error != 0)
	{
	  myprintf("error during middle step\n");
	  endrun(1005);
	}
      for(i = 0; i < matrixsize; i++)
	{
	  dy[i] = dy[i] + 2.0 * rhs[i];
	  yn[i] = yn[i] + dy[i];
	}
    }

  /* last step */
  sanitize(yn, nd);
  if(traj)
    network_solver_interpolate_traj_rhotemp(traj, traj->time + steps * h, yn, &rho, &temp);
  network_getrhs(rho, temp, yn, 0, nd, nw, rhs, deriv);
  for(i = 0; i < matrixsize; i++)
    rhs[i] = h * rhs[i] - dy[i];
  iparm[7] = refinement;	/* Max numbers of iterative refinement steps. */
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &matrixsize, values, rowstart, columns, &idum, &nrhs, iparm,
	  &msglvl, rhs, ddumarray, &error);
  if(error != 0)
    {
      myprintf("error during last step\n");
      endrun(1006);
    }
  for(i = 0; i < matrixsize; i++)
    {
      yn[i] = yn[i] + rhs[i];
    }

  /* release memory */
  phase = -1;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &matrixsize, values, rowstart, columns, &idum, &nrhs, iparm,
	  &msglvl, &ddum, &ddum, &error);
  free(ddumarray);
}

#else
static void network_solver_calc_dy(double temp, double rho, double dt, int steps, double *y, double *yn,
				   const double *first_rhs, const jacob_t * old_jacob,
				   const struct network_data *nd, struct network_workspace *nw,
				   struct network_solver_trajectory *traj)
{
  double *dy;
  network_var *deriv;
  double *rhs;
  jacob_t *jacob;
  int matrixsize, matrixsize2;
  int i, j, s;
  double h;
  gsl_matrix_view A;
  gsl_vector_view b;
  gsl_vector *x;
  gsl_permutation *p;
  struct network_solver_data *nsd = nw->nsd;

  /* take variables from global struct */
  dy = nsd->dy;
  rhs = nsd->rhs;
  deriv = nsd->rhs_deriv;
  jacob = &nsd->mod_jacob;
  matrixsize = nsd->matrixsize;
  matrixsize2 = nsd->matrixsize2;

  h = dt / steps;		/* subtimestep */
  memcpy(rhs, first_rhs, sizeof(*first_rhs) * matrixsize);
  prepare_matrix(jacob, old_jacob, matrixsize, h);

  A = gsl_matrix_view_array(jacob->matrix, matrixsize, matrixsize);
  x = gsl_vector_alloc(matrixsize);
  p = gsl_permutation_alloc(matrixsize);

  /* do decomposition only once */
  gsl_linalg_LU_decomp(&A.matrix, p, &s);

  /* first step */
  for(i = 0; i < matrixsize; i++)
    {
      rhs[i] = h * rhs[i];
    }
  b = gsl_vector_view_array(rhs, matrixsize);
  gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);
  for(i = 0; i < matrixsize; i++)
    {
      dy[i] = gsl_vector_get(x, i);
      yn[i] = y[i] + dy[i];
    }

  /* middle step */
  for(j = 1; j < steps; j++)
    {
      sanitize(yn, nd);
      if(traj)
	network_solver_interpolate_traj_rhotemp(traj, traj->time + j * h, yn, &rho, &temp);
      network_getrhs(rho, temp, yn, 0, nd, nw, rhs, deriv);
      for(i = 0; i < matrixsize; i++)
	{
	  rhs[i] = h * rhs[i] - dy[i];
	}
      b = gsl_vector_view_array(rhs, matrixsize);
      gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);
      for(i = 0; i < matrixsize; i++)
	{
	  dy[i] = dy[i] + 2.0 * gsl_vector_get(x, i);
	  yn[i] = yn[i] + dy[i];
	}
    }

  /* last step */
  sanitize(yn, nd);
  if(traj)
    network_solver_interpolate_traj_rhotemp(traj, traj->time + steps * h, yn, &rho, &temp);
  network_getrhs(rho, temp, yn, 0, nd, nw, rhs, deriv);
  for(i = 0; i < matrixsize; i++)
    {
      rhs[i] = h * rhs[i] - dy[i];
    }
  b = gsl_vector_view_array(rhs, matrixsize);
  gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);
  for(i = 0; i < matrixsize; i++)
    {
      yn[i] = yn[i] + gsl_vector_get(x, i);
    }

  gsl_vector_free(x);
  gsl_permutation_free(p);
}
#endif

void network_solver_extrapolate(int iter, int matrixsize, double *x, double *y, double *dy, double *qcol)
{
  /* This implements polynomial extrapolation as in Deuflhard 1983 eq. (1.3) */
  size_t vsize = matrixsize * sizeof(double);

  if(iter == 0)
    {
      memcpy(dy, y, vsize);
    }
  else
    {
      int k;

      for(k = 0; k < iter; k++)
	{
	  int i;
	  double fac = 1.0 / (x[iter - k - 1] / x[iter] - 1.0);

	  for(i = 0; i < matrixsize; i++)
	    {
	      /* change in y to reach the next column in the extrapolation tableau */
	      dy[i] = (y[i] - qcol[k * matrixsize + i]) * fac;

	      /* save the old approximation for the next step */
	      qcol[k * matrixsize + i] = y[i];

	      y[i] += dy[i];
	    }
	}
    }

  memcpy(&qcol[iter * matrixsize], y, vsize);
}

static void print_matrix(const int n, const double *matrix, const int *rowstart, const int *columns,
			 const double *rhs, const int offset)
{
  int i, j;
  const int *col;
  const double *a;
  FILE *out;

  if(!(out = fopen("matrix.m", "w")))
    {
      perror("cannot open matrix output file");
      return;
    }
  fprintf(out, "A = [ ");
  col = columns + *rowstart - offset;
  a = matrix + *rowstart - offset;
  for(i = 0; i < n; i++)
    {
      assert(col == columns + rowstart[i] - offset);
      assert(a == matrix + rowstart[i] - offset);
      for(j = 0; j < n; j++)
	{
	  if(col == columns + rowstart[i + 1] - offset)
	    break;
	  if(*col - offset == j)
	    {
	      fprintf(out, "%g ", *a);
	      col++;
	      a++;
	    }
	  else
	    fprintf(out, "%g ", 0.);
	}
      for(; j < n; j++) {fprintf(out, "%g ", 0.);}
      fprintf(out, "; ");
    }
  assert(col == columns + rowstart[i] - offset);
  assert(a == matrix + rowstart[i] - offset);
  fprintf(out, "]\n");
  fprintf(out, "b = [ ");
  for(i = 0; i < n; i++) {fprintf(out, "%g; ", rhs[i]);}
  fprintf(out, "]\n");
  fclose(out);
}

static void normalize(double *y, const struct network_solver_data *nsd)
{
  int i;
  double sum;

  sum = 0;
  for(i = 0; i < nsd->nelements; i++)
    {
      sum += y[i] * nsd->aion[i];
    }
  for(i = 0; i < nsd->nelements; i++)
    y[i] /= sum;
}

static void sanitize(double *y, const struct network_data *nd)
{
  double *yy, *yend = &y[nd->nuc_count];
  for(yy = y; yy != yend; yy++)
    *yy = fmin(1.0, fmax(1e-30, *yy));
#ifdef NETWORK_SEP_YZ
  y[nd->iYz] = fmax(1e-30, fmin(y[nd->iYz], nd->nucdata[nd->nuc_count - 1].nz));
#endif
#if NETWORK_VARIABLE
  y[nd->iTemp] = fmax(1e7, fmin(y[nd->iTemp], NETWORK_T_MAX));
#endif
#if NETWORK_VARIABLE == NETWORK_VAR_RHO_TEMP
  y[nd->iRho] = fmax(1e-6, fmin(y[nd->iRho], 1e11));
#endif
}

static void prepare_matrix(jacob_t * mod_jacob, const jacob_t * old_jacob, size_t n_matrix, double h)
{
  /* this computes (1 - h * Jacob) */
  int i;
#if NETWORK_SPARSE
#if defined(NETWORK_SUPERLU)
  const int offset = 0;		/* C-style indices */
#elif defined(NETWORK_PARDISO)
  const int offset = 1;		/* Fortran-style indices */
#else
#error "unknown sparse matrix library"
#endif

  mod_jacob->n_elements = old_jacob->n_elements;
  memcpy(mod_jacob->rowstart, old_jacob->rowstart, sizeof(*old_jacob->rowstart) * (n_matrix + 1));
  memcpy(mod_jacob->columns, old_jacob->columns, sizeof(*old_jacob->columns) * old_jacob->n_elements);

  for(i = 0; i < n_matrix; i++)
    {
      int j;

      for(j = old_jacob->rowstart[i] - offset; j < old_jacob->rowstart[i + 1] - offset; j++)
	{
	  mod_jacob->values[j] = -h * old_jacob->values[j];
	  if(mod_jacob->columns[j] - offset == i)
	    mod_jacob->values[j] += 1.0;
	}
    }
#else
  size_t n_matrix2 = n_matrix * n_matrix;
  const double *old_matrix = old_jacob->matrix;
  double *mod_matrix = mod_jacob->matrix;
  for(i = 0; i < n_matrix2; i++)
    {
      mod_matrix[i] = (-h) * old_matrix[i];
      if(i % (n_matrix + 1) == 0)
	mod_matrix[i] += 1.0;
    }
#endif /* NETWORK_SPARSE */
}

#endif /* NUCLEAR_NETWORK */
