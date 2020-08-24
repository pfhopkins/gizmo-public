#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#define GSLWORKSIZE 100000

/*! \file sidm_routines.c
 *  \brief Fuctions and routines needed for the calculations of dark matter self interactions
 *
 *  This file contains the functions and routines necesary for the computation of
 *  the self-interaction probabilities and the velocity kicks due to the interactios.
 *  Originally written by Miguel Rocha, rocham@uci.edu. Oct 2010. Updated on 2014 & re-written by PFH March 2018
 */

/*! This function calculates the interaction probability between two particles.
 *  It checks if comoving integration is on and does the necesary change of
 *  variables and units.
 */

#ifdef DM_SIDM


double prob_of_interaction(double mass, double r, double h_si, double dV[3], double dt)
{
    double dVmag = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]) / All.cf_atime; // velocity in physical
    double rho_eff = mass / (h_si*h_si*h_si) * All.cf_a3inv; // density in physical
    double cx_eff = All.DM_InteractionCrossSection * g_geo(r/h_si); // effective cross section (physical) scaled to cgs
    double units = UNIT_SURFDEN_IN_CGS; // needed to convert everything to cgs
    if(All.DM_InteractionVelocityScale>0) {double x=dVmag/All.DM_InteractionVelocityScale; cx_eff/=1+x*x*x*x;} // take velocity dependence
    return rho_eff * cx_eff * dVmag * dt * units; // dimensionless probability
}

/*! This routine sets the kicks for each particle after it has been decided that they will
 *  interact. It uses an algorithm tha conserves energy and momentum but picks a random direction so it does not conserves angular momentum. */
#if !defined(GRAIN_COLLISIONS) /* if using the 'grain collisions' module, these functions will be defined elsewhere [in the grains subroutines] */
void calculate_interact_kick(double dV[3], double kick[3], double m)
{
    double dVmag = (1-All.DM_DissipationFactor)*sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]);
    if(dVmag<0) {dVmag=0;}
    if(All.DM_KickPerCollision>0) {double v0=All.DM_KickPerCollision; dVmag=sqrt(dVmag*dVmag+v0*v0);}
    double cos_theta = 2.0*gsl_rng_uniform(random_generator)-1.0, sin_theta = sqrt(1.-cos_theta*cos_theta), phi = gsl_rng_uniform(random_generator)*2.0*M_PI;
    kick[0] = 0.5*(dV[0] + dVmag*sin_theta*cos(phi));
    kick[1] = 0.5*(dV[1] + dVmag*sin_theta*sin(phi));
    kick[2] = 0.5*(dV[2] + dVmag*cos_theta);
}
#endif


/*! This function returns the value of the geometrical factor needed for the calculation of the interaction probability. */
double g_geo(double r)
{
    double f, u; int i; u = r / 2.0 * GEOFACTOR_TABLE_LENGTH; i = (int) u;
    if(i >= GEOFACTOR_TABLE_LENGTH) {i = GEOFACTOR_TABLE_LENGTH - 1;}
    if(i <= 1) {f = 0.992318  + (GeoFactorTable[0] - 0.992318)*u;} else {f = GeoFactorTable[i - 1] + (GeoFactorTable[i] - GeoFactorTable[i - 1]) * (u - i);}
    return f;
}

/*! This routine initializes the table that will be used to get the geometrical factor
 *  as a function of the two particle separations. It populates a table with the results of the numerical integration */
void init_geofactor_table(void)
{
    int i; double result, abserr,r;
    gsl_function F; gsl_integration_workspace *workspace; workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);
    for(i = 0; i < GEOFACTOR_TABLE_LENGTH; i++)
    {
        r =  2.0/GEOFACTOR_TABLE_LENGTH * (i + 1);
        F.function = &geofactor_integ;
        F.params = &r;
        gsl_integration_qag(&F, 0.0, 1.0, 0, 1.0e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41,workspace, &result, &abserr);
        GeoFactorTable[i] = 2*M_PI*result;
    }
    gsl_integration_workspace_free(workspace);
}

/*! This function returns the integrand of the numerical integration done on init_geofactor_table(). */
double geofactor_integ(double x, void * params)
{
    double result, abserr, r, newparams[2];
    r = *(double *) params; newparams[0] = r; newparams[1] = x;
    gsl_function F; gsl_integration_workspace *workspace; workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);
    F.function = &geofactor_angle_integ; F.params = newparams;
    
    gsl_integration_qag(&F, -1.0, 1.0, 0, 1.0e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41,workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    /*! This function returns the value W(x). The values of the density kernel as a funtion of x=r/h */
    double wk=0; if(x<1) kernel_main(x, 1, 1, &wk, &wk, -1);
    return x*x*wk*result;
}

/*! This function returns the integrand of the angular part of the integral done on init_geofactor_table(). */
double geofactor_angle_integ(double u, void * params)
{
    double x,r,f;
    r = *(double *) params;
    x = *(double *) (params + sizeof(double));
    f = sqrt(x*x + r*r + 2*x*r*u);
    double wk=0; if(f<1) kernel_main(f, 1, 1, &wk, &wk, -1); /*! This function returns the value W(x). The values of the density kernel as a funtion of x=r/h */
    return wk;
}

/*! This function simply initializes some variables to prevent memory errors */
void init_self_interactions() {int i; for(i = 0; i < NumPart; i++) {P[i].dtime_sidm = 0; P[i].NInteractions = 0;}}

#endif
