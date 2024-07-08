#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "./cooling.h"

/*
 * This file contains the routines for optically-thin cooling (generally aimed towards simulations of the ISM,
 *   galaxy formation, and cosmology). A wide range of heating/cooling processes are included, including
 *   free-free, metal-line, Compton, collisional, photo-ionization and recombination, and more. Some of these
 *   are controlled by individual modules that need to be enabled or disabled explicitly.
 *
 * This file was originally part of the GADGET3 code developed by Volker Springel. The code has been modified heavily by
 *   Phil Hopkins (phopkins@caltech.edu) for GIZMO; essentially everything has been re-written at this point */


#ifdef COOLING

/* these are variables of the cooling tables. they are static but this shouldnt be a problem for shared-memory structure because
    they are only defined once in a global operation, then locked for particle-by-particle operations */
/* requires the cooling table TREECOOL, which is included in the GIZMO source in the cooling directory */
#define NCOOLTAB  2000 /* defines size of cooling table */

#if !defined(CHIMES)
static double Tmin = -1.0, Tmax = 9.0, deltaT; /* minimum/maximum temp, in log10(T/K) and temperature gridding: will be appropriately set in make_cooling_tables subroutine below */
static double *BetaH0, *BetaHep, *Betaff, *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp, *GammaeH0, *GammaeHe0, *GammaeHep; // UV background parameters
#ifdef COOL_METAL_LINES_BY_SPECIES
/* if this is enabled, the cooling table files should be in a folder named 'spcool_tables' in the run directory.
 cooling tables can be downloaded at: http://www.tapir.caltech.edu/~phopkins/public/spcool_tables.tgz or on the Bitbucket site (downloads section) */
static float *SpCoolTable0, *SpCoolTable1;
#endif
/* these are constants of the UV background at a given redshift: they are interpolated from TREECOOL but then not modified particle-by-particle */
static double J_UV = 0, gJH0 = 0, gJHep = 0, gJHe0 = 0, epsH0 = 0, epsHep = 0, epsHe0 = 0;
#endif

#if defined(CHIMES)
int ChimesEqmMode, ChimesUVBMode, ChimesInitIonState, N_chimes_full_output_freq, Chimes_incl_full_output = 1;
double chimes_rad_field_norm_factor, shielding_length_factor, cr_rate;
char ChimesDataPath[256], ChimesEqAbundanceTable[196], ChimesPhotoIonTable[196];
struct gasVariables *ChimesGasVars;
struct globalVariables ChimesGlobalVars;
#ifdef CHIMES_METAL_DEPLETION
struct Chimes_depletion_data_structure *ChimesDepletionData;
#endif
#endif



/* this is the 'parent' loop to do the cell cooling+chemistry. this is now openmp-parallelized, since the semi-implicit iteration can be a non-negligible cost */
void cooling_parent_routine(void)
{
    PRINT_STATUS("Cooling and Chemistry update");
    /* Determine indices of active particles. */
    int N_active=0, i, j, *active_indices; active_indices = (int *) malloc(N_gas * sizeof(int));
    for (i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 0) {continue;}
        if(P[i].Mass <= 0) {continue;}
#ifdef GALSF_EFFECTIVE_EQS
        if((SphP[i].Density*All.cf_a3inv > All.PhysDensThresh) && ((All.ComovingIntegrationOn==0) || (SphP[i].Density>=All.OverDensThresh))) {continue;} /* no cooling for effective-eos star-forming particles */
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        if(SphP[i].DelayTimeCoolingSNe > 0) {continue;} /* no cooling for particles marked in delayed cooling */
#endif
        active_indices[N_active] = i;
        N_active++;
	}

#ifdef _OPENMP
#pragma omp parallel private(i, j)
#endif
    { /* open parallel block */
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for(j=0;j<N_active;j++)
    {
        i=active_indices[j]; /* actual particle index */
        do_the_cooling_for_particle(i); /* do the actual cooling */
    }
    } /* close parallel block */
    free(active_indices); /* free memory */

#ifdef CHIMES /* CHIMES records some extra timing information here owing to large possible imbalances */
  CPU_Step[CPU_COOLINGSFR] += measure_time(); MPI_Barrier(MPI_COMM_WORLD);
  CPU_Step[CPU_COOLSFRIMBAL] += measure_time(); PRINT_STATUS("CHIMES chemistry and cooling finished");
#endif
}




/* subroutine which actually sends the particle data to the cooling routine and updates the entropies */
void do_the_cooling_for_particle(int i)
{
    double unew, dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i), ne_in, ne_out;

    if((dtime>0)&&(P[i].Mass>0)&&(P[i].Type==0))  // upon start-up, need to protect against dt==0 //
    {
        double uold = DMAX(All.MinEgySpec, SphP[i].InternalEnergy); int k; k=0; ne_in=0; ne_out=0;
        
#if defined(RADTRANSFER)
        for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[i].Lambda_RadiativeCooling_toRHDBins[k]=0;} // zero these out before cooling subroutine
#endif
        
#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction(i, 0.5*dtime*UNIT_TIME_IN_CGS); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif

#ifndef COOLING_OPERATOR_SPLIT
        double DtInternalEnergyEffCGS = SphP[i].DtInternalEnergy;
        /* do some prep operations on the hydro-step determined heating/cooling rates before passing to the cooling subroutine */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
        double grav_acc;
        for(k = 0; k < 3; k++)
        {
            grav_acc = All.cf_a2inv * P[i].GravAccel[k];
#ifdef PMGRID
            grav_acc += All.cf_a2inv * P[i].GravPM[k];
#endif
            DtInternalEnergyEffCGS -= SphP[i].GravWorkTerm[k] * All.cf_atime * grav_acc;
        }
#endif
        /* limit the magnitude of the hydro dtinternalenergy */
        if(DtInternalEnergyEffCGS < 0) {
            double qfac = DMIN(0,DMAX(DMAX(-0.9, exp(DtInternalEnergyEffCGS*dtime/SphP[i].InternalEnergy)-1.), All.MinEgySpec/SphP[i].InternalEnergy-1.)); // equivalent to saying this wouldn't lower internal energy to below 10% in one timestep
            DtInternalEnergyEffCGS = DMAX(DtInternalEnergyEffCGS , qfac*SphP[i].InternalEnergy/dtime );
            double u_gamma_minus_1 = (GAMMA(i)-1.) * SphP[i].InternalEnergy, rho = SphP[i].Density*All.cf_a3inv, pressure_thermalonly = u_gamma_minus_1 * rho;
            double vA = Get_Gas_Alfven_speed_i(i), pressure_total = 0.5*vA*vA*rho + SphP[i].Pressure*All.cf_a3inv;
            if(pressure_thermalonly < 0.05*pressure_total) {
                double DtInternalEnergyPdV = - u_gamma_minus_1 * (P[i].Particle_DivVel*All.cf_a2inv); /* change from expansion in PdV term */
                DtInternalEnergyEffCGS = DMAX(DtInternalEnergyEffCGS , DMIN(DtInternalEnergyPdV, 0)); /* limit to PdV expansion change in limit where the thermal energy is small compared to the total */
            }
        }
        DtInternalEnergyEffCGS = DMIN(DtInternalEnergyEffCGS ,  1.e4*SphP[i].InternalEnergy/dtime ); // equivalent to saying we cant massively enhance internal energy in a single timestep from the hydro work terms: should be big, since just numerical [shocks are real!]
        /* and convert to cgs before use in the cooling sub-routine */
        DtInternalEnergyEffCGS *= (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS_CGS/HYDROGEN_MASSFRAC);
        if(SphP[i].CoolingIsOperatorSplitThisTimestep==0) {SphP[i].DtInternalEnergy = DtInternalEnergyEffCGS;} // if unsplit, send this converted variable to cooling below
#endif


#if !defined(RT_COOLING_PHOTOHEATING_OLDFORMAT) /* standard behavior, call the actual cooling subroutine */
#if !defined(CHIMES)
        ne_in = SphP[i].Ne; ne_out = ne_in; /* this variable is not defined if chimes is on */
#endif
        unew = DoCooling(uold, SphP[i].Density * All.cf_a3inv, dtime, ne_in, &ne_out, i);
#if !defined(CHIMES)
        SphP[i].Ne = ne_out; /* update the free electron variable */
#endif
#else
        unew = uold + dtime * (rt_DoHeating(i, dtime) + rt_DoCooling(i, dtime));  /* call this special subroutine for cooling for the rt 'oldformat' [limited applications] */
#endif




#if defined(BH_THERMALFEEDBACK)
        if(SphP[i].Injected_BH_Energy) {unew += SphP[i].Injected_BH_Energy / P[i].Mass; SphP[i].Injected_BH_Energy = 0;}
#endif


#if defined(COSMIC_RAY_FLUID) && !defined(CRFLUID_ALT_DISABLE_LOSSES)
        CR_cooling_and_losses(i, SphP[i].Ne, SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, dtime*UNIT_TIME_IN_CGS );
#endif
        

#if defined(RADTRANSFER) /* account for cooling radiation which should, according to our modules, come out in certain bands */
        double nHcgs = HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * SphP[i].Density * All.cf_a3inv / PROTONMASS_CGS;    /* hydrogen number dens in cgs units */
        double ratefact = (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE) * nHcgs * nHcgs / (SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS) * (dtime*UNIT_TIME_IN_CGS) / (UNIT_SPECEGY_IN_CGS) * P[i].Mass; /* need to account for RSOL factors in emission/absorption rates */
        double de_u = (unew - SphP[i].InternalEnergy) * P[i].Mass; /* change in the total internal energy of the gas cell [integrating over everything] */
        double de_rad_tot_final = 0, de_rad_tot = 0; for(k=0;k<N_RT_FREQ_BINS;k++) {de_rad_tot += SphP[i].Lambda_RadiativeCooling_toRHDBins[k] * ratefact;} /* energy gained by gas needs to be subtracted from radiation. positive lambda means gas cooling (gas energy loss, so radiation energy gain, so positive here) */
        /* Removing the de_u * de_rad_tot > 0 limiter because this is not ruled out physically, e.g. if gas is being heated by PdV work while radiating away energy. 
        Can only do the limiter if the cooling function accounts for radiative processes only. */
	    //if(de_u * de_rad_tot > 0) {de_rad_tot = 0;} /* if radiation gains but gas net loses (or vice versa), could occur across different bands but don't do the routine below */
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            if((fabs(SphP[i].Lambda_RadiativeCooling_toRHDBins[k]) > MIN_REAL_NUMBER) && (fabs(de_rad_tot) > MIN_REAL_NUMBER))
            {
                double de_rad = SphP[i].Lambda_RadiativeCooling_toRHDBins[k] * ratefact; /* energy gained by gas needs to be subtracted from radiation. positive lambda means gas cooling (gas energy loss, so radiation energy gain, so positive here) */
		        //if(de_u * de_rad > 0) {de_rad = 0;} /* if radiation gains but gas net loses (or vice versa), could occur across different bands but don't do the routine below  - removed, see note above */
                if(fabs(de_rad) > MIN_REAL_NUMBER)
                {
                    double de_rad_min = DMIN(-0.99*SphP[i].Rad_E_gamma[k], 0); // don't let the radiation loss take all the radiation energy into negative, or more than the energy gained from cooling+heating
                    double de_rad_max = DMAX(100.*unew*P[i].Mass, 0); // don't let the radiation gain take more than some large factor times the current energy, or more than the energy lost from cooling+heating
                    de_rad = DMAX(DMIN(de_rad, de_rad_max), de_rad_min); // limit de_rad appropriately                    
                    if(fabs(de_rad) > MIN_REAL_NUMBER)
                    {
                        de_rad_tot_final += de_rad; // add to our running total                        
#ifdef RT_INFRARED  
                       if(k==RT_FREQ_BIN_INFRARED) {SphP[i].Radiation_Temperature = SphP[i].Radiation_Temperature_CoolingWeighted;}// need to also update the IR band temperature measure
#endif
                        double Rad_E_gamma_before = SphP[i].Rad_E_gamma[k]; // save for immediate use below
                        SphP[i].Rad_E_gamma[k] += de_rad; /* energy gained by gas is lost here (or vice versa if dust is acting as a net coolant) */
                        SphP[i].Rad_E_gamma_Pred[k] = SphP[i].Rad_E_gamma[k]; /* updated drifted */
#if defined(RT_EVOLVE_INTENSITIES)
                        int k_tmp; for(k_tmp=0;k_tmp<N_RT_INTENSITY_BINS;k_tmp++) {SphP[i].Rad_Intensity[k][k_tmp] += de_rad/RT_INTENSITY_BINS_DOMEGA; SphP[i].Rad_Intensity_Pred[k][k_tmp] += de_rad/RT_INTENSITY_BINS_DOMEGA;}
#endif
                        int kv; // add leading-order relativistic corrections here, accounting for gas motion in the addition/subtraction to the flux
#if defined(RT_EVOLVE_FLUX)
                        double corrfac = 0; if(Rad_E_gamma_before > 0 && SphP[i].Rad_E_gamma[k] > 0) {corrfac = SphP[i].Rad_E_gamma[k] / (MIN_REAL_NUMBER + Rad_E_gamma_before);}
                        for(kv=0;kv<3;kv++) {if(corrfac > 0) {SphP[i].Rad_Flux[k][kv] *= corrfac; SphP[i].Rad_Flux_Pred[k][kv] *= corrfac;} else {double fluxfac = RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS*SphP[i].VelPred[kv]/All.cf_atime * de_rad; SphP[i].Rad_Flux[k][kv] += fluxfac; SphP[i].Rad_Flux_Pred[k][kv] += fluxfac;}}
#endif
                        double momfac = 1. - de_rad / (P[i].Mass * C_LIGHT_CODE*C_LIGHT_CODE_REDUCED); // back-reaction on gas from emission [note peculiar units here, its b/c of how we fold in the existing value of v and tilde[u] in our derivation - one rsol factor in denominator needed]
                        for(kv=0;kv<3;kv++) {P[i].Vel[kv] *= momfac; SphP[i].VelPred[kv] *= momfac;}
                    }
                }
            }
        }
        
#endif // done with RHD-cooling block update
        

        /* InternalEnergy, InternalEnergyPred, Pressure, ne are now immediately updated; however, if COOLING_OPERATOR_SPLIT
         is set, then DtInternalEnergy carries information from the hydro loop which is only half-stepped here, so is -not- updated.
         if the flag is not set (default), then the full hydro-heating is accounted for in the cooling loop, so it should be re-zeroed here */
        SphP[i].InternalEnergy = unew;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        SphP[i].Pressure = get_pressure(i);
#ifndef COOLING_OPERATOR_SPLIT
        if(SphP[i].CoolingIsOperatorSplitThisTimestep==0) {SphP[i].DtInternalEnergy=0;} // if unsplit, zero the internal energy change here
#endif

#if defined(GALSF_ISMDUSTCHEM_MODEL)
        update_dust_acc_and_sput(i, dtime*UNIT_TIME_IN_MYR*0.001);
#endif

#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction(i, 0.5*dtime*UNIT_TIME_IN_CGS); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif
        

    } // closes if((dt>0)&&(P[i].Mass>0)&&(P[i].Type==0)) check
}


/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
double DoCooling(double u_old, double rho, double dt, double ne_guess, double *ne_eval, int target)
{
    double u, du; u=0; du=0;

#ifdef COOL_GRACKLE
#ifndef COOLING_OPERATOR_SPLIT
    /* because grackle uses a pre-defined set of libraries, we can't properly incorporate the hydro heating
     into the cooling subroutine. instead, we will use the approximate treatment below to split the step */
    if(SphP[target].CoolingIsOperatorSplitThisTimestep==0)
    {
        du = dt * SphP[target].DtInternalEnergy / ( (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS_CGS/HYDROGEN_MASSFRAC));
        u_old += 0.5*du;
        u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
        /* now we attempt to correct for what the solution would have been if we had included the remaining half-step heating
         term in the full implicit solution. The term "r" below represents the exact solution if the cooling function has
         the form d(u-u0)/dt ~ -a*(u-u0)  around some u0 which is close to the "ufinal" returned by the cooling routine,
         to which we then add the heating term from hydro and compute the solution over a full timestep */
        double r=u/u_old; if(r>1) {r=1/r;} if(fabs(r-1)>1.e-4) {r=(r-1)/log(r);} r=DMAX(0,DMIN(r,1));
        du *= 0.5*r; if(du<-0.5*u) {du=-0.5*u;} u+=du;
    }
#else
    /* with full operator splitting we just call grackle normally. note this is usually fine,
     but can lead to artificial noise at high densities and low temperatures, especially if something
     like artificial pressure (but not temperature) floors are used such that the temperature gets
     'contaminated' by the pressure terms */
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
#endif
    return DMAX(u,All.MinEgySpec);
#endif

#ifdef CHIMES
    chimes_update_gas_vars(target);

    /* Call CHIMES to evolve the chemistry and temperature over
     * the hydro timestep. */
    chimes_network(&(ChimesGasVars[target]), &ChimesGlobalVars);

    // Compute updated internal energy
    u = (double) ChimesGasVars[target].temperature * BOLTZMANN_CGS / ((GAMMA(target)-1) * PROTONMASS_CGS * calculate_mean_molecular_weight(&(ChimesGasVars[target]), &ChimesGlobalVars));
    u /= UNIT_SPECEGY_IN_CGS;  // code units

#ifdef CHIMES_TURB_DIFF_IONS 
    chimes_update_turbulent_abundances(target, 1); 
#endif 

    return DMAX(u, All.MinEgySpec);

#else // CHIMES

    int iter=0, iter_upper=0, iter_lower=0, iter_condition = 0, maxiter_uplo=10*MAXITER; double LambdaNet, ratefact, u_upper, u_lower;
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    double u_min = All.MinEgySpec * UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho;
    u = u_upper = u_lower =  u_old; /* initialize values */
    #define ROOTFIND_FUNCTION(du) du - ratefact * CoolingRateFromU(u_old+du, rho, ne_guess, ne_eval, target) * dt // control the *relative* error on the *change* in u
    double du_net = ROOTFIND_FUNCTION(u - u_old), du_net_upper = du_net, du_net_lower = du_net;

    /* bracketing */
    double u_step_fac = 1.1;
    int bracket_iter = 0, skip_rootfind = 0;
    while(du_net_upper * du_net_lower > 0 && bracket_iter<MAXITER){
        if((u_lower <= All.MinEgySpec) && (du_net_lower > 0) && (du_net_upper > 0)){skip_rootfind = 1; break;} // will never find the root because bouncing off the lower limit
        u_upper *= u_step_fac; 
        du_net_upper = ROOTFIND_FUNCTION(u_upper - u_old);
	if(du_net*du_net_upper < 0){u_lower = u; du_net_lower = du_net; break;} // let u_upper and u_old be the brackets
	u_lower = DMAX(u_lower/u_step_fac,All.MinEgySpec); // bound u_lower because we don't trust cooling function below this
        du_net_lower = ROOTFIND_FUNCTION(u_lower - u_old); 
	if(du_net*du_net_lower < 0){u_upper = u; du_net_upper = du_net; break;} // let u_lower an u_old be the brackets
        u_step_fac *= 1.1;
        bracket_iter++;
    }
    if(fabs(du_net_upper) < MIN_REAL_NUMBER) {u=u_upper;}
    else if(fabs(du_net_lower) < MIN_REAL_NUMBER) {u=u_lower;}
    else {
        if(!skip_rootfind){ // assuming we're not bouncing off the min temp
            if((du_net_upper * du_net_lower >= 0) || isnan(du_net_lower) || isnan(du_net_upper)) {PRINT_WARNING("Could not bracket cooling solution. ID=%lld u_min=%g u=%g u_lower=%g u_upper=%g f_lower=%g f_upper=%g\n", (long long)P[target].ID, u_min, u, u_lower,u_upper, du_net_lower, du_net_upper); endrun(10);}
            
            /* core iteration to convergence */
            double ROOTFIND_X_a = u_upper-u_old, ROOTFIND_X_b = u_lower-u_old, ROOTFUNC_a = du_net_upper, ROOTFUNC_b = du_net_lower, ROOTFIND_REL_X_tol = 1e-2, ROOTFIND_ABS_X_tol = 1e-15 * u_old;
            #include "../system/bracketed_rootfind.h"
            u = ROOTFIND_X_new + u_old;
            
            /* crash condition */
            if((ROOTFIND_ITER >= MAXITER) || isnan(u)) {
                printf("failed to converge in DoCooling(): ROOTFIND_X_new=%g ROOTFIND_X_a=%g ROOTFIND_X_b=%g ROOTFIND_X_error=%g u_in=%g u_upper=%g u_lower=%g rho_in=%g dt=%g ne_in=%g ne_out=%g target=%d ID=%ld \n",ROOTFIND_X_new, ROOTFIND_X_a, ROOTFIND_X_b,  ROOTFIND_X_error, u_old, u_upper, u_lower, rho,dt,ne_guess,*ne_eval,target, (long)P[target].ID); endrun(10);
            }
            u = DMAX(u_min,u);
        } else {u = All.MinEgySpec;}
    }

    double specific_energy_codeunits_toreturn = u / UNIT_SPECEGY_IN_CGS;    /* in internal units */
    SphP[target].Ne = *ne_eval;
#ifdef RT_CHEM_PHOTOION
    /* set variables used by RT routines; this must be set only -outside- of iteration, since this is the key chemistry update */
    double u_in=specific_energy_codeunits_toreturn, rho_in=SphP[target].Density*All.cf_a3inv, mu=1, temp, ne=1, nHI=SphP[target].HI, nHII=SphP[target].HII, nHeI=1, nHeII=0, nHeIII=0;
    temp = ThermalProperties(u_in, rho_in, target, &mu, &ne, &nHI, &nHII, &nHeI, &nHeII, &nHeIII);
    SphP[target].HI = nHI; SphP[target].HII = nHII;
#ifdef RT_CHEM_PHOTOION_HE
    SphP[target].HeI = nHeI; SphP[target].HeII = nHeII; SphP[target].HeIII = nHeIII;
#endif
#endif

    /* safe return */
    return specific_energy_codeunits_toreturn;
#endif // CHIMES
}



#ifndef CHIMES
/* returns cooling time.
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime(double u_old, double rho, double ne_guess, double *ne_eval, int target)
{
#if defined(COOL_GRACKLE) && !defined(GALSF_EFFECTIVE_EQS)
    double LambdaNet = CallGrackle(u_old, rho, 0.0, ne_guess, target, 1);
    if(LambdaNet >= 0) LambdaNet = 0.0;
    return LambdaNet / UNIT_TIME_IN_CGS;
#else
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    double LambdaNet = CoolingRateFromU(u_old, rho, ne_guess, ne_eval, target);
    if(LambdaNet >= 0) {return 0;} /* net heating due to UV background */
    return u_old / (-(nHcgs * nHcgs / rho) * LambdaNet) / UNIT_TIME_IN_CGS;
#endif
}


/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double ne_guess, double *ne_eval, int target)
{
    if(fac <= 0) {return 0.01*m_old;} /* the hot phase is actually colder than the cold reservoir! */
    double m, dm, m_lower, m_upper, ratefact, LambdaNet;
    int iter = 0;

    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u *= UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    fac /= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho * fac;
    m = m_old; m_lower = m; m_upper = m;
    LambdaNet = CoolingRateFromU(u, rho, ne_guess, ne_eval, target);

    /* bracketing */
    if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt < 0)	/* heating */
    {
        m_upper *= sqrt(1.1); m_lower /= sqrt(1.1);
        while(m_upper - m_old - m_upper * m_upper / m_old * ratefact * CoolingRateFromU(u, rho * m_upper / m_old, ne_guess, ne_eval, target) * dt < 0)
        {
            m_upper *= 1.1; m_lower *= 1.1;
        }
    }
    if(m - m_old - m_old * ratefact * LambdaNet * dt > 0)
    {
        m_lower /= sqrt(1.1); m_upper *= sqrt(1.1);
        while(m_lower - m_old - m_lower * m_lower / m_old * ratefact * CoolingRateFromU(u, rho * m_lower / m_old, ne_guess, ne_eval, target) * dt > 0)
        {
            m_upper /= 1.1; m_lower /= 1.1;
        }
    }

    do
    {
        m = 0.5 * (m_lower + m_upper);
        LambdaNet = CoolingRateFromU(u, rho * m / m_old, ne_guess, ne_eval, target);
        if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt > 0) {m_upper = m;} else {m_lower = m;}
        dm = m_upper - m_lower;
        iter++;
        if(iter >= (MAXITER - 10)) {printf("->m= %g\n", m);}
    }
    while(fabs(dm / m) > 1.0e-6 && iter < MAXITER);
    if(iter >= MAXITER) {printf("failed to converge in DoInstabilityCooling(): m_in=%g u_in=%g rho=%g dt=%g fac=%g ne_in=%g target=%d ID=%ld\n",m_old,u,rho,dt,fac,ne_guess,target,(long)P[target].ID); endrun(11);}
    return m;
}

#endif // !(CHIMES)






#ifdef CHIMES
/* This function converts thermal energy to temperature, using the mean molecular weight computed from the non-equilibrium CHIMES abundances. */
double chimes_convert_u_to_temp(double u, double rho, int target)
{
  return u * (GAMMA(target)-1) * PROTONMASS_CGS * ((double) calculate_mean_molecular_weight(&(ChimesGasVars[target]), &ChimesGlobalVars)) / BOLTZMANN_CGS;
}
// CHIMES
#elif defined(EOS_SUBSTELLAR_ISM)

#define NUM_SPECIES_IN_EOS 5
/*
Given the temperature, determine the internal energy assuming local thermodynamic equilibrium between the species present.

Returns the internal energy of the gas mixture per unit mass in erg/g:

u = Etot/Mtot = SUM_i(N_i E_i) / SUM_i(N_i m_i) over all species i with energy E_i, mass m_i, and relative abundance N_i
 */
double convert_temp_to_u(double temp, double rho, int target, double *ne, double *nH0, double *nHp, double *nHe0, double *nHep, double *nHepp, double *mu) {
    double dummy;
    find_abundances_and_rates(log10(temp), rho, target, -1, 0, ne, nH0, nHp, nHe0, nHep, nHepp, mu, &dummy, &dummy, &dummy,
                              &dummy); // all the thermo variables for this T
    double X = HYDROGEN_MASSFRAC, Y = 1. - X, Z = 0, fmol;
#ifdef METALS
    if (target >= 0) {
        Z = DMIN(0.25, P[target].Metallicity[0]);
        if (NUM_METAL_SPECIES >= 10) {
            Y = DMIN(0.35, P[target].Metallicity[1]);
        }
        X = 1. - (Y + Z);
    }
#endif
    double urad_from_uvb_in_G0 = MIN_REAL_NUMBER;                                         // pass this eventually
    fmol = Get_Gas_Molecular_Mass_Fraction(target, temp, *nH0, *ne, urad_from_uvb_in_G0); /* use our simple subroutine to estimate this, ignoring UVB and with clumping factor=1 */

    /* For full generality, make arrays of species' mean energy, masses, and abundances. Indices: 0: H_2 1: H 2: He 3: e 4: metals */
    double E_i[NUM_SPECIES_IN_EOS] = {0}, m_i[NUM_SPECIES_IN_EOS] = {0}, N_i[NUM_SPECIES_IN_EOS] = {0};
    double e_mono = 1.5 * BOLTZMANN_CGS * temp;
    E_i[1] = E_i[2] = E_i[3] = E_i[4] = e_mono;

    E_i[0] = e_mono;
    if (fmol > MIN_REAL_NUMBER) {
        E_i[0] = hydrogen_molecule_energy(temp);
    }
    m_i[0] = 2.;
    N_i[0] = 0.5 * X * fmol;

    // all mono species have 3/2 kB T
    m_i[1] = 1.;
    N_i[1] = X * (1 - fmol); // H
    m_i[2] = 4.;
    N_i[2] = 0.25 * Y; // He
    m_i[3] = ELECTRONMASS_CGS / PROTONMASS_CGS;
    N_i[3] = *ne * X; // e
    m_i[4] = 16. + 12. * fmol;
    N_i[4] = Z / (16. + 12. * fmol); // metals, consistent with mean molecular weight calculation

    int k;
    double sum_E = 0., sum_M = 0.;
    for (k = 0; k < NUM_SPECIES_IN_EOS; k++) {
        sum_E += N_i[k] * E_i[k];
        sum_M += N_i[k] * m_i[k];
    }
    return sum_E / sum_M / PROTONMASS_CGS;
}

/*
   EOS_SUBSTELLAR_ISM implementation of convert_u_to_temp: we call convert_temp_to_u
   and root-find to solve it = u
*/
double convert_u_to_temp(double u, double rho, int target, double *ne, double *nH0, double *nHp, double *nHe0, double *nHep, double *nHepp, double *mu) {
    double T_guess = u * PROTONMASS_CGS / BOLTZMANN_CGS;
    double T_max = 3 * u * PROTONMASS_CGS / BOLTZMANN_CGS, T_min = 0.5 * u * PROTONMASS_CGS / BOLTZMANN_CGS;
    if (All.Time == 0) {
        T_min = pow(10., Tmin);
        T_max = pow(10., Tmax);
    }
#define ROOTFIND_FUNCTION_INNER(temp) convert_temp_to_u(temp, rho, target, ne, nH0, nHp, nHe0, nHep, nHepp, mu) - u
    double ROOTFIND_X_b = T_guess, ROOTFUNC_b = ROOTFIND_FUNCTION_INNER(T_guess), ROOTFIND_X_a = T_min, ROOTFUNC_a = ROOTFIND_FUNCTION_INNER(T_min);
    if (ROOTFUNC_a * ROOTFUNC_b > 0) {
        ROOTFUNC_a = ROOTFIND_FUNCTION_INNER(T_max);
        ROOTFIND_X_a = T_max;
    } // make sure it's bracketed
    double ROOTFIND_REL_X_tol = 1e-4, ROOTFIND_ABS_X_tol = 0;
#include "../system/bracketed_rootfind.h"
    double temp = ROOTFIND_X_new, logtemp = log10(temp);
    if (temp <= 0) {
        return pow(10.0, Tmin);
    }
    if (logtemp < Tmin) {
        return pow(10.0, Tmin);
    }
    if (logtemp > Tmax) {
        return pow(10., Tmax);
    }
    return temp;
}
// elif defined(EOS_SUBSTELLAR_ISM)
#else 

/* this function determines the electron fraction, and hence the mean molecular weight. With it arrives at a self-consistent temperature.
 * Ionization abundances and the rates for the emission are also computed */
double convert_u_to_temp(double u, double rho, int target, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess)
{
    int iter = 0;
    double temp, temp_old, temp_old_old = 0, temp_new, prefac_fun_old, prefac_fun, fac, err_old, err_new, T_bracket_errneg = 0, T_bracket_errpos = 0, T_bracket_min = 0, T_bracket_max = 1.e20, bracket_sign = 0, Lambda_filler = 0; // double max = 0;
    double u_input = u, rho_input = rho, temp_guess;
    double T_0 = u * PROTONMASS_CGS / BOLTZMANN_CGS; // this is the dimensional temperature, which since u is fixed is -frozen- in this calculation: we can work dimensionlessly below
    temp_guess = (GAMMA(target)-1) * T_0; // begin assuming mu ~ 1
    *mu_guess = Get_Gas_Mean_Molecular_Weight_mu(temp_guess, rho, nH0_guess, ne_guess, 0., target); // get mu with that temp
    prefac_fun = (GAMMA(target)-1) * (*mu_guess); // dimensionless pre-factor determining the temperature
    err_new = prefac_fun - temp_guess / T_0; // define initial error from this iteration
    if(err_new < 0) {T_bracket_errneg = temp_guess;} else {T_bracket_errpos = temp_guess;}
    temp = prefac_fun * T_0; // re-calculate temp with the new mu

    do
    {
        //qfun_old = *ne_guess; // guess for ne
        //qfun_old = *mu_guess; // guess for mu
        prefac_fun_old = prefac_fun;
        err_old = err_new; // error from previous timestep
        find_abundances_and_rates(log10(temp), rho, target, -1, 0, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess, &Lambda_filler, &Lambda_filler, &Lambda_filler, &Lambda_filler); // all the thermo variables for this T
        prefac_fun = (GAMMA(target)-1) * (*mu_guess); // new value of the dimensionless pre-factor we need to solve
        temp_old = temp; // guess for T we just used
        temp_new = prefac_fun * T_0; // updated temp using the new values from the iteration of find_abundances_and_rates above
        err_new = (temp_new - temp_old) / T_0; // new error
        if(T_bracket_errpos == 0) {if(err_new > 0) {T_bracket_errpos = temp_old;} else {T_bracket_errneg = temp_old;}} // update the bracket values to the new T while its error still reflects here
        if(T_bracket_errneg == 0) {if(err_new < 0) {T_bracket_errneg = temp_old;} else {T_bracket_errpos = temp_old;}} // update the bracket values to the new T while its error still reflects here
        if(T_bracket_errneg > 0 && T_bracket_errpos > 0)
        {
            if(bracket_sign == 0) {if(T_bracket_errpos > T_bracket_errneg) {bracket_sign=1;} else {bracket_sign=-1;}}
            if(err_new > 0) {
                if(bracket_sign > 0) {T_bracket_errpos = DMIN(T_bracket_errpos, temp_old); /* Tpos>Tneg */} else {T_bracket_errpos = DMAX(T_bracket_errpos, temp_old); /* Tpos<Tneg */}
            } else {
                if(bracket_sign > 0) {T_bracket_errneg = DMAX(T_bracket_errneg, temp_old); /* Tpos>Tneg */} else {T_bracket_errneg = DMIN(T_bracket_errneg, temp_old); /* Tpos<Tneg */}
            } /* update bracket values if we can */
            if(bracket_sign > 0) {T_bracket_max=T_bracket_errpos; T_bracket_min=T_bracket_errneg;} else {T_bracket_max=T_bracket_errneg; T_bracket_min=T_bracket_errpos;}
        }

        //max = DMAX(max, temp_new * (*mu_guess) * HYDROGEN_MASSFRAC * fabs((*ne_guess - qfun_old) / (temp_new - temp_old + 1.0))); // old iteration: hardwired assumption that ne is only varying quanity in mu, and that Tmin ~ 1e4 or so
        //max = DMAX(max , temp_new / (*mu_guess) * fabs(*mu_guess - qfun_old) / (fabs(temp_new - temp_old) + 1.e-4*(All.MinGasTemp+0.1))); // newer - more flexible mu, and dimensionless T dependence
        //temp = temp_old + (temp_new - temp_old) / (1 + max);
        
        if((fabs(prefac_fun-prefac_fun_old) < 1.e-4) && (fabs(temp_new-temp_old)/(temp_new+temp_old) < 1.e-4)) {break;} // break pre-emptively if we'll trigger a nan below
        fac = (prefac_fun-prefac_fun_old)*T_0 / (temp_old-temp_old_old); // numerical derivative factor: want to use this to limit for convergence

        if(fac > 1) {fac = 1;} // don't allow us to move in the opposite direction from the new evaluation (should 'guess' in the direction of T_new-T_old) -- this tells us Newton-Raphson/Secant-type method fails here, so we simply follow the iteration to t_new
        if(fac > 0.9) {fac=0.9;} // don't allow us to 'jump' by a factor >10 times the temperature difference (arbitrary choice, slows convergence a bit but helps limit bad overshoot)
        if(fac < -9999.) {fac=-9999.;} // don't allow smaller step than 1e-4 times the temperature difference (since that's below our error tolerance anyways)

        temp = temp_old + (temp_new - temp_old) * 1./(1. - fac); // standard Newton-Raphson-type (technically Secant method) iteration
        if(temp < 0.5*temp_old) {temp = 0.5*temp_old;} // limiter to prevent un-physical overshoot before we have bracketing established
        if(temp > 3.0*temp_old) {temp = 3.0*temp_old;} // limiter to prevent un-physical overshoot before we have bracketing established
        
        if(T_bracket_errneg > 0 && T_bracket_errpos > 0) // if have bracketing and this wants to go outside brackets, revert to bisection
        {
            if(temp >= T_bracket_max || temp <= T_bracket_min) {temp = sqrt(T_bracket_min*T_bracket_max);} // bisect (in log-space)
        }
        /*
#ifndef RT_INFRARED        
        if(fabs(temp-temp_old_old)/(temp+temp_old_old) < 1.e-3) {double wt=get_random_number(12*iter+340*ThisTask+5435*target); temp=(wt*temp_old + (1.-wt)*temp_new);}
#endif
        */
        temp_old_old = temp_old;
        iter++;
        if(iter > (MAXITER - 10)) {printf("-> temp_next/new/old/oldold=%g/%g/%g/%g ne=%g mu=%g rho=%g iter=%d target=%d err_new/prev=%g/%g gamma_minus_1_mu_new/prev=%g/%g Brackets: Error_bracket_positive=%g Error_bracket_negative=%g T_bracket_Min/Max=%g/%g fac_for_SecantDT=%g \n", temp,temp_new,temp_old,temp_old_old,*ne_guess, (*mu_guess) ,rho,iter,target,err_new,err_old,prefac_fun,prefac_fun_old,T_bracket_errpos,T_bracket_errneg,T_bracket_min,T_bracket_max,fac); fflush(stdout);}
    }
    while(
#if defined(RT_INFRARED)
        (fabs(temp - temp_old) > 1e-3 * temp) && iter < MAXITER);
#else
          ((fabs(temp - temp_old) > 0.25 * temp) ||
           ((fabs(temp - temp_old) > 0.1 * temp) && (temp > 20.)) ||
           ((fabs(temp - temp_old) > 0.05 * temp) && (temp > 200.)) ||
           ((fabs(temp - temp_old) > 0.01 * temp) && (temp > 200.) && (iter<100)) ||
           ((fabs(temp - temp_old) > 1.0e-3 * temp) && (temp > 200.) && (iter<10))) && iter < MAXITER);
#endif
    if(iter >= MAXITER) {printf("failed to converge in convert_u_to_temp(): u_input= %g rho_input=%g n_elec_input=%g target=%d ID=%ld\n", u_input, rho_input, *ne_guess, target, (long)P[target].ID); endrun(12);}

    if(temp<=0) temp=pow(10.0,Tmin);
    if(log10(temp)<Tmin) temp=pow(10.0,Tmin);
    return temp;
}
#endif // CHIMES




#ifndef CHIMES
/* this function computes the actual ionization states, relative abundances, and returns the ionization/recombination rates if needed */
double find_abundances_and_rates(double logT, double rho, int target, double shieldfac, int return_cooling_mode,
                                 double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess,
                                 double *mu_guess, double *LambdaExc_return, double *LambdaIon_return, double *LambdaRec_return, double *LambdaFF_return)
{
    int j, jp, niter;
    double Tlow, Thi, flow, fhi, t, gJH0ne, gJHe0ne, gJHepne, logT_input, rho_input, ne_input, neold, nenew;
    double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep, EPSILON_SMALL=1.e-40;
    double n_elec, nH0, nHe0, nHp, nHep, nHepp; /* ionization states */
    logT_input = logT; rho_input = rho; ne_input = *ne_guess; /* save inputs (in case of failed convergence below) */
    if(!isfinite(logT)) {logT=Tmin;}    /* nan trap (just in case) */
    if(!isfinite(rho)) {logT=Tmin;}

    if(logT <= Tmin)		/* everything neutral */
    {
        nH0 = 1.0; nHe0 = yhelium(target); nHp = 0; nHep = 0; nHepp = 0; n_elec = 1.e-22;
        *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec;
        *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, 0, target);
        return 0;
    }
    if(logT >= Tmax)		/* everything is ionized */
    {
        nH0 = 0; nHe0 = 0; nHp = 1.0; nHep = 0; nHepp = yhelium(target); n_elec = nHp + 2.0 * nHepp;
        *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec;
        *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, 1.e3, target);
        return 0;
    }

    /* initialize quantities needed for iteration below */
    t = (logT - Tmin) / deltaT;
    j = (int) t;
    if(j<0) {j=0;}
    if(j>NCOOLTAB){
        PRINT_WARNING("j>=NCOOLTAB : j=%d t %g Tlow %g Thi %g logT %g Tmin %g deltaT %g \n",j,t,Tmin+deltaT*j,Tmin+deltaT*(j+1),logT,Tmin,deltaT);fflush(stdout);
        j=NCOOLTAB;
    }
    jp = j + 1;
    if(jp > NCOOLTAB) {jp=NCOOLTAB;}
    Tlow = Tmin + deltaT * j;
    Thi = Tlow + deltaT;
    fhi = t - j;
    flow = 1 - fhi;
    if(*ne_guess == 0) /* no guess provided, try to start from something sensible */
    {
        *ne_guess = 1.0;
        if(logT < 3.8) {*ne_guess = 0.1;}
        if(logT < 2) {*ne_guess = 1.e-10;}
    }
    /* CAFG: this is the density that we should use for UV background threshold */
    double local_gammamultiplier = return_local_gammamultiplier(target); // account for local UVB terms in some expressions below
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    if(shieldfac < 0) {shieldfac = return_uvb_shieldfac(target, local_gammamultiplier*gJH0/1.0e-12, nHcgs, logT);} // if < 0, that's a key to tell us this needs to be recalculated
    n_elec = *ne_guess; if(!isfinite(n_elec)) {n_elec=1;}
    neold = n_elec; niter = 0;
    double dt = 0, fac_noneq_cgs = 0, necgs = n_elec * nHcgs, ne_lower=0, ne_upper=2.; /* more initialized quantities */
    int bisection_mode=0; // 0 if doing the usual fixed-point iteration; 1 if switched to bisection method
    if(target >= 0) {dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target);} // dtime [code units]
    fac_noneq_cgs = (dt * UNIT_TIME_IN_CGS) * (necgs + 1.e-30*nHcgs); // factor needed below to asses whether timestep is larger/smaller than recombination time

#if defined(RT_CHEM_PHOTOION)
    double c_light_ne=0;
    if(target >= 0)
    {
        nH0 = SphP[target].HI; // need to initialize a value for the iteration below
#ifdef RT_CHEM_PHOTOION_HE
        nHe0 = SphP[target].HeI; nHep = SphP[target].HeII; // need to intialize a value for the iteration below
#endif
    }
#endif

    /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
    do
    {
        niter++;
        
        aHp = flow * AlphaHp[j] + fhi * AlphaHp[jp];
        aHep = flow * AlphaHep[j] + fhi * AlphaHep[jp];
        aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[jp];
        ad = flow * Alphad[j] + fhi * Alphad[jp];
        geH0 = flow * GammaeH0[j] + fhi * GammaeH0[jp];
        geH0 = DMAX(geH0, EPSILON_SMALL);
        geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[jp];
        geHe0 = DMAX(geHe0, EPSILON_SMALL);
        geHep = flow * GammaeHep[j] + fhi * GammaeHep[jp];
        geHep = DMAX(geHep, EPSILON_SMALL);
        fac_noneq_cgs = (dt * UNIT_TIME_IN_CGS) * (necgs + 1.e-30*nHcgs); // factor needed below to asses whether timestep is larger/smaller than recombination time
        if(necgs <= 1.e-25 || J_UV == 0)
        {
            gJH0ne = gJHe0ne = gJHepne = 0;
        }
        else
        {
            /* account for self-shielding in calculating UV background effects */
            gJH0ne = gJH0 * local_gammamultiplier / necgs * shieldfac; // check units, should be = c_light * n_photons_vol * rt_ion_sigma_HI[0] / necgs;
            gJH0ne = DMAX(gJH0ne, EPSILON_SMALL); if(!isfinite(gJH0ne)) {gJH0ne=0;} // need traps here b/c very small numbers assigned in some newer TREECOOL versions cause a nan underflow
            gJHe0ne = gJHe0 * local_gammamultiplier / necgs * shieldfac;
            gJHe0ne = DMAX(gJHe0ne, EPSILON_SMALL); if(!isfinite(gJHe0ne)) {gJHe0ne=0;}
            gJHepne = gJHep * local_gammamultiplier / necgs * shieldfac;
            gJHepne = DMAX(gJHepne, EPSILON_SMALL); if(!isfinite(gJHepne)) {gJHepne=0;}
        }
#if defined(RT_DISABLE_UV_BACKGROUND)
        gJH0ne = gJHe0ne = gJHepne = 0;
#endif
#if defined(RT_CHEM_PHOTOION)
        /* add in photons from explicit radiative transfer (on top of assumed background) */
        if(target >= 0)
        {
            int k;
            c_light_ne = C_LIGHT_CGS / ((MIN_REAL_NUMBER + necgs) * UNIT_LENGTH_IN_CGS); // want physical cgs units for quantities below
            double gJH0ne_0=gJH0 * local_gammamultiplier / (MIN_REAL_NUMBER + necgs), gJHe0ne_0=gJHe0 * local_gammamultiplier / (MIN_REAL_NUMBER + necgs), gJHepne_0=gJHep * local_gammamultiplier / (MIN_REAL_NUMBER + necgs); // need a baseline, so we don't over-shoot below
            gJH0ne = DMAX(gJH0ne, EPSILON_SMALL); if(!isfinite(gJH0ne)) {gJH0ne=0;} // need traps here b/c very small numbers assigned in some newer TREECOOL versions cause a nan underflow
            gJHe0ne = DMAX(gJHe0ne, EPSILON_SMALL); if(!isfinite(gJHe0ne)) {gJHe0ne=0;}
            gJHepne = DMAX(gJHepne, EPSILON_SMALL); if(!isfinite(gJHepne)) {gJHepne=0;}
#if defined(RT_DISABLE_UV_BACKGROUND)
            gJH0ne_0=gJHe0ne_0=gJHepne_0=MAX_REAL_NUMBER;
#endif
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if(RT_BAND_IS_IONIZING(k))
                {
                    double n_gamma_tot = rt_return_photon_number_density(target,k);
#ifdef RT_INFRARED
                    n_gamma_tot += rt_irband_egydensity_in_band(target,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]) / (DMAX(rt_nu_eff_eV[RT_FREQ_BIN_H0],SphP[target].Radiation_Temperature/2959.81)*ELECTRONVOLT_IN_ERGS/UNIT_ENERGY_IN_CGS);
#endif
                    double c_ne_time_n_photons_vol = c_light_ne * n_gamma_tot; // gives photon flux
                    double cross_section_ion, dummy, thold=1.0e20;
#ifdef GALSF
                    if(All.ComovingIntegrationOn) {thold=1.0e10;}
#endif
                    if(rt_ion_G_HI[k] > 0)
                    {
                        cross_section_ion = nH0 * rt_ion_sigma_HI[k];
                        dummy = rt_ion_sigma_HI[k] * c_ne_time_n_photons_vol;// egy per photon x cross section x photon flux (w attenuation factors already included in flux/energy update:) * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        //if(dummy > thold*gJH0ne_0) {dummy = thold*gJH0ne_0;}
                        gJH0ne += dummy;
                    }
#if 1 //def RT_CHEM_PHOTOION_HE
                    if(rt_ion_G_HeI[k] > 0)
                    {
                        cross_section_ion = nHe0 * rt_ion_sigma_HeI[k];
                        dummy = rt_ion_sigma_HeI[k] * c_ne_time_n_photons_vol;// * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        //if(dummy > thold*gJHe0ne_0) {dummy = thold*gJHe0ne_0;}
                        gJHe0ne += dummy;
                    }
                    if(rt_ion_G_HeII[k] > 0)
                    {
                        cross_section_ion = nHep * rt_ion_sigma_HeII[k];
                        dummy = rt_ion_sigma_HeII[k] * c_ne_time_n_photons_vol;// * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        //if(dummy > thold*gJHepne_0) {dummy = thold*gJHepne_0;}
                        gJHepne += dummy;
                    }
#endif
                }
            }
        }
#endif


        nH0 = aHp / (MIN_REAL_NUMBER + aHp + geH0 + gJH0ne);	/* eqn (33) */
#if defined(RT_CHEM_PHOTOION)
        if(target >= 0) {nH0 = (SphP[target].HI + fac_noneq_cgs * aHp) / (1 + fac_noneq_cgs * (aHp + geH0 + gJH0ne));} // slightly more general formulation that gives linear update but interpolates to equilibrium solution when dt >> dt_recombination
#endif // note that the above can produce inconsistencies if used without He ionization routines below, because it will still assume equilibrium ionization of He which can lead to a much higher-than-expected free-electron fraction relative to the ionized H fraction, so be careful
        nHp = 1.0 - nH0;		/* eqn (34) */

        if( ((gJHe0ne + geHe0) <= MIN_REAL_NUMBER) || (aHepp <= MIN_REAL_NUMBER) ) 	/* no ionization at all */
        {
            nHep = 0.0;
            nHepp = 0.0;
            nHe0 = yhelium(target);
        }
        else
        {
            nHep = yhelium(target) / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
            nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
            nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
        }
#if defined(RT_CHEM_PHOTOION) && defined(RT_CHEM_PHOTOION_HE)
        if(target >= 0)
        {
            double yHe = yhelium(target); // will use helium fraction below
            nHep = SphP[target].HeII + yHe * fac_noneq_cgs * (geHe0 + gJHe0ne) - SphP[target].HeIII * (fac_noneq_cgs*(geHe0 + gJHe0ne - aHepp) / (1.0 + fac_noneq_cgs*aHepp));
            nHep /= 1.0 + fac_noneq_cgs*(geHe0 + gJHe0ne + aHep + ad + geHep + gJHepne) + (fac_noneq_cgs*(geHe0 + gJHe0ne - aHepp) / (1.0 + fac_noneq_cgs*aHepp)) * fac_noneq_cgs*(geHep + gJHepne);
            if(nHep < 0) {nHep=0;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            if(nHep > yHe) {nHep=yHe;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            nHepp = (SphP[target].HeIII + SphP[target].HeII * fac_noneq_cgs*(geHep + gJHepne)) / (1. + fac_noneq_cgs*aHepp);
            if(nHepp < 0) {nHepp=0;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            if(nHepp > yHe-nHep) {nHepp=yHe-nHep;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            nHe0 = yHe - (nHep + nHepp); // remainder is neutral
        }
#endif
        if(!isfinite(n_elec)) {printf("target=%d niter=%d logT=%g n_elec/old=%g/%g nHp/nHep/nHepp=%g/%g/%g nHcgs=%g yHe=%g dt=%g shieldfac/local_gammamult=%g/%g aHp/aHep/aHepp=%g/%g/%g geH0/geHe0/geHep=%g/%g/%g gJH0ne/gJHe0ne/gJHepne=%g/%g/%g \n",target,niter,logT,n_elec,neold,nHp,nHep,nHepp,nHcgs,yhelium(target),dt,shieldfac,local_gammamultiplier,aHp,aHep,aHepp,geH0,geHe0,geHep,gJH0ne,gJHe0ne,gJHepne);}

        double error_old = fabs(n_elec - neold); // save the old error to compare with the new one, so we can check if we're converging
        neold = n_elec;
        n_elec = nHp + nHep + 2 * nHepp;	/* eqn (38) */
#ifdef COOL_LOW_TEMPERATURES
        n_elec += return_electron_fraction_from_heavy_ions(target, pow(10.,logT), rho, n_elec);
#endif       
	
	// keep track of these bounds in case we need to switch to bisection
	if(n_elec > neold) {ne_lower = DMAX(neold, ne_lower);}
	if(n_elec < neold) {ne_upper = DMIN(neold, ne_upper);}

        double nenew_tolmin = DMIN(1.0e-14, 0.01 * 1.e-3/nHcgs); // make tolerance 1e-14 (normally more than sufficient), unless super-dense where smaller tolerance needed owing to expectation [-guessed- here for physics] of even smaller ne needed
	if(bisection_mode) { // if fixed-point mode is not converging fast enough and we switched to bisection mode
	    if(n_elec < neold) {nenew = 0.5*(ne_lower + neold); ne_upper=neold;} // go to left midpoint and update the upper bound
	    else {nenew = 0.5*(ne_upper + neold); ne_lower = neold;} // go to right midpoint and update the lower bound
        } else { // otherwise we do the usual fixed-point iteration
	    nenew = 0.5 * (n_elec + neold);
	    if(niter>30 && (fabs(n_elec - neold) > 0.6 * error_old) && (nenew > nenew_tolmin)) { bisection_mode = 1;}  // if we're converging slower than bisection, just switch to bisection
	} 

        n_elec = nenew;
        if(!isfinite(n_elec)) {n_elec=1;}
        necgs = n_elec * nHcgs;

        double dneTHhold = DMAX(n_elec*0.01 , nenew_tolmin); // desired absolute tolerance for n_elec
        if(fabs(n_elec - neold) < dneTHhold) break;

        if(niter > (MAXITER - 10)) {printf("n_elec= %g/%g/%g yh=%g nHcgs=%g niter=%d\n", n_elec,neold,nenew, yhelium(target), nHcgs, niter);}
    }
    while(niter < MAXITER);

    if(niter >= MAXITER) {printf("failed to converge in find_abundances_and_rates(): logT_input=%g  rho_input=%g  ne_input=%g target=%d ID=%ld shieldfac=%g cooling_return=%d", logT_input, rho_input, ne_input, target, (long)P[target].ID, shieldfac, return_cooling_mode); endrun(13);}

    bH0 = flow * BetaH0[j] + fhi * BetaH0[jp];
    bHep = flow * BetaHep[j] + fhi * BetaHep[jp];
    bff = flow * Betaff[j] + fhi * Betaff[jp];
    *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec; /* write to send back */
    *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, sqrt(shieldfac)*(gJH0/2.29e-10), target);
    if(target >= 0) /* if this is a cell, update some of its thermodynamic stored quantities */
    {
#if defined(OUTPUT_MOLECULAR_FRACTION)
        SphP[target].MolecularMassFraction = Get_Gas_Molecular_Mass_Fraction(target, pow(10.,logT), nH0, n_elec, sqrt(shieldfac)*(gJH0/2.29e-10));
#endif
    }

    /* now check if we want to return the ionization/recombination heating/cooling rates calculated with all the above quantities */
    if(return_cooling_mode==1)
    {
        /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
        double LambdaExcH0 = bH0 * n_elec * nH0;
        double LambdaExcHep = bHep * n_elec * nHep;
        double LambdaExc = LambdaExcH0 + LambdaExcHep;	/* collisional excitation */

        double LambdaIonH0 = 2.18e-11 * geH0 * n_elec * nH0;
        double LambdaIonHe0 = 3.94e-11 * geHe0 * n_elec * nHe0;
        double LambdaIonHep = 8.72e-11 * geHep * n_elec * nHep;
        double LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* collisional ionization */

        double T_lin = pow(10.0, logT);
        double LambdaRecHp = 1.036e-16 * T_lin * n_elec * (aHp * nHp);
        double LambdaRecHep = 1.036e-16 * T_lin * n_elec * (aHep * nHep);
        double LambdaRecHepp = 1.036e-16 * T_lin * n_elec * (aHepp * nHepp);
        double LambdaRecHepd = 6.526e-11 * ad * n_elec * nHep;
        double LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd; /* recombination */

        double LambdaFF = bff * (nHp + nHep + 4 * nHepp) * n_elec; /* free-free (Bremsstrahlung) */
        double Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF; /* sum all of the above */

        *LambdaExc_return = LambdaExc; *LambdaIon_return = LambdaIon; *LambdaRec_return = LambdaRec; *LambdaFF_return = LambdaFF; /* send data back if needed */
        return Lambda; /* send it back */
    }
    return 0;
} // end of find_abundances_and_rates() //



/*  this function first computes the self-consistent temperature and abundance ratios, and then it calculates (heating rate-cooling rate)/n_h^2 in cgs units */
double CoolingRateFromU(double u, double rho, double ne_guess, double *ne_eval, int target)
{
    double nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu; nH0_guess = DMAX(0,DMIN(1,1.-ne_guess/1.2));
    double temp = convert_u_to_temp(u, rho, target, &ne_guess, &nH0_guess, &nHp_guess, &nHe0_guess, &nHep_guess, &nHepp_guess, &mu);
    double Lambda = CoolingRate(log10(temp), rho, ne_guess, ne_eval, target);
    return Lambda;
}


#endif // !(CHIMES)



extern FILE *fd;



#ifndef CHIMES
/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units
 */
double CoolingRate(double logT,  double rho, double n_elec_guess, double *n_elec_eval, int target)
{
    double n_elec=n_elec_guess, nH0, nHe0, nHp, nHep, nHepp, mu; /* ionization states [computed below] */
    double Lambda, Heat, LambdaFF, LambdaCompton, LambdaExc, LambdaExcH0, LambdaExcHep, LambdaIon, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
    double LambdaRec, LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd, T, T_cmb_radeff, shieldfac, LambdaMol, LambdaMetal, LambdaPElec, LambdaDust;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS;	/* hydrogen number dens in cgs units */
    Lambda=0; Heat=0; LambdaMol=0; LambdaFF=0; LambdaRec=0; LambdaExc=0; LambdaIon=0; LambdaMetal=0; LambdaCompton=0; LambdaPElec=0; LambdaDust=0; /* make sure these are all initialized to zero */
    if(logT <= Tmin) {logT = Tmin + 0.5 * deltaT;}	/* floor at Tmin */
    if(!isfinite(rho)) {return 0;}
    T = pow(10.0, logT);
    T_cmb_radeff = get_background_radiation_temperature_for_emission_corrections(target); /* CMB temperature, used below */

    /* some blocks below to define useful variables before calculation of cooling rates: */

#ifdef COOL_METAL_LINES_BY_SPECIES
    double *Z;
    if(target>=0)
    {
        Z = P[target].Metallicity;
#if defined(GALSF_ISMDUSTCHEM_MODEL) && !defined(GALSF_ISMDUSTCHEM_PASSIVE)
        int k; for(k=0;k<NUM_ISMDUSTCHEM_ELEMENTS;k++) {Z[k] = DMAX(0.,P[target].Metallicity[k]-SphP[target].ISMDustChem_Dust_Metal[k]);}
#endif
    } else { /* initialize dummy values here so the function doesn't crash, if called when there isn't a target particle */
        int k; double Zsol[NUM_METAL_SPECIES]; for(k=0;k<NUM_METAL_SPECIES;k++) {Zsol[k]=All.SolarAbundances[k];}
        Z = Zsol;
    }
#endif
    double local_gammamultiplier = return_local_gammamultiplier(target);
    shieldfac = return_uvb_shieldfac(target, local_gammamultiplier*gJH0/1.0e-12, nHcgs, logT);
    
#if defined(COOL_LOW_TEMPERATURES)
    double Tdust = 30.; /* set variables needed for dust heating/cooling. if dust cooling not calculated, default to 0 */
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) && (FLAG_NOT_IN_PUBLIC_CODE > 2)) || defined(SINGLE_STAR_SINK_DYNAMICS)
    Tdust = get_equilibrium_dust_temperature_estimate(target, shieldfac, T);
#endif
#endif


#if defined(RT_CHEM_PHOTOION) || defined(RT_PHOTOELECTRIC)
    double cx_to_kappa = HYDROGEN_MASSFRAC / PROTONMASS_CGS * UNIT_MASS_IN_CGS; // pre-factor for converting cross sections into opacities
#endif
    if(logT < Tmax)
    {
        /* get ionization states for H and He with associated ionization, collision, recombination, and free-free heating/cooling */
        Lambda += find_abundances_and_rates(logT, rho, target, shieldfac, 1, &n_elec, &nH0, &nHp, &nHe0, &nHep, &nHepp, &mu, &LambdaExc, &LambdaIon, &LambdaRec, &LambdaFF); /* adds all of these to our running total for cooling */
        *n_elec_eval = n_elec; /* save this value for the output cycle */
        LambdaCompton = evaluate_Compton_heating_cooling_rate(target,T,nHcgs,n_elec,shieldfac); /* note this can have either sign: heating or cooling */
        if(LambdaCompton > 0) {Lambda += LambdaCompton;}
        
#ifdef COOL_METAL_LINES_BY_SPECIES
        /* can restrict to low-densities where not self-shielded, but let shieldfac (in ne) take care of this self-consistently */
        if((J_UV != 0)&&(logT > 4.00))
        {
            /* cooling rates tabulated for each species from Wiersma, Schaye, & Smith tables (2008) */
            LambdaMetal = GetCoolingRateWSpecies(nHcgs, logT, Z); //* nHcgs*nHcgs;
            /* tables normalized so ne*ni/(nH*nH) included already, so just multiply by nH^2 */
            /* (sorry, -- dont -- multiply by nH^2 here b/c that's how everything is normalized in this function) */
            LambdaMetal *= n_elec;
            /* (modified now to correct out tabulated ne so that calculated ne can be inserted; ni not used b/c it should vary species-to-species */
            Lambda += LambdaMetal;
#if defined(OUTPUT_COOLRATE_DETAIL)
            if(target >= 0) {SphP[target].MetalCoolingRate = LambdaMetal;}
#endif
        }
#endif

#if defined(GALSF_ISMDUSTCHEM_HIGHTEMPDUSTCOOLING)
        Lambda += Lambda_Dust_HighTemperature_Gas_ISM(target,T,n_elec);
#endif

#ifdef COOL_LOW_TEMPERATURES
        if(logT <= 5.3)
        {
            /* approx to cooling function for solar metallicity and nH=1 cm^(-3) -- want to do something
             much better, definitely, but for now use this just to get some idea of system with cooling to very low-temp */
            LambdaMol = 2.8958629e-26/(pow(T/125.21547,-4.9201887)+pow(T/1349.8649,-1.7287826)+pow(T/6450.0636,-0.30749082));
            LambdaMol *= (1-shieldfac) / (1. + nHcgs/700.); // above the critical density, cooling rate suppressed by ~1/n; use critical density of CO[J(1-0)] as a proxy for this
            double Z_sol=1, truncation_factor=1; /* if don't have actual metallicities, we'll assume solar */
            if(logT>4.5) {double dx=(logT-4.5)/0.20; truncation_factor *= exp(-DMIN(dx*dx,40.));} /* continuous cutoff here just to avoid introducing artificial features in temperature-density */
#ifdef COOL_METAL_LINES_BY_SPECIES
            Z_sol = Z[0] / All.SolarAbundances[0]; /* use actual metallicity for this */
#endif
            LambdaMol *= (1+Z_sol)*(0.001 + 0.1*nHcgs/(1.+nHcgs) + 0.09*nHcgs/(1.+0.1*nHcgs) + Z_sol*Z_sol/(1.0+nHcgs)); // gives very crude estimate of metal-dependent terms //
#if defined(COOL_METAL_LINES_BY_SPECIES) && ((defined(FLAG_NOT_IN_PUBLIC_CODE) && (FLAG_NOT_IN_PUBLIC_CODE > 2)) || !defined(FLAG_NOT_IN_PUBLIC_CODE))
            double column = evaluate_NH_from_GradRho(SphP[target].Gradients.Density,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target) * UNIT_SURFDEN_IN_CGS; // converts to cgs            
            double Z_C = DMAX(1.e-6, Z[2]/All.SolarAbundances[2]), sqrt_T=sqrt(T), ncrit_CO=1.9e4*sqrt_T, Sigma_crit_CO=3.0e-5*T/Z_C, T3=T/1.e3, EXPmax=90.; // carbon abundance (relative to solar and 1/2 factor for original assumed 0.5 depletion), critical density and column
#if defined(GALSF_ISMDUSTCHEM_MODEL) && !defined(GALSF_ISMDUSTCHEM_PASSIVE)
            Z_C = DMAX(1.e-6, Z[2]/(0.5*All.SolarAbundances[2])); // gas-phase carbon abundance (relative to solar/2, usual assumption implicitly)
#endif
            double f_Cplus_CCO=1./(1.+nHcgs/3.e3), photoelec=0; // very crude estimate used to transition between C+ cooling curve and C/CO [nearly-identical] cooling curves above C+ critical density, where C+ rate rapidly declines
#ifdef GALSF_FB_FIRE_RT_UVHEATING
            photoelec = SphP[target].Rad_Flux_UV; if(gJH0>0 && shieldfac>0) {photoelec+=sqrt(shieldfac) * (gJH0/2.29e-10);} // uvb contribution //
#endif
#ifdef RT_PHOTOELECTRIC
            photoelec = SphP[target].Rad_E_gamma[RT_FREQ_BIN_PHOTOELECTRIC] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; photoelec=DMAX(DMIN(photoelec,1.e4),0); // convert to Habing field //
#endif
#if defined(RT_ISRF_BACKGROUND) && (!defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)) // latter flag decides whether we do treecol/sobolev here to get the background intensity
            photoelec += All.InterstellarRadiationFieldStrength * 1.7 * exp(-DMAX(P[target].Metallicity[0]/All.SolarAbundances[0],1e-4) * column * 500.); // RT_ISRF_BACKGROUND rescales the overal ISRF, factor of 1.7 gives Draine 1978 field in Habing units, extinction factor assumes the same FUV band-integrated dust opacity as RT module
#endif
#if !(defined(GALSF_FB_RT_UVHEATING) || defined(RT_PHOTOELECTRIC) || defined(RT_ISRF_BACKGROUND))
            photoelec = 1; // if no explicit modeling of FUV, just assume Habing
#endif
            f_Cplus_CCO = (nHcgs/(340.*DMAX(0.1,photoelec))); f_Cplus_CCO=1./(1.+f_Cplus_CCO*f_Cplus_CCO/sqrt_T); // fco/(1-fco) ~ 0.0022 * ((n/50 cm^-3)/G0)^2 * (100K/T)^(1/2) from Tielens
            double Lambda_Cplus = Z_C * (4.7e-28 * (pow(T,0.15) + 1.04e4*n_elec/sqrt_T) * exp(-DMIN(91.211/T,EXPmax)) + 2.08e-29*exp(-DMIN(23.6/T,EXPmax))); // fit from Barinovs et al., ApJ, 620, 537, 2005, and Wilson & Bell MNRAS 337 1027 2002; assuming factor of 0.5 depletion factor in ISM; rate per C+ relative to solar; + plus [CI]-609 µm line cooling from Hocuk⋆ et al. 2016MNRAS.456.2586H
            double Lambda_CCO = Z_C * T*sqrt_T * 2.73e-31 / (1. + (nHcgs/ncrit_CO)*(1.+1.*DMAX(column,0.017)/Sigma_crit_CO)); // fit from Hollenbach & McKee 1979 for CO (+CH/OH/HCN/OH/HCl/H20/etc., but those don't matter), with slight re-calibration of normalization (factor ~1.4 or so) to better fit the results from the full Glover+Clark network. As Glover+Clark show, if you shift gas out of CO into C+ and O, you have almost no effect on the integrated cooling rate, so this is a surprisingly good approximation without knowing anything about the detailed chemical/molecular state of the gas. uncertainties in e.g. ambient radiation are -much- larger. also note this rate is really carbon-dominated as the limiting abundance, so should probably use that.
	    
	    /* Large-velocity-gradient limiter for the CO cooling rate from Whitworth & Jaffa arXiv:1811.06814; compute Lambda_HI and use this as an upper bound */ 
	    double gradv_norm_kms_pc = velocity_gradient_norm(target) * UNIT_VEL_IN_KMS / UNIT_LENGTH_IN_PC; // velocity gradient in km/s/pc
	    double Lambda_CO_HI = 4.42e-28 * gradv_norm_kms_pc * pow(T,4) / (nHcgs * nHcgs);
	    double Lambda_CO = DMIN(Lambda_CO_HI, (1.-f_Cplus_CCO) * Lambda_CCO); // Let the LVG value be the limiter

            double Lambda_Metals = f_Cplus_CCO * Lambda_Cplus + Lambda_CO; // interpolate between both regimes //
            /* in the above Lambda_Metals expression, the column density expression attempts to account for the optically-thick correction in a slab. this is largely redundant (not exactly, b/c this is specific for CO-type molecules) with our optically-thick cooling module already included below, so we will not double-count it here [coefficient set to zero]. But it's included so you can easily turn it back on, if desired, instead of using the module below. */
            double Lambda_H2_thick = (6.7e-19*exp(-DMIN(5.86/T3,EXPmax)) + 1.6e-18*exp(-DMIN(11.7/T3,EXPmax)) + 3.e-24*exp(-DMIN(0.51/T3,EXPmax)) + 9.5e-22*pow(T3,3.76)*exp(-DMIN(0.0022/(T3*T3*T3),EXPmax))/(1.+0.12*pow(T3,2.1))) / nHcgs; // super-critical H2-H cooling rate [per H2 molecule]
            double Lambda_HD_thin = ((1.555e-25 + 1.272e-26*pow(T,0.77))*exp(-DMIN(128./T,EXPmax)) + (2.406e-25 + 1.232e-26*pow(T,0.92))*exp(-DMIN(255./T,EXPmax))) * exp(-DMIN(T3*T3/25.,EXPmax)); // optically-thin HD cooling rate [assuming all D locked into HD at temperatures where this is relevant], per molecule
            double f_molec = 0.5 * Get_Gas_Molecular_Mass_Fraction(target, T, nH0, n_elec, sqrt(shieldfac)*(gJH0/2.29e-10)); // [0.5*f_molec for H2/HD cooling b/c cooling rates above are per molecule, not per nucleon]

            double q = logT - 3., Y_Hefrac=DMAX(0.,DMIN(1.,Z[1])), X_Hfrac=DMAX(0.,DMIN(1.,1.-Y_Hefrac-Z[0])); // variable used below
            double Lambda_H2_thin = DMAX(nH0-2.*f_molec,0) * X_Hfrac * pow(10., DMAX(-103. + 97.59*logT - 48.05*logT*logT + 10.8*logT*logT*logT - 0.9032*logT*logT*logT*logT , -50.)); // sub-critical H2 cooling rate from H2-H collisions [per H2 molecule]; this from Galli & Palla 1998
            Lambda_H2_thin += Y_Hefrac * pow(10., DMAX(-23.6892 + 2.18924*q -0.815204*q*q + 0.290363*q*q*q -0.165962*q*q*q*q + 0.191914*q*q*q*q*q, -50.)); // H2-He; often more efficient than H2-H at very low temperatures (<100 K); this and other H2-x terms below from Glover & Abel 2008
            Lambda_H2_thin += f_molec * X_Hfrac * pow(10., DMAX(-23.9621 + 2.09434*q -0.771514*q*q + 0.436934*q*q*q -0.149132*q*q*q*q -0.0336383*q*q*q*q*q, -50.)); // H2-H2; can be more efficient than H2-H when H2 fraction is order-unity
            Lambda_H2_thin += nHp * X_Hfrac * pow(10., DMAX(-21.7167 + 1.38658*q -0.379153*q*q + 0.114537*q*q*q -0.232142*q*q*q*q + 0.0585389*q*q*q*q*q, -50.)); // H2-H+; very efficient if somehow appreciable H+ fraction remains
            double logLambdaH2_e = -34.2862 -48.5372*q -77.1212*q*q -51.3525*q*q*q -15.1692*q*q*q*q -0.981203*q*q*q*q*q; // H2-e [generally sub-dominant to H2-H+; can dominate if free e- largely from other sources (e.g. Mg, etc.), but in those conditions essentially impossible for H2 cooling to dominate
            if(logT>2.30103) {logLambdaH2_e = -22.1903 + 1.5729*q -0.213351*q*q + 0.961498*q*q*q -0.910232*q*q*q*q + 0.137497*q*q*q*q*q;}
            Lambda_H2_thin += n_elec * X_Hfrac * pow(10., DMAX(logLambdaH2_e, -50.));

            double f_HD = DMIN(0.00126*f_molec , 4.0e-5*nH0); // ratio of HD molecules to H2 molecules: in low limit, HD easier to form so saturates at about 0.13% of H2 molecules, following Galli & Palla 1998, but obviously cannot exceed the cosmic ratio of D/H=4e-5
            double nH_over_ncrit = Lambda_H2_thin / Lambda_H2_thick , Lambda_HD = f_HD * Lambda_HD_thin / (1. + (f_HD/(f_molec+MIN_REAL_NUMBER))*nH_over_ncrit), Lambda_H2 = f_molec * Lambda_H2_thin / (1. + nH_over_ncrit); // correct cooling rates for densities above critical
            double Lambda_Metals_Neutral = nH0 * Lambda_Metals; // finally note our metal terms here are all for atomic or molecular, not ionized (handled in tables above)
            if(!isfinite(Lambda_Metals_Neutral) || Lambda_Metals_Neutral < 0) {Lambda_Metals_Neutral=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(Lambda_H2) || Lambda_H2 < 0) {Lambda_H2=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(Lambda_HD) || Lambda_HD < 0) {Lambda_HD=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            LambdaMol = Lambda_Metals_Neutral + Lambda_H2 + Lambda_HD; // combine to get total cooling rate
#endif
            LambdaMol *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
            LambdaDust = gas_dust_heating_coeff(target,T,Tdust) * (T-Tdust);// Note our sign convention is such that positive lambda = gas cooling
#if !defined(GALSF_ISMDUSTCHEM_MODEL) || defined(GALSF_ISMDUSTCHEM_PASSIVE)
            if(T>3.e5) {double dx=(T-3.e5)/2.e5; LambdaDust *= exp(-DMIN(dx*dx,40.));} /* needs to truncate at high temperatures b/c of dust destruction (in some modules we solve for this explicitly - in that case can protect this more explicitly, but here, we will make a simple approximation, otherwise we run into problems. note this is not sublimation generally, but sputtering, that causes the destruction */
#endif
            LambdaDust *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
            if(!isfinite(LambdaDust)) {LambdaDust=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(!isfinite(LambdaMol)) {LambdaMol=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            LambdaMol *= ((T-T_cmb_radeff)/(T+T_cmb_radeff)); // account (approximately) for the CMB temperature 'bath' (could more accurately subtract Lambda(T[cmb]), but that's an approximation as well that can give some odd results owing to not treating the solve for molecules indepedently there, so we use this form instead, which is generally good)
            if(LambdaMol > 0) {Lambda += LambdaMol;}
        }
#endif

        Heat = 0;  /* Now, collect heating terms */

        if(J_UV != 0) {Heat += local_gammamultiplier * (nH0 * epsH0 + nHe0 * epsHe0 + nHep * epsHep) / nHcgs * shieldfac;} // shieldfac allows for self-shielding from background
#if defined(RT_DISABLE_UV_BACKGROUND)
        Heat = 0;
#endif
#if defined(RT_CHEM_PHOTOION)
        /* add in photons from explicit radiative transfer (on top of assumed background) */
        if((target >= 0) && (nHcgs > MIN_REAL_NUMBER))
        {
            int k; double c_light_nH = C_LIGHT_CGS / (nHcgs * UNIT_LENGTH_IN_CGS) * UNIT_ENERGY_IN_CGS; // want physical cgs units for quantities below
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if(RT_BAND_IS_IONIZING(k))
                {
                    double n_gamma_tot = rt_return_photon_number_density(target,k);
#ifdef RT_INFRARED
                    n_gamma_tot += rt_irband_egydensity_in_band(target,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]) / (DMAX(rt_nu_eff_eV[RT_FREQ_BIN_H0],SphP[target].Radiation_Temperature/2959.81)*ELECTRONVOLT_IN_ERGS/UNIT_ENERGY_IN_CGS);
#endif
                    double c_nH_time_n_photons_vol = c_light_nH * n_gamma_tot; // gives photon flux
                    double cross_section_ion, kappa_ion, dummy;
                    if(rt_ion_G_HI[k] > 0)
                    {
                        cross_section_ion = nH0 * rt_ion_sigma_HI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HI[k] * cross_section_ion * c_nH_time_n_photons_vol;// (egy per photon x cross section x photon flux) :: attenuation factors [already in flux/energy update]: * slab_averaging_function(kappa_ion * Sigma_particle); // egy per photon x cross section x photon flux (w attenuation factors) // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);  // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                    }
                    if(rt_ion_G_HeI[k] > 0)
                    {
                        cross_section_ion = nHe0 * rt_ion_sigma_HeI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HeI[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion * Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);  // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                    }
                    if(rt_ion_G_HeII[k] > 0)
                    {
                        cross_section_ion = nHep * rt_ion_sigma_HeII[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HeII[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion*Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt); // commented-out terms not appropriate here based on how we treat RSOL terms
                        Heat += dummy;
                    }
                }
            }
        }
#endif

        Heat += CR_gas_heating(target, n_elec, nH0, nHcgs); // CR hadronic+Coulomb+ionization heating //
#if defined(COOL_LOW_TEMPERATURES)
        if(LambdaMol<0) {Heat -= LambdaMol;} // Molecular line heating (Trad_mol_cooling_batch > Tgas) //
        if(LambdaMetal<0) {Heat -= LambdaMetal;} // potential net heating from low-temperature gas-phase metal line absorption //
#endif
        if(LambdaCompton<0) {Heat -= LambdaCompton;} /* Compton heating rather than cooling */

#if defined(GALSF_FB_FIRE_RT_UVHEATING) || defined(RT_PHOTOELECTRIC) || defined(RT_ISRF_BACKGROUND)
        /* Photoelectric heating following Bakes & Thielens 1994 (also Wolfire 1995); now with 'update' from Wolfire 2005 for PAH [fudge factor 0.5 below] */
        if((target >= 0) && (T < 1.0e6))
        {
            double photoelec = 0;
#ifdef GALSF_FB_FIRE_RT_UVHEATING
            photoelec += SphP[target].Rad_Flux_UV;
#ifdef COOL_UVB_SELFSHIELD_RAHMATI
            if(gJH0>0 && shieldfac>0) {photoelec += sqrt(shieldfac) * (gJH0 / 2.29e-10);} // uvb contribution //
#endif
#endif
#ifdef RT_PHOTOELECTRIC
            photoelec += SphP[target].Rad_E_gamma[RT_FREQ_BIN_PHOTOELECTRIC] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // convert to Habing field //
#endif
#if defined(RT_ISRF_BACKGROUND) && (!defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)) // latter flag decides whether we do treecol/sobolev here to get the background intensity // add a constant assumed FUV background, for isolated ISM simulations that don't get FUV from local sources self-consistently
            double column = evaluate_NH_from_GradRho(P[target].GradRho,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target) * UNIT_SURFDEN_IN_CGS; // converts to cgs            
            photoelec += All.InterstellarRadiationFieldStrength * 1.7 * exp(-DMAX(P[target].Metallicity[0]/All.SolarAbundances[0],1e-4) * column * 500.); // RT_ISRF_BACKGROUND rescales the overal ISRF, factor of 1.7 gives Draine 1978 field in Habing units, extinction factor assumes the same FUV band-integrated dust opacity as RT module
#endif
            if(photoelec > 0) {if(photoelec > 1.e4) {photoelec = 1.e4;}}

            if(photoelec > 0)
            {
                LambdaPElec = -1.3e-24 * photoelec / nHcgs * (P[target].Metallicity[0]/All.SolarAbundances[0]) * return_dust_to_metals_ratio_vs_solar(target,0); // negative sign for lambda b/c heating
                double x_photoelec = photoelec * sqrt(T) / (0.5 * (1.0e-12+n_elec) * nHcgs);
                LambdaPElec *= 0.049/(1+pow(x_photoelec/1925.,0.73)) + 0.037*pow(T/1.0e4,0.7)/(1+x_photoelec/5000.);
                Heat -= LambdaPElec;
            }
        }
#endif
    }
  else				/* here we're outside of tabulated rates, T>Tmax K */
    {
        /* at high T (fully ionized); only free-free and Compton cooling are present.  Assumes no heating. */
        Heat = LambdaExc = LambdaExcH0 = LambdaExcHep = LambdaIon = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRec = LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;
        nHp = 1.0; nHep = 0; nHepp = yhelium(target); n_elec = nHp + 2.0 * nHepp; /* very hot: H and He both fully ionized */
        *n_elec_eval = n_elec; /* save this value for the output cycle */

        LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * n_elec; // free-free
        LambdaCompton = evaluate_Compton_heating_cooling_rate(target,T,nHcgs,n_elec,shieldfac); // Compton
        Lambda = LambdaFF + DMAX(LambdaCompton,0);
    }

#if defined(RT_INFRARED)
    if(target >= 0) { /* attempt to account for gas-phase absorption self-opacity to cooling radiation here in limit where optically-thick regions are poorly resolved: only valid in some limits, so not really general here */
        double gas_self_absorption_opacity = rt_kappa_adaptive_IR_band(target,T,T,-1,-1), surface_density_fromcenter = 0.5 * (SphP[target].Density*All.cf_a3inv) * (Get_Particle_Size(target)*All.cf_atime);
        double tau_self = gas_self_absorption_opacity * surface_density_fromcenter, fcorr = 1./(1.+tau_self*tau_self);
        Heat*=fcorr; Lambda*=fcorr; LambdaMetal*=fcorr; LambdaExc*=fcorr; LambdaRec*=fcorr; LambdaIon*=fcorr; LambdaPElec*=fcorr; LambdaFF*=fcorr; LambdaMol*=fcorr; LambdaDust*=fcorr; LambdaCompton*=fcorr;
    }
#endif


#if defined(RT_NUV)
    double Lambda_rad_NUV = LambdaMetal; // most of LambdaMetal coming out in the NUV, as we define it
    Lambda_rad_NUV += LambdaExc + LambdaIon; // this represents gas kinetic energy lost to collisional ionization and excitation. but each is assumed to produce a recombination or cascade back to the ground state, which should re-emit. we're usually assuming case B recombination (UV emitted photons re-absorbed), so we'll assume a cascade for these into NUV/optical and other bands [otherwise should be added to photo-ionizing band]. still ignore LambdaRec, because otherwise this will double-count the UV background
#if !defined(RT_PHOTOELECTRIC) // if this module is active, these photons are accounted for explicitly in the photoelectric bands
    Lambda_rad_NUV += LambdaPElec; // otherwise lump it in here as well since it overlaps this band (should deplete it appropriately)
#endif
#if !defined(RT_FREEFREE) // if this module is active, these photons are accounted for explicitly in the free-free bands
    if(logT >= 5.) {Lambda_rad_NUV += LambdaFF;} // this can be coming at a wide range of wavelengths of course, but if not explicitly followed, it will go in one bin or another depending on the gas temperature
#endif
    if(target>=0) {
        int which_bin_to_cool_to = RT_FREQ_BIN_NUV; // default to NUV bin
#if defined(RT_INFRARED)
        SphP[target].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_NUV] = SphP[target].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_INFRARED] = 0; // set to nil before deciding if we will add radiation here, to avoid double-counting with previous-loop information
        if(SphP[target].Radiation_Temperature > 1.e4) {which_bin_to_cool_to = RT_FREQ_BIN_INFRARED;} // our more-accurate effective/adaptive IR band is already covering these wavelengths, rather than do the noisy step of cooling to NUV, re-absorbing with less accurate opacities and down-grading to the grey-band, just dump directly to the grey-band
#endif
        SphP[target].Lambda_RadiativeCooling_toRHDBins[which_bin_to_cool_to] += Lambda_rad_NUV;} // save this to be used later (include all misc terms that will appear in our NUV radiation umbrella)
#endif
    
#if defined(RT_INFRARED)
    double Lambda_rad_IR = LambdaMol; // include molecular line cooling, cold atomic cooling (in LambdaMol) and dust cooling here, since all coming out in the IR
    Lambda_rad_IR += LambdaCompton; // Compton cooling/heating influences things here as well. Need to be a bit careful because Compton cooling can, in very low-density regions, be off the CMB, which we aren't explicitly evolving. so may want to add a check here for this case. but in the limits of interest that term is small, so not so important here.
#if !defined(RT_FREEFREE) // if this module is active, these photons are accounted for explicitly in the free-free bands
    if(logT < 5.) {Lambda_rad_IR += LambdaFF;} // this can be coming at a wide range of wavelengths of course, but if not explicitly followed, it will go in one bin or another depending on the gas temperature
#endif
    if(target>=0) {SphP[target].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_INFRARED] += Lambda_rad_IR;} // save this to be used later (include all misc terms that will appear in our IR radiation umbrella)
#endif

#if defined(COOL_LOW_TEMPERATURES)
    /* now add the dust cooling/heating terms - for RT+cooling runs this must be performed after we know all other radiation cooling terms, 
     because we are performing a backward-Euler solve and need to determine the consistent IR radiation energy and Tdust at t+dt, solving
     the equations of gas+dust+radiation energy conservation with 0 dust heat capacity - note that Tdust is always in equilibrium in
     this system, and we should overwrite it with whatever we get here. */
#ifdef RT_INFRARED
    if(target >= 0) {LambdaDust = rt_ir_lambdadust(target, T);} /* This updates Dust_Temperature and Lambda_RadiativeCooling_toRHDBins */
#endif
    if(LambdaDust>0) {Lambda += LambdaDust;} /* add the -positive- Lambda-dust associated with cooling */
    if(LambdaDust<0) {Heat -= LambdaDust;} // Dust collisional heating (Tdust > Tgas) //
#endif
#ifdef RT_COOLING_DUST_ONLY
    if(LambdaDust > 0){Lambda = LambdaDust; Heat = 0;}
    if(LambdaDust <= 0){Heat = -LambdaDust; Lambda = 0;}
#endif
    double Q = Heat - Lambda;
#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target>=0) {SphP[target].CoolingRate = Lambda; SphP[target].HeatingRate = Heat;}
#endif
    
#if defined(COOL_LOW_TEMPERATURES) && !defined(COOL_LOWTEMP_THIN_ONLY)
    /* if we are in the optically thick limit, we need to modify the cooling/heating rates according to the appropriate limits;
        this flag does so by using a simple approximation. we consider the element as if it were a slab, with a column density
        calculated from the simulation properties and the Sobolev approximation. we then assume it develops an equilibrium internal
        temperature structure on a radiative diffusion timescale much faster than the dynamical time, and so the surface radiation
        from a photosphere can be simply related to the local density by the optical depth to infinity. the equations here follow
        Rafikov, 2007 (ApJ, 662, 642):
            denergy/dt/dArea = sigma*T^4 / fc(tau)
            fc(tau) = tau^eta + 1/tau (taking chi, phi~1; the second term describes the optically thin limit, which is calculated above
                more accurately anyways - that was just Kirchoff's Law; so we only need to worry about the first term)
            eta = 4*(gamma-1) / [gamma*(1+alpha+beta*(gamma-1)/gamma)], where gamma=real polytropic index, and alpha/beta follow
                an opacity law kappa=kappa_0 * P^alpha * T^beta. for almost all the regimes of interest, however, eta~1, which is also
                what is obtained for a convectively stable slab. so we will use this.
            now, this gives sigma*T^4/tau * Area_eff / nHcgs as the 'effective' cooling rate in our units of Heat or Lambda above.
                the nHcgs just puts it in the same volumetric terms. The Area_eff must be defined as ~m_particle/surface_density
                to have the same meaning for a slab as assumed in Rafikov (and to integrate correctly over all particles in the slab,
                if/when the slab is resolved). We estimate this in our usual fashion with the Sobolev-type column density
            tau = kappa * surface_density; we estimate kappa ~ 5 cm^2/g * (0.001+Z/Z_solar), as the frequency-integrated kappa for warm
                dust radiation (~150K), weighted by the dust-to-gas ratio (with a floor for molecular absorption). we could make this
                temperature-dependent, though, fairly easily - for this particular problem it won't make much difference
        This rate then acts as an upper limit to the net heating/cooling calculated above (restricts absolute value)
     */
    if( (nHcgs > 0.1) && (target >= 0) )  /* don't bother at very low densities, since youre not optically thick, and protect from target=-1 with GALSF_EFFECTIVE_EQS */
    {
        double surface_density = evaluate_NH_from_GradRho(SphP[target].Gradients.Density,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target);
        surface_density *=  UNIT_SURFDEN_IN_CGS; // converts to cgs
        double effective_area = 2.3 * PROTONMASS_CGS / surface_density; // since cooling rate is ultimately per-particle, need a particle-weight here
        double kappa_eff; // effective kappa, accounting for metal abundance, temperature, and density //
	
	    kappa_eff = rt_kappa_adaptive_IR_band(target,T,T,0,1) / UNIT_SURFDEN_IN_CGS; // will return simple opacity law kappa = 0.1cm^2/g (T/10K)^2, capped at 5 cm^2/g, in code units [convert to physical here]
        if(kappa_eff < 0.1) {kappa_eff=0.1;}
        if(T>1500.){
            /* this is an approximate result for high-temperature opacities, but provides a pretty good fit from 1.5e3 - 1.0e9 K */
            double k_electron = 0.2 * (1. + HYDROGEN_MASSFRAC); //0.167 * n_elec; /* Thompson scattering (non-relativistic) */
            double k_molecular = 0.1 * P[target].Metallicity[0]; /* molecular line opacities */

            double k_Hminus = 1.1e-25 * sqrt(P[target].Metallicity[0] * rho) * pow(T,7.7); /* negative H- ion opacity */
            double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (P[target].Metallicity[0]+0.001) * rho / (T*T*T*sqrt(T)); /* free-free, bound-free, bound-bound transitions */
            double k_radiative = k_molecular + 1./(1./k_Hminus + 1./(k_electron+k_Kramers)); /* approximate interpolation between the above opacities */
            double k_conductive = 2.6e-7 * n_elec * T*T/(rho*rho); //*(1+pow(rho/1.e6,0.67) /* e- thermal conductivity can dominate at low-T, high-rho, here it as expressed as opacity */
            kappa_eff += 1./(1./k_radiative + 1./k_conductive); /* effective opacity including both heat carriers (this is exact) */
        }
        double tau_eff = kappa_eff * surface_density;
        double Lambda_Thick_BlackBody = 5.67e-5 * (T*T*T*T) * effective_area / ((1.+tau_eff) * nHcgs);
        if(Q > 0) {if(Q > Lambda_Thick_BlackBody) {Q=Lambda_Thick_BlackBody;}} else {if(Q < -Lambda_Thick_BlackBody) {Q=-Lambda_Thick_BlackBody;}}
    }
#endif

#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target>=0) {SphP[target].NetHeatingRateQ = Q;}
#endif
#ifdef OUTPUT_MOLECULAR_FRACTION
    if(target>0) {SphP[target].MolecularMassFraction = Get_Gas_Molecular_Mass_Fraction(target, T, nH0, n_elec, sqrt(shieldfac)*(gJH0/2.29e-10));}
#endif

#ifndef COOLING_OPERATOR_SPLIT
    /* add the hydro energy change directly: this represents an additional heating/cooling term, to be accounted for in the semi-implicit solution determined here. this is more accurate when tcool << tdynamical */
    if(target >= 0) {if(SphP[target].CoolingIsOperatorSplitThisTimestep==0) {Q += SphP[target].DtInternalEnergy / nHcgs;}}
#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target >= 0) {SphP[target].HydroHeatingRate = SphP[target].DtInternalEnergy / nHcgs;}
#endif
#endif
    return Q;
} // ends CoolingRate



void InitCoolMemory(void)
{
    BetaH0 = (double *) mymalloc("BetaH0", (NCOOLTAB + 1) * sizeof(double));
    BetaHep = (double *) mymalloc("BetaHep", (NCOOLTAB + 1) * sizeof(double));
    AlphaHp = (double *) mymalloc("AlphaHp", (NCOOLTAB + 1) * sizeof(double));
    AlphaHep = (double *) mymalloc("AlphaHep", (NCOOLTAB + 1) * sizeof(double));
    Alphad = (double *) mymalloc("Alphad", (NCOOLTAB + 1) * sizeof(double));
    AlphaHepp = (double *) mymalloc("AlphaHepp", (NCOOLTAB + 1) * sizeof(double));
    GammaeH0 = (double *) mymalloc("GammaeH0", (NCOOLTAB + 1) * sizeof(double));
    GammaeHe0 = (double *) mymalloc("GammaeHe0", (NCOOLTAB + 1) * sizeof(double));
    GammaeHep = (double *) mymalloc("GammaeHep", (NCOOLTAB + 1) * sizeof(double));
    Betaff = (double *) mymalloc("Betaff", (NCOOLTAB + 1) * sizeof(double));

#ifdef COOL_METAL_LINES_BY_SPECIES
    long i_nH=41; long i_T=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES;
    SpCoolTable0 = (float *) mymalloc("SpCoolTable0",(kspecies*i_nH*i_T)*sizeof(float));
    if(All.ComovingIntegrationOn) {SpCoolTable1 = (float *) mymalloc("SpCoolTable1",(kspecies*i_nH*i_T)*sizeof(float));}
#endif
}



void MakeCoolingTable(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19
        Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
    int i; double T,Tfact;
    if(All.MinGasTemp > 0.0) {Tmin = log10(All.MinGasTemp);} else {Tmin=-1.0;} // set minimum temperature in this table to some very low value if zero, where none of the cooling approximations above make sense
    deltaT = (Tmax - Tmin) / NCOOLTAB;
    /* minimum internal energy for neutral gas */
    for(i = 0; i <= NCOOLTAB; i++)
    {
        BetaH0[i] = BetaHep[i] = Betaff[i] = AlphaHp[i] = AlphaHep[i] = AlphaHepp[i] = Alphad[i] = GammaeH0[i] = GammaeHe0[i] = GammaeHep[i] = 0;
        T = pow(10.0, Tmin + deltaT * i);
        Tfact = 1.0 / (1 + sqrt(T / 1.0e5));
        if(118348. / T < 70.) {BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;}
        if(473638. / T < 70.) {BetaHep[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;}

        Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));
        //AlphaHp[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);	/* old Cen92 fit */
        //AlphaHep[i] = 1.5e-10 * pow(T, -0.6353); /* old Cen92 fit */
        //AlphaHepp[i] = 4. * AlphaHp[i];	/* old Cen92 fit */
        AlphaHp[i] = 7.982e-11 / ( sqrt(T/3.148) * pow((1.0+sqrt(T/3.148)), 0.252) * pow((1.0+sqrt(T/7.036e5)), 1.748) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHep[i]= 9.356e-10 / ( sqrt(T/4.266e-2) * pow((1.0+sqrt(T/4.266e-2)), 0.2108) * pow((1.0+sqrt(T/3.676e7)), 1.7892) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHepp[i] = 2. * 7.982e-11 / ( sqrt(T/(4.*3.148)) * pow((1.0+sqrt(T/(4.*3.148))), 0.252) * pow((1.0+sqrt(T/(4.*7.036e5))), 1.748) ); /* Verner & Ferland (1996) : ~ Z*alphaHp[1,T/Z^2] */

        if(470000.0 / T < 70) {Alphad[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));}
        if(157809.1 / T < 70) {GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;}
        if(285335.4 / T < 70) {GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;}
        if(631515.0 / T < 70) {GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;}
    }
}


#ifdef COOL_METAL_LINES_BY_SPECIES

void LoadMultiSpeciesTables(void)
{
    if(All.ComovingIntegrationOn) {
        int i;
        double z;
        if(All.Time==All.TimeBegin) {
            All.SpeciesTableInUse=48;
            ReadMultiSpeciesTables(All.SpeciesTableInUse);
        }
        z=log10(1/All.Time)*48;
        i=(int)z;
        if(i<48) {
            if(i<All.SpeciesTableInUse) {
                All.SpeciesTableInUse=i;
                ReadMultiSpeciesTables(All.SpeciesTableInUse);
            }}
    } else {
        if(All.Time==All.TimeBegin) ReadMultiSpeciesTables(0);
    }
}

void ReadMultiSpeciesTables(int iT)
{
    /* read table w n,T for each species */
    long i_nH=41; long i_Temp=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES; long i,j,k,r;
    /* int i_He=7;  int l; */
    FILE *fdcool; char *fname;

    fname=GetMultiSpeciesFilename(iT,0);
    if(ThisTask == 0) printf(" ..opening Cooling Table %s \n",fname);
    if(!(fdcool = fopen(fname, "r"))) {
        printf(" Cannot read species cooling table in file `%s'\n", fname); endrun(456);}
    for(i=0;i<kspecies;i++) {
        for(j=0;j<i_nH;j++) {
            for(k=0;k<i_Temp;k++) {
                r=fread(&SpCoolTable0[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                if(r!=1) {printf(" Reached Cooling EOF! \n");
                }
            }}}
    fclose(fdcool);
    /*
     GetMultiSpeciesFilename(iT,&fname,1);
     if(!(fdcool = fopen(fname, "r"))) {
     printf(" Cannot read species (He) cooling table in file `%s'\n", fname); endrun(456);}
     for(i=0;i<2;i++)
     for(j=0;j<i_nH;j++)
     for(k=0;k<i_Temp;k++)
     for(l=0;l<i_He;l++)
     fread(&SpCoolTable0_He[i][j][k][l],sizeof(float),1,fdcool);
     fclose(fdcool);
     */
    if (All.ComovingIntegrationOn && i<48) {
        fname=GetMultiSpeciesFilename(iT+1,0);
        if(ThisTask == 0) printf(" ..opening (z+) Cooling Table %s \n",fname);
        if(!(fdcool = fopen(fname, "r"))) {
            printf(" Cannot read species 1 cooling table in file `%s'\n", fname); endrun(456);}
        for(i=0;i<kspecies;i++) {
            for(j=0;j<i_nH;j++) {
                for(k=0;k<i_Temp;k++) {
                    r=fread(&SpCoolTable1[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                    if(r!=1) {printf(" Reached Cooling EOF! \n");
                    }
                }}}
        fclose(fdcool);
        /*
         GetMultiSpeciesFilename(iT+1,&fname,1);
         if(!(fdcool = fopen(fname, "r"))) {
         printf(" Cannot read species 1 (He) cooling table in file `%s'\n", fname); endrun(456);}
         for(i=0;i<2;i++)
         for(j=0;j<i_nH;j++)
         for(k=0;k<i_Temp;k++)
         for(l=0;l<i_He;l++)
         fread(&SpCoolTable1_He[i][j][k][l],sizeof(float),1,fdcool);
         fclose(fdcool);
         */
    }
}

char *GetMultiSpeciesFilename(int i, int hk)
{
    static char fname[100];
    if(i<0) i=0; if(i>48) i=48;
    if(hk==0) {
        sprintf(fname,"./spcool_tables/spcool_%d",i);
    } else {
        sprintf(fname,"./spcool_tables/spcool_He_%d",i);
    }
    return fname;
}

#endif



/* table input (from file TREECOOL) for ionizing parameters */
#define JAMPL	1.0		/* amplitude factor relative to input table */
#define TABLESIZE 250		/* Max # of lines in TREECOOL -- needs to be at least one larger than number of non-zero lines */
static float inlogz[TABLESIZE];
static double gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE]; // upgrade from float to double, should read fine
static double eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE]; // upgrade from float to double, should read fine
static int nheattab;		/* length of table */


void ReadIonizeParams(char *fname)
{
    int i; FILE *fdcool;
    if(!(fdcool = fopen(fname, "r"))) {printf(" Cannot read ionization table in file `%s'. Make sure the correct TREECOOL file is placed in the code run-time directory, and that any leading comments (e.g. lines preceded by ##) are deleted from the file.\n", fname); endrun(456);}
    for(i=0; i<TABLESIZE; i++) {inlogz[i]=100; gH0[i]=0; gHe[i]=0; gHep[i]=0; eH0[i]=0; eHe[i]=0; eHep[i]=0;}
    for(i=0; i<TABLESIZE; i++) {if(fscanf(fdcool, "%g %lg %lg %lg %lg %lg %lg", &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF) {break;}}
    fclose(fdcool);
    for(i=0, nheattab=0; i<TABLESIZE; i++) {if(gH0[i] != 0.0) {nheattab++;} else {break;}} /*  nheattab is the number of entries in the table */
    if(ThisTask == 0) printf(" ..read ionization table [TREECOOL] with %d non-zero UVB entries in file `%s'. Make sure to cite the authors from which the UV background was compiled! (See user guide for the correct references).\n", nheattab, fname);
}


void IonizeParams(void)
{
    IonizeParamsTable();
}



void IonizeParamsTable(void)
{
    int i, ilow, ihi;
    double logz, dzlow, dzhi;
    double redshift;

    if(All.ComovingIntegrationOn)
        {redshift = 1 / All.Time - 1;}
    else
    {
        /* in non-cosmological mode, still use, but adopt z=0 background */
#ifdef RT_ISRF_BACKGROUND
	redshift = All.RadiationBackgroundRedshift;
#else
        redshift = 0;
#endif
        /*
         gJHe0 = gJHep = gJH0 = epsHe0 = epsHep = epsH0 = J_UV = 0;
         return;
         */
    }

    logz = log10(redshift + 1.0);
    ilow = 0;
    if(nheattab <= 0) {gJHe0 = gJHep = gJH0 = epsHe0 = epsHep = epsH0 = J_UV = 0; return;}
    for(i=0; i<nheattab; i++) {if(inlogz[i] < logz) {ilow = i;} else {break;}}
    ihi = i + 1;
    if(ilow >= nheattab) {ihi = i;}
    dzlow = logz - inlogz[ilow];
    dzhi = inlogz[ihi] - logz;
    if((logz > inlogz[nheattab - 1]) || (gH0[ilow] == 0) || (gH0[ihi] == 0) || (ilow > nheattab))
    {
        gJHe0 = gJHep = gJH0 = epsHe0 = epsHep = epsH0 = J_UV = 0;
        return;
    }
    else {J_UV = 1.e-21;}		/* irrelevant as long as it's not 0 */

    gJH0 = JAMPL * pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ihi])) / (dzlow + dzhi));
    gJHe0 = JAMPL * pow(10., (dzhi * log10(gHe[ilow]) + dzlow * log10(gHe[ihi])) / (dzlow + dzhi));
    gJHep = JAMPL * pow(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ihi])) / (dzlow + dzhi));
    epsH0 = JAMPL * pow(10., (dzhi * log10(eH0[ilow]) + dzlow * log10(eH0[ihi])) / (dzlow + dzhi));
    epsHe0 = JAMPL * pow(10., (dzhi * log10(eHe[ilow]) + dzlow * log10(eHe[ihi])) / (dzlow + dzhi));
    epsHep = JAMPL * pow(10., (dzhi * log10(eHep[ilow]) + dzlow * log10(eHep[ihi])) / (dzlow + dzhi));

    return;
}


void SetZeroIonization(void)
{
    gJHe0 = gJHep = gJH0 = 0; epsHe0 = epsHep = epsH0 = 0; J_UV = 0;
}


void IonizeParamsFunction(void)
{
    int i, nint;
    double a0, planck, ev, e0_H, e0_He, e0_Hep;
    double gint, eint, t, tinv, fac, eps;
    double at, beta, s;
    double pi;

#define UVALPHA         1.0
    double Jold = -1.0;
    double redshift;

    J_UV = 0.;
    gJHe0 = gJHep = gJH0 = 0.;
    epsHe0 = epsHep = epsH0 = 0.;


    if(All.ComovingIntegrationOn)	/* analytically compute params from power law J_nu */
    {
        redshift = 1 / All.Time - 1;

        if(redshift >= 6) {J_UV = 0.;}
        else
        {
            if(redshift >= 3) {J_UV = 4e-22 / (1 + redshift);}
            else
            {
                if(redshift >= 2) {J_UV = 1e-22;}
                else {J_UV = 1.e-22 * pow(3.0 / (1 + redshift), -3.0);}
            }
        }
        if(J_UV == Jold) {return;}
        Jold = J_UV;
        if(J_UV == 0) {return;}


        a0 = 6.30e-18;
        planck = 6.6262e-27;
        ev = 1.6022e-12;
        e0_H = 13.6058 * ev;
        e0_He = 24.59 * ev;
        e0_Hep = 54.4232 * ev;

        gint = 0.0;
        eint = 0.0;
        nint = 5000;
        at = 1. / ((double) nint);

        for(i = 1; i <= nint; i++)
        {
            t = (double) i;
            t = (t - 0.5) * at;
            tinv = 1. / t;
            eps = sqrt(tinv - 1.);
            fac = exp(4. - 4. * atan(eps) / eps) / (1. - exp(-2. * M_PI / eps)) * pow(t, UVALPHA + 3.);
            gint += fac * at;
            eint += fac * (tinv - 1.) * at;
        }

        gJH0 = a0 * gint / planck;
        epsH0 = a0 * eint * (e0_H / planck);
        gJHep = gJH0 * pow(e0_H / e0_Hep, UVALPHA) / 4.0;
        epsHep = epsH0 * pow((e0_H / e0_Hep), UVALPHA - 1.) / 4.0;

        at = 7.83e-18;
        beta = 1.66;
        s = 2.05;

        gJHe0 = (at / planck) * pow((e0_H / e0_He), UVALPHA) *
        (beta / (UVALPHA + s) + (1. - beta) / (UVALPHA + s + 1));
        epsHe0 = (e0_He / planck) * at * pow(e0_H / e0_He, UVALPHA) *
        (beta / (UVALPHA + s - 1) + (1 - 2 * beta) / (UVALPHA + s) - (1 - beta) / (UVALPHA + s + 1));

        pi = M_PI;
        gJH0 *= 4. * pi * J_UV;
        gJHep *= 4. * pi * J_UV;
        gJHe0 *= 4. * pi * J_UV;
        epsH0 *= 4. * pi * J_UV;
        epsHep *= 4. * pi * J_UV;
        epsHe0 *= 4. * pi * J_UV;
    }
}
#endif // !(CHIMES)




void InitCool(void)
{
    if(ThisTask == 0) printf("Initializing cooling ...\n");

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();

#ifdef COOL_GRACKLE
    InitGrackle();
#endif

#ifdef CHIMES
    sprintf(ChimesGlobalVars.MainDataTablePath, "%s/chimes_main_data.hdf5", ChimesDataPath);
    sprintf(ChimesGlobalVars.EqAbundanceTablePath, "%s/EqAbundancesTables/%s", ChimesDataPath, ChimesEqAbundanceTable);
    sprintf(ChimesGlobalVars.PhotoIonTablePath[0], "%s/%s", ChimesDataPath, ChimesPhotoIonTable);

    // By default, use 1 spectrum, unless
    // stellar fluxes are enabled.
    ChimesGlobalVars.N_spectra = 1;

#ifdef CHIMES_STELLAR_FLUXES
    ChimesGlobalVars.N_spectra += CHIMES_LOCAL_UV_NBINS;

    int spectrum_idx;
    for (spectrum_idx = 1; spectrum_idx < ChimesGlobalVars.N_spectra; spectrum_idx++)
      sprintf(ChimesGlobalVars.PhotoIonTablePath[spectrum_idx], "%s/starburstCrossSections/cross_sections_SB%d.hdf5", ChimesDataPath, spectrum_idx);
#endif

    if (ChimesUVBMode > 0)
      {
	ChimesGlobalVars.redshift_dependent_UVB_index = 0;
	if (ChimesUVBMode == 1)
	  ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0;
	else
	  ChimesGlobalVars.use_redshift_dependent_eqm_tables = 1;
      }
    else
      {
	ChimesGlobalVars.redshift_dependent_UVB_index = -1;
	ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0;
      }

    // Hybrid cooling has not yet been
    // implemented in Gizmo. Switch
    // it off for now.
    ChimesGlobalVars.hybrid_cooling_mode = 0;

    // Set the chimes_exit() function
    // to use the Gizmo-specific
    // chimes_gizmo_exit().
    chimes_exit = &chimes_gizmo_exit;

    // Initialise the CHIMES module.
    init_chimes(&ChimesGlobalVars);

#ifdef CHIMES_METAL_DEPLETION
    chimes_init_depletion_data();
#endif

#else // CHIMES
    InitCoolMemory();
    MakeCoolingTable();
    ReadIonizeParams("TREECOOL");
    IonizeParams();
#ifdef COOL_METAL_LINES_BY_SPECIES
    LoadMultiSpeciesTables();
#endif
#endif // CHIMES
}



#ifndef CHIMES
#ifdef COOL_METAL_LINES_BY_SPECIES
double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z)
{
    double ne_over_nh_tbl=1, Lambda=0;
    int k, N_species_active = (int)NUM_LIVE_SPECIES_FOR_COOLTABLES;

    /* pre-calculate the indices for density and temperature, then we just need to call the tables by species */
    int ixmax=40, iymax=175;
    int ix0, iy0, ix1, iy1;
    double dx, dy, dz, mdz;
    long i_T=iymax+1, inHT=i_T*(ixmax+1);
    if(All.ComovingIntegrationOn && All.SpeciesTableInUse<48) {dz=log10(1/All.Time)*48; dz=dz-(int)dz; mdz=1-dz;} else {dz=0; mdz=1;}

    dx = (log10(nHcgs)-(-8.0))/(0.0-(-8.0))*ixmax;
    dy = (logT-2.0)/(9.0-2.0)*iymax;
    if(dx<0) {dx=0;} else {if(dx>ixmax) {dx=ixmax;}}
    ix0=(int)dx; ix1=ix0+1; if(ix1>ixmax) {ix1=ixmax;}
    dx=dx-ix0;
    if(dy<0) {dy=0;} else {if(dy>iymax) {dy=iymax;}}
    iy0=(int)dy; iy1=iy0+1; if(iy1>iymax) {iy1=iymax;}
    dy=dy-iy0;
    long index_x0y0=iy0+ix0*i_T, index_x0y1=iy1+ix0*i_T, index_x1y0=iy0+ix1*i_T, index_x1y1=iy1+ix1*i_T;

    ne_over_nh_tbl = GetLambdaSpecies(0,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz);
    if(ne_over_nh_tbl > 0)
    {
        double zfac = 0.0127 / All.SolarAbundances[0];
        for (k=1; k<N_species_active; k++)
        {
            long k_index = k * inHT;
            Lambda += GetLambdaSpecies(k_index,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz) * Z[k+1]/(All.SolarAbundances[k+1]*zfac);
        }
        Lambda /= ne_over_nh_tbl;
    }
    return Lambda;
}


double GetLambdaSpecies(long k_index, long index_x0y0, long index_x0y1, long index_x1y0, long index_x1y1, double dx, double dy, double dz, double mdz)
{
    long x0y0 = index_x0y0 + k_index;
    long x0y1 = index_x0y1 + k_index;
    long x1y0 = index_x1y0 + k_index;
    long x1y1 = index_x1y1 + k_index;
    double i1, i2, j1, j2, w1, w2, u1;
    i1 = SpCoolTable0[x0y0];
    i2 = SpCoolTable0[x0y1];
    j1 = SpCoolTable0[x1y0];
    j2 = SpCoolTable0[x1y1];
    if(dz > 0)
    {
        i1 = mdz * i1 + dz * SpCoolTable1[x0y0];
        i2 = mdz * i2 + dz * SpCoolTable1[x0y1];
        j1 = mdz * j1 + dz * SpCoolTable1[x1y0];
        j2 = mdz * j2 + dz * SpCoolTable1[x1y1];
    }
    w1 = i1*(1-dy) + i2*dy;
    w2 = j1*(1-dy) + j2*dy;
    u1 = w1*(1-dx) + w2*dx;
    return u1;
}

#endif // COOL_METAL_LINES_BY_SPECIES
#endif // !(CHIMES)


#ifdef GALSF_FB_FIRE_RT_UVHEATING
void selfshield_local_incident_uv_flux(void)
{   /* include local self-shielding with the following */
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type==0)
        {
            if((SphP[i].Rad_Flux_UV>0) && (PPP[i].Hsml>0) && (SphP[i].Density>0) && (P[i].Mass>0) && (All.Time>0))
            {
                SphP[i].Rad_Flux_UV *= UNIT_FLUX_IN_CGS * 1276.19; SphP[i].Rad_Flux_EUV *= UNIT_FLUX_IN_CGS * 1276.19; // convert to Habing units (normalize strength to local MW field in this [narrow] band, so not the 'full' Habing flux)
                double surfdensity = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i); // in CGS
                double tau_nuv = rt_kappa(i,RT_FREQ_BIN_FIRE_UV) * surfdensity; // optical depth: this part is attenuated by dust //
                double tau_euv = 3.7e6 * surfdensity * UNIT_SURFDEN_IN_CGS; // optical depth: 912 angstrom kappa_euv: opacity from neutral gas //
                SphP[i].Rad_Flux_UV *= exp(-DMIN(tau_nuv,90.)); // attenuate [important in newer modules depending on UV flux to fully-attenuate down to << 1e-6 in dense gas]
                //SphP[i].Rad_Flux_UV *= 0.01 + 0.99/(1.0 + 0.8*tau_nuv + 0.85*tau_nuv*tau_nuv); // attenuate (for clumpy medium with hard-minimum 1% scattering)
                //SphP[i].Rad_Flux_EUV *= exp(-DMIN(tau_euv,90.)); // attenuate [important in newer modules depending on UV flux to fully-attenuate down to << 1e-6 in dense gas]
                SphP[i].Rad_Flux_EUV *= 0.01 + 0.99/(1.0 + 0.8*tau_euv + 0.85*tau_euv*tau_euv); // attenuate (for clumpy medium with hard-minimum 1% scattering) //
            } else {SphP[i].Rad_Flux_UV = SphP[i].Rad_Flux_EUV = 0;}
        }
    }
}
#endif




/* subroutine to update the molecular fraction using our implicit solver for a simple --single-species-- network (just H2) */
void update_explicit_molecular_fraction(int i, double dtime_cgs)
{
    if(dtime_cgs <= 0) {return;}
#ifdef COOL_MOLECFRAC_NONEQM
    // first define a number of environmental variables that are fixed over this update step
    double fH2_initial = SphP[i].MolecularMassFraction_perNeutralH; // initial molecular fraction per H atom, entering this subroutine, needed for update below
    double xn_e=1, nh0=0, nHe0, nHepp, nhp, nHep, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergy;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &xn_e, &nh0, &nhp, &nHe0, &nHep, &nHepp); // get thermodynamic properties [will assume fixed value of fH2 at previous update value]
    double T=1, f_dustgas_solar=1, urad_G0=1, xH0=1, x_e=0, nH_cgs=rho*UNIT_DENSITY_IN_NHCGS; // initialize definitions of some variables used below to prevent compiler warnings
    f_dustgas_solar=1; urad_G0=1; // initialize metal/dust and radiation fields. will assume solar-Z and spatially-uniform Habing field for incident FUV radiation unless reset below.
#ifdef RT_ISRF_BACKGROUND
    urad_G0 = All.InterstellarRadiationFieldStrength;
#endif
#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
    if(SphP[i].DelayTimeHII > 0) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} // force gas flagged as in HII regions to have zero molecular fraction
#endif
    if(temperature > 3.e5) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} else {T=temperature;} // approximations below not designed for high temperatures, should simply give null
    xH0 = DMIN(DMAX(nh0, 0.),1.); // get neutral fraction [given by call to this program]
    if(xH0 <= MIN_REAL_NUMBER) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} // no neutral gas, no molecules!
    x_e = DMIN(DMAX(xn_e, 0.),2.); // get free electron ratio [number per H nucleon]
    double log_T=log10(T), ln_T=log(T), gamma_12=return_local_gammamultiplier(i)*gJH0/1.0e-12, shieldfac=return_uvb_shieldfac(i,gamma_12,nH_cgs,log_T), urad_from_uvb_in_G0=sqrt(shieldfac)*(gJH0/2.29e-10); // estimate UVB contribution if we have partial shielding, to full photo-dissociation rates //
#ifdef METALS
    f_dustgas_solar=(P[i].Metallicity[0]/All.SolarAbundances[0])*return_dust_to_metals_ratio_vs_solar(i,0); // this is only used for the dust-phase formation rates below, so just the dust term here
#endif
    /* get incident radiation field from whatever module we are using to track it */
#ifdef GALSF_FB_FIRE_RT_UVHEATING
    urad_G0 = DMAX(SphP[i].Rad_Flux_UV, 1.e-10); // note this is ALREADY self-shielded by dust, so we need to be careful about 2x-counting the self-shielding approximation below; hence limit this to a rather sizeable value  //
#endif
#if defined(RT_PHOTOELECTRIC) || defined(RT_LYMAN_WERNER)
    int whichbin = RT_FREQ_BIN_LYMAN_WERNER;
#if !defined(RT_LYMAN_WERNER)
    whichbin = RT_FREQ_BIN_PHOTOELECTRIC; // use photo-electric bin as proxy (very close) if don't evolve LW explicitly
#endif
    urad_G0 = SphP[i].Rad_E_gamma[whichbin] * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // convert to Habing field //
#endif
    urad_G0 += urad_from_uvb_in_G0; // include whatever is contributed from the meta-galactic background, fed into this routine
    urad_G0 = DMIN(DMAX( urad_G0 , 1.e-10 ) , 1.e10 ); // limit values, because otherwise exponential self-shielding approximation easily artificially gives 0 incident field
#ifdef RT_INFRARED
    urad_G0 += rt_irband_egydensity_in_band(i,11.2,500.) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // add contribution from the adaptive band
#endif
    // define a number of variables needed in the shielding module
    double dx_cell = Get_Particle_Size(i) * All.cf_atime; // cell size
    double surface_density_H2_0 = 5.e14 * PROTONMASS_CGS, x_exp_fac=0.00085, w0=0.2; // characteristic cgs column for -molecular line- self-shielding
    w0 = 0.035; // actual calibration from Drain, Gnedin, Richings, others: 0.2 is more appropriate as a re-calibration for sims doing local eqm without ability to resolve shielding at higher columns
    //double surface_density_local = xH0 * SphP[i].Density * All.cf_a3inv * dx_cell * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth through the local cell/slab. note G0 is -already- attenuated in the pre-processing by dust.
    double surface_density_local = xH0 * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i) * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth to infinity with our Sobolev-type approximation. Note G0 is already attenuated by dust, but we need to include H2 self-shielding, for which it is appropriate to integrate to infinity.
    double v_thermal_rms = 0.111*sqrt(T); // sqrt(3*kB*T/2*mp), since want rms thermal speed of -molecular H2- in kms
    double gradv = velocity_gradient_norm(i);
    double dv_turb=gradv*dx_cell*UNIT_VEL_IN_KMS; // delta-velocity across cell
    double x00 = surface_density_local / surface_density_H2_0, x01 = x00 / (sqrt(1. + 3.*dv_turb*dv_turb/(v_thermal_rms*v_thermal_rms)) * sqrt(2.)*v_thermal_rms), y_ss, x_ss_1, x_ss_sqrt, fH2_tmp, fH2_max, fH2_min, Q_max, Q_min, Q_initial; // variable needed below. note the x01 term corrects following Gnedin+Draine 2014 for the velocity gradient at the sonic scale, assuming a Burgers-type spectrum [their Eq. 3]
    double b_time_Mach = 0.5 * dv_turb / (v_thermal_rms/sqrt(3.)); // cs_thermal for molecular [=rms v_thermal / sqrt(3)], dv_turb to full inside dx, assume "b" prefactor for compressive-to-solenoidal ratio corresponding to the 'natural mix' = 0.5. could further multiply by 1.58 if really needed to by extended dvturb to 2h = H, and vthermal from molecular to atomic for the generating field, but not as well-justified
    double clumping_factor = 1. + b_time_Mach*b_time_Mach; // this is the exact clumping factor for a standard lognormal PDF with S=ln[1+b^2 Mach^2] //
    double clumping_factor_3 = clumping_factor*clumping_factor*clumping_factor; // clumping factor N for <rho^n>/<rho>^n = clumping factor^(N*(N-1)/2) //
    
    /* evolve dot[nH2]/nH0 = d_dt[fH2[neutral]] = (1/nH0) * (a_Z*rho_dust*nHI [dust formation] + a_GP*nHI*ne [gas-phase formation] + b_3B*nHI*nHI*(nHI+nH2/8) [3-body collisional form] - b_H2HI*nHI*nH2 [collisional dissociation]
        - b_H2H2*nH2*nH2 [collisional mol-mol dissociation] - Gamma_H2^LW * nH2 [photodissociation] - Gamma_H2^+ [photoionization] - xi_H2*nH2 [CR ionization/dissociation] ) */
    double fH2=0, sqrt_T=sqrt(T), nH0=xH0*nH_cgs, EXPmax=90.; int iter=0; // define some variables for below, including neutral H number density, free electron number, etc.
    double x_p = DMIN(DMAX(nhp , x_e/10.), 2.); // get free H+ fraction [cap because irrelevant to below in very low regime //
    /* use interpolation function from Glover & Abel 2008 [GA08], section 2.1.3, for interpolating between ground state (v=0) and LTE assumptions for states for collisional dissociation rates */
    double XH=HYDROGEN_MASSFRAC, xH2_guess=XH*DMAX(DMIN(SphP[i].MolecularMassFraction,1.),0.), xH_guess=DMAX(XH-xH2_guess,0), xHe_guess=nHe0+nHep+nHepp;
    double logT4=log_T-4., ncr_H=pow(10.,3.0-0.416*logT4-0.327*logT4*logT4), ncr_H2=pow(10.,4.845-1.3*logT4+1.62*logT4*logT4), ncr_He=pow(10.,5.0792*(1.-1.23e-5*(T-2000.)));
    double ncrit = 1./(xH_guess/ncr_H + xH2_guess/ncr_H2 + xHe_guess/ncr_He), n_ncrit=nH_cgs/ncrit, f_v0_LTE=1./(1.+n_ncrit), f_LTE_v0 = 1.-f_v0_LTE;

    double b_H2Hp = DMAX(0., -3.3232183e-7 + 3.3735382e-7*ln_T -1.4491368e-7*ln_T*ln_T + 3.4172805e-8*ln_T*ln_T*ln_T -4.7813720e-9*ln_T*ln_T*ln_T*ln_T +3.9731542e-10*ln_T*ln_T*ln_T*ln_T*ln_T -1.8171411e-11*ln_T*ln_T*ln_T*ln_T*ln_T*ln_T +3.5311932e-13*ln_T*ln_T*ln_T*ln_T*ln_T*ln_T*ln_T) * exp(-DMIN(21237.15/T,EXPmax)) * (nhp*nH_cgs) * clumping_factor; // H2-H+ dissociation, GA08-TableA1-7; note their expression (GA08) has an error where they write log[T] but this gives unphysical values. it should be ln[T], as it is correctly written in the original Savin et al. 2004 paper from which this fitting function is taken
    double b_H2e_v0 = 4.49e-9 * pow(T,0.11) * exp(-DMIN(101858./T,EXPmax)), b_H2e_LTE = 1.91e-9 * pow(T,0.136) * exp(-DMIN(53407.1/T,EXPmax)), b_H2e = pow(10., f_v0_LTE*log10(b_H2e_v0) + f_LTE_v0*log10(b_H2e_LTE)) * (x_e*nH_cgs) * clumping_factor; // collisional H2-e- dissociation; GA08-A1-8
    double b_H2HI_v0 = 6.67e-12 * sqrt_T * exp(-DMIN(1.+63593./T,EXPmax)), b_H2HI_LTE = 3.52e-9 * exp(-DMIN(43900./T,EXPmax)), b_H2HI = pow(10., f_v0_LTE*log10(b_H2HI_v0) + f_LTE_v0*log10(b_H2HI_LTE)) * (xH0*nH_cgs) * clumping_factor; // collisional H2-H dissociation; GA08-A1-10
    double b_H2H2_v0 = 5.996e-30 * pow(T,4.1881) * exp(-DMIN(54657.4/T,EXPmax)) / pow(1. + 6.761e-6*T , 5.6881), b_H2H2_LTE = 1.3e-9 * exp(-DMIN(53300./T,EXPmax)), b_H2H2 = pow(10., f_v0_LTE*log10(b_H2H2_v0) + f_LTE_v0*log10(b_H2H2_LTE)) * (xH0*nH_cgs/2.) * clumping_factor; // collisional H2-H2 dissociation; GA08-A1-10
    double b_H2He_v0_log = -27.029 + 3.801*log_T - 29487./T, b_H2He_LTE_log = -2.729 - 1.75*log_T - 23474./T, b_H2He = pow(10., f_v0_LTE*b_H2He_v0_log + f_LTE_v0*b_H2He_LTE_log) * (nHe0*nH_cgs) * clumping_factor; // collisional H2-He dissociation, GA08-A1-11
    double b_H2Hep = (3.7e-14*exp(DMIN(35./T,EXPmax)) + 7.2e-15) * ((nHep+nHepp)*nH_cgs) * clumping_factor; // collisional H2-He+ dissociation, GA08-A1-24+25
    // D questionable - this will really just convert to HD, should exclude here
    double b_H2D; if(T<=2000.) {b_H2D=pow(10.,-56.4737 +5.88886*log_T +7.19692*log_T*log_T +2.25069*log_T*log_T*log_T -2.16903*log_T*log_T*log_T*log_T +0.317887*log_T*log_T*log_T*log_T*log_T);} else {b_H2D=3.17e-10 * exp(-DMIN(5207./T,EXPmax));}
    b_H2D *= (xH0*2.527e-5*nH_cgs) * clumping_factor; // collisional H2-D dissociation, GA08-A1-37, using D abundance from Cooke, Pettini,& Steidel 2018
    double b_H2Dp = 1.e-9 * (DMAX(0.417 + 0.846*log_T - 0.137*log_T*log_T,  0.)) * (nhp*2.527e-5*nH_cgs) * clumping_factor; // collisional H2-D+ dissociation, GA08-A1-39, using D abundance from Cooke, Pettini,& Steidel 2018
    b_H2D*=1.e-10; b_H2Dp*=1.e-10; // see note above on D: include some non-zero here as a 'safety factor', but the overwhelming fraction is going to HD which we account for implicitly in cooling functions so dont explicitly solve here
    double b_H2ext = b_H2Hp + b_H2e + b_H2He + b_H2Hep + b_H2Dp; b_H2HI += b_H2D; b_H2ext*=1./2.; b_H2HI*=1./2.; b_H2H2*=1./4.; // collect dissociation terms where the secondary (e.g. e- does -not- scale with fmol as we define it here, and those where it does to different powers; 1/2 here is to account for nH2 = (1/2) * fH2 * nH because we will solve for fH2 as a mass fraction, becomes 1/4 in H2-H2 equation
    
    double Tdust = 30.; // need to assume something about dust temperature for reaction rates below for dust-phase formation
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) && (FLAG_NOT_IN_PUBLIC_CODE > 2)) || defined(SINGLE_STAR_SINK_DYNAMICS)
    Tdust = get_equilibrium_dust_temperature_estimate(i, shieldfac, T);
#endif
    double a_Z = 3.e-18*sqrt_T / ((1. +4.e-2*sqrt(T+Tdust) +2.e-3*T +8.e-6*T*T )*(1. +1.e4/exp(DMIN(EXPmax,600./Tdust)))) * f_dustgas_solar * nH0 * clumping_factor; // dust surface formation (assuming dust-to-metals ratio is 0.5*(Z/solar)*dust-to-gas-relative-to-solar in all regions where this is significant), from Glover & Jappsen 2007

    //double a_GP = (1.833e-18 * pow(T,0.88)) / (1. + x_p*1846.*(1.+T/20000.)/sqrt(T)) * xH0 * x_e*nH_cgs * clumping_factor; // gas-phase formation [Glover & Abel 2008, using fitting functions slightly more convenient and assuming H-->H2 much more rapid than other reactions, from Krumholz & McKee 2010; denominator factor accounts for p+H- -> H + H, instead of H2]; note the Nickerson version of this expression omits the ne,-3 term replacing it with ne from Krumholz+McKee which makes it incorrect by a factor of ~1000 in normalization
    double k1,k2,k5,k15,k16,k17,lnTeV=ln_T-9.35915,R51_n=3.62e-17/nH_cgs; // R51_n is contribution from photo-diss assuming ISRF-like since very low-E threshold (=0.755) photons, serves only as a 'floor' here
    if(T<=6000.) {k1=-17.845 + 0.762*log_T + 0.1523*log_T*log_T - 0.03274*log_T*log_T*log_T;} else {k1=-16.420 + 0.1998*log_T*log_T - 5.447e-3*log_T*log_T*log_T*log_T + 4.0415e-5*log_T*log_T*log_T*log_T*log_T*log_T;}
    k1=pow(10.,DMAX(k1,-50.));
    if(T<=300.) {k2=1.5e-9;} else {k2=4.0e-9 * pow(T,-0.17);}
    k5=5.7e-6/sqrt_T + 6.3e-8 - 9.2e-11*sqrt_T + 4.4e-13*T;
    k15=exp(DMAX(-EXPmax, -1.801849334e1 + 2.36085220e0*lnTeV -2.827443e-1*lnTeV*lnTeV +1.62331664e-2*lnTeV*lnTeV*lnTeV
                 -3.36501203e-2*lnTeV*lnTeV*lnTeV*lnTeV +1.17832978e-2*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                 -1.65619470e-3*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV +1.06827520e-4*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                 -2.63128581e-6*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV));
    if(T<1160.45) {k16=1.46629e-16 * pow(T,1.78186);} else {
        k16=exp(DMAX(-EXPmax, -2.0372609e1 + 1.13944933e0*lnTeV
                     -1.4210135e-1*lnTeV*lnTeV +8.4644554e-3*lnTeV*lnTeV*lnTeV
                     -1.4327641e-3*lnTeV*lnTeV*lnTeV*lnTeV +2.0122503e-4*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                     +8.6639632e-5*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV -2.5850097e-5*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV
                     +2.4555012e-6*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV -8.0683825e-8*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV*lnTeV));}
    if(T<=8000.) {k17=6.9e-9 * pow(T,-0.35);} else {k17=9.6e-7 * pow(T,-0.90);}
    double x_Hminus = k1 * xH0 * x_e / ((k2+k16)*xH0 + (k5+k17)*x_p + k15*x_e + R51_n); // assuming equilibrium H-, using full set of terms from Glover & Jappsen 2007 // use for H- opacity

    double a_GP = k2 * x_Hminus * nH0 * clumping_factor; // actual rate for 2-body gas-phase formation, given x_Hminus calculated above from local equilibrium
    
    double b_3B = (6.0e-32/sqrt(sqrt_T) + 2.0e-31/sqrt_T) * nH0 * nH0 * xH0 * clumping_factor_3; // 3-body collisional formation
    double G_LW = 3.3e-11 * urad_G0 * (1./2.); // photo-dissociation (+ionization); note we're assuming a spectral shape identical to the MW background mean, scaling by G0, 1/2 here is to account for nH2 = (1/2) * fH2 * nH because we will solve for fH2 as a mass fraction
    double xi_cr_H2 = (7.525e-16) * (1./2.); // CR dissociation (+ionization), 1/2 here is to account for nH2 = (1/2) * fH2 * nH because we will solve for fH2 as a mass fraction. Using value favored by Indriolo et al for dense GMCs.
//#if defined(COSMIC_RAY_FLUID) || defined(COSMIC_RAY_SUBGRID_LEBRON) || defined(RT_ISRF_BACKGROUND) // scale ionization+dissociation rates with local CR energy density
    xi_cr_H2 = Get_CosmicRayIonizationRate_cgs(i) * (1./2.); // scales following Cummings et al. 2016 to 1.6e-17 per eV/cm^3 [this function now defined for everything; differs slightly from default assumed above but thats just being consistent about the baseline normalization here]
//#endif
    // want to solve the implicit equation: f_f = f_0 + g[f_f]*dt, where g[f_f] = df_dt evaluated at f=f_f, so root-find: dt*g[f_f] + f_0-f_f = 0
    // can write this as a quadtratic: x_a*f^2 - x_b_0*f - xb_LW*f + x_c = 0, where xb_LW is a non-linear function of f accounting for the H2 self-shielding terms
    double G_LW_dt_unshielded = G_LW * dtime_cgs; // LW term without shielding, multiplied by timestep for dimensions needed below
    double x_a_00 = (b_3B + b_H2HI - b_H2H2) * dtime_cgs; // terms quadratic in f -- this term can in principle be positive or negative, usually positive
    double x_b_00 = (a_GP + a_Z + 2.*b_3B + b_H2HI + b_H2ext + xi_cr_H2) * dtime_cgs; // terms linear in f [note sign, pulling the -sign out here] -- positive-definite, not including '1' for f term
    double x_c_00 = (a_GP + a_Z + b_3B) * dtime_cgs; // terms independent of f -- positive-definite, not including f0 term

    // use the previous-timestep value of fH2 to guess the shielding term and then compute the resulting fH2
    fH2_tmp=fH2_initial; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); // calculate the shielding term
    double x_a=x_a_00,x_c=x_c_00,y_a=0,x_b=0,y_b=0,z_a=0, x_b_0 = x_b_00 + y_ss*G_LW_dt_unshielded;
    double dfH2_linear = fH2_tmp*fH2_tmp*x_a_00 - fH2_tmp*x_b_0 + x_c_00; // linear derivative term for dependence: ffinal = fH2_tmp + dfH2_linear to linear, explicit order
    int change_in_fH2_is_small = 0; // key to know what to do below
    if(fabs(dfH2_linear) < 0.01*fH2_tmp) {change_in_fH2_is_small=1;} // this should be a sufficient criterion below
    
    if(change_in_fH2_is_small==1) // change is sufficiently small, can linearly approximate, essentially with aleading-order expansion of the non-linear solution for small timesteps
    {
        double order_2_corr = 2.*fH2_tmp*x_a_00 - x_b_0;
        if(fabs(order_2_corr) < 1.) {fH2 = fH2_tmp + dfH2_linear * (1. + order_2_corr);} else {fH2 = fH2_tmp + dfH2_linear;}
    } else { // we do a nonlinear solve
        x_b_0=x_b_00 + 1.; x_c=x_c_00 + fH2_initial; // x_c and x_b re-incorporate their constant terms in this limit to make the math easier
        y_a=x_a/(x_c + MIN_REAL_NUMBER); // convenient to convert to dimensionless variable needed for checking definite-ness
        #define ROOTFIND_FUNCTION(x) molecfrac_rootfind_function(x, x00, x01, x_b_0, x_c, y_a, G_LW_dt_unshielded); // want to find f_mol such that this is 0
	Q_initial = ROOTFIND_FUNCTION(fH2_initial); 

        x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
        z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
        double fH2_mid = fH2;
        double Q_mid = ROOTFIND_FUNCTION(fH2_mid);

        // OK now let's let the initial and previous-shielding values by the candidate bracket, and if that fails then find another value to bracket the other end
        double ROOTFIND_X_a, ROOTFIND_X_b, ROOTFUNC_a, ROOTFUNC_b;
        ROOTFIND_X_a = fH2_mid; ROOTFUNC_a = Q_mid; 
        ROOTFIND_X_b = fH2_initial; ROOTFUNC_b = Q_initial;	    
            //if not bracketing we must try other bounds
        if(ROOTFUNC_b * ROOTFUNC_a > 0){	
            // lower bound
            x_b=x_b_0+G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // recalculate all terms that depend on the shielding
            fH2_min = DMAX(0,DMIN(1,fH2)); // this serves as a lower-limit for fH2
            Q_min = ROOTFIND_FUNCTION(fH2_min); //molecfrac_rootfind_function(fH2_min, x00, x01, x_b_0, x_c, y_a, G_LW_dt_unshielded);
            ROOTFIND_X_b = fH2_min;  ROOTFUNC_b = Q_min;
            if(Q_min * Q_mid > 0){
            // upper bound
                fH2_tmp=1.; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding    	
                z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
                fH2_max = DMAX(0,DMIN(1,fH2)); // this serves as an upper-limit for fH2	
                Q_max = ROOTFIND_FUNCTION(fH2_max);
                ROOTFIND_X_b = fH2_max;  ROOTFUNC_b = Q_max;
                if(Q_max * Q_mid > 0){
                    if(Q_min*Q_max > 0){
                        ROOTFIND_X_a = 0; ROOTFIND_X_b = 1; ROOTFUNC_a = ROOTFIND_FUNCTION(0); ROOTFUNC_b = ROOTFIND_FUNCTION(1.);
                    } else { 
                        ROOTFIND_X_a = fH2_min; ROOTFIND_X_b = fH2_max; ROOTFUNC_a = Q_min; ROOTFUNC_b = Q_max;
                    }
                }
            }
        }
        
	if(ROOTFUNC_a * ROOTFUNC_b < 0){
        // specify desired relative error in fH2 and call the rootfinder
	    double ROOTFIND_REL_X_tol=1e-3, ROOTFIND_ABS_X_tol=0;
        #include "../system/bracketed_rootfind.h"
	    fH2 = ROOTFIND_X_new;
	    if(ROOTFIND_ITER > MAXITER){PRINT_WARNING("WARNING: Particle %lld did not converge to desired H_2 abundance tolerance\n",(long long)P[i].ID);}
        } else { // must be at 0 or 1 within machine precision of solution but not bracketing - choose the closer bracketing value of 0 or 1
	    if(fabs(ROOTFUNC_a) < fabs(ROOTFUNC_b)){fH2 = ROOTFIND_X_a;} else {fH2 = ROOTFIND_X_b;}
	}
    } // end nonlinear solve part

    if(!isfinite(fH2)) {fH2=0;} else {if(fH2>1) {fH2=1;} else if(fH2<0) {fH2=0;}} // check vs nans, valid values
    SphP[i].MolecularMassFraction_perNeutralH = fH2; // record variable -- this will be used for the next update, meanwhile the total fraction will be used in various routines through the code
    SphP[i].MolecularMassFraction = xH0 * SphP[i].MolecularMassFraction_perNeutralH; // record variable -- this is largely what is needed below
#endif
}

double molecfrac_rootfind_function(double fH2, double x00, double x01, double x_b_0, double x_c, double y_a, double G_LW_dt_unshielded){
    const double x_exp_fac=0.00085, w0 = 0.035, EXPmax=90.;
    double x_ss_1=1.+fH2*x01,
	x_ss_sqrt=sqrt(1.+fH2*x00),
	y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)), 
	x_b=x_b_0+y_ss*G_LW_dt_unshielded, 
	y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
    return 1 + y_a*fH2*fH2 - y_b*fH2; // value of the function we are trying to zero, with the updated value of fH2
}



/* simple subroutine to estimate the dust temperatures in our runs without detailed tracking of these individually [more detailed chemistry models do this] */
double get_equilibrium_dust_temperature_estimate(int i, double shielding_factor_for_exgalbg, double T)
{
#if defined(RT_INFRARED)
    if(i >= 0) {return SphP[i].Dust_Temperature;} // this is pre-computed -- simply return it
#endif
    /* simple three-component model [can do fancier] with cmb, dust, high-energy photons */
    double e_CMB=0.262*All.cf_a3inv/All.cf_atime, T_cmb=2.73/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    double e_IR=0.31, Tdust_ext=DMAX(30.,T_cmb); // Milky way ISRF from Draine (2011), assume peak of dust emission at ~100 microns
    double e_HiEgy=0.66, T_hiegy=5800.; // Milky way ISRF from Draine (2011), assume peak of stellar emission at ~0.6 microns [can still have hot dust, this effect is pretty weak]
#ifdef RT_ISRF_BACKGROUND
    e_IR *= All.InterstellarRadiationFieldStrength; e_HiEgy *= All.InterstellarRadiationFieldStrength; // need to re-scale the assumed ISRF components
    if(!All.ComovingIntegrationOn){
	e_CMB *= pow(1.+All.RadiationBackgroundRedshift,4);
 	T_cmb *= (1.+All.RadiationBackgroundRedshift);
    }
#endif

    if(i >= 0)
    {
#ifdef SINGLE_STAR_SINK_DYNAMICS // treatment using direct dust temperature solver accounting for absorption and gas-dust coupling - want this when capturing the dynamics of dense collapsing cores
	double absorption_rate=0, vol_inv = SphP[i].Density * All.cf_a3inv / P[i].Mass, fac_abs = C_LIGHT_CODE * SphP[i].Density * All.cf_a3inv;
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // we have information about individual radiation bands and their opacities; use these to compute dust absorption rate
	for(int k=0;k<N_RT_FREQ_BINS;k++){
	    if(RT_BAND_IS_IONIZING(k)) {continue;} // skip ionizing bands where the dust cross section is not accounted for
	    absorption_rate += fac_abs * rt_kappa(i,k) * SphP[i].Rad_E_gamma_Pred[k] * vol_inv;
	}
#endif
	absorption_rate += (e_CMB/UNIT_PRESSURE_IN_EV) * fac_abs * rt_kappa_adaptive_IR_band(i,T_cmb,T_cmb,0,1); // CMB absorption; assume cloud is optically-thin to the CMB
#if defined(RT_ISRF_BACKGROUND) // account for additional optical + IR radiation field with extinction
	double column = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i); // column density in code units
	double kappa_IR = rt_kappa_adaptive_IR_band(i,20.,20.,0,1); // assume Trad=20 for IR dust opacity
	double Zfac = 1.;
#ifdef METALS
	Zfac = P[i].Metallicity[0]/All.SolarAbundances[0];
#endif
	double kappa_opt = 180. * Zfac * UNIT_SURFDEN_IN_CGS;
	double tau_opt = kappa_opt*column;
	e_HiEgy += 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // extragalactic UV/optical background
	absorption_rate += fac_abs * kappa_opt * (e_HiEgy/UNIT_PRESSURE_IN_EV) * exp(DMAX(-tau_opt,-100));
	absorption_rate += fac_abs * kappa_IR * ((-0.5*expm1(DMAX(-tau_opt,-100)) * e_HiEgy + e_IR)/UNIT_PRESSURE_IN_EV); // this assumes absorbed ONIR photons are reradiated into IR, factor of 0.5 assumes 1/2 of reradiated IR photons do not go deeper into the cloud
#endif
	// OK now we have our dust absorption rate, let's call the solver
	double Tdust = rt_eqm_dust_temp(i, T, absorption_rate);
	return Tdust;
#endif // SINGLE_STAR_SINK_DYNAMICS

#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
        e_HiEgy=0; e_IR = 0; int k; double E_tot_to_evol_eVcgs = (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_EV;
        for(k=0;k<N_RT_FREQ_BINS;k++) {e_HiEgy+=SphP[i].Rad_E_gamma_Pred[k];}
#if defined(GALSF_FB_FIRE_RT_LONGRANGE)
        e_IR += SphP[i].Rad_E_gamma_Pred[RT_FREQ_BIN_FIRE_IR]; // note IR
#endif
#if defined(RT_INFRARED)
        e_IR += SphP[i].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED]; Tdust_ext = SphP[i].Radiation_Temperature; // note IR [irrelevant b/c of call above, but we'll keep this as a demo]
#endif
        e_HiEgy -= e_IR; // don't double-count the IR component flagged above //
        e_IR *= E_tot_to_evol_eVcgs; e_HiEgy *= E_tot_to_evol_eVcgs;
#endif
    }
    e_HiEgy += shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
    double Tdust_eqm = 10.; // arbitrary initial value //
    if(Tdust_ext*e_IR < 1.e-10 * (T_cmb*e_CMB + T_hiegy*e_HiEgy)) { // IR term is totally negligible [or zero exactly], use simpler expression assuming constant temperature for it to avoid sensitivity to floating-pt errors //
        Tdust_eqm = 2.92 * pow(Tdust_ext*e_IR + T_cmb*e_CMB + T_hiegy*e_HiEgy, 1./5.); // approximate equilibrium temp assuming Q~1/lambda [beta=1 opacity law], assuming background IR temp is a fixed constant [relevant in IR-thin limit, but we don't know T_rad, so this is a guess anyways]
    } else { // IR term is not vanishingly small. we will assume the IR radiation temperature is equal to the local Tdust. lacking any direct evolution of that field, this is a good proxy, and exact in the locally-IR-optically-thick limit. in the locally-IR-thin limit it slightly under-estimates Tdust, but usually in that limit the other terms dominate anyways, so this is pretty safe //
        double T0=2.92, q=pow(T0*e_IR,0.25), y=(T_cmb*e_CMB + T_hiegy*e_HiEgy)/(T0*e_IR*q); if(y<=1) {Tdust_eqm=T0*q*(0.8+sqrt(0.04+0.1*y));} else {double y5=pow(y,0.2), y5_3=y5*y5*y5, y5_4=y5_3*y5; Tdust_eqm=T0*q*(1.+15.*y5_4+sqrt(1.+30.*y5_4+25.*y5_4*y5_4))/(20.*y5_3);} // this gives an extremely accurate and exactly-joined solution to the full quintic equation assuming T_rad_IR=T_dust
    }
    return DMAX(DMIN(Tdust_eqm , 2000.) , 1.); // limit at sublimation temperature or some very low temp //
}



/* Calculates the coefficient for gas-dust collisional heat transfer, such that LambdaDust = gas_dust_heating_coeff * (T-Tdust) in erg cm^3 s^-1. */
double gas_dust_heating_coeff(int i, double T, double Tdust) 
{
    double Z_sol=1;
#ifdef METALS
    if(i>=0) {Z_sol = P[i].Metallicity[0]/All.SolarAbundances[0];}
#endif
    double fdust = return_dust_to_metals_ratio_vs_solar(i,Tdust); // accounting for dust destruction; we avoid calling the function for this because it can create a circular dependency
    return 1.116e-32 * sqrt(T)*(1.-0.8*exp(-75./T)) * Z_sol * fdust;  // Meijerink & Spaans 2005; Hollenbach & McKee 1979,1989. Assumes 10 Angstrom minimum grain size.
}



/* this function estimates the free electron fraction from heavy ions, assuming a simple mix of cold molecular gas, Mg, and dust, with the ions from singly-ionized Mg, to prevent artificially low free electron fractions */
double return_electron_fraction_from_heavy_ions(int target, double temperature, double density_cgs, double n_elec_HHe)
{
    if(All.ComovingIntegrationOn) {double rhofac=density_cgs/(1000.*COSMIC_BARYON_DENSITY_CGS); if(rhofac<0.2) {return 0;}} // ignore these reactions in the IGM
    double zeta_cr=1.0e-17, f_dustgas=0.01, n_ion_max=4.1533e-5, XH=HYDROGEN_MASSFRAC; // cosmic ray ionization rate (fixed as constant for non-CR runs) and dust-to-gas ratio
    if(target >= 0) {zeta_cr = Get_CosmicRayIonizationRate_cgs(target);} // convert to ionization rate, using models as in Cummings et al. 2016
#ifdef METALS
    if(target>=0) {f_dustgas=0.5*P[target].Metallicity[0]*return_dust_to_metals_ratio_vs_solar(target,0) + 1.e-15;} // constant dust-to-metals ratio [floor purely numerical here, needed to avoid a spurious divergence below]
#ifdef COOL_METAL_LINES_BY_SPECIES
    if(target>=0) {n_ion_max = (All.SolarAbundances[6]/24.3)/XH;} // limit, to avoid over-ionization at low metallicities
#endif
#endif
    /* Regime I: highly/photo-ionized, any contributions here would be negligible -- no need to continue */
    if(n_elec_HHe > 0.01) {return n_ion_max;} // contribute something negligible, doesn't matter here //
    double a_grain_micron=0.1, m_ion=24.305*PROTONMASS_CGS, mu_eff=2.38, m_neutrals=mu_eff*PROTONMASS_CGS, m_grain=4.189e-12*(2.4)*a_grain_micron*a_grain_micron*a_grain_micron, ngrain_ngas=(m_neutrals/m_grain)*f_dustgas; // effective size of grains that matter at these densities, and ions [here Mg] that dominate
    double k_ei=9.77e-8, y0=sqrt(m_ion/ELECTRONMASS_CGS);
    double y = exp(1.)*y0, ln_oneplusy=log(1.+y), psi_0 = 1.-ln_oneplusy + ln_oneplusy/(1.+ln_oneplusy) * log(ln_oneplusy*(1.+1./y)); // changed convention: using second-order expansion to solve for psi, psi should be < 0, and electron absorption should be suppressed for larger psi because of grain charge leading to coulomb repulsion
    double k_eg_00=0.0195*a_grain_micron*a_grain_micron*sqrt(temperature), k_eg_0=k_eg_00*exp(psi_0);
    double n_crit=k_ei*zeta_cr/(k_eg_0*ngrain_ngas*k_eg_0*ngrain_ngas), n_eff=density_cgs/m_neutrals;

    /* Regime II: CR-ionized with high enough ionization fraction that gas-phase recombinations dominate */
    if(n_eff < 0.01*n_crit) {return DMIN(n_ion_max , sqrt(zeta_cr / (k_ei * n_eff)) / (XH*mu_eff) );} // CR-dominated off gas with recombination via ions -- put into appropriate units here, but basically just what we would expect here. mu_eff factor transforms from n_e/n_neutrals to n_e/n_H because neutrals if molecular not as numerous
    if(n_eff < 100.*n_crit) {return DMIN(n_ion_max , 0.5 * (sqrt(4.*n_crit/n_eff + 1.) - 1.) * (k_eg_0*ngrain_ngas)/(k_ei*XH*mu_eff));} // interpolates between gas-recombination and dust-dominated regimes (just a simple interpolation function)
    
    /* Regime III: recombination dominated by dust, but dust has a 'normal' efficiency of absorbing grains, and most of the charge is still free+ion (n_ions ~ n_free_electrons >> Z_grains*n_grains, even if -rate- of grain absorption of free e- dominates over k_ei) */
    double psi_fac=16.71/(a_grain_micron*temperature), alpha=zeta_cr*psi_fac/(k_eg_00 * ngrain_ngas*ngrain_ngas * n_eff), alpha_min=0.02, alpha_max=10.; /* Z*psi_fac = Z*e^2/(a_grain*kB*T) to dimensionless [psi] units; alpha=prefactor in charge equation: psi = alpha*(exp[-psi] - y/(1-psi)); this alpha factor gives n_e/n_charge_summed_grains = (this)/(exp[psi]*psi) -- when this gets to smaller than ~10, can start to deplete all the charge onto grains, need to move to the next regime [note here psi only puts you more-safely into free e- regime, so dont need solution below if psi very small, hence ignore it here to be conservative */
    if(alpha > alpha_max) {return DMIN(n_ion_max , zeta_cr / (k_eg_0*ngrain_ngas * XH*mu_eff*n_eff ) );} // regime III limit above: here XH*mu_eff factor again just converts to the correct units. here balancing ionization [zeta_cr*n_neutrals] with absorption by dust grains [k_eg_0*ngrain_ngas*n_neutrals*x_e*n_neutrals]
    
    /* Regime IV: recombination dominated by dust and grains dominate the charge, strongly suppressing the free charges (n_ion ~ Z_grains*n_grains >> n_free_electrons) */
    if(alpha < 1.e-4) {return DMIN(n_ion_max , zeta_cr / (k_eg_00*ngrain_ngas * XH*mu_eff*n_eff) );} // alpha very small means that there is a very tiny number of free e- per dust grain to potentially absorb. this means the mean grain charge must be << 1, so psi->small, negligible correction here, just same expression as above with k_eg_0 (no psi correction, or k_eg calculated without any charge correction here. note in this limit, the ion fraction is larger than the free electron fraction by a ratio sqrt[m_ion/m_e] ~ k_eg/k_ig [see Keith+Wardle 2014, section 3.2, which adopts this limit
    if(alpha < alpha_min) {double psi=0.5*(1.-sqrt(1.+4.*y0*alpha)); return DMIN(n_ion_max , zeta_cr / (k_eg_00*exp(psi)*ngrain_ngas * XH*mu_eff*n_eff) );} // small-psi-limit. when alpha is sufficiently small, here in limit where [already above] can safely neglect k_ei, so solution for x_e is same as above/below, but with a different psi, here need to solve psi = exp[-psi]*alpha - y0*alpha/(1-psi). if alpha is small, psi should also be small, can drop the exp[-psi]*alpha term and this becomes a reasonable approximation here.
    double psi_xmin=0.5*(1.-sqrt(1.+4.*y0*alpha_min)), psi=psi_0 + (psi_xmin-psi_0)*2./(1.+alpha/alpha_min); // this interpolates between the asymptotic limmits at low and high alpha, where we can obtain high-accuracy solutions here. this is a completely ad-hoc fitting function for the numerical solutions in this intermediate range, which only relates to a rather modest range of alpha and is ensured to interpolate correctly to the extremes we already pulled out above.
    return DMIN(n_ion_max , zeta_cr / (k_eg_00*exp(psi)*ngrain_ngas * XH*mu_eff*n_eff));
}




/* this function evaluates Compton heating+cooling rates and synchrotron cooling for thermal gas populations, accounting for the
    explicitly-evolved radiation field if it is evolved (otherwise assuming a standard background), and B-fields if they
    are evolved, as well as the proper relativistic or non-relativistic effects and two-temperature plasma effects. */
double evaluate_Compton_heating_cooling_rate(int target, double T, double nHcgs, double n_elec, double shielding_factor_for_exgalbg)
{
    double Lambda = 0;
    double compton_prefac_eV = 2.16e-35 / nHcgs; // multiply field in eV/cm^3 by this and temperature difference to obtain rate

    double e_CMB_eV=0.262*All.cf_a3inv/All.cf_atime, T_cmb = 2.73/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    Lambda += compton_prefac_eV * n_elec * e_CMB_eV * (T-T_cmb);

    double e_UVB_eV = shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
    Lambda += compton_prefac_eV * n_elec * e_UVB_eV * (T-2.e4); // assume very crude approx Compton temp ~2e4 for UVB
    
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
    if(target >= 0)
    {
        int k; double E_tot_to_evol_eVcgs = (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_EV;
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            double e_tmp = SphP[target].Rad_E_gamma_Pred[k] * E_tot_to_evol_eVcgs, Teff = 0;
            
#if defined(GALSF_FB_FIRE_RT_LONGRANGE) /* three-band (UV, OPTICAL, IR) approximate spectra for stars as used in the FIRE (Hopkins et al.) models */
            if(k==RT_FREQ_BIN_FIRE_IR) {Teff=30.;}
            if(k==RT_FREQ_BIN_FIRE_OPT) {Teff=4000.;}
            if(k==RT_FREQ_BIN_FIRE_UV) {Teff=15000.;}
#endif
#if defined(RT_INFRARED) /* special mid-through-far infrared band, which includes IR radiation temperature evolution */
            if(k==RT_FREQ_BIN_INFRARED) {Teff=SphP[target].Radiation_Temperature;}
#endif
#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models; from 0.41-3.4 eV */
            if(k==RT_FREQ_BIN_OPTICAL_NIR) {Teff=2800.;}
#endif
#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models; from 3.4-8 eV */
            if(k==RT_FREQ_BIN_NUV) {Teff=12000.;}
#endif
#if defined(RT_PHOTOELECTRIC) /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
            if(k==RT_FREQ_BIN_PHOTOELECTRIC) {Teff=24400.;}
#endif
#if defined(RT_LYMAN_WERNER) /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
            if(k==RT_FREQ_BIN_LYMAN_WERNER) {Teff=28800.;}
#endif
#if defined(RT_CHEM_PHOTOION) /* Hydrogen and Helium ionizing bands: H0 here */
            if(k==RT_FREQ_BIN_H0) {Teff=2340.*rt_nu_eff_eV[k];}
#endif
#if defined(RT_PHOTOION_MULTIFREQUENCY) /* Hydrogen and Helium ionizing bands: He bands */
            if(k==RT_FREQ_BIN_He0 || k==RT_FREQ_BIN_He1 || k==RT_FREQ_BIN_He2) {Teff=2340.*rt_nu_eff_eV[k];}
#endif
#if defined(RT_SOFT_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
            if(k==RT_FREQ_BIN_SOFT_XRAY) {Teff=3.6e6;}
#endif
#if defined(RT_HARD_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
            if(k==RT_FREQ_BIN_HARD_XRAY) {Teff=1.7e7;}
#endif
            if(Teff < 3.e4) {e_tmp *= n_elec;} // low-energy radiation acts inefficiently on neutrals here
            Lambda += compton_prefac_eV * e_tmp * (T - Teff); // add to compton heating/cooling terms
        }
    }
#else // no explicit RHD terms evolved, so assume a MW-like ISRF instead
    double e_IR_eV=0.31, T_IR=DMAX(30.,T_cmb); // Milky way ISRF from Draine (2011), assume peak of dust emission at ~100 microns
    double e_OUV_eV=0.66, T_OUV=5800.; // Milky way ISRF from Draine (2011), assume peak of stellar emission at ~0.6 microns [can still have hot dust, this effect is pretty weak]
    Lambda += compton_prefac_eV * n_elec * (e_IR_eV*(T-T_IR) + e_OUV_eV*(T-T_OUV));
#endif

#ifdef BH_COMPTON_HEATING /* custom band to represent (non)relativistic X-ray compton cooling from an AGN source without full RHD */
    if(target >= 0)
    {
        double e_agn = (SphP[target].Rad_Flux_AGN * UNIT_FLUX_IN_CGS) / (C_LIGHT_CGS * ELECTRONVOLT_IN_ERGS), T_agn=2.e7; /* approximate from Sazonov et al. */
        Lambda += compton_prefac_eV * e_agn * (T-T_agn); // since the heating here is primarily hard X-rays, and the cooling only relevant for very high temps, do not have an n_elec here
    }
#endif

#ifdef MAGNETIC /* include sychrotron losses as well as long as we're here, since these scale more or less identically just using the magnetic instead of radiation energy */
    if(target >= 0)
    {
        double b_muG = get_cell_Bfield_in_microGauss(target), U_mag_ev=0.0248342*b_muG*b_muG, T_rad_background_at_emission = get_background_radiation_temperature_for_emission_corrections(target);
        Lambda += compton_prefac_eV * U_mag_ev * (T-T_rad_background_at_emission); // synchrotron losses proportional to temperature (non-relativistic here), as inverse compton, just here without needing to worry about "T-T_eff", as if T_eff->0
    }
#endif

    double T_eff_for_relativistic_corr = T; /* used below, but can be corrected */
    if(Lambda > 0) /* per CAFG's calculations, we should note that at very high temperatures, the rate-limiting step may be the Coulomb collisions moving energy from protons to e-; which if slow will prevent efficient e- cooling */
    {
        double Lambda_limiter_var = 1.483e34 * Lambda*Lambda*T; /* = (Lambda/2.6e-22)^2 * (T/1e9): if this >> 1, follow CAFG and cap at cooling rate assuming equilibrium e- temp from Coulomb exchange balancing compton */
        if(Lambda_limiter_var > 1) {Lambda_limiter_var = 1./pow(Lambda_limiter_var,0.2); Lambda*=Lambda_limiter_var; T_eff_for_relativistic_corr*=Lambda_limiter_var;}
    }
    if(T_eff_for_relativistic_corr > 3.e7) {Lambda *= (T_eff_for_relativistic_corr/1.5e9) / (1-exp(-T_eff_for_relativistic_corr/1.5e9));} /* relativistic correction term, becomes important at >1e9 K, enhancing rate */
    
    return Lambda;
}


/* this function defines an effective background radiation temperature for purposes of computing the emission corrections above */
double get_background_radiation_temperature_for_emission_corrections(int target)
{
    double T_cmb = 2.73/All.cf_atime;
#ifdef RT_ISRF_BACKGROUND
    if(!All.ComovingIntegrationOn){
	T_cmb *= 1.+All.RadiationBackgroundRedshift;
    }
#endif	
#if defined(RT_INFRARED)
    if(target >= 0)
    {
        double e_tmp_CMB = 0.262*All.cf_a3inv/All.cf_atime;
#ifdef RT_ISRF_BACKGROUND
	if(!All.ComovingIntegrationOn){
	    e_tmp_CMB *= pow(1.+All.RadiationBackgroundRedshift,4);
	}
#endif
	double e_tot_to_evol_eVcgs = (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_EV;
        double e_tmp_IR = SphP[target].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED] * e_tot_to_evol_eVcgs, T_tmp_IR = SphP[target].Radiation_Temperature;
        return (e_tmp_IR * T_tmp_IR + e_tmp_CMB * T_cmb) / (e_tmp_IR + e_tmp_CMB); // if evolving IR band, use sum of it plus cmb for background
    }
#endif
    return T_cmb; // default to assume total radiation temperature is cmb-dominated
}



/*  this function computes the self-consistent temperature and electron fraction */
double ThermalProperties(double u, double rho, int target, double *mu_guess, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess)
{
#if defined(CHIMES)
    int i = target; *ne_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_elec]]; *nH0_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HI]];
    *nHp_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HII]]; *nHe0_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HeI]];
    *nHep_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HeII]]; *nHepp_guess = ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HeIII]];
    double temp = ChimesGasVars[target].temperature;
    *mu_guess = Get_Gas_Mean_Molecular_Weight_mu(temp, rho, nH0_guess, ne_guess, 0, target);
    return temp;
#else
    if(target >= 0) {*ne_guess=SphP[target].Ne; *nH0_guess = DMAX(0,DMIN(1,1.-( *ne_guess / 1.2 )));} else {*ne_guess=1.; *nH0_guess=0.;}
    rho *= UNIT_DENSITY_IN_CGS; u *= UNIT_SPECEGY_IN_CGS;   /* convert to physical cgs units */
    double temp = convert_u_to_temp(u, rho, target, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess);
    return temp;
#endif
}



/* function to return the local multiplier relative to the UVB model to account in some local RHD models for local ionizing sources */
double return_local_gammamultiplier(int target)
{
#if defined(GALSF_FB_FIRE_RT_UVHEATING) && !defined(CHIMES)
    if((target >= 0) && (gJH0 > 0))
    {
        double local_gammamultiplier = SphP[target].Rad_Flux_EUV * 2.29e-10; // converts to GammaHI for typical SED (rad_uv normalized to Habing)
        local_gammamultiplier = 1.0 + local_gammamultiplier / gJH0; // this needs to live here in cooling.c where gJH0 is declared as a global shared variable!
        if(!isfinite(local_gammamultiplier)) {local_gammamultiplier=1;} // check for divide-by zero errors
        return DMAX(1., DMIN(2./gJH0, DMIN(1.e20, local_gammamultiplier))); // set a cap at large number and at value spectrum would have at extremely short distance (~1000 au) from a super-luminous O-star to prevent unphysically high values
    }
#endif
    return 1;
}


/* function to attenuate the UVB to model self-shielding in optically-thin simulations */
double return_uvb_shieldfac(int target, double gamma_12, double nHcgs, double logT)
{
#ifdef GALSF_EFFECTIVE_EQS
    return 1; // self-shielding is implicit in the sub-grid model already //
#endif
    double NH_SS_z, NH_SS = 0.0123; /* CAFG: H number density above which we assume no ionizing bkg (proper cm^-3): note this is a factor ~2 higher than Schaye and Rahmati normalization, owing to CAFG calibration (lower effective cross-section, owing to different UVB spectrum and potentially clumping factors, etc) */
    if(gamma_12>0) {NH_SS_z = NH_SS*pow(gamma_12,0.66)*pow(10.,0.173*(logT-4.));} else {NH_SS_z = NH_SS*pow(10.,0.173*(logT-4.));}
    double q_SS = nHcgs/NH_SS_z;
#ifdef COOL_UVB_SELFSHIELD_RAHMATI
    return 0.98 / pow(1 + pow(q_SS,1.64), 2.28) + 0.02 / pow(1 + q_SS*(1.+1.e-4*nHcgs*nHcgs*nHcgs*nHcgs), 0.84); // from Rahmati et al. 2012: gives gentler cutoff at high densities. but we need to modify it with the extra 1+(nHcgs/10)^4 denominator term since at very high nH, this cuts off much too-slowly (as nH^-0.84), which means UVB heating can be stronger than molecular cooling even at densities >> 1e4
#else
    return 1./(1.+q_SS*(1.+q_SS/2.*(1.+q_SS/3.*(1.+q_SS/4.*(1.+q_SS/5.*(1.+q_SS/6.*q_SS)))))); // this is exp(-q) down to ~1e-5, then a bit shallower, but more than sufficient approximation here //
#endif
}



#ifdef CHIMES
/* This routine updates the ChimesGasVars structure for particle target. */
void chimes_update_gas_vars(int target)
{
  double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target);
  double u_old_cgs = DMAX(All.MinEgySpec, SphP[target].InternalEnergy) * UNIT_SPECEGY_IN_CGS;
  double rho_cgs = SphP[target].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

#ifdef COOL_METAL_LINES_BY_SPECIES
  double H_mass_fraction = 1.0 - (P[target].Metallicity[0] + P[target].Metallicity[1]);
#else
  double H_mass_fraction = XH;
#endif

  ChimesGasVars[target].temperature = (ChimesFloat) chimes_convert_u_to_temp(u_old_cgs, rho_cgs, target);
  ChimesGasVars[target].nH_tot = (ChimesFloat) (H_mass_fraction * rho_cgs / PROTONMASS_CGS);
  ChimesGasVars[target].ThermEvolOn = All.ChimesThermEvolOn;

  // If there is an EoS, need to set TempFloor to that instead.
  ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(All.MinGasTemp, 10.1);
#if (FLAG_NOT_IN_PUBLIC_CODE <= 2) && defined(GALSF_FB_FIRE_RT_HIIHEATING)
    if (SphP[target].DelayTimeHII > 0) {ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(HIIRegion_Temp, 10.1);}
        else {ChimesGasVars[target].TempFloor = (ChimesFloat) DMAX(All.MinGasTemp, 10.1);}
#endif

  // Flag to control how the temperature
  // floor is implemented in CHIMES.
  ChimesGasVars[target].temp_floor_mode = 1;

  // Extragalactic UV background
  ChimesGasVars[target].isotropic_photon_density[0] = chimes_table_spectra.isotropic_photon_density[0];
  ChimesGasVars[target].isotropic_photon_density[0] *= chimes_rad_field_norm_factor;

  ChimesGasVars[target].G0_parameter[0] = chimes_table_spectra.G0_parameter[0];
  ChimesGasVars[target].H2_dissocJ[0] =  chimes_table_spectra.H2_dissocJ[0];

#ifdef CHIMES_STELLAR_FLUXES
  int kc;
  for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++)
    {
      ChimesGasVars[target].isotropic_photon_density[kc + 1] = (ChimesFloat) (SphP[target].Chimes_fluxPhotIon[kc] / C_LIGHT_CGS);

#ifdef CHIMES_HII_REGIONS
    if(SphP[target].DelayTimeHII > 0) {
        ChimesGasVars[target].isotropic_photon_density[kc + 1] += (ChimesFloat) (SphP[target].Chimes_fluxPhotIon_HII[kc] / C_LIGHT_CGS);
        ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) ((SphP[target].Chimes_G0[kc] + SphP[target].Chimes_G0_HII[kc]) / DMAX((SphP[target].Chimes_fluxPhotIon[kc] + SphP[target].Chimes_fluxPhotIon_HII[kc]), 1.0e-300));
	} else {ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) (SphP[target].Chimes_G0[kc] / DMAX(SphP[target].Chimes_fluxPhotIon[kc], 1.0e-300));}
    ChimesGasVars[target].H2_dissocJ[kc + 1] = (ChimesFloat) (ChimesGasVars[target].G0_parameter[kc + 1] * (chimes_table_spectra.H2_dissocJ[kc + 1] / chimes_table_spectra.G0_parameter[kc + 1]));
#else
      ChimesGasVars[target].G0_parameter[kc + 1] = (ChimesFloat) (SphP[target].Chimes_G0[kc] / DMAX(SphP[target].Chimes_fluxPhotIon[kc], 1.0e-300));
      ChimesGasVars[target].H2_dissocJ[kc + 1] = (ChimesFloat) (ChimesGasVars[target].G0_parameter[kc + 1] * (chimes_table_spectra.H2_dissocJ[kc + 1] / chimes_table_spectra.G0_parameter[kc + 1]));
#endif
    }
#endif

  ChimesGasVars[target].cr_rate = (ChimesFloat) cr_rate;  // For now, assume a constant cr_rate.
  ChimesGasVars[target].hydro_timestep = (ChimesFloat) (dt * UNIT_TIME_IN_CGS);

  ChimesGasVars[target].ForceEqOn = ChimesEqmMode;
  ChimesGasVars[target].divVel = (ChimesFloat) P[target].Particle_DivVel / UNIT_TIME_IN_CGS;
  if (All.ComovingIntegrationOn)
    {
      ChimesGasVars[target].divVel *= (ChimesFloat) All.cf_a2inv;
      ChimesGasVars[target].divVel += (ChimesFloat) (3 * All.cf_hubble_a / UNIT_TIME_IN_CGS);  /* Term due to Hubble expansion */
    }
  ChimesGasVars[target].divVel = (ChimesFloat) fabs(ChimesGasVars[target].divVel);

#ifndef COOLING_OPERATOR_SPLIT
  if(SphP[target].CoolingIsOperatorSplitThisTimestep==0) {ChimesGasVars[target].constant_heating_rate = ChimesGasVars[target].nH_tot * ((ChimesFloat) SphP[target].DtInternalEnergy);}
#else
  ChimesGasVars[target].constant_heating_rate = 0.0;
#endif

#ifdef CHIMES_SOBOLEV_SHIELDING
  double surface_density;
  surface_density = evaluate_NH_from_GradRho(P[target].GradRho,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target) * UNIT_SURFDEN_IN_CGS; // converts to cgs
  ChimesGasVars[target].cell_size = (ChimesFloat) (shielding_length_factor * surface_density / rho_cgs);
#else
  ChimesGasVars[target].cell_size = 1.0;
#endif

  ChimesGasVars[target].doppler_broad = 7.1;  // km/s. For now, just set this constant. Thermal broadening is also added within CHIMES.

  ChimesGasVars[target].InitIonState = ChimesInitIonState;

#ifdef CHIMES_HII_REGIONS
  if (SphP[target].DelayTimeHII > 0.0) {ChimesGasVars[target].cell_size = 1.0;} // switch of shielding in HII regions
#endif

#ifdef CHIMES_TURB_DIFF_IONS
  chimes_update_turbulent_abundances(target, 0);
#endif

#if defined(COOL_METAL_LINES_BY_SPECIES) && !defined(GALSF_FB_NOENRICHMENT)
  chimes_update_element_abundances(target);
#endif

  return;
}

#ifdef COOL_METAL_LINES_BY_SPECIES
/* This routine re-computes the element abundances from
 * the metallicity array and updates the individual ion
 * abundances accordingly. */
void chimes_update_element_abundances(int i)
{
  double H_mass_fraction = 1.0 - (P[i].Metallicity[0] + P[i].Metallicity[1]);

  /* Update the element abundances in ChimesGasVars. */
  ChimesGasVars[i].element_abundances[0] = (ChimesFloat) (P[i].Metallicity[1] / (4.0 * H_mass_fraction));   // He
  ChimesGasVars[i].element_abundances[1] = (ChimesFloat) (P[i].Metallicity[2] / (12.0 * H_mass_fraction));  // C
  ChimesGasVars[i].element_abundances[2] = (ChimesFloat) (P[i].Metallicity[3] / (14.0 * H_mass_fraction));  // N
  ChimesGasVars[i].element_abundances[3] = (ChimesFloat) (P[i].Metallicity[4] / (16.0 * H_mass_fraction));  // O
  ChimesGasVars[i].element_abundances[4] = (ChimesFloat) (P[i].Metallicity[5] / (20.0 * H_mass_fraction));  // Ne
  ChimesGasVars[i].element_abundances[5] = (ChimesFloat) (P[i].Metallicity[6] / (24.0 * H_mass_fraction));  // Mg
  ChimesGasVars[i].element_abundances[6] = (ChimesFloat) (P[i].Metallicity[7] / (28.0 * H_mass_fraction));  // Si
  ChimesGasVars[i].element_abundances[7] = (ChimesFloat) (P[i].Metallicity[8] / (32.0 * H_mass_fraction));  // S
  ChimesGasVars[i].element_abundances[8] = (ChimesFloat) (P[i].Metallicity[9] / (40.0 * H_mass_fraction));  // Ca
  ChimesGasVars[i].element_abundances[9] = (ChimesFloat) (P[i].Metallicity[10] / (56.0 * H_mass_fraction)); // Fe

  ChimesGasVars[i].metallicity = (ChimesFloat) (P[i].Metallicity[0] / 0.0129);  // In Zsol. CHIMES uses Zsol = 0.0129.
  ChimesGasVars[i].dust_ratio = ChimesGasVars[i].metallicity;

#ifdef CHIMES_METAL_DEPLETION
#ifdef _OPENMP
  int ThisThread = omp_get_thread_num();
#else
  int ThisThread = 0;
#endif
  chimes_compute_depletions(ChimesGasVars[i].nH_tot, ChimesGasVars[i].temperature, ThisThread);

  ChimesGasVars[i].element_abundances[1] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[0]; // C
  ChimesGasVars[i].element_abundances[2] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[1]; // N
  ChimesGasVars[i].element_abundances[3] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[2]; // O
  ChimesGasVars[i].element_abundances[5] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[3]; // Mg
  ChimesGasVars[i].element_abundances[6] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[4]; // Si
  ChimesGasVars[i].element_abundances[7] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[5]; // S
  ChimesGasVars[i].element_abundances[9] *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDepletionFactors[6]; // Fe

  ChimesGasVars[i].dust_ratio *= (ChimesFloat) ChimesDepletionData[ThisThread].ChimesDustRatio;
#endif // CHIMES_METAL_DEPLETION

  /* The element abundances may have changed, so use the check_constraint_equations()
   * routine to update the abundance arrays accordingly. If metals have been injected
   * into a particle, it is distributed across all of the metal's atomic/ionic/molecular
   * species, preserving the ion and molecule fractions of that element. */

  check_constraint_equations(&(ChimesGasVars[i]), &ChimesGlobalVars);
}
#endif // COOL_METAL_LINES_BY_SPECIES

#ifdef CHIMES_TURB_DIFF_IONS
/* mode == 0: re-compute the CHIMES abundance array from the ChimesNIons and ChimesHtot
 *            arrays that are used to track turbulent diffusion of ions and molecules.
 * mode == 1: update the ChimeSNIons array from the current CHIMES abundance array.
 */
void chimes_update_turbulent_abundances(int i, int mode)
{
  int k_species;
  double NHtot = (1.0 - (P[i].Metallicity[0] + P[i].Metallicity[1])) * (P[i].Mass * UNIT_MASS_IN_CGS) / PROTONMASS_CGS;

  if (mode == 0)
    {
      for (k_species = 0; k_species < ChimesGlobalVars.totalNumberOfSpecies; k_species++)
	ChimesGasVars[i].abundances[k_species] = (ChimesFloat) (SphP[i].ChimesNIons[k_species] / NHtot);
    }
  else
    {
      for (k_species = 0; k_species < ChimesGlobalVars.totalNumberOfSpecies; k_species++)
	SphP[i].ChimesNIons[k_species] = ((double) ChimesGasVars[i].abundances[k_species]) * NHtot;
    }
}
#endif // CHIMES_TURB_DIFF_IONS

#ifdef CHIMES_METAL_DEPLETION
void chimes_init_depletion_data(void)
{
  int i;

#ifdef _OPENMP
  ChimesDepletionData = (struct Chimes_depletion_data_structure *) malloc(maxThreads * sizeof(struct Chimes_depletion_data_structure));
#else
  ChimesDepletionData = (struct Chimes_depletion_data_structure *) malloc(sizeof(struct Chimes_depletion_data_structure));
#endif

  // Elements in Jenkins (2009) in the order
  // C, N, O, Mg, Si, P, S, Cl, Ti, Cr, Mn, Fe,
  // Ni, Cu, Zn, Ge, Kr
  // Solar abundances are as mass fractions,
  // taken from the Cloudy default values, as
  // used in CHIMES.
  double SolarAbund[DEPL_N_ELEM] = {2.07e-3, 8.36e-4, 5.49e-3, 5.91e-4, 6.83e-4, 7.01e-6, 4.09e-4, 4.72e-6, 3.56e-6, 1.72e-5, 1.12e-5, 1.1e-3, 7.42e-5, 7.32e-7, 1.83e-6, 2.58e-7, 1.36e-7};

  // Fit parameters, using equation 10 of Jenkins (2009).
  // Where possible, we take the updated fit parameters
  // A2 and B2 from De Cia et al. (2016), which we convert
  // to the Ax and Bx of J09 (with zx = 0). Otherwise, we
  // use the original J09 parameters. We list these in the
  // order Ax, Bx, zx.
  double DeplPars[DEPL_N_ELEM][3] = {{-0.101, -0.193, 0.803}, // C
				     {0.0, -0.109, 0.55},     // N
				     {-0.101, -0.172, 0.0},     // O
				     {-0.412, -0.648, 0.0},     // Mg
				     {-0.426, -0.669, 0.0},     // Si
				     {-0.068, -0.091, 0.0},      // P
				     {-0.189, -0.324, 0.0},     // S
				     {-1.242, -0.314, 0.609}, // Cl
				     {-2.048, -1.957, 0.43},  // Ti
				     {-0.892, -1.188, 0.0},      // Cr
				     {-0.642, -0.923, 0.0},      // Mn
				     {-0.851, -1.287, 0.0},     // Fe
				     {-1.490, -1.829, 0.599}, // Ni
				     {-0.710, -1.102, 0.711}, // Cu
				     {-0.182, -0.274, 0.0},       // Zn
				     {-0.615, -0.725, 0.69},  // Ge
				     {-0.166, -0.332, 0.684}}; // Kr

  int n_thread;
#ifdef _OPENMP
  n_thread = maxThreads;
#else
  n_thread = 1;
#endif

  for (i = 0; i < n_thread; i++)
    {
      memcpy(ChimesDepletionData[i].SolarAbund, SolarAbund, DEPL_N_ELEM * sizeof(double));
      memcpy(ChimesDepletionData[i].DeplPars, DeplPars, DEPL_N_ELEM * 3 * sizeof(double));
      // DustToGasSaturated is the dust to gas ratio when F_star = 1.0, i.e. at maximum depletion onto grains.
      ChimesDepletionData[i].DustToGasSaturated = 5.9688e-03;
    }
}

/* Computes the linear fits as in Jenkins (2009)
 * or De Cia et al. (2016). Note that this returns
 * log10(fraction in the gas phase) */
double chimes_depletion_linear_fit(double nH, double T, double Ax, double Bx, double zx)
{
  // First, compute the parameter F_star, using the
  // best-fit relation from Fig. 16 of Jenkins (2009).
  double F_star = 0.772 + (0.461 * log10(nH));
  double depletion;

  // Limit F_star to be no greater than unity
  if (F_star > 1.0)
    F_star = 1.0;

  // All metals are in the gas phase if T > 10^6 K
  if (T > 1.0e6)
    return 0.0;
  else
    {
      depletion = Bx + (Ax * (F_star - zx));

      // Limit depletion to no greater than 0.0 (remember: it is a log)
      if (depletion > 0.0)
	return 0.0;
      else
	return depletion;
    }
}

void chimes_compute_depletions(double nH, double T, int thread_id)
{
  int i;
  double pars[DEPL_N_ELEM][3];
  memcpy(pars, ChimesDepletionData[thread_id].DeplPars, DEPL_N_ELEM * 3 * sizeof(double));

  // ChimesDepletionFactors array is for the metals
  // in the order C, N, O, Mg, Si, S, Fe. The other
  // metals in CHIMES are not depleted.
  ChimesDepletionData[thread_id].ChimesDepletionFactors[0] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[0][0], pars[0][1], pars[0][2])); // C
  ChimesDepletionData[thread_id].ChimesDepletionFactors[1] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[1][0], pars[1][1], pars[1][2])); // N
  ChimesDepletionData[thread_id].ChimesDepletionFactors[2] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[2][0], pars[2][1], pars[2][2])); // O
  ChimesDepletionData[thread_id].ChimesDepletionFactors[3] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[3][0], pars[3][1], pars[3][2])); // Mg
  ChimesDepletionData[thread_id].ChimesDepletionFactors[4] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[4][0], pars[4][1], pars[4][2])); // Si
  ChimesDepletionData[thread_id].ChimesDepletionFactors[5] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[6][0], pars[6][1], pars[6][2])); // S
  ChimesDepletionData[thread_id].ChimesDepletionFactors[6] = pow(10.0, chimes_depletion_linear_fit(nH, T, pars[11][0], pars[11][1], pars[11][2])); // Fe

  // The dust abundance as used in CHIMES will be the
  // metallicity in solar units multiplied by ChimesDustRatio
  ChimesDepletionData[thread_id].ChimesDustRatio = 0.0;
  for (i = 0; i < DEPL_N_ELEM; i++)
    ChimesDepletionData[thread_id].ChimesDustRatio += ChimesDepletionData[thread_id].SolarAbund[i] * (1.0 - pow(10.0, chimes_depletion_linear_fit(nH, T, pars[i][0], pars[i][1], pars[i][2])));

  // The above sum gives the dust to gas mass ratio.
  // Now normalise it by the saturated value.
  ChimesDepletionData[thread_id].ChimesDustRatio /= ChimesDepletionData[thread_id].DustToGasSaturated;
}
#endif // CHIMES_METAL_DEPLETION


void chimes_gizmo_exit(void)
{
  endrun(56275362);
}
#endif // CHIMES
#endif
