#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/* this file handles the FIRE short-range radiation-pressure and
    photo-ionization terms. written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */



    
    
    
/* Routines for simple FIRE local photo-ionization heating feedback model. This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */

#ifdef CHIMES_HII_REGIONS 
/* This routine is based heavily on the HII_heating_singledomain() routine 
 * used in FIRE for HII heating. I have modified this to make use of the 
 * stellar luminosities used with the CHIMES routines, and it now only flags 
 * gas particles deemed to be within HII regions so that shielding in the CHIMES 
 * routines can be disabled for this particles. This routine does not actually 
 * heat and ionise these particles explicitly. */
void chimes_HII_regions_singledomain(void)
{
  if(All.Time<=0) 
    return;

  MyDouble *pos;
  int startnode, numngb, j, n, i, k;
  int do_ionize,dummy, n_iter_HII, age_bin;
  MyFloat h_i, dt, rho;
  double dx, dy, dz, r2, r, eps_cgs, prandom;
  double mionizable, mionized, RHII, RHIImax, RHIImin, R_search;
  double stellum, stellum_G0, prob, M_ionizing_emitted;
  double m_available, m_effective, RHIImultiplier;
  double stellar_age, stellar_mass, log_age_Myr;
  
  int max_n_iterations_HII = 5; 

  Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
    
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if((P[i].Type == 4) || ((All.ComovingIntegrationOn==0) && ((P[i].Type == 2) || (P[i].Type==3))))
	{
#ifndef WAKEUP
	  dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
	  dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
	  if(dt<=0) 
	    continue; // don't keep going with this loop

	  stellar_age = evaluate_stellar_age_Gyr(P[i].StellarAge); 
	  stellar_mass = P[i].Mass * All.UnitMass_in_g / (All.HubbleParam  * SOLAR_MASS); 
	  
	  // stellum is the number of H-ionising photons per second 
	  // produced by the star particle 
	  stellum = chimes_ion_luminosity(stellar_age * 1000.0, stellar_mass); 
	  if(stellum <= 0) 
	    continue;
	  
	  // Luminosity in the 6-13.6 eV band. 
	  stellum_G0 = chimes_G0_luminosity(stellar_age * 1000.0, stellar_mass); 

	  // Gravitational Softening (cgs units) 
	  eps_cgs = All.SofteningTable[P[i].Type] * All.cf_atime * All.UnitLength_in_cm / All.HubbleParam; 
	  
	  // Determine stellar age bin 
	  log_age_Myr = log10(stellar_age * 1000.0); 	  
	  if (log_age_Myr < CHIMES_LOCAL_UV_AGE_LOW) 
	    age_bin = 0; 
	  else if (log_age_Myr < CHIMES_LOCAL_UV_AGE_MID) 
	    age_bin = (int) floor(((log_age_Myr - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW) + 1); 
	  else 
	    { 
	      age_bin = (int) floor((((log_age_Myr - CHIMES_LOCAL_UV_AGE_MID) / CHIMES_LOCAL_UV_DELTA_AGE_HI) + ((CHIMES_LOCAL_UV_AGE_MID - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW)) + 1); 
	      if (age_bin > CHIMES_LOCAL_UV_NBINS - 1) 
		age_bin = CHIMES_LOCAL_UV_NBINS - 1; 
	    }
	  
	  pos = P[i].Pos;
	  rho = P[i].DensAroundStar;
	  h_i = PPP[i].Hsml;
	  
	  // Stromgren radius, RHII, computed using a case B recombination coefficient 
	  // at 10^4 K of 2.59e-13 cm^3 s^-1, as used in CHIMES, and assuming a 
	  // Hydrogen mass fraction XH = 0.7. 
	  RHII = 1.7376e-12 * pow(stellum, 0.33333) * pow(rho * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam, -0.66667);
	  
	  // Convert RHII from cm to code units 
	  RHII /= All.cf_atime*All.UnitLength_in_cm/All.HubbleParam;
	  
	  /* Impose a maximum RHII, to prevent the code trying to search 
	   * for neighbours too far away. Unlike the standard FIRE routines, 
	   * I do not base this on an estimate for where the flux falls below 
	   * the cosmic background. Instead, note that, for the maximum ionising 
	   * flux per Msol that we get from the Starburst99 models (which occurs 
	   * at a stellar age of 3.71 Myr), the ratio of ionisable gas mass to 
	   * stellar mass is 286 / nH. In other words, at nH = 1 cm^-3, a single 
	   * star particle can ionise 286 gas particles (assuming equal-mass 
	   * particles). The star particle's smoothing length h_i should contain
	   * DesNumNgb gas particles (typically 32). So if we set RHIImax to 
	   * 10 * h_i, this should be enough to handle HII regions down to 
	   * nH ~ 1 cm^-3. */ 
	  RHIImax = 10.0 * h_i; 
	  RHIImin = 0.5 * h_i; 
	  
	  // Ionizable gas mass in code units, based on the gas density 
	  // evaluated at the position of the star. Prefactor is 4pi/3. 
	  mionizable = 4.18879 * rho * pow(RHII, 3.0);  

	  // number of ionizing photons times proton mass, gives max mass ionized 
	  M_ionizing_emitted = PROTONMASS * stellum * (dt * All.UnitTime_in_s / All.HubbleParam);  // g
	  mionizable = DMIN(mionizable , M_ionizing_emitted/(All.UnitMass_in_g/All.HubbleParam)); 
	  
	  // Now limit RHII to be between the min and max defined above. 
	  if(RHII > RHIImax) 
	    RHII = RHIImax;

	  if(RHII < RHIImin) 
	    RHII = RHIImin;

	  /* Skip star particles that can ionise <10% of its own mass (this is  
	   * lower than 50% here, because there can be some variation between 
	   * particle masses, and in gas densities). */ 
	  if(mionizable / P[i].Mass > 0.1) 
	    {	      
	      prandom = get_random_number(P[i].ID + 7); 
	      mionized = 0.0;
	      startnode = All.MaxPart;     /* root node */
	      dummy = 0; 
	      n_iter_HII = 0;
	     
	      do {
		R_search = RHII;
		if(h_i > R_search) 
		  R_search = h_i;
		numngb = ngb_treefind_variable_threads(pos, R_search, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist);
		if(numngb>0)
		  {
		    for(n = 0; n < numngb; n++)
		      {
			j = Ngblist[n];
			if(P[j].Type == 0 && P[j].Mass > 0)
			  {
			    dx = pos[0] - P[j].Pos[0];
			    dy = pos[1] - P[j].Pos[1];
			    dz = pos[2] - P[j].Pos[2];
			    NEAREST_XYZ(dx, dy, dz, 1); /*  now find the closest image in the given box size  */
			    r2 = dx * dx + dy * dy + dz * dz;
			    r = sqrt(r2);
			   
			    /* If inside RHII and mionized<mionizeable and not already ionized, can be ionized! */
			    do_ionize=0; 
			    if((r <= RHII) && (SphP[j].DelayTimeHII <= 0) && (mionized < mionizable)) 
			      {
				m_effective = P[j].Mass * (SphP[j].Density / rho);
				// weight by density b/c of how the recomination rate in each particle scales 

				m_available = mionizable - mionized;
				if(m_effective <= m_available) 
				  {
				    // Enough photons to ionise the whole particle. 
				    do_ionize = 1;
				    mionized += m_effective; 
				  }
				else 
				  {
				    // Not enough to ionise a whole particle. 
				    // Use random number to determine whether 
				    // to ionise. 
				    prob = m_available/m_effective; 
				   
				    if(prandom < prob) 
				      do_ionize = 1;

				    mionized += prob * m_effective; 
				  } // if(m_effective<=m_available) 
			       
				if(do_ionize==1) 
				  {
				    SphP[j].DelayTimeHII = dt;
				   
				    for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++) 
				      {
					SphP[j].Chimes_fluxPhotIon_HII[k] = 0.0; 
					SphP[j].Chimes_G0_HII[k] = 0.0; 
				      }
				    
				    SphP[j].Chimes_fluxPhotIon_HII[age_bin] = (1.0 - All.Chimes_f_esc_ion) * stellum / (pow(r * All.cf_atime * All.UnitLength_in_cm / All.HubbleParam, 2.0) + pow(eps_cgs, 2.0)) ; 
				    SphP[j].Chimes_G0_HII[age_bin] = (1.0 - All.Chimes_f_esc_G0) * stellum_G0 / (pow(r * All.cf_atime * All.UnitLength_in_cm / All.HubbleParam, 2.0) + pow(eps_cgs, 2.0)); 
				  }
			      } // if((r <= RHII) && (SphP[j].DelayTimeHII <= 0) && (mionized<mionizable)) 
			  } // if(P[j].Type == 0 && P[j].Mass > 0)
		      } // for(n = 0; n < numngb; n++)
		  } // if(numngb>0)

		/* now check if we have ionized sufficient material, and if not, 
		   iterate with larger regions until we do */
		RHIImultiplier=1.10;
		if(mionized < 0.95 * mionizable) 
		  {
		    /* ok, this guy did not find enough gas to ionize, it needs to expand its search */
		    if((RHII >= RHIImax) || (n_iter_HII >= max_n_iterations_HII))
		      {
			/* we're done looping, this is just too big an HII region */
			mionized = 1.001*mionizable;
		      } 
		    else 
		      {
			/* in this case we're allowed to keep expanding RHII */
			if(mionized <= 0) 
			  RHIImultiplier = 2.0;
			else 
			  {
			    RHIImultiplier = pow(mionized / mionizable, -0.333);
			    if(RHIImultiplier > 5.0) 
			      RHIImultiplier=5.0;
			    if(RHIImultiplier < 1.26) 
			      RHIImultiplier=1.26;
			  } // if(mionized <= 0) 
		       
			RHII *= RHIImultiplier;
			if(RHII > 1.26*RHIImax) 
			  RHII=1.26*RHIImax;

			startnode=All.MaxPart; // this will trigger the while loop to continue
		      } // if((RHII>=RHIImax) || (n_iter_HII >= max_n_iterations_HII))
		  } // if(mionized < 0.95*mionizable) 
		n_iter_HII++;
	      } while(startnode >= 0);
	    } // if(mionizable / P[i].Mass > 0.1)
	} // if((P[i].Type == 4)||(P[i].Type == 2)||(P[i].Type == 3))
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  myfree(Ngblist);
} 
#endif // CHIMES_HII_REGIONS 
