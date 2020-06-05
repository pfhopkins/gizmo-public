#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

/*! \file rt_utilities.c
 *  \brief useful functions for radiation modules
 *
 *  This file contains a variety of useful functions having to do with radiation in different modules
 *    A number of the radiative transfer subroutines and more general mass-to-light ratio calculations
 *    will refer to these routines.
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */



/***********************************************************************************************************
 *
 * ROUTINES IN THIS BLOCK MUST BE MODIFIED FOR NEW MODULES USING DIFFERENT WAVEBANDS/PHYSICS
 *
 *  (these routines depend on compiler-time choices for which frequencies will be followed, and the 
 *    physics used to determine things like the types of source particles, source luminosities, 
 *    and how opacities are calculated)
 *
 ***********************************************************************************************************/



#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)

/***********************************************************************************************************/
/* routine which returns the luminosity [total volume/mass integrated] for the desired source particles in physical code units (energy/time),
    as a function of whatever the user desires, in the relevant bands. inpute here:
    'i' = index of target particle/cell for which the luminosity should be computed
    'mode' = flag for special behaviors. if <0 (e.g. -1), just returns whether or not a particle is 'active' (eligible as an RT source). if =0, normal behavior. if =1, then some bands have special behavior, for example self-shielding estimated -at the source-
    'lum' = pointer to vector of length N_RT_FREQ_BINS to hold luminosities for all bands
 */
/***********************************************************************************************************/
#ifdef CHIMES_STELLAR_FLUXES  
int rt_get_source_luminosity(int i, int mode, double *lum, double *chimes_lum_G0, double *chimes_lum_ion)
#else 
int rt_get_source_luminosity(int i, int mode, double *lum)
#endif 
{
    int active_check = 0;
    


    

    
#if defined(RT_INFRARED) /* can add direct infrared sources, but default to no direct IR (just re-emitted light) */
    if((1 << P[i].Type) & (RT_SOURCES))
    {
        lum[RT_FREQ_BIN_INFRARED] = 0.0; //default to no direct IR (just re-emitted light)
    }
#endif

    
#if defined(RT_NUV)
    /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models */
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        if( ((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3)))) && P[i].Mass>0 && PPP[i].Hsml>0 )
        {
            if(mode<0) {return 1;} active_check = 1;
            double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), f_op=0;
            if(star_age <= 0.0025) {f_op=0.09;} else {
                if(star_age <= 0.006) {f_op=0.09*(1+((star_age-0.0025)/0.004)*((star_age-0.0025)/0.004));
                } else {f_op=1-0.8410937/(1+sqrt((star_age-0.006)/0.3));}}
            double fac = P[i].Mass * UNIT_MASS_IN_SOLAR / UNIT_LUM_IN_SOLAR; // converts to code units
            lum[RT_FREQ_BIN_NUV] = (1-f_op) * fac * evaluate_light_to_mass_ratio(star_age, i);
        }
    }
#endif

    
#if defined(RT_OPTICAL_NIR)
    /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models */
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        if( ((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3)))) && P[i].Mass>0 && PPP[i].Hsml>0 )
        {
            if(mode<0) {return 1;} active_check = 1;
            double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), f_op=0;
            if(star_age <= 0.0025) {f_op=0.09;} else {
                if(star_age <= 0.006) {f_op=0.09*(1+((star_age-0.0025)/0.004)*((star_age-0.0025)/0.004));
                } else {f_op=1-0.8410937/(1+sqrt((star_age-0.006)/0.3));}}
            double fac = P[i].Mass * UNIT_MASS_IN_SOLAR / UNIT_LUM_IN_SOLAR; // converts to code units
            lum[RT_FREQ_BIN_OPTICAL_NIR] = f_op * fac * evaluate_light_to_mass_ratio(star_age, i);
        }
    }
#endif

    
#ifdef RT_PHOTOELECTRIC
    /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        if( ((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3)))) && P[i].Mass>0 && PPP[i].Hsml>0 )
        {
            if(mode<0) {return 1;} active_check = 1;
            double fac = (P[i].Mass * UNIT_MASS_IN_SOLAR) / UNIT_LUM_IN_CGS; // converts to code units
            //double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge); 
            //double l_band = 2.14e36 / sqrt(1. + pow(star_age/4.e-3,3.6)) * fac; // solar tracks, no nebular
            double l_band, x_age = evaluate_stellar_age_Gyr(P[i].StellarAge) / 3.4e-3;
            if(x_age <= 1) 
            { 
                l_band = 1.07e36 * (1.+x_age*x_age) * fac;
            } else {
                l_band = 2.14e36 / (x_age * sqrt(x_age)) * fac;
            } // 0.1 solar, with nebular. very weak metallicity dependence, with slightly slower decay in time for lower-metallicity pops; effect smaller than binaries
            lum[RT_FREQ_BIN_PHOTOELECTRIC] = l_band; // band luminosity //
        }
    }
#endif
    

#ifdef RT_LYMAN_WERNER
    /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        if( ((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3)))) && P[i].Mass>0 && PPP[i].Hsml>0 )
        {
            if(mode<0) {return 1;} active_check = 1;
            double fac = (P[i].Mass * UNIT_MASS_IN_SOLAR) / UNIT_LUM_IN_CGS; // converts to code units
            double l_band, x_age = evaluate_stellar_age_Gyr(P[i].StellarAge) / 3.4e-3;
            if(x_age <= 1) 
            { 
                l_band = 0.429e36 * (1.+x_age*x_age) * fac;
            } else {
                l_band = 0.962e36 * pow(x_age,-1.6) * exp(-x_age/117.6) * fac;
            } // 0.1 solar, with nebular. very weak metallicity dependence, with slightly slower decay in time for lower-metallicity pops; effect smaller than binaries
            lum[RT_FREQ_BIN_LYMAN_WERNER] = l_band; // band luminosity //
        }
    }
#endif
        
    
#if defined(RT_CHEM_PHOTOION)
    /* Hydrogen and Helium ionizing bands */
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        lum[RT_FREQ_BIN_H0] = 0; // begin from zero //
        double fac = 0;
#if defined(GALSF) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        /* calculate ionizing flux based on actual stellar or BH physics */
        if( ((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3)))) && P[i].Mass>0 && PPP[i].Hsml>0 )
        {
            if(mode<0) {return 1;} active_check=1;
            fac += particle_ionizing_luminosity_in_cgs(i) / UNIT_LUM_IN_CGS;
        }
#else
#ifdef RT_ILIEV_TEST1
        if(P[i].Type==4) {if(mode<0) {return 1;} active_check=1; fac += 5.0e48 * (13.6*ELECTRONVOLT_IN_ERGS) / UNIT_LUM_IN_CGS;} // 5e48 in ionizing photons per second //
#else
        if(P[i].Type==4) {if(mode<0) {return 1;} active_check=1; fac += All.IonizingLuminosityPerSolarMass_cgs * (P[i].Mass * UNIT_MASS_IN_SOLAR) / UNIT_LUM_IN_CGS;} // flux from star particles according to mass
#endif
#endif // GALSF else
#if defined(RT_PHOTOION_MULTIFREQUENCY)
        // we should have pre-tabulated how much luminosity gets assigned to each different waveband according to the following function //
        lum[RT_FREQ_BIN_He0]=lum[RT_FREQ_BIN_He1]=lum[RT_FREQ_BIN_He2]=0;
        int k; for(k=0;k<N_RT_FREQ_BINS;k++) {lum[k] += fac * precalc_stellar_luminosity_fraction[k];}
#else
        lum[RT_FREQ_BIN_H0] += fac;
#endif
    }
#endif // RT_CHEM_PHOTOION


#if defined(RT_HARD_XRAY) || defined(RT_SOFT_XRAY)
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
#if defined(RT_HARD_XRAY) 
            lum[RT_FREQ_BIN_HARD_XRAY] = 0; // LMXBs+HMXBs
#endif
#if defined(RT_SOFT_XRAY) 
            lum[RT_FREQ_BIN_SOFT_XRAY] = 0; // LMXBs+HMXBs
#endif
#if defined(BLACK_HOLES)
        if(P[i].Type == 5) 
        {
            if(mode<0) {return 1;} active_check=1;
            double lbol = bh_lum_bol(P[i].BH_Mdot,P[i].Mass,i); // luminosity in physical code units // 
            double lbol_lsun = lbol * UNIT_LUM_IN_SOLAR;
            double bol_corr = 0;
#if defined(RT_HARD_XRAY) 
            bol_corr = 0.43 * (10.83 * pow(lbol_lsun/1.e10,0.28) + 6.08 * pow(lbol_lsun/1.e10,-0.02)); // 0.5 for -ALL- hard-x-ray, 1.0 prefactor for just 2-10 keV
            lum[RT_FREQ_BIN_HARD_XRAY] = lbol / bol_corr; // typical bolometric correction from Hopkins, Richards, & Hernquist 2007 
#endif
#if defined(RT_SOFT_XRAY) 
            bol_corr = 17.87 * pow(lbol_lsun/1.e10,0.28) + 10.0 * pow(lbol_lsun/1.e10,-0.02);
            lum[RT_FREQ_BIN_SOFT_XRAY] = lbol / bol_corr; // typical bolometric correction from Hopkins, Richards, & Hernquist 2007 
#endif
        }
#endif
        if(P[i].Type == 4) 
        {
            if(mode<0) {return 1;} active_check=1;
            double fac = (P[i].Mass * UNIT_MASS_IN_SOLAR) / UNIT_LUM_IN_CGS; // converts to code units
            double L_HMXBs = 0.0; 
#ifdef GALSF
            double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
            if(star_age > 0.01) {L_HMXBs = 1.0e29 / (star_age*star_age);}
#endif
#if defined(RT_HARD_XRAY) 
            lum[RT_FREQ_BIN_HARD_XRAY] = (6.3e27 + 0.6*L_HMXBs) * fac; // LMXBs+HMXBs
#endif
#if defined(RT_SOFT_XRAY) 
            lum[RT_FREQ_BIN_SOFT_XRAY] = (8.2e27 + 0.4*L_HMXBs) * fac; // LMXBs+HMXBs
#endif
        }
    }
#endif // RT_HARD_XRAY

    
#if defined(RT_GENERIC_USER_FREQ)
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        if(P[i].Type == 4)
        {
            if(mode<0) {return 1;} active_check=1;
            lum[RT_FREQ_BIN_GENERIC_USER_FREQ] = 0;
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION /* assume special units for this problem, and that total mass of 'sources' is 1 */
            lum[RT_FREQ_BIN_GENERIC_USER_FREQ] = P[i].Mass * All.Vertical_Grain_Accel * C_LIGHT_CODE / (0.75*GRAIN_RDI_TESTPROBLEM_Q_AT_GRAIN_MAX/All.Grain_Size_Max); // special behavior for particular test of stratified boxes compared to explicit dust opacities
#endif
        }
    }
#endif
    
    
#ifdef RADTRANSFER
    /* generic sub routines for gas as a source term */
    if((1 << P[i].Type) & (RT_SOURCES)) // check if the particle falls into one of the allowed source types
    {
        if(P[i].Type == 0)
        {
            if(mode<0) {return 1;} active_check=1; // active //
            rt_get_lum_gas(i,lum); /* optionally re-distributes cooling flux as a blackbody */
            int k; for(k=0;k<N_RT_FREQ_BINS;k++) {lum[k] += SphP[i].Rad_Je[k];}        
        }
    }
#endif
    
    /* need to renormalize ALL sources for reduced speed of light */
    {int k; for(k=0;k<N_RT_FREQ_BINS;k++) {lum[k] *= RT_SPEEDOFLIGHT_REDUCTION;}}
    return active_check;
}





/***********************************************************************************************************/
/* calculate the opacity for use in radiation transport operations [in physical code units = Length^2/Mass]. this should
    be a total extinction opacity, i.e. kappa = kappa_scattering + kappa_absorption */
/***********************************************************************************************************/
double rt_kappa(int i, int k_freq)
{

#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION /* special test problem implementation */
    return 1*SphP[i].Interpolated_Opacity[k_freq] + 0.001 * All.Dust_to_Gas_Mass_Ratio * 0.75*GRAIN_RDI_TESTPROBLEM_Q_AT_GRAIN_MAX/All.Grain_Size_Max; /* enforce minimum */
#endif
    return MIN_REAL_NUMBER + SphP[i].Interpolated_Opacity[k_freq]; /* this is calculated in a different routine, just return it now */
#endif

#ifdef RT_CHEM_PHOTOION
    /* opacity to ionizing radiation for Petkova & Springel bands. note rt_update_chemistry is where ionization is actually calculated */
    double nH_over_Density = HYDROGEN_MASSFRAC / PROTONMASS * UNIT_MASS_IN_CGS;
    double kappa = nH_over_Density * (SphP[i].HI + MIN_REAL_NUMBER) * rt_sigma_HI[k_freq];
#if defined(RT_CHEM_PHOTOION_HE) && defined(RT_PHOTOION_MULTIFREQUENCY)
    kappa += nH_over_Density * ((SphP[i].HeI + MIN_REAL_NUMBER) * rt_sigma_HeI[k_freq] + (SphP[i].HeII + MIN_REAL_NUMBER) * rt_sigma_HeII[k_freq]);
    if(k_freq==RT_FREQ_BIN_He0)  {return kappa;}
    if(k_freq==RT_FREQ_BIN_He1)  {return kappa;}
    if(k_freq==RT_FREQ_BIN_He2)  {return kappa;}
#endif
    if(k_freq==RT_FREQ_BIN_H0)  {return kappa;}
#endif

#if defined(RT_HARD_XRAY) || defined(RT_SOFT_XRAY) || defined(RT_PHOTOELECTRIC) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(RT_NUV) || defined(RT_OPTICAL_NIR) || defined(RT_LYMAN_WERNER) || defined(RT_INFRARED) || defined(RT_FREEFREE)
    double fac = UNIT_SURFDEN_IN_CGS, Zfac; /* units */
    Zfac = 1.0; // assume solar metallicity 
#ifdef METALS
    Zfac = P[i].Metallicity[0]/All.SolarAbundances[0];
#endif
#ifdef RT_FREEFREE /* pure (grey, non-relativistic) Thompson scattering opacity + free-free absorption opacity */
    if(k_freq==RT_FREQ_BIN_FREEFREE)
    {
        double T_eff=0.59*(GAMMA(i)-1.)*U_TO_TEMP_UNITS*SphP[i].InternalEnergyPred, rho=SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS; // we're assuming fully-ionized gas with a simple equation-of-state here, nothing fancy, to get the temperature //
        double kappa_abs = 1.e30*rho*pow(T_eff,-3.5);
        return (0.35 + kappa_abs) * fac;
    }
#endif
#ifdef RT_HARD_XRAY
    /* opacity comes from H+He (Thompson) + metal ions */
    if(k_freq==RT_FREQ_BIN_HARD_XRAY) {return (0.53 + 0.27*Zfac) * fac;}
#endif
#ifdef RT_SOFT_XRAY
    /* opacity comes from H+He (Thompson) + metal ions */
    if(k_freq==RT_FREQ_BIN_SOFT_XRAY) {return (127. + 50.0*Zfac) * fac;}
#endif
#ifdef RT_PHOTOELECTRIC
    /* opacity comes primarily from dust (ignoring H2 molecular opacities here) */
    if(k_freq==RT_FREQ_BIN_PHOTOELECTRIC) {return 2000. * DMAX(1.e-4,Zfac) * fac;}
#endif
#ifdef RT_LYMAN_WERNER
    /* opacity from molecular H2 and dust (dominant at higher-metallicity) should be included */
    if(k_freq==RT_FREQ_BIN_LYMAN_WERNER) {return 2400.*Zfac * fac;} // just dust term for now
#endif
#ifdef RT_NUV
    /* opacity comes primarily from dust */
    if(k_freq==RT_FREQ_BIN_NUV) {return 1800.*Zfac * fac;}
#endif
#ifdef RT_OPTICAL_NIR
    /* opacity comes primarily from dust */
    if(k_freq==RT_FREQ_BIN_OPTICAL_NIR) {return 180.*Zfac * fac;}
#endif
#ifdef RT_INFRARED
    /* IR with dust opacity */
    double T_min = get_min_allowed_dustIRrad_temperature();
    if(k_freq==RT_FREQ_BIN_INFRARED)
    {
        if(isnan(SphP[i].Dust_Temperature) || SphP[i].Dust_Temperature<=T_min) {SphP[i].Dust_Temperature=T_min;} // reset baseline
        if(isnan(SphP[i].Radiation_Temperature) || SphP[i].Radiation_Temperature<=T_min) {SphP[i].Radiation_Temperature=T_min;} // reset baseline
        
        double T_dust_em = SphP[i].Dust_Temperature; // dust temperature in K //
        double Trad = SphP[i].Radiation_Temperature; // radiation temperature in K //
        if(Trad <= 0) {Trad = 5600.;}
        double kappa = 0.0; 
        /* opacities are from tables of Semenov et al 2003; we use their 'standard' 
            model, for each -dust- temperature range (which gives a different dust composition, 
            hence different wavelength-dependent specific opacity). We then integrate to 
            get the Rosseland-mean opacity for the given dust composition, assuming 
            the radiation is a blackbody with the specified -radiation- temperature. 
            We adopt their 'porous 5-layered sphere' model for dust composition. 
            We use simple fitting functions to the full tabulated data: however, note that
            (because the blackbody assumption smoothes fine structure in the opacities), 
            the deviations from the fit functions are much smaller than the deviations owing 
            to different grain composition choices (porous/non, composite/non, 5-layer/aggregated/etc) 
            in Semenov et al's paper */
        if(T_dust_em < 1500.) // < 1500 K, dust is present
        {
            double x = 4.*log10(Trad) - 8.; // needed for fitting functions to opacities (may come up with cheaper function later)
            double dx_excess=0; if(x > 7.) {dx_excess=x-7.; x=7.;} // cap for maximum temperatures at which fit-functions should be used //
            if(x < -4.) {x=-4.;} // cap for minimum temperatures at which fit functions below should be used //
            if(T_dust_em < 160.) // Tdust < 160 K (all dust constituents present)
            {
                kappa = exp(0.72819004 + 0.75142468*x - 0.07225763*x*x - 0.01159257*x*x*x + 0.00249064*x*x*x*x);
            } else if(T_dust_em < 275.) { // 160 < Tdust < 275 (no ice present)
                kappa = exp(0.16658241 + 0.70072926*x - 0.04230367*x*x - 0.01133852*x*x*x + 0.0021335*x*x*x*x);
            } else if(T_dust_em < 425.) { // 275 < Tdust < 425 (no ice or volatile organics present)
                kappa = exp(0.03583845 + 0.68374146*x - 0.03791989*x*x - 0.01135789*x*x*x + 0.00212918*x*x*x*x);        
            } else if(T_dust_em < 680.) { // 425 < Tdust < 680 (silicates, iron, & troilite present)
                kappa = exp(-0.76576135 + 0.57053532*x - 0.0122809*x*x - 0.01037311*x*x*x + 0.00197672*x*x*x*x);
            } else { // 680 < Tdust < 1500 (silicates & iron present)
                kappa = exp(-2.23863222 + 0.81223269*x + 0.08010633*x*x + 0.00862152*x*x*x - 0.00271909*x*x*x*x);
            }
            if(dx_excess > 0) {kappa *= exp(0.57*dx_excess);} // assumes kappa scales linearly with temperature (1/lambda) above maximum in fit; pretty good approximation //
            kappa *= Zfac; // the above are all dust opacities, so they scale with metallicity
        }
#ifdef COOLING
        kappa += 0.35 * SphP[i].Ne; // Thompson scattering
#endif
        return kappa * fac; // convert units and return
    }
#endif
#endif
    
    
    return 0;
}





/***********************************************************************************************************/
/* calculate absorbed fraction of opacity = 1-albedo = kappa_absorption / (kappa_scattering + kappa_absorption) needed for RT operations */
/***********************************************************************************************************/
double rt_absorb_frac_albedo(int i, int k_freq)
{
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
    return DMAX(1.e-6,DMIN(1.0-1.e-6,(1.0*GRAIN_RDI_TESTPROBLEM_SET_ABSFRAC)));
#endif
    return 0.5; /* appropriate for single-scattering (e.g. ISM dust at optical wavelengths) */
    //return 1.-1.e-6; /* appropriate for multiple-scattering at far-IR (wavelength much longer than dust size) */
#endif

#ifdef RT_CHEM_PHOTOION
    if(k_freq==RT_FREQ_BIN_H0)  {return 1.-1.e-6;} /* negligible scattering for ionizing radiation */
#if defined(RT_CHEM_PHOTOION_HE) && defined(RT_PHOTOION_MULTIFREQUENCY)
    if(k_freq==RT_FREQ_BIN_He0 || k_freq==RT_FREQ_BIN_He1 || k_freq==RT_FREQ_BIN_He2)  {return 1.-1.e-6;}
#endif
#endif

#if defined(RT_HARD_XRAY) || defined(RT_SOFT_XRAY) || defined(RT_INFRARED) /* these have mixed opacities from dust(assume albedo=1/2), ionization(albedo=0), and Thompson (albedo=1) */
    double fac; fac = UNIT_SURFDEN_IN_CGS; /* units */
#ifdef RT_HARD_XRAY /* opacity comes from H+He (Thompson) + metal ions -- assume 0 scattering from ions, 1 from Thompson */
    if(k_freq==RT_FREQ_BIN_HARD_XRAY) {return 1.-0.5*(0. + DMIN(1.,0.35*fac/rt_kappa(i,k_freq)));}
#endif
#ifdef RT_SOFT_XRAY /* opacity comes from H+He (Thompson) + metal ions -- assume 0 scattering from ions, 1 from Thompson */
    if(k_freq==RT_FREQ_BIN_SOFT_XRAY) {return 1.-0.5*(0. + DMIN(1.,0.35*fac/rt_kappa(i,k_freq)));}
#endif
#ifdef RT_INFRARED /* opacity comes from Thompson + dust -- assume 0.5/(1 + (Trad/725K)^(-2)) scattering from dust [Rayleigh, since we're in the long-wavelength limit by definition here], 1 from Thompson */
    if(k_freq==RT_FREQ_BIN_INFRARED)
    {
        double fA_tmp = (1.-0.5/(1.+((725.*725.)/(1.+SphP[i].Radiation_Temperature*SphP[i].Radiation_Temperature))));
#ifdef COOLING
        fA_tmp *= (1.-DMIN(1.,0.35*SphP[i].Ne*fac/rt_kappa(i,k_freq)));
#endif
        return fA_tmp;
    }
#endif
#endif
    
#ifdef RT_FREEFREE
    if(k_freq==RT_FREQ_BIN_FREEFREE)
    {
        double T_eff=0.59*(GAMMA(i)-1.)*U_TO_TEMP_UNITS*SphP[i].InternalEnergyPred, rho=SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS, kappa_abs = 1.e30*rho*pow(T_eff,-3.5);
        return kappa_abs / (0.35 + kappa_abs);
    }
#endif
    
    return 0.5; /* default to assuming kappa_scattering = kappa_absorption (pretty reasonable for dust at most wavelengths) */
}



#endif // #if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)






/***********************************************************************************************************
 *
 * ROUTINES WHICH DO NOT NEED TO BE MODIFIED SHOULD GO BELOW THIS BREAK
 *
 *  (these routines may depend on the RT solver or other numerical choices, but below here, place routines 
 *    which shouldn't needed to be hard-coded for different assumptions about the bands of the RT module, etc)
 *
 ***********************************************************************************************************/



/***********************************************************************************************************/
/* rate of photon absorption [absorptions per unit time per photon]: this, times the timestep dt, times the photon energy density E,
    gives the change in the energy density from absorptions (the sink term) */
/***********************************************************************************************************/
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
double rt_absorption_rate(int i, int k_freq)
{
    /* should be equal to (c * Kappa_opacity * rho) */
    return C_LIGHT_CODE_REDUCED * rt_absorb_frac_albedo(i, k_freq) * rt_kappa(i, k_freq) * SphP[i].Density*All.cf_a3inv;
}
#endif 




#ifdef RADTRANSFER

/***********************************************************************************************************/
/* returns the photon diffusion coefficient = fluxlimiter * speed_of_light[reduced] / (kappa_opacity * density)  [physical units] */
/***********************************************************************************************************/
double rt_diffusion_coefficient(int i, int k_freq)
{
    return return_flux_limiter(i,k_freq) * C_LIGHT_CODE_REDUCED / (1.e-45 + SphP[i].Rad_Kappa[k_freq] * SphP[i].Density*All.cf_a3inv);
}



/***********************************************************************************************************/
/* calculate the eddington tensor according to the different closure option[s] adopted */
/***********************************************************************************************************/
void rt_eddington_update_calculation(int j)
{
#ifdef RT_OTVET
    return; /* eddington tensor is calculated elsewhere [in the gravity subroutine]: don't mess with it here! */
#endif
#ifdef RT_M1
    /* calculate the eddington tensor with the M1 closure */
    int k_freq, k; double c_light = C_LIGHT_CODE_REDUCED, n_flux_j[3], fmag_j, V_j_inv = SphP[j].Density / P[j].Mass;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        n_flux_j[0]=n_flux_j[1]=n_flux_j[2]=0;
        double flux_vol[3]; for(k=0;k<3;k++) {flux_vol[k] = SphP[j].Rad_Flux[k_freq][k] * V_j_inv;}
        fmag_j = 0; for(k=0;k<3;k++) {fmag_j += flux_vol[k]*flux_vol[k];}
        if(fmag_j <= 0) {fmag_j=0;} else {fmag_j=sqrt(fmag_j); for(k=0;k<3;k++) {n_flux_j[k]=flux_vol[k]/fmag_j;}}
        double f_chifac = fmag_j / (MIN_REAL_NUMBER + c_light * SphP[j].Rad_E_gamma[k_freq] * V_j_inv);
        if(f_chifac < 0) {f_chifac=0;}
        if(fmag_j <= 0) {f_chifac = 0;}
        // restrict values of f_chifac to physical range.
        double f_min = 0.01, f_max = 0.99;
        if((f_chifac < f_min) || (isnan(f_chifac))) {f_chifac = f_min;}
        if(f_chifac > f_max) {f_chifac = f_max;}
        double chi_j = (3.+4.*f_chifac*f_chifac) / (5. + 2.*sqrt(4. - 3.*f_chifac*f_chifac));
        double chifac_iso_j = 0.5 * (1.-chi_j);
        double chifac_n_j = 0.5 * (3.*chi_j-1.);
        for(k=0;k<6;k++)
        {
            SphP[j].ET[k_freq][k] = 0;
            if(k<3)
            {
                SphP[j].ET[k_freq][k] = chifac_iso_j + chifac_n_j * n_flux_j[k]*n_flux_j[k];
            } else {
                if(k==3) {SphP[j].ET[k_freq][k] = chifac_n_j * n_flux_j[0]*n_flux_j[1];} // recall, for ET: 0=xx,1=yy,2=zz,3=xy,4=yz,5=xz
                if(k==4) {SphP[j].ET[k_freq][k] = chifac_n_j * n_flux_j[1]*n_flux_j[2];}
                if(k==5) {SphP[j].ET[k_freq][k] = chifac_n_j * n_flux_j[0]*n_flux_j[2];}
            }
        }
    }
    return;
#endif
#ifdef RT_FLUXLIMITEDDIFFUSION
    /* always assume the isotropic eddington tensor */
    int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {SphP[j].ET[k_freq][0]=SphP[j].ET[k_freq][1]=SphP[j].ET[k_freq][2]=1./3.; SphP[j].ET[k_freq][3]=SphP[j].ET[k_freq][4]=SphP[j].ET[k_freq][5]=0;}
    return;
#endif
    
    /* if nothing is set, default to guess the isotropic eddington tensor */
    {int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {SphP[j].ET[k_freq][0]=SphP[j].ET[k_freq][1]=SphP[j].ET[k_freq][2]=1./3.; SphP[j].ET[k_freq][3]=SphP[j].ET[k_freq][4]=SphP[j].ET[k_freq][5]=0;}}
    return;
}


/***********************************************************************************************************/
/*! simple subroutine to compute the dot product of the symmetric Eddington tensor ET=D with a
    vector v=vec_in, returning u=vec_out as u=D.v. Here u[0,1,2]=u[x,y,z], and
    D[0]=xx,D[1]=yy,D[2]=zz,D[3]=xy,D[4]=yz,D[5]=xz components of ET following our convention */
/***********************************************************************************************************/
void eddington_tensor_dot_vector(double ET[6], double vec_in[3], double vec_out[3])
{
    vec_out[0] = vec_in[0]*ET[0] + vec_in[1]*ET[3] + vec_in[2]*ET[5];
    vec_out[1] = vec_in[0]*ET[3] + vec_in[1]*ET[1] + vec_in[2]*ET[4];
    vec_out[2] = vec_in[0]*ET[5] + vec_in[1]*ET[4] + vec_in[2]*ET[2];
    return;
}




/***********************************************************************************************************/
/*! return the value of the flux-limiter function, as needed */
/***********************************************************************************************************/
double return_flux_limiter(int target, int k_freq)
{
#ifdef RT_FLUXLIMITER
    return SphP[target].Rad_Flux_Limiter[k_freq]; // apply flux-limiter
#endif
    return 1;
}



/***********************************************************************************************************/
/*
  routine which does the drift/kick operations on radiation quantities. separated here because we use a non-trivial
    update to deal with potentially stiff absorption terms (could be done more rigorously with something fully implicit in this 
    step, in fact). 
    mode = 0 == kick operation (update the conserved quantities)
    mode = 1 == predict/drift operation (update the predicted quantities)
 */
/***********************************************************************************************************/
void rt_update_driftkick(int i, double dt_entr, int mode)
{
#if defined(RT_EVOLVE_ENERGY) || defined(RT_EVOLVE_INTENSITIES)
    int kf, k_tmp; double total_erad_emission_minus_absorption = 0;
#if defined(RT_EVOLVE_INTENSITIES)
    for(kf=0;kf<N_RT_FREQ_BINS;kf++) {SphP[i].Rad_E_gamma[kf]=0; for(k_tmp=0;k_tmp<N_RT_INTENSITY_BINS;k_tmp++) {SphP[i].Rad_E_gamma[kf]+=RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity[kf][k_tmp];}}
    double de_emission_minus_absorption_saved[N_RT_FREQ_BINS][N_RT_INTENSITY_BINS]; // save this for use below
#endif
#ifdef RT_INFRARED
    double E_abs_tot = 0;/* energy absorbed in other bands is transfered to IR, by default: track it here */
    double Rad_E_gamma_tot = 0; // dust temperature defined by total radiation energy density //
    {int j; for(j=0;j<N_RT_FREQ_BINS;j++) {Rad_E_gamma_tot += SphP[i].Rad_E_gamma[j];}}
    double u_gamma = Rad_E_gamma_tot * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_CGS; // photon energy density in CGS //
    double Dust_Temperature_4 = C_LIGHT_CODE_REDUCED*UNIT_VEL_IN_CGS * u_gamma / (4. * 5.67e-5); // estimated effective temperature of local rad field in equilibrium with dust emission //
    SphP[i].Dust_Temperature = sqrt(sqrt(Dust_Temperature_4));
    double T_min = get_min_allowed_dustIRrad_temperature();
    if(SphP[i].Dust_Temperature <= T_min) {SphP[i].Dust_Temperature = T_min;} // dust temperature shouldn't be below CMB
#endif
    for(k_tmp=0; k_tmp<N_RT_FREQ_BINS; k_tmp++)
    {
#ifdef RT_INFRARED
        // need to do IR last after summing absorption from other bands //
        if(RT_FREQ_BIN_INFRARED < N_RT_FREQ_BINS-1) {if(kf == RT_FREQ_BIN_INFRARED) {kf = N_RT_FREQ_BINS-1;} if(kf == N_RT_FREQ_BINS-1) {kf = RT_FREQ_BIN_INFRARED;}}
#endif
#if defined(RT_EVOLVE_INTENSITIES)
        int k_angle; for(k_angle=0;k_angle<N_RT_INTENSITY_BINS;k_angle++)
#endif
        {
            kf = k_tmp; // normal loop
            double e0, dt_e_gamma_band=0, total_de_dt=0, a0 = -rt_absorption_rate(i,kf);
#if defined(RT_EVOLVE_INTENSITIES)
            if(mode==0) {e0 = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity[kf][k_angle];} else {e0 = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity_Pred[kf][k_angle];}
            dt_e_gamma_band = RT_INTENSITY_BINS_DOMEGA*SphP[i].Dt_Rad_Intensity[kf][k_angle];
#else
            if(mode==0) {e0 = SphP[i].Rad_E_gamma[kf];} else {e0 = SphP[i].Rad_E_gamma_Pred[kf];}
            dt_e_gamma_band = SphP[i].Dt_Rad_E_gamma[kf];
#endif
            total_de_dt = SphP[i].Rad_Je[kf] + dt_e_gamma_band;

#ifdef RT_INFRARED
            if(kf == RT_FREQ_BIN_INFRARED)
            {
                if((mode==0) && (dt_e_gamma_band!=0) && (dt_entr>0)) // only update temperatures on kick operations //
                {
                    // advected radiation changes temperature of radiation field, before absorption //
                    double dE_fac = dt_e_gamma_band * dt_entr; // change in energy from advection
                    double dTE_fac = SphP[i].Dt_Rad_E_gamma_T_weighted_IR * dt_entr; // T-weighted change from advection
                    double dE_abs = -e0 * (1. - exp(a0*dt_entr)); // change in energy from absorption
                    double rfac=1; if(dE_fac < -0.5*(e0+dE_abs)) {rfac=fabs(0.5*(e0+dE_abs))/fabs(dE_fac);} else {if(dE_fac > 0.5*e0) {rfac=0.5*e0/dE_fac;}}
                    dE_fac*=rfac; dTE_fac*=rfac; // limit temperature change from advection to prevent spurious divergences
                    
                    double T_max = DMAX(SphP[i].Radiation_Temperature , dE_fac / dTE_fac);
                    SphP[i].Radiation_Temperature = (e0 + dE_fac) / (MIN_REAL_NUMBER + DMAX(0., e0 / SphP[i].Radiation_Temperature + dTE_fac));
                    SphP[i].Radiation_Temperature = DMIN(SphP[i].Radiation_Temperature, T_max);
                    a0 = -rt_absorption_rate(i,kf); // update absorption rate using the new radiation temperature //
                }
                double total_emission_rate = E_abs_tot + fabs(a0)*e0 + SphP[i].Rad_Je[kf]; // add the summed absorption as emissivity here //
                total_de_dt = E_abs_tot + SphP[i].Rad_Je[kf] + dt_e_gamma_band;
                if(fabs(a0)>0)
                {
                    Dust_Temperature_4 = total_emission_rate * (SphP[i].Density*All.cf_a3inv/P[i].Mass) / (4. * (MIN_REAL_NUMBER + fabs(a0)) / C_LIGHT_CODE_REDUCED); // flux units
                    Dust_Temperature_4 *= UNIT_FLUX_IN_CGS / (5.67e-5); // convert to cgs
                    SphP[i].Dust_Temperature = sqrt(sqrt(Dust_Temperature_4));
                    if(SphP[i].Dust_Temperature < T_min) {SphP[i].Dust_Temperature = T_min;} // dust temperature shouldn't be below CMB
                }
                if(mode==0) // only update temperatures on kick operations //
                {
                    // dust absorption and re-emission brings T_rad towards T_dust: //
                    double dE_abs = -e0 * (1. - exp(a0*dt_entr)); // change in energy from absorption
                    double T_max = DMAX(SphP[i].Radiation_Temperature , SphP[i].Dust_Temperature); // should not exceed either initial temperature //
                    SphP[i].Radiation_Temperature = (e0 + dE_abs + total_emission_rate*dt_entr) / (MIN_REAL_NUMBER + (e0 + dE_abs) / SphP[i].Radiation_Temperature + total_emission_rate*dt_entr / SphP[i].Dust_Temperature);
                    SphP[i].Radiation_Temperature = DMIN(SphP[i].Radiation_Temperature, T_max);
                }
                if(SphP[i].Radiation_Temperature < T_min) {SphP[i].Radiation_Temperature = T_min;} // radiation temperature shouldn't be below CMB
            }
#endif
            
            /*---------------------------------------------------------------------------------------------------
             the following block is for absorption and special behavior where
             photons absorbed in one band are re-radiated [or up/down-scattered] into another.
             this must be hard-coded to maintain conservation (as opposed to treated as a source term)
             -----------------------------------------------------------------------------------------------------*/
            if(fabs(a0)*dt_entr > 50.) {a0 *= 50./(fabs(a0)*dt_entr);}
            double abs_0 = DMAX(0,fabs(a0)*dt_entr); double slabfac = slab_averaging_function(abs_0); double e_abs_0=exp(-abs_0); if(abs_0>100.) {e_abs_0=0;}
            /* since we're taking exponentials and inverses of some large numbers here, need to be careful not to let floating point errors cause nan's */
            if((dt_entr <= 0.)||(a0 >= 0.)||(abs_0 <= 0.)) {abs_0=0.; slabfac=e_abs_0=1.;} else {if(abs_0 < 1.e-5) {slabfac=1.-0.5*abs_0; e_abs_0 = 1.-abs_0;} else {if(abs_0 > 100.) {slabfac = 1./abs_0; e_abs_0 = 0.;}}}
            double e0_postabs = e0*e_abs_0, de_postabs = total_de_dt * dt_entr * slabfac, f_min = 0.01;
            if(e0_postabs+de_postabs < f_min*e0_postabs) {slabfac *= fabs((1.-f_min)*e0_postabs)/(fabs(de_postabs)+MIN_REAL_NUMBER);}
            
            double ef = e0 * e_abs_0 + total_de_dt * dt_entr * slabfac; // gives exact solution for dE/dt = -E*abs + de , the 'reduction factor' appropriately suppresses the source term //
            if((ef < 0)||(isnan(ef))) {ef=0;}
            double de_abs = e0 + total_de_dt * dt_entr - ef; // energy removed by absorption alone
            double de_emission_minus_absorption = (ef - (e0 + dt_e_gamma_band * dt_entr)); // total change, relative to what we would get with just advection (positive = net energy increase in the gas)
            if((dt_entr <= 0) || (de_abs <= 0)) {de_abs = 0;}
            
#if defined(RT_RAD_PRESSURE_FORCES) && defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_EVOLVE_FLUX) && !defined(RT_RADPRESSURE_IN_HYDRO)
            // for OTVET/FLD methods, need to apply radiation pressure term here so can limit this b/c just based on a gradient which is not flux-limited [as in hydro operators] //
            {
                double radacc[3]={0}, rmag=0, vel_i[3], L_particle = Get_Particle_Size(i)*All.cf_atime; // particle effective size/slab thickness
                double Sigma_particle = P[i].Mass / (M_PI*L_particle*L_particle); // effective surface density through particle
                double abs_per_kappa_dt = C_LIGHT_CODE_REDUCED * (SphP[i].Density*All.cf_a3inv) * dt_entr; // fractional absorption over timestep
                double f_kappa_abs = rt_absorb_frac_albedo(i,kf); // get albedo, we'll need this below
                double slabfac_rp = slab_averaging_function(f_kappa_abs*SphP[i].Rad_Kappa[kf]*Sigma_particle) * slab_averaging_function(f_kappa_abs*SphP[i].Rad_Kappa[kf]*abs_per_kappa_dt); // reduction factor for absorption over dt
                int kx; for(kx=0;kx<3;kx++)
                {
                    radacc[kx] = -dt_entr * slabfac_rp * return_flux_limiter(i,kf) * (SphP[i].Gradients.Rad_E_gamma_ET[kf][kx] / SphP[i].Density) / All.cf_atime; // naive radiation-pressure calc for FLD methods [physical units]
                    rmag += radacc[kx]*radacc[kx]; // compute magnitude
                    if(mode==0) {vel_i[kx]=P[i].Vel[kx]/All.cf_atime;} else {vel_i[kx]=SphP[i].VelPred[kx]/All.cf_atime;}
                }
                if(rmag > 0)
                {
                    rmag = sqrt(rmag); for(kx=0;kx<3;kx++) {radacc[kx] /= rmag;} // normalize
                    double rmag_max = de_abs / (P[i].Mass * C_LIGHT_CODE_REDUCED * (MIN_REAL_NUMBER + f_kappa_abs)); // limit magnitude to absorbed photon momentum
                    if(rmag > rmag_max) {rmag=rmag_max;}
#if defined(RT_ENABLE_R15_GRADIENTFIX)
                    rmag = rmag_max; // set to maximum (optically thin limit)
#endif
                    double work_band = 0;
                    for(kx=0;kx<3;kx++)
                    {
                        double radacc_eff = radacc[kx] * rmag; // re-normalize according to the criterion above
                        work_band += vel_i[kx] * radacc_eff * P[i].Mass; // PdV work done by photons [absorbed ones are fully-destroyed, so their loss of energy and momentum is already accounted for by their deletion in this limit //
                        if(mode==0) {P[i].Vel[kx] += radacc_eff * All.cf_atime;} else {SphP[i].VelPred[kx] += radacc_eff * All.cf_atime;}
                    }
                    double d_egy_rad = (2.*f_kappa_abs-1.)*work_band , d_egy_int = -2.*f_kappa_abs*work_band;
                    if(mode==0) {SphP[i].InternalEnergy += d_egy_int;} else {SphP[i].InternalEnergyPred += d_egy_int;}
#if defined(RT_EVOLVE_INTENSITIES)
                    {int k_q; for(k_q=0;k_q<N_RT_INTENSITY_BINS;k_q++) {if(mode==0) {SphP[i].Rad_Intensity[kf][k_q]+=d_egy_rad/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[kf][k_q]+=d_egy_rad/RT_INTENSITY_BINS_DOMEGA;}}}
#else
                    if(mode==0) {SphP[i].Rad_E_gamma[kf]+=d_egy_rad;} else {SphP[i].Rad_E_gamma_Pred[kf]+=d_egy_rad;}
#endif
                }
            }
#endif
            
            int donation_target_bin = -1; // frequency into which the photons will be deposited, if any //
#if defined(RT_CHEM_PHOTOION) && defined(RT_OPTICAL_NIR)
            /* assume absorbed ionizing photons are re-emitted via recombination into optical-NIR bins. valid if recombination time is fast.
             more accurately, this should be separately calculated in the cooling rates, and gas treated as a source */
            if(kf==RT_FREQ_BIN_H0) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
#ifdef RT_PHOTOION_MULTIFREQUENCY
            if(kf==RT_FREQ_BIN_He0) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
            if(kf==RT_FREQ_BIN_He1) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
            if(kf==RT_FREQ_BIN_He2) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
#endif
#endif
#if defined(RT_PHOTOELECTRIC) && defined(RT_INFRARED)
            if(kf==RT_FREQ_BIN_PHOTOELECTRIC) {donation_target_bin=RT_FREQ_BIN_INFRARED;} /* this is direct dust absorption, re-radiated in IR */
#endif
#if defined(RT_NUV) && defined(RT_INFRARED)
            if(kf==RT_FREQ_BIN_NUV) {donation_target_bin=RT_FREQ_BIN_INFRARED;} /* this is direct dust absorption, re-radiated in IR */
#endif
#if defined(RT_OPTICAL_NIR) && defined(RT_INFRARED)
            if(kf==RT_FREQ_BIN_OPTICAL_NIR) {donation_target_bin=RT_FREQ_BIN_INFRARED;} /* this is direct dust absorption, re-radiated in IR */
#endif
#ifdef RT_INFRARED
            if(donation_target_bin==RT_FREQ_BIN_INFRARED) {E_abs_tot += de_abs/(MIN_REAL_NUMBER + dt_entr);} /* donor bin is yourself in the IR - all self-absorption is re-emitted */
            if(kf==RT_FREQ_BIN_INFRARED) {ef = e0 + total_de_dt * dt_entr;} /* donor bin is yourself in the IR - all self-absorption is re-emitted */
#endif
            // isotropically re-emit the donated radiation into the target bin[s] //
#if defined(RT_EVOLVE_INTENSITIES)
            if(donation_target_bin >= 0) {int k_q; for(k_q=0;k_q<N_RT_INTENSITY_BINS;k_q++) {if(mode==0) {SphP[i].Rad_Intensity[donation_target_bin][k_q] += de_abs/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[donation_target_bin][k_q] += de_abs/RT_INTENSITY_BINS_DOMEGA;}}}
            if((ef < 0)||(isnan(ef))) {ef=0;}
            if(mode==0) {SphP[i].Rad_Intensity[kf][k_angle] = ef/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[kf][k_angle] = ef/RT_INTENSITY_BINS_DOMEGA;}
#else
            if(donation_target_bin >= 0) {if(mode==0) {SphP[i].Rad_E_gamma[donation_target_bin] += de_abs;} else {SphP[i].Rad_E_gamma_Pred[donation_target_bin] += de_abs;}}
            if((ef < 0)||(isnan(ef))) {ef=0;}
            if(mode==0) {SphP[i].Rad_E_gamma[kf] = ef;} else {SphP[i].Rad_E_gamma_Pred[kf] = ef;}
#endif

#if defined(RT_EVOLVE_FLUX)
            int k_dir; double f_mag=0, E_rad_forflux=0, vdot_h[3]={0}, vel_i[3]={0}, DeltaFluxEff[3]={0}, rho=SphP[i].Density*All.cf_a3inv; E_rad_forflux=0.5*(e0+ef); // use energy density averaged over this update for the operation below
            for(k_dir=0;k_dir<3;k_dir++) {if(mode==0) {vel_i[k_dir]=P[i].Vel[k_dir]/All.cf_atime;} else {vel_i[k_dir]=SphP[i].VelPred[k_dir]/All.cf_atime;}} // need gas velocity at this time
            double teqm_inv = SphP[i].Rad_Kappa[kf] * rho * C_LIGHT_CODE_REDUCED + MIN_REAL_NUMBER; // physical code units of 1/time, defines characteristic timescale for coming to equilibrium flux. see notes for CR second-order module for details. //
            eddington_tensor_dot_vector(SphP[i].ET[kf],vel_i,vdot_h); // calculate volume integral of scattering coefficient t_inv * (gas_vel . [e_rad*I + P_rad_tensor]), which gives an additional time-derivative term. this is the P term //
            for(k_dir=0;k_dir<3;k_dir++) {vdot_h[k_dir] = E_rad_forflux * (vel_i[k_dir] + vdot_h[k_dir]);} // and this is the eI term, multiply both by radiation energy to use in this step //
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR
            for(k_dir=0;k_dir<3;k_dir++) {DeltaFluxEff[k_dir] -= (P[i].Mass/rho) * (C_LIGHT_CODE_REDUCED*C_LIGHT_CODE_REDUCED/teqm_inv) * SphP[i].Gradients.Rad_E_gamma_ET[kf][k_dir];}
#else
            for(k_dir=0;k_dir<3;k_dir++) {DeltaFluxEff[k_dir] += (SphP[i].Dt_Rad_Flux[kf][k_dir]/teqm_inv);}
#endif
            for(k_dir=0;k_dir<3;k_dir++) {DeltaFluxEff[k_dir] += vdot_h[k_dir] * dt_entr;} // add the 'enthalpy advection' term here, vdot_h = Erad v.(e*I + P_rad)

            double tau=dt_entr*teqm_inv, f00=exp(-tau), f11=1.-f00;
            if(tau > 0 && isfinite(tau))
            {
                if(tau <= 0.04) {f11=tau-0.5*tau*tau+tau*tau*tau/6.; f00=1.-f11;} // some limits to prevent small/large number problems here
                if(!isfinite(f00) || !isfinite(f11)) {f00=1.; f11=0.;} // some limits to prevent small/large number problems here
                if(tau >= 20.) {f00=DMAX(0.,DMIN(1.e-11,f00)); f11=1.-f00;} // some limits to prevent small/large number problems here
                for(k_dir=0;k_dir<3;k_dir++)
                {
                    double flux_0; if(mode==0) {flux_0 = SphP[i].Rad_Flux[kf][k_dir];} else {flux_0 = SphP[i].Rad_Flux_Pred[kf][k_dir];}
                    flux_0 += vel_i[k_dir] * de_emission_minus_absorption; // add Lorentz term from net energy injected by absorption and re-emission (effectively, we operator-split this term and solve it -BEFORE- going to the next step)
                    double flux_f = flux_0 * f00 + DeltaFluxEff[k_dir] * f11; // exact solution for dE/dt = -E*abs + de , the 'reduction factor' appropriately suppresses the source term //
                    if(mode==0) {SphP[i].Rad_Flux[kf][k_dir] = flux_f;} else {SphP[i].Rad_Flux_Pred[kf][k_dir] = flux_f;}
                    f_mag += flux_f*flux_f; // magnitude of flux vector
                }
                if(f_mag > 0) // limit the flux according the physical (optically thin) maximum //
                {
                    f_mag=sqrt(f_mag); double fmag_max = 1.0 * C_LIGHT_CODE * ef; // maximum flux should be optically-thin limit: e_gamma/c: here allow some tolerance for numerical leapfrogging in timestepping. NOT the reduced RSOL here.
                    if(f_mag > fmag_max) {for(k_dir=0;k_dir<3;k_dir++) {if(mode==0) {SphP[i].Rad_Flux[kf][k_dir] *= fmag_max/f_mag;} else {SphP[i].Rad_Flux_Pred[kf][k_dir] *= fmag_max/f_mag;}}}
                }
            }
#endif
            total_erad_emission_minus_absorption += de_emission_minus_absorption; // add to cumulative sum for back-reaction to gas
#if defined(RT_EVOLVE_INTENSITIES)
            de_emission_minus_absorption_saved[k_tmp][k_angle] = de_emission_minus_absorption; // save this for use below
#endif
        } // clause for radiation angle [needed for evolving intensities]
    } // loop over frequencies
    
#if defined(RT_EVOLVE_INTENSITIES)
    for(kf=0;kf<N_RT_FREQ_BINS;kf++)
    {
        int k,k_om; double rho=SphP[i].Density*All.cf_a3inv, ceff=C_LIGHT_CODE_REDUCED, teq_inv=SphP[i].Rad_Kappa[kf]*rho*ceff, beta[3], f_a=rt_absorb_frac_albedo(i,kf), f_s=1.-f_a, b_dot_n[N_RT_INTENSITY_BINS]={0}, beta_2=0.;
        int n_iter = 1 + (int)(DMIN(DMAX(4. , dt_entr/teq_inv), 1000.)); // number of iterations to subcycle everything below //
        double dt=dt_entr/n_iter, tau=dt*teq_inv, i0[N_RT_INTENSITY_BINS]={0}, invfourpi=1./(4.*M_PI), J, b_dot_H, b2_dot_K; int i_iter;
        for(i_iter=0; i_iter<n_iter; i_iter++)
        {
            double egy_0=0,flux_0[3]={0},egy_f=0,flux_f[3]={0}; // compute total change over sub-cycle, to update gas properties
            // load all the gas and intensity properties we need [all can change on the subcycle so some re-computing here]
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if(mode==0) {i0[k_om] = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity[kf][k_om];} else {i0[k_om] = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity_Pred[kf][k_om];}}
            for(k=0;k<3;k++) {if(mode==0) {beta[k]=P[i].Vel[k]/(All.cf_atime*ceff);} else {beta[k]=SphP[i].VelPred[k]/(All.cf_atime*ceff);}} // need gas velocity at this time
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {b_dot_n[k_om]=0; for(k=0;k<3;k++) {b_dot_n[k_om]+=All.Rad_Intensity_Direction[k_om][k]*beta[k];}}
            beta_2=0; for(k=0;k<3;k++) {beta_2+=beta[k]*beta[k];}
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {egy_0+=i0[k_om]; for(k=0;k<3;k++) {flux_0[k]+=All.Rad_Intensity_Direction[k_om][k]*i0[k_om];}}
            J=0,b_dot_H=0,b2_dot_K=0; for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {J+=i0[k_om]*invfourpi; b_dot_H=b_dot_n[k_om]*i0[k_om]*invfourpi; b2_dot_K=b_dot_n[k_om]*b_dot_n[k_om]*i0[k_om]*invfourpi;}

            // isotropic terms that change total energy in bin (part of the 'work term' for the photon momentum)
            double work = ((f_s-f_a)*(beta_2*J + b2_dot_K) - 2.*f_s*b_dot_H) * tau; // will be shared isotropically.
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if((work>0) || (i0[k_om]<=0)) {i0[k_om]+=work;} else {i0[k_om]/=(1-work/i0[k_om]);}} // gaurantees linearized sum is still correct, and symmetric with positive changes, but can't get negative energies. shared isotropically.
            
            J=0; for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {J+=i0[k_om]*invfourpi;} // prepare to calculate isotropic scattering term
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {i0[k_om] = J + (i0[k_om]-J)*exp(-f_s*tau);} // isotropic scattering conserving total energy over step
            
            double fboost[N_RT_INTENSITY_BINS], fboost_avg=0, fboost_avg_p=0, fboost_avg_m=0; // calculate flux 'boost' terms
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {fboost[k_om] = 3.*b_dot_n[k_om] * ((de_emission_minus_absorption_saved[kf][k_om]/dt_entr) + (f_s*J + f_a*i0[k_om])); fboost_avg += fboost[k_om]/N_RT_INTENSITY_BINS;} // pre-calculate to get mean value, will divide out
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {work=(fboost[k_om]-fboost_avg)*tau; if((work>0) || (i0[k_om]<=0)) {fboost[k_om]=work; fboost_avg_p+=fboost[k_om];} else {fboost[k_om]=work/(1.-work/i0[k_om]); fboost_avg_m+=fboost[k_om];}} // zero total energy change at linear order ensured by subtracting out sum here; non-linearization ensures i0 cannot be negative, but does allow second-order dt work term to appear, that's ok for now
            if(fboost_avg_p>0 && fboost_avg_m<0) {double fc=-fboost_avg_m/fboost_avg_p; fboost_avg_m=(1.+fc)/(1.+fc*fc); fboost_avg_p=fc*fboost_avg_m;} else {fboost_avg_m=fboost_avg_p=0;} // // these re-weight to gaurantee the non-linear sum is identically zero while preserving positive-definite behavior
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if(fboost[k_om]>0) {i0[k_om]+=fboost_avg_p*fboost[k_om];} else {i0[k_om]+=fboost_avg_m*fboost[k_om];}} // alright done!
            
            // flux work term, allowed to both do work and be asymmetric so just need to ensure it retains positive-definite intensities
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {work=b_dot_n[k_om]*(f_a+f_s)*i0[k_om]; if((work>0) || (i0[k_om]<=0)) {i0[k_om]+=work;} else {i0[k_om]/=(1-work/i0[k_om]);}}

            // ok -now- calculate the net change in momentum and energy, for updating the gas quantities
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {egy_f+=i0[k_om]; for(k=0;k<3;k++) {flux_f[k]+=All.Rad_Intensity_Direction[k_om][k]*i0[k_om];}}
            double dv_gas[3]={0}, ke_gas_0=0, ke_gas_f=0, v0g=0, u0=0;
            for(k=0;k<3;k++) {dv_gas[k] = -(flux_f[k]-flux_0[k])/(ceff*P[i].Mass); v0g=ceff*beta[k]; ke_gas_0+=(v0g*v0g); ke_gas_f+=(v0g+dv_gas[k])*(v0g+dv_gas[k]);} // note everything is volume-integrated, accounted for above, and we defined flux for convience without the c, so just one power of c here.
            double d_ke_gas = 0.5*(ke_gas_f - ke_gas_0)*P[i].Mass, de_gas=-(egy_f-egy_0), de_gas_internal=(de_gas-d_ke_gas)/P[i].Mass;
            if(mode==0) {u0=SphP[i].InternalEnergy;} else {u0=SphP[i].InternalEnergyPred;} // for updating gas internal energy (work terms, after subtracting kinetic energy changes)
            if(de_gas_internal<=-0.9*u0) {de_gas_internal = DMIN(de_gas_internal/(1.-de_gas_internal/u0), -0.9*u0);} // just a catch to avoid negative energies (will break energy conservation if you are slamming into it, however!
            
            // assign everything back to the appropriate variables after update
            for(k=0;k<3;k++) {if(mode==0) {P[i].Vel[k] += dv_gas[k]*All.cf_atime;} else {SphP[i].VelPred[k] += dv_gas[k]*All.cf_atime;}} // update gas velocities (radiation pressure forces here)
            if(mode==0) {SphP[i].InternalEnergy += de_gas_internal;} else {SphP[i].InternalEnergyPred += de_gas_internal;} // update gas internal energy (work terms, after subtracting kinetic energy changes)
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if(mode==0) {SphP[i].Rad_Intensity[kf][k_om] = i0[k_om]/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[kf][k_om] = i0[k_om]/RT_INTENSITY_BINS_DOMEGA;}} // update intensities (all of the above)
            SphP[i].Rad_E_gamma[kf]=egy_f; // set this every time this subroutine is called, so it is accessible everywhere else //
        } // loop over iterations
    } // loop over frequencies
#else
    double mom_fac = 1. - total_erad_emission_minus_absorption / (P[i].Mass * C_LIGHT_CODE_REDUCED*C_LIGHT_CODE_REDUCED); // back-reaction on gas from emission
    {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {if(mode==0) {P[i].Vel[k_dir] *= mom_fac;} else {SphP[i].VelPred[k_dir] *= mom_fac;}}}
#endif

    if(mode > 0) {rt_eddington_update_calculation(i);} /* update the eddington tensor (if we calculate it) as well */
#endif
}



#endif





#ifdef RADTRANSFER
/***********************************************************************************************************/
/* this function initializes some of the variables we need */
/***********************************************************************************************************/
void rt_set_simple_inits(int RestartFlag)
{
    if(RestartFlag==1) return;
    
    int i; for(i = 0; i < NumPart; i++)
    {
        if(P[i].Type == 0)
        {
            int k;
#ifdef RT_INFRARED
            SphP[i].Dust_Temperature = All.InitGasTemp; //get_min_allowed_dustIRrad_temperature(); // in K, floor = CMB temperature or 10K
            SphP[i].Radiation_Temperature = All.InitGasTemp; //SphP[i].Dust_Temperature;
            SphP[i].Dt_Rad_E_gamma_T_weighted_IR = 0;
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            for(k=0;k<3;k++) {SphP[i].Rad_Accel[k]=0;}
#endif
#ifdef RT_CHEM_PHOTOION
            SphP[i].HII = MIN_REAL_NUMBER;
            SphP[i].HI = 1.0 - SphP[i].HII;
            SphP[i].Ne = SphP[i].HII;
#ifdef RT_CHEM_PHOTOION_HE
            double fac = (1-HYDROGEN_MASSFRAC) / (4.0 * HYDROGEN_MASSFRAC);
            SphP[i].HeIII = MIN_REAL_NUMBER * fac;
            SphP[i].HeII = MIN_REAL_NUMBER * fac;
            SphP[i].HeI = (1.0 - SphP[i].HeII - SphP[i].HeIII) * fac;
            SphP[i].Ne += SphP[i].HeII + 2.0 * SphP[i].HeIII;
#endif
#endif
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if(RestartFlag==0) {SphP[i].Rad_E_gamma[k] = MIN_REAL_NUMBER;}
                SphP[i].ET[k][0]=SphP[i].ET[k][1]=SphP[i].ET[k][2]=1./3.; SphP[i].ET[k][3]=SphP[i].ET[k][4]=SphP[i].ET[k][5]=0;
                SphP[i].Rad_Je[k] = 0; SphP[i].Rad_Kappa[k] = rt_kappa(i,k);
#ifdef RT_FLUXLIMITER
                SphP[i].Rad_Flux_Limiter[k] = 1;
#endif
#ifdef RT_INFRARED
                if(k==RT_FREQ_BIN_INFRARED) {SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED] = 5.67e-5 * 4 / (C_LIGHT * RT_SPEEDOFLIGHT_REDUCTION) * pow(All.InitGasTemp,4.) / UNIT_PRESSURE_IN_CGS * P[i].Mass / (SphP[i].Density*All.cf_a3inv);}
#endif
#ifdef RT_EVOLVE_ENERGY
                SphP[i].Rad_E_gamma_Pred[k] = SphP[i].Rad_E_gamma[k]; SphP[i].Dt_Rad_E_gamma[k] = 0;
#endif
#ifdef RT_EVOLVE_FLUX
                int k_dir; for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Rad_Flux_Pred[k][k_dir] = SphP[i].Rad_Flux[k][k_dir] = SphP[i].Dt_Rad_Flux[k][k_dir] = 0;}
#endif
#ifdef RT_EVOLVE_INTENSITIES
                int k_dir; for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++) {SphP[i].Rad_Intensity_Pred[k][k_dir] = SphP[i].Rad_Intensity[k][k_dir] = MIN_REAL_NUMBER; SphP[i].Dt_Rad_Intensity[k][k_dir] = 0;}
#endif
                
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                double q_a=0.75*GRAIN_RDI_TESTPROBLEM_Q_AT_GRAIN_MAX/All.Grain_Size_Max, e0=All.Vertical_Grain_Accel/q_a, kappa0=All.Dust_to_Gas_Mass_Ratio*q_a;
                e0 *= (P[i].Mass/SphP[i].Density) * exp(-kappa0*(1.-exp(-P[i].Pos[2]))); // attenuate according to equilibrium expectation, if we're using single-scattering radiation pressure [otherwise comment this line out] //
                SphP[i].Rad_E_gamma_Pred[k]=SphP[i].Rad_E_gamma[k]=e0;
#if defined(RT_EVOLVE_FLUX)
                SphP[i].Rad_Flux_Pred[k][2]=SphP[i].Rad_Flux[k][2] = e0*C_LIGHT_CODE_REDUCED;
                SphP[i].Rad_Flux[k][0]=SphP[i].Rad_Flux[k][1]=SphP[i].Rad_Flux_Pred[k][0]=SphP[i].Rad_Flux_Pred[k][1]=0;
#endif
#endif
                
            }
        }
    }
}
#endif



#if defined(RT_EVOLVE_INTENSITIES)
/***********************************************************************************************************/
/* routine to initialize the distribution of ray angles along which the intensities are explicitly evolved.
    this follows Bruhls et al. 1999 [Appendix B]. the rays uniformly sample the unit sphere, are distributed
    isotropically, and satisfy the numerical quadtratures through second order
    J=SUM[dOmega*I]/4PI=I, F=SUM[dOmega*I*Omega_hat]=0, K=SUM[dOmega*I*Omega_hat.x.Omega_hat]/4PI=(1/3)*Identity*J
    to machine accuracy for an isotropic (I=constant along all discrete directions) radiation field.
    The user specifies the number of independent polar angles Np=RT_LOCALRAYGRID (any integer). The total number of rays
    per octant of the unit sphere is then Np*(Np+1)/2, and the total number of rays
    altogether is 4*Np*(Np+1). The exact numerical quadrature assumption is used throughout.
 */
/***********************************************************************************************************/
void rt_init_intensity_directions(void)
{
    int n_polar = RT_LOCALRAYGRID;
    if(n_polar < 1) {printf("Number of rays is invalid (<1). Terminating.\n"); endrun(5346343);}

    double mu[n_polar]; int i,j,k,l,n=0,n_oct=n_polar*(n_polar+1)/2;
    double Rad_Intensity_Direction_tmp[n_oct][3];
    for(j=0;j<n_polar;j++) {mu[j] = sqrt( (j + 1./6.) / (n_polar - 1./2.) );}
    
    for(i=0;i<n_polar;i++)
    {
        for(j=0;j<n_polar-i;j++)
        {
            k=n_polar-1-i-j;
            Rad_Intensity_Direction_tmp[n][0]=mu[i]; Rad_Intensity_Direction_tmp[n][1]=mu[j]; Rad_Intensity_Direction_tmp[n][2]=mu[k];
            n++;
        }
    }
    n=0;
    for(i=0;i<2;i++)
    {
        double sign_x = 1 - 2*i;
        for(j=0;j<2;j++)
        {
            double sign_y = 1 - 2*j;
            for(k=0;k<2;k++)
            {
                double sign_z = 1 - 2*k;
                for(l=0;l<n_oct;l++)
                {
                    All.Rad_Intensity_Direction[n][0] = Rad_Intensity_Direction_tmp[l][0] * sign_x;
                    All.Rad_Intensity_Direction[n][1] = Rad_Intensity_Direction_tmp[l][1] * sign_y;
                    All.Rad_Intensity_Direction[n][2] = Rad_Intensity_Direction_tmp[l][2] * sign_z;
                    n++;
                }
            }
        }
    }
}
#endif


/***********************************************************************************************************/
/* optional routine to distribute cooling or other sources from gas: currently empty */
/***********************************************************************************************************/
void rt_get_lum_gas(int target, double *je)
{
#ifdef RT_FREEFREE
    int k = RT_FREQ_BIN_FREEFREE;
    double t_eff = 0.59 * (GAMMA(target)-1.) * U_TO_TEMP_UNITS * SphP[target].InternalEnergyPred; // we're assuming fully-ionized gas with a simple equation-of-state here, nothing fancy, to get the temperature //
    je[k] = rt_absorb_frac_albedo(target,k) * rt_kappa(target,k) * P[target].Mass * ((4. * 5.67e-5) * t_eff*t_eff*t_eff*t_eff) / UNIT_FLUX_IN_CGS; // blackbody emissivity (Kirchoff's law): account for albedo [absorption opacity], and units //
#endif
}



/***********************************************************************************************************/
/* below returns (1-exp(-x))/x , which is needed for averages through slabs and over time for radiative quantities. here for convenience */
/***********************************************************************************************************/
double slab_averaging_function(double x)
{
    /* this fitting function is accurate to ~0.1% at all x, and extrapolates to the correct limits without the pathological behaviors of (1-exp(-x))/x for very small/large x */
    return (1.00000000000 + x*(0.21772719088733913 + x*(0.047076512011644776 + x*0.005068307557496351))) /
           (1.00000000000 + x*(0.71772719088733920 + x*(0.239273440788647680 + x*(0.046750496137263675 + x*0.005068307557496351))));
}



#ifdef RT_INFRARED

/* return the minimum user-specified dust temperature. note there is nothing physical about this, just a convenience function since we enforce a minimum -gas- temperature */
double get_min_allowed_dustIRrad_temperature(void)
{
#if defined(GALSF)
    return DMAX(All.MinGasTemp, 2.73/All.cf_atime);
#endif
    return MIN_REAL_NUMBER;
}

/* return LambdaDust, the dust heating/cooling rate (>0 is heating, <0 is cooling) */
double get_rt_ir_lambdadust_effective(double T, double rho, double *nH0_guess, double *ne_guess, int target)
{
#ifdef COOLING
    double Tdust_0 = SphP[target].Dust_Temperature; // dust temperature estimate from previous loop over radiation operators
    double LambdaDust_initial_guess = 1.116e-32 * (Tdust_0-T) * sqrt(T)*(1.-0.8*exp(-75./T)) * (P[target].Metallicity[0]/All.SolarAbundances[0]); // guess value based on the -current- values of T, Tdust //
        
    double egy_therm = SphP[target].InternalEnergyPred*P[target].Mass; // true internal energy (before this cooling loop)
    double egy_rad = SphP[target].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED]; // true radiation field energy (before this cooling loop)
    double egy_tot = egy_rad + egy_therm; // true total energy [in code units]
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS;    // effective hydrogen number dens in cgs units (just for normalization convention)
    double volume = (P[target].Mass / (SphP[target].Density*All.cf_a3inv)); // particle volume in code units
    double ratefact = (nHcgs*nHcgs) * volume / (UNIT_PRESSURE_IN_CGS /UNIT_TIME_IN_CGS); // conversion b/t Lambda and du used by code
    double Erad_to_T4_fac = RT_SPEEDOFLIGHT_REDUCTION * 1.32e14 * UNIT_PRESSURE_IN_CGS / volume; // conversion from absolute rad energy to T^4 units, used multiple places below, coefficient = cL_reduced/(4*sigma_B)
    double Teff = Get_Gas_Mean_Molecular_Weight_mu(T, rho, nH0_guess, ne_guess, 0, target) * (GAMMA(target)-1.) * U_TO_TEMP_UNITS * (egy_tot / P[target].Mass); // convert from internal energy to temperature units for factor below

    double xf, a = Teff*Teff*Teff*Teff / (Erad_to_T4_fac*egy_tot); // dimensionless factors needed to solve for the equilibrium Tdust-Tgas relation
    if(a<0.2138) {xf=(1+19*a+132*a*a+418*a*a*a+580*a*a*a*a+243*a*a*a*a*a)/(1+20*a+148*a*a+508*a*a*a+796*a*a*a*a+432*a*a*a*a*a);} // eqm solution (power series approx)
     else {double a0=pow(a,0.25); xf=(-704-1045*a0+128*a0*a0*a0*(39+32*a0*(4+7*a0+64*a0*a0*a0*(-1+8*a0*(-1+4*a0)))))/(8388608.*a*a*a0*a0);} // eqm solution (power series approx)

    double L0_abs = fabs(LambdaDust_initial_guess); // absolute value of the initially-computed guess for the cooling/heating rate
    double Edot0 = L0_abs * ratefact; // now this is an absolute Edot in code units
    double efinal_minus_einitial = egy_tot*xf - egy_therm; // change in energy if we went all the way to equilibrium
    double t_cooling_eff = fabs(efinal_minus_einitial) / Edot0; // effective cooling time at the initially-estimated rate here
    double sign_term=1.; if(efinal_minus_einitial < 0.) {sign_term=-1.;} // sign of the cooling/heating (to keep for below)
    double dt = (P[target].TimeBin ? (((integertime) 1) << P[target].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; // timestep being taken [code units]
    double tau = dt/t_cooling_eff, xfac=(1.-exp(-tau))/tau; if(tau<0.05) {xfac=1.-0.5*tau+tau*tau/6.;} else {if(tau>20.) {xfac=1./tau;}} // correct rate to asymptote to equilibrium
    double lambda_eff = sign_term * L0_abs * xfac; // final effective cooling/heating rate

    SphP[target].Dust_Temperature = DMAX(pow(Erad_to_T4_fac*DMAX( 0., egy_rad - lambda_eff*ratefact*dt ), 0.25), get_min_allowed_dustIRrad_temperature()); // update dust temperature guess //
    return lambda_eff;
#endif
    return 0;
}

#endif




#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE)
#ifdef CHIMES_STELLAR_FLUXES
/* The following routines are fitting functions that are used to
 * obtain the luminosities in the 6-13.6 eV energy band (i.e. G0)
 * and the >13.6 eV band (i.e. H-ionising), which will be used
 * by CHIMES. These functions were fit to Starburst99 models
 * that used the Geneva 2012/13 tracks with v=0.4 rotation
 * and Z=0.014 metallicity. */

double chimes_G0_luminosity(double stellar_age, double stellar_mass)
{
  // stellar_age in Myr.
  // stellar_mass (current, not initial) in Msol.
  // return value in Habing units * cm^2.
  double zeta = 6.5006802e29;
  if (stellar_age < 4.07)
    return stellar_mass * exp(89.67 + (0.172 * pow(stellar_age, 0.916)));
  else
    return stellar_mass * zeta * pow(1773082.52 / stellar_age, 1.667) * pow(1.0 + pow(stellar_age / 1773082.52, 28.164), 1.64824);
}

double chimes_ion_luminosity(double stellar_age, double stellar_mass)
{
  // stellar_age in Myr.
  // stellar_mass (current, not initial) in Msol.
  // return value in s^-1.
  double zeta = 3.2758118e21;
  if (stellar_age < 3.71)
    return stellar_mass * exp(107.21 + (0.111 * pow(stellar_age, 0.974)));
  else
    return stellar_mass * zeta * pow(688952.27 / stellar_age, 4.788) * pow(1.0 + pow(stellar_age / 688952.27, 1.124), -17017.50356);
}
#endif
#endif
