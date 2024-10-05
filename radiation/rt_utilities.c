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

#define SET_ACTIVE_RT_CHECK() if(mode<0) {return 1;} else {active_check=1;}


/***********************************************************************************************************/
/* routine which returns the luminosity [total volume/mass integrated] for the desired source particles in physical code units (energy/time),
    as a function of whatever the user desires, in the relevant bands. inpute here:
    'i' = index of target particle/cell for which the luminosity should be computed
    'mode' = flag for special behaviors. if <0 (e.g. -1), just returns whether or not a particle is 'active' (eligible as an RT source). if =0, normal behavior. if =1, then some bands have special behavior, for example self-shielding estimated -at the source-
    'lum' = pointer to vector of length N_RT_FREQ_BINS to hold luminosities for all bands.
    Note that for a number of the bands below where 'sources' are stars or star-like objects,
        there are two default 'versions' implemented: one assuming the sources are -individual- stars/protostars/compact objects, the other
        assuming the sources represent stellar -populations-. Make sure you implement the correct assumptions with appropriate flags for the simulations you wish to run.
 */
/***********************************************************************************************************/
int rt_get_source_luminosity(int i, int mode, double *lum)
{
    if(!((1 << P[i].Type) & (RT_SOURCES))) {return 0;}; // boolean test of whether i is a source or not - end if not a valid source particle
    if(P[i].Mass <= 0) {return 0;} // reject invalid particles scheduled for deletion
    int active_check = 0; // default to inactive //
    
#if defined(GALSF)
#if defined(SINGLE_STAR_SINK_DYNAMICS)
    active_check += rt_get_lum_band_singlestar(i,mode,lum); // get luminosities for individual star/sink particles assuming they are protostars or stars
#else
    active_check += rt_get_lum_band_stellarpopulation(i,mode,lum); // get luminosities for star particles assuming they represent IMF-averaged populations
#if defined(BLACK_HOLES)
    active_check += rt_get_lum_band_agn(i,mode,lum); // get luminosities for BH/sink particles assuming they represent AGN
#endif
#endif
#endif
    if(mode < 0 && active_check) {return 1;} // if got a positive answer already, that's all we are checking here, we are done

    
#if defined(RT_CHEM_PHOTOION) && !defined(GALSF) /* Hydrogen and Helium ionizing bands; this is an idealized test-problem version implementation */
    if(P[i].Type==4)
    {
        SET_ACTIVE_RT_CHECK(); double l_ion=All.IonizingLuminosityPerSolarMass_cgs * (P[i].Mass * UNIT_MASS_IN_SOLAR) / UNIT_LUM_IN_CGS; // flux from star particles according to mass
#ifdef RT_ILIEV_TEST1
        l_ion = 5.0e48 * (13.6*ELECTRONVOLT_IN_ERGS) / UNIT_LUM_IN_CGS; // 5e48 in ionizing photons per second -- constant for idealized test problem //
#endif
        lum[RT_FREQ_BIN_H0] = l_ion; // default to all flux into single-band
#if defined(RT_PHOTOION_MULTIFREQUENCY)
        int i_vec[4] = {RT_FREQ_BIN_H0, RT_FREQ_BIN_He0, RT_FREQ_BIN_He1, RT_FREQ_BIN_He2}; // these will all be the same if not using multi-frequency module //
        int k; for(k=0;k<4;k++) {lum[i_vec[k]] = l_ion * rt_ion_precalc_stellar_luminosity_fraction[i_vec[k]];} // assign flux appropriately according to pre-tabulated result //
#endif
    }
#endif

    
#if defined(RT_GENERIC_USER_FREQ)   /* example code to be modified as-needed for custom RT problems */
    if(P[i].Type == 4) // set this to whichever type you want to use for the specific sources in this band
    {
        SET_ACTIVE_RT_CHECK(); // flag that tells the code that indeed this particle should be active!
        lum[RT_FREQ_BIN_GENERIC_USER_FREQ] = 0; // set the actual luminosity here for your test problem!
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION /* assume special units for this problem, and that total mass of 'sources' is 1 */
        double m_total_expected = 1; // assume total mass of sources is 1, and we want to weight such that fractional emission per source is equal to their mass fraction
        double A_base = boxSize_X * boxSize_X; // area of the base of the box used for scaling to get the desired flux
#if (NUMDIMS == 3)
        A_base = boxSize_X * boxSize_Y;
#endif
        lum[RT_FREQ_BIN_GENERIC_USER_FREQ] = (P[i].Mass/1.) * All.Vertical_Grain_Accel * C_LIGHT_CODE * ((All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS)*(All.Grain_Size_Max/UNIT_LENGTH_IN_CGS)) * A_base / (0.75*All.Grain_Q_at_MaxGrainSize); // special behavior for particular test of stratified boxes compared to explicit dust opacities
#endif
    }
#endif
    
    
#ifdef RADTRANSFER
    if(P[i].Type == 0) /* generic sub routines for gas as a source term, should go at the very end of this routine */
    {
        SET_ACTIVE_RT_CHECK(); rt_get_lum_gas(i,lum); // optionally re-distributes cooling flux as a blackbody; but also where bands like free-free reside //
        int k; for(k=0;k<N_RT_FREQ_BINS;k++) {lum[k] += SphP[i].Rad_Je[k];}
    }
#endif
    
    /* need to renormalize ALL sources for reduced speed of light */
    {int k; for(k=0;k<N_RT_FREQ_BINS;k++) {lum[k] *= (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE);}}
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
    return SphP[i].Interpolated_Opacity[k_freq] + 1.e-3 * All.Dust_to_Gas_Mass_Ratio*0.75*All.Grain_Q_at_MaxGrainSize/((All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS)*(All.Grain_Size_Max/UNIT_LENGTH_IN_CGS)); /* enforce minimum; note kappa in code units here so need to convert appropriately */
#endif
    return MIN_REAL_NUMBER + SphP[i].Interpolated_Opacity[k_freq]; /* this is calculated in a different routine, just return it now */
#endif

#ifdef RT_CHEM_PHOTOION
    /* opacity to ionizing radiation for Petkova & Springel bands. note cooling.c or rt_update_chemistry is where ionization is actually calculated */
    double nH_over_Density = HYDROGEN_MASSFRAC / PROTONMASS_CGS * UNIT_MASS_IN_CGS;
    double kappa = nH_over_Density * (SphP[i].HI + MIN_REAL_NUMBER) * rt_ion_sigma_HI[k_freq]; // note this is designed for specific applications: does not include dust, or free-free, or free-electron scattering contributions here, all of which can be important.
#if defined(RT_CHEM_PHOTOION_HE) && defined(RT_PHOTOION_MULTIFREQUENCY)
    kappa += nH_over_Density * ((SphP[i].HeI + MIN_REAL_NUMBER) * rt_ion_sigma_HeI[k_freq] + (SphP[i].HeII + MIN_REAL_NUMBER) * rt_ion_sigma_HeII[k_freq]);
    if(k_freq==RT_FREQ_BIN_He0)  {return kappa;}
    if(k_freq==RT_FREQ_BIN_He1)  {return kappa;}
    if(k_freq==RT_FREQ_BIN_He2)  {return kappa;}
#endif
    if(k_freq==RT_FREQ_BIN_H0)  {return kappa;}
#endif

#if defined(RT_HARD_XRAY) || defined(RT_SOFT_XRAY) || defined(RT_PHOTOELECTRIC) || defined (GALSF_FB_FIRE_RT_LONGRANGE) || defined(RT_NUV) || defined(RT_OPTICAL_NIR) || defined(RT_LYMAN_WERNER) || defined(RT_INFRARED) || defined(RT_FREEFREE)
    double fac = UNIT_SURFDEN_IN_CGS, Zfac, dust_to_metals_vs_standard, kappa_HHe; /* units */
    Zfac = 1.0; kappa_HHe=0.35; // assume solar metallicity, simple Thompson cross-section limit for various processes below
    dust_to_metals_vs_standard = return_dust_to_metals_ratio_vs_solar(i,0); // many of the dust opacities below will need this as the dimensionless dust-to-metals ratio normalized to the canonical Solar value of ~1/2
#ifdef METALS
    if(i>=0) {Zfac = P[i].Metallicity[0]/All.SolarAbundances[0];}
#endif
#ifdef COOLING
    if(i>=0) {kappa_HHe=0.02 + 0.35*SphP[i].Ne;}
#endif

    
#ifdef RT_FREEFREE /* pure (grey, non-relativistic) Thompson scattering opacity + free-free absorption opacity. standard expressions here from Rybicki & Lightman. */
    if(k_freq==RT_FREQ_BIN_FREEFREE)
    {
        double T_eff=0.59*(GAMMA(i)-1.)*U_TO_TEMP_UNITS*SphP[i].InternalEnergyPred, rho=SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS; // we're assuming fully-ionized gas with a simple equation-of-state here, nothing fancy, to get the temperature //
        double kappa_abs = 1.e30*rho*pow(T_eff,-3.5);
        return (0.35 + kappa_abs) * fac;
    }
#endif
#ifdef RT_HARD_XRAY
    /* opacity comes from H+He (Thompson) + metal ions. expressions here for metal ions come from integrating over the standard Morrison & McCammmon 1983 metal cross-sections for a standard slope gamma in the band */
    if(k_freq==RT_FREQ_BIN_HARD_XRAY) {return (0.53 + 0.27*Zfac) * fac;}
#endif
#ifdef RT_SOFT_XRAY
    /* opacity comes from H+He (Thompson) + metal ions. expressions here for metal ions come from integrating over the standard Morrison & McCammmon 1983 metal cross-sections for a standard slope gamma in the band */
    if(k_freq==RT_FREQ_BIN_SOFT_XRAY) {return (127. + 50.0*Zfac) * fac;}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    /* three-band (UV, OPTICAL, IR) approximate spectra for stars as used in the FIRE (Hopkins et al.) models. mean opacities here come from integrating over the Hopkins 2004 (Pei 1992) opacities versus wavelength for the large bands here, using a luminosity-weighted mean stellar spectrum from the same starburst99 models used to compute the stellar feedback */
    if(k_freq==RT_FREQ_BIN_FIRE_UV)  {return (1800.) * fac;}
    if(k_freq==RT_FREQ_BIN_FIRE_OPT) {return (180.)  * fac;} /* note this is roughly equivalent to the specific extinction at R-band [bit higher in UBV, lower in IJHK] */
    if(k_freq==RT_FREQ_BIN_FIRE_IR)  {return (10.) * fac * (0.1 + Zfac);}
#endif
#ifdef RT_PHOTOELECTRIC
    /* opacity comes primarily from dust (ignoring H2 molecular opacities here). band is 8-13.6 eV [0.091-0.155 micron] (note overlap with LW band) */
    if(k_freq==RT_FREQ_BIN_PHOTOELECTRIC) {return DMAX(kappa_HHe, 720.*DMAX(1.e-4,Zfac*dust_to_metals_vs_standard)) * fac;} // depends rather sensitively on assumed input spectrum+dust composition (e.g. MW vs SMC-like), using 2007+ Draine+Li MW-like dust results here, weighted over a flat spectrum with 1/2 of the weight for the portion of the band which overlaps with LW.
#endif
#ifdef RT_LYMAN_WERNER
    /* opacity from molecular H2 and dust (dominant at higher-metallicity) should be included. band is 11.2-13.6 eV [0.111-0.155 micron] */
    if(k_freq==RT_FREQ_BIN_LYMAN_WERNER) {return DMAX(kappa_HHe, 900.*Zfac*dust_to_metals_vs_standard) * fac;} // just dust term for now, depends rather sensitively on assumed input spectrum+dust composition (e.g. MW vs SMC-like). using 2007+ Draine+Li MW-like dust results here, weighted by the LW cross-section and a flat spectrum in the range.
#endif
#ifdef RT_NUV
    /* opacity comes primarily from dust. effective waveband is 0.16-0.36 micron [~3.444-8 eV] */
    if(k_freq==RT_FREQ_BIN_NUV) {return DMAX(kappa_HHe, 480.*Zfac*dust_to_metals_vs_standard) * fac;} // depends rather sensitively on assumed input spectrum+dust composition (e.g. MW vs SMC-like), using 2007+ Draine+Li MW-like dust results here
#endif
#ifdef RT_OPTICAL_NIR
    /* opacity comes primarily from dust. effective waveband is 0.36-3 microns [~0.4133-3.444 eV], e.g. between U-K+ band  */
    if(k_freq==RT_FREQ_BIN_OPTICAL_NIR) {return DMAX(kappa_HHe, 180.*Zfac*dust_to_metals_vs_standard) * fac;} // this is close to the specific opacity at R-band, which you can treat very crudely as a sort of 'effective wavelength' for this
#endif
#ifdef RT_INFRARED
    /* IR with dust opacity */
    double T_min = get_min_allowed_dustIRrad_temperature();
    if(k_freq==RT_FREQ_BIN_INFRARED)
    {
        if(isnan(SphP[i].Dust_Temperature)) {PRINT_WARNING("\n NaN dust temperature for cell-ID=%llu  \n", (unsigned long long) P[i].ID); SphP[i].Dust_Temperature = 1.e4;}
        if(isnan(SphP[i].Radiation_Temperature)) {PRINT_WARNING("\n NaN gas temperature for cell-ID=%llu  \n", (unsigned long long) P[i].ID);}
        if(SphP[i].Dust_Temperature<=T_min) {SphP[i].Dust_Temperature=T_min;} // reset baseline
        if(SphP[i].Radiation_Temperature<=T_min) {SphP[i].Radiation_Temperature=T_min;} // reset baseline
        double T_dust_em = SphP[i].Dust_Temperature; // dust temperature in K //
        double Trad = SphP[i].Radiation_Temperature; // radiation temperature in K //
        if((Trad <= 0) || (T_dust_em<=0)) {PRINT_WARNING("\n Cell-ID=%llu  has T_rad=%g and T_dust=%g\n", (unsigned long long) P[i].ID, Trad, T_dust_em);}
        return rt_kappa_adaptive_IR_band(i,T_dust_em,Trad,0,0); // < 1500 K, dust is present; here first flag 0 uses the radiation temperature because we want to know the Planck-mean *absorption* opacity. Second flag 0 says to include both dust and gas opacities. In the subroutine, divide by fac because the function outputs in code units but we're working in CGS here
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
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS) && defined(RT_GENERIC_USER_FREQ)
    if(k_freq==RT_FREQ_BIN_GENERIC_USER_FREQ) {return DMAX(1.e-6, DMIN(1.0 - 1.e-6, All.Grain_Absorbed_Fraction_vs_Total_Extinction));}
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
        double fA_tmp = (1.-0.5/(1.+((725.*725.)/(1.+SphP[i].Radiation_Temperature*SphP[i].Radiation_Temperature)))); // rough interpolation depending on the radiation temperature: high Trad, this is 1/2, low Trad, gets closer to unity; need to revise for sublimated cases here ??? */
#ifdef COOLING
		if(rt_kappa(i,k_freq)>0) {fA_tmp *= (1.-DMIN(1.,0.35*SphP[i].Ne*fac/rt_kappa(i,k_freq)));} else {return 1.0;} // the value should not matter if rt_kappa=0 // ??? correct this for Klein-Nishina as well?
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




/* subroutine for 'rt_get_source_luminosity', with identical variables, for cases where the radiation
    represents IMF-averaged stellar populations, i.e. the sort of thing which would be used in galaxy simulations.
 */
int rt_get_lum_band_stellarpopulation(int i, int mode, double *lum)
{
    if(!((P[i].Type == 4) || ((All.ComovingIntegrationOn==0)&&((P[i].Type==2)||(P[i].Type==3))))) {return 0;} // only star-type particles act in this subroutine //
    int active_check = 0; // default to inactive //
#if defined(GALSF) /* basically none of these modules make sense without the GALSF module active */
    double star_age = evaluate_stellar_age_Gyr(i), m_sol = P[i].Mass * UNIT_MASS_IN_SOLAR;
    if((star_age<=0) || isnan(star_age)) {return 0;} // calculate stellar age, will be used below, and catch for bad values

    
#if defined(GALSF_FB_FIRE_RT_LONGRANGE) /* three-band (UV, OPTICAL, IR) approximate spectra for stars as used in the FIRE (Hopkins et al.) models */
    SET_ACTIVE_RT_CHECK();
    double f_uv=All.PhotonMomentum_fUV, f_op=All.PhotonMomentum_fOPT;
    double L = evaluate_light_to_mass_ratio(star_age, i) * m_sol / UNIT_LUM_IN_SOLAR; if(L<=0 || isnan(L)) {L=0;}
#ifndef RT_FIRE_FIX_SPECTRAL_SHAPE
    double sigma_eff = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,P[i].DensAroundStar,PPP[i].NumNgb,0,i); if((sigma_eff <= 0)||(isnan(sigma_eff))) {sigma_eff=0;} // sigma here is in code units
    if(star_age <= 0.0025) {f_op=0.09;} else {if(star_age <= 0.006) {f_op=0.09*(1+((star_age-0.0025)/0.004)*((star_age-0.0025)/0.004));} else {f_op=1-0.8410937/(1+sqrt((star_age-0.006)/0.3));}}
    /* note that the metallicity doing attenuation is the -gas- opacity around the star, while here we only know the stellar metallicity,
        so we use this as a guess, but this could substantially under-estimate opacities for old stars in MW-like galaxies. But for young stars (which dominate) this is generally ok. */
    double tau_uv = sigma_eff*rt_kappa(i,RT_FREQ_BIN_FIRE_UV); double tau_op = sigma_eff*rt_kappa(i,RT_FREQ_BIN_FIRE_OPT); // kappa returned in code units
    f_uv = (1-f_op)*(All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV)/(1+0.8*tau_uv+0.85*tau_uv*tau_uv));
    f_op *= All.PhotonMomentum_fOPT + (1-All.PhotonMomentum_fOPT)/(1+0.8*tau_op+0.85*tau_op*tau_op); /* this is a fitting function for tau_disp~0.22 'tail' w. exp(-tau) 'core', removes expensive functions [f_uv = (1-f_op)*(All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV)*exp(-tau_uv)); f_op *= All.PhotonMomentum_fOPT + (1-All.PhotonMomentum_fOPT)*exp(-tau_op);]
     :: accounting for leakage for P(tau) ~ exp(-|logtau/tau0|/sig), following Hopkins et al. 2011, we have: [f_uv = (1-f_op)*(All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV)/ (1 + pow(tau_uv,1./(4.*tau_disp))/(3.*tau_disp) + pow(2.*tau_disp*tau_uv,1./tau_disp))); f_op *= All.PhotonMomentum_fOPT + (1-All.PhotonMomentum_fOPT)/ (1 + pow(tau_op,1./(4.*tau_disp))/(3.*tau_disp) + pow(2.*tau_disp*tau_op,1./tau_disp));] */
#endif
    lum[RT_FREQ_BIN_FIRE_UV]  = L * f_uv;
    lum[RT_FREQ_BIN_FIRE_OPT] = L * f_op;
    lum[RT_FREQ_BIN_FIRE_IR]  = L * (1-f_uv-f_op);
#endif
    
#if defined(RT_INFRARED) /* can add direct infrared sources, but default to no direct IR (just re-emitted light) */
    lum[RT_FREQ_BIN_INFRARED] = 0; //default to no direct IR (just re-emitted light)
#endif

#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models */
    SET_ACTIVE_RT_CHECK();
    double f_op=0; if(star_age <= 0.0025) {f_op=0.09;} else {
        if(star_age <= 0.006) {f_op=0.09*(1+((star_age-0.0025)/0.004)*((star_age-0.0025)/0.004));} else {f_op=1-0.8410937/(1+sqrt((star_age-0.006)/0.3));}}
    lum[RT_FREQ_BIN_OPTICAL_NIR] = f_op * evaluate_light_to_mass_ratio(star_age, i) * m_sol / UNIT_LUM_IN_SOLAR;
#endif

#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models */
    SET_ACTIVE_RT_CHECK();
#if !defined(RT_OPTICAL_NIR)
    double f_op=0; if(star_age <= 0.0025) {f_op=0.09;} else {
        if(star_age <= 0.006) {f_op=0.09*(1+((star_age-0.0025)/0.004)*((star_age-0.0025)/0.004));} else {f_op=1-0.8410937/(1+sqrt((star_age-0.006)/0.3));}}
#endif
    lum[RT_FREQ_BIN_NUV] = (1-f_op) * evaluate_light_to_mass_ratio(star_age, i) * m_sol / UNIT_LUM_IN_SOLAR;
#endif

#if defined(RT_PHOTOELECTRIC) /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK();
    double l_band_pe, x_age_pe = star_age / 3.4e-3; // converts to code units, and defines age relative to convenient break time
    if(x_age_pe <= 1) {l_band_pe = 1.07e36 * (1.+x_age_pe*x_age_pe) * m_sol / UNIT_LUM_IN_CGS;}
        else {l_band_pe = 2.14e36 / (x_age_pe * sqrt(x_age_pe)) * m_sol / UNIT_LUM_IN_CGS;} // 0.1 solar, with nebular. very weak metallicity dependence, with slightly slower decay in time for lower-metallicity pops; effect smaller than binaries
    lum[RT_FREQ_BIN_PHOTOELECTRIC] = l_band_pe; // band luminosity //
#endif
    
#if defined(RT_LYMAN_WERNER)  /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK();
    double l_band_lw, x_age_lw = star_age / 3.4e-3; // converts to code units, and defines age relative to convenient break time
    if(x_age_lw <= 1) {l_band_lw = 0.429e36 * (1.+x_age_lw*x_age_lw) * m_sol / UNIT_LUM_IN_CGS;}
        else {l_band_lw = 0.962e36 * pow(x_age_lw,-1.6) * exp(-x_age_lw/117.6) * m_sol / UNIT_LUM_IN_CGS;} // 0.1 solar, with nebular. very weak metallicity dependence, with slightly slower decay in time for lower-metallicity pops; effect smaller than binaries
    lum[RT_FREQ_BIN_LYMAN_WERNER] = l_band_lw; // band luminosity //
#endif

#if defined(RT_CHEM_PHOTOION)   /* Hydrogen and Helium ionizing bands */
    SET_ACTIVE_RT_CHECK();
    double l_ion = particle_ionizing_luminosity_in_cgs(i) / UNIT_LUM_IN_CGS; /* calculate ionizing flux based on actual stellar or BH physics */
    lum[RT_FREQ_BIN_H0] = l_ion; // default to putting everything into a single band //
#if defined(RT_PHOTOION_MULTIFREQUENCY)
    int i_vec[4] = {RT_FREQ_BIN_H0, RT_FREQ_BIN_He0, RT_FREQ_BIN_He1, RT_FREQ_BIN_He2}; // these will all be the same if not using multi-frequency module //
    int k; for(k=0;k<4;k++) {lum[i_vec[k]] = l_ion * rt_ion_precalc_stellar_luminosity_fraction[i_vec[k]];} // assign flux appropriately according to pre-tabulated result //
#endif
#endif

#if defined(RT_HARD_XRAY) || defined(RT_SOFT_XRAY) /* soft and hard X-rays for e.g. Compton heating by X-ray binaries */
    SET_ACTIVE_RT_CHECK();
    double L_HMXBs=0; if(star_age > 0.01) {L_HMXBs = 1.0e29 / (star_age*star_age);}
#if defined(RT_SOFT_XRAY)
    lum[RT_FREQ_BIN_SOFT_XRAY] = (8.2e27 + 0.4*L_HMXBs) * m_sol / UNIT_LUM_IN_CGS; // LMXBs+HMXBs
#endif
#if defined(RT_HARD_XRAY)
    lum[RT_FREQ_BIN_HARD_XRAY] = (6.3e27 + 0.6*L_HMXBs) * m_sol / UNIT_LUM_IN_CGS; // LMXBs+HMXBs
#endif
#endif
    
#endif
    return active_check;
}



/* subroutine for 'rt_get_source_luminosity', with identical variables, for cases where the radiation
   represents flux from sink particles representing AGN, i.e. the sort of thing which would be used in galaxy simulations.
*/
int rt_get_lum_band_agn(int i, int mode, double *lum)
{
    if(P[i].Type != 5) {return 0;} // only go forward for BH-type particles
    int active_check = 0; // default to inactive //
#if defined(BLACK_HOLES)
    double l_bol = bh_lum_bol(P[i].BH_Mdot,P[i].Mass,i); if(l_bol <= 0) {return 0;} // no accretion luminosity -- no point in going further!
    // corrections below follow  Shen, PFH, et al. 2020 to account for alpha-ox and template spectrum to get AGN set in different bands as a function of bolometric luminosity. functional form very similar to Hopkins, Richards, & Hernquist 2007, but updated values. //
    double lbol_lsun = l_bol * UNIT_LUM_IN_SOLAR, R_opt_xr; // luminosity in physical code units //
    double f_xr_0=0.0461795, R_xr_opt = pow(lbol_lsun/1.e10,0.026) / (0.0455713 + 0.140974*pow(lbol_lsun/1.e10,0.304)), Rfxr=R_xr_opt*f_xr_0; // x-ray to optical ratio normalized to its value at Lbol=1e13 solar
    if(Rfxr > 0.5) {R_xr_opt /= pow(1.+Rfxr*Rfxr*Rfxr*Rfxr, 0.25);} // this just prevents unphysical divergences
    R_opt_xr = (1.-R_xr_opt*f_xr_0) / (1.-f_xr_0); // this corrects the IR/optical/UV portion of the spectrum
    
#if defined(RT_INFRARED) /* special mid-through-far infrared band, which includes IR radiation temperature evolution */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_INFRARED] = 0.273 * R_opt_xr * l_bol;
#endif
#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models; from 0.41-3.4 eV */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_OPTICAL_NIR] = 0.181 * R_opt_xr * l_bol;
#endif
#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models; from 3.4-8 eV */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_NUV] = 0.141 * R_opt_xr * l_bol;
#endif
#ifdef RT_PHOTOELECTRIC /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_PHOTOELECTRIC] = 0.117 * R_opt_xr * l_bol; // broad band here [note can 2x-count with LW because that is a sub-band, but include it b/c need to total for dust PE heating
#endif
#ifdef RT_LYMAN_WERNER  /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_LYMAN_WERNER] = 0.0443 * R_opt_xr * l_bol;
#endif
#if defined(RT_CHEM_PHOTOION)   /* Hydrogen and Helium ionizing bands */
    SET_ACTIVE_RT_CHECK();
#if defined(RT_PHOTOION_MULTIFREQUENCY)
    lum[RT_FREQ_BIN_H0]  = 0.1130 * R_opt_xr * l_bol; // total ionizing flux: 13.6-24.6
    lum[RT_FREQ_BIN_He0] = 0.0820 * R_opt_xr * l_bol; // total ionizing flux: 24.6-54.5
    lum[RT_FREQ_BIN_He1] = 0.0111 * R_opt_xr * l_bol; // total ionizing flux: 54.5-70
    lum[RT_FREQ_BIN_He2] = 0.0243 * R_opt_xr * l_bol; // total ionizing flux: 70-500
#else
    lum[RT_FREQ_BIN_H0] = 0.230 * R_opt_xr * l_bol; // total ionizing flux
#endif
#endif
#if defined(RT_SOFT_XRAY) /* soft x-ray 0.5-2 keV band, for compton heating */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_SOFT_XRAY] = 0.00803 * R_xr_opt * l_bol;
#endif
#if defined(RT_HARD_XRAY) /* hard x-ray 2-10+ keV band, for compton heating; since used for that we include some higher-frequence radiation as well */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_HARD_XRAY] = 2.33 * 0.0113 * R_xr_opt * l_bol; // [2.33 factor is extrapolating to include -ultra-hard- X-rays beyond 10keV, useful for some Compton heating estimates]
#endif
    /* note: once account for 2x-counting of LW and PE bands above, this adds up to almost the entire Lbol, but fraction ~ 0.0236 * R_xr_opt * l_bol remains, divided between radio [radio-quiet agn here] and gamma-rays */


#endif
    return active_check;
}



/* subroutine for 'rt_get_source_luminosity', with identical variables, for cases where the radiation
   represents flux from individual sink or star particles representing individual stars or protostars,
   i.e. the sort of thing which would be used in star or planet formation simulations.
*/
int rt_get_lum_band_singlestar(int i, int mode, double *lum)
{
    if(P[i].Type < 4) {return 0;} // only go forward with star or sink-type particles
    int active_check = 0, k; // default to inactive //
    
#if defined(RT_INFRARED) /* special mid-through-far infrared band, which includes IR radiation temperature evolution */
    SET_ACTIVE_RT_CHECK(); k=RT_FREQ_BIN_INFRARED; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]);
#endif
#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models; from 0.41-3.4 eV */
    SET_ACTIVE_RT_CHECK(); k=RT_FREQ_BIN_OPTICAL_NIR; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]);
#endif
#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models; from 3.4-8 eV */
    SET_ACTIVE_RT_CHECK(); k=RT_FREQ_BIN_NUV; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]);
#endif
#ifdef RT_PHOTOELECTRIC /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK(); k=RT_FREQ_BIN_PHOTOELECTRIC; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]); // broad band here [note can 2x-count with LW because that is a sub-band, but include it b/c need to total for dust PE heating
#endif
#ifdef RT_LYMAN_WERNER  /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK(); k=RT_FREQ_BIN_LYMAN_WERNER; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]);
#endif
#if defined(GALSF_FB_FIRE_RT_LONGRANGE) && defined(RADTRANSFER) /* set of FIRE default bands, if used here for stars as well, though currently not cross-linked with some of the other physics */
    SET_ACTIVE_RT_CHECK(); k=RT_FREQ_BIN_FIRE_UV; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]);
    k=RT_FREQ_BIN_FIRE_OPT; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]);
    k=RT_FREQ_BIN_FIRE_IR; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]);
#endif
#if defined(RT_CHEM_PHOTOION)   /* Hydrogen and Helium ionizing bands */
    SET_ACTIVE_RT_CHECK();
#if defined(RT_PHOTOION_MULTIFREQUENCY)
    int i_vec[4] = {RT_FREQ_BIN_H0, RT_FREQ_BIN_He0, RT_FREQ_BIN_He1, RT_FREQ_BIN_He2}; // these will all be the same if not using multi-frequency module //
    for(k=0;k<4;k++) {lum[i_vec[k]] = stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[i_vec[k]],All.RHD_bins_nu_max_ev[i_vec[k]]);} // integrate between band boundaries, defined in global 'nu' in eV
#else
    SET_ACTIVE_RT_CHECK(); k=RT_FREQ_BIN_H0; lum[k]=stellar_lum_in_band(i,All.RHD_bins_nu_min_ev[k],All.RHD_bins_nu_max_ev[k]); // total ionizing flux
#ifdef RT_STARBENCH_TEST
    lum[RT_FREQ_BIN_H0] = 1e49 * (rt_nu_eff_eV[RT_FREQ_BIN_H0]*ELECTRONVOLT_IN_ERGS) / UNIT_LUM_IN_CGS;
#endif
#endif
#endif

#ifdef RT_SOFT_XRAY
    /* currently assume zero here, need to add function here if desired from XRBs or coronal activity, b/c model assumes a thermal spectrum which will give null */
#endif
#ifdef RT_HARD_XRAY
    /* currently assume zero here, need to add function here if desired from XRBs or coronal activity, b/c model assumes a thermal spectrum which will give null */
#endif
#ifdef RT_FREEFREE
    /* negligible free-free emissivity from stars here */
#endif
#ifdef RT_GENERIC_USER_FREQ
    /* code whatever is desired */
#endif

    return active_check;
}




/* this initializes the list of the effective min and max frequencies associated with each waveband, to be used throughout the code */
void rt_define_effective_frequencies_in_bands(void)
{
    /* initialize the table that contains the effective wavelengths of all the different bansd we are actually evolving */
    int k; double rhd_bins_nu_min_ev[N_RT_FREQ_BINS], rhd_bins_nu_max_ev[N_RT_FREQ_BINS]; for(k=0;k<N_RT_FREQ_BINS;k++) {rhd_bins_nu_min_ev[k]=0; rhd_bins_nu_max_ev[k]=MAX_REAL_NUMBER;}
#ifdef RT_CHEM_PHOTOION
#if defined(RT_PHOTOION_MULTIFREQUENCY)
    int i_vec[4] = {RT_FREQ_BIN_H0, RT_FREQ_BIN_He0, RT_FREQ_BIN_He1, RT_FREQ_BIN_He2};
    rhd_bins_nu_min_ev[i_vec[3]]=rt_ion_nu_min[i_vec[3]]; rhd_bins_nu_max_ev[i_vec[3]]=500; for(k=0;k<3;k++) {rhd_bins_nu_min_ev[i_vec[k]]=rt_ion_nu_min[i_vec[k]]; rhd_bins_nu_max_ev[i_vec[k]]=rt_ion_nu_min[i_vec[k+1]];}
#else
    k=RT_FREQ_BIN_H0; rhd_bins_nu_min_ev[k]=13.6; rhd_bins_nu_max_ev[k]=500;
#endif
#endif
#ifdef RT_SOFT_XRAY
    k=RT_FREQ_BIN_SOFT_XRAY; rhd_bins_nu_min_ev[k]=500; rhd_bins_nu_max_ev[k]=2000;
#endif
#ifdef RT_HARD_XRAY
    k=RT_FREQ_BIN_HARD_XRAY; rhd_bins_nu_min_ev[k]=2000; rhd_bins_nu_max_ev[k]=10000;
#endif
#ifdef RT_PHOTOELECTRIC
    k=RT_FREQ_BIN_PHOTOELECTRIC; rhd_bins_nu_min_ev[k]=8; rhd_bins_nu_max_ev[k]=13.6;
#endif
#ifdef RT_LYMAN_WERNER
    k=RT_FREQ_BIN_LYMAN_WERNER; rhd_bins_nu_min_ev[k]=11.2; rhd_bins_nu_max_ev[k]=13.6;
#endif
#ifdef RT_NUV
    k=RT_FREQ_BIN_NUV; rhd_bins_nu_min_ev[k]=3.444; rhd_bins_nu_max_ev[k]=8.;
#endif
#ifdef RT_OPTICAL_NIR
    k=RT_FREQ_BIN_OPTICAL_NIR; rhd_bins_nu_min_ev[k]=0.4133; rhd_bins_nu_max_ev[k]=3.444;
#endif
#ifdef RT_GENERIC_USER_FREQ
    k=RT_FREQ_BIN_GENERIC_USER_FREQ; rhd_bins_nu_min_ev[k]=0; rhd_bins_nu_max_ev[k]=MAX_REAL_NUMBER;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    k=RT_FREQ_BIN_FIRE_UV; rhd_bins_nu_min_ev[k]=3.444; rhd_bins_nu_max_ev[k]=13.6;
    k=RT_FREQ_BIN_FIRE_OPT; rhd_bins_nu_min_ev[k]=0.365; rhd_bins_nu_max_ev[k]=3.444;
    k=RT_FREQ_BIN_FIRE_IR; rhd_bins_nu_min_ev[k]=0.01; rhd_bins_nu_max_ev[k]=0.365;
#endif
#ifdef RT_INFRARED
    k=RT_FREQ_BIN_INFRARED; rhd_bins_nu_min_ev[k]=0.001; rhd_bins_nu_max_ev[k]=0.4133;
#endif
#ifdef RT_FREEFREE
    k=RT_FREQ_BIN_FREEFREE; rhd_bins_nu_min_ev[k]=0; rhd_bins_nu_max_ev[k]=MAX_REAL_NUMBER;
#endif
    for(k=0;k<N_RT_FREQ_BINS;k++) {All.RHD_bins_nu_min_ev[k]=rhd_bins_nu_min_ev[k]; All.RHD_bins_nu_max_ev[k]=rhd_bins_nu_max_ev[k];}
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
    /* should be equal to (c_reduced * Kappa_opacity * rho) */
    return (C_LIGHT_CODE_REDUCED) * rt_absorb_frac_albedo(i,k_freq) * (rt_kappa(i,k_freq) * SphP[i].Density*All.cf_a3inv);
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
    int k_freq, k; double n_flux_j[3], fmag_j, V_j_inv = SphP[j].Density / P[j].Mass;
    for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        n_flux_j[0]=n_flux_j[1]=n_flux_j[2]=0;
        double flux_vol[3]; for(k=0;k<3;k++) {flux_vol[k] = SphP[j].Rad_Flux[k_freq][k] * V_j_inv;}
        fmag_j = 0; for(k=0;k<3;k++) {fmag_j += flux_vol[k]*flux_vol[k];}
        if(fmag_j <= 0) {fmag_j=0;} else {fmag_j=sqrt(fmag_j); for(k=0;k<3;k++) {n_flux_j[k]=flux_vol[k]/fmag_j;}}
        double f_chifac = fmag_j / (MIN_REAL_NUMBER + C_LIGHT_CODE_REDUCED * SphP[j].Rad_E_gamma[k_freq] * V_j_inv);
        if(f_chifac < 0) {f_chifac=0;}
        if(fmag_j <= 0) {f_chifac = 0;}
        // restrict values of f_chifac to physical range.
        double f_min = 0.01, f_max = 0.9999;
        if(f_chifac < f_min) {f_chifac = f_min;}
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
    double E_abs_tot_toIR = 0;/* energy absorbed in other bands is transfered to IR, by default: track it here */
    double Rad_E_gamma_tot = 0; // dust temperature defined by total radiation energy density //
    {int j; for(j=0;j<N_RT_FREQ_BINS;j++) {Rad_E_gamma_tot += SphP[i].Rad_E_gamma[j];}}
    double a_rad_inverse=C_LIGHT_CGS/(4.*5.67e-5), vol_inv_phys=(SphP[i].Density*All.cf_a3inv/P[i].Mass), u_gamma = Rad_E_gamma_tot * vol_inv_phys * UNIT_PRESSURE_IN_CGS; // photon energy density in CGS //
    double Dust_Temperature_4 = u_gamma * a_rad_inverse; // estimated effective temperature of local rad field in equilibrium with dust emission. note that for our definitions, rad energy density has its 'normal' value independent of RSOL, so Tdust should as well; emission -and- absorption are both lower by a factor of RSOL, but these cancel in the Tdust4 here //
#if !defined(COOLING) // if cooling is active, don't reset this here, because it needs to include the gas coupling term which will be self-consistently calculated there
    SphP[i].Dust_Temperature = sqrt(sqrt(Dust_Temperature_4)); // just set this to the local radiation equilibrium temperature
#endif
    double T_min = get_min_allowed_dustIRrad_temperature();
    if(SphP[i].Dust_Temperature <= T_min) {SphP[i].Dust_Temperature = T_min;} // dust temperature shouldn't be below CMB
    double IRBand_opacity_fraction_from_gas_absorption = 0; // needed below to know what fraction is immediately re-radiated or not
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
            double e0, dt_e_gamma_band=0, total_de_dt=0, a0_abs = -rt_absorption_rate(i,kf);
#if defined(RT_EVOLVE_INTENSITIES)
            if(mode==0) {e0 = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity[kf][k_angle];} else {e0 = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity_Pred[kf][k_angle];}
            dt_e_gamma_band = RT_INTENSITY_BINS_DOMEGA*SphP[i].Dt_Rad_Intensity[kf][k_angle];
#else
            if(mode==0) {e0 = SphP[i].Rad_E_gamma[kf];} else {e0 = SphP[i].Rad_E_gamma_Pred[kf];}
            dt_e_gamma_band = SphP[i].Dt_Rad_E_gamma[kf];
#endif
#ifdef RT_COMOVING
            double ET_dotdot_GradVcom = SphP[i].ET[kf][0]*SphP[i].Gradients.Velocity[0][0] + SphP[i].ET[kf][1]*SphP[i].Gradients.Velocity[1][1] + SphP[i].ET[kf][2]*SphP[i].Gradients.Velocity[2][2]
                + SphP[i].ET[kf][3]*(SphP[i].Gradients.Velocity[0][1]+SphP[i].Gradients.Velocity[1][0]) + SphP[i].ET[kf][4]*(SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2]) + SphP[i].ET[kf][5]*(SphP[i].Gradients.Velocity[0][2]+SphP[i].Gradients.Velocity[2][0]);
            double VolP_dotdot_GradV = e0 * ET_dotdot_GradVcom * All.cf_a2inv; // convert to physical units and multiply by radiation energy density to get into appropriate units
            dt_e_gamma_band += (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE) * (-VolP_dotdot_GradV); // account for RSOL term here as usual
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
                    double dE_abs = -e0 * (1. - exp(a0_abs*dt_entr)); // change in energy from absorption
                    double T_max = DMAX(SphP[i].Radiation_Temperature , dE_fac / dTE_fac); 
                    double rfac=1; if(dE_fac < -0.5*(e0+dE_abs)) {rfac=fabs(0.5*(e0+dE_abs))/fabs(dE_fac);} else {if(dE_fac > 0.5*e0) {rfac=0.5*e0/dE_fac;}}
                    dE_fac*=rfac; dTE_fac*=rfac; // limit temperature change from advection to prevent spurious divergences
                    
                    SphP[i].Radiation_Temperature = (e0 + dE_fac) / (MIN_REAL_NUMBER + DMAX(0., e0 / SphP[i].Radiation_Temperature + dTE_fac));
                    SphP[i].Radiation_Temperature = DMIN(SphP[i].Radiation_Temperature, T_max);
                    a0_abs = -rt_absorption_rate(i,kf); // update absorption rate using the new radiation temperature //
                }
                double total_absorption_rate = E_abs_tot_toIR + fabs(a0_abs)*e0; // add the summed absorption and equate to dust emission //
#ifdef COOLING  // we account for gas-dust coupling as an additional heat source to be radiated away
		        double temp = get_temperature(i);
		        double nHcgs = HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * SphP[i].Density * All.cf_a3inv / PROTONMASS_CGS;

                /* 
                Here we are splitting the re-emission by dust (performed here) and the gas-dust-radiation energy transfer handled in the cooling solver;
                need to know the dust temperature to get the re-radiated radiation temperature. OK as long as the opacity doesn't
                change dramatically in a single timestep, otherwise have to couple all matter-radiation source terms implicitly. */    
		        SphP[i].Dust_Temperature = rt_eqm_dust_temp(i, temp, total_absorption_rate * vol_inv_phys / RT_SPEEDOFLIGHT_REDUCTION);
#else
                SphP[i].Dust_Temperature = rt_eqm_dust_temp(i, 0, total_absorption_rate * vol_inv_phys / RT_SPEEDOFLIGHT_REDUCTION); // Calling with T=0 will account for dust absorption only
#endif
                if(SphP[i].Dust_Temperature < T_min){SphP[i].Dust_Temperature = T_min;}
                double Tdust_eff = SphP[i].Dust_Temperature, Trad_eff = SphP[i].Radiation_Temperature;
                double kappa_gas = rt_kappa_adaptive_IR_band(i,Tdust_eff,Trad_eff,-1,-1), kappa_total = rt_kappa_adaptive_IR_band(i,Tdust_eff,Trad_eff,0,0);
                IRBand_opacity_fraction_from_gas_absorption = kappa_gas / (kappa_total + MIN_REAL_NUMBER); /* gas absorption opacity only, relative to total opacity (all sources+scattering) */
                double total_emission_rate = total_absorption_rate * (1.-IRBand_opacity_fraction_from_gas_absorption) + SphP[i].Rad_Je[kf]; /* we will re-radiate this much because the component due to gas-dust coupling is accounted for in the cooling loop */
                total_de_dt = E_abs_tot_toIR + SphP[i].Rad_Je[kf] + dt_e_gamma_band;
                if((mode==0) && (Tdust_eff < MAX_DUST_TEMP)) // only update temperatures on kick operations and Tdust is meaningful //
                {
                    /* dust absorption and re-emission brings T_rad towards T_dust: */
                    double dE_abs_IR = -e0 * (1. - exp(a0_abs*dt_entr)); /* change in energy from absorption */
                    double T_max = DMAX(SphP[i].Radiation_Temperature , Tdust_eff); /* should not exceed either initial temperature */
                    SphP[i].Radiation_Temperature = (e0 + dE_abs_IR + total_emission_rate*dt_entr) / (MIN_REAL_NUMBER + (e0 + dE_abs_IR) / SphP[i].Radiation_Temperature + total_emission_rate*dt_entr / Tdust_eff);
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
            if(fabs(a0_abs)*dt_entr > 50.) {a0_abs *= 50./(fabs(a0_abs)*dt_entr);}
            double abs_0 = DMAX(0,fabs(a0_abs)*dt_entr); double slabfac = slab_averaging_function(abs_0); double e_abs_0=exp(-abs_0); if(abs_0>100.) {e_abs_0=0;}
            /* since we're taking exponentials and inverses of some large numbers here, need to be careful not to let floating point errors cause nan's */
            if((dt_entr <= 0.)||(a0_abs >= 0.)||(abs_0 <= 0.)) {abs_0=0.; slabfac=e_abs_0=1.;} else {if(abs_0 < 1.e-5) {slabfac=1.-0.5*abs_0; e_abs_0 = 1.-abs_0;} else {if(abs_0 > 100.) {slabfac = 1./abs_0; e_abs_0 = 0.;}}}
            double e0_postabs = e0*e_abs_0, de_postabs = total_de_dt * dt_entr * slabfac, f_min = 0.01;
            if(e0_postabs+de_postabs < f_min*e0_postabs) {slabfac *= fabs((1.-f_min)*e0_postabs)/(fabs(de_postabs)+MIN_REAL_NUMBER);}
            
            double ef = e0 * e_abs_0 + total_de_dt * dt_entr * slabfac; // gives exact solution for dE/dt = -E*abs + de , the 'reduction factor' appropriately suppresses the source term //
#ifdef RT_INFRARED
            if(isnan(ef)) {PRINT_WARNING("\n ef energy prediction is NaN for cell-ID=%llu, e0=%g e_abs_0=%g abs_0=%g a0_abs=%g total_de_dt=%g dt_entr=%g slabfac=%g Trad=%g Tdust=%g\n", (unsigned long long) P[i].ID,e0, e_abs_0,abs_0, a0_abs, total_de_dt,dt_entr,slabfac,SphP[i].Radiation_Temperature,SphP[i].Dust_Temperature);}
#else
            if(isnan(ef)) {PRINT_WARNING("\n ef energy prediction is NaN for cell-ID=%llu, e0=%g e_abs_0=%g abs_0=%g a0_abs=%g total_de_dt=%g dt_entr=%g slabfac=%g\n", (unsigned long long) P[i].ID,e0, e_abs_0,abs_0, a0_abs, total_de_dt,dt_entr,slabfac);}
#endif
            if(ef < 0) {ef=0;}
            double de_abs = e0 + total_de_dt * dt_entr - ef; // energy removed by absorption alone
            double de_emission_minus_absorption = (ef - DMAX(0, (e0 + dt_e_gamma_band * dt_entr * slabfac))); // total change, relative to what we would get with just advection (positive = net energy increase in the gas)
            if((dt_entr <= 0) || (de_abs <= 0)) {de_abs = 0;}
            
#if defined(RT_RAD_PRESSURE_FORCES) && defined(RT_COMPGRAD_EDDINGTON_TENSOR) && !defined(RT_EVOLVE_FLUX) && !defined(RT_RADPRESSURE_IN_HYDRO)
            // for OTVET/FLD methods, need to apply radiation pressure term here so can limit this b/c just based on a gradient which is not flux-limited [as in hydro operators] //
            {
                double radacc[3]={0}, rmag=0, vel_i[3], L_particle = Get_Particle_Size(i)*All.cf_atime; // particle effective size/slab thickness
                double Sigma_particle = P[i].Mass / (M_PI*L_particle*L_particle), abs_per_kappa_dt = C_LIGHT_CODE_REDUCED * (SphP[i].Density*All.cf_a3inv) * dt_entr; // effective surface density through particle & fractional absorption over timestep
                double f_kappa_abs = rt_absorb_frac_albedo(i,kf); // get albedo, we'll need this below
                double slabfac_rp=1; if(check_if_absorbed_photons_can_be_reemitted_into_same_band(kf)==0) {slabfac_rp=slab_averaging_function(f_kappa_abs*SphP[i].Rad_Kappa[kf]*Sigma_particle) * slab_averaging_function(f_kappa_abs*SphP[i].Rad_Kappa[kf]*abs_per_kappa_dt);} // reduction factor for absorption over dt
                int kx; for(kx=0;kx<3;kx++)
                {
                    radacc[kx] = -dt_entr * slabfac_rp * return_flux_limiter(i,kf) * (SphP[i].Gradients.Rad_E_gamma_ET[kf][kx] / SphP[i].Density) / All.cf_atime; // naive radiation-pressure calc for FLD methods [physical units]
                    rmag += radacc[kx]*radacc[kx]; // compute magnitude
                    if(mode==0) {vel_i[kx]=(C_LIGHT_CODE_REDUCED/C_LIGHT_CODE)*P[i].Vel[kx]/All.cf_atime;} else {vel_i[kx]=(C_LIGHT_CODE_REDUCED/C_LIGHT_CODE)*SphP[i].VelPred[kx]/All.cf_atime;} // [for comoving] note this is the 'effective' u appearing in the RHD equations for an RSOL, care needed with these factors!
                }
                if(rmag > 0)
                {
                    rmag = sqrt(rmag); for(kx=0;kx<3;kx++) {radacc[kx] /= rmag;} // normalize
                    double rmag_max = de_abs / (P[i].Mass * C_LIGHT_CODE_REDUCED * (MIN_REAL_NUMBER + f_kappa_abs)); // limit magnitude to absorbed photon momentum
                    if(check_if_absorbed_photons_can_be_reemitted_into_same_band(kf)==0) {if(rmag > rmag_max) {rmag=rmag_max;}}
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
                    double d_egy_rad = (2.*f_kappa_abs-1.)*work_band , d_egy_int = -2.*f_kappa_abs*work_band * (C_LIGHT_CODE/C_LIGHT_CODE_REDUCED); // correct for rsol factor above which reduced vel_i by rsol; -only- add back this term for gas
                    if(mode==0) {SphP[i].InternalEnergy += d_egy_int;} else {SphP[i].InternalEnergyPred += d_egy_int;}
#if defined(RT_EVOLVE_INTENSITIES)
                    {int k_q; for(k_q=0;k_q<N_RT_INTENSITY_BINS;k_q++) {if(mode==0) {SphP[i].Rad_Intensity[kf][k_q]+=d_egy_rad/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[kf][k_q]+=d_egy_rad/RT_INTENSITY_BINS_DOMEGA;}}}
#else
                    if(mode==0) {SphP[i].Rad_E_gamma[kf]+=d_egy_rad;} else {SphP[i].Rad_E_gamma_Pred[kf]+=d_egy_rad;}
#endif
                }
            }
#endif
            
            int donation_target_bin = rt_get_donation_target_bin(kf); // frequency into which the photons will be deposited, if any //
#ifdef RT_INFRARED
            if((donation_target_bin == RT_FREQ_BIN_INFRARED) && (kf != RT_FREQ_BIN_INFRARED)) {E_abs_tot_toIR += de_abs/(MIN_REAL_NUMBER + dt_entr);} /* donor bin is yourself in the IR - some self-absorption is re-emitted, but this is handled explicitly below, so don't need to include it in sum here */
            if(kf==RT_FREQ_BIN_INFRARED) {
#ifdef COOLING
                ef += de_abs*(1.-IRBand_opacity_fraction_from_gas_absorption); /* update: assume a fraction de_abs * IRBand_opacity_fraction_from_gas_absorption is absorbed by the gas, which will not be instantly re-emitted here, but later in the cooling subroutines */
                if(mode==0) {SphP[i].DtInternalEnergy += (de_abs * IRBand_opacity_fraction_from_gas_absorption) / ((MIN_REAL_NUMBER + dt_entr) * P[i].Mass);} /* this fraction absorbed by gas goes into a heating rate which can be balanced implicitly in the cooling function later */
#else
                ef = e0 + total_de_dt * dt_entr; // previous version: assumes all self-absorption is re-emitted
#endif
            } /* donor bin is yourself in the IR - just need to decide what to do with the photons */
#endif
            // isotropically re-emit the donated radiation into the target bin[s] //
#if defined(RT_EVOLVE_INTENSITIES)
            // this is the leading-order (isotropic) emission-absorption step, i.e. the psi_a * (j_e - I) term in the intensity equation. solved by the methods above to deal generically with stiff emission-absorption problems, re-used below if needed //
            if(donation_target_bin >= 0) {int k_q; for(k_q=0;k_q<N_RT_INTENSITY_BINS;k_q++) {if(mode==0) {SphP[i].Rad_Intensity[donation_target_bin][k_q] += de_abs/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[donation_target_bin][k_q] += de_abs/RT_INTENSITY_BINS_DOMEGA;}}}
            if(ef < 0) {ef=0;}
            if(mode==0) {SphP[i].Rad_Intensity[kf][k_angle] = ef/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[kf][k_angle] = ef/RT_INTENSITY_BINS_DOMEGA;}
#else
            if(donation_target_bin >= 0) {if(mode==0) {SphP[i].Rad_E_gamma[donation_target_bin] += de_abs;} else {SphP[i].Rad_E_gamma_Pred[donation_target_bin] += de_abs;}}
            if(ef < 0) {ef=0;}

            if(mode==0) {SphP[i].Rad_E_gamma[kf] = ef;} else {SphP[i].Rad_E_gamma_Pred[kf] = ef;}
#endif

#if defined(RT_EVOLVE_FLUX)
            int k_dir; double f_mag=0, E_rad_forflux=0, vdot_h[3]={0}, vel_i[3]={0}, DeltaFluxEff[3]={0}, rho=SphP[i].Density*All.cf_a3inv; E_rad_forflux=0.5*(e0+ef); // use energy density averaged over this update for the operation below
            for(k_dir=0;k_dir<3;k_dir++) {if(mode==0) {vel_i[k_dir]=RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS*P[i].Vel[k_dir]/All.cf_atime;} else {vel_i[k_dir]=RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS*SphP[i].VelPred[k_dir]/All.cf_atime;}} // need gas velocity at this time [effective v - note RSOL terms]
            double teqm_inv = SphP[i].Rad_Kappa[kf] * rho * C_LIGHT_CODE_REDUCED + MIN_REAL_NUMBER; // physical code units of 1/time, defines characteristic timescale for coming to equilibrium flux. see notes for CR second-order module for details. //
            eddington_tensor_dot_vector(SphP[i].ET[kf], vel_i, vdot_h); // calculate volume integral of scattering coefficient t_inv * (gas_vel . [e_rad*I + P_rad_tensor]), which gives an additional time-derivative term. this is the P term //
            for(k_dir=0;k_dir<3;k_dir++) {vdot_h[k_dir] = E_rad_forflux * (vel_i[k_dir] + vdot_h[k_dir]);} // and this is the eI term, multiply both by radiation energy to use in this step //
#ifdef RT_COMPGRAD_EDDINGTON_TENSOR // definitely favor this for greater accuracy and reduced noise //
            for(k_dir=0;k_dir<3;k_dir++) {DeltaFluxEff[k_dir] -= (P[i].Mass/rho) * (C_LIGHT_CODE_REDUCED*C_LIGHT_CODE_REDUCED/teqm_inv) * SphP[i].Gradients.Rad_E_gamma_ET[kf][k_dir]*All.cf_a3inv/All.cf_atime;} // here we compute the nabla.pressure_gradient_tensor term from gradients directly, and use this in the next step after multiplying the flux equation by (tilde[c]^2/dt_eqm_inv) and working in dimensionless time units
#else
            for(k_dir=0;k_dir<3;k_dir++) {DeltaFluxEff[k_dir] += (SphP[i].Dt_Rad_Flux[kf][k_dir]/teqm_inv);} // the nabla.pressure_gradient_tensor is computed in the finite-volume solver, here
#endif
            for(k_dir=0;k_dir<3;k_dir++) {DeltaFluxEff[k_dir] += vdot_h[k_dir];} // add the 'enthalpy advection' term here, vdot_h = Erad v.(e*I + P_rad)

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
                    f_mag=sqrt(f_mag); double fmag_max = C_LIGHT_CODE_REDUCED * ef; // maximum flux should be optically-thin limit: e_gamma*c: here allow some tolerance for numerical leapfrogging in timestepping. should be the RSOL here, although in principle equations can allow exceeding this if we have reached equilibrium, it really violates the M1 closure assumptions. see discussion in Skinner+Ostriker 2013 or Levermore et al. 1984
                    if(f_mag > fmag_max) {for(k_dir=0;k_dir<3;k_dir++) {if(mode==0) {SphP[i].Rad_Flux[kf][k_dir] *= fmag_max/f_mag;} else {SphP[i].Rad_Flux_Pred[kf][k_dir] *= fmag_max/f_mag;}}}
#if defined(GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION)
                    if(P[i].Pos[2]<=0.1) {if(mode==0) {SphP[i].Rad_Flux[kf][0]=SphP[i].Rad_Flux[kf][1]=0; SphP[i].Rad_Flux[kf][2]=fmag_max;} else {SphP[i].Rad_Flux_Pred[kf][0]=SphP[i].Rad_Flux_Pred[kf][1]=0; SphP[i].Rad_Flux_Pred[kf][2]=fmag_max;}}
#endif
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
    if(dt_entr > 0) { // none of this is worth doing if we don't have a finite timestep here
    for(kf=0;kf<N_RT_FREQ_BINS;kf++)
    {
        int k,k_om; double rho=SphP[i].Density*All.cf_a3inv, ceff=C_LIGHT_CODE_REDUCED, ctrue=C_LIGHT_CODE, teq_inv=SphP[i].Rad_Kappa[kf]*rho*ceff, beta[3], f_a=rt_absorb_frac_albedo(i,kf), f_s=1.-f_a, b_dot_n[N_RT_INTENSITY_BINS]={0}, beta_2=0.;
        int n_iter = 1 + (int)(DMIN(DMAX(4. , dt_entr/teq_inv), 1000.)); // number of iterations to subcycle everything below //
        double dt=dt_entr/n_iter, tau=dt*teq_inv, i0[N_RT_INTENSITY_BINS]={0}, invfourpi=1./(4.*M_PI), J, b_dot_H, b2_dot_K; int i_iter;
        for(i_iter=0; i_iter<n_iter; i_iter++)
        {
            double egy_0=0,flux_0[3]={0},egy_f=0,flux_f[3]={0}; // compute total change over sub-cycle, to update gas properties
            // load all the gas and intensity properties we need [all can change on the subcycle so some re-computing here]
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if(mode==0) {i0[k_om] = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity[kf][k_om];} else {i0[k_om] = RT_INTENSITY_BINS_DOMEGA*SphP[i].Rad_Intensity_Pred[kf][k_om];}}
            for(k=0;k<3;k++) {if(mode==0) {beta[k]=P[i].Vel[k]/(All.cf_atime*ctrue);} else {beta[k]=SphP[i].VelPred[k]/(All.cf_atime*ctrue);}} // need gas velocity at this time; with equations written this way, the 'beta' term is the -true- beta, so we have to use the true SOL
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {b_dot_n[k_om]=0; for(k=0;k<3;k++) {b_dot_n[k_om]+=All.Rad_Intensity_Direction[k_om][k]*beta[k];}}
            beta_2=0; for(k=0;k<3;k++) {beta_2+=beta[k]*beta[k];}
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {egy_0+=i0[k_om]; for(k=0;k<3;k++) {flux_0[k]+=All.Rad_Intensity_Direction[k_om][k]*i0[k_om];}}
            J=0,b_dot_H=0,b2_dot_K=0; for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {J+=i0[k_om]*invfourpi; b_dot_H+=b_dot_n[k_om]*i0[k_om]*invfourpi; b2_dot_K+=b_dot_n[k_om]*b_dot_n[k_om]*i0[k_om]*invfourpi;}

            // isotropic terms that change total energy in bin (part of the 'work term' for the photon momentum): this includes the beta.beta*(J+K) and beta.H terms
            double work = (1. * (f_s-f_a)*(beta_2*J + b2_dot_K) - 2.*f_s*b_dot_H) * tau; // will be shared isotropically.
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if((work>0) || (i0[k_om]<=0)) {i0[k_om]+=work;} else {i0[k_om]/=(1-work/i0[k_om]);}} // gaurantees linearized sum is still correct, and symmetric with positive changes, but can't get negative energies. shared isotropically.

            // isotropic scattering term [scattering * (J - I) term in the intensity equation] [recall, our general update for the 'energy term' for absorption and emission above already took care of the psi_a*(j_em - I) term in the intensity equation
            J=0; for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {J+=i0[k_om]*invfourpi;} // prepare to calculate isotropic scattering term
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {i0[k_om] = J + (i0[k_om]-J)*exp(-f_s*tau);} // isotropic scattering conserving total energy over step
            
            // flux 'boost' and 'beaming' terms (go as n.beta). note we replace je -> je-I + I, and use the fact that we have solved already for the time-integral of (psi_a*(je-I)*dt) = de_emission_minus_absorption_saved, which can be re-used here, in average form <psi_a*(je-I)> = dE/dt --> just make sure the units are correct! because we're working in dimensionless units below, we should divide by tau, to be working in the same delta-units here: these are the 3 n.beta * [ psi_a*(j_em-J_nu)*(creduced/c)^2 + (psi_a+psi_s)*J_nu) ] in the intensity equation
            double fboost[N_RT_INTENSITY_BINS], fboost_avg=0, fboost_avg_p=0, fboost_avg_m=0; // calculate flux 'boost' terms
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {fboost[k_om] = 3.*b_dot_n[k_om] * ((de_emission_minus_absorption_saved[kf][k_om]/tau) + ((f_a+f_s)*J)); fboost_avg += fboost[k_om]/N_RT_INTENSITY_BINS;} // pre-calculate to get mean value, will divide out
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {work=(fboost[k_om]-fboost_avg)*tau; if((work>0) || (i0[k_om]<=0)) {fboost[k_om]=work; fboost_avg_p+=fboost[k_om];} else {fboost[k_om]=work/(1.-work/i0[k_om]); fboost_avg_m+=fboost[k_om];}} // zero total energy change at linear order ensured by subtracting out sum here; non-linearization ensures i0 cannot be negative, but does allow second-order dt work term to appear, that's ok for now
            if(fboost_avg_p>0 && fboost_avg_m<0) {double fc=-fboost_avg_m/fboost_avg_p; fboost_avg_m=(1.+fc)/(1.+fc*fc); fboost_avg_p=fc*fboost_avg_m;} else {fboost_avg_m=fboost_avg_p=0;} // // these re-weight to gaurantee the non-linear sum is identically zero while preserving positive-definite behavior
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if(fboost[k_om]>0) {i0[k_om]+=fboost_avg_p*fboost[k_om];} else {i0[k_om]+=fboost_avg_m*fboost[k_om];}} // alright done!
            
            // flux work term, allowed to both do work and be asymmetric so just need to ensure it retains positive-definite intensities: n.beta * (psi_a+psi_s) * I term in intensity equation
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {work=b_dot_n[k_om]*(f_a+f_s)*i0[k_om] * tau; if((work>0) || (i0[k_om]<=0)) {i0[k_om]+=work;} else {i0[k_om]/=(1-work/i0[k_om]);}}

            // ok -now- calculate the net change in momentum and energy, for updating the gas quantities
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {egy_f+=i0[k_om]; for(k=0;k<3;k++) {flux_f[k]+=All.Rad_Intensity_Direction[k_om][k]*i0[k_om];}}
            double dv_gas[3]={0}, ke_gas_0=0, ke_gas_f=0, v0g=0, u0=0;
            for(k=0;k<3;k++) {dv_gas[k] = -(flux_f[k]-flux_0[k])/(ceff*P[i].Mass); v0g=ctrue*beta[k]; ke_gas_0+=(v0g*v0g); ke_gas_f+=(v0g+dv_gas[k])*(v0g+dv_gas[k]);} // note everything is volume-integrated, accounted for above, and we defined flux for convience without the c, so just one power of c here.
            double d_ke_gas = 0.5*(ke_gas_f - ke_gas_0)*P[i].Mass, de_gas=-(ctrue/ceff)*(egy_f-egy_0), de_gas_internal=(de_gas-d_ke_gas)/P[i].Mass; // note ctrue/ceff factor here, accounting for rsol difference in gas heating/cooling rates vs RHD
            if(mode==0) {u0=SphP[i].InternalEnergy;} else {u0=SphP[i].InternalEnergyPred;} // for updating gas internal energy (work terms, after subtracting kinetic energy changes)
            if(de_gas_internal<=-0.9*u0) {de_gas_internal = DMIN(de_gas_internal/(1.-de_gas_internal/u0), -0.9*u0);} // just a catch to avoid negative energies (will break energy conservation if you are slamming into it, however!
            
            // assign everything back to the appropriate variables after update
            for(k=0;k<3;k++) {if(mode==0) {P[i].Vel[k] += dv_gas[k]*All.cf_atime;} else {SphP[i].VelPred[k] += dv_gas[k]*All.cf_atime;}} // update gas velocities (radiation pressure forces here)
            if(mode==0) {SphP[i].InternalEnergy += de_gas_internal;} else {SphP[i].InternalEnergyPred += de_gas_internal;} // update gas internal energy (work terms, after subtracting kinetic energy changes)
            for(k_om=0;k_om<N_RT_INTENSITY_BINS;k_om++) {if(mode==0) {SphP[i].Rad_Intensity[kf][k_om] = i0[k_om]/RT_INTENSITY_BINS_DOMEGA;} else {SphP[i].Rad_Intensity_Pred[kf][k_om] = i0[k_om]/RT_INTENSITY_BINS_DOMEGA;}} // update intensities (all of the above)
            SphP[i].Rad_E_gamma[kf]=egy_f; // set this every time this subroutine is called, so it is accessible everywhere else //
        } // loop over iterations
    } // loop over frequencies
    } // finite timestep requirement
#else
    double mom_fac = 1. - RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS * total_erad_emission_minus_absorption / (P[i].Mass * C_LIGHT_CODE_REDUCED*C_LIGHT_CODE_REDUCED); // back-reaction on gas from emission, which is isotropic in the fluid frame but anisotropic in the lab frame. this effect is only important in actually semi-relativistic problems so we use "real" C here, not a RSOL, and match the corresponding term above in the radiation flux equation (if that is evolved explicitly). careful checking-through gives the single termm here, not both
    if(fabs(mom_fac - 1) > 0.1) {printf("WARNING: Large radiation backreaction for cell %d (mom_fac=%g), check the RT solver stability if this is not a relativistic problem.\n",i,mom_fac);}
    {int k_dir; for(k_dir=0;k_dir<3;k_dir++) {if(mode==0) {P[i].Vel[k_dir] *= mom_fac;} else {SphP[i].VelPred[k_dir] *= mom_fac;}}}
#endif

    if(mode > 0) {rt_eddington_update_calculation(i);} /* update the eddington tensor (if we calculate it) as well */

#ifdef RT_ISRF_BACKGROUND
    if(mode==0) {rt_apply_boundary_conditions(i);} /* if we have any special boundary conditions (e.g. fixed ISRF at box edge) apply this here */
#endif    
#endif
}



#endif

#if defined(RT_ISRF_BACKGROUND) && defined(RADTRANSFER)
void rt_apply_boundary_conditions(int i)
{
    double urad[N_RT_FREQ_BINS]; int k, k_dir;
    get_background_isrf_urad(i, urad);
    // if we are within 10% of the box length of the edge:
    if(DMAX(DMAX(P[i].Pos[0],P[i].Pos[1]),P[i].Pos[2]) > 0.9*All.BoxSize || DMIN(DMIN(P[i].Pos[0],P[i].Pos[1]),P[i].Pos[2]) < 0.1*All.BoxSize)
    {
        for(k = 0; k < N_RT_FREQ_BINS; k++)
        {
            SphP[i].Rad_E_gamma[k] = urad[k] * P[i].Mass/(SphP[i].Density * All.cf_a3inv);
#ifdef RT_EVOLVE_FLUX
            for(k_dir = 0; k_dir < 3; k_dir++){SphP[i].Rad_Flux[k][k_dir] = 0;}
#endif
#ifdef RT_INFRARED
            if(k==RT_FREQ_BIN_INFRARED) {
                SphP[i].Radiation_Temperature = background_isrf_cmb_Teff();
                SphP[i].Dust_Temperature = DMIN(All.InitGasTemp,100.);
            }
#endif
        }
    } else {
        for(k = 0; k < N_RT_FREQ_BINS; k++){SphP[i].Rad_E_gamma[k] = DMAX(SphP[i].Rad_E_gamma[k], MIN_REAL_NUMBER);}
    }
}

// computes the background ISRF energy density in code units for in each band, in code units
void get_background_isrf_urad(int i, double *urad){
    int k;
    for(k = 0; k < N_RT_FREQ_BINS; k++)
    {
        urad[k] = MIN_REAL_NUMBER;
#ifdef RT_INFRARED
	double fac_uCMB = 1.;
	fac_uCMB = pow(1+All.RadiationBackgroundRedshift, 4);
        if(k==RT_FREQ_BIN_INFRARED){urad[k] = (All.InterstellarRadiationFieldStrength * 0.39 + fac_uCMB * 0.26) * ELECTRONVOLT_IN_ERGS / UNIT_PRESSURE_IN_CGS;} // 0.33 eV/cm^3 is dust emission peak, 0.26 is CMB - note how this bin actually lumps the two together
#endif
#ifdef RT_OPTICAL_NIR
        if(k==RT_FREQ_BIN_OPTICAL_NIR){urad[k] = All.InterstellarRadiationFieldStrength * 0.54 * ELECTRONVOLT_IN_ERGS / UNIT_PRESSURE_IN_CGS;} // stellar emission
#endif
#ifdef RT_NUV
        if(k==RT_FREQ_BIN_NUV){urad[k] = All.InterstellarRadiationFieldStrength * 0.024 * ELECTRONVOLT_IN_ERGS / UNIT_PRESSURE_IN_CGS;} // stellar emission
#endif
#ifdef RT_PHOTOELECTRIC
        if(k==RT_FREQ_BIN_PHOTOELECTRIC){urad[k] = All.InterstellarRadiationFieldStrength * 1.7 * 3.9e-14 / UNIT_PRESSURE_IN_CGS;} // Draine 1978 value = 1.7 Habing
#endif
    }
}

double background_isrf_cmb_Teff(){
    // Returns the energy-weighted effective temperature of the background ISRF that has equivalent average photon energy to the sum of the ISRF and CMB
    // Necessary because current IR band treatment lumps both radiation fields together
    double urad_ISRF_CGS_eV = All.InterstellarRadiationFieldStrength * 0.39, Trad_ISRF = DMIN(All.InitGasTemp,100.);
    double fac_TCMB= 1.+All.RadiationBackgroundRedshift, fac_uCMB = pow(fac_TCMB,4);
    double urad_CMB_CGS_eV = fac_uCMB * 0.262, Trad_CMB = 2.73 * fac_TCMB;
    return (urad_ISRF_CGS_eV * Trad_ISRF + urad_CMB_CGS_eV * Trad_CMB) / (urad_ISRF_CGS_eV + urad_CMB_CGS_eV); // weighting by SED energy
}

#endif



#ifdef RADTRANSFER
/***********************************************************************************************************/
/* this function initializes some of the variables we need */
/***********************************************************************************************************/
void rt_set_simple_inits(int RestartFlag)
{
    if(RestartFlag==1) return;
    int flag_to_reset_values_on_startup = 0;
    if(RestartFlag==0) {flag_to_reset_values_on_startup = 1;}
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL) && defined(SINGLE_STAR_RESTART_FROM_FIRESIM)
    if(RestartFlag==2) {flag_to_reset_values_on_startup = 1;}
#endif
    
    int i; for(i = 0; i < NumPart; i++)
    {
        if(P[i].Type == 0)
        {
            int k;
#ifdef RT_INFRARED
            if(flag_to_reset_values_on_startup) {SphP[i].Radiation_Temperature = SphP[i].Dust_Temperature = DMIN(All.InitGasTemp,100.);} //get_min_allowed_dustIRrad_temperature(); // in K, floor = CMB temperature or 10K
#ifdef RT_ISRF_BACKGROUND
            if(flag_to_reset_values_on_startup) {SphP[i].Radiation_Temperature = background_isrf_cmb_Teff();} //SphP[i].Dust_Temperature;
#endif
            SphP[i].Dt_Rad_E_gamma_T_weighted_IR = 0;
#endif
#ifdef RT_CHEM_PHOTOION
            if(flag_to_reset_values_on_startup)
            {
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
            } else {
#ifdef RT_CHEM_PHOTOION_HE
                SphP[i].HeIII = DMIN(1.0, DMAX(MIN_REAL_NUMBER, 1.0 - SphP[i].HeII - SphP[i].HeI * (4.0 * HYDROGEN_MASSFRAC)/(1-HYDROGEN_MASSFRAC))); // not read in since for this subroutine follows from others
#endif
            }
#endif
#ifdef RT_ISRF_BACKGROUND
            double urad[N_RT_FREQ_BINS];
            get_background_isrf_urad(i, urad);
#endif
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if(flag_to_reset_values_on_startup) {SphP[i].Rad_E_gamma[k] = MIN_REAL_NUMBER;}
                int flag_to_reset_values_on_startup_et = flag_to_reset_values_on_startup;
#if !defined(OUTPUT_EDDINGTON_TENSOR)
                flag_to_reset_values_on_startup_et = 1;
#endif
                if(flag_to_reset_values_on_startup_et) {SphP[i].ET[k][0]=SphP[i].ET[k][1]=SphP[i].ET[k][2]=1./3.; SphP[i].ET[k][3]=SphP[i].ET[k][4]=SphP[i].ET[k][5]=0;}
                SphP[i].Rad_Je[k] = 0;
#ifdef RT_FLUXLIMITER
                SphP[i].Rad_Flux_Limiter[k] = 1;
#endif

#ifdef RT_INFRARED
                if(flag_to_reset_values_on_startup && k==RT_FREQ_BIN_INFRARED) { // only initialize the IR energy if starting a new run, otherwise use what's in the snapshot
                    SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED] = (4.*5.67e-5 / C_LIGHT_CGS) * pow(DMIN(All.InitGasTemp,100.),4.) / UNIT_PRESSURE_IN_CGS * P[i].Mass / (SphP[i].Density*All.cf_a3inv);
                }
#endif
#ifdef RT_ISRF_BACKGROUND
                if(flag_to_reset_values_on_startup) {SphP[i].Rad_E_gamma[k] = urad[k] * P[i].Mass / (SphP[i].Density*All.cf_a3inv);}
#endif
#ifdef RT_EVOLVE_ENERGY
                SphP[i].Rad_E_gamma_Pred[k] = SphP[i].Rad_E_gamma[k]; 
                SphP[i].Dt_Rad_E_gamma[k] = 0;
#endif
#ifdef RT_EVOLVE_FLUX
                int k_dir, flag_to_reset_values_on_startup_flux = flag_to_reset_values_on_startup;
#if !defined(OUTPUT_RT_RAD_FLUX)
                flag_to_reset_values_on_startup_flux = 1;
#endif
                if(flag_to_reset_values_on_startup_flux) {
                    for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Rad_Flux_Pred[k][k_dir] = 0;}
                } else {
                    for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Rad_Flux_Pred[k][k_dir] *= P[i].Mass/(SphP[i].Density*All.cf_a3inv);} // need to correct the units here before using
                }
                for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Rad_Flux[k][k_dir] = SphP[i].Rad_Flux_Pred[k][k_dir];}
                for(k_dir=0;k_dir<3;k_dir++) {SphP[i].Dt_Rad_Flux[k][k_dir] = 0;}
#endif
#ifdef RT_EVOLVE_INTENSITIES
                int k_dir; for(k_dir=0;k_dir<N_RT_INTENSITY_BINS;k_dir++) {SphP[i].Rad_Intensity_Pred[k][k_dir] = SphP[i].Rad_Intensity[k][k_dir] = MIN_REAL_NUMBER; SphP[i].Dt_Rad_Intensity[k][k_dir] = 0;}
#endif
                
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                double q_a = (0.75*All.Grain_Q_at_MaxGrainSize) / ((All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS)*(All.Grain_Size_Max/UNIT_LENGTH_IN_CGS)), e0 = All.Vertical_Grain_Accel / q_a, kappa_0 = All.Grain_Absorbed_Fraction_vs_Total_Extinction * q_a * All.Dust_to_Gas_Mass_Ratio, cell_vol = (P[i].Mass/SphP[i].Density);
                double rho_base_setup = 1., H_scale_setup = 1.*boxSize_X; // define in code units the -assumed- initial scaling of the base gas density and vertical scale-length (PROBLEM SPECIFIC HERE!)
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
                kappa_0 *= sqrt(All.Grain_Size_Max / All.Grain_Size_Min); // opacity must be corrected for dependence of Q on grainsize or lack thereof
#endif
                double E_cell_thin = cell_vol * e0 * exp(-kappa_0*rho_base_setup*H_scale_setup*(1.-exp(-P[i].Pos[2]/H_scale_setup))), E_cell=E_cell_thin; // attenuate according to equilibrium expectation, if we're using single-scattering radiation pressure [otherwise comment this line out] //
                double tau_tot = q_a * All.Dust_to_Gas_Mass_Ratio * rho_base_setup*H_scale_setup; if(tau_tot>1) {E_cell = cell_vol * (3.*All.Vertical_Grain_Accel*All.Dust_to_Gas_Mass_Ratio*rho_base_setup*H_scale_setup) * (exp(-P[i].Pos[2]/H_scale_setup) + 1./tau_tot);} // attenuate according to approximate optically-thick expression with free-streaming from the 'photosphere' when optically thin
                SphP[i].Rad_E_gamma_Pred[k] = SphP[i].Rad_E_gamma[k] = E_cell;
#if defined(RT_EVOLVE_FLUX)
                SphP[i].Rad_Flux_Pred[k][2]=SphP[i].Rad_Flux[k][2] = E_cell_thin*C_LIGHT_CODE_REDUCED;
                SphP[i].Rad_Flux[k][0]=SphP[i].Rad_Flux[k][1]=SphP[i].Rad_Flux_Pred[k][0]=SphP[i].Rad_Flux_Pred[k][1]=0;
#endif
#endif
                SphP[i].Rad_Kappa[k] = rt_kappa(i,k); // let everything else get reset first before calling this //
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
    if(n_polar < 1) {printf("Number of rays is invalid (<1). Terminating.\n"); fflush(stdout); endrun(5346343);}

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
/* subroutine specific to some bands where, rather than modeling absorption and emission explicitly (with some emissivity
    associated with e.g. gas or dust), we assume an instantaneous exact local radiative equilibrium re-emission from
    band A into band B, in the absorption step essentially transferring photon energy instantly between bins */
/***********************************************************************************************************/
int rt_get_donation_target_bin(int bin)
{
    int donation_target_bin = -1; // default here is to assume no 'target bin' -- meaning the absorbed radiation does not get 'transferred' but simply absorbed (emission handled separately)
#if defined(RT_CHEM_PHOTOION) && defined(RT_OPTICAL_NIR)
    /* in these modules, as typically applied to galaxy and star-formation simulations, we will
       assume absorbed ionizing photons are instantly re-emitted via recombination into optical-NIR bins. valid if recombination time is
       fast compared to all our timesteps. more accurately, this should be separately calculated in the cooling rates, and gas treated as a source, but for now we use this module */
    if(bin==RT_FREQ_BIN_H0) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
#ifdef RT_PHOTOION_MULTIFREQUENCY
    if(bin==RT_FREQ_BIN_He0) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
    if(bin==RT_FREQ_BIN_He1) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
    if(bin==RT_FREQ_BIN_He2) {donation_target_bin=RT_FREQ_BIN_OPTICAL_NIR;}
#endif
#endif
#if defined(RT_PHOTOELECTRIC) && defined(RT_INFRARED)
    if(bin==RT_FREQ_BIN_PHOTOELECTRIC) {donation_target_bin=RT_FREQ_BIN_INFRARED;} /* this is direct dust absorption, re-radiated in IR if we aren't explicitly modeling dust thermal physics (module here) */
#endif
#if defined(RT_NUV) && defined(RT_INFRARED)
    if(bin==RT_FREQ_BIN_NUV) {donation_target_bin=RT_FREQ_BIN_INFRARED;} /* this is direct dust absorption, re-radiated in IR if we aren't explicitly modeling dust thermal physics (module here) */
#endif
#if defined(RT_OPTICAL_NIR) && defined(RT_INFRARED)
    if(bin==RT_FREQ_BIN_OPTICAL_NIR) {donation_target_bin=RT_FREQ_BIN_INFRARED;} /* this is direct dust absorption, re-radiated in IR if we aren't explicitly modeling dust thermal physics (module here) */
#endif
    return donation_target_bin;
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



/* simple code to make special exceptions for certain bands to allow radiation pressure forces to exceed absorbed photon momentum as needed */
int check_if_absorbed_photons_can_be_reemitted_into_same_band(int kfreq)
{
    int checker = 0; // default to no, but this isn't always true
#if defined(GALSF_FB_FIRE_RT_LONGRANGE)
    if(kfreq==RT_FREQ_BIN_FIRE_IR) {checker=1;} // skip
#endif
#if defined(RT_FREEFREE)
    if(kfreq==RT_FREQ_BIN_FREEFREE) {checker=1;} // skip
#endif
#if defined(RT_GENERIC_USER_FREQ)
    if(kfreq==RT_FREQ_BIN_GENERIC_USER_FREQ) {checker=1;} // skip
#endif
#if defined(RT_INFRARED)
    if(kfreq==RT_FREQ_BIN_INFRARED) {checker=1;} // skip
#endif
    return checker;
}



/*****************************************************************************
Routines specifically for handling thermal and radiative coupling between gas,
dust, and the IR radiation field component (RT_INFRARED)
*****************************************************************************/

#ifdef RT_INFRARED

/* return the minimum user-specified dust temperature. note there is nothing physical about this, just a convenience function since we enforce a minimum -gas- temperature */
double get_min_allowed_dustIRrad_temperature(void)
{
#if defined(GALSF)
    return DMAX(All.MinGasTemp, 2.73/All.cf_atime);
#endif
    return MIN_REAL_NUMBER;
}

/* dust_dE_cooling
Returns the derivative of the dust energy, to be root-solved to 0 to obtain
the equilibrium dust temperature assuming 0 dust heat capacity. Accounts for
dust emission/absorption to obtain a solution consistent with a backward-Euler
implicit cooling solution.

Parameters
----------
i - index of particle in particle list
Tgas - Assumed final gas temperature at the end of the timestep
Tdust - Dust temperature (generally the unknown quantity we need to solve for)
Tdust_fixedpoint_1 - Stores one possible fixed-point iterate for approximating the root (can use as a guess)
Tdust_fixedpoint_2 - Stores another possible fixed-point iterate for approximating the root 

Returns
-------
dE - net dust heating (=0 for dust in equilibrium)
*/
double dust_dE_cooling(int i, double Tgas, double Tdust, double* Tdust_fixedpoint_1, double* Tdust_fixedpoint_2){
    double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
    double nHcgs = HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * SphP[i].Density * All.cf_a3inv / PROTONMASS_CGS;
    double lambda_to_dErad = (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE) * nHcgs * nHcgs * (dt*UNIT_TIME_IN_CGS) / (SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS) / (UNIT_SPECEGY_IN_CGS) * P[i].Mass; /* need to account for RSOL factors in emission/absorption rates */
    
    double dust_absorption_nonIR = 0;
    for(int k=0; k < N_RT_FREQ_BINS; k++){
        if(RT_BAND_IS_IONIZING(k)){continue;} /* gas-phase absorption */
        if(k==RT_FREQ_BIN_INFRARED){continue;} /* this is only counting up non-IR contributions, e.g. nebular NUV */
        double e_final = SphP[i].Rad_E_gamma[k] + SphP[i].Lambda_RadiativeCooling_toRHDBins[k] * lambda_to_dErad;
        e_final = DMAX(0,e_final); // check against overshoot into negative values
        double absrate_k = rt_absorption_rate(i, k) * dt; // this needs to be positive to sensible behavior here
        if(absrate_k > 0) {dust_absorption_nonIR += e_final * fabs(expm1(-absrate_k));}
    }
    double alpha_gd = gas_dust_heating_coeff(i,Tgas,Tdust);
    double LambdaDust = alpha_gd * (Tgas-Tdust);
    double de_IR_dust = LambdaDust * lambda_to_dErad; // equates to *net* emission of radiation by dust (emission - absorption)
    double LambdaIR_gas = SphP[i].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_INFRARED];
    double de_IR_gas = LambdaIR_gas * lambda_to_dErad; // net emission by gas
    
    double kappa_dust_emission = rt_kappa_adaptive_IR_band(i, Tdust, Tdust, 1,1);
    double fac_emission = 4.*5.67e-5/(UNIT_PRESSURE_IN_CGS*UNIT_VEL_IN_CGS)*P[i].Mass*RT_SPEEDOFLIGHT_REDUCTION*dt;
    double dust_emission = fac_emission*kappa_dust_emission*pow(Tdust,4); // *total* dust emission
    
    double T_IR_0 = SphP[i].Radiation_Temperature;
    double Tmax = DMAX(DMAX(Tgas, Tdust),T_IR_0), Tmin = DMIN(DMIN(Tgas,Tdust), T_IR_0);
    double e_IR_0 = SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED];
    double e_IR_final = DMAX(e_IR_0 + de_IR_dust + de_IR_gas, 0);
    /* below we have some approximate scalings for the temperature modification, since the different photon-gas interchange operations can behave differently in terms of photon number changes (N_emitted - N_absorbed); we interpolate between regimes to prevent unphysical behavior in the strong-coupling limit */
    //double T_IR_final = (e_IR_final/e_IR_0) * T_IR_0; /* begin by guessing the value we would have if these operations conserved photon number */
    //if(T_IR_final > Tmax || T_IR_final < Tmin) { /* this gives an unphysical evolution in radiation temperature */
    //if(de_IR_gas + de_IR_dust < 0) {T_IR_final = T_IR_0;} else {  /* only absorption occurs, radiation temperature unmodified */
    double T_IR_final = e_IR_final / (DMAX(e_IR_0/T_IR_0,0) + DMAX(de_IR_gas / Tgas,0) + DMAX(de_IR_dust / Tdust,0) + MIN_REAL_NUMBER); // assume gas-phase IR component emitted at Tgas. /* if new emission occurs, radiation temperature updated assuming all energy gained comes from emission with radiation temperature = gas temperature; if absorption occurs, radiation temperature unmodified  */
    T_IR_final = DMAX(Tmin, DMIN(T_IR_final, Tmax)); // limiters for edge-cases above //
    //}}
    SphP[i].Radiation_Temperature_CoolingWeighted = T_IR_final; // PFH: this had been protected by a Tdust < max, but that will give the wrong radiation temperature when radiation-gas (non dust) interactions like Compton or Kramers are important, needs to be set regardless of whether the dust itself is present

    double dE_dust = 0; // now count up the energy changes in the dust for us to solve for 0
    double dust_absorption = dust_absorption_nonIR;
    dust_absorption += e_IR_final * C_LIGHT_CODE_REDUCED * rt_kappa_adaptive_IR_band(i, Tdust, T_IR_final,-1,1) * SphP[i].Density*All.cf_a3inv * dt;
    double result = LambdaDust * lambda_to_dErad + dust_absorption - dust_emission;

    double Tdust_fixed1_tmp = Tgas + (dust_absorption - dust_emission)/(alpha_gd*lambda_to_dErad + MIN_REAL_NUMBER); // make sure to include term in denominator to protect vs nans
    double Tdust_fixed2_tmp = sqrt(sqrt(DMAX(0,LambdaDust * lambda_to_dErad + dust_absorption)/(fac_emission * kappa_dust_emission + MIN_REAL_NUMBER))); // make sure to include term in denominator and MAX in numerator to protect vs nans
    // PFH: better behavior for 2 above is to bracket with solution for fixed dust absorption, which gives T2 = (T0 == above with Tdust = 0 in LambdaDust, positive-definite) x , where x solves x^4+b*x=1, b=((alpha_gd * lambda_to_dErad)/(fac_emission * kappa_dust_emission * T0^3)), should be dimensionless, for b << 1, x=1, for b >> 1, x->1/b; the above essentially jumps to b->infinity in the latter regime
    
    *Tdust_fixedpoint_1 = DMAX(Tdust_fixed1_tmp, 0); // check against overshoot into negative values
    *Tdust_fixedpoint_2 = DMAX(Tdust_fixed2_tmp, 0); // check against overshoot into negative values
    return result;
}

/* Returns the dust cooling rate LambdaDust in erg s^-1 cm^3, accounting for 
dust emission
*/
double rt_ir_lambdadust(int i, double T){
    double Tdust, T_lower, T_upper, dE, dE1, dE2, dE_lower, dE_upper, dE_guess, dTdust_tol=1e-6;
    double Tdust_fixedpoint_1, Tdust_fixedpoint_2, dummy;
    // define ROOTFIND_FUNCTION_INNER because this gets called nested inside the cooling solver, and needs to be def'd distinctly from the overlying ROOTFIND_FUNCTION
    #define ROOTFIND_FUNCTION_INNER(dTdust) dust_dE_cooling(i, T, T+dTdust, &Tdust_fixedpoint_1, &Tdust_fixedpoint_2)
    if((All.Time==0 )|| (!isfinite(SphP[i].Dust_Temperature))){Tdust=T;} else {Tdust = DMIN(1e3,SphP[i].Dust_Temperature);}

    dE = dE_guess = ROOTFIND_FUNCTION_INNER(Tdust-T);
   
    if(Tdust_fixedpoint_1 > 0 && Tdust_fixedpoint_1 < MAX_DUST_TEMP) {dE1 =  dust_dE_cooling(i, T, Tdust_fixedpoint_1, &dummy, &dummy);} else {dE1 = MAX_REAL_NUMBER;}
    if(Tdust_fixedpoint_2 > 0 && Tdust_fixedpoint_2 < MAX_DUST_TEMP) {dE2 =  dust_dE_cooling(i, T, Tdust_fixedpoint_2, &dummy, &dummy);} else {dE2 = MAX_REAL_NUMBER;}

    // error estimate of the fixed-point guesses, used for bracketing below
    double fixedpoint_error = DMIN(fabs(Tdust-Tdust_fixedpoint_2), fabs(Tdust-Tdust_fixedpoint_1))/Tdust; 
    if(fabs(dE1) < fabs(dE)){Tdust = Tdust_fixedpoint_1; dE=dE_guess=dE1;}
    if(fabs(dE2) < fabs(dE)){Tdust = Tdust_fixedpoint_2; dE=dE_guess=dE2;}
    
    /* bracketing the dust temperature */
    int n_iter = 0;
    if(dE < 0)
    {
        double scalefac = DMAX(0.9, 1-fixedpoint_error);
        T_upper = Tdust;
        dE_upper = dE_guess; 
        while(dE < 0) {
            Tdust *= scalefac; 
            dE = ROOTFIND_FUNCTION_INNER(Tdust-T);
            if(dE==0){break;}
            scalefac *= 0.9; 
            n_iter++;
        }
        T_lower = Tdust, dE_lower = dE;
    } else {
        T_lower = Tdust, dE_lower = dE_guess;
        double scalefac = DMIN(1.1, 1+fixedpoint_error);
        while(dE > 0 && Tdust < MAX_DUST_TEMP) {
            Tdust *= scalefac; Tdust = DMIN(Tdust,MAX_DUST_TEMP);
            dE = ROOTFIND_FUNCTION_INNER(Tdust-T); 
            if(dE==0){break;}
            scalefac *= 1.1; 
            n_iter++;
        }
        T_upper = Tdust, dE_upper = dE;
    }     
    if(T_upper>=MAX_DUST_TEMP && dE_upper > 0){SphP[i].Dust_Temperature = MAX_DUST_TEMP; return 0;}

    if(dE_lower * dE_upper > 0){ PRINT_WARNING("Failed to bracket Tdust solution for ID=%lld T=%g T_lower=%g T_upper=%g dE_lower=%g dE_upper=%g\n", P[i].ID, T, T_lower,T_upper, dE_lower, dE_upper);}

    if(dE!=0){  // root-solve for Tdust
        double ROOTFIND_X_a = T_lower-T, ROOTFIND_X_b = T_upper-T;
        double ROOTFUNC_a = dE_lower; double ROOTFUNC_b = dE_upper;
        double ROOTFIND_REL_X_tol = dTdust_tol, ROOTFIND_ABS_X_tol=0.;
        #include "../system/bracketed_rootfind.h"
        Tdust = ROOTFIND_X_new+T;
    }
    double LambdaDust = gas_dust_heating_coeff(i,T,Tdust) * (T-Tdust);
    SphP[i].Lambda_RadiativeCooling_toRHDBins[RT_FREQ_BIN_INFRARED] += LambdaDust;
    SphP[i].Dust_Temperature = Tdust;
    return LambdaDust;
}

#endif

/******************************************************************************************************
This returns the volumetric quantity de/dt = heat transfer from gas + photon absorption - emission 
for dust, which we will root-find to determine the dust temperature. 

dust_absorption_rate must be passed as the dust photon absorption rate per unit volume in code units,
correcting for the reduced speed of light if applicable.
******************************************************************************************************/
double dust_dEdt(int i, double T, double Tdust, double dust_absorption_rate)
{
    double nHcgs = HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * SphP[i].Density * All.cf_a3inv / PROTONMASS_CGS;    /* hydrogen number dens in cgs units */
    double fac_emission = 4.*5.67e-5/(UNIT_PRESSURE_IN_CGS*UNIT_VEL_IN_CGS)*SphP[i].Density*All.cf_a3inv; // in code units
    double LambdaDust_fac = 0;
#ifdef COOLING
    if(T>0) {LambdaDust_fac = gas_dust_heating_coeff(i,T,Tdust) * nHcgs * nHcgs /(UNIT_PRESSURE_IN_CGS/UNIT_TIME_IN_CGS);}
#endif    
    double kappa_emission = rt_kappa_adaptive_IR_band(i, Tdust, Tdust, 1, 1);
    double dust_emission = fac_emission * kappa_emission * pow(Tdust,4);
#if defined(COOLING) && !defined(RT_INFRARED) // if we aren't doing RT self-consistently, approximate outward radiative transport rate in optically-thick regime
    double column = evaluate_NH_from_GradRho(SphP[i].Gradients.Density,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i);
    double tau = column * kappa_emission;
    dust_emission /= (1 + tau*tau); // e.g. Masunaha & Inutsuka 1999, Rafikov 2007
#endif
    return LambdaDust_fac * (T-Tdust) + dust_absorption_rate - dust_emission;
}

/***********************************************************************************************************
Returns the equilibrium dust temperature as a function of gas temperature and dust absorption rate.
dust_absorption_rate must be passed as the photon energy absorption rate per unit volume in code units,
correcting for the reduced speed of light if applicable. If T=0 then gas-dust coupling is neglected
and we only solve for equilibrium between emission and absorption.
************************************************************************************************************/
double rt_eqm_dust_temp(int i, double T, double dust_absorption_rate)
{
    double T_old, T_lower=0, T_upper=MAX_REAL_NUMBER, T_secant, Tdust_guess, Tdust, dEdt, dEdt_upper, dEdt_lower, fac, dEdt_guess, scalefac;
    double Tmax=1e10; // upper-bound dust temperature above which we definitely don't believe our detailed (tiny) dust abundance
    Tmax = MAX_DUST_TEMP; // this is now a global variable
    /* First we come up with a reasonable guess for the dust temp based on available info */
#ifdef RT_INFRARED
    Tdust_guess = DMIN(DMAX(SphP[i].Dust_Temperature,1.),1e3); // previous dust temperature should be a good guess
#else // case where we don't have a pre-computed dust temp, use asymptotic limits to get a good guess
    double Zfac = 1.0;
#ifdef METALS
    if(i>=0) {Zfac = P[i].Metallicity[0]/All.SolarAbundances[0];}
#endif
    double rho_c_arad_fac = (4.*5.67e-5)/(UNIT_VEL_IN_CGS*UNIT_PRESSURE_IN_CGS)*SphP[i].Density*All.cf_a3inv; // a c rho in code units
    Tdust_guess = sqrt(cbrt(100 * dust_absorption_rate/(rho_c_arad_fac * (0.1*UNIT_SURFDEN_IN_CGS) * Zfac)));  // guess neglecting gas-dust coupling term and assuming a beta=2 emission opacity law kappa = 0.1 cm^2/g Z (T/10K)^2
    Tdust_guess = DMAX(Tdust_guess, sqrt(sqrt(dust_absorption_rate / (rho_c_arad_fac * (5.*UNIT_SURFDEN_IN_CGS) * Zfac)))); // account for how opacity tops out around 5 Z cm^2/g
#ifdef COOLING // account for gas-dust coupling
    double nHcgs = HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * SphP[i].Density * All.cf_a3inv / PROTONMASS_CGS;    /* hydrogen number dens in cgs units */
    double LambdaDust_fac = gas_dust_heating_coeff(i,T,Tdust_guess) * nHcgs * nHcgs /(UNIT_PRESSURE_IN_CGS/UNIT_TIME_IN_CGS);
    double Tdust_coupled = T - rho_c_arad_fac * rt_kappa_adaptive_IR_band(i,T,T,1,0) * pow(T,4) / (LambdaDust_fac+MIN_REAL_NUMBER); // bound for the gas-dust coupled regime assuming T ~ Td
    Tdust_guess = DMAX(Tdust_coupled, Tdust_guess);
#endif
#endif // end non-RT case for guess
    /* We now have our initial guess */    
    if(T==0) {return Tdust_guess;} // if just calling for a rough estimate this is good enough

    Tdust = Tdust_guess;
    int n_iter=0;    
    dEdt_guess = dEdt = dust_dEdt(i,T,Tdust_guess,dust_absorption_rate);
    
    if(dEdt==0){return Tdust_guess;}
    /* bracketing the dust temperature */
    if(dEdt < 0)
    {
	scalefac = 0.9;
	T_upper = DMIN(Tmax,Tdust), dEdt_upper = dEdt_guess; 
	while(dEdt<0) {
	    Tdust *= scalefac; 
	    dEdt = dust_dEdt(i,T,Tdust,dust_absorption_rate); 
        if(dEdt==0){return Tdust;}
	    scalefac *= 0.9; 
	    n_iter++;
	}
	T_lower = Tdust, dEdt_lower = dEdt;
    } else {
	T_lower = Tdust, dEdt_lower = dEdt_guess;
	scalefac = 1.1;
	while(dEdt>0 && Tdust < Tmax) {
	    Tdust *= scalefac; Tdust = DMIN(Tdust,Tmax);
	    dEdt = dust_dEdt(i,T,Tdust,dust_absorption_rate); 
        if(dEdt==0){return Tdust;}
	    scalefac *= 1.1; 
	    n_iter++;
	    }
	    T_upper = Tdust, dEdt_upper = dEdt;
    }     
    if(T_upper==Tmax && dEdt_upper > 0) {return Tmax;}
    if(T_lower>=Tmax) {return Tmax;}

#if 1  // PFH: still testing which option is better, but the new rootfind struggles here, in hyper-zoom-in runs when given dust close to max temperature (raising max temp resolves the failure to converge or Nan's but then jumps to very high solutions somewhat randomly, where it shouldnt. The old secant routine below appears stable and more robust in this particular instance for now.
    #define ROOTFIND_FUNCTION(dTdust) dust_dEdt(i,T,T+dTdust,dust_absorption_rate); // here we want to converge on a relative tolerance for Tdust-Tgas
    double ROOTFIND_X_a = T_upper-T, ROOTFIND_X_b = T_lower-T, ROOTFUNC_a = dEdt_upper, ROOTFUNC_b = dEdt_lower, ROOTFIND_REL_X_tol = 1e-6, ROOTFIND_ABS_X_tol=0.;
    #include "../system/bracketed_rootfind.h"
    Tdust = ROOTFIND_X_new + T;
    if(ROOTFIND_ITER > MAXITER || isnan(Tdust)){PRINT_WARNING("WARNING: Particle %lld did not converge to desired Tdust tolerance (iter=%d, Tdust=%g, Tgas=%g)\n",(long long)P[i].ID,ROOTFIND_ITER,Tdust,T);}
    
#else

    T_old = Tdust; double dEdt_old = dEdt; Tdust = Tdust_guess; dEdt = dEdt_guess; // For our second guess we take the backeting value opposite of the initial guess.
    double dT_dustgas = T-Tdust;
    do  // secant method iterations with bisection as a backup; usually converges to machine epsilon in 4-5 iterations
    {
        dT_dustgas = T - Tdust;
        T_secant = Tdust - dEdt * (Tdust - T_old) / (dEdt - dEdt_old);
        T_secant = DMAX(DMIN(T_secant,T_upper),T_lower);
        dEdt_old = dEdt;
        dEdt = dust_dEdt(i,T,T_secant,dust_absorption_rate);
        fac = fabs(T_secant - Tdust)/(MIN_REAL_NUMBER+fabs(Tdust-T_old)); //fabs(dEdt)/(MIN_REAL_NUMBER+fabs(dEdt_old));
        if(fac < 0.5) { // accept the secant iteration if it is converging more rapidly
            T_old=Tdust;
            Tdust=T_secant;
        } else { // if secant isn't working do bisection iteration instead; guaranteed to reduce the error
            T_old = Tdust;
            Tdust = sqrt(T_lower*T_upper);
            dEdt = dust_dEdt(i,T,Tdust,dust_absorption_rate);
            fac = 0.5;
        }
        if(dEdt>0) {T_lower=Tdust;} else {T_upper=Tdust;} // either way, update upper and lower bounds
        n_iter++;
        if(n_iter > MAXITER-10) {
            PRINT_WARNING("Warning: Dust temperature iteration converging slowly: ID=%lld iter=%d T=%g Tdust=%g Tdust_guess=%g T_upper=%g T_lower=%g dEdt=%g fac=%g.\n",(long long)P[i].ID,n_iter,T,Tdust,Tdust_guess, T_upper, T_lower,dEdt, fac);
            if(n_iter > MAXITER){break;}
        }
    } while(fabs(dT_dustgas - (T-Tdust)) > 1.e-3 * fabs(T-Tdust)); // sufficient to converge dust cooling to 10^-3 tolerance, at this point uncertainties in dust properties will dominate the error budget

#endif
    
    return Tdust;
}

/***********************************************************************************************************/
/* returns the fraction of a blackbody SED in a given photon energy band - accurate to <1% over all wavelengths
   E_lower - lower end of the energy band in eV
   E_upper - upper end of the energy band in eV
   T_eff - effective blackbody temperature of the SED, in K
*/
/***********************************************************************************************************/
double blackbody_lum_frac(double E_lower, double E_upper, double T_eff)
{
    double k_B = 8.617e-5; // Boltzmann constant in eV/K
    double x1 = E_lower / (k_B * T_eff), x2 = E_upper / (k_B * T_eff), f_lower, f_upper;
    if(x1 < 3.40309){
      f_lower = (131.4045728599595*x1*x1*x1)/(2560. + x1*(960. + x1*(232. + 39.*x1))); // approximation of integral of Planck function from 0 to x1, valid for x1 << 1
    } else {
      f_lower = 1 - (0.15398973382026504*(6. + x1*(6. + x1*(3. + x1))))*exp(-DMIN(x1,40.)); // approximation of Planck integral for large x
    }
    if(x2 < 3.40309){
      f_upper = (131.4045728599595*x2*x2*x2)/(2560. + x2*(960. + x2*(232. + 39.*x2))); // approximation of integral of Planck function from 0 to x2, valid for x2 << 1
    } else {
      f_upper = 1 - (0.15398973382026504*(6. + x2*(6. + x2*(3. + x2))))*exp(-DMIN(x2,40.)); // approximation of Planck integral for large x
    }
    double df = f_upper - f_lower;
    if(df<=0) {if(x2<=x1) {return 0;} else {if(x1>4.) {if(x1<120.) {return 0.15398973382026504*(6.+x1*(6.+x1*(3.+x1)))*exp(-DMIN(x1,120.));} else {return 2.e-47;}}}}
    return DMAX(df, 0);
}

/* subroutine to return the photon energy density [in physical code units] in a given band range [i - index of star particle, E_lower - lower end of the energy band in eV, E_upper - upper end of the energy band in eV] */
double rt_irband_egydensity_in_band(int i, double E_lower, double E_upper)
{
#if defined(RT_INFRARED)
    double u_gamma = SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED] * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * blackbody_lum_frac(E_lower, E_upper, SphP[i].Radiation_Temperature);
    if(!isfinite(u_gamma) || (u_gamma<0)) {u_gamma = 0;}
    return u_gamma;
#else
    return 0;
#endif
}

/***********************************************************************************************************/
/* returns the fraction of a star's SED (approximated as a blackbody) in a given photon energy band - accurate to <1% over all wavelengths
   i - index of star particle
   E_lower - lower end of the energy band in eV
   E_upper - upper end of the energy band in eV
*/
/***********************************************************************************************************/
double stellar_lum_in_band(int i, double E_lower, double E_upper)
{
#if   defined(SINGLE_STAR_SINK_DYNAMICS) // use generic fits based on mass
    double l_sol=bh_lum_bol(0,P[i].Mass,i)*UNIT_LUM_IN_SOLAR, m_sol=P[i].Mass*UNIT_MASS_IN_SOLAR, r_sol=pow(m_sol,0.738); // L/Lsun, M/Msun, R/Rsun
#else
    double l_sol=1., r_sol=1.; // nothing usefully defined for the above - default to solar-type stars //
#endif
    double T_eff = 5780. * pow(l_sol/(r_sol*r_sol), 0.25);
    double f = blackbody_lum_frac(E_lower, E_upper, T_eff);
    return f * l_sol / UNIT_LUM_IN_SOLAR;
}




#if defined(CHIMES_STELLAR_FLUXES) && (defined(RADTRANSFER) || defined(RT_USE_GRAVTREE))
/* The following routines are fitting functions that are used to obtain the luminosities in the 6-13.6 eV energy band (i.e. G0)
 * and the >13.6 eV band (i.e. H-ionising), which will be used by CHIMES. These functions were fit to Starburst99 models
 * that used the Geneva 2012/13 tracks with v=0.4 rotation and Z=0.014 metallicity. These are separated because CHIMES uses its special age-bins isntead of freq-bins */

double chimes_G0_luminosity(double stellar_age, double stellar_mass) // age in Myr, mass in Msol, return value in Habing units * cm^2
{
  double zeta = 6.5006802e29;
  if (stellar_age < 4.07) {return stellar_mass * exp(89.67 + (0.172 * pow(stellar_age, 0.916)));}
    else {return stellar_mass * zeta * pow(1773082.52 / stellar_age, 1.667) * pow(1.0 + pow(stellar_age / 1773082.52, 28.164), 1.64824);}
}

double chimes_ion_luminosity(double stellar_age, double stellar_mass) // age in Myr, mass in Msol, return value in s^-1
{
  double zeta = 3.2758118e21;
  if (stellar_age < 3.71) {return stellar_mass * exp(107.21 + (0.111 * pow(stellar_age, 0.974)));}
    else {return stellar_mass * zeta * pow(688952.27 / stellar_age, 4.788) * pow(1.0 + pow(stellar_age / 688952.27, 1.124), -17017.50356);}
}

int rt_get_source_luminosity_chimes(int i, int mode, double *lum, double *chimes_lum_G0, double *chimes_lum_ion)
{
    int value_to_return = 0;
    value_to_return = rt_get_source_luminosity(i, mode, lum); // call routine as normal for all bands, before adding chimes-specific details
    if( ((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3)))) && (P[i].Mass>0) && (PPP[i].Hsml>0) )
    {
        int age_bin, j; double age_Myr=1000.*evaluate_stellar_age_Gyr(i), log_age_Myr=log10(age_Myr), stellar_mass=P[i].Mass*UNIT_MASS_IN_SOLAR;
        if(log_age_Myr < CHIMES_LOCAL_UV_AGE_LOW) {age_bin = 0;} else if (log_age_Myr < CHIMES_LOCAL_UV_AGE_MID) {age_bin = (int) floor(((log_age_Myr - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW) + 1);} else {
            age_bin = (int) floor((((log_age_Myr - CHIMES_LOCAL_UV_AGE_MID) / CHIMES_LOCAL_UV_DELTA_AGE_HI) + ((CHIMES_LOCAL_UV_AGE_MID - CHIMES_LOCAL_UV_AGE_LOW) / CHIMES_LOCAL_UV_DELTA_AGE_LOW)) + 1);
            if (age_bin > CHIMES_LOCAL_UV_NBINS - 1) {age_bin = CHIMES_LOCAL_UV_NBINS - 1;}}
    
        for(j=0;j<CHIMES_LOCAL_UV_NBINS;j++) {chimes_lum_G0[j]=0; chimes_lum_ion[j]=0;}
        chimes_lum_G0[age_bin] = chimes_G0_luminosity(age_Myr,stellar_mass) * All.Chimes_f_esc_G0;
        chimes_lum_ion[age_bin] = chimes_ion_luminosity(age_Myr,stellar_mass) * All.Chimes_f_esc_ion;
    }
    return value_to_return;
}
#endif 



/*--------------------------------------------------------------------
  calculate the IR dust opacity [in physical code units = Length^2/Mass].
   NOTE: The flag do_emission_absorption_scattering_opacity toggles special behaviour.
   -1: returns the absorption opacity only, using the radiation temperature
    0: returns the scattering+absorption opacity using the radiation temperature (usually want this)
    1: returns the *emission* opacity, assuming the dust+gas radiates as a blackbody (depends only on T_dust)
   likewise, dust_or_gas_opacity_only_flag toggles different behaviors:
    0: total IR-band opacity,
    1: opacity -only- from dust,
   -1: opacity -only- from non-dust 
--------------------------------------------------------------------*/
double rt_kappa_adaptive_IR_band(int i, double T_dust, double Trad, int do_emission_absorption_scattering_opacity, int dust_or_gas_opacity_only_flag)
{
    if(do_emission_absorption_scattering_opacity==1) {Trad = T_dust;} // if we want the emissivity then we assume radiation emitted at T_dust
    double fac=UNIT_SURFDEN_IN_CGS, x = 4.*log10(Trad) - 8., kappa=0, T_dust_opacitytable = T_dust; // needed for fitting functions to opacities (may come up with cheaper function later)
    double dx_excess=0; if(x > 7.) {dx_excess=x-7.; x=7.;} // cap for maximum temperatures at which fit-functions should be used //
    //if(x < -4.) {x=-4.;} // cap for minimum temperatures at which fit functions below should be used //
    double Zfac = 1.0, dust_to_metals_vs_standard = return_dust_to_metals_ratio_vs_solar(i,T_dust); // avoid call to return_dust_to_metals_ratio_vs_solar to avert circular dependency
#ifdef METALS
    if(i>=0) {Zfac = P[i].Metallicity[0]/All.SolarAbundances[0];}
#endif

    if(dust_or_gas_opacity_only_flag >= 0) // dust opacities
    {
#ifdef RT_INFRARED // use fancy detailed fit with composition varying by dust temperature
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
        
#if defined(RT_INFRARED) || defined(COOL_LOW_TEMPERATURES)
        T_dust_opacitytable = DMIN(T_dust , 1499.9); // limit to <1500 so always use opacities for 'capped' value at 1500 below, but don't ignore, because we're assuming the dust destruction above 1500K is accounted for in the self-consistent calculation of the dust-to-metals ratio, NOT in the opacities here //
#endif
        if(T_dust_opacitytable < 160.) // Tdust < 160 K (all dust constituents present)
        {
            kappa = exp(0.72819004 + 0.75142468*x - 0.07225763*x*x - 0.01159257*x*x*x + 0.00249064*x*x*x*x);
        } else if(T_dust_opacitytable < 275.) { // 160 < Tdust < 275 (no ice present)
            kappa = exp(0.16658241 + 0.70072926*x - 0.04230367*x*x - 0.01133852*x*x*x + 0.0021335*x*x*x*x);
        } else if(T_dust_opacitytable < 425.) { // 275 < Tdust < 425 (no ice or volatile organics present)
            kappa = exp(0.03583845 + 0.68374146*x - 0.03791989*x*x - 0.01135789*x*x*x + 0.00212918*x*x*x*x);
        } else if(T_dust_opacitytable < 680.) { // 425 < Tdust < 680 (silicates, iron, & troilite present)
            kappa = exp(-0.76576135 + 0.57053532*x - 0.0122809*x*x - 0.01037311*x*x*x + 0.00197672*x*x*x*x);
        } else if(T_dust < MAX_DUST_TEMP) { // 680 < Tdust < 1500 (silicates & iron present)
            kappa = exp(-2.23863222 + 0.81223269*x + 0.08010633*x*x + 0.00862152*x*x*x - 0.00271909*x*x*x*x);
        } else {
            kappa = MIN_REAL_NUMBER; // dust completely absent above MAX_DUST_TEMP
        }
        if(dx_excess > 0) {kappa *= exp(0.57*dx_excess);} // assumes kappa scales linearly with temperature (1/lambda) above maximum in fit; pretty good approximation //
    	kappa = DMIN(1.e-3 * Trad * Trad, kappa); // ensure that we extrapolate to low temperatures with a beta=2 law, like in the S03 paper fiducial model
#else
        kappa = DMIN(1.e-3 * Trad * Trad, 5.); // beta=2 law capped at 5 cm^2/g, rough approximation of Semenov model neglecting jumps in composition
#endif
#ifdef RADTRANSFER
        if((do_emission_absorption_scattering_opacity==1) || (do_emission_absorption_scattering_opacity==-1)) {
            kappa *= (1.-0.5/(1.+((725.*725.)/(1.+Trad*Trad)))); /* rough interpolation for dust depending on the radiation temperature: high Trad, this is 1/2, low Trad, gets closer to unity */
        } /* multiply by (1-albedo) because absorption depends only on albedo, and emission cross section depends only on kappa_absorption */
#endif
        kappa *= Zfac*dust_to_metals_vs_standard; // the above are all dust opacities, so they scale with dust content per our usual expressions
    }
    
    if(dust_or_gas_opacity_only_flag <= 0) // non-dust (e.g. gas-phase) IR opacities
    {
        /* this is an approximate result for a wide range of low-to-high-temperature opacities -not- from the dust phase, but provides a pretty good fit from 1.5e3 - 1.0e9 K, and valid at O(1) level down to <10 K, with updates from PFH in Sept 2022 */
        double x_elec = 1., zmetals = 0.014;
#ifdef COOLING
        x_elec = SphP[i].Ne; // actual free electron fraction
#endif
#ifdef METALS
        zmetals = P[i].Metallicity[0];
#endif
        double f_neutral_approx = DMAX(0., 1.-x_elec); /* approximate neutral fraction (good enough for us for what we need below) */
        double f_free_metals_approx = zmetals * DMAX(0, 1.-0.5*dust_to_metals_vs_standard); /* metal mass fraction times the free (not locked in dust abundance), assuming the default solar scaling is 1/2 */
        double Tgas=1. + 0.59*(GAMMA(i)-1.)*U_TO_TEMP_UNITS*SphP[i].InternalEnergyPred, rho_cgs = SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS; /* crude estimate of gas temperature to use with scalings below, and gas density in cgs */
        double k_electron = 0.4 * HYDROGEN_MASSFRAC * x_elec / ((1. + 2.7e11*rho_cgs/(Tgas*Tgas)) * (1. + pow(Trad/4.5e8, 0.86))); /* Thompson scattering (non-relativistic), scaling with free electron fraction [remembering that in our units, x_elec is n_e/n_H_nuclei, not scaled to total nuclear number]; includes corrections for partial degeneracy at low gas temperatures from Buchler et al. 1976, and Klein-Nishina terms at high radiation temperatures >1e9 */
        double k_molecular = 0.1 * (f_free_metals_approx + 3.e-9) * f_neutral_approx; /* molecular line opacities, which should only dominate at low-temperatures in the fits below, but are not really assumed to extrapolate to the very low densities we apply this to here; this works ok comparing e.g. Lenzuni, Chernoff & Salpeter 1991 ApJS 76 759L [opacities for metal free gases], using the 3e-9 to represent the H2 molecular opacity (really low, only here for completeness) */
#if defined(COOL_MOLECFRAC_NONEQM)
        k_molecular *= SphP[i].MolecularMassFraction;
#endif
        double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (f_free_metals_approx * exp(-DMIN(1.5e5/Trad,40.)) + 0.001*x_elec) * rho_cgs / (Trad*Trad*Trad*sqrt(Tgas)); /* free-free, bound-free, bound-bound transitions. bound-bound is small except at discrete wavelengths, so in a mean for a broad-band like we have here, is negligible. the 0.001 term is free-free, independent of metallicity, but note the power of the free electron fraction. the bound-free depends on metal ions here by assumption, specifically those not locked in dust, being ionized -- hence the exponential suppression at low radiation temperatures where the bound states cannot be ionized. the overall Tgas dependence here comes from the sound speed, the Trad from the wavelength (1/nu^3) dependence of the opacity */
        double k_effective_Fe = 1.5e20 * f_free_metals_approx * rho_cgs / (Trad*Trad) * exp(-DMIN(pow(0.8e4/Trad,4),40.)) * exp(-DMIN(pow(Trad/0.7e6,2),40.)); /* crude approximation to the iron line-blanketing opacity calculations from Jiang et al. 2015+2016 */
        k_Kramers += k_effective_Fe;
        double k_Rayleigh = f_neutral_approx * DMIN(5.e-19 * pow(Trad,4) , 0.2*(1.+ HYDROGEN_MASSFRAC)); /* rayleigh scattering from atomic gas [caps at thompson, much lower at low-T here] */
#ifdef COOLING
#ifdef RT_CHEM_PHOTOION
        double x_Hp = SphP[i].HII, x_H0 = SphP[i].HI;
#else
        double u_in=SphP[i].InternalEnergy, rho_in=SphP[i].Density*All.cf_a3inv, mu=1, ne=1, nHI=0, nHII=0, nHeI=1, nHeII=0, nHeIII=0;
        double temp = ThermalProperties(u_in, rho_in, i, &mu, &ne, &nHI, &nHII, &nHeI, &nHeII, &nHeIII);
        double x_Hp = nHII, x_H0 = nHI;
#endif
        double x_Hminus = 4.e-10 * Tgas * x_elec * x_H0 / ((1. + x_Hp*300. + x_elec*1000.*(Tgas/1.3e4)*(Tgas/1.3e4)/(1.+(Tgas/1.3e4)*(Tgas/1.3e4)) + 4.e-17*1.) * (1. + Tgas/3.e4)); /* H- abundance: see series of equations in our non-equilbrium molecular solver (from e.g. Glover and Jappsen 2007 and other sources), with simple but accurate enough for our purposes replacements to make it quick to compute these to the needed accuracy for our purposes. note we need the free-electron fraction, neutral fraction, and free proton fraction. these denominator terms quantify differences from the idealized scaling assumed here, which assumes an idealized scaling of xH0~1~constant and near-vanishing xHp and x_e, for lower temperatures. last term assumes a constant photon-to-baryon ratio for scaling to different environments */
        double k_Hminus_bf = 4.2e7 * pow(8760./Trad, 1.5) * exp(-DMIN(8760./Trad,40.)); /* bound-free H- opacity, from using the fitting functions in John 1988 [A&A, 193, 189], integrating over the Planck function for a flux-mean opacity (Rosseland mean ill-defined here because need all components since this vanishes outside certain ranges) */
        double phi_hm = DMIN(Tgas/5040.,2.), k_Hminus_ff = 1.9e6 * pow(8760./Trad, 2) * exp(-DMIN(8760./Trad,40.)) * (0.6-2.5*sqrt(phi_hm)+2.5*phi_hm+2.7*phi_hm*sqrt(phi_hm)); /* free-free H- opacity, mixing the fits from John and references in Lenzuni, Chernoff, & Salpeter, but re-calculated for arbitrary radiation vs gas temperature. note this will appear to give differences from their opacities, the main difference comes not from this expression (which is simplified) but from the different x_H- and x_e, which owes to a very different chain of expressions, which give a quite different result in the end. */
        double k_Hminus = x_Hminus * (k_Hminus_bf + k_Hminus_ff); /* add both together */
#else
        double k_Hminus = 1.1e-25 * sqrt((zmetals + 1.e-5) * rho_cgs) * pow(Tgas,7.7) * exp(-DMIN(8760./Trad,40.)); /* negative H- ion opacity (this is a fit for stellar atmospheres, which has a very strong temp dependence because of implicit free-electron and H- scaling with T, but that's not as useful for us since we're tracking the chemistry we need here) */
#endif
        double k_radiative = k_molecular + k_Kramers + k_Hminus + k_electron + k_Rayleigh; /* we don't want a rosseland mean here given our band divisions (already kramers and H- and molecular are rosseleand-mean-ized in fact within themselves), but here different sources should add linearly for a flux-mean */
        if((do_emission_absorption_scattering_opacity==1) || (do_emission_absorption_scattering_opacity==-1)) {k_radiative -= k_electron;} /* here we just want absorption/emission, not scattering opacity, so we do not include the free electron term */
        kappa += k_radiative; /* add this to the dust opacity we have already calculated above */
    }
    
    return kappa * fac;
}
