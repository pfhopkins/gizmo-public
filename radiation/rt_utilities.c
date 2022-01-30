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
        lum[RT_FREQ_BIN_GENERIC_USER_FREQ] = (P[i].Mass/1.) * All.Vertical_Grain_Accel * C_LIGHT_CODE * (All.Grain_Internal_Density*All.Grain_Size_Max) * A_base / (0.75*All.Grain_Q_at_MaxGrainSize); // special behavior for particular test of stratified boxes compared to explicit dust opacities
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
    return SphP[i].Interpolated_Opacity[k_freq] + 1.e-3 * All.Dust_to_Gas_Mass_Ratio*0.75*All.Grain_Q_at_MaxGrainSize/(All.Grain_Internal_Density*All.Grain_Size_Max); /* enforce minimum */
#endif
    return MIN_REAL_NUMBER + SphP[i].Interpolated_Opacity[k_freq]; /* this is calculated in a different routine, just return it now */
#endif

#ifdef RT_CHEM_PHOTOION
    /* opacity to ionizing radiation for Petkova & Springel bands. note rt_update_chemistry is where ionization is actually calculated */
    double nH_over_Density = HYDROGEN_MASSFRAC / PROTONMASS_CGS * UNIT_MASS_IN_CGS;
    double kappa = nH_over_Density * (SphP[i].HI + MIN_REAL_NUMBER) * rt_ion_sigma_HI[k_freq];
#if defined(RT_CHEM_PHOTOION_HE) && defined(RT_PHOTOION_MULTIFREQUENCY)
    kappa += nH_over_Density * ((SphP[i].HeI + MIN_REAL_NUMBER) * rt_ion_sigma_HeI[k_freq] + (SphP[i].HeII + MIN_REAL_NUMBER) * rt_ion_sigma_HeII[k_freq]);
    if(k_freq==RT_FREQ_BIN_He0)  {return kappa;}
    if(k_freq==RT_FREQ_BIN_He1)  {return kappa;}
    if(k_freq==RT_FREQ_BIN_He2)  {return kappa;}
#endif
    if(k_freq==RT_FREQ_BIN_H0)  {return kappa;}
#endif

#if defined(RT_HARD_XRAY) || defined(RT_SOFT_XRAY) || defined(RT_PHOTOELECTRIC) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(RT_NUV) || defined(RT_OPTICAL_NIR) || defined(RT_LYMAN_WERNER) || defined(RT_INFRARED) || defined(RT_FREEFREE)
    double fac = UNIT_SURFDEN_IN_CGS, Zfac, dust_to_metals_vs_standard, kappa_HHe; /* units */
    Zfac = 1.0; kappa_HHe=0.35; // assume solar metallicity, simple Thompson cross-section limit for various processes below
    dust_to_metals_vs_standard = return_dust_to_metals_ratio_vs_solar(i); // many of the dust opacities below will need this as the dimensionless dust-to-metals ratio normalized to the canonical Solar value of ~1/2
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
        if(isnan(SphP[i].Dust_Temperature)) {PRINT_WARNING("\n NaN dust temperature for cell-ID=%llu  \n", (unsigned long long) P[i].ID);}
        if(isnan(SphP[i].Radiation_Temperature)) {PRINT_WARNING("\n NaN gas temperature for cell-ID=%llu  \n", (unsigned long long) P[i].ID);}
        if(SphP[i].Dust_Temperature<=T_min) {SphP[i].Dust_Temperature=T_min;} // reset baseline
        if(SphP[i].Radiation_Temperature<=T_min) {SphP[i].Radiation_Temperature=T_min;} // reset baseline
        double T_dust_em = SphP[i].Dust_Temperature; // dust temperature in K //
        double Trad = SphP[i].Radiation_Temperature; // radiation temperature in K //
        if(Trad <= 0) {Trad = 5600.;}
        double kappa = 0.0, x_elec = 1.;
#ifdef COOLING
        x_elec = SphP[i].Ne; // actual free electron fraction
#endif
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
            kappa *= Zfac*dust_to_metals_vs_standard; // the above are all dust opacities, so they scale with dust content per our usual expressions
            kappa += 0.35 * x_elec; // plus Thompson scattering
        } else {
            /* this is an approximate result for high-temperature opacities -not- from the dust phase, but provides a pretty good fit from 1.5e3 - 1.0e9 K */
            double Tgas=10. + 0.59*(GAMMA(i)-1.)*U_TO_TEMP_UNITS*SphP[i].InternalEnergyPred, rho_cgs = SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS; // crude estimate of gas temperature to use with scalings below, and gas density in cgs
            double k_electron = 0.2 * (1. + HYDROGEN_MASSFRAC) * x_elec; /* Thompson scattering (non-relativistic), scaling with free electron fraction */
            double k_molecular = 0.1 * P[i].Metallicity[0] * DMAX(0,1.-x_elec); /* molecular line opacities, which should only dominate at low-temperatures in the fits below, but are not really assumed to extrapolate to the very low densities we apply this to here */
#if defined(COOL_MOLECFRAC_NONEQM)
            return k_molecular *= SphP[i].MolecularMassFraction;
#endif
            double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (P[i].Metallicity[0]+0.001) * rho_cgs / (Trad*Trad*Trad*sqrt(Tgas)); /* free-free, bound-free, bound-bound transitions */
            double k_Hminus = 1.1e-25 * sqrt((P[i].Metallicity[0] + 1.e-5) * rho_cgs) * pow(Tgas,7.7); /* negative H- ion opacity */
            double k_radiative = k_molecular + 1./(1./k_Hminus + 1./(k_electron+k_Kramers)); /* approximate interpolation between the above opacities */
            kappa = DMAX(k_radiative, k_electron); /* set floor at Thompson here */
        }
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
        double fA_tmp = (1.-0.5/(1.+((725.*725.)/(1.+SphP[i].Radiation_Temperature*SphP[i].Radiation_Temperature)))); // rough interpolation depending on the radiation temperature: high Trad, this is 1/2, low Trad, gets closer to unity */
#ifdef COOLING
		if (rt_kappa(i,k_freq)>0) {fA_tmp *= (1.-DMIN(1.,0.35*SphP[i].Ne*fac/rt_kappa(i,k_freq)));} else {return 1.0;} // the value should not matter if rt_kappa=0
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
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge), m_sol = P[i].Mass * UNIT_MASS_IN_SOLAR;
    if((star_age<=0) || isnan(star_age)) {return 0;} // calculate stellar age, will be used below, and catch for bad values

    
    
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
    int active_check = 0; // default to inactive //
    
#if defined(RT_INFRARED) /* special mid-through-far infrared band, which includes IR radiation temperature evolution */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_INFRARED] = stellar_lum_in_band(i, 0, 0.4133);
#endif
#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models; from 0.41-3.4 eV */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_OPTICAL_NIR] = stellar_lum_in_band(i, 0.4133, 3.444);
#endif
#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models; from 3.4-8 eV */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_NUV] = stellar_lum_in_band(i, 3.444, 8.);
#endif
#ifdef RT_PHOTOELECTRIC /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_PHOTOELECTRIC] = stellar_lum_in_band(i, 8, 13.6); // broad band here [note can 2x-count with LW because that is a sub-band, but include it b/c need to total for dust PE heating
#endif
#ifdef RT_LYMAN_WERNER  /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
    SET_ACTIVE_RT_CHECK(); lum[RT_FREQ_BIN_LYMAN_WERNER] = stellar_lum_in_band(i, 11.2, 13.6);
#endif
#if defined(RT_CHEM_PHOTOION)   /* Hydrogen and Helium ionizing bands */
    SET_ACTIVE_RT_CHECK();
#if defined(RT_PHOTOION_MULTIFREQUENCY)
    int i_vec[4] = {RT_FREQ_BIN_H0, RT_FREQ_BIN_He0, RT_FREQ_BIN_He1, RT_FREQ_BIN_He2}; // these will all be the same if not using multi-frequency module //    
    int k; for(k=0;k<3;k++) {lum[i_vec[k]] = stellar_lum_in_band(i, rt_ion_nu_min[i_vec[k]], rt_ion_nu_min[i_vec[k+1]]);} // integrate between band boundaries, defined in global 'nu' in eV
    lum[i_vec[3]] = stellar_lum_in_band(i, rt_ion_nu_min[i_vec[3]], 500.); // integrate to end of spectrum for the last band
#else
    lum[RT_FREQ_BIN_H0] = stellar_lum_in_band(i, 13.6, 500.); // total ionizing flux
#ifdef RT_STARBENCH_TEST
    lum[RT_FREQ_BIN_H0] = 1e49 * (rt_nu_eff_eV[RT_FREQ_BIN_H0]*ELECTRONVOLT_IN_ERGS) / UNIT_LUM_IN_CGS;
#endif    
#endif
#endif

    return active_check;
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
    double E_abs_tot = 0;/* energy absorbed in other bands is transfered to IR, by default: track it here */
    double Rad_E_gamma_tot = 0; // dust temperature defined by total radiation energy density //
    {int j; for(j=0;j<N_RT_FREQ_BINS;j++) {Rad_E_gamma_tot += SphP[i].Rad_E_gamma[j];}}
    double a_rad_inverse=C_LIGHT_CGS/(4.*5.67e-5), vol_inv_phys=(SphP[i].Density*All.cf_a3inv/P[i].Mass), u_gamma = Rad_E_gamma_tot * vol_inv_phys * UNIT_PRESSURE_IN_CGS; // photon energy density in CGS //
    double Dust_Temperature_4 = u_gamma * a_rad_inverse; // estimated effective temperature of local rad field in equilibrium with dust emission. note that for our definitions, rad energy density has its 'normal' value independent of RSOL, so Tdust should as well; emission -and- absorption are both lower by a factor of RSOL, but these cancel in the Tdust4 here //
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
                    double dE_abs = -e0 * (1. - exp(a0*dt_entr)); // change in energy from absorption
                    double T_max = DMAX(SphP[i].Radiation_Temperature , dE_fac / dTE_fac); 
                    double rfac=1; if(dE_fac < -0.5*(e0+dE_abs)) {rfac=fabs(0.5*(e0+dE_abs))/fabs(dE_fac);} else {if(dE_fac > 0.5*e0) {rfac=0.5*e0/dE_fac;}}
                    dE_fac*=rfac; dTE_fac*=rfac; // limit temperature change from advection to prevent spurious divergences
                    
                    SphP[i].Radiation_Temperature = (e0 + dE_fac) / (MIN_REAL_NUMBER + DMAX(0., e0 / SphP[i].Radiation_Temperature + dTE_fac));
                    SphP[i].Radiation_Temperature = DMIN(SphP[i].Radiation_Temperature, T_max);
                    a0 = -rt_absorption_rate(i,kf); // update absorption rate using the new radiation temperature //
                }
                double total_emission_rate = E_abs_tot + fabs(a0)*e0 + SphP[i].Rad_Je[kf]; // add the summed absorption as emissivity here //
                total_de_dt = E_abs_tot + SphP[i].Rad_Je[kf] + dt_e_gamma_band;
                if(fabs(a0)>0)
                {
                    Dust_Temperature_4 = total_emission_rate * vol_inv_phys / (MIN_REAL_NUMBER + fabs(a0)); // physical energy density units, both emission and absorption rates reduced by one power of RSOL so cancel
                    Dust_Temperature_4 *= UNIT_PRESSURE_IN_CGS * a_rad_inverse; // convert to cgs; physical a_r here
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
#ifdef RT_INFRARED
            if(isnan(ef)) {PRINT_WARNING("\n ef energy prediction is NaN for cell-ID=%llu, e0=%g e_abs_0=%g abs_0=%g a0=%g total_de_dt=%g dt_entr=%g slabfac=%g Trad=%g Tdust=%g\n", (unsigned long long) P[i].ID,e0, e_abs_0,abs_0, a0, total_de_dt,dt_entr,slabfac,SphP[i].Radiation_Temperature,SphP[i].Dust_Temperature);}
#else
            if(isnan(ef)) {PRINT_WARNING("\n ef energy prediction is NaN for cell-ID=%llu, e0=%g e_abs_0=%g abs_0=%g a0=%g total_de_dt=%g dt_entr=%g slabfac=%g\n", (unsigned long long) P[i].ID,e0, e_abs_0,abs_0, a0, total_de_dt,dt_entr,slabfac);}
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
            if(donation_target_bin==RT_FREQ_BIN_INFRARED) {E_abs_tot += de_abs/(MIN_REAL_NUMBER + dt_entr);} /* donor bin is yourself in the IR - all self-absorption is re-emitted */
            if(kf==RT_FREQ_BIN_INFRARED) {ef = e0 + total_de_dt * dt_entr;} /* donor bin is yourself in the IR - all self-absorption is re-emitted */
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
            for(k_dir=0;k_dir<3;k_dir++) {DeltaFluxEff[k_dir] -= (P[i].Mass/rho) * (C_LIGHT_CODE_REDUCED*C_LIGHT_CODE_REDUCED/teqm_inv) * SphP[i].Gradients.Rad_E_gamma_ET[kf][k_dir];} // here we compute the nabla.pressure_gradient_tensor term from gradients directly, and use this in the next step after multiplying the flux equation by (tilde[c]^2/dt_eqm_inv) and working in dimensionless time units
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
            if(k==RT_FREQ_BIN_INFRARED){SphP[i].Radiation_Temperature = DMIN(All.InitGasTemp,100.);}
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
        if(k==RT_FREQ_BIN_INFRARED){urad[k] = (RT_ISRF_BACKGROUND * 0.39 + 0.26) * ELECTRONVOLT_IN_ERGS / UNIT_PRESSURE_IN_CGS;} // 0.33 eV/cm^3 is dust emission peak, 0.26 is CMB - note how this bin actually lumps the two together
#endif
#ifdef RT_OPTICAL_NIR
        if(k==RT_FREQ_BIN_OPTICAL_NIR){urad[k] = RT_ISRF_BACKGROUND * 0.54 * ELECTRONVOLT_IN_ERGS / UNIT_PRESSURE_IN_CGS;} // stellar emission
#endif
#ifdef RT_NUV
        if(k==RT_FREQ_BIN_NUV){urad[k] = RT_ISRF_BACKGROUND * 0.024 * ELECTRONVOLT_IN_ERGS / UNIT_PRESSURE_IN_CGS;} // stellar emission
#endif
#ifdef RT_PHOTOELECTRIC
        if(k==RT_FREQ_BIN_PHOTOELECTRIC){urad[k] = RT_ISRF_BACKGROUND * 1.7 * 3.9e-14 / UNIT_PRESSURE_IN_CGS;} // Draine 1978 value = 1.7 Habing
#endif
    }
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
            SphP[i].Dust_Temperature = DMIN(All.InitGasTemp,100.); //get_min_allowed_dustIRrad_temperature(); // in K, floor = CMB temperature or 10K
            SphP[i].Radiation_Temperature = DMIN(All.InitGasTemp,100.); //SphP[i].Dust_Temperature;
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
#ifdef RT_ISRF_BACKGROUND
            double urad[N_RT_FREQ_BINS];
            get_background_isrf_urad(i, urad);
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
                if(RestartFlag==0 && k==RT_FREQ_BIN_INFRARED){ // only initialize the IR energy if starting a new run, otherwise use what's in the snapshot
                    SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED] = (4.*5.67e-5 / C_LIGHT_CGS) * pow(DMIN(All.InitGasTemp,100.),4.) / UNIT_PRESSURE_IN_CGS * P[i].Mass / (SphP[i].Density*All.cf_a3inv);
                }
#endif
#ifdef RT_ISRF_BACKGROUND
                if(RestartFlag == 0) {SphP[i].Rad_E_gamma[k] = urad[k] * P[i].Mass / (SphP[i].Density*All.cf_a3inv);}
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
                double q_a = (0.75*All.Grain_Q_at_MaxGrainSize) / (All.Grain_Internal_Density*All.Grain_Size_Max), e0 = All.Vertical_Grain_Accel / q_a, kappa_0 = All.Grain_Absorbed_Fraction_vs_Total_Extinction * q_a * All.Dust_to_Gas_Mass_Ratio, cell_vol = (P[i].Mass/SphP[i].Density);
                double rho_base_setup = 1., H_scale_setup = 1.; // define in code units the -assumed- initial scaling of the base gas density and vertical scale-length (PROBLEM SPECIFIC HERE!)
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
double get_rt_ir_lambdadust_effective(double T, double rho, double *nH0_guess, double *ne_guess, int target, int update_Tdust)
{
#ifdef COOLING
    double Tdust_0 = SphP[target].Dust_Temperature; // dust temperature estimate from previous loop over radiation operators
    double Tdust = Tdust_0;
    double egy_therm = SphP[target].InternalEnergyPred*P[target].Mass; // true internal energy (before this cooling loop)
    double egy_rad = (C_LIGHT_CODE/C_LIGHT_CODE_REDUCED) * SphP[target].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED]; // effective radiation field energy (before this cooling loop), accounting for the effects of an RSOL on the difference between the gas and radiation field energy equations
    double egy_tot = egy_rad + egy_therm; // true total energy [in code units]
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS_CGS; // effective hydrogen number dens in cgs units (just for normalization convention)
    double volume = (P[target].Mass / (SphP[target].Density*All.cf_a3inv)); // particle volume in code units
    double ratefact = (nHcgs*nHcgs) * volume / (UNIT_PRESSURE_IN_CGS / UNIT_TIME_IN_CGS); // conversion b/t Lambda and du used by code
    double Erad_to_T4_fac = (C_LIGHT_CODE_REDUCED/C_LIGHT_CODE) * (C_LIGHT_CGS/(4. * 5.67e-5)) * UNIT_PRESSURE_IN_CGS / volume; // conversion from absolute rad energy to T^4 units, used multiple places below, coefficient = cL_reduced/(4*sigma_B); note RSOL power above cancels here b/c of the cancellation in absorption+emission coefficients
    double Teff = Get_Gas_Mean_Molecular_Weight_mu(T, rho, nH0_guess, ne_guess, 0, target) * (GAMMA(target)-1.) * U_TO_TEMP_UNITS * (egy_tot / P[target].Mass); // convert from internal energy to temperature units for factor below

    double xf, a = Teff*Teff*Teff*Teff / (Erad_to_T4_fac*egy_tot); // dimensionless factors needed to solve for the equilibrium Tdust-Tgas relation; this assumes total 'effective' energy between gas and radiation is conserved, and Td=Tgas, and solves for what the equilibrium temperature would be in terms of xf = Teqm / Teff
    if(a<0.2138) {xf=(1+19*a+132*a*a+418*a*a*a+580*a*a*a*a+243*a*a*a*a*a)/(1+20*a+148*a*a+508*a*a*a+796*a*a*a*a+432*a*a*a*a*a);} // eqm solution (power series approx)
     else {double a0=pow(a,0.25); xf=(-704-1045*a0+128*a0*a0*a0*(39+32*a0*(4+7*a0+64*a0*a0*a0*(-1+8*a0*(-1+4*a0)))))/(8388608.*a*a*a0*a0);} // eqm solution (power series approx)

    double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target); // timestep being taken [code units]
    double LambdaDust_initial_guess, lambda_eff, L0_abs, Edot0, efinal_minus_einitial, t_cooling_eff, sign_term, tau, xfac, lambda_fac=1.116e-32 * sqrt(T)*(1.-0.8*exp(-75./T)) * (P[target].Metallicity[0]/All.SolarAbundances[0]) * return_dust_to_metals_ratio_vs_solar(target); int iter=0;
    do // LambdaDust implicitly depends nonlinearly on Tdust the way we have this set up here, so we do fixed-point iteration to convergence - typically only a few iters needed
    {
      LambdaDust_initial_guess = lambda_fac * (Tdust-T); // guess value based on the -current- values of T, Tdust //
      L0_abs = fabs(LambdaDust_initial_guess); // absolute value of the initially-computed guess for the cooling/heating rate of the gas
      Edot0 = L0_abs * ratefact; // now this is an absolute Edot in code units, for the gas loss/gain from dust
      efinal_minus_einitial = egy_tot*xf - egy_therm; // change in gas thermal energy if we went all the way to equilibrium
      t_cooling_eff = fabs(efinal_minus_einitial) / Edot0; // effective cooling time at the initially-estimated rate here: we'll use this to stably interpolate
      sign_term=1.; if(efinal_minus_einitial < 0.) {sign_term=-1.;} // sign of the cooling/heating (to keep for below)
      tau = dt/t_cooling_eff;
      xfac=(1.-exp(-tau))/tau; if(tau<0.05) {xfac=1.-0.5*tau+tau*tau/6.;} else {if(tau>20.) {xfac=1./tau;}} // correct rate to asymptote to equilibrium
      lambda_eff = sign_term * L0_abs * xfac; // final effective gas cooling/heating rate
      Tdust_0 = Tdust;
      Tdust = DMAX(pow(Erad_to_T4_fac*DMAX( 0., egy_rad - lambda_eff*ratefact*dt ), 0.25), get_min_allowed_dustIRrad_temperature());
      iter += 1;
    } while ((fabs(Tdust - Tdust_0) > 1e-14 * Tdust) && (iter<MAXITER));

    if(update_Tdust) {SphP[target].Dust_Temperature = Tdust;} //DMAX(pow(Erad_to_T4_fac*DMAX( 0., egy_rad - lambda_eff*ratefact*dt ), 0.25), get_min_allowed_dustIRrad_temperature());} // update dust temperature guess //
    return lambda_eff;
#endif
    return 0;
}

#endif



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
      f_lower = 1 - (0.15398973382026504*(6. + x1*(6. + x1*(3. + x1))))*exp(-x1); // approximation of Planck integral for large x
    }
    if(x2 < 3.40309){
      f_upper = (131.4045728599595*x2*x2*x2)/(2560. + x2*(960. + x2*(232. + 39.*x2))); // approximation of integral of Planck function from 0 to x2, valid for x2 << 1
    } else {
      f_upper = 1 - (0.15398973382026504*(6. + x2*(6. + x2*(3. + x2))))*exp(-x2); // approximation of Planck integral for large x
    }
    return DMAX(f_upper - f_lower, 0);
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
        int age_bin, j; double age_Myr=1000.*evaluate_stellar_age_Gyr(P[i].StellarAge), log_age_Myr=log10(age_Myr), stellar_mass=P[i].Mass*UNIT_MASS_IN_SOLAR;
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
