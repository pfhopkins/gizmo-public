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

/* Routines for models that require stellar evolution: luminosities, mass loss, SNe rates, etc. 
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifdef GALSF


/* return the (solar-scaled) light-to-mass ratio of an SSP with a given age; used throughout the code */
double evaluate_light_to_mass_ratio(double stellar_age_in_gyr, int i)
{
    double lum=1; if(stellar_age_in_gyr < 0.01) {lum=1000;} // default to a dumb imf-averaged 'young/high-mass' vs 'old/low-mass' distinction 
#ifdef SINGLE_STAR_SINK_DYNAMICS // calculate single-star luminosity (and convert to solar luminosity-to-mass ratio, which this output assumes) 
    lum=calculate_individual_stellar_luminosity(0, P[i].BH_Mass, i) / P[i].BH_Mass * (All.UnitEnergy_in_cgs / (All.UnitTime_in_s * SOLAR_LUM)) / (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS));
#endif
    if(stellar_age_in_gyr<0.033) {lum*=calculate_relative_light_to_mass_ratio_from_imf(stellar_age_in_gyr,i);} // account for IMF variation model [if used]
    return lum;
}


/* subroutine to calculate luminosity of an individual star, according to accretion rate, 
    mass, age, etc. Modify your assumptions about main-sequence evolution here. */
double calculate_individual_stellar_luminosity(double mdot, double mass, long i)
{
    double lum = 0;
#ifdef SINGLE_STAR_SINK_DYNAMICS
    double c_code = C_LIGHT_CODE;
    double m_solar = mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
    /* if below the deuterium burning limit, just use the potential energy efficiency at the surface of a jupiter-density object */
    double rad_eff_protostar = 5.0e-7;
    if(m_solar < 0.012) {rad_eff_protostar = 5.e-8 * pow(m_solar/0.00095,2./3.);}
    lum = rad_eff_protostar * mdot * c_code*c_code;
    /* now for pre-main sequence, need to also check the mass-luminosity relation */
    double lum_sol = 0;
    lum_sol *= SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s);
    lum += lum_sol;
#endif    
    return lum;
}


/* return the light-to-mass ratio, for the IMF of a given particle, relative to the Chabrier/Kroupa IMF */
double calculate_relative_light_to_mass_ratio_from_imf(double stellar_age_in_gyr, int i)
{
#ifdef GALSF_SFR_IMF_VARIATION // fitting function from David Guszejnov's IMF calculations (ok for Mturnover in range 0.01-100) for how mass-to-light ratio varies with IMF shape/effective turnover mass 
    double log_mimf = log10(P[i].IMF_Mturnover);
    return (0.051+0.042*(log_mimf+2)+0.031*(log_mimf+2)*(log_mimf+2)) / 0.31;
#endif
#ifdef GALSF_SFR_IMF_SAMPLING // account for IMF sampling model if not evolving individual stars
    double mu = 0.01 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam / (1.989e33); // 1 O-star per 100 Msun
    if(stellar_age_in_gyr > 0.003) {mu *= 0.326 * (0.003 / stellar_age_in_gyr);} // expectation value is declining with time, so 'effective multiplier' is larger
    return P[i].IMF_NumMassiveStars / mu;
#endif
    return 1; // Chabrier or Kroupa IMF //
}


#if defined(FLAG_NOT_IN_PUBLIC_CODE) || (defined(RT_CHEM_PHOTOION) && defined(GALSF))
double particle_ionizing_luminosity_in_cgs(long i)
{
    double lm_ssp=0;
    if(P[i].Type != 5)
    {
        /* use updated SB99 tracks: including rotation, new mass-loss tracks, etc. */
        double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if(star_age >= 0.02) return 0; // skip since old stars don't contribute
        if(star_age < 0.0035) {lm_ssp=500.;} else {double log_age=log10(star_age/0.0035); lm_ssp=470.*pow(10.,-2.24*log_age-4.2*log_age*log_age) + 60.*pow(10.,-3.6*log_age);}
        lm_ssp *= calculate_relative_light_to_mass_ratio_from_imf(star_age, i);
#ifdef SINGLE_STAR_SINK_DYNAMICS /* use effective temperature as a function of stellar mass and size to get ionizing photon production */
        double l_sol = bh_lum_bol(0,P[i].Mass,i) * (All.UnitEnergy_in_cgs / (All.UnitTime_in_s * SOLAR_LUM)); // L/Lsun
        double m_sol = P[i].Mass * (All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS)); // M/Msun
        double r_sol = pow(m_sol, 0.738); // R/Rsun
        double T_eff = 5780. * pow(l_sol/(r_sol*r_sol), 0.25); // ZAMS effective temperature
        double x0 = 157800./T_eff; // h*nu/kT for nu>13.6 eV
        double fion = 0.0; // fraction of blackbody emitted above x0
        if(x0 < 30.) {double q=18./(x0*x0) + 1./(8. + x0 + 20.*exp(-x0/10.)); fion = exp(-1./q);} // accurate to <10% for a Planck spectrum to x0>30, well into vanishing flux //
        lm_ssp = fion * l_sol / m_sol; // just needs to be multiplied by the actual stellar luminosity to get luminosity to mass ratio
#endif
    } // (P[i].Type != 5)
    lm_ssp *= (1.95*P[i].Mass*All.UnitMass_in_g/All.HubbleParam); // convert to luminosity from L/M
    if((lm_ssp <= 0) || (!isfinite(lm_ssp))) {lm_ssp=0;} // trap for negative values and nans (shouldnt happen)
    return lm_ssp;
}
#endif









/* this routine tells the feedback algorithms what to 'inject' when a stellar feedback event occurs.
    you must define the mass, velocity (which defines the momentum and energy), and metal content (yields)
    of the ejecta for the event[s] of interest. Mass [Msne] and velocity [SNe_v_ejecta] should
    be in code units. yields[k] should be defined for all metal species [k], and in dimensionless units
    (mass fraction of the ejecta in that species). */
void particle2in_addFB_fromstars(struct addFB_evaluate_data_in_ *in, int i, int fb_loop_iteration)
{
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
    if(P[i].SNe_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    // 'dummy' example model assumes all SNe are identical with IMF-averaged properties from the AGORA model (Kim et al., 2016 ApJ, 833, 202)
    in->Msne = P[i].SNe_ThisTimeStep * (14.8*SOLAR_MASS)/(All.UnitMass_in_g/All.HubbleParam); // assume every SNe carries 14.8 solar masses (IMF-average)
    in->SNe_v_ejecta = 2.607e8 / All.UnitVelocity_in_cm_per_s; // assume ejecta are ~2607 km/s [KE=1e51 erg, for M=14.8 Msun]
#ifdef METALS
    int k; for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.178*All.SolarAbundances[k]/All.SolarAbundances[0];} // assume a universal solar-type yield with ~2.63 Msun of metals
    if(NUM_METAL_SPECIES>=10) {in->yields[1] = 0.4;} // (catch for Helium, which the above scaling would give bad values for)
#endif
#endif
}


/* this routine calculates the event rates for different types of mechanical/thermal feedback
    algorithms. things like SNe rates, which determine when energy/momentum/mass are injected, should go here.
    you can easily modify this to accomodate any form of thermal or mechanical feedback/injection of various
    quantities from stars. */
double mechanical_fb_calculate_eventrates(int i, double dt)
{
    double RSNe = 0;
#ifdef GALSF_FB_THERMAL
    // pure thermal feedback: assumes AGORA model (Kim et al., 2016 ApJ, 833, 202) where everything occurs at 5Myr exactly //
    if(P[i].SNe_ThisTimeStep != 0) {P[i].SNe_ThisTimeStep=-1; return 0;} // already had an event, so this particle is "done"
    if(evaluate_stellar_age_Gyr(P[i].StellarAge) < 0.005) {return 0;} // enforce age limit of 5 Myr
    P[i].SNe_ThisTimeStep = P[i].Mass * (All.UnitMass_in_g / All.HubbleParam) / (91. * SOLAR_MASS); // 1 event per 91 solar masses
    return 1;
#endif
#ifdef GALSF_FB_MECHANICAL
    // mechanical feedback: 'dummy' example model below assumes a constant SNe rate for t < 30 Myr, then nothing. experiment! //
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    if(star_age < 0.03)
    {
        RSNe = 3.e-4; // assume a constant rate ~ 3e-4 SNe/Myr/solar mass for t = 0-30 Myr //
        double p = dt * RSNe * P[i].Mass * (All.UnitTime_in_Megayears/All.HubbleParam) * (All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS; // unit conversion factor
        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
        P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
    }
#endif
    return RSNe;
}




    
#endif /* GALSF */
