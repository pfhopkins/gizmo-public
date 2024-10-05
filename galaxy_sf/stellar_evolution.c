#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
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


/* function to say whether a given particle should be treated as a "single star" or as -eligible to become- a single star (if it is gas), or
    whether it should be treated as a stellar population with some total mass, with IMF-integrated properties */
int is_particle_single_star_eligible(long i)
{
#if defined(SINGLE_STAR_SINK_DYNAMICS)
    if(P[i].Type == 0 || P[i].Type == 5) // only type=0 or type=5 (sinks) are eligible in our models here to be 'single star' candidates
    {
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL) // here's the interesting regime, where we have some criterion for deciding which cells are eligible for 'single-star' status
        if(P[i].Type == 5) {return 1;} // all type-5 elements are assumed sinks
        if(P[i].Type == 0) {if(P[i].Mass*UNIT_MASS_IN_SOLAR > (SINGLE_STAR_AND_SSP_HYBRID_MODEL)) {return 0;} else {return 1;}} // use a simple mass threshold to decide which model we will use, specified by using this as a compile-time parameter
#else
        return 1; // no hybrid model, so all particles satisfying these criteria are automatically single-star eligible
#endif
    }
#endif
    return 0; // catch - default to non-single-star (SSP), unless satisfy some of the criteria above
}


/* return the light-to-mass ratio [in units of Lsun/Msun] of a star or stellar population with a given age; used throughout the code below */
double evaluate_light_to_mass_ratio(double stellar_age_in_gyr, int i)
{
    if(is_particle_single_star_eligible(i)) // SINGLE-STAR VERSION: calculate single-star luminosity (and convert to solar luminosity-to-mass ratio, which this output assumes)
    {
#ifdef SINGLE_STAR_SINK_DYNAMICS
        double m0=P[i].Mass; if(P[i].Type == 5) {m0=P[i].BH_Mass;}
        return calculate_individual_stellar_luminosity(0, m0, i) / m0 * (UNIT_LUM_IN_SOLAR) / (UNIT_MASS_IN_SOLAR);
#endif
    }
    else // STELLAR-POPULATION VERSION: compute integrated mass-to-light ratio of an SSP
    {
        double lum=1; if(stellar_age_in_gyr < 0.01) {lum=1000;} // default to a dumb imf-averaged 'young/high-mass' vs 'old/low-mass' distinction
#ifdef GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF
        lum *= calculate_relative_light_to_mass_ratio_from_imf(stellar_age_in_gyr,i,-1); // account for IMF variation model [if used; currently must be custom set as desired for modules]
#endif
        return lum;
    }
    return 0; // catch
}


/* subroutine to calculate luminosity of an individual star, according to accretion rate,
    mass, age, etc. Modify your assumptions about main-sequence evolution here. ONLY relevant for SINGLE-STAR inputs. */
double calculate_individual_stellar_luminosity(double mdot, double mass, long i)
{
#if !defined(SINGLE_STAR_SINK_DYNAMICS)
    return 0; /* not defined */
#endif
    /* if above flags not defined, estimate accretion + main-sequence luminosity as simply as possible */
    double lum=0, lum_sol=0, c_code = C_LIGHT_CODE, m_solar = mass * UNIT_MASS_IN_SOLAR;
    double rad_eff_protostar = 5.0e-7; /* if below the deuterium burning limit, just use the potential energy efficiency at the surface of a jupiter-density object */
    if(m_solar < 0.012) {rad_eff_protostar = 5.e-8 * pow(m_solar/0.00095,2./3.);}
    lum = rad_eff_protostar * mdot * c_code*c_code;
    if(m_solar >= 0.012) /* now for pre-main sequence and main sequence, need to also check the mass-luminosity relation */
    {
        if(m_solar < 0.43) {lum_sol = 0.185 * m_solar*m_solar;}
        else if(m_solar < 2.) {lum_sol = m_solar*m_solar*m_solar*m_solar;}
        else if(m_solar < 53.9) {lum_sol = 1.5 * m_solar*m_solar*m_solar * sqrt(m_solar);}
        else {lum_sol = 32000. * m_solar;}
    }
    return lum;

}


/* return the light-to-mass ratio, for the IMF of a given particle, relative to the Chabrier/Kroupa IMF.
    ONLY relevant for STELLAR POPULATION integrated inputs. "mode" denotes if we're interested in very massive stars (or other special behaviors): for now -1=bolometric, +1=very massive */
double calculate_relative_light_to_mass_ratio_from_imf(double stellar_age_in_gyr, int i, int mode)
{
#ifdef GALSF_SFR_IMF_VARIATION // fitting function from David Guszejnov's IMF calculations (ok for Mturnover in range 0.01-100) for how mass-to-light ratio varies with IMF shape/effective turnover mass
    double log_mimf = log10(P[i].IMF_Mturnover);
    return (0.051+0.042*(log_mimf+2)+0.031*(log_mimf+2)*(log_mimf+2)) / 0.31;
#endif
#ifdef GALSF_SFR_IMF_SAMPLING // account for IMF sampling model if not evolving individual stars
    double mu = 0.0115 * P[i].Mass * UNIT_MASS_IN_SOLAR; // 1 O-star per 100 Msun [more exactly calculated here as number of stars per solar mass with mass > 8 Msun, from our adopted Kroupa IMF from 0.01-100 Msun]
    double t=stellar_age_in_gyr*1000.,t1=3.7,t2=7.,t3=44.,a0=0.13,mu_min=1.e-3*mu;
    if(t>t3) {mu*=0;} else {if(t>t2) {mu*=(1.-a0)*(1.-(t-t2)/(t3-t2));} else {if(t>t1) {mu*=1.-a0*(t-t1)/(t2-t1);}}} // expectation value is declining with time, so 'effective multiplier' is larger
    if(mode > 0) {
        return P[i].IMF_NumMassiveStars / DMAX(mu,mu_min); // scales just with the number of massive stars
    } else {
        if(t>=t3) {return 1;} else {return 0.01 + P[i].IMF_NumMassiveStars / DMAX(mu,mu_min);} // scales with a minimum and a long-timescale baseline
    }
#endif
    return 1; // Chabrier or Kroupa IMF //
}


#if defined(GALSF_FB_FIRE_RT_HIIHEATING) || (defined(RT_CHEM_PHOTOION) && defined(GALSF))
/* routine to compute the -ionizing- luminosity coming from either individual stars or an SSP */
double particle_ionizing_luminosity_in_cgs(long i)
{
    if(is_particle_single_star_eligible(i)) /* SINGLE STAR VERSION: use effective temperature as a function of stellar mass and size to get ionizing photon production */
    {
#ifdef SINGLE_STAR_SINK_DYNAMICS
        double l_sol=bh_lum_bol(0,P[i].Mass,i)*(UNIT_LUM_IN_SOLAR), m_sol=P[i].Mass*UNIT_MASS_IN_SOLAR, r_sol=pow(m_sol,0.738); // L/Lsun, M/Msun, R/Rsun
        double T_eff=5780.*pow(l_sol/(r_sol*r_sol),0.25), x0=157800./T_eff, fion=0; // ZAMS effective temperature; x0=h*nu/kT for nu>13.6 eV; fion=fraction of blackbody emitted above x0
        if(x0 < 30.) {double q=18./(x0*x0) + 1./(8. + x0 + 20.*exp(-x0/10.)); fion = exp(-1./q);} // accurate to <10% for a Planck spectrum to x0>30, well into vanishing flux //
        return fion * l_sol * SOLAR_LUM_CGS; // return value in cgs, as desired for this routine [l_sol is in L_sun, by definition above] //
#endif
    }
    else /* STELLAR POPULATION VERSION: use updated SB99 tracks: including rotation, new mass-loss tracks, etc. */
    {
        if(P[i].Type != 5)
        {
            double lm_ssp=0, star_age=evaluate_stellar_age_Gyr(i), t0=0.0035, tmax=0.02;
            if(star_age < t0) {lm_ssp=500.;} else {double log_age=log10(star_age/t0); lm_ssp=470.*pow(10.,-2.24*log_age-4.2*log_age*log_age) + 60.*pow(10.,-3.6*log_age);}
            if(star_age < 0.033) {lm_ssp *= 1.e-4 + calculate_relative_light_to_mass_ratio_from_imf(star_age,i,1);}
            if(star_age >= tmax) {return 0;} // skip since old stars don't contribute
            return lm_ssp * SOLAR_LUM_CGS * (P[i].Mass*UNIT_MASS_IN_SOLAR); // converts to cgs luminosity [lm_ssp is in Lsun/Msun, here]
        } // (P[i].Type != 5)
#ifdef BH_HII_HEATING /* AGN template: light-to-mass ratio L(>13.6ev)/Mparticle in Lsun/Msun, above is dNion/dt = 5.5e54 s^-1 (Lbol/1e45 erg/s) */
        if(P[i].Type == 5) {return 0.18 * bh_lum_bol(P[i].BH_Mdot,P[i].Mass,i) * UNIT_LUM_IN_CGS;}
#endif
    }
    return 0; // catch
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
    in->Msne = P[i].SNe_ThisTimeStep * (14.8/UNIT_MASS_IN_SOLAR); // assume every SNe carries 14.8 solar masses (IMF-average)
    in->SNe_v_ejecta = 2607. / UNIT_VEL_IN_KMS; // assume ejecta are ~2607 km/s [KE=1e51 erg, for M=14.8 Msun], which is IMF-averaged
#ifdef SINGLE_STAR_SINK_DYNAMICS // if single-star exploding or returning mass, use its actual mass & assumed energy to obtain the velocity
    if(is_particle_single_star_eligible(i))
    {
        in->Msne = DMIN(1.,P[i].SNe_ThisTimeStep) * P[i].Mass; // mass fraction of star being returned this timestep
        in->SNe_v_ejecta = sqrt(2.*(1.e51/UNIT_ENERGY_IN_CGS)/P[i].Mass); // for SNe [total return], simple v=sqrt(2E/m)should be fine without relativistic corrections
        if(P[i].SNe_ThisTimeStep<1) {double m_msun=P[i].Mass*UNIT_MASS_IN_SOLAR; in->SNe_v_ejecta = (616. * sqrt((1.+0.1125*m_msun)/(1.+0.0125*m_msun)) * pow(m_msun,0.131)) / UNIT_VEL_IN_KMS;} // scaling from size-mass relation+eddington factor, assuming line-driven winds //
    }
#endif
#ifdef METALS
    int k; for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.178*All.SolarAbundances[k]/All.SolarAbundances[0];} // assume a universal solar-type yield with ~2.63 Msun of metals
    if(NUM_LIVE_SPECIES_FOR_COOLTABLES>=10) {in->yields[1] = 0.4;} // (catch for Helium, which the above scaling would give bad values for)
#endif
#endif
}


/* this routine calculates the event rates for different types of mechanical/thermal feedback
    algorithms. things like SNe rates, which determine when energy/momentum/mass are injected, should go here.
    you can easily modify this to accomodate any form of thermal or mechanical feedback/injection of various
    quantities from stars. */
double mechanical_fb_calculate_eventrates(int i, double dt)
{

#if defined(SINGLE_STAR_SINK_DYNAMICS) && !defined(FLAG_NOT_IN_PUBLIC_CODE) /* SINGLE-STAR version: simple implementation of single-star wind mass-loss and SNe rates */
    if(is_particle_single_star_eligible(i))
    {
        double m_sol,l_sol; m_sol=P[i].Mass*UNIT_MASS_IN_SOLAR; l_sol=bh_lum_bol(0,P[i].Mass,i)*UNIT_LUM_IN_SOLAR;
#ifdef SINGLE_STAR_FB_WINDS
        double gam=DMIN(0.5,3.2e-5*l_sol/m_sol), alpha=0.5+0.4/(1.+16./m_sol), q0=(1.-alpha)*gam/(1.-gam), k0=1./30.; // Eddington factor (~L/Ledd for winds), capped at 1/2 for sanity reasons, approximate scaling for alpha factor with stellar type (weak effect)
        P[i].SNe_ThisTimeStep = DMIN(0.5, (2.338 * alpha * pow(l_sol,7./8.) * pow(m_sol,0.1845) * (1./q0) * pow(q0*k0,1./alpha) / m_sol) * (dt*UNIT_TIME_IN_GYR)); // Castor, Abbot, & Klein scaling
#endif
#ifdef SINGLE_STAR_FB_SNE
        double t_lifetime_Gyr = 10.*(m_sol/l_sol) + 0.003; /* crude estimate of main-sequence lifetime, capped at 3 Myr*/
        if(evaluate_stellar_age_Gyr(i) >= t_lifetime_Gyr) {P[i].SNe_ThisTimeStep=1;}
#endif
        return 1;
    }
#endif

#ifdef GALSF_FB_THERMAL /* STELLAR-POPULATION version: pure thermal feedback: assumes AGORA model (Kim et al., 2016 ApJ, 833, 202) where everything occurs at 5Myr exactly */
    if(P[i].SNe_ThisTimeStep != 0) {P[i].SNe_ThisTimeStep=-1; return 0;} // already had an event, so this particle is "done"
    if(evaluate_stellar_age_Gyr(i) < 0.005) {return 0;} // enforce age limit of 5 Myr
    P[i].SNe_ThisTimeStep = P[i].Mass*UNIT_MASS_IN_SOLAR / 91.; // 1 event per 91 solar masses
    return 1;
#endif

#ifdef GALSF_FB_MECHANICAL /* STELLAR-POPULATION version: mechanical feedback: 'dummy' example model below assumes a constant SNe rate for t < 30 Myr, then nothing. experiment! */
    double star_age = evaluate_stellar_age_Gyr(i);
    if(star_age < 0.03)
    {
        double RSNe = 3.e-4; // assume a constant rate ~ 3e-4 SNe/Myr/solar mass for t = 0-30 Myr //
        double p = RSNe * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
        double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
        P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
        return RSNe;
    }
#endif

    return 0;
}






double Z_for_stellar_evol(int i)
{
    if(i<0) {return 1;}
#ifdef METALS
    double Z_solar = P[i].Metallicity[0]/All.SolarAbundances[0]; // use total metallicity
    return DMIN(DMAX(Z_solar,0.01),3.); // stellar evolution tables here are not arbitrarily extrapolable, so this is bounded //
#else
    return 1; // metals not evolved, return unity
#endif
}


#if defined(GALSF_SFR_IMF_SAMPLING_DISTRIBUTE_SF)
void update_stellarnumber_and_timedistribofstarformation(void)
{
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 4)
        {
            double dt_since_form_code = evaluate_time_since_t_initial_in_Gyr(P[i].StellarAge) / UNIT_TIME_IN_GYR; // time since spawn in code [physical] units
            if((dt_since_form_code < P[i].TimeDistribOfStarFormation) && (P[i].TimeDistribOfStarFormation > 0)) // candidate for 'spawning'
            {
                double dt_remaining = P[i].TimeDistribOfStarFormation - dt_since_form_code; // time remaining to spawn
                double dt_timestep = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); // timestep being taken [code units]
                double d_tau = DMAX(DMIN( dt_timestep , dt_remaining) , 0) / P[i].TimeDistribOfStarFormation; // effective step size in dimensionless units
                double f_highz=0.0115, f_m; f_m = f_highz; // 1/mass in solar per O-star number
#ifdef FIRE_SNE_ENERGY_METAL_DEPENDENCE_EXPERIMENT // SIMPLE_POPTHREE_MODEL
                double f_lowz=2., zmin=-7, zmax=-5, z0=log10(P[i].Metallicity[0]/0.01 + 1.e-15);
                if(z0<=zmin) {f_m=f_lowz;} else if(z0<zmax) {f_m=f_lowz + (f_highz-f_lowz)*(z0-zmin)/(zmax-zmin);}
#endif
                double n_expected_total = f_m * P[i].Mass * UNIT_MASS_IN_SOLAR; // default f_m above is 1 O-star per 100 Msun [more exactly calculated here as number of stars per solar mass with mass > 8 Msun, from our adopted Kroupa IMF from 0.01-100 Msun]
                double dn_expected_in_dt = n_expected_total * d_tau; // expected number to spawn in this particular interval
                gsl_rng *random_generator_for_massivestars; // allocate rng
                random_generator_for_massivestars = gsl_rng_alloc(gsl_rng_ranlxd1); // setup generator
                gsl_rng_set(random_generator_for_massivestars, P[i].ID + 49531 + All.NumCurrentTiStep); // seed random number
                unsigned int n_to_add = gsl_ran_poisson(random_generator_for_massivestars, dn_expected_in_dt); // actually draw the number
                if(n_to_add > 0) // hey, new stars!
                {
                    P[i].IMF_WeightedMeanStellarFormationTime = (P[i].IMF_NumMassiveStars * P[i].StellarAge + n_to_add * All.Time) / (P[i].IMF_NumMassiveStars + n_to_add); // update the effective stellar age (so e.g. if we lose stars, this can keep updating)
                    P[i].IMF_NumMassiveStars += n_to_add; // add this to the number of massive stars that we are tracking explicitly
                }
            }
        }
    }
}
#endif


#ifdef METALS
void get_jet_yields(double *yields, int i) {
    int k; for(k=0;k<NUM_METAL_SPECIES;k++) {yields[k]=P[i].Metallicity[k];} /* return surface abundances, to leading order */
#ifdef SNE_NONSINK_SPAWN
    if(P[i].Type==4) {double Msne; get_SNe_yields(yields,i,stellar_lifetime_in_Gyr(i),0,&Msne);}
#endif
#ifdef STARFORGE_FEEDBACK_TRACERS
    for(k=0;k<NUM_STARFORGE_FEEDBACK_TRACERS;k++) {yields[NUM_METAL_SPECIES-NUM_STARFORGE_FEEDBACK_TRACERS+k]=0;} yields[NUM_METAL_SPECIES-NUM_STARFORGE_FEEDBACK_TRACERS+0]=1; // this is 'fully' jet material, so mark as such here, so it is noted for all wind routines [whichever form of the wind subroutine we actually use, otherwise it would only appear in the jet version]
#endif
}
#endif



#ifdef SINGLE_STAR_FB_JETS
/* Calculates the launch velocity of jets, in physical (not comoving) units */
double single_star_jet_velocity(int n)
{
    return (All.BAL_f_launch_v * sqrt(All.G * P[n].BH_Mass / (10. / UNIT_LENGTH_IN_SOLAR))); // we use the flag as a multiplier times the Kepler velocity at the protostellar radius. Really we'd want v_kick = v_kep * m_accreted / m_kicked to get the right momentum; without a better guess, assume fiducial protostellar radius of 10*Rsun, as in Federrath 2014
}
#endif










#endif /* GALSF */
