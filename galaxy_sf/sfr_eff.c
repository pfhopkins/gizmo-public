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

/*!
 *  routines for star formation in cosmological/galaxy/single-star/black hole simulations
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel,
 *   but the physical modules for star formation and feedback have been
 *   replaced, and the algorithm is mostly new to GIZMO. Many additional modules
 *   added since, with significant contributions from Mike Grudic.
 */


#ifdef GALSF // top-level switch for compiling the routines below //


#if defined(GALSF_SFR_IMF_VARIATION) || defined(GALSF_SFR_IMF_SAMPLING)
/* function to determine what the IMF of a new star particle will be, based
    on the gas properties of the particle out of which it forms */
void assign_imf_properties_from_starforming_gas(int i)
{
#ifdef GALSF_SFR_IMF_VARIATION
    double h = Get_Particle_Size(i) * All.cf_atime;
    double cs = Get_Gas_effective_soundspeed_i(i) * All.cf_afac3; // actual sound speed in the simulation: might be unphysically high for SF conditions!
    cs = 0.2 / UNIT_VEL_IN_KMS; // set to a minimum cooling temperature, for the actual star-forming conditions. for now, just use a constant //
    double dv2_abs = 0; /* calculate local velocity dispersion (including hubble-flow correction) in physical units */
    // squared norm of the trace-free symmetric [shear] component of the velocity gradient tensor //
    dv2_abs = ((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1])*(SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) +
                        (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2])*(SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) +
                        (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])*(SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
               (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] +
                         SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] +
                         SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) -
                        (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] +
                         SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] +
                         SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2]))) * All.cf_a2inv*All.cf_a2inv;
    double M_sonic = cs*cs*cs*cs / (All.G * dv2_abs * h) * UNIT_MASS_IN_SOLAR; // sonic mass in solar units //
    P[i].IMF_Mturnover = DMAX(0.01,DMIN(M_sonic,100.));
    P[i].IMF_Mturnover = 2.0; // 'normal' IMF in our definitions


    /* now we need to record all the properties we care to save about the star-forming gas, for the sake of later use: */
    int j,k;
    double NH = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i);
    double dv2abs_tot = 0; /* calculate complete velocity dispersion (including hubble-flow correction) in physical units */
    for(j=0;j<3;j++)
    {
        for(k=0;k<3;k++)
        {
            double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
            if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
            dv2abs_tot += vt*vt;
        }
    }
    double acc=0,vel=0;
    for(k=0;k<3;k++)
    {
        double acc_tmp = P[i].GravAccel[k];
#ifdef PMGRID
        acc_tmp += P[i].GravPM[k];
#endif
        acc_tmp *= All.cf_a2inv;
        acc += acc_tmp * acc_tmp;
        vel += SphP[i].VelPred[k]*SphP[i].VelPred[k];
    }
    double b_mag = 0;
#ifdef MAGNETIC
    double gizmo2gauss_2 = UNIT_B_IN_GAUSS*UNIT_B_IN_GAUSS;
    for(k=0;k<3;k++) {b_mag += Get_Gas_BField(i,k)*Get_Gas_BField(i,k) * gizmo2gauss_2;}
#endif
    double rad_flux_uv = 1;
    double cr_energy_density = 0;
#ifdef COSMIC_RAY_FLUID
    int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {cr_energy_density += SphP[i].CosmicRayEnergyPred[k_CRegy] * SphP[i].Density * All.cf_a3inv / P[i].Mass;}
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
    P[i].IMF_FormProps[0] = P[i].min_dist_to_bh; // min distance to nearest sink particle
#else
    P[i].IMF_FormProps[0] = P[i].IMF_Mturnover; // IMF turnover mass as defined above
#endif
    P[i].IMF_FormProps[1] = SphP[i].Density * All.cf_a3inv; // density
    P[i].IMF_FormProps[2] = SphP[i].InternalEnergyPred; // thermal internal energy (use to calculate temperature)
    P[i].IMF_FormProps[3] = Get_Gas_effective_soundspeed_i(i) * All.cf_afac3; // sound speed (not trivially related to temperature if CRs, etc included)
    P[i].IMF_FormProps[4] = sqrt(dv2_abs); // shear velocity gradient (norm of shear gradient tensor)
    P[i].IMF_FormProps[5] = h; // particle length/size (inter-particle spacing)
    P[i].IMF_FormProps[6] = NH; // local gas surface density (our usual estimator) in the cloud where the particle formed
    P[i].IMF_FormProps[7] = sqrt(dv2abs_tot) * h; // total rms/turbulent velocity dispersion
    P[i].IMF_FormProps[8] = sqrt(acc); // gravitational acceleration
    P[i].IMF_FormProps[9] = sqrt(vel); // total velocity (use with acceleration to estimate shear omega, etc)
    P[i].IMF_FormProps[10] = sqrt(b_mag) * All.cf_a2inv; // magnetic field strength |B|
    P[i].IMF_FormProps[11] = rad_flux_uv; // incident UV flux normalized to MW 'canonical' (Habing) field value
    P[i].IMF_FormProps[12] = cr_energy_density; // cosmic ray energy density (if CRs are enabled)

#endif

#ifdef GALSF_SFR_IMF_SAMPLING
    gsl_rng *random_generator_for_massivestars;
    random_generator_for_massivestars = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(random_generator_for_massivestars, P[i].ID+121);
    double mu = 0.01 * P[i].Mass * UNIT_MASS_IN_SOLAR; // 1 O-star per 100 Msun
    unsigned int kk = gsl_ran_poisson(random_generator_for_massivestars, mu);
    P[i].IMF_NumMassiveStars = (double)kk;
#endif
}
#endif



/* return the stellar age in Gyr for a given labeled age, needed throughout for stellar feedback */
double evaluate_stellar_age_Gyr(double stellar_tform)
{
    double age,a0,a1,a2,x0,x1,x2;
    if(All.ComovingIntegrationOn)
    {
        a0 = stellar_tform;
        a2 = All.Time;
        if(fabs(1-(All.OmegaMatter+All.OmegaLambda))<=0.01)
        {
            /* use exact solution for flat universe, ignoring the radiation-dominated epoch [no stars forming then] */
            x0 = (All.OmegaMatter/(1-All.OmegaMatter))/(a0*a0*a0);
            x2 = (All.OmegaMatter/(1-All.OmegaMatter))/(a2*a2*a2);
            age = (2./(3.*sqrt(1-All.OmegaMatter)))*log(sqrt(x0*x2)/((sqrt(1+x2)-1)*(sqrt(1+x0)+1)));
            age *= 1./All.Hubble_H0_CodeUnits;
        } else {
            /* use simple trap rule integration */
            a1 = 0.5*(a0+a2);
            x0 = 1./(a0*hubble_function(a0));
            x1 = 1./(a1*hubble_function(a1));
            x2 = 1./(a2*hubble_function(a2));
            age = (a2-a0)*(x0+4.*x1+x2)/6.;
        }
    } else {
        /* time variable is simple time, when not in comoving coordinates */
        age=All.Time-stellar_tform;
    }
    age *= UNIT_TIME_IN_GYR; // convert to absolute Gyr
    if((age<=1.e-5)||(isnan(age))) {age=1.e-5;}
    return age;
}



/* simple routine to determine density thresholds and other common units for SF routines */
void set_units_sfr(void)
{
    All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
    All.PhysDensThresh = All.CritPhysDensity / UNIT_DENSITY_IN_NHCGS;
#ifdef GALSF_EFFECTIVE_EQS
    All.EgySpecCold = All.TempClouds / ((4 / (1 + 3 * HYDROGEN_MASSFRAC)) * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS); /* note: assuming fully-neutral atomic H+He primordial mixture */
    All.EgySpecSN = All.TempSupernova / ((4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS); /* note: assuming fully-ionized H+He primordial mixture */
#endif
}



/* function which takes properties of a gas particle 'i' and returns probability of its turning into a BH seed particle */
double return_probability_of_this_forming_bh_from_seed_model(int i)
{
    double p=0;
#ifdef BH_SEED_FROM_LOCALGAS
    double Z_threshold_solar = 0.01, surfacedensity_threshold_cgs = 1.0; /* metallicity below which, and density above which, seed BH formation is efficient */
    /* note surface dens in g/cm^2; threshold for bound cluster formation in our experiments is ~0.2-2 g/cm^2 (10^3 - 10^4 M_sun/pc^2) */
    if(All.ComovingIntegrationOn) {if(All.Time > 1/(1+All.SeedBlackHoleMinRedshift)) {return 0;}} /* outside allowed redshift */
    if(SphP[i].Density*All.cf_a3inv < All.PhysDensThresh) {return 0;} /* must be above SF density threshold */
    double Z_in_solar = P[i].Metallicity[0]/All.SolarAbundances[0], surfacedensity = MIN_REAL_NUMBER;
    /* now calculate probability of forming a BH seed particle */
    p = P[i].Mass / All.SeedBlackHolePerUnitMass; /* probability of forming a seed per unit mass [in code units] */
#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    double Rcrit = PPP[i].Hsml;
    Z_threshold_solar = 0.1; /* based on Linhao's paper, we need to allow formation at somewhat higher metallicity or we tail to get BHs in the central density concentrations when they form */
#if !defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORALL)
    Rcrit = All.ForceSoftening[0]; /* search radius is not h, in this case, but the force softening, but this is really not the case we want to study */
#endif
    Rcrit = DMAX( Rcrit , 0.1/(UNIT_LENGTH_IN_KPC*All.cf_atime)); /* set a baseline Rcrit_min, otherwise we get statistics that are very noisy */
#ifdef BH_CALC_DISTANCES
    if(P[i].min_dist_to_bh < 10.*Rcrit) {return 0;} /* don't allow formation if there is already a sink nearby, akin to SF sink rules */
#endif
    surfacedensity = P[i].MencInRcrit / (M_PI*Rcrit*Rcrit) * UNIT_SURFDEN_IN_CGS * All.cf_a2inv; /* this is the -total- mass density inside the critical kernel radius Rcrit, evaluated within the tree walk */
    double Z_u = Z_in_solar/Z_threshold_solar, S_u = surfacedensity / surfacedensity_threshold_cgs;
    if(!isfinite(Z_u) || !isfinite(S_u) || (S_u < 0.01)) {return 0;}
    if(S_u < 3.5) {p *= 1 - exp(-S_u*S_u);} // quadratic cutoff at low densities: probability drops as S^(2), saturates at 1
    p /= 1 + Z_u + 0.5*Z_u*Z_u; // quadratic expansion of exponential cutoff: probability drops as Z^(-2) rather than exp(-Z), saturates at 1
#else
    surfacedensity = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i) * UNIT_SURFDEN_IN_CGS; /* this gives the Sobolev-estimated column density of -gas- alone */
    if(surfacedensity>0.1*surfacedensity_threshold_cgs) {p *= (1-exp(-surfacedensity/surfacedensity_threshold_cgs)) * exp(-Z_in_solar/Z_threshold_solar);} else {p=0;} /* apply threshold metallicity and density cutoff */
#endif
#endif
    if(p > 12.) {p=1;} else {if(p > 1.e-4) {p=1-exp(-p);}}
    return p;
}



/* Routine to actually determine the SFR assigned to an individual gas particle at each time. i is the index. mode is normal [0] or for snapshot output [1], where the latter is included to make sure certain flags dont give misleading outputs */
double get_starformation_rate(int i, int mode)
{
    double rateOfSF,tsfr,y; y=0; int flag=1, j, k; /* flag to proceed to SFR calc */
    if(P[i].Mass <= 0 || SphP[i].Density <= 0) {flag=0;} /* zero-mass elements [for deletion] not eligible for SF */
#ifdef GALSF_SUBGRID_WINDS
    if(SphP[i].DelayTime > 0) {flag=0;} /* 'decoupled' wind elements not eligible for SF */
#endif
#ifdef BH_WIND_SPAWN
    if(P[i].ID == All.AGNWindID) {flag=0;} /* spawned hyper-resolution elements not eligible for SF */
#endif
    if(All.ComovingIntegrationOn && SphP[i].Density < All.OverDensThresh) {flag=0;} /* below overdensity threshold required for SF */
    if(SphP[i].Density*All.cf_a3inv < All.PhysDensThresh) {flag=0;} /* below physical density threshold */
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED) && !defined(GALSF_SFR_VIRIAL_CONTINUOUS)
    if(flag==0) {SphP[i].AlphaVirial_SF_TimeSmoothed=0;} /* for time-smoothed virial param, reset to nil if fall below threshold */
#endif
    tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * All.cf_a3inv)) * All.MaxSfrTimescale; /* set default SFR timescale to scale appropriately with the gas dynamical time */
    rateOfSF = P[i].Mass / tsfr; /* 'normal' sfr from density law above */
    if(tsfr<=0 || rateOfSF<=0 || flag==0) {return 0;} /* nonsense here, return 0 */
    
#if defined(SINGLE_STAR_SINK_DYNAMICS) && defined(SINGLE_STAR_SINK_FORMATION)
    int cell_can_be_singlestar = is_particle_single_star_eligible(i); // call function to determine if we're actually eligible to be a true single-star element
#endif

#ifdef GALSF_EFFECTIVE_EQS /* do the SFR calc for the Springel-Hernquist EOS, before any 'fancy' sf criteria, when above-threshold, or else risk incorrect entropies */
    double factorEVP = pow(SphP[i].Density * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP; /* evaporation factor */
    double egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold; /* specific energy of hot [volume-filling] phase gas */
    double tcool = GetCoolingTime(egyhot, SphP[i].Density * All.cf_a3inv, SphP[i].Ne, i); /* cooling time of two-phase mix */
    y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold); /* parameter */
    double cloudmass_fraction = (1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y))); /* quasi-equilibrium mass in cold phase */
    rateOfSF = (1 - All.FactorSN) * cloudmass_fraction * P[i].Mass / tsfr; /* SFR given by cold mass (less SNe-entrainment fraction) divided by tSFR */
    update_internalenergy_for_galsf_effective_eos(i,tcool,tsfr,cloudmass_fraction,rateOfSF); /* updates entropies for the effective equation-of-state */
#endif

    int exceeds_force_softening_threshold; exceeds_force_softening_threshold = 0; /* flag that notes if the density is so high such that gravity is non-Keplerian [inside of smallest force-softening limits] */
#if (SINGLE_STAR_SINK_FORMATION & 1024)
    if(DMIN(PPP[i].Hsml, 2.*Get_Particle_Size(i)) <= DMAX(All.MinHsml, 2.*All.ForceSoftening[0])) {exceeds_force_softening_threshold=1;}
    if(mode == 0) {if(exceeds_force_softening_threshold) {return 10. * rateOfSF;}}
#endif

    /* compute various velocity-gradient terms which are potentially used in the various criteria below */
    double dv2abs=0, dv2abs_0=0, divv=0, gradv[9]={0}, cs_eff=0, vA=0, v_fast=0; /* calculate local velocity dispersion (including hubble-flow correction) in physical units */
    cs_eff=Get_Gas_thermal_soundspeed_i(i); vA=Get_Gas_Alfven_speed_i(i); /* specifically get the -isothermal- soundspeed and Alfven speed  (since we're doing a local Jeans analysis using these terms) [dont include terms like radiation pressure or cosmic ray pressure in the relevant speeds here] */
    v_fast=sqrt(cs_eff*cs_eff + vA*vA); /* calculate fast magnetosonic speed for use below */
    for(j=0;j<3;j++) {
        for(k=0;k<3;k++) {
            double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
            if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
            gradv[3*j + k]=vt; dv2abs+=vt*vt; if(j==k) {divv+=vt;} // save for possible use below
        }}
    dv2abs_0=dv2abs; // save for possible use below
#if defined(SINGLE_STAR_SINK_DYNAMICS) && defined(SINGLE_STAR_SINK_FORMATION) && ((defined(COOLING) && !defined(COOL_LOWTEMP_THIN_ONLY)) || defined(EOS_GMC_BAROTROPIC)) // if we have to deal with optically-thick thermo
    double nHcgs = HYDROGEN_MASSFRAC * (SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS);
    if(nHcgs > 1e13) {v_fast = DMIN(v_fast, 0.2/UNIT_VEL_IN_KMS);} // limiter to permit sink formation in simulations that really resolve the opacity limit and bog down when an optically-thick core forms. Modify this if you want to follow first collapse more/less - scale as c_s ~ n^(1/5)
#endif

#if defined(GALSF_SFR_VIRIAL_SF_CRITERION) /* apply standard virial-parameter criteria here. note that our definition of virial parameter here is ratio of [Kinetic+Internal Energy]/|Gravitational Energy| -- so <1=bound, >1=unbound, 1/2=bound-and-virialized, etc. */
    double k_cs = 1. * v_fast / (Get_Particle_Size(i)*All.cf_atime), alpha_crit; alpha_crit = 1.0; /* effective wavenumber for thermal+B-field+CR+whatever internal energy support, and threshold virial parameter */
#if defined(SINGLE_STAR_SINK_DYNAMICS)
    if(cell_can_be_singlestar) {k_cs *= M_PI;} // use a stricter version here, because the relevant pre-factor depends on whether we expect Jeans collapse at the thermal limit to be resolved or un-resolved
#endif
#if (GALSF_SFR_VIRIAL_SF_CRITERION > 0)
    if(divv < 0) {dv2abs -= divv*divv/3.;} // this is mathematically identical to taking the Frobenius norm of the divergence-free (trace-free) part of the shear tensor instead of the whole tensor. when using the stricter criterion, if the gas is in inflow (divv<0), don't want to count the inflow motion itself against the virial criterion, since this could if close to free-fall artificially bias us against recognizing real star formation
#endif
    double Mach_eff_2=0, cs2_contrib=2.*k_cs*k_cs; Mach_eff_2=dv2abs/cs2_contrib; dv2abs+=2.*k_cs*k_cs; // account for thermal+magnetic pressure with standard Jeans criterion (k^2*cs^2 vs 4pi*G*rho) //
    double alpha_vir = dv2abs / (8.*M_PI * All.G * SphP[i].Density * All.cf_a3inv); // coefficient comes from different density profiles, assuming a constant velocity gradient tensor: 22.6=constant-density cube, 8pi[approximate]=constant-density sphere, e.g. rho~exp(-r^n) n={4,8,16,32,64}->{17.1,22.1,24.1,24.9,25.1,25.15} [approaches uniform-density sphere as n->infinity]
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED) /* compute and prepare to use our time-rolling average virial criterion */
    double dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); /* the physical time-step */
    double alpha_0=1./(1.+alpha_vir), dtau=DMIN(1.,DMAX(0.,exp(-(DMIN(DMAX(8.*dtime/tsfr,0.),20.))))); /* dimensionless units for below */
    SphP[i].AlphaVirial_SF_TimeSmoothed = DMIN(DMAX(SphP[i].AlphaVirial_SF_TimeSmoothed * dtau + alpha_0 * (1.-dtau) , 1.e-10), 1.); /* update rolling time-averaged virial parameter */
    alpha_vir = 1./SphP[i].AlphaVirial_SF_TimeSmoothed - 1.; /* use the rolling average below */
#endif
    if(exceeds_force_softening_threshold) {alpha_vir /= 10.;} /* account for gravitational softening effects here, making this threshold less steep */
#if (GALSF_SFR_VIRIAL_SF_CRITERION <= 0) && !defined(GALSF_SFR_VIRIAL_CONTINUOUS) && !(SINGLE_STAR_SINK_FORMATION & 512) /* 'weakest' mode: reduce [do not zero] SFR if above alpha_crit, and not -too- dense */
    if((alpha_vir>alpha_crit) && (SphP[i].Density*All.cf_a3inv<100.*All.PhysDensThresh)) {rateOfSF *= 0.0015;} /* PFH: note the 100x threshold limit here is an arbitrary choice currently set -by hand- to prevent runaway densities from this prescription! */
#endif
#if (GALSF_SFR_VIRIAL_SF_CRITERION > 0) && !defined(GALSF_SFR_VIRIAL_CONTINUOUS) && !(SINGLE_STAR_SINK_FORMATION & 512) /* 'normal' mode: zero SF if don't meet virial threshold */
    if(alpha_vir>alpha_crit) {rateOfSF=0;} /* simple 'hard' threshold here */
#endif
#if (SINGLE_STAR_SINK_FORMATION & 512) || defined(GALSF_SFR_VIRIAL_CONTINUOUS) /* semi-continuous SF as a function of alpha_vir */
    double sf_eff_multiplier = 1;
    if(alpha_vir>alpha_crit) {sf_eff_multiplier = 0.01;}
#if (GALSF_SFR_VIRIAL_CONTINUOUS == 1)
    sf_eff_multiplier = exp(-1.4 * DMIN(DMIN(DMAX(sqrt(DMAX(MIN_REAL_NUMBER,alpha_vir)), 1.e-4), 1.e10),22.)); /* continuous cutoff of rateOfSF with increasing virial parameter as ~exp[-1.4*sqrt(alpha_vir)], following fitting function from Padoan 2012 [limit the values of sqrt(alpha_vir) here since we'll take an exponential so don't want a nan] */
#endif
#if (GALSF_SFR_VIRIAL_CONTINUOUS == 2) || (SINGLE_STAR_SINK_FORMATION & 512)
    Mach_eff_2 = DMIN(DMAX(1.e-5, Mach_eff_2/3.), 1.e4); if(!isfinite(Mach_eff_2)) {Mach_eff_2=1.e4;}
    double S_ln=log(1.+Mach_eff_2/4.), S_crit=log(alpha_vir*(1.+2.*Mach_eff_2*Mach_eff_2/(1.+Mach_eff_2*Mach_eff_2))); // Mach_eff_2 is determined by the ratio of the kinetic to the thermal terms in the virial parameter, corrected to the 1D dispersion here
    sf_eff_multiplier *= 0.5 * exp(3.*S_ln/8.) * (1. + erf((S_ln-S_crit)/sqrt(2.*S_ln))); // multi-free-fall model, as in e.g. Federrath+Klessen 2012/2013 ApJ 761,156; 763,51 (similar to that implemented in e.g. Kretschmer+Teyssier 2020), based on the analytic models in Hopkins MNRAS 2013, 430 1653, with correct virial parameter [K+T used a definition which gives the wrong value for thermally-supported clouds]
#endif
    rateOfSF *= sf_eff_multiplier;
#endif
#endif

#if (SINGLE_STAR_SINK_FORMATION & 256) || defined(GALSF_SFR_MOLECULAR_CRITERION) /* scale SFR to fraction of 'molecular' gas in cell */
    double ne=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred; // pull various known thermal properties, prepare to extract others //
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties, like neutral fraction, temperature, etc, that we will use below //
    rateOfSF *= Get_Gas_Molecular_Mass_Fraction(i, temperature, nh0, ne, 0.);
#endif

#if (SINGLE_STAR_SINK_FORMATION & 2) || (GALSF_SFR_VIRIAL_SF_CRITERION >= 4) /* restrict to convergent flows */
#if !defined(GALSF_SFR_VIRIAL_SF_CRITERION) || defined(SINGLE_STAR_SINK_DYNAMICS)
    if(divv >= 0) {rateOfSF=0;} /* diverging flow, no SF [simplest version of this] */
#else
    if(divv < 0) { // here need to more carefully account for the fact that we are not trying to capture laminar collapse, but whether a patch should or should not -fragment-
        /* double t_compression = -1./divv; if(t_compression < tsfr) {rateOfSF *= tsfr/t_compression;} // if you are collapsing or being externally compressed faster than the local dynamical time, allow SF to characteristically follow -that- timescale, to prevent runaway densities from appearing spuriously (will have no effect in single star runs since sf is deterministic and this normalization is set to infinity at the end) */
    } else {
        double dv_turb_est = sqrt(DMAX(dv2abs_0 - divv*divv/3., MIN_REAL_NUMBER)); // magnitude of Frobenius norm of the trace-free shear tensor. if we had perfectly isotropic turbulence, the rms value of this would be ~sqrt(8)*q, where q = <|dvi/dxj|^2>^(1/2). meanwhile the rms value of the divergence would be ~sqrt[3]*x, or ~sqrt[3/8]~0.6 times this value. so this gives us a reasonable guess of when divv is smaller than we would expect of rms fluctuations in an isotropic turbulent random field.
        if(P[i].Particle_DivVel>0 && 1./divv<tsfr && divv > 0.3*dv_turb_est) {rateOfSF=0;} // take a threshold which is ~(1/2) of the rms, corresponding to divergence being ~1% (0.2) or 3% (0.3) of the contribution to the Frobenius norm above: very low threshold to consider this, could easily raise it. also check against noise by confirming that the P[i].Particle_DivVel variable has the same sign. And require that the characteristic timescale for the velocity divergence to influence the density field is actually shorter than the gravitational collapse timescale.
    }
#endif
#endif

#if (SINGLE_STAR_SINK_FORMATION & 128) || (GALSF_SFR_VIRIAL_SF_CRITERION >= 4) /* check that the velocity gradient is negative-definite, ie. converging along all principal axes, which is much stricter than div v < 0 */
    for(j=0;j<3;j++){ // symmetrize the velocity gradient
      for(k=0;k<j;k++){double temp = gradv[3*j + k]; gradv[3*j + k] = 0.5*(gradv[3*j + k] + gradv[3*k + j]); gradv[3*k + j] = 0.5*(temp + gradv[3*k + j]);}}
    gsl_matrix_view M = gsl_matrix_view_array(gradv, 3, 3); gsl_vector *eval1 = gsl_vector_alloc(3);
    gsl_eigen_symm_workspace *v = gsl_eigen_symm_alloc(3); gsl_eigen_symm(&M.matrix, eval1,  v);
    if(exceeds_force_softening_threshold==0) {for(k=0;k<3;k++) if(gsl_vector_get(eval1,k) >= 0) {rateOfSF=0;}} /* cannot apply this criterion when we exceed the limits where gravity is treated as fully-Newtonian, it will severely suppress 'true' collapse */
    gsl_eigen_symm_free(v); gsl_vector_free(eval1);
#endif

#if (SINGLE_STAR_SINK_FORMATION & 64) || (GALSF_SFR_VIRIAL_SF_CRITERION >= 3) /* check if Jeans mass is low enough for conceivable formation of 'stars' */
    double cs_touse=cs_eff, MJ_crit=DMAX(DMIN(1.e6, 1.*P[i].Mass*UNIT_MASS_IN_SOLAR), 100.); /* for galaxy-scale SF, default to large ~1000 Msun threshold */
    if(exceeds_force_softening_threshold) {MJ_crit = DMAX(1.e6 , 100.*P[i].Mass*UNIT_MASS_IN_SOLAR);}
#if defined(SINGLE_STAR_SINK_DYNAMICS) && defined(SINGLE_STAR_SINK_FORMATION)
    if(cell_can_be_singlestar) {cs_touse=v_fast; MJ_crit=DMIN(1.e4, DMAX(1.e-3 , 100.*P[i].Mass*UNIT_MASS_IN_SOLAR));} /* for single-star formation use un-resolved Jeans mass criterion, with B+thermal pressure */
#endif
    double MJ_solar = 2.*pow(cs_touse*UNIT_VEL_IN_KMS/0.2,3)/sqrt(SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS / (HYDROGEN_MASSFRAC*1.0e3));
    if(MJ_solar > MJ_crit) {rateOfSF=0;} /* if too massive Jeans mass, go no further */
#endif

#ifdef SINGLE_STAR_AND_SSP_HYBRID_MODEL
    if(cell_can_be_singlestar==0) {return rateOfSF;} // special call for our hybrid model - it's important that "sink-only" SF criteria fall below this line, otherwise wont be able to form stars in the more diffuse medium as you should
#endif
    
#ifdef GALSF_SFR_TIDAL_HILL_CRITERION /* check that the tidal tensor is negative-definite, ie. converging along all principal axes, indicating that we're dominating our environment gravitationally and are living in our own Hill sphere */
    if(exceeds_force_softening_threshold==0) {for(k=0;k<3;k++) {if(P[i].tidal_tensorps[k][k] >= 0) {rateOfSF=0;}}} /* we've already diagonized this in gravtree.c, so this is a straightforward check. again should only be applied where force calculation is fully-reliable */
#endif

#if (SINGLE_STAR_SINK_FORMATION & 4) /* restrict to local density/potential maxima */
    if(SphP[i].Density_Relative_Maximum_in_Kernel > 0) {rateOfSF=0;}
#endif

#if (SINGLE_STAR_SINK_FORMATION & 8) /* restrict to cell which neither 'sees' or 'is seen by' a sink too close */
    if(P[i].BH_Ngb_Flag) {rateOfSF=0;} /* cell cannot be 'seen' by -any- sink as a potential interacting neighbor */
    if(P[i].min_dist_to_bh < P[i].Hsml){rateOfSF=0;} /* cell does not overlap with a sink */
#endif

#if (SINGLE_STAR_SINK_FORMATION & 16) /* restrict to cells which have a local SF time or free-fall time shorter than their free-fall time onto the nearest sink */
    if(DMIN(P[i].min_bh_approach_time, P[i].min_bh_freefall_time) < tsfr) {rateOfSF=0;}
#endif



#if defined(SINGLE_STAR_SINK_DYNAMICS) && defined(SINGLE_STAR_SINK_FORMATION)
    rateOfSF *= 1.0e20; /* make sink formation guaranteed to happen, where it can, by setting rate super-high if non-zero */
#endif
    return rateOfSF; /* finally, we have a SFR! */
}



#ifdef GALSF_EFFECTIVE_EQS
/* compute the 'effective eos' cooling/heating, including thermal feedback sources, here */
void update_internalenergy_for_galsf_effective_eos(int i, double tcool, double tsfr, double cloudmass_fraction, double rateOfSF)
{
    double dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); /*  the actual time-step */
    double x = cloudmass_fraction, factorEVP = pow(SphP[i].Density * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP, trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
    double egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold, egyeff = egyhot * (1 - x) + All.EgySpecCold * x, egycurrent = SphP[i].InternalEnergy, ne;
    ne=1.0;

#if defined(BH_THERMALFEEDBACK)
    if((SphP[i].Injected_BH_Energy > 0) && (P[i].Mass>0))
    {
        egycurrent += SphP[i].Injected_BH_Energy / P[i].Mass;
        if(egycurrent > egyeff)
        {
            tcool = GetCoolingTime(egycurrent, SphP[i].Density * All.cf_a3inv, ne, i);
            if(tcool < trelax && tcool > 0) trelax = tcool;
        }
        SphP[i].Injected_BH_Energy = 0;
    }
#endif

    /* now update the thermal variables */
    SphP[i].InternalEnergy = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));
    SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
    SphP[i].Pressure = get_pressure(i);
    //SphP[i].dInternalEnergy = 0;
    SphP[i].DtInternalEnergy = 0; /* HERE, it's ok, b/c effective EOS is designed to model new pressure even under compressions,
                                 (since we're zero'ing the second-half-step from the hydro step) */
}
#endif // GALSF_EFFECTIVE_EQS //




/* parent routine for star formation. for 'effective equation of state' models for star-forming gas, this also updates their effective EOS parameters */
void star_formation_parent_routine(void)
{
    int i, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
    unsigned int bits; double dtime, mass_of_star, p, prob, rate_in_msunperyear, sfrrate, totsfrrate, sum_sm, total_sm, sm=0, rate, sum_mass_stars, total_sum_mass_stars;
#if defined(BH_SEED_FROM_LOCALGAS) || defined(SINGLE_STAR_SINK_DYNAMICS)
    int num_bhformed=0, tot_bhformed=0;
#endif
    for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinActive[bin]) {TimeBinSfr[bin] = 0;}}
    stars_spawned = stars_converted = 0; sum_sm = sum_mass_stars = 0;
    for(bits = 0; GALSF_GENERATIONS > (1 << bits); bits++);

    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if((P[i].Type == 0)&&(P[i].Mass>0))
        {
            SphP[i].Sfr = 0; flag = 1; /* will be reset below if flag==0, but default to flag = 1 (non-eligible) */
            dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); /*  the actual time-step */
            
            /* check whether an initial (not fully-complete!) conditions for star formation are fulfilled for a given particle */
            if(SphP[i].Density * All.cf_a3inv >= All.PhysDensThresh) {flag = 0;} // if sufficiently dense, go forward into SF routine //
            if(All.ComovingIntegrationOn) {if(SphP[i].Density < All.OverDensThresh) {flag = 1;}} // (additional density check for cosmological runs) //
            
#ifdef GALSF_SUBGRID_WINDS
            if(SphP[i].DelayTime > 0) {flag=1; SphP[i].DelayTime -= dtime;} /* no star formation for particles in the wind; update our wind delay-time calculations */
            if((SphP[i].DelayTime < 0) || (SphP[i].Density*All.cf_a3inv < All.WindFreeTravelDensFac*All.PhysDensThresh)) {SphP[i].DelayTime=0;}
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
            if(SphP[i].DelayTimeCoolingSNe > 0) {flag=1; SphP[i].DelayTimeCoolingSNe -= dtime;} /* no star formation for particles in the wind; update our wind delay-time calculations */
#endif
            
            if((flag == 0)&&(dtime>0))        /* active star formation (upon start-up, we need to protect against dt==0) */
            {
                /* Alright, now we consider the actual gas-to-star particle conversion and associated steps */
                sm = get_starformation_rate(i, 0) * dtime; // expected stellar mass formed this timestep (this also updates entropies for the effective equation-of-state model) //
                p = sm / P[i].Mass;
                double pfac = 1.-exp(-p);
                if(p<0.1) {pfac=p*(1.-0.5*p*(1.-p/3.));} // avoids this accidentally setting to 0 if p is small, from floating-point errors
                sum_sm += P[i].Mass * pfac;
                
                /* the upper bits of the gas particle ID store how many stars this gas particle gas already generated */
                if(bits == 0) {number_of_stars_generated = 0;} else {number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));}
                
                mass_of_star = P[i].Mass / (GALSF_GENERATIONS - number_of_stars_generated);
                if(number_of_stars_generated >= GALSF_GENERATIONS-1) mass_of_star=P[i].Mass;
                
                SphP[i].Sfr = sm / dtime * UNIT_MASS_IN_SOLAR / UNIT_TIME_IN_YR;
                if(dtime>0) {TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;}
                
                prob = P[i].Mass / mass_of_star * pfac;
                
#if defined(METALS) && defined(GALSF_EFFECTIVE_EQS) // does instantaneous enrichment //
                double w = get_random_number(P[i].ID);
                P[i].Metallicity[0] += w * All.SolarAbundances[0] * pfac;
                if(NUM_METAL_SPECIES>=10) {int k; for(k=1;k<NUM_METAL_SPECIES;k++) {P[i].Metallicity[k] += w * All.SolarAbundances[k] * pfac;}}
#endif
                
                if(get_random_number(P[i].ID + 1) < prob)    /* ok, make a star */
                {
                    
#ifdef BH_SEED_FROM_LOCALGAS
                    /* before making a star, assess whether or not we can instead make a BH seed particle */
                    p = return_probability_of_this_forming_bh_from_seed_model(i);
                    if(get_random_number(P[i].ID + 2) < p)
                    {
                        /* make a BH particle */
                        P[i].Type = 5;
                        TimeBinCountSph[P[i].TimeBin]--;
                        num_bhformed++;
                        Stars_converted++;
                        stars_converted++;
                        P[i].StellarAge = All.Time;
                        
                        P[i].BH_Mass = All.SeedBlackHoleMass;
                        if(All.SeedBlackHoleMassSigma > 0)
                        {
                            gsl_rng *random_generator_forbh; /* generate gaussian random number for random BH seed mass */
                            random_generator_forbh = gsl_rng_alloc(gsl_rng_ranlxd1); gsl_rng_set(random_generator_forbh, P[i].ID+17);
                            P[i].BH_Mass = pow( 10., log10(All.SeedBlackHoleMass) + gsl_ran_gaussian(random_generator_forbh, All.SeedBlackHoleMassSigma) );
                        }
                        
                        if(p>1) P[i].BH_Mass *= p; /* assume multiple seeds in particle merge */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        P[i].Mass = SphP[i].MassTrue + SphP[i].dMass;
#endif
#ifdef BH_INCREASE_DYNAMIC_MASS
                        P[i].Mass *= BH_INCREASE_DYNAMIC_MASS;
#endif
#ifdef BH_ALPHADISK_ACCRETION
                        P[i].BH_Mass_AlphaDisk = All.SeedAlphaDiskMass;
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                        double bh_mu=2.0*get_random_number(P[i].ID+3)-1.0, bh_phi=2*M_PI*get_random_number(P[i].ID+4), bh_sin=sqrt(1-bh_mu*bh_mu);
                        double spin_prefac = All.G * P[i].BH_Mass / C_LIGHT_CODE; // assume initially maximally-spinning BH with random orientation
                        P[i].BH_Specific_AngMom[0]=spin_prefac * bh_sin*cos(bh_phi); P[i].BH_Specific_AngMom[1]=spin_prefac * bh_sin*sin(bh_phi); P[i].BH_Specific_AngMom[2]=spin_prefac * bh_mu;
#endif
#ifdef BH_WIND_SPAWN
                        P[i].unspawned_wind_mass = 0;
#endif
#ifdef BH_COUNTPROGS
                        P[i].BH_CountProgs = 1;
#endif
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                        P[i].SinkRadius = All.ForceSoftening[P[i].Type];
#endif
                        P[i].BH_Mdot = 0;
                        P[i].DensAroundStar = SphP[i].Density;
                    } else {
#endif /* closes ifdef(BH_SEED_FROM_LOCALGAS) */
                        
                        /* ok, we're going to make a star! */
#if defined(GALSF_SFR_IMF_VARIATION) || defined(GALSF_SFR_IMF_SAMPLING)
                        /* if we're allowing for a variable IMF, this is where we will calculate the IMF properties produced from the gas forming stars */
                        assign_imf_properties_from_starforming_gas(i);
#endif
                        
                        if(number_of_stars_generated == (GALSF_GENERATIONS - 1))
                        {
                            /* here we turn the gas particle itself into a star */
                            Stars_converted++;
                            stars_converted++;
                            sum_mass_stars += P[i].Mass;
                            
                            TimeBinCountSph[P[i].TimeBin]--;
                            TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;
                            
                            P[i].StellarAge = All.Time;
                            
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
                            P[i].DensAroundStar = SphP[i].Density;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                            P[i].Mass = SphP[i].MassTrue + SphP[i].dMass;
#endif
                            
#ifdef SINGLE_STAR_SINK_DYNAMICS
                            if(is_particle_single_star_eligible(i))
                            {
                                P[i].Type = 5;
                                num_bhformed++;
                                P[i].BH_Mass = All.SeedBlackHoleMass; // if desired to make this appreciable fraction of particle mass, please do so in params file
#ifdef HERMITE_INTEGRATION
                                P[i].AccretedThisTimestep = 0;
#endif
#ifdef GRAIN_FLUID
                                P[i].BH_Dust_Mass = 0;
#endif
#ifdef BH_RETURN_BFLUX
                                P[i].B[0] = P[i].B[1] = P[i].B[2] = 0;
#endif
                                TreeReconstructFlag = 1;
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                                P[i].SinkRadius = All.ForceSoftening[5];
                                double cs = 0.2 / UNIT_VEL_IN_KMS;
#if (defined(COOLING) && !defined(COOL_LOWTEMP_THIN_ONLY)) || defined(EOS_GMC_BAROTROPIC)
                                double nHcgs = HYDROGEN_MASSFRAC * (SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_NHCGS);
                                if(nHcgs > 1e10) cs *= pow(nHcgs/1e10, 1./5); // if we're getting opacity-limited then we can set a smaller sink radius, since cs ~ n^1/5
#endif
                                P[i].SinkRadius = DMAX(0.79 * P[i].Mass * All.G / (cs * cs), All.ForceSoftening[5]); // volume-equivalent particle radius R= (3V/(4PI))^(1/3) at the density where cell length = Jeans length/2
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
                                P[i].min_bh_t_orbital=MAX_REAL_NUMBER; P[i].comp_dx[0]=P[i].comp_dx[1]=P[i].comp_dx[2]=P[i].comp_dv[0]=P[i].comp_dv[1]=P[i].comp_dv[2]=P[i].is_in_a_binary = 0;
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
                                P[i].SuperTimestepFlag=P[i].COM_GravAccel[0]=P[i].COM_GravAccel[1]=P[i].COM_GravAccel[2]=P[i].comp_Mass=P[i].COM_dt_tidal=0;
#endif
#ifdef BH_ALPHADISK_ACCRETION
                                P[i].BH_Mass_AlphaDisk = DMAX(DMAX(0, P[i].Mass-P[i].BH_Mass), All.SeedAlphaDiskMass);
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                                double bh_mu=2.0*get_random_number(P[i].ID+3)-1.0, bh_phi=2*M_PI*get_random_number(P[i].ID+4), bh_sin=sqrt(1-bh_mu*bh_mu);
                                double spin_prefac = All.G * P[i].BH_Mass / C_LIGHT_CODE; // assume initially maximally-spinning BH with random orientation
                                P[i].BH_Specific_AngMom[0]=spin_prefac*bh_sin*cos(bh_phi); P[i].BH_Specific_AngMom[1]= spin_prefac * bh_sin*sin(bh_phi); P[i].BH_Specific_AngMom[2]=spin_prefac * bh_mu;
#endif
#ifdef BH_COUNTPROGS
                                P[i].BH_CountProgs = 1;
#endif
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
                                P[i].dt_since_last_gas_search = 0;
                                P[i].do_gas_search_this_timestep = 1;
                                P[i].BH_TimeBinGasNeighbor = P[i].TimeBin;
#endif
                                P[i].BH_Mdot = 0;
                                P[i].DensAroundStar = SphP[i].Density;
                                
#ifdef OUTPUT_SINK_FORMATION_PROPS //save the at-formation properties of sink particles
                                MyDouble tempB[3]={0,0,0}; double NH = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i); double dv2_abs = 0; /* calculate local velocity dispersion (including hubble-flow correction) in physical units */
#ifdef MAGNETIC
                                tempB[0]=SphP[i].B[0];tempB[1]=SphP[i].B[1];tempB[2]=SphP[i].B[2]; //use particle magnetic field
#endif
                                dv2_abs = ((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1])*(SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) // squared norm of the trace-free symmetric [shear] component of the velocity gradient tensor //
                                                    + (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2])*(SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) + (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])*(SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
                                           (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] + SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] + SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) - (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] + SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] + SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2]))) * All.cf_a2inv*All.cf_a2inv;
                                //Saves at formation sink properties in a table: 0:Time 1:ID 2:Mass 3-5:Position 6-8:Velocity 9-11:Magnetic field 12:Internal energy 13:Density 14:cs_effective 15:particle size 16:local surface density 17:local velocity dispersion 18: distance to closest BH
                                fprintf(FdBhFormationDetails,"%g %llu %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n", All.Time, (unsigned long long)P[i].ID, P[i].Mass, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],  P[i].Vel[0], P[i].Vel[1],P[i].Vel[2], tempB[0], tempB[1], tempB[2], SphP[i].InternalEnergyPred, SphP[i].Density * All.cf_a3inv, Get_Gas_effective_soundspeed_i(i) * All.cf_afac3, Get_Particle_Size(i) * All.cf_atime, NH, dv2_abs, P[i].min_dist_to_bh ); fflush(FdBhFormationDetails);
#endif
                            }
#endif // SINGLE_STAR_SINK_DYNAMICS			   
			    if(P[i].Type != 5) {P[i].Type = 4;} // if we didn't set to type 5 above, default to type 4
                        } /* closes final generation from original gas particle */
                        else
                        {
                            /* here we spawn a new star particle */
                            if(NumPart + stars_spawned >= All.MaxPart)
                            {
                                PRINT_WARNING("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)",ThisTask, NumPart, stars_spawned, All.MaxPart);
                                fflush(stdout);
                                endrun(8888);
                            }
                            
                            P[NumPart + stars_spawned] = P[i];
                            P[NumPart + stars_spawned].Type = 4;
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
                            P[NumPart + stars_spawned].DensAroundStar = SphP[i].Density;
#endif
                            NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
                            FirstActiveParticle = NumPart + stars_spawned;
                            NumForceUpdate++;
                            
                            TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;
                            PrevInTimeBin[NumPart + stars_spawned] = i;
                            NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
                            if(NextInTimeBin[i] >= 0) {PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;}
                            NextInTimeBin[i] = NumPart + stars_spawned;
                            if(LastInTimeBin[P[i].TimeBin] == i) {LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;}
                            
                            P[i].ID += ((MyIDType) 1 << (sizeof(MyIDType) * 8 - bits));
                            
                            P[NumPart + stars_spawned].Mass = mass_of_star;
                            P[i].Mass -= P[NumPart + stars_spawned].Mass;
                            if(P[i].Mass<0) P[i].Mass=0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                            SphP[i].MassTrue -= P[NumPart + stars_spawned].Mass;
                            if(SphP[i].MassTrue<0) SphP[i].MassTrue=0;
#endif
                            sum_mass_stars += P[NumPart + stars_spawned].Mass;
                            P[NumPart + stars_spawned].StellarAge = All.Time;
                            
                            force_add_star_to_tree(i, NumPart + stars_spawned);
                            
                            stars_spawned++;
                        }
#ifdef BH_SEED_FROM_LOCALGAS
                    } /* closes else for decision to make a BH particle */
#endif
                }
                
#if defined(METALS) && defined(GALSF_EFFECTIVE_EQS) // does instantaneous enrichment //
                if(P[i].Type == 0)    /* to protect using a particle that has been turned into a star */
                {
                    P[i].Metallicity[0] += (1 - w) * All.SolarAbundances[0] * pfac;
                    if(NUM_METAL_SPECIES>=10) {int k; for(k=1;k<NUM_METAL_SPECIES;k++) {P[i].Metallicity[k] += (1-w) * All.SolarAbundances[k] * pfac;}}
                }
#endif
            } // closes check of flag==0 for star-formation operation
            
#if defined(GALSF_SUBGRID_WINDS)
            if( (flag==0 || All.ComovingIntegrationOn==0) && (P[i].Mass>0) && (P[i].Type==0) && (dtime>0) && (All.Time>0) )
            {
                double pvtau_return[4];
                assign_wind_kick_from_sf_routine(i,sm,dtime,pvtau_return);
            }
#endif
            
        } /* End of If Type = 0 */
    } /* end of main loop over active particles, huzzah! */
    
    
    
    
#if defined(BH_SEED_FROM_LOCALGAS) || defined(SINGLE_STAR_SINK_DYNAMICS)
    MPI_Allreduce(&num_bhformed, &tot_bhformed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if( (ThisTask==0) && (tot_bhformed > 0) )
    {
        printf("BH/Sink formation: %d gas particles converted into BHs\n",tot_bhformed);
        All.TotBHs += tot_bhformed;
    } // if(tot_bhformed > 0)
#endif
    
    MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(tot_spawned > 0 || tot_converted > 0)
    {
        if(ThisTask==0) printf("SFR: spawned %d stars, converted %d gas particles into stars\n", tot_spawned, tot_converted);
        All.TotNumPart += tot_spawned;
        All.TotN_gas -= tot_converted;
        NumPart += stars_spawned;
        /* Note: N_gas is only reduced once rearrange_particle_sequence is called */
        /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    } //(tot_spawned > 0 || tot_converted > 0)
    
    for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++) {if(TimeBinCount[bin]) {sfrrate += TimeBinSfr[bin];}}
    
#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    {
        MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(ThisTask == 0)
        {
            if(All.TimeStep > 0) {rate = total_sm / (All.TimeStep / (All.cf_atime*All.cf_hubble_a));} else {rate = 0;}
            /* convert to solar masses per yr */
            rate_in_msunperyear = rate * UNIT_MASS_IN_SOLAR / UNIT_TIME_IN_YR;
            fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars);
            fflush(FdSfr); // can flush it, because only occuring on domain-level steps anyways
        } // thistask==0
    }
    CPU_Step[CPU_COOLINGSFR] += measure_time();
} /* end of main sfr_cooling routine!!! */





#if defined(GALSF_SUBGRID_WINDS)
void assign_wind_kick_from_sf_routine(int i, double sm, double dtime, double pvtau_return[4])
{
    int j; double v,p,prob, norm, dir[3];

#if (GALSF_SUBGRID_WIND_SCALING == 0)
    /* this is the simple, old standard wind model, with constant velocity & loading with SFR */
    p = All.WindEfficiency * sm / P[i].Mass;
    v = sqrt(2 * All.WindEnergyFraction*All.FactorSN*All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency);
    prob = 1 - exp(-p);
#endif

#if (GALSF_SUBGRID_WIND_SCALING == 1)
    /* wind model where launching scales with halo/galaxy bulk properties (as in Romeel's simulations) */
    if(SphP[i].HostHaloMass > 0 && sm > 0)
    {
        double HaloConcentrationNorm = 9.;  /* concentration c0 of a halo of unit mass */
        double HaloConcentrationSlope = -0.15;  /* slope n of mass concentration relation, namely c = c0 * M_200,crit^n */

        double r200c, v_esc, c_halo, wind_energy, wind_momentum, wind_mass;
        double rhocrit = 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
        rhocrit *= All.OmegaRadiation*All.cf_a2inv*All.cf_a2inv + All.OmegaMatter*All.cf_a3inv + (1-All.OmegaMatter-All.OmegaLambda-All.OmegaRadiation)*All.cf_a2inv + All.OmegaLambda; /* physical critical density at redshift z */

        r200c = pow(SphP[i].HostHaloMass / (4 * M_PI / 3.0 * 200 * rhocrit), 1.0 / 3.0);	/* physical r_200,crit value, assuming FoF mass = M_200,crit */
        v_esc = sqrt(All.G * SphP[i].HostHaloMass / r200c);	/* physical circular velocity at r_200,crit */
        c_halo = HaloConcentrationNorm * pow(SphP[i].HostHaloMass, HaloConcentrationSlope);
        v_esc *= sqrt(2 * c_halo / (log(1 + c_halo) - c_halo / (1 + c_halo)));	/* physical escape velocity of halo */
        v = All.VariableWindVelFactor * v_esc;	/* physical wind velocity */

        wind_momentum = sm * All.VariableWindSpecMomentum;
        wind_energy = sm * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN);

        wind_mass = (wind_energy + sqrt(wind_energy * wind_energy + v * v * wind_momentum * wind_momentum)) / (v * v);
        /* wind mass for this particle, assuming the wind is first given the energy wind_energy and then the momentum wind_momentum */
        p = wind_mass / P[i].Mass;
    }
    else
    {
        v = 0;
        p = 0;
    }
    prob = 1 - exp(-p);
#endif

#if (GALSF_SUBGRID_WIND_SCALING == 2)
    /* wind model where launching scales with halo/galaxy bulk properties (as in Vogelsberger's simulations) */
    if(SphP[i].DM_VelDisp > 0 && sm > 0)
    {
        double wind_energy, wind_momentum, wind_mass;
        v = All.VariableWindVelFactor * SphP[i].DM_VelDisp;  /* physical wind velocity */
        //      if(v < 50.0) v = 50.0;
        wind_momentum = sm * All.VariableWindSpecMomentum;
        wind_energy = sm * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN);
        wind_mass = (wind_energy + sqrt(wind_energy * wind_energy + v * v * wind_momentum * wind_momentum)) / (v * v);
        /* wind mass for this particle, assuming the wind is first given the energy wind_energy and then the momentum wind_momentum */
        p = wind_mass / P[i].Mass;
    }
    else
    {
        v = 0;
        p = 0;
    }
    prob = 1 - exp(-p);
#endif

    if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
    {
#if !defined(GALSF_WINDS_ORIENTATION)
#define GALSF_WINDS_ORIENTATION 0   // determine the wind acceleration orientation //
#endif

#if (GALSF_WINDS_ORIENTATION==0) // random wind direction
        double theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
        double phi = 2 * M_PI * get_random_number(P[i].ID + 4);
        dir[0] = sin(theta) * cos(phi); dir[1] = sin(theta) * sin(phi); dir[2] = cos(theta);
        if(get_random_number(P[i].ID + 5) < 0.5) {for(j=0;j<3;j++) dir[j]=-dir[j];}
#endif
#if (GALSF_WINDS_ORIENTATION==1) // polar wind (defined by accel.cross.vel)
        dir[0] = P[i].GravAccel[1] * P[i].Vel[2] - P[i].GravAccel[2] * P[i].Vel[1];
        dir[1] = P[i].GravAccel[2] * P[i].Vel[0] - P[i].GravAccel[0] * P[i].Vel[2];
        dir[2] = P[i].GravAccel[0] * P[i].Vel[1] - P[i].GravAccel[1] * P[i].Vel[0];
        if(get_random_number(P[i].ID + 5) < 0.5) {for(j=0;j<3;j++) dir[j]=-dir[j];}
#endif
#if (GALSF_WINDS_ORIENTATION==2) // along density gradient //
        for(j=0;j<3;j++) dir[j]=-P[i].GradRho[j];
#endif

        // now actually do the kick for the wind //
        for(j=0,norm=0;j<3;j++) norm+=dir[j]*dir[j];
        if(norm>0) {norm=sqrt(norm);} else {dir[0]=dir[1]=0; dir[2]=norm=1;}
        for(j = 0; j < 3; j++)
        {
            P[i].Vel[j] += v * All.cf_atime * dir[j]/norm;
            SphP[i].VelPred[j] += v * All.cf_atime * dir[j]/norm;
        }
        SphP[i].DelayTime = All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
    } /* if(get_random_number(P[i].ID + 2) < prob) */
}
#endif // defined(GALSF_SUBGRID_WINDS)




#if defined(GALSF_EFFECTIVE_EQS)
/* Routine to initialize quantities needed for the Spingel & Hernquist effective equation of state */
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight, gamma_minus1_eff;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, thresholdStarburst;
#ifdef COOL_METAL_LINES_BY_SPECIES
  int k; double Z[NUM_METAL_SPECIES]; for(k=0;k<NUM_METAL_SPECIES;k++) Z[k]=0;
#endif

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;
      egyhot = All.EgySpecSN / A0;
      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
      u4 = 1.0e4 / (meanweight * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS);
      dens = 1.0e6 * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
      tcool = GetCoolingTime(egyhot, dens, ne, -1);
      coolrate = egyhot / tcool / dens;
      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh = x / pow(1 - x, 2) * (All.FactorSN * All.EgySpecSN - (1 -
						      All.FactorSN) * All.EgySpecCold) / (All.MaxSfrTimescale * coolrate);

      if(ThisTask == 0)
	{
      printf("\n Springel-Hernquist EOS model: A0= %g  \n", A0);
	  printf(" ..computed: PhysDensThresh= %g  (int units)         %g cm^-3\n", All.PhysDensThresh, All.PhysDensThresh * UNIT_DENSITY_IN_NHCGS);
	  printf(" ..expected fraction of cold gas at threshold = %g\n", x);
	  printf(" ..tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);
	}

      dens = All.PhysDensThresh * 10;
      do
	{
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
      tcool = GetCoolingTime(egyhot, dens, ne, -1);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;
      gamma_minus1_eff = (GAMMA_DEFAULT-1);
	  peff = gamma_minus1_eff * dens * egyeff;
	  fac = 1 / (log(dens * 1.025) - log(dens));
	  dens *= 1.025;
	  neff = -log(peff) * fac;
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
      tcool = GetCoolingTime(egyhot, dens, ne, -1);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;
	  peff = gamma_minus1_eff * dens * egyeff;
	  neff += log(peff) * fac;
	}
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;

      if(ThisTask == 0)
	{
	  printf("Run-away sets in for dens=%g\n", thresholdStarburst);
	  printf("Dynamic range for quiescent star formation= %g\n", thresholdStarburst / All.PhysDensThresh);
	  fflush(stdout);
	}

      if(All.ComovingIntegrationOn)
	{
	  All.Time = All.TimeBegin;
	  set_cosmo_factors_for_current_time();
	  IonizeParams();
	}

#if defined(GALSF_SUBGRID_WINDS)
        if(All.WindEfficiency > 0) {if(ThisTask == 0) {printf("Windspeed: %g\n", sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency));}}
#endif
    }
}
#endif // GALSF_EFFECTIVE_EQS //



#endif // GALSF
