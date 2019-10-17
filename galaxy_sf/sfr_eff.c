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


#ifdef GALSF // master switch for compiling the routines below //


#if defined(GALSF_SFR_IMF_VARIATION) || defined(GALSF_SFR_IMF_SAMPLING)
/* function to determine what the IMF of a new star particle will be, based 
    on the gas properties of the particle out of which it forms */
void assign_imf_properties_from_starforming_gas(int i)
{
#ifdef GALSF_SFR_IMF_VARIATION
    double h = Get_Particle_Size(i) * All.cf_atime;
    double cs = Particle_effective_soundspeed_i(i) * All.cf_afac3; // actual sound speed in the simulation: might be unphysically high for SF conditions!
    cs = (1.9e4 / All.UnitVelocity_in_cm_per_s); // set to a minimum cooling temperature, for the actual star-forming conditions. for now, just use a constant //
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
    double M_sonic = cs*cs*cs*cs / (All.G * dv2_abs * h);
    M_sonic *= All.UnitMass_in_g / All.HubbleParam / (1.989e33); // sonic mass in solar units //
    P[i].IMF_Mturnover = DMAX(0.01,DMIN(M_sonic,100.));
    P[i].IMF_Mturnover = 2.0; // 'normal' IMF in our definitions
    
    
    /* now we need to record all the properties we care to save about the star-forming gas, for the sake of later use: */
    int j,k;
    double NH = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1);
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
    double gizmo2gauss = 4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam;
    for(k=0;k<3;k++) {b_mag += Get_Particle_BField(i,k)*Get_Particle_BField(i,k) * gizmo2gauss;}
#endif
    double rad_flux_uv = 1;
    double cr_energy_density = 0;

    P[i].IMF_FormProps[0] = P[i].IMF_Mturnover; // IMF turnover mass as defined above
    P[i].IMF_FormProps[1] = SphP[i].Density * All.cf_a3inv; // density
    P[i].IMF_FormProps[2] = SphP[i].InternalEnergyPred; // thermal internal energy (use to calculate temperature)
    P[i].IMF_FormProps[3] = Particle_effective_soundspeed_i(i) * All.cf_afac3; // sound speed (not trivially related to temperature if CRs, etc included)
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
    double mu = 0.01 * P[i].Mass * All.UnitMass_in_g / All.HubbleParam / (1.989e33); // 1 O-star per 100 Msun
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
        if(fabs(1-(All.Omega0+All.OmegaLambda))<=0.01)
        {
            /* use exact solution for flat universe */
            x0 = (All.Omega0/(1-All.Omega0))/(a0*a0*a0);
            x2 = (All.Omega0/(1-All.Omega0))/(a2*a2*a2);
            age = (2./(3.*sqrt(1-All.Omega0)))*log(sqrt(x0*x2)/((sqrt(1+x2)-1)*(sqrt(1+x0)+1)));
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
    age *= 0.001*All.UnitTime_in_Megayears/All.HubbleParam; // convert to absolute Gyr
    if((age<=1.e-5)||(isnan(age))) {age=1.e-5;}
    return age;
}



/* simple routine to determine density thresholds and other common units for SF routines */
void set_units_sfr(void)
{
    All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
    All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
#ifdef GALSF_EFFECTIVE_EQS
    double meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */
    All.EgySpecCold = All.TempClouds / (meanweight * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS);
    meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
    All.EgySpecSN = All.TempSupernova / (meanweight * (GAMMA_DEFAULT-1) * U_TO_TEMP_UNITS);
#endif // GALSF_EFFECTIVE_EQS
}



/* function which takes properties of a gas particle 'i' and returns probability of its turning into a BH seed particle */
double return_probability_of_this_forming_bh_from_seed_model(int i)
{
    double p=0;
#ifdef BH_SEED_FROM_LOCALGAS
    if(All.Time < 1.0/(1.0+All.SeedBlackHoleMinRedshift)) /* within the allowed redshift range for forming seeds */
    if(SphP[i].Density*All.cf_a3inv > All.PhysDensThresh) /* require it be above the SF density threshold */
    if(P[i].Metallicity[0]/All.SolarAbundances[0] < 0.1) /* and below some metallicity */
    {
        double GradRho = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1) * All.UnitDensity_in_cgs * All.UnitLength_in_cm * All.HubbleParam; /* this gives the Sobolev-estimated column density */
        /* surface dens in g/cm^2; threshold for bound cluster formation in our experiments is ~2 g/cm^2 (10^4 M_sun/pc^2) */
        if (GradRho > 0.1)
        {
            /* now calculate probability of forming a BH seed particle */
            p = P[i].Mass / All.SeedBlackHolePerUnitMass; /* probability of forming a seed per unit mass [in code units] */
            if(p > 1.e-4) {p = 1.-exp(-p);}
            p *= (1-exp(-GradRho/1.0)) * exp(-(P[i].Metallicity[0]/All.SolarAbundances[0])/0.01); /* apply threshold metallicity and density cutoff */
            /* want to add factors to control this probability in zoom-in runs */
        }
    }
#endif
    return p;
}



/* Routine to actually determine the SFR assigned to an individual gas particle at each time */
double get_starformation_rate(int i)
{
    double rateOfSF,tsfr,y; y=0;
    int flag;
#ifdef GALSF_EFFECTIVE_EQS
    double factorEVP, egyhot, ne, tcool, x, cloudmass;
#endif
#ifdef GALSF_SUBGRID_WINDS
    if(SphP[i].DelayTime > 0) return 0;
#endif
    
#ifdef BH_WIND_SPAWN
    if(P[i].ID == All.AGNWindID) return 0;
#endif

    flag = 1;			/* default is normal cooling */
    if(SphP[i].Density*All.cf_a3inv >= All.PhysDensThresh) {flag = 0;}
#if (GALSF_SFR_VIRIAL_SF_CRITERION>=3)
    else {SphP[i].AlphaVirial_SF_TimeSmoothed = 0.;}
#endif
    if(All.ComovingIntegrationOn)
    if(SphP[i].Density < All.OverDensThresh)
    flag = 1;
    if((flag == 1)||(P[i].Mass<=0))
    return 0;
#if (GALSF_SFR_VIRIAL_SF_CRITERION>=3)
    double dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval;
    double dtime = dt / All.cf_hubble_a; /*  the actual time-step */
#endif
    tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * All.cf_a3inv)) * All.MaxSfrTimescale;
    if(tsfr<=0) return 0;
    
    
#ifndef GALSF_EFFECTIVE_EQS
    /* 'normal' sfr from density law above */
    rateOfSF = P[i].Mass / tsfr;
#else
    factorEVP = pow(SphP[i].Density * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;
    egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
    ne = SphP[i].Ne;
    tcool = GetCoolingTime(egyhot, SphP[i].Density * All.cf_a3inv, ne, i);
    y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
    x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
    cloudmass = x * P[i].Mass;
    rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;
    
    update_internalenergy_for_galsf_effective_eos(i,tcool,tsfr,x,rateOfSF); // updates entropies for the effective equation-of-state //
#endif // GALSF_EFFECTIVE_EQS
    
    
#ifdef GALSF_SFR_MOLECULAR_CRITERION
    /* Krumholz & Gnedin fitting function for f_H2 as a function of local properties */
    double tau_fmol = evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1);
    tau_fmol *= (0.1 + P[i].Metallicity[0]/All.SolarAbundances[0]);
    if(tau_fmol>0) {
        tau_fmol *= 434.78 * All.UnitDensity_in_cgs * All.UnitLength_in_cm * All.HubbleParam;
        y = 0.756 * (1 + 3.1*pow(P[i].Metallicity[0]/All.SolarAbundances[0],0.365));
        y = log(1 + 0.6*y + 0.01*y*y) / (0.6*tau_fmol);
        y = 1 - 0.75*y/(1 + 0.25*y);
        if(y<0) y=0; if(y>1) y=1;
        rateOfSF *= y;
    } // if(tau_fmol>0)
#endif // GALSF_SFR_MOLECULAR_CRITERION

    
#ifdef CHIMES_SFR_MOLECULAR_CRITERION 
    /* This is similar to GALSF_SFR_MOLECULAR_CRITERION, except that 
     * the H2 fraction is taken from the CHIMES network. */
    y = ChimesGasVars[i].abundances[H2] * 2.0; 
    if (y < 0) 
      y = 0.0; 
    if (y > 1) 
      y = 1.0; 
    rateOfSF *= y; 
#endif 
    
    
#ifdef GALSF_SFR_VIRIAL_SF_CRITERION
    int j,k; double dv2abs=0, divv=0, gradv[9]={0}; /* calculate local velocity dispersion (including hubble-flow correction) in physical units */
    for(j=0;j<3;j++)
    {
        for(k=0;k<3;k++)
        {
            double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
            if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
            gradv[3*j + k] = vt; // save for possible use below
            if(j==k) {divv += vt;} // save for possible use below
            dv2abs += vt*vt; // save for possible use below
        }
    }
    /* add thermal support, although it is almost always irrelevant on large scales */
    double cs_eff = Particle_effective_soundspeed_i(i);    
    double k_cs = cs_eff / (Get_Particle_Size(i)*All.cf_atime);
    
#ifdef SINGLE_STAR_SINK_FORMATION
#if (defined(COOLING) && !defined(COOL_LOWTEMP_THIN_ONLY)) || defined(EOS_GMC_BAROTROPIC) // if we have to deal with optically-thick thermo
    double nHcgs = HYDROGEN_MASSFRAC * (SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) / PROTONMASS;
    if(nHcgs > 1e13) cs_eff=DMIN(cs_eff, 2e4/All.UnitVelocity_in_cm_per_s); //1.62e5/All.UnitVelocity_in_cm_per_s); // limiter to permit sink formation in simulations that really resolve the opacity limit and bog down when an optically-thick core forms. Modify this if you want to follow first collapse more/less - scale as c_s ~ n^(1/5)
#endif
#ifdef MAGNETIC
    double bmag=0; for(k=0;k<3;k++) {bmag+=Get_Particle_BField(i,k)*Get_Particle_BField(i,k);}
    cs_eff = sqrt(cs_eff*cs_eff + bmag/SphP[i].Density);
#endif
    k_cs = M_PI * cs_eff / (Get_Particle_Size(i)*All.cf_atime);
#endif
                                            
    dv2abs += 2.*k_cs*k_cs; // account for thermal pressure with standard Jeans criterion (k^2*cs^2 vs 4pi*G*rho) //
    double alpha_vir = dv2abs / (8. * M_PI * All.G * SphP[i].Density * All.cf_a3inv); // 1/4 or 1/8 -- going more careful here //
#if (GALSF_SFR_VIRIAL_SF_CRITERION > 0)
    if(alpha_vir < 1.0)
    {   /* check if Jeans mass is remotely close to solar; if not, dont allow it to form 'stars' */
        double q = cs_eff * All.UnitVelocity_in_cm_per_s / (0.2e5);
        double q2 = SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / (HYDROGEN_MASSFRAC*1.0e3*PROTONMASS);
        double MJ_solar = 2.*q*q*q/sqrt(q2);
        double MJ_crit = 1000.;
#ifdef SINGLE_STAR_SINK_FORMATION
        MJ_crit = DMIN(1.e4, DMAX(1.e-3 , 100.*P[i].Mass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS)));
#endif
        if(MJ_solar > MJ_crit) {alpha_vir = 100.;}
    }
#endif
#if (GALSF_SFR_VIRIAL_SF_CRITERION >= 3)
    SphP[i].AlphaVirial_SF_TimeSmoothed += 8.*(1./(1+alpha_vir) - SphP[i].AlphaVirial_SF_TimeSmoothed) * dtime/tsfr;
    if (SphP[i].AlphaVirial_SF_TimeSmoothed < 0.5 || divv >= 0) rateOfSF *= 0.0;
#if (GALSF_SFR_VIRIAL_SF_CRITERION >= 4) 
    // we check that the velocity gradient is negative-definite, ie. converging along all principal axes, which is much stricter than div v < 0
    gsl_matrix_view M = gsl_matrix_view_array (gradv, 3, 3);
    gsl_vector *eval1 = gsl_vector_alloc (3);
    gsl_eigen_symm_workspace *v = gsl_eigen_symm_alloc (3);
    gsl_eigen_symm(&M.matrix, eval1,  v);
    for(k=0; k<3; k++) if (gsl_vector_get(eval1,k) >= 0) rateOfSF = 0; // check each eigenvalue
    gsl_eigen_symm_free (v);
    gsl_vector_free (eval1);
#endif
#elif (GALSF_SFR_VIRIAL_SF_CRITERION > 1)
    if(alpha_vir >= 1.0) {rateOfSF *= 0.0;}
#endif
#if (GALSF_SFR_VIRIAL_SF_CRITERION<3)
    if((alpha_vir<1.0)||(SphP[i].Density*All.cf_a3inv>100.*All.PhysDensThresh)) {rateOfSF *= 1.0;} else {rateOfSF *= 0.0015;} // PFH: note the latter flag is an arbitrary choice currently set -by hand- to prevent runaway densities from this prescription! //
#endif
#endif // GALSF_SFR_VIRIAL_SF_CRITERION

    
#ifdef GALSF_SFR_TIDAL_HILL_CRITERION // we check that the tidal tensor is negative-definite, ie. converging along all principal axes, indicating that we're dominating our environment gravitationally and are living in our own Hill sphere
    {int k; for(k=0;k<3;k++) {if(P[i].tidal_tensorps[k][k] >= 0) {rateOfSF = 0;}}} // we've already diagonized this bad boy in gravtree.c - MYG
#endif

    
#ifdef SINGLE_STAR_SINK_FORMATION
    rateOfSF *= 1.0e10; // make sink formation guaranteed to happen, where it can
#endif
#if (SINGLE_STAR_SINK_FORMATION & 2) // restrict to convergent flows //
    {int k; double divv=0; for(k=0;k<3;k++) {divv += SphP[i].Gradients.Velocity[k][k] * All.cf_a2inv;}
        if(All.ComovingIntegrationOn) {divv += 3.*All.cf_hubble_a;}
        if(divv >= 0) {rateOfSF=0;}}
#endif
#if (SINGLE_STAR_SINK_FORMATION & 4)
    if(SphP[i].Density_Relative_Maximum_in_Kernel > 0) {rateOfSF=0;} // restrict to local density/potential maxima //
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
#ifndef SLOPE2_SINKS
    if(P[i].BH_Ngb_Flag) {rateOfSF=0;} // particle cannot be 'seen' by -any- sink as a potential interacting neighbor //
#endif
    if(P[i].min_dist_to_bh < 1.24*Get_Particle_Size(i)) {rateOfSF=0;} // particle does not see a sink within a volume = 8x=2^3 times its cell volume [set coefficient =1.86 for 27x=3^3 its cell volume] //
#endif
#if (SINGLE_STAR_SINK_FORMATION & 16)
    if(DMIN(P[i].min_bh_approach_time, P[i].min_bh_freefall_time) < tsfr) {rateOfSF = 0;} // probably not about to get gobbled up by a sink before it can collapse //
#endif


    
    return rateOfSF;
}



#ifdef GALSF_EFFECTIVE_EQS
/* compute the 'effective eos' cooling/heating, including thermal feedback sources, here */
void update_internalenergy_for_galsf_effective_eos(int i, double tcool, double tsfr, double x, double rateOfSF)
{
    double dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval, dtime = dt / All.cf_hubble_a; /*  the actual time-step */
    double factorEVP = pow(SphP[i].Density * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP, trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
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
#endif // defined(BH_THERMALFEEDBACK)

    /* now update the thermal variables */
    SphP[i].InternalEnergy = (egyeff + (egycurrent - egyeff) * exp(-dtime / trelax));
    SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
    SphP[i].Pressure = get_pressure(i);
    //SphP[i].dInternalEnergy = 0;
    SphP[i].DtInternalEnergy = 0; /* HERE, it's ok, b/c effective EOS is designed to model new pressure even under compressions, 
                                 (since we're zero'ing the second-half-step from the hydro step) */
}
#endif // GALSF_EFFECTIVE_EQS //




/* master routine for star formation. for 'effective equation of state' models for star-forming gas, this also updates their effective EOS parameters */
void star_formation_parent_routine(void)
{
  int i, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  unsigned int bits;
  double dtime, mass_of_star, p, prob, rate_in_msunperyear, sfrrate, totsfrrate;
  double sum_sm, total_sm, sm=0, rate, sum_mass_stars, total_sum_mass_stars;
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
        dtime = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; /*  the actual time-step */
        
        /* check whether an initial (not fully-complete!) conditions for star formation are fulfilled for a given particle */
        if(SphP[i].Density * All.cf_a3inv >= All.PhysDensThresh) {flag = 0;} // if sufficiently dense, go forward into SF routine //
        if(All.ComovingIntegrationOn) {if(SphP[i].Density < All.OverDensThresh) flag = 1;} // (additional density check for cosmological runs) //

#ifdef GALSF_SUBGRID_WINDS
        if(SphP[i].DelayTime > 0) {flag=1; SphP[i].DelayTime -= dtime;} /* no star formation for particles in the wind; update our wind delay-time calculations */
        if((SphP[i].DelayTime<0) || (SphP[i].Density*All.cf_a3inv < All.WindFreeTravelDensFac*All.PhysDensThresh)) {SphP[i].DelayTime=0;}
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        if(SphP[i].DelayTimeCoolingSNe > 0) {flag=1; SphP[i].DelayTimeCoolingSNe -= dtime;} /* no star formation for particles in the wind; update our wind delay-time calculations */
#endif
        
        
    if((flag == 0)&&(dtime>0)&&(P[i].TimeBin))		/* active star formation (upon start-up, we need to protect against dt==0) */
	    {
          sm = get_starformation_rate(i) * dtime; // expected stellar mass formed this timestep
            // (this also updates entropies for the effective equation-of-state model) //
	      p = sm / P[i].Mass;
	      sum_sm += P[i].Mass * (1 - exp(-p));
            

        /* Alright, now we consider the actual gas-to-star particle conversion and associated steps */

	      /* the upper bits of the gas particle ID store how many stars this gas particle gas already generated */
	      if(bits == 0)
            number_of_stars_generated = 0;
	      else
            number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));

	      mass_of_star = P[i].Mass / (GALSF_GENERATIONS - number_of_stars_generated);
            if(number_of_stars_generated >= GALSF_GENERATIONS-1) mass_of_star=P[i].Mass;

          SphP[i].Sfr = sm / dtime *
            (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
	      if(dtime>0) TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;

          prob = P[i].Mass / mass_of_star * (1 - exp(-p));
        
#if defined(METALS) && defined(GALSF_EFFECTIVE_EQS) // does instantaneous enrichment //
            double w = get_random_number(P[i].ID);
            P[i].Metallicity[0] += w * All.SolarAbundances[0] * (1 - exp(-p));
            if(NUM_METAL_SPECIES>=10)
            {
                int k;
                for(k=1;k<NUM_METAL_SPECIES;k++) {P[i].Metallicity[k] += w * All.SolarAbundances[k] * (1 - exp(-p));}
            }
#endif
            
        if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
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
                P[i].BH_Mdot = 0;
                P[i].DensAroundStar = SphP[i].Density;
            } else {
#endif /* closes ifdef(BH_SEED_FROM_LOCALGAS) */ 

            /* ok, we're going to make a star! */
#if defined(GALSF_SFR_IMF_VARIATION) || defined(GALSF_SFR_IMF_SAMPLING)
            /* if we're allowing for a variable IMF, this is where we will 
                calculate the IMF properties produced from the gas forming stars */
            assign_imf_properties_from_starforming_gas(i);
#endif
                
            if(number_of_stars_generated == (GALSF_GENERATIONS - 1))
		    {
		      /* here we turn the gas particle itself into a star */
		      Stars_converted++;
		      stars_converted++;
		      sum_mass_stars += P[i].Mass;

		      P[i].Type = 4;
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
                P[i].Type = 5;
                num_bhformed++;
                P[i].BH_Mass = All.SeedBlackHoleMass; // if desired to make this appreciable fraction of particle mass, please do so in params file
                TreeReconstructFlag = 1;
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                P[i].SinkRadius = All.ForceSoftening[5];
#ifdef SINGLE_STAR_SINK_DYNAMICS
                double cs = 2e4 / All.UnitVelocity_in_cm_per_s;
#if (defined(COOLING) && !defined(COOL_LOWTEMP_THIN_ONLY)) || defined(EOS_GMC_BAROTROPIC)
                double nHcgs = HYDROGEN_MASSFRAC * (SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) / PROTONMASS;
                if(nHcgs > 1e10) cs *= pow(nHcgs/1e10, 1./5); // if we're getting opacity-limited then we can set a smaller sink radius, since cs ~ n^1/5
#endif
                P[i].SinkRadius = DMAX(3 * P[i].Mass * All.G / (M_PI * cs * cs), All.ForceSoftening[5]); // volume-equivalent particle radius R= (3V/(4PI))^(1/3) at the density where M_Jeans = particle mass
#endif	
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
#ifdef SLOPE2_SINKS
                spin_prefac = sqrt(All.G * P[i].Mass * All.ForceSoftening[5]); // assume material is initially in a circular orbit at the resolution limit
#endif
                P[i].BH_Specific_AngMom[0]=spin_prefac*bh_sin*cos(bh_phi); P[i].BH_Specific_AngMom[1]= spin_prefac * bh_sin*sin(bh_phi); P[i].BH_Specific_AngMom[2]=spin_prefac * bh_mu;
#endif
#ifdef BH_COUNTPROGS
                P[i].BH_CountProgs = 1;
#endif
                P[i].BH_Mdot = 0;
                P[i].DensAroundStar = SphP[i].Density;
#endif // SINGLE_STAR_SINK_DYNAMICS
                
		    } /* closes final generation from original gas particle */
		  else
		    {
		      /* here we spawn a new star particle */

		      if(NumPart + stars_spawned >= All.MaxPart)
			{
			  printf
			    ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			     ThisTask, NumPart, stars_spawned, All.MaxPart);
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
		      if(NextInTimeBin[i] >= 0)
                  PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
		      NextInTimeBin[i] = NumPart + stars_spawned;
		      if(LastInTimeBin[P[i].TimeBin] == i)
                  LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;

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
	    if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
        {
            P[i].Metallicity[0] += (1 - w) * All.SolarAbundances[0] * (1 - exp(-p));
            if(NUM_METAL_SPECIES>=10)
            {
                int k;
                for(k=1;k<NUM_METAL_SPECIES;k++) {P[i].Metallicity[k] += (1-w) * All.SolarAbundances[k] * (1 - exp(-p));}
            }
        }
#endif
        } // closes check of flag==0 for star-formation operation

#if defined(GALSF_SUBGRID_WINDS)
        if( (flag==0 || All.ComovingIntegrationOn==0) &&
           (P[i].Mass>0) && (P[i].Type==0) && (dtime>0) && (All.Time>0) )
        {
            double pvtau_return[4];
            assign_wind_kick_from_sf_routine(i,sm,dtime,pvtau_return);
        }
#endif

	} /* End of If Type = 0 */
    } /* end of main loop over active particles, huzzah! */



    
#if defined(BH_SEED_FROM_LOCALGAS) || defined(SINGLE_STAR_SINK_DYNAMICS)
  MPI_Allreduce(&num_bhformed, &tot_bhformed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_bhformed > 0)
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

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    {
        MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(ThisTask == 0)
        {
            if(All.TimeStep > 0)
                rate = total_sm / (All.TimeStep / (All.cf_atime*All.cf_hubble_a));
            else
                rate = 0;
            /* convert to solar masses per yr */
            rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
            fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars);
            fflush(FdSfr); // can flush it, because only occuring on master steps anyways
        } // thistask==0
    }

    // TO: Don't call rearrange_particle_sequence(). This makes the cell array inconsistent with the tree
    //if(tot_converted+tot_spawned > 0) {rearrange_particle_sequence();}

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
        rhocrit *= All.Omega0/All.cf_a3inv + (1-All.Omega0-All.OmegaLambda)/All.cf_a2inv + All.OmegaLambda; /* physical critical density at redshift z */

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
	  printf(" ..computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
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


