#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! Routines for gas equation-of-state terms (collects things like calculation of gas pressure)
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

/*! this pair of functions: 'return_user_desired_target_density' and 'return_user_desired_target_pressure' should be used
 together with 'HYDRO_GENERATE_TARGET_MESH'. This will attempt to move the mesh and mass
 towards the 'target' pressure profile. Use this to build your ICs.
 The 'desired' pressure and density as a function of particle properties (most commonly, position) should be provided in the function below */
double return_user_desired_target_density(int i)
{
    return 1; // uniform density everywhere -- will try to generate a glass //
    /*
     // this example would initialize a constant-density (density=rho_0) spherical cloud (radius=r_cloud) with a smooth density 'edge' (width=interp_width) surrounded by an ambient medium of density =rho_0/rho_contrast //
     double dx=P[i].Pos[0]-boxHalf_X, dy=P[i].Pos[1]-boxHalf_Y, dz=P[i].Pos[2]-boxHalf_Z, r=sqrt(dx*dx+dy*dy+dz*dz);
     double rho_0=1, r_cloud=0.5*boxHalf_X, interp_width=0.1*r_cloud, rho_contrast=10.;
     return rho_0 * ((1.-1./rho_contrast)*0.5*erfc(2.*(r-r_cloud)/interp_width) + 1./rho_contrast);
     */
}
double return_user_desired_target_pressure(int i)
{
    return 1; // uniform pressure everywhere -- will try to generate a constant-pressure medium //
    /*
     // this example would initialize a radial pressure gradient corresponding to a self-gravitating, spherically-symmetric, infinite power-law
     //   density profile rho ~ r^(-b) -- note to do this right, you need to actually set that power-law for density, too, in 'return_user_desired_target_density' above
     double dx=P[i].Pos[0]-boxHalf_X, dy=P[i].Pos[1]-boxHalf_Y, dz=P[i].Pos[2]-boxHalf_Z, r=sqrt(dx*dx+dy*dy+dz*dz);
     double b = 2.; return 2.*M_PI/fabs((3.-b)*(1.-b)) * pow(return_user_desired_target_density(i),2) * r*r;
     */
}

/*!
Simple getter for the Pressure attribute - will calculate it on-the-fly if EOS quantities are not being cached
 */
double get_pressure(int i) {
#ifndef EOS_PRECOMPUTE
    set_eos_pressure(i);
#endif
    return SphP[i].Pressure;
}

/*!
    Updates the thermodynamic quantities determined by the current internal
    energy, density, and chemistry: pressure, and potentially adiabatic index
    and temperature if those are getting cached.

    this subroutine needs to set the value of the 'press' variable (pressure), which you can see from the
    templates below can follow an arbitrary equation-of-state. for more general equations-of-state you want to specifically set the soundspeed
    variable as well. */
void set_eos_pressure(int i)
{
    double soundspeed, press=0, gamma_eos_index = GAMMA(i); soundspeed=0; /* get effective adiabatic index */
    press = (gamma_eos_index-1) * SphP[i].InternalEnergyPred * Get_Gas_density_for_energy_i(i); /* ideal gas EOS (will get over-written it more complex EOS assumed) */

#if (defined(EOS_PRECOMPUTE) && defined(EOS_CARRIES_TEMPERATURE)) || defined(EOS_SUBSTELLAR_ISM) // will do temperature here
    double ne=1, nh0=0, nHe0, nHepp, nhp, nHeII, temp, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred;
    temp = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
#endif
#ifdef EOS_PRECOMPUTE
#ifdef EOS_CARRIES_TEMPERATURE
    SphP[i].Temperature = temp; // cache the tempature
#endif
#ifdef EOS_CARRIES_GAMMA
    SphP[i].Gamma = gamma_eos(i); // cache the adiabatic index; this will reuse the pre-computed SphP[i].Temperature assigned above
#endif
#endif

#ifdef EOS_SUBSTELLAR_ISM
    press = Get_Gas_density_for_energy_i(i) * BOLTZMANN_CGS * temp / UNIT_ENERGY_IN_CGS / (mu_meanwt * PROTONMASS_CGS / UNIT_MASS_IN_CGS);
#endif
    
#ifdef GALSF_EFFECTIVE_EQS /* modify pressure to 'interpolate' between effective EOS and isothermal, with the Springel & Hernquist 2003 'effective' EOS */
    if(SphP[i].Density*All.cf_a3inv >= All.PhysDensThresh) {press = All.FactorForSofterEQS * press + (1 - All.FactorForSofterEQS) * All.cf_afac1 * (gamma_eos_index-1) * SphP[i].Density * All.InitGasU;}
#endif    
    
#ifdef EOS_HELMHOLTZ /* pass the necessary quantities to wrappers for the Timms EOS */
    struct eos_input eos_in;
    struct eos_output eos_out;
    eos_in.rho  = SphP[i].Density;
    eos_in.eps  = SphP[i].InternalEnergyPred;
    eos_in.Ye   = SphP[i].Ye;
    eos_in.Abar = SphP[i].Abar;
    eos_in.temp = SphP[i].Temperature;
    int ierr = eos_compute(&eos_in, &eos_out);
    assert(!ierr);
    press      = eos_out.press;
    soundspeed = eos_out.csound;
    SphP[i].Temperature = eos_out.temp;
#endif
    
#ifdef EOS_TILLOTSON
    press = calculate_eos_tillotson(i); soundspeed = SphP[i].SoundSpeed; /* done in subroutine, save for below */
#endif
    
#ifdef EOS_ENFORCE_ADIABAT
    press = EOS_ENFORCE_ADIABAT * pow(SphP[i].Density, gamma_eos_index);
#ifdef TURB_DRIVING
    SphP[i].EgyDiss += (SphP[i].InternalEnergy - press / (SphP[i].Density * (gamma_eos_index-1.))); /* save the change in energy imprinted by this enforced equation of state here */
#endif
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = press / (SphP[i].Density * (gamma_eos_index-1.)); /* reset internal energy: particles live -exactly- along this relation */
#endif
    
#ifdef EOS_GMC_BAROTROPIC // barytropic EOS calibrated to Masunaga & Inutsuka 2000, eq. 4 in Federrath 2014 Apj 790. Reasonable over the range of densitites relevant to some small-scale star formation problems
    gamma_eos_index=7./5.; double rho=Get_Gas_density_for_energy_i(i), nH_cgs=rho*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS;
    if(nH_cgs > 2.30181e16) {gamma_eos_index=5./3.;} /* dissociates to atomic if dense enough (hot) */
    if (nH_cgs < 1.49468e8) {press = 6.60677e-16 * nH_cgs;} // isothermal below ~10^8 cm^-3 (adiabatic gamma=5/3 for soundspeed, etc, but this assumes effective eos from cooling, etc
    else if (nH_cgs < 2.30181e11) {press = 1.00585e-16 * pow(nH_cgs, 1.1);} // 'transition' region
    else if (nH_cgs < 2.30181e16) {press = 3.92567e-20 * pow(nH_cgs, gamma_eos_index);} // adiabatic molecular
    else if (nH_cgs < 2.30181e21) {press = 3.1783e-15 * pow(nH_cgs, 1.1);} // 'transition' region
    else {press = 2.49841e-27 * pow(nH_cgs, gamma_eos_index);} // adiabatic atomic
#if CHECK_IF_PREPROCESSOR_HAS_NUMERICAL_VALUE_(EOS_GMC_BAROTROPIC)
#if (EOS_GMC_BAROTROPIC==1) // EOS used in Bate Bonnell & Bromm 2003 and related works - isothermal below 6e10 cm^-3, adiabatic above. Assumes c_s = 200m/s at low density
    if (nH_cgs < 6e10) {press = 6.60677e-16 * nH_cgs;} // isothermal below 6e10 cm^-3 (adiabatic gamma=5/3 for soundspeed, etc, but this assumes effective eos from cooling, etc
    else press = 3.964062e-5 * pow(nH_cgs/6e10,1.4);
#endif
#endif
    press /= UNIT_PRESSURE_IN_CGS;
    /* in this case the EOS is modeling cooling, etc, so -does not- allow shocks or deviations from adiabat, so we reset the internal energy every time this is checked */
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = press / (rho * (gamma_eos_index-1.));
#endif
    
#ifdef COSMIC_RAY_FLUID /* compute the CR contribution to the total pressure and effective soundspeed here */
    double soundspeed2 = gamma_eos_index*(gamma_eos_index-1) * SphP[i].InternalEnergyPred;
    int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        press += Get_Gas_CosmicRayPressure(i,k_CRegy);
        soundspeed2 += GAMMA_COSMICRAY(k_CRegy) * (GAMMA_COSMICRAY(k_CRegy)-1.) * SphP[i].CosmicRayEnergyPred[k_CRegy] / P[i].Mass;
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES // using effective gamma of the alfven component = 3/2
        press += (1.5-1) * SphP[i].Density * (SphP[i].CosmicRayAlfvenEnergy[k_CRegy][0]+SphP[i].CosmicRayAlfvenEnergy[k_CRegy][1]);
        soundspeed2 += 1.5*(1.5-1)*(SphP[i].CosmicRayAlfvenEnergy[k_CRegy][0]+SphP[i].CosmicRayAlfvenEnergy[k_CRegy][1]) / P[i].Mass;
#endif
    }
    soundspeed = sqrt(soundspeed2);
#endif
    
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    soundspeed = sqrt(gamma_eos_index*(gamma_eos_index-1) * SphP[i].InternalEnergyPred + (4./3.)*(1./3.)*SphP[i].SubGrid_CosmicRayEnergyDensity/SphP[i].Density);
    press += (1./3.) * SphP[i].SubGrid_CosmicRayEnergyDensity;
#endif
    
#ifdef RT_RADPRESSURE_IN_HYDRO /* add radiation pressure in the Riemann problem directly */
    int k_freq; double gamma_rad=4./3., fluxlim=1; double soundspeed2 = gamma_eos_index*(gamma_eos_index-1) * SphP[i].InternalEnergyPred;
    if(P[i].Mass>0 && SphP[i].Density>0) {for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
    {
        press += (gamma_rad-1.) * return_flux_limiter(i,k_freq) * SphP[i].Rad_E_gamma_Pred[k_freq] * SphP[i].Density / P[i].Mass;
        soundspeed2 += gamma_rad*(gamma_rad-1.) * SphP[i].Rad_E_gamma_Pred[k_freq] / P[i].Mass;
    }}
    soundspeed = sqrt(soundspeed2);
#endif
    
#if defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE)
    /* add an artificial pressure term to suppress fragmentation at/below the explicit resolution scale */
    double h_eff = DMAX(Get_Particle_Size(i), KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER*ForceSoftening_KernelRadius(i)); /* need to include latter to account for inter-particle spacing << grav soft cases */
    /* standard finite-volume formulation of this (note there is some geometric ambiguity about whether there should be a "pi" in the equation below, but this 
        can be completely folded into the (already arbitrary) definition of NJeans, so we just use the latter parameter */
    double NJeans = 4; // set so that resolution = lambda_Jeans/NJeans -- fragmentation with Jeans/Toomre scales below this will be artificially suppressed now
    double xJeans = (NJeans * NJeans / gamma_eos_index) * All.G * h_eff*h_eff * SphP[i].Density * SphP[i].Density * All.cf_afac1/All.cf_atime;
    if(xJeans>press) press=xJeans;
#endif
    
#if defined(HYDRO_GENERATE_TARGET_MESH)
    press = return_user_desired_target_pressure(i) * (SphP[i].Density / return_user_desired_target_density(i)); // define pressure by reference to 'desired' fluid quantities //
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = return_user_desired_target_pressure(i) / ((gamma_eos_index-1) * SphP[i].Density);
#endif
    
#ifdef EOS_GENERAL /* need to be sure soundspeed variable is set: if not defined above, set it to the default which is given by the effective gamma */
    if(soundspeed == 0) {SphP[i].SoundSpeed = sqrt(gamma_eos_index * press / Get_Gas_density_for_energy_i(i));} else {SphP[i].SoundSpeed = soundspeed;}
#endif

    /* Finally, set the pressure as advertised */
    SphP[i].Pressure = press;
}


/*! this function allows the user to specify an arbitrarily complex adiabatic index. note that for pure adiabatic evolution, one can simply set the pressure to obey some barytropic equation-of-state and use EOS_GENERAL to tell the code to deal with it appropriately.
      but for more general functionality, we want this index here to be appropriately variable. */
double gamma_eos(int i)
{
#if defined(COOL_MOLECFRAC_NONEQM)
    if(i >= 0) {
        if(P[i].Type==0) {
            double fH = HYDROGEN_MASSFRAC, f = SphP[i].MolecularMassFraction, xe = SphP[i].Ne; // use the variables below to update the EOS as needed
            double f_mono = fH*(xe + 1.-f) + (1.-fH)/4., f_di = fH*f/2., gamma_mono=5./3., gamma_di=7./5.; // sum e-, H or p, He, which act monotomic, and molecular, by number
#ifdef EOS_SUBSTELLAR_ISM
            gamma_di = hydrogen_molecule_gamma(get_temperature(i));
#endif
            return 1. + (f_mono + f_di) / (f_mono/(gamma_mono-1.) + f_di/(gamma_di-1.)); // weighted sum by number to compute effective EOS
        }
    }
#endif
    return GAMMA_DEFAULT; /* default to universal constant here */
}

/* Returns the temperature, either pre-computed or calling the routine to re-compute it*/
double get_temperature(int i){
#if defined(EOS_PRECOMPUTE) && defined(EOS_CARRIES_TEMPERATURE)
    if(All.Time > 0 || SphP[i].Temperature > 0) {return SphP[i].Temperature;}
#endif
    return compute_temperature(i);
}

/* Simple wrapper for calling ThermalProperties for temperature only - should only be called by get_temperature() above */
double compute_temperature(int i){
    double ne=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);
    return temperature;
}


/* trivial function to check if particle falls below the minimum allowed temperature */
void check_particle_for_temperature_minimum(int i)
{
    if(All.MinEgySpec)
    {
        if(SphP[i].InternalEnergy < All.MinEgySpec)
        {
            SphP[i].InternalEnergy = All.MinEgySpec;
            SphP[i].DtInternalEnergy = 0;
        }
    }
}



double INLINE_FUNC Get_Gas_density_for_energy_i(int i)
{
#ifdef HYDRO_PRESSURE_SPH
    return SphP[i].EgyWtDensity;
#endif
    return SphP[i].Density;
}


/* calculate the 'total' effective sound speed in a given element [including e.g. cosmic ray pressure and other forms of pressure, if present] */
double INLINE_FUNC Get_Gas_effective_soundspeed_i(int i)
{
#ifdef EOS_GENERAL
    return SphP[i].SoundSpeed;
#else
    /* if nothing above triggers, then we resort to good old-fashioned ideal gas */
    return sqrt(GAMMA(i) * SphP[i].Pressure / Get_Gas_density_for_energy_i(i));
#endif
}


/* calculate the thermal sound speed (using just the InternalEnergy variable) in a given element */
double INLINE_FUNC Get_Gas_thermal_soundspeed_i(int i)
{
    return sqrt(convert_internalenergy_soundspeed2(i,SphP[i].InternalEnergyPred));
}


/* calculate the Alfven speed in a given element */
double Get_Gas_Alfven_speed_i(int i)
{
#if defined(MAGNETIC)
    int k; double bmag=0; for(k=0;k<3;k++) {bmag+=Get_Gas_BField(i,k)*All.cf_a2inv * Get_Gas_BField(i,k)*All.cf_a2inv;}
    if(bmag > 0) {return sqrt(bmag / (MIN_REAL_NUMBER + SphP[i].Density*All.cf_a3inv));}
#endif
    return 0;
}


/* calculate the fast MHD wave speed in a given element */
double Get_Gas_Fast_MHD_wavespeed_i(int i)
{
    double cs = Get_Gas_thermal_soundspeed_i(i), vA = Get_Gas_Alfven_speed_i(i);
    return sqrt(cs*cs + vA*vA);    
}


/* calculate and return the actual B Field of a cell (in comoving units, so multiply by All.cf_a2inv to get physical) */
double INLINE_FUNC Get_Gas_BField(int i_particle_id, int k_vector_component)
{
#if defined(MAGNETIC)
    return SphP[i_particle_id].BPred[k_vector_component] * SphP[i_particle_id].Density / P[i_particle_id].Mass;
#endif
    return 0;
}


/* handy function that just returns the B-field magnitude in microGauss, physical units. purely here to save us time re-writing this */
double get_cell_Bfield_in_microGauss(int i)
{
    double Bmag=0;
#ifdef MAGNETIC
    int k; for(k=0;k<3;k++) {double B=Get_Gas_BField(i,k)*All.cf_a2inv; Bmag+=B*B;} // actual B-field in code units
#else
    Bmag=2.*SphP[i].Pressure*All.cf_a3inv; // assume equipartition
#endif
    return UNIT_B_IN_GAUSS * sqrt(DMAX(Bmag,0)) * 1.e6; // return B in microGauss
}


/* returns the conversion factor to go -approximately- (for really quick estimation) in code units, from internal energy to soundspeed */
double INLINE_FUNC convert_internalenergy_soundspeed2(int i, double u)
{
    double gamma_eos_touse = GAMMA(i);
    return gamma_eos_touse * (gamma_eos_touse-1) * u;
}


/* returns the ionized fraction of gas, meant as a reference for runs outside of the cooling routine which include cooling+other physics */
double Get_Gas_Ionized_Fraction(int i)
{
#ifdef COOLING
#ifdef CHIMES 
  return (double) ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_HII]]; 
#else 
    double ne=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
    double f_ion = DMIN(DMAX(DMAX(DMAX(1-nh0, nhp), ne/1.2), 1.e-8), 1.); // account for different measures above (assuming primordial composition)
    if((!isfinite(f_ion)) || (f_ion<0)) {f_ion=0;}
    return f_ion;
#endif
#endif
    return 1;
}


/* returns the dust-to-metals ratio normalized to the canonical solar value of 1/2: i.e. for 'standard' conditions, should = 1, but if e.g. all dust is sublimated, = 0;
    explicit dust formation-destruction modules should plug in here, to communicate to relevant opacity and other routines in heating/cooling/molecular chemistry cross-code
 T_dust_manual_override here designates whether we want to pass a specific dust temp or allow the code to call one itself here, to avoid creating circular dependencies and to
    allow for self-consistent coupling to various other dust dynamics/formation/chemistry modules. if you just want 'default' behavior and aren't worried about this, call with this parameter set to 0
 */
double return_dust_to_metals_ratio_vs_solar(int i, double T_dust_manual_override)
{
    if(i<0 || P[i].Type!=0) {return 1;}
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS) /* note since the applications of this module really want a -surface area per unit mass- ratio to scale off of, we actually want to scale this by the geometric opacity, relative to what 'typical' solar conditions would give */
    double kappa_interp_geo_cgs = SphP[i].InterpolatedGeometricDustCrossSection / UNIT_SURFDEN_IN_CGS; // this should be in cgs, in cm^2/g
    double kappa_solar_geo_cgs = 3300.; // this is a rough estimate of what one would get as a 'reference' opacity if one assumed a maximum grain size of 0.1 micron, grain density of 2.25 g/cm^3, as in e.g. Weingartner & Draine 2001 (roughly what their models would give, for the size range we usually model, to compare)
    double Z_scaled = P[i].Metallicity[0]/All.SolarAbundances[0]; // metallicity of the particle in solar
    return (kappa_interp_geo_cgs / kappa_solar_geo_cgs) / (Z_scaled); // will be multiplied by metallicity to convert later
#endif
#if defined(GALSF_ISMDUSTCHEM_MODEL) && !defined(GALSF_ISMDUSTCHEM_PASSIVE)
    if(P[i].Metallicity[0]>0) {return (SphP[i].ISMDustChem_Dust_Metal[0]/P[i].Metallicity[0])/0.5;} else {return 0;} // use total amount of dust from 'live' dust evolution models
#endif
#if defined(RT_INFRARED)
    double T_evap = 1500.; // 2e3 * pow(SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS, 0.0195); // latter function from Kuiper 2010 eqs 21-22; sublimation temeprature from Isella & Natta 2005, fit to Pollack 1994. works for protostellar environments, but extrapolates poorly to diffuse ISM and/or stellar/AGN atmosphere environments, so for now use a simpler 1500 K which is a rough median between these, with more sophisticated dust modules required to fit all different parameter regimes.
    double T_dust = T_dust_manual_override; if(T_dust == 0) {T_dust = SphP[i].Dust_Temperature;} // use this iff the dust temp sent is nil
    double Tdust_Tsub = T_dust / T_evap; // ratio for below
    double fdust = sigmoid_sqrt(9.*(1.-Tdust_Tsub)) * exp(-DMIN(40.,Tdust_Tsub*Tdust_Tsub/9.)); // crudely don't bother accounting for size spectrum, just adopt an exponential cutoff above the sublimation temperature
    return DMAX(fdust,1e-25); // floor at value too small to influence physical dust processes, just so dust temp root-finders have something finite and continuous to work with
#endif
#if defined(COOL_LOW_TEMPERATURES) && !defined(SINGLE_STAR_SINK_DYNAMICS) // skip this and assume fdust=1 if SINGLE_STAR_SINK_DYNAMICS on because it uses the fancy dust temp solver whose result depends implicitly on the dust fraction - if sublimation is important then we should be running full RT anyway
    double Tdust = T_dust_manual_override; if(Tdust == 0) {Tdust = get_equilibrium_dust_temperature_estimate(i,0,0);} // call this iff the dust temp sent is nil
    if(Tdust >= 2000.) {return 1.e-4;} else {return exp(-pow(Tdust/1000.,3));} // this hit the maximum allowed temperature in the routine if it gets >2000; for lower temps, let it smoothly cut off
#endif
    return 1; // default behavior
}


/* return an estimate of the Hydrogen molecular fraction of gas, intended for simulations of e.g. molecular clouds, galaxies, and star formation */
double Get_Gas_Molecular_Mass_Fraction(int i, double temperature, double neutral_fraction, double free_electron_ratio, double urad_from_uvb_in_G0)
{
    /* if tracking chemistry explicitly, return the explicitly-evolved H2 fraction */
#ifdef CHIMES // use the CHIMES molecular network for H2
    return DMIN(1,DMAX(0, ChimesGasVars[i].abundances[ChimesGlobalVars.speciesIndices[sp_H2]] * 2.0)); // factor 2 converts to mass fraction in molecular gas, as desired
#endif
    
#if (COOL_GRACKLE_CHEMISTRY >= 2) // Use GRACKLE explicitly-tracked H2 [using the molecular network if this is valid]
    return DMIN(1,DMAX(0, SphP[i].grH2I + SphP[i].grH2II)); // include both states of H2 tracked
#endif
    
#if defined(COOL_MOLECFRAC_NONEQM) // use our simple 1-species network for explicitly-evolved H2 fraction
    return DMIN(1, DMAX(0, SphP[i].MolecularMassFraction));
#endif

#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
    if(SphP[i].DelayTimeHII > 0) {return 0;} // force gas flagged as in HII regions to have zero molecular fraction
#endif

    
#if defined(COOL_MOLECFRAC_LOCALEQM) || defined(COOL_MOLECFRAC_KMT) || defined(COOL_MOLECFRAC_GD) // here are some of the 'fancy' molecular fraction estimators which need various additional properties
    double T=1, nH_cgs=1, Z_Zsol=1, urad_G0=1, xH0=1, x_e=0; // initialize definitions of some variables used below to prevent compiler warnings
    if(temperature > 3.e5) {return 0;} else {T=temperature;} // approximations below not designed for high temperatures, should simply give null
    xH0 = DMIN(DMAX(neutral_fraction,0.),1.); // get neutral fraction [given by call to this program]
    x_e = DMIN(DMAX(free_electron_ratio,0.),2.); // get free electron ratio [number per H nucleon]
    nH_cgs = SphP[i].Density*All.cf_a3inv * UNIT_DENSITY_IN_NHCGS; // get nH defined as number of nucleons per cm^3
    Z_Zsol=1; urad_G0=1; // initialize metal and radiation fields. will assume solar-Z and spatially-uniform Habing field for incident FUV radiation unless reset below.
#ifdef METALS
    Z_Zsol = P[i].Metallicity[0]/All.SolarAbundances[0]; // metallicity in solar units [scale to total Z, since this mixes dust and C opacity], and enforce a low-Z floor to prevent totally unphysical behaviors at super-low Z [where there is still finite opacity in reality; e.g. Kramer's type and other opacities enforce floor around ~1e-3]
#endif
    /* get incident radiation field from whatever module we are using to track it */
#ifdef GALSF_FB_FIRE_RT_UVHEATING
    urad_G0 = DMAX(SphP[i].Rad_Flux_UV, 1.e-10); // note this is ALREADY self-shielded, so we need to be careful about 2x-counting the self-shielding approximation below; hence limit this to a rather sizeable value  //
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
#endif
        
    
#if defined(COOL_MOLECFRAC_LOCALEQM) // ??? -- update to match noneqm fancier cooling functions --
    /* estimate local equilibrium molecular fraction actually using the real formation and destruction rates. expressions for the different rate terms
        as used here are collected in Nickerson, Teyssier, & Rosdahl et al. 2018. Expression for the line self-shielding here
        including turbulent and cell line blanketing terms comes from Gnedin & Draine 2014. below solves this all exactly, using the temperature, metallicity,
        density, ionization states, FUV incident radiation field, and column densities in the simulations. */
    /* take eqm of dot[nH2] = a_H2*rho_dust*nHI [dust formation] + a_GP*nHI*ne [gas-phase formation] + b_3B*nHI*nHI*(nHI+nH2/8) [3-body collisional form] - b_H2HI*nHI*nH2 [collisional dissociation]
        - b_H2H2*nH2*nH2 [collisional mol-mol dissociation] - Gamma_H2^LW * nH2 [photodissociation] - Gamma_H2^+ [photoionization] - xi_H2*nH2 [CR ionization/dissociation] */
    double fH2=0, sqrt_T=sqrt(T), nH0=xH0*nH_cgs, n_e=x_e*nH_cgs, EXPmax=40., clumping_factor=1; // define some variables for below, including neutral H number density, free electron number, etc.
    double f_dustgas_solar = 0.5*Z_Zsol*return_dust_to_metals_ratio_vs_solar(i,0); // dust-to-gas ratio locally
    double a_Z = 3.e-18*sqrt_T / ((1. +4.e-2*sqrt(T+Tdust) +2.e-3*T +8.e-6*T*T )*(1. +1.e4/exp(DMIN(EXPmax,600./Tdust)))) * f_dustgas_solar * nH_cgs * nH0 * clumping_factor; // dust surface formation (assuming dust-to-metals ratio is 0.5*(Z/solar)*dust-to-gas-relative-to-solar in all regions where this is significant), from Glover & Jappsen 2007
    double Tdust = 30.; // need to assume something about dust temperature for reaction rates below for dust-phase formation
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) && (FLAG_NOT_IN_PUBLIC_CODE > 2)) || defined(SINGLE_STAR_SINK_DYNAMICS)
    Tdust = get_equilibrium_dust_temperature_estimate(i, 1, T);
#endif
    //double a_GP = (1.833e-21 * pow(T,0.88)) * nH0 * n_e; // gas-phase formation [old form, from Nickerson et al., appears to be a significant typo in their expression compared to the sources from which they extracted it]
    double a_GP = (1.833e-18 * pow(T,0.88)) * nH0 * n_e / (1. + x_e*1846.*(1.+T/20000.)/sqrt(T)); // gas-phase formation [Glover & Abel 2008, using fitting functions slightly more convenient and assuming H-->H2 much more rapid than other reactions, from Krumholz & McKee 2010; denominator factor accounts for p+H- -> H + H, instead of H2]
    double b_3B = (6.0e-32/sqrt(sqrt_T) + 2.0e-31/sqrt_T) * nH0 * nH0 * nH0; // 3-body collisional formation
    double b_H2HI = (7.073e-19 * pow(T,2.012) * exp(-DMIN(5.179e4/T,EXPmax)) / pow(1. + 2.130e-5*T , 3.512)) * nH0 * (nH0/2.); // collisional dissociation
    b_H2HI += 4.49e-9 * pow(T,0.11) * exp(-DMIN(101858./T,EXPmax)) * (n_e) * (nH0/2.); // collisional H2-e- dissociation [note assuming ground-state optically thin dissociation here as thats where this is most relevant, see Glover+Abel 2008)
    double b_H2H2 = (5.996e-30 * pow(T,4.1881) * exp(-DMIN(5.466e4/T,EXPmax)) / pow(1. + 6.761e-6*T , 5.6881)) * (nH0/2.) * (nH0/2.); // collisional mol-mol dissociation
    double G_LW = 3.3e-11 * urad_G0 * (nH0/2.); // photo-dissociation (+ionization); note we're assuming a spectral shape identical to the MW background mean, scaling by G0
    double xi_cr_H2 = (7.525e-16) * (nH0/2.); // CR dissociation (+ionization)
    // can write this as a quadtratic: 0 = x_a*f^2 - x_b*f + x_c, with f = molec mass fraction
    double x_a = (b_3B + b_H2HI - b_H2H2); // terms quadratic in f -- this term can in principle be positive or negative, usually positive
    double x_b = (a_GP + a_Z + 2.*b_3B + b_H2HI + G_LW + xi_cr_H2); // terms linear in f [note sign, pulling the -sign out here] -- positive-definite
    double x_c = (a_GP + a_Z + b_3B); // terms independent of f -- positive-definite
    double y_a = x_a / (x_c + MIN_REAL_NUMBER), y_b = x_b / (x_c + MIN_REAL_NUMBER), z_a = 4. * y_a / (y_b*y_b + MIN_REAL_NUMBER); // convenient to convert to dimensionless variable needed for checking definite-ness
    if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // checking limits of terms for accuracy

    /* now comes the tricky bit -- need to account for the -molecular- self-shielding [depends on fH2, not just the dust external shielding already accounted for */
    double xb0 = a_GP + a_Z + 2.*b_3B + b_H2HI + xi_cr_H2;
    if(fH2 > 1.e-10 && fH2 < 0.99 && G_LW > 0.1*xb0) // fH2 is non-trivial, and the radiation term is significant, so we need to think about molecular self-shielding
    {
        double fH2_min = fH2; // we have just calculated fH2 with -no- molecular self-shielding, so this number can only go up from here
        // calculate a bundle of variables we will need below, to account for the velocity-gradient Sobolev approximation and slab attenuation of G0 //
        double dx_cell = Get_Particle_Size(i) * All.cf_atime; // cell size
        double surface_density_H2_0 = 5.e14 * PROTONMASS_CGS, x_exp_fac=0.00085, w0=0.2; // characteristic cgs column for -molecular line- self-shielding
        double surface_density_local = xH0 * SphP[i].Density * All.cf_a3inv * dx_cell * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth through the local cell/slab. that's closer to what we want here, since G0 is -already- attenuated in the pre-processing step!
        double v_thermal_rms = 0.111*sqrt(T); // sqrt(3*kB*T/2*mp), since want rms thermal speed of -molecular H2- in kms
        double gradv=velocity_gradient_norm(i);
        double dv_turb=gradv*dx_cell*UNIT_VEL_IN_KMS; // delta-velocity across cell
        double x00 = surface_density_local / surface_density_H2_0, x01 = x00 / (sqrt(1. + 3.*dv_turb*dv_turb/(v_thermal_rms*v_thermal_rms)) * sqrt(2.)*v_thermal_rms), y_ss, x_ss_1, x_ss_sqrt, fH2_tmp, fH2_max, Qmax, Qmin; // variable needed below. note the x01 term corrects following Gnedin+Draine 2014 for the velocity gradient at the sonic scale, assuming a Burgers-type spectrum [their Eq. 3]

        fH2_tmp = 1.; // now consider the maximally shielded case, if you had fmol = 1 in the shielding terms
        x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
        z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
        fH2_max = DMAX(0,DMIN(1,fH2)); // this serves as an upper-limit for f
        
        if(fH2_max > 1.1*fH2_min)
        {
            fH2_tmp = fH2_max; // re-calculate the maximally-shielded case
            x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
            z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
            fH2_tmp=fH2; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
            fH2_max = fH2; Qmax = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // set the new max fH2, from this, and set the corresponding value of the function we are trying to root-find for

            fH2_tmp = fH2_min; // re-calculate the minimally-shielded case
            x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
            z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
            fH2_tmp=fH2; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
            fH2_min = fH2; Qmin = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // set the new min fH2, from this, and set the corresponding value of the function we are trying to root-find for

            fH2 = exp( (log(fH2_min)*Qmax - log(fH2_max)*Qmin) / (Qmax-Qmin) ); // do a Newton-Raphson step in log[f_H2] space now that we have good initial brackets
            if((fH2_max > 1.5*fH2_min) && (Qmax*Qmin < 0) && (fH2_max > 1.1*fH2)) // have a big enough dynamic range, and bracketing Qmin/max, to make further iteration meaningful
            {
                double f_p=fH2_min, Q_p=Qmin, Q, fH2_new; int iter=0; // define variables for iteration below
                while(1)
                {
                    x_ss_1=1.+fH2*x01; x_ss_sqrt=sqrt(1.+fH2*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
                    Q = 1 + y_a*fH2*fH2 - y_b*fH2; // update the value of the function we are trying to zero
                    if(iter==0) {if(Q*Q_p>=0) {f_p=fH2_max; Q_p=Qmax;}} // check in case we attempted to bracket from the 'wrong side'
                    if(Q*Q_p >= 0) {break;} // no longer bracketing, end while loop
                    fH2_new = exp( (log(f_p)*Q - log(fH2)*Q_p) / (Q-Q_p) ); f_p=fH2; fH2=fH2_new; Q_p=Q; // update guess and previous values //
                    iter++; // count iterations
                    if(fabs(fH2-f_p) < 0.1*0.5*(f_p+fH2)) {break;} // converged well enough for our purposes!
                    if((y_ss > 0.85) || (y_ss*G_LW < xb0)) {break;} // negligible shielding, or converged to point where external LW is not dominant dissociator so no further iteration needed
                    if((fH2 > 0.95*fH2_max) || (fH2 > 0.99) || (fH2 < 1.e-10) || (fH2 < 1.05*fH2_min) || (iter > 10)) {break;} // approached physical limits or bounds of validity, or end of iteration cycle
                } // end of convergence iteration to find solution for fmol
            } // opening condition for iteration requiring large enough dynamic range, valid bracketing
        } // opening condition for even checking iteration with fmax > 1.5*fmin
    } // opening condition for considering any molecular self-shielding terms at all
    if(!isfinite(fH2)) {fH2=0;} else {if(fH2>1) {fH2=1;} else if(fH2<0) {fH2=0;}} // check vs nans, valid values
    return xH0 * fH2; // return answer
#endif

    
#if defined(COOL_MOLECFRAC_KMT)
    /* use the simpler Kumholz, McKee, & Tumlinson 2009 sub-grid model for molecular fractions in equilibrium, which is a function modeling spherical clouds
        of internally uniform properties exposed to incident radiation. Depends on column density, metallicity, and incident FUV field. */
    /* get estimate of mass column density integrated away from this location for self-shielding */
    double surface_density_Msun_pc2_infty = 0.05 * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i) * UNIT_SURFDEN_IN_CGS / 0.000208854; // approximate column density with Sobolev or Treecol methods as appropriate; converts to M_solar/pc^2
    /* 0.05 above is in testing, based on calculations by Laura Keating: represents a plausible re-scaling of the shielding length for sub-grid clumping */
    double surface_density_Msun_pc2_local = SphP[i].Density * Get_Particle_Size(i) * All.cf_a2inv * UNIT_SURFDEN_IN_CGS / 0.000208854; // this is -just- the depth through the local cell/slab. that's closer to what we want here, since G0 is -already- attenuated in the pre-processing step!
    double surface_density_Msun_pc2 = DMIN( surface_density_Msun_pc2_local, surface_density_Msun_pc2_infty);
    //double surface_density_Msun_pc2 = surface_density_Msun_pc2_local;
    /* now actually do the relevant calculation with the KMT fitting functions */
    double clumping_factor_for_unresolved_densities = 1; // Gnedin et al. add a large clumping factor to account for inability to resolve high-densities, here go with what is resolved
    double chi = 0.766 * (1. + 3.1*pow(Z_Zsol, 0.365)); // KMT estimate of chi, useful if we do -not- know anything about actual radiation field
    if(urad_G0 >= 0) {chi = 71. * urad_G0 / (clumping_factor_for_unresolved_densities * nH_cgs);} // their actual fiducial value including radiation information
    double psi = chi * (1.+0.4*chi)/(1.+1.08731*chi); // slightly-transformed chi variable
    double s = (Z_Zsol + 1.e-3) * surface_density_Msun_pc2 / (MIN_REAL_NUMBER + psi); // key variable controlling shielding in the KMT approximaton
    double q = s * (125. + s) / (11. * (96. + s)); // convert to more useful form from their Eq. 37
    double fH2 = 1. - pow(1.+q*q*q , -1./3.); // full KMT expression [unlike log-approximation, this extrapolates physically at low-q]
    if(q<0.2) {fH2 = q*q*q * (1. - 2.*q*q*q/3.)/3.;} // catch low-q limit more accurately [prevent roundoff error problems]
    if(q>10.) {fH2 = 1. - 1./q;} // catch high-q limit more accurately [prevent roundoff error problems]
    fH2 = DMIN(1,DMAX(0, fH2)); // multiple by neutral fraction, as this is ultimately the fraction of the -neutral- gas in H2
    return xH0 * fH2;
#endif
    

#if defined(COOL_MOLECFRAC_GD)
    /* use the sub-grid final expression calibrated to ~60pc resolution simulations with equilibrium molecular chemistry and post-processing radiative
        transfer from Gnedin & Draine 2014 (Eqs. 5-7) */
    double S_slab = Get_Particle_Size(i) * All.cf_atime * UNIT_LENGTH_IN_PC / 100.; // slab size in units of 100 pc
    double D_star = 0.17 * (2. + S_slab*S_slab*S_slab*S_slab*S_slab) / (1. + S_slab*S_slab*S_slab*S_slab*S_slab); // intermediate variable
    double U_star = 9. * D_star / S_slab, n_star = 14. * sqrt(D_star) / S_slab; // intermediate variables
    double g_eff = sqrt(D_star*D_star + Z_Zsol*Z_Zsol); // intermediate variable parameterizing the dust-to-gas ratio here [assuming the dust-to-gas ratio relative to solar scales linearly with metallicity, giving Z_Zsol = their D_MW parameter]
    double Lambda_incident = log(1. + pow(0.05/g_eff + urad_G0, 2./3.) * pow(g_eff, 1./3.) / U_star); // intermediate variable parameterizing the incident radiation, takes input UV radiation field relative to MW
    double nHalf = n_star * Lambda_incident / g_eff; // intermediate variable
    double w_x = 0.8 + sqrt(Lambda_incident) / pow(S_slab, 1./3.); // intermediate variable
    double x_f = w_x * log(nH_cgs / nHalf); // intermediate variable
    double fH2_gd = 1./(1. + exp(-x_f*(1.-0.02*x_f+0.001*x_f*x_f)));
    return xH0 * fH2_gd;
#endif
    
    
#if (SINGLE_STAR_SINK_FORMATION & 256) || defined(GALSF_SFR_MOLECULAR_CRITERION) || defined(COOL_MOLECFRAC_KG) /* estimate f_H2 with Krumholz & Gnedin 2010 fitting function, assuming simple scalings of radiation field, clumping, and other factors with basic gas properties so function only of surface density and metallicity, truncated at low values (or else it gives non-sensical answers) */
    double clumping_factor=1, fH2_kg=0, tau_fmol = (0.1 + P[i].Metallicity[0]/All.SolarAbundances[0]) * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i) * 434.78 * UNIT_SURFDEN_IN_CGS; // convert units for surface density. also limit to Z>=0.1, where their fits were actually good, or else get unphysically low molecular fractions
    if(tau_fmol>0) {double y = 0.756 * (1 + 3.1*pow(P[i].Metallicity[0]/All.SolarAbundances[0],0.365)) / clumping_factor; // this assumes all the equilibrium scalings of radiation field, density, SFR, etc, to get a trivial expression
        y = log(1 + 0.6*y + 0.01*y*y) / (0.6*tau_fmol); y = 1 - 0.75*y/(1 + 0.25*y); fH2_kg=DMIN(1,DMAX(0,y));}
    return fH2_kg * neutral_fraction;
#endif
    
    
#if defined(COOLING) || defined(COOL_MOLECFRAC_GC) /* if none of the above is set, default to a wildly-oversimplified scaling set by fits to the temperature below which gas at a given density becomes molecular from cloud simulations in Glover+Clark 2012 */
    double T_mol = DMAX(1.,DMIN(8000., SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS));
    return neutral_fraction / (1. + temperature*temperature/(T_mol*T_mol));
#endif
    
    return 0; // catch //
}


/* return helium -number- fraction, not mass fraction */
double yhelium(int target)
{
#ifdef COOL_METAL_LINES_BY_SPECIES
    if(target >= 0) {double ytmp=DMIN(0.5,P[target].Metallicity[1]); return 0.25*ytmp/(1.-ytmp);} else {return ((1.-HYDROGEN_MASSFRAC)/(4.*HYDROGEN_MASSFRAC));}
#else
    return ((1.-HYDROGEN_MASSFRAC)/(4.*HYDROGEN_MASSFRAC)); // assume uniform H-He gas
#endif
}


/* return mean molecular weight, appropriate for the approximations of the user-selected chemical network[s] */
double Get_Gas_Mean_Molecular_Weight_mu(double T_guess, double rho, double *xH0, double *ne_guess, double urad_from_uvb_in_G0, int target)
{
#if defined(CHIMES)
    return calculate_mean_molecular_weight(&(ChimesGasVars[target]), &ChimesGlobalVars);
#elif defined(COOLING)
    double X=HYDROGEN_MASSFRAC, Y=1.-X, Z=0, fmol;
#ifdef METALS
    if(target >= 0)
    {
        Z = DMIN(0.25,P[target].Metallicity[0]); if(NUM_METAL_SPECIES>=10) {Y = DMIN(0.35,P[target].Metallicity[1]);}
        X = 1. - (Y+Z);
    }
#endif
    fmol = Get_Gas_Molecular_Mass_Fraction(target, T_guess, *xH0, *ne_guess, urad_from_uvb_in_G0); /* use our simple subroutine to estimate this, ignoring UVB and with clumping factor=1 */
    return 1. / ( X*(1-0.5*fmol) + Y/4. + *ne_guess*HYDROGEN_MASSFRAC + Z/(16.+12.*fmol) ); // since our ne is defined in some routines with He, should multiply by universal
#else
    return 4./(3.+5.*HYDROGEN_MASSFRAC); // fully-ionized H-He plasma
#endif
}



/* subroutine to calculate the conductivities and assign the various non-ideal MHD coefficients (Ohmic, Hall, Ambipolar). can modify or expand chemistry or assumptions about grains herein. */
void calculate_and_assign_nonideal_mhd_coefficients(int i)
{
#ifdef MHD_NON_IDEAL
#ifdef COOLING // only try to get self-consistent resistivities after the first timestep, when we have calculated the self-consistent ionization state - on the first timestep default to 0 (=ideal MHD)
    if(All.Time <= All.TimeBegin) {SphP[i].Eta_MHD_OhmicResistivity_Coeff = SphP[i].Eta_MHD_HallEffect_Coeff = SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = 0; return;} // =0 on the first timestep, since we don't know the ionization yet
#endif
    /* calculations below follow Wardle 2007 and Keith & Wardle 2014, for the equation sets */
    double mean_molecular_weight = 2.38; // molecular H2, +He with solar mass fractions and metals
    double a_grain_micron = 0.1, f_dustgas = 0.01; // effective size of grains that matter at these densities
    double m_ion = 24.3; // Mg dominates ions in dense gas [where this is relevant]; this is ion mass in units of proton mass
    double zeta_cr = Get_CosmicRayIonizationRate_cgs(i); // cosmic ray ionization rate (fixed as constant for non-CR runs)
#ifdef COOLING
    double T_eff_atomic = 1.23 * (5./3.-1.) * U_TO_TEMP_UNITS * SphP[i].InternalEnergyPred; /* we'll use this to make a quick approximation to the actual mean molecular weight here */
    double nH_cgs = SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, T_transition=DMIN(8000.,nH_cgs), f_mol=1./(1. + T_eff_atomic*T_eff_atomic/(T_transition*T_transition));
    mean_molecular_weight = 4. / (1. + (3. + 4.*SphP[i].Ne - 2.*f_mol) * HYDROGEN_MASSFRAC);
#endif
#ifdef METALS
    f_dustgas = 0.5 * P[i].Metallicity[0] * return_dust_to_metals_ratio_vs_solar(i,0); // appropriate dust-to-metals ratio
#endif
    double temperature = mean_molecular_weight * (GAMMA(i)-1.) * U_TO_TEMP_UNITS * SphP[i].InternalEnergyPred; // will use appropriate EOS to estimate temperature
    // now everything should be fully-determined (given the inputs above and the known properties of the gas) //
    double m_neutral = mean_molecular_weight; // in units of the proton mass
    double ag01 = a_grain_micron/0.1, m_grain = 7.51e9 * ag01*ag01*ag01; // grain mass [internal density =3 g/cm^3]
    double rho = SphP[i].Density*All.cf_a3inv * UNIT_DENSITY_IN_CGS, n_eff = rho / PROTONMASS_CGS; // density in cgs
    // calculate ionization fraction in dense gas; use rate coefficients k to estimate grain charge
    double k0 = 1.95e-4 * ag01*ag01 * sqrt(temperature); // prefactor for rate coefficient for electron-grain collisions
    double ngr_ngas = (m_neutral/m_grain) * f_dustgas; // number of grains per neutral
    double psi_prefac = 167.1 / (ag01 * temperature); // e*e/(a_grain*k_boltzmann*T): Z_grain = psi/psi_prefac where psi is constant determines charge
    double alpha = zeta_cr * psi_prefac / (ngr_ngas*ngr_ngas * k0 * (n_eff/m_neutral)); // coefficient for equation that determines Z_grain
    // psi solves the equation: psi = alpha * (exp[psi] - y/(1+psi)) where y=sqrt(m_ion/m_electron); note the solution for small alpha is independent of m_ion, only large alpha
    //   (where the non-ideal effects are weak, generally) produces a difference: at very high-T, appropriate m_ion should be hydrogen+helium, but in this limit our cooling
    //    routines will already correctly determine the ionization states. so we can safely adopt Mg as our ion of consideration
    double y=sqrt(m_ion*PROTONMASS_CGS/ELECTRONMASS_CGS), psi_0 = 0.5188025-0.804386*log(y), psi=psi_0; // solution for large alpha [>~10]
    if(alpha<0.002) {psi=alpha*(1.-y)/(1.+alpha*(1.+y));} else if(alpha<10.) {psi=psi_0/(1.+0.027/alpha);} // accurate approximation for intermediate values we can use here
    double k_e = k0 * exp(psi); // e-grain collision rate coefficient
    double k_i = k0 * sqrt(ELECTRONMASS_CGS / (m_ion*PROTONMASS_CGS)) * (1 - psi); // i-grain collision rate coefficient
    double n_elec = zeta_cr / (ngr_ngas * k_e); // electron number density
    double n_ion = zeta_cr / (ngr_ngas * k_i); // ion number density
    double Z_grain = psi / psi_prefac; // mean grain charge (note this is signed, will be negative)
#ifdef COOLING
    double mu_eff=2.38, x_elec=DMAX(1.e-18, SphP[i].Ne*HYDROGEN_MASSFRAC*mu_eff), R=x_elec*psi_prefac/ngr_ngas; psi_0=-3.787124454911839; n_elec=x_elec*n_eff/mu_eff; // R is essentially the ratio of negative charge in e- to dust: determines which regime we're in to set quantities below
    if(R > 100.) {psi=psi_0;} else if(R < 0.002) {psi=R*(1.-y)/(1.+2.*y*R);} else {psi=psi_0/(1.+pow(R/0.18967,-0.5646));} // simple set of functions to solve for psi, given R above, using the same equations used to determine low-temp ion fractions
    n_ion = n_elec * y * exp(psi)/(1.-psi); Z_grain = psi / psi_prefac; // we can immediately now calculate these from the above
#endif
    // now define more variables we will need below //
    double gizmo2gauss = UNIT_B_IN_GAUSS; // convert to B-field to gauss (units)
    double B_Gauss = 0; int k; for(k=0;k<3;k++) {B_Gauss += Get_Gas_BField(i,k)*Get_Gas_BField(i,k);} // get magnitude of B //
    if(B_Gauss<=0) {B_Gauss=0;} else {B_Gauss = sqrt(B_Gauss) * All.cf_a2inv * gizmo2gauss;} // B-field magnitude in Gauss
    double xe = n_elec / n_eff;
    double xi = n_ion / n_eff;
    double xg = ngr_ngas;
    // get collision rates/cross sections for different species //
    double nu_g = 7.90e-6 * ag01*ag01 * sqrt(temperature/m_neutral) / (m_neutral+m_grain); // Pinto & Galli 2008
    double nu_ei = 51.*xe*pow(temperature,-1.5); // Pandey & Wardle 2008 (e-ion)
    double nu_e = nu_ei + 6.21e-9*pow(temperature/100.,0.65)/m_neutral; // Pinto & Galli 2008 for latter (e-neutral)
    double nu_ie = ((ELECTRONMASS_CGS*xe)/(m_ion*PROTONMASS_CGS*xi))*nu_ei; // Pandey & Wardle 2008 for former (e-ion)
    double nu_i = nu_ie + 1.57e-9/(m_neutral+m_ion); // Pandey & Wardle 2008 for former (e-ion), Pinto & Galli 2008 for latter (i-neutral)
    // use the cross sections to determine the hall parameters and conductivities //
    double beta_prefac = ELECTRONCHARGE_CGS * B_Gauss / (PROTONMASS_CGS * C_LIGHT_CGS * n_eff);
    double beta_i = beta_prefac / (m_ion * nu_i); // standard beta factors (Hall parameters)
    double beta_e = beta_prefac / (ELECTRONMASS_CGS/PROTONMASS_CGS * nu_e);
    double beta_g = beta_prefac / (m_grain * nu_g) * fabs(Z_grain);
    double be_inv = 1/(1 + beta_e*beta_e), bi_inv = 1/(1 + beta_i*beta_i), bg_inv = 1/(1 + beta_g*beta_g);
    double sigma_O = xe*beta_e + xi*beta_i + xg*fabs(Z_grain)*beta_g; // ohmic conductivity
    double sigma_H = -xe*be_inv + xi*bi_inv + xg*Z_grain*bg_inv; // hall conductivity
    double sigma_P = xe*beta_e*be_inv + xi*beta_i*bi_inv + xg*fabs(Z_grain)*beta_g*bg_inv; // pedersen conductivity
    double sign_Zgrain = Z_grain/fabs(Z_grain); if(Z_grain==0) {sign_Zgrain=0;}
    double sigma_A2 = (xe*beta_e*be_inv)*(xi*beta_i*bi_inv)*pow(beta_i+beta_e,2) +
    (xe*beta_e*be_inv)*(xg*fabs(Z_grain)*beta_g*bg_inv)*pow(sign_Zgrain*beta_g+beta_e,2) +
    (xi*beta_i*bi_inv)*(xg*fabs(Z_grain)*beta_g*bg_inv)*pow(sign_Zgrain*beta_g-beta_i,2); // alternative formulation which is automatically positive-definite
    // now we can finally calculate the diffusivities //
    double eta_prefac = B_Gauss * C_LIGHT_CGS / (4 * M_PI * ELECTRONCHARGE_CGS * n_eff );
    double eta_O = eta_prefac / sigma_O;
    double sigma_perp2 = sigma_H*sigma_H + sigma_P*sigma_P;
    double eta_H = eta_prefac * sigma_H / sigma_perp2;
    double eta_A = eta_prefac * (sigma_A2)/(sigma_O*sigma_perp2);
    eta_O = fabs(eta_O); eta_A = fabs(eta_A); // these depend on the absolute values and should be written as such, so eta is always positive [not true for eta_h]
    // convert units to code units
    double x_neutral = DMAX(0., 1-m_ion*xi); // neutral fraction of mass
    double units_cgs_to_code = UNIT_TIME_IN_CGS / (UNIT_LENGTH_IN_CGS * UNIT_LENGTH_IN_CGS); // convert coefficients (L^2/t) to code units [physical]
    double eta_ohmic = eta_O*units_cgs_to_code, eta_hall = eta_H*units_cgs_to_code, eta_ad = x_neutral * eta_A*units_cgs_to_code;
//#define MHD_NON_IDEAL_CORRECTIONTERMS 1
#ifdef MHD_NON_IDEAL_CORRECTIONTERMS /* account for unphysical or not internally self-consistent drift/slip speeds (PFH+Squire 24) */
    double gradbmag2=0,gradbmag=0,btmp=0,L_B=MAX_REAL_NUMBER;
    int j; for(k=0;k<3;k++) {for(j=0;j<3;j++) {btmp=SphP[i].Gradients.B[k][j]; gradbmag2+=btmp*btmp;}} // need to get magnitude of B gradient for below
    if(gradbmag2>0) {gradbmag=sqrt(gradbmag2)*(All.cf_a2inv/All.cf_atime)*gizmo2gauss; L_B=(B_Gauss/gradbmag)*UNIT_LENGTH_IN_CGS;} // L_B is gradient length in cgs
    double xi_AbsZi_eff = (xe*beta_e + (xi+1.e-25)*beta_i + xg*fabs(Z_grain)*beta_g) / (beta_e + beta_i + beta_g); // weighted mean xi_qi to use
    double psi_n = (xe/(1.+beta_e) + xi/(1.+beta_i) + xg*fabs(Z_grain)/(1.+beta_g)) / (xe+xi+xg*fabs(Z_grain) + 1.e-20); // coupling parameter for weight of neutrals in effective speed
    if(xe+xi+xg < 1.e-20) {psi_n = 1.;}
    double m_carrier_weighted = PROTONMASS_CGS * (xe*ELECTRONMASS_CGS/PROTONMASS_CGS + xi*m_ion + xg*fabs(Z_grain)*m_grain + psi_n*m_neutral); // effective weight of the dragged carriers for speeds below
    double vT_crit = sqrt( BOLTZMANN_CGS*temperature * (xe + xi + xg*fabs(Z_grain) + psi_n) / m_carrier_weighted ); // salient thermal speed for superthermal drift
    double vA_crit = B_Gauss / sqrt(4.*M_PI*m_carrier_weighted*n_eff); // effective Alfven speed to compare as well
    double vT_an = DMIN(vT_crit, vA_crit); // speed for anomalous term to kick on
    double vdrift_mag = ((B_Gauss / L_B) * C_LIGHT_CGS) / (4.*M_PI*ELECTRONCHARGE_CGS * n_eff * xi_AbsZi_eff), vslip_mag = eta_A/L_B, epstein_corr = sqrt(1. + (vdrift_mag*vdrift_mag + vslip_mag*vslip_mag)/(vT_crit*vT_crit));
    double eta_an = (eta_prefac * units_cgs_to_code / xi_AbsZi_eff) * sqrt(1. + 4.*M_PI*C_LIGHT_CGS*C_LIGHT_CGS*ELECTRONMASS_CGS*n_eff*xe/(B_Gauss*B_Gauss)); // anomalous resistivity, set to max of either gyro frequency or electron plasma frequency (dust plasma frequency should be lower unless in dusty plasma regime, where behavior is much more complicated
    eta_ohmic*=epstein_corr; eta_ad/=sqrt(epstein_corr); // epstein-type correction for superthermal drift or slip
    if(vdrift_mag > vT_an/40.) {eta_ohmic += eta_an * exp(-vT_an/vdrift_mag) * sqrt(1. + (vdrift_mag*vdrift_mag)/(vT_an*vT_an));} // factor to represent boosted ohmic to diffuse when exceed critical limit
    double eta_max = 1.e24;
    if(eta_ohmic > eta_max) {eta_ohmic = eta_max;}
#endif
    SphP[i].Eta_MHD_OhmicResistivity_Coeff = eta_ohmic;     /*!< Ohmic resistivity coefficient [physical units of L^2/t] */
    SphP[i].Eta_MHD_HallEffect_Coeff = eta_hall;            /*!< Hall effect coefficient [physical units of L^2/t] */
    SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = eta_ad;      /*!< Hall effect coefficient [physical units of L^2/t] */
#endif
}


/* subroutine to calculate and assign the Spitzer-Braginskii conduction and viscosity coefficients, with appropriate limiters for high-beta plasmas, saturated behaviors, and optional practical numerical limiters when in extreme regimes */
void calculate_and_assign_conduction_and_viscosity_coefficients(int i)
{
#if defined(CONDUCTION) || defined(VISCOSITY)
    
#if (defined(CONDUCTION_SPITZER) || defined(VISCOSITY_BRAGINSKII))
    double ion_frac, beta_i; int k; k=0; ion_frac=1; beta_i=1.e20; double rho=SphP[i].Density*All.cf_a3inv, u_int=SphP[i].InternalEnergyPred;
#if defined(COOLING) /* get the ionized fraction. NOTE we CANNOT call 'ThermalProperties' or functions like 'Get_Ionized_Fraction' here in gradients.c, as we have not done self-shielding steps yet and most modules will yield unphysical answers! */
    ion_frac = SphP[i].Ne / (1. + 2.*yhelium(i)); /* quick estimator. this is actually what we need for conduction since its the free electrons conducting, and we want number relative to fully-ionized gas */
#endif
    double vf_lim,cs,cs_therm; cs=Get_Gas_effective_soundspeed_i(i); cs_therm=Get_Gas_thermal_soundspeed_i(i); vf_lim = cs*All.cf_afac3;
#ifdef MAGNETIC
    double bhat[3]={0},bmag=0,double_dot_dv=0; for(k=0;k<3;k++) {bhat[k]=Get_Gas_BField(i,k); bmag+=bhat[k]*bhat[k];}
    if(bmag>0) {bmag = sqrt(bmag); for(k=0;k<3;k++) {bhat[k]/=bmag;}}
    beta_i = bmag*bmag * All.cf_afac1 * All.cf_a3inv / (All.cf_atime * rho * cs_therm * cs_therm);
    vf_lim *= DMIN(1.e4 , sqrt(1.+beta_i));
#endif
#endif
    
#ifdef CONDUCTION
    SphP[i].Kappa_Conduction = All.ConductionCoeff;
#ifdef CONDUCTION_SPITZER
    /* calculate the thermal conductivities: use the Spitzer formula */
    SphP[i].Kappa_Conduction *= ion_frac * pow(u_int, 2.5);
    /* account for saturation (when the mean free path of electrons is large): estimate whether we're in that limit with the gradients */
    double electron_free_path = All.ElectronFreePathFactor * u_int * u_int / rho;
    double du_conduction=0; for(k=0;k<3;k++) {du_conduction += SphP[i].Gradients.InternalEnergy[k] * SphP[i].Gradients.InternalEnergy[k];}
    double temp_scale_length = u_int / sqrt(du_conduction) * All.cf_atime;
    /* also the Whistler instability limits the heat flux at high-beta; Komarov et al., arXiv:1711.11462 (2017) */
    SphP[i].Kappa_Conduction /= (1 + (4.2 + 1./(3.*beta_i)) * electron_free_path / temp_scale_length); /* should be in physical units */
#ifdef DIFFUSION_OPTIMIZERS
    SphP[i].Kappa_Conduction = DMIN(SphP[i].Kappa_Conduction , 42.85 * rho * vf_lim * DMIN(20.*Get_Particle_Size(i)*All.cf_atime , temp_scale_length));
#endif
#endif
#endif
    
#ifdef VISCOSITY
    SphP[i].Eta_ShearViscosity = All.ShearViscosityCoeff; SphP[i].Zeta_BulkViscosity = All.BulkViscosityCoeff;
#ifdef VISCOSITY_BRAGINSKII
    /* calculate the viscosity coefficients: use the Braginskii shear tensor formulation expanded to first order */
    SphP[i].Eta_ShearViscosity *= ion_frac * pow(u_int, 2.5); SphP[i].Zeta_BulkViscosity = 0;
    /* again need to account for possible saturation (when the mean free path of ions is large): estimate whether we're in that limit with the gradients */
    double ion_free_path = All.ElectronFreePathFactor * u_int * u_int / rho; double dv_magnitude=0, v_magnitude=0;
    /* need an estimate of the internal energy gradient scale length, which we get by d(P/rho) = P/rho * (dP/P - drho/rho) */
    for(k=0;k<3;k++) {int k1;
        for(k1=0;k1<3;k1++) {
            dv_magnitude += SphP[i].Gradients.Velocity[k][k1]*SphP[i].Gradients.Velocity[k][k1];
#ifdef MAGNETIC
            double_dot_dv += SphP[i].Gradients.Velocity[k][k1] * bhat[k]*bhat[k1] * All.cf_a2inv; // physical units
#endif
        }
        v_magnitude += SphP[i].VelPred[k]*SphP[i].VelPred[k];
    }
    double vel_scale_length = sqrt( v_magnitude / dv_magnitude ) * All.cf_atime;
    /* also limit to saturation magnitude ~ signal_speed / lambda_MFP^2 */
    SphP[i].Eta_ShearViscosity /= (1 + 4.2 * ion_free_path / vel_scale_length); /* should be in physical units */
#ifdef MAGNETIC
    // following Jono Squire's notes, the mirror and firehose instabilities limit pressure anisotropies [which scale as the viscous term inside the gradient: nu_braginskii*(bhat.bhat:grad.v)] to >-2*P_magnetic and <1*P_magnetic
    double P_effective_visc = SphP[i].Eta_ShearViscosity * double_dot_dv;
    double P_magnetic = 0.5 * (bmag*All.cf_a2inv) * (bmag*All.cf_a2inv);
    if(P_effective_visc < -2.*P_magnetic) {SphP[i].Eta_ShearViscosity = 2.*P_magnetic / fabs(double_dot_dv);}
    if(P_effective_visc > P_magnetic) {SphP[i].Eta_ShearViscosity = P_magnetic / fabs(double_dot_dv);}
#endif
#ifdef DIFFUSION_OPTIMIZERS
    double eta_sat = rho * vf_lim / (ion_free_path * (1 + 4.2 * ion_free_path / vel_scale_length));
    if(eta_sat <= 0) {SphP[i].Eta_ShearViscosity=0;}
    if(SphP[i].Eta_ShearViscosity>0) {SphP[i].Eta_ShearViscosity = 1. / (1./SphP[i].Eta_ShearViscosity + 1./eta_sat);} // again, all physical units //
#endif
#endif
#endif
    
#endif
}


#ifdef TURB_DIFFUSION
void calculate_and_assign_turbulent_diffusion_coefficients(int i)
{
    /* estimate local turbulent diffusion coefficient from velocity gradients using Smagorinsky mixing model:
     we do this after slope-limiting to prevent the estimated velocity gradients from being unphysically large */
    double h_turb = Get_Particle_Size(i) * All.cf_atime; // physical
    if(h_turb > 0)
    {
        // overall normalization //
        double C_Smagorinsky_Lilly = 0.15; // this is the standard Smagorinsky-Lilly constant, calculated from Kolmogorov theory: should be 0.1-0.2 //
        double turb_prefactor = 0.25 * All.TurbDiffusion_Coefficient * C_Smagorinsky_Lilly*C_Smagorinsky_Lilly * sqrt(2.0);
        // then scale with inter-particle spacing //
        turb_prefactor *= h_turb*h_turb;
        // calculate frobenius norm of symmetric shear velocity gradient tensor //
        double shear_factor = sqrt((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) *
                                            (SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) +
                                            (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) *
                                            (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) +
                                            (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2]) *
                                            (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
                                   (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] +
                                             SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] +
                                             SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) -
                                            (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] +
                                             SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] +
                                             SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2])));
        // slope-limit and convert to physical units //
        double shearfac_max = 0.5 * sqrt(SphP[i].VelPred[0]*SphP[i].VelPred[0]+SphP[i].VelPred[1]*SphP[i].VelPred[1]+SphP[i].VelPred[2]*SphP[i].VelPred[2]) / h_turb;
        shear_factor = DMIN(shear_factor , shearfac_max * All.cf_atime) * All.cf_a2inv; // physical
#ifdef TURB_DIFF_DYNAMIC
        int u, v; double trace = 0;
        shearfac_max = 0.5 * sqrt(SphP[i].Velocity_bar[0] * SphP[i].Velocity_bar[0] + SphP[i].Velocity_bar[1] * SphP[i].Velocity_bar[1]+SphP[i].Velocity_bar[2] * SphP[i].Velocity_bar[2]) * All.cf_atime / h_turb;
        for (u = 0; u < 3; u++) {
            for (v = 0; v < 3; v++) {
                if (SphP[i].VelShear_bar[u][v] < 0) {SphP[i].VelShear_bar[u][v] = DMAX(SphP[i].VelShear_bar[u][v], -shearfac_max);}
                else {SphP[i].VelShear_bar[u][v] = DMIN(SphP[i].VelShear_bar[u][v], shearfac_max);}
                if (u == v) {trace += SphP[i].VelShear_bar[u][u];}}}
        /* If it was already trace-free, don't zero out the diagonal components */
        if (trace != 0 && NUMDIMS > 1) {for (u = 0; u < NUMDIMS; u++) {SphP[i].VelShear_bar[u][u] -= 1.0 / NUMDIMS * trace;}}
        for (u = 0; u < 3; u++) { /* Don't want to recalculate these a bunch later on, so save them */
            SphP[i].Velocity_hat[u] *= All.TurbDynamicDiffSmoothing;
            for (v = 0; v < 3; v++) {SphP[i].MagShear_bar += SphP[i].VelShear_bar[u][v] * SphP[i].VelShear_bar[u][v];}}
        SphP[i].MagShear = sqrt(2.0) * shear_factor / All.cf_a2inv; // Don't want this physical
        SphP[i].MagShear_bar = DMIN(sqrt(2.0 * SphP[i].MagShear_bar), shearfac_max); turb_prefactor /= 0.25;
#endif
        // ok, combine to get the diffusion coefficient //
        SphP[i].TD_DiffCoeff = turb_prefactor * shear_factor; // physical
    } else {
        SphP[i].TD_DiffCoeff = 0;
    }
#ifdef TURB_DIFF_ENERGY
#if !defined(CONDUCTION_SPITZER)
    SphP[i].Kappa_Conduction = 0;
#endif
    SphP[i].Kappa_Conduction += All.ConductionCoeff * SphP[i].TD_DiffCoeff * SphP[i].Density * All.cf_a3inv; // physical
#endif
#ifdef TURB_DIFF_VELOCITY
#if !defined(VISCOSITY_BRAGINSKII)
    SphP[i].Eta_ShearViscosity = 0; SphP[i].Zeta_BulkViscosity = 0;
#endif
    SphP[i].Eta_ShearViscosity += All.ShearViscosityCoeff * SphP[i].TD_DiffCoeff * SphP[i].Density * All.cf_a3inv; // physical
    SphP[i].Zeta_BulkViscosity += All.BulkViscosityCoeff * SphP[i].TD_DiffCoeff * SphP[i].Density * All.cf_a3inv; // physical
#endif
}
#endif
