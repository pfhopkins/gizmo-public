#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"

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




/*! return the pressure of particle i: this subroutine needs to  set the value of the 'press' variable (pressure), which you can see from the
    templates below can follow an arbitrary equation-of-state. for more general equations-of-state you want to specifically set the soundspeed
    variable as well. */
double get_pressure(int i)
{
    double soundspeed=0, press=0, gamma_eos_index = GAMMA(i); /* get effective adiabatic index */
    press = (gamma_eos_index-1) * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i); /* ideal gas EOS (will get over-written it more complex EOS assumed) */
    
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
    gamma_eos_index=7./5.; double rho=Particle_density_for_energy_i(i), nH_cgs=rho*All.cf_a3inv * ( All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam ) / PROTONMASS;
    if(nH_cgs > 2.30181e16) {gamma_eos_index=5./3.;} /* dissociates to atomic if dense enough (hot) */
#if (EOS_GMC_BAROTROPIC==0) // exact EOS used in Bate Bonnell & Bromm 2003 and related works - isothermal below 6e10 cm^-3, adiabatic above
    if (nH_cgs < 6e10) {press = 5.663e-16 * nH_cgs;} // isothermal below 6e10 cm^-3 (adiabatic gamma=5/3 for soundspeed, etc, but this assumes effective eos from cooling, etc
    else press = 3.4e-5 * pow(nH_cgs/6e10,1.4);
#else    
    if (nH_cgs < 1.49468e8) {press = 6.60677e-16 * nH_cgs;} // isothermal below ~10^8 cm^-3 (adiabatic gamma=5/3 for soundspeed, etc, but this assumes effective eos from cooling, etc
    else if (nH_cgs < 2.30181e11) {press = 1.00585e-16 * pow(nH_cgs, 1.1);} // 'transition' region
    else if (nH_cgs < 2.30181e16) {press = 3.92567e-20 * pow(nH_cgs, gamma_eos_index);} // adiabatic molecular
    else if (nH_cgs < 2.30181e21) {press = 3.1783e-15 * pow(nH_cgs, 1.1);} // 'transition' region
    else {press = 2.49841e-27 * pow(nH_cgs, gamma_eos_index);} // adiabatic atomic
#endif
    press /= All.UnitPressure_in_cgs;
    /* in this case the EOS is modeling cooling, etc, so -does not- allow shocks or deviations from adiabat, so we reset the internal energy every time this is checked */
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = press / (rho * (gamma_eos_index-1.));
#endif
    
    
    
#if defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE)
    /* add an artificial pressure term to suppress fragmentation at/below the explicit resolution scale */
    double h_eff = DMAX(Get_Particle_Size(i), All.ForceSoftening[0]/2.8); /* need to include latter to account for inter-particle spacing << grav soft cases */
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
    if(soundspeed == 0) {SphP[i].SoundSpeed = sqrt(gamma_eos_index * press / Particle_density_for_energy_i(i));} else {SphP[i].SoundSpeed = soundspeed;}
#endif
    return press;
}




/*! this function allows the user to specify an arbitrarily complex adiabatic index. note that for pure adiabatic evolution, one can simply set the pressure to obey some barytropic equation-of-state and use EOS_GENERAL to tell the code to deal with it appropriately.
      but for more general functionality, we want this index here to be appropriately variable. */
double gamma_eos(int i)
{
#ifdef EOS_SUBSTELLAR_ISM
    if(i>=0) {
        if(P[i].Type==0) {
            double T_eff_atomic = 1.23 * (5./3.-1.) * U_TO_TEMP_UNITS * SphP[i].InternalEnergyPred;
            double nH_cgs = SphP[i].Density*All.cf_a3inv * ( All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam ) / PROTONMASS;
            double T_transition=DMIN(8000.,nH_cgs), f_mol=1./(1. + T_eff_atomic*T_eff_atomic/(T_transition*T_transition));
            /* double gamma_mol_atom = (29.-8./(2.-f_mol))/15.; // interpolates between 5/3 (fmol=0) and 7/5 (fmol=1) */
            /* return gamma_mol_atom + (5./3.-gamma_mol_atom) / (1 + T_eff_atomic*T_eff_atomic/(40.*40.)); // interpolates back up to 5/3 when temps fall below ~30K [cant excite upper states] */
            
            /* We take a detailed fit from Vaidya et al. A&A 580, A110 (2015) for n_H ~ 10^7, which accounts for collisional dissociation at 2000K and ionization at 10^4K,
               and take the fmol-weighted average with 5./3 at the end to interpolate between atomic/not self-shielding and molecular/self-shielding. Gamma should technically
               really come from calculating the species number-weighted specific heats, but fmol is very approximate so this should be OK */
            double gamma_mol = 5./3, logT = log10(T_eff_atomic);
            gamma_mol -= 0.381374640 * sigmoid_sqrt(5.946*(logT-1.248)); // going down from 5./3 at 10K to the dip at ~1.2
            gamma_mol += 0.220724233 * sigmoid_sqrt(6.176*(logT-1.889)); // peak at ~ 80K
            gamma_mol -= 0.067922267 * sigmoid_sqrt(10.26*(logT-2.235)); // plateau at ~1.4
            gamma_mol -= 0.418671231 * sigmoid_sqrt(7.714*(logT-3.134)); // collisional dissociation, down to ~1.1
            gamma_mol += 0.6472439052 * sigmoid_sqrt(98.87*(logT-4.277)); // back to 5/3 once we're fully dissociated
            // comment out the above line and uncomment the two lines below if you want the exact version from Vaidya+15, which rolls the heat of ionization into the EOS - note that this should NOT be used with the standard cooling module
//            gamma_mol += 0.659888854 / (1 + (logT-4.277)*(logT-4.277)/0.176); // peak at ~5./3 for atomic H after dissoc but before ionization
//            gamma_mol += 0.6472439052 * sigmoid_sqrt(98.87*(logT-5077)); // ionization at 10^4K (note this happens at logT ~ 5 because we're just adopting a simple conversion factor from u to T
            return gamma_mol*f_mol + (1-f_mol)*5./3;
        }
    }
#endif
    return GAMMA_DEFAULT; /* default to universal constant here */
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



double INLINE_FUNC Particle_density_for_energy_i(int i)
{
#ifdef HYDRO_PRESSURE_SPH
    return SphP[i].EgyWtDensity;
#endif
    return SphP[i].Density;
}



double INLINE_FUNC Particle_effective_soundspeed_i(int i)
{
#ifdef EOS_GENERAL
    return SphP[i].SoundSpeed;
#else
    /* if nothing above triggers, then we resort to good old-fashioned ideal gas */
    return sqrt(GAMMA(i) * SphP[i].Pressure / Particle_density_for_energy_i(i));
#endif
}



/* returns the conversion factor to go -approximately- (for really quick estimation) in code units, from internal energy to soundspeed */
double INLINE_FUNC convert_internalenergy_soundspeed2(int i, double u)
{
    double gamma_eos_touse = GAMMA(i);
    return gamma_eos_touse * (gamma_eos_touse-1) * u;
}





