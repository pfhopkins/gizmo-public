#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../../allvars.h"
#include "../../proto.h"

/*! Routines for cosmic ray 'fluid' modules (as opposed to the explicit CR-PIC methods, which are in the grain+particles section of the code)
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */






/* return total CR energy density associated with a cell */
double INLINE_FUNC Get_CosmicRayEnergyDensity_cgs(int i)
{
    if(i<=0) {return 0;}
    return 1.6e-12; // eV/cm-3, approximate from Cummings et al. 2016 V1 data
}


/* return total CR ionization rate zeta_cr in s^-1 associated with a cell */
double Get_CosmicRayIonizationRate_cgs(int i)
{
    return 1.e-5 * Get_CosmicRayEnergyDensity_cgs(i); // scales following Cummings et al. 2016 to 1.6e-17 per eV/cm^3
}


/* cosmic ray heating of gas, from Guo & Oh 2008, following Mannheim & Schlickeiser 1994.
 We assume all the electron losses go into radiation (and ionization), as the electron coulomb losses into gas are lower than protons by factor of energy and me/mp.
 For protons, we assume 1/6 of the hadronic losses (based on branching ratios) and all of the Coulomb losses thermalize.
 Do want to make sure that the rates we assume here are consistent with those used in the CR cooling routine above. */
double CR_gas_heating(int target, double n_elec, double nH0, double nHcgs)
{
#if defined(CRFLUID_ALT_DISABLE_LOSSES)
    return 0;
#endif
    double a_hadronic, b_coulomb_ion_per_GeV, f_heat_hadronic;
    a_hadronic = 6.37e-16; b_coulomb_ion_per_GeV = 3.09e-16*(n_elec + 0.57*nH0)*HYDROGEN_MASSFRAC; f_heat_hadronic=1./6.; /* some coefficients; a_hadronic is the default coefficient, b_coulomb_ion_per_GeV the default divided by GeV, b/c we need to divide the energy per CR. note there is an extra factor in principle for the ionization term here compared to its version in the CR losses module above: this represents the fraction of CR energy going into the thermal energy of the gas, as opposed to ionization energy, but this is close to unity */
#if   defined(COOL_LOW_TEMPERATURES) // no CR module, but low-temperature cooling is on, we account for the CRs as a heating source, assuming a MW-like background scaled cosmologically to avoid over-heating IGM at high redshifts //
    double prefac_CR=1.; if(All.ComovingIntegrationOn) {double rhofac = (PROTONMASS*nHcgs/HYDROGEN_MASSFRAC) / (1000.*COSMIC_BARYON_DENSITY_CGS); if(rhofac < 0.2) {prefac_CR=0;} else {if(rhofac > 200.) {prefac_CR=1;} else {prefac_CR=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off CR heating for any gas with density unless it's >1000 times the cosmic mean density
    return (0.87*f_heat_hadronic*a_hadronic + 0.53*b_coulomb_ion_per_GeV) * (1.6e-12*prefac_CR) / (1.e-2 + nHcgs); // assume MW-like CR background modulated by above factor (1.6e-12*prefac_CR)=eCR_cgs here //
#else
    return 0;
#endif
}


/* subroutine to calculate which part of the adiabatic PdV work from the RP gets assigned to the CRs vs the gas; since the CRs are always smooth by definition under this operation this follows simply from the local cell divergence and the effective CR eos */
double CR_calculate_adiabatic_gasCR_exchange_term(int i, double dt_entr, double gamma_minus_eCR_tmp, int mode)
{
    double u0, d_CR; if(mode==0) {u0=SphP[i].InternalEnergy;} else {u0=SphP[i].InternalEnergyPred;} // initial energy
    if(u0<All.MinEgySpec) {u0=All.MinEgySpec;} // enforced throughout code
    
    double divv_p=-dt_entr*P[i].Particle_DivVel*All.cf_a2inv, divv_f=divv_p, divv_u=0; // get locally-estimated gas velocity divergence for cells - if using non-Lagrangian method, need to modify. take negative of this [for sign of change to energy] and multiply by timestep
    if(All.ComovingIntegrationOn) {double divv_h=-dt_entr*(3.*All.cf_hubble_a); divv_p+=divv_h; divv_f+=divv_h;} // include hubble-flow terms
    double P_cr = gamma_minus_eCR_tmp * SphP[i].Density * All.cf_a3inv / P[i].Mass, P_tot = SphP[i].Pressure * All.cf_a3inv; // define the pressure from CRs and total pressure (physical units)
#ifdef MAGNETIC
    double B2=0; int k; for(k=0;k<3;k++) {double B=Get_Gas_BField(i,k)*All.cf_a2inv; B2+=B*B;}
    P_tot += 0.5*B2; // add magnetic pressure [B^2/2], in physical code units, since it contributes to the PdV work but not included in 'pressure' total above
#endif
    double fac_P = DMAX(0, DMIN(1, P_cr/(P_tot + 1.e-10*P_cr + MIN_REAL_NUMBER))); // fraction of total pressure from CRs
    double Ui = u0 * P[i].Mass; // factor for multiplication below, and initial thermal energy
    double dtI_hydro = SphP[i].DtInternalEnergy * P[i].Mass * dt_entr; // change given by hydro-step computed delta_InternalEnergy
    double min_IEgy = P[i].Mass * All.MinEgySpec; // minimum internal energy - in total units -
    
    if(divv_p*dtI_hydro > 0 || divv_f*dtI_hydro > 0) // same sign from hydro and from smooth-flow-estimator, suggests we are in a smooth flow, so we'll use stronger assumptions about the effective 'entropy' here
    {
        if(divv_p*dtI_hydro <= 0) {divv_u=divv_f;} // if divv_f agrees in sign here, use it
        if(divv_f*dtI_hydro <= 0) {divv_u=divv_p;} // if divv_p agrees in sign here, use it
        if(divv_p*divv_f > 0) {if(fabs(divv_p) > fabs(divv_f)) {divv_u=divv_p;} else {divv_u=divv_f;}} // if both agree in sign here, use -larger- since more accurately captures CR-dominated limit
        d_CR = gamma_minus_eCR_tmp * divv_u; // expected PdV CR energy change
        if(fabs(d_CR) > fabs(dtI_hydro)) {d_CR = dtI_hydro;} // do not allow this to exceed the sum (since all terms have the same sign here, in a well-ordered smooth flow)
        if(fabs(d_CR) < fac_P*fabs(dtI_hydro)) {d_CR = fac_P*dtI_hydro;} // but also do not allow CR term to be -below- CR pressure fraction times total term, since that should be attributed to the CR (as this is all a quasi-adiabatic term)
    } else { // both divv terms agree with each other, but dis-agree with the sign of the total change. can't assume anything about smoothness-of-the-flow
        if(fabs(divv_p) > fabs(divv_f)) {divv_u=divv_f;} else {divv_u=divv_p;} // pick the divv estimator with the smaller absolute magnitude, since it deviates
        d_CR = gamma_minus_eCR_tmp * divv_u; // expected PdV CR energy change
        double f_limiter, fac_test=fabs(d_CR)/fabs(dtI_hydro); if(fac_test>fac_P) {d_CR*=fac_P/fac_test;} // don't let CR change exceed their pressure fraction
        if(d_CR > 0) {if(Ui <= min_IEgy) {f_limiter = 1.e-20;} else {f_limiter=0.5;} // gas will be 'cooled', limit so don't overshoot when Pcr is large
            if(d_CR > f_limiter*(Ui-min_IEgy)) {d_CR = f_limiter*(Ui-min_IEgy);} // limit fractional loss to gas
        } else {f_limiter = 1000.; if(fabs(d_CR)>f_limiter*Ui) {d_CR=-f_limiter*Ui;}} // gas will be heated, limit fractional gain
    }
    return d_CR; // return final value
}

