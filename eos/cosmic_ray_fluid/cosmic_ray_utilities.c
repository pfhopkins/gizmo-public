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

#ifdef COSMIC_RAY_FLUID




/* routine which determines the fraction of injected CR energy per 'bin' of CR energy.
    so it can be used, we send the CR bin, the source type [0=SNe; 1=stellar winds; 2-4=unused now; 5=sink/AGN;],
    the target gas cell, the shock velocity (used to scale if desired),
    and a flag indicating the option to return the spectral slope/index */
double CR_energy_spectrum_injection_fraction(int k_CRegy, int source_type, double shock_vel, int return_index_in_bin, int target)
{
    double f_bin = 1./N_CR_PARTICLE_BINS; /* uniformly distributed */
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
#if (N_CR_PARTICLE_BINS == 2) /* one-bin protons, one electrons */
    double f_bin_v[2]={0.95 , 0.05}; f_bin=f_bin_v[k_CRegy]; // 5% of injection into e-, roughly motivated by observed spectra and nearby SNRs
#endif
#if (N_CR_PARTICLE_BINS > 2) /* multi-bin spectrum for p and e-: inset assumptions about injection spectrum here! */
    double f_elec = 0.02; // fraction of the energy to put into e- as opposed to p+ at injection [early experiments with 'observed'  fraction ~ 1% give lower e-/p+ actually observed in the end, so tentative favoring closer to equal at injection? but not run to z=0, so U_rad high from CMB; still experimenting here]
    double inj_slope = 4.25; // injection slope with j(p) ~ p^(-inj_slope), so dN/dp ~ p^(2-inj_slope)
    double R_break_e = 1.0; // location of spectral break for injection e- spectrum, in GV
    double inj_slope_lowE_e = 4.2; // injection slope with j(p) ~ p^(-inj_slope), so dN/dp ~ p^(2-inj_slope), for electrons below R_break_e
#if !defined(CRFLUID_ALT_RSOL_FORM)
    inj_slope_lowE_e=4.25; // slightly better fit with this scheme, though a bit marginal compared to using the above defaults (also tried R_break_e=2.0; inj_slope=4.30; , similar, but no significant improvement)
#endif
    double R=return_CRbin_CR_rigidity_in_GV(-1,k_CRegy); int species=return_CRbin_CR_species_ID(k_CRegy); // get bin-centered R and species type
    //if(species < 0 && R < R_break_e) {inj_slope = inj_slope_lowE_e;} // follow model injection spectra favored in Strong et al. 2011 (A+A, 534, A54), who argue the low-energy e- injection spectrum must break to a lower slope by ~1 independent of propagation and re-acceleration model
    if(species > -200 && R < R_break_e) {inj_slope = inj_slope_lowE_e;} // follow model injection spectra favored in Strong et al. 2011 (A+A, 534, A54), who argue the low-energy e- injection spectrum must break to a lower slope by ~1 independent of propagation and re-acceleration model
    double EGeV = return_CRbin_kinetic_energy_in_GeV_binvalsNRR(k_CRegy); // get bin-centered E_GeV for normalizing total energy in bin
    f_bin = EGeV * pow(R/R_break_e , 3.-inj_slope) * log(CR_global_max_rigidity_in_bin[k_CRegy] / CR_global_min_rigidity_in_bin[k_CRegy]); // normalize accounting for slope, isotropic spectrum, logarithmic bin width [which can vary], and energy per N

    if(return_index_in_bin) {return 2.-inj_slope;} // this is the index corresponding to our dN/dp ~ p^gamma
    double f_norm = 1.e-20; // default to very little energy
    if(species == -1) {f_norm = f_elec;} // e-
    if(species == +1) {f_norm = 1.-f_elec;} // p
    if(species == -2) {f_norm = 1.e-10 * f_elec;} // e+ (assuming negligible e+ injection to start)
    if(species > 1 && species != 7) // heavy elements need to scale injection rates appropriately
    {
        double Zfac=0, Zfac_ISM=P[target].Metallicity[0]/All.SolarAbundances[0], mu_wt=return_CRbin_CRmass_in_mp(-1,k_CRegy), Z_cr=fabs(return_CRbin_CR_charge_in_e(-1,k_CRegy)), Mism_over_Mej=1; // scale heavier elements to the metallicity of the gas into which CRs are being accelerated
        Zfac = Zfac_ISM; // assume abundance of ejecta is identical to ambient ISM into which its being ejected
        if(source_type == 1) {Zfac = (Mism_over_Mej*Zfac_ISM + 1.4*DMIN(Zfac_ISM,1.))/(Mism_over_Mej*HYDROGEN_MASSFRAC + HYDROGEN_MASSFRAC);} // stellar outflows. using FIRE-3 yields this is exact after IMF-integrating for Z_ism=Z_star, for CNO; basically get slight enhancement, but not much, b/c these are OB winds; even including AGB, would only move factor to 3.4 from 1.4, which is halved, so not much effect at all.
        if(source_type == 0) {if(shock_vel>5000./UNIT_VEL_IN_KMS) {Zfac=(Mism_over_Mej*Zfac_ISM + 9.70)/(Mism_over_Mej*HYDROGEN_MASSFRAC + 0.025);} else {Zfac=(Mism_over_Mej*Zfac_ISM + 13.645)/(Mism_over_Mej*HYDROGEN_MASSFRAC + 0.441);}} // shock_vel here tells us if its a 1a [faster] or CCSNe. for CCSNe, use Iwamoto 1999 and Nomoto 2006 to get abundances of ejecta in CNO, integrated over IMF and species of interest. for both, assume acceleration efficiency is maximized at highest mach numbers after shock actually develops, so swept-up ISM mass is ~ejecta mass
        // now scale to mass fraction for solar abundances which gives the units we work with above
        if(species == 2) {Zfac *= 3.7e-9;} // B (for standard elements initialize to solar ratios assuming similar energy/nucleon)
        if(species == 3) {Zfac *= 2.4e-3;} // C
        if(species == 4) {Zfac *= 1.4e-10;} // Be7+9 (stable)
        if(species == 5) {Zfac *= 1.4e-20;} // Be10 (radioactive)
        if(species == 6) {Zfac *= 0.0094;} // CNO (combined bin)
        f_norm = Zfac * pow(mu_wt/Z_cr , inj_slope-3.) / mu_wt; // approximate injection factor for a constant-beta distribution at a given R_GV needed below
    }
    f_bin *= f_norm; // normalize injection depending on the species (e- or p+, etc)
#endif
#endif
    if(return_index_in_bin) {return 0;}
    return f_bin;
}


/* routine which gives diffusion coefficient as a function of energy for the 'constant diffusion coefficient' models:
    current default: -extremely- simple power-law, assuming diffusion coefficient increases with CR energy per unit charge as (E/Z)^(1/2) */
double diffusion_coefficient_constant(int target, int k_CRegy)
{
    double dimensionless_kappa_relative_to_GV_protons = 1;
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
#ifdef CRFLUID_DIFFUSION_CORRECTION_TERMS
    int target_bin_centering_for_CR_quantities = -1; // the correction terms depend on these being evaluated at their bin-centered locations
    dimensionless_kappa_relative_to_GV_protons = return_CRbin_beta_factor(target_bin_centering_for_CR_quantities,k_CRegy) * pow( CR_global_min_rigidity_in_bin[k_CRegy]*CR_global_max_rigidity_in_bin[k_CRegy] , 0.5 * 0.6 ); // assume a quasi-empirical scaling here, and for these correction terms its important that the 'bin center' being used for the zero point here is the geometric mean of the bin edges, hence the 0.5 term b/c geometric mean is sqrt[min*max] //
#else
    int target_bin_centering_for_CR_quantities = target; // if this = target, evaluate quantities like R_GV at the CR-energy weighted mean of the bin, if =-1, evaluate them at the bin center instead: important for some subtle effects especially if using numerical derivatives for correction terms
    dimensionless_kappa_relative_to_GV_protons = return_CRbin_beta_factor(target_bin_centering_for_CR_quantities,k_CRegy) * pow( return_CRbin_CR_rigidity_in_GV(-1,k_CRegy) , 0.6 ); // assume a quasi-empirical scaling here //
#endif
#endif
    return All.CosmicRayDiffusionCoeff * dimensionless_kappa_relative_to_GV_protons;
}


                                                                                                                              
/* routine which gives diffusion coefficient as a function of CR bin for the self-confinement models [in local equilibrium]. mode sets what we assume about the 'sub-grid'
    parameters f_QLT (rescales quasi-linear theory) or f_cas (rescales turbulence strength)
      <=0: fQLT=1 [most naive quasi-linear theory, ruled out by observations],  fcas=1 [standard Goldreich-Shridar cascade]
        1: fQLT=100, fcas=1
        2: fQLT=1, fcas=100
        3: fQLT=1, fcas-K41 from Hopkins et al. 2020 paper, for pure-Kolmogorov isotropic spectrum
        4: fQLT=1, fcas-IK, IK spectrum instead of GS
   if set mode < 0, will also ignore the dust-damping contribution from Squire et al. 2020.
   coefficient is returned in cgs units
 */
#ifndef CRFLUID_SET_SC_MODEL
#define CRFLUID_SET_SC_MODEL 1 /* set which mode to return from the SC subroutine here, of the various choices for how to e.g. model fCas, fQLT */
#endif
double diffusion_coefficient_self_confinement(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG,
    double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion)
{
    double vol_inv = SphP[target].Density*All.cf_a3inv / P[target].Mass, fturb_multiplier=1, f_QLT=1, R_CR_GV, Z_charge_CR, M_cr_mp, b0[3]={0}, p0[3]={0};
    int target_bin_centering_for_CR_quantities = target; // if this = target, evaluate quantities like R_GV at the CR-energy weighted mean of the bin, if =-1, evaluate them at the bin center instead: important for some subtle effects especially if using numerical derivatives for correction terms
#ifdef CRFLUID_DIFFUSION_CORRECTION_TERMS
    target_bin_centering_for_CR_quantities = -1; // the correction terms depend on these being evaluated at their bin-centered locations
#endif
    R_CR_GV=return_CRbin_CR_rigidity_in_GV(target_bin_centering_for_CR_quantities,k_CRegy); Z_charge_CR=return_CRbin_CR_charge_in_e(target,k_CRegy); M_cr_mp=return_CRbin_CRmass_in_mp(target,k_CRegy);
    int k; double n_cgs=rho_cgs/PROTONMASS_CGS, EPSILON_SMALL=1.e-50, e_CR=0, e_B=0, bhat_dot_CR_Pgrad=0, B2=0;
#ifdef MAGNETIC
    for(k=0;k<3;k++) {b0[k]=SphP[target].BPred[k]*vol_inv*All.cf_a2inv; B2+=b0[k]*b0[k];}
    e_B=0.5*B2; B2=1./sqrt(B2+MIN_REAL_NUMBER); for(k=0;k<3;k++) {b0[k]*=B2;} // calculate B-field energy and bhat vector
#else
    for(k=0;k<3;k++) {b0[k]=SphP[target].Gradients.CosmicRayPressure[k_CRegy][k]; B2+=b0[k]*b0[k];} // just assume equilibrium here, convert to physical units, pick arbitrary direction here //
    e_B=SphP[target].Pressure*All.cf_a3inv; for(k=0;k<3;k++) {b0[k]/=sqrt(B2+MIN_REAL_NUMBER);} // calculate B-field energy and bhat vector
#endif
    e_CR=SphP[target].CosmicRayEnergyPred[k_CRegy]*vol_inv; for(k=0;k<3;k++) {p0[k]=SphP[target].Gradients.CosmicRayPressure[k_CRegy][k]*All.cf_a3inv/All.cf_atime;}
    for(k=0;k<3;k++) {bhat_dot_CR_Pgrad += b0[k]*p0[k];} // dot product of bhat and CR pressure gradient, summed over relevant bins
    double beta=return_CRbin_beta_factor(target_bin_centering_for_CR_quantities,k_CRegy), Omega_gyro=beta*(0.00898734*b_muG/R_CR_GV) * UNIT_TIME_IN_CGS, r_L=beta*C_LIGHT_CODE/Omega_gyro, kappa_0=r_L*beta*C_LIGHT_CODE; /* all in physical -code- units */
    double x_LL = DMAX( r_L / L_scale, EPSILON_SMALL ), vA_code=Get_Gas_ion_Alfven_speed_i(target), k_turb=1./L_scale, k_L=1./r_L;

    if(mode==1) {f_QLT = 100;} // multiplier to account for arbitrary deviation from QLT, applies to all damping mechanisms [100 = favored value in our study; or could use fcas = 100]
    fturb_multiplier = pow(M_A,3./2.); // multiplier to account for different turbulent cascade models (fcas = 1)
    if(mode==2) {fturb_multiplier *= 100.;} // arbitrary multiplier (fcas = 100, here)
    if(mode==3) {fturb_multiplier = pow(M_A,3./2.) * 1./(pow(M_A,1./2.)*pow(x_LL,1./6.));} // pure-Kolmogorov (fcas-K41)
    if(mode==4) {fturb_multiplier = pow(M_A,3./2.) / pow(x_LL,1./10.);} // GS anisotropic but perp cascade is IK (fcas-IK) /

    /* ok now we finally have all the terms needed to calculate the various damping rates that determine the equilibrium diffusivity */
    double U0bar_grain=3., rhograin_int_cgs=1., fac_grain=R_CR_GV*sqrt(n_cgs)*U0bar_grain/(b_muG*rhograin_int_cgs), f_grainsize = DMAX(8.e-4*pow(fac_grain*(temperature/1.e4),0.25), 3.e-3*sqrt(fac_grain)), Z_sol=1.; // b=2, uniform logarithmic grain spectrum over a factor of ~100 in grain size; f_grainsize = 0.07*pow(sqrt(fion*n1)*EcrGeV*T4/BmuG,0.25); // MRN size spectrum
#ifdef METALS
    Z_sol = P[target].Metallicity[0]/0.014;
#endif
    double G_dust = vA_code*k_L * Z_sol * f_grainsize; // also can increase by up to a factor of 2 for regimes where charge collisionally saturated, though this is unlikely to be realized
    if(mode<0) {G_dust = 0;} // for this choice, neglect the dust-damping term 
    double G_ion_neutral = (5.77e-11 * n_cgs * (0.97*nh0 + 0.03*nHe0) * sqrt(temperature)) * UNIT_TIME_IN_CGS / sqrt(M_cr_mp); // ion-neutral damping: need to get thermodynamic quantities [neutral fraction, temperature in Kelvin] to compute here -- // G_ion_neutral = (xiH + xiHe); // xiH = nH * siH * sqrt[(32/9pi) *kB*T*mH/(mi*(mi+mH))]
    double fac_turb = sqrt(k_turb*k_L) * fturb_multiplier; // factor to use below
    double G_turb_plus_linear_landau = (vA_noion + sqrt(M_PI)*cs_thermal/4.) * fac_turb; // linear Landau + turbulent (both have same form, assume k_turb from cascade above)

#ifdef CRFLUID_ALT_FLUX_FORM_JOCH
    double G0 = G_ion_neutral + G_turb_plus_linear_landau + G_dust; // linear terms all add into single G0 term
    double CRPressureGradScaleLength=Get_CosmicRayGradientLength(target,k_CRegy), x_EB_ECR=(e_B+EPSILON_SMALL)/(e_CR+EPSILON_SMALL); // get scale length used below
    double Gamma_effective = G0, phi_0 = fabs(bhat_dot_CR_Pgrad) * ((sqrt(M_PI)*cs_thermal*vA_code*k_L)/(2.*e_B*G0*G0 + EPSILON_SMALL)); // parameter which determines whether NLL dominates
    if(isfinite(phi_0) && (phi_0>0.01)) {Gamma_effective *= phi_0/(2.*(sqrt(1.+phi_0)-1.));} // this accounts exactly for the steady-state solution for the Thomas+Pfrommer formulation, including both the linear [Landau,turbulent,ion-neutral] + non-linear terms. can estimate (G_nonlinear_landau_effective = Gamma_effective - G0)
    /* with damping rates above, equilibrium transport is equivalent to pure streaming, with v_stream = vA + (diffusive equilibrium part) give by the solution below, proportional to Gamma_effective and valid to O(v^2/c^2) */
    double v_st_eff = vA_code * (1. + f_QLT * 4. * kappa_0 * Gamma_effective * x_EB_ECR * (1. + 2.*vA_code*vA_code/(C_LIGHT_CODE*C_LIGHT_CODE)) / (M_PI*vA_code*vA_code + EPSILON_SMALL)); // effective equilibrium streaming speed for all terms accounted
    return (GAMMA_COSMICRAY(k_CRegy) * v_st_eff * CRPressureGradScaleLength) * (UNIT_VEL_IN_CGS*UNIT_LENGTH_IN_CGS); // convert to effective diffusivity from 'streaming speed' [this introduces the gamma and gradient length], and convert to CGS from code units
#else
    double G_adiabatic = 0.5*P[target].Particle_DivVel*All.cf_a2inv; // adiabatic term [signed like the other linear terms
    double Gamma_NLL = (sqrt(M_PI)/8.) * cs_thermal / r_L; // NLL prefactor: Gamma_NLL is this times (e_A/e_B)
    double f_cas_ET = 7. * (vA_noion / (beta * C_LIGHT_CODE)) * log(1.+L_scale/r_L); /* Alfvenic turbulence following an anisotropic Goldreich-Shridar cascade, per Chandran 2000 */
    double S_ext_turb = f_cas_ET * vA_noion * fac_turb * M_A*M_A * (r_L/L_scale); // extrinsic turbulence cascade term. note consistency means multiplying by pitch-angle and gyro-averaging factors, and fcas above; expression here assumes whatever you do its a balanced cascade (defauly to GS95 scalings)
    double S_ext_gri  = vA_code * fabs(bhat_dot_CR_Pgrad) / (e_B + MIN_REAL_NUMBER); // Flux-steady-state value of the GRI term, normalized to e_B. note unless we rewrite this to the 5th-order polynomial version, assumption here is that return_CRbin_nuplusminus_asymmetry(i,k_CRegy) -> 1 for the term here in steady state
    double Gamma_LIN = -(G_ion_neutral + G_turb_plus_linear_landau + G_dust + G_adiabatic); // sum of all the linear damping/growth terms
    double S_ext = S_ext_turb + S_ext_gri; // total driving term, for the flux-steady assumption

    if(mode==5) {S_ext=S_ext_turb + S_ext_gri; Gamma_LIN=-DMAX(DMIN(fturb_multiplier,1.),100.)*vA_noion*k_turb*(0.+1.*pow(k_L/k_turb,0.25))*((e_CR+EPSILON_SMALL)/(e_B+EPSILON_SMALL));} // resolve fundamental issues with SC+ET models by invoking alternative damping, following Hopkins et al. 2021
    if(mode==6) {double S_lin = 9.0e-19*UNIT_LENGTH_IN_CGS * (1.+M_A) * sqrt(vA_noion*vA_noion + cs_thermal*cs_thermal) * pow(r_L*UNIT_LENGTH_IN_CGS/1.5e12 , -0.66);
        S_ext = 0.1*S_ext_turb + 0.01*S_ext_gri; Gamma_LIN=S_lin - (G_ion_neutral+G_adiabatic+1.e-10*G_dust); Gamma_NLL += vA_noion/r_L;} // resolve fundamental issues with SC+ET models by invoking alternative linear driving, following Hopkins et al. 2021
    if(mode==7) {f_cas_ET=vA_noion/(0.007 * C_LIGHT_CODE); S_ext_turb=f_cas_ET*vA_noion*fac_turb*M_A*M_A*pow(r_L/L_scale,2./3.); S_ext=S_ext_turb+S_ext_gri;} // resolve fundamental issues with SC+ET models by invoking alternative constant driving, following Hopkins et al. 2021

    double fac=0, f0 = Gamma_LIN/(2.*Gamma_NLL + MIN_REAL_NUMBER), f1 = 4.*Gamma_NLL*S_ext / (Gamma_LIN*Gamma_LIN + MIN_REAL_NUMBER);
    if(f0>0) {fac=f0*(1.+sqrt(1.+f1));} else {if(f1>0.1) {fac=f0*(1.-sqrt(1.+f1));} else {fac=(S_ext/(fabs(Gamma_LIN)+MIN_REAL_NUMBER))*(1.-f1/4.);}}
    double gyro_avg_factor = 3./4.; // weighting factor from pitch-angle averaging over nu, gives ~3/4 for nu~|mu-vA/c|^2, etc.
    return (kappa_0 / fac) * (4./(3.*M_PI*gyro_avg_factor)) * (UNIT_VEL_IN_CGS*UNIT_LENGTH_IN_CGS);
#endif
}


                                                                                                                              
/* routine which gives diffusion coefficient [in cgs] for extrinsic turbulence models. 'mode' sets whether we assume Alfven modes (mode<0), Fast-mode scattering (mode>0), or both (=0),
     0: 'default' Alfven + Fast modes (both, summing scattering rates linearly)
    -1: 'default' Alfven modes: correctly accounting for an anisotropic Goldreich-Shridar cascade, per Chandran 2000
    -2: Alfven modes in pure Goldreich-Shridhar cascade, ignoring anisotropic effects [*much* higher scattering rate, artificially]
     1: 'default' Fast modes: following Yan & Lazarian 2002, accounting for damping from viscous, ion-neutral, and other effects, and suppression if beta>1
     2: Fast modes following a pure isotropic Kolmogorov cascade down to gyro radius [*much* higher scattering rate, artificially]
 */
double diffusion_coefficient_extrinsic_turbulence(int mode, int target, int k_CRegy, double M_A, double L_scale, double b_muG,
    double vA_noion, double rho_cgs, double temperature, double cs_thermal, double nh0, double nHe0, double f_ion)
{
    double f_cas_ET=MAX_REAL_NUMBER, EPSILON_SMALL=1.e-50, h0_kpc=L_scale*UNIT_LENGTH_IN_KPC;
    if(mode <= 0) /* Alfvenic turbulence [default here following Chandran 200, including anisotropy effects] */
    {
        f_cas_ET = 0.007 * C_LIGHT_CODE / vA_noion; /* damping in Alfvenic turbulence following an anisotropic Goldreich-Shridar cascade, per Chandran 2000 */
        if(mode==-2) {f_cas_ET = 1;} /* pure Goldreich-Shridhar cascade, ignoring anisotropic effects per Chandran-00 */ //if(M_A < 1.) {f_cas_ET = pow(1./DMAX(M_A*M_A , x_LL), 1./3.);} /* Lazarian '16 modification for weak cascade in sub-Alfvenic turbulence */
    }
    if(mode >= 0) /* Fast modes [default here following Yan & Lazarian 2002, including damping effects] */
    {
        int target_bin_centering_for_CR_quantities = target; // if this = target, evaluate quantities like R_GV at the CR-energy weighted mean of the bin, if =-1, evaluate them at the bin center instead: important for some subtle effects especially if using numerical derivatives for correction terms
#ifdef CRFLUID_DIFFUSION_CORRECTION_TERMS
        target_bin_centering_for_CR_quantities = -1; // the correction terms depend on these being evaluated at their bin-centered locations
#endif
        double R_CR_GV=return_CRbin_CR_rigidity_in_GV(target_bin_centering_for_CR_quantities,k_CRegy);
        double n1=rho_cgs/PROTONMASS_CGS, T4=temperature/1.e4, fcasET_colless = 0.04*cs_thermal/vA_noion; /* collisionless [Landau] damping of fast modes */
        double fcasET_viscBrg = 0.03*pow(EPSILON_SMALL + M_A,4./3.)*T4/pow(EPSILON_SMALL + b_muG*h0_kpc*n1*R_CR_GV*T4,1./6.); /* Spitzer/Braginski viscous damping of fast modes */
        double fcasET_viscMol = 0.41*pow(EPSILON_SMALL + M_A,4./3.)*nh0/pow(EPSILON_SMALL + b_muG*h0_kpc*n1*R_CR_GV/(EPSILON_SMALL + T4),1./6.); /* atomic/molecular collisional damping of fast modes */
        double f_cas_ET_fast = fcasET_colless + fcasET_viscBrg + fcasET_viscMol; /* fast modes, accounting for damping, following Yan+Lazarian 2005 */
        double fast_gyrores_dampingsuppression = 1.; // term to account for the fact that small pitch angles become unscattered when neutral fraction is large or beta >~ 1, making kappa blow up rapidly */
        double beta_half = cs_thermal / vA_noion; if(beta_half > 1.) {fast_gyrores_dampingsuppression=0;} else {fast_gyrores_dampingsuppression*=exp(-beta_half*beta_half*beta_half);} // parallel modes strongly damped if beta >~ 1
        double f_neutral_crit = 0.001 * pow(T4,0.25) / (pow(n1*beta_half*beta_half,0.75) * sqrt(h0_kpc)); // neutral fraction above which the parallel modes are strongly damped
        if(nh0 > 2.*f_neutral_crit) {fast_gyrores_dampingsuppression=0;} else {fast_gyrores_dampingsuppression*=exp(-(nh0*nh0*nh0*nh0)/(f_neutral_crit*f_neutral_crit*f_neutral_crit*f_neutral_crit));} // suppression very rapid, as exp(-[fn/f0]^4)
        if(mode==2) {fast_gyrores_dampingsuppression=0; f_cas_ET_fast = 0.0009 * pow(R_CR_GV*h0_kpc*h0_kpc/b_muG,1./3.)/(M_A*M_A);} /* kappa~9e28 * (l_alfven/kpc)^(2/3) * RGV^(1/3) * (B/muG)^(-1/3),  follows Jokipii 1966, with our corrections for spectral shape */
        if(mode==3) {fast_gyrores_dampingsuppression=0; f_cas_ET_fast = 0.003 * (h0_kpc + 0.1*pow(R_CR_GV*h0_kpc*h0_kpc/b_muG,1./3.)/(M_A*M_A) + 2.4e-7*R_CR_GV/b_muG);} /* Snodin et al. 2016 -- different expression for extrinsic MHD-turb diffusivity, using the proper definition of L_Alfven for scaling to the correct limit */
        f_cas_ET = 1./(EPSILON_SMALL + 1./(EPSILON_SMALL+f_cas_ET) + fast_gyrores_dampingsuppression / (EPSILON_SMALL+f_cas_ET_fast)); /* combine fast-mode and Alfvenic scattering */
    }
    return 1.e32 * h0_kpc / (EPSILON_SMALL + M_A*M_A) * f_cas_ET;
}

                                                                                                                              

/*!----------------------------------------------------------------------------------------------------------------------------------------------------
 routines below are more general and/or numerical: they generally do NOT need to be modified even if you are changing the
   physical assumptions, energies, or other properties of the CRs
 ----------------------------------------------------------------------------------------------------------------------------------------------------*/


/* cosmic ray interactions affecting the -thermal- temperature of the gas are included in the actual cooling/heating functions;
    they are solved implicitly above. however we need to account for energy losses of the actual cosmic ray fluid, here. The
    timescale for this is reasonably long, so we can treat it semi-explicitly, as we do here.
    -- We use the estimate for combined hadronic + Coulomb losses from Volk 1996, Ensslin 1997, as updated in Guo & Oh 2008: */
void CR_cooling_and_losses(int target, double n_elec, double nHcgs, double dtime_cgs)
{

    if(dtime_cgs <= 0) {return;} /* catch */
    int k_CRegy; double f_ion=DMAX(DMIN(Get_Gas_Ionized_Fraction(target),1.),0.);
    double a_hadronic = 6.37e-16, b_coulomb_ion_per_GeV = 3.09e-16*(n_elec + 0.57*(1.-f_ion))*HYDROGEN_MASSFRAC; /* some coefficients; a_hadronic is the default coefficient, b_coulomb_ion_per_GeV the default Coulomb+ionization (the two scale nearly-identically) normalization divided by GeV, b/c we need to divide the energy per CR  */
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double CR_coolrate,Z,species_ID; CR_coolrate=0; Z=fabs(return_CRbin_CR_charge_in_e(target,k_CRegy)); species_ID=return_CRbin_CR_species_ID(k_CRegy);
        if(species_ID > 0) /* protons and nuclei here */
        {
#if (N_CR_PARTICLE_BINS > 2) /* note these are currently energy-loss expressions; for truly multi-bin, probably better to work with dp/dt, instead of dE/dt */
            double E_GeV=return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), beta=return_CRbin_beta_factor(target,k_CRegy);
            CR_coolrate += b_coulomb_ion_per_GeV * ((Z*Z)/(beta*E_GeV)) * nHcgs; // all protons Coulomb-interact, can be rapid for low-E
            if(E_GeV>=0.28) {CR_coolrate += a_hadronic * nHcgs;} // only GeV CRs or higher trigger above threshold for collisions
#else
            CR_coolrate = (0.87*a_hadronic + 0.53*b_coulomb_ion_per_GeV) * nHcgs; /* for N<=2, assume a universal spectral shape, the factor here corrects for the fraction above-threshold for hadronic interactions, and 0.53 likewise for averaging  */
#endif
        } else { /* electrons here: note for electrons and positrons, always in the relativistic limit, don't need to worry about beta << 1 limits */
            /* bremsstrahlung [following Blumenthal & Gould, 1970]: dEkin/dt=4*alpha_finestruct*r_classical_elec^2*c * SUM[n_Z,ion * Z * (Z+1) * (ln[2*gamma_elec]-1/3) * E_kin */
            double E_GeV=return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), E_rest=0.000511, gamma=1.+E_GeV/E_rest;
            CR_coolrate += n_elec * nHcgs * 1.39e-16 * DMAX(log(2.*gamma)-0.33,0);
            /* synchrotron and inverse compton scale as dE/dt=(4/3)*sigma_Thompson*c*gamma_elec^2*(U_mag+U_rad), where U_mag and U_rad are the magnetic and radiation energy densities, respectively. Ignoring Klein-Nishina corrections here, as they are negligible at <40 GeV and only a ~15% correction up to ~1e5 GeV */
            double b_muG = get_cell_Bfield_in_microGauss(target), U_mag_ev=0.0248342*b_muG*b_muG, U_rad_ev = get_cell_Urad_in_eVcm3(target);
            CR_coolrate += 5.2e-20 * gamma * (U_mag_ev + U_rad_ev); // U_mag_ev=(B^2/8pi)/(eV/cm^(-3)), here; U_rad=U_rad/(eV/cm^-3) //
        }
        
        /* here, cooling is being treated as energy loss 'within the bin'. with denser bins, should allow for movement -between- bins, using the spectral method below */
        CR_coolrate *= CosmicRayFluid_RSOL_Corrfac(k_CRegy); // account for RSOL terms as needed
        double q_CR_cool = exp(-CR_coolrate * dtime_cgs); if(CR_coolrate * dtime_cgs > 20.) {q_CR_cool = 0;}
        SphP[target].CosmicRayEnergyPred[k_CRegy] *= q_CR_cool; SphP[target].CosmicRayEnergy[k_CRegy] *= q_CR_cool;
#ifdef CRFLUID_M1
        int k; for(k=0;k<3;k++) {SphP[target].CosmicRayFlux[k_CRegy][k] *= q_CR_cool; SphP[target].CosmicRayFluxPred[k_CRegy][k] *= q_CR_cool;}
#endif
    }
    return;
}

                                                                                                                              

/* utility to estimate -locally- (without multi-pass filtering) the local Alfven Mach number */
double Get_AlfvenMachNumber_Local(int i, double vA_idealMHD_codeunits, int use_shear_corrected_vturb_flag)
{
    int i1,i2; double v2_t=0,dv2_t=0,b2_t=0,db2_t=0,M_A,h0,EPSILON_SMALL=1.e-50; // factor which will represent which cascade model we are going to use
    for(i1=0;i1<3;i1++)
    {
        v2_t += SphP[i].VelPred[i1]*SphP[i].VelPred[i1];
        for(i2=0;i2<3;i2++) {dv2_t += SphP[i].Gradients.Velocity[i1][i2]*SphP[i].Gradients.Velocity[i1][i2];}
#ifdef MAGNETIC
        b2_t += Get_Gas_BField(i,i1) * Get_Gas_BField(i,i1);
        for(i2=0;i2<3;i2++) {db2_t += SphP[i].Gradients.B[i1][i2]*SphP[i].Gradients.B[i1][i2];}
#endif
    }
    v2_t=sqrt(v2_t); b2_t=sqrt(b2_t); dv2_t=sqrt(dv2_t); db2_t=sqrt(db2_t); dv2_t/=All.cf_atime; db2_t/=All.cf_atime; b2_t*=All.cf_a2inv; db2_t*=All.cf_a2inv; v2_t/=All.cf_atime; dv2_t/=All.cf_atime;
    h0=Get_Particle_Size(i)*All.cf_atime; // physical units

    if(use_shear_corrected_vturb_flag == 1)
    {
        dv2_t = sqrt((1./2.)*((SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) *
            (SphP[i].Gradients.Velocity[1][0]+SphP[i].Gradients.Velocity[0][1]) + (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) *
            (SphP[i].Gradients.Velocity[2][0]+SphP[i].Gradients.Velocity[0][2]) + (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2]) * (SphP[i].Gradients.Velocity[2][1]+SphP[i].Gradients.Velocity[1][2])) +
            (2./3.)*((SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[0][0] + SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[1][1] +
            SphP[i].Gradients.Velocity[2][2]*SphP[i].Gradients.Velocity[2][2]) - (SphP[i].Gradients.Velocity[1][1]*SphP[i].Gradients.Velocity[2][2] + SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[1][1] +
            SphP[i].Gradients.Velocity[0][0]*SphP[i].Gradients.Velocity[2][2]))) * All.cf_a2inv;
#ifdef MAGNETIC
        db2_t = sqrt((1./2.)*((SphP[i].Gradients.B[1][0]+SphP[i].Gradients.B[0][1]) * (SphP[i].Gradients.B[1][0]+SphP[i].Gradients.B[0][1]) +
            (SphP[i].Gradients.B[2][0]+SphP[i].Gradients.B[0][2]) * (SphP[i].Gradients.B[2][0]+SphP[i].Gradients.B[0][2]) +
            (SphP[i].Gradients.B[2][1]+SphP[i].Gradients.B[1][2]) * (SphP[i].Gradients.B[2][1]+SphP[i].Gradients.B[1][2])) +
            (2./3.)*((SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[0][0] + SphP[i].Gradients.B[1][1]*SphP[i].Gradients.B[1][1] +
            SphP[i].Gradients.B[2][2]*SphP[i].Gradients.B[2][2]) - (SphP[i].Gradients.B[1][1]*SphP[i].Gradients.B[2][2] +
            SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[1][1] + SphP[i].Gradients.B[0][0]*SphP[i].Gradients.B[2][2]))) * All.cf_a3inv;
#endif
        double db_v_equiv = h0 * db2_t * vA_idealMHD_codeunits / (EPSILON_SMALL + b2_t); /* effective delta-v corresponding to delta-B fluctuations */
        double dv_e = sqrt((h0*dv2_t)*(h0*dv2_t) + db_v_equiv*db_v_equiv + EPSILON_SMALL); /* total effective delta-velocity */
        double vA_eff = sqrt(vA_idealMHD_codeunits*vA_idealMHD_codeunits + db_v_equiv*db_v_equiv + EPSILON_SMALL); /* effective Alfven speed including fluctuation-B */
        M_A = (EPSILON_SMALL + dv_e) / (EPSILON_SMALL + vA_eff);
    } else {
        M_A = h0*(EPSILON_SMALL + dv2_t) / (EPSILON_SMALL + vA_idealMHD_codeunits); /* velocity fluctuation-inferred Mach number */
        M_A = DMAX(M_A , h0*(EPSILON_SMALL + db2_t) / (EPSILON_SMALL + b2_t)); /* B-field fluctuation-inferred Mach number [in incompressible B-turb, this will 'catch' where dv is locally low instanteously] */
    }
    M_A = DMAX( EPSILON_SMALL , M_A ); // proper calculation of the local Alfven Mach number
    return M_A;
}

                                                                                                                           


/* parent routine to assign diffusion coefficients. for the most relevant physical models, we do a lot of utility here but do the more interesting
    (and uncertain) physical calculation in the relevant sub-routines above, so you don't need to modify all of this in most cases */
void CalculateAndAssign_CosmicRay_DiffusionAndStreamingCoefficients(int i)
{
    /* first define some very general variables, and calculate some useful quantities that will be used for any model */
    int k_CRegy; double DiffusionCoeff, CR_kappa_streaming, CRPressureGradScaleLength, v_streaming;
#if (CRFLUID_DIFFUSION_MODEL > 0)
    double cs_thermal,M_A,L_scale,vA_code,vA_noion,gizmo2gauss,Omega_per_GeV_ifveqc,Bmag,unit_kappa_code,b_muG,E_B,f_ion,temperature,EPSILON_SMALL; int k; k=0;
    unit_kappa_code=UNIT_VEL_IN_CGS*UNIT_LENGTH_IN_CGS; gizmo2gauss=UNIT_B_IN_GAUSS; f_ion=1; temperature=0; EPSILON_SMALL=1.e-50;
    Bmag=2.*SphP[i].Pressure*All.cf_a3inv; cs_thermal=sqrt(convert_internalenergy_soundspeed2(i,SphP[i].InternalEnergyPred)); /* quick thermal pressure properties (we'll assume beta=1 if MHD not enabled) */
#ifdef MAGNETIC /* get actual B-field */
    double B[3]={0}; Bmag=0; for(k=0;k<3;k++) {B[k]=Get_Gas_BField(i,k)*All.cf_a2inv; Bmag+=B[k]*B[k];} // B-field in code units (physical)
#endif
    Bmag=sqrt(DMAX(Bmag,0)); b_muG=Bmag*gizmo2gauss/1.e-6; b_muG=sqrt(b_muG*b_muG + 1.e-6); vA_code=sqrt(Bmag*Bmag/(SphP[i].Density*All.cf_a3inv)); vA_noion=vA_code; E_B=0.5*Bmag*Bmag*(P[i].Mass/(SphP[i].Density*All.cf_a3inv));
    Omega_per_GeV_ifveqc=(0.00898734*b_muG) * UNIT_TIME_IN_CGS; /* B-field in units of physical microGauss; set a floor at nanoGauss level. convert to physical code units */
#ifdef COOLING
    double ne=1, nh0=0, nHe0=0, nHepp, nhp, nHeII, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, rho_cgs, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); rho_cgs=rho*UNIT_DENSITY_IN_CGS; // get thermodynamic properties
    f_ion = DMIN(DMAX(DMAX(DMAX(1-nh0, nhp), ne/1.2), 1.e-8), 1.); // account for different measures above (assuming primordial composition)
#endif
    M_A = Get_AlfvenMachNumber_Local(i,vA_noion,0); /* get turbulent Alfven Mach number estimate. 0 or 1 to turn on shear-correction */
    L_scale = Get_Particle_Size(i)*All.cf_atime; /* define turbulent scales [estimation of M_A defined by reference to this scale */
#endif
    
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        v_streaming=Get_CosmicRayStreamingVelocity(i,k_CRegy);
        DiffusionCoeff=0; CR_kappa_streaming=0; CRPressureGradScaleLength=Get_CosmicRayGradientLength(i,k_CRegy); /* set these for the bin as we get started */
#ifndef CRFLUID_ALT_DISABLE_STREAMING /* self-consistently calculate the diffusion coefficients for cosmic ray fluids; first the streaming part of this (kappa~v_stream*L_CR_grad) following e.g. Wentzel 1968, Skilling 1971, 1975, Holman 1979, as updated in Kulsrud 2005, Yan & Lazarian 2008, Ensslin 2011 */
        CR_kappa_streaming = GAMMA_COSMICRAY(k_CRegy) * v_streaming * CRPressureGradScaleLength; /* the diffusivity is now just the product of these two coefficients (all physical units) */
#endif
#if (CRFLUID_DIFFUSION_MODEL == 0) /* set diffusivity to a universal power-law scaling (constant per-bin)  */
        DiffusionCoeff = diffusion_coefficient_constant(i,k_CRegy); //  this is the input value of the diffusivity, for constant-kappa models
#endif
#if (CRFLUID_DIFFUSION_MODEL < 0) /* disable CR diffusion, specifically */
        DiffusionCoeff = 0; // no diffusion (but -can- allow streaming)
#endif
#if (CRFLUID_DIFFUSION_MODEL == 3) /* Farber et al. 2018 -- higher coeff in neutral gas, lower in ionized gas */
        DiffusionCoeff = (3.e29/unit_kappa_code) * (1.-f_ion + f_ion/30.); // 30x lower in neutral (note use f_ion directly here, not temperature as they do)
#endif
#if (CRFLUID_DIFFUSION_MODEL == 4) /* Wiener et al. 2017 style pure-streaming but with larger streaming speeds and limited losses, using their scaling for assumption that turbulent+non-linear Landau only dominate damping */
        double ni_m3=f_ion*(rho_cgs/PROTONMASS_CGS)/1.e-3, T6=temperature/1.e6, Lturbkpc=L_scale*UNIT_LENGTH_IN_KPC, Lgradkpc=CRPressureGradScaleLength*UNIT_LENGTH_IN_KPC, h0_fac=0.1*Get_Particle_Size(i)*All.cf_atime*All.cf_a2inv*UNIT_VEL_IN_KMS, dv2_10=0; for(k=0;k<3;k++) {int j; for(j=0;j<3;j++) {dv2_10 += SphP[i].Gradients.Velocity[j][k]*SphP[i].Gradients.Velocity[j][k]*h0_fac*h0_fac;}}
        double ecr_14 = SphP[i].CosmicRayEnergyPred[k_CRegy] * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_CGS / 1.0e-14; // CR energy density in CGS units //
        CR_kappa_streaming = GAMMA_COSMICRAY(k_CRegy) * CRPressureGradScaleLength * (v_streaming + (1./UNIT_VEL_IN_KMS)*(4.1*pow(MIN_REAL_NUMBER+ni_m3*T6,0.25)/pow(MIN_REAL_NUMBER+ecr_14*Lgradkpc,0.5) + 1.2*pow(MIN_REAL_NUMBER+dv2_10*ni_m3,0.75)/(MIN_REAL_NUMBER+ecr_14*sqrt(Lturbkpc)))); // convert to effective diffusivity
#endif
#if (CRFLUID_DIFFUSION_MODEL == 5) /* streaming at fast MHD wavespeed [just to see what it does] */
        CR_kappa_streaming = GAMMA_COSMICRAY(k_CRegy) * sqrt(v_streaming*v_streaming + cs_thermal*cs_thermal) * CRPressureGradScaleLength;
#endif
#if (CRFLUID_DIFFUSION_MODEL == 1) || (CRFLUID_DIFFUSION_MODEL == 2) || (CRFLUID_DIFFUSION_MODEL == 7) /* textbook extrinsic turbulence model: kappa~v_CR*r_gyro * B_bulk^2/(B_random[scale~r_gyro]^2) v_CR~c, r_gyro~p*c/(Z*e*B)~1e12 cm * RGV *(3 muG/B)  (RGV~1 is the magnetic rigidity). assuming a Kolmogorov spectrum */
        int scatter_modes = 0; /* default to using both Alfven+damped-fast modes */
#if (CRFLUID_DIFFUSION_MODEL==1)
        scatter_modes = -1; /* Alfven modes only*/
#endif
#if (CRFLUID_DIFFUSION_MODEL==2)
        scatter_modes = 1; /* Fast modes only*/
#endif
#if defined(CRFLUID_SET_ET_MODEL)
        scatter_modes = CRFLUID_SET_ET_MODEL; /* set to user-defined value */
#endif
        DiffusionCoeff = diffusion_coefficient_extrinsic_turbulence(scatter_modes,i,k_CRegy,M_A,L_scale,b_muG,vA_noion,rho_cgs,temperature,cs_thermal,nh0,nHe0,f_ion) / unit_kappa_code;
#endif
#if (CRFLUID_DIFFUSION_MODEL == 6) || (CRFLUID_DIFFUSION_MODEL == 7) /* self-confinement-based diffusivity */
        int target_bin_centering_for_CR_quantities = i; // if this = i, evaluate quantities like R_GV at the CR-energy weighted mean of the bin, if =-1, evaluate them at the bin center instead: important for some subtle effects especially if using numerical derivatives for correction terms
#ifdef CRFLUID_DIFFUSION_CORRECTION_TERMS
        target_bin_centering_for_CR_quantities = -1; // the correction terms depend on these being evaluated at their bin-centered locations
#endif
        double Omega_gyro_ifveqc=(0.00898734*b_muG/return_CRbin_CR_rigidity_in_GV(target_bin_centering_for_CR_quantities,k_CRegy)) * UNIT_TIME_IN_CGS, r_L=C_LIGHT_CODE/Omega_gyro_ifveqc, kappa_0=r_L*C_LIGHT_CODE; // some handy numbers for limiting extreme-kappa below. all in -physical- code units //
        CR_kappa_streaming = diffusion_coefficient_self_confinement(CRFLUID_SET_SC_MODEL,i,k_CRegy,M_A,L_scale,b_muG,vA_noion,rho_cgs,temperature,cs_thermal,nh0,nHe0,f_ion) / unit_kappa_code;
        if(!isfinite(CR_kappa_streaming)) {CR_kappa_streaming = 1.e30/unit_kappa_code;} /* apply some limiters since its very easy for the routine above to give wildly-large-or-small diffusivity, which wont make a difference compared to just 'small' or 'large', but will mess things up numerically */
        CR_kappa_streaming = DMIN( DMAX( DMIN(DMAX(CR_kappa_streaming,kappa_0) , 1.0e10*GAMMA_COSMICRAY(k_CRegy) * CRPressureGradScaleLength*CRFLUID_REDUCED_C_CODE(k_CRegy)) , 1.e25/unit_kappa_code ) , 1.e34/unit_kappa_code );
#endif

        /* -- ok, we've done what we came to do -- everything below here is pure-numerical, not physics, and should generally not be modified -- */

#if (CRFLUID_DIFFUSION_MODEL == 7) /* 'combined' extrinsic turbulence + self-confinement model: add scattering rates linearly (plus lots of checks to prevent unphysical bounds) */
        CR_kappa_streaming = 1. / (EPSILON_SMALL +  1./(CR_kappa_streaming+EPSILON_SMALL) + 1./(DiffusionCoeff+EPSILON_SMALL) ); DiffusionCoeff=0; if(!isfinite(CR_kappa_streaming)) {CR_kappa_streaming = 1.e30/unit_kappa_code;} // if scattering rates add linearly, this is a rough approximation to the total transport (essentially, smaller of the two dominates)
        CR_kappa_streaming = DMIN( DMAX( CR_kappa_streaming , kappa_0 ) , 1.0e10*GAMMA_COSMICRAY(k_CRegy) * CRPressureGradScaleLength*CRFLUID_REDUCED_C_CODE(k_CRegy) ); CR_kappa_streaming = DMIN( DMAX( CR_kappa_streaming , 1.e25/unit_kappa_code ) , 1.e34/unit_kappa_code );
#endif
        DiffusionCoeff = DiffusionCoeff + CR_kappa_streaming; //  add 'diffusion' and 'streaming' terms since enter numerically the same way
            
#ifndef CRFLUID_M1 /* now we apply a limiter to prevent the coefficient from becoming too large: cosmic rays cannot stream/diffuse with v_diff > c */
        // [all of this only applies if we are using the pure-diffusion description: the M1-type description should -not- use a limiter here, or negative kappa]
        double diffusion_velocity_limit=C_LIGHT_CODE, L_eff=DMAX(Get_Particle_Size(i)*All.cf_atime,CRPressureGradScaleLength); /* maximum diffusion velocity (set <c if desired) */
        double kappa_diff_vel = DiffusionCoeff * (GAMMA_COSMICRAY(k_CRegy)-1.) / L_eff; DiffusionCoeff *= 1 / (1 + kappa_diff_vel/diffusion_velocity_limit); /* caps maximum here */
#ifdef GALSF /* for multi-physics problems, we suppress diffusion [in the FLD-like limit] where it is irrelevant for timestepping-sake */
        DiffusionCoeff *= 1 / (1 + kappa_diff_vel/(0.01*diffusion_velocity_limit)); /* caps maximum here */
        double P_cr_Ratio = Get_Gas_CosmicRayPressure(i,k_CRegy) / (MIN_REAL_NUMBER + SphP[i].Pressure), P_min=1.0e-4; if(P_cr_Ratio < P_min) {DiffusionCoeff *= pow(P_cr_Ratio/P_min,2);}
        P_min = 1.0e-6; if(P_cr_Ratio < P_min) {DiffusionCoeff *= pow(P_cr_Ratio/P_min,2);}
#endif
        DiffusionCoeff /= (GAMMA_COSMICRAY(k_CRegy)-1.); // ensure correct units for subsequent operations //
#endif
        if((DiffusionCoeff<=0)||(isnan(DiffusionCoeff))) {DiffusionCoeff=0;} /* nan check */
        SphP[i].CosmicRayDiffusionCoeff[k_CRegy] = DiffusionCoeff; /* final assignment! */
    } // end CR bin loop
}



/* utility routine which handles the numerically-necessary parts of the CR 'injection' for you */
void inject_cosmic_rays(double CR_energy_to_inject, double injection_velocity, int source_type, int target, double *dir)
{
    if(CR_energy_to_inject <= 0) {return;}
    double f_injected[N_CR_PARTICLE_BINS]; f_injected[0]=1; int k_CRegy;;
#if (N_CR_PARTICLE_BINS > 1) /* add a couple steps to make sure injected energy is always normalized properly! */
    double sum_in=0.0; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {f_injected[k_CRegy]=CR_energy_spectrum_injection_fraction(k_CRegy,source_type,injection_velocity,0,target); sum_in+=f_injected[k_CRegy];}
    if(sum_in>0.0) {for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {f_injected[k_CRegy]/=sum_in;}} else {for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) {f_injected[k_CRegy]=1./N_CR_PARTICLE_BINS;}}
#endif
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double dEcr = evaluate_cr_transport_reductionfactor(target, k_CRegy, 0) * CR_energy_to_inject * f_injected[k_CRegy]; // normalized properly to sum to unity, and account for RSOL in injection rate [akin to RHD treatment]
        if(dEcr <= 0) {continue;}
        #pragma omp atomic
        SphP[target].CosmicRayEnergy[k_CRegy] += dEcr; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
        #pragma omp atomic
        SphP[target].CosmicRayEnergyPred[k_CRegy] += dEcr; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
#ifdef CRFLUID_M1
        double dir_mag=0, flux_mag=dEcr * CRFLUID_REDUCED_C_CODE(k_CRegy), dir_to_use[3]={0}; int k;
#ifdef MAGNETIC
        double B_dot_dir=0, Bdir[3]={0}; for(k=0;k<3;k++) {Bdir[k]=SphP[target].BPred[k]; B_dot_dir+=dir[k]*Bdir[k];} // the 'default' direction is projected onto B
        for(k=0;k<3;k++) {dir_to_use[k]=B_dot_dir*Bdir[k];} // launch -along- B, projected [with sign determined] by the intially-desired direction
#else
        for(k=0;k<3;k++) {dir_to_use[k]=dir[k];} // launch in the 'default' direction
#endif
        for(k=0;k<3;k++) {dir_mag += dir_to_use[k]*dir_to_use[k];}
        if(dir_mag <= 0) {dir_to_use[0]=0; dir_to_use[1]=0; dir_to_use[2]=1; dir_mag=1;}
        for(k=0;k<3;k++) {
            double dflux=flux_mag*dir_to_use[k]/sqrt(dir_mag);
            #pragma omp atomic
            SphP[target].CosmicRayFlux[k_CRegy][k]+=dflux; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
            #pragma omp atomic
            SphP[target].CosmicRayFluxPred[k_CRegy][k]+=dflux; // update injected CR energy. needs to be done thread-safely, but since the above routines dont depend on this, it should be safe to do here.
        }
#endif
    }
    return;
}



/* return CR pressure within a given bin */
double INLINE_FUNC Get_Gas_CosmicRayPressure(int i, int k_CRegy)
{
    if((P[i].Mass > 0) && (SphP[i].Density>0) && (SphP[i].CosmicRayEnergyPred[k_CRegy] > 0))
    {
        return (GAMMA_COSMICRAY(k_CRegy)-1.) * (SphP[i].CosmicRayEnergyPred[k_CRegy] * SphP[i].Density) / P[i].Mass; // cosmic ray pressure = (4/3-1) * e_cr = 1/3 * (E_cr/Vol) //
    } else {return 0;}
}



/* return CR gradient scale length, with various physical limiters applied: not intended for pure numerical gradient-length calculations (where units dont matter), but for preventing some unphysical situations */
double Get_CosmicRayGradientLength(int i, int k_CRegy)
{
    /* now we need the -parallel- cosmic ray pressure or energy density scale length */
    double CRPressureGradMag = 0.0;
    int k; for(k=0;k<3;k++) {CRPressureGradMag += SphP[i].Gradients.CosmicRayPressure[k_CRegy][k]*SphP[i].Gradients.CosmicRayPressure[k_CRegy][k];}
    CRPressureGradMag = sqrt(1.e-46 + CRPressureGradMag); // sqrt to make absolute value
#ifdef MAGNETIC /* with anisotropic transport, we really want the -parallel- gradient scale-length, so need another factor here */
    double B2_tot=0.0, b=0; CRPressureGradMag=0; for(k=0;k<3;k++) {b=Get_Gas_BField(i,k); B2_tot+=b*b; CRPressureGradMag+=b*SphP[i].Gradients.CosmicRayPressure[k_CRegy][k];} // note, this is signed!
    CRPressureGradMag = sqrt((1.e-40 + CRPressureGradMag*CRPressureGradMag) / (1.e-46 + B2_tot)); // divide B-magnitude to get scalar magnitude, and take sqrt[(G.P)^2] to get absolute value
#endif
    
    /* limit the scale length: if too sharp, need a slope limiter at around the particle size */
    double L_gradient_min = Get_Particle_Size(i) * All.cf_atime;
    /* limit this scale length; if the gradient is too shallow, there is no information beyond a few smoothing lengths, so we can't let streaming go that far */
    double L_gradient_max = DMAX(1000.*L_gradient_min, 500.0*PPP[i].Hsml*All.cf_atime);

    /* also, physically, cosmic rays cannot stream/diffuse with a faster coefficient than ~v_max*L_mean_free_path, where L_mean_free_path ~ 2.e20 * (cm^-3/n) [collisional here] */
    double nH_cgs = SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_NHCGS;
    double L_mean_free_path = (3.e25 / nH_cgs) / UNIT_LENGTH_IN_CGS;
    L_gradient_max = DMIN(L_gradient_max, L_mean_free_path);
    
    double CRPressureGradScaleLength = Get_Gas_CosmicRayPressure(i,k_CRegy) / CRPressureGradMag * All.cf_atime;
    if(CRPressureGradScaleLength > 0) {CRPressureGradScaleLength = 1.0/(1.0/CRPressureGradScaleLength + 1.0/L_gradient_max);} else {CRPressureGradScaleLength=0;}
    CRPressureGradScaleLength = sqrt(L_gradient_min*L_gradient_min + CRPressureGradScaleLength*CRPressureGradScaleLength);
    return CRPressureGradScaleLength; /* this is returned in -physical- units */
}



/* return the effective CR 'streaming' velocity for sub-grid [unresolved] models with streaming velocity set by e.g. the Alfven speed along the gradient of the CR pressure */
double Get_CosmicRayStreamingVelocity(int i, int k_CRegy)
{
#if (defined(CRFLUID_M1) && !defined(CRFLUID_ALT_FLUX_FORM_JOCH)) || (defined(CRFLUID_ALT_DISABLE_STREAMING))
    return 0; // with this option, the streaming is included by default in the flux equation, different from what we do below where we include it as an 'effective diffusivity'
#endif
    double v_streaming = Get_Gas_ion_Alfven_speed_i(i); // limit to Alfven speed, but put some limiters for extreme cases //
#if defined(CRFLUID_M1)
    v_streaming = DMIN(v_streaming , CRFLUID_REDUCED_C_CODE(k_CRegy)); // limit to maximum transport speed //
#endif
    v_streaming *= All.cf_afac3; // converts to physical units and rescales according to chosen coefficient //
    return v_streaming;
}


/* routine to quickly estimate the atomic mass of the CR particles (in proton masses): making assumptions here [e.g. no positrons] -- could always hard-code some choice here for more complicated scenarios */
double return_CRbin_CRmass_in_mp(int target, int k_CRegy)
{
    int species = return_CRbin_CR_species_ID(k_CRegy);
    if(species < 0) {return 0.000544618;} // electrons or positrons
        else {
            if(species == 1) {return  1.0;} // p
            if(species == 7) {return  1.0;} // pbar (anti-p)
            if(species == 2) {return 10.8;} // B: mostly 11, but non-negligible 10, with 10B/11B ~ 0.58, so should really account for both. won't make huge difference here but some care needed in propagation models and abundance models and corrections
            if(species == 3) {return 12.0;} // C
            if(species == 4) {return  9.0;} // Be7-9: for now going with 9 as the reference, but could do 7 as well...
            if(species == 5) {return  10.;} // Be10
            if(species == 6) {return 14.8;} // CNO bin
            /* if don't have a rule for the species, default to 2x Z for heavy species */
            double Z_abs = fabs(return_CRbin_CR_charge_in_e(target,k_CRegy));
            if(Z_abs > 1.5) {return 2.*Z_abs;} else {return 1;} // nuclei, making a simple assumption about the nuclear weight [can modify easily with custom flag!]
    }
}


/* routine which returns the dimensionless velocity beta = v/c of the CRs (this is not the RSOL, but true c, for e.g. cooling, energies, etc) */
double return_CRbin_beta_factor(int target, int k_CRegy)
{
    double m_cr_mp=return_CRbin_CRmass_in_mp(target,k_CRegy); // mass in proton masses
    double q = return_CRbin_CR_rigidity_in_GV(target,k_CRegy)*1.06579*fabs(return_CRbin_CR_charge_in_e(target,k_CRegy))/m_cr_mp; // dimensionless factor to convert from R to beta
    double gamma = sqrt(1.+q*q), beta = q/gamma;
    return beta;
}


/* routine which returns the dimensionless lorentz factor gamma=1/sqrt[1-beta^2] of the CRs (this is not the RSOL, but true c, for e.g. cooling, energies, etc) */
double return_CRbin_gamma_factor(int target, int k_CRegy)
{
    double m_cr_mp=return_CRbin_CRmass_in_mp(target,k_CRegy); // mass in proton masses
    double q = return_CRbin_CR_rigidity_in_GV(target,k_CRegy)*1.06579*fabs(return_CRbin_CR_charge_in_e(target,k_CRegy))/m_cr_mp; // dimensionless factor to convert from R to beta
    return sqrt(1.+q*q);
}


/* routine which returns the effective adiabatic index of CRs, P_cr = (gamma_eos - 1) * e_kinetic */
double gamma_eos_of_crs_in_bin(int k_CRegy)
{
    return (4. + 1./return_CRbin_gamma_factor(-1,k_CRegy)) / 3.;
}

/* routine which returns the CR kinetic energy in GeV, for a given rigidity, etc. */
double return_CRbin_kinetic_energy_in_GeV(int target, int k_CRegy)
{
    double m_cr_mp=return_CRbin_CRmass_in_mp(target,k_CRegy); // mass in proton masses
    double R_GV = return_CRbin_CR_rigidity_in_GV(target,k_CRegy), Z = fabs(return_CRbin_CR_charge_in_e(target,k_CRegy));
    double q = R_GV*Z*1.06579/m_cr_mp; // dimensionless factor to convert from R to beta
    double gamma = sqrt(1.+q*q), beta = q/gamma;
    double KE_fac = DMAX(0.,(1.-sqrt(DMAX(0.,1.-beta*beta)))) / DMAX(beta,MIN_REAL_NUMBER); if(beta < 0.01) {KE_fac = 0.5*beta;} // non-relativistic expansion used when gamma very close to one, to prevent numerical errors
    return R_GV * Z * KE_fac;
}


/* routine which returns the CR number density at a given bin in cm^-3 */
double return_CRbin_numberdensity_in_cgs(int target, int k_CRegy)
{
    double e_CR_tot = SphP[target].CosmicRayEnergyPred[k_CRegy] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * (UNIT_PRESSURE_IN_CGS);
    double ECR_per  = return_CRbin_kinetic_energy_in_GeV(target,k_CRegy) * 0.00160218; /* converts to energy in erg */
    return e_CR_tot / ECR_per;
}


/* handy function that just returns the radiation energy density in eV/cm^-3, physical units. purely here to save us time re-writing this */
double get_cell_Urad_in_eVcm3(int i)
{
    double erad = 0.26*All.cf_a3inv/All.cf_atime; // default with the CMB energy density, which we assume here is a baseline minimum
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
    int kfreq; double e_units = (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_EV;
    for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad+=SphP[i].Rad_E_gamma_Pred[kfreq]*e_units;}
#else
    double uRad_MW = 0.31 + 0.66; /* dust (0.31) and stars (0.66) for Milky way ISRF from Draine (2011); want this to decay as we approach the IGM (where CMB totally dominates) */
    double prefac_rad=1, rho_cgs=SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS; if(All.ComovingIntegrationOn) {double rhofac = rho_cgs / (1000.*COSMIC_BARYON_DENSITY_CGS);
        if(rhofac < 0.2) {prefac_rad=0;} else {if(rhofac > 200.) {prefac_rad=1;} else {prefac_rad=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off stellar background for any gas with density unless it's >1000 times the cosmic mean density
    prefac_rad *= rho_cgs/(0.01*PROTONMASS_CGS + rho_cgs); // cut off below low densities, ~1e-2
    erad += uRad_MW * prefac_rad;
#endif
    return erad;
}



/* return pre-factor for CR streaming losses, such that loss rate dE/dt = -E * streamfac */
double CR_get_streaming_loss_rate_coefficient(int target, int k_CRegy)
{
    double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target), streamfac = 0;
    if((SphP[target].CosmicRayEnergy[k_CRegy] <= MIN_REAL_NUMBER) || (SphP[target].CosmicRayEnergy[k_CRegy] <= MIN_REAL_NUMBER) || (dt <= 0)) {return 0;}
#if !defined(CRFLUID_ALT_DISABLE_STREAMING)
    double vstream_0 = Get_CosmicRayStreamingVelocity(target,k_CRegy), vA=Get_Gas_ion_Alfven_speed_i(target); /* define naive streaming and Alfven speeds */

#if defined(CRFLUID_M1) && !defined(CRFLUID_ALT_FLUX_FORM_JOCH)
    double v_flux_eff=0; int k; for(k=0;k<3;k++) {v_flux_eff += SphP[target].CosmicRayFluxPred[k_CRegy][k] * SphP[target].CosmicRayFluxPred[k_CRegy][k];} // need magnitude of flux vector
    if(v_flux_eff > 0) {v_flux_eff=sqrt(v_flux_eff) / (MIN_REAL_NUMBER + SphP[target].CosmicRayEnergyPred[k_CRegy]);} else {v_flux_eff=0;} // effective speed of CRs = |F|/E
    int target_for_cr_gamma = target; // if this = -1, use the gamma factor at the bin-center for evaluating this, if this = target, use the mean gamma of the bin, weighted by the CR energy -- won't give exactly the same result here
    double gamma_0=return_CRbin_gamma_factor(target_for_cr_gamma,k_CRegy), gamma_fac=gamma_0/(gamma_0-1.), beta_fac=return_CRbin_beta_factor(target_for_cr_gamma,k_CRegy); // lorentz factor here, needed in next line, because the loss term here scales with -total- energy, not kinetic energy
    if(beta_fac<0.1) {gamma_fac=2./(beta_fac*beta_fac) -0.5 - 0.125*beta_fac*beta_fac;} // avoid accidental nan
    streamfac = (vA * (beta_fac*beta_fac) / fabs(3.*SphP[target].CosmicRayDiffusionCoeff[k_CRegy])) * ((gamma_fac) * return_CRbin_nuplusminus_asymmetry(target,k_CRegy) * v_flux_eff/CosmicRayFluid_RSOL_Corrfac(k_CRegy) - (3.*(GAMMA_COSMICRAY(k_CRegy)-1.) + (gamma_fac)) * vA * (2./3.) * return_cosmic_ray_anisotropic_closure_function_threechi(target,k_CRegy)); // this is (vA/[3kappa])*(F - 2*chifac*vA*(ecr+3*Pcr))/ecr, using the 'full F' [corrected back from rsol, b/c rsol correction moves outside this for loss terms]
    double sfac_max = fabs(0.1 * SphP[target].CosmicRayEnergy[k_CRegy] / (MIN_REAL_NUMBER + dt));
    if(fabs(streamfac) > sfac_max) {streamfac *= sfac_max / fabs(streamfac);}
    return streamfac; // probably want to limit to make sure above doesn't take on too extreme a value... also above, initially only had positive term since this removes energy from CRs when streaming super-Alfvenically, but when streaming sub-Alfvenically, could this become a source term with energy going into CRs? seems problematic if vA very high, but then scattering would work inefficiently... so plausible, but really need to be careful again about magnitude...
#endif
    
    if(vA>0) {vstream_0 = DMIN(vA, vstream_0);} /* account for the fact that the loss term is always [or below] the Alfven speed, regardless of the bulk streaming speed */
    double L_cr = Get_CosmicRayGradientLength(target,k_CRegy), v_st_eff = SphP[target].CosmicRayDiffusionCoeff[k_CRegy] / (GAMMA_COSMICRAY(k_CRegy) * L_cr + MIN_REAL_NUMBER); // maximum possible streaming speed from combined diffusivity
    streamfac = fabs((1./3.) * DMIN(v_st_eff, vstream_0) / L_cr); // if upper-limit to streaming is less than nominal 'default' v_stream/loss term, this should be lower too
#endif
    return streamfac;
}



/* routine to do the drift/kick operations for CRs: mode=0 is kick, mode=1 is drift */
#if !defined(CRFLUID_EVOLVE_SCATTERINGWAVES)
double CosmicRay_Update_DriftKick(int i, double dt_entr, int mode)
{
    if(dt_entr <= 0) {return 0;} // no update

    int k_CRegy;
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        int k; double eCR, u0; k=0; if(mode==0) {eCR=SphP[i].CosmicRayEnergy[k_CRegy]; u0=SphP[i].InternalEnergy;} else {eCR=SphP[i].CosmicRayEnergyPred[k_CRegy]; u0=SphP[i].InternalEnergyPred;} // initial energy
        if(u0<All.MinEgySpec) {u0=All.MinEgySpec;} // enforced throughout code
        if(eCR < 0) {eCR=0;} // limit to physical values
        double closure_f1, closure_f2; closure_f1=1, closure_f2=0; // prefactors for below
#if defined(CRFLUID_M1) && !defined(CRFLUID_ALT_FLUX_FORM_JOCH)
        double three_chi = return_cosmic_ray_anisotropic_closure_function_threechi(i,k_CRegy); // 3*chi = 3*(1-<mu^2>)/2 closure function //
        closure_f1 = 3.-2.*three_chi; closure_f2 = 1.-three_chi; // prefactors for both terms below //
#endif
#if defined(CRFLUID_M1) /* CR FLUX VECTOR UPDATE */
        // this is the exact solution for the CR flux-update equation over a finite timestep dt: it needs to be solved this way [implicitly] as opposed to explicitly for dt because in the limit of dt_cr_dimless being large, the problem exactly approaches the diffusive solution
        double DtCosmicRayFlux[3]={0}, flux[3]={0}, CR_veff[3]={0}, CR_vmag=0, q_cr=0, cr_speed=CRFLUID_REDUCED_C_CODE(k_CRegy), rsol_correction_factor=CosmicRayFluid_RSOL_Corrfac(k_CRegy), V_i=P[i].Mass/SphP[i].Density, P0_cr, fac_for_DtCosmicRayFlux; P0_cr=Get_Gas_CosmicRayPressure(i,k_CRegy);
        //cr_speed = DMAX(All.cf_afac3*SphP[i].MaxSignalVel , CRFLUID_REDUCED_C_CODE(k_CRegy)); // may give slightly improved performance with modified rsol formulation
        cr_speed = DMAX(All.cf_afac3*SphP[i].MaxSignalVel , DMIN(CRFLUID_REDUCED_C_CODE(k_CRegy) , 10.*fabs(SphP[i].CosmicRayDiffusionCoeff[k_CRegy])/(Get_Particle_Size(i)*All.cf_atime)));
        fac_for_DtCosmicRayFlux = -rsol_correction_factor * fabs(SphP[i].CosmicRayDiffusionCoeff[k_CRegy]) * V_i / (GAMMA_COSMICRAY(k_CRegy)-1.);
        for(k=0;k<3;k++) {DtCosmicRayFlux[k] = SphP[i].Gradients.CosmicRayPressure[k_CRegy][k];}
#ifdef MAGNETIC // do projection onto field lines
        double bhat[3]={0}, B0[3]={0}, Bmag2=0, Bmag=0, bbGB=0, DtCRDotBhat=0;
        for(k=0;k<3;k++)
        {
            if(mode==0) {B0[k]=SphP[i].B[k]/V_i;} else {B0[k]=SphP[i].BPred[k]/V_i;}
            DtCRDotBhat += DtCosmicRayFlux[k] * B0[k]; Bmag2 += B0[k]*B0[k]; bhat[k]=B0[k];
        }
        if(Bmag2 > 0) {Bmag=sqrt(Bmag2); for(k=0;k<3;k++) {bhat[k]/=Bmag;}}
        for(k=0;k<3;k++) {int k2; for(k2=0;k2<3;k2++) {bbGB += -bhat[k]*bhat[k2]*SphP[i].Gradients.B[k][k2]/Bmag;}}
        for(k=0;k<3;k++) {
            DtCosmicRayFlux[k] = fac_for_DtCosmicRayFlux * (closure_f1*DtCRDotBhat*B0[k]/Bmag2 + closure_f2*P0_cr*bbGB);
        }
#endif
#if defined(CRFLUID_M1) && !defined(CRFLUID_ALT_FLUX_FORM_JOCH) && !defined(CRFLUID_ALT_DISABLE_STREAMING)
        double v_Alfven = three_chi * Get_Gas_ion_Alfven_speed_i(i) * return_CRbin_nuplusminus_asymmetry(i,k_CRegy); /* define naive streaming and Alfven speeds */
        double dt_f_m=0; for(k=0;k<3;k++) {dt_f_m+=DtCosmicRayFlux[k]*DtCosmicRayFlux[k];}
        if(dt_f_m>0) {for(k=0;k<3;k++) {DtCosmicRayFlux[k] += rsol_correction_factor * (DtCosmicRayFlux[k]/sqrt(dt_f_m)) * v_Alfven * (GAMMA_COSMICRAY(k_CRegy) * eCR);}} // (tilde[c]/c) * v_a * (ecr+Pcr), in same direction as gradient wants to 'push' naturally [natural direction of F]
#endif
        if(mode==0) {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFlux[k_CRegy][k];}} else {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFluxPred[k_CRegy][k];}}
#ifdef MAGNETIC // do projection onto field lines
        double fluxmag=0, fluxdot=0; for(k=0;k<3;k++) {fluxmag+=flux[k]*flux[k]; fluxdot+=flux[k]*B0[k];}
        if(fluxmag>0) {fluxmag=sqrt(fluxmag);} else {fluxmag=0;}
        if(fluxdot<0) {fluxmag*=-1;} // points down-field
        // before acting on the 'stiff' sub-system, account for the 'extra' advection term that accounts for 'twisting' of B: note more careful derivation shows this is sub-leading order in v/c, should not be included here
        //double fac_bv=0; for(k=0;k<3;k++) {fac_bv += All.cf_a2inv * bhat[k] * (bhat[0]*SphP[i].Gradients.Velocity[k][0] + bhat[1]*SphP[i].Gradients.Velocity[k][1] + bhat[2]*SphP[i].Gradients.Velocity[k][2]);}
        //if(All.ComovingIntegrationOn) {fac_bv += All.cf_hubble_a;} // adds cosmological/hubble flow term here [not included in peculiar velocity gradient]
        //fluxmag *= exp(-DMAX(-2.,DMIN(2.,rsol_correction_factor*fac_bv*dt_entr))); // limit factor for change here, should be small given Courant factor, then update flux term accordingly, before next step -- acts like a mod of the divv term //
        if(Bmag2>0) {for(k=0;k<3;k++) {flux[k] = fluxmag * B0[k] / sqrt(Bmag2);}} // re-assign to be along field
#endif
        int target_for_CR_beta_factor = i; // if this =1, use energy-weighted mean value in bin for CR beta, otherwise if =-1, use median point of bin
#ifdef CRFLUID_DIFFUSION_CORRECTION_TERMS
        target_for_CR_beta_factor = -1;
#endif
        double beta_fac = return_CRbin_beta_factor(target_for_CR_beta_factor,k_CRegy); // velocity beta, to account for non-relativistic CRs
        double dt_cr_dimless = dt_entr * beta_fac*beta_fac * cr_speed*cr_speed * (1./3.) / (MIN_REAL_NUMBER + fabs(SphP[i].CosmicRayDiffusionCoeff[k_CRegy] * rsol_correction_factor));
        dt_cr_dimless = DMIN(dt_cr_dimless , 0.1); // arbitrary limiter here for some additional numerical stability
        if((dt_cr_dimless > 0)&&(dt_cr_dimless < 20.)) {q_cr = exp(-dt_cr_dimless);} // factor for CR interpolation
        for(k=0;k<3;k++) {flux[k] = q_cr*flux[k] + (1.-q_cr)*DtCosmicRayFlux[k];} // updated flux
        for(k=0;k<3;k++) {CR_veff[k]=flux[k]/(eCR+MIN_REAL_NUMBER); CR_vmag+=CR_veff[k]*CR_veff[k];} // effective streaming speed
        if((CR_vmag <= 0) || (isnan(CR_vmag))) // check for valid numbers
        {
            for(k=0;k<3;k++) {flux[k]=0; for(k=0;k<3;k++) {CR_veff[k]=0;}} // zero if invalid
        } else {
            //double CR_vmax = CRFLUID_REDUCED_C_CODE(k_CRegy); // enforce a hard upper limit here, though shouldn't be needed with modern formulation
            double CR_vmax = CRFLUID_M1; // [use stricter limit here, for timestep concordance] enforce a hard upper limit here, though shouldn't be needed with modern formulation
            CR_vmag = sqrt(CR_vmag); if(CR_vmag > CR_vmax) {for(k=0;k<3;k++) {flux[k]*=CR_vmax/CR_vmag; CR_veff[k]*=CR_vmax/CR_vmag;}} // limit flux to free-streaming speed [as with RT]
        }
        if(mode==0) {for(k=0;k<3;k++) {SphP[i].CosmicRayFlux[k_CRegy][k]=flux[k];}} else {for(k=0;k<3;k++) {SphP[i].CosmicRayFluxPred[k_CRegy][k]=flux[k];}}
#endif
    
        /* update scalar CR energy. first update the CR energies from fluxes. since this is positive-definite, some additional care is needed */
        double dCR_dt = SphP[i].DtCosmicRayEnergy[k_CRegy], eCR_tmp = eCR;
        double dCR = dCR_dt*dt_entr, dCRmax = 1.e10*(eCR_tmp+MIN_REAL_NUMBER);
#if defined(GALSF)
        dCRmax = DMAX(2.0*eCR_tmp , 0.1*u0*P[i].Mass);
#endif
        if(dCR > dCRmax) {dCR=dCRmax;} // don't allow excessively large values
        if(dCR < -eCR_tmp) {dCR=-eCR_tmp;} // don't allow it to go negative
        double eCR_0, eCR_00; eCR_00 = eCR_tmp; eCR_tmp += dCR; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy]=eCR_tmp;} else {SphP[i].CosmicRayEnergyPred[k_CRegy]=eCR_tmp;} // updated energy
        eCR_0 = eCR_tmp; // save this value for below
        
        
#if defined(COOLING_OPERATOR_SPLIT)
        /* now need to account for the adiabatic heating/cooling of the 'fluid', here, with gamma=GAMMA_COSMICRAY(k_CRegy) */
        double dCR_div = CR_calculate_adiabatic_gasCR_exchange_term(i, dt_entr, (GAMMA_COSMICRAY(k_CRegy)-1.)*eCR_tmp, mode); // this will handle the update below - separate subroutine b/c we want to allow it to appear in a couple different places
        double uf = DMAX(u0 - dCR_div/P[i].Mass , All.MinEgySpec); // final updated value of internal energy per above
        if(mode==0) {SphP[i].InternalEnergy = uf;} else {SphP[i].InternalEnergyPred = uf;} // update gas
        if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy] += dCR_div;} else {SphP[i].CosmicRayEnergyPred[k_CRegy] += dCR_div;} // update CRs: note if explicitly evolving spectrum, this is done separately below //
#endif

    } // loop over CR bins complete
    return 1;
}
#endif







/* optional code to allow the RSOL to depend on bin energy, still testing this */
double return_CRbin_M1speed(int k_CRegy)
{
#if defined(CRFLUID_M1)
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
#if defined(CRFLUID_ALT_VARIABLE_RSOL) // experimental block here //
    double R = CR_global_rigidity_at_bin_center[k_CRegy];
    double f = All.CosmicRayDiffusionCoeff * UNIT_LENGTH_IN_KPC * pow(R , 0.8);
    if(f > CRFLUID_M1) {return f;}
#endif
#endif
    return CRFLUID_M1;
#else
    return C_LIGHT_CODE;
#endif
}


/* estimate amount by which flux of CRs has been reduced relative to solution with c_reduced = c_true, for RSOL with M1 */
double evaluate_cr_transport_reductionfactor(int target, int k_CRegy, int mode)
{
#if defined(CRFLUID_M1)
#if defined(CRFLUID_ALT_RSOL_FORM)
    return CosmicRayFluid_RSOL_Corrfac(k_CRegy); // uniform reduction factor for all terms
#else
    double kappa = SphP[target].CosmicRayDiffusionCoeff[k_CRegy]; /* diffusion coefficient [physical units] */
    double fluxmag=0, Bmag=0, gradmag=0, Lgrad=0, veff=0, P0=Get_Gas_CosmicRayPressure(target,k_CRegy); int m;
    for(m=0;m<3;m++) {
        double f_0=SphP[target].CosmicRayFluxPred[k_CRegy][m], g_0=SphP[target].Gradients.CosmicRayPressure[k_CRegy][m], B_0=f_0;
#ifdef MAGNETIC
        B_0 = Get_Gas_BField(target,m);
#endif
        Bmag += B_0*B_0; fluxmag += f_0*B_0; gradmag += g_0*B_0;
    }
    if(Bmag>0) {fluxmag=fabs(fluxmag)/sqrt(Bmag); gradmag=fabs(gradmag)/sqrt(Bmag);}
    if(gradmag>0) {Lgrad = All.cf_atime * P0 / gradmag;}
    if(fluxmag>0 && SphP[target].CosmicRayEnergyPred[k_CRegy] > MIN_REAL_NUMBER) {veff = fluxmag / SphP[target].CosmicRayEnergyPred[k_CRegy];}
    //if(mode==0) {veff = CRFLUID_REDUCED_C_CODE(k);} // we're injecting, so the relevant speed here is just the injection speed
    int use_injectionmod=0; if(mode==0) {use_injectionmod=1;} else {if(return_CRbin_CR_species_ID(k_CRegy) < 0) {use_injectionmod=1;}} // for injection, or leptons, where loss=injection at high-RGV (high-diffcoeff, so high veff), more accurate to match suppression factors this way
    if(use_injectionmod) {veff = CRFLUID_M1;} else {veff = DMIN(veff, CRFLUID_M1);} // we're injecting, so the relevant speed here is just the injection speed (note we speed-limit flux to this for timestepping reasons)
    if(use_injectionmod) {Lgrad = 5./UNIT_LENGTH_IN_KPC;} // set initial gradient length to a constant to reduce noise [looking at newer tests, don't -really- need this, but doesn't hurt either]
    double v_max = DMIN( C_LIGHT_CODE , kappa / (MIN_REAL_NUMBER + Lgrad) ); // attempt at a limiter function here to determine if being flux-limited in the equations below //
    double RSOL_over_v_desired = veff / (MIN_REAL_NUMBER + v_max);
    if(isfinite(RSOL_over_v_desired) && (RSOL_over_v_desired > 0) && (RSOL_over_v_desired < MAX_REAL_NUMBER)) {if(RSOL_over_v_desired < 1) {return RSOL_over_v_desired;}}
    return 1;
#endif
#else
    return 1;
#endif
}


/*<! closure function needed for arbitrarily anisotropic CR distribution function, from Hopkins '21 */
double return_cosmic_ray_anisotropic_closure_function_threechi(int target, int k_CRegy)
{
#if !defined(CRFLUID_ALT_M1_ISO_CLOSURE) && defined(CRFLUID_M1)
    double fluxmag2=0,ecr,ecrv,f,mu1_2,mu2; int k; ecr=SphP[target].CosmicRayEnergyPred[k_CRegy];
    for(k=0;k<3;k++) {f=SphP[target].CosmicRayFluxPred[k_CRegy][k]; fluxmag2+=f*f;}
#if defined(CRFLUID_ALT_RSOL_FORM)
    double v_eff_cr = CRFLUID_REDUCED_C_CODE(k_CRegy); // universal reduction factor
#else
    double kappa=SphP[target].CosmicRayDiffusionCoeff[k_CRegy], Lgrad_inv=0, P=0; for(k=0;k<3;k++) {P=SphP[target].Gradients.CosmicRayPressure[k_CRegy][k]; Lgrad_inv+=P*P;}
    Lgrad_inv = (sqrt(Lgrad_inv) / Get_Gas_CosmicRayPressure(target, k_CRegy)) / All.cf_atime;
    double v_eff_cr = DMIN(DMAX(CRFLUID_REDUCED_C_CODE(k_CRegy) , kappa*Lgrad_inv) , C_LIGHT_CODE); // more complicated factor b/c transport not universally slowed-down
#endif
    ecrv=ecr*v_eff_cr; f=fluxmag2/(ecrv*ecrv); mu1_2=DMAX(0,DMIN(1,f));
    if(!isfinite(mu1_2) || fluxmag2 < MIN_REAL_NUMBER || ecrv < MIN_REAL_NUMBER || !isfinite(f) || f < MIN_REAL_NUMBER || !isfinite(ecrv) || !isfinite(fluxmag2)) {mu1_2=0;} else {if(isfinite(f) && f > 1) {mu1_2=1;}} // safety checks for initialization where may have 0's, etc.
    mu1_2 = DMAX(0,DMIN(1,mu1_2)); // one more safety check
    mu2 = (3.+4.*mu1_2) / (5.+2.*sqrt(4.-3.*mu1_2)); // actual closure relation
    return 3.*(1.-mu2)/2.; // definition of chi variable
#endif
    return 1;
}


/* routine which returns the typical absolute value of the rigidity of a given CR [in GV] in a given 'bin' of our multi-bin approximation */
double return_CRbin_CR_rigidity_in_GV(int target, int k_CRegy)
{
    double R = 1;
#if (N_CR_PARTICLE_BINS > 1)    /* insert physics here */
#if (N_CR_PARTICLE_BINS == 2) /* one-bin protons, one electrons */
    double Rv[2]={1.8 , 0.6}; R=Rv[k_CRegy]; // approximate peak energies of each from Cummings et al. 2016 Fig 15
#endif
#if (N_CR_PARTICLE_BINS > 2) /* arbitrary numbers of bins for CR spectra, assumes binning same in e- and p */
    if(target >= 0) {R=CR_return_mean_rigidity_in_bin_in_GV(target,k_CRegy);} else {R=CR_global_rigidity_at_bin_center[k_CRegy];} // this is pre-defined globally for this bin list
#endif
#endif
    return R;
}

/* routine which returns the typical charge of CRs in a given 'bin' of our multi-bin approximation */
double return_CRbin_CR_charge_in_e(int target, int k_CRegy)
{
#if (N_CR_PARTICLE_BINS > 2)
    return CR_global_charge_in_bin[k_CRegy]; // this is pre-defined globally for this bin list
#endif
#if (N_CR_PARTICLE_BINS == 2) /* one-bin protons, one electrons */
    if(k_CRegy==0) {return 1} else {return -1;} // proton, then e- bin
#endif
    return 1; // default one-bin model is only protons
}

/* routine which returns the species ID of the CRs allowing for different code modes */
int return_CRbin_CR_species_ID(int k_CRegy)
{
#if (N_CR_PARTICLE_BINS == 2) /* one-bin protons, one electrons */
    if(k_CRegy==0) {return 1} else {return -1;}; // proton, then e- bin
#endif
    return 1; // default to protons
}



/* calculate the -ion- Alfven speed in a given element, relevant for very short-wavelength modes with frequency larger than the ion-neutral collision timescale (relevant for CRs in particular) */
double Get_Gas_ion_Alfven_speed_i(int i)
{
#if defined(MAGNETIC)
    return Get_Gas_thermal_soundspeed_i(i); // if no B-fields, just assume Alfven speed equal to thermal sound speed
#endif
    double vA = Get_Gas_Alfven_speed_i(i); // normal ideal-MHD Alfven speed
#ifdef CRFLUID_ION_ALFVEN_SPEED
    vA /= sqrt(1.e-10 + Get_Gas_Ionized_Fraction(i)); // Alfven speed of interest is that of the ions alone, not the ideal MHD Alfven speed //
#endif
    return vA;
}


/* calculate |nu_+ - nu_-|/|nu_+ + nu_-|, i.e. the asymmetry between scattering by forward-vs-backward propagating waves. in SC limit, this is 1, in ET limit with fully-isotropic scattering, this is 0. relevant here because we want in principle to be able to interpolate between the two, and that's important in particular in the strong ion-neutral damping limit. */
double return_CRbin_nuplusminus_asymmetry(int i, int k_CRegy)
{
    return 1; // default to unity; basically always true with reasonable SC growth rates, even for extremely low ionization fractions //
}


#endif // closes block for entire file for COSMIC_RAY_FLUID





/* return total CR energy density associated with a cell */
double INLINE_FUNC Get_CosmicRayEnergyDensity_cgs(int i)
{
    if(i<=0) {return 0;}
#ifdef COSMIC_RAY_FLUID
    double u_cr=0; int k; for(k=0;k<N_CR_PARTICLE_BINS;k++) {u_cr += SphP[i].CosmicRayEnergyPred[k];}
    return u_cr * (SphP[i].Density*All.cf_a3inv / P[i].Mass) * UNIT_PRESSURE_IN_CGS;
#endif
    return 1.6e-12; // eV/cm-3, approximate from Cummings et al. 2016 V1 data
}


/* return total CR ionization rate zeta_cr in s^-1 associated with a cell */
double Get_CosmicRayIonizationRate_cgs(int i)
{
#if defined(COSMIC_RAY_FLUID) && (N_CR_PARTICLE_BINS > 2)
    double ecr_units=(SphP[i].Density*All.cf_a3inv/P[i].Mass)*UNIT_PRESSURE_IN_CGS, zeta_cr=0; int k;
    for(k=0;k<N_CR_PARTICLE_BINS;k++)
    {
        double T_GeV=return_CRbin_kinetic_energy_in_GeV(-1,k), beta=return_CRbin_beta_factor(-1,k), Z=return_CRbin_CR_charge_in_e(-1,k), gamma=return_CRbin_gamma_factor(-1,k);
        zeta_cr += 3.43e-18 * (Z*Z/T_GeV) * ((1.-0.069*beta*beta+0.14*log(beta*gamma))/beta) * (SphP[i].CosmicRayEnergyPred[k]*ecr_units); // cross sections from standard Bethe-Blocke formulation, valid at all CR energies we consider explicitly
    }
    return zeta_cr;
#else
    return 1.e-5 * Get_CosmicRayEnergyDensity_cgs(i); // scales following Cummings et al. 2016 to 1.6e-17 per eV/cm^3
#endif
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
#if defined(COSMIC_RAY_FLUID) || defined(FLAG_NOT_IN_PUBLIC_CODE)
#if (N_CR_PARTICLE_BINS > 2)
    double e_heat=0, e_CR_units_0=(SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_CGS / nHcgs; int k_CRegy;
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        double e_cr_units = SphP[target].CosmicRayEnergyPred[k_CRegy] * e_CR_units_0;
        if(return_CRbin_CR_species_ID(k_CRegy) > 0)
        {
            double E_GeV = return_CRbin_kinetic_energy_in_GeV(target,k_CRegy), beta = return_CRbin_beta_factor(target,k_CRegy), Z=fabs(return_CRbin_CR_charge_in_e(target,k_CRegy));
            double T_eff_fullion = 0.59*(5./3.-1.)*U_TO_TEMP_UNITS*SphP[target].InternalEnergyPred, xm = 0.0286*sqrt(T_eff_fullion/2.e6);
            e_heat += b_coulomb_ion_per_GeV * ((Z*Z*beta*beta)/((beta*beta*beta+xm*xm*xm)*E_GeV)) * e_cr_units; // all protons Coulomb-heat, can be rapid for low-E
            if(E_GeV>=0.28) {e_heat += f_heat_hadronic * a_hadronic * e_cr_units;} // only GeV CRs or higher trigger above threshold for collisions
        }
    }
    return e_heat;
#else
    return (0.87*f_heat_hadronic*a_hadronic + 0.53*b_coulomb_ion_per_GeV) * Get_CosmicRayEnergyDensity_cgs(target) / nHcgs; /* for N<=2, assume a universal spectral shape, the factor here corrects for the fraction above-threshold for hadronic interactions, and 0.53 likewise for averaging  */
#endif
#elif defined(COOL_LOW_TEMPERATURES) // no CR module, but low-temperature cooling is on, we account for the CRs as a heating source, assuming a MW-like background scaled cosmologically to avoid over-heating IGM at high redshifts //
    double prefac_CR=1.; if(All.ComovingIntegrationOn) {double rhofac = (PROTONMASS_CGS*nHcgs/HYDROGEN_MASSFRAC) / (1000.*COSMIC_BARYON_DENSITY_CGS); if(rhofac < 0.2) {prefac_CR=0;} else {if(rhofac > 200.) {prefac_CR=1;} else {prefac_CR=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off CR heating for any gas with density unless it's >1000 times the cosmic mean density
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
#ifdef COSMIC_RAY_FLUID
    divv_f=-dt_entr*SphP[i].Face_DivVel_ForAdOps;
#endif
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

