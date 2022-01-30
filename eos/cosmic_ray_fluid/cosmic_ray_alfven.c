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

#if defined(CRFLUID_EVOLVE_SCATTERINGWAVES)

/*! To Do: with new RSOL scheme, needs some RSOL factors more carefully placed, here.
    generally can be updated to be a bit more flexible
    to handle newer damping rate estimates more modularly, and to deal with extinsic turbulence as well
    (which just appears as a source term in the e_A equations) */

/* routine to do the drift/kick operations for CRs: mode=0 is kick, mode=1 is drift */
double CosmicRay_Update_DriftKick(int i, double dt_entr, int mode)
{
    int k_CRegy;
    if(dt_entr <= 0) {return 0;} // no update
    for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++)
    {
        
    int k; double eCR, u0; k=0; if(mode==0) {eCR=SphP[i].CosmicRayEnergy[k_CRegy]; u0=SphP[i].InternalEnergy;} else {eCR=SphP[i].CosmicRayEnergyPred[k_CRegy]; u0=SphP[i].InternalEnergyPred;} // initial energy
    if(u0<All.MinEgySpec) {u0=All.MinEgySpec;} // enforced throughout code
    if(eCR < 0) {eCR=0;} // limit to physical values

    // now update all scalar fields (CR energies and Alfvenic energies, if those are followed) from fluxes and adiabatic terms //
    int q_whichupdate, q_N_updates = 3; // update the Alfvenic energy terms from (0) advection with gas [already solved], (1) their fluxes, and (2) their adiabatic terms. this should be basically identical to the CR term.
    for(q_whichupdate=0; q_whichupdate<q_N_updates; q_whichupdate++)
    {
        // first update the CR energies from fluxes. since this is positive-definite, some additional care is needed //
        double dCR_dt = SphP[i].DtCosmicRayEnergy[k_CRegy], gamma_eff = GAMMA_COSMICRAY(k_CRegy), eCR_tmp = eCR;
        if(q_whichupdate>0) {dCR_dt=SphP[i].DtCosmicRayAlfvenEnergy[k_CRegy][q_whichupdate-1]; gamma_eff=(3./2.); if(mode==0) {eCR_tmp=SphP[i].CosmicRayAlfvenEnergy[k_CRegy][q_whichupdate-1];} else {eCR_tmp=SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][q_whichupdate-1];}}
        double dCR = dCR_dt*dt_entr, dCRmax = 1.e10*(eCR_tmp+MIN_REAL_NUMBER);
        if(dCR > dCRmax) {dCR=dCRmax;} // don't allow excessively large values
        if(dCR < -eCR_tmp) {dCR=-eCR_tmp;} // don't allow it to go negative
	    double eCR_00 = eCR_tmp; eCR_tmp += dCR; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        if(q_whichupdate==0) {if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy]=eCR_tmp;} else {SphP[i].CosmicRayEnergyPred[k_CRegy]=eCR_tmp;}} // updated energy
        if(q_whichupdate>0) {if(mode==0) {SphP[i].CosmicRayAlfvenEnergy[k_CRegy][q_whichupdate-1]=eCR_tmp;} else {SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][q_whichupdate-1]=eCR_tmp;}} // updated energy
        // now need to account for the adiabatic heating/cooling of the 'fluid', here, with gamma=gamma_eff //
        double eCR_0 = eCR_tmp, d_div = (-(gamma_eff-1.) * SphP[i].Face_DivVel_ForAdOps*All.cf_a2inv) * dt_entr;
        if(All.ComovingIntegrationOn) {d_div += (-3.*(gamma_eff-1.) * All.cf_hubble_a) * dt_entr;} /* adiabatic term from Hubble expansion (needed for cosmological integrations */
        double dCR_div = DMIN(eCR_tmp*d_div , 0.5*u0*P[i].Mass); // limit so don't take away all the gas internal energy [to negative values]
        if(dCR_div + eCR_tmp < 0) {dCR_div = -eCR_tmp;} // check against energy going negative
        eCR_tmp += dCR_div; if((eCR_tmp<0)||(isnan(eCR_tmp))) {eCR_tmp=0;} // check against energy going negative or nan
        dCR_div = eCR_tmp - eCR_0; // actual change that is going to be applied
        if(dCR_div < -0.5*P[i].Mass*u0) {dCR_div=-0.5*P[i].Mass*u0;} // before re-coupling, ensure this will not cause negative energies
        if(dCR_div < -0.9*eCR_00) {dCR_div=-0.9*eCR_00;} // before re-coupling, ensure this will not cause negative energies 
        if(q_whichupdate==0) {if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy] += dCR_div; SphP[i].InternalEnergy -= dCR_div/P[i].Mass;} else {SphP[i].CosmicRayEnergyPred[k_CRegy] += dCR_div; SphP[i].InternalEnergyPred -= dCR_div/P[i].Mass;}}
        if(q_whichupdate>0) {if(mode==0) {SphP[i].CosmicRayAlfvenEnergy[k_CRegy][q_whichupdate-1] += dCR_div; SphP[i].InternalEnergy -= dCR_div/P[i].Mass;} else {SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][q_whichupdate-1] += dCR_div; SphP[i].InternalEnergyPred -= dCR_div/P[i].Mass;}}
    }

    int target_bin_centering_for_CR_quantities = i; // if this = i, evaluate quantities like R_GV at the CR-energy weighted mean of the bin, if =-1, evaluate them at the bin center instead: important for some subtle effects especially if using numerical derivatives for correction terms
    double E_CRs_Gev=return_CRbin_CR_rigidity_in_GV(target_bin_centering_for_CR_quantities,k_CRegy), Z_charge_CR=fabs(return_CRbin_CR_charge_in_e(i,k_CRegy)), M_cr_mp=return_CRbin_CRmass_in_mp(i,k_CRegy); // charge and energy and resonant Alfven wavenumber (in gyro units) of the CR population we're evolving
        
    // ok, the updates from [0] advection w gas, [1] fluxes, [2] adiabatic, [-] catastrophic (in cooling.c) are all set, just need exchange terms b/t CR and Alfven //
    double EPSILON_SMALL = 1.e-77; // want a very small number here
    double bhat[3], Bmag=0, Bmag_Gauss, clight_code=C_LIGHT_CODE, Omega_gyro, eA[2], vA_code, vA2_c2, E_B, fac, flux_G, fac_Omega, flux[3], f_CR, f_CR_dot_B, cs_thermal, r_turb_driving, G_ion_neutral=0, G_turb_plus_linear_landau=0, G_nonlinear_landau_prefix=0;
    double ne=1, f_ion=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, rho_cgs=rho*UNIT_DENSITY_IN_CGS;
#ifdef COOLING 
    temperature = ThermalProperties(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
    f_ion = DMIN(DMAX(DMAX(DMAX(1-nh0, nhp), ne/1.2), 1.e-8), 1.); // account for different measures above (assuming primordial composition)
#endif
    for(k=0;k<3;k++) {if(mode==0) {bhat[k]=SphP[i].B[k];} else {bhat[k]=SphP[i].BPred[k];}} // grab whichever B field we need for our mode
    if(mode==0) {eCR=SphP[i].CosmicRayEnergy[k_CRegy]; u0=SphP[i].InternalEnergy;} else {eCR=SphP[i].CosmicRayEnergyPred[k_CRegy]; u0=SphP[i].InternalEnergyPred;} // initial energy
    if(u0<All.MinEgySpec) {u0=All.MinEgySpec;} // enforce the usual minimum thermal energy the code requires
    for(k=0;k<2;k++) {if(mode==0) {eA[k]=SphP[i].CosmicRayAlfvenEnergy[k_CRegy][k];} else {eA[k]=SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][k];}} // Alfven energy
    if(mode==0) {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFlux[k_CRegy][k];}} else {for(k=0;k<3;k++) {flux[k]=SphP[i].CosmicRayFluxPred[k_CRegy][k];}} // load flux
    f_CR=0; f_CR_dot_B=0; for(k=0;k<3;k++) {f_CR+=flux[k]*flux[k]; f_CR_dot_B+=bhat[k]*flux[k];} // compute the magnitude of the flux density
    f_CR=sqrt(f_CR); if(f_CR_dot_B<0) {f_CR*=-1;} // initialize the flux density variable from the previous timestep, appropriately signed with respect to the b-field
    for(k=0;k<3;k++) {Bmag+=bhat[k]*bhat[k];} // compute magnitude
    Bmag = sqrt(Bmag); for(k=0;k<3;k++) {bhat[k]/=(EPSILON_SMALL+Bmag);} // now it's bhat we have here
    Bmag *= SphP[i].Density/P[i].Mass * All.cf_a2inv; // convert to actual B in physical units
    E_B = 0.5*Bmag*Bmag * (P[i].Mass/(SphP[i].Density*All.cf_a3inv)); // B-field energy (energy density times volume, for ratios with energies above)
    double Eth_0 = EPSILON_SMALL + 1.e-8 * P[i].Mass*u0; // set minimum magnetic energy relative to thermal (maximum plasma beta ~ 1e8) to prevent nasty divergences
    if(E_B < Eth_0) {Bmag = sqrt(2.*Eth_0/((P[i].Mass/(SphP[i].Density*All.cf_a3inv))));} // enforce this maximum beta for purposes of "B" to insert below 
    E_B = 0.5*Bmag*Bmag * (P[i].Mass/(SphP[i].Density*All.cf_a3inv)); // B-field energy (energy density times volume, for ratios with energies above)
    Bmag_Gauss = Bmag * UNIT_B_IN_GAUSS; // turn it into Gauss
    Omega_gyro = (8987.34 * Bmag_Gauss * (Z_charge_CR/E_CRs_Gev)) * UNIT_TIME_IN_CGS; // gyro frequency of the CR population we're evolving, converted to physical code units //
    double vA_noion = Get_Gas_Alfven_speed_i(i); // Alfven speed in code units [recall B units such that there is no 4pi here]
    vA_code = Get_Gas_ion_Alfven_speed_i(i); // include ionization appropriately for small-scale modes
    cs_thermal = sqrt(convert_internalenergy_soundspeed2(i,u0)); // thermal sound speed at appropriate drift-time [in code units, physical]
    vA2_c2 = vA_code*vA_code / (clight_code*clight_code); // Alfven speed vs speed of light
    fac_Omega = (3.*M_PI/16.) * Omega_gyro * (1.+2.*vA2_c2); // factor which will be used heavily below
    /* for turbulent (anisotropic and linear landau) damping terms: need to know the turbulent driving scale: assume a cascade with a driving length equal to the pressure gradient scale length */
    r_turb_driving = 0; for(k=0;k<3;k++) {r_turb_driving += SphP[i].Gradients.Pressure[k]*SphP[i].Gradients.Pressure[k];} // compute gradient magnitude
    r_turb_driving = DMAX( SphP[i].Pressure / (EPSILON_SMALL + sqrt(r_turb_driving)) , Get_Particle_Size(i) ) * All.cf_atime; // maximum of gradient scale length or resolution scale
    double k_turb = 1./r_turb_driving, k_L = Omega_gyro / clight_code;
    
    // before acting on the 'stiff' sub-system, account for the 'extra' advection term that accounts for 'twisting' of B:
    fac=0; for(k=0;k<3;k++) {fac += All.cf_a2inv * bhat[k] * (bhat[0]*SphP[i].Gradients.Velocity[k][0] + bhat[1]*SphP[i].Gradients.Velocity[k][1] + bhat[2]*SphP[i].Gradients.Velocity[k][2]);}
    if(All.ComovingIntegrationOn) {fac += All.cf_hubble_a;} // adds cosmological/hubble flow term here [not included in peculiar velocity gradient]
    fac *= -dt_entr; if(!isfinite(fac)) {fac=0;} else {if(fac>2.) {fac=2.;} else {if(fac<-2.) {fac=-2.;}}} // limit factor for change here, should be small given Courant factor
    f_CR *= exp(fac); // update flux term accordingly, before next step //

    // because the equations below will very much try to take things to far-too-small values for numerical precision, we need to define a bunch of sensible bounds for values to allow, to prevent divergences, but also enforce conservation
    // calculate minimum eA,eCR to enforce; needed because if eA is identically zero, nothing can get amplified, and it will always be zero. but for large enough seed to amplify, results should not depend on seed //
    eA[0]=DMAX(eA[0],0); eA[1]=DMAX(eA[1],0); eCR=DMAX(eCR,0); // enforce non-negative energies 
    double Min_Egy=0, e_tot=0, e_tot_new=0, fmax=0; e_tot = eCR + eA[0] + eA[1] + EPSILON_SMALL; // sum total energy, enforce positive-definite: will use this to ensure total energy conservation when enforcing minima below
    { 
        double h=Get_Particle_Size(i)*All.cf_atime; int k2; for(k=0;k<3;k++) {for(k2=0;k2<3;k2++) {Min_Egy+=SphP[i].Gradients.B[k][k2]*SphP[i].Gradients.B[k][k2];}}
        Min_Egy=h*sqrt(Min_Egy/9.)*All.cf_a2inv; Min_Egy=DMIN(Min_Egy,Bmag); r_turb_driving=DMAX(h,r_turb_driving); Min_Egy=DMIN(Min_Egy,Bmag*pow(h/r_turb_driving,1./3.)); Min_Egy=Min_Egy*pow(DMIN(clight_code/Omega_gyro,DMIN(h,r_turb_driving))/h,1./3.); // Min_Egy is now magnetic field extrap to r_gyro
        Min_Egy = 0.5 * (Min_Egy*Min_Egy) * P[i].Mass/(SphP[i].Density*All.cf_a3inv); // magnetic energy at this scale, from the above //
        double epsilon = 1.e-15; Min_Egy *= epsilon; // minimum energy is a tiny fraction of B at the dissipation scale
        if(Min_Egy <= 0 || !isfinite(Min_Egy)) {Min_Egy = 1.e-15*eCR;} // if this minimum-energy calculation failed, enforce a tiny fraction of the CR energy
        if(Min_Egy <= 0 || !isfinite(Min_Egy)) {Min_Egy = 1.e-15*P[i].Mass*u0;} // if this minimum-energy calculation failed, enforce a tiny fraction of the thermal energy
        if(Min_Egy <= 0) {Min_Egy = EPSILON_SMALL;} // if this still failed, simply enforce a tiny positive-definite value
    }
    eCR=DMAX(eCR,Min_Egy); eA[0]=DMAX(eA[0],Min_Egy); eA[1]=DMAX(eA[1],Min_Egy); // enforce

    // ok, now all the advection and adiabatic operations should be complete. they are split above. 
    //  what remains is the stiff, coupled subsystem of wave growth+damping, which needs to be treated 
    //  more carefully or else we get very large over/under-shoots

    // first define some convenient units and dimensionless quantities, and enforce limits on values of input quantities
    double cr_speed = CRFLUID_REDUCED_C_CODE(k_CRegy);
    double eCR_0 = 1.e-6*(E_B + P[i].Mass*u0) + eCR + eA[0] + eA[1] + fabs(f_CR/cr_speed); // this can be anything, just need a normalization for the characteristic energy scale of the problem //
    double ceff2_va2=(cr_speed*cr_speed)/(vA_code*vA_code), t0=1./(fac_Omega*(eCR_0/E_B)*vA2_c2), gammCR=GAMMA_COSMICRAY(k_CRegy), f_unit=vA_code*eCR_0, volume=P[i].Mass/(SphP[i].Density*All.cf_a3inv); // factors used below , and for units
    double x_e=eCR/eCR_0, x_f=f_CR/f_unit, x_up=eA[0]/eCR_0, x_um=eA[1]/eCR_0, dtau=dt_entr/t0; e_tot/=eCR_0; Min_Egy/=eCR_0; // initial values in relevant units
    Min_Egy=DMAX(DMIN(Min_Egy,DMIN(x_e,DMIN(x_up,x_um))),EPSILON_SMALL); if(!isfinite(Min_Egy)) {Min_Egy=EPSILON_SMALL;} // enforce positive-definite-ness
    // we can more robustly define a minimum and maximum e_A by reference to a minimum and maximum 'effective diffusivity' over which it is physically meaningful, and numerically possible to evolve them
    double ref_diffusivity = 4.4e26 / (UNIT_VEL_IN_CGS * UNIT_LENGTH_IN_CGS); // define a unit diffusivity in code units for reference below
    double xkappa_min = DMAX(vA_code*vA_code*t0/(3.e8*ref_diffusivity) , EPSILON_SMALL); // maximum diffusivity ~1e35, but be non-zero
    double xkappa_max = DMAX(DMIN(vA_code*vA_code*t0/(3.e-8*ref_diffusivity) , 0.5*E_B/eCR_0), xkappa_min); // minimum diffusivity at ~1e19, but cannot have more energy in eAp+eAm than total magnetic energy! (equations below assume -small- fraction of E_B in eA!, or growth rates non-linearly modified)
    if(e_tot < Min_Egy || !isfinite(e_tot)) {e_tot = Min_Egy;} // enforce minima/maxima
    if(x_e   < Min_Egy || !isfinite(x_e)  ) {x_e   = Min_Egy;} // enforce minima/maxima
    if(x_um<EPSILON_SMALL || !isfinite(x_um)) {x_um=EPSILON_SMALL;} else {if(x_um>xkappa_max) {x_um=xkappa_max;}} // enforce minima/maxima
    if(x_up<EPSILON_SMALL || !isfinite(x_up)) {x_up=EPSILON_SMALL;} else {if(x_up>xkappa_max) {x_up=xkappa_max;}} // enforce minima/maxima
    if(x_um+x_up<xkappa_min) {fac=xkappa_min/(x_um+x_up); x_um*=fac; x_up*=fac;} // only want to enforce -sum- having effective diffusivity, not both
    e_tot_new=x_e+x_um+x_up; x_e*=e_tot/e_tot_new; x_up*=e_tot/e_tot_new; x_um*=e_tot/e_tot_new; // check energy after limit-enforcement
    fmax = x_e*sqrt(ceff2_va2); if(!isfinite(x_f)) {x_f=0;} else {if(x_f>fmax) {x_f=fmax;} else {if(x_f<-fmax) {x_f=-fmax;}}} // check for flux maximum/minimum

    // calculate the dimensionless flux source term for the stiff part of the equations
    flux_G=0; for(k=0;k<3;k++) {flux_G += bhat[k] * SphP[i].Gradients.CosmicRayPressure[k_CRegy][k];} // b.gradient[P] -- flux source term
    double psifac = flux_G * (vA_code*t0) / (eCR/volume); // this gives the strength of the gradient source term, should remain fixed over stiff part of loop

    // calculate the wave-damping rates (again in appropriate dimensionless units)
    /* ion-neutral damping: need thermodynamic information (neutral fractions, etc) to compute self-consistently */
    G_ion_neutral = (5.77e-11 * (rho_cgs/PROTONMASS_CGS) * nh0 * sqrt(temperature)) * UNIT_TIME_IN_CGS / sqrt(M_cr_mp); // need to get thermodynamic quantities [neutral fraction, temperature in Kelvin] to compute here -- // G_ion_neutral = (xiH + xiHe); // xiH = nH * siH * sqrt[(32/9pi) *kB*T*mH/(mi*(mi+mH))]. converted to -physical- code units

    int i1,i2; double v2_t=0,dv2_t=0,b2_t=0,db2_t=0,x_LL,M_A,h0,fturb_multiplier=1; // factor which will represent which cascade model we are going to use
    for(i1=0;i1<3;i1++)
    {
        v2_t += SphP[i].VelPred[i1]*SphP[i].VelPred[i1]; b2_t += Get_Gas_BField(i,i1) * Get_Gas_BField(i,i1);
        for(i2=0;i2<3;i2++) {dv2_t += SphP[i].Gradients.Velocity[i1][i2]*SphP[i].Gradients.Velocity[i1][i2]; db2_t += SphP[i].Gradients.B[i1][i2]*SphP[i].Gradients.B[i1][i2];}
    }
    v2_t=sqrt(v2_t); b2_t=sqrt(b2_t); dv2_t=sqrt(dv2_t); db2_t=sqrt(db2_t); dv2_t/=All.cf_atime; db2_t/=All.cf_atime; b2_t*=All.cf_a2inv; db2_t*=All.cf_a2inv; v2_t/=All.cf_atime; dv2_t/=All.cf_atime; h0=Get_Particle_Size(i)*All.cf_atime; // physical units
    M_A = h0*(EPSILON_SMALL + dv2_t) / (EPSILON_SMALL + vA_noion); M_A = DMAX(M_A , h0*(EPSILON_SMALL + db2_t) / (EPSILON_SMALL + b2_t)); M_A = DMAX( EPSILON_SMALL , M_A ); // proper calculation of the local Alfven Mach number
    x_LL = clight_code / (Omega_gyro * h0); x_LL=DMAX(x_LL,EPSILON_SMALL); k_turb = 1./h0; // scale at which turbulence is being measured here //
    fturb_multiplier = pow(M_A,3./2.); // corrects to Alfven scale, for correct estimate according to Farmer and Goldreich, Lazarian, etc.
    if(M_A<1.) {fturb_multiplier*=DMIN(sqrt(M_A),pow(M_A,7./6.)/pow(x_LL,1./6.));} else {fturb_multiplier*=DMIN(1.,1./(pow(M_A,1./2.)*pow(x_LL,1./6.)));} /* Lazarian+ 2016 multi-part model for where the resolved scales lie on the cascade */
    G_turb_plus_linear_landau = (vA_noion + sqrt(M_PI/16.)*cs_thermal) * sqrt(k_turb*k_L) * fturb_multiplier; // linear Landau + turbulent (both have same form, assume k_turb from cascade above)

    G_nonlinear_landau_prefix = (sqrt(M_PI)/8.) * (1./E_B) * (cs_thermal*k_L); // non-linear Landau damping (will be multiplied by eA)
    double gamma_in_t_ll = (G_ion_neutral + G_turb_plus_linear_landau) * t0; // dimensionless now and appropriate code units
    double gamma_nll = G_nonlinear_landau_prefix * eCR_0 * t0; // dimensionless now and appropriate code units


    // now we are ready to actually integrate these equations, in a numerically-stable manner, with protection from over/under-shooting
    double dtau_max = 1.e-5;
    double dx_e,dx_f,dx_up,dx_um,x_e_0=x_e,x_f_0=x_f,x_up_0=x_up,x_um_0=x_um,dtaux=0.,efmax=50.,expfac;
    double x_e_prev,x_f_prev,x_up_prev,x_um_prev,n_eqm_loops=1.; fmax=1./EPSILON_SMALL; // (need to set initial fmax to large value)
    long n_iter=0, n_iter_max=100000; // sets the maximum number of sub-cycles which we will allow below for any sub-process
    while(1)
    {
        /* here's the actual set of remaining stiff equations to be solved
            dx_e  = gammCR*(x_up+x_um)*x_e + (x_um-x_up)*x_f;                     // deCR_dt
            dx_f  = -ceff2_va2*(psifac + (x_um-x_up)*x_e + (x_up+x_um)*x_f);      // df_dt
            dx_up = -x_up*(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_up - x_f);    // deAp_dt
            dx_um = -x_um*(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_um + x_f);    // deAm_dt
        */
        x_e_prev=x_e; x_f_prev=x_f; x_up_prev=x_up; x_um_prev=x_um; // reset values at the beginning of the loop (these will be cycled multiple times below)

        // for eqm: if psi>0, f<0, um->grows, up->pure-damping //
        double f_eqm, up_eqm, um_eqm, tinv_u, tinv_f;
        double q_tmp = (gammCR-1.)*x_e + gamma_in_t_ll, q_inner = (4.*gamma_nll*fabs(psifac)) / (q_tmp*q_tmp);
        if(q_inner < 1.e-4) {q_inner=q_inner/2.;} else {q_inner=sqrt(1.+q_inner)-1.;}
        double x_nonzero = q_tmp * q_inner / (2.*gamma_nll); if(fabs(gamma_nll) < EPSILON_SMALL) {x_nonzero = fabs(psifac)/q_tmp;}    
        double x_f_magnitude = x_e + fabs(psifac) / x_nonzero;
        if(psifac > 0)
        {
            up_eqm=xkappa_min; um_eqm=x_nonzero; f_eqm=-x_f_magnitude;
            tinv_u = EPSILON_SMALL + fabs(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_um + x_f);
        } else {
            um_eqm=xkappa_min; up_eqm=x_nonzero; f_eqm=+x_f_magnitude;
            tinv_u = EPSILON_SMALL + fabs(gammCR*x_e + gamma_in_t_ll + gamma_nll*x_up - x_f);
        }
        tinv_f = fabs( ceff2_va2*(psifac + (x_um-x_up)*x_e + (x_up+x_um)*x_f) ) * (1./(EPSILON_SMALL + fabs(f_eqm)) + 1./(EPSILON_SMALL+fabs(x_f)));
        double t_eqm = 1./(tinv_u + tinv_f); // timescale to approach equilibrium solution
        
        // set timestep (steadily  growing from initial conservative value ) //
        dtaux = dtau; if(dtaux > dtau) {dtaux=dtau;}
        if(dtaux > dtau_max) {dtaux = dtau_max;}
        dtau_max *= 2.; if(dtaux > 10.) {dtaux=10.;}
        
        double jump_fac = 0.5; // fraction towards equilibrium to 'jump' each time
        //if(dtaux >= 0.33*jump_fac*t_eqm)
        if(dtau >= jump_fac*t_eqm)
        {
            // timestep is larger than the timescale to approach the equilibrium solution, 
            //  so move the systems towards equilibrium, strictly
            //
            if(dtaux > jump_fac*t_eqm) {dtaux = jump_fac*t_eqm;} else {jump_fac = dtaux/t_eqm;} // initial 'step' is small fraction of equilibrium
            dtaux = n_eqm_loops * t_eqm; jump_fac = 1. + (jump_fac-1.)/n_eqm_loops; n_eqm_loops*=1.1; // each sub-cycle consecutively in eqm, allow longer step
            if((x_f<=-fmax && f_eqm<=-fmax) || (x_f>=+fmax && f_eqm>=+fmax)) {t_eqm=1./EPSILON_SMALL; jump_fac=1.; dtaux=dtau;} // if slamming into limits, terminate cycle with big step
            if(f_eqm > 0)
            {
                x_up = exp( log(x_up)*(1.-jump_fac) + log(up_eqm)*jump_fac );
                if(x_f > 0) {x_f = +exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );} else {
                    if(fabs(x_f)<10.*fabs(f_eqm)) {x_f=x_f*(1.-jump_fac)+f_eqm*jump_fac;} else {
                        x_f = -exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );}}
            } else {
                x_um = exp( log(x_um)*(1.-jump_fac) + log(um_eqm)*jump_fac );
                if(x_f < 0) {x_f = -exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );} else {
                    if(fabs(x_f)<10.*fabs(f_eqm)) {x_f=x_f*(1.-jump_fac)+f_eqm*jump_fac;} else {
                        x_f = +exp( log(fabs(x_f))*(1.-jump_fac) + log(fabs(f_eqm))*jump_fac );}}
            }

        } else {

            // timestep is smaller than the timescale to approach equilibrium, so integrate directly, 
            //  but we will still use a fully implicit backwards-Euler type scheme for the two 'stiffest' 
            //  components of the system (namely, the flux and eA term corresponding to the multiplicative direction)
            //  [the other terms, e.g. the damped energy change and the CR energy change, can be dealt with after]
            //
            double x_dum=0, x_out=0; n_eqm_loops=1.; // (if we enter this,  we need to terminate the parent loop above)
            if(f_eqm>0) {x_dum=x_um;} else {x_dum=x_up;}
            double q0 = 1.+ceff2_va2*dtaux*x_dum, g00 = gammCR*x_e + gamma_in_t_ll, psi00 = psifac + x_dum*x_e;
            double a_m1 = -x_up_prev/dtaux, a_0 = g00 + 1./dtaux , a_1 = gamma_nll, c2dt = ceff2_va2*dtaux;
            if(f_eqm<0) {psi00=psifac-x_dum*x_e; a_1=-gamma_nll; a_0=-(g00 + 1./dtaux); a_m1=x_um_prev/dtaux;}
            double d0 = -a_m1*q0, c0 = x_f_prev - a_0*q0 - c2dt*(a_m1 + psi00), b0 = -a_1*q0 - c2dt*(a_0-x_e), a0 = -a_1*c2dt; 
            if(f_eqm<0) {b0 = -a_1*q0 - c2dt*(a_0+x_e);;}
            if(fabs(a0) < EPSILON_SMALL)
            {
                if(fabs(b0) < EPSILON_SMALL)
                {
                    x_out = fabs(d0/c0); // linear solve
                } else {
                    d0/=c0; b0/=c0; c0=fabs(4.*b0*d0); if(c0<1.e-4) {c0*=0.5;} else {c0=sqrt(1.+c0)-1.;}
                    x_out = c0/(2.*fabs(b0)); // quadratic solve
                }    
            } else {
                // cubic solve
                double p0=-b0/(3.*a0), q0=p0*p0*p0 + (b0*c0-3.*a0*d0)/(6.*a0*a0), r0=c0/(3.*a0), f0=r0-p0*p0, g0=q0*q0 + f0*f0*f0; 
                if(g0 >= 0.)
                {
                    g0=sqrt(g0); a0=q0+g0; b0=q0-g0; 
                    q0=pow(fabs(a0),1./3.); if(a0<0) {q0*=-1.;}
                    r0=pow(fabs(b0),1./3.); if(b0<0) {r0*=-1.;}
                } else {
                    g0=sqrt(-g0); a0=sqrt(q0*q0+g0*g0); b0=atan(g0/q0); r0=0.; q0=2.*a0*cos(b0/3.);
                }
                x_out = fabs(p0 + q0 + r0);
            }
            x_f = a_m1/x_out + a_0 + a_1*x_out;
            if(f_eqm>0) {x_up=x_out;} else {x_um=x_out;}

        }

        // now deal with the non-stiff part of the equations, namely the other Alfven-energy component + CR energy
        if(f_eqm > 0) // do the evolution for the eA term -not- involved in the stiff part of the equations //
        {
            double g0 = gammCR*x_e + gamma_in_t_ll + x_f; // pure-damping for um
            if(g0 > 0) {
                expfac=g0*dtaux; if(expfac>efmax) {expfac=efmax;} 
                if(expfac>1.e-6) {expfac=exp(expfac)-1.;}
                x_um /= (1. + (1.+gamma_nll*x_um/g0)*expfac);
            } else {x_um -= x_um*(g0 + gamma_nll*x_up)*dtaux;} // (linear if x_f hasn't behaved yet)
        } else {
            double g0 = gammCR*x_e + gamma_in_t_ll - x_f; // pure-damping for up
            if(g0 > 0) {
                expfac=g0*dtaux; if(expfac>efmax) {expfac=efmax;} 
                if(expfac>1.e-6) {expfac=exp(expfac)-1.;}
                x_up /= (1. + (1.+gamma_nll*x_up/g0)*expfac);
            } else {x_up -= x_up*(g0 + gamma_nll*x_up)*dtaux;} // (linear if x_f hasn't behaved yet)
        }
        
        // calculate total-energy damping (needed for deriving change in e_cr, which is then given by energy conservation) //
        double x_um_eff=0.5*(x_um+x_um_prev), x_up_eff=0.5*(x_up+x_up_prev), x_f_eff=0.5*(x_f+x_f_prev); // effective values for use in damping rates below
        expfac=gamma_in_t_ll*dtaux; if(expfac>efmax) {expfac=efmax;} 
        if(expfac>1.e-6) {expfac=exp(expfac)-1.;}
        double de_damp = x_um_eff/(1.+1./(expfac*(1.+gamma_nll*x_um_eff/gamma_in_t_ll))) + 
                         x_up_eff/(1.+1./(expfac*(1.+gamma_nll*x_up_eff/gamma_in_t_ll))); // energy loss to thermalized wave-damping
        if(!isfinite(de_damp)) {de_damp=0;}
        double e_tot = DMAX(x_um_prev,0) + DMAX(x_up_prev,0) + DMAX(x_e_prev,0) - de_damp; // total energy (less damping) before step
        double x_e_egycon = DMAX(e_tot-(x_up+x_um), Min_Egy);
        expfac = gammCR*(x_up_eff+x_um_eff)*dtaux; double x_numer = x_e_prev + dtaux*(x_um_eff-x_up_eff)*x_f_eff;
        if(expfac<0.9 && x_numer>0.) {x_e=x_numer/(1.-expfac);} else {
            if(expfac*x_e_prev+x_numer > 0.01*x_e) {x_e=expfac*x_e_prev+x_numer;} else {x_e*=0.01;}}
        if(fabs(x_e_egycon-x_e_prev) < fabs(x_e-x_e_prev)) {x_e=x_e_egycon;}
        if(e_tot < Min_Egy || !isfinite(e_tot)) {e_tot = Min_Egy;} // enforce minima/maxima
        if(x_e   < Min_Egy || !isfinite(x_e)  ) {x_e   = Min_Egy;} // enforce minima/maxima
        if(x_um<EPSILON_SMALL || !isfinite(x_um)) {x_um=EPSILON_SMALL;} else {if(x_um>xkappa_max) {x_um=xkappa_max;}} // enforce minima/maxima
        if(x_up<EPSILON_SMALL || !isfinite(x_up)) {x_up=EPSILON_SMALL;} else {if(x_up>xkappa_max) {x_up=xkappa_max;}} // enforce minima/maxima
        if(x_um+x_up<xkappa_min) {expfac=xkappa_min/(x_um+x_up); x_um*=expfac; x_up*=expfac;} // only want to enforce -sum- having effective diffusivity, not both
        double e_tot_new=x_e+x_um+x_up; x_e*=e_tot/e_tot_new; x_up*=e_tot/e_tot_new; x_um*=e_tot/e_tot_new; // check energy after limit-enforcement
	    fmax = x_e*sqrt(ceff2_va2); if(!isfinite(x_f)) {x_f=0;} else {if(x_f>fmax) {x_f=fmax;} else {if(x_f<-fmax) {x_f=-fmax;}}} // check for flux maximum/minimum
        
        // calculate change in parameters to potentially break the cycle here
        dx_e  = (x_e - x_e_prev) / (EPSILON_SMALL + x_e + x_e_prev);
        dx_up = (x_up - x_up_prev) / (EPSILON_SMALL + x_up + x_up_prev);
        dx_um = (x_um - x_um_prev) / (EPSILON_SMALL + x_um + x_um_prev);
        dx_f  = (x_f - x_f_prev) / (EPSILON_SMALL + fabs(x_f) + fabs(x_f_prev));
        double dx_max = sqrt(dx_e*dx_e + dx_up*dx_up + dx_um*dx_um + dx_f*dx_f); // sum in quadrature
        if(!isfinite(dx_max)) {dx_max=1;} // enforce validity for check below 

        if((n_iter > 0) && (n_iter % 10000 == 0)) // print diagnostics if the convergence is happening slowly
        {
            printf("niter/max=%ld/%ld dtau/step=%g/%g d_params=%g (init/previous/now) xeCR=%g/%g/%g xeAp=%g/%g/%g xeAm=%g/%g/%g xflux=%g/%g/%g ceff2_va2=%g damp_g_intll=%g damp_g_nll=%g psi_gradientfac=%g heat_term=%g xkappa_min/max=%g/%g egy_min=%g flux_max=%g \n",
                n_iter,n_iter_max,dtau,dtaux,dx_max,x_e_0,x_e_prev,x_e,x_up_0,x_up_prev,x_up,x_um_0,x_um_prev,x_um,x_f_0,x_f_prev,x_f,ceff2_va2,gamma_in_t_ll,gamma_nll,psifac,(x_e_0+x_up_0+x_um_0)-(x_e+x_up+x_um),xkappa_min,xkappa_max,Min_Egy,fmax);
            fflush(stdout);
        }
        dtau -= dtaux; // subtract the time we've already integrated from the total timestep
        n_iter++; // count the number of iterations
        if(dtau <= 0.) break; // we have reached the end of the integration time for our sub-stepping. we are done!
        if(dx_max <= 1.e-3*dtaux/dtau) break; // the values of -all- the parameters are changing by much less than the floating-point errors. we are done!
        if(n_iter > n_iter_max) break; // we have reached the maximum allowed number of iterations. we give up!
    }

    // ok! done with the main integration/sub-cycle loop, now just do various clean-up operations
    double thermal_heating = eCR_0 * ((x_e_0+x_up_0+x_um_0)-(x_e+x_up+x_um)); // net thermalized energy from damping terms
    if(thermal_heating < 0 || !isfinite(thermal_heating)) // if this is less than zero (from residual floating-point error), then the energy goes up, which shouldnt happen: set to zero
    {
        e_tot=x_e_0+x_up_0+x_um_0; e_tot_new = x_e+x_up+x_um; // initial and final energies should be equal in this case
        x_e *= e_tot/e_tot_new; x_up *= e_tot/e_tot_new; x_um *= e_tot/e_tot_new; // enforce that equality
        thermal_heating=0; // set thermal change to nil
    }
    eCR=eCR_0*x_e; eA[0]=eCR_0*x_up; eA[1]=eCR_0*x_um; f_CR=f_unit*x_f; Min_Egy*=eCR_0; xkappa_min*=eCR_0; xkappa_max*=eCR_0; // re-assign dimensional quantities

    // assign the updated values back to the resolution elements, finally!
    if(mode==0) {SphP[i].CosmicRayEnergy[k_CRegy]=eCR;} else {SphP[i].CosmicRayEnergyPred[k_CRegy]=eCR;} // CR energy
    for(k=0;k<2;k++) {if(mode==0) {SphP[i].CosmicRayAlfvenEnergy[k_CRegy][k]=eA[k];} else {SphP[i].CosmicRayAlfvenEnergyPred[k_CRegy][k]=eA[k];}} // Alfven energy
    if(mode==0) {for(k=0;k<3;k++) {SphP[i].CosmicRayFlux[k_CRegy][k]=f_CR*bhat[k];}} else {for(k=0;k<3;k++) {SphP[i].CosmicRayFluxPred[k_CRegy][k]=f_CR*bhat[k];}} // assign to flux vector
    if(mode==0) {SphP[i].InternalEnergy+=thermal_heating/P[i].Mass;} else {SphP[i].InternalEnergyPred+=thermal_heating/P[i].Mass;} // heating term from damping
    SphP[i].CosmicRayDiffusionCoeff[k_CRegy] = 1. / (fac_Omega*((eA[0]+eA[1])/E_B)/(clight_code*clight_code)); // effective diffusion coefficient in code units

        
    } // complete loop over CR bins
    return 1; // exit
}




#endif // closes block for entire file


