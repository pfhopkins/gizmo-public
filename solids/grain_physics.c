#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*

 This module contains the self-contained sub-routines needed for
 grain-specific physics in proto-planetary/proto-stellar/planetary cases,
 GMC and ISM/CGM/IGM dust dynamics, dust dynamics in cool-star atmospheres,
 winds, and SNe remnants, as well as terrestrial turbulence and
 particulate-laden turbulence. Anywhere where particles coupled to gas
 via coulomb, aerodynamic, or lorentz forces are interesting.

 This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.

 */



#ifdef GRAIN_FLUID

#if defined(GRAIN_BACKREACTION)
void grain_backrx(void);
#endif

/* function to apply the drag on the grains from surrounding gas properties */
void apply_grain_dragforce(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    int i, k; PRINT_STATUS("Beginning particulate/grain/PIC force evaluation.");
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) /* loop over active particles */
    {
        if(!((1 << P[i].Type) & (GRAIN_PTYPES))) {P[i].Grain_AccelTimeMin = MAX_REAL_NUMBER;} /* for active elements, set this large to re-set below */
#ifdef BOX_BND_PARTICLES
        if(P[i].ID > 0) /* 'frozen' particles are excluded */
#endif
        if(((1 << P[i].Type) & (GRAIN_PTYPES)) && (P[i].Mass>0)) /* only particles of designated type[s] are eligible for this routine */
        {
#if defined(GRAIN_BACKREACTION)
            for(k=0;k<3;k++) {P[i].Grain_DeltaMomentum[k]=0;} /* reset momentum to couple back to gas (or else would diverge) */
#endif
            double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
            double vgas_mag = 0.0; for(k=0;k<3;k++) {vgas_mag+=(P[i].Gas_Velocity[k]-P[i].Vel[k])*(P[i].Gas_Velocity[k]-P[i].Vel[k]);}
            vgas_mag = sqrt(vgas_mag) / All.cf_atime; /* convert to physical units */
            int grain_subtype = 1; /* default assumption about particulate sub-type for operations below */
#if defined(PIC_MHD)
            grain_subtype = P[i].MHD_PIC_SubType;
#endif
            if((grain_subtype <= 2) && (dt > 0) && (P[i].Gas_Density>0) && (vgas_mag > 0)) /* only bother with particles moving wrt gas with finite gas density and timestep */
            {
                double gamma_eff = GAMMA_DEFAULT; // adiabatic index to use below
                double cs = sqrt( (gamma_eff*(gamma_eff-1)) * P[i].Gas_InternalEnergy);
                double R_grain_cgs = P[i].Grain_Size, R_grain_code = R_grain_cgs / UNIT_LENGTH_IN_CGS;
                double rho_gas = P[i].Gas_Density * All.cf_a3inv, rho_grain_physical = All.Grain_Internal_Density, rho_grain_code = rho_grain_physical / UNIT_DENSITY_IN_CGS; // rho_grain in cgs and code units //
                double x0 = 0.469993*sqrt(gamma_eff) * vgas_mag/cs; // (3/8)*sqrt[pi/2]*|vgas-vgrain|/cs //
                double tstop_inv = 1.59577/sqrt(gamma_eff) * rho_gas * cs / (R_grain_code * rho_grain_code); // 2*sqrt[2/pi] * 1/tstop //
#ifdef GRAIN_LORENTZFORCE /* calculate the grain charge following Draine & Sutin */
                double cs_cgs = cs * UNIT_VEL_IN_CGS;
                double tau_draine_sutin = R_grain_cgs * (2.3*PROTONMASS_CGS) * (cs_cgs*cs_cgs) / (gamma_eff * ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS);
                double Z_grain = -DMAX( 1./(1. + sqrt(1.0e-3/tau_draine_sutin)) , 2.5*tau_draine_sutin ); /* note: if grains moving super-sonically with respect to gas, and charge equilibration time is much shorter than the streaming/dynamical timescales, then the charge is slightly reduced, because the ion collision rate is increased while the electron collision rate is increased less (since electrons are moving much faster, we assume the grain is still sub-sonic relative to the electron sound speed. in this case, for the large-grain limit, the Draine & Sutin results can be generalized; the full expressions are messy but can be -approximated- fairly well for Mach numbers ~3-30 by simply suppressing the equilibrium grain charge by a power ~exp[-0.04*mach]  (weak effect, though can be significant for mach>10) */
                if(isnan(Z_grain)||(Z_grain>=0)) {Z_grain=0;}
#endif
#ifdef GRAIN_EPSTEIN_STOKES
                if(grain_subtype == 0 || grain_subtype == 1)
                {
                    double mu = 2.3*PROTONMASS_CGS, temperature = (mu/PROTONMASS_CGS) * (1.4-1.) * U_TO_TEMP_UNITS * P[i].Gas_InternalEnergy; // assume molecular gas (as its the only regime where this is relevant) with gamma=1.4
                    double cross_section = GRAIN_EPSTEIN_STOKES * 2.0e-15 * (1. + 70./temperature);
                    cross_section /= UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS;
                    double n_mol = rho_gas * UNIT_MASS_IN_CGS / mu, mean_free_path = 1 / (n_mol * cross_section); // should be in code units now //
                    double corr_mfp = R_grain_code / ((9./4.) * mean_free_path);
                    if(corr_mfp > 1) {tstop_inv /= corr_mfp;}
                }
#endif
#if defined(GRAIN_LORENTZFORCE) && defined(GRAIN_EPSTEIN_STOKES) /* also have charged grains, so we will calculate Coulomb forces as well */
                if(grain_subtype == 1)
                {
                    double a_Coulomb = sqrt(2.*gamma_eff*gamma_eff*gamma_eff/(9.*M_PI));
                    double tstop_Coulomb_inv = 0.797885/sqrt(gamma_eff) * rho_gas * cs / (R_grain_code * rho_grain_code); // base normalization //
                    tstop_Coulomb_inv /= (1. + a_Coulomb *(vgas_mag/cs)*(vgas_mag/cs)*(vgas_mag/cs)) * sqrt(1.+x0*x0); // velocity dependence (force becomes weak when super-sonic)
                    tstop_Coulomb_inv *= (Z_grain/tau_draine_sutin) * (Z_grain/tau_draine_sutin) / 17.; // coulomb attraction terms, assuming ions have charge ~1, and Coulomb logarithm is 17
                    // don't need super-accuration gas ionization states, just need approximate estimate, which we can make based on temperature //
                    double T_Kelvin = (2.3*PROTONMASS_CGS) * (cs_cgs*cs_cgs) / (1.3807e-16 * gamma_eff), f_ion_to_use = 0; // temperature in K
                    if(T_Kelvin > 1000.) {f_ion_to_use = exp(-15000./T_Kelvin);} /* default to a simple approximate guess for ionization, without cooling active */
#ifdef COOLING  // in this case, have the ability to calculate more accurate ionization fraction
                    double u_tmp, ne_tmp = 1, nh0_tmp = 0, mu_tmp = 1, temp_tmp, nHeII_tmp, nhp_tmp, nHe0_tmp, nHepp_tmp;
                    u_tmp = T_Kelvin / (2.3 * U_TO_TEMP_UNITS); // needs to be in code units; 2.3 for mean molecular weight factor and gamma_eos factor //
                    temp_tmp = ThermalProperties(u_tmp, rho_gas, -1, &mu_tmp, &ne_tmp, &nh0_tmp, &nhp_tmp, &nHe0_tmp, &nHeII_tmp, &nHepp_tmp);
                    f_ion_to_use = DMIN(ne_tmp , 1.);
#endif
                    tstop_Coulomb_inv *= f_ion_to_use; // correct for ionization fraction
                    tstop_inv += tstop_Coulomb_inv; // add both forces
                }
#endif // LORENTZ + EPSTEIN/STOKES force (Coulomb computation)


                /* this external_forcing parameter includes additional grain-specific forces. note that -anything- which imparts an
                 identical acceleration onto gas and dust will cancel in the terms in t_stop, and just act like a 'normal' acceleration
                 on the dust. for this reason the gravitational acceleration doesn't need to enter our 'external_forcing' parameter */
                double external_forcing[3]={0}, eps=MIN_REAL_NUMBER; P[i].Grain_AccelTimeMin = DMAX(1./(eps+tstop_inv) , sqrt(Get_Particle_Size(i)*All.cf_atime/(eps+vgas_mag*tstop_inv)));
#ifdef GRAIN_LORENTZFORCE
                if(grain_subtype == 1)
                {
                    /* Lorentz force on a grain = Z*e/c * ([v_grain-v_gas] x B) :: this comes from E-field E0 = -v_gas x B,
                        with force per particle Fp = (1 - R) * (np*E0 + J_p/c x B) / np = (1-R)*(E0 + v_p x B);
                        we ignore the Hall effect setting R=0 (ignore current carried by the particles themselves in induction) */
                    double grain_mass = (4.*M_PI/3.) * R_grain_code*R_grain_code*R_grain_code * rho_grain_code; // code units
                    double lorentz_units = UNIT_B_IN_GAUSS; // code B to Gauss
                    lorentz_units *= (ELECTRONCHARGE_CGS/C_LIGHT_CGS) * UNIT_VEL_IN_CGS / UNIT_MASS_IN_CGS; // converts acceleration to cgs
                    lorentz_units /= UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS; // converts it to code-units acceleration

                    double bhat[3], bmag=0, efield[3]={0}, efield_coeff=0, dv[3]; /* define unit vectors and B for evolving the lorentz force */
                    for(k=0;k<3;k++) {bhat[k]=P[i].Gas_B[k]*All.cf_a2inv; bmag+=bhat[k]*bhat[k]; dv[k]=(P[i].Vel[k]-P[i].Gas_Velocity[k])/All.cf_atime;}
                    if(bmag>0) {bmag=sqrt(bmag); for(k=0;k<3;k++) {bhat[k]/=bmag;}} else {bmag=0;}
                    double grain_charge_cinv = Z_grain / grain_mass * lorentz_units;
#ifdef GRAIN_RDI_TESTPROBLEM
                    if(All.Grain_Charge_Parameter != 0) {grain_charge_cinv = -All.Grain_Charge_Parameter/All.Grain_Size_Max * pow(All.Grain_Size_Max/P[i].Grain_Size,2);} // set charge manually //
                    //if(fabs(grain_charge_cinv)>0) {grain_charge_cinv /= 1.e-3 + P[i].Gas_Density;} /* this is the 'photoelectric' scaling for isothermal gas; modify for your charge law */
#endif
                    /* now apply the boris integrator */
                    double lorentz_coeff = (0.5*dt) * bmag * grain_charge_cinv; // dimensionless half-timestep term for boris integrator //
                    double v_m[3]={0}, v_t[3]={0}, v_p[3]={0}, vcrosst[3]={0};
                    for(k=0;k<3;k++) {v_m[k] = dv[k] + 0.5*efield_coeff*efield[k];} // half-step from E-field
                    /* cross-product for rotation */
                    vcrosst[0] = v_m[1]*bhat[2] - v_m[2]*bhat[1]; vcrosst[1] = v_m[2]*bhat[0] - v_m[0]*bhat[2]; vcrosst[2] = v_m[0]*bhat[1] - v_m[1]*bhat[0];
                    double tL=1./(eps+0.5*bmag*fabs(grain_charge_cinv)), vgasXB_mag=0; for(k=0;k<3;k++) {vgasXB_mag+=vcrosst[k]*vcrosst[k];}
                    P[i].Grain_AccelTimeMin = DMIN(P[i].Grain_AccelTimeMin, DMAX(tL , sqrt(Get_Particle_Size(i)*All.cf_atime/(eps+sqrt(vgasXB_mag)/tL))));
                    for(k=0;k<3;k++) {v_t[k] = v_m[k] + lorentz_coeff * vcrosst[k];} // first half-rotation
                    vcrosst[0] = v_t[1]*bhat[2] - v_t[2]*bhat[1]; vcrosst[1] = v_t[2]*bhat[0] - v_t[0]*bhat[2]; vcrosst[2] = v_t[0]*bhat[1] - v_t[1]*bhat[0];
                    for(k=0;k<3;k++) {v_p[k] = v_m[k] + (2.*lorentz_coeff/(1.+lorentz_coeff*lorentz_coeff)) * vcrosst[k];} // second half-rotation
                    for(k=0;k<3;k++) {v_p[k] += 0.5*efield_coeff*efield[k];} // half-step from E-field
                    /* calculate effective acceleration from discrete step in velocity */
                    for(k=0;k<3;k++) {external_forcing[k] += (v_p[k] - dv[k]) / dt;} // boris integrator
                }
#endif

                double delta_egy=0, delta_mom[3]={0}, xf=0, dt_tinv=dt*tstop_inv, C1=-(1+sqrt(1+x0*x0))/x0;
                if(dt_tinv < 100.) {double C2 = C1*exp(dt_tinv); xf=-2*C2/(C2*C2-1);}
                double slow_fac = 1 - xf/x0; /* note that, with an external (gravitational) acceleration, we can still solve this equation for the relevant update */
                for(k=0; k<3; k++)
                {
                    /* measure the imparted energy and momentum as if there were no external acceleration */
                    double v_init = P[i].Vel[k] / All.cf_atime, v_gas_i = P[i].Gas_Velocity[k] / All.cf_atime; /* physical units */
                    double vel_new = v_init + slow_fac * (v_gas_i - v_init);
                    /* now calculate the updated velocity accounting for any external, non-standard accelerations */
                    double vdrift = 0, dv[3];
                    if(tstop_inv > 0) {vdrift = external_forcing[k] / (tstop_inv * sqrt(1+x0*x0));}
                    dv[k] = slow_fac * (v_gas_i - v_init + vdrift);
                    if(isnan(vdrift)||isnan(slow_fac)) {dv[k] = 0;}

                    vel_new = v_init + dv[k];
                    delta_mom[k] = P[i].Mass * (vel_new - v_init);
                    delta_egy += 0.5*P[i].Mass * (vel_new*vel_new - v_init*v_init);
#ifdef GRAIN_BACKREACTION
                    P[i].Grain_DeltaMomentum[k] = delta_mom[k] * All.cf_atime; /* converted back to code units here */
#endif
                    /* note, we can directly apply this by taking P[i].Vel[k] += dv[k]; but this is not as accurate as our
                     normal leapfrog integration scheme. we can also account for the -gas- acceleration, by including it like vdrift;
                     for a constant t_stop, the gas acceleration term appears as P[i].Vel[l] += Gas_Accel[k] * dt + slow_fac * (Gas-Accel[k] / tstop_inv) */
                    P[i].GravAccel[k] += dv[k] / (dt * All.cf_a2inv); /* we solve the equations with an external acceleration already (external_forcing above): therefore add to forces like gravity that are acting on the gas and dust in the same manner (in terms of acceleration). put into code units */
                }
            } // closes check for gas density, dt, vmag > 0, subtype valid


#ifdef PIC_MHD
#ifndef PIC_SPEEDOFLIGHT_REDUCTION
#define PIC_SPEEDOFLIGHT_REDUCTION (1)
#endif            
            if((grain_subtype >= 3) && (dt > 0) && (P[i].Gas_Density>0) && (vgas_mag > 0)) /* only bother with particles moving wrt gas with finite gas density and timestep */
            {
                double reduced_C = PIC_SPEEDOFLIGHT_REDUCTION * C_LIGHT_CODE; /* effective speed of light for this part of the code */
                double charge_to_mass_ratio_dimensionless = All.PIC_Charge_to_Mass_Ratio; /* dimensionless q/m in units of e/mp */
                //if(grain_subtype==4) {charge_to_mass_ratio_dimensionless = -1836.15; /* electrons */

                double lorentz_units = UNIT_B_IN_GAUSS * UNIT_VEL_IN_CGS * (ELECTRONCHARGE_CGS/(PROTONMASS_CGS*C_LIGHT_CGS)) / (UNIT_VEL_IN_CGS/UNIT_TIME_IN_CGS); // code velocity to CGS and B to Gauss, times base units e/(mp*c), then convert 'back' to code-units acceleration
#ifdef PIC_MHD_NEW_RSOL_METHOD
                lorentz_units *= PIC_SPEEDOFLIGHT_REDUCTION; // the rsol enters by slowing down the forces here, acts as a unit shift for time
#endif
                double efield[3], bhat[3]={0}, bmag=0, v_g[3]; /* define unit vectors and B for evolving the lorentz force */
                for(k=0;k<3;k++) {bhat[k]=P[i].Gas_B[k]*All.cf_a2inv; bmag+=bhat[k]*bhat[k]; v_g[k]=P[i].Gas_Velocity[k]/(All.cf_atime*reduced_C);} /* get magnitude and unit vector for B, and vector beta [-true- beta here] */
                if(bmag>0) {bmag=sqrt(bmag); for(k=0;k<3;k++) {bhat[k]/=bmag;}} else {bmag=0;} /* take it correctly assuming its non-zero */
                double efield_coeff = (0.5*dt) * charge_to_mass_ratio_dimensionless * bmag * lorentz_units; // dimensionless half-timestep term for boris integrator //
                efield[0] = -v_g[1]*bhat[2] + v_g[2]*bhat[1]; efield[1] = -v_g[2]*bhat[0] + v_g[0]*bhat[2]; efield[2] = -v_g[0]*bhat[1] + v_g[1]*bhat[0]; /* efield term, but with magnitude of B factored out for units above */
                double v_0[3],v0[3],vf[3],v2=0; for(k=0;k<3;k++) {v0[k]=P[i].Vel[k]/All.cf_atime; v2+=v0[k]*v0[k];} // magnitude of velocity [this is reduced from c]
                if(v2 >= reduced_C*reduced_C) {PRINT_WARNING("VELOCITY HAS EXCEEDED THE SPEED OF LIGHT. BAD.");} // check against reduced c
                double gamma_0=1/sqrt(1-v2/(reduced_C*reduced_C)); for(k=0;k<3;k++) {v_0[k]=v0[k]*gamma_0/reduced_C;} // calculate true gamma, convert to the momentum term ~gamma*beta (this times mc is true scalar momentum)

                /* now apply the boris integrator */
                double v_m[3]={0}, v_t[3]={0}, v_p[3]={0}, vcrosst[3]={0}, lorentz_coeff=efield_coeff;
                for(k=0;k<3;k++) {v_m[k] = v_0[k] + efield_coeff*efield[k];} // first half-step from E-field
                lorentz_coeff /= sqrt(1+v_m[0]*v_m[0]+v_m[1]*v_m[1]+v_m[2]*v_m[2]); // lorentz factor at this mid-point jump (recall v_m is a 'beta~v/c' here) is needed to correct the factor for the B-field term
                vcrosst[0] = v_m[1]*bhat[2] - v_m[2]*bhat[1]; vcrosst[1] = v_m[2]*bhat[0] - v_m[0]*bhat[2]; vcrosst[2] = v_m[0]*bhat[1] - v_m[1]*bhat[0]; // cross-product for rotation
                for(k=0;k<3;k++) {v_t[k] = v_m[k] + lorentz_coeff * vcrosst[k];} // first half-rotation
                vcrosst[0] = v_t[1]*bhat[2] - v_t[2]*bhat[1]; vcrosst[1] = v_t[2]*bhat[0] - v_t[0]*bhat[2]; vcrosst[2] = v_t[0]*bhat[1] - v_t[1]*bhat[0];
                for(k=0;k<3;k++) {v_p[k] = v_m[k] + (2.*lorentz_coeff/(1.+lorentz_coeff*lorentz_coeff)) * vcrosst[k];} // second half-rotation
                for(k=0;k<3;k++) {v_p[k] += efield_coeff*efield[k];} // second half-step from E-field. v_p now contains the final scalar momentum in dimensionless units, i.e. gamma*beta. so this divided by gamma gives final beta
                double vp2=v_p[0]*v_p[0]+v_p[1]*v_p[1]+v_p[2]*v_p[2], gamma_f=sqrt(1+vp2); for(k=0;k<3;k++) {vf[k]=reduced_C*v_p[k]/gamma_f;} // convert back to a velocity 'vf' which is always <= reduced_C - this is now the 'effective' velocity with which CRs will propagate

                for(k=0;k<3;k++)
                {
#ifdef GRAIN_BACKREACTION
                    double delta_momentum = P[i].Mass * (vf[k]*gamma_f - v0[k]*gamma_0) * All.cf_atime; /* account for lorentz factor in calculating the discrete momentum change here [put into code units] */
#ifdef PIC_MHD_NEW_RSOL_METHOD
                    delta_momentum /= PIC_SPEEDOFLIGHT_REDUCTION*PIC_SPEEDOFLIGHT_REDUCTION; // the real force back on the gas is the difference in the conserved quantity, Delta[(c/tilde[c]*gamma*beta_vector*mc], which requires multiplying the above by (c/RSOL)^2
#endif
                    P[i].Grain_DeltaMomentum[k] += delta_momentum; // save to couple back to gas in loop below
#endif
                    P[i].GravAccel[k] += (vf[k]-v0[k]) / (dt * All.cf_a2inv); /* update acceleration with the kick from the full boris push above [put into code units] */
                }
            } // closes check for gas density, dt, vmag > 0, subtype valid
#endif

        } // closes check for particle type, id
    } // closes main particle loop (loop over active particles)
#if defined(GRAIN_BACKREACTION)
    grain_backrx(); /* call parent routine to assign the back-reaction force among neighbors */
#endif
    PRINT_STATUS(" ..particulate/grain/PIC force evaluation done.");
    CPU_Step[CPU_DRAGFORCE] += measure_time();
}




/* this is a template for fully-automated parallel (hybrid MPI+OpenMP/Pthreads) neighbor communication
 written in a completely modular fashion. this works as long as what you are trying to do isn't too
 complicated (from a communication point-of-view). You specify a few key variables, and then
 define the variables that need to be passed, and write the actual sub-routine that does the actual
 'work' between neighbors, but all of the parallelization, looping, communication blocks,
 etc, are all handled for you. */

#define CORE_FUNCTION_NAME grain_backrx_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(((1 << P[i].Type) & (GRAIN_PTYPES))&&(P[i].TimeBin>=0)&&(P[i].Mass>0)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */


/* this structure defines the variables that need to be sent -from- the 'searching' element */
struct INPUT_STRUCT_NAME
{
    int NodeList[NODELISTLENGTH]; MyDouble Pos[3], Hsml; /* these must always be defined */
#if defined(GRAIN_BACKREACTION)
    double Grain_DeltaMomentum[3], Gas_Density, Grain_AccelTimeMin;
#endif
}
*DATAIN_NAME, *DATAGET_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{   /* "i" is the particle from which data will be assigned, to structure "in" */
    int k; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];} /* good example - always needed */
    in->Hsml = PPP[i].Hsml; /* also always needed for search (can change the radius "PPP[i].Hsml" but in->Hsml must be defined */
#if defined(GRAIN_BACKREACTION)
    for(k=0;k<3;k++) {in->Grain_DeltaMomentum[k]=P[i].Grain_DeltaMomentum[k];}
    in->Gas_Density = P[i].Gas_Density;
    in->Grain_AccelTimeMin = P[i].Grain_AccelTimeMin;
#endif
}

/* this structure defines the variables that need to be sent -back to- the 'searching' element */
struct OUTPUT_STRUCT_NAME
{ /* define variables below as e.g. "double X;" */
}
*DATARESULT_NAME, *DATAOUT_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{  /* "i" is the particle to which data from structure "out" will be assigned. mode=0 for local communication,
    =1 for data sent back from other processors. you must account for this. */
    /* example: ASSIGN_ADD(P[i].X,out->X,mode); which is short for: if(mode==0) {P[i].X=out->X;} else {P[i].X+=out->X;} */
}


/* this subroutine does the actual neighbor-element calculations (this is the 'core' of the loop, essentially) */
/*!   -- this subroutine writes to shared memory [updating the neighbor values]: need to protect these writes for openmp below. modified values for the minimum timestep are read, so both read and write need to be protected. */
int grain_backrx_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, n; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* define variables and zero memory and import data for local target*/
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* imports the data to the correct place and names */

    if(local.Hsml <= 0) {return 0;} /* don't bother doing a loop if this isnt going to do anything */
    int kernel_shared_BITFLAG = 1; /* grains 'see' gas in this loop */

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, kernel_shared_BITFLAG);
            if(numngb_inbox < 0) {return -2;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if((P[j].Mass <= 0)||(P[j].Hsml <= 0)) {continue;} /* make sure neighbor is valid */
                int k; double dp[3]; for(k=0;k<3;k++) {dp[k]=local.Pos[k]-P[j].Pos[k];} /* position offset */
                NEAREST_XYZ(dp[0],dp[1],dp[2],1); double r2=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2]; /* box-wrap appropriately and calculate distance */
#ifdef BOX_BND_PARTICLES
                if(P[j].ID > 0) {r2 = -1;} /* ignore frozen boundary particles */
#endif
                if((r2>0)&&(r2<local.Hsml*local.Hsml)) /* only keep elements inside search radius */
                {
                    double hinv,hinv3,hinv4,wk_i=0,dwk_i=0,r=sqrt(r2); kernel_hinv(local.Hsml,&hinv,&hinv3,&hinv4);
                    kernel_main(r*hinv, hinv3, hinv4, &wk_i, &dwk_i, 0); /* kernel quantities that may be needed */
#if defined(GRAIN_BACKREACTION)
                    double wt = -wk_i / local.Gas_Density, dv2=0; /* degy=wt*delta_egy; */
                    for(k=0;k<3;k++) {
                        double dv = wt*local.Grain_DeltaMomentum[k]; // momentum to be sent to this neighbor element
                        dv2+=dv*dv; // save squared sum
                        #pragma omp atomic
                        P[j].Vel[k] += dv; // add the velocity (checking for thread safety in doing so!)
                        #pragma omp atomic
                        SphP[j].VelPred[k] += dv; // add the velocity (checking for thread safety in doing so!)
                    }
                    
                    double taccel_min_prev = 0, taccel_min_new = 0;
                    #pragma omp atomic read
                    taccel_min_prev = P[j].Grain_AccelTimeMin; // this can be modified below so needs to be done in a thread-safe manner here //
                    taccel_min_new = DMIN(DMIN(2.*All.ErrTolIntAccuracy*Get_Particle_Size(j)*All.cf_atime*All.cf_atime/sqrt(dv2+MIN_REAL_NUMBER) , 4.*local.Grain_AccelTimeMin), taccel_min_prev);
                    if(taccel_min_new < taccel_min_prev)
                    {
                        #pragma omp atomic write
                        P[j].Grain_AccelTimeMin = taccel_min_new;
                    }
                    /* for(k=0;k<3;k++) {degy-=P[j].Mass*dv*(VelPred_j[k]+0.5*dv);} SphP[j].InternalEnergy += degy; */ // ignoring these terms -- if re-add them be sure to do so thread-safely //
#endif
                }
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


void grain_backrx(void)
{
    PRINT_STATUS(" ..assigning grain back-reaction to gas");
     //grain_backrx_initial_operations_preloop(); /* do initial pre-processing operations as needed before main loop [nothing needed here] */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    //grain_backrx_final_operations_and_cleanup(); /* do final operations on results [nothing needed here] */
    CPU_Step[CPU_DRAGFORCE] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */







#ifdef GRAIN_COLLISIONS /* these routines [combined with the DM_SIDM modules] allow for arbitrary grain-grain scattering-sticking-destruction-type interactions */

/*! This function returns the cross-section per unit mass of the grains represented by a single 'super-particle'. Modify this appropriately for your assumed grain
    physics. Here we are simply assuming hard-sphere scattering. */
double return_grain_cross_section_per_unit_mass(int i)
{   /* All.DM_InteractionCrossSection just serves here as a renormalization parameter, optionally */
    return All.DM_InteractionCrossSection * 0.75 / (P[i].Grain_Size * All.Grain_Internal_Density); // this will be in cgs units, in our standard convention //
}

/*! This routine determines the probability of a grain-grain interaction within the kernel. Note that the function 'g_geo' does all the accounting for the kernel
    structure and finite size -- the rest of this can be initialized as a normal scattering rate (~density*cross_section_per_unit_mass*velocity)
    where here 'All.DM_InteractionCrossSection' is the cross-section read in from the params file, and other params like DM_InteractionVelocityScale
    allow the user to control the collision velocity dependence. This function should be appropriately modified to the actual grain physics being represented.
    Here, the default assumption is simple hard-sphere scattering with a constant cross section per unit grain mass, set by the grain size */
double prob_of_grain_interaction(double cx_per_unitmass, double mass, double r, double h_si, double dV[3], double dt, int j_ngb)
{
    double dVmag = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]) / All.cf_atime; // velocity in physical
    double rho_eff = 0.5*(mass + P[j_ngb].Mass) / (h_si*h_si*h_si) * All.cf_a3inv; // density in physical
    double cx_eff = g_geo(r/h_si) * (mass*cx_per_unitmass + P[j_ngb].Mass*return_grain_cross_section_per_unit_mass(j_ngb)) / (mass + P[j_ngb].Mass); // mass-weighted effective cross section (physical) scaled to cgs
    double units = UNIT_SURFDEN_IN_CGS; // needed to convert everything to cgs
    if(All.DM_InteractionVelocityScale>0) {double x=dVmag/All.DM_InteractionVelocityScale; cx_eff/=1+x*x*x*x;} // take velocity dependence
    return rho_eff * cx_eff * dVmag * dt * units; // dimensionless probability
}

/*! This routine sets the kicks for each grain after it has been decided that they will interact. By default at present this results only in velocity 'kicks',
    but one can modify this function to allow other types of interactions. By default, it will assume elastic collisions,
    and an algorithm that conserves energy and momentum but picks a random direction so it does not conserves angular momentum. */
void calculate_interact_kick(double dV[3], double kick[3], double m)
{
    double dVmag = (1-All.DM_DissipationFactor)*sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]);
    if(dVmag<0) {dVmag=0;}
    if(All.DM_KickPerCollision>0) {double v0=All.DM_KickPerCollision; dVmag=sqrt(dVmag*dVmag+v0*v0);}
    double cos_theta = 2.0*gsl_rng_uniform(random_generator)-1.0, sin_theta = sqrt(1.-cos_theta*cos_theta), phi = gsl_rng_uniform(random_generator)*2.0*M_PI;
    kick[0] = 0.5*(dV[0] + dVmag*sin_theta*cos(phi));
    kick[1] = 0.5*(dV[1] + dVmag*sin_theta*sin(phi));
    kick[2] = 0.5*(dV[2] + dVmag*cos_theta);
}

#endif






#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)

#define CORE_FUNCTION_NAME interpolate_fluxes_opacities_gasgrains_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(((1 << P[i].Type) & (GRAIN_PTYPES+1))&&(P[i].TimeBin>=0)&&(P[i].Mass>0)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/* this structure defines the variables that need to be sent -from- the 'searching' element */
struct INPUT_STRUCT_NAME {
    int NodeList[NODELISTLENGTH], Type; MyDouble Mass, Hsml, Pos[3], Vel[3], Grain_Size, Grain_Abs_Coeff[N_RT_FREQ_BINS]; /* these must always be defined */
} *DATAIN_NAME, *DATAGET_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration) {   /* "i" is the particle from which data will be assigned, to structure "in" */
    in->Type=P[i].Type; in->Mass=P[i].Mass; in->Hsml=PPP[i].Hsml; int k; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k]; in->Vel[k]=P[i].Vel[k];}
    if((1 << P[i].Type) & (GRAIN_PTYPES))
    {
        in->Grain_Size=P[i].Grain_Size; int k_freq;
        double R_grain_code=P[i].Grain_Size/UNIT_LENGTH_IN_CGS, rho_grain_code=All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS, rho_gas_code=P[i].Gas_Density*All.cf_a3inv; /* internal grain density in code units */
        for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
        {
            double Q_abs_eff = return_grain_extinction_efficiency_Q(i, k_freq); /* need this to calculate the absorption efficiency in each band */
            in->Grain_Abs_Coeff[k_freq] = Q_abs_eff * 3. / (4. * C_LIGHT_CODE_REDUCED * rho_grain_code * R_grain_code * rho_gas_code);
        }
    }
}

/* this structure defines the variables that need to be sent -back to- the 'searching' element */
struct OUTPUT_STRUCT_NAME { /* define variables below as e.g. "double X;" */
    MyDouble Interpolated_Radiation_Acceleration[3]; /* flux values to return to grains */
    MyDouble Interpolated_Opacity[N_RT_FREQ_BINS]; /* opacity values interpolated to gas positions */
} *DATARESULT_NAME, *DATAOUT_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration) {  /* "i" is the particle to which data from structure "out" will be assigned. mode=0 for local communication, =1 for data sent back from other processors. you must account for this. */
    int k,k_freq;
    if(P[i].Type==0) {for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {ASSIGN_ADD(SphP[i].Interpolated_Opacity[k_freq],out->Interpolated_Opacity[k_freq],mode);}}
    if((1 << P[i].Type) & (GRAIN_PTYPES)) {for(k=0;k<3;k++) {P[i].GravAccel[k] += out->Interpolated_Radiation_Acceleration[k]/All.cf_a2inv;}} /* this simply adds to the 'gravitational' acceleration for kicks */
}

/* this subroutine does the actual neighbor-element calculations (this is the 'core' of the loop, essentially) */
int interpolate_fluxes_opacities_gasgrains_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, n; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* define variables and zero memory and import data for local target*/
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* imports the data to the correct place and names */

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            if(local.Type == 0)
            {
                numngb_inbox = ngb_treefind_pairs_threads_targeted(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, GRAIN_PTYPES); /* gas searches for grains which can see -it- */
            } else {
                numngb_inbox = ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist); /* grains search for gas -they- can see */
            }
            if(numngb_inbox < 0) {return -2;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if((P[j].Mass <= 0)||(PPP[j].Hsml <= 0)) {continue;} /* make sure neighbor is valid */
                int k,k_freq; double dp[3],h_to_use; for(k=0;k<3;k++) {dp[k]=local.Pos[k]-P[j].Pos[k];} /* position offset */
                NEAREST_XYZ(dp[0],dp[1],dp[2],1); double r2=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2]; /* box-wrap appropriately and calculate distance */
                if(local.Type == 0) {h_to_use = PPP[j].Hsml;} else {h_to_use = local.Hsml;}
                if((r2>0)&&(r2<h_to_use*h_to_use)) /* only keep elements inside search radius */
                {
                    double wt=0,hinv,hinv3,hinv4,wk_i=0,dwk_i=0,r=sqrt(r2); kernel_hinv(h_to_use,&hinv,&hinv3,&hinv4);
                    kernel_main(r*hinv, hinv3, hinv4, &wk_i, &dwk_i, 0); /* kernel quantities that may be needed */

                    if(local.Type==0) /* sitting on a -gas- element, want to interpolate opacity to it */
                    {
                        wt = P[j].Mass * (wk_i / P[j].Gas_Density); /* dimensionless weight of this gas element as 'seen' by the grain: = (grain_part_mass/gas_part_mass) * (gas_part_mass * Wk / gas_density [=sum gas_part_mass * Wk]) */
                        double R_grain_code=P[j].Grain_Size/UNIT_LENGTH_IN_CGS, rho_grain_code=All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS; /* internal grain density in code units */
                        for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
                        {
                            double Q_abs_eff = return_grain_extinction_efficiency_Q(j, k_freq); /* need this to calculate the absorption efficiency in each band */
                            out.Interpolated_Opacity[k_freq] += wt * Q_abs_eff * 3. / (4. * rho_grain_code * R_grain_code);
                        }
                    } else { /* sitting on a -grain- element, want to interpolate flux to it and calculate radiation pressure force */
                        wt = SphP[j].Density*All.cf_a3inv * wk_i; /* weight of element to 'i, with appropriate coefficient from above */
                        double radacc[3]={0},vel_i[3]={0},dtEgamma_work_done=0; for(k=0;k<3;k++) {vel_i[k]=RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS*local.Vel[k]/All.cf_atime;} /* velocity of interest here is the grain velocity (radiation in lab frame) */
                        for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
                        {
                            double f_kappa_abs=0.5,vdot_h[3]={0},flux_i[3]={0},flux_mag=0,erad_i=0,flux_corr=1;
                            f_kappa_abs = 0.5; // rt_absorb_frac_albedo(i,k_freq); -- this is set to 1/2 anyways but would require extra passing, ignore for now //
#if defined(RT_EVOLVE_FLUX) || (defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX) && defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY))
                            erad_i = SphP[j].Rad_E_gamma_Pred[k_freq];
                            eddington_tensor_dot_vector(SphP[j].ET[k_freq],vel_i,vdot_h); for(k=0;k<3;k++) {vdot_h[k] = erad_i * (vel_i[k] + vdot_h[k]);} // calculate volume integral of scattering coefficient t_inv * (gas_vel . [e_rad*I + P_rad_tensor]), which gives an additional time-derivative term. this is the P term //
                            for(k=0;k<3;k++) {flux_i[k]=SphP[j].Rad_Flux_Pred[k_freq][k]; flux_mag+=flux_i[k]*flux_i[k];}
                            if(flux_mag>0 && isfinite(flux_mag)) {flux_mag=sqrt(flux_mag);} else {flux_mag=MIN_REAL_NUMBER; flux_i[0]=flux_i[1]=0; flux_i[2]=flux_mag;}
#elif (defined(RT_OTVET) || defined(RT_FLUXLIMITEDDIFFUSION))
                            erad_i = SphP[j].Rad_E_gamma_Pred[k_freq];
                            for(k=0;k<3;k++) {flux_i[k] = -SphP[j].Gradients.Rad_E_gamma_ET[k_freq][k]; flux_mag+=flux_i[k]*flux_i[k];}
                            if(flux_mag>0) {for(k=0;k<3;k++) {flux_i[k]/=sqrt(flux_mag);}} else {flux_i[0]=0;flux_i[1]=0;flux_i[2]=1;}
                            flux_mag = erad_i*C_LIGHT_CODE_REDUCED; for(k=0;k<3;k++) {flux_i[k]*=flux_mag;}
#endif
                            if(!isfinite(flux_mag) || flux_mag<=MIN_REAL_NUMBER) {flux_mag=MIN_REAL_NUMBER; flux_i[0]=flux_i[1]=0; flux_i[2]=flux_mag;}
                            double flux_thin = erad_i * C_LIGHT_CODE_REDUCED; if(!isfinite(flux_thin) || flux_thin<=0) {flux_thin=0;}
                            flux_corr = DMIN(1., 100.*flux_thin/flux_mag);
                            for(k=0;k<3;k++)
                            {
                                radacc[k] += local.Grain_Abs_Coeff[k_freq] * wt * (flux_corr*flux_i[k] - vdot_h[k]); // note these 'vdoth' terms shouldn't be included in FLD, since its really assuming the entire right-hand-side of the flux equation reaches equilibrium with the pressure tensor, which gives the expression in rt_utilities
                                dtEgamma_work_done += (2.*f_kappa_abs-1.) * radacc[k] * vel_i[k] * local.Mass; // PdV work done by photons [absorbed ones are fully-destroyed, so their loss of energy and momentum is already accounted for by their deletion in this limit -- note that we have to be careful about the RSOL factors here! //
                            }
                        }
                        for(k=0;k<3;k++) {out.Interpolated_Radiation_Acceleration[k] += radacc[k];} /* prepare to send back to grain */
                    }
                }
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}

void interpolate_fluxes_opacities_gasgrains(void)
{
    PRINT_STATUS(" ..assigning opacities to gas from the grain distribution, and interpolating radiation fields to grains");
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    CPU_Step[CPU_DRAGFORCE] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */



double return_grain_extinction_efficiency_Q(int i, int k_freq)
{
    double Q = 1; /* default to geometric opacity */
#if defined(GRAIN_RDI_TESTPROBLEM)
    Q *= All.Grain_Q_at_MaxGrainSize; // this needs to be set by-hand, Q for the maximum sized grains. irrelevant for the scale-free problem (degenerate with flux), but important here */
#if !defined(GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE)
    Q *= P[i].Grain_Size / All.Grain_Size_Max;
#endif
#else
    /* INSERT PHYSICS HERE -- this is where you want to specify the optical properties of grains relative to the frequency bins being evolved. could code up something for -ALL- the bins we do, but that's a lot, so we'll do these as-needed, for runs with different frequencies */
    if(ThisTask==0) {PRINT_WARNING("Code does not have entered grain absorption efficiency/optical properties for your specific wavelength being evolved. Please enter that information in the routine 'return_grain_extinction_efficiency_Q'. For now will assume geometric absorption (Q=1). \n");}
#endif
    return Q;
}



#endif //defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)




#endif
