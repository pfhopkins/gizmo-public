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
 winds, and SNe remnants, as well as terrestrial turbulance and
 particulate-laden turbulence. Anywhere where particles coupled to gas
 via coulomb, aerodynamic, or lorentz forces are interesting.
 
 This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 
 */



#ifdef GRAIN_FLUID

#ifdef GRAIN_BACKREACTION
void grain_backrx(void);
#endif

/* function to apply the drag on the grains from surrounding gas properties */
void apply_grain_dragforce(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    int i, k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef GRAIN_BACKREACTION
        for(k=0;k<3;k++) {P[i].Grain_DeltaMomentum[k]=0;}
#endif
        if((P[i].Type != 0)&&(P[i].Type != 4))
        {
#ifdef BOX_BND_PARTICLES
            if(P[i].ID > 0)
#endif
                if(P[i].Gas_Density > 0)
                {
                    double dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
                    if(dt > 0)
                    {
                        double cs = sqrt( GAMMA * GAMMA_MINUS1 * P[i].Gas_InternalEnergy);
                        double R_grain_cgs = P[i].Grain_Size;
                        double R_grain_code = R_grain_cgs / (All.UnitLength_in_cm / All.HubbleParam);
                        double rho_gas = P[i].Gas_Density * All.cf_a3inv;
                        double rho_grain_physical = All.Grain_Internal_Density; // cgs units //
                        double rho_grain_code = rho_grain_physical / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam); // code units //
                        double vgas_mag = 0.0;
                        for(k=0;k<3;k++) {vgas_mag+=(P[i].Gas_Velocity[k]-P[i].Vel[k])*(P[i].Gas_Velocity[k]-P[i].Vel[k]);}
                        
                        if(vgas_mag > 0)
                        {
                            vgas_mag = sqrt(vgas_mag) / All.cf_atime;
                            double x0 = 0.469993*sqrt(GAMMA) * vgas_mag/cs; // (3/8)*sqrt[pi/2]*|vgas-vgrain|/cs //
                            double tstop_inv = 1.59577/sqrt(GAMMA) * rho_gas * cs / (R_grain_code * rho_grain_code); // 2*sqrt[2/pi] * 1/tstop //
#ifdef GRAIN_LORENTZFORCE
                            /* calculate the grain charge following Draine & Sutin */
                            double cs_cgs = cs * All.UnitVelocity_in_cm_per_s;
                            double tau_draine_sutin = R_grain_cgs * (2.3*PROTONMASS) * (cs_cgs*cs_cgs) / (GAMMA * ELECTRONCHARGE*ELECTRONCHARGE);
                            double Z_grain = -DMAX( 1./(1. + sqrt(1.0e-3/tau_draine_sutin)) , 2.5*tau_draine_sutin );
                            if(isnan(Z_grain)||(Z_grain>=0)) {Z_grain=0;}
#endif
#ifdef GRAIN_EPSTEIN_STOKES
                            double mu = 2.3 * PROTONMASS;
                            double temperature = mu * (P[i].Gas_InternalEnergy*All.UnitEnergy_in_cgs*All.HubbleParam/All.UnitMass_in_g) / BOLTZMANN;
                            double cross_section = GRAIN_EPSTEIN_STOKES * 2.0e-15 * (1. + 70./temperature);
                            cross_section /= (All.UnitLength_in_cm * All.UnitLength_in_cm / (All.HubbleParam*All.HubbleParam));
                            double n_mol = rho_gas / (mu * All.HubbleParam/All.UnitMass_in_g);
                            double mean_free_path = 1 / (n_mol * cross_section); // should be in code units now //
                            double corr_mfp = R_grain_code / ((9./4.) * mean_free_path);
                            if(corr_mfp > 1) {tstop_inv /= corr_mfp;}
#ifdef GRAIN_LORENTZFORCE
                            /* also have charged grains, so we will calculate Coulomb forces as well */
                            double a_Coulomb = sqrt(2.*GAMMA*GAMMA*GAMMA/(9.*M_PI));
                            double tstop_Coulomb_inv = 0.797885/sqrt(GAMMA) * rho_gas * cs / (R_grain_code * rho_grain_code); // base normalization //
                            tstop_Coulomb_inv /= (1. + a_Coulomb *(vgas_mag/cs)*(vgas_mag/cs)*(vgas_mag/cs)) * sqrt(1.+x0*x0); // velocity dependence (force becomes weak when super-sonic)
                            tstop_Coulomb_inv *= (Z_grain/tau_draine_sutin) * (Z_grain/tau_draine_sutin) / 17.; // coulomb attraction terms, assuming ions have charge ~1, and Coulomb logarithm is 17
                            // don't need super-accuration gas ionization states, just need approximate estimate, which we can make based on temperature //
                            double T_Kelvin = (2.3*PROTONMASS) * (cs_cgs*cs_cgs) / (1.3807e-16 * GAMMA), f_ion_to_use = 0; // temperature in K
#ifdef COOLING  // in this case, have the ability to calculate more accurate ionization fraction
                            {
                                double u_tmp, ne_tmp = 1, nh0_tmp = 0, mu_tmp = 1, temp_tmp, nHeII_tmp, nhp_tmp, nHe0_tmp, nHepp_tmp;
                                u_tmp = 1.3807e-16 * T_Kelvin / (2.3*PROTONMASS) * (All.UnitMass_in_g/All.UnitEnergy_in_cgs); // needs to be in code units
                                temp_tmp = ThermalProperties(u_tmp, rho_gas, -1, &mu_tmp, &ne_tmp, &nh0_tmp, &nhp_tmp, &nHe0_tmp, &nHeII_tmp, &nHepp_tmp);
                                f_ion_to_use = DMIN(ne_tmp , 1.);
                            }
#else           // without cooling active, use a simple approximate guess for ionization
                            if(T_Kelvin > 1000.) {f_ion_to_use = exp(-15000./T_Kelvin);}
#endif
                            tstop_Coulomb_inv *= f_ion_to_use; // correct for ionization fraction
                            tstop_inv += tstop_Coulomb_inv; // add both forces
#endif // LORENTZ force (Coulomb computation)
#endif
                            double C1 = (-1-sqrt(1+x0*x0)) / x0;
                            double xf = 0.0;
                            double dt_tinv = dt * tstop_inv;
                            if(dt_tinv < 100.)
                            {
                                double C2 = C1 * exp( dt_tinv );
                                xf = -2 * C2 / (C2*C2 -1);
                            }
                            double slow_fac = 1 - xf / x0;
                            // note that, with an external (gravitational) acceleration, we can still solve this equation for the relevant update //
                            
                            double dv[3]={0}, external_forcing[3]={0};
                            for(k=0;k<3;k++) {external_forcing[k] = 0;}
                            /* this external_forcing parameter includes additional grain-specific forces. note that -anything- which imparts an
                             identical acceleration onto gas and dust will cancel in the terms in t_stop, and just act like a 'normal' acceleration
                             on the dust. for this reason the gravitational acceleration doesn't need to enter our 'external_forcing' parameter */
#ifdef GRAIN_LORENTZFORCE
                            /* Lorentz force on a grain = Z*e/c * ([v_grain-v_gas] x B) */
                            double v_cross_B[3];
                            v_cross_B[0] = (P[i].Vel[1]-P[i].Gas_Velocity[1])*P[i].Gas_B[2] - (P[i].Vel[2]-P[i].Gas_Velocity[2])*P[i].Gas_B[1];
                            v_cross_B[1] = (P[i].Vel[2]-P[i].Gas_Velocity[2])*P[i].Gas_B[0] - (P[i].Vel[0]-P[i].Gas_Velocity[0])*P[i].Gas_B[2];
                            v_cross_B[2] = (P[i].Vel[0]-P[i].Gas_Velocity[0])*P[i].Gas_B[1] - (P[i].Vel[1]-P[i].Gas_Velocity[1])*P[i].Gas_B[0];
                            
                            double grain_mass = (4.*M_PI/3.) * R_grain_code*R_grain_code*R_grain_code * rho_grain_code; // code units
                            double lorentz_units = sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam); // code B to Gauss
                            lorentz_units *= (ELECTRONCHARGE/C) * All.UnitVelocity_in_cm_per_s / (All.UnitMass_in_g / All.HubbleParam); // converts acceleration to cgs
                            lorentz_units /= All.UnitVelocity_in_cm_per_s / (All.UnitTime_in_s / All.HubbleParam); // converts it to code-units acceleration
                            
                            /* define unit vectors and B for evolving the lorentz force */
                            double bhat[3]={0}, bmag=0, efield[3]={0}, efield_coeff=0;
                            for(k=0;k<3;k++) {bhat[k]=P[i].Gas_B[k]; bmag+=bhat[k]*bhat[k]; dv[k]=P[i].Vel[k]-P[i].Gas_Velocity[k];}
                            if(bmag>0) {bmag=sqrt(bmag); for(k=0;k<3;k++) {bhat[k]/=bmag;}} else {bmag=0;}
                            double grain_charge_cinv = Z_grain / grain_mass * lorentz_units;
#ifdef GRAIN_RDI_TESTPROBLEM
                            if(All.Grain_Charge_Parameter != 0) {grain_charge_cinv = -All.Grain_Charge_Parameter/All.Grain_Size_Max * pow(All.Grain_Size_Max/P[i].Grain_Size,2);} // set charge manually //
#endif
                            /* now apply the boris integrator */
                            double lorentz_coeff = (0.5*dt) * bmag * grain_charge_cinv; // dimensionless half-timestep term for boris integrator //
                            double v_m[3]={0}, v_t[3]={0}, v_p[3]={0}, vcrosst[3]={0};
                            for(k=0;k<3;k++) {v_m[k] = dv[k] + 0.5*efield_coeff*efield[k];} // half-step from E-field
                            /* cross-product for rotation */
                            vcrosst[0] = v_m[1]*bhat[2] - v_m[2]*bhat[1]; vcrosst[1] = v_m[2]*bhat[0] - v_m[0]*bhat[2]; vcrosst[2] = v_m[0]*bhat[1] - v_m[1]*bhat[0];
                            for(k=0;k<3;k++) {v_t[k] = v_m[k] + lorentz_coeff * vcrosst[k];} // first half-rotation
                            vcrosst[0] = v_t[1]*bhat[2] - v_t[2]*bhat[1]; vcrosst[1] = v_t[2]*bhat[0] - v_t[0]*bhat[2]; vcrosst[2] = v_t[0]*bhat[1] - v_t[1]*bhat[0];
                            for(k=0;k<3;k++) {v_p[k] = v_m[k] + (2.*lorentz_coeff/(1.+lorentz_coeff*lorentz_coeff)) * vcrosst[k];} // second half-rotation
                            for(k=0;k<3;k++) {v_p[k] += 0.5*efield_coeff*efield[k];} // half-step from E-field
                            /* calculate effective acceleration from discrete step in velocity */
                            for(k=0;k<3;k++) {external_forcing[k] += (v_p[k] - dv[k]) / dt;} // boris integrator
                            //for(k=0;k<3;k++) {external_forcing[k] += grain_charge_cinv * v_cross_B[k];} // standard explicit integrator
                            
                            /* note: if grains moving super-sonically with respect to gas, and charge equilibration time is much shorter than the
                             streaming/dynamical timescales, then the charge is slightly reduced, because the ion collision rate is increased while the
                             electron collision rate is increased less (since electrons are moving much faster, we assume the grain is still sub-sonic
                             relative to the electron sound speed. in this case, for the large-grain limit, the Draine & Sutin results can be generalized;
                             the full expressions are messy but can be -approximated- fairly well for Mach numbers ~3-30 by simply
                             suppressing the equilibrium grain charge by a power ~exp[-0.04*mach]  (weak effect, though can be significant for mach>10) */
#endif
                            double delta_egy = 0;
                            double delta_mom[3];
                            for(k=0; k<3; k++)
                            {
                                /* measure the imparted energy and momentum as if there were no external acceleration */
                                double v_init = P[i].Vel[k];
                                double vel_new = v_init + slow_fac * (P[i].Gas_Velocity[k]-v_init);
                                /* now calculate the updated velocity accounting for any external, non-standard accelerations */
                                double vdrift = 0;
                                if(tstop_inv > 0) {vdrift = external_forcing[k] / (tstop_inv * sqrt(1+x0*x0));}
                                dv[k] = slow_fac * (P[i].Gas_Velocity[k] - v_init + vdrift);
                                if(isnan(vdrift)||isnan(slow_fac)) {dv[k] = 0;}
                                
                                vel_new = v_init + dv[k];
                                delta_mom[k] = P[i].Mass * (vel_new - v_init);
                                delta_egy += 0.5*P[i].Mass * (vel_new*vel_new - v_init*v_init);
#ifdef GRAIN_BACKREACTION
                                P[i].Grain_DeltaMomentum[k] = delta_mom[k];
#endif
                                /* note, we can directly apply this by taking P[i].Vel[k] += dv[k]; but this is not as accurate as our
                                 normal leapfrog integration scheme.
                                 we can also account for the -gas- acceleration, by including it like vdrift;
                                 for a constant t_stop, the gas acceleration term appears as
                                 P[i].Vel[l] += Gas_Accel[k] * dt + slow_fac * (Gas-Accel[k] / tstop_inv) */
                                /* note that we solve the equations with an external acceleration already (external_forcing above): therefore add to forces
                                 like gravity that are acting on the gas and dust in the same manner (in terms of acceleration) */
                                P[i].GravAccel[k] += dv[k] / dt;
                                //P[i].Vel[k] += dv[k];
                            }
                        } // closes check for if(v_mag > 0)
                    } // closes check for if(dt > 0)
                } // closes check for if(P[i].Gas_Density > 0)
        } // closes check for if(P[i].Type != 0)
    } // closes main particle loop
    
#ifdef GRAIN_BACKREACTION
    grain_backrx(); /* call master routine to assign the back-reaction force among neighbors */
#endif
    CPU_Step[CPU_DRAGFORCE] += measure_time();
}




/* this is a template for fully-automated parallel (hybrid MPI+OpenMP/Pthreads) neighbor communication
 written in a completely modular fashion. this works as long as what you are trying to do isn't too
 complicated (from a communication point-of-view). You specify a few key variables, and then
 define the variables that need to be passed, and write the actual sub-routine that does the actual
 'work' between neighbors, but all of the parallelization, looping, communication blocks,
 etc, are all handled for you. */

#define MASTER_FUNCTION_NAME grain_backrx_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int MASTER_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if((P[i].Type==3)&&(P[i].TimeBin>=0)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */


/* this structure defines the variables that need to be sent -from- the 'searching' element */
struct INPUT_STRUCT_NAME
{
    int NodeList[NODELISTLENGTH]; MyDouble Pos[3], Hsml; /* these must always be defined */
#ifdef GRAIN_BACKREACTION
    double Grain_DeltaMomentum[3], Gas_Density;
#endif
}
*DATAIN_NAME, *DATAGET_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{   /* "i" is the particle from which data will be assigned, to structure "in" */
    int k; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];} /* good example - always needed */
    in->Hsml = PPP[i].Hsml; /* also always needed for search (can change the radius "PPP[i].Hsml" but in->Hsml must be defined */
#ifdef GRAIN_BACKREACTION
    for(k=0;k<3;k++) {in->Grain_DeltaMomentum[k]=P[i].Grain_DeltaMomentum[k];}
    in->Gas_Density = P[i].Gas_Density;
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
            if(numngb_inbox < 0) {return -1;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; if((P[j].Mass <= 0)||(P[j].Hsml <= 0)) {continue;} /* make sure neighbor is valid */
                int k; double dp[3]; for(k=0;k<3;k++) {dp[k]=local.Pos[k]-P[j].Pos[k];} /* position offset */
                NEAREST_XYZ(dp[0],dp[1],dp[2],1); double r2=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2]; /* box-wrap appropriately and calculate distance */
#ifdef BOX_BND_PARTICLES
                if(P[j].ID > 0) {r2 = -1;} /* ignore frozen boundary particles */
#endif
                if((r2>0)&&(r2<local.Hsml*local.Hsml)) /* only keep elements inside search radius */
                {
                    double hinv,hinv3,hinv4,wk_i=0,dwk_i=0,r=sqrt(r2); kernel_hinv(local.Hsml,&hinv,&hinv3,&hinv4);
                    kernel_main(r*hinv, hinv3, hinv4, &wk_i, &dwk_i, 0); /* kernel quantities that may be needed */
#ifdef GRAIN_BACKREACTION
                    double wt = -wk_i / local.Gas_Density; /* degy=wt*delta_egy; */
                    for(k=0;k<3;k++) {double dv=wt*local.Grain_DeltaMomentum[k]; P[j].Vel[k]+=dv; SphP[j].VelPred[k]+=dv;}
                    /* for(k=0;k<3;k++) {degy-=P[j].Mass*dv*(VelPred_j[k]+0.5*dv);} SphP[j].InternalEnergy += degy; */
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
    PRINT_STATUS(" ..assigning grain back-reaction to gas\n");
     //grain_backrx_initial_operations_preloop(); /* do initial pre-processing operations as needed before main loop [nothing needed here] */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    //grain_backrx_final_operations_and_cleanup(); /* do final operations on results [nothing needed here] */
    CPU_Step[CPU_DRAGFORCE] += timeall; /* collect timing information [here lumping it all together] */
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
double prob_of_grain_interaction(double cx_per_unitmass, double mass, double r, double h_si, double dV[3], integertime dt_step, int j_ngb)
{
    double dVmag = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]) / All.cf_atime; // velocity in physical
    double dt = dt_step * All.Timebase_interval / All.cf_hubble_a; // time in physical
    double rho_eff = 0.5*(mass + P[j_ngb].Mass) / (h_si*h_si*h_si) * All.cf_a3inv; // density in physical
    double cx_eff = g_geo(r/h_si) * (mass*cx_per_unitmass + P[j_ngb].Mass*return_grain_cross_section_per_unit_mass(j_ngb)) / (mass + P[j_ngb].Mass); // mass-weighted effective cross section (physical) scaled to cgs
    double units = All.UnitDensity_in_cgs * All.UnitLength_in_cm * All.HubbleParam; // needed to convert everything to cgs
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






#endif


