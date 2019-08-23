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
void GrainLoop_Master(void);
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
    GrainLoop_Master(); /* call master routine to assign the back-reaction force among neighbors */
#endif
    CPU_Step[CPU_DRAGFORCE] += measure_time();
}




/* this is a template for fully-automated parallel (hybrid MPI+OpenMP/Pthreads) neighbor communication
 written in a completely modular fashion. this works as long as what you are trying to do isn't too
 complicated (from a communication point-of-view). You specify a few key variables, and then
 define the variables that need to be passed, and write the actual sub-routine that does the actual
 'work' between neighbors, but all of the parallelization, looping, communication blocks,
 etc, are all handled for you. */

#define MASTER_FUNCTION_NAME GrainLoop_Master /* this function -must- be defined somewhere as "void MASTER_FUNCTION_NAME(void);" in order to actually be call-able by other parent routines! */
#define KERNEL_BITFLAG_DEFINITION_LOCAL 1 /* set to bitflag or name of sub-routine which returns bitflag for valid -neighbor- types */
#define CPU_COST_CODE_NAME CPU_DRAGFORCE /* cpu 'bin' name to charge costs to (for displaying/load-balancing) */
# include "../system/code_predefs_for_standard_neighbor_loops.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together */

/* this structure defines the variables that need to be sent -from- the 'searching' particle */
struct INPUT_STRUCT_NAME
{
    int NodeList[NODELISTLENGTH]; MyDouble Pos[3], Hsml; /* these must always be defined */
#ifdef GRAIN_BACKREACTION
    double Grain_DeltaMomentum[3], Gas_Density;
#endif
}
*DATAIN_NAME, *DATAGET_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' particle */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i);
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i)
{   /* "i" is the particle from which data will be assigned, to structure "in" */
    int k; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];} /* good example - always needed */
    in->Hsml = PPP[i].Hsml; /* also always needed for search (can change the radius "PPP[i].Hsml" but in->Hsml must be defined */
#ifdef GRAIN_BACKREACTION
    for(k=0;k<3;k++) {in->Grain_DeltaMomentum[k]=P[i].Grain_DeltaMomentum[k];}
    in->Gas_Density = P[i].Gas_Density;
#endif
}

/* this structure defines the variables that need to be sent -back to- the 'searching' particle */
struct OUTPUT_STRUCT_NAME
{ /* define variables below as e.g. "double X;" */
}
*DATARESULT_NAME, *DATAOUT_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' particle */
static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration);
static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{  /* "i" is the particle to which data from structure "out" will be assigned. mode=0 for local communication,
    =1 for data sent back from other processors. you must account for this. */
    /* example: ASSIGN_ADD(P[i].X,out->X,mode); which is short for: if(mode==0) {P[i].X=out->X;} else {P[i].X+=out->X;} */
}

/* this subroutine defines the conditions (true/false) for a particle to be considered 'active' and summon the evaluation subroutines */
int CONDITIONFUNCTION_FOR_EVALUATION(int i);
int CONDITIONFUNCTION_FOR_EVALUATION(int i) /* "i" is the index of the particle for which we evaluate */
{
    if(P[i].TimeBin<0) {return 0;}
    if(P[i].Type==3) {return 1;}
    return 0;
}

/* this subroutine does the actual neighbor-element calculations (this is the 'core' of the loop, essentially) */
static inline void NEIGHBOROPS_FUNCTION_NAME(struct INPUT_STRUCT_NAME *local, struct OUTPUT_STRUCT_NAME *out, int j, int loop_iteration);
static inline void NEIGHBOROPS_FUNCTION_NAME(struct INPUT_STRUCT_NAME *local, struct OUTPUT_STRUCT_NAME *out, int j, int loop_iteration)
{
    int k; double dp[3]; for(k=0;k<3;k++) {dp[k]=local->Pos[k]-P[j].Pos[k];} /* position */
#ifdef BOX_PERIODIC /* box-wrap appropriately  */
    NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif
    double r2=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2];
#ifdef BOX_BND_PARTICLES
    if(P[j].ID > 0) {r2 = -1;}
#endif
    if((r2>0)&&(r2<local->Hsml*local->Hsml)) /* only keep elements inside search radius */
    {
        double hinv,hinv3,hinv4,wk_i=0,dwk_i=0,r=sqrt(r2); kernel_hinv(local->Hsml,&hinv,&hinv3,&hinv4);
        kernel_main(r*hinv, hinv3, hinv4, &wk_i, &dwk_i, 0); /* kernel quantities that may be needed */
        /*
         double VelPred_j[3]; for(k=0;k<3;k++) {VelPred_j[k]=SphP[j].VelPred[k];}
         #ifdef BOX_SHEARING
         if(local->Pos[0]-P[j].Pos[0] > +boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
         if(local->Pos[0]-P[j].Pos[0] < -boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
         #endif
         */
#ifdef GRAIN_BACKREACTION
        double wt = -wk_i / local->Gas_Density;/*, degy = wt*delta_egy; */
        for(k=0;k<3;k++)
        {
            double dv=wt*local->Grain_DeltaMomentum[k]; P[j].Vel[k]+=dv; SphP[j].VelPred[k]+=dv;
            /*            degy-=P[j].Mass*dv*(VelPred_j[k]+0.5*dv); */
        }
        /*      SphP[j].InternalEnergy += degy; */
#endif
    }
}

/* this subroutine computes any final operations on the particles after the main loop is completed */
static inline void FINAL_OPERATIONS_FUNCTION_NAME(int i, int loop_iteration);
static inline void FINAL_OPERATIONS_FUNCTION_NAME(int i, int loop_iteration)
{ /* here "i" is the particle for which we've returned data -- do what you like to it! */
}

/* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
#include "../system/code_block_for_standard_neighbor_loops.h"







#ifdef GRAIN_COLLISIONS

void grain_collisions(void)
{
    int i,k;
    double dt, dvel, degy, soundspeed, R_grain, t_stop, slow_fac, vel_new, delta_mom[3], delta_egy;
    int N_MAX_KERNEL,N_MIN_KERNEL,MAXITER_FB,NITER,startnode,dummy,numngb_inbox,jnearest,i,j,k,n;
    double *pos,h,h2,hinv,hinv3,r2,rho,u,wk;
    Ngblist = (int *) mymalloc(NumPart * sizeof(int));
    
    CPU_Step[CPU_MISC] += measure_time();
#ifndef IO_REDUCED_MODE
    if(ThisTask == 0)
    {
        printf("Beginning Grain Collisions & Interactions \n");
    }
#endif
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 3)
        {
            if(P[i].Grain_Density > 0)
            {
#ifndef WAKEUP
                dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
                dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
                if(dt > 0)
                {
                    soundspeed = Particle_effective_soundspeed_i(i);
                    R_grain = P[i].Grain_Size; // grain size in --cgs-- units //
                } // closes check for if(dt > 0)
            } // closes check for if(P[i].Gas_Density > 0)
        } // closes check for if(P[i].Type != 0)
    } // closes main particle loop
    myfree(Ngblist);
    CPU_Step[CPU_DRAGFORCE] += measure_time();
} /* closes grain_collisions routine */




/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct grain_densdata_in
{
    MyDouble Pos[3];
    MyFloat Vel[3];
    MyFloat Hsml;
    int NodeList[NODELISTLENGTH];
}
*GrnDensDataIn, *GrnDensDataGet;

static struct grain_densdata_out
{
    MyLongDouble RhoGrains;
    MyLongDouble GrainVel[3];
}
*GrnDensDataResult, *GrnDensDataOut;


/*
 calculates density and velocity of surrounding grain particles;
 currently is the full routine allowing for MPI communications, etc;
 but have commented out the iteration to solve for the neighbor number (we just use Hsml
 as defined by the SPH; could easily make this free to iterate itself, at small cost)
 */
void grain_density(void)
{
    MyFloat *Left, *Right;
    int i, j, ndone, ndone_flag, npleft, dummy, iter = 0;
    integertime  dt_step;
    int ngrp, sendTask, recvTask, place, nexport, nimport;
    long long ntot;
    double fac;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
    double timecomp, timecomm, timewait;
    double tstart, tend, t0, t1;
    double desnumngb, valuenorm;
    
    CPU_Step[CPU_DENSMISC] += measure_time();
    
    Left = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc(NumPart * sizeof(MyFloat));
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(grain_density_isactive(i))
        {
            Left[i] = Right[i] = 0;
        }
    }
    
    /* allocate buffers to arrange communication */
    Ngblist = (int *) mymalloc(NumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                           sizeof(struct grain_densdata_in) + sizeof(struct grain_densdata_out) +
                                                           sizemax(sizeof(struct grain_densdata_in),sizeof(struct grain_densdata_out))));
    DataIndexTable = (struct data_index *) mymalloc(All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc(All.BunchSize * sizeof(struct data_nodelist));
    
    t0 = my_second();
    desnumngb = All.DesNumNgb;
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        i = FirstActiveParticle;    /* begin with this index */
        do
        {
            for(j = 0; j < NTask; j++)
            {
                Send_count[j] = 0;
                Exportflag[j] = -1;
            }
            /* do local particles and prepare export list */
            tstart = my_second();
            for(nexport = 0; i >= 0; i = NextActiveParticle[i])
            {
                if(grain_density_isactive(i))
                {
                    if(grain_density_evaluate(i, 0, &nexport, Send_count) < 0)
                        break;
                }
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
            MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
            
            tstart = my_second();
            MPI_Allgather(Send_count, NTask, MPI_INT, Sendcount_matrix, NTask, MPI_INT, MPI_COMM_WORLD);
            tend = my_second();
            timewait1 += timediff(tstart, tend);
            
            for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
            {
                Recv_count[j] = Sendcount_matrix[j * NTask + ThisTask];
                nimport += Recv_count[j];
                if(j > 0)
                {
                    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                    Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }
            GrnDensDataGet = (struct grain_densdata_in *) mymalloc(nimport * sizeof(struct grain_densdata_in));
            GrnDensDataIn = (struct grain_densdata_in *) mymalloc(nexport * sizeof(struct grain_densdata_in));
            
            /* prepare particle data for export */
            for(j = 0; j < nexport; j++)
            {
                place = DataIndexTable[j].Index;
                GrnDensDataIn[j].Pos[0] = P[place].Pos[0];
                GrnDensDataIn[j].Pos[1] = P[place].Pos[1];
                GrnDensDataIn[j].Pos[2] = P[place].Pos[2];
                GrnDensDataIn[j].Hsml = PPP[place].Hsml;
                memcpy(GrnDensDataIn[j].NodeList,
                       DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
            }
            /* exchange particle data */
            tstart = my_second();
            for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;
                
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                        /* get the particles */
                        MPI_Sendrecv(&GrnDensDataIn[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct grain_densdata_in), MPI_BYTE,
                                     recvTask, TAG_GRDENS_A,
                                     &GrnDensDataGet[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct grain_densdata_in), MPI_BYTE,
                                     recvTask, TAG_GRDENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second();
            timecommsumm1 += timediff(tstart, tend);
            myfree(GrnDensDataIn);
            GrnDensDataResult = (struct grain_densdata_out *) mymalloc(nimport * sizeof(struct grain_densdata_out));
            GrnDensDataOut = (struct grain_densdata_out *) mymalloc(nexport * sizeof(struct grain_densdata_out));
            
            /* now do the particles that were sent to us */
            tstart = my_second();
            for(j = 0; j < nimport; j++)
                grain_density_evaluate(j, 1, &dummy, &dummy);
            tend = my_second();
            timecomp2 += timediff(tstart, tend);
            
            if(i < 0)
                ndone_flag = 1;
            else
                ndone_flag = 0;
            
            tstart = my_second();
            MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            tend = my_second();
            timewait2 += timediff(tstart, tend);
            
            /* get the result */
            tstart = my_second();
            for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
                sendTask = ThisTask;
                recvTask = ThisTask ^ ngrp;
                if(recvTask < NTask)
                {
                    if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                        /* send the results */
                        MPI_Sendrecv(&GrnDensDataResult[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct grain_densdata_out),
                                     MPI_BYTE, recvTask, TAG_GRDENS_B,
                                     &GrnDensDataOut[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct grain_densdata_out),
                                     MPI_BYTE, recvTask, TAG_GRDENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second();
            timecommsumm2 += timediff(tstart, tend);
            
            /* add the result to the local particles */
            tstart = my_second();
            for(j = 0; j < nexport; j++)
            {
                place = DataIndexTable[j].Index;
                if(P[place].Type == 3)
                {
                    P[place].Grain_Density += GrnDensDataOut[j].RhoGrains;
                    P[place].Grain_Velocity[0] += GrnDensDataOut[j].GrainVel[0];
                    P[place].Grain_Velocity[1] += GrnDensDataOut[j].GrainVel[1];
                    P[place].Grain_Velocity[2] += GrnDensDataOut[j].GrainVel[2];
                }
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
            myfree(GrnDensDataOut);
            myfree(GrnDensDataResult);
            myfree(GrnDensDataGet);
        }
        while(ndone < NTask);
        
        
        /* do final operations on results */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(grain_density_isactive(i))
            {
                if(P[i].Grain_Density > 0)
                {
                    P[i].Grain_Velocity[0] /= P[i].Grain_Density;
                    P[i].Grain_Velocity[1] /= P[i].Grain_Density;
                    P[i].Grain_Velocity[2] /= P[i].Grain_Density;
                }
                
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            sumup_large_ints(1, &npleft, &ntot);
        }
    } while(ntot > 0); /* closes the check for particles that need iteration */
    
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    //myfree(Right);
    //myfree(Left);
    
    /* mark as active again */
    /*
     for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
     if(P[i].TimeBin < 0)
     P[i].TimeBin = -P[i].TimeBin - 1;
     */
    
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    
    CPU_Step[CPU_DENSCOMPUTE] += timecomp;
    CPU_Step[CPU_DENSWAIT] += timewait;
    CPU_Step[CPU_DENSCOMM] += timecomm;
    CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
    
} // done with routine!





/*! core of sph density computation, adapted to search for grains now */
int grain_density_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
    int j, n;
    int startnode, numngb, numngb_inbox, listindex = 0;
    double h, h2, fac, hinv, hinv3, hinv4, wk, dwk;
    double dx, dy, dz, r, r2, u, mass_j;
    MyLongDouble sum_variable;
    MyLongDouble rho;
    MyLongDouble weighted_numngb;
    MyLongDouble gasvel[3];
    MyDouble *pos;
    MyFloat *vel;
    gasvel[0] = gasvel[1] = gasvel[2] = 0;
    rho = weighted_numngb = 0;
    
    if(mode == 0)
    {
        pos = P[target].Pos;
        h = PPP[target].Hsml;
        vel = P[target].Vel;
    }
    else
    {
        pos = GrnDensDataGet[target].Pos;
        vel = GrnDensDataGet[target].Vel;
        h = GrnDensDataGet[target].Hsml;
    }
    
    h2 = h * h;
    hinv = 1.0 / h;
#if (NUMDIMS==1)
    hinv3 = hinv / (boxSize_Y * boxSize_Z);
#endif
#if (NUMDIMS==2)
    hinv3 = hinv * hinv / boxSize_Z;
#endif
#if (NUMDIMS==3)
    hinv3 = hinv * hinv * hinv;
#endif
    hinv4 = hinv3 * hinv;
    
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = GrnDensDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }
    
    numngb = 0;
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_targeted(pos, h, target, &startnode, mode, nexport, nsend_local, 8); // 8 = 2^3: search for particles of type=3
            if(numngb_inbox < 0) return -1;
            
            for(n = 0; n < numngb_inbox; n++)
            {
                j = Ngblist[n];
                if(P[j].Mass == 0) continue;
                dx = pos[0] - P[j].Pos[0];
                dy = pos[1] - P[j].Pos[1];
                dz = pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC            /*  now find the closest image in the given box size  */
                NEAREST_XYZ(dx,dy,dz,1);
#endif
                r2 = dx*dx + dy*dy + dz*dz;
                if(r2 < h2)
                {
                    numngb++;
                    r = sqrt(r2);
                    u = r * hinv;
                    kernel_main(u,hinv3,hinv4,&wk,&dwk,0);
                    mass_j = P[j].Mass;
                    rho += FLT(mass_j * wk);
                    weighted_numngb += FLT(NORM_COEFF * wk / hinv3);
                    MyDouble VelPred_j[3];
                    for(k=0;k<3;k++) {VelPred_j[k]=SphP[j].VelPred[k];}
#ifdef BOX_SHEARING
                    if(pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                    if(pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
                    gasvel[0] += FLT(mass_j * wk * VelPred_j[0]);
                    gasvel[1] += FLT(mass_j * wk * VelPred_j[1]);
                    gasvel[2] += FLT(mass_j * wk * VelPred_j[2]);
                }
            }
        }
    }
    
    if(mode == 1)
    {
        listindex++;
        if(listindex < NODELISTLENGTH)
        {
            startnode = GrnDensDataGet[target].NodeList[listindex];
            if(startnode >= 0)
                startnode = Nodes[startnode].u.d.nextnode;    /* open it */
        }
    }
    if(mode == 0)
    {
        P[target].Grain_Density = rho;
        P[target].Grain_Velocity[0] = 0;
        P[target].Grain_Velocity[1] = 0;
        P[target].Grain_Velocity[2] = 0;
    }
    else
    {
        GrnDensDataResult[target].RhoGrains = rho;
        GrnDensDataResult[target].GrainVel[0] = 0;
        GrnDensDataResult[target].GrainVel[1] = 0;
        GrnDensDataResult[target].GrainVel[2] = 0;
    }
    
    return 0;
}



/* code to tell the grain density routine which particles to use */
int grain_density_isactive(int n)
{
    if(P[n].TimeBin < 0) return 0;
    if(P[n].Type == 3) return 1;
    return 0;
}



#endif // closes GRAIN_COLLISIONS



#endif


