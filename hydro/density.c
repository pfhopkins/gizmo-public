#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#include "../mesh_motion.h"

/*! \file density.c
 *  \brief hydro kernel size and neighbor determination, volumetric quantities calculated
 *
 *  This file contains the "first hydro loop", where the gas densities and some
 *  auxiliary quantities are computed.  There is also functionality that corrects the kernel length if needed.
 */
/*!
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * The code has been modified substantially (condensed, different criteria for kernel lengths, optimizatins,
 * rewritten parallelism, new physics included, new variable/memory conventions added) by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

struct kernel_density /*! defines a number of useful variables we will use below */
{
  double dp[3],dv[3],r, wk, dwk, hinv, hinv3, hinv4, mj_wk, mj_dwk_r;
};


/*! routine to determine if a given element is actually going to be active in the density subroutines below */
int density_isactive(int n)
{
    /* first check our 'marker' for particles which have finished iterating to an Hsml solution (if they have, dont do them again) */
    if(P[n].TimeBin < 0) {return 0;}

#if defined(GRAIN_FLUID)
    if((1 << P[n].Type) & (GRAIN_PTYPES)) {return 1;} /* any of the particle types flagged as a valid grain-type is active here */
#endif

#if defined(BH_INTERACT_ON_GAS_TIMESTEP)
    if(P[n].Type == 5){if(!P[n].do_gas_search_this_timestep && All.Ti_Current > 0) return 0;} /* not enough time has elapsed since the last gas interaction */
#endif
#if defined(RT_SOURCE_INJECTION)
    if((1 << P[n].Type) & (RT_SOURCES))
    {
#if defined(GALSF)
       if(((P[n].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[n].Type == 2)||(P[n].Type==3))))&&(P[n].Mass>0))
        {
            double star_age = evaluate_stellar_age_Gyr(P[n].StellarAge);
            if((star_age < 0.1)&&(star_age > 0)&&(!isnan(star_age))) return 1;
        }
#else
        if(Flag_FullStep) {return 1;} // only do on full timesteps
#endif
    }
#endif

#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    if(((P[n].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[n].Type == 2)||(P[n].Type==3))))&&(P[n].Mass>0))
    {
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
        /* check if there is going to be a SNe this timestep, in which case, we want the density info! */
        if(P[n].SNe_ThisTimeStep>0) return 1;
#endif
        
#if defined(GALSF)
        if(P[n].DensAroundStar<=0) return 1;
        if(All.ComovingIntegrationOn==0) // only do stellar age evaluation if we have to //
        {
            double star_age = evaluate_stellar_age_Gyr(P[n].StellarAge);
            if(star_age < 0.035) return 1;
        }
#endif
#if (defined(GRAIN_FLUID) || defined(RADTRANSFER)) && (!defined(GALSF) && !(defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)))
        return 1;
#endif
    }
#endif

#ifdef BLACK_HOLES
    if(P[n].Type == 5) return 1;
#endif

    if(P[n].Type == 0 && P[n].Mass > 0) return 1;
    return 0; /* default to 0 if no check passed */
}




#define CORE_FUNCTION_NAME density_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME hydrokerneldensity_particle2in    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME hydrokerneldensity_out2particle  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/*! this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
  MyDouble Pos[3];
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
  MyFloat Accel[3];
#endif
  MyFloat Vel[3];
  MyFloat Hsml;
#ifdef GALSF_SUBGRID_WINDS
  MyFloat DelayTime;
#endif
  int NodeList[NODELISTLENGTH];
  int Type;
}
 *DATAIN_NAME, *DATAGET_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
void hydrokerneldensity_particle2in(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k;
    in->Type = P[i].Type;
    in->Hsml = PPP[i].Hsml;
    for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    for(k=0;k<3;k++) {if(P[i].Type==0) {in->Vel[k]=SphP[i].VelPred[k];} else {in->Vel[k]=P[i].Vel[k];}}
    if(P[i].Type == 0)
    {
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        for(k=0;k<3;k++) {in->Accel[k] = All.cf_a2inv*P[i].GravAccel[k] + SphP[i].HydroAccel[k];} // PHYSICAL units //
#endif
#ifdef GALSF_SUBGRID_WINDS
        in->DelayTime = SphP[i].DelayTime;
#endif
    }
}

/*! this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble Ngb;
    MyLongDouble Rho;
    MyLongDouble DhsmlNgb;
    MyLongDouble Particle_DivVel;
    MyLongDouble NV_T[3][3];
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
    MyDouble ParticleVel[3];
#endif
#ifdef HYDRO_SPH
    MyLongDouble DhsmlHydroSumFactor;
#endif
#ifdef RT_SOURCE_INJECTION
    MyLongDouble KernelSum_Around_RT_Source;
#endif
#ifdef HYDRO_PRESSURE_SPH
    MyLongDouble EgyRho;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
    MyFloat AGS_zeta;
#endif
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    MyFloat NV_D[3][3];
    MyFloat NV_A[3][3];
#endif
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    MyFloat GradRho[3];
#endif
#if defined(BLACK_HOLES)
    int BH_TimeBinGasNeighbor;
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
    MyDouble BH_dr_to_NearestGasNeighbor;
#endif
#endif
#if defined(TURB_DRIVING) || defined(GRAIN_FLUID)
    MyDouble GasVel[3];
#endif
#if defined(GRAIN_FLUID)
    MyDouble Gas_InternalEnergy;
#if defined(GRAIN_LORENTZFORCE)
    MyDouble Gas_B[3];
#endif
#endif
}
 *DATARESULT_NAME, *DATAOUT_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
void hydrokerneldensity_out2particle(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    int j,k;
    ASSIGN_ADD(PPP[i].NumNgb, out->Ngb, mode);
    ASSIGN_ADD(PPP[i].DhsmlNgbFactor, out->DhsmlNgb, mode);
    ASSIGN_ADD(P[i].Particle_DivVel, out->Particle_DivVel,   mode);

    if(P[i].Type == 0)
    {
        ASSIGN_ADD(SphP[i].Density, out->Rho, mode);
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
        for(k=0;k<3;k++) ASSIGN_ADD(SphP[i].ParticleVel[k], out->ParticleVel[k],   mode);
#endif
        for(k=0;k<3;k++) {for(j=0;j<3;j++) {ASSIGN_ADD(SphP[i].NV_T[k][j], out->NV_T[k][j], mode);}}

#ifdef HYDRO_SPH
        ASSIGN_ADD(SphP[i].DhsmlHydroSumFactor, out->DhsmlHydroSumFactor, mode);
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
        ASSIGN_ADD(PPPZ[i].AGS_zeta, out->AGS_zeta,   mode);
#endif

#ifdef HYDRO_PRESSURE_SPH
        ASSIGN_ADD(SphP[i].EgyWtDensity,   out->EgyRho,   mode);
#endif

#if defined(TURB_DRIVING)
        for(k = 0; k < 3; k++) {ASSIGN_ADD(SphP[i].SmoothedVel[k], out->GasVel[k], mode);}
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        for(k = 0; k < 3; k++)
            for(j = 0; j < 3; j++)
            {
                ASSIGN_ADD(SphP[i].NV_D[k][j], out->NV_D[k][j], mode);
                ASSIGN_ADD(SphP[i].NV_A[k][j], out->NV_A[k][j], mode);
            }
#endif
    } // P[i].Type == 0 //

#if defined(GRAIN_FLUID)
    if((1 << P[i].Type) & (GRAIN_PTYPES))
    {
        ASSIGN_ADD(P[i].Gas_Density, out->Rho, mode);
        ASSIGN_ADD(P[i].Gas_InternalEnergy, out->Gas_InternalEnergy, mode);
        for(k = 0; k<3; k++) {ASSIGN_ADD(P[i].Gas_Velocity[k], out->GasVel[k], mode);}
#if defined(GRAIN_LORENTZFORCE)
        for(k = 0; k<3; k++) {ASSIGN_ADD(P[i].Gas_B[k], out->Gas_B[k], mode);}
#endif
    }
#endif

#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    ASSIGN_ADD(P[i].DensAroundStar, out->Rho, mode);
    for(k = 0; k<3; k++) {ASSIGN_ADD(P[i].GradRho[k], out->GradRho[k], mode);}
#endif

#if defined(RT_SOURCE_INJECTION)
#if defined(RT_BH_ANGLEWEIGHT_PHOTON_INJECTION)
    if(All.TimeStep == 0) // we only do this on the 0'th timestep, since we haven't done a BH loop yet to get the angle weights we'll use normally
#endif    
    if((1 << P[i].Type) & (RT_SOURCES)) {ASSIGN_ADD(P[i].KernelSum_Around_RT_Source, out->KernelSum_Around_RT_Source, mode);}
#endif

#ifdef BLACK_HOLES
    if(P[i].Type == 5)
    {
        if(mode == 0) {BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;} else {if(BPP(i).BH_TimeBinGasNeighbor > out->BH_TimeBinGasNeighbor) {BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;}}
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
        if(mode == 0) {BPP(i).BH_dr_to_NearestGasNeighbor = out->BH_dr_to_NearestGasNeighbor;} else {if(BPP(i).BH_dr_to_NearestGasNeighbor > out->BH_dr_to_NearestGasNeighbor) {BPP(i).BH_dr_to_NearestGasNeighbor = out->BH_dr_to_NearestGasNeighbor;}}
#endif
    } /* if(P[i].Type == 5) */
#endif
}

/*! declare this utility function here now that the relevant structures it uses have been defined */
void density_evaluate_extra_physics_gas(struct INPUT_STRUCT_NAME *local, struct OUTPUT_STRUCT_NAME *out, struct kernel_density *kernel, int j);


/*! This function represents the core of the initial hydro kernel-identification and volume computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine should in general contain no writes to shared memory. for optimization reasons, a couple of such writes have been included here in the sub-code for some sink routines -- those need to be handled with special care, both for thread safety and because of iteration. in general writes to shared memory in density.c are strongly discouraged -- */
int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, n, startnode, numngb_inbox, listindex = 0; double r2, h2, u, mass_j, wk;
    struct kernel_density kernel; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    if(mode == 0) {hydrokerneldensity_particle2in(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    h2 = local.Hsml * local.Hsml; kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
#if defined(BLACK_HOLES)
    out.BH_TimeBinGasNeighbor = TIMEBINS;
#ifdef BH_ACCRETE_NEARESTFIRST
    out.BH_dr_to_NearestGasNeighbor = MAX_REAL_NUMBER;
#endif
#endif
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
#ifdef GALSF_SUBGRID_WINDS /* check if partner is a wind particle: if I'm not wind, then ignore the wind particle */
                if(SphP[j].DelayTime > 0) {if(!(local.DelayTime > 0)) {continue;}}
#endif
                if(P[j].Mass <= 0) continue;
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                if(r2 < h2) /* this loop is only considering particles inside local.Hsml, i.e. seen-by-main */
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                    mass_j = P[j].Mass;
                    kernel.mj_wk = FLT(mass_j * kernel.wk);

                    out.Ngb += kernel.wk;
                    out.Rho += kernel.mj_wk;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                    if(local.Type == 0 && kernel.r==0) {int kv; for(kv=0;kv<3;kv++) {out.ParticleVel[kv] += kernel.mj_wk * SphP[j].VelPred[kv];}} // just the self-contribution //
#endif
#if defined(RT_SOURCE_INJECTION)
#if defined(RT_BH_ANGLEWEIGHT_PHOTON_INJECTION)
                    if(All.TimeStep == 0) // we only do this on the 0'th timestep, since we haven't done a BH loop yet to get the angle weights we'll use normally
#endif                        
                    if((1 << local.Type) & (RT_SOURCES)) {out.KernelSum_Around_RT_Source += 1.-u*u;}
#endif
                    out.DhsmlNgb += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#ifdef HYDRO_SPH
                    double mass_eff = mass_j;
#ifdef HYDRO_PRESSURE_SPH
                    mass_eff *= SphP[j].InternalEnergyPred;
                    out.EgyRho += kernel.wk * mass_eff;
#endif
                    out.DhsmlHydroSumFactor += -mass_eff * (NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
                    if(local.Type == 0) {out.AGS_zeta += mass_j * kernel_gravity(u, kernel.hinv, kernel.hinv3, 0);}
#endif
                    /* for everything below, we do NOT include the particle self-contribution! */
                    if(kernel.r > 0)
                    {
                        if(local.Type == 0)
                        {
                            wk = kernel.wk; /* MAKE SURE THIS MATCHES CHOICE IN GRADIENTS.c!!! */
                            /* the weights for the MLS tensor used for gradient estimation */
                            out.NV_T[0][0] +=  wk * kernel.dp[0] * kernel.dp[0];
                            out.NV_T[0][1] +=  wk * kernel.dp[0] * kernel.dp[1];
                            out.NV_T[0][2] +=  wk * kernel.dp[0] * kernel.dp[2];
                            out.NV_T[1][1] +=  wk * kernel.dp[1] * kernel.dp[1];
                            out.NV_T[1][2] +=  wk * kernel.dp[1] * kernel.dp[2];
                            out.NV_T[2][2] +=  wk * kernel.dp[2] * kernel.dp[2];
                            /* these will temporarily hold the 'face area' terms */
                            out.NV_T[1][0] += wk * kernel.dp[0];
                            out.NV_T[2][0] += wk * kernel.dp[1];
                            out.NV_T[2][1] += wk * kernel.dp[2];
                        }
                        kernel.dv[0] = local.Vel[0] - SphP[j].VelPred[0];
                        kernel.dv[1] = local.Vel[1] - SphP[j].VelPred[1];
                        kernel.dv[2] = local.Vel[2] - SphP[j].VelPred[2];
                        NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,kernel.dv,1); /* wrap velocities for shearing boxes if needed */
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                        // do neighbor contribution to smoothed particle velocity here, after wrap, so can account for shearing boxes correctly //
                        {int kv; for(kv=0;kv<3;kv++) {out.ParticleVel[kv] += kernel.mj_wk * (local.Vel[kv] - kernel.dv[kv]);}}
#endif
                        out.Particle_DivVel -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how Hsml will evolve (particle drift) */

                        density_evaluate_extra_physics_gas(&local, &out, &kernel, j);
                    } // kernel.r > 0
                } // if(r2 < h2)
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {hydrokerneldensity_out2particle(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}



/*! this is an extra function to simplify additional computations within the kernel that need to be done as part of the evaluation above */
void density_evaluate_extra_physics_gas(struct INPUT_STRUCT_NAME *local, struct OUTPUT_STRUCT_NAME *out, struct kernel_density *kernel, int j)
{
    kernel->mj_dwk_r = P[j].Mass * kernel->dwk / kernel->r;

    if(local->Type != 0)
    {
#if defined(GRAIN_FLUID)
        if((1 << local->Type) & (GRAIN_PTYPES))
        {
            out->Gas_InternalEnergy += kernel->mj_wk * SphP[j].InternalEnergyPred;
            out->GasVel[0] += kernel->mj_wk * (local->Vel[0]-kernel->dv[0]);
            out->GasVel[1] += kernel->mj_wk * (local->Vel[1]-kernel->dv[1]);
            out->GasVel[2] += kernel->mj_wk * (local->Vel[2]-kernel->dv[2]);
#if defined(GRAIN_LORENTZFORCE)
            out->Gas_B[0] += kernel->wk * SphP[j].BPred[0];
            out->Gas_B[1] += kernel->wk * SphP[j].BPred[1];
            out->Gas_B[2] += kernel->wk * SphP[j].BPred[2];
#endif
        }
#endif

#if defined(BLACK_HOLES)
        /* note, we will have some writes to shared memory here for some initializations of 'j' quantities. fortunately these do not depend on previous values of those quantities, so can be done thread-safely with minor edits using the constructs below */
        if(local->Type == 5)
        {
            #pragma omp atomic write
            P[j].SwallowID = 0;  // this way we don't have to do a global loop over local particles in blackhole_accretion() to reset these quantities...
#ifdef SINGLE_STAR_SINK_DYNAMICS
            #pragma omp atomic write
            P[j].SwallowTime = MAX_REAL_NUMBER; // initialize as a large number before looking
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
            if(kernel->r < DMAX(P[j].Hsml, All.ForceSoftening[5])) {
                #pragma omp atomic write
                P[j].BH_Ngb_Flag = 1; // note that this particle is inside of a BH's kernel function
            }
#endif
            short int TimeBin_j = P[j].TimeBin; if(TimeBin_j < 0) {TimeBin_j = -TimeBin_j - 1;} // need to make sure we correct for the fact that TimeBin is used as a 'switch' here to determine if a particle is active for iteration, otherwise this gives nonsense!
            if(out->BH_TimeBinGasNeighbor > TimeBin_j) {out->BH_TimeBinGasNeighbor = TimeBin_j;}
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
            double dr_eff_wtd = Get_Particle_Size(j);
            dr_eff_wtd=sqrt(dr_eff_wtd*dr_eff_wtd + (kernel->r)*(kernel->r)); /* effective distance for Gaussian-type kernel, weighted by density */
            if((dr_eff_wtd < out->BH_dr_to_NearestGasNeighbor) && (P[j].Mass > 0)) {out->BH_dr_to_NearestGasNeighbor = dr_eff_wtd;}
#endif
        }
#endif // BLACK_HOLES
        
    } else { /* local.Type == 0 */

#if defined(TURB_DRIVING)
        out->GasVel[0] += kernel->mj_wk * (local->Vel[0]-kernel->dv[0]);
        out->GasVel[1] += kernel->mj_wk * (local->Vel[1]-kernel->dv[1]);
        out->GasVel[2] += kernel->mj_wk * (local->Vel[2]-kernel->dv[2]);
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        double wk = kernel->wk;
        out->NV_A[0][0] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dp[0] * wk;
        out->NV_A[0][1] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dp[1] * wk;
        out->NV_A[0][2] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dp[2] * wk;
        out->NV_A[1][0] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dp[0] * wk;
        out->NV_A[1][1] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dp[1] * wk;
        out->NV_A[1][2] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dp[2] * wk;
        out->NV_A[2][0] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dp[0] * wk;
        out->NV_A[2][1] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dp[1] * wk;
        out->NV_A[2][2] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dp[2] * wk;

        out->NV_D[0][0] += kernel->dv[0] * kernel->dp[0] * wk;
        out->NV_D[0][1] += kernel->dv[0] * kernel->dp[1] * wk;
        out->NV_D[0][2] += kernel->dv[0] * kernel->dp[2] * wk;
        out->NV_D[1][0] += kernel->dv[1] * kernel->dp[0] * wk;
        out->NV_D[1][1] += kernel->dv[1] * kernel->dp[1] * wk;
        out->NV_D[1][2] += kernel->dv[1] * kernel->dp[2] * wk;
        out->NV_D[2][0] += kernel->dv[2] * kernel->dp[0] * wk;
        out->NV_D[2][1] += kernel->dv[2] * kernel->dp[1] * wk;
        out->NV_D[2][2] += kernel->dv[2] * kernel->dp[2] * wk;
#endif

    } // Type = 0 check

#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    /* this is here because for the models of BH growth and self-shielding of stars, we
     just need a quick-and-dirty, single-pass approximation for the gradients (the error from
     using this as opposed to the higher-order gradient estimators is small compared to the
     Sobolev approximation): use only for -non-gas- particles */
    out->GradRho[0] += kernel->mj_dwk_r * kernel->dp[0];
    out->GradRho[1] += kernel->mj_dwk_r * kernel->dp[1];
    out->GradRho[2] += kernel->mj_dwk_r * kernel->dp[2];
#endif

}




/*! This function computes the local neighbor kernel for each active hydro element, the number of neighbours in the current kernel radius, and the divergence
 * and rotation of the velocity field.  This is used then to compute the effective volume of the element in MFM/MFV/SPH-type methods, which is then used to
 * update volumetric quantities like density and pressure. The routine iterates to attempt to find a target kernel size set adaptively -- see code user guide for details
 */
void density(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second(); MyFloat *Left, *Right; double fac, fac_lim, desnumngb, desnumngbdev; long long ntot;
    int i, npleft, iter=0, redo_particle, particle_set_to_minhsml_flag = 0, particle_set_to_maxhsml_flag = 0;
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
    
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES /* define a variable for below to know which stellar types qualify here */
    int valid_stellar_types = 2+4+8+16, invalid_stellar_types = 1+32; // allow types 1,2,3,4 here //
#if (defined(GRAIN_FLUID) || defined(RADTRANSFER)) && (!defined(GALSF) && !(defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)))
    valid_stellar_types = 16; invalid_stellar_types = 1+2+4+8+32; // -only- type-4 sources in these special problems
#ifdef RADTRANSFER
    invalid_stellar_types = 64; valid_stellar_types = RT_SOURCES; // any valid 'injection' source is allowed
#endif
#ifdef GRAIN_FLUID
    invalid_stellar_types = GRAIN_PTYPES;
#endif
#endif
#endif
    
    /* initialize anything we need to about the active particles before their loop */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if(density_isactive(i)) {
            Left[i] = Right[i] = 0;
#ifdef BLACK_HOLES
            P[i].SwallowID = 0;
#ifdef SINGLE_STAR_SINK_DYNAMICS
            P[i].SwallowTime = MAX_REAL_NUMBER;
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
            P[i].BH_Ngb_Flag = 0;
#endif
#endif
            double maxsoft = All.MaxHsml; /* before the first pass, need to ensure the particles do not exceed the maximum Hsml allowed */
#if defined(DO_DENSITY_AROUND_STAR_PARTICLES) && defined(GALSF)
            if( ((1 << P[i].Type) & (valid_stellar_types)) && !((1 << P[i].Type) & (invalid_stellar_types)) ) {maxsoft = 2.0 / (UNIT_LENGTH_IN_KPC*All.cf_atime);}
#endif
#ifdef BLACK_HOLES
            if(P[i].Type == 5) {maxsoft = All.BlackHoleMaxAccretionRadius / All.cf_atime;}  // MaxAccretionRadius is now defined in params.txt in PHYSICAL units
#endif
            if((PPP[i].Hsml < 0) || !isfinite(PPP[i].Hsml) || (PPP[i].Hsml > 0.99*maxsoft)) {PPP[i].Hsml = 0.99*maxsoft;} /* don't set to exactly maxsoft because our looping below won't treat this correctly */
            
        }} /* done with intial zero-out loop */
    desnumngb = All.DesNumNgb; desnumngbdev = All.MaxNumNgbDeviation;
    /* in the initial timestep and iteration, use a much more strict tolerance for the neighbor number */
    if(All.Time==All.TimeBegin) {if(All.MaxNumNgbDeviation > 0.05) desnumngbdev=0.05;}
    MyLongDouble desnumngbdev_0 = desnumngbdev, Tinv[3][3], ConditionNumber=0; int k,k1,k2; k=0;

    /* allocate buffers to arrange communication */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */

        /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        double tstart = my_second(), tend;
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(density_isactive(i))
            {
                if(PPP[i].NumNgb > 0)
                {
                    PPP[i].DhsmlNgbFactor *= PPP[i].Hsml / (NUMDIMS * PPP[i].NumNgb);
                    P[i].Particle_DivVel /= PPP[i].NumNgb;
                    /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                    PPP[i].NumNgb *= NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS);
                } else {
                    PPP[i].NumNgb = PPP[i].DhsmlNgbFactor = P[i].Particle_DivVel = 0;
                }
#if defined(ADAPTIVE_GRAVSOFT_FORALL) /* if particle is AGS-active and non-gas, set DivVel to zero because it will be reset in ags_hsml routine */
                if(ags_density_isactive(i) && (P[i].Type > 0)) {PPP[i].Particle_DivVel = 0;}
#endif

                // inverse of SPH volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(PPP[i].DhsmlNgbFactor > -0.9) {PPP[i].DhsmlNgbFactor = 1 / (1 + PPP[i].DhsmlNgbFactor);} else {PPP[i].DhsmlNgbFactor = 1;} /* note: this would be -1 if only a single particle at zero lag is found */
                P[i].Particle_DivVel *= PPP[i].DhsmlNgbFactor;

                double dimless_face_leak=0; MyLongDouble NV_T_prev[6]; NV_T_prev[0]=SphP[i].NV_T[0][0]; NV_T_prev[1]=SphP[i].NV_T[1][1]; NV_T_prev[2]=SphP[i].NV_T[2][2]; NV_T_prev[3]=SphP[i].NV_T[0][1]; NV_T_prev[4]=SphP[i].NV_T[0][2]; NV_T_prev[5]=SphP[i].NV_T[1][2];
                if(P[i].Type == 0) /* invert the NV_T matrix we just measured */
                {
                    /* use the single-moment terms of NV_T to construct the faces one would have if the system were perfectly symmetric in reconstruction 'from both sides' */
                    double V_i = NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) / PPP[i].NumNgb, dx_i = pow(V_i , 1./NUMDIMS); // this is the effective volume which will be used below
                    dx_i = sqrt(V_i * (SphP[i].NV_T[0][0] + SphP[i].NV_T[1][1] + SphP[i].NV_T[2][2])); // this is the sqrt of the weighted sum of (w*r^2)
                    double Face_Area_OneSided_Estimator_in[3]={0}, Face_Area_OneSided_Estimator_out[3]={0}; Face_Area_OneSided_Estimator_in[0]=SphP[i].NV_T[1][0]; Face_Area_OneSided_Estimator_in[1]=SphP[i].NV_T[2][0]; Face_Area_OneSided_Estimator_in[2]=SphP[i].NV_T[2][1];
                    /* now fill in the missing elements of NV_T (it's symmetric, so we saved time not computing these directly) */
                    SphP[i].NV_T[1][0]=SphP[i].NV_T[0][1]; SphP[i].NV_T[2][0]=SphP[i].NV_T[0][2]; SphP[i].NV_T[2][1]=SphP[i].NV_T[1][2];
                    double dimensional_NV_T_normalizer = pow( PPP[i].Hsml , 2-NUMDIMS ); /* this has the same dimensions as NV_T here */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {SphP[i].NV_T[k1][k2] /= dimensional_NV_T_normalizer;}} /* now NV_T should be dimensionless */
                    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
                        this will tell us how robust our procedure is (and let us know if we need to expand the neighbor number */
                    double ConditionNumber_threshold = 10. * CONDITION_NUMBER_DANGER; /* set a threshold condition number - above this we will 'pre-condition' the matrix for better behavior */
                    double trace_initial = SphP[i].NV_T[0][0] + SphP[i].NV_T[1][1] + SphP[i].NV_T[2][2]; /* initial trace of this symmetric, positive-definite matrix; used below as a characteristic value for adding the identity */
                    double conditioning_term_to_add = 1.05 * (trace_initial / NUMDIMS) / ConditionNumber_threshold; /* this will be added as a test value if the code does not reach the desired condition number */
                    while(1)
                    {
                        ConditionNumber = matrix_invert_ndims(SphP[i].NV_T, Tinv);
                        if(ConditionNumber < ConditionNumber_threshold) {break;}
                        for(k1=0;k1<NUMDIMS;k1++) {SphP[i].NV_T[k1][k1] += conditioning_term_to_add;} /* add the conditioning term which should make the matrix better-conditioned for subsequent use */
                        conditioning_term_to_add *= 1.2; /* multiply the conditioning term so it will grow and eventually satisfy our criteria */
                    }
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {SphP[i].NV_T[k1][k2] = Tinv[k1][k2] / dimensional_NV_T_normalizer;}} /* re-insert normalization correctly */
                    /* now NV_T holds the inverted matrix elements, for use in hydro */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Face_Area_OneSided_Estimator_out[k1] += 2.*V_i*SphP[i].NV_T[k1][k2]*Face_Area_OneSided_Estimator_in[k2];}} /* calculate mfm/mfv areas that we would have by default, if both sides of reconstruction were symmetric */
                    for(k1=0;k1<3;k1++) {dimless_face_leak += fabs(Face_Area_OneSided_Estimator_out[k1]) / NUMDIMS;} // average of absolute values
#ifdef HYDRO_KERNEL_SURFACE_VOLCORR
                    double closure_asymm=0; for(k1=0;k1<3;k1++) {closure_asymm += Face_Area_OneSided_Estimator_in[k1]*Face_Area_OneSided_Estimator_in[k1];}
                    double particle_inverse_volume = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
                    closure_asymm = sqrt(closure_asymm) / (PPP[i].Hsml * particle_inverse_volume); // dimensionnless measure of asymmetry in kernel
                    SphP[i].FaceClosureError = DMIN(DMAX(1.0259-2.52444*closure_asymm,0.344301),1.); // correction factor for 'missing' volume assuming a wendland C2 kernel and a sharp surface from Reinhardt & Stadel 2017 (arXiv:1701.08296)
#else
                    SphP[i].FaceClosureError = dimless_face_leak / (2.*NUMDIMS*pow(dx_i,NUMDIMS-1));
#endif
                } // P[i].Type == 0 //

                /* now check whether we had enough neighbours */
                double ncorr_ngb = 1.0;
                double cn=1;
                double c0 = 0.1 * (double)CONDITION_NUMBER_DANGER;
                if(P[i].Type==0)
                {
                    /* use the previous timestep condition number to correct how many neighbors we should use for stability */
                    if((iter==0)&&(ConditionNumber>SphP[i].ConditionNumber))
                    {
                        /* if we find ourselves with a sudden increase in condition number - check if we have a reasonable
                            neighbor number for the previous iteration, and if so, use the new (larger) correction */
                        ncorr_ngb=1; cn=SphP[i].ConditionNumber; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
                        double dn_ngb = fabs(PPP[i].NumNgb-All.DesNumNgb*ncorr_ngb)/(desnumngbdev_0*ncorr_ngb);
                        ncorr_ngb=1; cn=ConditionNumber; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
                        double dn_ngb_alt = fabs(PPP[i].NumNgb-All.DesNumNgb*ncorr_ngb)/(desnumngbdev_0*ncorr_ngb);
                        dn_ngb = DMIN(dn_ngb,dn_ngb_alt);
                        if(dn_ngb < 10.0) SphP[i].ConditionNumber = ConditionNumber;
                    }
                    ncorr_ngb=1; cn=SphP[i].ConditionNumber; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
#if !defined(HYDRO_KERNEL_SURFACE_VOLCORR)
                    double d00=0.35; if(SphP[i].FaceClosureError > d00) {ncorr_ngb = DMAX(ncorr_ngb , DMIN(SphP[i].FaceClosureError/d00 , 2.));}
#endif
                }
                desnumngb = All.DesNumNgb * ncorr_ngb;
                desnumngbdev = desnumngbdev_0 * ncorr_ngb;
                /* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
#if !defined(EOS_ELASTIC)
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}
#endif

#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                {
                    desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;
#ifdef SINGLE_STAR_SINK_DYNAMICS
                    desnumngbdev = (All.BlackHoleNgbFactor+1);
#else
                    desnumngbdev = 4 * (All.BlackHoleNgbFactor+1);
#endif
                }
#endif

#ifdef GRAIN_FLUID /* for the grains, we only need to estimate neighboring gas properties, we don't need to worry about condition numbers or conserving an exact neighbor number */
                if((1 << P[i].Type) & (GRAIN_PTYPES))
                {
                    desnumngb = All.DesNumNgb; desnumngbdev = All.DesNumNgb / 4;
#if defined(GRAIN_BACKREACTION)
                    desnumngbdev = desnumngbdev_0;
#endif
                }
#endif

                double minsoft = All.MinHsml;
                double maxsoft = All.MaxHsml;

#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
                /* use a much looser check for N_neighbors when the central point is a star particle,
                 since the accuracy is limited anyways to the coupling efficiency -- the routines use their
                 own estimators+neighbor loops, anyways, so this is just to get some nearby particles */
                if( ((1 << P[i].Type) & (valid_stellar_types)) && !((1 << P[i].Type) & (invalid_stellar_types)) )
                {
                    desnumngb = All.DesNumNgb;
#if defined(RT_SOURCE_INJECTION)
                    if(desnumngb < 64.0) {desnumngb = 64.0;} // we do want a decent number to ensure the area around the particle is 'covered'
#endif
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                    desnumngb = DMAX(desnumngb , 128); // we do want a decent number to ensure the area around the particle is 'covered'
                    if(KERNEL_FUNCTION > 3) {desnumngb = DMAX(desnumngb, 256);}
#endif
#ifdef GALSF
                    if(desnumngb < 64.0) {desnumngb = 64.0;} // we do want a decent number to ensure the area around the particle is 'covered'
                    // if we're finding this for feedback routines, there isn't any good reason to search beyond a modest physical radius //
                    double unitlength_in_kpc=UNIT_LENGTH_IN_KPC*All.cf_atime;
                    maxsoft = 2.0 / unitlength_in_kpc;
#endif
                    desnumngbdev = desnumngb / 2; // enforcing exact number not important
                }
#endif

#ifdef BLACK_HOLES
                if(P[i].Type == 5) {maxsoft = All.BlackHoleMaxAccretionRadius / All.cf_atime;}  // MaxAccretionRadius is now defined in params.txt in PHYSICAL units
#ifdef SINGLE_STAR_SINK_DYNAMICS
		        if(P[i].Type == 5) {minsoft = All.ForceSoftening[5] / All.cf_atime;} // we should always find all neighbours within the softening kernel/accretion radius, which is a lower bound on the accretion radius
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                if(P[i].Type == 5) {minsoft = DMAX(minsoft, P[i].SinkRadius);}
#endif
#endif
#endif

                redo_particle = 0;

                /* check if we are in the 'normal' range between the max/min allowed values */
                if((PPP[i].NumNgb < (desnumngb - desnumngbdev) && PPP[i].Hsml < 0.999*maxsoft) ||
                   (PPP[i].NumNgb > (desnumngb + desnumngbdev) && PPP[i].Hsml > 1.001*minsoft))
                    {redo_particle = 1;}

                /* check maximum kernel size allowed */
                particle_set_to_maxhsml_flag = 0;
                if((PPP[i].Hsml >= 0.999*maxsoft) && (PPP[i].NumNgb < (desnumngb - desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].Hsml == maxsoft)
                    {
                        /* iteration at the maximum value is already complete */
                        particle_set_to_maxhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        PPP[i].Hsml = maxsoft;
                        particle_set_to_maxhsml_flag = 1;
                    }
                }

                /* check minimum kernel size allowed */
                particle_set_to_minhsml_flag = 0;
                if((PPP[i].Hsml <= 1.001*minsoft) && (PPP[i].NumNgb > (desnumngb + desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].Hsml == minsoft)
                    {
                        /* this means we've already done an iteration with the MinHsml value, so the
                         neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
                        particle_set_to_minhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        PPP[i].Hsml = minsoft;
                        particle_set_to_minhsml_flag = 1;
                    }
                }

#ifdef GALSF
                if((All.ComovingIntegrationOn)&&(All.Time>All.TimeBegin))
                {
                    if((P[i].Type==4)&&(iter>1)&&(PPP[i].NumNgb>4)&&(PPP[i].NumNgb<100)&&(redo_particle==1)) {redo_particle=0;}
                }
#endif

                if((redo_particle==0)&&(P[i].Type == 0))
                {
                    /* ok we have reached the desired number of neighbors: save the condition number for next timestep */
                    if(ConditionNumber > 1e6 * (double)CONDITION_NUMBER_DANGER) {
                        PRINT_WARNING("Condition number=%g CNum_prevtimestep=%g CNum_danger=%g iter=%d Num_Ngb=%g desnumngb=%g Hsml=%g Hsml_min=%g Hsml_max=%g \n i=%d task=%d ID=%llu Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)  \n NVT=%.17g/%.17g/%.17g %.17g/%.17g/%.17g %.17g/%.17g/%.17g NVT_inv=%.17g/%.17g/%.17g %.17g/%.17g/%.17g %.17g/%.17g/%.17g ",
                               ConditionNumber,SphP[i].ConditionNumber,CONDITION_NUMBER_DANGER,iter,PPP[i].NumNgb,desnumngb,PPP[i].Hsml,All.MinHsml,All.MaxHsml, i, ThisTask,
                               (unsigned long long) P[i].ID, P[i].Type, PPP[i].Hsml, PPP[i].DhsmlNgbFactor, Left[i], Right[i],
                               (float) PPP[i].NumNgb, Right[i] - Left[i], particle_set_to_maxhsml_flag, particle_set_to_minhsml_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
                               SphP[i].NV_T[0][0],SphP[i].NV_T[0][1],SphP[i].NV_T[0][2],SphP[i].NV_T[1][0],SphP[i].NV_T[1][1],SphP[i].NV_T[1][2],SphP[i].NV_T[2][0],SphP[i].NV_T[2][1],SphP[i].NV_T[2][2],
                               NV_T_prev[0],NV_T_prev[3],NV_T_prev[4],NV_T_prev[3],NV_T_prev[1],NV_T_prev[5],NV_T_prev[4],NV_T_prev[5],NV_T_prev[2]);}
                    SphP[i].ConditionNumber = ConditionNumber;
                }

                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
                        PRINT_WARNING("i=%d task=%d ID=%llu iter=%d Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)",
                               i, ThisTask, (unsigned long long) P[i].ID, iter, P[i].Type, PPP[i].Hsml, PPP[i].DhsmlNgbFactor, Left[i], Right[i],
                               (float) PPP[i].NumNgb, Right[i] - Left[i], particle_set_to_maxhsml_flag, particle_set_to_minhsml_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                    }

                    /* need to redo this particle */
                    npleft++;

                    if(Left[i] > 0 && Right[i] > 0)
                        if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                            SphP[i].ConditionNumber = ConditionNumber;
                            continue;
                        }

                    if((particle_set_to_maxhsml_flag==0)&&(particle_set_to_minhsml_flag==0))
                    {
                        if(PPP[i].NumNgb < (desnumngb - desnumngbdev)) {Left[i] = DMAX(PPP[i].Hsml, Left[i]);}
                        else
                        {
                            if(Right[i] != 0) {if(PPP[i].Hsml < Right[i]) {Right[i] = PPP[i].Hsml;}} else {Right[i] = PPP[i].Hsml;}
                        }

                        // right/left define upper/lower bounds from previous iterations
                        if(Right[i] > 0 && Left[i] > 0)
                        {
                            // geometric interpolation between right/left //
                            double maxjump=0;
                            if(iter>1) {maxjump = 0.2*log(Right[i]/Left[i]);}
                            if(PPP[i].NumNgb > 1)
                            {
                                double jumpvar = PPP[i].DhsmlNgbFactor * log( desnumngb / PPP[i].NumNgb ) / NUMDIMS;
                                if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
                                PPP[i].Hsml *= exp(jumpvar);
                            } else {
                                PPP[i].Hsml *= 2.0;
                            }
                            if((PPP[i].Hsml<Right[i])&&(PPP[i].Hsml>Left[i]))
                            {
                                if(iter > 1)
                                {
                                    double hfac = exp(maxjump);
                                    if(PPP[i].Hsml > Right[i] / hfac) {PPP[i].Hsml = Right[i] / hfac;}
                                    if(PPP[i].Hsml < Left[i] * hfac) {PPP[i].Hsml = Left[i] * hfac;}
                                }
                            } else {
                                if(PPP[i].Hsml>Right[i]) PPP[i].Hsml=Right[i];
                                if(PPP[i].Hsml<Left[i]) PPP[i].Hsml=Left[i];
                                PPP[i].Hsml = pow(PPP[i].Hsml * Left[i] * Right[i] , 1.0/3.0);
                            }
                        }
                        else
                        {
                            if(Right[i] == 0 && Left[i] == 0)
                            {
                                char buf[1000]; sprintf(buf, "Right[i] == 0 && Left[i] == 0 && PPP[i].Hsml=%g\n", PPP[i].Hsml); terminate(buf);
                            }

                            if(Right[i] == 0 && Left[i] > 0)
                            {
                                if (PPP[i].NumNgb > 1)
                                    {fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS;} // this would give desnumgb if constant density (+0.231=2x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(PPP[i].DhsmlNgbFactor==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac < fac_lim+0.231)
                                    {
                                        PPP[i].Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                    {
                                        PPP[i].Hsml *= exp(fac_lim+0.231);
                                        // fac~0.26 leads to expected doubling of number if density is constant,
                                        //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                    }
                                }
                                else
                                    {PPP[i].Hsml *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }

                            if(Right[i] > 0 && Left[i] == 0)
                            {
                                if(PPP[i].NumNgb > 1)
                                    {fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS;} // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if(fac_lim < -1.535) {fac_lim = -1.535;} // decreasing N_ngb by factor ~100

                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(PPP[i].DhsmlNgbFactor==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac > fac_lim-0.231)
                                    {
                                        PPP[i].Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                        {PPP[i].Hsml *= exp(fac_lim-0.231);} // limiter to prevent --too-- far a jump in a single iteration
                                }
                                else
                                    {PPP[i].Hsml *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                        } // closes if[particle_set_to_max/minhsml_flag]
                    } // closes redo_particle
                    /* resets for max/min values */
                    if(PPP[i].Hsml < minsoft) {PPP[i].Hsml = minsoft;}
                    if(particle_set_to_minhsml_flag==1) {PPP[i].Hsml = minsoft;}
                    if(PPP[i].Hsml > maxsoft) {PPP[i].Hsml = maxsoft;}
                    if(particle_set_to_maxhsml_flag==1) {PPP[i].Hsml = maxsoft;}
                }
                else {P[i].TimeBin = -P[i].TimeBin - 1;}	/* Mark as inactive */
            } //  if(density_isactive(i))
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])

        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 10) {PRINT_STATUS("ngb iteration %d: need to repeat for %d%09d particles", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("failed to converge in neighbour iteration in density()\n"); fflush(stdout); endrun(1155);}
        }
    }
    while(ntot > 0);

    /* iteration is done - de-malloc everything now */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    myfree(Right); myfree(Left);

    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0) {P[i].TimeBin = -P[i].TimeBin - 1;}
    }


    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results
     ( any quantities that only need to be evaluated once, on the final iteration --
     won't save much b/c the real cost is in the neighbor loop for each particle, but it's something )
     -- also, some results (for example, viscosity suppression below) should not be calculated unless
     the quantities are 'stabilized' at their final values -- */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(density_isactive(i))
        {
            if(P[i].Type == 0 && P[i].Mass > 0)
            {
                if(SphP[i].Density > 0)
                {
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
                    /* set motion of the mesh-generating points */
#if (HYDRO_FIX_MESH_MOTION==4)
                    set_mesh_motion(i); // use user-specified analytic function to define mesh motions //
#elif ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                    double eps_pvel = 0.3; // normalization for how much 'weight' to give to neighbors (unstable if >=0.5)
                    for(k=0;k<3;k++) {SphP[i].ParticleVel[k] = (1.-eps_pvel)*SphP[i].VelPred[k] + eps_pvel*SphP[i].ParticleVel[k]/SphP[i].Density;} // assign mixture velocity
#elif (HYDRO_FIX_MESH_MOTION==7)
                    for(k=0;k<3;k++) {SphP[i].ParticleVel[k] = SphP[i].VelPred[k];} // move with fluid
#endif
#endif

#ifdef HYDRO_SPH
#ifdef HYDRO_PRESSURE_SPH
                    if(SphP[i].InternalEnergyPred > 0)
                    {
                        SphP[i].EgyWtDensity /= SphP[i].InternalEnergyPred;
                    } else {
                        SphP[i].EgyWtDensity = 0;
                    }
#endif
                    /* need to divide by the sum of x_tilde=1, i.e. numden_ngb */
                    if((PPP[i].Hsml > 0)&&(PPP[i].NumNgb > 0))
                    {
                        double numden_ngb = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
                        SphP[i].DhsmlHydroSumFactor *= PPP[i].Hsml / (NUMDIMS * numden_ngb);
                        SphP[i].DhsmlHydroSumFactor *= -PPP[i].DhsmlNgbFactor; /* now this is ready to be called in hydro routine */
                    } else {
                        SphP[i].DhsmlHydroSumFactor = 0;
                    }
#endif


#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
                    for(k1 = 0; k1 < 3; k1++)
                        for(k2 = 0; k2 < 3; k2++)
                        {
                            SphP[i].NV_D[k2][k1] *= All.cf_a2inv; // converts to physical velocity/length
                            SphP[i].NV_A[k2][k1] /= All.cf_atime; // converts to physical accel/length
                        }
                    // all quantities below in this block should now be in proper PHYSICAL units, for subsequent operations //
                    double dtDV[3][3], A[3][3], V[3][3], S[3][3];
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            V[k1][k2] = SphP[i].NV_D[k1][0]*SphP[i].NV_T[0][k2] + SphP[i].NV_D[k1][1]*SphP[i].NV_T[1][k2] + SphP[i].NV_D[k1][2]*SphP[i].NV_T[2][k2];
                            A[k1][k2] = SphP[i].NV_A[k1][0]*SphP[i].NV_T[0][k2] + SphP[i].NV_A[k1][1]*SphP[i].NV_T[1][k2] + SphP[i].NV_A[k1][2]*SphP[i].NV_T[2][k2];
                        }
                    SphP[i].NV_DivVel = V[0][0] + V[1][1] + V[2][2];
                    SphP[i].NV_trSSt = 0;
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            dtDV[k1][k2] = A[k1][k2] - (V[k1][0]*V[0][k2] + V[k1][1]*V[1][k2] + V[k1][2]*V[2][k2]);
                            /* S = 0.5*(V+V_transpose) - delta_ij*div_v/3 */
                            S[k1][k2] = 0.5 * (V[k1][k2] + V[k2][k1]);
                            if(k2==k1) S[k1][k2] -= SphP[i].NV_DivVel / NUMDIMS;
                            /* Trace[S*S_transpose] = SSt[0][0]+SSt[1][1]+SSt[2][2] = |S|^2 = sum(Sij^2) */
                            SphP[i].NV_trSSt += S[k1][k2]*S[k1][k2];
                        }
                    SphP[i].NV_dt_DivVel = dtDV[0][0] + dtDV[1][1] + dtDV[2][2];
#endif


#if defined(TURB_DRIVING)
                    if(SphP[i].Density > 0)
                    {
                        SphP[i].SmoothedVel[0] /= SphP[i].Density;
                        SphP[i].SmoothedVel[1] /= SphP[i].Density;
                        SphP[i].SmoothedVel[2] /= SphP[i].Density;
                    } else {
                        SphP[i].SmoothedVel[0] = SphP[i].SmoothedVel[1] = SphP[i].SmoothedVel[2] = 0;
                    }
#endif
                }

#ifndef HYDRO_SPH
                if((PPP[i].Hsml > 0)&&(PPP[i].NumNgb > 0))
                {
                    SphP[i].Density = P[i].Mass * PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) ); // divide mass by volume
                } else {
                    if(PPP[i].Hsml <= 0)
                    {
                        SphP[i].Density = 0; // in this case, give up, no meaningful volume
                    } else {
                        SphP[i].Density = P[i].Mass / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) ); // divide mass (lone particle) by volume
                    }
                }
#endif
                double Volume_0; Volume_0 = P[i].Mass / SphP[i].Density; // save for potential later use
#if defined(HYDRO_KERNEL_SURFACE_VOLCORR)
                SphP[i].Density /= SphP[i].FaceClosureError; // correct volume of the cell based on the free surface correction above
                SphP[i].FaceClosureError = Volume_0;
#endif
#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
                Volume_0 = P[i].Mass / SphP[i].Density;
                if(All.Time == All.TimeBegin) {SphP[i].Density_ExplicitInt = SphP[i].Density;} // set initial value to density calculated above
                    else {SphP[i].Density = SphP[i].Density_ExplicitInt;} // set to explicitly-evolved density field
                SphP[i].FaceClosureError = Volume_0;
#endif
#ifdef HYDRO_VOLUME_CORRECTIONS
                SphP[i].Volume_1 = SphP[i].Volume_0 = Volume_0; // initialize this value for use in the correction loop, and in case this is not set in the subsequent loop because of inactivity, set this first to the zeroth-order estimator
#endif
                SphP[i].Pressure = get_pressure(i);		// should account for density independent pressure

            } // P[i].Type == 0


#if defined(GRAIN_FLUID)
            if((1 << P[i].Type) & (GRAIN_PTYPES))
            {
                if(P[i].Gas_Density > 0)
                {
                    P[i].Gas_InternalEnergy /= P[i].Gas_Density;
                    for(k = 0; k<3; k++) {P[i].Gas_Velocity[k] /= P[i].Gas_Density;}
                } else {
                    P[i].Gas_InternalEnergy = 0;
                    for(k = 0; k<3; k++) {P[i].Gas_Velocity[k] = 0;}
#if defined(GRAIN_LORENTZFORCE)
                    for(k = 0; k<3; k++) {P[i].Gas_B[k] = 0;}
#endif
                }
            }
#endif


#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
            /* non-gas particles are handled separately, in the ags_hsml routine */
            if(P[i].Type==0)
            {
                PPPZ[i].AGS_zeta = 0;
                double zeta_0 = 0; // 2.0 * P[i].Mass*P[i].Mass * PPP[i].Hsml*PPP[i].Hsml; // self-value of zeta if no neighbors are found //
                if((PPP[i].Hsml > 0)&&(PPP[i].NumNgb > 0))
                {
                    /* the zeta terms ONLY control errors if we maintain the 'correct' neighbor number: for boundary
                        particles, it can actually be worse. so we need to check whether we should use it or not */
                    if((PPP[i].Hsml > 1.001*All.MinHsml) && (PPP[i].Hsml < 0.999*All.MaxHsml) &&
                        (fabs(PPP[i].NumNgb-All.DesNumNgb)/All.DesNumNgb < 0.05))
                    {
                        double ndenNGB = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
                        PPPZ[i].AGS_zeta *= 0.5 * P[i].Mass * PPP[i].Hsml / (NUMDIMS * ndenNGB) * PPP[i].DhsmlNgbFactor;
                    } else {
                        PPPZ[i].AGS_zeta = zeta_0;
                    }
                } else {
                    PPPZ[i].AGS_zeta = zeta_0;
                }
            }
#endif
            apply_pm_hires_region_clipping_selection(i);

         /* finally, convert NGB to the more useful format, NumNgb^(1/NDIMS),
            which we can use to obtain the corrected particle sizes. Because of how this number is used above, we --must-- make
            sure that this operation is the last in the loop here */
            if(PPP[i].NumNgb > 0) {PPP[i].NumNgb=pow(PPP[i].NumNgb,1./NUMDIMS);} else {PPP[i].NumNgb=0;}

        } // density_isactive(i)

#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) /* re-assign magnetic fields after getting the correct density for newly-spawned cells when these options are enabled */
        if(P[i].Type==0) {if(P[i].ID==All.AGNWindID && SphP[i].IniDen<0) {SphP[i].IniDen=SphP[i].Density; int k; for(k=0;k<3;k++) {SphP[i].BPred[k]=SphP[i].B[k]=SphP[i].IniB[k]*(All.UnitMagneticField_in_gauss/UNIT_B_IN_GAUSS)*(P[i].Mass/(All.cf_a2inv*SphP[i].Density));}}}
#endif

    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])

    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_DENSCOMPUTE] += timecomp; CPU_Step[CPU_DENSWAIT] += timewait;
    CPU_Step[CPU_DENSCOMM] += timecomm; CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */



/* Routines for a loop after the iterative density loop needed to find neighbors, etc, once all have converged, to apply additional correction terms to the cell volumes and faces (for those needed -before- the gradients loop because they alter primitive quantities needed for gradients, such as particle densities, pressures, etc.)
    This was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */
#ifdef HYDRO_VOLUME_CORRECTIONS

#define CORE_FUNCTION_NAME cellcorrections_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_cellcorrections    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_cellcorrections  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(GasGrad_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. for current implementation, only cells eligible for gradients and hydro should be called */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/* define structures to use below */
struct INPUT_STRUCT_NAME {MyDouble Pos[3], Hsml, Volume_0; int NodeList[NODELISTLENGTH];} *DATAIN_NAME, *DATAGET_NAME;

/* define properties to be sent to nodes */
void particle2in_cellcorrections(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{in->Volume_0=SphP[i].Volume_0; in->Hsml=PPP[i].Hsml; int k; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];}}

/* define output structure to use below */
struct OUTPUT_STRUCT_NAME {MyFloat Volume_1;} *DATARESULT_NAME, *DATAOUT_NAME;

/* define properties to be collected from nodes */
void out2particle_cellcorrections(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{ASSIGN_ADD(SphP[i].Volume_1, out->Volume_1, mode);}

/* core subroutine. this does not write to shared memory. */
int cellcorrections_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, n, k, startnode, numngb_inbox, listindex = 0; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* generic variables we always use, and set initial memory */
    if(mode == 0) {particle2in_cellcorrections(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode; /* open it */} /* start usual neighbor tree search */
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}
            for(n=0; n<numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                double dp[3]; for(k=0;k<3;k++) {dp[k]=local.Pos[k]-P[j].Pos[k];}
                NEAREST_XYZ(dp[0],dp[1],dp[2],1); // find the closest image in the given box size  //
                double r2=0; for(k=0;k<3;k++) {r2+=dp[k]*dp[k];} // distance
                if(r2 >= PPP[j].Hsml*PPP[j].Hsml) {continue;} // need to be inside of 'j's kernel search
                double u,hinv,hinv3,hinv4,wk,dwk; kernel_hinv(PPP[j].Hsml, &hinv, &hinv3, &hinv4); u=sqrt(r2)*hinv; wk=0; dwk=0; // define kernel-needed variables
                kernel_main(u, hinv3, hinv4, &wk, &dwk, -1); // calculate the normal kernel weight 'wk'
                out.Volume_1 += SphP[j].Volume_0 * SphP[j].Volume_0 * wk; // this is the next-order correction to the cell volume quadrature
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode=DATAGET_NAME[target].NodeList[listindex]; if(startnode>=0) {startnode=Nodes[startnode].u.d.nextnode;}}} // handle opening nodes
    } // closes while(startnode >= 0)
    if(mode == 0) {out2particle_cellcorrections(&out, target, 0, 0);} else {DATARESULT_NAME[target] = out;} /* collect the result at the right place */
    return 0; /* done */
}

/* final operations for after the updates are computed */
void cellcorrections_final_operations_and_cleanup(void)
{
    int i; for(i=FirstActiveParticle; i>=0; i=NextActiveParticle[i]) { /* check all active elements */
        CONDITIONFUNCTION_FOR_EVALUATION /* ensures only the ones which met our criteria above are actually treated here */
        {
            if(SphP[i].Volume_1 > 0) {SphP[i].Density = P[i].Mass / SphP[i].Volume_1;} else {SphP[i].Volume_1 = SphP[i].Volume_0;} // set the updated density. other variables that need volumes will all scale off this, so we can rely on it to inform everything else [if bad value here, revert to the 0th-order volume quadrature]
            SphP[i].Pressure = get_pressure(i);
        }}
}

/* parent routine which calls the work loop above */
void cellcorrections_calc(void)
{
    CPU_Step[CPU_DENSMISC] += measure_time(); double t00_truestart = my_second();
    PRINT_STATUS(" ..calculating first-order corrections to cell sizes/faces");
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    cellcorrections_final_operations_and_cleanup(); /* do final operations on results */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_DENSCOMPUTE] += timecomp; CPU_Step[CPU_DENSWAIT] += timewait; CPU_Step[CPU_DENSCOMM] += timecomm;
    CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

#endif // parent if statement for all code in the HYDRO_VOLUME_CORRECTIONS block
