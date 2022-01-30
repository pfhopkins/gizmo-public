/*! \file blackhole_swallow_and_kick.c
*  \brief routines for gas accretion onto black holes, and black hole mergers
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
/*
* This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
* see notes in blackhole.c for details on code history.
*/


#ifdef BLACK_HOLES // top-level flag [needs to be here to prevent compiler breaking when this is not active] //


static int N_gas_swallowed, N_star_swallowed, N_dm_swallowed, N_BH_swallowed;

#ifdef BH_ALPHADISK_ACCRETION
#define out_accreted_BH_Mass_alphaornot out.accreted_BH_Mass_alphadisk
#else
#define out_accreted_BH_Mass_alphaornot out.accreted_BH_Mass
#endif


#define CORE_FUNCTION_NAME blackhole_swallow_and_kick_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(bhsink_isactive(i) && P[i].SwallowID==0) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */


/* this structure defines the variables that need to be sent -from- the 'searching' element */
struct INPUT_STRUCT_NAME
{
    int NodeList[NODELISTLENGTH]; MyDouble Pos[3]; MyFloat Vel[3], Hsml, Mass, BH_Mass, Dt, Mdot; MyIDType ID;
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS) || defined(BH_WIND_KICK)
    MyFloat Jgas_in_Kernel[3];
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat BH_Mass_AlphaDisk;
#endif
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    MyFloat BH_angle_weighted_kernel_sum;
#endif
#if defined(BH_RETURN_ANGMOM_TO_GAS)
    MyFloat BH_Specific_AngMom[3], angmom_norm_topass_in_swallowloop;
#endif
#if defined(BH_RETURN_BFLUX)
    MyFloat B[3];
    MyFloat kernel_norm_topass_in_swallowloop;
#endif
#ifdef SINGLE_STAR_FB_LOCAL_RP
    MyFloat Luminosity;
#endif    
}
*DATAIN_NAME, *DATAGET_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k, j_tempinfo; j_tempinfo = P[i].IndexMapToTempStruc; /* link to the location in the shared structure where this is stored */
    for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k]; in->Vel[k]=P[i].Vel[k];} /* good example - always needed */
    in->Hsml = PPP[i].Hsml; in->Mass = P[i].Mass; in->BH_Mass = BPP(i).BH_Mass; in->ID = P[i].ID; in->Mdot = BPP(i).BH_Mdot;
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS) || defined(BH_WIND_KICK)
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
    for(k=0;k<3;k++) {in->Jgas_in_Kernel[k] = P[i].BH_Specific_AngMom[k];}
#else
    for(k=0;k<3;k++) {in->Jgas_in_Kernel[k] = BlackholeTempInfo[j_tempinfo].Jgas_in_Kernel[k];}
#endif
#endif
#ifdef BH_ALPHADISK_ACCRETION
    in->BH_Mass_AlphaDisk = BPP(i).BH_Mass_AlphaDisk;
#endif
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    in->BH_angle_weighted_kernel_sum = BlackholeTempInfo[j_tempinfo].BH_angle_weighted_kernel_sum;
#endif
    in->Dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
    if(P[i].Type == 5){in->Dt = P[i].dt_since_last_gas_search;}
#endif
#if defined(BH_RETURN_ANGMOM_TO_GAS)
    for(k=0;k<3;k++) {in->BH_Specific_AngMom[k] = BPP(i).BH_Specific_AngMom[k];}
    in->angmom_norm_topass_in_swallowloop = BlackholeTempInfo[j_tempinfo].angmom_norm_topass_in_swallowloop;
#endif
#if defined(BH_RETURN_BFLUX)
    for(k=0;k<3;k++) {in->B[k] = BPP(i).B[k];}
    in->kernel_norm_topass_in_swallowloop = BlackholeTempInfo[j_tempinfo].kernel_norm_topass_in_swallowloop;
#endif
#ifdef SINGLE_STAR_FB_LOCAL_RP
    in->Luminosity = bh_lum_bol(in->Mdot, in->BH_Mass, i);
#endif    
}


/* this structure defines the variables that need to be sent -back to- the 'searching' element */
struct OUTPUT_STRUCT_NAME
{ /* define variables below as e.g. "double X;" */
    MyDouble accreted_Mass, accreted_BH_Mass, accreted_BH_Mass_alphadisk;
#ifdef GRAIN_FLUID
    MyDouble accreted_dust_Mass;
#endif
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
    MyDouble accreted_momentum[3];
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
    MyDouble accreted_centerofmass[3];
#endif
#if defined(BH_RETURN_BFLUX)
    MyDouble accreted_B[3];
//    MyDouble accreted_Phi;
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
    MyDouble accreted_J[3];
#endif
#ifdef BH_COUNTPROGS
    int BH_CountProgs;
#endif
#ifdef GALSF
    MyFloat Accreted_Age;
#endif
}
*DATARESULT_NAME, *DATAOUT_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

#define ASSIGN_ADD_PRESET(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    int k, target = P[i].IndexMapToTempStruc; k=0;
    ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_Mass, out->accreted_Mass, mode);
    ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_BH_Mass, out->accreted_BH_Mass, mode);
    ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_BH_Mass_alphadisk, out->accreted_BH_Mass_alphadisk, mode);
#ifdef GRAIN_FLUID
    ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_dust_Mass, out->accreted_dust_Mass, mode);
#endif
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
    for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_momentum[k], out->accreted_momentum[k], mode);}
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
    for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_centerofmass[k], out->accreted_centerofmass[k], mode);}
#endif
#if defined(BH_RETURN_BFLUX)
    for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_B[k], out->accreted_B[k], mode);}
//    ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_Phi, out->accreted_Phi, mode);
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
    for(k=0;k<3;k++) {ASSIGN_ADD_PRESET(BlackholeTempInfo[target].accreted_J[k], out->accreted_J[k], mode);}
#endif
#ifdef BH_COUNTPROGS
    BPP(i).BH_CountProgs += out->BH_CountProgs;
#endif
#ifdef GALSF
    if(P[i].StellarAge > out->Accreted_Age) P[i].StellarAge = out->Accreted_Age;
#endif
}


/* subroutine that handles the 'final' swallowing and kicking operations for coupling feedback or accretion between sinks and neighboring elements */
/*!   -- this subroutine writes to shared memory [updating the neighbor values]: need to protect these writes for openmp below. there's a lot of writing, so this needs to be done carefully. */
int blackhole_swallow_and_kick_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration);
int blackhole_swallow_and_kick_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb, listindex = 0, j, k, n, bin; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* define variables and zero memory and import data for local target*/
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* imports the data to the correct place and names */
    double h_i=local.Hsml, hinv=1/h_i, hinv3, f_accreted; hinv3=hinv*hinv*hinv; f_accreted=0;
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    double kernel_zero,dwk; kernel_main(0.0,1.0,1.0,&kernel_zero,&dwk,-1); dwk=0;
#ifdef SINGLE_STAR_FB_LOCAL_RP
    double mom = local.Luminosity * local.Dt / C_LIGHT_CODE, mom_wt = 0;
#else    
    double mom = bh_lum_bol(local.Mdot, local.BH_Mass, -1) * local.Dt / C_LIGHT_CODE, mom_wt = 0;
#endif
#endif
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS) || defined(BH_WIND_KICK)
    double J_dir[3]; for(k=0;k<3;k++) {J_dir[k] = local.Jgas_in_Kernel[k];}
    double norm=0; for(k=0;k<3;k++) {norm+=J_dir[k]*J_dir[k];}
    if(norm>0) {norm=1/sqrt(norm); for(k=0;k<3;k++) {J_dir[k]*=norm;}} else {J_dir[0]=J_dir[1]=0; J_dir[2]=1;}
#endif
#if defined(BH_WIND_KICK)
    double bh_mass_withdisk=local.BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += local.BH_Mass_AlphaDisk;
#endif
#endif
#ifdef GALSF
    out.Accreted_Age = MAX_REAL_NUMBER;
#endif

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb = ngb_treefind_pairs_threads_targeted(local.Pos, h_i, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, BH_NEIGHBOR_BITFLAG);
            if(numngb < 0) {return -2;}
            for(n = 0; n < numngb; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                MyIDType OriginallyMarkedSwallowID; OriginallyMarkedSwallowID = P[j].SwallowID; // record this to help prevent double-counting below
                double Mass_j, Vel_j[3], InternalEnergy_j; // velocity kicks to apply to 'j's, accumulate here to apply at bottom of code
                #pragma omp atomic read
                Mass_j = P[j].Mass; // this can get modified -a lot- below, so we need to read it carefully right now
                if(P[j].Type==0)
                {
                    #pragma omp atomic read
                    InternalEnergy_j = SphP[j].InternalEnergy; // this can get modified -a lot- below, so we need to read it carefully right now
                }
                for(k=0;k<3;k++)
                {
                    #pragma omp atomic read
                    Vel_j[k] = P[j].Vel[k]; // this can get modified -a lot- below, so we need to read it carefully right now
                }
                double Mass_j_0 = Mass_j, InternalEnergy_j_0 = InternalEnergy_j, Vel_j_0[3]; for(k=0;k<3;k++) {Vel_j_0[k]=Vel_j[k];} // save initial values to know if we need to update neighbor values below
                
                double dpos[3]={0},dvel[3]={0}; for(k=0;k<3;k++) {dpos[k]=P[j].Pos[k]-local.Pos[k]; dvel[k]=Vel_j[k]-local.Vel[k];}
                NEAREST_XYZ(dpos[0],dpos[1],dpos[2],-1); /*  find the closest image in the given box size  */
                NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,dvel,-1); /* wrap velocities for shearing boxes if needed */

                
#if defined(BH_RETURN_ANGMOM_TO_GAS) || defined(BH_RETURN_BFLUX)
                double wk, dwk, u=0;
                if(P[j].Type == 0){
                    for(k=0;k<3;k++) {u+=dpos[k]*dpos[k];}
                    u=sqrt(u)/DMAX(h_i, P[j].Hsml); if(u<1) { kernel_main(u,1., 1.,&wk,&dwk,-1); } else {wk=dwk=0;}
                }
#endif
#if defined(BH_RETURN_ANGMOM_TO_GAS) /* this should go here [right before the loop that accretes it back onto the BH] */
                if(P[j].Type == 0)
                {
                    double dlv[3]; dlv[0]=local.BH_Specific_AngMom[1]*dpos[2]-local.BH_Specific_AngMom[2]*dpos[1]; dlv[1]=local.BH_Specific_AngMom[2]*dpos[0]-local.BH_Specific_AngMom[0]*dpos[2]; dlv[2]=local.BH_Specific_AngMom[0]*dpos[1]-local.BH_Specific_AngMom[1]*dpos[0];
                    for(k=0;k<3;k++) {
                        dlv[k] *= wk * local.angmom_norm_topass_in_swallowloop;
                        Vel_j[k] += dlv[k];
                        out.accreted_momentum[k]-=Mass_j*dlv[k];
                    }
                    out.accreted_J[0]-=Mass_j*(dpos[1]*dlv[2] - dpos[2]*dlv[1]); out.accreted_J[1]-=Mass_j*(dpos[2]*dlv[0] - dpos[0]*dlv[2]); out.accreted_J[2]-=Mass_j*(dpos[0]*dlv[1] - dpos[1]*dlv[0]);
                }
#endif
#if defined(BH_RETURN_BFLUX) // do a kernel-weighted redistribution of the magnetic flux in the sink into surrounding particles
                if((P[j].Type == 0) && (local.kernel_norm_topass_in_swallowloop > 0)){
                    double dB, b_fraction_toreturn = DMIN(0.1, local.Dt / (local.BH_Mass_AlphaDisk / local.Mdot)) * wk / local.kernel_norm_topass_in_swallowloop; // return a fraction dt/t_accretion of the total flux, with simple kernel weighting for each particle
                    for(k=0; k<3;k++) {
                        dB = b_fraction_toreturn * local.B[k];
                        #pragma omp atomic
                        SphP[j].B[k] += dB;
                        #pragma omp atomic
                        SphP[j].BPred[k] +=  dB;
                        out.accreted_B[k] -= dB;
                    }
                }
#endif
                
                
                /* we've found a particle to be swallowed.  This could be a BH merger, DM particle, or baryon w/ feedback */
                if(P[j].SwallowID == local.ID && Mass_j > 0)
                {   /* accreted quantities to be added [regardless of particle type] */
                    f_accreted = 1; /* default to accreting entire particle */
#ifdef BH_WIND_KICK
                    if(P[j].Type == 0)
                    {
                        f_accreted = All.BAL_f_accretion; /* if particle is gas, only a fraction gets accreted in these particular modules */
#ifndef BH_GRAVCAPTURE_GAS
                        if((All.BlackHoleFeedbackFactor > 0) && (All.BlackHoleFeedbackFactor != 1.)) {f_accreted /= All.BlackHoleFeedbackFactor;} else {if(All.BAL_v_outflow > 0) f_accreted = 1./(1. + fabs(1.*BH_WIND_KICK)*All.BlackHoleRadiativeEfficiency*C_LIGHT_CODE/(All.BAL_v_outflow));}
                        if((bh_mass_withdisk - local.Mass) <= 0) {f_accreted=0;} // DAA: no need to accrete gas particle to enforce mass conservation (we will simply kick),  note that here the particle mass P.Mass is larger than the physical BH mass P.BH_Mass
#endif
                    }
#endif

                    
                    /* handle accretion/conservation of certain conserved quantities, depending on whether we are intending our sub-grid model to follow them */
                    double mcount_for_conserve; mcount_for_conserve = f_accreted * Mass_j;
#if (BH_FOLLOW_ACCRETED_ANGMOM == 1) /* in this case we are only counting this if its coming from BH particles */
                    if(P[j].Type != 5) {mcount_for_conserve=0;} else {mcount_for_conserve=BPP(j).BH_Mass;}
#ifdef BH_ALPHADISK_ACCRETION
                    if(P[j].Type == 5) {mcount_for_conserve += BPP(j).BH_Mass_AlphaDisk;}
#endif
#endif
#ifdef GRAIN_FLUID
                    if((1<<P[j].Type) & GRAIN_PTYPES) {out.accreted_dust_Mass += FLT(Mass_j);}
#endif
#if defined(BH_FOLLOW_ACCRETED_MOMENTUM)
                    for(k=0;k<3;k++) {out.accreted_momentum[k] += FLT( mcount_for_conserve * dvel[k]);}
#endif
#if defined(BH_FOLLOW_ACCRETED_COM)
                    for(k=0;k<3;k++) {out.accreted_centerofmass[k] += FLT(mcount_for_conserve * dpos[k]);}
#endif
#ifdef BH_RETURN_BFLUX
                    for(k=0;k<3;k++) {out.accreted_B[k] += FLT(SphP[j].BPred[k]);}
#endif
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                    out.accreted_J[0] += FLT(mcount_for_conserve * ( dpos[1]*dvel[2] - dpos[2]*dvel[1] ));
                    out.accreted_J[1] += FLT(mcount_for_conserve * ( dpos[2]*dvel[0] - dpos[0]*dvel[2] ));
                    out.accreted_J[2] += FLT(mcount_for_conserve * ( dpos[0]*dvel[1] - dpos[1]*dvel[0] ));
                    if(P[j].Type == 5) {for(k=0;k<3;k++) {out.accreted_J[k] += FLT(mcount_for_conserve * BPP(j).BH_Specific_AngMom[k]);}}
#endif

                    
                    
                    if(P[j].Type == 5)  /* this is a BH-BH merger */
                    {
#ifdef BH_OUTPUT_MOREINFO
                        fprintf(FdBhMergerDetails,"%g  %u %g %2.7f %2.7f %2.7f  %u %g %2.7f %2.7f %2.7f\n", All.Time,  local.ID,local.BH_Mass,local.Pos[0],local.Pos[1],local.Pos[2],  P[j].ID,BPP(j).BH_Mass,P[j].Pos[0],P[j].Pos[1],P[j].Pos[2]); fflush(FdBhMergerDetails);
#elif !defined(IO_REDUCED_MODE)
                        fprintf(FdBlackHolesDetails,"ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n", ThisTask, All.Time, local.ID, P[j].ID, local.BH_Mass, BPP(j).BH_Mass); fflush(FdBlackHolesDetails);
#endif
#ifdef BH_INCREASE_DYNAMIC_MASS
                        /* the true dynamical mass of the merging BH is Mass_j/BH_INCREASE_DYNAMIC_MASS unless exceeded by physical growth
                         - in the limit BPP(j).BH_Mass > BH_INCREASE_DYNAMIC_MASS x m_b, then bh_mass=Mass_j on average and we are good as well  */
                        out.accreted_Mass    += FLT( DMAX(BPP(j).BH_Mass, Mass_j/BH_INCREASE_DYNAMIC_MASS) );
#else
                        out.accreted_Mass    += FLT(Mass_j);
#endif
                        out.accreted_BH_Mass += FLT(BPP(j).BH_Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        out.accreted_BH_Mass_alphadisk += FLT(BPP(j).BH_Mass_AlphaDisk);
#endif
#ifdef BH_WIND_SPAWN
                        out_accreted_BH_Mass_alphaornot += FLT(BPP(j).unspawned_wind_mass);
#endif
#ifdef BH_COUNTPROGS
                        out.BH_CountProgs += BPP(j).BH_CountProgs;
#endif
                        bin = P[j].TimeBin;
                        #pragma omp atomic
                        TimeBin_BH_mass[bin] -= BPP(j).BH_Mass;
                        #pragma omp atomic
                        TimeBin_BH_dynamicalmass[bin] -= Mass_j;
                        #pragma omp atomic
                        TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
                        if(BPP(j).BH_Mass > 0) {
                            #pragma omp atomic
                            TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;
                        }
                        Mass_j = 0;
#ifdef GALSF
                        out.Accreted_Age = P[j].StellarAge;
#endif
                        #pragma omp atomic
                        N_BH_swallowed++;
                    } // if(P[j].Type == 5) -- BH + BH merger


                    
                    
#ifdef BH_GRAVCAPTURE_NONGAS /* DM and star particles can only be accreted ifdef BH_GRAVCAPTURE_NONGAS */
                    if((P[j].Type == 1) || (All.ComovingIntegrationOn && (P[j].Type==2||P[j].Type==3)) )
                    {   /* this is a DM particle: In this case, no kick, so just zero out the mass and 'get rid of' the particle (preferably by putting it somewhere irrelevant) */
#ifdef BH_OUTPUT_MOREINFO
                        printf(" ..BH_swallow_DM: j %d Type(j) %d  M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g \n", j,P[j].Type,Mass_j,Vel_j[0],Vel_j[1],Vel_j[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],local.Pos[0],local.Pos[1],local.Pos[2]);
#endif
                        out.accreted_Mass += FLT(Mass_j); out.accreted_BH_Mass += FLT(Mass_j);
                        Mass_j = 0;
                        #pragma omp atomic
                        N_dm_swallowed++;
                    }
                    if((P[j].Type==4) || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn) ))
                    {   /* this is a star particle: If there is an alpha-disk, we let them go to the disk. If there is no alpha-disk, stars go to the BH directly and won't affect feedback. (Can be simply modified if we need something different.) */
                        out.accreted_Mass += FLT(Mass_j);
                        out_accreted_BH_Mass_alphaornot += FLT(Mass_j);
                        Mass_j = 0;
                        #pragma omp atomic
                        N_star_swallowed++;
                    }
#endif // #ifdef BH_GRAVCAPTURE_NONGAS -- BH + DM or Star merger


                    
                    
                    /* this is a gas particle: DAA: we need to see if the gas particle has to be accreted in full or not, depending on BH_WIND_KICK
                     the only difference with BH_ALPHADISK_ACCRETION should be that the mass goes first to the alphadisk */
                    if(P[j].Type == 0)
                    {
                        out.accreted_Mass += FLT(f_accreted*Mass_j);
#ifdef BH_GRAVCAPTURE_GAS
                        out_accreted_BH_Mass_alphaornot += FLT(f_accreted*Mass_j);
#endif
			double Mass_initial = Mass_j; // save this for possible IO below
                        Mass_j *= (1-f_accreted);
#ifdef BH_WIND_KICK     /* BAL kicking operations. NOTE: we have two separate BAL wind models, particle kicking and smooth wind model. This is where we do the particle kicking BAL model. This should also work when there is alpha-disk. */
                        double v_kick=All.BAL_v_outflow, dir[3]; for(k=0;k<3;k++) {dir[k]=dpos[k];} // DAA: default direction is radially outwards
#if defined(BH_COSMIC_RAYS) /* inject cosmic rays alongside wind injection */
                        double dEcr = All.BH_CosmicRay_Injection_Efficiency * Mass_j * (All.BAL_f_accretion/(1.-All.BAL_f_accretion)) * C_LIGHT_CODE*C_LIGHT_CODE;
                        inject_cosmic_rays(dEcr,All.BAL_v_outflow,5,j,dir);
#endif
#if (BH_WIND_KICK < 0)  /* DAA: along polar axis defined by angular momentum within Kernel (we could add finite opening angle) work out the geometry w/r to the plane of the disk */
                        if((dir[0]*J_dir[0] + dir[1]*J_dir[1] + dir[2]*J_dir[2]) > 0){for(k=0;k<3;k++) {dir[k]=J_dir[k];}} else {for(k=0;k<3;k++) {dir[k]=-J_dir[k];}}
#endif
                        for(k=0,norm=0;k<3;k++) {norm+=dir[k]*dir[k];} if(norm<=0) {dir[0]=0;dir[1]=0;dir[2]=1;norm=1;} else {norm=sqrt(norm); dir[0]/=norm;dir[1]/=norm;dir[2]/=norm;}
                        for(k=0;k<3;k++) {Vel_j[k]+=v_kick*All.cf_atime*dir[k];}
#ifdef GALSF_SUBGRID_WINDS // if sub-grid galactic winds are decoupled from the hydro, we decouple the BH kick winds as well
                        #pragma omp atomic write
                        SphP[j].DelayTime = All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
#endif
#ifdef BH_OUTPUT_MOREINFO
                        printf(" ..BAL kick: P[j].ID %llu ID %llu Type(j) %d f_acc %g M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g v_out %g \n",(unsigned long long) P[j].ID, (unsigned long long) P[j].SwallowID,P[j].Type, All.BAL_f_accretion,Mass_j,Vel_j[0],Vel_j[1],Vel_j[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],local.Pos[0],local.Pos[1],local.Pos[2],v_kick);
                        fprintf(FdBhWindDetails,"%g  %llu %g  %2.7f %2.7f %2.7f  %2.7f %2.7f %2.7f %g %g %g %llu  %2.7f %2.7f %2.7f\n",All.Time, (unsigned long long)P[j].ID, Mass_j, P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],  Vel_j[0],Vel_j[1],Vel_j[2],dir[0]/norm,dir[1]/norm,dir[2]/norm, (unsigned long long)local.ID, local.Pos[0],local.Pos[1],local.Pos[2]); fflush(FdBhWindDetails);
#endif
#endif // #ifdef BH_WIND_KICK
                        #pragma omp atomic
                        N_gas_swallowed++;
#ifdef OUTPUT_SINK_ACCRETION_HIST
                        MyDouble tempB[3]={0,0,0};
#ifdef MAGNETIC
                        for(k=0;k<3;k++) {tempB[k]=Get_Gas_BField(j,k);} //use particle magnetic field
#endif
                        fprintf(FdBhSwallowDetails,"%g %llu %g %2.7f %2.7f %2.7f %llu %g %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f\n", All.Time, (unsigned long long)local.ID,local.Mass,local.Pos[0],local.Pos[1],local.Pos[2],  (unsigned long long)P[j].ID, Mass_initial, (P[j].Pos[0]-local.Pos[0]),(P[j].Pos[1]-local.Pos[1]),(P[j].Pos[2]-local.Pos[2]), (Vel_j[0]-local.Vel[0]),(Vel_j[1]-local.Vel[1]),(Vel_j[2]-local.Vel[2]), SphP[j].InternalEnergy, tempB[0], tempB[1], tempB[2], SphP[j].Density); fflush(FdBhSwallowDetails);
#endif
                    }  // if(P[j].Type == 0)
                    //P[j].SwallowID = 0; /* DAA: make sure it is not accreted (or ejected) by the same BH again if inactive in the next timestep [PFH: no longer necessary with the new way we re-initialize the SwallowIDs] */
                } // if(P[j].SwallowID == id)  -- particles being entirely or partially swallowed!!!

                
                
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
                /* now, do any other feedback "kick" operations (which used the previous loops to calculate weights) */
                if(mom>0 && local.Dt>0 && OriginallyMarkedSwallowID==0 && P[j].SwallowID==0 && Mass_j>0 && P[j].Type==0) // particles NOT being swallowed!
                {
                    double r=0, dir[3]; for(k=0;k<3;k++) {dir[k]=dpos[k]; r+=dir[k]*dir[k];} // should be away from BH
                    if(r>0)
                    {
                        r=sqrt(r); for(k=0;k<3;k++) {dir[k]/=r;} /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                        for(norm=0,k=0;k<3;k++) {norm+=dir[k]*J_dir[k];}
                        mom_wt = bh_angleweight_localcoupling(j,norm,r,h_i) / local.BH_angle_weighted_kernel_sum;
                        if(local.BH_angle_weighted_kernel_sum<=0) {mom_wt=0;}

#ifdef BH_PHOTONMOMENTUM /* inject radiation pressure: add initial L/c optical/UV coupling to the gas at the dust sublimation radius */
                        double v_kick = All.BH_Rad_MomentumFactor * mom_wt * mom / Mass_j;
                        for(k=0;k<3;k++) {Vel_j[k]+=v_kick*All.cf_atime*dir[k];}
#endif
#if defined(BH_COSMIC_RAYS) && defined(BH_WIND_CONTINUOUS) /* inject cosmic rays alongside continuous wind injection */
                        double dEcr = (evaluate_blackhole_cosmicray_efficiency(local.Mdot,local.BH_Mass,-1)) * mom_wt * C_LIGHT_CODE*C_LIGHT_CODE * local.Mdot*local.Dt;
                        inject_cosmic_rays(dEcr,All.BAL_v_outflow,5,j,dir);
#endif
#if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK) /* inject BAL winds, this is the more standard smooth feedback model */
                        double m_wind = mom_wt * (1-All.BAL_f_accretion)/(All.BAL_f_accretion) * local.Mdot*local.Dt; /* mass to couple */
                        if(local.BH_angle_weighted_kernel_sum<=0) {m_wind=0;}
                        //1. check if (Vw-V0)*rhat <= 0   [ equivalently, check if   |Vw| <= V0*rhat ]
                        //2. if (1) is False, the wind will catch the particle, couple mass, momentum, energy, according to the equations above
                        //3. if (1) is True, the wind will not catch the particle, or will only asymptotically catch it. For the sake of mass conservation in the disk, I think it is easiest to treat this like the 'marginal' case where the wind barely catches the particle. In this case, add the mass normally, but no momentum, and no energy, giving:
                        //dm = m_wind, dV = 0, du = -mu*u0   [decrease the thermal energy slightly to account for adding more 'cold' material to it]
                        double dvr_gas_to_bh, dr_gas_to_bh;
                        for(dvr_gas_to_bh=dr_gas_to_bh=0, k=0;k<3;k++) {dvr_gas_to_bh += dvel[k]*dpos[k]; dr_gas_to_bh += dpos[k]*dpos[k];}
                        dvr_gas_to_bh /= dr_gas_to_bh ;

                        /* add wind mass to particle. here we ignore the density correction to the particle from the updated mass, relying on the next-step to catch this. if injected mass is very large, this can cause problems! */
                        Mass_j += m_wind;

                        /* now add wind momentum to particle */
                        if(dvr_gas_to_bh < All.BAL_v_outflow)   // gas moving away from BH at v < BAL speed
                        {
                            double e_wind = 0;
                            for(k=0;k<3;k++)
                            {
                                norm = All.cf_atime*All.BAL_v_outflow*dir[k] - dvel[k]; // relative wind-particle velocity (in code units) including BH-particle motion;
                                Vel_j[k] += All.BlackHoleFeedbackFactor * norm * m_wind/Mass_j; // momentum conservation gives updated velocity
                                e_wind += (norm/All.cf_atime)*(norm/All.cf_atime); // -specific- shocked wind energy
                            }
                            e_wind *= 0.5*m_wind/Mass_j; // make total wind energy, add to particle as specific energy of -particle-
                            InternalEnergy_j += e_wind;
                        } else {    // gas moving away from BH at wind speed (or faster) already.
                            if(InternalEnergy_j * ( Mass_j - m_wind ) / Mass_j > 0) {InternalEnergy_j = InternalEnergy_j * ( Mass_j - m_wind ) / Mass_j;}
                        }
#endif // if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)
                    } // r > 0
                } // (check if valid gas neighbor of interest)
#endif // defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
                
                /* ok, now deal with coupling back terms to the actual neighbor 'j'. to do this we need to be careful to ensure thread-safety */
                if(Mass_j != Mass_j_0 || Vel_j[0]-Vel_j_0[0] != 0 || Vel_j[1]-Vel_j_0[1] != 0 || Vel_j[2]-Vel_j_0[2] != 0 || InternalEnergy_j != InternalEnergy_j_0)
                {
                    /* we updated variables that need to get assigned to element 'j' -- let's do it */
                    if(Mass_j > 0)
                    {
                        // do one more check here to make sure we don't accidentally 'resurrect' an element with zero mass set by a different process. technically it is conceivable something could come in the meantime but this makes it vanishingly unlikely [we hope]
                        double Mass_j_0_n = 0;
                        #pragma omp atomic read
                        Mass_j_0_n = P[j].Mass;
                        if(Mass_j_0_n > 0)
                        {
                            #pragma omp atomic
                            P[j].Mass += Mass_j - Mass_j_0; // finite mass update [delta difference added here, allowing for another element to update in the meantime]
                        }
                    } else {
                        #pragma omp atomic write
                        P[j].Mass = 0; // make sure the mass is -actually- zero'd here
                    }
                    
                    for(k=0;k<3;k++) {
                        #pragma omp atomic
                        P[j].Vel[k] += Vel_j[k] - Vel_j_0[k]; // discrete momentum change
                        #pragma omp atomic
                        P[j].dp[k] += Mass_j * Vel_j[k] - Mass_j_0*Vel_j_0[k]; // discrete momentum change
                    }

                    if(P[j].Type==0) // gas variable, need to update the 'pred' and truemass and related variables as well
                    {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        if(Mass_j > 0)
                        {
                            #pragma omp atomic
                            SphP[j].MassTrue += Mass_j - Mass_j_0; // finite mass update
                        } else {
                            #pragma omp atomic write
                            SphP[j].MassTrue = 0; // make sure the mass is -actually- zero'd here
                        }
#endif
                        for(k=0;k<3;k++) {
                            #pragma omp atomic
                            SphP[j].VelPred[k] += Vel_j[k] - Vel_j_0[k]; // delta-update
                        }
                        #pragma omp atomic
                        SphP[j].InternalEnergy += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                        #pragma omp atomic
                        SphP[j].InternalEnergyPred += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                    }
                }
                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    } // while(startnode >= 0) (outer of the double-loop)
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
} /* closes bh_evaluate_swallow routine */



void blackhole_swallow_and_kick_loop(void)
{
    N_gas_swallowed = N_star_swallowed = N_dm_swallowed = N_BH_swallowed = 0;
    #include "../../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    /* collect and print results on any swallow operations in this pass */
    int Ntot_gas_swallowed=0, Ntot_star_swallowed=0, Ntot_dm_swallowed=0, Ntot_BH_swallowed=0;
    MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_star_swallowed, &Ntot_star_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&N_dm_swallowed, &Ntot_dm_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if((ThisTask == 0)&&(Ntot_gas_swallowed+Ntot_star_swallowed+Ntot_dm_swallowed+Ntot_BH_swallowed>0))
    {
        printf("Accretion done: swallowed %d gas, %d star, %d dm, and %d BH particles\n",
               Ntot_gas_swallowed, Ntot_star_swallowed, Ntot_dm_swallowed, Ntot_BH_swallowed);
    }
    CPU_Step[CPU_BLACKHOLES] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */





#ifdef BH_WIND_SPAWN
void spawn_bh_wind_feedback(void)
{
    int i, n_particles_split = 0, MPI_n_particles_split, dummy_gas_tag=0;
    for(i = 0; i < NumPart; i++)
        if(P[i].Type==0)
        {
            dummy_gas_tag=i;
            break;
        }

    /* don't loop or go forward if there are no gas particles in the domain, or the code will crash */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        long nmax = (int)(0.99*All.MaxPart); if(All.MaxPart-20 < nmax) nmax=All.MaxPart-20;
        if((NumPart+n_particles_split+(int)(2.*(BH_WIND_SPAWN+0.1)) < nmax) && (P[i].Type == 5)) // basic condition: particle is a 'spawner' (sink), and code can handle the event safely without crashing.
        {
            int sink_eligible_to_spawn = 0; // flag to check eligibility for spawning
            if(BPP(i).unspawned_wind_mass >= (BH_WIND_SPAWN)*All.BAL_wind_particle_mass) {sink_eligible_to_spawn=1;} // have 'enough' mass to spawn
#if defined(SINGLE_STAR_SINK_DYNAMICS)
            if((P[i].Mass <= 3.5*All.MeanGasParticleMass) || (P[i].BH_Mass*UNIT_MASS_IN_SOLAR < 0.01)) {sink_eligible_to_spawn=0;}  // spawning causes problems in these modules for low-mass sinks, so arbitrarily restrict to this, since it's roughly a criterion on the minimum particle mass. and for <0.01 Msun, in pre-collapse phase, no jets
#endif
            if(sink_eligible_to_spawn)
            {
                int j; dummy_gas_tag=-1; double r2=MAX_REAL_NUMBER;
                for(j=0; j<N_gas; j++) /* find the closest gas particle on the domain to act as the dummy */
                {
                    if(P[j].Type==0)
                    {
                        if((P[j].Mass>0) && (SphP[j].Density>0))
                        {
                            double dx2=(P[j].Pos[0]-P[i].Pos[0])*(P[j].Pos[0]-P[i].Pos[0]) + (P[j].Pos[1]-P[i].Pos[1])*(P[j].Pos[1]-P[i].Pos[1]) + (P[j].Pos[2]-P[i].Pos[2])*(P[j].Pos[2]-P[i].Pos[2]);
                            if(dx2 < r2) {r2=dx2; dummy_gas_tag=j;}
                        }
                    }
                }
                if(dummy_gas_tag >= 0)
                {
                    n_particles_split += blackhole_spawn_particle_wind_shell( i , dummy_gas_tag, n_particles_split);
                }
            }
        }
    }
    MPI_Allreduce(&n_particles_split, &MPI_n_particles_split, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(MPI_n_particles_split>0){
#ifdef MAINTAIN_TREE_IN_REARRANGE
        All.NumForcesSinceLastDomainDecomp +=  0.0001 * All.TreeDomainUpdateFrequency * All.TotNumPart; // we can insert spawned particles in the tree, but still a good idea to rebuild the tree every now and then, so we make the next domain+treebuild come a bit sooner; additional cost should be small
#else
        TreeReconstructFlag = 1; // otherwise just wipe and rebuild the tree next chance you get - more expensive but more accurate
#endif
        if(ThisTask == 0) {printf(" ..BH-Spawn Event: %d particles spawned \n", MPI_n_particles_split);}
    }

    /* rearrange_particle_sequence -must- be called immediately after this routine! */
    All.TotNumPart += (long long)MPI_n_particles_split;
    All.TotN_gas   += (long long)MPI_n_particles_split;
    Gas_split       = n_particles_split;                    // specific to the local processor //
}


void get_random_orthonormal_basis(int seed, double *nx, double *ny, double *nz)
{
    double phi, cos_theta, sin_theta, sin_phi, cos_phi;
    phi=2.*M_PI*get_random_number(seed+1+ThisTask), cos_theta=2.*(get_random_number(seed+3+2*ThisTask)-0.5); sin_theta=sqrt(1-cos_theta*cos_theta), sin_phi=sin(phi), cos_phi=cos(phi);
    /* velocities (determined by wind velocity direction) */
    nz[0]=sin_theta*cos_phi; nz[1]=sin_theta*sin_phi; nz[2]=cos_theta; // random z axis

    double dot_product=0, norm=0; int k;
    while(norm==0){ // necessary in case ny is parallel to nz - believe it or not this happened once!
        phi=2.*M_PI*get_random_number(seed+4+ThisTask), cos_theta=2.*(get_random_number(seed+5+2*ThisTask)-0.5); sin_theta=sqrt(1-cos_theta*cos_theta), sin_phi=sin(phi), cos_phi=cos(phi);
        ny[0]=sin_theta*cos_phi; ny[1]=sin_theta*sin_phi; ny[2]=cos_theta; // random y axis, needs to have its z component deprojected
        // do Gram-Schmidt to get an orthonormal basis
        for(k=0; k<3; k++) {dot_product += ny[k] * nz[k];}
        for(k=0; k<3; k++) {ny[k] -= dot_product * nz[k]; norm += ny[k]*ny[k];} // deproject component along z
        if(norm==0) continue;
        norm = 1./sqrt(norm); for(k=0; k<3; k++) {ny[k] *= norm;}
    }
    nx[0] = ny[1]*nz[2] - ny[2]*nz[1];
    nx[1] = ny[2]*nz[0] - ny[0]*nz[2];
    nx[2] = ny[0]*nz[1] - ny[1]*nz[0];
    return;
}

/* Convenience function to compute the direction to launch a wind particle                                      */
/*                                                                                                              */
/* i - index of particle doing the spawning                                                                     */
/* num_spawned_this_call - how many we have already spawned in this call of blackhole_spawn_particle_wind_shell */
/* mode - 0 for random, 1 for collimated, 2 for isotropized random, 3 for angular grid                          */
/* ny, nz - shape (3,) arrays containing 2 vectors in the fixed orthonormal basis - for collimated winds, nz    */
/*          is the axis                                                                                         */
/* dir - shape (3,) array containing the direction - pass as an input to remember the previous direction        */

void get_wind_spawn_direction(int i, int num_spawned_this_call, int mode, double *ny, double *nz, double *veldir, double *dpdir)
{
    int k;
    if((mode != 3) && (num_spawned_this_call % 2)) { // every second particle is spawned in the opposite direction to the last, conserving momentum and COM
        for(k=0; k<3;k++) {veldir[k] = -veldir[k]; dpdir[k]=-dpdir[k];}
        return; // we're done
    }
    double nx[3] = {ny[1]*nz[2] - ny[2]*nz[1], ny[2]*nz[0] - ny[0]*nz[2], ny[0]*nz[1] - ny[1]*nz[0]};
    // now do the actual direction based on the mode we're in
    double phi, cos_theta, sin_theta, sin_phi, cos_phi;
    if(mode==0){ // fully random
        phi=2.*M_PI*get_random_number(num_spawned_this_call+1+ThisTask), cos_theta=2.*(get_random_number(num_spawned_this_call+3+2*ThisTask)-0.5); sin_theta=sqrt(1-cos_theta*cos_theta), sin_phi=sin(phi), cos_phi=cos(phi);
        veldir[0]=sin_theta*cos_phi; veldir[1]=sin_theta*sin_phi; veldir[2]=cos_theta;dpdir[0]=veldir[0];dpdir[1]=veldir[1];dpdir[2]=veldir[2];
    } else if (mode==1){ // collimated according to a conical velocity field
        double theta0=0.01, thetamax=30.*(M_PI/180.); // "flattening parameter" and max opening angle of jet velocity distribution from Matzner & McKee 1999, sets the collimation of the jets
        double theta=atan(theta0*tan(get_random_number(num_spawned_this_call+7+5*ThisTask)*atan(sqrt(1+theta0*theta0)*tan(thetamax)/theta0))/sqrt(1+theta0*theta0)); // biased sampling to get collimation
        phi=2.*M_PI*get_random_number(num_spawned_this_call+1+ThisTask);
        cos_theta = cos(theta), sin_theta=sin(theta), sin_phi=sin(phi), cos_phi=cos(phi);
        for(k=0;k<3;k++) {veldir[k] = sin_theta*cos_phi*nx[k] + sin_theta*sin_phi*ny[k] + cos_theta*nz[k]; dpdir[k]=veldir[k];} //converted from angular momentum relative to into standard coordinates
    }
#if defined(BH_DEBUG_SPAWN_JET_TEST)
    else if (mode==4) { // old-style of jet being initialized as a cylinder around the BH
        phi=2.*M_PI*get_random_number(num_spawned_this_call+1+ThisTask), cos_theta=2.*(get_random_number(num_spawned_this_call+3+2*ThisTask)-0.5); sin_theta=sqrt(1-cos_theta*cos_theta), sin_phi=sin(phi), cos_phi=cos(phi);
        if(P[i].ID_child_number % 2 == 0) {cos_theta=fabs(cos_theta);} else {cos_theta=-1.0*fabs(cos_theta);} // balance vertical directions
        double ct_v=1.-(1-cos((BH_DEBUG_SPAWN_JET_TEST)*(M_PI/180.)))*(1.-fabs(cos_theta)), st_v=sqrt(1-ct_v*ct_v); if(cos_theta<0) {ct_v*=-1;}
        for(k=0;k<3;k++) {dpdir[k] = sin_theta*cos_phi*nx[k] + sin_theta*sin_phi*ny[k] + cos_theta*nz[k]; veldir[k] = st_v*cos_phi*nx[k] + st_v*sin_phi*ny[k] + ct_v*nz[k];}
    }
#endif
    return;
}

#ifdef MAGNETIC
void get_wind_spawn_magnetic_field(int j, int mode, double *ny, double *nz,  double *dpdir, double d_r)
{
    int k; SphP[j].divB = 0;
#ifdef DIVBCLEANING_DEDNER
    SphP[j].DtPhi = SphP[j].PhiPred = SphP[j].Phi = 0;
#endif

#ifdef BH_WIND_SPAWN_SET_BFIELD_POLTOR /* user manually sets the poloidal and toroidal components here */
    double inj_scale=All.BH_spawn_rinj/All.cf_atime, Bfield[3]={0}, nx[3]={ny[1]*nz[2]-ny[2]*nz[1], ny[2]*nz[0]-ny[0]*nz[2], ny[0]*nz[1]-ny[1]*nz[0]};
    double cos_theta=nz[0]*dpdir[0]+nz[1]*dpdir[1]+nz[2]*dpdir[2], sin_theta=sqrt(1-cos_theta*cos_theta), cos_phi=(nx[0]*dpdir[0]+nx[1]*dpdir[1]+nx[2]*dpdir[2])/sin_theta, sin_phi=(ny[0]*dpdir[0]+ny[1]*dpdir[1]+ny[2]*dpdir[2])/sin_theta;
    /* initialize poloidal component, in the nx,ny,nz coordinate frame */
    Bfield[0]+= All.B_spawn_pol*d_r*cos_theta*d_r*sin_theta/inj_scale/inj_scale*cos_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    Bfield[1]+= All.B_spawn_pol*d_r*cos_theta*d_r*sin_theta/inj_scale/inj_scale*sin_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    Bfield[2]+= All.B_spawn_pol*(1-d_r*cos_theta*d_r*sin_theta/inj_scale/inj_scale)    *exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    /* initialize toroidal component, in the nx,ny,nz coordinate frame */
    Bfield[0]+= -1*All.B_spawn_tor*(d_r/inj_scale)*sin_theta*sin_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    Bfield[1]+=    All.B_spawn_tor*(d_r/inj_scale)*sin_theta*cos_phi*exp(-1.0*d_r*d_r/inj_scale/inj_scale)/exp(-1.0);
    /* assign it back to the actual evolved B in the lab/simulation coordinate frame */
    for(k=0;k<3;k++) {SphP[j].IniB[k] = Bfield[0]*nx[k] + Bfield[1]*ny[k] + Bfield[2]*nz[k]; SphP[j].DtB[k]=0;
        SphP[j].BPred[k]=SphP[j].B[k]=SphP[j].IniB[k]*(All.UnitMagneticField_in_gauss/UNIT_B_IN_GAUSS)*(P[j].Mass/(All.cf_a2inv*SphP[j].Density));}
    
#else /* set B-fields to be weak relative to local ISM values */

    double Bmag=0, Bmag_0=0;
    for(k=0;k<3;k++) {double B=SphP[j].B[k]*SphP[j].Density/P[j].Mass*All.cf_a2inv; Bmag+=B*B; Bmag_0+=SphP[j].B[k]*SphP[j].B[k];} // get actual Bfield
    double Bmag_low_rel_to_progenitor = 1.e-10 * sqrt(Bmag); // set to some extremely low value relative to cloned element
    double u_internal_new_cell = All.BAL_internal_temperature / (  0.59 * (5./3.-1.) * U_TO_TEMP_UNITS ); // internal energy of new wind cell
    double Bmag_low_rel_to_pressure = 1.e-3 * sqrt(2.*SphP[j].Density*All.cf_a3inv * u_internal_new_cell); // set to beta = 1e6
    Bmag = DMAX(Bmag_low_rel_to_progenitor , Bmag_low_rel_to_pressure); // pick the larger of these (still small) B-field values
#ifdef MHD_B_SET_IN_PARAMS
    double Bmag_IC = sqrt(All.BiniX*All.BiniX + All.BiniY*All.BiniY + All.BiniZ*All.BiniZ) * All.UnitMagneticField_in_gauss / UNIT_B_IN_GAUSS; // IC B-field sets floor as well
    Bmag = DMAX(Bmag , 0.1 * Bmag_IC);
#endif
#if defined(SINGLE_STAR_FB_SNE)
    if(P[j].ProtoStellarStage == 6) {Bmag=0;} // No need to have flux in SN ejecta - note that this assumes we inherited this attribute from the spawning sink before calling this routine
#endif
    Bmag = DMAX(Bmag, MIN_REAL_NUMBER); // floor to prevent underflow errors
    /* add magnetic flux here to 'Bmag' if desired */
    Bmag *= P[j].Mass / (All.cf_a2inv * SphP[j].Density); // convert back to code units
    for(k=0;k<3;k++) {if(Bmag_0>0) {SphP[j].B[k]*=Bmag/sqrt(Bmag_0);} else {SphP[j].B[k]=Bmag;}} // assign if valid values
    for(k=0;k<3;k++) {SphP[j].BPred[k]=SphP[j].B[k]; SphP[j].DtB[k]=0;} // set predicted = actual, derivative to null
    
#endif

    return;
}
#endif


/*! this code copies what was used in merge_split.c for the gas particle split case */
int blackhole_spawn_particle_wind_shell( int i, int dummy_sph_i_to_clone, int num_already_spawned )
{
    double total_mass_in_winds = BPP(i).unspawned_wind_mass;
    int n_particles_split   = (int) floor( total_mass_in_winds / All.BAL_wind_particle_mass ); /* if we set BH_WIND_SPAWN we presumably wanted to do this in an exactly-conservative manner, which means we want to have an even number here. */
    int k=0; long j;


    if((((int)BH_WIND_SPAWN) % 2) == 0) {if(( n_particles_split % 2 ) != 0) {n_particles_split -= 1;}} /* n_particles_split was not even. we'll wait to spawn this last particle, to keep an even number, rather than do it right now and break momentum conservation */
    if( (n_particles_split == 0) || (n_particles_split < 1) ) {return 0;}
    int n0max = DMAX(20 , (int)(3.*(BH_WIND_SPAWN)+0.1)); if((n0max % 2) != 0) {n0max += 1;} // should ensure n0max is always an even number //
    if(n_particles_split > n0max) {n_particles_split = n0max;}


    /* here is where the details of the split are coded, the rest is bookkeeping */
    //double mass_of_new_particle = total_mass_in_winds / n_particles_split; /* don't do this, as can produce particles with extremely large masses; instead wait to spawn */
    double mass_of_new_particle = All.BAL_wind_particle_mass;
    printf("Task %d wants to create %g mass in wind with %d new particles each of mass %g \n .. splitting BH %d using hydro element %d\n", ThisTask,total_mass_in_winds, n_particles_split, mass_of_new_particle, i, dummy_sph_i_to_clone);

    if(NumPart + num_already_spawned + n_particles_split >= All.MaxPart)
    {
        printf("On Task=%d with NumPart=%d (+N_spawned=%d) we tried to split a particle, but there is no space left...(All.MaxPart=%d). Try using more nodes, or raising PartAllocFac, or changing the split conditions to avoid this.\n", ThisTask, NumPart, num_already_spawned, All.MaxPart);
        fflush(stdout); endrun(8888);
    }
    double d_r = 0.25 * KERNEL_CORE_SIZE*PPP[i].Hsml; // needs to be epsilon*Hsml where epsilon<<1, to maintain stability //
    double r2=0; for(k=0;k<3;k++) {r2+=(P[dummy_sph_i_to_clone].Pos[k]-P[i].Pos[k])*(P[dummy_sph_i_to_clone].Pos[k]-P[i].Pos[k]);}
    d_r = DMIN(d_r, 0.5*sqrt(r2));
#ifndef SELFGRAVITY_OFF
    d_r = DMAX(d_r , 2.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0]);
#endif
#ifdef BH_WIND_SPAWN_SET_BFIELD_POLTOR
    d_r = DMIN(d_r , All.BH_spawn_rinj/All.cf_atime); /* KYSu: sets spawn scale manually */
#endif
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
    d_r = DMIN(P[i].SinkRadius, d_r); //launch close to the sink
#endif
    long bin, bin_0; for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinCount[bin] > 0) break;} /* gives minimum active timebin of any particle */
    bin_0 = bin; int i0 = i; /* save minimum timebin, also save ID of BH particle for use below */
    bin = P[i0].TimeBin; /* make this particle active on the BH/star timestep */
#ifdef BH_DEBUG_SPAWN_JET_TEST
    bin = bin_0; i0 = dummy_sph_i_to_clone; /* make this particle active on the minimum timestep, and order with respect to the cloned particle */
#endif

    double veldir[3], dpdir[3]; // velocity direction to spawn in - declare outside the loop so we remember it from the last iteration
    int mode = 0; // 0 if doing totally random directions, 1 if collimated, 2 for 3-axis isotropized, and 3 if using an angular grid,  4 old collimatation script, position isotropic velicity coliminated within certain open angle (might be useful to still keep this option owing to the free open angle choice and better sampling the magnetic field geometry)
#if defined(SINGLE_STAR_FB_JETS) || defined(JET_DIRECTION_FROM_KERNEL_AND_SINK) || defined(BH_FB_COLLIMATED)
    mode = 1; // collimated mode
#endif
#if defined(BH_DEBUG_SPAWN_JET_TEST)
    mode = 4;
#endif

    // based on the mode we're in, let's pick a fixed orthonormal basis that all spawned elements are aware of
    double jz[3]={0,0,1},jy[3]={0,1,0},jx[3]={1,0,0};  /* set up a coordinate system [xyz if we don't have any other information */
#ifdef BH_FOLLOW_ACCRETED_ANGMOM  /* use local angular momentum to estimate preferred directions/coordinates for spawning */
    if(mode==1){ // set up so that the z axis is the angular momentum vector
#ifdef JET_DIRECTION_FROM_KERNEL_AND_SINK // Jgas stores total angmom in COM frame of sink-gas system; use this for direction
        double Jtot=0; for(k=0;k<3;k++) {Jtot+=P[i].Jgas_in_Kernel[k]*P[i].Jgas_in_Kernel[k];}
        if(Jtot>0) {Jtot=1/sqrt(Jtot); for(k=0;k<3;k++) {jz[k]=P[i].Jgas_in_Kernel[k]*Jtot;}}
#else
        double Jtot=0; for(k=0;k<3;k++) {Jtot+=P[i].BH_Specific_AngMom[k]*P[i].BH_Specific_AngMom[k];}
        if(Jtot>0) {Jtot=1/sqrt(Jtot); for(k=0;k<3;k++) {jz[k]=P[i].BH_Specific_AngMom[k]*Jtot;}}
#endif
        Jtot=jz[1]*jz[1]+jz[2]*jz[2]; if(Jtot>0) {Jtot=1/sqrt(Jtot); jy[1]=jz[2]*Jtot; jy[2]=-jz[1]*Jtot;}
        jx[0]=jz[1]*jy[2]-jz[2]*jy[1]; jx[1]=jz[2]*jy[0]-jz[0]*jy[2]; jx[2]=jz[0]*jy[1]-jz[1]*jy[0];
    }
#endif
    if(mode == 3){ // if doing an angular grid, need some fixed coordinates to orient it, but want to switch em up each time to avoid artifacts
        get_random_orthonormal_basis(P[i].ID_child_number, jx, jy, jz);
    }
#ifdef BH_WIND_SPAWN_SET_JET_PRECESSION /* rotate the jet angle according to the explicitly-included precession parameters */
    double degree = All.BH_jet_precess_degree, period = All.BH_jet_precess_period/UNIT_TIME_IN_GYR, new_dir[3];
    new_dir[0]= jx[0]*cos(degree/180.*M_PI)-jx[2]*sin(degree/180.*M_PI); new_dir[1]= 1.0*jx[1]; new_dir[2]= jx[0]*sin(degree/180.*M_PI)+jx[2]*cos(degree/180.*M_PI);
    jx[0]= new_dir[0]*cos(2.*M_PI/period*All.Time)-new_dir[1]*sin(2.*M_PI/period*All.Time); jx[1]= new_dir[0]*sin(2.*M_PI/period*All.Time)+new_dir[1]*cos(2.*M_PI/period*All.Time); jx[2]= new_dir[2];

    new_dir[0]= jy[0]*cos(degree/180.*M_PI)-jy[2]*sin(degree/180.*M_PI); new_dir[1]= 1.0*jy[1]; new_dir[2]= jy[0]*sin(degree/180.*M_PI)+jy[2]*cos(degree/180.*M_PI);
    jy[0]= new_dir[0]*cos(2.*M_PI/period*All.Time)-new_dir[1]*sin(2.*M_PI/period*All.Time); jy[1]= new_dir[0]*sin(2.*M_PI/period*All.Time)+new_dir[1]*cos(2.*M_PI/period*All.Time); jy[2]= new_dir[2];

    new_dir[0]= jz[0]*cos(degree/180.*M_PI)-jz[2]*sin(degree/180.*M_PI); new_dir[1]= 1.0*jz[1]; new_dir[2]= jz[0]*sin(degree/180.*M_PI)+jz[2]*cos(degree/180.*M_PI);
    jz[0]= new_dir[0]*cos(2.*M_PI/period*All.Time)-new_dir[1]*sin(2.*M_PI/period*All.Time); jz[1]= new_dir[0]*sin(2.*M_PI/period*All.Time)+new_dir[1]*cos(2.*M_PI/period*All.Time); jz[2]= new_dir[2];
#endif    

    /* create the  new particles to be added to the end of the particle list :
        i is the BH particle tag, j is the new "spawed" particle's location, dummy_sph_i_to_clone is a dummy SPH particle's tag to be used to init the wind particle */
    for(j = NumPart + num_already_spawned; j < NumPart + num_already_spawned + n_particles_split; j++)
    {   /* first, clone the 'dummy' particle so various fields are set appropriately */
        P[j] = P[dummy_sph_i_to_clone]; SphP[j] = SphP[dummy_sph_i_to_clone]; /* set the pointers equal to one another -- all quantities get copied, we only have to modify what needs changing */

        /* now we need to make sure everything is correctly placed in timebins for the tree */
        P[j].TimeBin = bin; // get the timebin, and put this particle into the appropriate timebin
        NextActiveParticle[j] = FirstActiveParticle; FirstActiveParticle = j; NumForceUpdate++;
        TimeBinCount[bin]++; TimeBinCountSph[bin]++; PrevInTimeBin[j] = i0; /* likewise add it to the counters that register how many particles are in each timebin */
#ifndef BH_DEBUG_SPAWN_JET_TEST
        NextInTimeBin[j] = NextInTimeBin[i0]; if(NextInTimeBin[i0] >= 0) {PrevInTimeBin[NextInTimeBin[i0]] = j;} NextInTimeBin[i0] = j; if(LastInTimeBin[bin] == i0) {LastInTimeBin[bin] = j;}
#else
        if(FirstInTimeBin[bin] < 0) {FirstInTimeBin[bin]=j; LastInTimeBin[bin]=j; NextInTimeBin[j]=-1; PrevInTimeBin[j]=-1;} /* only particle in this time bin on this task */
            else {NextInTimeBin[j]=FirstInTimeBin[bin]; PrevInTimeBin[j]=-1; PrevInTimeBin[FirstInTimeBin[bin]]=j; FirstInTimeBin[bin]=j;} /* there is already at least one particle; add this one "to the front" of the list */
#endif
        P[j].Ti_begstep = All.Ti_Current; P[j].Ti_current = All.Ti_Current;
#ifdef WAKEUP /* note - you basically MUST have this flag on for this routine to work at all -- */
        P[j].dt_step = GET_INTEGERTIME_FROM_TIMEBIN(bin);
        PPPZ[j].wakeup = 1;
        NeedToWakeupParticles_local = 1;
#endif
        /* this is a giant pile of variables to zero out. dont need everything here because we cloned a valid particle, but handy anyways */
        P[j].Particle_DivVel = 0; SphP[j].DtInternalEnergy = 0; for(k=0;k<3;k++) {SphP[j].HydroAccel[k] = 0; P[j].GravAccel[k] = 0;}
        P[j].NumNgb=cbrt(All.DesNumNgb); // this gets cube rooted at the end of the density loop, so take cbrt here
#ifdef PMGRID
        for(k=0;k<3;k++) {P[j].GravPM[k] = 0;}
#endif
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        SphP[j].MaxKineticEnergyNgb = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].dMass = 0; SphP[j].DtMass = 0; SphP[j].MassTrue = P[j].Mass; for(k=0;k<3;k++) {SphP[j].GravWorkTerm[k] = 0;}
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        PPPZ[j].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        PPP[j].AGS_Hsml = PPP[j].Hsml;
#endif
#endif
#ifdef ADAPTIVE_TREEFORCE_UPDATE
        P[j].tdyn_step_for_treeforce = 0;
        P[j].time_since_last_treeforce = MAX_REAL_NUMBER; // make sure we get a new tree force right off the bat
#endif
#ifdef CONDUCTION
        SphP[j].Kappa_Conduction = 0;
#endif
#ifdef MHD_NON_IDEAL
        SphP[j].Eta_MHD_OhmicResistivity_Coeff = 0; SphP[j].Eta_MHD_HallEffect_Coeff = 0; SphP[j].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
#ifdef VISCOSITY
        SphP[j].Eta_ShearViscosity = 0; SphP[j].Zeta_BulkViscosity = 0;
#endif
#ifdef TURB_DIFFUSION
        SphP[j].TD_DiffCoeff = 0;
#endif
#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==1)
        SphP[j].HostHaloMass = 0;
#endif
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        SphP[j].DelayTimeCoolingSNe = 0;
#endif
#ifdef GALSF
        SphP[j].Sfr = 0;
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        SphP[j].alpha = 0.0;
#endif
#if defined(BH_THERMALFEEDBACK)
        SphP[j].Injected_BH_Energy = 0;
#endif
#ifdef RADTRANSFER
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            SphP[j].Rad_E_gamma[k] = 0;
#if defined(RT_EVOLVE_ENERGY)
            SphP[j].Rad_E_gamma_Pred[k] = 0; SphP[j].Dt_Rad_E_gamma[k] = 0;
#endif
#if defined(RT_EVOLVE_FLUX)
            int kdir; for(kdir=0;kdir<3;kdir++){SphP[j].Rad_Flux[k][kdir] = 0;}
#endif
        }
#endif
#if defined(GALSF)
        P[j].StellarAge = All.Time; // use this attibute to save the gas cell's formation time for possible subsequent checks for special behavior on its first timestep
#endif
        

        /* now set the real hydro variables. */
        /* set the particle ID */ // unsigned int bits; int SPLIT_GENERATIONS = 4; for(bits = 0; SPLIT_GENERATIONS > (1 << bits); bits++); /* the particle needs an ID: we give it a bit-flip from the original particle to signify the split */
        P[j].ID = All.AGNWindID; /* update:  We are using a fixed wind ID, to allow for trivial wind particle identification */
#if defined(SINGLE_STAR_SINK_DYNAMICS)
        if(mass_of_new_particle >= 0.5*All.MeanGasParticleMass) {P[j].ID = All.AGNWindID + 1;} // this just has the nominal mass resolution, so no special treatment - this avoids the P[i].ID == All.AGNWindID checks throughout the code
#endif

        P[j].ID_child_number = P[i].ID_child_number; P[i].ID_child_number +=1; P[j].ID_generation = P[i].ID; // this allows us to track spawned particles by giving them unique sub-IDs

        P[j].Mass = mass_of_new_particle; /* assign masses to both particles (so they sum correctly) */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].MassTrue = P[j].Mass;
#endif
#ifndef BH_DEBUG_FIX_MDOT_MBH
        P[i].Mass -= P[j].Mass; /* make sure the operation is mass conserving! */
#endif
        BPP(i).unspawned_wind_mass -= P[j].Mass; /* remove the mass successfully spawned, to update the remaining unspawned mass */

        double v_magnitude = All.BAL_v_outflow * All.cf_atime; // velocity of the jet: default mode is to set this manually to a specific value
#ifdef SINGLE_STAR_FB_JETS
        v_magnitude = single_star_jet_velocity(i); // get velocity from our more detailed function
#endif
        

#if defined(METALS) && (defined(SINGLE_STAR_FB_JETS) || defined(SINGLE_STAR_FB_WINDS) || defined(SINGLE_STAR_FB_SNE))
        double yields[NUM_METAL_SPECIES]={0}; get_jet_yields(yields,i); // default to jet-type
        for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k]=yields[k];}  // update metallicity of spawned cell modules
#endif
        
        // actually lay down position and velocities using coordinate basis
        get_wind_spawn_direction(i, j - (NumPart + num_already_spawned), mode, jy, jz, veldir, dpdir);
        for(k=0;k<3;k++) {P[j].Pos[k]=P[i].Pos[k] + dpdir[k]*d_r; P[j].Vel[k]=P[i].Vel[k] + veldir[k]*v_magnitude; SphP[j].VelPred[k]=P[j].Vel[k];}

        /* condition number, smoothing length, and density */
        SphP[j].ConditionNumber *= 100.0; /* boost the condition number to be conservative, so we don't trigger madness in the kernel */
#if defined(SINGLE_STAR_SINK_DYNAMICS) 
        SphP[j].MaxSignalVel = 2*DMAX(v_magnitude, SphP[j].MaxSignalVel);// need this to satisfy the Courant condition in the first timestep after spawn
#endif
        /* note, if you want to use this routine to inject magnetic flux or cosmic rays, do this below */
#ifdef BH_DEBUG_SPAWN_JET_TEST
        PPP[j].Hsml=5.*d_r; SphP[j].Density=mass_of_new_particle/pow(KERNEL_CORE_SIZE*PPP[j].Hsml,NUMDIMS); /* PFH: need to write this in a way that does not make assumptions about units/problem structure */
#endif
#if defined(BH_WIND_SPAWN_SET_BFIELD_POLTOR) 
        SphP[j].IniDen = -1. * SphP[j].Density; /* this is essentially acting like a bitflag, to signal to the code that the density needs to be recalculated because a spawn event just occurred */
#endif
#ifdef MAGNETIC
        get_wind_spawn_magnetic_field(j, mode, jy, jz, dpdir, d_r);
#endif
#ifdef COSMIC_RAY_FLUID
        int k_CRegy; for(k_CRegy=0;k_CRegy<N_CR_PARTICLE_BINS;k_CRegy++) /* initialize CR energy and other related terms to nil */
        {
            SphP[j].CosmicRayEnergyPred[k_CRegy]=SphP[j].CosmicRayEnergy[k_CRegy]=SphP[j].DtCosmicRayEnergy[k_CRegy]=0;
#ifdef CRFLUID_M1
            for(k=0;k<3;k++) {SphP[j].CosmicRayFlux[k_CRegy][k]=SphP[j].CosmicRayFluxPred[k_CRegy][k]=0;}
#endif
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(k=0;k<3;k++) {SphP[j].CosmicRayAlfvenEnergy[k_CRegy][k]=SphP[j].CosmicRayAlfvenEnergyPred[k_CRegy][k]=SphP[j].DtCosmicRayAlfvenEnergy[k_CRegy][k]=0;}
#endif        
        } /* complete CR initialization to null */
#endif
        SphP[j].InternalEnergy = All.BAL_internal_temperature / (  0.59 * (5./3.-1.) * U_TO_TEMP_UNITS ); /* internal energy, determined by desired wind temperature (assume fully ionized primordial gas with gamma=5/3) */
        SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;

#if defined(BH_COSMIC_RAYS) /* inject cosmic rays alongside wind injection */
        double dEcr = evaluate_blackhole_cosmicray_efficiency(BPP(i).BH_Mdot,BPP(i).BH_Mass,i) * P[j].Mass * (All.BAL_f_accretion/(1.-All.BAL_f_accretion)) * C_LIGHT_CODE*C_LIGHT_CODE;
        inject_cosmic_rays(dEcr,All.BAL_v_outflow,5,j,veldir);
#endif
        /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
        force_add_star_to_tree(i0, j);// (buggy) /* we solve this by only calling the merge/split algorithm when we're doing the new domain decomposition */
    }
    if(BPP(i).unspawned_wind_mass < 0) {BPP(i).unspawned_wind_mass=0;}
    return n_particles_split;
}
#endif


#endif // top-level flag
