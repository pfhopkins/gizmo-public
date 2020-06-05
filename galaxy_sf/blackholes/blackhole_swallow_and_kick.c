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


#ifdef BLACK_HOLES // master flag [needs to be here to prevent compiler breaking when this is not active] //


static int N_gas_swallowed, N_star_swallowed, N_dm_swallowed, N_BH_swallowed;

#ifdef BH_ALPHADISK_ACCRETION
#define out_accreted_BH_Mass_alphaornot out.accreted_BH_Mass_alphadisk
#else
#define out_accreted_BH_Mass_alphaornot out.accreted_BH_Mass
#endif


#define MASTER_FUNCTION_NAME blackhole_swallow_and_kick_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int MASTER_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(P[i].Type==5 && P[i].SwallowID==0) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
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
#ifndef WAKEUP
    in->Dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
    in->Dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
#if defined(BH_RETURN_ANGMOM_TO_GAS)
    for(k=0;k<3;k++) {in->BH_Specific_AngMom[k] = BPP(i).BH_Specific_AngMom[k];}
    in->angmom_norm_topass_in_swallowloop = BlackholeTempInfo[j_tempinfo].angmom_norm_topass_in_swallowloop;
#endif
#if defined(BH_RETURN_BFLUX)
    for(k=0;k<3;k++) {in->B[k] = BPP(i).B[k];}
    in->kernel_norm_topass_in_swallowloop = BlackholeTempInfo[j_tempinfo].kernel_norm_topass_in_swallowloop;
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



int blackhole_swallow_and_kick_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration);
int blackhole_swallow_and_kick_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb, listindex = 0, j, k, n, bin; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* define variables and zero memory and import data for local target*/
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* imports the data to the correct place and names */
    double h_i=local.Hsml, hinv=1/h_i, hinv3, f_accreted; hinv3=hinv*hinv*hinv; f_accreted=0;
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    double kernel_zero,dwk; kernel_main(0.0,1.0,1.0,&kernel_zero,&dwk,-1); dwk=0;
    double mom = bh_lum_bol(local.Mdot, local.BH_Mass, -1) * local.Dt / C_LIGHT_CODE, mom_wt = 0;
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
            if(numngb < 0) return -1;
            for(n = 0; n < numngb; n++)
            {
                j = ngblist[n]; MyIDType OriginallyMarkedSwallowID; OriginallyMarkedSwallowID = P[j].SwallowID; // record this to help prevent double-counting below
                double dpos[3]={0},dvel[3]={0}; for(k=0;k<3;k++) {dpos[k]=P[j].Pos[k]-local.Pos[k]; dvel[k]=P[j].Vel[k]-local.Vel[k];}
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
                    for(k=0;k<3;k++) {dlv[k] *= wk * local.angmom_norm_topass_in_swallowloop; P[j].Vel[k]+=dlv[k]; SphP[j].VelPred[k]+=dlv[k]; P[j].dp[k]+=P[j].Mass*dlv[k]; out.accreted_momentum[k]-=P[j].Mass*dlv[k];}
                    out.accreted_J[0]-=P[j].Mass*(dpos[1]*dlv[2] - dpos[2]*dlv[1]); out.accreted_J[1]-=P[j].Mass*(dpos[2]*dlv[0] - dpos[0]*dlv[2]); out.accreted_J[2]-=P[j].Mass*(dpos[0]*dlv[1] - dpos[1]*dlv[0]);
                }
#endif
#if defined(BH_RETURN_BFLUX) // do a kernel-weighted redistribution of the magnetic flux in the sink into surrounding particles
                if((P[j].Type == 0) && (local.kernel_norm_topass_in_swallowloop > 0)){                    
                    double dB, b_fraction_toreturn = DMIN(0.1, local.Dt / (local.BH_Mass_AlphaDisk / local.Mdot)) * wk / local.kernel_norm_topass_in_swallowloop; // return a fraction dt/t_accretion of the total flux, with simple kernel weighting for each particle
                    for(k=0; k<3;k++) {
                        dB = b_fraction_toreturn * local.B[k];
                        SphP[j].B[k] +=  dB;
                        SphP[j].BPred[k] +=  dB;
                        out.accreted_B[k] -= dB;                        
                    }
                }
#endif                
                /* we've found a particle to be swallowed.  This could be a BH merger, DM particle, or baryon w/ feedback */
                if(P[j].SwallowID == local.ID && P[j].Mass > 0)
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
                    double mcount_for_conserve; mcount_for_conserve = f_accreted * P[j].Mass;
#if (BH_FOLLOW_ACCRETED_ANGMOM == 1) /* in this case we are only counting this if its coming from BH particles */
                    if(P[j].Type!=5) {mcount_for_conserve=0;} else {mcount_for_conserve=BPP(j).BH_Mass;}
#ifdef BH_ALPHADISK_ACCRETION
                    if(P[j].Type==5) {mcount_for_conserve += BPP(j).BH_Mass_AlphaDisk;}
#endif
#endif
#ifdef GRAIN_FLUID
                    if((1<<P[j].Type) & GRAIN_PTYPES) {out.accreted_dust_Mass += FLT(P[j].Mass);}
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
                    if(P[j].Type==5) {for(k=0;k<3;k++) {out.accreted_J[k] += FLT(mcount_for_conserve * BPP(j).BH_Specific_AngMom[k]);}}
#endif

                    if(P[j].Type == 5)  /* this is a BH-BH merger */
                    {
#ifdef BH_OUTPUT_MOREINFO
                        fprintf(FdBhMergerDetails,"%g  %u %g %2.7f %2.7f %2.7f  %u %g %2.7f %2.7f %2.7f\n", All.Time,  local.ID,local.BH_Mass,local.Pos[0],local.Pos[1],local.Pos[2],  P[j].ID,BPP(j).BH_Mass,P[j].Pos[0],P[j].Pos[1],P[j].Pos[2]);
#else
#ifndef IO_REDUCED_MODE
                        fprintf(FdBlackHolesDetails,"ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n", ThisTask, All.Time, local.ID, P[j].ID, local.BH_Mass, BPP(j).BH_Mass);
#endif
#endif
#ifdef BH_INCREASE_DYNAMIC_MASS
                        /* the true dynamical mass of the merging BH is P[j].Mass/BH_INCREASE_DYNAMIC_MASS unless exceeded by physical growth
                         - in the limit BPP(j).BH_Mass > BH_INCREASE_DYNAMIC_MASS x m_b, then bh_mass=P[j].Mass on average and we are good as well  */
                        out.accreted_Mass    += FLT( DMAX(BPP(j).BH_Mass, P[j].Mass/BH_INCREASE_DYNAMIC_MASS) );
#else
                        out.accreted_Mass    += FLT(P[j].Mass);
#endif
                        out.accreted_BH_Mass += FLT(BPP(j).BH_Mass);
#ifdef BH_ALPHADISK_ACCRETION
                        out.accreted_BH_Mass_alphadisk += FLT(BPP(j).BH_Mass_AlphaDisk);
#endif
#ifdef BH_COUNTPROGS
                        out.BH_CountProgs += BPP(j).BH_CountProgs;
#endif
                        bin = P[j].TimeBin; TimeBin_BH_mass[bin] -= BPP(j).BH_Mass; TimeBin_BH_dynamicalmass[bin] -= P[j].Mass; TimeBin_BH_Mdot[bin] -= BPP(j).BH_Mdot;
                        if(BPP(j).BH_Mass > 0) {TimeBin_BH_Medd[bin] -= BPP(j).BH_Mdot / BPP(j).BH_Mass;}
                        P[j].Mass = 0; BPP(j).BH_Mass = 0; BPP(j).BH_Mdot = 0; 
#ifdef GALSF
                        out.Accreted_Age = P[j].StellarAge;
#endif
                        N_BH_swallowed++;
                    } // if(P[j].Type == 5) -- BH + BH merger


#ifdef BH_GRAVCAPTURE_NONGAS /* DM and star particles can only be accreted ifdef BH_GRAVCAPTURE_NONGAS */
                    if((P[j].Type == 1) || (All.ComovingIntegrationOn && (P[j].Type==2||P[j].Type==3)) )
                    {   /* this is a DM particle: In this case, no kick, so just zero out the mass and 'get rid of' the particle (preferably by putting it somewhere irrelevant) */
#ifdef BH_OUTPUT_MOREINFO
                        printf(" ..BH_swallow_DM: j %d Type(j) %d  M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g \n", j,P[j].Type,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],local.Pos[0],local.Pos[1],local.Pos[2]);
#endif
                        out.accreted_Mass += FLT(P[j].Mass); out.accreted_BH_Mass += FLT(P[j].Mass); P[j].Mass = 0;
                        N_dm_swallowed++;
                    }
                    if((P[j].Type==4) || ((P[j].Type==2||P[j].Type==3) && !(All.ComovingIntegrationOn) ))
                    {   /* this is a star particle: If there is an alpha-disk, we let them go to the disk. If there is no alpha-disk, stars go to the BH directly and won't affect feedback. (Can be simply modified if we need something different.) */
                        out.accreted_Mass += FLT(P[j].Mass); out_accreted_BH_Mass_alphaornot += FLT(P[j].Mass); P[j].Mass = 0;
                        N_star_swallowed++;
                    }
#endif // #ifdef BH_GRAVCAPTURE_NONGAS -- BH + DM or Star merger


                    /* this is a gas particle: DAA: we need to see if the gas particle has to be accreted in full or not, depending on BH_WIND_KICK
                     the only difference with BH_ALPHADISK_ACCRETION should be that the mass goes first to the alphadisk */
                    if(P[j].Type == 0)                    
                    {
                        out.accreted_Mass += FLT(f_accreted*P[j].Mass);
#ifdef BH_GRAVCAPTURE_GAS
                        out_accreted_BH_Mass_alphaornot += FLT(f_accreted*P[j].Mass);
#endif
                        P[j].Mass *= (1-f_accreted);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue *= (1-f_accreted);
#endif
                        for(k=0; k<3; k++) {P[j].dp[k] *= (1-f_accreted);}

#ifdef BH_WIND_KICK     /* BAL kicking operations. NOTE: we have two separate BAL wind models, particle kicking and smooth wind model. This is where we do the particle kicking BAL model. This should also work when there is alpha-disk. */
                        double v_kick=All.BAL_v_outflow, dir[3]; for(k=0;k<3;k++) {dir[k]=dpos[k];} // DAA: default direction is radially outwards
#if defined(BH_COSMIC_RAYS) /* inject cosmic rays alongside wind injection */
                        double dEcr = All.BH_CosmicRay_Injection_Efficiency * P[j].Mass * (All.BAL_f_accretion/(1.-All.BAL_f_accretion)) * C_LIGHT_CODE*C_LIGHT_CODE;
                        inject_cosmic_rays(dEcr,All.BAL_v_outflow,5,j,dir);
#endif
#if (BH_WIND_KICK < 0)  /* DAA: along polar axis defined by angular momentum within Kernel (we could add finite opening angle) work out the geometry w/r to the plane of the disk */
                        if((dir[0]*J_dir[0] + dir[1]*J_dir[1] + dir[2]*J_dir[2]) > 0){for(k=0;k<3;k++) {dir[k]=J_dir[k];}} else {for(k=0;k<3;k++) {dir[k]=-J_dir[k];}}
#endif
                        for(k=0,norm=0;k<3;k++) {norm+=dir[k]*dir[k];} if(norm<=0) {dir[0]=0;dir[1]=0;dir[2]=1;norm=1;} else {norm=sqrt(norm); dir[0]/=norm;dir[1]/=norm;dir[2]/=norm;}
                        for(k=0;k<3;k++) {P[j].Vel[k]+=v_kick*All.cf_atime*dir[k]; SphP[j].VelPred[k]+=v_kick*All.cf_atime*dir[k];}
#ifdef GALSF_SUBGRID_WINDS // if sub-grid galactic winds are decoupled from the hydro, we decouple the BH kick winds as well
                        SphP[j].DelayTime = All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a;
#endif
#ifdef BH_OUTPUT_MOREINFO
                        printf(" ..BAL kick: P[j].ID %llu ID %llu Type(j) %d f_acc %g M(j) %g V(j).xyz %g/%g/%g P(j).xyz %g/%g/%g p(i).xyz %g/%g/%g v_out %g \n",(unsigned long long) P[j].ID, (unsigned long long) P[j].SwallowID,P[j].Type, All.BAL_f_accretion,P[j].Mass,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],local.Pos[0],local.Pos[1],local.Pos[2],v_kick);
                        fprintf(FdBhWindDetails,"%g  %llu %g  %2.7f %2.7f %2.7f  %2.7f %2.7f %2.7f %g %g %g %llu  %2.7f %2.7f %2.7f\n",All.Time, P[j].ID, P[j].Mass, P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],  P[j].Vel[0],P[j].Vel[1],P[j].Vel[2],dir[0]/norm,dir[1]/norm,dir[2]/norm, local.ID, local.Pos[0],local.Pos[1],local.Pos[2]);
#endif
#endif // #ifdef BH_WIND_KICK
                        N_gas_swallowed++;
#ifdef BH_OUTPUT_GASSWALLOW
                        MyDouble tempB[3]={0,0,0};
#ifdef MAGNETIC
                        for(k=0;k<3;k++) {tempB[k]=Get_Gas_BField(j,k);} //use particle magnetic field
#endif
                        fprintf(FdBhSwallowDetails,"%g %llu %g %2.7f %2.7f %2.7f %llu %g %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f %2.7f\n", All.Time, local.ID,local.Mass,local.Pos[0],local.Pos[1],local.Pos[2],  P[j].ID, P[j].Mass, (P[j].Pos[0]-local.Pos[0]),(P[j].Pos[1]-local.Pos[1]),(P[j].Pos[2]-local.Pos[2]), (P[j].Vel[0]-local.Vel[0]),(P[j].Vel[1]-local.Vel[1]),(P[j].Vel[2]-local.Vel[2]), SphP[j].InternalEnergy, tempB[0], tempB[1], tempB[2], SphP[j].Density);
#endif
                    }  // if(P[j].Type == 0)

                    P[j].SwallowID = 0; /* DAA: make sure it is not accreted (or ejected) by the same BH again if inactive in the next timestep */
                } // if(P[j].SwallowID == id)  -- particles being entirely or partially swallowed!!!

#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
                /* now, do any other feedback "kick" operations (which used the previous loops to calculate weights) */
                if(mom>0 && local.Mdot>0 && local.Dt>0 && OriginallyMarkedSwallowID==0 && P[j].SwallowID==0 && P[j].Mass>0 && P[j].Type==0) // particles NOT being swallowed!
                {
                    double r=0, dir[3]; for(k=0;k<3;k++) {dir[k]=dpos[k]; r+=dir[k]*dir[k];} // should be away from BH
                    if(r>0)
                    {
                        r=sqrt(r); for(k=0;k<3;k++) {dir[k]/=r;} /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                        for(norm=0,k=0;k<3;k++) {norm+=dir[k]*J_dir[k];}
                        mom_wt = bh_angleweight_localcoupling(j,norm,r,h_i) / local.BH_angle_weighted_kernel_sum;
                        if(local.BH_angle_weighted_kernel_sum<=0) mom_wt=0;
                                
#if defined(BH_COSMIC_RAYS) && defined(BH_WIND_CONTINUOUS) /* inject cosmic rays alongside continuous wind injection */
                        double dEcr = All.BH_CosmicRay_Injection_Efficiency * mom_wt * C_LIGHT_CODE*C_LIGHT_CODE * local.Mdot*local.Dt;
                        inject_cosmic_rays(dEcr,All.BAL_v_outflow,5,j,dir);
#endif
#if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK) /* inject BAL winds, this is the more standard smooth feedback model */
                        double m_wind = mom_wt * (1-All.BAL_f_accretion)/(All.BAL_f_accretion) * local.Mdot*local.Dt; /* mass to couple */
                        if(local.BH_angle_weighted_kernel_sum<=0) m_wind=0;
                        //1. check if (Vw-V0)*rhat <= 0   [ equivalently, check if   |Vw| <= V0*rhat ]
                        //2. if (1) is False, the wind will catch the particle, couple mass, momentum, energy, according to the equations above
                        //3. if (1) is True, the wind will not catch the particle, or will only asymptotically catch it. For the sake of mass conservation in the disk, I think it is easiest to treat this like the 'marginal' case where the wind barely catches the particle. In this case, add the mass normally, but no momentum, and no energy, giving:
                        //dm = m_wind, dV = 0, du = -mu*u0   [decrease the thermal energy slightly to account for adding more 'cold' material to it]
                        double dvr_gas_to_bh, dr_gas_to_bh;
                        for(dvr_gas_to_bh=dr_gas_to_bh=0, k=0;k<3;k++) {dvr_gas_to_bh += dvel[k]*dpos[k]; dr_gas_to_bh  += dpos[k]*dpos[k];}
                        dvr_gas_to_bh /= dr_gas_to_bh ;
                        
                        /* add wind mass to particle, correcting density as needed */
                        if(P[j].Hsml<=0)
                        {
                            if(SphP[j].Density>0) {SphP[j].Density*=(1+m_wind/P[j].Mass);} else {SphP[j].Density=m_wind*hinv3;}
                        } else {
                            SphP[j].Density += kernel_zero * m_wind/(P[j].Hsml*P[j].Hsml*P[j].Hsml);
                        }
                        P[j].Mass += m_wind;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                        SphP[j].MassTrue += m_wind;
#endif
                        /* now add wind momentum to particle */
                        if(dvr_gas_to_bh < All.BAL_v_outflow)   // gas moving away from BH at v < BAL speed
                        {
                            double e_wind = 0;
                            for(k=0;k<3;k++)
                            {
                                norm = All.cf_atime*All.BAL_v_outflow*dir[k] - dvel[k]; // relative wind-particle velocity (in code units) including BH-particle motion;
                                P[j].Vel[k] += All.BlackHoleFeedbackFactor * norm * m_wind/P[j].Mass; // momentum conservation gives updated velocity
                                SphP[j].VelPred[k] += All.BlackHoleFeedbackFactor * norm * m_wind/P[j].Mass;
                                e_wind += (norm/All.cf_atime)*(norm/All.cf_atime); // -specific- shocked wind energy
                            }
                            e_wind *= 0.5*m_wind/P[j].Mass; // make total wind energy, add to particle as specific energy of -particle-
                            SphP[j].InternalEnergy += e_wind; SphP[j].InternalEnergyPred += e_wind;
                        } else {    // gas moving away from BH at wind speed (or faster) already.
                            if(SphP[j].InternalEnergy * ( P[j].Mass - m_wind ) / P[j].Mass > 0) {SphP[j].InternalEnergy = SphP[j].InternalEnergy * ( P[j].Mass - m_wind ) / P[j].Mass;}
                        }
#endif // if defined(BH_WIND_CONTINUOUS) && !defined(BH_WIND_KICK)
                    } // r > 0
                } // (check if valid gas neighbor of interest)
#endif // defined(BH_CALC_LOCAL_ANGLEWEIGHTS)                
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







#endif // master flag
