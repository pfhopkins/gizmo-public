/*! \file blackhole_feed.c
*  \brief This is where particles are marked for gas accretion.
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


#define CORE_FUNCTION_NAME blackhole_feed_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(bhsink_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */


/* this structure defines the variables that need to be sent -from- the 'searching' element */
struct INPUT_STRUCT_NAME
{
    int NodeList[NODELISTLENGTH]; MyDouble Pos[3]; MyFloat Vel[3], Hsml, Mass, BH_Mass, Dt, Density, Mdot; MyIDType ID;
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    MyFloat Jgas_in_Kernel[3];
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    MyFloat mass_to_swallow_edd;
#endif
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
    MyFloat SinkRadius;
#endif
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
    MyFloat AGS_Hsml;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    MyFloat BH_Mass_AlphaDisk;
#endif
#ifdef BH_ACCRETE_NEARESTFIRST
    MyFloat BH_dr_to_NearestGasNeighbor;
#endif
}
*DATAIN_NAME, *DATAGET_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k, j_tempinfo; j_tempinfo=P[i].IndexMapToTempStruc; /* link to the location in the shared structure where this is stored */
    for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k]; in->Vel[k]=P[i].Vel[k];} /* good example - always needed */
    in->Hsml = PPP[i].Hsml; in->Mass = P[i].Mass; in->BH_Mass = BPP(i).BH_Mass; in->ID = P[i].ID; in->Density = BPP(i).DensAroundStar; in->Mdot = BPP(i).BH_Mdot;
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
    in->SinkRadius = PPP[i].SinkRadius;
#endif
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
    in->AGS_Hsml = PPP[i].AGS_Hsml;
#endif
#ifdef BH_ALPHADISK_ACCRETION
    in->BH_Mass_AlphaDisk = BPP(i).BH_Mass_AlphaDisk;
#endif
    in->Dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
    if(P[i].Type == 5){in->Dt = P[i].dt_since_last_gas_search;}
#endif
#ifdef BH_ACCRETE_NEARESTFIRST
    in->BH_dr_to_NearestGasNeighbor = P[i].BH_dr_to_NearestGasNeighbor;
#endif
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
    for(k=0;k<3;k++) {in->Jgas_in_Kernel[k] = P[i].BH_Specific_AngMom[k];}
#else
    for(k=0;k<3;k++) {in->Jgas_in_Kernel[k] = BlackholeTempInfo[j_tempinfo].Jgas_in_Kernel[k];}
#endif
#endif
#if defined(BH_GRAVCAPTURE_GAS)
    in->mass_to_swallow_edd = BlackholeTempInfo[j_tempinfo].mass_to_swallow_edd;
#endif
}


/* this structure defines the variables that need to be sent -back to- the 'searching' element */
struct OUTPUT_STRUCT_NAME
{ /* define variables below as e.g. "double X;" */
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    double BH_angle_weighted_kernel_sum;
#endif
#ifdef BH_REPOSITION_ON_POTMIN
    double BH_MinPot, BH_MinPotPos[3];
#endif
}
*DATARESULT_NAME, *DATAOUT_NAME; /* dont mess with these names, they get filled-in by your definitions automatically */

#define ASSIGN_ADD_PRESET(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    int k, target; k=0; target = P[i].IndexMapToTempStruc;
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    ASSIGN_ADD_PRESET(BlackholeTempInfo[target].BH_angle_weighted_kernel_sum, out->BH_angle_weighted_kernel_sum, mode);
#endif
#ifdef BH_REPOSITION_ON_POTMIN
    if(mode==0) {BPP(i).BH_MinPot=out->BH_MinPot; for(k=0;k<3;k++) {BPP(i).BH_MinPotPos[k]=out->BH_MinPotPos[k];}
        } else {if(out->BH_MinPot < BPP(i).BH_MinPot) {BPP(i).BH_MinPot=out->BH_MinPot; for(k=0;k<3;k++) {BPP(i).BH_MinPotPos[k]=out->BH_MinPotPos[k];}}}
#endif
}



/* do loop over neighbors to get quantities for accretion */
/*!   -- this subroutine writes to shared memory [updating the neighbor values]: need to protect these writes for openmp below. none of the modified values are read, so only the write block is protected. */
int blackhole_feed_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration);
int blackhole_feed_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    /* initialize variables before loop is started */
    int startnode, numngb, listindex = 0, j, k, n; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* define variables and zero memory and import data for local target*/
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* imports the data to the correct place and names */
    double h_i = local.Hsml, wk, dwk, vrel, vesc, dpos[3], dvel[3], f_accreted; f_accreted=1;
    if((local.Mass<0)||(h_i<=0)) {return 0;}
    double w, p, r2, r, u, sink_radius=All.ForceSoftening[5], h_i2 = h_i * h_i, hinv = 1 / h_i, hinv3 = hinv * hinv * hinv, ags_h_i = All.ForceSoftening[5]; p=0; w=0;
#ifdef BH_REPOSITION_ON_POTMIN
    out.BH_MinPot = BHPOTVALUEINIT;
#endif
#if (ADAPTIVE_GRAVSOFT_FORALL & 32)
    ags_h_i = local.AGS_Hsml;
#endif
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    double J_dir[3]; for(k=0;k<3;k++) {J_dir[k] = local.Jgas_in_Kernel[k];}
#endif
#if defined(BH_GRAVCAPTURE_GAS) && defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
    double meddington = bh_eddington_mdot(local.BH_Mass), medd_max_accretable = All.BlackHoleEddingtonFactor * meddington * local.Dt, eddington_factor = local.mass_to_swallow_edd / medd_max_accretable;   /* if <1 no problem, if >1, need to not set some swallowIDs */
#endif
#if defined(BH_SWALLOWGAS)
    double mass_markedswallow,bh_mass_withdisk; mass_markedswallow=0; bh_mass_withdisk=local.BH_Mass;
#ifdef BH_ALPHADISK_ACCRETION
    bh_mass_withdisk += local.BH_Mass_AlphaDisk;
#endif
#endif
#if defined(BH_WIND_KICK) && !defined(BH_GRAVCAPTURE_GAS) /* DAA: increase the effective mass-loading of BAL winds to reach the desired momentum flux given the outflow velocity "All.BAL_v_outflow" chosen --> appropriate for cosmological simulations where particles are effectively kicked from ~kpc scales (i.e. we need lower velocity and higher mass outflow rates compared to accretion disk scales) - */
    f_accreted = All.BAL_f_accretion; if((All.BlackHoleFeedbackFactor > 0) && (All.BlackHoleFeedbackFactor != 1.)) {f_accreted /= All.BlackHoleFeedbackFactor;} else {if(All.BAL_v_outflow > 0) {f_accreted = 1./(1. + fabs(1.*BH_WIND_KICK)*All.BlackHoleRadiativeEfficiency*C_LIGHT_CODE/(All.BAL_v_outflow));}}
#endif
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS)
    double norm=0; for(k=0;k<3;k++) {norm+=J_dir[k]*J_dir[k];}
    if(norm>0) {norm=1/sqrt(norm); for(k=0;k<3;k++) {J_dir[k]*=norm;}} else {J_dir[0]=J_dir[1]=0; J_dir[2]=1;}
#endif
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
    sink_radius = local.SinkRadius;
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
                if(P[j].Mass > 0)
                {
                    for(k=0;k<3;k++) {dpos[k] = P[j].Pos[k] - local.Pos[k]; dvel[k]=P[j].Vel[k]-local.Vel[k];}
                    NEAREST_XYZ(dpos[0],dpos[1],dpos[2],-1); r2=0; for(k=0;k<3;k++) {r2 += dpos[k]*dpos[k];}
                    NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,P[j].Pos,dvel,-1); /* wrap velocities for shearing boxes if needed */
                    if(r2 < h_i2 || r2 < PPP[j].Hsml*PPP[j].Hsml)
                    {
                        vrel=0; for(k=0;k<3;k++) {vrel += dvel[k]*dvel[k];}
                        r=sqrt(r2); vrel=sqrt(vrel)/All.cf_atime;  /* relative velocity in physical units. do this once and use below */
                        vesc=bh_vesc(j, local.Mass, r, ags_h_i);
                        
                        /* note that SwallowID is both read and potentially re-written below: we need to make sure this is done in a thread-safe manner */
                        MyIDType SwallowID_j;
                        #pragma omp atomic read
                        SwallowID_j = P[j].SwallowID; // ok got a clean read. -not- gauranteed two threads won't see this at the same time and compete over it [both think they get it here]. but only one will -actually- get it, and that's ok.
                        
#ifdef BH_REPOSITION_ON_POTMIN
                        /* check if we've found a new potential minimum which is not moving too fast to 'jump' to */
                        double boundedness_function, potential_function; boundedness_function = P[j].Potential + 0.5 * vrel*vrel * All.cf_atime; potential_function = P[j].Potential;
#if (BH_REPOSITION_ON_POTMIN == 2)
                        if( boundedness_function < 0 )
                        {
                            double wt_rsoft = r / (3.*All.ForceSoftening[5]); // normalization arbitrary here, just using for convenience for function below
                            boundedness_function *= 1./(1. + wt_rsoft*wt_rsoft); // this down-weights particles which are very far away, relative to the user-defined force softening scale, which should define some 'confidence radius' of resolution around the BH particle
                        }
                        potential_function = boundedness_function; // jumps based on -most bound- particle, not just deepest potential (down-weights fast-movers)
#endif
                        if(potential_function < out.BH_MinPot)
#if (BH_REPOSITION_ON_POTMIN == 1)
                        if( P[j].Type == 4 && vrel <= vesc )   // DAA: only if it is a star particle & bound
#endif
#if (BH_REPOSITION_ON_POTMIN == 2)
                        if( (P[j].Type != 0) && (P[j].Type != 5) )   // allow stars or dark matter but exclude gas, it's too messy! also exclude BHs, since we don't want to over-merge them
#endif
                        {
                            out.BH_MinPot=potential_function; for(k=0;k<3;k++) {out.BH_MinPotPos[k] = P[j].Pos[k];}
                        }
#endif // BH_REPOSITION_ON_POTMIN
			
                        
                        
                        /* check_for_bh_merger.  Easy.  No Edd limit, just a pos and vel criteria. */
#if !defined(BH_DEBUG_DISABLE_MERGERS)
                        if(P[j].Type == 5)  /* we may have a black hole merger -- check below if allowed */
                            if((local.ID != P[j].ID) && (SwallowID_j == 0) && (BPP(j).BH_Mass < local.BH_Mass)) /* we'll assume most massive BH swallows the other - simplifies analysis and ensures unique results */
#ifdef SINGLE_STAR_SINK_DYNAMICS
                            if((r < 1.0001*P[j].min_dist_to_bh) && (r < PPP[j].Hsml) && (r < sink_radius) && (P[j].Mass < local.Mass) && (P[j].Mass < 10*All.MeanGasParticleMass)) /* only merge away stuff that is within the softening radius, and is no more massive that a few gas particles */
#endif
                            {
                                if((vrel < vesc) && (bh_check_boundedness(j,vrel,vesc,r,sink_radius)==1))
                                {
                                    printf(" ..BH-BH Merger: P[j.]ID=%llu to be swallowed by id=%llu \n", (unsigned long long) P[j].ID, (unsigned long long) local.ID);
                                    SwallowID_j = local.ID;
                                } else {
#if defined(BH_OUTPUT_MOREINFO)     // DAA: BH merger info will be saved in a separate output file
                                    printf(" ..ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n", ThisTask, All.Time, local.ID, P[j].ID, vrel, vesc);
#elif !defined(IO_REDUCED_MODE)
                                    fprintf(FdBlackHolesDetails, "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g vesc=%g\n", ThisTask, All.Time, local.ID, P[j].ID, vrel, vesc); fflush(FdBlackHolesDetails);
#endif
                                }
                            } // if eligible for bh-bh mergers //
#endif // BH_DEBUG_DISABLE_MERGERS
                        
                        
                        
                        /* This is a similar loop to what we already did in blackhole_environment, but here we stochastically
                         reduce GRAVCAPT events in order to (statistically) obey the eddington limit */
#if defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVCAPTURE_NONGAS)
                        if((P[j].Type != 5) && (SwallowID_j < local.ID)) // we have a particle not already marked to swallow
                        {
#ifdef SINGLE_STAR_SINK_DYNAMICS
                            double eps = DMAX( r , DMAX(P[j].Hsml , ags_h_i) * KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER); // plummer-equivalent
			                if(eps*eps*eps /(P[j].Mass + local.Mass) <= P[j].SwallowTime)
#endif
#if defined(BH_ALPHADISK_ACCRETION)
                            if(local.BH_Mass_AlphaDisk < BH_ALPHADISK_ACCRETION*local.BH_Mass)
#endif
#if defined(BH_ACCRETE_NEARESTFIRST)
                            if((P[j].Type != 0) || (r<=1.0001*local.BH_dr_to_NearestGasNeighbor))
#endif
                            if((vrel < vesc)) // && (particles_swallowed_this_bh_this_process < particles_swallowed_this_bh_this_process_max))
                            { /* bound */
#ifdef BH_GRAVCAPTURE_FIXEDSINKRADIUS
                                double spec_mom=0; for(k=0;k<3;k++) {spec_mom += dvel[k]*dpos[k];} // delta_x.delta_v
                                spec_mom = (r2*vrel*vrel - spec_mom*spec_mom*All.cf_a2inv); // specific angular momentum^2 = r^2(delta_v)^2 - (delta_v.delta_x)^2;
				                if(spec_mom < All.G * (local.Mass + P[j].Mass) * sink_radius)  // check Bate 1995 angular momentum criterion (in addition to bounded-ness)
#endif
                                if( bh_check_boundedness(j,vrel,vesc,r,sink_radius)==1 ) /* bound and apocenter within target distance */
                                {
#ifdef BH_GRAVCAPTURE_NONGAS        /* simply swallow non-gas particle if BH_GRAVCAPTURE_NONGAS enabled */
                                    if(P[j].Type != 0) {SwallowID_j = local.ID;}
#endif
#if defined(BH_GRAVCAPTURE_GAS)     /* now deal with gas */
                                    if(P[j].Type == 0)
                                    {
#if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION) /* if Eddington-limited and NO alpha-disk, do this stochastically */
                                        p = 1. / eddington_factor;
#if defined(BH_WIND_CONTINUOUS) || defined(BH_WIND_KICK)
                                        p /= All.BAL_f_accretion; // we need to accrete more, then remove the mass in winds
#endif
                                        w = get_random_number(P[j].ID);
                                        if(w < p)
                                        {
#ifdef BH_OUTPUT_MOREINFO
                                            printf(" ..BH-Food Marked: P[j.]ID=%llu to be swallowed by id=%llu \n", (unsigned long long) P[j].ID, (unsigned long long) local.ID);
#endif
                                            SwallowID_j = local.ID;
                                        }
#else //if defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                        SwallowID_j = local.ID; /* in other cases, just swallow the particle */  //particles_swallowed_this_bh_this_process++;
#endif //else defined(BH_ENFORCE_EDDINGTON_LIMIT) && !defined(BH_ALPHADISK_ACCRETION)
                                    } //if (P[j].Type == 0)
#endif //ifdef BH_GRAVCAPTURE_GAS
                                } // if( apocenter in tolerance range )
                            } // if(vrel < vesc)
                        } //if(P[j].Type != 5)
#endif // defined(BH_GRAVCAPTURE_GAS) || defined(BH_GRAVCAPTURE_NONGAS)
                        
                        
                        
                        /* now is the more standard accretion only of gas, according to the mdot calculated before */
                        if(P[j].Type == 0) /* here we have a gas particle */
                        {
                            u=r*hinv; if(u<1) {kernel_main(u,hinv3,hinv*hinv3,&wk,&dwk,-1);} else {wk=dwk=0;}
#if defined(BH_SWALLOWGAS) && !defined(BH_GRAVCAPTURE_GAS) /* compute accretion probability, this below is only meaningful if !defined(BH_GRAVCAPTURE_GAS)... */
                            if(SwallowID_j < local.ID)
                            {
                                double dm_toacc = bh_mass_withdisk - (local.Mass + mass_markedswallow); if(dm_toacc>0) {p=dm_toacc*wk/local.Density;} else {p=0;}
#ifdef BH_WIND_KICK /* DAA: for stochastic winds (BH_WIND_KICK) we remove a fraction of mass from gas particles prior to kicking --> need to increase the probability here to balance black hole growth */
                                if(f_accreted>0) {p /= f_accreted; if((bh_mass_withdisk - local.Mass) < 0) {p = ( (1-f_accreted)/f_accreted ) * local.Mdot * local.Dt * wk / local.Density;}} /* DAA: compute outflow probability when "bh_mass_withdisk < mass" - we don't need to enforce mass conservation in this case, relevant only in low-res sims where the BH seed mass is much lower than the gas particle mass */
#endif
#ifdef BH_ACCRETE_NEARESTFIRST /* put all the weight on the single nearest gas particle, instead of spreading it in a kernel-weighted fashion */
                                p=0; if(dm_toacc>0 && P[j].Mass>0 && r<1.0001*local.BH_dr_to_NearestGasNeighbor) {p=dm_toacc/P[j].Mass;}
#endif
                                w = get_random_number(P[j].ID);
                                if(w < p)
                                {
#ifdef BH_OUTPUT_MOREINFO
                                    printf(" ..BH-Food Marked: j %d w %g p %g TO_BE_SWALLOWED \n",j,w,p);
#endif
                                    SwallowID_j = local.ID;
                                    mass_markedswallow += P[j].Mass*f_accreted;
                                } // if(w < p)
                            } // swallowID < localID
#endif // BH_SWALLOWGAS
#if defined(BH_CALC_LOCAL_ANGLEWEIGHTS) /* calculate the angle-weighting for the photon momentum */
                            if((local.Dt>0)&&(r>0)&&(SwallowID_j==0)&&(P[j].Mass>0)&&(P[j].Type==0))
                            { /* cos_theta with respect to disk of BH is given by dot product of r and Jgas */
                                norm=0; for(k=0;k<3;k++) {norm+=(dpos[k]/r)*J_dir[k];}
                                out.BH_angle_weighted_kernel_sum += bh_angleweight_localcoupling(j,norm,r,h_i);
                            }
#endif
#ifdef BH_THERMALFEEDBACK
                            double energy = bh_lum_bol(local.Mdot, local.BH_Mass, -1) * local.Dt;
                            if(local.Density > 0) {
                                #pragma omp atomic
                                SphP[j].Injected_BH_Energy += (wk/local.Density) * energy * P[j].Mass;
                            }
#endif                            
                        } // if(P[j].Type == 0)
                        
                        
                        /* ok, before exiting this loop, need to mark whether or not we actually designated a particle for accretion! */
                        if(SwallowID_j > 0)
                        {
                            #pragma omp atomic write
                            P[j].SwallowID = SwallowID_j;  // ok got a clean write. -not- gauranteed two threads won't see this at the same time and compete over it [both think they get it here]. but only one will -actually- get it, and that's ok.
                        }
                        
                    } // if(r2 < h_i2)
                } // if(P[j].Mass > 0)
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    } // while(startnode >= 0) (outer of the double-loop)
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
} /* closes bh_evaluate routine */



void blackhole_feed_loop(void)
{
#include "../../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
#include "../../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
#include "../../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
CPU_Step[CPU_BLACKHOLES] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */


#endif // top-level flag
