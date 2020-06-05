#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"


/* Routines for pure thermal/scalar feedback/enrichment models: these are intended to represent
    extremely simplified models for stellar feedback manifest as a "pure thermal energy dump"
    (potentially with some cooling turnoff). This is -not- a model for mechanical feedback
    (which, critically, must include the actual momentum and solve for wind/SNe shock properties
    at the interface with the ISM). Those physics are included in the mechanical_fb.c file and
    algorithms therein. This also uses an extremely simple kernel-weighting, rather than a
    more self-consistent area weighting. */

/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#if defined(GALSF_FB_THERMAL)

/* routine that evaluates whether a FB event occurs in a given particle, in a given timestep */
void determine_where_addthermalFB_events_occur(void)
{
    int i; double check = 0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 4) continue;
        if(P[i].Mass <= 0) continue;
        check += mechanical_fb_calculate_eventrates(i,1); // this should do the calculation and add to number of SNe as needed //
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //
}

struct kernel_addthermalFB {double dp[3], r, wk, dwk, hinv, hinv3, hinv4;};


#define MASTER_FUNCTION_NAME addthermalFB_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int MASTER_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_addthermalFB    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_addthermalFB  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(addthermalFB_evaluate_active_check(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */


/* define structures to use below */
struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3], Hsml, Msne, Esne, wt_sum;
#ifdef METALS
    MyDouble yields[NUM_METAL_SPECIES];
#endif
    int NodeList[NODELISTLENGTH];
}
*DATAIN_NAME, *DATAGET_NAME;

/* define properties to be injected. these must be scalar-only -- the simple routine below will not conserve vector inputs/ejecta (e.g. momentum) */
void particle2in_addthermalFB(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    if((P[i].SNe_ThisTimeStep<=0)||(P[i].DensAroundStar<=0)) {in->Msne=0; return;} // trap for no sne
    int k; in->Hsml=PPP[i].Hsml; in->wt_sum=P[i].DensAroundStar; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];} // simple kernel-weighted deposition
    struct addFB_evaluate_data_in_ local; particle2in_addFB_fromstars(&local,i,0); // get feedback properties from generic routine //
    in->Msne = local.Msne; in->Esne = 0.5 * local.Msne * local.SNe_v_ejecta*local.SNe_v_ejecta; // assign mass and energy to be used below
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=local.yields[k];} // assign yields //
#endif
}


struct OUTPUT_STRUCT_NAME
{
    MyFloat M_coupled;
}
*DATARESULT_NAME, *DATAOUT_NAME;

void out2particle_addthermalFB(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    P[i].Mass -= out->M_coupled;
    if((P[i].Mass<0)||(isnan(P[i].Mass))) {P[i].Mass=0;}
}


int addthermalFB_evaluate_active_check(int i);
int addthermalFB_evaluate_active_check(int i)
{
    if(P[i].Type != 4) return 0;
    if(P[i].Mass <= 0) return 0;
    if(PPP[i].Hsml <= 0) return 0;
    if(PPP[i].NumNgb <= 0) return 0;
    if(P[i].SNe_ThisTimeStep>0) {return 1;}
    return 0;
}


int addthermalFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double u,r2,h2,kernel_zero,wk;
    struct kernel_addthermalFB kernel;
    struct INPUT_STRUCT_NAME local;
    struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1); wk=0;
    
    /* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_addthermalFB(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(local.Msne<=0) return 0; // no SNe for the master particle! nothing to do here //
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml;
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    
    /* Now start the actual FB computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = DATAGET_NAME[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) return -1;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) continue; // same particle //
                if(r2>=h2) continue; // outside kernel //
                // calculate kernel quantities //
                kernel.r = sqrt(r2);
                if(kernel.r <= 0) continue;
                u = kernel.r * kernel.hinv;
                if(u<1) {kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);} else {kernel.wk=kernel.dwk=0;}
                if((kernel.wk <= 0)||(isnan(kernel.wk))) continue;
                wk = P[j].Mass * kernel.wk / local.wt_sum; // normalized weight function
                
                /* inject mass */
                double dM_ejecta_in = wk * local.Msne;
                if(P[j].Hsml<=0) {if(SphP[j].Density>0){SphP[j].Density*=(1+dM_ejecta_in/P[j].Mass);} else {SphP[j].Density=dM_ejecta_in*kernel.hinv3;}} else {SphP[j].Density+=kernel_zero*dM_ejecta_in/(P[j].Hsml*P[j].Hsml*P[j].Hsml);}
                SphP[j].Density *= 1 + dM_ejecta_in/P[j].Mass; // inject mass at constant particle volume //
                P[j].Mass += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[j].MassTrue += dM_ejecta_in;
#endif
#ifdef METALS
                /* inject metals */
                for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k]=(1-dM_ejecta_in/P[j].Mass)*P[j].Metallicity[k] + dM_ejecta_in/P[j].Mass*local.yields[k];}
#endif
                /* inject energy */
                SphP[j].InternalEnergy += wk * local.Esne / P[j].Mass;
                SphP[j].InternalEnergyPred += wk * local.Esne / P[j].Mass;
#ifdef GALSF_FB_TURNOFF_COOLING
                /* if the sub-grid 'cooling turnoff' model is enabled, turn off cooling for the 'blastwave timescale',
                 which is physically the timescale for the blastwave to be completely stopped by ISM ram-pressure
                 (much longer than the actual cooling time of the blastwave) */
                double Esne51 = local.Esne * UNIT_ENERGY_IN_CGS/1.e51;
                double density_to_n = All.cf_a3inv*UNIT_DENSITY_IN_NHCGS;
                double pressure_to_p4 = (1/All.cf_afac1)*density_to_n*U_TO_TEMP_UNITS / 1.0e4; 
                double dt_ram = 7.08 * pow(Esne51*SphP[j].Density*density_to_n,0.34) * pow(SphP[j].Pressure*pressure_to_p4,-0.70) / (UNIT_TIME_IN_MYR);
                if(dt_ram > SphP[j].DelayTimeCoolingSNe) SphP[j].DelayTimeCoolingSNe = dt_ram;
#endif
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
{
                startnode = DATAGET_NAME[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}    /* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addthermalFB(&out, target, 0, 0);} else {DATARESULT_NAME[target] = out;}
    return 0;
} // int addthermalFB_evaluate



void thermal_fb_calc(void)
{
    PRINT_STATUS(" ..depositing thermal feedback to gas");
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    CPU_Step[CPU_SNIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */




#endif /* GALSF_FB_THERMAL */

