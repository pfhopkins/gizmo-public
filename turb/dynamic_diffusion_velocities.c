#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"


/*! \file dynamic_diffusion_velocities.c
 *  \brief need filtered velocity information to calculate filtered gradients
 *
 */
/*
 * This file was rewritten by Doug Rennehan (douglas.rennehan@gmail.com) for GIZMO, and was
 * copied with modifications from gradients.c, which was written by Phil Hopkins
 * (phopkins@caltech.edu) for GIZMO.
 */

#ifdef TURB_DIFF_DYNAMIC

#define ASSIGN_ADD_PRESET(x,y,mode) (x+=y)
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#define SHOULD_I_USE_SPH_GRADIENTS(condition_number) ((condition_number > CONDITION_NUMBER_DANGER) ? (1):(0))
#define NV_MYSIGN(x) (( x > 0 ) - ( x < 0 ))

struct kernel_DiffFilter {
    double dp[3], r, wk_i, wk_j, dwk_i, dwk_j, h_i;
};


#define CORE_FUNCTION_NAME DiffFilter_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_DiffFilter    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_DiffFilter  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if((P[i].Type==0)&&(P[i].TimeBin>=0)&&(P[i].Mass>0)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

struct INPUT_STRUCT_NAME {
    MyDouble Pos[3];
    MyFloat Mass;
    MyFloat Hsml;
    MyDouble Density;
#ifndef DONOTUSENODELIST
    int NodeList[NODELISTLENGTH];
#endif
#ifdef GALSF_SUBGRID_WINDS
    MyFloat DelayTime;
#endif
    MyDouble VelPred[3];
}
*DATAIN_NAME, *DATAGET_NAME;

static inline void particle2in_DiffFilter(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration) {
    int k; for(k = 0; k < 3; k++) {in->Pos[k] = P[i].Pos[k]; in->VelPred[k] = SphP[i].VelPred[k];}
    in->Density = SphP[i].Density;
    in->Hsml = PPP[i].Hsml;
    in->Mass = DMAX(0,P[i].Mass);
#ifdef GALSF_SUBGRID_WINDS
    in->DelayTime = SphP[i].DelayTime;
#endif
}


struct OUTPUT_STRUCT_NAME {
    MyDouble Norm_hat;
    MyDouble Velocity_bar[3];
    MyFloat FilterWidth_bar;
    MyFloat MaxDistance_for_grad;
}
*DATARESULT_NAME, *DATAOUT_NAME;

#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))
static inline void out2particle_DiffFilter(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration) {
    int k; for (k = 0; k < 3; k++) {ASSIGN_ADD_PRESET(SphP[i].Velocity_bar[k], out->Velocity_bar[k], mode);}
    MAX_ADD(SphP[i].FilterWidth_bar, out->FilterWidth_bar, mode);
    MAX_ADD(SphP[i].MaxDistance_for_grad, out->MaxDistance_for_grad, mode);
    ASSIGN_ADD_PRESET(SphP[i].Norm_hat, out->Norm_hat, mode);
}

/* operations that need to be performed before entering the main loop */
void dynamic_diff_vel_calc_initial_operations_preloop(void);
void dynamic_diff_vel_calc_initial_operations_preloop(void)
{
    /* Because of the smoothing operation, need to set bar quantity to current SPH value first */
    int i;
    for (i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if (P[i].Type == 0) {
            SphP[i].Norm_hat = 0;
            SphP[i].h_turb = Get_Particle_Size(i); // All.cf_atime unnecessary, will multiply later
            SphP[i].FilterWidth_bar = 0;
            SphP[i].MaxDistance_for_grad = 0;
            int k; for (k = 0; k < 3; k++) {SphP[i].Velocity_bar[k] = SphP[i].VelPred[k] / All.TurbDynamicDiffSmoothing;}
        }
    }
}


/**
 *
 *  Computes the smoothed velocity field Velocity_bar according to eq. 2.17 in
 *  Monaghan 2011 (turbulence for SPH).
 *  - D. Rennehan
 *
 */
/*!   -- this subroutine contains no writes to shared memory -- */
int DiffFilter_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration) {
    /* initialize and check if we should bother doing a neighbor loop */
    int startnode, numngb, listindex = 0, j, k, n;
    double hinv, hinv3, hinv4, r2, u;
    struct kernel_DiffFilter kernel;
    struct INPUT_STRUCT_NAME local;
    struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    memset(&kernel, 0, sizeof(struct kernel_DiffFilter));
    if (mode == 0) {particle2in_DiffFilter(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if (local.Hsml <= 0) {return 0;}
    if (local.Mass == 0) {return 0;}
    if (local.Density <= 0) {return 0;}

    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.Hsml;
    double h2_i = kernel.h_i * kernel.h_i;
    int kernel_mode_i = 0;
    kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);

    /* Now start the actual neighbor computation for this particle */
    if (mode == 0) {startnode = All.MaxPart;}	/* root node */
        else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;}	/* open it */

    while (startnode >= 0) {
        while (startnode >= 0) {
            numngb = ngb_treefind_pairs_threads(local.Pos, All.TurbDynamicDiffFac * kernel.h_i, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb < 0) {return -2;}

            for (n = 0; n < numngb; n++) {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
#ifdef GALSF_SUBGRID_WINDS
                if (local.DelayTime == 0 && SphP[j].DelayTime > 0) {continue;}
                if (local.DelayTime > 0 && SphP[j].DelayTime == 0) {continue;}
#endif
                if (P[j].Mass <= 0) {continue;}
                if (SphP[j].Density <= 0) {continue;}

                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); /*  now find the closest image in the given box size  */
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                double mean_weight = 0.5 * (local.Density + SphP[j].Density) / (local.Density * SphP[j].Density);
                double h_j = PPP[j].Hsml;
                double V_j = P[j].Mass * mean_weight;
                if (r2 <= 0) {continue;}

                double h_avg = 0.5 * (kernel.h_i + h_j);
                double hhat = All.TurbDynamicDiffFac * kernel.h_i;
                double hhat_j = All.TurbDynamicDiffFac * h_j;
                double hhat2 = hhat * hhat;
                double hhatj2 = hhat_j * hhat_j;
                if ((r2 >= hhat2) && (r2 >= hhatj2)) {continue;}

                double hhatinv, hhatinv3, hhatinv4, uhat, wkhat, dwkhat, rhat;
                rhat = sqrt(r2);
                if (rhat < hhat) {
                    kernel_hinv(hhat, &hhatinv, &hhatinv3, &hhatinv4);
                    uhat = DMIN(rhat * hhatinv, 1.0);
                    kernel_main(uhat, hhatinv3, hhatinv4, &wkhat, &dwkhat, kernel_mode_i);
                } else {wkhat = dwkhat = 0;}
                if (rhat < hhat) {out.Norm_hat += P[j].Mass * wkhat;}
                if ((r2 >= h2_i) && (r2 >= (h_j * h_j))) {continue;}
                kernel.r = sqrt(r2);
                if (kernel.r > out.FilterWidth_bar) {out.FilterWidth_bar = kernel.r;}

                // This is very, very important for supersonic flows, or any flow with highly varying smoothing lengths. The FilterWidth_bar (h_bar) *must* extend out to the maximum interaction distance.
                if (kernel.r > out.MaxDistance_for_grad) {out.MaxDistance_for_grad = kernel.r;}
                kernel_hinv(h_avg, &hinv, &hinv3, &hinv4);
                u = DMIN(kernel.r * hinv, 1.0);
                if(u<1) {kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, kernel_mode_i);} else {kernel.wk_i=kernel.dwk_i=0;}
                double weight_i = kernel.wk_i * V_j;
                double VelPred_diff[3];
                /* Because we are using the average h value between i,j: W_ij = W_ji */
                if (kernel.r < h_avg) {
                    for (k = 0; k < 3; k++) {
                        VelPred_diff[k] = SphP[j].VelPred[k] - local.VelPred[k];
                        out.Velocity_bar[k] += VelPred_diff[k] * weight_i;
                    }
                }
            } // numngb loop
        } // while(startnode)

#ifndef DONOTUSENODELIST
        if (mode == 1) {listindex++;
            if (listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if (startnode >= 0) startnode = Nodes[startnode].u.d.nextnode;}}	/* open it */
#endif
    }
    if(mode == 0) {out2particle_DiffFilter(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* Now collect the result at the right place */
    return 0;
}


/**
 * primary routine being called for this calculation
 */
void dynamic_diff_vel_calc(void) {
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    PRINT_STATUS("Start velocity smoothing computation...");
    dynamic_diff_vel_calc_initial_operations_preloop(); /* any initial operations */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    PRINT_STATUS(" ..velocity smoothing done.");
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_IMPROVDIFFCOMPUTE] += timecomp; CPU_Step[CPU_IMPROVDIFFWAIT] += timewait; CPU_Step[CPU_IMPROVDIFFCOMM] += timecomm;
    CPU_Step[CPU_IMPROVDIFFMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */


#endif /* End TURB_DIFF_DYNAMIC */
