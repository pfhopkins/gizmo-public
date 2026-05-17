#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "difffilter_loop_api.h"   /* difffilter_vel_calc_gpu_toplevel */


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

/* dynamic_diff_vel_calc() below drives the DiffFilterSpec runner loop via
   difffilter_vel_calc_gpu_toplevel() (turb/difffilter_loop.cc). Legacy
   CPU-tree scaffolding (code_block_xchange_initialize.h, INPUT_STRUCT_NAME,
   OUTPUT_STRUCT_NAME, particle2in_DiffFilter, out2particle_DiffFilter,
   DiffFilter_evaluate) was retired in Step 5 Phase D1; the legacy GPU shim
   difffilter_gpu.cc was retired by DiffFilterSpec (Wave 5). */

/* operations that need to be performed before entering the main loop */
void dynamic_diff_vel_calc_initial_operations_preloop(void);
void dynamic_diff_vel_calc_initial_operations_preloop(void)
{
    /* Because of the smoothing operation, need to set bar quantity to current fluid value first */
    for (int i : ActiveParticleList) {
        if (P[i].Type == 0) {
            CellP[i].Norm_hat = 0;
            CellP[i].h_turb = P[i].Get_Particle_Size(); // All.cf_atime unnecessary, will multiply later
            CellP[i].FilterWidth_bar = 0;
            CellP[i].MaxDistance_for_grad = 0;
            CellP[i].Velocity_bar = CellP[i].VelPred / All.TurbDynamicDiffSmoothing;
        }
    }
}


/**
 * primary routine being called for this calculation
 */
void dynamic_diff_vel_calc(void) {
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    PRINT_STATUS("Start velocity smoothing computation...");
    dynamic_diff_vel_calc_initial_operations_preloop(); /* any initial operations */

    /* Runner-template path (DiffFilterSpec): scaled-symmetric gas-gas
       velocity-smoothing loop. The toplevel builds the active gas list,
       drives run_neighbor_loop, and scatters the per-active result straight
       into CellP[] via apply_active_writeback. Replaces the legacy
       difffilter_evaluate_gpu() call + host scatter. */
    difffilter_vel_calc_gpu_toplevel();

    PRINT_STATUS(" ..velocity smoothing done.");
    double t1; t1 = WallclockTime = my_second();
    double timeall_local = timediff(t00_truestart, t1);
    CPU_Step[CPU_IMPROVDIFFCOMPUTE] += timeall_local;
}


#endif /* End TURB_DIFF_DYNAMIC */
