#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified (re-arranged, consolidated, and a number of additional
 * sub-loops and other structures for e.g. feedback, gradients, neighbor operations on non-gas, etc,
 * added) by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void compute_grav_accelerations(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  PRINT_STATUS("Start gravity force computation...");

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      long_range_force();
      CPU_Step[CPU_MESH] += measure_time();
    }
#endif

  gravity_tree();		/* computes gravity accel. */

  /* For the first timestep, we redo it to allow usage of relative opening criterion for consistent accuracy */
  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0) {gravity_tree();}

  PRINT_STATUS(" ..gravity force computation done");
}



void compute_hydro_densities_and_forces(void)
{
  if(All.TotN_gas > 0)
    {
        PRINT_STATUS("Start hydrodynamics computation...");
        density();		/* computes density, and pressure */
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
        ags_density();
#endif
        force_update_hmax();	/* update kernel lengths in tree */
        /*! This function updates the hmax-values in tree nodes that hold SPH
         *  particles. These values are needed to find all neighbors in the
         *  hydro-force computation.  Since the Hsml-values are potentially changed
         *  in the SPH-denity computation, force_update_hmax() should be carried
         *  out before the hydrodynamical SPH forces are computed, i.e. after
         *  density().
         */

        PRINT_STATUS(" ..density & tree-update computation done...");

#ifdef HYDRO_VOLUME_CORRECTIONS
        cellcorrections_calc(); /* must be called after density, and after the update of hmax in the tree [because it depends on bi-directional search], but before gradients where quantities dependent on volumetric elements such as density are needed */
#endif
#ifdef TURB_DIFF_DYNAMIC
        dynamic_diff_vel_calc(); /* This must be called between density and gradient calculations */
#endif
#if defined(RADTRANSFER) && defined(GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION)
        rt_source_injection(); /* doing source injection here (just before interpolation and hydro gradients) is slightly more accurate for this setup, but not possible in total generality owing to dependence of some injection modules on quantities calculated below */
#endif
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
        interpolate_fluxes_opacities_gasgrains(); /* this must be called here for the computation of opacities and radiative quantity gradients below to be correct */
#endif
#ifdef GALSF /* PFH set of feedback routines; here because for e.g. strong SNe, obtain better stability if they are coupled discretely just -before- the hydro force is computed */
        compute_stellar_feedback();
#endif

        hydro_gradient_calc(); /* calculates the gradients of hydrodynamical quantities  */
        PRINT_STATUS(" ..gradient computation done.");

#ifdef TURB_DIFF_DYNAMIC
        dynamic_diff_calc(); /* This MUST be called immediately following gradient calculations */
#endif
        hydro_force();		/* adds hydrodynamical accelerations and computes du/dt  */
        compute_additional_forces_for_all_particles(); /* other accelerations that need to be computed are done here */
        PRINT_STATUS(" ..hydro force computation done.");

    } else {
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
        ags_density(); // if there are no gas particles but ags-all is active, still need to enter this loop //
        force_update_hmax();    /* update kernel lengths in tree */
#endif
        compute_additional_forces_for_all_particles();
    }
}



void compute_additional_forces_for_all_particles(void)
{
#ifdef DM_FUZZY
    DMGrad_gradient_calc();
#endif
#if defined(DM_FUZZY) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(DM_SIDM)
    AGSForce_calc();
#endif
#ifdef GRAIN_FLUID
    apply_grain_dragforce(); /* if we are solving a coupled set of grains via aerodynamic drag, this is where their acceleration should be calculated */
#endif
}




#ifdef GALSF
void compute_stellar_feedback(void)
{
    CPU_Step[CPU_MISC] += measure_time();

#ifdef GALSF_FB_MECHANICAL /* check the mechanical sources of feedback */
    mechanical_fb_calc_toplevel();  /* call the parent loop for the different mechanical fb sub-loops */
    MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_SNIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */
#endif

#ifdef GALSF_FB_THERMAL
    thermal_fb_calc(); /* thermal feedback */
    MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_SNIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */
#endif
    
    
    
    CPU_Step[CPU_MISC] += measure_time();
}
#endif // GALSF //
