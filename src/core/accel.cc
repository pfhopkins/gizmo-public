#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../declarations/allvars.h"
#include "../core/proto.h"

/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active gas/fluid cells is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */
/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been rewritten and extensively modified (re-arranged, consolidated, and a number of additional
 * sub-loops and other structures for e.g. feedback, gradients, neighbor operations on non-gas, etc,
 * added) by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void compute_grav_accelerations(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  PRINT_STATUS("Start gravity force computation...");

#ifdef PMGRID
  {
    if(All.PM_Ti_endstep == All.Ti_Current)
      {
        long_range_force();
        CPU_Step[CPU_MESH] += measure_time();
      }
  }
#endif

  {
    gravity_tree();		/* computes gravity accel. */

    /* For the first timestep, we redo it to allow usage of relative opening criterion for consistent accuracy.
       The estimator this seeds is refreshed only when the tree is built, and the import is pruned against it,
       so the repeat pass has to rebuild -- otherwise it walks the first pass's tree and seeds nothing. */
    if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0) {TreeReconstructFlag = 1; gravity_tree();}
  }

  PRINT_STATUS(" ..gravity force computation done");
}



void compute_hydro_densities_and_forces(void)
{
  if(All.TotN_gas > 0)
    {
        PRINT_STATUS("Start hydrodynamics computation...");
        /* Hydro corridor entry: decide Mode A vs Mode B once for the whole
         * corridor span (density → cellcorrections → gradients → hydro_force),
         * so every consumer in the span dispatches coherently. */
        gizmo_hydro_corridor_decide_mode();
        /* density() internally handles ghost_exchange prep + redo (neighbor-list path)
           via the ghost_symlist_lifecycle helpers. Same for gradients and hydro_force. */
        double t_bench_density_start = my_second(); double child0_density = CPU_ChildCharged;
        density();		/* computes density, and pressure */
        double t_bench_density_only = cpu_minus_children(timediff(t_bench_density_start, my_second()), child0_density);
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
        ags_density();
#endif
        force_update_hmax();	/* update kernel lengths in tree */
        /*! This function updates the hmax-values in tree nodes that hold gas. These values are needed to find all neighbors in the hydro-force computation.  Since the KernelRadius-values are potentially changed in the gas-denity computation, force_update_hmax() should be carried out before the hydrodynamical forces are computed, i.e. after density(). */

        PRINT_STATUS(" ..density & tree-update computation done...");

#if defined(RADTRANSFER) && defined(GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION)
        rt_source_injection(); /* doing source injection here (just before interpolation and hydro gradients) is slightly more accurate for this setup, but not possible in total generality owing to dependence of some injection modules on quantities calculated below */
#endif
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
        interpolate_fluxes_opacities_gasgrains(); /* this must be called here for the computation of opacities and radiative quantity gradients below to be correct */
        CPU_Step[CPU_RTNONFLUXOPS] += measure_time();
#endif
#ifdef GALSF /* PFH set of feedback routines; here because for e.g. strong SNe, obtain better stability if they are coupled discretely just -before- the hydro force is computed */
        compute_stellar_feedback();
#endif

        /* Hydro corridor begin: in Mode A (any rank count) imports the gas
         * ghost pool and builds the shared active list + symmetric gas CSR
         * consumed by cellcorrections + gradients + hydro_force via
         * gizmo_hydro_corridor_external_csr(). In Mode B builds only the
         * shared active list (request-driven, no CSR). See hydro_corridor.h
         * for the ownership contract and mid-span refresh points.
         *
         * It sits below the injection/feedback block above because each of
         * those routines runs its own neighbor loop, and a neighbor loop's
         * ghost import first tears down the single global ghost pool. With
         * them above the corridor they can no longer tear its pool, so the
         * refresh before gradients stays on the in-place path instead of
         * re-importing and rebuilding the CSR for a neighbor list that never
         * changed. This does NOT make the span tear-free: a loop inside the
         * span that still runs its own ghost import tears it just the same.
         * Under TURB_DIFF_DYNAMIC, dynamic_diff_calc between gradients and
         * hydro_force does exactly that and costs one corridor re-import per
         * step. Consuming this CSR, the way cellcorrections and gradients do,
         * is what removes such a tear. */
        gizmo_hydro_corridor_begin();
        CPU_Step[CPU_DENSMISC] += measure_time();

#ifdef HYDRO_VOLUME_CORRECTIONS
        /* Running below feedback is NOT neutral for this routine, by intent.
         * Its volume correction reads only positions, kernel radii and the
         * zeroth-order volumes, none of which feedback touches, so that part
         * is unchanged. But its closing Density = Mass/Volume_1 and pressure
         * update read mass and internal energy, which feedback DOES write, so
         * in a feedback build those primitives are now derived from the
         * post-feedback state — the same state the gradients below consume,
         * which is the self-consistent choice. */
        cellcorrections_calc(); /* must be called after density, and after the update of hmax in the tree [because it depends on bi-directional search], but before gradients where quantities dependent on volumetric elements such as density are needed */
#endif
#ifdef TURB_DIFF_DYNAMIC
        /* Between density and gradients, and below feedback: the diffusion
         * coefficient this filter estimates is applied to the post-feedback
         * state, so measuring it on the post-feedback velocity field is the
         * self-consistent choice. */
        dynamic_diff_vel_calc();
#endif

        double t_bench_grad_start = my_second(); double child0_grad = CPU_ChildCharged;
        hydro_gradient_calc(); /* calculates the gradients of hydrodynamical quantities  */
#ifdef MHD_MODIFIED_GRADIENT
        {   /* determine whether the active gas fraction is large enough to justify the global MG solve */
            long long ngas_active_local = 0;
            for(int tbin = 0; tbin < TIMEBINS; tbin++) {if(TimeBinActive[tbin]) {ngas_active_local += TimeBinCountGas[tbin];}}
            long long ngas_active_global; MPI_Allreduce(&ngas_active_local, &ngas_active_global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            double active_gas_fraction = (All.TotN_gas > 0) ? ((double)ngas_active_global / (double)All.TotN_gas) : 1.0;
            if(active_gas_fraction >= All.ActiveFractionForMGSweep) {
                All.Flag_SkipMGSolve = 0;
                const double t_mg_solve = my_second();
                mg_gradient_correction_calc(); /* MG method: global sparse-matrix solve for exact div(B)=0 correction coefficients (Tu et al. 2026) */
                /* The solve is charged to its own row rather than left in the gradient
                 * residual.  It has to be charged SOMEWHERE: measure_time() discounts only
                 * what a child bucket claims, so an unclaimed span of this size lands both
                 * in the enclosing bracket and in the next residual, and is counted twice. */
                cpu_charge_child(CPU_MHD_MG, timediff(t_mg_solve, my_second()));
            } else {
                All.Flag_SkipMGSolve = 1;
                if(ThisTask == 0) {PRINT_STATUS("Skipping MG global solve (active gas fraction %g < %g), using CG fallback for active cells", active_gas_fraction, All.ActiveFractionForMGSweep);}
            }
        }
#endif
#if defined(COOLING) && defined(GALSF_FB_FIRE_RT_LONGRANGE)
        selfshield_local_incident_uv_flux(); /* needs to be called after gravity tree (where raw flux is calculated) and the local gradient calculation (GradRho) to properly self-shield the particles that had this calculated */
#endif
        double t_bench_grad = cpu_minus_children(timediff(t_bench_grad_start, my_second()), child0_grad);
        PRINT_STATUS(" ..gradient computation done.");

#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES >= 4)
        special_rt_feedback_injection(); /* do before proper hydro loop */
#endif

#ifdef TURB_DIFF_DYNAMIC
        dynamic_diff_calc(); /* This MUST be called immediately following gradient calculations */
#endif
        double t_bench_hydro_start = my_second(); double child0_hydro = CPU_ChildCharged;
        hydro_force();		/* adds hydrodynamical accelerations and computes du/dt  */
        double t_bench_hydro = cpu_minus_children(timediff(t_bench_hydro_start, my_second()), child0_hydro);
        /* Hydro corridor exit: reset mode to UNSET. Any unexpected
         * gizmo_hydro_corridor_get_mode() outside the corridor span will
         * read UNSET and is therefore detectable. */
        gizmo_hydro_corridor_end();
        compute_additional_forces_for_all_particles(); /* other accelerations that need to be computed are done here */
        /* Feed GPU-path kernel timings into CPU_Step accumulators for cpu.txt output.
           drift+ghost+redo feeds are done inside density()/gradients()/hydro_force() via
           the ghost_symlist_lifecycle helpers. */
        CPU_Step[CPU_DENSCOMPUTE] += t_bench_density_only;
        CPU_Step[CPU_DENSWAIT] += t_bench_grad;
        CPU_Step[CPU_HYDCOMPUTE] += t_bench_hydro;
        PRINT_STATUS(" ..hydro force computation done.");

    } else {
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
        ags_density(); // if there are no gas particles but ags-all is active, still need to enter this loop //
        force_update_hmax();    /* update kernel lengths in tree */
#endif
        compute_additional_forces_for_all_particles();
    }
}



void compute_additional_forces_for_all_particles(void)
{
#if defined(DM_FUZZY)
    DMGrad_gradient_calc();
#endif
#if defined(CBE_INTEGRATOR_WITHGRADIENTS)
    /* Persistent CBE basis-moment gradients refreshed for AGSForce-active
     * particles each call; inactive particles retain prior-step gradients
     * (hydro semantics). Independent of DMGrad_gradient_calc (different
     * Spec, different storage); both fire when both flags are on. */
    CBEGrad_gradient_calc();
#endif
#if defined(DM_FUZZY) || defined(CBE_INTEGRATOR) || defined(DM_SIDM)
    AGSForce_calc();
#endif
#if defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION)
    apply_grain_dragforce(); /* if we are solving a coupled set of grains via aerodynamic drag, this is where their acceleration should be calculated */
#endif
}




#ifdef GALSF
void compute_stellar_feedback(void)
{
    CPU_Step[CPU_MISC] += measure_time();
#ifdef GALSF_LIMIT_FBTIMESTEPS_FROM_BELOW
    if(All.Dt_Since_LastFBCalc_Gyr > All.Dt_Min_Between_FBCalc_Gyr) {All.Dt_Since_LastFBCalc_Gyr = 0;}
    All.Dt_Since_LastFBCalc_Gyr += All.TimeStep / All.cf_hubble_a * UNIT_TIME_IN_GYR; // augment by timestep
    if(All.Dt_Since_LastFBCalc_Gyr < All.Dt_Min_Between_FBCalc_Gyr) {return;}
#endif

#ifdef GALSF_FB_MECHANICAL /* check the mechanical sources of feedback */
    mechanical_fb_calc_toplevel();  /* call the parent loop for the different mechanical fb sub-loops */
    gizmo_exit_bad_stop_if_requested("accel:after_mechanical_fb"); CPU_Step[CPU_SNIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */
#endif

#ifdef GALSF_FB_THERMAL
    thermal_fb_calc(); /* thermal feedback */
    gizmo_exit_bad_stop_if_requested("accel:after_thermal_fb"); CPU_Step[CPU_SNIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */
#endif

#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
    HII_heating_singledomain(); /* local photo-ionization heating */
    gizmo_exit_bad_stop_if_requested("accel:after_hii_heating"); CPU_Step[CPU_HIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */
#endif

#ifdef GALSF_FB_FIRE_RT_LOCALRP
    radiation_pressure_winds_consolidated(); /* local radiation pressure */
    gizmo_exit_bad_stop_if_requested("accel:after_local_rp"); CPU_Step[CPU_LOCALWIND] += measure_time(); /* collect timings and reset clock for next timing */
#endif
    
    CPU_Step[CPU_MISC] += measure_time();
}
#endif // GALSF //
