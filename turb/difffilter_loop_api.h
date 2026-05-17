/* turb/difffilter_loop_api.h — POD API surface for the DiffFilter runner port.
 *
 * Plain-old-data only: NO Kokkos, NO runner template headers. Included by the
 * host-only physics TUs (dynamic_diffusion.cc, dynamic_diffusion_velocities.cc)
 * so they can drive the GPU loops + name temporary_data_dyndiff WITHOUT being
 * pulled through nvcc_wrapper. The Kokkos/Spec definitions live in the
 * companion turb/difffilter_loop.h (included only by difffilter_loop.cc and
 * mesh/neighbor_loop_runner.cc).
 *
 * See OPEN_3d_difffilter_design.md §4.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) and Claude for GIZMO.
 */
#ifndef DIFFFILTER_LOOP_API_H
#define DIFFFILTER_LOOP_API_H

#include "../declarations/allvars.h"

#ifdef TURB_DIFF_DYNAMIC

/* Per-particle smoothed-gradient min/max bracket. Migrated verbatim from
 * dynamic_diffusion.cc (SSOT — the definition there is deleted). */
struct Quantities_for_Smooth_Gradients {
    double Velocity_hat[3];
};

/* Per-gas-particle accumulator for the dynamic Smagorinsky solve. Migrated
 * verbatim from dynamic_diffusion.cc; dynamic_diff_calc() owns one array of
 * these (DynamicDiffDataPasser, length N_gas) and dynamicdiff_gpu_toplevel
 * scatters the GPU loop output into it. Field set/types unchanged. */
struct temporary_data_dyndiff {
    struct Quantities_for_Smooth_Gradients Maxima;
    struct Quantities_for_Smooth_Gradients Minima;
    MyFloat  FilterWidth_hat;
    MyDouble Dynamic_numerator_hat;
    MyDouble Dynamic_denominator_hat;
    MyDouble GradVelocity_hat[3][3];
    MyDouble dynamic_fac[3][3];
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
    MyDouble dynamic_fac_const[3][3];
#endif
    MyDouble ProductVelocity_hat[3][3];
};

/* ----------------------------------------------------------------------------
 * Toplevel callers — defined in turb/difffilter_loop.cc.
 *
 * difffilter_vel_calc_gpu_toplevel:
 *   Builds the active gas list (DiffFilterSpec::is_active) and runs the
 *   DiffFilter velocity-smoothing loop. Output scattered directly into CellP[]
 *   by the Spec's apply_active_writeback. Replaces difffilter_evaluate_gpu().
 *
 * dynamicdiff_gpu_toplevel:
 *   One runner pass for the given external dynamic_iteration. Builds the active
 *   gas list (DynDiffSpec::is_active) and scatters the loop output into
 *   dddp[] (the caller's DynamicDiffDataPasser, length N_gas, indexed by
 *   particle index). Replaces dynamicdiff_evaluate_gpu().
 * -------------------------------------------------------------------------- */
void difffilter_vel_calc_gpu_toplevel(void);
void dynamicdiff_gpu_toplevel(int dynamic_iteration,
                              struct temporary_data_dyndiff *dddp);

#endif /* TURB_DIFF_DYNAMIC */
#endif /* DIFFFILTER_LOOP_API_H */
