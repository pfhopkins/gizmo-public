/* gravity/ags_density_loop.cc — host-only hooks + ghost-writeback manifest
 * for AgsDensitySpec.
 *
 * KOKKOS_INLINE_FUNCTION hooks (load_active, load_neighbor, pair_kernel,
 * zero_accum, active_subgroup_key) and the inline pair body
 * (ags_density_pair_kernel_body) live in gravity/ags_density_loop.h. This
 * TU holds:
 *   - host-side per-active radius + per-call scalar capture
 *   - populate / cleanup of the device context (single sticky-call-scope
 *     need_wakeup_uvm scratch — no per-iter reset, since it accumulates
 *     across all iters of the call)
 *   - apply_active_writeback (no-op for AgsDensitySpec — after_iter writes
 *     post-processed values to P[i] directly on host)
 *   - merge_accum (Mode B remote peer merge — manifest pattern)
 *   - after_iter (legacy convergence test + bisection; also syncs
 *     P[i].AGS_KernelRadius every iter)
 *   - after_iter_global (iter > 10 print only)
 *   - ghost-writeback bundle manifest (PARTICLE_MAX wakeup) + lifecycle
 *
 * Replaces the now-deleted gravity/ags_density_gpu.cc (retired in the
 * cleanup after two-binary parity passed).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"  /* MUST precede allvars.h: installs device-pass `#define All AllDeviceMirror` redirect before cell_data.h is parsed */
#include "../declarations/allvars.h"
#include "../declarations/gpu_numeric_macros.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"               /* MUST precede ags_density_loop.h */
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_writeback_ops.h"
#include "../mesh/ghost_symlist_lifecycle.h"  /* gizmo_ghost_safety_factor */
#include "ags_density_loop.h"

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE

/* AGS_DSOFT_TOL (softening growth tolerance per step) is defined globally in
 * precompiler_logic.h so after_iter's minsoft/maxsoft clamp matches the AGS
 * radius solver in gravity/ags_rkern.cc. */

/* Targeted diagnostic trace.
 *
 * Compile with -DGIZMO_NLR_AGS_DEBUG_ID_TRACE=<particle_ID> to enable
 * per-particle per-iter state dump for the specified ID. Used to diagnose
 * a first-Vista-validation regression (runner crashes endrun(888) on
 * type-5 sink ID=2190205 with soft=0; legacy passes same config).
 *
 * Production builds without the macro carry no overhead and no print
 * sites. The legacy mirror in ags_rkern.cc was deleted with the rest of
 * the legacy body in the cleanup. */
#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
static inline void ags_debug_trace_id_state(const char *label, int i)
{
    if((long long)P[i].ID != (long long)GIZMO_NLR_AGS_DEBUG_ID_TRACE) return;
    printf("[AGS_TRACE rank=%d %s] ID=%lld Type=%d Mass=%g TimeBin=%d "
           "AGS_KernelRadius=%g KernelRadius=%g ForceSoftening=%g "
           "AGS_zeta=%g NumNgb=%g AGS_vsig=%g\n",
           ThisTask, label, (long long)P[i].ID, (int)P[i].Type,
           (double)P[i].Mass, (int)P[i].TimeBin,
           (double)P[i].AGS_KernelRadius, (double)P[i].KernelRadius,
           (double)P[i].ForceSoftening,
           (double)P[i].AGS_zeta, (double)P[i].NumNgb,
           (double)P[i].AGS_vsig);
    fflush(stdout);
}
#endif

/* ============================================================================
 * PHYSICS HOOKS
 * ========================================================================== */

double AgsDensitySpec::search_radius(const neighbor_loop_args& args,
                                      int /*active_slot*/, int i)
{
    return (double)args.P[i].AGS_KernelRadius;
}

/* Per-call scalars — snapshot TimeBinActive[] + AGS limits + cosmology
 * factors once per outer call. The kernel reads scalars.TimeBinActive[]
 * by value (captured into the device lambda via the runner's
 * trivially-copyable CallScalars staging).
 *
 * Uses nlr_host_all_ptr() for explicit host-snapshot intent — the underlying
 * macro redirect is now device-pass-only, but the accessor convention
 * documents which side of host/device this read lives on. A bug here
 * originally surfaced via debug-trace instrumentation (AGS_DesNumNgb=0
 * → bisection inversion → radius collapse → endrun(888)). */
AgsDensitySpec::CallScalars
AgsDensitySpec::populate_call_scalars(const neighbor_loop_args& /*args*/)
{
    const struct global_data_all_processes *h = nlr_host_all_ptr();
    CallScalars scalars;
    scalars.common = nlr_common_scalars_from_all();
    for(int k = 0; k < TIMEBINS; k++) scalars.TimeBinActive[k] = TimeBinActive[k];
    scalars.Time                   = h->Time;
    scalars.TimeBegin              = h->TimeBegin;
    scalars.AGS_DesNumNgb          = h->AGS_DesNumNgb;
    scalars.AGS_MaxNumNgbDeviation = h->AGS_MaxNumNgbDeviation;
    scalars.fac_mu                 = -3.0 / h->cf_atime;
    return scalars;
}

/* ============================================================================
 * DEVICE CONTEXT LIFECYCLE
 *
 * populate allocates + zeros ctx.need_wakeup_uvm (single UVM int).
 * cleanup propagates ctx.need_wakeup_uvm into NeedToWakeupParticles_local
 * (legacy ags_density_gpu.cc:210 path) and frees the UVM.
 *
 * STICKY ACROSS ITERS — a per-iter reset would lose non-final-iter
 * wakeups (iter 1 writes, iter N doesn't, cleanup sees 0). The flag
 * accumulates across all iters of the iterative call;
 * legacy ran one evaluator-call per outer iter so the legacy
 * *d_need_wakeup was per-call too. There is NO reset_per_iter_device_context
 * hook for AgsDensitySpec.
 *
 * Runs unconditionally via NlrDeviceContextCleanupGuard at runner exit on
 * every dispatch path.
 * ========================================================================== */

void AgsDensitySpec::populate_device_context(const neighbor_loop_args& /*args*/,
                                              DeviceContext& ctx)
{

    /* Single-int UVM scratch — pair_kernel atomic_or-s a 1 here on any
     * wakeup write. Accumulates across all iters; cleanup_device_context
     * propagates to NeedToWakeupParticles_local once at runner exit. */
    ctx.need_wakeup_uvm = NULL;
    ctx.wakeup_dirty_base = WakeupDirty;   /* global UVM sidecar base; kernel marks WakeupDirty[j] on wakeup */
    int *nw = (int *) gizmo_gpu_alloc_shared(sizeof(int), NULL);
    if(!nw) {
        ctx.populate_failed = 1;
        gizmo_request_controlled_stop(7724,
            "adaptive-softening density: could not stage the wakeup flag; "
            "the kernel radii are not updated",
            __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    *nw = 0;
    ctx.need_wakeup_uvm = nw;
}

void AgsDensitySpec::cleanup_device_context(const neighbor_loop_args& /*args*/,
                                             DeviceContext& ctx)
{
    if(ctx.need_wakeup_uvm) {
        if(*ctx.need_wakeup_uvm) {
            /* Path 1 of design v0.4.2 §6: active rank's local-side wakeup
             * write raises the global flag. Path 2 (home rank receiving
             * cross-rank delta) is handled inside ghost_writeback.cc when
             * wakeups_applied > 0. */
            NeedToWakeupParticles_local = 1;
        }
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(ctx.need_wakeup_uvm);
        ctx.need_wakeup_uvm = nullptr;
    }
}

/* ============================================================================
 * APPLY_ACTIVE_WRITEBACK — NO-OP for AgsDensitySpec.
 *
 * The original per_active_accum[active_slot] write was multi-subgroup-
 * collision-prone (active_slot is subgroup-local). The cleaner design is:
 * after_iter writes the post-processed values
 * directly to P[i] on host (legacy semantics — rkern.cc:179-196 does the
 * same). The caller's post-runner finalize pass reads P[i] for the
 * final-final operations (rkern.cc:431-453). No host aux buffer needed.
 *
 * This hook stays declared (Spec contract requirement) but does nothing.
 * Runner invokes it but the call collapses to no work.
 * ========================================================================== */
void AgsDensitySpec::apply_active_writeback(const neighbor_loop_args& /*args*/,
                                             int                       /*active_slot*/,
                                             int                       /*i*/,
                                             const AccumData&          /*accum*/)
{
    /* Intentionally empty — see banner. */
}


/* ============================================================================
 * AFTER_ITER — per-active per-iter post-process + convergence test.
 *
 * Legacy equivalent: gravity/ags_rkern.cc:174-405 (the inner block of the
 * outer do-while loop). Verbatim port with these changes:
 *
 *   - Per-iter accumulator is the runner's AccumData `accum` (zero_accum
 *     resets per iter) — replaces direct P[ii].NumNgb mutation that legacy
 *     uses as cross-bm-subgroup scratch. Each active belongs to exactly
 *     one bm subgroup (actives_partition_by_subgroup=true), so no cross-
 *     subgroup merge is needed; accum is the single source.
 *
 *   - Iterative state (Left, Right, AGS_Prev, set_to_max/minrkern flags)
 *     lives in IterScratch (ctx.scratch), not in standalone arrays as in
 *     legacy.
 *
 *   - Radius mutation is via IterResult::AdjustRadius (runner-owned); we
 *     never write to P[i].AGS_KernelRadius directly in the loop. The
 *     final converged radius lands in P[i].AGS_KernelRadius via the
 *     post-runner caller scatter.
 *
 *   - Post-processed per-iter values (normalized NumNgb, inverted
 *     DrkernNgbFactor, Particle_DivVel correction, etc.) are written to
 *     P[i] here on the host as a side-effect — matches legacy lines
 *     183-196 which do the same. These values are read by the runner's
 *     next iter's pair_kernel ONLY via the runner accum path (kernel does
 *     not read P[i].NumNgb on subsequent iters), so the side-effect is
 *     contained.
 *
 *   - TimeBin negation (legacy line 279 / 403) is NOT done here — design
 *     v0.4.2 §10.5 lock: TimeBin marker stays caller-side post-runner.
 *
 *   - AGS_FACE_CALCULATION_IS_ACTIVE NV_T inversion IS done here (below):
 *     the accumulated raw moment matrix is inverted into the symmetric face
 *     operator P[i].NV_T via ags_invert_nvt_for_faces, and its condition
 *     number drives the neighbor-count expansion in the convergence block.
 * ========================================================================== */
IterResult AgsDensitySpec::after_iter(const AfterIterContext<AgsDensitySpec>& ctx,
                                       const AccumData&                       accum)
{
    const int    i                 = ctx.i;
    const int    iter              = ctx.iter_index;
    const double current_h         = ctx.h_search_current;
    AgsDensityIterScratch& scratch = ctx.scratch;
    const AgsDensityCallScalars& scalars = ctx.scalars;

    /* (0) Iter-0 initialization of IterScratch. The runner zero-memsets
     * scratch at the start of every iterative call (see runner.cc:2370),
     * so on iter 0 all fields are 0. AGS_Prev needs to capture the
     * ENTRY-time radius for the AGS_DSOFT_TOL minsoft/maxsoft clamp
     * (rkern.cc:90 + :205-206). Left = Right = 0 and the set_to_*rkern
     * flags = 0 are already correct from the zero-memset; AGS_Prev is
     * the only field needing explicit iter-0 setup. */
    if(iter == 0) {
        scratch.AGS_Prev = current_h;
        /* Caller resets P[i].wakeup = 0 + AGS_vsig = 0 pre-runner
         * (legacy rkern.cc:96-97 semantics, kept caller-side per
         * design v0.4.2 §8). */
    }

#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
    if((long long)P[i].ID == (long long)GIZMO_NLR_AGS_DEBUG_ID_TRACE) {
        printf("[AGS_TRACE rank=%d after_iter ENTER iter=%d] ID=%lld "
               "current_h=%g AGS_Prev=%g accum.Ngb=%g accum.AGS_vsig_max=%g "
               "Left=%g Right=%g\n",
               ThisTask, iter, (long long)P[i].ID,
               current_h, scratch.AGS_Prev,
               accum.Ngb, accum.AGS_vsig_max,
               scratch.Left, scratch.Right);
        fflush(stdout);
    }
#endif

    /* Final-radius writeback: the runner owns
     * its own per-active radii_uvm[sg][slot] (mutated via AdjustRadius
     * returns), but P[i].AGS_KernelRadius is what the post-runner
     * caller finalize pass + downstream physics reads. Sync P[i] to
     * the runner's current view EVERY iter so converged actives have
     * their final radius land in P[i] without needing a separate
     * writeback channel. Updated below at each return path. */
    P[i].AGS_KernelRadius = (MyFloat)current_h;

    /* Mode A per-iter state visibility: the
     * P[i].* writes below are visible to subsequent Mode A pair_kernel
     * walks because of UVM-canonical particle aliasing (system/
     * gpu_particles_arena.cc: arena.P / arena.CellP
     * are pure pointer aliases of host P / CellP under SharedSpace).
     * Mode B request-driven walks rebuild the per-call slab from host
     * P each call, so the host writes are also visible there. No
     * separate arena-sync needed. */

    /* (1) Post-process per-iter accumulator into normalized values. Writes
     * into P[i] mirror legacy rkern.cc:179-196. */
    double NumNgb_raw          = accum.Ngb;
    double DrkernNgbFactor_raw = accum.DrkernNgb;
    double Particle_DivVel_raw = accum.Particle_DivVel;
    double AGS_zeta_raw        = accum.AGS_zeta;
    double AGS_vsig            = accum.AGS_vsig_max;

#ifdef DM_FUZZY
    /* AGS_Density tracks Mass * NumNgb (raw weighted count). Legacy
     * rkern.cc:179. */
    P[i].AGS_Density = (MyFloat)((double)P[i].Mass * NumNgb_raw);
#endif

    double NumNgb_norm          = 0;
    double DrkernNgbFactor_norm = 0;
    double Particle_DivVel_norm = 0;
    if(NumNgb_raw > 0) {
        DrkernNgbFactor_norm = DrkernNgbFactor_raw * current_h / (NUMDIMS * NumNgb_raw);
        Particle_DivVel_norm = Particle_DivVel_raw / NumNgb_raw;
        NumNgb_norm          = NumNgb_raw * VOLUME_NORM_COEFF_FOR_NDIMS * pow(current_h, NUMDIMS);
    }

    /* Inverse-of-volume-element (legacy lines 191-195). */
    if(DrkernNgbFactor_norm > -0.9) DrkernNgbFactor_norm = 1.0 / (1.0 + DrkernNgbFactor_norm);
    else                            DrkernNgbFactor_norm = 1.0;
    Particle_DivVel_norm *= DrkernNgbFactor_norm;

    P[i].NumNgb          = (MyFloat)NumNgb_norm;
    P[i].DrkernNgbFactor = (MyFloat)DrkernNgbFactor_norm;
    P[i].Particle_DivVel = (MyFloat)Particle_DivVel_norm;
    P[i].AGS_zeta        = (MyFloat)AGS_zeta_raw;
    P[i].AGS_vsig        = (MyFloat)AGS_vsig;
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    /* Invert the accumulated raw moment matrix into the symmetric, fully-
     * populated face operator P[i].NV_T (the inverted tensor that the
     * CBE/dm_fuzzy face construction and the DMGrad gradient estimator both
     * consume). Storing the raw moment matrix here -- the prior behavior --
     * produced unphysical faces (wrong magnitude, sometimes wrong sign) with
     * no conditioning control. The condition number captured here drives the
     * neighbor-count expansion in the convergence block below. Mirrors the
     * hydro density path (hydro/density_loop.cc:508-633). */
    {
        const double sym6[6] = { (double)accum.NV_T_00, (double)accum.NV_T_01, (double)accum.NV_T_02,
                                 (double)accum.NV_T_11, (double)accum.NV_T_12, (double)accum.NV_T_22 };
        double Tinv[3][3];
        AgsNvtInversion inv = ags_invert_nvt_for_faces(sym6, current_h, Tinv);
        for(int a = 0; a < 3; a++) { for(int b = 0; b < 3; b++) { P[i].NV_T[a][b] = (MyDouble)Tinv[a][b]; } }
        if(inv.cn_expansion > scratch.condition_number_max) { scratch.condition_number_max = inv.cn_expansion; }
    }
#endif

    /* (2) Convergence test. minsoft / maxsoft clamping from
     * scratch.AGS_Prev (rkern.cc:201-207). */
    double minsoft = ags_return_minsoft(i);
    double maxsoft = ags_return_maxsoft(i);
    if(scalars.Time > scalars.TimeBegin) {
        if(scratch.AGS_Prev * AGS_DSOFT_TOL > minsoft) minsoft = scratch.AGS_Prev * AGS_DSOFT_TOL;
        if(scratch.AGS_Prev / AGS_DSOFT_TOL < maxsoft) maxsoft = scratch.AGS_Prev / AGS_DSOFT_TOL;
    }
    double desnumngb    = scalars.AGS_DesNumNgb;
    double desnumngbdev = scalars.AGS_MaxNumNgbDeviation;

#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    /* Condition-number-driven neighbor expansion (mirrors hydro
     * density_loop.cc:619-633): an ill-conditioned face matrix inflates the
     * desired neighbor count, expanding the search radius next iteration,
     * which lowers the condition number. Applied to BOTH target and deviation,
     * and BEFORE the late-iter deviation relaxation below so that relaxation
     * cannot silently undo the expanded target. */
    {
        const double c0 = 0.1 * (double)CONDITION_NUMBER_DANGER;
        const double cn = scratch.condition_number_max;
        double ncorr_ngb = 1.0;
        if(cn > c0) { ncorr_ngb = sqrt(1.0 + (cn - c0) / (double)CONDITION_NUMBER_DANGER); }
        if(ncorr_ngb > 2.0) { ncorr_ngb = 2.0; }
        desnumngb    *= ncorr_ngb;
        desnumngbdev *= ncorr_ngb;
    }

    /* AGS_FACE relaxes desnumngbdev later (iter > 10) since the face
     * inversion is more sensitive (legacy rkern.cc:213-214). */
    if(iter > 10) {
        const double relax = exp(0.1 * log(desnumngb / (16.0 * desnumngbdev)) * ((double)iter - 9.0));
        const double cap   = 0.25 * desnumngb;
        const double bumped = desnumngbdev * relax;
        desnumngbdev = (cap < bumped) ? cap : bumped;
    }
#else
    if(iter > 4) {
        const double relax = exp(0.1 * log(desnumngb / (16.0 * desnumngbdev)) * ((double)iter - 3.0));
        const double cap   = 0.25 * desnumngb;
        const double bumped = desnumngbdev * relax;
        desnumngbdev = (cap < bumped) ? cap : bumped;
    }
#endif
    if(scalars.Time <= scalars.TimeBegin) {
        if(desnumngbdev > 0.0005) desnumngbdev = 0.0005;
        if(iter > 50) {
            const double relax = exp(0.1 * log(desnumngb / (16.0 * desnumngbdev)) * ((double)iter - 49.0));
            const double cap   = 0.25 * desnumngb;
            const double bumped = desnumngbdev * relax;
            desnumngbdev = (cap < bumped) ? cap : bumped;
        }
    }

    /* "Normal" range check (rkern.cc:221-223). */
    int redo_particle = 0;
    if((NumNgb_norm < (desnumngb - desnumngbdev) && current_h < 0.999 * maxsoft) ||
       (NumNgb_norm > (desnumngb + desnumngbdev) && current_h > 1.001 * minsoft)) {
        redo_particle = 1;
    }

    /* Max-kernel cap (rkern.cc:227-241). */
    double new_h = current_h;
    scratch.set_to_maxrkern = 0;
    if((current_h >= 0.999 * maxsoft) && (NumNgb_norm < (desnumngb - desnumngbdev))) {
        redo_particle = 0;
        if(current_h == maxsoft) {
            scratch.set_to_maxrkern = 0;
        } else {
            redo_particle           = 1;
            new_h                   = maxsoft;
            scratch.set_to_maxrkern = 1;
        }
    }

    /* Min-kernel cap (rkern.cc:244-259). */
    scratch.set_to_minrkern = 0;
    if((current_h <= 1.001 * minsoft) && (NumNgb_norm > (desnumngb + desnumngbdev))) {
        redo_particle = 0;
        if(current_h == minsoft) {
            scratch.set_to_minrkern = 0;
        } else {
            redo_particle           = 1;
            new_h                   = minsoft;
            scratch.set_to_minrkern = 1;
        }
    }

    if(!redo_particle) {
#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
        if((long long)P[i].ID == (long long)GIZMO_NLR_AGS_DEBUG_ID_TRACE) {
            printf("[AGS_TRACE rank=%d after_iter Converged] ID=%lld iter=%d "
                   "current_h=%g AGS_Prev=%g NumNgb_raw=%g NumNgb_norm=%g "
                   "Left=%g Right=%g minsoft=%g maxsoft=%g "
                   "set_to_maxrkern=%d set_to_minrkern=%d\n",
                   ThisTask, (long long)P[i].ID, iter,
                   current_h, scratch.AGS_Prev, NumNgb_raw, NumNgb_norm,
                   scratch.Left, scratch.Right, minsoft, maxsoft,
                   scratch.set_to_maxrkern, scratch.set_to_minrkern);
            fflush(stdout);
        }
#endif
        /* Converged — runner removes from active_set. Caller-side final
         * pass picks up P[i].AGS_KernelRadius at its current value. */
        return IterResult{IterStatus::Converged, 0.0};
    }

    /* (3) Redo path: maxiter warning (rkern.cc:263-269). Match legacy
     * print format verbatim so log grep parity is preserved. */
    if(iter >= MAXITER - 10) {
        PRINT_WARNING("AGS: i=%d task=%d ID=%llu Type=%d KernelRadius=%g Drkern=%g Left=%g "
                      "Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g "
                      "maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)\n",
                      i, ThisTask, (unsigned long long)P[i].ID, (int)P[i].Type,
                      current_h, P[i].DrkernNgbFactor, scratch.Left, scratch.Right,
                      (float)P[i].NumNgb, scratch.Right - scratch.Left,
                      scratch.set_to_maxrkern, scratch.set_to_minrkern,
                      minsoft, maxsoft, desnumngb, desnumngbdev, redo_particle,
                      P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
    }

    /* Early-converged: Left/Right brackets are tight (rkern.cc:274-281).
     * Legacy decrements npleft then sets TimeBin marker; we just return
     * Converged here (TimeBin marker handled caller-side). */
    if(scratch.Left > 0 && scratch.Right > 0 &&
       (scratch.Right - scratch.Left) < 1.0e-3 * scratch.Left) {
        return IterResult{IterStatus::Converged, 0.0};
    }

    if(scratch.set_to_maxrkern == 0 && scratch.set_to_minrkern == 0) {
        /* Standard bisection update (rkern.cc:284-393). */
        if(NumNgb_norm < (desnumngb - desnumngbdev)) {
#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
            if((long long)P[i].ID == (long long)GIZMO_NLR_AGS_DEBUG_ID_TRACE) {
                printf("[AGS_TRACE rank=%d after_iter BISECT_LEFT iter=%d] ID=%lld "
                       "NumNgb_norm=%g desnumngb=%g desnumngbdev=%g "
                       "current_h=%g Left_before=%g (will set Left=current_h)\n",
                       ThisTask, iter, (long long)P[i].ID,
                       NumNgb_norm, desnumngb, desnumngbdev,
                       current_h, scratch.Left);
                fflush(stdout);
            }
#endif
            if(current_h > scratch.Left) scratch.Left = current_h;
        } else {
#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
            if((long long)P[i].ID == (long long)GIZMO_NLR_AGS_DEBUG_ID_TRACE) {
                printf("[AGS_TRACE rank=%d after_iter BISECT_RIGHT iter=%d] ID=%lld "
                       "NumNgb_norm=%g desnumngb=%g desnumngbdev=%g "
                       "current_h=%g Right_before=%g (will set Right=current_h)\n",
                       ThisTask, iter, (long long)P[i].ID,
                       NumNgb_norm, desnumngb, desnumngbdev,
                       current_h, scratch.Right);
                fflush(stdout);
            }
#endif
            if(scratch.Right != 0) {
                if(current_h < scratch.Right) scratch.Right = current_h;
            } else {
                scratch.Right = current_h;
            }
        }

        if(scratch.Right > 0 && scratch.Left > 0) {
            /* Geometric interpolation between Left/Right (rkern.cc:302-326). */
            double maxjump = 0;
            if(iter > 1) maxjump = 0.2 * log(scratch.Right / scratch.Left);
            if(NumNgb_norm > 1) {
                double jumpvar = P[i].DrkernNgbFactor * log(desnumngb / NumNgb_norm) / NUMDIMS;
                if(iter > 1) {
                    if(fabs(jumpvar) < maxjump) {
                        jumpvar = (jumpvar < 0) ? -maxjump : maxjump;
                    }
                }
                new_h = current_h * exp(jumpvar);
            } else {
                new_h = current_h * 2.0;
            }
            if(new_h < scratch.Right && new_h > scratch.Left) {
                if(iter > 1) {
                    const double hfac = exp(maxjump);
                    if(new_h > scratch.Right / hfac) new_h = scratch.Right / hfac;
                    if(new_h < scratch.Left  * hfac) new_h = scratch.Left  * hfac;
                }
            } else {
                if(new_h > scratch.Right) new_h = scratch.Right;
                if(new_h < scratch.Left)  new_h = scratch.Left;
                new_h = pow(new_h * scratch.Left * scratch.Right, 1.0/3.0);
            }
        } else if(scratch.Right == 0 && scratch.Left == 0) {
            printf("AGS: Right==0 && Left==0 && P[%d].AGS_KernelRadius=%g (task=%d)\n", i, current_h, ThisTask); fflush(stdout);
            endrun(90001010);
            return IterResult{IterStatus::Converged, current_h};   /* graceful: bad-stop set; stop iterating this active with the last valid AGS radius; drains at runner completion -> phase poll (per-active: NO immediate collective) */
        } else if(scratch.Right == 0 && scratch.Left > 0) {
            /* No upper bound — geometric extrapolation (rkern.cc:336-364). */
            double fac_lim;
            if(NumNgb_norm > 1) fac_lim = log(desnumngb / NumNgb_norm) / NUMDIMS;
            else                 fac_lim = 1.4;
            if((NumNgb_norm < 2.0 * desnumngb) && (NumNgb_norm > 0.1 * desnumngb)) {
                double slope = P[i].DrkernNgbFactor;
                if(iter > 2 && slope < 1) slope = 0.5 * (slope + 1);
                double fac = fac_lim * slope;
                if(iter >= 4 && P[i].DrkernNgbFactor == 1.0) fac *= 10;
                if(fac < fac_lim + 0.231) new_h = current_h * exp(fac);
                else                       new_h = current_h * exp(fac_lim + 0.231);
            } else {
                new_h = current_h * exp(fac_lim);
            }
        } else /* Right > 0 && Left == 0 */ {
            /* No lower bound — geometric extrapolation (rkern.cc:366-393). */
            double fac_lim;
            if(NumNgb_norm > 1) fac_lim = log(desnumngb / NumNgb_norm) / NUMDIMS;
            else                 fac_lim = 1.4;
            if(fac_lim < -1.535) fac_lim = -1.535;
            if((NumNgb_norm < 2.0 * desnumngb) && (NumNgb_norm > 0.1 * desnumngb)) {
                double slope = P[i].DrkernNgbFactor;
                if(iter > 2 && slope < 1) slope = 0.5 * (slope + 1);
                double fac = fac_lim * slope;
                if(iter >= 10 && P[i].DrkernNgbFactor == 1.0) fac *= 10;
                if(fac > fac_lim - 0.231) new_h = current_h * exp(fac);
                else                       new_h = current_h * exp(fac_lim - 0.231);
            } else {
                new_h = current_h * exp(fac_lim);
            }
        }
    }
    /* else: set_to_max/minrkern path already set new_h above (maxsoft / minsoft) */

    /* Min/max clamps (rkern.cc:397-400). */
    if(new_h < minsoft) new_h = minsoft;
    if(scratch.set_to_minrkern == 1) new_h = minsoft;
    if(new_h > maxsoft) new_h = maxsoft;
    if(scratch.set_to_maxrkern == 1) new_h = maxsoft;

    /* Sync P[i] to the radius the runner will use next iter.
     * Without this, P[i].AGS_KernelRadius would hold the
     * iter-K value while radii_uvm[sg][slot] holds iter-(K+1) — and the
     * post-runner finalize pass would see the stale K value. */
    P[i].AGS_KernelRadius = (MyFloat)new_h;

#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
    if((long long)P[i].ID == (long long)GIZMO_NLR_AGS_DEBUG_ID_TRACE) {
        printf("[AGS_TRACE rank=%d after_iter AdjustRadius] ID=%lld iter=%d "
               "current_h=%g new_h=%g NumNgb_raw=%g NumNgb_norm=%g "
               "Left=%g Right=%g minsoft=%g maxsoft=%g\n",
               ThisTask, (long long)P[i].ID, iter,
               current_h, new_h, NumNgb_raw, NumNgb_norm,
               scratch.Left, scratch.Right, minsoft, maxsoft);
        fflush(stdout);
    }
#endif

    return IterResult{IterStatus::AdjustRadius, new_h};
}

/* after_iter_global — host-only, no MPI (TRAP 7).
 *
 * Legacy `ags_density` printed "AGS-ngb iteration N: need to repeat for K
 * particles" at iter > 10 (rkern.cc:413). The runner calls
 * after_iter_global BEFORE the per-subgroup MPI_Allreduce,
 * so `drv.global_active_total` here is the PREVIOUS iter's count, not the
 * current iter's. Rather than print a misleading number or duplicate
 * the Allreduce here (would be a TRAP-7 violation), we drop the count
 * from the message and just emit a "still iterating" marker. The legacy
 * exact count is lost from the AGS port's per-iter log; not load-bearing
 * for physics. */
void AgsDensitySpec::after_iter_global(const neighbor_loop_args& /*args*/,
                                        const NlrIterDriver<AgsDensitySpec>& drv)
{
    const int iter = drv.iter_index;
    if(iter > 10 && ThisTask == 0) {
        printf("AGS-ngb iteration %d: still iterating.\n", iter);
    }
}

/* ============================================================================
 * GHOST WRITEBACK MANIFEST (first iterative ghost-writeback user).
 *
 * Single op: PARTICLE_MAX(wakeup). Operates on positive-int wakeup values
 * (TimeBin+1 hydro convention) after the writer fix —
 * cross-rank propagation now correctly delivers the largest TimeBin+1 to
 * the home rank's P[j].wakeup.
 *
 * The bundle only applies P[j].wakeup on the home rank (PARTICLE_MAX merge);
 * it does NOT raise NeedToWakeupParticles_local. The global wakeup flag is
 * raised from the sticky need_wakeup_uvm counter on any rank that generated
 * a wakeup inside the pair body (see cleanup_device_context), and
 * timestep.cc's MPI_Allreduce then makes all ranks process the reverse-comm'd
 * P[j].wakeup. (Previously the now-deleted ghost_writeback_wakeup helper set
 * the global flag itself; the runner's bundle does not have that side
 * effect.)
 * ========================================================================== */
GHOST_WRITEBACK_BUNDLE_BEGIN(ags_density)
    GHOST_WRITEBACK_PARTICLE_MAX(wakeup)
GHOST_WRITEBACK_BUNDLE_END(ags_density)

/* ghost_write_detector_begin/end: runner default (loop_name = "ags_density"). */

void AgsDensitySpec::ghost_writeback_begin(const neighbor_loop_args& /*args*/,
                                            const NeighborLoopPlan& /*plan*/)
{
    ghost_writeback_begin_bundle(ags_density_ghost_writeback_bundle_ptr());
}

void AgsDensitySpec::ghost_writeback_end(const neighbor_loop_args& /*args*/,
                                          const NeighborLoopPlan& /*plan*/)
{
    ghost_writeback_end_bundle(ags_density_ghost_writeback_bundle_ptr());
}

/* AGS field parity was validated by the retired two-binary route: the runner
 * build against a post-fix legacy build. */

/* ============================================================================
 * ags_density() — runner-driven caller surface.
 *
 * Replaces the legacy ags_density() body that used to live in
 * gravity/ags_rkern.cc; the legacy body and its GIZMO_NLR_AGSDENSITY_USE_LEGACY
 * two-binary-parity compile gate were retired after parity passed on Vista.
 *
 *   1. Per-active pre-loop init (legacy rkern.cc:94-98): for each
 *      ags-active i, capture entry-time AGS_KernelRadius into the
 *      AGS_Prev[] array for the final-final pass's minsoft/maxsoft
 *      clamp; reset P[i].wakeup=0 and P[i].AGS_vsig=0.
 *
 *   2. Build per-bm subgroup partition (legacy rkern.cc:112-128):
 *      ags_gravity_kernel_shared_BITFLAG(P[i].Type) is a pure function of
 *      Type, so each active belongs to exactly one bm group (the runner's
 *      partition assertion checks this every iter under DEBUG /
 *      GIZMO_NLR_ASSERT_PARTITION). Allreduce the bm-presence mask so all
 *      ranks see the same subgroup ordering (empty-on-this-rank subgroups
 *      get nullptr active_indices + num_active_local=0).
 *
 *   3. Ghost lifecycle is the runner's: its Mode A import is request-driven
 *      for this loop's own sources and types, with per-iter regrow inside
 *      `rebuild_mode_a_arena_and_ctx_for_current_active_union`; Mode B and
 *      Mode D import nothing.
 *
 *   4. run_neighbor_loop_iterative<AgsDensitySpec>(args) — owns the iter
 *      loop, per-iter pair_kernel dispatch, per-iter after_iter call which
 *      writes post-processed values + new AGS_KernelRadius to P[i].
 *
 *   5. Final-final pass (legacy rkern.cc:431-453, verbatim): AGS_zeta
 *      normalization + NumNgb cube-root. Uses AGS_Prev[i] for the
 *      minsoft/maxsoft clamp; reads/writes P[i] directly.
 *
 *   6. Timing accounting (legacy rkern.cc:456-458).
 *
 * No outer do-while, no per-iter ghost_writeback_zero_wakeup /
 * ghost_writeback_wakeup ladder (the runner's bundle + Spec ghost_writeback
 * hooks subsume it; the legacy wakeup helpers were retired in the cleanup).
 * ========================================================================== */

void ags_density(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    double t00_truestart = my_second(); double child0_span = CPU_ChildCharged;

    /* Canonical host All accessor — stylistic intent-tag for host-snapshot
     * reads in this scope. The redirect is device-pass-only,
     * so bare All.* would also be correct here; the accessor documents that
     * these reads run on the host pass. Used at three sites below outside
     * Spec::populate_call_scalars. */
    const struct global_data_all_processes *host_all = nlr_host_all_ptr();

    /* (1) AGS_Prev[] alloc + per-active pre-loop init. AGS_Prev[] is used
     * by the final-final pass; runner's per-active IterScratch.AGS_Prev
     * handles the per-iter clamp internally. */
    MyFloat *AGS_Prev = (MyFloat *) mymalloc("AGS_Prev", NumPart * sizeof(MyFloat));
    for (int i : ActiveParticleList) {
        if(ags_density_isactive(i)) {
            AGS_Prev[i]     = P[i].AGS_KernelRadius;
            P[i].AGS_vsig   = 0;
            P[i].wakeup     = 0;
#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
            ags_debug_trace_id_state("pre_solve_init", i);
#endif
        }
    }

    /* Pre-solve hard guard. If any ags-active particle
     * enters with AGS_KernelRadius <= 0, the convergence test will
     * spuriously "converge" at h=0 (maxsoft clamp via AGS_Prev*AGS_DSOFT_TOL
     * collapses), leaving the particle with soft=0 → endrun(888) in
     * get_timestep. Either init.cc didn't seed properly, or some upstream
     * path corrupted the radius. Fail loudly with full state to surface
     * the root class. Debug-gated; remove (or convert to assert) once
     * the regression is fixed. */
#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
    for (int i : ActiveParticleList) {
        if(ags_density_isactive(i) && !(P[i].AGS_KernelRadius > 0)) {
            fprintf(stderr,
                "[AGS_GUARD rank=%d] FATAL: ags-active particle entered "
                "ags_density() with AGS_KernelRadius <= 0. "
                "ID=%lld Type=%d Mass=%g TimeBin=%d "
                "AGS_KernelRadius=%g KernelRadius=%g ForceSoftening=%g "
                "minsoft=%g maxsoft=%g\n",
                ThisTask, (long long)P[i].ID, (int)P[i].Type,
                (double)P[i].Mass, (int)P[i].TimeBin,
                (double)P[i].AGS_KernelRadius, (double)P[i].KernelRadius,
                (double)P[i].ForceSoftening,
                ags_return_minsoft(i), ags_return_maxsoft(i));
            fflush(stderr);
            endrun(81360);
        }
    }
#endif

    /* (2) Build per-bm subgroup partition. */
    std::map<int, std::vector<int>> bm_groups_host;
    uint64_t local_bm_presence = 0;
    for (int i : ActiveParticleList) {
        if(ags_density_isactive(i)) {
            int bm = ags_gravity_kernel_shared_BITFLAG(P[i].Type);
            if(bm > 0 && bm < 64) {
                bm_groups_host[bm].push_back(i);
                local_bm_presence |= (1ULL << bm);
            }
        }
    }
    uint64_t global_bm_presence = local_bm_presence;
    if(NTask > 1) {
        MPI_Allreduce(&local_bm_presence, &global_bm_presence, 1,
                      MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
    }

    /* Build NlrSubgroup array. ALL ranks enter with the same subgroup
     * ordering (per global_bm_presence); ranks with no actives for a
     * particular bm contribute num_active_local=0 + null active_indices. */
    std::vector<NlrSubgroup> subgroups;
    std::vector<std::vector<int>> subgroup_actives_storage;
    int total_active_local = 0;
    for(int bm = 1; bm < 64; bm++) {
        if(!(global_bm_presence & (1ULL << bm))) continue;
        std::vector<int> actives_for_bm;
        auto it = bm_groups_host.find(bm);
        if(it != bm_groups_host.end()) actives_for_bm = std::move(it->second);
        subgroup_actives_storage.push_back(std::move(actives_for_bm));
        NlrSubgroup sg{};
        sg.j_type_bitmask   = (unsigned int)bm;
        sg.num_active_local = (int)subgroup_actives_storage.back().size();
        sg.active_indices   = (sg.num_active_local > 0)
                              ? subgroup_actives_storage.back().data()
                              : nullptr;
        subgroups.push_back(sg);
        total_active_local += sg.num_active_local;
    }

    if(subgroups.empty()) {
        /* No globally-active AGS particles. Skip the runner call entirely
         * (runner asserts num_subgroups >= 1 at entry). */
        myfree(AGS_Prev);
        double t1 = my_second(); cpu_chain_sync(t1);
        CPU_Step[CPU_AGSDENSMISC] += cpu_minus_children(timediff(t00_truestart, t1), child0_span);
        return;
    }

    /* Concatenated active list across subgroups (runner's args.active_list).
     * Each active particle appears in exactly one subgroup, so the
     * concatenation has no duplicates. */
    std::vector<int> active_list_concat;
    active_list_concat.reserve(total_active_local);
    for(const auto& sg_actives : subgroup_actives_storage) {
        for(int i : sg_actives) active_list_concat.push_back(i);
    }

    /* (3) The runner owns the ghost lifecycle: on a path that imports, it imports
       exactly this loop's neighbours for exactly these sources, and the other paths
       import nothing.  An all-types import here would only be torn down again. */
    double ags_ghost_safety = gizmo_ghost_safety_factor();

    /* (4) Build iterative args + drive the runner. */
    AgsDensitySpec::Aux aux{};                    /* empty per design v0.4.3 */

    neighbor_loop_args_iterative args{};
    static_cast<neighbor_loop_args&>(args) = nlr_default_args();
    args.P                   = P;
    args.CellP               = (host_all->TotN_gas > 0) ? CellP : nullptr;
    args.num_total           = NumPart;
    args.active_list         = (total_active_local > 0)
                               ? active_list_concat.data() : nullptr;
    args.num_active          = total_active_local;
    args.aux                 = &aux;
    args.num_subgroups       = (int)subgroups.size();
    args.subgroups           = subgroups.data();
    args.ghost_safety_factor = ags_ghost_safety;

    /* The runner span is the loop's compute; everything around it is the
     * caller's own work and stays in the misc bucket. Same split density()
     * reports. */
    const double t_runner_start = my_second(); const double child0_runner = CPU_ChildCharged;
    run_neighbor_loop_iterative<AgsDensitySpec>(args);
    const double t_runner = cpu_minus_children(timediff(t_runner_start, my_second()), child0_runner);
    CPU_Step[CPU_AGSDENSCOMPUTE] += t_runner;

#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
    for (int i : ActiveParticleList) ags_debug_trace_id_state("post_runner_pre_final", i);
#endif

    /* Post-runner cleanup. */
    if(NTask > 1) ghost_exchange_cleanup();

    /* (5) Final-final pass — AGS_zeta normalization + NumNgb cube-root.
     * Verbatim from legacy rkern.cc:431-453. */
    for (int i : ActiveParticleList) {
        if(ags_density_isactive(i)) {
            if((P[i].Mass > 0) && (P[i].AGS_KernelRadius > 0) && (P[i].NumNgb > 0)) {
                double minsoft = ags_return_minsoft(i);
                double maxsoft = ags_return_maxsoft(i);
                minsoft = DMAX(minsoft, AGS_Prev[i] * AGS_DSOFT_TOL);
                maxsoft = DMIN(maxsoft, AGS_Prev[i] / AGS_DSOFT_TOL);
                if(P[i].AGS_KernelRadius >= maxsoft) { P[i].AGS_zeta = 0; }

                double z0 = 0.5 * P[i].AGS_zeta * P[i].AGS_KernelRadius /
                            (NUMDIMS * P[i].Mass * P[i].NumNgb /
                             (VOLUME_NORM_COEFF_FOR_NDIMS *
                              pow(P[i].AGS_KernelRadius, NUMDIMS)));
                double h_eff = 2.0 * (KERNEL_CORE_SIZE *
                                      host_all->ForceSoftening[P[i].Type]);
                double Prho = 0 * h_eff * h_eff / 2.0;
                if(P[i].Particle_DivVel > 0) Prho = -Prho;
                P[i].AGS_zeta = P[i].Mass * P[i].Mass * P[i].DrkernNgbFactor *
                                (z0 + Prho);
                P[i].NumNgb = pow(P[i].NumNgb, 1.0 / NUMDIMS);
            } else {
                P[i].AGS_zeta         = 0;
                P[i].NumNgb           = 0;
                P[i].AGS_KernelRadius = host_all->ForceSoftening[P[i].Type];
            }
        }
    }
#ifdef GIZMO_NLR_AGS_DEBUG_ID_TRACE
    for (int i : ActiveParticleList) ags_debug_trace_id_state("post_final_final", i);
#endif

    myfree(AGS_Prev);

    /* (6) Timing accounting. */
    double t1 = my_second(); cpu_chain_sync(t1);
    double timeall = cpu_minus_children(timediff(t00_truestart, t1), child0_span);
    /* The caller's own share: the whole call less the runner span charged above. */
    CPU_Step[CPU_AGSDENSMISC] += timeall - t_runner;
}

#endif /* AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE */
