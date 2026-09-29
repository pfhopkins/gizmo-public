/* cbe_integrator_gpu.cc — GPU/OMP per-particle EP dispatch for the CBE
 * integrator: drift-kick (called from core/kicks.cc each half-step) and
 * post-gravity finalization (called from gravity/ags_force_loop.cc after the
 * AGSForce iterative neighbor loop closes).
 *
 * Both entries follow the same shape as solids/grain_drag.cc:
 *   - tiny-N (num_active < GPU_MIN_PARTICLES_FOR_OFFLOAD): direct host call
 *     in an OMP parallel-for loop. No arena, no Kokkos allocation, no
 *     compact gather. Each active index is unique → thread-safe.
 *   - large-N: compact gather of particle_data[num_active] in
 *     Kokkos shared-space, single kernel launch, narrow scatter of only
 *     the fields the kernel writes, then kokkos_free.
 *
 * The narrow scatter is load-bearing: it prevents the kernel from
 * inadvertently overwriting unrelated P[i] fields touched by other code
 * between the gather and the scatter.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"   /* MUST precede allvars.h */
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"

#include "../declarations/gpu_numeric_macros.h"
#include "../declarations/gpu_error_check.h"
#include "../declarations/gpu_dispatch_templates.h"
#include "../declarations/macros.h"
#include "../declarations/constants.h"
#include "../system/gpu_particles_arena.h"
#include "sidm_gpu_decls.h"

#if defined(CBE_INTEGRATOR)

#include "cbe_integrator_functions.h"

/* ----------------------------------------------------------------------------
 * cbe_drift_kick_evaluate_gpu — first/second half-step CBE drift-kick.
 *
 * Kernel: do_cbe_drift_kick_kernel(pi, dt, dT_out).
 * Writes: pi.CBE_basis_moments[NBASIS][NMOMENTS] AND pi.Vel[0..2], plus the
 *        predictor reset's pi.CBE_basis_moments_pred / pi.CBE_VelPred (pred =
 *        post-kick conserved state). The absolute-update round-trip derives
 *        V_new from the post-update mass-weighted-mean basis velocity and
 *        writes it to pi.Vel directly. The GPU narrow scatter therefore must
 *        copy all of these fields back from compact_P.
 * Reads: All.Time, All.TimeBegin, All.Ti_Current (via the kernel's basis-
 *        resplit branch); All-mirror handles the device read.
 *
 * The kernel returns a per-particle SPD-repair trace increment (sum over
 * basis of trace_after - trace_before) via *dT_out. Compact per-active
 * scratch dT_scratch[a] is host-summed and
 * forwarded to cbe_step_diagnostics_observe_repair (col-7/8). The whole
 * accumulator path is guarded on OUTPUT_ADDITIONAL_RUNINFO ||
 * CBE_INTEGRATOR_OUTPUT_MOREINFO so production builds pay zero overhead
 * (kernel gets nullptr; repair still happens, only the diagnostic
 * bookkeeping disappears). dP is identically 0 (SPD touches only stress).
 * --------------------------------------------------------------------------*/
void cbe_drift_kick_evaluate_gpu(struct particle_data *P_host,
                                 const int *active_indices, int num_active,
                                 const double *dt_host)
{
    if(num_active <= 0) return;

    if(num_active < GPU_MIN_PARTICLES_FOR_OFFLOAD)
    {
        /* Tiny-N OMP path: pure host calls, no arena, no Kokkos. Each
         * active index writes its own particle (unique by construction). */
        PRINT_STATUS("  CBE drift-kick (OMP): %d active", num_active);
#if defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(CBE_INTEGRATOR_OUTPUT_MOREINFO)
        double *dT_scratch = (double *) malloc(num_active * sizeof(double));
        if(!dT_scratch) {
            fprintf(stderr, "[task %d] cbe_drift_kick_evaluate_gpu: malloc(%zu) failed for dT_scratch (num_active=%d)\n",
                    ThisTask, (size_t)(num_active * sizeof(double)), num_active);
            endrun(91501);
        }
#else
        double *dT_scratch = nullptr;
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for(int a = 0; a < num_active; a++) {
            double dT_local = 0.0;
            do_cbe_drift_kick_kernel(P_host[active_indices[a]], dt_host[a],
                                     dT_scratch ? &dT_local : (double*)nullptr);
            if(dT_scratch) dT_scratch[a] = dT_local;
        }
#if defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(CBE_INTEGRATOR_OUTPUT_MOREINFO)
        if(dT_scratch) {   /* skip the diagnostic sum/free if the scratch alloc soft-failed above (endrun 91501) */
            double dT_sum = 0.0;
            for(int a = 0; a < num_active; a++) dT_sum += dT_scratch[a];
            cbe_step_diagnostics_observe_repair(/* dP */ 0.0, dT_sum);
            free(dT_scratch);
        }
#endif
        return;
    }

    /* Large-N GPU path: compact gather → kernel → narrow scatter. */
    GIZMO_GPU_ENSURE_ALL_FRESH();

    const size_t cbe_stage_bytes = (size_t) num_active * (sizeof(struct particle_data)
                                                          + sizeof(int) + sizeof(double));
    struct particle_data *compact_P = (struct particle_data *) gizmo_gpu_alloc_shared((size_t) num_active * sizeof(struct particle_data), NULL);
    int    *d_active = (int *)    gizmo_gpu_alloc_shared((size_t) num_active * sizeof(int), NULL);
    double *d_dt     = (double *) gizmo_gpu_alloc_shared((size_t) num_active * sizeof(double), NULL);
    if(!compact_P || !d_active || !d_dt) {
        if(d_dt)      {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_dt);}
        if(d_active)  {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_active);}
        if(compact_P) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(compact_P);}
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "cbe drift-kick: could not stage %d active particles (%.1f MB); "
                 "their distribution-function moments are not advanced",
                 num_active, (double) cbe_stage_bytes / (1024.0 * 1024.0));
        gizmo_request_controlled_stop(7716, msg, __FILE__, __LINE__, __FUNCTION__);
        return;
    }
#if defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(CBE_INTEGRATOR_OUTPUT_MOREINFO)
    /* Diagnostic only: the drift-kick runs without it, as the host path already
     * allows (its own scratch is optional there too). */
    double *dT_scratch = (double *) gizmo_gpu_alloc_shared((size_t) num_active * sizeof(double), NULL);
#else
    double *dT_scratch = nullptr;
#endif
    for(int a = 0; a < num_active; a++) compact_P[a] = P_host[active_indices[a]];
    memcpy(d_active, active_indices, num_active * sizeof(int));
    memcpy(d_dt,     dt_host,        num_active * sizeof(double));

    PRINT_STATUS("  CBE drift-kick (GPU): %d active", num_active);
    {
        struct particle_data *kp = compact_P;
        double *kdt = d_dt;
        double *kdT = dT_scratch;
        gizmo_gpu_kernel_launch("cbe_drift_kick", num_active, KOKKOS_LAMBDA(int a) {
            double dT_local = 0.0;
            do_cbe_drift_kick_kernel(kp[a], kdt[a],
                                     kdT ? &dT_local : (double*)nullptr);
            if(kdT) kdT[a] = dT_local;
        });
    }
    /* gizmo_gpu_kernel_launch fences before returning (the narrow scatter
     * below already relies on this); host reads of dT_scratch[a] follow. */

    /* Narrow scatter: kernel writes CBE_basis_moments AND pi.Vel (the latter
     * from the absolute-update round-trip's V_new = MMV derivation) with its
     * momentum change in pi.dp, AND
     * the predictor reset writes CBE_basis_moments_pred /
     * CBE_VelPred (pred = post-kick conserved state). All must be scattered
     * back to P_host — the OMP path updates in place, but the GPU path runs on
     * the compact copy, so dropping the pred fields here would leave host pred
     * stale on the large-N path. */
    for(int a = 0; a < num_active; a++) {
        int ii = active_indices[a];
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++)
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
                P_host[ii].CBE_basis_moments[m][k]      = compact_P[a].CBE_basis_moments[m][k];
                P_host[ii].CBE_basis_moments_pred[m][k] = compact_P[a].CBE_basis_moments_pred[m][k];
            }
        for(int k = 0; k < 3; k++) {
            P_host[ii].Vel[k]         = compact_P[a].Vel[k];
            P_host[ii].dp[k]          = compact_P[a].dp[k];
            P_host[ii].CBE_VelPred[k] = compact_P[a].CBE_VelPred[k];
        }
    }

#if defined(OUTPUT_ADDITIONAL_RUNINFO) || defined(CBE_INTEGRATOR_OUTPUT_MOREINFO)
    if(dT_scratch) {   /* absent when the diagnostic scratch could not be had */
        double dT_sum = 0.0;
        for(int a = 0; a < num_active; a++) dT_sum += dT_scratch[a];
        cbe_step_diagnostics_observe_repair(/* dP */ 0.0, dT_sum);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(dT_scratch);
    }
#endif
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_dt);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_active);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(compact_P);
}


/* ----------------------------------------------------------------------------
 * cbe_postgravity_evaluate_gpu — per-active post-AGSForce finalization.
 *
 * Kernel: do_cbe_postgravity_kernel(pi).
 * Writes: pi.CBE_basis_moments_dt[basis][0] only (mass-closure safety net).
 *         The CBE bulk-velocity injection into pi.GravAccel is gone in the
 *         current design — bulk V is now derived from the post-update
 *         mass-weighted-mean basis velocity inside the drift-kick absolute
 *         round-trip.
 * Reads: pi.Mass, pi.CBE_basis_moments_dt[basis][0], pi.CBE_basis_moments[basis][0].
 *
 * Called with active_indices = ActiveParticleList.data(),
 * num_active = ActiveParticleList.size() (no extra gating — matches the
 * legacy CPU loop in gravity/ags_force_loop.cc).
 * --------------------------------------------------------------------------*/
void cbe_postgravity_evaluate_gpu(struct particle_data *P_host,
                                  const int *active_indices, int num_active)
{
    if(num_active <= 0) return;

    if(num_active < GPU_MIN_PARTICLES_FOR_OFFLOAD)
    {
        PRINT_STATUS("  CBE postgravity (OMP): %d active", num_active);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for(int a = 0; a < num_active; a++)
            do_cbe_postgravity_kernel(P_host[active_indices[a]]);
        return;
    }

    GIZMO_GPU_ENSURE_ALL_FRESH();

    const size_t postgravity_stage_bytes = (size_t) num_active * sizeof(struct particle_data);
    struct particle_data *compact_P = (struct particle_data *) gizmo_gpu_alloc_shared(postgravity_stage_bytes, NULL);
    if(!compact_P) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "cbe post-gravity: could not stage %d active particles (%.1f MB); "
                 "the closure correction is not applied",
                 num_active, (double) postgravity_stage_bytes / (1024.0 * 1024.0));
        gizmo_request_controlled_stop(7716, msg, __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    for(int a = 0; a < num_active; a++) compact_P[a] = P_host[active_indices[a]];

    PRINT_STATUS("  CBE postgravity (GPU): %d active", num_active);
    {
        struct particle_data *kp = compact_P;
        gizmo_gpu_kernel_launch("cbe_postgravity", num_active, KOKKOS_LAMBDA(int a) {
            do_cbe_postgravity_kernel(kp[a]);
        });
    }

    /* Narrow scatter: kernel writes only the mass slot of the dt-accumulator
     * (closure-enforcing subtraction). Scatter that slot back; pi.GravAccel
     * is NOT touched anymore in the post-2026-06-04 design. */
    for(int a = 0; a < num_active; a++) {
        int ii = active_indices[a];
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++)
            P_host[ii].CBE_basis_moments_dt[m][0] = compact_P[a].CBE_basis_moments_dt[m][0];
    }

    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(compact_P);
}


#else /* stubs for builds without CBE_INTEGRATOR */

void cbe_drift_kick_evaluate_gpu(struct particle_data *, const int *, int, const double *) {}
void cbe_postgravity_evaluate_gpu(struct particle_data *, const int *, int) {}

#endif /* CBE_INTEGRATOR */
