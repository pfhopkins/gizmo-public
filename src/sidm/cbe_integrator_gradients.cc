/*
 * sidm/cbe_integrator_gradients.cc
 *
 * CBE pre-force gradient module. Out-of-line host hooks for
 * CBEGradSpec plus the toplevel driver CBEGrad_gradient_calc(). Mirrors
 * sidm/dm_fuzzy_loop.cc (DMGradSpec / DMGrad_gradient_calc) — same shape,
 * two passes orchestrated at the toplevel via Aux::loop_iteration.
 *
 * Persistent gradient storage lives on P[i].Gradients_CBE_basis_moments
 * (declarations/particle_data.h, gated on CBE_INTEGRATOR_WITHGRADIENTS).
 * Standard P[]-driven ghost import (gizmo_request_filtered_ghost_import_fresh)
 * carries the field naturally; there is no scratch UVM array and no custom
 * Alltoallv. The earlier scratch architecture (CbeGradScratch,
 * CbeGradScratchOwner, cbe_grad_import_ghosts, CbeGradientsSpec,
 * CbeBjLimiterSpec) is gone.
 *
 * Pass-0 writeback solves M^{-1} . B with the same Tikhonov-style
 * ill-conditioning guard used by hydro NV_T inversion
 * (hydro/density_loop.cc:514-525). Pass-1 writeback rescales the
 * already-persistent field in place by phi[m][k] — it MUST NOT re-zero or
 * re-compute the field.
 *
 * Primitive-gradient swap: the stored field now holds primitive gradients
 * (∂ρ, ∂v, ∂S) rather than moment gradients (∂m, ∂p, ∂T), packed in the
 * same slot layout. The writeback math is unchanged (M^{-1} · B, then
 * per-(m,k) rescale by phi) — only the upstream LSQ accumulator B is now
 * built from primitive deltas (see cbe_integrator_gradients.h pair
 * bodies). The field name `Gradients_CBE_basis_moments` is now STALE;
 * rename out of scope.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"   /* MUST precede allvars.h */
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "cbe_integrator_gradients.h"

#if defined(CBE_INTEGRATOR_WITHGRADIENTS)

/* ============================================================================
 * Out-of-line CBEGradSpec hooks.
 * ========================================================================== */

/* search_radius — mirrors AgsForceSpec so the symmetric face set the pass-1
 * limiter walks is a SUPERSET of the AgsForce flux pair set. The pass-0 LSQ
 * pair body internally re-narrows acceptance to r < h_i. */
double CBEGradSpec::search_radius(const neighbor_loop_args& /*args*/,
                                   int /*active_slot*/, int i)
{
    return (double)P[i].AGS_KernelRadius;
}

CBEGradSpec::CallScalars
CBEGradSpec::populate_call_scalars(const neighbor_loop_args& args)
{
    CallScalars s;
    s.common         = nlr_common_scalars_from_all();
    s.loop_iteration = static_cast<const Aux*>(args.aux)->loop_iteration;
    return s;
}

/* apply_active_writeback — pass-gated:
 *   pass 0 (LSQ):       solve M^{-1} . B with Tikhonov-style diagonal-loading
 *                       guard; write P[i].Gradients_CBE_basis_moments.
 *   pass 1 (BJ-style):  rescale the already-persistent gradient by
 *                       phi[m][k] IN PLACE (no re-zero, no recompute). */
void CBEGradSpec::apply_active_writeback(const neighbor_loop_args& args,
                                          int active_slot, int i,
                                          const AccumData& accum)
{
    /* These should be impossible under the runner contract — args.aux is
     * set by CBEGrad_gradient_calc before every per-bm dispatch, and the
     * runner only invokes apply_active_writeback for indices it built into
     * the active list (which is bounds-clamped in toplevel). A failure
     * here is a runner / toplevel bug, not a recoverable physics state;
     * loudly abort rather than silently skipping a gradient. */
    const Aux *aux = static_cast<const Aux*>(args.aux);
    if(aux == nullptr) {
        fprintf(stderr,
                "[CBEGradSpec::apply_active_writeback] args.aux == nullptr "
                "(active_slot=%d, i=%d, NumPart=%d, rank=%d) — runner contract violation\n",
                active_slot, i, NumPart, ThisTask);
        endrun(91420);
        return;   /* soft stop: aux is NULL; returning avoids the aux->loop_iteration deref below */
    }
    if(i < 0 || i >= NumPart) {
        fprintf(stderr,
                "[CBEGradSpec::apply_active_writeback] i=%d out of [0,%d) "
                "(active_slot=%d, rank=%d) — runner handed out-of-range active index\n",
                i, NumPart, active_slot, ThisTask);
        endrun(91421);
        return;   /* soft stop: i is out of range; returning avoids the args.P[i] OOB reference below */
    }

    struct particle_data& Pi = args.P[i];

    if(aux->loop_iteration <= 0) {
        /* ------------------------------------------------------------------
         * Pass 0 — LSQ inversion with Tikhonov-style diagonal-loading guard.
         *
         * Same pattern as hydro NV_T inversion (hydro/density_loop.cc:514-525):
         *   threshold = 10 * CONDITION_NUMBER_DANGER
         *   while cond > threshold:
         *     M += (1.05 * trace/NUMDIMS / threshold) * I on diag
         *     cond_term *= 1.2
         * 50-iter hard cap; if it falls through, zero this row so commit #5's
         * flux body cannot ingest an arbitrarily-large gradient.
         * ------------------------------------------------------------------ */
        double M[3][3];
        double M_inv[3][3];
        for(int a = 0; a < 3; a++)
            for(int b = 0; b < 3; b++)
                M[a][b] = accum.M[a][b];

        const double cond_threshold = 10.0 * CONDITION_NUMBER_DANGER;
        const double trace_initial  = M[0][0] + M[1][1] + M[2][2];
        double       cond_term      = 1.05 * (trace_initial / 3.0) / cond_threshold;
        bool         singular       = false;

        for(int it = 0; it < 50; it++) {
            const double cond = matrix_invert_ndims(M, M_inv);
            if(isfinite(cond) && cond < cond_threshold) {
                bool any_bad = false;
                for(int a = 0; a < 3 && !any_bad; a++)
                    for(int b = 0; b < 3 && !any_bad; b++)
                        if(!isfinite(M_inv[a][b])) any_bad = true;
                if(!any_bad) break;
            }
            if(!isfinite(cond_term) || cond_term <= 0) cond_term = MIN_REAL_NUMBER;
            for(int a = 0; a < 3; a++) M[a][a] += cond_term;
            cond_term *= 1.2;
            if(it == 49) singular = true;
        }

        if(singular) {
            for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++)
                for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++)
                    for(int d = 0; d < 3; d++)
                        Pi.Gradients_CBE_basis_moments[m][k][d] = 0.0;
            return;
        }

        /* grad[m][k][d] = sum_e M_inv[d][e] * B[m][k][e] */
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
                for(int d = 0; d < 3; d++) {
                    double g = 0.0;
                    for(int e = 0; e < 3; e++) {
                        g += M_inv[d][e] * accum.B[m][k][e];
                    }
                    Pi.Gradients_CBE_basis_moments[m][k][d] = g;
                }
            }
        }
    } else {
        /* ------------------------------------------------------------------
         * Pass 1 — pairwise BJ-style limiter rescale, IN PLACE.
         *
         * Pi.Gradients_CBE_basis_moments was written by pass 0 (and
         * already imported on ghosts via standard P[] ghost transport when
         * the pass-1 runner call ran). We multiply each (m,k,d) component
         * by phi[m][k] clipped to [0, 1] with non-finite -> 0. Also sanitize
         * the persistent value (non-finite -> 0) so a bad pass-0 entry does
         * not propagate into commit #5's flux body.
         * ------------------------------------------------------------------ */
        for(int m = 0; m < CBE_INTEGRATOR_NBASIS; m++) {
            for(int k = 0; k < CBE_INTEGRATOR_NMOMENTS; k++) {
                double phi = accum.phi[m][k];
                if(!isfinite(phi)) phi = 0.0;
                else if(phi < 0.0) phi = 0.0;
                else if(phi > 1.0) phi = 1.0;
                for(int d = 0; d < 3; d++) {
                    const double g_old = Pi.Gradients_CBE_basis_moments[m][k][d];
                    double g_new = (isfinite(g_old) && isfinite(phi))
                                   ? (g_old * phi) : 0.0;
                    if(!isfinite(g_new)) g_new = 0.0;
                    Pi.Gradients_CBE_basis_moments[m][k][d] = g_new;
                }
            }
        }
    }
}


/* ============================================================================
 * CBEGrad_gradient_calc — toplevel. Two passes through the runner, separated
 * so pass 1's standard P[] ghost import sees the persistent
 * P[].Gradients_CBE_basis_moments written by pass 0.
 *
 * Active set: AGSForce_isactive (same predicate as the force consumer).
 * Subgroup partitioning: ags_gravity_kernel_shared_BITFLAG (same as
 *   AGSForce_calc, DMGrad_gradient_calc).
 *
 * Pre-zero: pass 0 only (LSQ accumulator scatter); pass 1 rescales in place.
 * Called from core/accel.cc in the pre-force phase, paralleling
 * DMGrad_gradient_calc.
 * ========================================================================== */
void CBEGrad_gradient_calc(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    const double t00 = my_second();  const double child0_span = CPU_ChildCharged;
    PRINT_STATUS(" ..calculating CBE basis-moment gradients\n");

    /* Partition local actives by shared AGS neighbor-type bitmask. bm==0
     * particles (not AGS-kernel-sharing) are dropped, matching DMGrad. */
    std::map<int, std::vector<int>> bm_groups;
    for(int i : ActiveParticleList) {
        if(AGSForce_isactive(i)) {
            const int bm = ags_gravity_kernel_shared_BITFLAG(P[i].Type);
            if(bm > 0 && bm < 64) { bm_groups[bm].push_back(i); }
        }
    }

    /* Global bm-presence union — every rank must iterate the same bm set in
     * the same order so the per-bm run_neighbor_loop calls (collective) stay
     * synchronized. Mirrors DMGrad_gradient_calc / ags_density_loop.cc. */
    uint64_t local_bm_presence = 0;
    for(const auto& kv : bm_groups) { local_bm_presence |= (1ULL << kv.first); }
    uint64_t global_bm_presence = local_bm_presence;
    if(NTask > 1) {
        MPI_Allreduce(&local_bm_presence, &global_bm_presence, 1,
                      MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
    }
    if(global_bm_presence == 0) {
        { const double t_end = my_second();
          CPU_Step[CPU_AGSDENSMISC] += cpu_minus_children(timediff(t00, t_end), child0_span);
          cpu_chain_sync(t_end); }
        return;   /* no AGSForce-active particles anywhere */
    }
    /* (CPU accounting bucket: parallel DMGrad_gradient_calc — both paths
     * land in CPU_AGSDENSMISC.) */

    /* Local union active list — used for the pass-0 pre-zero. */
    std::vector<int> union_actives;
    for(const auto& kv : bm_groups) {
        union_actives.insert(union_actives.end(), kv.second.begin(), kv.second.end());
    }

    /* Two passes: 0 = raw LSQ -> writes P[].Gradients_CBE_basis_moments;
     *             1 = pairwise BJ-style limiter -> rescales the same field
     *                 in place. Each pass is a fresh set of runner calls,
     *                 so pass 1's standard P[] ghost import sees pass-0's
     *                 freshly-written gradients on ghosts. */
    for(int pass = 0; pass < 2; pass++) {
        /* Pre-zero only on pass 0; pass 1 rescales the persistent field. */
        if(pass == 0) {
            for(int i : union_actives) {
                /* Snapshot last step's final gradient as the 'prev' (stale) row
                 * before the pre-zero wipes it — the density-continuity matching
                 * cost reconstructs face densities from this prior-step gradient. */
                std::memcpy(P[i].Gradients_CBE_basis_moments_prev,
                            P[i].Gradients_CBE_basis_moments,
                            sizeof(P[i].Gradients_CBE_basis_moments));
                std::memset(P[i].Gradients_CBE_basis_moments, 0,
                            sizeof(P[i].Gradients_CBE_basis_moments));
            }
        }

        /* One runner call per bm in the GLOBAL union, ascending order. A
         * rank with no local actives for a bm still enters with
         * num_active=0 — run_neighbor_loop is collective. */
        for(int bm = 1; bm < 64; bm++) {
            if(!(global_bm_presence & (1ULL << bm))) { continue; }
            auto it = bm_groups.find(bm);
            std::vector<int>* lst = (it != bm_groups.end()) ? &it->second : nullptr;

            CBEGradSpec::Aux aux;
            aux.loop_iteration = pass;

            neighbor_loop_args args = nlr_default_args();
            args.active_list = (lst && !lst->empty()) ? lst->data() : nullptr;
            args.num_active  = lst ? (int)lst->size() : 0;
            args.aux         = &aux;
            args.neighbor_type_mask_override = (unsigned int)bm;
            run_neighbor_loop<CBEGradSpec>(args);
        }
    }

    { const double t_end = my_second();
          CPU_Step[CPU_AGSDENSMISC] += cpu_minus_children(timediff(t00, t_end), child0_span);
          cpu_chain_sync(t_end); }
}

#endif /* CBE_INTEGRATOR_WITHGRADIENTS */
