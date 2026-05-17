/* mechanical_fb_gpu.cc — GPU neighbor-loop port for addFB_evaluate (B8).
 *
 * Dispatches all 6 modes (-2, -1, 0, 1, 2, 3) of the default-scheme
 * mechanical_fb on the GPU, sharing a single neighbor list across modes.
 * Per-gas deltas atomically accumulated into a MechFBGasDelta buffer; after
 * all modes complete, ghost deltas are reduced to home ranks via
 * ghost_writeback_mechfb; home-cell deltas are scattered by the existing
 * host-side verify_and_assign_local_mechfb_integrals (conservation-enforced).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../system/gpu_particles_arena.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"

#include "../declarations/gpu_numeric_macros.h"
#include "../declarations/gpu_error_check.h"
#include "../declarations/gpu_dispatch_templates.h"
#include "../declarations/macros.h"
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "../mesh/ghost_exchange_spec.h"  /* mech_fb_v1 inline spec literal */
#include "galsf_gpu_decls.h"

#if defined(GALSF_FB_MECHANICAL)

/* Legacy GPU evaluator: pulls in mechanical_fb_functions.h (MechFBCallScalars,
 * MechFBLocalIn/Out, mech_fb_* helper prototypes, inject_cosmic_rays_into_delta).
 * Does NOT include mechfb_loop.h — MechFBSpec / runner-template types are only
 * needed by the runner-port path. */
#include "mechanical_fb_functions.h"


/* mech_fb_local_fill, mech_fb_apply_aws_out, mech_fb_apply_source_mass_out
 * were promoted from `static` (TU-local) to public symbols in mechfb_loop.cc
 * (Phase 4 / Wave 3 / 3e.1 milestone 3) so the MechFBSpec runner-port's
 * after_iter_global + reset_per_iter_device_context can share the legacy
 * per-source pack + per-mode source-side host writes — single source of
 * truth across legacy + runner-port. Prototypes in mechanical_fb_functions.h. */


/* ================================================================
   GPU mechanical-feedback evaluator
   ----------------------------------------------------------------
   Kernel writes (host-visible, by index):
     i-side: out_arr[aa] (per-source MechFBOut: M_coupled, Area_weighted_sum)
     j-side (gas neighbors, indexed by j = neighbors[nn]; only Type==0):
       d_gas[j].N_injected, m_injected, TE_injected, KE_injected,
       d_gas[j].Z_injected[k], p_injected[k], Mass_Where_Dust_Shocked
                                           — Kokkos::atomic_add (all)
     P_gpu / CellP_gpu: NOT WRITTEN by the kernel.

   Sparse-scatter target: gas_delta_host[j] over touched j (walk neighbors[],
   filter Type==0). The full memcpy(P_host, P_gpu, num_all*...) and
   memcpy(CellP_host, CellP_gpu, num_all*...) currently at lines ~292-293
   are CARGO-CULT (kernel doesn't touch P/CellP). Source-side mass loss is
   already applied explicitly via mech_fb_apply_source_mass_out for active
   sources only. Delete the full memcpys; replace gas_delta full-N copy
   with sparse scatter.
   ================================================================ */
void mechanical_fb_evaluate_gpu(struct particle_data *P_host,
                                 struct gas_cell_data *CellP_host,
                                 int num_total,
                                 struct MechFBGasDelta *gas_delta_host,
                                 int n_gas,
                                 int *i_active_host, int num_active,
                                 const double *src_radii_host,
                                 int *n_couplings_out)
{
    ghost_write_detector_begin("mechanical_fb");
    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* Caller-supplies-active-list rule: i_active_host[] holds LOCAL indices
       built from ActiveParticleList + per-mode active-check (superset across
       modes). Never iterate num_total here — ghost imports would double-deposit.
       Multi-rank correctness: gizmo_density_prep_ghosts is an MPI collective so
       every rank must agree on whether to call it. MPI_Allreduce num_src first;
       if no rank has any source, every rank skips ghost_prep + the rest. */
    int num_src = num_active;
    int global_num_src = num_src;
    MPI_Allreduce(&num_src, &global_num_src, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(global_num_src == 0) { ghost_write_detector_end(); return; } /* nothing anywhere; ghost_writeback would write nothing */

    int imported_ghosts = 0;
    {   /* Collective import decision: active-aware exchange can leave rank A with 0
         * ghosts while rank B has some (e.g. after density redo on tiny-N).  Deciding
         * per-rank whether to call ghost_exchange causes a collective desync.  Use an
         * Allreduce so all ranks agree: if any rank needs ghosts, all re-import. */
        int need_import_local = (ghost_get_num_ghosts() <= 0) ? 1 : 0;
        int need_import = 0;
        MPI_Allreduce(&need_import_local, &need_import, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(need_import) {
            if(ghost_get_num_ghosts() > 0) ghost_exchange_cleanup();
            /* mech_fb_v1 caller-built spec literal (Bucket 1.2 PoC). The literal
             * here IS the SSOT for this loop's physics; ghost_exchange.cc has no
             * mech_fb-specific knowledge.
             *
             *   sources : caller-filtered i_active_host[] (per-mode active-check
             *             superset; type filter happens in addFB_evaluate_active_check
             *             upstream — no request_type_mask filter inside ghost_exchange)
             *   supply  : Type 0 cells only (kernel deposits energy/mass/momentum
             *             into the continuum/cell population)
             *   mode    : SYMMETRIC (paired r_ij<max(h_i,h_j) — matches the legacy
             *             ngb_treefind_pairs_threads in old GIZMO)
             *   h_q     : src_radii_host[] (caller-supplied per-source radii) */
            double t_drift = my_second();
            /* Codex 2026-05-12 v2: use the OUT-OF-LINE host accessor
             * `gizmo_host_ti_current()` from predict.cc — the inline
             * gizmo_gpu_host_all_ptr() in gpu_all_mirror.h proved
             * unreliable under nvcc_wrapper (dbg27). See
             * feedback_all_dev_trap_host_side.md. */
            move_particles(gizmo_host_ti_current());  /* lazy drift — was inside gizmo_density_prep_ghosts */
            CPU_Step[CPU_DENSMISC] += timediff(t_drift, my_second());
            int nq = num_active;
            double (*qpos)[3] = (double(*)[3]) malloc((size_t)(nq > 0 ? nq : 1) * sizeof(double[3]));
            double  *qh       = (double *)     malloc((size_t)(nq > 0 ? nq : 1) * sizeof(double));
            for(int a = 0; a < nq; a++) {
                int i = i_active_host[a];
                qpos[a][0] = P[i].Pos[0]; qpos[a][1] = P[i].Pos[1]; qpos[a][2] = P[i].Pos[2];
                qh[a]      = src_radii_host[a];
            }
            struct ghost_exchange_spec_t sp = {
                /* request_type_mask  */ 0u,
                /* supply_type_mask   */ GHOST_TYPE_0,
                /* search_mode        */ NGB_SEARCH_SYMMETRIC,
                /* safety_factor      */ gizmo_ghost_safety_factor(),
                /* caller_name        */ "mech_fb_v1",
                /* n_queries          */ nq,
                /* query_pos          */ (const double (*)[3]) qpos,
                /* query_h            */ qh
            };
            double t_ghost = my_second();
            ghost_exchange_run(&sp);
            CPU_Step[CPU_DENSCOMM] += timediff(t_ghost, my_second());
            free(qpos); free(qh);
            imported_ghosts = 1;
        }
    }

    int num_all = ghost_get_num_local() + ghost_get_num_ghosts();
    if(num_all <= 0) num_all = num_total;

    /* Wrapper fast-path: with no source particles, no kernel writes occur, so
     * we can skip arena_acquire (~1.4s slow-path memcpy when invalidated),
     * d_gas allocation+memset (~2.5 GB SharedSpace for num_all=12.4M), per-source
     * allocations, the gpu_ngb_list_build call, the per-mode kernel dispatch,
     * and the unconditional 17 GB host scatter at the end of the function.  The
     * ghost_writeback_mechfb collective is a no-op when NTask<=1 or num_ghosts<=0
     * (it returns at its own guard); on multi-rank with imported ghosts we
     * participate via a calloc'd zero buffer to keep MPI_Alltoallv consistent. */
    if(num_src == 0) {
        /* MULTI-RANK COLLECTIVE CORRECTNESS: ghost_writeback_mechfb is collective
         * (MPI_Alltoall + Alltoallv inside). Every rank that passes the global-
         * active-zero gate at line ~152 must enter the collective in lock-step,
         * EVEN if this rank has 0 local ghosts. With narrow-supply imports
         * (mech_fb_v1 supply=Type 0), ghost imports are asymmetric: rank-with-
         * sources gets imports back, rank-without doesn't. Skipping the call
         * on num_ghosts==0 alone deadlocked the rank-with-imports.  Drop the
         * `&& ghost_get_num_ghosts() > 0` from the gate; ghost_writeback_mechfb
         * now participates with empty buffers when num_ghosts==0 (line-876 fix). */
        if(NTask > 1) {
            struct MechFBGasDelta *zero_gas = (struct MechFBGasDelta *)
                calloc(num_all > 0 ? num_all : 1, sizeof(struct MechFBGasDelta));
            ghost_writeback_mechfb(zero_gas, gas_delta_host, n_gas);
            free(zero_gas);
        }
        ghost_write_detector_end();
        if(n_couplings_out) *n_couplings_out = 0;
        if(imported_ghosts) ghost_exchange_cleanup();
        return;
    }

    /* Step 13 Phase 1 arena. */
    gpu_particles_arena_set_site("mechanical_fb_evaluate_gpu");
    gpu_particles_arena_acquire(num_all, P_host, CellP_host);
    struct particle_data *P_gpu = gpu_particles_arena_P();
    struct gas_cell_data *CellP_gpu = gpu_particles_arena_CellP();

    /* Per-source input / output. 1-element backstop when num_src==0. */
    int alloc_n = (num_src > 0) ? num_src : 1;
    struct MechFBLocalIn *d_local = (struct MechFBLocalIn *)
        Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(struct MechFBLocalIn));
    struct MechFBOut *d_out = (struct MechFBOut *)
        Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(struct MechFBOut));
    int *d_active = (int *)
        Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(int));
    if(num_src > 0) memcpy(d_active, i_active_host, num_src * sizeof(int));

    /* Per-gas delta accumulator, sized for home + ghost cells. Zero before kernels. */
    struct MechFBGasDelta *d_gas = (struct MechFBGasDelta *)
        Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(num_all * sizeof(struct MechFBGasDelta));
    memset(d_gas, 0, num_all * sizeof(struct MechFBGasDelta));

    /* Build ONE neighbor list; reused across all 6 modes. */
    gpu_neighbor_list_t gnl;
    gpu_ngb_list_build(P_gpu, num_all,
                       i_active_host, num_src,
                       NGB_SEARCH_SYMMETRIC, 1 /* Type 0 only */,
                       &gnl, gpu_step_sidx_ptr(), 1.0, src_radii_host, NULL, "mech_fb");

    PRINT_STATUS("  GPU mech_fb: %d sources, %lld pairs", num_src, (long long)gnl.total_pairs);

    /* Codex r6 fix 2026-05-14: build the per-call MechFBCallScalars once
     * host-side via nlr_host_all_ptr(); pass into the lambda capture so the
     * refactored mechanical_fb_per_source_setup / pair_kernel /
     * inject_cosmic_rays_into_delta read scalars.* instead of bare All.*
     * (the per-TU All_dev mirror in this GPU TU is not a reliable substitute). */
    MechFBCallScalars mechfb_scalars;
    mechfb_fill_call_scalars(&mechfb_scalars);

    /* Per-mode dispatch. Modes beyond -2,-1,0 are only meaningful under
       GALSF_FB_FIRE_STELLAREVOLUTION (mass return, r-process, age tracers). */
    int modes[6] = {-2, -1, 0, 1, 2, 3};
    int num_modes = 3;
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION)
    num_modes = 4;  /* adds mode 1 */
#if defined(GALSF_FB_FIRE_RPROCESS)
    num_modes = 5;  /* adds mode 2 */
#endif
#if defined(GALSF_FB_FIRE_AGE_TRACERS)
    num_modes = 6;  /* adds mode 3 */
#endif
#endif

    for(int mi = 0; mi < num_modes; mi++) {
        int loop_iteration = modes[mi];

        /* Host: pack per-source input for this mode (reads P[i].Area_weighted_sum
           which was updated by the previous weighting pass via mech_fb_apply_aws_out). */
        for(int aa = 0; aa < num_src; aa++) {
            mech_fb_local_fill(i_active_host[aa], loop_iteration, &d_local[aa]);
        }

        /* Zero per-source output. */
        memset(d_out, 0, alloc_n * sizeof(struct MechFBOut));

        /* Launch kernel. */
        {
            int64_t *offsets   = gnl.offsets;
            int     *neighbors = gnl.neighbors;
            struct MechFBLocalIn *local_arr = d_local;
            struct MechFBOut     *out_arr   = d_out;
            int                  *active_arr = d_active;
            struct particle_data *kp = P_gpu;
            struct gas_cell_data *kc = CellP_gpu;
            struct MechFBGasDelta *kg = d_gas;

            Kokkos::parallel_for("mech_fb", num_src, KOKKOS_LAMBDA(int aa) {
                int i = active_arr[aa];
                /* per-mode active mask — superset list contains stars active in
                   ANY mode; skip those not active in this specific mode. */
                if(!mechanical_fb_star_active_check(i, loop_iteration, kp)) return;

                const struct MechFBLocalIn& loc = local_arr[aa];
                if(loc.Msne <= 0) return;
                if(loc.KernelRadius <= 0) return;

                struct MechFBSourceMode mstate;
                mechanical_fb_per_source_setup(loc, loop_iteration, mechfb_scalars, mstate);

                struct MechFBOut myout;
                myout.M_coupled = 0;
                for(int k = 0; k < AREA_WEIGHTED_SUM_ELEMENTS; k++) myout.Area_weighted_sum[k] = 0;

                int64_t start = offsets[aa], end = offsets[aa+1];
                for(int64_t nn = start; nn < end; nn++) {
                    int j = neighbors[nn];
                    if(kp[j].Type != 0) continue;
                    if(kp[j].Mass <= 0) continue;
                    Vec3<double> dp = loc.Pos - kp[j].Pos;
                    nearest_xyz(dp);
                    double r2 = dp.norm_sq();
                    if(r2 <= 0) continue;
                    double h2j = (double)kp[j].KernelRadius * (double)kp[j].KernelRadius;
                    if((r2 > mstate.h2) && (r2 > h2j)) continue;
                    if(r2 > mstate.r2max_phys) continue;
                    /* Legacy GPU evaluator: single-buffer d_gas covers home + ghost
                     * (num_local_gas = num_local_particles = num_all). Ghost ptr =
                     * nullptr (unreachable: j < num_local_gas always holds, so the
                     * first branch in mechfb_target_gas_delta fires for every j).
                     * oracle_dry_run = false (legacy has no oracle gate). scalars
                     * threaded by-value through the Kokkos lambda capture.
                     *
                     * num_local_particles added in milestone 3.5: in the runner
                     * port it splits non-gas locals from imported ghosts; here
                     * setting it to num_all alongside num_local_gas keeps legacy
                     * semantics identical (both thresholds equal → ghost branch
                     * unreachable). */
                    mechanical_fb_pair_kernel(loc, mstate, loop_iteration, j,
                                              kp, kc,
                                              mechfb_scalars,
                                              /*home_gas_delta     */ kg,
                                              /*ghost_gas_delta    */ (struct MechFBGasDelta*)nullptr,
                                              /*num_local_gas      */ num_all,
                                              /*num_local_particles*/ num_all,
                                              /*oracle_dry_run     */ false,
                                              dp, r2, myout);
                }
                out_arr[aa] = myout;
            });
        }
        Kokkos::fence();
        gizmo_gpu_check_last_error("mech_fb", num_src);

        /* Host: consume per-source output. For weighting modes (-2, -1), write
           Area_weighted_sum back into P[i] so the NEXT mode's local_fill sees
           the updated value. For coupling modes (>= 0), immediately apply
           source mass loss so later modes pack from the updated source state. */
        if(loop_iteration < 0) {
            for(int aa = 0; aa < num_src; aa++) {
                mech_fb_apply_aws_out(&d_out[aa], i_active_host[aa], loop_iteration);
            }
        } else {
            for(int aa = 0; aa < num_src; aa++) {
                int i = i_active_host[aa];
                MyFloat M_coupled = d_out[aa].M_coupled;
                if(M_coupled <= 0) continue;
                mech_fb_apply_source_mass_out(P_host, CellP_host, i, M_coupled);
                if(P_gpu != P_host || CellP_gpu != CellP_host) {
                    mech_fb_apply_source_mass_out(P_gpu, CellP_gpu, i, M_coupled);
                }
            }
        }
    } /* per-mode loop */

    /* Row 2 of arena-scope sweep: the former full memcpy(P_host, P_gpu, num_all*...)
     * and memcpy(CellP_host, CellP_gpu, num_all*...) here were CARGO-CULT — per
     * the kernel-writes audit (commit a06e30ca) the mech_fb pair kernel writes
     * only d_gas[j], never P/CellP. Source-side mass loss is already applied
     * explicitly to BOTH P_host and P_gpu via mech_fb_apply_source_mass_out
     * (lines ~301-304) for active sources only, so the host already has the
     * correct state. The full memcpys also masked a latent correctness bug:
     * they could overwrite newer host-side mutations with stale arena state.
     * DELETED. */

    /* Sparse gas_delta scatter: kernel writes d_gas[j] (atomic_add) only for
     * j's appearing in gnl.neighbors. Walk the CSR; copy the home-cell subset
     * (j < n_gas) into caller's gas_delta_host. Ghost-cell deltas (j >= n_gas)
     * stay in d_gas for ghost_writeback_mechfb to ship to the home rank.
     * Idempotent assignment makes duplicate j's correct.
     *
     * gnl.neighbors lives in DEVICE_SPACE (CudaSpace) — not host-readable on
     * GH200. Deep-copy to host once, then iterate the host buffer. */
    if(gnl.total_pairs > 0) {
        std::vector<int> gnl_neighbors_host((size_t)gnl.total_pairs);
        gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
        for(int64_t idx = 0; idx < gnl.total_pairs; idx++) {
            int j = gnl_neighbors_host[idx];
            if(j < n_gas) gas_delta_host[j] = d_gas[j];
        }
    }

    /* Ghost writeback: reduce ghost-side deltas into home-rank home_buf entries.
       d_gas is in SharedSpace (host-accessible), passed as the source; caller's
       gas_delta_host is the home destination. */
    ghost_writeback_mechfb(d_gas, gas_delta_host, n_gas);
    ghost_write_detector_end();

    /* Count gas cells that received coupling, matching N_Gas_Couplings_ThisTask.
     * TODO(arena-scope): this remains an O(n_gas) scan because ghost_writeback_mechfb
     * may merge remote contributions into home cells that weren't in our local
     * touched-list. A future change would have writeback report its touched
     * indices so this count can be sparse. The cost is small relative to the
     * deleted full-N P/CellP memcpy (~5GB), but it is the next refinement. */
    int n_coup = 0;
    for(int j = 0; j < n_gas; j++) if(gas_delta_host[j].N_injected > 0) n_coup++;
    if(n_couplings_out) *n_couplings_out = n_coup;

    /* No full P/CellP scatter remains; arena state is consistent with host
     * because (a) kernel doesn't write P/CellP and (b) per-source mass-loss
     * was applied to BOTH host and arena above. ghost_writeback_mechfb may
     * have merged into gas_delta_host (caller's buffer), which doesn't touch
     * the arena. */
    gpu_ngb_list_free(&gnl, gpu_step_sidx_ptr());
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_gas);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_active);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_out);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_local);
    /* Per-mode source updates above mutate host particle state directly, so
     * invalidate the arena before any later GPU consumer can assume freshness. */
    gpu_particles_arena_invalidate();

    if(imported_ghosts) ghost_exchange_cleanup();
}


#else

void mechanical_fb_evaluate_gpu(struct particle_data *p,
                                 struct gas_cell_data *cp,
                                 int num_total,
                                 struct MechFBGasDelta *gas_delta_host,
                                 int n_gas,
                                 int *i_active_host, int num_active,
                                 const double *src_radii_host,
                                 int *n_couplings_out)
{
    (void)p; (void)cp; (void)num_total; (void)gas_delta_host; (void)n_gas;
    (void)i_active_host; (void)num_active; (void)src_radii_host;
    if(n_couplings_out) *n_couplings_out = 0;
}

#endif /* GALSF_FB_MECHANICAL */
