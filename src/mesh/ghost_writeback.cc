/* ghost_writeback.cc — reverse communication of ghost particle modifications.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "ghost_writeback.h"
#include "gpu_neighbor_list.h"   /* the motion-target set the apply loop marks into */
#include "../system/gpu_particles_arena.h"
#include "../system/mpi_alltoallv_typed.h"


/* ============================================================================
 * Generic ghost-writeback scaffold (bundle-level exchange).
 *
 * begin_bundle: per-callback snapshot when num_ghosts>0; exactly-once
 *               arena_invalidate when num_ghosts>0. Strict no-op when
 *               n_callbacks==0 (no snapshot, no invalidate, no MPI, no
 *               cleanup).
 *
 * end_bundle:   exactly-once register_writeback (always, when bundle
 *               non-empty); BUNDLE-LEVEL single-pass count + pack +
 *               ONE Alltoall + ONE Alltoallv + apply (NTask>1 only);
 *               per-callback cleanup; exactly-once arena_invalidate
 *               (NTask>1 only). Strict no-op when n_callbacks==0.
 *
 * MPI-call invariant: total MPI calls per non-trivial end_bundle is O(1)
 * in the number of callbacks, NOT O(n_callbacks). All callbacks in a bundle
 * share one Alltoall (count exchange) and one gizmo_mpi_alltoallv_typed
 * (data exchange); the receiver demultiplexes by canonical callback order.
 * This makes it safe to put dozens-to-hundreds of writeback fields in one
 * Spec's manifest (mechFB, chemistry, etc.) without ballooning MPI cost.
 *
 * Per-rank send-buffer layout (rank-major, callback-minor):
 *
 *     [rank 0 segment][rank 1 segment]...[rank NTask-1 segment]
 *
 *     each rank-t segment = concatenation in canonical bundle order:
 *       [ cb0 records to t | cb1 records to t | ... | cb(N-1) records to t ]
 *
 * Both sender and receiver derive offsets from the per-callback per-rank
 * count matrix (exchanged by the count Alltoall). Canonical bundle order
 * is identical across ranks because every rank assembles the same
 * manifest; the per-process bundle pointer addresses differ but the
 * ordered list of callbacks they describe is the same.
 *
 * Mode B paths (local + remote): bundle hooks are NOT called by the runner
 * — j-side writes happen on the rank that owns j (Mode B local: same rank
 * as walker; Mode B remote: peer rank running the walker on its own pool),
 * so no reverse-comm is needed. The runner's nlr_path_uses_imported_ghosts
 * predicate gates the begin_bundle/end_bundle calls to Mode A only.
 *
 * Legacy invalidation timing preserved exactly:
 *   begin: arena_invalidate only when num_ghosts > 0
 *   end:   arena_invalidate only when NTask > 1
 *   end:   register_writeback always (when bundle non-empty)
 * ========================================================================== */

void ghost_writeback_begin_bundle(const struct ghost_writeback_bundle *bundle)
{
    if (!bundle || bundle->n_callbacks == 0) return;  /* strict no-op */
    int num_ghosts = ghost_get_num_ghosts();
    int num_local  = ghost_get_num_local();
    if (num_ghosts <= 0) return;                       /* legacy gate */
    for (int c = 0; c < bundle->n_callbacks; c++) {
        const ghost_writeback_callback *cb = bundle->callbacks[c];
        cb->snapshot(cb->ctx, num_ghosts, num_local);
    }
    gpu_particles_arena_invalidate();
}

void ghost_writeback_end_bundle(const struct ghost_writeback_bundle *bundle)
{
    if (!bundle || bundle->n_callbacks == 0) return;  /* strict no-op */

    /* Always register the writeback when the bundle is non-empty, even on
     * the NTask<=1 short-circuit, so the ghost-write detector counts it. */
    ghost_write_detector_register_writeback();

    int num_ghosts = ghost_get_num_ghosts();
    int num_local  = ghost_get_num_local();
    const int N = bundle->n_callbacks;

    if (NTask <= 1) {
        /* No reverse-comm needed; still let each callback free its snap. */
        for (int c = 0; c < N; c++) {
            const ghost_writeback_callback *cb = bundle->callbacks[c];
            cb->cleanup(cb->ctx);
        }
        return;  /* NO arena_invalidate on NTask<=1 (matches legacy). */
    }

    int *home_rank  = ghost_get_home_rank();
    int *home_index = ghost_get_home_index();
    (void)home_index;  /* used by callbacks; some need it indirectly */

    /* ---- Per-callback per-rank count matrix.
     * send_count[c * NTask + t] = number of cb-c records to send to rank t. */
    int *send_count = (int *) calloc((size_t)N * NTask, sizeof(int));

    /* ---- SINGLE COUNT PASS over ghosts (not per-callback). For each ghost
     * and each callback, ask whether it contributes; accumulate into the
     * count matrix. Hot path: N predicate calls per ghost. */
    for (int g = 0; g < num_ghosts; g++) {
        int t = home_rank[g];
        for (int c = 0; c < N; c++) {
            const ghost_writeback_callback *cb = bundle->callbacks[c];
            if (cb->delta_for_ghost(cb->ctx, g, num_local)) {
                send_count[c * NTask + t]++;
            }
        }
    }

    /* Local record count = total records this rank will send. Captured here
     * (before send_count is freed below) for the optional rank-0 traffic
     * print after the apply pass. Global SUM == global received (matter
     * conservation in the alltoallv), so one MPI_Reduce gives both numbers. */
    long long local_records_sent = 0;
    for (int c = 0; c < N; c++) {
        for (int t = 0; t < NTask; t++) {
            local_records_sent += send_count[c * NTask + t];
        }
    }

    /* ---- Per-rank send byte counts and displacements. */
    int *send_bytes     = (int *) calloc(NTask, sizeof(int));
    for (int c = 0; c < N; c++) {
        const size_t ds = bundle->callbacks[c]->delta_size;
        for (int t = 0; t < NTask; t++) {
            send_bytes[t] += send_count[c * NTask + t] * (int)ds;
        }
    }
    int *send_byte_disp = (int *) malloc(NTask * sizeof(int));
    send_byte_disp[0] = 0;
    for (int t = 1; t < NTask; t++) send_byte_disp[t] = send_byte_disp[t-1] + send_bytes[t-1];
    int total_send_bytes = send_byte_disp[NTask-1] + send_bytes[NTask-1];

    /* ---- Per-rank per-callback OFFSET within the rank's segment.
     * cb_off[c * NTask + t] = start of cb-c records inside rank-t's segment
     * (relative to send_byte_disp[t]). */
    int *cb_off = (int *) calloc((size_t)N * NTask, sizeof(int));
    for (int t = 0; t < NTask; t++) {
        int running = 0;
        for (int c = 0; c < N; c++) {
            cb_off[c * NTask + t] = running;
            running += send_count[c * NTask + t] * (int)bundle->callbacks[c]->delta_size;
        }
    }

    /* ---- Single send buffer, packed in canonical (rank, callback) order. */
    char *send_buf = (char *) malloc((size_t)(total_send_bytes > 0 ? total_send_bytes : 1));

    /* ---- SINGLE PACK PASS over ghosts. Per-(c, t) progress counters track
     * how many cb-c records have already been packed for rank t inside its
     * segment. */
    int *pack_progress = (int *) calloc((size_t)N * NTask, sizeof(int));
    for (int g = 0; g < num_ghosts; g++) {
        int t = home_rank[g];
        for (int c = 0; c < N; c++) {
            const ghost_writeback_callback *cb = bundle->callbacks[c];
            if (cb->delta_for_ghost(cb->ctx, g, num_local)) {
                const size_t ds = cb->delta_size;
                int p = pack_progress[c * NTask + t]++;
                size_t offset = (size_t)send_byte_disp[t]
                              + (size_t)cb_off[c * NTask + t]
                              + (size_t)p * ds;
                cb->pack_delta(cb->ctx, g, num_local, send_buf + offset);
            }
        }
    }
    free(pack_progress);
    free(cb_off);

    /* ---- ONE MPI_Alltoall on the count matrix. Each rank sends N integers
     * (one per callback) to every other rank, packaged as "what I'm going
     * to send to you, broken down by callback." Receiver gets the per-cb
     * counts of records arriving FROM each rank. */
    int *alltoall_send = (int *) malloc((size_t)NTask * N * sizeof(int));
    int *alltoall_recv = (int *) malloc((size_t)NTask * N * sizeof(int));
    for (int t = 0; t < NTask; t++) {
        for (int c = 0; c < N; c++) {
            alltoall_send[t * N + c] = send_count[c * NTask + t];
        }
    }
    MPI_Alltoall(alltoall_send, N, MPI_INT, alltoall_recv, N, MPI_INT, MPI_COMM_WORLD);
    /* Untranspose to recv_count[c * NTask + t] = # cb-c records arriving from rank t. */
    int *recv_count = (int *) calloc((size_t)N * NTask, sizeof(int));
    for (int t = 0; t < NTask; t++) {
        for (int c = 0; c < N; c++) {
            recv_count[c * NTask + t] = alltoall_recv[t * N + c];
        }
    }
    free(alltoall_send); free(alltoall_recv);

    /* ---- Per-rank recv byte counts and displacements. */
    int *recv_bytes     = (int *) calloc(NTask, sizeof(int));
    for (int c = 0; c < N; c++) {
        const size_t ds = bundle->callbacks[c]->delta_size;
        for (int t = 0; t < NTask; t++) {
            recv_bytes[t] += recv_count[c * NTask + t] * (int)ds;
        }
    }
    int *recv_byte_disp = (int *) malloc(NTask * sizeof(int));
    recv_byte_disp[0] = 0;
    for (int t = 1; t < NTask; t++) recv_byte_disp[t] = recv_byte_disp[t-1] + recv_bytes[t-1];
    int total_recv_bytes = recv_byte_disp[NTask-1] + recv_bytes[NTask-1];

    char *recv_buf = (char *) malloc((size_t)(total_recv_bytes > 0 ? total_recv_bytes : 1));

    /* ---- ONE gizmo_mpi_alltoallv_typed on the byte stream (size = 1). */
    gizmo_mpi_alltoallv_typed(send_buf, send_bytes, send_byte_disp,
                              recv_buf, recv_bytes, recv_byte_disp,
                              1, MPI_COMM_WORLD);
    free(send_buf); free(send_bytes); free(send_byte_disp); free(send_count);

    /* ---- SINGLE APPLY PASS over received bytes. For each rank-t segment
     * (in canonical rank order), demultiplex by callback in canonical
     * bundle order using the per-cb-per-rank recv_count. */
    for (int t = 0; t < NTask; t++) {
        size_t cb_base = (size_t)recv_byte_disp[t];
        for (int c = 0; c < N; c++) {
            const ghost_writeback_callback *cb = bundle->callbacks[c];
            const size_t ds = cb->delta_size;
            int n_records = recv_count[c * NTask + t];
            for (int r = 0; r < n_records; r++) {
                const void *rec = recv_buf + cb_base + (size_t)r * ds;
                cb->apply_delta(cb->ctx, rec);
                /* A loop that changes its neighbours' velocities has the runner
                 * arm the motion-target set for its writeback: the owner this
                 * delta landed on is then a particle whose bound must be raised.
                 * Every delta record begins with its home index (the callbacks'
                 * contract, ghost_writeback.h). */
                if(gx_motion_target_armed()) {gx_motion_target_mark_host(*(const int *)rec);}
            }
            cb_base += (size_t)n_records * ds;
        }
    }

    free(recv_buf); free(recv_bytes); free(recv_byte_disp); free(recv_count);

    /* Per-callback cleanup, then exactly-once arena_invalidate. */
    for (int c = 0; c < N; c++) {
        const ghost_writeback_callback *cb = bundle->callbacks[c];
        cb->cleanup(cb->ctx);
    }
    gpu_particles_arena_invalidate();

    /* Optional rank-0 traffic print. Backwards-compatible: bundles without a
     * loop_name (manually-constructed bundles) stay silent. Uses MPI-summed
     * global record count so rank 0 reports total traffic, not just its own
     * share — avoids false-fail on multi-rank MA-N validation when rank 0
     * is purely a receiver. */
    if (bundle->loop_name != nullptr) {
        long long global_records = 0;
        MPI_Reduce(&local_records_sent, &global_records, 1, MPI_LONG_LONG,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (ThisTask == 0 && global_records > 0) {
            printf("  Ghost writeback (%s): %lld records exchanged (global, sent==received)\n",
                   bundle->loop_name, global_records);
            fflush(stdout);
        }
    }
}


/* Compact delta struct for hydro j-writes.
 * Contains only the fields that hydro_force writes to j-particles. */
struct ghost_delta_hydro_t {
    int home_index;         /* P[]/CellP[] index on the home rank */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    MyDouble dMass;         /* additive: CellP[j].dMass */
#endif
    short int wakeup;       /* max: P[j].wakeup */
};


/* --- AGSForce variant (RETIRED) -----------------------------------------
 * The standalone ghost_writeback_{zero_,}agsforce path was retired when
 * AgsForceSpec migrated to the generic ghost_writeback bundle scaffold
 * (PARTICLE_ADD_VEC3 on Vel/dp, PARTICLE_ADD on NInteractions, PARTICLE_MAX
 * on wakeup) wired up directly in gravity/ags_force_loop.cc. Deleted once
 * migrated; no callers remained. */




/* --- SwallowTime variant (RETIRED) --------------------------------------
 * The ghost_writeback_{zero_,}swallowtime compatibility wrappers were
 * retired once SinkEnv1Spec migrated to the generic bundle scaffold:
 * sinks/sink_env1_loop.cc reverse-communicates SwallowTime via a
 * PARTICLE_MIN(particle_data::SwallowTime) manifest op. No callers
 * remained; the singleton-bundle wrappers and their ghost_writeback_ops.h
 * include were deleted with them. */


/* ThermalFB snapshot-based ghost_writeback_thermalfb / ghost_writeback_zero_thermalfb
 * retired; replaced by generic bundle in thermal_fb_loop.cc. */




/* --- MechFB variant (RETIRED) -------------------------------------------
 * The ghost_writeback_mechfb compatibility wrapper was retired once the
 * legacy galaxy_sf/mechanical_fb_gpu.cc evaluator was deleted: MechFBSpec
 * in galaxy_sf/mechfb_loop.cc owns the per-gas MechFBGasDelta reverse-comm
 * via its mechfb_writeback_detail custom generic-bundle callback. No
 * callers remained. */


    /* --- Grain backreaction variant (RETIRED) --------------------------------
     * The ghost_writeback_{zero_,}grainbackrx compatibility wrappers were
     * retired once GrainBackrxSpec migrated to the generic bundle scaffold:
     * solids/grain_physics_loop.cc reverse-communicates the grain_backrx
     * j-writes (Vel/VelPred/dp + Grain_AccelTimeMin min-update, plus the
     * GRAIN_EVOLUTION COND/SUBL fields) via the generic ghost-writeback
     * bundle. The multifluid GRAIN_FLUID-gate removal is captured upstream
     * via DO_FLUID_ALTSPECIES_DRAG_CALCULATION on grain_physics_loop.{cc,h}. */



/* RadFBRP variant retired.
 * RadFBRPSpec in galaxy_sf/radfb_rp_loop.{h,cc} now owns the j-side
 * reverse-comm via the generic ghost-writeback bundle (PARTICLE_ADD_VEC3
 * on Vel/dp + GAS_ADD_VEC3 on VelPred). */

/* RT source injection variant retired.
 * RtSrcInjectionSpec in
 * radiation/rt_source_injection_loop.{h,cc} now owns the j-side
 * reverse-comm via the generic ghost-writeback bundle. Three new
 * generic ops added at the same time to mesh/ghost_writeback_ops.h
 * (GAS_ADD_ARRAY, GAS_ADD_2D, GAS_ADD_VEC3_ARRAY) cover the field
 * shapes the legacy bespoke snapshot-diff handled: Rad_Je / Rad_E_gamma /
 * Rad_E_gamma_Pred (gas-side 1D), Rad_Intensity[/Pred] (gas-side 2D),
 * Rad_Flux[/Pred] (gas-side array-of-Vec3). */


/* --- Ghost-write detector (debug builds only) ----------------------------- */

#ifdef GIZMO_GPU_ARENA_DEBUG

static struct particle_data *gwd_P_snap = NULL;
static struct gas_cell_data *gwd_CellP_snap = NULL;
static int gwd_n_snap = 0;
static int gwd_local_at_snap = 0;
static int gwd_wb_count_at_snap = 0;
static int gwd_wb_count = 0;
static const char *gwd_kernel = "(none)";
static int gwd_active = 0;

void ghost_write_detector_register_writeback(void) { gwd_wb_count++; }

/* Re-snapshot ghost P/CellP at the current state. Called by gpu_ngb_list_build
 * after its lazy-drift loop updates Ti_current + Pos on touched ghost slots —
 * those updates are NOT a "kernel write" the detector should flag (they're
 * caller-side predicted-state setup), so we move the comparison baseline to
 * the post-drift state. The remaining detector window (kernel + scatter) then
 * catches any genuine kernel-side writes that need a writeback.
 *
 * CellP may be NULL when this fires for a DM/sidm/CBE-only run with no gas
 * particles allocated (AgsDensity opted into the detector via the default
 * runner hook, and AgsDensity is the first detector
 * client that can run with CellP unallocated). Skip the CellP snapshot in
 * that case; the gwd_CellP_snap pointer stays NULL from begin(). */
void ghost_write_detector_resnapshot_after_lazy_drift(void)
{
    if(!gwd_active) return;
    if(gwd_n_snap <= 0) return;
    int local = ghost_get_num_local();
    int n = ghost_get_num_ghosts();
    /* If layout changed mid-window we already can't compare — leave the
     * snapshot alone; the original alarm logic handles it. */
    if(local != gwd_local_at_snap || n != gwd_n_snap) return;
    memcpy(gwd_P_snap,    P     + local, n * sizeof(struct particle_data));
    if(CellP != NULL && gwd_CellP_snap != NULL) {
        memcpy(gwd_CellP_snap, CellP + local, n * sizeof(struct gas_cell_data));
    }
}

void ghost_write_detector_begin(const char *kernel_name)
{
    if(NTask <= 1) return;
    int n = ghost_get_num_ghosts();
    if(n <= 0) return;
    int local = ghost_get_num_local();
    if(gwd_active) {
        /* Nested begin without an end — programmer error. */
        fprintf(stderr, "[ghost_write_detector] nested begin('%s') while '%s' still active\n",
                kernel_name ? kernel_name : "(null)", gwd_kernel);
        endrun(91234);
    }
    /* new[], not malloc: particle_data is over-aligned (32 bytes) and malloc
     * cannot honour that. These cross function boundaries, so the LIFO arena is
     * not an option; new[] aligns an over-aligned type on its own. */
    gwd_P_snap = new struct particle_data[n];
    memcpy(gwd_P_snap, P + local, n * sizeof(struct particle_data));
    /* CellP=NULL guard: DM-only / sidm-only / CBE-only runs have no gas
     * particles and no CellP backing. The detector was added via the
     * default runner hook for AgsDensity, which is the first
     * detector client that can fire with CellP unallocated. Leave
     * gwd_CellP_snap=NULL in that case; end() skips the CellP memcmp. */
    if(CellP != NULL) {
        gwd_CellP_snap = new struct gas_cell_data[n];
        memcpy(gwd_CellP_snap, CellP + local, n * sizeof(struct gas_cell_data));
    } else {
        gwd_CellP_snap = NULL;
    }
    gwd_n_snap = n;
    gwd_local_at_snap = local;
    gwd_wb_count_at_snap = gwd_wb_count;
    gwd_kernel = kernel_name ? kernel_name : "(unnamed)";
    gwd_active = 1;
}

void ghost_write_detector_end(void)
{
    if(!gwd_active) return;
    int wb_ran = (gwd_wb_count > gwd_wb_count_at_snap);
    if(!wb_ran) {
        /* Use the local count captured at begin(); if it changed, ghost layout
         * was rebuilt mid-window — treat as a write-detection failure. */
        int local_now = ghost_get_num_local();
        int n_now = ghost_get_num_ghosts();
        int n_diff = 0, first_diff = -1;
        if(local_now == gwd_local_at_snap && n_now == gwd_n_snap) {
            /* CellP comparison skipped when CellP wasn't allocated at begin()
             * (DM/sidm/CBE-only runs); the gwd_CellP_snap pointer is NULL
             * in that case and there is nothing to diff. */
            int cellp_ok = (CellP != NULL && gwd_CellP_snap != NULL);
            for(int g = 0; g < gwd_n_snap; g++) {
                int p_diff = (memcmp(&P[gwd_local_at_snap + g], &gwd_P_snap[g], sizeof(struct particle_data)) != 0);
                int c_diff = cellp_ok &&
                             (memcmp(&CellP[gwd_local_at_snap + g], &gwd_CellP_snap[g], sizeof(struct gas_cell_data)) != 0);
                if(p_diff || c_diff) { n_diff++; if(first_diff < 0) first_diff = g; }
            }
        } else {
            n_diff = -1; /* layout changed — can't compare */
        }
        if(n_diff != 0) {
            fprintf(stderr, "[ghost_write_detector] kernel '%s' modified ghost particles "
                    "but no ghost_writeback_* was called.\n"
                    "  ghost_count_at_begin=%d, local_at_begin=%d, ghost_count_now=%d, local_now=%d, "
                    "first_differing_ghost_index=%d, n_diff=%d\n"
                    "  → either add a ghost_writeback call, or confirm this kernel is read-only on j.\n",
                    gwd_kernel, gwd_n_snap, gwd_local_at_snap, n_now, local_now,
                    first_diff, n_diff);
            /* Byte-level diff for the first differing slot — narrows down
             * which struct field the kernel mutated. */
            if(first_diff >= 0) {
                const unsigned char *cur_P = (const unsigned char *) &P[gwd_local_at_snap + first_diff];
                const unsigned char *snap_P = (const unsigned char *) &gwd_P_snap[first_diff];
                size_t sz = sizeof(struct particle_data);
                size_t first_byte = (size_t)-1, last_byte = 0, diff_bytes = 0;
                for(size_t b = 0; b < sz; b++) {
                    if(cur_P[b] != snap_P[b]) {
                        if(first_byte == (size_t)-1) first_byte = b;
                        last_byte = b; diff_bytes++;
                    }
                }
                fprintf(stderr, "  P[]: diff at bytes [%zu..%zu], %zu bytes total within sizeof(particle_data)=%zu\n",
                        first_byte, last_byte, diff_bytes, sz);
                if(first_byte != (size_t)-1) {
                    fprintf(stderr, "  P[].first_byte = 0x%02x → 0x%02x (snap → now)\n",
                            snap_P[first_byte], cur_P[first_byte]);
                }
                if(CellP != NULL && gwd_CellP_snap != NULL) {
                    const unsigned char *cur_C = (const unsigned char *) &CellP[gwd_local_at_snap + first_diff];
                    const unsigned char *snap_C = (const unsigned char *) &gwd_CellP_snap[first_diff];
                    size_t sc = sizeof(struct gas_cell_data);
                    size_t first_byteC = (size_t)-1, last_byteC = 0, diff_bytesC = 0;
                    for(size_t b = 0; b < sc; b++) {
                        if(cur_C[b] != snap_C[b]) {
                            if(first_byteC == (size_t)-1) first_byteC = b;
                            last_byteC = b; diff_bytesC++;
                        }
                    }
                    fprintf(stderr, "  CellP[]: diff at bytes [%zu..%zu], %zu bytes total within sizeof(gas_cell_data)=%zu\n",
                            first_byteC, last_byteC, diff_bytesC, sc);
                }
            }
            endrun(91234);
        }
    }
    delete[] gwd_P_snap;     gwd_P_snap = NULL;
    delete[] gwd_CellP_snap; gwd_CellP_snap = NULL;
    gwd_n_snap = 0;
    gwd_active = 0;
}

#endif /* GIZMO_GPU_ARENA_DEBUG */
