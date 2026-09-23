/* gpu_neighbor_list.cc — GPU-accelerated neighbor list construction.
 *
 * Extracted from hydro/density_gpu.cc for reuse by any code that needs
 * a GPU-built CSR neighbor list (density, symmetric list, future loops).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <Kokkos_Core.hpp>

/* GPU All mirror: per-TU managed pointer to shared UVM allocation. */
#include "../declarations/gpu_all_mirror.h"
#include "../declarations/allvars.h"
#include "../system/gpu_particles_arena.h"
#include "../core/timestep_functions.h"   /* DriftKickTableView for the walk widening */

/* The drift-factor interpolator the walk's widening uses.
 *
 * Refreshed once per call rather than per node: it is the same table the node
 * sweep mirrors (gpu_force_drift.cc), and a walk that built its own would be a
 * second interpolator answering the same question. Null on failure, which simply
 * disables widening -- the walk then opens on the stored length, and the caller's
 * certification is what makes that legal. */
static struct DriftKickTableView g_walk_drift_tables;
static double *g_walk_drift_storage = NULL;
static int g_walk_drift_tables_ok = 0;

/* Refresh once per preparation. Shares drift_kick_table_mirror_refresh with the
   particle drift and the node sweep -- same table, same units, one owner. */
static void gx_walk_drift_tables_refresh(void)
{
    g_walk_drift_tables_ok =
        (drift_kick_table_mirror_refresh(&g_walk_drift_storage, &g_walk_drift_tables) == 0);
}
#include "../core/proto.h"
#include "../declarations/gpu_error_check.h"

#include "sfc_tiles.h"
#include "sfc_tiles_functions.h"
#include "gpu_neighbor_list.h"
#include "gpu_dirty_tracker.h"
#include "neighbor_list.h"
#include "ghost_writeback.h"  /* ghost_write_detector_resnapshot_after_lazy_drift, ghost_get_num_local */
#include "ghost_exchange_functions.h"   /* the shared accept predicates (canonical box wrap) */
#include "device_tree_walk.h"          /* the one device tree traversal */
#include "../gravity/forcetree.h"       /* BITFLAG_TOPLEVEL, force_host_lazy_drift_ti */
#include "../gravity/gpu_gravity_tree.h" /* node SoA + drift certification */

/* TILE_PERIODIC_X/Y/Z defined in sfc_tiles.h (included via gpu_neighbor_list.h) */

/* Persistent gas-only SIDX shared across density+symlist within a step.
 * Lifetime: built lazily on first gas (type_bitmask=1) ngb_list_build, reused
 * for all subsequent gas builds, freed by gpu_step_sidx_invalidate() after drift. */
static gpu_spatial_index_t g_step_sidx{};
static gpu_spatial_index_t g_step_sidx_alltypes{};

/* Lazy-drift h-slack: under lazy drift, neighbor j's KernelRadius in P[] may
 * be at j's old Ti_current, not at time1. drift_particle's gas-extras block
 * grows h via exp(divv_fac/NDIMS) per drift (capped at exp(divv_fac_max=±0.3
 * /NDIMS) ≈ ±10% per drift), accumulating over multiple deferred drifts.
 * compact_xyzh[j*4+3] inherits the same staleness when populated from P[].
 *
 * The BVH walk reads compact_xyzh[j*4+3] for tile-overlap decisions in
 * symmetric mode. Stale-small h_j → BVH overlap test underestimates →
 * legitimate neighbors missed.
 *
 * Mitigation: at compact_xyzh population time, multiply h by (1 + slack) so
 * the BVH over-includes tiles, absorbing accumulated h-growth from
 * undrifted particles. Per-pair r² acceptance inside kernels reads the
 * REAL P[j].KernelRadius (UVM) — over-inclusion is wasted work, never a
 * silent miss. 0.5 (50% inflation) covers the per-decomp-interval h-growth
 * with margin. */
static constexpr double SIDX_H_SLACK = 0.5;

gpu_spatial_index_t *gpu_step_sidx_ptr(void) { return &g_step_sidx; }
gpu_spatial_index_t *gpu_step_sidx_alltypes_ptr(void) { return &g_step_sidx_alltypes; }


/* Dirty-index tracking for compact_xyzh.h field.
 *
 * Pre-tracker: a single global g_dirty_list/g_dirty_all pair was shared by
 * both g_step_sidx (gas-only) and g_step_sidx_alltypes. Both caches persist
 * across ghost-import-only changes, so one cache consuming and clearing the
 * shared global state would silently leave the other stale -- a physics-
 * correctness hole. Now: per-cache state via gpu_dirty_tracker.
 *
 * Caches register their dense particle-index range [base, base+count) on
 * build, unregister on free. Marks route to ALL caches whose range covers
 * the j (each cache has its own bitset). Refresh consumes only its own
 * cache's bitset.
 *
 * mark_h_dirty_all preserves global semantics: it sets all_dirty on every
 * registered cache (matching the old "unknown-scope mutation"). Per-cache
 * promote-to-all still fires when one cache's popcount exceeds threshold
 * inside the tracker. */

void gpu_compact_xyzh_mark_h_dirty_all(void)
{
    gpu_dirty_tracker_mark_all_global();
}

void gpu_compact_xyzh_mark_h_dirty_idx(int i)
{
    if(i < 0) return;
    int idx_arr[1] = { i };
    gpu_dirty_tracker_mark_indices(idx_arr, 1);
}

void gpu_compact_xyzh_mark_h_dirty_range(int start, int end)
{
    gpu_dirty_tracker_mark_range(start, end);
}

void gpu_compact_xyzh_mark_h_dirty_indices(const int *indices, int n)
{
    gpu_dirty_tracker_mark_indices(indices, n);
}

/* Backwards-compat: any caller that doesn't know which indices it dirtied
 * conservatively forces a full-pool refresh on every cache. */
void gpu_compact_xyzh_mark_h_dirty(void) { gpu_compact_xyzh_mark_h_dirty_all(); }

/* SSOT mark helpers — see header for design.  Route to BOTH the GPU SIDX
 * dirty tracker AND the host glt cache dirty tracker. Adding a further cache
 * later requires only register/unregister inside that cache's lifetime —
 * these helpers automatically include it. */
void gizmo_mark_kernel_radius_dirty_indices(const int *indices, int n)
{
    if(!indices || n <= 0) return;
    gpu_dirty_tracker_mark_indices(indices, n);
    ghost_exchange_local_tree_mark_h_dirty_indices(indices, n);
}
void gizmo_mark_kernel_radius_dirty_range(int start, int end)
{
    if(end <= start) return;
    gpu_dirty_tracker_mark_range(start, end);
    ghost_exchange_local_tree_mark_h_dirty_range(start, end);
}

/* SIDX lifecycle epoch counters, bumped by the notify hooks below on every
 * ghost import, ghost cleanup and pool membership change. gpu_ngb_list_build
 * stamps them into each index at build and requires them to still match before
 * reusing it, so a cached index cannot survive a change of ghost or pool
 * identity that happens to preserve the particle count. */
static uint64_t g_sidx_ghost_epoch = 0;
static uint64_t g_sidx_pool_epoch  = 0;

void gpu_sidx_notify_ghost_imported(int start, int count)
{
    /* Contract requires unconditional call on every rank, including count==0. */
    (void)start; (void)count;
    g_sidx_ghost_epoch++;
}

void gpu_sidx_notify_ghost_cleanup(void)
{
    /* Called BEFORE NumPart shrinks. */
    g_sidx_ghost_epoch++;
}

void gpu_sidx_notify_pool_changed(void)
{
    /* Type/Mass/membership change in the home pool — invalidate any pool-
     * dependent cache on next access. */
    g_sidx_pool_epoch++;
}


/* Drift-time SIDX refresh (incremental rebuild).
 *
 * The unconditional full rebuild post-drift in the prior code cost
 * ~1.3s/step on fire_m11i tiny-N (the dominant tiny-N bucket post-
 * UVM-canonical). Of that, ~1s is build_sfc_tiles streaming 12.4M
 * particles to re-derive pool membership and re-tile it — work that's
 * only needed when particle layout actually changes (i.e.,
 * domain_decomp). There is no sort to skip: P[] arrives Peano-Hilbert
 * ordered from the domain decomposition and tiles are consecutive runs
 * of it, so the cost is the passes over P[], not any ordering step.
 *
 * Between domain decomps, particles drift but their pool/tile
 * assignments are still meaningful: each tile still references the
 * same particles. Their positions just changed slightly. So the
 * refresh path is:
 *
 *   1. For each tile, recompute lo/hi/hmax from the current particle
 *      positions in its pool slice (host OMP — random P[].Pos reads
 *      go to host memory directly under UVM-canonical, no device
 *      page-fault storm).
 *   2. Re-fit the BVH from the updated tile bboxes via build_tile_bvh
 *      (CPU, O(ntiles) ≈ ms for ~70k tiles; structurally identical
 *      to the original BVH because ntiles and tile order are unchanged).
 *   3. Stage updated tiles + BVH to device (deep_copy, ms-scale).
 *   4. Refresh compact_xyzh[i*4+0..2] device-side (existing parallel_for,
 *      ms-scale). NOTE: drift_particle DOES change KernelRadius (predict.cc
 *      lines ~160 and ~229: P[i].KernelRadius *= exp(divv_fac/NUMDIMS)).
 *      The drift-time refresh here updates positions only; the h component
 *      is refreshed separately via mark_h_dirty machinery, which the lazy-
 *      drift path inside gpu_ngb_list_build populates with the j's it
 *      drifted (so compact_xyzh[j*4+3] gets re-read from P_shared on the
 *      next build).
 *
 * Correctness invariant: each tile's bbox covers all current positions
 * of the particles in its pool. Tile assignments are frozen, so as
 * particles "wander" spatially they accumulate into multiple tiles'
 * bbox regions — BVH queries may visit a couple extra tiles per query
 * (inefficiency, never a missed neighbor; each tile still iterates its
 * pool, particles' actual positions are checked).
 *
 * Reset boundary: gpu_step_sidx_invalidate_full() at every domain_decomp
 * frees the SIDX, forcing a fresh rebuild of pool and tiles. Domain
 * decomp is already a heavy step, so the marginal cost is small. This
 * naturally bounds bbox dispersion within a decomp interval and resets
 * tile assignments to current spatial layout.
 *
 */

/* Recomputes tile bboxes from each pool member's "virtual at-time1 position"
 * AND fills the host-side position-staging buffer (idx->h_pos_buf) with the
 * same values.
 *
 * "Virtual at-time1 position" = P[j].Pos + P[j].Vel * get_drift_factor(
 *     P[j].Ti_current, All.Ti_Current, j, 0). This matches what
 * drift_particle's Pos update would produce IF / WHEN the particle is
 * lazily drifted by a downstream consumer. Under the current full-drift
 * regime (move_particles iterates every NumPart particle), Ti_current ==
 * All.Ti_Current for all j, dt = 0, virt_pos == P[j].Pos — i.e. this code
 * is a no-op in absolute value, just exercising the threadsafe drift-factor
 * code path so it's already wired for a future active-only iteration mode
 * in move_particles.
 *
 * Correctness invariant: the bbox covers each particle's actual location
 * at time1, regardless of whether the particle has been drifted yet. BVH
 * queries against the bbox find every potentially-relevant pool member;
 * the consumer then drifts the particle on first read. */
static void sidx_refresh_tile_bboxes_host(gpu_spatial_index_t *idx,
                                             struct particle_data *P_shared)
{
    sfc_tile_t *h_tiles = idx->h_tiles;
    int *h_pool = idx->h_pool;
    int ntiles = idx->ntiles;
    double *pos_buf = idx->h_pos_buf;
    /* SSOT per-particle reach under the cached policy (LEGACY for non-runner
     * callers → byte-equivalent to legacy P[j].KernelRadius aggregation). */
    const mode_b_radius_policy_t policy_capture = idx->cache_radius_policy;
    /* Out-of-line host accessor for the host-side drift-factor input. */
    integertime time1 = gizmo_host_ti_current();

    #pragma omp parallel for schedule(static)
    for(int t = 0; t < ntiles; t++) {
        sfc_tile_t *tile = &h_tiles[t];
        if(tile->count <= 0) continue;
        int j0 = h_pool[tile->first];
        double dt0 = get_drift_factor(P_shared[j0].Ti_current, time1, j0, 0);
        double x0 = P_shared[j0].Pos[0] + P_shared[j0].Vel[0] * dt0;
        double y0 = P_shared[j0].Pos[1] + P_shared[j0].Vel[1] * dt0;
        double z0 = P_shared[j0].Pos[2] + P_shared[j0].Vel[2] * dt0;
        double lo0 = x0, hi0 = x0;
        double lo1 = y0, hi1 = y0;
        double lo2 = z0, hi2 = z0;
        double hmax = nlr_particle_symmetric_radius(P_shared[j0], policy_capture);
        /* Per-type bands recomputed alongside scalar hmax: under the new
         * invariant the bands are policy-aware and would otherwise stay
         * frozen at build-time values, contradicting the conservative-
         * upper-bound rule.  Restart from 0 and aggregate from current
         * particles. */
        double hbt[TILE_NUM_PTYPES] = {0};
        {
            int t0 = (int)P_shared[j0].Type;
            if(t0 >= 0 && t0 < TILE_NUM_PTYPES && hmax > hbt[t0]) hbt[t0] = hmax;
        }
        if(pos_buf) { pos_buf[j0*3+0] = x0; pos_buf[j0*3+1] = y0; pos_buf[j0*3+2] = z0; }
        for(int s = 1; s < tile->count; s++) {
            int j = h_pool[tile->first + s];
            double dt = get_drift_factor(P_shared[j].Ti_current, time1, j, 0);
            double x = P_shared[j].Pos[0] + P_shared[j].Vel[0] * dt;
            double y = P_shared[j].Pos[1] + P_shared[j].Vel[1] * dt;
            double z = P_shared[j].Pos[2] + P_shared[j].Vel[2] * dt;
            if(x < lo0) lo0 = x; else if(x > hi0) hi0 = x;
            if(y < lo1) lo1 = y; else if(y > hi1) hi1 = y;
            if(z < lo2) lo2 = z; else if(z > hi2) hi2 = z;
            double h = nlr_particle_symmetric_radius(P_shared[j], policy_capture);
            if(h > hmax) hmax = h;
            int tj = (int)P_shared[j].Type;
            if(tj >= 0 && tj < TILE_NUM_PTYPES && h > hbt[tj]) hbt[tj] = h;
            if(pos_buf) { pos_buf[j*3+0] = x; pos_buf[j*3+1] = y; pos_buf[j*3+2] = z; }
        }
        tile->lo[0] = lo0; tile->hi[0] = hi0;
        tile->lo[1] = lo1; tile->hi[1] = hi1;
        tile->lo[2] = lo2; tile->hi[2] = hi2;
        tile->hmax = hmax;
        for(int tt = 0; tt < TILE_NUM_PTYPES; tt++) tile->hmax_by_type[tt] = hbt[tt];

    }
}

/* Refresh compact_xyzh positions on device via host-staged buffer.
 *
 * Avoids the UVM-fault-storm cost of a parallel_for that reads P_shared.Pos
 * directly on device (which costs ~1.25s/step on fire_m11i 12.4M after host
 * drift just wrote those pages — every page faults migration to GPU on first
 * access). Instead: positions were already filled into idx->h_pos_buf by the
 * host bbox-recompute loop; bulk deep_copy to d_pos_buf (~200MB at NVLink
 * ~600GB/s = sub-ms), then a small device kernel scatters into the
 * interleaved d_compact_xyzh array. h field intentionally untouched here —
 * NOTE: drift_particle DOES change KernelRadius (predict.cc:160,229), but
 * THIS function is the position-only fast path used at drift-time; the h
 * component is refreshed via the mark_h_dirty machinery on the next
 * gpu_ngb_list_build, which consumes the dirty list populated by every
 * h-writer (lazy drift, density iter, etc.). Splitting pos and h refresh
 * lets us amortize the position update across the whole step while only
 * the touched h slots get refreshed per-build. */
static void sidx_refresh_compact_positions_device(gpu_spatial_index_t *idx)
{
    int num_total = idx->num_total;
    if(!idx->h_pos_buf || !idx->d_pos_buf) return;
    /* Bulk host->device copy of the staged position buffer. Single linear
     * cudaMemcpy under the hood; pages are migrated as one transfer rather
     * than fault-by-fault on first device read. */
    using UV = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
    Kokkos::View<double*, Kokkos::HostSpace, UV>            h_v(idx->h_pos_buf, 3 * num_total);
    Kokkos::View<double*, GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_v(idx->d_pos_buf, 3 * num_total);
    Kokkos::deep_copy(d_v, h_v);

    double *compact = idx->d_compact_xyzh;
    double *pos_buf = idx->d_pos_buf;
    Kokkos::parallel_for("compact_xyzh_pos_scatter", num_total, KOKKOS_LAMBDA(int i) {
        compact[i*4+0] = pos_buf[i*3+0];
        compact[i*4+1] = pos_buf[i*3+1];
        compact[i*4+2] = pos_buf[i*3+2];
    });
    Kokkos::fence();
    gizmo_gpu_check_last_error("compact_xyzh_pos_scatter", num_total);
}

/* Re-fit BVH from updated tile bboxes (structural rebuild — left/right links
 * are recomputed identically since tile order is unchanged, but bboxes/hmax
 * propagate from the refreshed tiles). Calls build_tile_bvh (which mymalloc's
 * a fresh BVH), copies into the persistent HostSpace h_bvh buffer, then frees
 * the mymalloc'd transient. The HostSpace buffer was allocated to size
 * 2*ntiles-1 at build time; ntiles is unchanged across drifts, so the buffer
 * always has room. */
static void sidx_rebuild_bvh_inplace(gpu_spatial_index_t *idx)
{
    tile_bvh_node_t *h_bvh_tmp = NULL;
    int new_nnodes = build_tile_bvh(idx->h_tiles, idx->ntiles, &h_bvh_tmp);
    if(new_nnodes > 0 && h_bvh_tmp && idx->h_bvh) {
        memcpy(idx->h_bvh, h_bvh_tmp, new_nnodes * sizeof(tile_bvh_node_t));
    }
    if(h_bvh_tmp) myfree(h_bvh_tmp);
    idx->h_bvh_nnodes = new_nnodes;
    idx->bvh_root = new_nnodes - 1;
}

/* Stage updated host h_tiles + h_bvh to device d_tiles + d_bvh.
 * h_pool / d_pool unchanged across drifts (tile assignments frozen). */
static void sidx_stage_to_device(gpu_spatial_index_t *idx)
{
    using UV = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
    Kokkos::View<sfc_tile_t*,      Kokkos::HostSpace, UV>            h_tiles_v(idx->h_tiles, idx->ntiles);
    Kokkos::View<sfc_tile_t*,      GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_tiles_v(idx->d_tiles, idx->ntiles);
    Kokkos::View<tile_bvh_node_t*, Kokkos::HostSpace, UV>            h_bvh_v(idx->h_bvh, idx->h_bvh_nnodes);
    Kokkos::View<tile_bvh_node_t*, GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_bvh_v(idx->d_bvh, idx->h_bvh_nnodes);
    Kokkos::deep_copy(d_tiles_v, h_tiles_v);
    Kokkos::deep_copy(d_bvh_v,   h_bvh_v);
}

/* Forward decl */
void gpu_step_sidx_invalidate_full(void);

/* Drift-time refresh: reuse pool/tile membership, recompute bboxes/BVH, refresh compact_xyzh. */
static void sidx_refresh_after_drift(gpu_spatial_index_t *idx,
                                      struct particle_data *P_shared)
{
    sidx_refresh_tile_bboxes_host(idx, P_shared);
    sidx_rebuild_bvh_inplace(idx);
    sidx_stage_to_device(idx);
    sidx_refresh_compact_positions_device(idx);
    idx->positions_stale_after_drift = 0;
}

void gpu_step_sidx_invalidate(void)
{
    /* No All-mirror belt here. Since the position refresh moved to the point of
     * reuse, no path in this function launches a device kernel: the common path
     * only sets a flag, and the rest free buffers. The kernels this file owns
     * are reached through gpu_spatial_index_build and gpu_ngb_list_build, which
     * carry their own belts, and the refresh now runs inside the
     * latter. The leading fence the belt also provided before frees now lives
     * inside gpu_spatial_index_free, with the release it protects. */
    struct particle_data *P_shared = gpu_particles_arena_P();
    if(!P_shared) {
        /* No arena -> no canonical particle storage; fall back to full free. */
        gpu_step_sidx_invalidate_full();
        return;
    }

    /* Current scope: refresh the gas SIDX (the hot-path 1.3s/step bucket
     * from density_sidx_prebuild); alltypes goes through full-free since it
     * only builds on sink-active steps and isn't on the dominant tiny-N
     * path. The refresh path could be extended to alltypes once gas is
     * validated. */
    if(g_step_sidx_alltypes.valid) gpu_spatial_index_free(&g_step_sidx_alltypes);

    gpu_spatial_index_t *idx = &g_step_sidx;
    if(idx->valid) {
        if(!idx->h_tiles || !idx->h_pool || !idx->d_compact_xyzh || idx->ntiles <= 0) {
            /* Defensive: incomplete state -> free, fall back to full rebuild. */
            gpu_spatial_index_free(idx);
            gpu_compact_xyzh_mark_h_dirty_all();
        } else {
            /* Mark the positions stale rather than refreshing them here. The
             * refresh happens at the point of reuse, in gpu_ngb_list_build, so
             * a consumer that invalidates the index on count or epoch — and so
             * rebuilds from current positions regardless — never pays for a
             * refresh whose result it discards. Correctness is unchanged: the
             * only paths that walk this index go through that build, which
             * refreshes first and hard-aborts if it ever sees a stale one.
             * The h component was already handled exactly this way. */
            idx->positions_stale_after_drift = 1;
            /* h-dirty state intentionally left intact. drift_particle DOES
             * change KernelRadius (predict.cc:160,229) — those h updates are
             * marked into the per-cache dirty tracker (gpu_dirty_tracker) by
             * the lazy-drift loop in the previous step's gpu_ngb_list_build,
             * by move_particles/gizmo_full_drift_to (predict.cc:307,351), and by other
             * h-writers like density iter. The next gpu_ngb_list_build
             * consume()s this cache's bits and refreshes compact_xyzh[*4+3]
             * from current P_shared.KernelRadius before walking. */
        }
    }

}

void gpu_step_sidx_invalidate_full(void)
{
    if(g_step_sidx_alltypes.valid) gpu_spatial_index_free(&g_step_sidx_alltypes);
    if(g_step_sidx.valid) gpu_spatial_index_free(&g_step_sidx);
    /* Force full-mode dirty so next build's compact_xyzh seeds correctly. */
    gpu_compact_xyzh_mark_h_dirty_all();
}


/* An index that could not be built must not be left looking usable: the walk reads
 * its tiles, BVH and compact positions directly. Release what was built and leave
 * the index marked invalid, which is the signal every consumer already tests, so
 * the caller can hand back an empty neighbour list. gpu_spatial_index_free fences
 * before releasing device memory, so nothing is released under a running kernel.
 * The three host build buffers are arena allocations, released in reverse order
 * exactly as the success path does. */
static void sidx_build_leave_invalid(gpu_spatial_index_t *idx, int num_total,
                                     tile_bvh_node_t *h_bvh, sfc_tile_t *h_tiles, int *h_pool,
                                     const char *what, size_t bytes)
{
    myfree(h_bvh);
    myfree(h_tiles);
    myfree(h_pool);
    gpu_spatial_index_free(idx);
    char msg[256];
    snprintf(msg, sizeof(msg),
             "gpu_spatial_index_build: could not allocate %s (%.1f MB) for %d particles; "
             "spatial index left unbuilt",
             what, (double) bytes / (1024.0 * 1024.0), num_total);
    gizmo_request_controlled_stop(7712, msg, __FILE__, __LINE__, __FUNCTION__);
}

void gpu_spatial_index_build(struct particle_data *P_shared, int num_total,
                             int type_bitmask, gpu_spatial_index_t *idx,
                             const char *caller_label,
                             mode_b_radius_policy_t radius_policy)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

#if defined(BOX_PERIODIC)
    /* Hard guard on the box lengths the wrap depends on.  The device predicates
     * image through the canonical macros, which read the boxSize_ and boxHalf_
     * lengths from
     * this TU's AllDeviceMirror; if that mirror is unsynced they read zero,
     * every wrapped separation collapses, the per-axis prunes go dead and the
     * BVH silently degenerates into an exhaustive scan.  Fail loud at the call
     * site so a regression in the mirror sync can never quietly cost the index.
     * (Previously this checked the per-axis flags/box_sizes copies; those are
     * gone with that API, but the failure mode belongs to the globals and is
     * unchanged, so the guard now reads them directly.) */
    {
        const double box_len[3] = { boxSize_X, boxSize_Y, boxSize_Z };
        const int    wraps[3]   = { TILE_PERIODIC_X, TILE_PERIODIC_Y, TILE_PERIODIC_Z };
        for(int k = 0; k < 3; k++) {
            if(wraps[k] && !(box_len[k] > 0.0)) {
                printf("gpu_spatial_index_build: axis %d is periodic but its box length is %g "
                       "(caller='%s'). Likely cause: this TU's AllDeviceMirror not synced from "
                       "host All. Confirm gizmo_gpu_sync_all() has run for this timestep.\n",
                       k, box_len[k], caller_label ? caller_label : "?");
                fflush(stdout);
                endrun(913004);
            }
        }
    }
#endif

    /* Build SFC tiles + BVH on CPU */
    sfc_tile_t *h_tiles;
    int *h_pool;
    int num_pool;
    int ntiles = build_sfc_tiles(P_shared, num_total, type_bitmask, TILE_TARGET_SIZE,
                                 &h_tiles, &h_pool, &num_pool, radius_policy);
    idx->ntiles = ntiles;

    tile_bvh_node_t *h_bvh;
    int bvh_nnodes = build_tile_bvh(h_tiles, ntiles, &h_bvh);
    idx->bvh_root = bvh_nnodes - 1;

    /* Allocate kernel-read-path arrays in DEVICE_SPACE (CudaSpace HBM on GPU
     * builds, falls back to SharedSpace elsewhere).  This eliminates HMM/TLB-
     * miss overhead on the small-N kernel hot path where one thread does
     * ~1000s of scattered reads through bvh/tiles/pool/compact_xyzh — that
     * scattered-UVM-access pattern is the suspected source of the residual
     * 1.4s "fused_fnc" floor on 1-active-particle calls.  CPU host arrays
     * h_tiles/h_bvh/h_pool are transferred via Kokkos::deep_copy through
     * unmanaged-View wrappers (cudaMemcpy under the hood on CUDA builds);
     * compact_xyzh is built directly on-device by a parallel_for that reads
     * from UVM-backed P_shared and writes to DEVICE_SPACE compact_xyzh. */
    int bvh_size = (2 * ntiles - 1);
    if(bvh_size < 1) bvh_size = 1;
    int pool_size = (num_pool > 0) ? num_pool : 1;
    /* Sized to at least one element, as bvh_size and pool_size already are: a type
     * mask that matches nothing would otherwise ask for zero bytes, and a request
     * for none of something cannot be distinguished from a refusal. */
    size_t sidx_tiles_bytes = (size_t)((ntiles > 0) ? ntiles : 1) * sizeof(sfc_tile_t);
    size_t sidx_bvh_bytes   = (size_t) bvh_size * sizeof(tile_bvh_node_t);
    size_t sidx_pool_bytes  = (size_t) pool_size * sizeof(int);
    idx->d_tiles = (sfc_tile_t *) gizmo_gpu_alloc_device(sidx_tiles_bytes, "ngl_sidx_dev_tiles");
    if(!idx->d_tiles) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the device tile array", sidx_tiles_bytes); return;}
    idx->d_bvh = (tile_bvh_node_t *) gizmo_gpu_alloc_device(sidx_bvh_bytes, "ngl_sidx_dev_bvh");
    if(!idx->d_bvh) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the device tile BVH", sidx_bvh_bytes); return;}
    idx->d_pool = (int *) gizmo_gpu_alloc_device(sidx_pool_bytes, "ngl_sidx_dev_pool");
    if(!idx->d_pool) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the device tile membership pool", sidx_pool_bytes); return;}

    /* Stage host buffers into device memory.  On non-CUDA builds DEVICE_SPACE
     * == SharedSpace and Kokkos::deep_copy reduces to a memcpy. */
    {
        using UV = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
        Kokkos::View<sfc_tile_t*,        Kokkos::HostSpace, UV>            h_tiles_v(h_tiles, ntiles);
        Kokkos::View<sfc_tile_t*,        GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_tiles_v(idx->d_tiles, ntiles);
        Kokkos::View<tile_bvh_node_t*,   Kokkos::HostSpace, UV>            h_bvh_v(h_bvh, bvh_nnodes);
        Kokkos::View<tile_bvh_node_t*,   GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_bvh_v(idx->d_bvh, bvh_nnodes);
        Kokkos::View<int*,               Kokkos::HostSpace, UV>            h_pool_v(h_pool, num_pool);
        Kokkos::View<int*,               GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_pool_v(idx->d_pool, num_pool);
        Kokkos::deep_copy(d_tiles_v, h_tiles_v);
        Kokkos::deep_copy(d_bvh_v,   h_bvh_v);
        Kokkos::deep_copy(d_pool_v,  h_pool_v);
    }

    /* Build compact double4 position+h array for cache-efficient GPU BVH traversal.
       DOUBLE positions: float ABSOLUTE positions are invalid for GIZMO's ~1e11
       dynamic range (see §37/§38) — they must NOT decide neighbour inclusion.
       ~64MB for 2M particles vs 800MB for full P_shared. h (slot 3) is a relative
       reach; kept double so the leaf accept stays consistent with the double opener.
       In DEVICE_SPACE so the BVH-walk kernels read from HBM directly. */
    size_t sidx_compact_bytes = (size_t)((num_total > 0) ? num_total : 1) * 4 * sizeof(double);
    idx->d_compact_xyzh = (double *) gizmo_gpu_alloc_device(sidx_compact_bytes, "ngl_sidx_dev_compact_xyzh");
    if(!idx->d_compact_xyzh) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the device compact position array", sidx_compact_bytes); return;}
    {
        double *compact = idx->d_compact_xyzh;
        double h_inflate = 1.0 + SIDX_H_SLACK; /* lazy-drift slack: see SIDX_H_SLACK comment */
        const mode_b_radius_policy_t policy_capture = radius_policy;
        Kokkos::parallel_for("compact_xyzh_build", num_total, KOKKOS_LAMBDA(int i) {
            compact[i*4+0] = P_shared[i].Pos[0];
            compact[i*4+1] = P_shared[i].Pos[1];
            compact[i*4+2] = P_shared[i].Pos[2];
            /* SSOT per-particle reach under policy_capture — see nlr_radius_policy.h.
             * LEGACY default policy recovers raw P[j].KernelRadius for every type. */
            double h_j = nlr_particle_symmetric_radius(P_shared[i], policy_capture);
            compact[i*4+3] = h_j * h_inflate;
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("compact_xyzh_build", num_total);
    }

    /* Keep host-side persistent copies of tiles/pool/BVH alive across drifts
     * so the incremental refresh path (gpu_step_sidx_invalidate ->
     * sidx_refresh_after_drift) can recompute tile bboxes from current
     * particle positions on host without re-deriving pool membership or
     * re-tiling it.
     *
     * MUST use Kokkos::HostSpace allocator (heap-based) NOT mymalloc — the
     * latter is LIFO-stack-disciplined and persistent SIDX buffers would
     * sit on top of any subsequent transient allocation (e.g. density's
     * per-step mymalloc'd Left/Right arrays), preventing those transients
     * from being freed in LIFO order. The transient mymalloc'd h_tiles /
     * h_pool / h_bvh from build_sfc_tiles + build_tile_bvh are copied
     * out then myfree'd in proper LIFO order below. */
    idx->h_tiles = (sfc_tile_t *) gizmo_gpu_alloc_host(sidx_tiles_bytes, "ngl_sidx_host_tiles");
    if(!idx->h_tiles) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the retained host tile array", sidx_tiles_bytes); return;}
    memcpy(idx->h_tiles, h_tiles, ntiles * sizeof(sfc_tile_t));
    idx->h_pool = (int *) gizmo_gpu_alloc_host(sidx_pool_bytes, "ngl_sidx_host_pool");
    if(!idx->h_pool) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the retained host tile membership pool", sidx_pool_bytes); return;}
    memcpy(idx->h_pool, h_pool, num_pool * sizeof(int));
    idx->h_bvh = (tile_bvh_node_t *) gizmo_gpu_alloc_host(sidx_bvh_bytes, "ngl_sidx_host_bvh");
    if(!idx->h_bvh) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the retained host tile BVH", sidx_bvh_bytes); return;}
    memcpy(idx->h_bvh, h_bvh, bvh_nnodes * sizeof(tile_bvh_node_t));
    idx->h_bvh_nnodes = bvh_nnodes;
    idx->num_pool = num_pool;

    /* Drift-refresh staging buffers: h_pos_buf is filled per-pool-member by the
     * host bbox-recompute loop; bulk deep_copy to d_pos_buf; device scatter
     * into d_compact_xyzh. Sized 3*num_total doubles so each pool index can
     * write directly to h_pos_buf[j*3+0..2] without remapping. Non-pool
     * entries stay uninitialized in h_pos_buf and their stale d_compact_xyzh
     * positions are never read by BVH queries (BVH only visits tiles, tiles
     * only contain pool members). */
    int pos_buf_count = (num_total > 0 ? num_total : 1);
    size_t sidx_pos_buf_bytes = (size_t) 3 * pos_buf_count * sizeof(double);
    idx->h_pos_buf = (double *) gizmo_gpu_alloc_host(sidx_pos_buf_bytes, "ngl_sidx_host_pos_buf");
    if(!idx->h_pos_buf) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the host position staging buffer", sidx_pos_buf_bytes); return;}
    idx->d_pos_buf = (double *) gizmo_gpu_alloc_device(sidx_pos_buf_bytes, "ngl_sidx_dev_pos_buf");
    if(!idx->d_pos_buf) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_tiles, h_pool, "the device position staging buffer", sidx_pos_buf_bytes); return;}


    /* Free the transient mymalloc'd build buffers in proper LIFO order
     * (build_tile_bvh allocated h_bvh last; build_sfc_tiles allocated
     * h_pool then h_tiles). */
    myfree(h_bvh);
    myfree(h_tiles);
    myfree(h_pool);

    idx->num_total = num_total;
    idx->cache_tbm = type_bitmask;
    idx->cache_radius_policy = radius_policy;
    idx->ghost_epoch_when_built = g_sidx_ghost_epoch;
    idx->pool_epoch_when_built  = g_sidx_pool_epoch;
    idx->valid = 1;
    /* Built from current positions, so no refresh is outstanding for this index.
     * Explicit because the struct may be a cached one being rebuilt in place. */
    idx->positions_stale_after_drift = 0;
    /* Register this cache with the dirty tracker over [0, num_total). The
     * compact_xyzh build above wrote every row's h from the live P[] under this
     * cache's radius policy, and nothing between there and here can mutate it,
     * so the range starts clean: the first refresh would recompute values it
     * already holds, over the whole pool. */
    if(idx->dirty_handle >= 0) gpu_dirty_tracker_unregister(idx->dirty_handle);
    idx->dirty_handle = gpu_dirty_tracker_register(0, num_total, 1);

}

void gpu_spatial_index_free(gpu_spatial_index_t *idx)
{
    /* Ordering, not cleanup. kokkos_free does not synchronize, so releasing a
     * device allocation while a kernel may still be reading it is a
     * use-after-free. The fence lives HERE, with the release, so that no caller
     * can omit it -- this function is reached from the step loop, the
     * decomposition boundary and the cached-index staleness guard, and putting
     * the rule in any one of those leaves the next caller free to reintroduce
     * the hazard. Skipped when there is nothing device-side to release. */
    if(idx->d_compact_xyzh || idx->d_pool || idx->d_bvh || idx->d_tiles || idx->d_pos_buf) {
        Kokkos::fence();
    }
    if(idx->d_compact_xyzh) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_compact_xyzh); idx->d_compact_xyzh = NULL;}
    if(idx->d_pool) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_pool); idx->d_pool = NULL;}
    if(idx->d_bvh) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_bvh); idx->d_bvh = NULL;}
    if(idx->d_tiles) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_tiles); idx->d_tiles = NULL;}
    /* Free host-side persistent buffers kept alive across drifts.
     * mymalloc uses LIFO stack discipline; free in reverse-allocation order:
     * h_pool -> h_tiles -> h_bvh
     * (h_bvh allocated first by build_tile_bvh, but build_tile_bvh may have
     * been re-called via sidx_rebuild_bvh_inplace which freed-then-allocated,
     * so h_bvh is on top of the stack at this point in normal flow). */
    /* Persistent host-side buffers live in Kokkos::HostSpace (heap-allocated,
     * not mymalloc) so they don't pin the LIFO stack across other transient
     * mymalloc'd state. Free order doesn't matter. */
    if(idx->h_bvh)   { Kokkos::kokkos_free<Kokkos::HostSpace>(idx->h_bvh);   idx->h_bvh = NULL; }
    if(idx->h_tiles) { Kokkos::kokkos_free<Kokkos::HostSpace>(idx->h_tiles); idx->h_tiles = NULL; }
    if(idx->h_pool)  { Kokkos::kokkos_free<Kokkos::HostSpace>(idx->h_pool);  idx->h_pool = NULL; }
    if(idx->h_pos_buf) { Kokkos::kokkos_free<Kokkos::HostSpace>(idx->h_pos_buf); idx->h_pos_buf = NULL; }
    if(idx->d_pos_buf) { Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_pos_buf); idx->d_pos_buf = NULL; }
    idx->h_bvh_nnodes = 0;
    idx->num_pool = 0;
    idx->num_total = 0;
    idx->cache_tbm = -1;
    idx->cache_radius_policy = MODE_B_RADIUS_DEFAULT;
    idx->valid = 0;
    idx->positions_stale_after_drift = 0;  /* nothing left to refresh */
    if(idx->dirty_handle >= 0) {
        gpu_dirty_tracker_unregister(idx->dirty_handle);
        idx->dirty_handle = -1;
    }
}


/* Exhaustion is reported by returning NULL, so the checks below are live code and
 * the run stops cleanly instead of aborting mid-flight. A zero-byte request is
 * treated as nothing to allocate, which the callers here read as a refusal. */
static void *ngl_alloc_shared(size_t bytes, const char *label)
{
    if(bytes == 0) {return NULL;}
    return gizmo_gpu_alloc_shared(bytes, label);
}
static void *ngl_alloc_device(size_t bytes, const char *label)
{
    if(bytes == 0) {return NULL;}
    return gizmo_gpu_alloc_device(bytes, label);
}

/* The per-active and per-pair arrays are the largest transients this loop asks for --
 * the walk scratchpad alone is 512 int slots per active particle. When one cannot be
 * had, the node is out of memory at the size this loop needs. Release the transients,
 * leave a VALID EMPTY list (every active with zero neighbours) so that whatever runs
 * between here and the stop draining walks nothing rather than reading a half-built
 * list, and name the buffer that could not be had. */
static void ngl_build_leave_empty(gpu_neighbor_list_t *gnl, int num_active,
                                  double *d_radii, double *d_source_pos,
                                  int *d_scratch, int *d_counts,
                                  const char *what, size_t bytes)
{
    if(d_counts)     {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_counts);}
    if(d_scratch)    {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_scratch);}
    if(d_source_pos) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_source_pos);}
    if(d_radii)      {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);}
    /* The row offsets are what makes the empty list readable: every consumer walks
     * offsets[aa]..offsets[aa+1] unconditionally, and only reaches the neighbour
     * array when total_pairs > 0. So when the build fails before the offsets exist,
     * this allocates them here rather than handing back a null row index. The
     * request is (num_active+1) 8-byte slots -- orders of magnitude below the walk
     * scratchpad and the pair list whose failure brings us here -- so it is served
     * from what those releases just returned. Should even that fail, the consumer
     * reads a null row index and dies where it would have died anyway; there is no
     * smaller allocation left to fall back to. */
    if(!gnl->offsets) {
        gnl->offsets = (int64_t *) ngl_alloc_shared((size_t)(num_active + 1) * sizeof(int64_t),
                                                    "ngl_pairs_offsets");
    }
    if(gnl->offsets) {for(int aa = 0; aa <= num_active; aa++) {gnl->offsets[aa] = 0;}}
    gnl->total_pairs = 0;
    char msg[256];
    snprintf(msg, sizeof(msg),
             "gpu_ngb_list_build: could not allocate %s (%.1f MB) for %d active particles; "
             "neighbour list left empty",
             what, (double) bytes / (1024.0 * 1024.0), num_active);
    gizmo_request_controlled_stop(7711, msg, __FILE__, __LINE__, __FUNCTION__);
}

void gpu_ngb_list_build(struct particle_data *P_shared, int num_total,
                        int *active_indices_host, int num_active,
                        int search_mode, int type_bitmask,
                        gpu_neighbor_list_t *gnl,
                        gpu_spatial_index_t *cached_idx,
                        double search_radius_factor,
                        const double *search_radii_host,
                        const double *source_positions_host,
                        const char *caller_label,
                        double j_kernel_radius_scale,
                        mode_b_radius_policy_t radius_policy)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    gnl->num_active = num_active;
    /* Everything this function hands back starts null, so a failure part-way through
     * leaves a list whose unbuilt parts are recognisable rather than whatever the
     * caller's struct happened to contain. Several callers declare it uninitialised,
     * and the free path releases the spatial-index mirrors too, so those must be
     * cleared here as well and not only where they are copied from a built index. */
    gnl->d_active = NULL; gnl->offsets = NULL; gnl->neighbors = NULL; gnl->total_pairs = 0;
    gnl->d_tiles = NULL; gnl->d_bvh = NULL; gnl->d_pool = NULL; gnl->d_compact_xyzh = NULL;
    gnl->ntiles = 0; gnl->bvh_root = 0;
    double t_entry = my_second(); /* DIAG: entry */
    const double cpu_rows_child0 = CPU_ChildCharged;
    /* Defensive guard: a cached SIDX's compact_xyzh / pool only contains the
     * originally-built types. If the caller's type_bitmask differs from the
     * cache's, the walker would return neighbors of types outside the
     * caller's mask (e.g., DM neighbors leaking into a gas-only density
     * walk → lazy-drift attempts on Type=1 → drift_particle abort
     * 'no prediction into past allowed'). HARD-ABORT with a clear message
     * so any future Spec-author who routes a tbm-mismatched cache fails at
     * the right layer, not via downstream nonsense. The companion warning
     * in gpu_neighbor_list.h on gpu_step_sidx_alltypes_ptr documents this
     * invariant. */
    if (cached_idx && cached_idx->valid && cached_idx->cache_tbm >= 0 &&
        cached_idx->cache_tbm != type_bitmask) {
        fprintf(stderr,
            "gpu_ngb_list_build FATAL: caller='%s' type_bitmask=0x%x but cached "
            "SIDX was built with tbm=0x%x. The cached compact_xyzh / pool only "
            "contains the originally-built types — walking it under a different "
            "mask returns wrong-type neighbors and triggers downstream aborts "
            "(e.g. drift_particle 'no prediction into past allowed' on Type=1 "
            "during a gas-only density walk). Spec author: declare "
            "sidx_cache_kind matching your neighbor_type_mask, or rebuild "
            "the cache. See gpu_neighbor_list.h docstring on "
            "gpu_step_sidx_alltypes_ptr.\n",
            caller_label ? caller_label : "?", type_bitmask, cached_idx->cache_tbm);
        fflush(stderr);
        endrun(913005);
    }
    /* Companion HARD-ABORT for radius_policy mismatch.  The
     * cached compact_xyzh[j*4+3] reflects the build-time policy's per-j reach;
     * walking it under a different Spec policy would silently use the wrong
     * leaf-side h_j and miss valid pairs.  Same shape as the cache_tbm gate. */
    if (cached_idx && cached_idx->valid &&
        cached_idx->cache_radius_policy != radius_policy) {
        fprintf(stderr,
            "gpu_ngb_list_build FATAL: caller='%s' radius_policy=0x%x but cached "
            "SIDX was built with radius_policy=0x%x. The cached compact_xyzh[j*4+3] "
            "encodes the build-time policy's per-particle pair-search reach; "
            "walking it under a different policy returns wrong leaf-h_j values "
            "and silently misses valid pairs.  Spec author: ensure Spec::radius_policy "
            "matches the cache's owner (or declare SidxCacheKind::None to rebuild "
            "per call).\n",
            caller_label ? caller_label : "?",
            (unsigned)radius_policy, (unsigned)cached_idx->cache_radius_policy);
        fflush(stderr);
        endrun(913006);
    }
    /* A caller that reuses a cached index without supplying source geometry
     * takes each active's position and radius from compact_xyzh[active_index],
     * which holds POOL MEMBERS ONLY.  INVARIANT: that is correct exactly while
     * every active is itself a pool member — on a reused cache the incremental
     * refresh updates pool members, so a non-pool active would read an
     * unrefreshed position and get a wrong neighbour set, silently.  A caller
     * that cannot guarantee it passes explicit source positions instead.  The
     * runner's compile-time check covers the Specs it dispatches; the symmetric
     * list path below satisfies the invariant by being gas-only into a gas-only
     * pool.  The call form is what should make the unsafe combination
     * unrepresentable, rather than a scan looking for it after the fact. */

    /* Early-out: with no active particles there is nothing to search.
     * Skip the SIDX build/refresh AND all kernel launches.  Allocate 1-element
     * stubs so the caller's gpu_ngb_list_free path is well-defined (it always
     * frees neighbors/offsets/d_active).  Use cached SIDX pointers if available
     * so the free path's "is this from cache?" comparison still works. */
    if(num_active == 0) {
        gpu_spatial_index_t *idx_for_stubs = NULL;
        if(cached_idx && cached_idx->valid && cached_idx->num_total == num_total) {
            idx_for_stubs = cached_idx;
        }
        if(idx_for_stubs) {
            gnl->d_tiles  = idx_for_stubs->d_tiles;
            gnl->d_bvh    = idx_for_stubs->d_bvh;
            gnl->d_pool   = idx_for_stubs->d_pool;
            gnl->d_compact_xyzh = idx_for_stubs->d_compact_xyzh;
            gnl->ntiles   = idx_for_stubs->ntiles;
            gnl->bvh_root = idx_for_stubs->bvh_root;
        } else {
            gnl->d_tiles = NULL; gnl->d_bvh = NULL; gnl->d_pool = NULL; gnl->d_compact_xyzh = NULL;
            gnl->ntiles = 0; gnl->bvh_root = 0;
        }
        gnl->d_active  = (int *) ngl_alloc_shared(sizeof(int), "ngl_pairs_active_stub");
        gnl->offsets   = (int64_t *) ngl_alloc_shared(sizeof(int64_t), "ngl_pairs_offsets_stub");
        gnl->neighbors = (int *) ngl_alloc_device(sizeof(int), "ngl_pairs_neighbors_stub");
        /* Single elements: failing these means the node has no memory left at all.
         * The empty list is the same one every other exhausted allocation here
         * hands back, so the callers need no separate case. */
        if(!gnl->d_active || !gnl->offsets || !gnl->neighbors) {
            ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL,
                                  "the empty-list placeholders", sizeof(int64_t));
            cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
            return;
        }
        gnl->offsets[0] = 0;
        gnl->total_pairs = 0;
        cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
        return;
    }

    /* Use cached spatial index if available, otherwise build fresh.
     * If caller provided a cached_idx but it's not yet built, populate it (this
     * enables persistent caching across calls — caller controls invalidation). */
    /* Value-initialise: the positional form listed fewer entries than the struct has members and
     * put them in the wrong slots, so it overrode the members that carry a non-zero default --
     * dirty_handle (-1 = not registered), cache_tbm (-1 = no build yet) and cache_radius_policy.
     * A zero dirty_handle passes the `>= 0` test that guards gpu_dirty_tracker_unregister, so an
     * index that had never registered anything would unregister handle 0, which belongs to
     * whoever did register first. */
    gpu_spatial_index_t local_idx = {};
    gpu_spatial_index_t *idx;
    /* Invalidate the cached SIDX unless it still describes the same particles.
     * num_total: the compact_xyzh and pool arrays were sized for the old count,
     * so accessing beyond them is UB (ghost exchange redo, particle creation).
     * Epochs: a cleanup-and-reimport can land the SAME ghost count with
     * different ghost contents, which no count test can see. The index would
     * then hold stale positions, tile bounds and BVH, because a ghost import
     * marks h-dirty only and the refresh kernels rewrite compact_xyzh[i*4+3]
     * alone. Drift does not bump either epoch -- membership is unchanged there,
     * and the drift refresh path handles moved positions -- so this costs no
     * rebuild on the common path. */
    if(cached_idx && cached_idx->valid &&
       (cached_idx->num_total          != num_total          ||
        cached_idx->ghost_epoch_when_built != g_sidx_ghost_epoch ||
        cached_idx->pool_epoch_when_built  != g_sidx_pool_epoch)) {
        gpu_spatial_index_free(cached_idx);
    }
    if(cached_idx && cached_idx->valid) {
        /* Reuse: settle the outstanding position refresh BEFORE the index becomes
         * consumer-visible below. This is the only path on which a drift-time
         * refresh is actually needed, which is why it waits until here. */
        if(cached_idx->positions_stale_after_drift) {
            /* Charged to the refresh bucket even though it runs here, so the
             * bucket keeps naming the work it holds rather than the place the
             * work happens; the enclosing list-build charge deducts it as a
             * child. Without this a refresh running here would silently inflate the
             * list-build row instead. */
            const double t_refresh_start = my_second();
            const double child0_refresh = CPU_ChildCharged;
            sidx_refresh_after_drift(cached_idx, P_shared);
            cpu_charge_child(CPU_SIDX_REFRESH,
                             cpu_minus_children(timediff(t_refresh_start, my_second()), child0_refresh));
        }
        idx = cached_idx;
    } else if(cached_idx) {
        gpu_spatial_index_build(P_shared, num_total, type_bitmask, cached_idx, caller_label, radius_policy);
        idx = cached_idx;
    } else {
        gpu_spatial_index_build(P_shared, num_total, type_bitmask, &local_idx, caller_label, radius_policy);
        idx = &local_idx;
    }
    /* A build that ran out of memory leaves the index invalid and has already asked
     * for the stop, naming the buffer. There is nothing to walk, so hand back the
     * same empty list any other exhausted allocation here produces. */
    if(!idx->valid) {
        ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL, "the spatial index", 0);
        return;
    }

    /* Coordinate-staleness invariant, same class as the cache_tbm and
     * radius_policy guards above: a walk over an index whose tile bboxes, BVH
     * and compact positions predate the last drift silently misses neighbours.
     * Every path reaching here has either rebuilt the index from current
     * positions or refreshed it, so this can only fire if a future caller
     * introduces a third path. Fail loudly at the right layer. */
    if(idx->positions_stale_after_drift) {
        fprintf(stderr,
            "gpu_ngb_list_build FATAL: caller='%s' is about to walk a spatial index "
            "whose positions predate the last drift. The tile bounding boxes, BVH and "
            "compact positions describe where particles WERE, so the walk would miss "
            "genuine neighbours without any error. Every consumer must either rebuild "
            "the index or refresh it before use; a path that does neither has been "
            "added.\n", caller_label ? caller_label : "?");
        fflush(stderr);
        endrun(913007);
    }

    /* Copy spatial index pointers to gnl for use by free */
    gnl->d_tiles = idx->d_tiles;
    gnl->d_bvh = idx->d_bvh;
    gnl->d_pool = idx->d_pool;
    gnl->d_compact_xyzh = idx->d_compact_xyzh;
    gnl->ntiles = idx->ntiles;
    gnl->bvh_root = idx->bvh_root;

    /* Refresh the h component of the compact array. Two modes (driven by the
     * per-cache gpu_dirty_tracker):
     *  - all-dirty: full-pool parallel_for(num_total, ...) — pays ~1.1-1.2s
     *    per call on fire_m11i 12.4M pool (UVM fault latency on P_shared
     *    KernelRadius reads).
     *  - bitset drain: parallel_for(n_dirty, ...) reading indices staged from
     *    this cache's bitset to a device buffer. ~ms per call when n_dirty is
     *    the active set. Per-cache state means consuming-and-clearing this
     *    cache's bits leaves the other registered caches' bitsets untouched.
     * Skipped entirely when this cache has neither all_dirty nor any set bits. */
    /* Per-cache dirty tracker query: this cache's bitset is independent of
     * other caches' state. consume() iterates set bits, populates d_dirty,
     * then clears bitset+all_dirty for THIS cache only. */
    int do_refresh = 0, refresh_all = 0;
    int handle = cached_idx ? cached_idx->dirty_handle : -1;
    if(cached_idx && cached_idx->valid && handle >= 0) {
        if(gpu_dirty_tracker_is_all_dirty(handle)) { do_refresh = 1; refresh_all = 1; }
        else if(gpu_dirty_tracker_popcount(handle) > 0) { do_refresh = 1; refresh_all = 0; }
    }
    if(do_refresh) {
        double *compact = idx->d_compact_xyzh;
        double h_inflate = 1.0 + SIDX_H_SLACK; /* see SIDX_H_SLACK — lazy-drift over-search slack */
        /* SSOT per-j reach under the cached policy.  HARD-ABORT above guarantees
         * cached_idx->cache_radius_policy == caller's radius_policy. */
        const mode_b_radius_policy_t policy_capture = idx->cache_radius_policy;
        if(refresh_all) {
            Kokkos::parallel_for("compact_h_refresh_all", num_total, KOKKOS_LAMBDA(int i) {
                double h_j = nlr_particle_symmetric_radius(P_shared[i], policy_capture);
                compact[i*4+3] = h_j * h_inflate;
            });
            /* Drain the bitset (no kernel use; just clear it). */
            struct {} dummy;
            gpu_dirty_tracker_consume(handle,
                [](int j, void *ud){ (void)j; (void)ud; },
                &dummy);
        } else {
            /* Drain bitset → host vector → device buffer → kernel. */
            std::vector<int> dirty_host;
            dirty_host.reserve(gpu_dirty_tracker_popcount(handle));
            gpu_dirty_tracker_consume(handle,
                [](int j, void *ud){ ((std::vector<int> *)ud)->push_back(j); },
                &dirty_host);
            int n_dirty = (int)dirty_host.size();
            size_t dirty_bytes = (size_t)((n_dirty > 0) ? n_dirty : 1) * sizeof(int);
            int *d_dirty = (int *) ngl_alloc_device(dirty_bytes, "ngl_dirty");
            /* The refresh cannot be skipped: the bits have already been consumed, so
             * leaving them unwritten would let a later walk read stale reaches. Stop
             * instead, and leave the list empty so nothing walks the index meanwhile. */
            if(!d_dirty) {
                ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL,
                                      "the refreshed-particle index list", dirty_bytes);
                cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
                return;
            }
            {
                Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                    hv(dirty_host.data(), n_dirty);
                Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                    dv(d_dirty, n_dirty);
                Kokkos::deep_copy(dv, hv);
            }
            Kokkos::parallel_for("compact_h_refresh_idx", n_dirty, KOKKOS_LAMBDA(int k) {
                int i = d_dirty[k];
                double h_j = nlr_particle_symmetric_radius(P_shared[i], policy_capture);
                compact[i*4+3] = h_j * h_inflate;
            });
            Kokkos::fence();
            gizmo_gpu_check_last_error("compact_h_refresh_idx", n_dirty);
            Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_dirty);
        }
        Kokkos::fence();
        gizmo_gpu_check_last_error("compact_h_refresh", num_total);
    }

    /* Active indices: always re-uploaded (changes per call) */
    size_t active_bytes = (size_t)((num_active > 0) ? num_active : 1) * sizeof(int);
    gnl->d_active = (int *) ngl_alloc_shared(active_bytes, "ngl_pairs_active");
    if(!gnl->d_active) {ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL, "the active-index list", active_bytes); return;}
    memcpy(gnl->d_active, active_indices_host, num_active * sizeof(int));

    /* Optional explicit per-active search radii (for loops with a different kernel
       than P[i].KernelRadius, e.g. KernelRadiusDM or AGS_Hsml). NULL → use P[i].KernelRadius. */
    double *d_radii = NULL;
    if(search_radii_host) {
        size_t radii_bytes = (size_t)((num_active > 0) ? num_active : 1) * sizeof(double);
        d_radii = (double *) ngl_alloc_shared(radii_bytes, "ngl_pairs_radii");
        if(!d_radii) {ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL, "the per-active search radii", radii_bytes); return;}
        memcpy(d_radii, search_radii_host, num_active * sizeof(double));
    }

    /* Optional explicit per-active source positions (for sources not backed by
       P[] entries, e.g. arbitrary grid cells). NULL → read pos from P[active[aa]].
       Layout in caller's array: source_positions_host[aa*3 + k] for axis k. */
    double *d_source_pos = NULL;
    if(source_positions_host) {
        size_t srcpos_bytes = (size_t)((num_active > 0) ? num_active : 1) * 3 * sizeof(double);
        d_source_pos = (double *) ngl_alloc_shared(srcpos_bytes, "ngl_pairs_source_pos");
        if(!d_source_pos) {ngl_build_leave_empty(gnl, num_active, d_radii, NULL, NULL, NULL, "the per-active source positions", srcpos_bytes); return;}
        memcpy(d_source_pos, source_positions_host, num_active * 3 * sizeof(double));
    }

    /* Allocate CSR offsets (64-bit row pointers) */
    size_t offsets_bytes = (size_t)(num_active + 1) * sizeof(int64_t);
    gnl->offsets = (int64_t *) ngl_alloc_shared(offsets_bytes, "ngl_pairs_offsets");
    if(!gnl->offsets) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, NULL, NULL, "the CSR row offsets", offsets_bytes); return;}

    /* Per-particle scratchpad for fused single-pass build. Each active particle
     * gets a fixed stride (NGL_SCRATCH_STRIDE) of int slots in d_scratch; the BVH
     * walk emits j-indices directly there, with the count tracked in d_counts.
     * After scan + compact we transcribe into the dense CSR neighbors[] array.
     * Memory: stride * num_active * 4 bytes (e.g. 256 * 2M * 4 = 2GB). */
    constexpr int NGL_SCRATCH_STRIDE = 512;
    /* size_t cast required: int * int overflows for num_active > ~4.19M (e.g. fire_m11i
     * gas-per-rank), wrapping to negative int → ~UINT64_MAX after promotion to size_t. */
    size_t na_safe = (size_t)((num_active > 0) ? num_active : 1);
    size_t scratch_bytes = na_safe * (size_t)NGL_SCRATCH_STRIDE * sizeof(int);
    int *d_scratch = (int *) ngl_alloc_device(scratch_bytes, "ngl_pairs_scratch");
    if(!d_scratch) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, NULL, NULL, "the neighbour-walk scratchpad", scratch_bytes); return;}
    size_t counts_bytes = na_safe * sizeof(int);
    int *d_counts  = (int *) ngl_alloc_device(counts_bytes, "ngl_pairs_counts");
    if(!d_counts) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, d_scratch, NULL, "the per-active neighbour counts", counts_bytes); return;}

    /* Drain any prior async GPU work before the passes below. */
    Kokkos::fence();
    /* Fused single pass: BVH walk + write neighbors into per-particle scratchpad */
    {
        sfc_tile_t *tiles = gnl->d_tiles;
        tile_bvh_node_t *bvh = gnl->d_bvh;
        int *pool = gnl->d_pool;
        int *active = gnl->d_active;
        int *scratch = d_scratch;
        int *counts = d_counts;
        int ntiles = gnl->ntiles;
        int bvh_root = gnl->bvh_root;
        int smode = search_mode;

        double sr_fac = search_radius_factor;
        double j_rad_scale = j_kernel_radius_scale;
        const double *radii = d_radii;
        const double *src_pos = d_source_pos;
        const double *compact_xyzh = gnl->d_compact_xyzh;
        Kokkos::parallel_for("ngb_fused", num_active, KOKKOS_LAMBDA(int aa) {
            int i = active[aa];
            double h_i = (radii ? radii[aa] : (double)compact_xyzh[i*4+3]) * sr_fac;
            double pos_i[3];
            if(src_pos) { pos_i[0] = src_pos[aa*3+0]; pos_i[1] = src_pos[aa*3+1]; pos_i[2] = src_pos[aa*3+2]; }
            else        { pos_i[0] = (double)compact_xyzh[i*4+0]; pos_i[1] = (double)compact_xyzh[i*4+1]; pos_i[2] = (double)compact_xyzh[i*4+2]; }
            int cnt = search_neighbors_sfc_gpu(compact_xyzh, pos_i, h_i, j_rad_scale,
                                               tiles, ntiles, pool, smode,
                                               bvh, bvh_root,
                                               &scratch[(size_t)aa * NGL_SCRATCH_STRIDE],
                                               NGL_SCRATCH_STRIDE);
            counts[aa] = cnt;
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("ngb_fused", num_active);

    }

    /* Count overflow particles (count > stride: they need a re-walk in compact phase) */
    int overflow_count = 0;
    {
        int *counts = d_counts;
        Kokkos::parallel_reduce("ngb_overflow_check", num_active,
            KOKKOS_LAMBDA(int aa, int &local) {
                if(counts[aa] > NGL_SCRATCH_STRIDE) local++;
            }, overflow_count);
        Kokkos::fence();
    }

    /* GPU exclusive prefix scan: counts → offsets, returning total.
       Counts are correct even for overflow particles (search_neighbors_sfc_gpu
       returns the true count regardless of bounded write). */
    int64_t total_ll = 0;
    {
        int *counts = d_counts;
        int64_t *offsets = gnl->offsets;
        Kokkos::parallel_scan("ngb_offsets_scan", num_active,
            KOKKOS_LAMBDA(int aa, int64_t &update, const bool final) {
                int64_t v = (int64_t)counts[aa];
                if(final) offsets[aa] = update;
                update += v;
            }, total_ll);
        Kokkos::fence();
    }
    /* Sanity guard: the CSR index (gnl->total_pairs / gnl->offsets) is 64-bit,
     * so INT_MAX is no longer a ceiling. Still abort cleanly on a nonsensical
     * total_pairs (negative => prefix-scan corruption; or a count so large the
     * neighbors allocation byte size would overflow size_t) rather than letting
     * a bad value reach kokkos_malloc / the device compact pass. */
    if(total_ll < 0 || (uint64_t)total_ll > (uint64_t)(SIZE_MAX / sizeof(int))) {
        fprintf(stderr,
            "[NGL FATAL rank=%d] CSR total_pairs implausible: caller=%s num_active=%d "
            "total_pairs=%lld  search_radius_factor=%g j_kernel_radius_scale=%g\n",
            ThisTask, caller_label ? caller_label : "?", num_active, (long long)total_ll,
            search_radius_factor, j_kernel_radius_scale);
        fflush(stderr);
        endrun(915100);
    }
    int64_t total = total_ll;
    gnl->offsets[num_active] = total;
    gnl->total_pairs = total;

    /* Allocate CSR neighbors array (length is 64-bit; element type stays int) */
    size_t neighbors_bytes = (size_t)((total > 0) ? total : 1) * sizeof(int);
    gnl->neighbors = (int *) ngl_alloc_device(neighbors_bytes, "ngl_pairs_neighbors");
    if(!gnl->neighbors) {ngl_build_leave_empty(gnl, num_active, d_radii, d_source_pos, d_scratch, d_counts, "the CSR neighbour list", neighbors_bytes); return;}

    /* Compact: copy from per-particle scratchpad into dense CSR neighbors[]. */
    {
        sfc_tile_t *tiles = gnl->d_tiles;
        tile_bvh_node_t *bvh = gnl->d_bvh;
        int *pool = gnl->d_pool;
        int *active = gnl->d_active;
        int *scratch = d_scratch;
        int *counts = d_counts;
        int64_t *offsets = gnl->offsets;
        int *neighbors = gnl->neighbors;
        int ntiles = gnl->ntiles;
        int bvh_root = gnl->bvh_root;
        int smode = search_mode;
        double sr_fac = search_radius_factor;
        double j_rad_scale = j_kernel_radius_scale;
        const double *radii = d_radii;
        const double *src_pos = d_source_pos;
        const double *compact_xyzh = gnl->d_compact_xyzh;
        Kokkos::parallel_for("ngb_compact", num_active, KOKKOS_LAMBDA(int aa) {
            int n = counts[aa];
            int64_t dst = offsets[aa];
            if(n <= NGL_SCRATCH_STRIDE) {
                size_t src = (size_t)aa * NGL_SCRATCH_STRIDE;
                for(int k = 0; k < n; k++) neighbors[dst + k] = scratch[src + k];
            } else {
                /* Overflow path: re-walk BVH writing directly into neighbors[] */
                int i = active[aa];
                double h_i = (radii ? radii[aa] : (double)compact_xyzh[i*4+3]) * sr_fac;
                double pos_i[3];
                if(src_pos) { pos_i[0] = src_pos[aa*3+0]; pos_i[1] = src_pos[aa*3+1]; pos_i[2] = src_pos[aa*3+2]; }
                else        { pos_i[0] = (double)compact_xyzh[i*4+0]; pos_i[1] = (double)compact_xyzh[i*4+1]; pos_i[2] = (double)compact_xyzh[i*4+2]; }
                search_neighbors_sfc_gpu(compact_xyzh, pos_i, h_i, j_rad_scale,
                                         tiles, ntiles, pool, smode,
                                         bvh, bvh_root,
                                         &neighbors[dst], 0x7fffffff);
            }
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("ngb_compact", num_active);
    }

    /* Free temporaries */
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_scratch);
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_counts);
    if(d_radii) Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);
    if(d_source_pos) Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_source_pos);


    /* Lazy-drift hook (Attack C): drift each neighbor in the freshly-built
     * CSR list to time1 on host. Active particle i is already drifted
     * (move_particles iterated ActiveParticleList). Each pool member j touched
     * by this kernel needs its predicted state (CellP[j].VelPred / Density /
     * InternalEnergyPred / KernelRadius) at time1 before the kernel reads it;
     * drift_particle(j, time1) handles all of that with a single call.
     *
     * drift_particle's "if(time1 == time0) return" early-exit dedups: a j
     * already drifted (e.g. it was in another active i's neighbor list
     * earlier this step, or it IS an active particle) is a fast no-op.
     *
     * Marks h_dirty for the touched j's so that the NEXT gpu_ngb_list_build
     * call's compact_h_refresh updates compact_xyzh[j*4+3] from the freshly-
     * drifted KernelRadius. The CURRENT call's compact_xyzh h field is
     * stale by up to one drift step; the SIDX_H_SLACK inflation in
     * compact_xyzh write paths absorbs that staleness in the BVH tile-overlap
     * test. Per-pair r² acceptance reads the actual P[j].KernelRadius (now
     * freshly drifted), so correctness is preserved. */
    /* The sweep below is dead work when the pool is already current: every
     * drift_particle call returns immediately, and the device->host copy and the
     * h-dirty marking exist only to serve those calls.
     *
     * Two proofs are needed, one per segment of the pool, and both are compared
     * against the time THIS call needs, so anything short of a match falls
     * through to the full sweep. gizmo_full_drift_ti() covers the local
     * particles: move_particles deliberately does not advance it, because it
     * drifts only the active set, which is what makes it a proof rather than a
     * convention. ghost_pool_current_ti() covers the imported segment,
     * established when the owners advanced their particles before packing them.
     *
     * It needs no separate invalidation. A new timestep advances All.Ti_Current,
     * so a certificate from an earlier time simply stops matching.
     *
     * The h-dirty marking stays covered in exactly this case: a full-N drift
     * marks the whole local range as it goes, and ghost slots are marked when
     * they are installed. MEASURED on a production run: on fulldrift steps this
     * is 92.0 billion pool visits with zero members behind, at ~97 s per rank. */
    const integertime t_pool = gizmo_host_ti_current();
    const int ghost_segment_current = (ghost_get_num_ghosts() == 0) ||
                                      (ghost_pool_current_ti() == t_pool);
    const int pool_already_current = (gizmo_full_drift_ti() == t_pool) && ghost_segment_current;
    if(gnl->total_pairs > 0 && gnl->neighbors && !pool_already_current) {
        std::vector<int> ngb_host((size_t)gnl->total_pairs);
        gpu_ngb_copy_neighbors_to_host(gnl, ngb_host.data());
        /* Out-of-line host accessor. Lazy-drift target for CSR
         * neighbors — host-side drift_particle calls. */
        integertime time1 = t_pool;
        /* Ghosts imported for this step were advanced to the current time by
         * their owners before being packed, so the whole imported segment is
         * already current and there is nothing to confirm per ghost. The pool's
         * stamp is compared against the time THIS call needs rather than trusted
         * on its own, so a pool carried over from an earlier time still gets
         * checked particle by particle. */
        const int ghost_start = num_total - ghost_get_num_ghosts();
        const int ghosts_certified = (ghost_pool_current_ti() == time1);
        /* Collect the distinct members that are behind, then advance them in one
         * threaded pass, rather than calling drift_particle once per visit.
         *
         * The drift is real per-particle work -- it runs the implicit
         * thermochemistry solve through set_eos_pressure -- so it must not be
         * strictly serial. A member appears once per PAIR, and two threads
         * testing one particle's Ti_current before either writes would advance
         * it twice, so the distinct set is established first. The stamp is
         * generation-counted and never needs clearing between calls.
         *
         * MEASURED on a production run: this leaves the loop's own rank skew at
         * a tenth of what the per-visit form generated, and that skew was being
         * absorbed by the convergence barrier downstream. */
        static std::vector<unsigned int> pool_seen;
        static unsigned int pool_seen_gen = 0;
        static std::vector<int> pool_behind;
        if((int)pool_seen.size() < num_total) {pool_seen.assign((size_t)num_total, 0u);}
        if(++pool_seen_gen == 0u) {std::fill(pool_seen.begin(), pool_seen.end(), 0u); pool_seen_gen = 1u;}
        pool_behind.clear();
        for(int64_t idx_n = 0; idx_n < gnl->total_pairs; idx_n++) {
            int j = ngb_host[idx_n];
            if(j < 0 || j >= num_total) continue;
            if(ghosts_certified && j >= ghost_start) continue;
            if(pool_seen[(size_t)j] == pool_seen_gen) continue;
            pool_seen[(size_t)j] = pool_seen_gen;
            if(P[j].Ti_current != time1) {pool_behind.push_back(j);}
        }
        {
            const int n_behind = (int)pool_behind.size();
            const int *behind_idx = pool_behind.data();
            drift_particles_batch(behind_idx, n_behind, time1);
        }
        /* Lazy drift just called drift_particle on each j in ngb_host.
         * drift_particle mutates Ti_current, Pos, AND KernelRadius
         * (predict.cc:160,229 — *= exp(divv_fac/N)). Mark h-dirty for both
         * GPU SIDX tracker and host glt cache via the SSOT helper.
         * Promote-to-all kicks in per cache if any cache's bitset popcount
         * exceeds threshold. */
        if(gnl->total_pairs <= (int64_t)INT_MAX) {
            gizmo_mark_kernel_radius_dirty_indices(ngb_host.data(), (int)gnl->total_pairs);
        } else {
            /* >2^31 neighbor pairs: an index list this large is hugely
             * redundant (num_total < 2^31), so truncating n would mark a
             * wrong subset. Escalate to a full-pool mark across both caches. */
            gizmo_mark_kernel_radius_dirty_range(0, num_total);
        }
        /* Move detector baseline past the lazy drift's Ti_current/Pos updates
         * — those are predicted-state setup, not kernel writes that need
         * writeback. Subsequent kernel-side writes to ghost particles will
         * still be flagged by ghost_write_detector_end(). No-op when
         * GIZMO_GPU_ARENA_DEBUG is undefined or detector is inactive. */
        ghost_write_detector_resnapshot_after_lazy_drift();
    }

    /* An index built for this call alone is rebuilt from scratch on the next one, so
     * once the list holds its four device arrays nothing reads the rest of it again:
     * the host tile/pool/BVH mirrors only stage the build, and the position buffers
     * only ever serve the drift refresh of a CACHED index. Hand those four to the list
     * by clearing them here, then release the remainder through the index's own free,
     * so a buffer added to the index later is covered without anyone recalling this
     * site. Without this the host mirrors and position buffers outlive every caller
     * that passes no cached index, which is what made the neighbour-list pool climb
     * without bound across a long run. */
    if(idx == &local_idx) {
        local_idx.d_tiles = NULL; local_idx.d_bvh = NULL;
        local_idx.d_pool  = NULL; local_idx.d_compact_xyzh = NULL;
        gpu_spatial_index_free(&local_idx);
    }

    /* Charge list-build wall, less any kernel time already charged inside it,
     * so the two rows never overlap. */
    cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
}


void gpu_ngb_list_free(gpu_neighbor_list_t *gnl, gpu_spatial_index_t *cached_idx)
{
    if(gnl->neighbors) Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->neighbors);
    if(gnl->offsets)   Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->offsets);
    if(gnl->d_active)  Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->d_active);
    /* Only free tiles/BVH/pool/compact_xyzh if they were NOT from the cached index.
     * Pointers may also be NULL (early-out path with no cache); guard each.
     * d_compact_xyzh in particular was previously leaked here for every
     * non-cached call (~199 MB for an all-types pool, ~73 MB for gas-only).
     * The callers that pass no cached index are merge_and_split_particles, the
     * turbulent power spectra and the two-point correlation. */
    if(!cached_idx || !cached_idx->valid ||
       gnl->d_tiles != cached_idx->d_tiles) {
        if(gnl->d_compact_xyzh) Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->d_compact_xyzh);
        if(gnl->d_pool)  Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->d_pool);
        if(gnl->d_bvh)   Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->d_bvh);
        if(gnl->d_tiles) Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->d_tiles);
    }
}

/* Copy gnl->neighbors (DEVICE_SPACE) into a host buffer. Caller owns host_dest.
   For host-side per-source loops (radfb_local, merge_split, density.cc:868
   HYDRO_VOLUME_CORRECTIONS path, turb_powerspectra, twopoint) that index
   neighbors[] from CPU code. */
void gpu_ngb_copy_neighbors_to_host(const gpu_neighbor_list_t *gnl, int *host_dest)
{
    if(!host_dest || !gnl || gnl->total_pairs <= 0 || !gnl->neighbors) {return;}
    Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        h(host_dest, (size_t)gnl->total_pairs);
    Kokkos::View<const int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        d(gnl->neighbors, (size_t)gnl->total_pairs);
    Kokkos::deep_copy(h, d);
}


/* High-level wrapper: build a symmetric neighbor list on GPU and return it
   in the mymalloc-based neighbor_list_t format expected by gradient/hydro.
   Called from accel.cc (which is NOT compiled by nvcc).
   search_radius_factor: multiplier on KernelRadius (default 1.0; >1 for TURB_DIFF_DYNAMIC). */
void gpu_build_symmetric_neighbor_list(struct particle_data *P_host, int num_total,
                                       int *active_indices, int num_active,
                                       neighbor_list_t *out,
                                       double search_radius_factor)
{
    /* Use the per-step particle arena instead of a dedicated full-NumPart memcpy.
     * The arena's fast path is a no-op when valid (e.g. when gradient/hydro
     * already populated it earlier in the step), avoiding ~2.3s of redundant
     * P-copy on small-N symlist invocations.  Pass the global CellP so the
     * arena's "valid" state remains consistent across mixed P/CellP consumers. */
    gpu_particles_arena_set_site("gpu_build_symmetric_neighbor_list");
    gpu_particles_arena_acquire(num_total, P_host, CellP);
    struct particle_data *P_shared = gpu_particles_arena_P();

    /* Build GPU CSR — share gas-only SIDX with density via the step-persistent cache.
     *
     * search_radius_factor is applied to BOTH the i-side query radius and the
     * j-side kernel radius — a genuinely symmetric scaled search. This repairs
     * the j-side under-search in the shared symlist (hydro-gradient Velocity_hat
     * wide filter under TURB_DIFF_DYNAMIC).
     *
     * RADIUS SEMANTICS: pass EXPLICIT raw per-active radii (P[i].KernelRadius)
     * so search_radius_factor multiplies the RAW kernel radius. With NULL
     * radii the builder would derive h_i from compact_xyzh[i*4+3], which is
     * already slack-inflated (P.KernelRadius * (1+SIDX_H_SLACK)) — compounding
     * the slack into the physics widening factor (e.g. fac=2 -> effective ~3h
     * search, ~27x neighbor volume, CSR-overflow / kernel stall). The runner
     * Spec path already passes explicit fac*raw radii; this matches it. */
    std::vector<double> symlist_raw_radii((num_active > 0) ? (size_t)num_active : 1);
    for(int aa = 0; aa < num_active; aa++) {
        symlist_raw_radii[aa] = (double) P_shared[active_indices[aa]].KernelRadius;
    }
    gpu_neighbor_list_t gpu_nl;
    gpu_ngb_list_build(P_shared, num_total, active_indices, num_active,
                       NGB_SEARCH_SYMMETRIC, 1 /* gas only */, &gpu_nl, gpu_step_sidx_ptr(),
                       search_radius_factor, symlist_raw_radii.data(), NULL, "symlist",
                       search_radius_factor /* j_kernel_radius_scale */);

    /* Copy CSR into mymalloc neighbor_list_t */
    out->num_active = num_active;
    out->total_pairs = gpu_nl.total_pairs;
    out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
    out->neighbors = (int *) mymalloc("ngb_neighbors", (size_t)(gpu_nl.total_pairs > 0 ? gpu_nl.total_pairs : 1) * sizeof(int));
    /* gpu_nl.offsets is SharedSpace (UVM) → host memcpy is fine.
     * gpu_nl.neighbors is DEVICE_SPACE (CudaSpace) → must use deep_copy, not host memcpy. */
    memcpy(out->offsets, gpu_nl.offsets, (size_t)(num_active + 1) * sizeof(int64_t));
    if(gpu_nl.total_pairs > 0) {
        Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            h_neighbors(out->neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::View<const int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            d_neighbors(gpu_nl.neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::deep_copy(h_neighbors, d_neighbors);
    }

    /* Free GPU temporaries (keep tiles/BVH alive — owned by g_step_sidx).
     * Arena is intentionally not released — subsequent gradient/hydro callers
     * benefit from the fast-path skip. */
    gpu_ngb_list_free(&gpu_nl, gpu_step_sidx_ptr());

}


/* Cross-type variant: i-side is the caller's active_indices (any type), j-side
   is filtered by j_type_bitmask. Caller supplies per-active search radii (so
   e.g. KernelRadiusDM or AGS_Hsml can be used instead of P[i].KernelRadius).
   Returns a neighbor_list_t in the mymalloc format, same as the symmetric
   variant — consumers look identical. */
void gpu_build_cross_type_neighbor_list(struct particle_data *P_host, int num_total,
                                        int *i_active_indices, int num_active,
                                        const double *i_search_radii_host,
                                        int j_type_bitmask, int search_mode,
                                        neighbor_list_t *out)
{
    /* Use the per-step particle arena to avoid a redundant full-NumPart memcpy
     * (see gpu_build_symmetric_neighbor_list for rationale). */
    gpu_particles_arena_set_site("gpu_build_cross_type_neighbor_list");
    gpu_particles_arena_acquire(num_total, P_host, CellP);
    struct particle_data *P_shared = gpu_particles_arena_P();

    /* Build GPU CSR with explicit per-i radii and j-side type filter */
    gpu_neighbor_list_t gpu_nl;
    gpu_ngb_list_build(P_shared, num_total, i_active_indices, num_active,
                       search_mode, j_type_bitmask, &gpu_nl, NULL,
                       1.0 /* search_radius_factor */, i_search_radii_host, NULL, "xtype");

    /* Copy CSR into mymalloc neighbor_list_t */
    out->num_active = num_active;
    out->total_pairs = gpu_nl.total_pairs;
    out->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(num_active + 1) * sizeof(int64_t));
    out->neighbors = (int *) mymalloc("ngb_neighbors", (size_t)(gpu_nl.total_pairs > 0 ? gpu_nl.total_pairs : 1) * sizeof(int));
    /* See gpu_build_symmetric_neighbor_list for why neighbors needs deep_copy. */
    memcpy(out->offsets, gpu_nl.offsets, (size_t)(num_active + 1) * sizeof(int64_t));
    if(gpu_nl.total_pairs > 0) {
        Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            h_neighbors(out->neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::View<const int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            d_neighbors(gpu_nl.neighbors, (size_t)gpu_nl.total_pairs);
        Kokkos::deep_copy(h_neighbors, d_neighbors);
    }

    /* Free GPU temporaries.  Arena is intentionally retained for subsequent callers. */
    gpu_ngb_list_free(&gpu_nl, NULL);
}


/* ===================================================================== */
/* Device receiver traversal for request-driven ghost discovery.          */
/* ===================================================================== */
/* The supply-rank half: each received envelope carries a peer's query plus the
 * start nodes that peer's walk reached in THIS rank's tree, and the answer is
 * the set of local particles the query admits.  The host does the same work in
 * mode_b_walk_from_start_nodes + the accept loop that follows it; this is the
 * device form of exactly that, and the two are required to produce the same set.
 *
 * Node geometry comes from the SoA mirror, never the managed Nodes[]/Extnodes[]
 * arrays: streaming those from a kernel is memory-bound to the point of erasing
 * the win.  Leaf fields are staged into a compact array first, for the same
 * reason and because the walk-export path does not build the tile geometry that
 * other device searches read.
 *
 * The device walk cannot drift a stale node the way the host walk does (that
 * needs a lock).  It does not have to: it only runs when the node sweep has
 * certified the whole tree current at this time AND no host lazy drift has
 * happened since, which is the same pair of conditions the device gravity walk
 * uses.  The certification is checked immediately before launch, because the
 * discovery walks themselves arm the lazy-drift latch — a sender walk earlier in
 * the same exchange can withdraw device legality.  Ordering violations therefore
 * make the test fail and route to the host, never corrupt a result. */

/* Per-leaf staged record.  `type` folds in the host's Mass > 0 test (negative =
 * the host would have rejected this particle) and `pool` folds in both the
 * NumPart-when-built bound and the pool-membership lookup (negative = not a
 * supply candidate), so the kernel tests two integers where the host tests four
 * conditions across two arrays. */
struct gx_recv_leaf_t {
    double pos[3];
    int    type;
    int    pool;
};

/* One envelope's traversal.  Returns the number of accepted pool slots, writing
 * the first `cap` of them to `out` (the caller re-runs with the true count when
 * a row overflows its scratch slot).  Mirrors mode_b_walk_impl's three index
 * classes exactly; see mesh/mode_b_local_walker.cc for the host original.
 *
 * `anomaly` reports the one state the host treats as fatal: an index in the gap
 * between the particle slots and the node base, which belongs to neither and
 * means the tree is malformed.  The host stops the run there, so the device
 * cannot simply stop walking -- that would silently truncate an envelope.  It
 * records the state and the caller reproduces the host's stop. */
/* The leaf half of the receiver walk: decide whether a locally-owned particle
 * is a supply candidate this query admits, and record its pool slot.
 *
 * The host filters a leaf twice -- once while walking, with NEAREST_XYZ, and
 * again in the accept pass, with NGB_PERIODIC_BOX_LONG_*.  The two macro
 * families differ only in that the first keeps the sign of the wrapped
 * separation and the second takes its magnitude, so for any separation the
 * search can actually admit they give the same squared distance and the accept
 * form alone reproduces the pair.  The equivalence gate is what establishes
 * that, per pair: it expects the candidate sets to match exactly, and any
 * disagreement has to be shown to be a boundary case rather than assumed to be
 * one.
 *
 * Accepted slots past `cap` are counted but not written, so the caller can
 * re-run the row against a buffer sized to the true count. */
struct GxRecvEmitPairs {
    const struct gx_recv_leaf_t *leaves;
    unsigned int supply_mask;
    int         *out;
    int          cap;
    int          n_found;

    KOKKOS_INLINE_FUNCTION
    void visit(int j, double qx, double qy, double qz, double reach)
    {
        const struct gx_recv_leaf_t &lf = leaves[j];
        if(lf.type < 0 || lf.pool < 0) {return;}
        if(!(supply_mask & (1u << lf.type))) {return;}
        if(gx_pair_accept_wrap_and_test(qx - lf.pos[0], qy - lf.pos[1], qz - lf.pos[2],
                                        reach, 0.0, NGB_SEARCH_ONEWAY)) {
            if(n_found < cap) {out[n_found] = lf.pool;}
            n_found++;
        }
    }
};

KOKKOS_INLINE_FUNCTION
static int gx_recv_walk_one(const struct gx_export_envelope_t &env,
                            unsigned int supply_mask,
                            const Vec3<MyFloat> *node_center,
                            const MyFloat *node_len,
                            const int *node_sibling,
                            const int *node_nextnode,
                            const unsigned int *node_bitflags,
                            const int *nextnode_aux,
                            const struct gx_recv_leaf_t *leaves,
                            int tree_base, int tree_slots, int node_capacity,
                            int foreign_base, int pseudo_start, int num_local,
                            int *anomaly, int *out, int cap)
{
    GxDeviceTreeView tree;
    tree.node_center    = node_center;
    tree.node_len       = node_len;
    tree.node_sibling   = node_sibling;
    tree.node_nextnode  = node_nextnode;
    tree.node_bitflags  = node_bitflags;
    tree.nextnode_aux   = nextnode_aux;
    tree.node_base      = tree_base;
    tree.particle_slots = tree_slots;
    tree.local_particle_slots = num_local;
    tree.node_capacity  = node_capacity;
    tree.foreign_base   = foreign_base;
    tree.pseudo_start   = pseudo_start;

    GxRecvEmitPairs emit;
    emit.leaves      = leaves;
    emit.supply_mask = supply_mask;
    emit.out         = out;
    emit.cap         = cap;
    emit.n_found     = 0;

    gx_device_tree_walk(env, tree, emit, anomaly);
    return emit.n_found;
}

/* Envelopes are processed in fixed-size batches, and each batch's accepted pairs
 * are emitted through one buffer of a fixed size, so neither the scratch nor the
 * output footprint is set by however many envelopes arrived or how many pairs
 * they admit.  Rows are independent and the answer is a set, so splitting
 * changes nothing about the result. */
static const long GX_RECV_BATCH  = 4096;   /* envelopes per walk launch */
static const int  GX_RECV_STRIDE = 512;    /* per-row scratch slots before a re-walk */
/* One budget serves the scratch and the emission buffer; a batch whose pairs do
 * not fit is emitted in several passes rather than abandoned, because declining
 * on a dense batch would send exactly the large-N case this exists to serve back
 * to the host. */
static const long GX_RECV_PAIR_FLOOR = GX_RECV_BATCH * (long)GX_RECV_STRIDE;

/* Exhaustion is reported by returning NULL, so the caller's NULL check decides
   what to do. */
template <class T>
static T *gx_recv_alloc(const char *label, size_t count)
{
    return (T *) gizmo_gpu_alloc_device(count * sizeof(T), label);
}


/* Describe this rank's tree to the device walk, or say why it cannot be
 * described.  Returns 0 with `out` filled, or 1 with `out` untouched.
 *
 * Every caller of the device traversal needs the same answer to the same two
 * questions -- is there a mirror, and does it cover everything a walk can reach
 * -- and the index arithmetic that separates the three index classes is the same
 * arithmetic in every case.  Written once, because a second copy would be a
 * second place for the class boundaries to be got wrong, and a walk that reads
 * one boundary wrong does not fail, it answers short.
 *
 * `local_particle_slots` is the CALLER's, and it is the one thing here that is
 * genuinely per-caller: a walk answering for this rank alone passes the owned
 * count, while one meant to see imported ghosts too passes the full slot count.
 * See mesh/device_tree_walk.h.
 *
 * What this does NOT decide is whether the geometry is CURRENT.  That is a
 * policy question -- who is expected to have drifted the nodes, and what to do
 * when nobody has -- and the answer differs by caller, so each one settles it at
 * its own site rather than inheriting a rule written for another. */
int gx_device_tree_view_build(struct GxDeviceTreeView *out, int local_particle_slots,
                              const char *caller)
{
    /* A tree mirror is required, and it exists whenever a tree does: the build
     * pipeline that fills it runs inside force_treebuild, which every
     * configuration performs because neighbour search needs the tree. */
    const struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa || !soa->center || !soa->len || !soa->sibling || !soa->nextnode ||
       !soa->bitflags || !soa->nextnode_aux ||
       All.TreeNodeIndexBase <= 0 || Numnodestree <= 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d has no usable tree mirror; answering on the host\n",
                   caller, ThisTask);
            fflush(stdout);
        }
        return 1;
    }

    const int tree_base     = All.TreeNodeIndexBase;
    const int tree_slots    = All.TreeParticleSlots;
    const int node_capacity = gpu_gravity_tree_capacity();
    const int foreign_base  = tree_base + MaxNodes;
    const int pseudo_start  = tree_base + MaxNodes + MaxForeignNodes;
    /* The walk may reach any foreign node that was installed, so the mirror has to
     * cover the foreign slots that have storage behind them -- AllocatedForeignNodes,
     * this rank's actual import, which is what the mirror is sized to.  NOT
     * MaxForeignNodes: that is the shared INDEX ceiling used above to place the
     * pseudo-particle region, it is the worst rank's import rather than this one's,
     * and no node is ever installed in the gap between the two, so nothing points
     * there.  Testing the ceiling would decline on every rank whose import is
     * smaller than the largest, which is nearly all of them.  A short mirror is a
     * precondition failure, not a malformed tree: decline, and the host answers --
     * but say so, because a run that quietly answered everything on the host would
     * otherwise look exactly like a run where the device did the work.  This says
     * nothing about whether the geometry is current; that is the caller's. */
    if(node_capacity < MaxNodes + AllocatedForeignNodes ||
       soa->nextnode_aux_size < tree_slots + NTopleaves) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d tree mirror covers %d nodes and %d particle links, short of the %d nodes and %d links the walk can reach; answering on the host\n",
                   caller, ThisTask, node_capacity, soa->nextnode_aux_size,
                   MaxNodes + AllocatedForeignNodes, tree_slots + NTopleaves);
            fflush(stdout);
        }
        return 1;
    }

    /* WIDEN-ON-OPEN inputs (landing 4). The walk re-bounds each node it opens from
       the pair (len, node_ti) plus vmax, instead of requiring a sweep to have
       advanced every node first. Left null if the mirror does not carry them, in
       which case the walk opens on the stored length alone. */
    out->node_vmax            = soa->vmax;
    out->node_ti              = soa->node_ti;
    out->ti_now               = All.Ti_Current;
    gx_walk_drift_tables_refresh();
    out->drift_tables_ok      = g_walk_drift_tables_ok;
    if(g_walk_drift_tables_ok) {out->drift_tables = g_walk_drift_tables;}
    out->node_center          = soa->center;
    out->node_len             = soa->len;
    out->node_sibling         = soa->sibling;
    out->node_nextnode        = soa->nextnode;
    out->node_bitflags        = soa->bitflags;
    out->nextnode_aux         = soa->nextnode_aux;
    out->node_base            = tree_base;
    out->particle_slots       = tree_slots;
    out->local_particle_slots = local_particle_slots;
    out->node_capacity        = node_capacity;
    out->foreign_base         = foreign_base;
    out->pseudo_start         = pseudo_start;
    return 0;
}

/* ---- The touched-set workspace ------------------------------------------
 *
 * Rank-local and persistent, so a call that reaches a thousand leaves pays for a
 * thousand rather than for the whole rank.  The contract is in the header; what
 * follows is the storage and the three operations that use it.
 *
 * The generation stamp is the same one the Mode A neighbour-list hook uses a few
 * hundred lines above (`pool_seen` / `pool_seen_gen`) -- allocated once, never
 * cleared, a wrap re-zeroing it -- moved onto the device because the recorder is
 * a kernel.  Nothing scans it: entries are touched only for leaves a walk
 * actually reaches, which is what keeps this inside the rule that a step with a
 * handful of active particles does no work proportional to the rank. */
static struct GxTouchedSet g_touched_set;

int gx_touched_set_ensure(int local_particle_slots)
{
    if(local_particle_slots <= 0) {return 1;}
    if(g_touched_set.capacity >= local_particle_slots && g_touched_set.seen) {return 0;}

    /* Grown, not resized in place: the stamps say which slots a PREVIOUS
     * generation claimed, and the slots have been re-indexed underneath them by
     * whatever grew the rank.  Re-zeroing and restarting the generation is the
     * honest response; carrying stamps across would let a stale one suppress a
     * particle that genuinely needs drifting. */
    /* Ordering, not cleanup -- the same rule gpu_spatial_index_free states above:
       kokkos_free does not synchronize, so releasing storage a kernel may still
       be reading is a use-after-free. The recording kernels are fenced by
       drift_and_mark before this is ever reached, but the fence belongs WITH the
       release so no future caller has to know that. Skipped when there is
       nothing to release. */
    if(g_touched_set.seen || g_touched_set.list || g_touched_set.counter) {Kokkos::fence();}
    if(g_touched_set.seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.seen);}
    if(g_touched_set.list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.list);}
    if(g_touched_set.counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.counter);}
    g_touched_set = GxTouchedSet{};

    unsigned int *seen    = (unsigned int *) ngl_alloc_shared((size_t)local_particle_slots * sizeof(unsigned int), "touched_set_seen");
    int          *list    = (int *)          ngl_alloc_shared((size_t)local_particle_slots * sizeof(int),          "touched_set_list");
    int          *counter = (int *)          ngl_alloc_shared(sizeof(int),                                         "touched_set_counter");
    if(!seen || !list || !counter) {
        if(seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(seen);}
        if(list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(list);}
        if(counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(counter);}
        return 1;
    }
    for(int k = 0; k < local_particle_slots; k++) {seen[k] = 0u;}
    *counter = 0;

    g_touched_set.seen     = seen;
    g_touched_set.list     = list;
    g_touched_set.counter  = counter;
    g_touched_set.capacity = local_particle_slots;
    g_touched_set.gen      = 0u;
    return 0;
}

/* Lifecycle faults, counted per rank.  A caller that is permanently declining to its own safe
 * route, or a retire that found work nobody drained, is otherwise indistinguishable from the
 * mechanism simply working. */
static long long g_touched_refused_epochs  = 0;
static long long g_touched_retire_faults   = 0;

static void touched_report_once_(const char *what)
{
    static int reported = 0;
    if(reported) {return;}
    reported = 1;
    printf("touched set: task %d: %s\n", ThisTask, what);
    fflush(stdout);
}

/* Open an epoch for `owner`.  Returns 0 when the recorder is now this owner's, nonzero when it
 * could not be handed over -- and a refusal leaves the generation, the cursor and the owner
 * exactly as they were, because the whole point is not to disturb whoever still holds them.
 *
 * THE ADMISSION CONDITION IS THAT NOBODY HOLDS THE RECORDER, and nothing else.  The cursor
 * cannot serve: a consume resets it while the generation deliberately lives on across the
 * remaining passes of the same call, so a zero cursor mid-call means "this pass is drained",
 * never "the recorder is free".  Admitting on a zero cursor would let a second owner take the
 * generation between two passes of a live call, and the first owner's next claim would then
 * arrive out of phase -- which is why the holder, not the cursor, is the question.
 *
 * ⚠ A refusal is NOT self-healing.  The generation is deliberately not advanced, so the
 * previous epoch's stamps still stand and the particles carrying them can no longer re-enter
 * the list.  Every caller therefore has to ACT on a refusal -- the fused walk declines the
 * whole call collectively, the gravity pre-walk takes the full drift -- and none of them may
 * proceed as if the epoch had opened. */
int gx_touched_set_begin_epoch_owned(int owner)
{
    if(!g_touched_set.seen || !g_touched_set.list || !g_touched_set.counter) {
        g_touched_refused_epochs++;
        return 1;
    }
    if(g_touched_set.owner != GX_TOUCHED_OWNER_NONE) {
        g_touched_refused_epochs++;
        touched_report_once_("an epoch was requested while another phase still held the recorder; "
                             "the requester takes its own fallback");
        return 1;
    }
    /* Zero is the never-claimed value, so a wrap has to skip it AND clear the
     * stamps -- otherwise a slot still carrying the old maximum would read as
     * claimed by the new generation and its particle would silently go
     * undrifted. */
    if(++g_touched_set.gen == 0u) {
        for(int k = 0; k < g_touched_set.capacity; k++) {g_touched_set.seen[k] = 0u;}
        g_touched_set.gen = 1u;
    }
    *g_touched_set.counter = 0;
    g_touched_set.owner    = owner;
    Kokkos::memory_fence();   /* the epoch is open before any claim in it is visible */
    return 0;
}

/* The fused walk's entry point.  It holds the recorder for its whole call, across every
 * discovery pass, and retires at the end of that call.  It reports rather than swallows,
 * because its caller (the walk preparation) already has the collective decline that a
 * refusal needs. */
int gx_touched_set_begin_call(void)
{
    return gx_touched_set_begin_epoch_owned(GX_TOUCHED_OWNER_FUSED_WALK);
}

/* End an epoch.  MUST be called after the owner's last consume, which is also what fences the
 * recording kernel -- retiring with a kernel still in flight would leave claims arriving into
 * an epoch that no longer exists.
 *
 * Both faults are loud.  Retiring someone else's epoch means two phases disagree about who is
 * running, and returning quietly would leave the real holder's epoch closed underneath it.  A
 * nonzero cursor means claims were recorded and never drained, i.e. particles the walk reached
 * were never drifted -- the silent staleness this recorder exists to prevent, so it is
 * reported rather than cleared quietly. */
void gx_touched_set_retire(int owner)
{
    if(g_touched_set.owner != owner) {
        g_touched_retire_faults++;
        touched_report_once_("a phase tried to retire a recorder epoch it does not hold");
        return;
    }
    /* Wait for the recording kernel before reading the cursor.  On the ordinary path the
     * owner's last consume already did this, so there is nothing outstanding and the wait
     * costs nothing -- but the check below exists for the path where NO consume ran, and on
     * that path nothing else has waited for the claims.  Reading the cursor there without
     * this would report the epoch as clean whenever the kernel's writes were simply not
     * visible yet, i.e. the detector would be blind in exactly the case it is for.  It must
     * be a full fence rather than a memory fence: the kernel may still be running, so
     * ordering the accesses is not enough, they have to have happened. */
    Kokkos::fence();
    if(g_touched_set.counter && *g_touched_set.counter != 0) {
        g_touched_retire_faults++;
        touched_report_once_("a recorder epoch retired with claims nobody consumed; "
                             "the particles they name were never drifted");
        *g_touched_set.counter = 0;
    }
    g_touched_set.owner = GX_TOUCHED_OWNER_NONE;
    Kokkos::memory_fence();   /* the epoch is closed before the next one can open */
}

long long gx_touched_set_refused_epochs(void) {return g_touched_refused_epochs;}
long long gx_touched_set_retire_faults(void)  {return g_touched_retire_faults;}

struct GxTouchedSet gx_touched_set_view(void) {return g_touched_set;}

/* Give the workspace back.  Called once, at shutdown, before Kokkos is
 * finalized -- Kokkos must not be torn down while an allocation it is tracking
 * is still owned.  Deliberately NOT hung off the tree epoch or the domain
 * decomposition: this storage is persistent on purpose, and freeing it there
 * would turn a once-per-run allocation into per-rebuild churn, which is the
 * cost the generation stamp exists to avoid. */
void gx_touched_set_release(void)
{
    /* Ordering, not cleanup (gpu_spatial_index_free states the rule): kokkos_free
       does not synchronize. This runs on the controlled-stop path as well as the
       normal one, where a kernel may well still be in flight, so incidental
       completion is not an ownership contract. Fence once if anything is held. */
    if(g_touched_set.seen || g_touched_set.list || g_touched_set.counter) {Kokkos::fence();}
    if(g_touched_set.seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.seen);}
    if(g_touched_set.list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.list);}
    if(g_touched_set.counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_touched_set.counter);}
    g_touched_set = GxTouchedSet{};
}

void gx_touched_set_drift_and_mark(integertime time1)
{
    if(!g_touched_set.counter || !g_touched_set.list) {return;}

    /* The recorder is a device kernel writing into shared space, so its writes
     * are not visible to the loop below until it has finished.  Without this the
     * count reads as whatever was there when the launch returned. */
    Kokkos::fence();

    /* The cursor counts every claim, including any the list had no room for, so
     * it is what the append attempted rather than what the list holds.  The two
     * agree unless the generation stamp and the list disagree about how many
     * owned slots there are, which the ensure makes impossible -- but the
     * consequence of trusting the cursor if they ever did is a read past the
     * allocation, here, several kernel launches before the anomaly that records
     * the overflow is ever looked at.  So the allocation bounds the read, and the
     * anomaly stays the thing that reports it. */
    const int claimed = *g_touched_set.counter;
    const int n = (claimed < g_touched_set.capacity) ? claimed : g_touched_set.capacity;
    /* First claims for THIS pass. Reported separately from the pass and call
       counts so a reduction factor is read rather than inferred from a quotient
       whose denominator has to be guessed. */
    if(n <= 0) {return;}

    /* Which of these are actually behind is drift_particles_batch's question and
     * it already answers it -- deciding it here as well would be a second place
     * owning the same test.  So it is asked to hand its compaction back, in
     * place: the recorded list becomes the advanced list, and `n_drifted` is how
     * much of it is live.  The status it returns says only that a controlled stop
     * is pending somewhere, so nothing branches on it. */
    int n_drifted = 0;
    (void) drift_particles_batch(g_touched_set.list, n, time1,
                                 g_touched_set.list, &n_drifted);

    /* drift_particle rescales KernelRadius, so a particle it ADVANCED owes a
     * dirty mark -- and only those.  Marking everything the walk recorded would
     * mark the already-current majority too: measured on a mixed-timebin vehicle,
     * only 19.3% of recorded particles were behind, so that is several times the
     * cache invalidation the work actually justifies, and it can push the dirty
     * tracker over its promote-to-all threshold for nothing. */
    if(n_drifted > 0) {gizmo_mark_kernel_radius_dirty_indices(g_touched_set.list, n_drifted);}

    *g_touched_set.counter = 0;
}


/* Put this rank into the state a fused device walk needs, and describe its tree.
 * Returns 0 with `out` filled, or 1 with the walk declined and the host to answer.
 *
 * A fused walk evaluates the pair kernel at the leaf it just reached, so unlike
 * a discovery walk it reads particle fields, and unlike the host walker it
 * cannot drift a stale one when it gets there -- that needs a lock.  So both the
 * particles and the node geometry have to be current BEFORE the launch, and this
 * is where that is arranged, once, ahead of any discovery round.
 *
 * The particle half is NOT arranged here.  It used to be a full-rank drift, on
 * the argument that discovering the reached set first would mean doing the walk
 * twice and buy nothing.  Both halves of that were wrong at small active counts:
 * a call below ten thousand actives reaches on the order of a thousand leaves
 * while the drift advanced essentially the entire local pool, and the second
 * traversal costs a small fraction of the drift it removes.  Discovery is now
 * what decides which particles are brought current, per pass, at the three
 * evaluation sites in mesh/neighbor_loop_runner.cc -- the two shapes this rank's
 * own self walk takes, and the queries its peers sent.
 *
 * What this preparation still owes the walk is everything that cannot be
 * discovered: the node geometry, which decides where the walk goes and so cannot
 * be repaired from what it reached.  It sweeps the whole tree when nothing else
 * has.  The receiver walk declines instead of sweeping when gravity is compiled
 * in, on the ground that it would drift the entire tree on behalf of a walk that
 * touches part of it; that reasoning applies here too, and the only thing
 * standing against it is that a fused walk cannot take the lock a host walk uses
 * to drift a node when it arrives.  Closing that gap is a separate piece of
 * work.  The one state a sweep cannot repair is a host lazy drift that has
 * already advanced nodes at this time: a sweep skips nodes that are current, so
 * their mirrors stay behind, and the host has to answer. */
int gx_device_fused_walk_prepare(struct GxDeviceTreeView *out, const char *caller)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    const int num_local = ghost_get_num_local();
    if(num_local <= 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d owns no local particles, so it cannot answer on the device; this call falls back for every rank\n",
                   caller, ThisTask);
            fflush(stdout);
        }
        return 1;
    }

    if(gx_device_tree_view_build(out, num_local, caller) != 0) {return 1;}

    /* A walk from the root dereferences the root node before it tests anything,
     * so an empty or half-built tree is not a walk that returns nothing -- it is
     * an unbounded read.  The node capacity cannot catch it, because that is the
     * mirror's allocation rather than the number of nodes actually built, so a
     * built count of zero passes every bound test in the traversal.  Test the
     * built count itself, here, where declining is still free. */
    if(Numnodestree <= 0 || out->node_base >= out->pseudo_start) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d has %d built tree nodes; a walk from the root has nothing to enter, answering on the host\n",
                   caller, ThisTask, Numnodestree);
            fflush(stdout);
        }
        return 1;
    }

    /* ---- Everything above is a CAPABILITY question, answered from state that
     * already exists: does this rank own particles, is there a mirror, does it
     * cover what a walk can reach, is there a root to enter.  Everything below
     * REPAIRS state and is the expensive half.
     *
     * The order matters and is not incidental.  A rank that cannot describe its
     * tree declines, and under the collective vote that decision pulls the WHOLE
     * call back to Mode A -- so any repair performed before the checks is work
     * paid for an answer that is then thrown away, on every rank.  Ask first,
     * repair second.
     *
     * Nothing below invalidates the view built above: the drift neither creates
     * nor destroys particles, and the node sweep writes through the mirror
     * arrays the view already points at rather than reallocating them. */

    /* The walk reads particle fields and cannot drift a stale one when it gets
     * there.  It does NOT follow that the whole rank has to be current: what the
     * walk reads is the leaves it reaches, and a call below ten thousand actives
     * reaches on the order of a thousand of them out of half a million.  So the
     * particles are brought current per discovery pass, against the set that pass
     * actually recorded, at mesh/neighbor_loop_runner.cc's three evaluation sites --
     * which is the same three-stage shape the host backend at those sites has
     * always had, and the same one move_particles, the neighbour-list hook, the
     * Mode B walker and the ghost send-list certify already use.
     *
     * Nothing is voted on here for the particles, and that is deliberate.  Every
     * rank carries the same obligation and discharges it the same way; the drift
     * enters no collective and its only failure report is that a controlled stop
     * is already pending, which is a property of the run rather than of this
     * rank.  A vote would create the divergence it was meant to prevent.
     *
     * The workspace the recording writes into IS a capability, and is arranged
     * here for exactly that reason: it can fail rank-locally, so it belongs where
     * declining is still collective. */
    if(gx_touched_set_ensure(out->local_particle_slots) != 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d could not hold the touched-set workspace for %d local particles; this call falls back for every rank\n",
                   caller, ThisTask, out->local_particle_slots);
            fflush(stdout);
        }
        return 1;
    }
    /* The epoch is a capability too, and it is refused only when another phase still holds the
     * recorder.  Declining here routes it through the same collective vote as a workspace that
     * could not be allocated, so every rank falls back together -- and a refusal must never be
     * walked past, because it leaves the previous epoch's stamps standing. */
    if(gx_touched_set_begin_call() != 0) {
        static int reported = 0;
        if(!reported) {
            reported = 1;
            printf("%s: task %d could not open a touched-set epoch; this call falls back for every rank\n",
                   caller, ThisTask);
            fflush(stdout);
        }
        return 1;
    }

    const int nodes_already_current = gpu_gravity_tree_nodes_current_at(All.Ti_Current) ? 1 : 0;
    /* Nothing dirtied them, so the span the census measures starts here. */
    if(!nodes_already_current) {
        /* A host lazy drift already advanced nodes at this time.  The sweep skips
         * nodes that are current, so their mirrors would stay behind it, and that
         * state cannot be repaired here.  Expected on some calls rather than
         * exceptional -- but reported, because the fallback is collective and an
         * unexplained absence of the device path is indistinguishable from a
         * device path that ran. */
        /* A host lazy drift may have advanced nodes at this time without writing
         * their mirrors, and the ordinary sweep skips exactly those.  Mode-D used
         * to answer that by forcing the variant that rewrites EVERY mirror --
         * ~1.99M of them, to serve a walk that opens a few thousand nodes.
         *
         * It no longer needs to, because the walk now carries its own answer to
         * both classes of staleness:
         *   (a) advanced-but-unmirrored -- the node dirty set recorded exactly
         *       which nodes those are, and repairing them is O(Ndirty), measured
         *       at ~9,200 per rank per span against ~1.99M mirrors;
         *   (b) never-advanced -- widen-on-open re-bounds the node from
         *       (len, node_ti, vmax), so it needs no advancement at all.
         * The two are one design: (b) is what makes (a) sufficient.
         *
         * ⛔ This is NOT "delete the sweep".  Mode-D stops FORCING one; the sweep
         * survives for its other three callers, and the sibling device receiver
         * walk already declines for this exact reason and says so in code
         * (gpu_neighbor_list.cc, "Gravity owns the sweep").  Mode-D is adopting
         * the policy its sibling already has.
         *
         * Any doubt falls back to the old behaviour: an unsafe or overflowed
         * dirty epoch, or a mirror that cannot carry the widening inputs, makes
         * gpu_node_dirty_repair() return nonzero and the full sweep runs. */
        int sweep_rc = 0;
        int _sweep_needed  = 1;
        /* ⛔ Widening needs its inputs. Without them the walk would open on the
         * stored length alone, which is only legal when something else certified
         * the geometry -- so a missing table or mirror DECLINES to the old sweep
         * rather than quietly narrowing the bound. */
        const int widen_armed = (out->drift_tables_ok && out->node_ti && out->node_vmax);
        if(gpu_gravity_tree_nodes_current_at(All.Ti_Current)) {
            _sweep_needed = 0;                       /* fully certified: widening not needed */
        } else if(!widen_armed) {
            _sweep_needed = 1;                       /* cannot widen -> must sweep */
        } else if(gpu_gravity_tree_oneway_safe_at(All.Ti_Current)) {
            _sweep_needed = 0;                       /* already safe to walk */
        } else if(gpu_node_dirty_repair(All.Ti_Current) == 0) {
            _sweep_needed = 0;                       /* O(Ndirty) repair sufficed */
        }
        if(_sweep_needed) {
            { sweep_rc = gpu_force_drift_nodes_ex(All.Ti_Current, /*refresh_mirrors_already_current=*/1); }

            if(sweep_rc != 0) {
                static int reported = 0;
                if(!reported) {
                    reported = 1;
                    printf("%s: task %d could not sweep the node geometry current; this call falls back for every rank\n",
                           caller, ThisTask);
                    fflush(stdout);
                }
                return 1;
            }
            /* The sweep has just rewritten every mirror, so whatever the lazy drift
             * had accumulated is answered and the census span restarts. ⛔ Only on
             * the path that actually swept: the repair path leaves the span alone,
             * because it answered the dirty set rather than the whole tree. */
        }
    }

    return 0;
}

int gx_device_receiver_walk(const struct gx_export_envelope_t *envelopes, long n_env,
                            const int *envelope_peer,
                            unsigned int supply_mask, int search_mode,
                            mode_b_radius_policy_t radius_policy, double j_reach_scale,
                            const int *j_to_pool, int npart_bound,
                            int num_pool, char *matched)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();
    (void)radius_policy; (void)j_reach_scale;

    /* Symmetric search needs the per-type node bands, which are not mirrored to
     * the device yet; that is the next stage of this work.  Until then the host
     * walk answers those callers. */
    if(search_mode != NGB_SEARCH_ONEWAY) {return 1;}
    if(n_env <= 0 || num_pool <= 0 || !matched) {return 1;}

    const int num_local = ghost_get_num_local();
    if(num_local <= 0) {return 1;}

    /* DISPATCH FLOOR.  Staging the local leaves costs O(num_local) whatever the
     * envelopes ask for, while the traversal it enables scales with the envelope
     * count -- so on a step with little to answer the staging is pure overhead,
     * and the corridor already loses the small and middle bins to fixed per-call
     * costs of exactly this kind.  Requiring each envelope to cover a fixed
     * number of staged particles keeps the device path out of that regime.
     * Structural and rank-local: no caller identity, no tuning knob.
     * The ratio is provisional and is what the pricing arm sets. */
    const long GX_RECV_LEAVES_PER_ENVELOPE = 32;
    if(n_env * GX_RECV_LEAVES_PER_ENVELOPE < (long)num_local) {return 1;}

    struct GxDeviceTreeView tree_view;
    if(gx_device_tree_view_build(&tree_view, num_local, "gx_device_receiver_walk") != 0) {return 1;}
    const int tree_base     = tree_view.node_base;
    const int tree_slots    = tree_view.particle_slots;
    const int node_capacity = tree_view.node_capacity;
    const int foreign_base  = tree_view.foreign_base;
    const int pseudo_start  = tree_view.pseudo_start;

    /* The traversal reads node geometry, which drifts, so the nodes have to be
     * current before it runs -- the host walk achieves that by drifting each
     * stale node as it reaches it, under a lock, which a kernel cannot do.
     *
     * Requiring some earlier caller to have swept them is not enough: the sweep
     * lives on the gravity path, so with self-gravity disabled nothing would
     * ever perform it and this traversal could never run at all, on problems
     * whose tree and mirror are perfectly valid.  Perform the sweep here when
     * nothing else has.  It is the same work the host walk would do node by
     * node, done once for the whole tree, and it is rank-local.
     *
     * The one case it cannot repair is a host lazy drift that has already
     * advanced nodes to this time: the sweep skips nodes that are already
     * current and so would leave their mirrors behind. Then, and only then, the
     * host answers. */
    if(!gpu_gravity_tree_nodes_current_at(All.Ti_Current)) {
        /* The geometry is current by neither route -- no sweep certified it and
         * this tree was not built at this time.  A host lazy drift at this time
         * advanced nodes without their mirrors, and a sweep skips already-current
         * nodes, so that state is unrepairable here and the host answers.
         * A tree BUILT after such a drift never reaches this branch: the build
         * rewrites every node and every mirror, and no later host walk can re-arm
         * the latch at this time because force_drift_node returns early on a node
         * that is already current. */
        if(force_host_lazy_drift_ti() == All.Ti_Current) {return 1;}
#ifdef SELFGRAVITY_OFF
        /* No gravity walk exists to sweep the nodes in this build, so without
         * this the traversal could never run at all, however valid the tree and
         * its mirror are. */
        if(gpu_force_drift_nodes(All.Ti_Current) != 0) {return 1;}
#else
        /* Gravity owns the sweep. Sweeping here instead would drift the whole
         * tree eagerly where the walks drift only what they touch, so leave the
         * geometry alone and let the host answer. */
        return 1;
#endif
    }

    using DevSp = GIZMO_KOKKOS_DEVICE_SPACE;

    /* One envelope's accepted pairs are distinct local particles, so a row can
     * never exceed num_local.  Sizing the emission buffer to at least that lets
     * every row fit whole, so no batch is ever abandoned for being dense and the
     * decline paths all sit before anything is written. */
    const long pair_cap = (GX_RECV_PAIR_FLOOR > (long)num_local) ? GX_RECV_PAIR_FLOOR : (long)num_local;

    /* Device buffers.  Allocated once for the whole call and reused across every
     * batch.  A failure here is reported as a decline, not an abort: the caller
     * runs the host walk, and this window holds no collectives so no other rank
     * needs to agree. */
    /* Host buffers first: std::vector throws rather than returning null, and a
     * throw escaping this function would take the rank down between the envelope
     * exchange and the caller's reduction.  Reserving them here turns that into
     * the same decline every other resource failure produces. */
    std::vector<struct gx_recv_leaf_t>     leaf_h;
    std::vector<struct gx_export_envelope_t> env_h;
    std::vector<int>     counts_h;
    std::vector<int64_t> offsets_h;
    std::vector<int>     pairs_h;
    try {
        leaf_h.resize((size_t)num_local);
        env_h.resize((size_t)GX_RECV_BATCH);
        counts_h.resize((size_t)GX_RECV_BATCH);
        offsets_h.resize((size_t)GX_RECV_BATCH);
        pairs_h.resize((size_t)pair_cap);
    } catch(const std::bad_alloc &) {
        printf("gx_device_receiver_walk: task %d could not reserve host staging for %d local leaves; answering on the host\n",
               ThisTask, num_local);
        fflush(stdout);
        return 1;
    }

    struct gx_recv_leaf_t *leaf_d = gx_recv_alloc<struct gx_recv_leaf_t>("gx_recv_leaf", (size_t)num_local);
    struct gx_export_envelope_t *env_d = gx_recv_alloc<struct gx_export_envelope_t>("gx_recv_env", (size_t)GX_RECV_BATCH);
    int     *scratch_d = gx_recv_alloc<int>("gx_recv_scratch", (size_t)GX_RECV_BATCH * GX_RECV_STRIDE);
    int     *counts_d  = gx_recv_alloc<int>("gx_recv_counts",  (size_t)GX_RECV_BATCH);
    int64_t *offsets_d = gx_recv_alloc<int64_t>("gx_recv_offsets", (size_t)GX_RECV_BATCH);
    int     *pairs_d   = gx_recv_alloc<int>("gx_recv_pairs",   (size_t)pair_cap);
    int     *anomaly_d = gx_recv_alloc<int>("gx_recv_anomaly", 1);

    if(!leaf_d || !env_d || !scratch_d || !counts_d || !offsets_d || !pairs_d || !anomaly_d) {
        printf("gx_device_receiver_walk: task %d could not reserve device buffers for %d local leaves; answering on the host\n",
               ThisTask, num_local);
        fflush(stdout);
        if(leaf_d)    {Kokkos::kokkos_free<DevSp>(leaf_d);}
        if(env_d)     {Kokkos::kokkos_free<DevSp>(env_d);}
        if(scratch_d) {Kokkos::kokkos_free<DevSp>(scratch_d);}
        if(counts_d)  {Kokkos::kokkos_free<DevSp>(counts_d);}
        if(offsets_d) {Kokkos::kokkos_free<DevSp>(offsets_d);}
        if(pairs_d)   {Kokkos::kokkos_free<DevSp>(pairs_d);}
        if(anomaly_d) {Kokkos::kokkos_free<DevSp>(anomaly_d);}
        return 1;
    }

    using UmHostI   = Kokkos::View<int*,     Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmHostI64 = Kokkos::View<int64_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmDevI    = Kokkos::View<int*,     DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmDevI64  = Kokkos::View<int64_t*, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    /* Stage the leaf fields the walk reads.  One pass over the local particles,
     * host-side, into a compact record; the AoS is never touched from device. */
    {
#pragma omp parallel for schedule(static)
        for(int j = 0; j < num_local; j++) {
            struct gx_recv_leaf_t rec;
            rec.pos[0] = (double)P[j].Pos[0];
            rec.pos[1] = (double)P[j].Pos[1];
            rec.pos[2] = (double)P[j].Pos[2];
            rec.type   = (P[j].Mass > 0) ? (int)P[j].Type : -1;
            int pool = -1;
            if(j_to_pool && j < npart_bound) {
                const int pp = j_to_pool[j];
                if(pp >= 0 && pp < num_pool) {pool = pp;}
            }
            rec.pool = pool;
            leaf_h[(size_t)j] = rec;
        }
        Kokkos::View<struct gx_recv_leaf_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            lh(leaf_h.data(), (size_t)num_local);
        Kokkos::View<struct gx_recv_leaf_t*, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            ld(leaf_d, (size_t)num_local);
        Kokkos::deep_copy(ld, lh);
    }

    {   /* clear the malformed-tree report */
        int zero = 0;
        Kokkos::View<const int, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> zh(&zero);
        Kokkos::View<int, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>> zd(anomaly_d);
        Kokkos::deep_copy(zd, zh);
    }

    const Vec3<MyFloat>  *node_center   = tree_view.node_center;
    const MyFloat        *node_len      = tree_view.node_len;
    const int            *node_sibling  = tree_view.node_sibling;
    const int            *node_nextnode = tree_view.node_nextnode;
    const unsigned int   *node_bitflags = tree_view.node_bitflags;
    const int            *nextnode_aux  = tree_view.nextnode_aux;

    int status = 0;

    for(long base = 0; base < n_env && status == 0; base += GX_RECV_BATCH) {
        const int nb = (int)((n_env - base < GX_RECV_BATCH) ? (n_env - base) : GX_RECV_BATCH);
        for(int b = 0; b < nb; b++) {env_h[(size_t)b] = envelopes[base + b];}
        {
            Kokkos::View<struct gx_export_envelope_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                eh(env_h.data(), (size_t)nb);
            Kokkos::View<struct gx_export_envelope_t*, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                ed(env_d, (size_t)nb);
            Kokkos::deep_copy(ed, eh);
        }

        /* Walk: bounded write into the row's scratch slot, true count returned. */
        {
            const struct gx_export_envelope_t *env_v = env_d;
            const struct gx_recv_leaf_t *leaf_v = leaf_d;
            int *scratch_v = scratch_d; int *counts_v = counts_d; int *anom_v = anomaly_d;
            const int stride = GX_RECV_STRIDE;
            Kokkos::parallel_for("gx_recv_walk", nb, KOKKOS_LAMBDA(int b) {
                counts_v[b] = gx_recv_walk_one(env_v[b], supply_mask,
                                               node_center, node_len, node_sibling,
                                               node_nextnode, node_bitflags, nextnode_aux,
                                               leaf_v,
                                               tree_base, tree_slots, node_capacity,
                                               foreign_base, pseudo_start, num_local,
                                               anom_v,
                                               &scratch_v[(size_t)b * stride], stride);
            });
            Kokkos::fence();
            gizmo_gpu_check_last_error("gx_recv_walk", nb);
        }

        int64_t batch_total = 0;
        {
            int *counts_v = counts_d; int64_t *offsets_v = offsets_d;
            Kokkos::parallel_scan("gx_recv_scan", nb,
                KOKKOS_LAMBDA(int b, int64_t &update, const bool final) {
                    const int64_t v = (int64_t)counts_v[b];
                    if(final) {offsets_v[b] = update;}
                    update += v;
                }, batch_total);
            Kokkos::fence();
        }
        if(batch_total < 0) {
            printf("gx_device_receiver_walk: task %d prefix sum over %d envelopes produced a negative pair count (%lld)\n",
                   ThisTask, nb, (long long)batch_total);
            fflush(stdout);
            endrun(90001025);
            status = 1; break;
        }
        if(batch_total == 0) {continue;}

        Kokkos::deep_copy(UmHostI(counts_h.data(), (size_t)nb),  UmDevI(counts_d, (size_t)nb));
        Kokkos::deep_copy(UmHostI64(offsets_h.data(), (size_t)nb), UmDevI64(offsets_d, (size_t)nb));

        /* Emit in as many passes as the fixed output buffer needs.  Split points
         * are row boundaries chosen from the counts already in hand, so every
         * pass fits by construction and no batch is ever abandoned for being
         * dense. */
        int r0 = 0;
        while(r0 < nb && status == 0) {
            int r1 = r0; int64_t sub_total = 0;
            while(r1 < nb && sub_total + (int64_t)counts_h[(size_t)r1] <= pair_cap) {
                sub_total += (int64_t)counts_h[(size_t)r1];
                r1++;
            }
            if(r1 == r0) {
                /* The buffer holds num_local, and one envelope cannot admit a
                 * local particle twice -- the exported subtrees are disjoint --
                 * so a row this large means the traversal emitted a duplicate. */
                printf("gx_device_receiver_walk: task %d envelope %ld admits %d pairs against %d local particles; the traversal has emitted a duplicate\n",
                       ThisTask, base + r0, counts_h[(size_t)r0], num_local);
                fflush(stdout);
                endrun(90001026);
                status = 1; break;
            }
            const int64_t sub_base = offsets_h[(size_t)r0];
            {
                const struct gx_export_envelope_t *env_v = env_d;
                const struct gx_recv_leaf_t *leaf_v = leaf_d;
                int *scratch_v = scratch_d; int *counts_v = counts_d;
                int64_t *offsets_v = offsets_d; int *pairs_v = pairs_d; int *anom_v = anomaly_d;
                const int stride = GX_RECV_STRIDE;
                const int rr0 = r0;
                Kokkos::parallel_for("gx_recv_compact", r1 - r0, KOKKOS_LAMBDA(int t) {
                    const int b = rr0 + t;
                    const int n = counts_v[b];
                    if(n <= 0) {return;}
                    const int64_t dst = offsets_v[b] - sub_base;
                    if(n <= stride) {
                        const size_t src = (size_t)b * stride;
                        for(int c = 0; c < n; c++) {pairs_v[dst + c] = scratch_v[src + c];}
                    } else {
                        /* Overflowed its scratch slot: re-walk straight into the
                         * final position. */
                        gx_recv_walk_one(env_v[b], supply_mask,
                                         node_center, node_len, node_sibling,
                                         node_nextnode, node_bitflags, nextnode_aux,
                                         leaf_v,
                                         tree_base, tree_slots, node_capacity,
                                         foreign_base, pseudo_start, num_local,
                                         anom_v, &pairs_v[dst], n);
                    }
                });
                Kokkos::fence();
                gizmo_gpu_check_last_error("gx_recv_compact", r1 - r0);
            }

            Kokkos::deep_copy(UmHostI(pairs_h.data(), (size_t)sub_total),
                              UmDevI(pairs_d, (size_t)sub_total));

            /* Scatter into the bitmap.  It is MPI-facing host memory, and setting
             * a bit twice is a no-op, so the host does it. */
            for(int b = r0; b < r1; b++) {
                const int t = envelope_peer[base + b];
                if(t < 0 || t >= NTask || t == ThisTask) {continue;}
                char *mf = matched + (size_t)t * (size_t)num_pool;
                const int64_t off = offsets_h[(size_t)b] - sub_base;
                const int n = counts_h[(size_t)b];
                for(int c = 0; c < n; c++) {mf[pairs_h[(size_t)(off + c)]] = 1;}
            }
            r0 = r1;
        }
    }

    int anomaly = 0;
    {
        Kokkos::View<int, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ah(&anomaly);
        Kokkos::View<const int, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ad(anomaly_d);
        Kokkos::deep_copy(ah, ad);
    }

    Kokkos::kokkos_free<DevSp>(leaf_d);
    Kokkos::kokkos_free<DevSp>(env_d);
    Kokkos::kokkos_free<DevSp>(scratch_d);
    Kokkos::kokkos_free<DevSp>(counts_d);
    Kokkos::kokkos_free<DevSp>(offsets_d);
    Kokkos::kokkos_free<DevSp>(pairs_d);
    Kokkos::kokkos_free<DevSp>(anomaly_d);

    /* The host walk stops the run on this state, so reaching it here means the
     * same thing.  Falling back would only hide a malformed tree: the host walk
     * would meet it too. */
    if(anomaly) {
        printf("gx_device_receiver_walk: task %d walked into the index gap between the particle slots and the node base, or was handed an incompletely filled tree view; either way the walk cannot answer\n",
               ThisTask);
        fflush(stdout);
        endrun(90001024);
        return 1;
    }
    /* Bits already set when a pass gave up need no undo: the caller falls back to
     * the host walk, which sets exactly the same bits, and setting a bit twice is
     * a no-op. */
    return status;
}


/* Per-TU init function: sets this TU's All_ptr to the shared UVM allocation */
