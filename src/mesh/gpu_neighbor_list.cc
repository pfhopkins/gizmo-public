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

/* SSOT mark helpers — see header for design.  The ghost-exchange supply cache
 * holds membership only and does not track kernel radii, so the GPU SIDX dirty
 * tracker is the sole consumer; a further cache registers here. */
void gizmo_mark_kernel_radius_dirty_indices(const int *indices, int n)
{
    if(!indices || n <= 0) return;
    gpu_dirty_tracker_mark_indices(indices, n);
}
void gizmo_mark_kernel_radius_dirty_range(int start, int end)
{
    if(end <= start) return;
    gpu_dirty_tracker_mark_range(start, end);
}

/* The owned particles' epoch: bumped whenever an index's rows over the owned
 * particles stop describing them in a way the particle count cannot show (see
 * gpu_sidx_notify_owned_changed). gpu_ngb_list_build stamps it into each index at
 * build and requires it to still match before reusing it. The imported particles
 * need no counter here: the ghost exchange already records whether a pool is live
 * and which import produced it (ghost_pool_is_live, ghost_provenance_epoch). */
static uint64_t g_sidx_owned_epoch = 0;

void gpu_sidx_notify_owned_changed(void)
{
    g_sidx_owned_epoch++;
}


/* A new sync point.  The gas index is KEPT: it describes its members as of its reference time, and
 * every walk reads it at the time of the search (sfc_tiles.h), so a drift needs nothing here.  The
 * all-types index is not kept -- only the gas index has its bounds raised when its members are kicked
 * or their motion is written -- so it is released, and the first sink call of the sync point rebuilds it. */
void gpu_step_sidx_invalidate(void)
{
    if(g_step_sidx_alltypes.valid) gpu_spatial_index_free(&g_step_sidx_alltypes);
}

void gpu_step_sidx_invalidate_full(void)
{
    if(g_step_sidx_alltypes.valid) gpu_spatial_index_free(&g_step_sidx_alltypes);
    if(g_step_sidx.valid) gpu_spatial_index_free(&g_step_sidx);
}


/* An index that could not be built must not be left looking usable: the walk reads
 * its tiles, BVH and rows directly. Release what was built and leave the index
 * marked invalid, which is the signal every consumer already tests, so the caller
 * can hand back an empty neighbour list. gpu_spatial_index_free fences before
 * releasing device memory, so nothing is released under a running kernel. The
 * four host build buffers are arena allocations, released in reverse order
 * exactly as the success path does. */
static void sidx_build_leave_invalid(gpu_spatial_index_t *idx, int num_total,
                                     tile_bvh_node_t *h_bvh, double *h_rows, sfc_tile_t *h_tiles, int *h_pool,
                                     const char *what, size_t bytes)
{
    myfree(h_bvh);
    myfree(h_rows);
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

/* How far a tile's box can grow per unit drift interval, relative to its own size or reach.  Once that
 * growth times the interval since the index was built passes SIDX_MAX_BOX_GROWTH, searching the kept
 * index costs more than rebuilding it.  Performance only: a value that lags a raise delays a rebuild and
 * never admits a wrong answer. */
static constexpr double SIDX_MAX_BOX_GROWTH = 1.0;

KOKKOS_INLINE_FUNCTION
double sidx_tile_looseness(const sfc_tile_t &tile)
{
    double w = 0;
    for(int k = 0; k < 3; k++) {
        const double growth = (tile.u_max[k] - tile.u_min[k]) + 2.0 * tile.rho;
        if(!(growth > 0)) {continue;}
        double size = tile.hi[k] - tile.lo[k];
        if(tile.hmax > size) {size = tile.hmax;}
        const double r = (size > 0) ? growth / size : MAX_REAL_NUMBER;
        if(r > w) {w = r;}
    }
    return w;
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

    /* Build SFC tiles + BVH on CPU.  Every member is described as of now, the index's reference time,
     * where the drift puts it then (sfc_tiles.h); the build drifts nothing. */
    const integertime ti_ref = gizmo_host_ti_current();
    const struct DriftKickTableView host_tables = drift_kick_table_view_host();
    sfc_tile_t *h_tiles;
    int *h_pool;
    int num_pool;
    double *h_rows;
    int all_current = 1;
    int ntiles = build_sfc_tiles(P_shared, num_total, type_bitmask, TILE_TARGET_SIZE,
                                 &h_tiles, &h_pool, &num_pool, &h_rows,
                                 ti_ref, &host_tables, &all_current, radius_policy);
    if(ntiles < 0) {
        gpu_spatial_index_free(idx);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "gpu_spatial_index_build (caller '%s'): a member's clock is outside [0, now] or its position "
                 "or velocity is not finite, so no search can bound where it is; spatial index left unbuilt",
                 caller_label ? caller_label : "?");
        gizmo_request_controlled_stop(7740, msg, __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    idx->ntiles = ntiles;

    tile_bvh_node_t *h_bvh;
    int bvh_nnodes = build_tile_bvh(h_tiles, ntiles, &h_bvh);
    idx->bvh_root = bvh_nnodes - 1;

    /* The order a raise that touches most of the index widens it in: nodes by height, leaves first.
     * Every node's index is above its children's, so one forward pass sets the heights. */
    std::vector<int> height((size_t)(bvh_nnodes > 0 ? bvh_nnodes : 1), 0);
    int nlevels = 0;
    for(int nd = 0; nd < bvh_nnodes; nd++) {
        height[nd] = (h_bvh[nd].left < 0) ? 0 : 1 + std::max(height[h_bvh[nd].left], height[h_bvh[nd].right]);
        if(height[nd] + 1 > nlevels) {nlevels = height[nd] + 1;}
    }
    std::vector<int> level_offsets((size_t)nlevels + 1, 0), level_nodes((size_t)(bvh_nnodes > 0 ? bvh_nnodes : 1), 0);
    for(int nd = 0; nd < bvh_nnodes; nd++) {level_offsets[height[nd] + 1]++;}
    for(int L = 0; L < nlevels; L++) {level_offsets[L + 1] += level_offsets[L];}
    {
        std::vector<int> cursor(level_offsets.begin(), level_offsets.end() - 1);
        for(int nd = 0; nd < bvh_nnodes; nd++) {level_nodes[cursor[height[nd]]++] = nd;}
    }
    double looseness = 0;
    for(int t = 0; t < ntiles; t++) {const double w = sidx_tile_looseness(h_tiles[t]); if(w > looseness) {looseness = w;}}

    /* Allocate kernel-read-path arrays in DEVICE_SPACE (CudaSpace HBM on GPU
     * builds, falls back to SharedSpace elsewhere).  This eliminates HMM/TLB-
     * miss overhead on the small-N kernel hot path where one thread does
     * ~1000s of scattered reads through bvh/tiles/pool/rows — that
     * scattered-UVM-access pattern is the suspected source of the residual
     * 1.4s "fused_fnc" floor on 1-active-particle calls.  The host build arrays
     * are transferred via Kokkos::deep_copy through unmanaged-View wrappers
     * (cudaMemcpy under the hood on CUDA builds). */
    int bvh_size = (2 * ntiles - 1);
    if(bvh_size < 1) bvh_size = 1;
    int pool_size = (num_pool > 0) ? num_pool : 1;
    /* Sized to at least one element, as bvh_size and pool_size already are: a type
     * mask that matches nothing would otherwise ask for zero bytes, and a request
     * for none of something cannot be distinguished from a refusal. */
    size_t sidx_tiles_bytes = (size_t)((ntiles > 0) ? ntiles : 1) * sizeof(sfc_tile_t);
    size_t sidx_bvh_bytes   = (size_t) bvh_size * sizeof(tile_bvh_node_t);
    size_t sidx_pool_bytes  = (size_t) pool_size * sizeof(int);
    /* Rows are DOUBLE (positions x,y,z; reach in slot 3): float ABSOLUTE positions are invalid for
     * GIZMO's ~1e11 dynamic range and must NOT decide neighbour inclusion (see §37/§38). */
    size_t sidx_rows_bytes  = (size_t) pool_size * SIDX_ROW_WIDTH * sizeof(double);
    size_t sidx_slot_bytes  = (size_t)((num_total > 0) ? num_total : 1) * sizeof(int);
    size_t sidx_level_bytes = (size_t) bvh_size * sizeof(int);
    size_t sidx_offs_bytes  = (size_t)(nlevels + 1) * sizeof(int);
    idx->d_tiles = (sfc_tile_t *) gizmo_gpu_alloc_device(sidx_tiles_bytes, "ngl_sidx_dev_tiles");
    if(!idx->d_tiles) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the device tile array", sidx_tiles_bytes); return;}
    idx->d_bvh = (tile_bvh_node_t *) gizmo_gpu_alloc_device(sidx_bvh_bytes, "ngl_sidx_dev_bvh");
    if(!idx->d_bvh) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the device tile BVH", sidx_bvh_bytes); return;}
    idx->d_pool = (int *) gizmo_gpu_alloc_device(sidx_pool_bytes, "ngl_sidx_dev_pool");
    if(!idx->d_pool) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the device tile membership pool", sidx_pool_bytes); return;}
    idx->d_compact_xyzh = (double *) gizmo_gpu_alloc_device(sidx_rows_bytes, "ngl_sidx_dev_rows");
    if(!idx->d_compact_xyzh) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the device member rows", sidx_rows_bytes); return;}
    idx->d_slot_of = (int *) gizmo_gpu_alloc_device(sidx_slot_bytes, "ngl_sidx_dev_slot_of");
    if(!idx->d_slot_of) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the device particle-to-slot map", sidx_slot_bytes); return;}
    idx->d_level_nodes = (int *) gizmo_gpu_alloc_device(sidx_level_bytes, "ngl_sidx_dev_level_nodes");
    if(!idx->d_level_nodes) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the device BVH level schedule", sidx_level_bytes); return;}
    idx->h_level_offsets = (int *) gizmo_gpu_alloc_host(sidx_offs_bytes, "ngl_sidx_host_level_offsets");
    if(!idx->h_level_offsets) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the BVH level offsets", sidx_offs_bytes); return;}
    idx->looseness = (double *) gizmo_gpu_alloc_shared(sizeof(double), "ngl_sidx_looseness");
    if(!idx->looseness) {sidx_build_leave_invalid(idx, num_total, h_bvh, h_rows, h_tiles, h_pool, "the looseness scalar", sizeof(double)); return;}

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
        Kokkos::View<double*,            Kokkos::HostSpace, UV>            h_rows_v(h_rows, (size_t) num_pool * SIDX_ROW_WIDTH);
        Kokkos::View<double*,            GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_rows_v(idx->d_compact_xyzh, (size_t) num_pool * SIDX_ROW_WIDTH);
        Kokkos::View<int*,               Kokkos::HostSpace, UV>            h_level_v(level_nodes.data(), bvh_nnodes);
        Kokkos::View<int*,               GIZMO_KOKKOS_DEVICE_SPACE, UV>    d_level_v(idx->d_level_nodes, bvh_nnodes);
        Kokkos::deep_copy(d_tiles_v, h_tiles_v);
        Kokkos::deep_copy(d_bvh_v,   h_bvh_v);
        Kokkos::deep_copy(d_pool_v,  h_pool_v);
        Kokkos::deep_copy(d_rows_v,  h_rows_v);
        Kokkos::deep_copy(d_level_v, h_level_v);
    }
    memcpy(idx->h_level_offsets, level_offsets.data(), (size_t)(nlevels + 1) * sizeof(int));
    *idx->looseness = looseness;

    /* Each particle's pool slot, or -1: how a raise or a reach refresh handed particle indices finds the
     * member's row and tile. */
    {
        int *slot_of = idx->d_slot_of;
        const int *pool = idx->d_pool;
        Kokkos::parallel_for("sidx_slot_of_clear", num_total, KOKKOS_LAMBDA(int i) {slot_of[i] = -1;});
        Kokkos::parallel_for("sidx_slot_of_fill", num_pool, KOKKOS_LAMBDA(int s) {slot_of[pool[s]] = s;});
        Kokkos::fence();
        gizmo_gpu_check_last_error("sidx_slot_of", num_total);
    }

    /* Free the transient mymalloc'd build buffers in proper LIFO order
     * (build_tile_bvh allocated h_bvh last; build_sfc_tiles allocated
     * h_pool, then h_tiles, then h_rows). */
    myfree(h_bvh);
    myfree(h_rows);
    myfree(h_tiles);
    myfree(h_pool);

    idx->bvh_nnodes = bvh_nnodes;
    idx->nlevels = nlevels;
    idx->num_pool = num_pool;
    idx->num_total = num_total;
    idx->cache_tbm = type_bitmask;
    idx->cache_radius_policy = radius_policy;
    idx->ghost_live_when_built       = ghost_pool_is_live();
    idx->ghost_provenance_when_built = ghost_provenance_epoch();
    idx->owned_epoch_when_built      = g_sidx_owned_epoch;
    idx->ti_ref = ti_ref;
    idx->rows_are_positions = all_current;
    idx->rebuild_needed = 0;
    idx->valid = 1;
    /* Register this cache with the dirty tracker over [0, num_total). The rows
     * were written from the live P[] under this cache's radius policy, and
     * nothing between there and here can mutate it, so the range starts clean:
     * the first refresh would recompute values it already holds. */
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
    if(idx->d_compact_xyzh || idx->d_pool || idx->d_bvh || idx->d_tiles || idx->d_slot_of || idx->d_level_nodes ||
       idx->looseness) {
        Kokkos::fence();
    }
    if(idx->d_compact_xyzh) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_compact_xyzh); idx->d_compact_xyzh = NULL;}
    if(idx->d_pool) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_pool); idx->d_pool = NULL;}
    if(idx->d_bvh) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_bvh); idx->d_bvh = NULL;}
    if(idx->d_tiles) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_tiles); idx->d_tiles = NULL;}
    if(idx->d_slot_of) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_slot_of); idx->d_slot_of = NULL;}
    if(idx->d_level_nodes) {Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(idx->d_level_nodes); idx->d_level_nodes = NULL;}
    if(idx->h_level_offsets) {Kokkos::kokkos_free<Kokkos::HostSpace>(idx->h_level_offsets); idx->h_level_offsets = NULL;}
    if(idx->looseness) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(idx->looseness); idx->looseness = NULL;}
    idx->bvh_nnodes = 0;
    idx->nlevels = 0;
    idx->num_pool = 0;
    idx->num_total = 0;
    idx->cache_tbm = -1;
    idx->cache_radius_policy = MODE_B_RADIUS_DEFAULT;
    idx->ti_ref = 0;
    idx->rows_are_positions = 0;
    idx->rebuild_needed = 0;
    idx->valid = 0;
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


/* ---- Keeping the kept index's bounds true (sfc_tiles.h) ----------------------------------------
 * A member's velocity range is re-read when its velocity changes, and its reach when its radius does;
 * the tile and the nodes above it are widened to cover the new value.  Bounds only ever widen.  A
 * raise costs the members it is handed and the paths from their tiles to the root -- or, when those
 * paths together cost more than one sweep of the index, one sweep of it, level by level. */

/* What one member contributes to the bounds of its tile and of the nodes above it. */
struct SidxRaise {
    double u_min[3], u_max[3];
    double rho;
    double hmax;
    double hmax_by_type[TILE_NUM_PTYPES];
};

/* Widen the bounds of a tile or node to cover `m` (a member, a tile or a node).  Called concurrently
 * -- members of one tile, walks up shared ancestors -- so each field is raised atomically; the plain
 * read first skips the atomic when the bound already covers the value.  Returns 1 if a bound moved. */
template <class Box, class Src>
KOKKOS_INLINE_FUNCTION
int sidx_widen(Box *b, const Src &m)
{
    int moved = 0;
    for(int k = 0; k < 3; k++) {
        if(m.u_min[k] < b->u_min[k] && Kokkos::atomic_fetch_min(&b->u_min[k], m.u_min[k]) > m.u_min[k]) {moved = 1;}
        if(m.u_max[k] > b->u_max[k] && Kokkos::atomic_fetch_max(&b->u_max[k], m.u_max[k]) < m.u_max[k]) {moved = 1;}
    }
    if(m.rho  > b->rho  && Kokkos::atomic_fetch_max(&b->rho,  m.rho)  < m.rho)  {moved = 1;}
    if(m.hmax > b->hmax && Kokkos::atomic_fetch_max(&b->hmax, m.hmax) < m.hmax) {moved = 1;}
    for(int t = 0; t < TILE_NUM_PTYPES; t++) {
        if(m.hmax_by_type[t] > b->hmax_by_type[t] &&
           Kokkos::atomic_fetch_max(&b->hmax_by_type[t], m.hmax_by_type[t]) < m.hmax_by_type[t]) {moved = 1;}
    }
    return moved;
}

/* Re-establish, level by level from the leaves, that every node covers its children. */
static void sidx_widen_all_levels(gpu_spatial_index_t *idx)
{
    const sfc_tile_t *tiles = idx->d_tiles;
    tile_bvh_node_t *bvh = idx->d_bvh;
    const int *level_nodes = idx->d_level_nodes;
    for(int L = 0; L < idx->nlevels; L++) {
        /* each level reads the one below it, so it waits for it */
        Kokkos::parallel_for("sidx_widen_level",
                             Kokkos::RangePolicy<>(idx->h_level_offsets[L], idx->h_level_offsets[L + 1]),
                             KOKKOS_LAMBDA(int q) {
            tile_bvh_node_t *node = &bvh[level_nodes[q]];
            if(node->left < 0) {sidx_widen(node, tiles[-(node->left + 1)]);}
            else {sidx_widen(node, bvh[node->left]); sidx_widen(node, bvh[node->right]);}
        });
        Kokkos::fence();
    }
    gizmo_gpu_check_last_error("sidx_widen_level", idx->nlevels);
}

enum { SIDX_RAISE_MOTION = 1, SIDX_RAISE_REACH = 2 };

/* The ratio of the two ways to widen the nodes above n touched members: walking each member's path, or
 * one sweep of the index.  The walk is taken while it costs at most this times the sweep. */
static constexpr double SIDX_PATH_WALK_FACTOR = 1.0;

/* What one particle contributes to the index.  Computed on the host, which reads P[] and CellP[] where
 * their pages live -- a device kernel reading them faults the pages across -- and applied on the device,
 * which owns the bounds.  j < 0 marks a particle that contributes nothing; type < 0, no reach. */
struct SidxRaiseRecord {
    int j, type;
    double r, r_drifted;
    double u_min[3], u_max[3], rho;
};

/* Raise the kept index over particles list[0..n) (host indices; every particle 0..n-1 when list is null).
 * MOTION re-reads each member's velocity range, after its velocity changed; REACH rewrites its row's reaches
 * from its current radius and raises the reach bands, after its radius changed.  A particle that is not a
 * member now -- another type (the index is then rebuilt), or no mass (no pair kernel takes it) -- is
 * skipped.
 * Returns 0; 1 when a member's motion is not finite, so nothing can bound it (the run is stopped); 2 when
 * the raise could not be staged.  Either way the index no longer bounds its members and is marked to be
 * rebuilt by the next list build -- never released here, since lists built from it may still be in use. */
static int sidx_raise_members(gpu_spatial_index_t *idx, const int *list, int n, int what)
{
    if(n <= 0 || !idx->valid) {return 0;}
    const int num_total = idx->num_total, type_bitmask = idx->cache_tbm;
    const mode_b_radius_policy_t policy = idx->cache_radius_policy;
    std::vector<SidxRaiseRecord> rec((size_t)n);
    int fault = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(|:fault)
#endif
    for(int k = 0; k < n; k++) {
        SidxRaiseRecord &r = rec[(size_t)k];
        const int j = list ? list[k] : k;
        r.j = -1; r.type = -1; r.r = 0; r.r_drifted = 0; r.rho = 0;
        for(int d = 0; d < 3; d++) {r.u_min[d] = MAX_REAL_NUMBER; r.u_max[d] = -MAX_REAL_NUMBER;}
        if(j < 0 || j >= num_total) {continue;}
        const int type = (int)P[j].Type;
        if(type < 0 || type >= TILE_NUM_PTYPES || !((1 << type) & type_bitmask) || !(P[j].Mass > 0)) {continue;}
        if((what & SIDX_RAISE_MOTION) && sfc_member_motion_range(j, P, CellP, r.u_min, r.u_max, &r.rho)) {fault |= 1; continue;}
        if(what & SIDX_RAISE_REACH) {
            r.r = nlr_particle_symmetric_radius(P[j], policy);
            r.r_drifted = nlr_particle_symmetric_radius_after_drift(j, P, policy);
            if(r.r_drifted < r.r) {r.r_drifted = r.r;}
            r.type = type;
        }
        r.j = j;
    }
    if(fault) {
        idx->rebuild_needed = 1;
        gizmo_request_controlled_stop(7740, "a particle kept in the neighbour index has a velocity that is not finite, "
                                      "so no search can bound where it is", __FILE__, __LINE__, __FUNCTION__);
        return 1;
    }
    SidxRaiseRecord *d_rec = (SidxRaiseRecord *) ngl_alloc_device((size_t) n * sizeof(SidxRaiseRecord), "sidx_raise_records");
    if(!d_rec) {idx->rebuild_needed = 1; return 2;}
    {
        Kokkos::View<const SidxRaiseRecord*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> hv(rec.data(), (size_t) n);
        Kokkos::View<SidxRaiseRecord*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>> dv(d_rec, (size_t) n);
        Kokkos::deep_copy(dv, hv);
    }
    const double touched_tiles = (n < idx->ntiles) ? (double) n : (double) idx->ntiles;
    const int walk_paths = (touched_tiles * idx->nlevels <= SIDX_PATH_WALK_FACTOR * (double)(idx->ntiles + idx->bvh_nnodes));
    sfc_tile_t *tiles = idx->d_tiles;
    tile_bvh_node_t *bvh = idx->d_bvh;
    const int *slot_of = idx->d_slot_of;
    double *rows = idx->d_compact_xyzh;
    double *looseness = idx->looseness;
    Kokkos::parallel_for("sidx_raise_members", n, KOKKOS_LAMBDA(int k) {
        const SidxRaiseRecord &r = d_rec[k];
        if(r.j < 0) {return;}
        const int slot = slot_of[r.j];
        if(slot < 0) {return;}
        struct SidxRaise m;
        for(int d = 0; d < 3; d++) {m.u_min[d] = r.u_min[d]; m.u_max[d] = r.u_max[d];}
        m.rho = r.rho; m.hmax = 0;
        for(int t = 0; t < TILE_NUM_PTYPES; t++) {m.hmax_by_type[t] = 0;}
        if(r.type >= 0) {
            /* the row holds the reaches themselves, which may fall; the bands only rise */
            rows[(size_t)slot * SIDX_ROW_WIDTH + 3] = r.r;
            rows[(size_t)slot * SIDX_ROW_WIDTH + 4] = r.r_drifted;
            m.hmax = r.r_drifted; m.hmax_by_type[r.type] = r.r_drifted;
        }
        sfc_tile_t *tile = &tiles[slot / TILE_TARGET_SIZE];
        /* A tile that already covers the member needs nothing above it either: whoever raised it is
         * raising its ancestors too, or they were covering it already. */
        if(!sidx_widen(tile, m)) {return;}
        Kokkos::atomic_max(looseness, sidx_tile_looseness(*tile));
        if(!walk_paths) {return;}
        for(int node = tile->bvh_leaf; node >= 0 && sidx_widen(&bvh[node], m); node = bvh[node].parent) {}
    });
    Kokkos::fence();
    gizmo_gpu_check_last_error("sidx_raise_members", n);
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_rec);
    if(!walk_paths) {sidx_widen_all_levels(idx);}
    return 0;
}

/* The kept gas index follows the velocities of its members: particles idx_host[0..n) (host indices)
 * just had their velocity changed.  Called after the kick of the active set, and by every other writer
 * of a particle's velocity through gizmo_motion_bound_raise.  Costs the list and the paths above it. */
void gpu_step_sidx_raise_motion(const int *idx_host, int n)
{
    gpu_spatial_index_t *idx = &g_step_sidx;
    if(!idx->valid || n <= 0 || !idx_host) {return;}
    GIZMO_GPU_ENSURE_ALL_FRESH();
    sidx_raise_members(idx, idx_host, n, SIDX_RAISE_MOTION);
}

/* Hand back a list with no pairs, after a stop has already been requested.  It
 * requests nothing itself, so the first error stays the one reported.  d_active
 * may be left null: the runner reads that as "a stop is pending, return". */
static void ngl_leave_csr_empty(gpu_neighbor_list_t *gnl, int num_active)
{
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
    ngl_leave_csr_empty(gnl, num_active);
    char msg[256];
    snprintf(msg, sizeof(msg),
             "gpu_ngb_list_build: could not allocate %s (%.1f MB) for %d active particles; "
             "neighbour list left empty",
             what, (double) bytes / (1024.0 * 1024.0), num_active);
    gizmo_request_controlled_stop(7711, msg, __FILE__, __LINE__, __FUNCTION__);
}

/* Keep, in each row, exactly the candidates the loop's pair test accepts: a member of the list's types
 * within R_i of query aa for ONEWAY, within max(R_i, j_radius_scale * r_j) for
 * SYMMETRIC, r_j being the member's reach under the loop's radius policy.  Every candidate must already
 * be current.  Rows keep their order and their members' order.  The compaction runs in place in the host
 * copy `ngb`: each row is trimmed on its own, then the kept prefixes move down in row order, which never
 * overwrites a row that has not moved yet; the result replaces the front of the device list. */
static void ngl_trim_rows_to_exact(gpu_neighbor_list_t *gnl, std::vector<int> &ngb, int num_active,
                                   const double *q_pos, const double *q_radius, double radius_factor,
                                   int search_mode, int type_bitmask, mode_b_radius_policy_t radius_policy,
                                   double j_radius_scale, const struct particle_data *Pp)
{
    int64_t *off = gnl->offsets;
    std::vector<int64_t> kept((size_t)num_active);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
#endif
    for(int aa = 0; aa < num_active; aa++) {
        const double R = q_radius[aa] * radius_factor;
        int64_t w = off[aa];
        for(int64_t n = off[aa]; n < off[aa + 1]; n++) {
            const int j = ngb[(size_t)n];
            const struct particle_data &pj = Pp[j];
            if(!((1 << pj.Type) & type_bitmask)) {continue;}
            const double h_j = (search_mode == NGB_SEARCH_ONEWAY) ? 0.0
                             : nlr_particle_symmetric_radius(pj, radius_policy) * j_radius_scale;
            if(!gx_pair_accept_wrap_and_test(q_pos[aa*3+0] - (double)pj.Pos[0], q_pos[aa*3+1] - (double)pj.Pos[1],
                                             q_pos[aa*3+2] - (double)pj.Pos[2], R, h_j, search_mode)) {continue;}
            ngb[(size_t)w++] = j;
        }
        kept[aa] = w - off[aa];
    }
    int64_t total = 0;
    for(int aa = 0; aa < num_active; aa++) {
        const int64_t start = off[aa];
        off[aa] = total;
        if(total != start && kept[aa] > 0) {memmove(&ngb[(size_t)total], &ngb[(size_t)start], (size_t)kept[aa] * sizeof(int));}
        total += kept[aa];
    }
    off[num_active] = total;
    gnl->total_pairs = total;
    if(total > 0) {
        Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(ngb.data(), (size_t)total);
        Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>> d(gnl->neighbors, (size_t)total);
        Kokkos::deep_copy(d, h);
    }
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
     * caller's struct happened to contain. Several callers declare it uninitialised. */
    gnl->d_active = NULL; gnl->offsets = NULL; gnl->neighbors = NULL; gnl->total_pairs = 0;
    double t_entry = my_second(); /* DIAG: entry */
    const double cpu_rows_child0 = CPU_ChildCharged;
    /* A cached index holds only the types it was built for. Walked under a
     * different type_bitmask it would return neighbours outside the caller's
     * mask (e.g. DM in a gas-only density walk, whose lazy drift then stops on
     * Type=1). Refuse it here with a clear message, so a Spec routed to the
     * wrong cache fails at the right layer: request the stop, hand back the
     * empty list and return without walking. The companion warning in
     * gpu_neighbor_list.h on gpu_step_sidx_alltypes_ptr documents this
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
        ngl_leave_csr_empty(gnl, num_active);
        return;
    }
    /* The same refusal for a radius_policy mismatch: a cached row's reach is the
     * build-time policy's, so walking it under another policy would use the
     * wrong leaf-side reach and miss valid pairs. */
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
        ngl_leave_csr_empty(gnl, num_active);
        return;
    }
    /* Early-out: with no active particles there is nothing to search.
     * Skip the SIDX build/refresh AND all kernel launches.  Allocate 1-element
     * stubs so the caller's gpu_ngb_list_free path is well-defined (it always
     * frees neighbors/offsets/d_active). */
    if(num_active == 0) {
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
    /* An index built for this call alone is released on every way out of it; the
     * list it produced holds no part of it. Freeing an index that was never built
     * (the cached route) does nothing. */
    struct ngl_call_index_release {
        gpu_spatial_index_t *owned;
        ~ngl_call_index_release() {gpu_spatial_index_free(owned);}
    } local_idx_release = {&local_idx};
    gpu_spatial_index_t *idx;
    /* Invalidate the cached SIDX unless it still describes the same particles.
     * num_total: the slot map was sized for the old count, so accessing beyond it
     * is UB (ghost exchange redo, particle creation).  The ghost pool's liveness
     * and import: a cleanup-and-reimport can land the SAME ghost count with
     * different ghost contents, which no count test can see.  The owned epoch: a
     * change of membership or a position written outside a drift likewise leaves
     * rows that no longer describe the members.  A drift changes none of these: the kept index is read at the
     * time of the search (sfc_tiles.h), so reuse across a drift is the common path. */
    const integertime t_now = gizmo_host_ti_current();
    const struct DriftKickTableView host_tables = drift_kick_table_view_host();
    /* Whether every pool member is already at the time of this search.  The local
     * particles by the full-drift certificate (move_particles drifts only the
     * active set, so it does not advance it, which is what makes it a proof rather
     * than a convention); the imported segment by its owners having advanced it
     * before packing.  A new timestep advances All.Ti_Current, so a certificate
     * from an earlier time simply stops matching. */
    const int ghost_segment_current = (ghost_get_num_ghosts() == 0) ||
                                      (ghost_pool_current_ti() == t_now);
    const int pool_current = (gizmo_full_drift_ti() == t_now) && ghost_segment_current;
    if(cached_idx && cached_idx->valid &&
       (cached_idx->num_total          != num_total          ||
        cached_idx->ghost_live_when_built       != ghost_pool_is_live()     ||
        cached_idx->ghost_provenance_when_built != ghost_provenance_epoch() ||
        cached_idx->owned_epoch_when_built      != g_sidx_owned_epoch       ||
        cached_idx->ti_ref > t_now || cached_idx->rebuild_needed)) {
        gpu_spatial_index_free(cached_idx);
    }
    /* A kept index is rebuilt instead when a fresh one is the better search: every
     * member is current, so a fresh index reads exactly and its list needs no trim;
     * or its boxes can have grown by more than their own size. */
    if(cached_idx && cached_idx->valid && cached_idx->ti_ref < t_now) {
        const double D_kept = get_drift_factor_impl(cached_idx->ti_ref, t_now, 1.0, &host_tables);
        if(pool_current || !(*cached_idx->looseness * D_kept <= SIDX_MAX_BOX_GROWTH)) {
            gpu_spatial_index_free(cached_idx);
        }
    }
    /* Bring a kept index's reaches up to date.  Every particle whose radius changed
     * since was marked dirty (gizmo_mark_kernel_radius_dirty_*); its row takes its
     * current reach and the bands above it are raised to cover it.  All-dirty, or
     * a refresh that cannot be staged, rebuilds the index instead: a fresh index
     * starts current.  Per-cache state means consuming-and-clearing this cache's
     * bits leaves the other registered caches' bitsets untouched. */
    if(cached_idx && cached_idx->valid && cached_idx->dirty_handle >= 0) {
        const int handle = cached_idx->dirty_handle;
        if(gpu_dirty_tracker_is_all_dirty(handle)) {
            gpu_spatial_index_free(cached_idx);
        } else if(gpu_dirty_tracker_popcount(handle) > 0) {
            /* Drain bitset -> host list -> raise. */
            std::vector<int> dirty_host;
            dirty_host.reserve(gpu_dirty_tracker_popcount(handle));
            gpu_dirty_tracker_consume(handle,
                [](int j, void *ud){ ((std::vector<int> *)ud)->push_back(j); },
                &dirty_host);
            if(sidx_raise_members(cached_idx, dirty_host.data(), (int)dirty_host.size(), SIDX_RAISE_REACH)) {
                gpu_spatial_index_free(cached_idx);
            }
        }
    }
    if(cached_idx && cached_idx->valid) {
        idx = cached_idx;
    } else if(cached_idx) {
        gpu_spatial_index_build(P_shared, num_total, type_bitmask, cached_idx, caller_label, radius_policy);
        idx = cached_idx;
    } else {
        gpu_spatial_index_build(P_shared, num_total, type_bitmask, &local_idx, caller_label, radius_policy);
        idx = &local_idx;
    }
    /* A build that ran out of memory, or met a member no search can bound, leaves
     * the index invalid and has already asked for the stop, naming the cause.
     * There is nothing to walk, so hand back the same empty list any other
     * exhausted allocation here produces. */
    if(!idx->valid) {
        ngl_leave_csr_empty(gnl, num_active);
        return;
    }

    /* How this search reads the index (sfc_walk_frame).  Exact when nothing can
     * have moved since the rows were written: the index was built at this time
     * from members that were all current then, and every member is current now.
     * Every member's reach is its stored one when the pool is current (the rows
     * were brought up to date above); otherwise it is its drifted reach. */
    struct sfc_walk_frame frame;
    frame.D = (idx->ti_ref < t_now) ? get_drift_factor_impl(idx->ti_ref, t_now, 1.0, &host_tables) : 0.0;
    frame.exact = (idx->ti_ref == t_now) && pool_current && idx->rows_are_positions;
    frame.reach_current = pool_current;
    if(!(frame.D >= 0.0 && frame.D < 1.0e30)) {
        char msg[256];
        snprintf(msg, sizeof(msg), "gpu_ngb_list_build (caller '%s'): the drift interval since the neighbour index "
                 "was built is %g, which bounds nothing; neighbour list left empty", caller_label ? caller_label : "?", frame.D);
        gizmo_request_controlled_stop(7740, msg, __FILE__, __LINE__, __FUNCTION__);
        ngl_leave_csr_empty(gnl, num_active);
        return;
    }

    /* Active indices: always re-uploaded (changes per call) */
    size_t active_bytes = (size_t)((num_active > 0) ? num_active : 1) * sizeof(int);
    gnl->d_active = (int *) ngl_alloc_shared(active_bytes, "ngl_pairs_active");
    if(!gnl->d_active) {ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL, "the active-index list", active_bytes); return;}
    memcpy(gnl->d_active, active_indices_host, num_active * sizeof(int));

    /* Each query's radius and position.  A caller may supply either (a loop with a
     * different kernel than the particle's own, e.g. KernelRadiusDM or AGS_Hsml; a
     * source not backed by a P[] entry, e.g. a grid cell).  Otherwise they are the
     * query particle's own at the time of the search -- its reach under the loop's
     * policy, its position once drifted to now (a query is normally an active
     * particle and current already) -- and never the index's rows, which describe
     * members at the index's reference time.  Positions: [aa*3 + k] for axis k. */
    size_t radii_bytes = (size_t)((num_active > 0) ? num_active : 1) * sizeof(double);
    double *d_radii = (double *) ngl_alloc_shared(radii_bytes, "ngl_pairs_radii");
    if(!d_radii) {ngl_build_leave_empty(gnl, num_active, NULL, NULL, NULL, NULL, "the per-active search radii", radii_bytes); return;}
    size_t srcpos_bytes = (size_t)((num_active > 0) ? num_active : 1) * 3 * sizeof(double);
    double *d_source_pos = (double *) ngl_alloc_shared(srcpos_bytes, "ngl_pairs_source_pos");
    if(!d_source_pos) {ngl_build_leave_empty(gnl, num_active, d_radii, NULL, NULL, NULL, "the per-active source positions", srcpos_bytes); return;}
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for(int aa = 0; aa < num_active; aa++) {
        const int i = active_indices_host[aa];
        d_radii[aa] = search_radii_host ? search_radii_host[aa] : nlr_particle_symmetric_radius(P_shared[i], radius_policy);
        if(source_positions_host) {
            for(int k = 0; k < 3; k++) {d_source_pos[aa*3 + k] = source_positions_host[aa*3 + k];}
        } else {
            double c[3], hw = 0.0;
            if(particle_motion_envelope(i, P_shared, CellP, t_now, &host_tables, c, &hw) == PARTICLE_MOTION_UNBOUNDED) {
                for(int k = 0; k < 3; k++) {c[k] = (double)P_shared[i].Pos[k];}
            }
            for(int k = 0; k < 3; k++) {d_source_pos[aa*3 + k] = c[k];}
        }
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
        sfc_tile_t *tiles = idx->d_tiles;
        tile_bvh_node_t *bvh = idx->d_bvh;
        int *pool = idx->d_pool;
        int *scratch = d_scratch;
        int *counts = d_counts;
        int ntiles = idx->ntiles;
        int bvh_root = idx->bvh_root;
        int smode = search_mode;

        double sr_fac = search_radius_factor;
        double j_rad_scale = j_kernel_radius_scale;
        const double *radii = d_radii;
        const double *src_pos = d_source_pos;
        const double *rows = idx->d_compact_xyzh;
        const struct sfc_walk_frame walk_frame = frame;
        Kokkos::parallel_for("ngb_fused", num_active, KOKKOS_LAMBDA(int aa) {
            double h_i = radii[aa] * sr_fac;
            double pos_i[3] = {src_pos[aa*3+0], src_pos[aa*3+1], src_pos[aa*3+2]};
            int cnt = search_neighbors_sfc_gpu(rows, pos_i, h_i, j_rad_scale,
                                               tiles, ntiles, pool, smode,
                                               bvh, bvh_root, walk_frame,
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
        sfc_tile_t *tiles = idx->d_tiles;
        tile_bvh_node_t *bvh = idx->d_bvh;
        int *pool = idx->d_pool;
        int *scratch = d_scratch;
        int *counts = d_counts;
        int64_t *offsets = gnl->offsets;
        int *neighbors = gnl->neighbors;
        int ntiles = idx->ntiles;
        int bvh_root = idx->bvh_root;
        int smode = search_mode;
        double sr_fac = search_radius_factor;
        double j_rad_scale = j_kernel_radius_scale;
        const double *radii = d_radii;
        const double *src_pos = d_source_pos;
        const double *rows = idx->d_compact_xyzh;
        const struct sfc_walk_frame walk_frame = frame;
        Kokkos::parallel_for("ngb_compact", num_active, KOKKOS_LAMBDA(int aa) {
            int n = counts[aa];
            int64_t dst = offsets[aa];
            if(n <= NGL_SCRATCH_STRIDE) {
                size_t src = (size_t)aa * NGL_SCRATCH_STRIDE;
                for(int k = 0; k < n; k++) neighbors[dst + k] = scratch[src + k];
            } else {
                /* Overflow path: re-walk BVH writing directly into neighbors[] */
                double h_i = radii[aa] * sr_fac;
                double pos_i[3] = {src_pos[aa*3+0], src_pos[aa*3+1], src_pos[aa*3+2]};
                search_neighbors_sfc_gpu(rows, pos_i, h_i, j_rad_scale,
                                         tiles, ntiles, pool, smode,
                                         bvh, bvh_root, walk_frame,
                                         &neighbors[dst], 0x7fffffff);
            }
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("ngb_compact", num_active);
    }

    /* Free temporaries.  The staged query radii and positions are kept for the trim below. */
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_scratch);
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(d_counts);


    /* Two things may remain to be done to the walk's list, both on one host copy of it.
     *
     * The members it found that are behind the time of this search are drifted to it
     * (Attack C: lazy drift).  A kernel reads each neighbour's predicted state --
     * CellP[j].VelPred / Density / InternalEnergyPred / KernelRadius -- at this time,
     * and drift_particle(j, time1) provides all of it.  Not needed when the pool is
     * already current (pool_current above): every drift would return at once.
     * MEASURED on a production run: on fulldrift steps this sweep was 92.0 billion
     * pool visits with zero members behind, at ~97 s per rank.
     *
     * Then, unless the walk read the index exactly, its list is a superset -- every
     * member the pair test can accept once drifted -- and each row is trimmed to the
     * members that test accepts at their current positions and reaches. */
    const integertime time1 = t_now;
    const int need_drift = !pool_current;
    const int need_trim = !frame.exact;
    if(gnl->total_pairs > 0 && gnl->neighbors && (need_drift || need_trim)) {
        std::vector<int> ngb_host((size_t)gnl->total_pairs);
        gpu_ngb_copy_neighbors_to_host(gnl, ngb_host.data());
        if(need_drift) {
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
            /* The drift moved the KernelRadius of exactly the members it advanced, so
             * those are the radii to mark dirty; a member it left alone kept its own. */
            gizmo_mark_kernel_radius_dirty_indices(pool_behind.data(), (int)pool_behind.size());
            /* Move detector baseline past the lazy drift's Ti_current/Pos updates
             * — those are predicted-state setup, not kernel writes that need
             * writeback. Subsequent kernel-side writes to ghost particles will
             * still be flagged by ghost_write_detector_end(). No-op when
             * GIZMO_GPU_ARENA_DEBUG is undefined or detector is inactive. */
            ghost_write_detector_resnapshot_after_lazy_drift();
        }
        if(need_trim) {
            ngl_trim_rows_to_exact(gnl, ngb_host, num_active, d_source_pos, d_radii, search_radius_factor,
                                   search_mode, type_bitmask, radius_policy, j_kernel_radius_scale, P_shared);
        }
    }
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_source_pos);

    /* Release this call's own index now, so its fence and frees are charged to this
     * build (the release at scope exit then finds nothing to do). */
    gpu_spatial_index_free(&local_idx);

    /* Charge list-build wall, less any kernel time already charged inside it,
     * so the two rows never overlap. */
    cpu_charge_child(CPU_NGB_BUILD, cpu_minus_children(timediff(t_entry, my_second()), cpu_rows_child0));
}


void gpu_ngb_list_free(gpu_neighbor_list_t *gnl)
{
    if(gnl->neighbors) Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->neighbors);
    if(gnl->offsets)   Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->offsets);
    if(gnl->d_active)  Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->d_active);
    gnl->neighbors = NULL;
    gnl->offsets = NULL;
    gnl->d_active = NULL;
    gnl->num_active = 0;
    gnl->total_pairs = 0;
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
     * so search_radius_factor multiplies the RAW kernel radius, as the runner
     * Spec path does with its explicit fac*raw radii. */
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


    /* Free the GPU list (the kept gas index is not the list's).
     * Arena is intentionally not released — subsequent gradient/hydro callers
     * benefit from the fast-path skip. */
    gpu_ngb_list_free(&gpu_nl);

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
    gpu_ngb_list_free(&gpu_nl);
}


/* ===================================================================== */
/* Device receiver traversal for request-driven ghost discovery.          */
/* ===================================================================== */
/* The supply-rank half: each received envelope carries a peer's query plus the
 * start nodes that peer's walk reached in THIS rank's tree, and the answer is
 * the set of local particles the query admits.  This is the device form of
 * mode_b_walk_from_start_nodes: it finds the particles that may be neighbours once
 * drifted, and both backends hand those to the same accept
 * (gx_send_set_accept_rows), which drifts them and keeps the exact set.
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

/* Per-leaf staged record.  `type` is negative for a particle the receiver can never
 * send -- no mass, or not in the supply pool -- so the kernel tests one integer where
 * the host tests three conditions across two arrays.  `pos` is where the particle is
 * once drifted to the current time, and `half_width` how far from there it can be
 * along any axis (particle_motion_envelope), negative when its motion cannot be
 * bounded.  The leaf test adds its own rounding allowance (motion_envelope_test_slack),
 * so the device's candidates are a superset of what the host's exact test accepts. */
struct gx_recv_leaf_t {
    double pos[3];
    double half_width;
    int    type;
};

/* The leaf half of the receiver walk: record a locally-owned particle that may be a
 * neighbour of the query once drifted to the current time: kept if the box round the
 * position the drift will give it reaches the query sphere; one whose motion cannot be bounded is kept.  This is
 * discovery only -- the candidates go to gx_send_set_accept_rows, which drifts them
 * and makes the exact decision -- so the walk may over-include but never decides.
 *
 * Candidates past `cap` are counted but not written, so the caller can re-run the
 * row against a buffer sized to the true count. */
struct GxRecvCandidates {
    const struct gx_recv_leaf_t *leaves;
    unsigned int supply_mask;
    int         *out;
    int          cap;
    int          n_found;

    KOKKOS_INLINE_FUNCTION
    void visit(int j, double qx, double qy, double qz, double reach)
    {
        const struct gx_recv_leaf_t &lf = leaves[j];
        if(lf.type < 0) {return;}
        if(!(supply_mask & (1u << lf.type))) {return;}
        int keep = 1;
        if(lf.half_width >= 0.0) {
            const double q[3] = {qx, qy, qz};
            const double hw = lf.half_width + motion_envelope_test_slack(lf.pos, q, reach);
            keep = gx_boxpair_overlap_wrap_and_test(lf.pos[0] - qx, lf.pos[1] - qy, lf.pos[2] - qz,
                                                    hw, hw, hw, reach, reach * reach);
        }
        if(keep) {
            if(n_found < cap) {out[n_found] = j;}
            n_found++;
        }
    }
};

/* One envelope's traversal.  Returns the number of candidates, writing the first
 * `cap` of them to `out` (the caller re-runs with the true count when a row
 * overflows its scratch slot).  Mirrors mode_b_walk_impl's three index classes
 * exactly; see mesh/mode_b_local_walker.cc for the host original.
 *
 * `anomaly` reports the one state the host treats as fatal: an index in the gap
 * between the particle slots and the node base, which belongs to neither and
 * means the tree is malformed.  The host stops the run there, so the device
 * cannot simply stop walking -- that would silently truncate an envelope.  It
 * records the state and the caller reproduces the host's stop. */
KOKKOS_INLINE_FUNCTION
static int gx_recv_walk_one(const struct gx_export_envelope_t &env,
                            unsigned int supply_mask,
                            int type_mask_trusted,
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
    tree.type_mask_trusted = type_mask_trusted;   /* read on the host; a device read of the
                                                   * host global would be silently wrong */
    tree.nextnode_aux   = nextnode_aux;
    tree.node_base      = tree_base;
    tree.particle_slots = tree_slots;
    tree.local_particle_slots = num_local;
    tree.node_capacity  = node_capacity;
    tree.foreign_base   = foreign_base;
    tree.pseudo_start   = pseudo_start;

    GxRecvCandidates emit;
    emit.leaves      = leaves;
    emit.supply_mask = supply_mask;
    emit.out         = out;
    emit.cap         = cap;
    emit.n_found     = 0;

    gx_device_tree_walk(env, tree, emit, anomaly, supply_mask);
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
    out->type_mask_trusted    = TypePresenceMaskTrusted;
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

void gx_touched_set_begin_call(void)
{
    /* Zero is the never-claimed value, so a wrap has to skip it AND clear the
     * stamps -- otherwise a slot still carrying the old maximum would read as
     * claimed by the new generation and its particle would silently go
     * undrifted. */
    if(++g_touched_set.gen == 0u) {
        for(int k = 0; k < g_touched_set.capacity; k++) {g_touched_set.seen[k] = 0u;}
        g_touched_set.gen = 1u;
    }
    if(g_touched_set.counter) {*g_touched_set.counter = 0;}
}

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

/* ============================================================================
 * The motion-target set and the raise it feeds (declared in gpu_neighbor_list.h).
 * ========================================================================== */
static struct GxMotionTargetSet g_motion_targets;
static int g_motion_targets_armed = 0;

int gx_motion_target_ensure(int local_particle_slots)
{
    if(local_particle_slots <= 0) {return 1;}
    if(g_motion_targets.capacity >= local_particle_slots && g_motion_targets.seen) {return 0;}
    gx_motion_target_release();
    unsigned int *seen    = (unsigned int *) ngl_alloc_shared((size_t)local_particle_slots * sizeof(unsigned int), "motion_target_seen");
    int          *list    = (int *)          ngl_alloc_shared((size_t)local_particle_slots * sizeof(int),          "motion_target_list");
    int          *counter = (int *)          ngl_alloc_shared(sizeof(int),                                         "motion_target_counter");
    if(!seen || !list || !counter) {
        if(seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(seen);}
        if(list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(list);}
        if(counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(counter);}
        return 1;
    }
    for(int k = 0; k < local_particle_slots; k++) {seen[k] = 0u;}
    *counter = 0;
    g_motion_targets.seen = seen; g_motion_targets.list = list; g_motion_targets.counter = counter;
    g_motion_targets.capacity = local_particle_slots; g_motion_targets.gen = 0u;
    return 0;
}

void gx_motion_target_begin_call(void)
{
    if(++g_motion_targets.gen == 0u) {
        for(int k = 0; k < g_motion_targets.capacity; k++) {g_motion_targets.seen[k] = 0u;}
        g_motion_targets.gen = 1u;
    }
    if(g_motion_targets.counter) {*g_motion_targets.counter = 0;}
}

struct GxMotionTargetSet gx_motion_target_view(void) {return g_motion_targets;}

void gx_motion_target_mark_host(int j)
{
    /* Serial host callers only (the writeback apply loop, a host module's own
     * loop); the device mark is the atomic form. */
    struct GxMotionTargetSet &ts = g_motion_targets;
    if(!ts.seen || j < 0 || j >= ts.capacity) {return;}
    if(ts.seen[j] == ts.gen) {return;}
    ts.seen[j] = ts.gen;
    const int slot = (*ts.counter)++;
    if(slot < ts.capacity) {ts.list[slot] = j;}
}

void gx_motion_target_set_armed(int armed) {g_motion_targets_armed = armed;}
int  gx_motion_target_armed(void) {return g_motion_targets_armed;}

void gx_motion_target_consume(void)
{
    if(!g_motion_targets.counter || !g_motion_targets.list) {return;}
    Kokkos::fence();   /* the markers may be device kernels */
    const int claimed = *g_motion_targets.counter;
    const int n = (claimed < g_motion_targets.capacity) ? claimed : g_motion_targets.capacity;
    if(n > 0) {gizmo_motion_bound_raise(g_motion_targets.list, n);}
    *g_motion_targets.counter = 0;
}

void gx_motion_target_release(void)
{
    if(g_motion_targets.seen || g_motion_targets.list || g_motion_targets.counter) {Kokkos::fence();}
    if(g_motion_targets.seen)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_motion_targets.seen);}
    if(g_motion_targets.list)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_motion_targets.list);}
    if(g_motion_targets.counter) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(g_motion_targets.counter);}
    g_motion_targets = GxMotionTargetSet{};
}

void gizmo_motion_bound_raise(const int *idx, int n)
{
    if(n <= 0 || !idx) {return;}
    /* The tree records which top-level nodes changed, for the exchange at the
     * next tree-update phase. */
    gravity_note_motion_bound(idx, n);
    /* The kept gas neighbour index follows the same velocities. */
    gpu_step_sidx_raise_motion(idx, n);
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
    gx_touched_set_begin_call();

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
                            mode_b_radius_policy_t radius_policy,
                            double j_radius_scale, double safety_factor,
                            const int *j_to_pool, int npart_bound,
                            int num_pool, struct ghost_send_set *send_set)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* ONEWAY only.  A symmetric search accepts a pair on the NEIGHBOUR's radius as
     * well as the query's, so it has to bound, per node, how far the particles of
     * each type below it reach -- Extnodes[].hmax_per_type[].  Those bands are host
     * AoS only.  The per-node type-presence bits the device walk does carry answer
     * whether a type is there, which is enough to skip a node but not to decide
     * whether one of its particles reaches back, so they do not close this gap.
     * Symmetric callers are answered by the host walk. */
    if(search_mode != NGB_SEARCH_ONEWAY) {return GX_RECEIVER_DECLINED;}
    if(n_env <= 0 || num_pool <= 0 || !send_set) {return GX_RECEIVER_DECLINED;}

    const int num_local = ghost_get_num_local();
    if(num_local <= 0) {return GX_RECEIVER_DECLINED;}

    /* DISPATCH FLOOR.  Staging the local leaves costs O(num_local) whatever the
     * envelopes ask for, while the traversal it enables scales with the envelope
     * count -- so on a step with little to answer the staging is pure overhead,
     * and the corridor already loses the small and middle bins to fixed per-call
     * costs of exactly this kind.  Requiring each envelope to cover a fixed
     * number of staged particles keeps the device path out of that regime.
     * Structural and rank-local: no caller identity, no tuning knob.
     * The ratio is provisional and is what the pricing arm sets. */
    const long GX_RECV_LEAVES_PER_ENVELOPE = 32;
    if(n_env * GX_RECV_LEAVES_PER_ENVELOPE < (long)num_local) {return GX_RECEIVER_DECLINED;}

    struct GxDeviceTreeView tree_view;
    if(gx_device_tree_view_build(&tree_view, num_local, "gx_device_receiver_walk") != 0) {return GX_RECEIVER_DECLINED;}
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
        if(force_host_lazy_drift_ti() == All.Ti_Current) {return GX_RECEIVER_DECLINED;}
#ifdef SELFGRAVITY_OFF
        /* No gravity walk exists to sweep the nodes in this build, so without
         * this the traversal could never run at all, however valid the tree and
         * its mirror are. */
        if(gpu_force_drift_nodes(All.Ti_Current) != 0) {return GX_RECEIVER_DECLINED;}
#else
        /* Gravity owns the sweep. Sweeping here instead would drift the whole
         * tree eagerly where the walks drift only what they touch, so leave the
         * geometry alone and let the host answer. */
        return GX_RECEIVER_DECLINED;
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
    std::vector<struct gx_candidate_row> rows_h;
    try {
        leaf_h.resize((size_t)num_local);
        env_h.resize((size_t)GX_RECV_BATCH);
        counts_h.resize((size_t)GX_RECV_BATCH);
        offsets_h.resize((size_t)GX_RECV_BATCH);
        pairs_h.resize((size_t)pair_cap);
        rows_h.resize((size_t)GX_RECV_BATCH);
    } catch(const std::bad_alloc &) {
        printf("gx_device_receiver_walk: task %d could not reserve host staging for %d local leaves; answering on the host\n",
               ThisTask, num_local);
        fflush(stdout);
        return GX_RECEIVER_DECLINED;
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
        return GX_RECEIVER_DECLINED;
    }

    using UmHostI   = Kokkos::View<int*,     Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmHostI64 = Kokkos::View<int64_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmDevI    = Kokkos::View<int*,     DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UmDevI64  = Kokkos::View<int64_t*, DevSp, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    int status_staging_failed = 0;
    /* Stage the leaf fields the walk reads.  One pass over the local particles,
     * host-side, into a compact record; the AoS is never touched from device.  The
     * position a particle that is behind will be drifted to is taken here, from the host
     * tables, so the kernel reads the record instead of the particle state and the tables. */
    {
        const integertime ti_now = All.Ti_Current;
        const struct DriftKickTableView drift_tables = drift_kick_table_view_host();
        long n_unmapped = 0;
#pragma omp parallel for schedule(static) reduction(+:n_unmapped)
        for(int j = 0; j < num_local; j++) {
            struct gx_recv_leaf_t rec;
            int in_pool = 0;
            if(j_to_pool && j < npart_bound) {
                const int pp = j_to_pool[j];
                in_pool = (pp >= 0 && pp < num_pool);
            }
            /* The pool holds every particle of positive mass (gx_send_set_accept_rows), so one
               missing from it is a stale or corrupt map, stopped below rather than hidden. */
            if(P[j].Mass > 0 && !in_pool) {n_unmapped++;}
            rec.type = (P[j].Mass > 0 && in_pool) ? (int)P[j].Type : -1;
            rec.pos[0] = (double)P[j].Pos[0]; rec.pos[1] = (double)P[j].Pos[1]; rec.pos[2] = (double)P[j].Pos[2];
            rec.half_width = 0.0;
            if(rec.type >= 0) {
                double hw = 0.0;
                const int motion = particle_motion_envelope(j, P, CellP, ti_now, &drift_tables, rec.pos, &hw);
                rec.half_width = (motion == PARTICLE_MOTION_UNBOUNDED) ? -1.0 : hw;
            }
            leaf_h[(size_t)j] = rec;
        }
        if(n_unmapped > 0) {
            printf("gx_device_receiver_walk: task %d has %ld particles of positive mass with no supply-pool slot\n",
                   ThisTask, n_unmapped);
            fflush(stdout);
            gizmo_request_controlled_stop(7738, "gx_device_receiver_walk: supply map missing particles of positive mass",
                                          __FILE__, __LINE__, __FUNCTION__);
            status_staging_failed = 1;
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

    int status = status_staging_failed ? GX_RECEIVER_FAILED : GX_RECEIVER_COMPLETED;

    for(long base = 0; base < n_env && status == GX_RECEIVER_COMPLETED; base += GX_RECV_BATCH) {
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
            const int type_mask_trusted_v = TypePresenceMaskTrusted;   /* captured by value */
            Kokkos::parallel_for("gx_recv_walk", nb, KOKKOS_LAMBDA(int b) {
                counts_v[b] = gx_recv_walk_one(env_v[b], supply_mask, type_mask_trusted_v,
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
            status = GX_RECEIVER_FAILED; break;
        }
        if(batch_total == 0) {continue;}

        Kokkos::deep_copy(UmHostI(counts_h.data(), (size_t)nb),  UmDevI(counts_d, (size_t)nb));
        Kokkos::deep_copy(UmHostI64(offsets_h.data(), (size_t)nb), UmDevI64(offsets_d, (size_t)nb));

        /* Emit in as many passes as the fixed output buffer needs.  Split points
         * are row boundaries chosen from the counts already in hand, so every
         * pass fits by construction and no batch is ever abandoned for being
         * dense. */
        int r0 = 0;
        while(r0 < nb && status == GX_RECEIVER_COMPLETED) {
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
                status = GX_RECEIVER_FAILED; break;
            }
            const int64_t sub_base = offsets_h[(size_t)r0];
            {
                const struct gx_export_envelope_t *env_v = env_d;
                const struct gx_recv_leaf_t *leaf_v = leaf_d;
                int *scratch_v = scratch_d; int *counts_v = counts_d;
                int64_t *offsets_v = offsets_d; int *pairs_v = pairs_d; int *anom_v = anomaly_d;
                const int stride = GX_RECV_STRIDE;
                const int rr0 = r0;
                const int type_mask_trusted_v = TypePresenceMaskTrusted;   /* captured by value */
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
                        gx_recv_walk_one(env_v[b], supply_mask, type_mask_trusted_v,
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
            /* Hand each row's candidates to the shared accept, which drifts them and
             * keeps the exact set.  Envelopes arrive grouped by peer in ascending
             * order, so a peer that spans a pass or a batch boundary simply carries
             * on in the next call. */
            long n_rows = 0;
            for(int b = r0; b < r1; b++) {
                const int t = envelope_peer[base + b];
                if(t < 0 || t >= NTask) {
                    /* Every envelope's sender is known; one that is not would have
                     * its candidates silently dropped, so stop rather than under-include. */
                    printf("gx_device_receiver_walk: task %d envelope %ld names sender %d, outside 0..%d\n",
                           ThisTask, base + b, t, NTask - 1);
                    fflush(stdout);
                    gizmo_request_controlled_stop(7737, "gx_device_receiver_walk: envelope with no valid sender",
                                                  __FILE__, __LINE__, __FUNCTION__);
                    status = GX_RECEIVER_FAILED;
                    break;
                }
                if(t == ThisTask) {continue;}
                const int64_t off = offsets_h[(size_t)b] - sub_base;
                struct gx_candidate_row row = {&env_h[(size_t)b], t, &pairs_h[(size_t)off], counts_h[(size_t)b]};
                rows_h[(size_t)n_rows++] = row;
            }
            if(status == GX_RECEIVER_COMPLETED &&
               gx_send_set_accept_rows(send_set, rows_h.data(), n_rows, search_mode,
                                       radius_policy, j_radius_scale, safety_factor) != 0) {
                status = GX_RECEIVER_FAILED;
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
        return GX_RECEIVER_FAILED;
    }
    /* A pass that gave up reports FAILED rather than declining: what it already
     * handed to the send set cannot be taken back by a host rerun. */
    return status;
}


/* Per-TU init function: sets this TU's All_ptr to the shared UVM allocation */
