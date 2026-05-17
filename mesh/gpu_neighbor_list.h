/* gpu_neighbor_list.h — GPU-accelerated neighbor list construction.
 *
 * Provides gpu_ngb_list_build(): builds a CSR neighbor list using
 * SFC tiles + BVH spatial index with GPU-parallel neighbor search.
 * The tile/BVH construction runs on CPU; the per-particle search
 * runs on GPU via Kokkos parallel_for.
 *
 * Reusable for any neighbor list type: density (one-way, gas-only),
 * symmetric (max(h_i,h_j)), or cross-type (different type_bitmask).
 * The spatial index can optionally be cached and reused across calls
 * with the same particle positions.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef GPU_NEIGHBOR_LIST_H
#define GPU_NEIGHBOR_LIST_H

#include <stdint.h>
#include "sfc_tiles.h"

/* GPU-resident neighbor list: CSR arrays.
   offsets: SharedSpace (host writes offsets[num_active]=total after scan).
   neighbors: DeviceSpace (GPU HBM on CUDA; never host-accessed directly). */
struct gpu_neighbor_list_t {
    int *offsets;       /* [num_active+1] in SharedSpace */
    int *neighbors;     /* [total_pairs] in DeviceSpace (GPU HBM — no UVM fault) */
    int num_active;
    int total_pairs;

    /* Device-resident copies of spatial index data */
    sfc_tile_t *d_tiles;
    tile_bvh_node_t *d_bvh;
    int *d_pool;
    int *d_active;
    int ntiles;
    int bvh_root;

    /* Periodicity parameters (copied to avoid global access in kernels) */
    int periodic_flags[3];
    double box_sizes[3];
    double box_halves[3];

    /* Compact position+h array for BVH traversal cache efficiency (points into
       spatial index memory; do NOT free from gnl — owned by gpu_spatial_index_t). */
    float *d_compact_xyzh;
};


/* Cached spatial index: tiles + BVH in SharedSpace.
   Depends only on particle positions (not KernelRadius), so can be
   reused across multiple neighbor list builds with different search radii. */
struct gpu_spatial_index_t {
    sfc_tile_t *d_tiles;
    tile_bvh_node_t *d_bvh;
    int *d_pool;
    int ntiles;
    int bvh_root;
    int periodic_flags[3];
    double box_sizes[3];
    double box_halves[3];
    /* Compact position+h array: d_compact_xyzh[i*4+0..3] = x,y,z,h for particle i.
       Fits ~32MB for 2M particles vs ~800MB for full P[], dramatically improving
       GPU BVH traversal cache efficiency. Built in gpu_spatial_index_build. */
    float *d_compact_xyzh;
    int num_total;  /* particle count when built; mismatch → invalidate */
    int valid;  /* 1 if built and usable */
    int dirty_handle = -1; /* gpu_dirty_tracker handle; -1 when not registered */
    /* Type bitmask used at the most recent build. The cached compact_xyzh /
     * pool only includes the originally-built types; later callers MUST pass
     * the same tbm to gpu_ngb_list_build. Codex 2026-05-12: density port
     * exposed a bug where a gas-only caller (tbm=1) used the all-types cache
     * (tbm=0x3f), causing the walker to return DM neighbors and lazy-drift
     * to abort on Type=1 particles. gpu_ngb_list_build now hard-aborts on
     * mismatch instead of silently mis-walking. Default -1 = "no build yet". */
    int cache_tbm = -1;
    /* Host-side persistent copies kept alive across drifts to support
     * gpu_step_sidx_invalidate's incremental refresh path: tile bboxes
     * are recomputed in place from current particle positions, the BVH
     * is re-fitted from the updated tile bboxes (build_tile_bvh re-call,
     * O(ntiles)), then re-staged to device. Skips the expensive SFC
     * sort (build_sfc_tiles ~1s on 12.4M particles) which is the
     * dominant cost in a fresh build. Only freed at full invalidate
     * (post-domain_decomp boundary). */
    sfc_tile_t *h_tiles;        /* [ntiles] */
    int *h_pool;                /* [num_pool] */
    tile_bvh_node_t *h_bvh;     /* [2*ntiles-1] */
    int h_bvh_nnodes;
    int num_pool;
    /* Per-tile original extent (max axis range at last full rebuild). Used to
     * monitor cumulative bbox dispersion across drifts; if it grows pathologically
     * the next gpu_step_sidx_invalidate_full() resets it. */
    double *h_tile_orig_max_extent; /* [ntiles] */
    /* Host + device position-staging buffers used by the drift refresh path to
     * bypass the UVM-fault storm of a device-side parallel_for over P_shared.Pos.
     * The bbox-recompute pass on host fills h_pos_buf for every pool member;
     * h_pos_buf -> d_pos_buf is one bulk deep_copy (~200MB at NVLink ~600GB/s
     * = sub-ms); then a tiny device-side scatter kernel writes into the
     * interleaved d_compact_xyzh[i*4+0..2]. Non-pool entries are unused by the
     * BVH walk so their stale h_pos_buf values are harmless. */
    float *h_pos_buf;           /* [3*num_total] in Kokkos::HostSpace */
    float *d_pos_buf;           /* [3*num_total] in DEVICE_SPACE */
};


/* Build spatial index (tiles + BVH) on CPU, copy to SharedSpace.
   P_shared must be in SharedSpace (managed memory).
   caller_label: short tag printed in DIAG_SIDX so we can attribute rebuilds. */
void gpu_spatial_index_build(struct particle_data *P_shared, int num_total,
                             int type_bitmask, gpu_spatial_index_t *idx,
                             const char *caller_label = "?");

/* Free spatial index SharedSpace memory. */
void gpu_spatial_index_free(gpu_spatial_index_t *idx);

/* Module-level persistent SIDX for gas-only (type_bitmask=1) neighbor builds.
 * Shared across density rounds + symlist within a single step so the BVH +
 * compact_xyzh build only happens once per step instead of 4× per step.
 * Must be invalidated after drift (positions change) — caller (run.cc) calls
 * gpu_step_sidx_invalidate() after find_next_sync_point_and_drift(). */
gpu_spatial_index_t *gpu_step_sidx_ptr(void);

/* Module-level persistent SIDX for all-types (type_bitmask=0x3f) builds.
 * Specifically for the SINK_PARTICLE codepath (sink_env1, sink_feed,
 * sink_swk all use the same all-types pool with the same num_total).
 * First sink call within a step builds it; subsequent sink calls hit the
 * cached BVH+compact_xyzh, saving ~1.5s × 2 per step on sink-active steps.
 * IMPORTANT: callers MUST pass tbm=0x3f when using this cache. Mixing
 * type bitmasks against a shared cache will produce wrong answers (the
 * cached compact_xyzh / pool only includes the originally-built types).
 * Invalidated alongside the gas-only SIDX by gpu_step_sidx_invalidate(). */
gpu_spatial_index_t *gpu_step_sidx_alltypes_ptr(void);

/* Drift-time refresh: incremental bbox + BVH update without an SFC re-sort.
 * Called from run.cc after find_next_sync_point_and_drift(). Recomputes
 * each tile's bbox from current particle positions, re-fits the BVH,
 * refreshes compact_xyzh[i*4+0..2]. Tile assignments stay frozen (so
 * particles can wander into other tiles' bbox regions — inefficiency,
 * not correctness loss; each tile's pool still references its original
 * particles whose actual current positions are inside the recomputed
 * bbox). Reset to a fresh full rebuild at domain_decomp via
 * gpu_step_sidx_invalidate_full(). */
void gpu_step_sidx_invalidate(void);

/* Full invalidate: free SIDXes so the next gpu_ngb_list_build does a
 * complete rebuild including the SFC sort. Called from run.cc after
 * any domain_decomp variant, since decomp shuffles particle indices
 * and pool/tile assignments become stale. */
void gpu_step_sidx_invalidate_full(void);

/* Dirty-index API for compact_xyzh h-field tracking.
 *
 * compact_xyzh[j*4+3] must equal arena[j].KernelRadius for every j that any
 * cached neighbor search will read as a candidate (BVH-walk reads compact for
 * pruning + the leaf check_tile_particles reads compact[j*4+3] for h_j in
 * SYMMETRIC mode).  Whenever code mutates arena[j].KernelRadius (or imports
 * a ghost slot that overwrites it), it must register j as dirty so the next
 * cached gpu_ngb_list_build's compact_h_refresh covers it.
 *
 * Three primitives:
 *   _idx(i)               — single index dirty
 *   _range(start, end)    — half-open range dirty (e.g. ghost import)
 *   _indices(arr, n)      — vector of indices dirty (e.g. density h-iter sync)
 *   _all()                — full-pool dirty (fresh arena alloc, fallback)
 *
 * Internal state auto-promotes to "all dirty" when the dirty list grows past
 * a memory-budget threshold (refreshing a few million indices via list is no
 * faster than the full-pool refresh and uses more bookkeeping).
 *
 * Multi-rank guarantee: ghost imports MUST register the imported range so
 * symmetric searches that read h_j for ghost candidates see fresh values.
 *
 * Cleared by: full SIDX rebuild (gpu_spatial_index_build), refresh fire
 * inside gpu_ngb_list_build, or explicit call to _all(). */
void gpu_compact_xyzh_mark_h_dirty_idx(int i);
void gpu_compact_xyzh_mark_h_dirty_range(int start, int end);
void gpu_compact_xyzh_mark_h_dirty_indices(const int *indices, int n);
void gpu_compact_xyzh_mark_h_dirty_all(void);
/* Filter the dirty list, dropping any indices >= threshold. Used by
 * ghost_exchange_cleanup to scrub ghost-slot indices when ghost slots leave
 * scope, instead of escalating to mark_h_dirty_all (which forces a 1.5s
 * full-pool refresh on every cached build until next clear). No-op when
 * dirty_all is already set or when the list is empty. */
void gpu_compact_xyzh_dirty_drop_above(int threshold);

/* SIDX lifecycle notification hooks. ghost_exchange owns the import/cleanup
 * lifecycle; SIDX owns the spatial-index representation. These are the
 * lifecycle signals SIDX consumes (currently just bumps epochs + diagnostic
 * counters; consumed by the segmented SIDX in a later commit).
 *
 * Contract:
 *  - ghost_imported(start, count): MUST be called on every rank at the end of
 *    every ghost_exchange_*_impl path, INCLUDING the count==0 / no-receive
 *    case. count==0 invalidates any cached ghost segment from a prior import
 *    (so a ghost->no-ghost transition can never leave stale ghost data behind).
 *  - ghost_cleanup(): called from ghost_exchange_cleanup BEFORE NumPart shrinks.
 *    Frees any cached ghost segment synchronously (memory safety).
 *  - pool_changed(): called whenever the home pool's membership may have
 *    changed in ways that aren't NumPart_local (Type write, Mass<=0, particle
 *    creation/deletion, merge_split). Bumps pool_epoch so home segment is
 *    rebuilt on next NGL build. */
void gpu_sidx_notify_ghost_imported(int start, int count);
void gpu_sidx_notify_ghost_cleanup(void);
void gpu_sidx_notify_pool_changed(void);

/* Diagnostic accessors for the SIDX lifecycle counters. Used by GX_RD_CACHE /
 * DIAG_NGL prints under the verbose-diag gate. */
#ifdef __cplusplus
extern "C" {
#endif
uint64_t gpu_sidx_ghost_epoch(void);
uint64_t gpu_sidx_pool_epoch(void);
int      gpu_sidx_last_ghost_start(void);
int      gpu_sidx_last_ghost_count(void);
#ifdef __cplusplus
}
#endif

/* Backwards-compat alias for callers that haven't been updated yet (treats
 * any unknown mutation as "all dirty"). New code should use the index-aware
 * variants above. */
void gpu_compact_xyzh_mark_h_dirty(void);

/* SSOT helpers for "P[i].KernelRadius was just written" — call ONE function
 * after any KernelRadius mutation and BOTH the GPU SIDX dirty tracker AND
 * the host glt cache dirty tracker get marked. Future caches added to either
 * layer pick this up automatically.
 *
 * Order between the two underlying marks is irrelevant — neither is consumed
 * until the next gpu_ngb_list_build (GPU side) or next ghost_exchange_run
 * (host side); both happen before return. */
void gizmo_mark_kernel_radius_dirty_indices(const int *indices, int n);
void gizmo_mark_kernel_radius_dirty_range(int start, int end);

/* Build GPU-accelerated CSR neighbor list.
   If cached_idx is non-NULL and valid, reuses its tiles+BVH.
   Otherwise builds a fresh spatial index internally.
   P_shared must be accessible from GPU (SharedSpace or managed memory).
   active_indices_host: host-side array of source identifiers, size num_active.
     Default mode (source_positions_host == NULL): these are P[] indices and
     source positions are read as P[active[aa]].Pos.
     Override mode (source_positions_host != NULL): these are caller-defined
     opaque IDs; pass any sentinel (e.g. 0..num_active-1) since the kernel
     reads positions from source_positions_host instead.
   search_mode: NGB_SEARCH_ONEWAY or NGB_SEARCH_SYMMETRIC.
   type_bitmask: which particle types to include in the search pool (j-side).
   search_radius_factor: multiplier on per-source radius (default 1.0).
   search_radii_host: optional per-active-source explicit search radii
     (size num_active). NULL → use P[active[aa]].KernelRadius * search_radius_factor.
     REQUIRED when source_positions_host is non-NULL (override sources have no
     P[] entry to fall back to).
   source_positions_host: optional per-active-source position array (size
     num_active * 3, doubles, layout pos[aa*3+k] for axis k). NULL → use
     P[active[aa]].Pos (current behavior). Non-NULL → arbitrary source
     positions decoupled from any P[] index (e.g. TURB_DRIVING_SPECTRUMGRID
     grid cell centers). */
void gpu_ngb_list_build(struct particle_data *P_shared, int num_total,
                        int *active_indices_host, int num_active,
                        int search_mode, int type_bitmask,
                        gpu_neighbor_list_t *gnl,
                        gpu_spatial_index_t *cached_idx,
                        double search_radius_factor = 1.0,
                        const double *search_radii_host = NULL,
                        const double *source_positions_host = NULL,
                        const char *caller_label = "?",
                        /* j_kernel_radius_scale: multiplier on the j-side kernel
                         * radius in SYMMETRIC mode (1.0 = legacy). Set to
                         * All.TurbDynamicDiffFac for the TURB_DIFF_DYNAMIC
                         * wide-filter loops so the pair reach is the genuinely
                         * symmetric max(fac*h_i, fac*h_j). Applied at query time;
                         * cached compact_xyzh / SIDX hmax stay keyed on raw radii.
                         * See OPEN_3d_difffilter_design.md §3. */
                        double j_kernel_radius_scale = 1.0);

/* Free CSR arrays + active indices. Does NOT free tiles/BVH/pool if they
   belong to the cached spatial index (use gpu_spatial_index_free for those). */
void gpu_ngb_list_free(gpu_neighbor_list_t *gnl, gpu_spatial_index_t *cached_idx);

/* Copy gnl->neighbors (DEVICE_SPACE / CudaSpace) into a caller-allocated host
   buffer. Use when host code needs to index gnl.neighbors[] directly (e.g.
   per-source CPU loops in radfb_local, merge_split). host_dest must hold at
   least gnl->total_pairs ints; no-op when total_pairs <= 0. */
void gpu_ngb_copy_neighbors_to_host(const gpu_neighbor_list_t *gnl, int *host_dest);

/* Cross-type high-level wrapper: i-list is caller-supplied active indices of
   any type(s); j-side is filtered by j_type_bitmask. Caller supplies explicit
   per-active search radii (for loops whose kernel isn't P[i].KernelRadius —
   e.g. KernelRadiusDM, AGS_Hsml). Returns a neighbor_list_t in the mymalloc
   format used by the existing symlist API. */
struct neighbor_list_t; /* forward decl from mesh/neighbor_list.h */
void gpu_build_cross_type_neighbor_list(struct particle_data *P_host, int num_total,
                                        int *i_active_indices, int num_active,
                                        const double *i_search_radii_host,
                                        int j_type_bitmask, int search_mode,
                                        neighbor_list_t *out);

#endif /* GPU_NEIGHBOR_LIST_H */
