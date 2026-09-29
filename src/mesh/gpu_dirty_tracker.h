/* gpu_dirty_tracker.h — per-cache h-dirty bitset registry.
 *
 * Replaces the global g_dirty_list/g_dirty_all pair in gpu_neighbor_list.cc
 * with per-cache state, so multiple persistent SIDX caches (gas-only,
 * alltypes, future home/ghost segments) can each track their own pending
 * h refresh state.
 *
 * Design (locked in plan v5):
 *  - Each cache registers a dense range [base, base+count) (segment-keyed
 *    by particle index, as the marks are).
 *    Bitset sized count bits, indexed by (j - base).
 *  - mark_indices(j[], n) walks every registered cache, range-routes per-j,
 *    sets bits in ALL caches whose range covers j. Same for mark_range.
 *  - mark_all_global() sets all_dirty on every registered cache (preserves
 *    the global-scope semantics of the old gpu_compact_xyzh_mark_h_dirty_all
 *    for unknown-scope mutations).
 *  - consume(handle, callback) iterates set bits via __builtin_ctzll,
 *    invokes callback(j) for each, clears that cache's bitset + all_dirty
 *    flag. Other caches' state is untouched.
 *  - popcount(handle) returns the unique-index count for promote-threshold
 *    decisions.
 *
 * Threading: marks fire from single-threaded host paths (density iter
 * post-loop, lazy drift post-walk, ghost exchange post-import, etc.).
 * No cross-thread concurrency on the tracker itself.
 *
 * Multi-rank: tracker state is per-rank. No MPI involved. */

#ifndef GPU_DIRTY_TRACKER_H
#define GPU_DIRTY_TRACKER_H

#include <stdint.h>

typedef int gpu_dirty_handle_t; /* opaque; -1 == invalid */

/* Register a dense particle-index range. Allocates a bitset sized to count
 * bits (rounded up to 64-bit word boundary). Returns a handle, or -1 on
 * out-of-slots. The caller stores the handle and unregisters it on cache
 * free.
 *
 * start_clean tells the tracker whether the cache's h values are already
 * current. Pass 1 when the caller has just written every row from the live
 * P[] (nothing can have changed since), so the first consume does no work.
 * Pass 0 when the backing store is uninitialised or of unknown age, which
 * makes the first consume refresh the whole range. */
gpu_dirty_handle_t gpu_dirty_tracker_register(int base, int count, int start_clean);

/* Free the bitset and free the slot for reuse. Safe to call with -1. */
void gpu_dirty_tracker_unregister(gpu_dirty_handle_t handle);

/* Mark a list of particle indices dirty across all registered caches.
 * Each j is range-checked against each cache; bit set in every cache that
 * covers it. Promote-to-all fires per cache when popcount > threshold.
 * O(n * n_caches) — n_caches is typically 2 (gas + alltypes).  */
void gpu_dirty_tracker_mark_indices(const int *indices, int n);

/* Mark a contiguous range [start, end). Faster than mark_indices when
 * marking a large region (e.g. full-drift over [0, NumPart)) because
 * each cache's bitset can be set in O((end-start)/64) word writes. */
void gpu_dirty_tracker_mark_range(int start, int end);

/* Mark every registered cache's all_dirty flag (e.g. for unknown-scope
 * mutations). Preserves the global semantics of the old
 * gpu_compact_xyzh_mark_h_dirty_all. */
void gpu_dirty_tracker_mark_all_global(void);

/* Iterate set bits for one cache's bitset, calling callback(j) for each.
 * If the cache's all_dirty flag is set, callback fires for every j in
 * [base, base+count). Either way, the bitset and all_dirty are cleared at
 * the end. */
void gpu_dirty_tracker_consume(gpu_dirty_handle_t handle,
                               void (*callback)(int j, void *userdata),
                               void *userdata);

/* Number of unique dirty indices for one cache. Returns count if
 * all_dirty is set; otherwise the popcount of the bitset. Used by callers
 * to decide whether to take a fast all-refresh path. */
int gpu_dirty_tracker_popcount(gpu_dirty_handle_t handle);

/* Reports whether a cache's all_dirty flag is set (without consuming).
 * Useful for diagnostic prints. */
int gpu_dirty_tracker_is_all_dirty(gpu_dirty_handle_t handle);

#endif /* GPU_DIRTY_TRACKER_H */
