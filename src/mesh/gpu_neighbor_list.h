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
#include "nlr_radius_policy.h"  /* mode_b_radius_policy_t + MODE_B_RADIUS_* + SSOT wrappers */

/* GPU-resident neighbor list: CSR arrays.
   offsets: SharedSpace (host writes offsets[num_active]=total after scan).
   neighbors: DeviceSpace (GPU HBM on CUDA; never host-accessed directly).
   The list owns only these three arrays; it holds no part of the spatial index
   it was searched from. */
struct gpu_neighbor_list_t {
    int64_t *offsets;   /* [num_active+1] in SharedSpace (CSR row pointers; 64-bit) */
    int *neighbors;     /* [total_pairs] in DeviceSpace (GPU HBM — no UVM fault).
                           values are particle indices (< num_total < 2^31); only
                           the array length is 64-bit. */
    int *d_active;      /* [num_active] in SharedSpace: the active particle of each row */
    int num_active;
    int64_t total_pairs;
};


/* One segment of a spatial index: tiles + BVH over the members of one particle range (sfc_tiles.h).
   It describes its members as of its reference time ti_ref, the time it was built, and every walk
   reads it at the time of the search, so it can be KEPT while its members drift.  The pool holds the
   members' particle indices; the slot map, where a segment keeps one, is indexed from source_base. */
struct gpu_index_segment_t {
    sfc_tile_t *d_tiles = nullptr;
    tile_bvh_node_t *d_bvh = nullptr;
    int *d_pool = nullptr;
    int ntiles = 0;
    int bvh_root = -1;         /* -1: the segment has no members to walk */
    int bvh_nnodes = 0;
    /* Member rows: d_compact_xyzh[slot*SIDX_ROW_WIDTH ...] for the member in pool slot
       `slot` -- its position at ti_ref, its reach at its own clock, and its reach once
       drifted (sfc_tiles.h).  DOUBLE: float absolute positions are invalid for GIZMO's
       dynamic range (see §37/§38). */
    double *d_compact_xyzh = nullptr;
    int num_pool = 0;          /* slots: the tiles' runs of TILE_TARGET_SIZE, including the unused end of a partial tile */
    int source_base = 0;       /* the members were drawn from particles [source_base, source_base + source_count) */
    int source_count = 0;
    integertime ti_ref = 0;    /* the reference time */
    int rows_are_positions = 0;/* every member was current at ti_ref: the rows are positions, not predictions */
    int valid = 0;             /* 1 if built and usable */
    /* What keeps a segment's bounds true while it is kept: present only on a segment that is raised
       (the owned one).  The imported segment is exact at its build time and never outlives it. */
    int *d_slot_of = nullptr;  /* [source_count] the pool slot of particle source_base + o, -1 if not a member */
    int *d_level_nodes = nullptr;   /* the BVH nodes by height, leaves first */
    int *d_shear_folds = nullptr;   /* [2 per slot] the x folds (up, down) that took each member to its primary-box
                                       image, for its velocity range there; BOX_SHEARING > 1 only, else null */
    int *h_level_offsets = nullptr; /* [nlevels+1] where each height starts in d_level_nodes (host) */
    int nlevels = 0;
    double *looseness = nullptr;    /* largest tile growth per unit drift interval (SharedSpace scalar) */
    int rebuild_needed = 0;    /* a raise could not be applied: the next list build rebuilds the segment */
    int dirty_handle = -1;     /* gpu_dirty_tracker handle over the source range; -1 when not registered */
    /* What the segment was built over, beyond its range: the owned segment the owned epoch (a change of
       membership gained or a position written outside a drift, which no count can show; a membership lost is
       retired slot by slot through the dirty tracker, gpu_sidx_notify_member_lost); the imported segment
       the ghost exchange's import (a cleanup-and-reimport can land the same count with other contents). */
    uint64_t owned_epoch_when_built = 0;
    uint64_t member_loss_epoch_when_built = 0;   /* the all-types owned segment: the member-loss count it was built at */
    unsigned long long ghost_provenance_when_built = 0;
};

/* Spatial index over the rank's particles: an OWNED segment over its own particles and a GHOST
   segment over the particles imported for the current neighbour exchange.  The owned segment of the
   gas index is kept across sync points and imports, its bounds raised when a member's velocity or
   radius changes; the ghost segment is built, exact, by the first list build that needs it after an
   import and released when the imported particles are (gpu_sidx_ghost_pool_cleanup).  A walk
   visits the owned segment, then the ghost segment, of one index.  Reused across list builds with
   different search radii. */
struct gpu_spatial_index_t {
    gpu_index_segment_t owned;
    gpu_index_segment_t ghost;
    /* Type bitmask used at the most recent build. The cached rows / pool only
     * include the originally-built types; later callers MUST pass the same tbm to
     * gpu_ngb_list_build. A mismatch (e.g. a gas-only caller with tbm=1 reusing an
     * all-types cache built with tbm=0x3f) causes the walker to return neighbors of
     * the wrong types and lazy-drift to abort on unexpected particle types.
     * gpu_ngb_list_build refuses a mismatch (controlled stop, empty list) instead
     * of walking it.
     * Default -1 = "no build yet". */
    int cache_tbm = -1;
    /* Radius-policy the cached rows / tile bands were built under. Mirrors
     * cache_tbm semantics: gpu_ngb_list_build refuses the index if the caller's Spec
     * radius_policy differs from this cached value (a row's reach reflects the
     * build-time policy, so a mismatched caller would see wrong leaf-h_j and miss
     * valid pairs).  Default MODE_B_RADIUS_DEFAULT matches the policy hydro Specs
     * use; runner Mode A passes Spec::radius_policy explicitly. */
    mode_b_radius_policy_t cache_radius_policy = MODE_B_RADIUS_DEFAULT;
};


/* Release both segments of an index. */
void gpu_spatial_index_free(gpu_spatial_index_t *idx);

/* Module-level persistent SIDX for gas-only (type_bitmask=1) neighbor builds.
 * Shared across density rounds + symlist.  Its owned segment is KEPT across sync
 * points and imports: a drift needs nothing (see gpu_spatial_index_t).  Rebuilt when
 * its range or epoch no longer match, and when a fresh segment is the better search
 * (gpu_ngb_list_build). */
gpu_spatial_index_t *gpu_step_sidx_ptr(void);

/* Module-level persistent SIDX for all-types (type_bitmask=0x3f) builds.
 * Specifically for the SINK_PARTICLE codepath (sink_env1, sink_feed,
 * sink_swk all use the same all-types pool).
 * First sink call within a step builds it; subsequent sink calls hit the
 * cached BVH and rows, saving ~1.5s × 2 per step on sink-active steps.
 * IMPORTANT: callers MUST pass tbm=0x3f when using this cache. Mixing
 * type bitmasks against a shared cache will produce wrong answers (the
 * cached rows / pool only include the originally-built types).
 * Released at every sync point by gpu_step_sidx_invalidate(): its owned segment
 * is not raised, so it is kept only while nothing can move (Ti_Current advances
 * only between sync points). */
gpu_spatial_index_t *gpu_step_sidx_alltypes_ptr(void);

/* A new sync point: releases the all-types index (the gas index is kept).
 * Called from run.cc after find_next_sync_point_and_drift(). */
void gpu_step_sidx_invalidate(void);

/* The kept gas index follows its members' velocities: particles idx[0..n) just
 * had their velocity changed (host indices; non-members are skipped).  Called
 * after the kick of the active set and from gizmo_motion_bound_raise. */
void gpu_step_sidx_raise_motion(const int *idx, int n);

/* Full invalidate: free SIDXes so the next gpu_ngb_list_build does a
 * complete rebuild, re-deriving pool membership and tiles. Called from
 * run.cc after
 * any domain_decomp variant, since decomp shuffles particle indices
 * and pool/tile assignments become stale. */
void gpu_step_sidx_invalidate_full(void);

/* The imported particles are about to be released (ghost_exchange_cleanup).  The gas index releases
 * its ghost segment and keeps its owned one; the all-types index, whose owned segment is rebuilt at
 * every import (gpu_ngb_list_build), is released whole, so it is never reused past the import it was
 * built under. */
void gpu_sidx_ghost_pool_cleanup(void);

/* The owned particles changed in a way a kept index's rows cannot show and the
 * particle count does not reveal: the decomposition re-laid them out
 * (domain_particle_layout_changed), particles were rearranged (merge/split), a
 * particle became gas in place (force_tree_note_type_presence: grain promotion),
 * or a member's position was written outside a drift (the Hermite corrector).  Bumps the owned epoch, so an index
 * built before it is rebuilt on its next use.  A member that loses its mass is not
 * announced: the pair kernel owns the Mass > 0 test, and
 * rearrange_particle_sequence, which removes the slot, announces it then.
 * Ghost import and cleanup need no call: the index keys its imported particles on
 * the ghost exchange's own record (ghost_pool_is_live, ghost_provenance_epoch). */
void gpu_sidx_notify_owned_changed(void);
/* Particle `particle` changed type in place to a type other than gas (star or sink formation from gas; a
 * star becoming a sink is announced the same way and costs only a no-op slot lookup).
 * The gas index keeps the segment: the particle is marked like a changed radius, and the next maintenance
 * retires its slot, so no walk returns it, while the bounds it widened stay conservative until a rebuild.
 * The all-types index, whose per-type reach bands it no longer matches, is rebuilt on its next use. */
void gpu_sidx_notify_member_lost(int particle);
/* The same layout change, made by a full decomposition that also ordered the particles along the
 * space-filling curve: a gas segment built over them at this time, before anything changes them, takes
 * that order instead of sorting. */
void gpu_sidx_notify_owned_reordered(void);

#ifdef __cplusplus
extern "C" {
#endif
#ifdef __cplusplus
}
#endif

/* "P[i].KernelRadius (or another input of its reach) was just written": call after any
 * such write.  A kept owned segment holds each member's reach in its row and the largest
 * reach below each tile and node; the next gpu_ngb_list_build over it rewrites the marked
 * rows and raises the bands above them before it walks (a mark past a few million
 * particles rebuilds the segment instead).  The imported segment needs no marks: an
 * imported particle's radius does not change while it is imported. */
void gizmo_mark_kernel_radius_dirty_indices(const int *indices, int n);
void gizmo_mark_kernel_radius_dirty_range(int start, int end);

/* Build GPU-accelerated CSR neighbor list.
   If cached_idx is non-NULL and valid, reuses its tiles+BVH.
   Otherwise builds a fresh spatial index internally.
   P_shared must be accessible from GPU (SharedSpace or managed memory).
   active_indices_host: host-side array of source identifiers, size num_active.
     Default mode (source_positions_host == NULL): these are P[] indices, and each
     source's position is the particle's own Pos, which must be at the time of the
     search (a source that is not stops the run, controlled stop 7743) -- never the
     index's rows, which describe members at the index's reference time.
     Override mode (source_positions_host != NULL): these are caller-defined
     opaque IDs; pass any sentinel (e.g. 0..num_active-1) since the kernel
     reads positions from source_positions_host instead.
   search_mode: NGB_SEARCH_ONEWAY or NGB_SEARCH_SYMMETRIC.
   type_bitmask: which particle types to include in the search pool (j-side).
   search_radius_factor: multiplier on per-source radius (default 1.0).
   search_radii_host: optional per-active-source explicit search radii
     (size num_active). NULL → the source particle's own reach under
     radius_policy, times search_radius_factor. REQUIRED when
     source_positions_host is non-NULL (override sources have no P[] entry to
     fall back to).
   source_positions_host: optional per-active-source position array (size
     num_active * 3, doubles, layout pos[aa*3+k] for axis k); arbitrary source
     positions decoupled from any P[] index (e.g. TURB_DRIVING_SPECTRUMGRID grid
     cell centers).
   The list is EXACT for the kept pool: each row holds the members within R_i of
   the source (ONEWAY) or within max(R_i, j_kernel_radius_scale * r_j) (SYMMETRIC),
   with every member drifted to the time of the search.  A member whose mass went
   to zero stays until rearrange_particle_sequence removes it: the pair kernel owns
   the Mass > 0 test. */
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
                         * cached rows / SIDX bands stay keyed on raw radii. */
                        double j_kernel_radius_scale = 1.0,
                        /* radius_policy: Spec::radius_policy from the runner;
                         * controls per-particle reach used to (re)populate
                         * the rows' reaches; a cached index built under another
                         * policy is refused (controlled stop, empty list).  Default
                         * MODE_B_RADIUS_DEFAULT (= GAS_KERNEL) matches the
                         * policy hydro Specs use, so non-runner callers
                         * (merge_split / radfb_local / twopoint / turb_powerspectra
                         * / legacy symlist) sharing the gas-only step-persistent
                         * SIDX are not refused.
                         * Multi-type non-runner callers (twopoint with 0xFF)
                         * pass cached_idx=NULL → build local SIDX → policy
                         * doesn't matter for caching. */
                        mode_b_radius_policy_t radius_policy = MODE_B_RADIUS_DEFAULT);

/* Free the list's three arrays and clear them.  Safe on a list that was never
   built, or built empty.  The spatial index is not the list's to free. */
void gpu_ngb_list_free(gpu_neighbor_list_t *gnl);

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

/* Which local supply-pool slots this rank sends to which peer (owned by
   ghost_exchange.cc).  A receiver walk hands it the pool slots it accepts:
   peers in ascending order, each peer's slots in one or more consecutive calls,
   repeats allowed -- the set keeps each (peer, slot) once.  Returns 0, or nonzero
   when the set could not take the slots; it is then failed, and the caller stops
   emitting and reports the failure rather than trying another backend. */
struct ghost_send_set;
int gx_send_set_emit(struct ghost_send_set *send_set, int peer, const int *pool_slots, int n);

/* One received envelope's candidates: the local particles (P[] indices) its walk
   found may be neighbours once drifted to the current time. */
struct gx_export_envelope_t;
struct gx_candidate_row {
    const struct gx_export_envelope_t *envelope;
    int  peer;           /* the task that sent the envelope */
    int *local_index;    /* candidates; overwritten in place with accepted pool slots */
    int  count;
};

/* The receiver backends' shared accept: drifts the distinct candidates of all rows
   to the current time, keeps exactly those each envelope's query accepts at their
   current position and radius, and hands their pool slots to the send set.  Rows
   must be in ascending peer order.  Returns 0, or nonzero with the send set failed
   (nothing emitted if the drift or its storage failed); the caller then reports the
   failure rather than trying another backend. */
int gx_send_set_accept_rows(struct ghost_send_set *send_set, struct gx_candidate_row *rows, long n_rows,
                            int search_mode, mode_b_radius_policy_t radius_policy,
                            double j_radius_scale, double safety_factor);

/* What a receiver backend did with the envelopes it was handed. */
enum {
    GX_RECEIVER_COMPLETED = 0,   /* every accepted pair is in the send set */
    GX_RECEIVER_DECLINED  = 1,   /* nothing emitted; the host walk must answer */
    GX_RECEIVER_FAILED    = 2    /* stopped after emitting; the send set is unusable */
};

/* Device traversal of received export envelopes (the supply-rank half of
   request-driven ghost discovery).  Resumes a bounded subtree walk from each
   envelope's start nodes, records at every local leaf the particles that may be
   neighbours once drifted, and hands those candidates to gx_send_set_accept_rows --
   the same accept the host walk's candidates go through, so the two backends are
   interchangeable.

   Returns GX_RECEIVER_DECLINED when it declined, in which case the caller must run
   the host walk instead.  Every decline is decided before anything is emitted, so
   the two backends never interleave in the send set.  Declining is rank-local and
   safe: the window this runs in contains no collectives.  Reasons to decline are a
   search mode not yet supported here, too little work to be worth staging the
   local leaves, the node geometry not being certified current on the device, a
   tree mirror that does not cover the range the walk may reach, and allocation
   failure.

   Returns GX_RECEIVER_FAILED when it stopped after it had begun emitting: the send
   set could not grow, or one of two states that also request a controlled stop,
   because neither can arise unless the tree or the traversal is already wrong and
   continuing would answer with a silently truncated or duplicated neighbour set --
   an index in the gap between the particle slots and the node base (what the host
   walk stops on), and a single envelope admitting more particles than the rank
   owns.  A host rerun cannot repair a half-filled set, so the caller fails too.

   envelope_peer[k] is the task that sent envelope k, non-decreasing in k.
   j_to_pool / npart_bound map a local particle index to its slot in the supply
   pool (negative = not in the pool).  radius_policy, j_radius_scale and
   safety_factor are the caller's symmetric-search reach, passed unchanged to the
   shared accept; ONEWAY ignores them. */
int gx_device_receiver_walk(const struct gx_export_envelope_t *envelopes, long n_env,
                            const int *envelope_peer,
                            unsigned int supply_mask, int search_mode,
                            mode_b_radius_policy_t radius_policy,
                            double j_radius_scale, double safety_factor,
                            const int *j_to_pool, int npart_bound,
                            int num_pool, struct ghost_send_set *send_set);

/* Describe this rank's tree to the device walk, or say why it cannot be
   described.  Returns 0 with *out filled, 1 with *out untouched and a one-shot
   line naming `caller`.  Every user of the device traversal asks the same two
   questions -- is there a mirror, does it cover what a walk can reach -- and
   derives the same three index-class boundaries, so that is written once here.

   `local_particle_slots` is the caller's own choice and the walk's ownership
   line: the owned count for a walk answering for this rank alone, the full slot
   count for one meant to see imported ghosts.  See mesh/device_tree_walk.h.

   This does NOT establish that the geometry is CURRENT.  Who is expected to
   have drifted the nodes, and what to do when nobody has, differs by caller;
   each settles it at its own site. */
struct GxDeviceTreeView;
int gx_device_tree_view_build(struct GxDeviceTreeView *out, int local_particle_slots,
                              const char *caller);

/* Put this rank into the state a FUSED device walk needs -- one that evaluates a
   pair kernel where it lands, so it reads particle fields and cannot drift a
   stale one when it arrives -- and describe its tree.  Returns 0 with *out
   filled and the rank current, or 1 with the walk declined and the host to
   answer.  Call once, before any discovery round: after it, nothing on the rank
   goes stale again within the call.

   Sweeps the node geometry when nothing else has, holds the touched-set
   workspace the discovery passes record into, sets the walk's ownership line to
   the owned count, and refuses a tree with no built nodes -- which a walk from
   the root would otherwise enter as an unbounded read rather than as an empty
   answer. It does NOT drift the particles: which ones the walk needs is decided
   by discovering them, per pass, at the evaluation sites. */
int gx_device_fused_walk_prepare(struct GxDeviceTreeView *out, const char *caller);


/* The touched-set workspace a fused walk records into.  See GxTouchedSet.
 *
 * Allocation is a CAPABILITY, not a repair: it is arranged in the preparation
 * above, where a rank that cannot have it declines and the collective readiness
 * vote pulls every rank to the host path together.  Doing it later, inside a
 * discovery round, would let one rank answer differently from its peers in the
 * middle of an exchange.
 *
 * `begin_call` opens a generation and is called once per fused call, not once
 * per pass: a particle drifted for the first pass is still current for the
 * later ones, so a generation per call is what stops the round loop re-examining
 * it.  `view` hands back a copy for a kernel to capture by value.
 *
 * `drift_and_mark` closes one pass: it waits for the recording kernel, advances
 * the particles it recorded, marks their kernel radii dirty, and reopens the
 * cursor for the next pass.  The three belong together because the fence is the
 * step that is invisible when it is missing.  It returns nothing deliberately --
 * the drift reports only whether a controlled stop is already pending, which is
 * a property of the run rather than of these particles, and a caller that
 * branched on it would take a different path through a collective exchange than
 * its peers.
 *
 * `begin_epoch_owned` is the general form and `begin_call` is the fused walk's wrapper for it.
 * It REFUSES (nonzero, recorder untouched) while claims are still outstanding, so a second
 * caller can never advance the generation over another's live stamps; a refused caller takes
 * its own safe route.
 * `retire` closes an epoch so that any later claim reports itself as out of phase.  The claim
 * itself is device-callable and lives in declarations/gpu_recorder_claim.h, because its body
 * needs Kokkos and this header is read by host-only units. */
int  gx_touched_set_ensure(int local_particle_slots);
int  gx_touched_set_begin_call(void);
int  gx_touched_set_begin_epoch_owned(int owner);
void gx_touched_set_retire(int owner);
long long gx_touched_set_refused_epochs(void);
long long gx_touched_set_retire_faults(void);
struct GxTouchedSet gx_touched_set_view(void);
void gx_touched_set_drift_and_mark(integertime time1);
/* Released once at shutdown, before Kokkos is finalized. Not on the tree or
   domain epoch: the storage is persistent by design. */
void gx_touched_set_release(void);

/* Raise the motion bound of particles idx[0..n) in the gravity tree's nodes,
 * with the top-level part carried to the other ranks at the next tree-update
 * phase.  Called by whoever changed a particle's velocity outside the kick,
 * with the list of particles it wrote; the cost is the list, never the rank. */
void gizmo_motion_bound_raise(const int *idx, int n);

/* The motion-target set (GxMotionTargetSet).  The runner opens a generation
 * per call for a loop that writes neighbour velocities, the kernels mark into
 * it, and `consume` raises the marked set and closes the generation.  The
 * host mark serves units without Kokkos (the reverse writeback landing on the
 * owner; a host module's own loop); the device mark is below.  `armed` is set
 * by the runner around a loop's reverse writeback so the apply loop knows to
 * mark the deltas' targets. */
int  gx_motion_target_ensure(int local_particle_slots);
void gx_motion_target_begin_call(void);
struct GxMotionTargetSet gx_motion_target_view(void);
void gx_motion_target_mark_host(int j);
void gx_motion_target_consume(void);
void gx_motion_target_set_armed(int armed);
int  gx_motion_target_armed(void);
void gx_motion_target_release(void);

#if defined(KOKKOS_VERSION)   /* device atomics: only a unit that carries Kokkos can compile this */
/* Mark owned particle j as a motion target for this call.  Exactly one marker
 * appends it, however many pairs reach it; an index outside the owned range
 * (a ghost copy) is not this rank's to raise and is ignored. */
KOKKOS_INLINE_FUNCTION
void gx_motion_target_mark(const struct GxMotionTargetSet &ts, int j)
{
    if(!ts.seen || j < 0 || j >= ts.capacity) {return;}
    if(Kokkos::atomic_exchange(&ts.seen[j], ts.gen) == ts.gen) {return;}
    const int slot = Kokkos::atomic_fetch_add(ts.counter, 1);
    if(slot < ts.capacity) {ts.list[slot] = j;}
}
#endif

#endif /* GPU_NEIGHBOR_LIST_H */
