/* ghost_writeback.h — reverse communication of ghost particle modifications.
 *
 * After a GPU neighbor loop writes to ghost j-particles (via Kokkos atomics),
 * ghost_writeback communicates those modifications back to the home MPI ranks.
 *
 * Runner-template loops reverse-communicate j-particle writes through the
 * generic bundle scaffold below (ghost_writeback_begin_bundle /
 * ghost_writeback_end_bundle, with manifest ops from ghost_writeback_ops.h).
 * Usage pattern for each loop that modifies j-particles:
 *   1. ghost_writeback_begin_bundle()  — snapshot ghost-region values
 *   2. [run GPU kernel with atomic j-writes]
 *   3. ghost_writeback_end_bundle()    — diff, pack, reverse MPI, apply
 * Each Spec declares its bundle via a GHOST_WRITEBACK_*(field) manifest.
 *
 * The ghost provenance map (home rank + index per ghost) is built during
 * ghost_exchange() and freed in ghost_exchange_cleanup(). Accessors are
 * in ghost_exchange.cc.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef GHOST_WRITEBACK_H
#define GHOST_WRITEBACK_H

#include <stddef.h>

/* ============================================================================
 * Generic ghost-writeback scaffold (Pass B.iv + B.iv.1 bundle-level exchange).
 *
 * One bundle owns N callbacks; each callback is one (op, field) pair.
 * The scaffold owns: ghost-count + home-rank/home-index lookups, snapshot
 * dispatch, single count pass, single pack pass, ONE Alltoall (count
 * matrix) + ONE gizmo_mpi_alltoallv_typed (data bytes), single apply
 * pass, exactly-once detector-register, exactly-once arena_invalidate,
 * cleanup dispatch.
 *
 * MPI-call invariant: total MPI calls per non-trivial end_bundle is O(1)
 * regardless of the number of callbacks. All callbacks in one bundle
 * share one count exchange and one byte exchange; the receiver
 * demultiplexes by canonical bundle order. This makes it safe to put
 * dozens-to-hundreds of writeback fields in one Spec's manifest (mechFB,
 * chemistry, etc.) without ballooning MPI cost.
 *
 * Strict no-op semantics: empty bundle (n_callbacks==0) returns from
 * begin_bundle and end_bundle at line 1, before any snapshot / register /
 * invalidate / MPI / cleanup call. Default builds without ghost-writeback
 * physics flags behave identically to "writeback absent."
 *
 * Legacy invalidation timing preserved exactly:
 *   begin: gpu_particles_arena_invalidate() only when num_ghosts > 0
 *   end:   gpu_particles_arena_invalidate() only when NTask > 1
 *   end:   register_writeback always (when bundle non-empty)
 *
 * Mode B paths (local + remote): bundle hooks are NOT called by the
 * neighbor-loop runner — j-side writes happen on the rank that owns j
 * (Mode B local: same rank as walker; Mode B remote: peer rank running
 * the walker on its own pool), so no reverse-comm is needed. The runner's
 * nlr_path_uses_imported_ghosts predicate gates begin_bundle/end_bundle
 * to Mode A only.
 *
 * Physics-author surface: see mesh/ghost_writeback_ops.h for the manifest
 * macros (GHOST_WRITEBACK_BUNDLE_BEGIN, GHOST_WRITEBACK_PARTICLE_MIN, ...)
 * that hide the callback / context plumbing.
 * ========================================================================== */

struct ghost_writeback_callback {
    /* Fixed transport-element size in bytes; must match what pack_delta
     * writes and apply_delta reads. */
    size_t delta_size;

    /* Pre-kernel snapshot. Called once per begin_bundle on each rank when
     * num_ghosts > 0. The callback is responsible for allocating any
     * dynamic state inside *ctx. */
    void (*snapshot)(void *ctx, int num_ghosts, int num_local);

    /* Returns 1 if ghost g (g in [0, num_ghosts)) needs reverse-comm; 0
     * otherwise. Hot path: called twice per ghost (count + pack passes). */
    int  (*delta_for_ghost)(void *ctx, int g, int num_local);

    /* Write delta_size bytes describing ghost g's change to *out_delta. */
    void (*pack_delta)(void *ctx, int g, int num_local, void *out_delta);

    /* Apply received delta to home particle. The delta encodes home_index
     * (callback's contract), and it is the record's FIRST field: the apply
     * loop reads it there to mark the owner as a motion target when the loop
     * writes neighbour velocities. */
    void (*apply_delta)(void *ctx, const void *in_delta);

    /* Free per-bundle dynamic state. Called once during end_bundle on every
     * rank (including NTask<=1 short-circuit), when bundle is non-empty.
     * Must be idempotent (legacy callbacks no-op when their snap was never
     * allocated, e.g. num_ghosts==0). */
    void (*cleanup)(void *ctx);

    /* Per-callback opaque state, stable across begin/end. Callback owns
     * any pointers stored inside *ctx. */
    void *ctx;
};

struct ghost_writeback_bundle {
    /* Array of callback pointers; bundle assembly macros may include
     * trailing nullptr entries which the count must exclude. */
    const struct ghost_writeback_callback *const *callbacks;
    int                                            n_callbacks;
    /* Optional loop name for the rank-0 traffic print in _end_bundle.
     * Set automatically by GHOST_WRITEBACK_BUNDLE_BEGIN(LOOP) to #LOOP.
     * nullptr for manually-constructed bundles → silent (backwards-compatible). */
    const char                                    *loop_name;
};

/* Scaffold entry points. Begin must be paired with end on the same bundle.
 *
 *   begin_bundle: per-callback snapshot (when num_ghosts>0); exactly-once
 *                 arena_invalidate (when num_ghosts>0). Strict no-op when
 *                 n_callbacks==0.
 *   end_bundle:   exactly-once register_writeback; per-callback count +
 *                 pack + alltoallv + apply (NTask>1) or skipped (NTask<=1);
 *                 per-callback cleanup; exactly-once arena_invalidate
 *                 (NTask>1 only). Strict no-op when n_callbacks==0.
 */
void ghost_writeback_begin_bundle(const struct ghost_writeback_bundle *bundle);
void ghost_writeback_end_bundle  (const struct ghost_writeback_bundle *bundle);

/* Per-physics compatibility wrappers retired:
 *   - ghost_writeback_{zero_,}hydro      — hydro corridor commit 8b retired
 *     these; hydro_force_loop reverse-communicates via the bundle scaffold.
 *   - ghost_writeback_{zero_,}swallowtime — audit J4 retired these;
 *     SinkEnv1Spec in sinks/sink_env1_loop.cc reverse-communicates SwallowTime
 *     via the bundle scaffold (PARTICLE_MIN on particle_data::SwallowTime). */


/* sinkfeed retired: replaced by the bundle scaffold in
 * sinks/sink_feed_loop.cc (PARTICLE_MAX(SwallowID) -- a claim token, not an ID, so the max picks
 * one claimant even when two sinks share an identifier +
 * GAS_ADD(Injected_Sink_Energy) ops). */

/* MechFB variant retired: the ghost_writeback_mechfb compatibility wrapper
 * had no remaining callers — MechFBSpec in galaxy_sf/mechfb_loop.cc owns the
 * per-gas MechFBGasDelta reverse-comm via its mechfb_writeback_detail custom
 * generic-bundle callback. */

/* Grain backreaction retired: the ghost_writeback_{zero_,}grainbackrx
 * compatibility wrappers had no remaining callers — GrainBackrxSpec in
 * solids/grain_physics_loop.cc reverse-communicates the grain_backrx
 * j-writes via the generic ghost-writeback bundle. The multifluid intent
 * to drop the GRAIN_FLUID half of the gate (so GRAIN_BACKREACTION compiles
 * under HYDRO_MULTIFLUID-only builds) is captured by the wider
 * DO_FLUID_ALTSPECIES_DRAG_CALCULATION outer guard on grain_physics_loop.{cc,h}. */

/* sinkswallow retired: replaced by the bundle
 * scaffold's compound callback in sinks/sink_swk_loop.cc (snapshot +
 * CR/B-aware predicate + clamped apply). */

/* RadFBRP runner-template port: the
 * legacy ghost_writeback_zero_radfbrp + ghost_writeback_radfbrp snapshot-
 * based reverse-comm path is retired; RadFBRPSpec in
 * galaxy_sf/radfb_rp_loop.{h,cc} now uses the generic ghost-writeback
 * bundle (PARTICLE_ADD_VEC3 on Vel/dp + GAS_ADD_VEC3 on VelPred). */

/* RT source injection variant retired. RtSrcInjectionSpec in
 * radiation/rt_source_injection_loop.{h,cc} now owns the j-side
 * reverse-comm via the generic ghost-writeback bundle. */

/* --- Ghost-write detector ---------------------------------------------------
 * Debug-build instrumentation that catches "GPU kernel touched ghost particles
 * but no ghost_writeback_* was called". Active only under GIZMO_GPU_ARENA_DEBUG;
 * compiles to no-ops otherwise. Pattern:
 *
 *   ghost_write_detector_begin("kernel_name");
 *   ... kernel + memcpy back to host ...
 *   ghost_writeback_foo();   // increments writeback counter
 *   ghost_write_detector_end();
 *
 * If end() sees ghost-byte differences AND no writeback ran since begin(), it
 * aborts with a diagnostic naming the kernel. Each writeback function calls
 * ghost_write_detector_register_writeback() at entry to mark itself accounted.
 * Catches the "forgot to call ghost_writeback_X entirely" failure mode that
 * silently dropped multi-rank RT source injection and HII heating physics. */
#ifdef GIZMO_GPU_ARENA_DEBUG
void ghost_write_detector_begin(const char *kernel_name);
void ghost_write_detector_end(void);
void ghost_write_detector_register_writeback(void);
/* Re-snapshot ghost state at the current point. Used by gpu_ngb_list_build's
 * lazy-drift loop to move the detector baseline past the per-ghost drift
 * (Ti_current + Pos updates), so the kernel-window comparison flags only
 * genuine kernel-side writes that need a writeback. */
void ghost_write_detector_resnapshot_after_lazy_drift(void);
#else
static inline void ghost_write_detector_begin(const char *k) { (void)k; }
static inline void ghost_write_detector_end(void) {}
static inline void ghost_write_detector_register_writeback(void) {}
static inline void ghost_write_detector_resnapshot_after_lazy_drift(void) {}
#endif

/* Accessors from ghost_exchange.cc */
int ghost_get_num_ghosts(void);
int ghost_get_num_local(void);
int *ghost_get_home_rank(void);
int *ghost_get_home_index(void);
int *ghost_get_wb_recv_count(void);
int *ghost_get_wb_recv_disp(void);
int *ghost_get_wb_send_count(void);
int *ghost_get_wb_send_disp(void);

#endif /* GHOST_WRITEBACK_H */
