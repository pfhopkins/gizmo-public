/*! \file ghost_exchange.cc
 *  \brief Ghost particle exchange for GPU-ready neighbor finding.
 *
 *  Replaces the pseudo-particle export mechanism with an upfront "import-the-neighbors"
 *  pattern: before any neighbor loop, exchange boundary particles between MPI ranks so
 *  that all neighbors are local. Subsequent kernels iterate over local + ghost particles
 *  without secondary MPI phases.
 *
 *  Ghost particles are appended to P[] and CellP[] arrays at indices >= NumPart_before_ghost.
 *  After all neighbor operations complete, ghost_exchange_cleanup() resets NumPart/N_gas.
 *
 *  This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *
 *  KNOWN LIMITATIONS (to be optimized):
 *  - Sends full P[i]/CellP[i] structs per ghost (~2-5 KB/particle). Should use compact
 *    struct with only fields needed by the active kernel set (~200 bytes).
 *  - O(NTopleaves^2) overlap check between local and remote leaves. Should use spatial
 *    sorting or tree-based pruning for simulations with many top-level leaves.
 *  - Global MPI_Allreduce on need_leaf (NTopleaves ints). Could use point-to-point for
 *    sparse communication patterns.
 *  - Per-task routing is approximate: uses MPI_Allreduce(MPI_MAX) on need_leaf, so a leaf
 *    requested by ANY task is sent to ALL requesting tasks. Per-task leaf request lists
 *    would reduce traffic.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include "../declarations/allvars.h"
#include "../declarations/lifecycle_counters.h"
#include "../core/proto.h"
#include "../system/mpi_alltoallv_typed.h"
#include "gpu_neighbor_list.h" /* gpu_compact_xyzh_mark_h_dirty_range */
#include "sfc_tiles.h"           /* build_sfc_tiles, build_tile_bvh, sfc_tile_t, tile_bvh_node_t */
#include "neighbor_list.h"       /* NGB_SEARCH_ONEWAY, NGB_SEARCH_SYMMETRIC */
#include "ghost_exchange_functions.h" /* gx_pair_accept_wrap_and_test: shared accept, wraps via the canonical macros */
#include "ghost_writeback.h"     /* ghost_get_num_local (bounded fine-tree walk) */
#include "ghost_exchange_spec.h"
#include "mode_b_local_walker.h"
#ifdef _OPENMP
#include <omp.h>                 /* threaded sender export + receiver walk below */
#endif

/*
 * ============================================================================
 * COMPACT GHOST STRUCT FIELD REQUIREMENTS (for future optimization)
 *
 * Currently sends full P[i] + CellP[i] structs. When optimizing, the minimum
 * fields needed per kernel are:
 *
 * ALL kernels need from P: Pos[3], Mass, Type, KernelRadius, NumNgb, TimeBin
 *
 * density_evaluate additionally needs:
 *   P: Vel[3]
 *   CellP: VelPred[3], InternalEnergyPred, Density
 *   + #ifdef MAGNETIC: CellP.BPred[3]
 *   + #ifdef COSMIC_RAY_FLUID: CellP.CosmicRayEnergyPred[N_CR_PARTICLE_BINS]
 *   + #ifdef RADTRANSFER: CellP.Rad_E_gamma[N_RT_FREQ_BINS], Rad_E_gamma_Pred
 *
 * hydro_gradient_calc additionally needs:
 *   All of density fields, plus:
 *   CellP: Pressure, MaxSignalVel
 *
 * hydro_force_evaluate additionally needs:
 *   All of gradient fields, plus:
 *   CellP: DtInternalEnergy, HydroAccel[3], SoundSpeed, Gradients
 *   P: GravAccel[3]
 *   + #ifdef DIVBCLEANING_DEDNER: CellP.PhiPred
 *   + all MHD/RT/CR evolved quantities
 * ============================================================================
 */

/* saved state for cleanup */
static int NumPart_before_ghost = -1;
static int N_gas_before_ghost = -1;
static int NumGhostParticles = 0;
/* The largest ghost import this rank has completed, over the epoch running now and the one before
 * it.  The capacity a rank needs is set by its worst import, not its most recent one, and ghost
 * demand is strongly uneven between ranks, so this is kept per rank and is the measured term the
 * epoch sizing asks for.  Two epochs rather than one because a single quiet epoch would otherwise
 * be enough to justify releasing storage that the next one immediately asks for again, and paying
 * two migrations to save memory for one epoch is a poor trade.  Neither value is carried across a
 * restart: a restored capacity is already whatever the run had grown to, and the sizing refuses to
 * lower it until it has observed an import. */
static int GhostEpochHighWater = 0;
static int GhostPreviousEpochHighWater = 0;

/* Ghost provenance map: for each ghost particle, the home MPI rank and index.
   Used by ghost_writeback to reverse-communicate j-particle modifications.
   Allocated with malloc (not mymalloc) to avoid stack ordering issues. */
static int *ghost_home_rank_map = NULL;     /* [NumGhostParticles] home MPI rank */
static int *ghost_home_index_map = NULL;    /* [NumGhostParticles] home P[]/CellP[] index */
static int *ghost_wb_recv_count = NULL;     /* [NTask] ghosts received from each rank */
static int *ghost_wb_recv_disp = NULL;      /* [NTask] displacement by source rank */
static int *ghost_wb_send_count = NULL;     /* [NTask] ghosts we sent to each rank */
static int *ghost_wb_send_disp = NULL;      /* [NTask] displacement for what each rank got from us */

/* Send-side provenance for ghost_refresh_values(): the ordered list of LOCAL
   indices this rank exported at the last import (grouped by ghost_wb_send_*),
   so a value-only refresh can re-pack current owner P/CellP without re-running
   discovery. Preserved unconditionally at import (ownership taken from the
   per-import send_home_idx buffer), freed in ghost_exchange_cleanup(). The actual
   refresh guard is (non-NULL ghost_send_home_idx) + (send/recv totals match the
   live pool): a non-NULL pointer implies "an import happened and no cleanup since"
   (cleanup NULLs it). The monotonic epoch below counts completed imports; it has
   two consumers: (a) the refresh diagnostic harness asserts a value-refresh does
   NO reimport (epoch unchanged) while a full cleanup+reimport bumps it; (b) the
   hydro corridor records the epoch its published CSR was built from and permits
   the value-refresh fast path ONLY on an epoch match — a live pool from some
   OTHER import could pass the count checks by coincidence while the CSR still
   indexes the old slot layout. */
static int *ghost_send_home_idx = NULL;     /* [ghost_send_home_count] exported local indices, send order */
static int  ghost_send_home_count = 0;      /* == total_send at last import */
static unsigned long long g_ghost_provenance_epoch = 0; /* import counter (see above) */

/* Persistent supply-pool cache for the request-driven ghost exchange.  Within
 * a step the local pool [0..NumPart_local) is stable across the several
 * ghost_exchange calls a step makes, so deriving membership once and reusing it
 * saves N-1 scans of P[] per step.
 *
 * Membership only: which particles are eligible supply, and the reverse map
 * from a particle to its pool slot.  Nothing position-dependent is cached --
 * the routed producer reads live positions -- so the entry stays valid as the
 * particles move, and is keyed on the supply-identity epoch that every event
 * changing pool membership already bumps.  NumPart and the eligible type mask
 * are checked at use time as a defensive cross-check.
 *
 * Allocated with plain malloc/free: the entry outlives the function frame, which
 * would violate the mymalloc stack's LIFO ordering. */

/* Membership/order epoch for the supply pool.  Rank-local and compared only for
 * equality: it answers "is the pool I cached still the same set, in the same
 * order?", nothing more.  Bumped ONLY where particles are created, eliminated,
 * or moved between slots (rearrange_particle_sequence); NEVER by drift, h, or a
 * radius policy, none of which can change membership.  64-bit so wraparound is
 * not a case anyone has to reason about. */
static long long g_supply_identity_epoch = 0;

extern "C" void ghost_exchange_supply_identity_changed(const char *reason)
{
    (void)reason;   /* named at the call site so the reason is greppable there */
    g_supply_identity_epoch++;
}

struct ghost_local_tree_cache_t {
    int valid;
    int NumPart_when_built;
    /* Identity generation this entry's pool/j_to_pool were built against. */
    long long identity_epoch_when_built;
    unsigned int eligible_type_mask_when_built;
    int num_pool;
    int *pool;                     /* [num_pool] malloc */
    int *j_to_pool;                /* [NumPart_when_built] malloc, j -> pool_pos or -1 */
    unsigned char *mark;           /* [num_pool] calloc, the send set's per-peer dedup marker;
                                      all zero between calls */
};
/* Membership only: type mask + positive mass, so an entry stays valid as the
 * particles move.  The routed producer reads live positions for everything
 * position-dependent, which is why no geometry is cached beside this. */
static struct ghost_local_tree_cache_t g_glt_cache = {
    .valid = 0,
    .NumPart_when_built = -1,
    .identity_epoch_when_built = -1,
    .eligible_type_mask_when_built = GHOST_TYPE_ALL,
    .num_pool = 0,
    .pool = NULL,
    .j_to_pool = NULL,
    .mark = NULL,
};

/* gx_policy_scaled_h — SSOT supply-side reach inside ghost_exchange.cc.
 *
 * The j-side reach this call searches with, for one particle, under the spec's
 * radius policy and scale.  Used by the routed producer on both sides of the
 * exchange, so sender and receiver agree on how far a supply particle reaches.
 *
 * Compile-flag gating on AGS_KernelRadius lives in nlr_radius_policy.h's
 * wrapper; never replicate the `#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE`
 * here. */
static inline double gx_policy_scaled_h(int j,
                                        mode_b_radius_policy_t radius_policy,
                                        double j_radius_scale,
                                        double safety_factor)
{
    return nlr_particle_symmetric_radius(P[j], radius_policy)
           * j_radius_scale * safety_factor;
}

/* Diagnostic counters. */
static long g_glt_cache_hits = 0;
static long g_glt_cache_misses = 0;

static void glt_cache_free(void)
{
    if(g_glt_cache.pool)      { free(g_glt_cache.pool);      g_glt_cache.pool = NULL; }
    if(g_glt_cache.j_to_pool) { free(g_glt_cache.j_to_pool); g_glt_cache.j_to_pool = NULL; }
    if(g_glt_cache.mark)      { free(g_glt_cache.mark);      g_glt_cache.mark = NULL; }
    g_glt_cache.valid = 0;
    g_glt_cache.NumPart_when_built = -1;
    g_glt_cache.identity_epoch_when_built = -1;
    g_glt_cache.eligible_type_mask_when_built = GHOST_TYPE_ALL;
    g_glt_cache.num_pool = 0;
}

/* The send set: which local pool slots this rank sends to which peer.
 *
 * The receiver walks hand it the slots they accept, peer by peer in ascending
 * order, a peer possibly over several calls.  It keeps each (peer, slot) once,
 * and when a peer is finished sorts that peer's slots, so the set reads out in
 * send order -- destination rank ascending, then pool index ascending -- with
 * no pass over the pool.  Work and memory are proportional to what is sent, not
 * to NTask x num_pool, which on a step sending a handful of particles is the
 * difference between reading a few slots and reading hundreds of millions.
 *
 * Duplicates are caught by a marker over the pool that is set for the peer being
 * filled and cleared again when that peer is finished, so it only ever holds one
 * peer's slots and is all zero between calls.  It is kept with the pool cache and
 * never cleared wholesale: clearing it per call would touch the whole pool on
 * every call, however little is sent.
 *
 * Emission never communicates, so ranks may emit in different numbers of calls
 * and from different receiver backends.  A failure is rank-local until the caller
 * agrees it across ranks. */
struct ghost_send_set {
    int  *slots;              /* unique matches, each peer's run contiguous */
    long  capacity, used;
    long *start;              /* [NTask] where peer t's run begins in slots */
    int  *count;              /* [NTask] how many slots peer t is sent */
    int   ntask_allocated;
    unsigned char *mark;      /* [num_pool] owned by the pool cache */
    int   num_pool;
    long  marks_outstanding;  /* marks set and not yet cleared; zero between calls */
    int   peer;               /* peer being filled, -1 when none is */
    int   last_peer;          /* highest peer begun this call, -1 before the first */
    int   failed;
    const int *j_to_pool;     /* [npart_bound] local particle -> pool slot, negative when not in the pool */
    int   npart_bound;
};
static struct ghost_send_set g_send_set = {NULL, 0, 0, NULL, NULL, 0, NULL, 0, 0, -1, -1, 0, NULL, 0};

/* Sort the finished peer's run into pool order and clear its marks. */
static void gx_send_set_finish_peer(struct ghost_send_set *s)
{
    if(s->peer < 0) {return;}
    int *run = s->slots + s->start[s->peer];
    const int n = s->count[s->peer];
    std::sort(run, run + n);
    for(int k = 0; k < n; k++) {s->mark[run[k]] = 0;}
    s->marks_outstanding -= n;
    s->peer = -1;
}

/* Clear whatever marks the unfinished peer holds and drop the output.  Finished
 * peers have cleared their own marks already, so this restores the all-zero
 * marker whenever it is called. */
static void gx_send_set_abandon(struct ghost_send_set *s)
{
    if(s->peer >= 0) {
        for(long k = s->start[s->peer]; k < s->used; k++) {s->mark[s->slots[k]] = 0;}
        s->marks_outstanding -= s->count[s->peer];
        s->peer = -1;
    }
    free(s->slots);
    s->slots = NULL;
    s->capacity = s->used = 0;
    s->failed = 1;
}

/* Returns 0, or nonzero with the set failed and nothing to abandon but its output. */
static int gx_send_set_begin(struct ghost_send_set *s, unsigned char *mark, int num_pool,
                             const int *j_to_pool, int npart_bound)
{
    s->failed = 1;
    s->peer = -1;
    s->last_peer = -1;
    s->used = 0;
    if(s->marks_outstanding != 0 || s->slots != NULL) {
        printf("ERROR: ghost send set on task %d was not left clean by the previous call: %ld marks "
               "still set, output buffer %s\n", ThisTask, s->marks_outstanding,
               s->slots ? "never taken or dropped" : "released");
        fflush(stdout);
        gizmo_request_controlled_stop(7735, "ghost_exchange: send set not clear at the start of a call",
                                      __FILE__, __LINE__, __FUNCTION__);
        return 1;
    }
    if(!mark) {return 1;}
    if(s->ntask_allocated < NTask) {
        free(s->start); free(s->count);
        s->start = (long *) malloc((size_t)NTask * sizeof(long));
        s->count = (int *)  malloc((size_t)NTask * sizeof(int));
        if(!s->start || !s->count) {
            free(s->start); free(s->count);
            s->start = NULL; s->count = NULL; s->ntask_allocated = 0;
            return 1;
        }
        s->ntask_allocated = NTask;
    }
    for(int t = 0; t < NTask; t++) {s->start[t] = 0; s->count[t] = 0;}
    s->capacity = 256;
    s->slots = (int *) malloc((size_t)s->capacity * sizeof(int));
    if(!s->slots) {s->capacity = 0; return 1;}
    s->mark = mark;
    s->num_pool = num_pool;
    s->j_to_pool = j_to_pool;
    s->npart_bound = (j_to_pool != NULL) ? npart_bound : 0;
    s->failed = 0;
    return 0;
}

static int gx_send_set_grow(struct ghost_send_set *s)
{
    if(s->capacity > LONG_MAX / 2) {return 1;}
    const long new_capacity = 2 * s->capacity;
    if((unsigned long)new_capacity > (unsigned long)(SIZE_MAX / sizeof(int))) {return 1;}
    int *grown = (int *) realloc(s->slots, (size_t)new_capacity * sizeof(int));
    if(!grown) {return 1;}
    s->slots = grown;
    s->capacity = new_capacity;
    return 0;
}

int gx_send_set_emit(struct ghost_send_set *s, int peer, const int *slots, int n)
{
    if(s->failed) {return 1;}
    /* A peer out of range, arriving out of order or a second time, or a slot
     * outside the pool, can only come from a broken receiver: its slots would be
     * sent out of order or split, or written outside the marker, which outlives
     * the call.  Stop rather than build on it. */
    int bad = (peer < 0 || peer >= NTask || n < 0 || (peer != s->peer && peer <= s->last_peer));
    for(int k = 0; k < n && !bad; k++) {bad = ((unsigned)slots[k] >= (unsigned)s->num_pool);}
    if(bad) {
        printf("ERROR: ghost send set on task %d was handed %d slots for peer %d after peer %d (pool %d); "
               "peers must arrive in ascending order, each in one run, with slots inside the pool\n",
               ThisTask, n, peer, s->last_peer, s->num_pool);
        fflush(stdout);
        gizmo_request_controlled_stop(7734, "ghost_exchange: send set handed an invalid emission",
                                      __FILE__, __LINE__, __FUNCTION__);
        s->failed = 1;
        return 1;
    }
    if(peer != s->peer) {
        gx_send_set_finish_peer(s);
        s->peer = peer;
        s->last_peer = peer;
        s->start[peer] = s->used;
    }
    for(int k = 0; k < n; k++) {
        const int p = slots[k];
        if(s->mark[p]) {continue;}
        if(s->used == s->capacity && gx_send_set_grow(s) != 0) {s->failed = 1; return 1;}
        s->mark[p] = 1;
        s->marks_outstanding++;
        s->slots[s->used++] = p;
        s->count[peer]++;
    }
    return 0;
}

/* ---- Utility: walk TopNodes to find which leaf a particle belongs to ---- */
static inline int ghost_toptree_leaf(peanokey key)
{
    int no = 0;
    peanokey mask = ((peanokey)7) << (3 * (BITS_PER_DIMENSION - 1));
    int shift = 3 * (BITS_PER_DIMENSION - 1);
    while(TopNodes[no].Daughter >= 0)
    {
        no = TopNodes[no].Daughter + (int)((key & mask) >> shift);
        mask >>= 3;
        shift -= 3;
    }
    return TopNodes[no].Leaf;
}


/*!
 * \brief Main ghost exchange routine. Call before neighbor operations.
 *
 * For each remote top-level leaf whose bounding region overlaps any local
 * particle's search sphere, imports all particles from that leaf.
 * Ghost particles are appended to P[]/CellP[] starting at NumPart.
 *
 * After all neighbor operations, call ghost_exchange_cleanup() to remove ghosts.
 *
 * safety_factor: multiplier on search_radius for the overlap criterion.
 *   1.0 = normal (previous-step hmax is accurate).
 *   >1.0 = inflate search radius to account for h-growth during density iteration
 *          (e.g. 2.0 on first timestep when densities are just guesses).
 */
static inline int ghost_type_passes(int ptype, unsigned int mask) { return (mask & (1u << (unsigned)ptype)) != 0u; }

/* Forward decls. */
/* Result of a single-backend ghost-exchange discovery attempt. The dispatcher
 * (ghost_exchange_impl) owns admission policy: a tile attempt that cannot fit
 * particle slots, or whose counts overflow the int MPI transport representation,
 * returns WITHOUT materialising ghosts (clean rollback), and the dispatcher
 * falls back to exact request-driven discovery. */
enum ghost_exchange_result {
    GHOST_EXCHANGE_COMPLETED = 0,
    GHOST_EXCHANGE_PARTICLE_CAPACITY_EXCEEDED,
    GHOST_EXCHANGE_COUNT_RANGE_EXCEEDED
};

static ghost_exchange_result ghost_exchange_request_driven_impl(const struct ghost_exchange_spec_t *spec);

/* Every spec uses the walk-export routed producer (sender fine-tree export +
 * bounded receiver walk).  Kept as a named predicate because the two search
 * modes reach that conclusion for different reasons, recorded here.
 *
 * ONEWAY: its reach is the query's own radius, which the sender knows exactly,
 * so routing needs nothing from the supply side.  The traversal AND export
 * opener are both R_open = h_q (mode_b_local_walker.cc, the R_open branch and
 * the topleaf export re-test), and the accept ignores h_j entirely
 * (ghost_exchange_functions.h gx_pair_accept_wrap_and_test, ONEWAY branch).
 * Query h already carries the spec safety factor, so a widened query cannot
 * outgrow the opener.  KEEP THOSE THREE IN STEP: if ONEWAY accept ever gains an
 * h_j term, or the opener stops using h_q, this must be re-derived.
 *
 * SYMMETRIC: the supply-side reach is bounded by the per-type node band the
 * sender opener walks against.  The band is seeded from the conservative union
 * of every radius source a leaf policy can select
 * (force_hmax_per_type_particle_radius, MODE_B_RADIUS_ALL_SOURCES), so it
 * dominates any radius_policy a spec can declare -- which is why this is a
 * property of the BAND, not of the loop.  A new radius-policy bit therefore
 * extends that union; it never needs a per-loop claim here.  A safety factor
 * above 1 (TURB_DIFF_DYNAMIC) is not a disqualifier either: both walks scale
 * the j-side reach they search with by the spec's safety factor, so their reach
 * equals the accept's at any safety.
 * Rank-uniform -- search_mode is a spec constant, identical on every rank, so
 * this never splits ranks across a collective. */
static inline int gx_walk_export_eligible(const struct ghost_exchange_spec_t *spec)
{
    return spec != NULL;
}

static void ghost_exchange_impl(const struct ghost_exchange_spec_t *spec)
{
    /* Tiny-N corridor counter: increments on API entry, before any
     * dispatch. Mode B paths in run_neighbor_loop must NOT enter this
     * function. Counts even single-rank early-out cases by design. See
     * declarations/lifecycle_counters.h. */
    g_ghost_import_counter++;

    /* Every request is answered by request-driven discovery, whatever its search
     * mode: ONEWAY opens on the query's own radius, SYMMETRIC on the per-type
     * node band that bounds any radius policy.  Neither is a property of a
     * specific loop, so there is nothing here to key on. */
    ghost_exchange_request_driven_impl(spec);
}

/* Public entry for new-style callers that build their own spec literal at
 * the call site (mech_fb_v1 onward). The literal IS the single source of
 * truth for that loop's physics — to flip mode / supply_mask / query list,
 * edit the literal at the caller. Dispatch keys only on spec fields
 * (explicit query list, or search_mode == NGB_SEARCH_ONEWAY). */
extern "C" void ghost_exchange_run(const struct ghost_exchange_spec_t *spec)
{
    ghost_exchange_impl(spec);
}

static inline int ghost_particle_slots_fit(long long required)
{
    if(required > (long long)All.MaxPart) {return 0;}
    if(All.TotN_gas > 0 && required > (long long)All.MaxPartGas) {return 0;}
    return 1;
}

/* SSOT for "what a send slot contains": one exported particle's P (+ gas CellP,
   zeroed for non-gas). Used by BOTH import pack loops and ghost_refresh_values()
   so the refresh cannot drift from import. src_P/src_CellP are the value source
   (production import + refresh pass P/CellP; the refresh harness passes a copy). */
static inline void gx_pack_send_slot(const struct particle_data *src_P,
                                     const struct gas_cell_data *src_CellP,
                                     int j,
                                     struct particle_data *dst_P,
                                     struct gas_cell_data *dst_CellP)
{
    *dst_P = src_P[j];
    if(src_P[j].Type == 0 && j < N_gas) *dst_CellP = src_CellP[j];
    else                                memset(dst_CellP, 0, sizeof(struct gas_cell_data));
}

/* Time every particle in the live ghost pool is current to, or -1 when there is
 * no live certified pool. Set once per import/refresh, after every exported slot
 * has been advanced to that time, and cleared wherever the pool it describes
 * stops existing.
 *
 * Consumers compare it against the time they need rather than trusting it
 * outright, so a stale or cleared stamp costs a per-particle check and never a
 * wrong answer. Nothing can corrupt it in the other direction, because no path
 * moves a particle backwards -- drift_particle refuses that outright. */
static integertime g_ghost_pool_current_ti = -1;
integertime ghost_pool_current_ti(void) { return g_ghost_pool_current_ti; }

/* Advance every particle about to be exported to the current time, so nothing
 * goes on the wire behind the time its receiver will read it at.
 *
 * A ghost shipped behind the current time has to be advanced again by every rank
 * that receives it, in neighbour-list order, inside the barrier-sensitive part of
 * the dispatch -- work the owner can do once, in pool order, here. It also lets a
 * receiver rely on the pool's currency instead of re-deriving it per particle,
 * and it is a prerequisite for ever shipping a compact ghost record: the fields
 * the drift needs are exactly the ones a compaction would drop.
 *
 * Threaded, because the drift is real per-particle work -- it runs the implicit
 * thermochemistry solve through set_eos_pressure. That requires each particle be
 * visited exactly once: a particle exported to several ranks appears once per
 * destination, and two threads testing its Ti_current before either writes would
 * advance it twice. Establishing the distinct set first is what keeps the
 * parallel pass safe, at one stamp compare per slot; the stamp is
 * generation-counted so it never needs clearing between calls. */
/* Returns 0 when every exported slot stands at t_now, nonzero when the drift that
 * would have advanced them did not complete -- see drift_particles_batch. The pool
 * stamp the callers publish rests on this having succeeded on EVERY rank, so a
 * caller that stamps regardless would assert a currency the pool does not have,
 * and the scan that would notice is the one the stamp suppresses. */
/* The distinct particles behind a time, gathered from one or more index lists and
 * then drifted together: begin, add each list, drift.  The batch drift threads over
 * the particles it is handed, so each must appear once -- two threads advancing the
 * same particle would drift it twice.  The stamp is generation-counted, so it never
 * needs clearing between calls, and it is sized once to the particle count.
 *
 * The work list holds only the particles actually behind, so it is both the drift
 * list and the h-dirty list.  drift_particle rescales KernelRadius, so a mark is owed
 * for a particle that moved -- but marking every listed particle would drive the
 * dirty tracker toward a full-pool refresh on the many steps where nothing was
 * behind.
 *
 * begin and add return nonzero when their storage cannot grow; nothing has been
 * drifted then, and a controlled stop is requested, because the particles these
 * lists hold would otherwise be used behind the current time. */
static std::vector<unsigned int> gx_behind_seen;
static unsigned int gx_behind_seen_gen = 0;
static std::vector<int> gx_behind_list;

static int gx_drift_distinct_storage_failed(void)
{
    printf("ghost_exchange: task %d could not grow the list of particles to drift\n", ThisTask);
    fflush(stdout);
    gizmo_request_controlled_stop(7739, "ghost_exchange: could not allocate the drift work list",
                                  __FILE__, __LINE__, __FUNCTION__);
    return 1;
}

static int gx_drift_distinct_begin(void)
{
    try {
        if(gx_behind_seen.size() < (size_t)NumPart) {gx_behind_seen.assign((size_t)NumPart, 0u);}
    } catch(const std::bad_alloc &) {return gx_drift_distinct_storage_failed();}
    if(++gx_behind_seen_gen == 0u) {std::fill(gx_behind_seen.begin(), gx_behind_seen.end(), 0u); gx_behind_seen_gen = 1u;}
    gx_behind_list.clear();
    return 0;
}

static int gx_drift_distinct_add(const int *idx, int n, integertime t_now)
{
    try {
        for(int k = 0; k < n; k++) {
            const int j = idx[k];
            if(j < 0 || j >= NumPart) {continue;}
            if(gx_behind_seen[(size_t)j] == gx_behind_seen_gen) {continue;}
            gx_behind_seen[(size_t)j] = gx_behind_seen_gen;
            if(P[j].Ti_current != t_now) {gx_behind_list.push_back(j);}
        }
    } catch(const std::bad_alloc &) {return gx_drift_distinct_storage_failed();}
    return 0;
}

/* Returns 0 when every gathered particle stands at t_now, nonzero when the drift that
 * would have advanced them did not complete -- see drift_particles_batch. */
static int gx_drift_distinct_to_time(integertime t_now)
{
    if(gx_behind_list.empty()) {return 0;}
    const int n_behind = (int)gx_behind_list.size();
    const int *behind_idx = gx_behind_list.data();
    const int drift_status = drift_particles_batch(behind_idx, n_behind, t_now);
    gizmo_mark_kernel_radius_dirty_indices(behind_idx, n_behind);
    return drift_status;
}

static int gx_certify_send_list_current(const int *home_idx, int n_slots, integertime t_now)
{
    if(!home_idx || n_slots <= 0) {return 0;}
    if(gx_drift_distinct_begin() != 0 || gx_drift_distinct_add(home_idx, n_slots, t_now) != 0) {return 1;}
    return gx_drift_distinct_to_time(t_now);
}

/* The receiver's half of discovery, shared by both backends.  A backend walks its
 * envelopes and records, per envelope, the local particles that may be neighbours
 * once drifted to the current time -- a superset, because a particle it did not
 * drift may have moved since its position was stored.  This settles each row:
 * drift the distinct candidates of the whole batch to the current time, keep
 * exactly those the query accepts at their current position and radius, and hand
 * their pool slots to the send set.  Accepting on stored positions instead would
 * both miss particles that moved into reach and send ones that moved out of it.
 *
 * Nothing is emitted unless every candidate was drifted: storage that cannot grow,
 * or a drift that did not complete, fails the batch before its first emission.
 * Rows must arrive with peers in ascending order (the send set's contract); a row's
 * candidates are compacted in place, so no second buffer is needed.  The accept
 * pass and the emission pass are separate loops, so the first can be threaded on
 * its own if it ever costs enough to matter. */
int gx_send_set_accept_rows(struct ghost_send_set *s, struct gx_candidate_row *rows, long n_rows,
                            int search_mode, mode_b_radius_policy_t radius_policy,
                            double j_radius_scale, double safety_factor)
{
    if(s->failed) {return 1;}
    const integertime t_now = All.Ti_Current;
    /* The pool is reused only while it was built over exactly this particle set, and holds every
       particle of positive mass; a candidate the walks could return is therefore always in it.  One
       that is not, or a map of another size, means the cache is stale or the walk is malformed, and
       dropping the candidate would silently narrow the answer -- so the run stops instead. */
    if(!s->j_to_pool || s->npart_bound != NumPart) {
        printf("ghost_exchange: task %d supply map covers %d particles, the rank holds %d\n", ThisTask, s->npart_bound, NumPart);
        fflush(stdout);
        gizmo_request_controlled_stop(7738, "ghost_exchange: receiver supply map does not match the particle set",
                                      __FILE__, __LINE__, __FUNCTION__);
        s->failed = 1;
        return 1;
    }

    /* Every candidate must be a particle of this set with a supply-pool slot (above). */
    for(long r = 0; r < n_rows; r++) {
        const struct gx_candidate_row *row = &rows[r];
        for(int c = 0; c < row->count; c++) {
            const int j = row->local_index[c];
            const int pp = (j >= 0 && j < NumPart) ? s->j_to_pool[j] : -1;
            if(pp < 0 || pp >= s->num_pool) {
                printf("ghost_exchange: task %d receiver candidate %d (of 0..%d) has no supply-pool slot\n", ThisTask, j, NumPart - 1);
                fflush(stdout);
                gizmo_request_controlled_stop(7738, "ghost_exchange: receiver candidate outside the supply pool",
                                              __FILE__, __LINE__, __FUNCTION__);
                s->failed = 1;
                return 1;
            }
        }
    }

    if(gx_drift_distinct_begin() != 0) {s->failed = 1; return 1;}
    for(long r = 0; r < n_rows; r++) {
        if(gx_drift_distinct_add(rows[r].local_index, rows[r].count, t_now) != 0) {s->failed = 1; return 1;}
    }
    if(gx_drift_distinct_to_time(t_now) != 0) {s->failed = 1; return 1;}

    /* Exact acceptance at current positions; accepted pool slots replace the
     * candidates at the front of each row.  Mass is tested again because the drift
     * can remove a particle (an outflow boundary zeroes it).  Reach comes from THIS caller's spec, not
     * from whatever policy the cached pool happened to be built under: pool
     * membership is reused across differing radius policies, and only the tile/BVH/
     * compact geometry may use the build-time policy, because its leaf h was baked
     * with it. */
    for(long r = 0; r < n_rows; r++) {
        struct gx_candidate_row *row = &rows[r];
        const struct gx_export_envelope_t *e = row->envelope;
        int n_accepted = 0;
        for(int c = 0; c < row->count; c++) {
            const int j = row->local_index[c];
            if(!(P[j].Mass > 0)) {continue;}
            const double hj = gx_policy_scaled_h(j, radius_policy, j_radius_scale, safety_factor);
            if(gx_pair_accept_wrap_and_test(e->pos[0] - (double)P[j].Pos[0],
                                            e->pos[1] - (double)P[j].Pos[1],
                                            e->pos[2] - (double)P[j].Pos[2],
                                            e->h, hj, search_mode)) {
                row->local_index[n_accepted++] = s->j_to_pool[j];
            }
        }
        row->count = n_accepted;
    }

    for(long r = 0; r < n_rows; r++) {
        if(rows[r].count > 0 &&
           gx_send_set_emit(s, rows[r].peer, rows[r].local_index, rows[r].count) != 0) {return 1;}
    }
    return 0;
}

/* SSOT for the forward particle+cell transport: pack the exported slots and carry
   them with two typed Alltoallv calls (P always; CellP only when gas exists
   globally) with element-unit counts. Used by BOTH import impls (materialising
   ghosts at &P[NumPart]) and ghost_refresh_values() (overwriting existing ghost
   slots at &P[NumPart_before_ghost]). Verbose per-impl diagnostics stay at the call
   sites; only the pack and the transport are factored here.

   The payload is carried in rounds rather than staged whole. A crowded step exports
   a large multiple of a rank's own particles -- a near-all-active FIF step asked for
   1110 MB of particle_data + gas_cell_data at once, more than the whole working pool
   -- and staging that in one piece made this transport something that has to fit,
   which is exactly what the pool is not sized for: it counts only what cannot be
   broken up, and expects everything else to take more rounds when there is less room
   (see arena_megabytes_from_tenants). Rounds put this transport back under that rule,
   so a bigger pool buys fewer rounds and never decides whether the step works.

   A round carries, for every peer, the next slice of that peer's run. Both ends
   apply the same rule to the same per-peer counts, which were exchanged before this
   is called, so each rank already knows what it will receive in each round and no
   further count exchange is needed. Slices land at their final offsets, so the
   delivered pool is identical to what a single round produced. */
static void gx_pack_and_forward_particle_exchange(const int *send_home_idx,
                                                  const int *send_count, const int *send_disp,
                                                  struct particle_data *dst_P,
                                                  struct gas_cell_data *dst_CellP,
                                                  const int *recv_count, const int *recv_disp)
{
    const size_t slot_bytes = sizeof(struct particle_data) + sizeof(struct gas_cell_data);

    /* How many slots per peer a round may carry. Sized from the room this rank has
       now, shared out over the peers a round touches, and reduced to what the
       tightest rank can manage so every rank cuts its runs at the same place. At
       least one slot per peer: a round that cannot be afforded is reported below
       rather than skipped, or the transport would never finish.

       Half the free room, not all of it: the two staging buffers are taken while
       the callers still hold their own working arrays, and a round sized to the
       last free byte would be refused by the check below and turn a run that fits
       into a stop. */
    int slots_per_peer = 1;
    long long tightest_free = 0;
    {
        const size_t per_peer = slot_bytes * (size_t) ((NTask > 0) ? NTask : 1);
        long long affordable = (per_peer > 0) ? (long long) (((size_t) FreeBytes / 2) / per_peer) : 1;
        if(affordable < 1) {affordable = 1;}
        if(affordable > INT_MAX) {affordable = INT_MAX;}
        /* The room left on the tightest rank rides along with the agreement, so
           reporting it costs no collective of its own. */
        long long want[2] = {affordable, (long long) FreeBytes}, least[2] = {1, 0};
        MPI_Allreduce(want, least, 2, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
        slots_per_peer = (least[0] > 0) ? (int) least[0] : 1;
        tightest_free  = least[1];
    }

    /* How many rounds cover the longest run anywhere. A rank with nothing to send
       still has to enter every round its peers send in, so the longest run in either
       direction on any rank sets the count for all of them.

       The first round is the largest -- every peer contributes a full slice, and a
       later round's slice can only be shorter -- and each round gives its buffers
       back before the next takes any, so asking once whether the first round fits
       answers the question for all of them. */
    int rounds = 1;
    long long heaviest_rank_bytes = 0;
    size_t staging_bytes = 0;
    int short_local = 0, short_any = 0;
    {
        long long longest_run = 0, carried = 0, first_round = 0;
        for(int t = 0; t < NTask; t++) {
            if(send_count[t] > longest_run) {longest_run = send_count[t];}
            if(recv_count[t] > longest_run) {longest_run = recv_count[t];}
            carried += send_count[t];
            first_round += (send_count[t] < slots_per_peer) ? send_count[t] : slots_per_peer;
        }
        const size_t first_slots = (size_t) (first_round > 0 ? first_round : 1);
        staging_bytes = gizmo_mymalloc_rounded_size(first_slots * sizeof(struct particle_data))
                      + gizmo_mymalloc_rounded_size(first_slots * sizeof(struct gas_cell_data));
        short_local = gizmo_alloc_fits_this_rank(staging_bytes, 2) ? 0 : 1;
        /* What this rank puts on the wire, and whether it can afford a round, travel
           with the same agreement: the figure worth reporting belongs to whichever
           rank carries most, and that is not reliably rank 0. One rank short stops
           them all, together, before a byte is taken. */
        long long mine[3] = {longest_run, carried * (long long) slot_bytes, short_local};
        long long most[3] = {0, 0, 0};
        MPI_Allreduce(mine, most, 3, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        rounds = (int) ((most[0] + slots_per_peer - 1) / slots_per_peer);
        if(rounds < 1) {rounds = 1;}
        heaviest_rank_bytes = most[1];
        short_any = (most[2] != 0);
    }

    if(short_any)
    {
        if(short_local) {
            /* Two different shortages reach here: not enough bytes free, or no room
               left in the arena's table of live blocks. Say which, because raising
               the pool answers the first and does nothing at all for the second. */
            if((long long) staging_bytes > (long long) FreeBytes) {
                printf("Ghost exchange: rank %d needs %g MB to carry one round of exported particles but "
                       "has %g MB of working memory left. Raise Working_Mem_Pool_Per_Task_in_MB, or spread "
                       "the run over more ranks. Stopping.\n",
                       ThisTask, (double) staging_bytes / (1024.0 * 1024.0),
                       (double) FreeBytes / (1024.0 * 1024.0));
            } else {
                printf("Ghost exchange: rank %d cannot take the 2 staging buffers for one round: the "
                       "working memory has %g MB free, so this is not a shortage of memory but of its "
                       "table of live blocks. Something above is holding an unusual number of "
                       "allocations. Stopping.\n",
                       ThisTask, (double) FreeBytes / (1024.0 * 1024.0));
            }
            fflush(stdout);
        }
        gizmo_request_controlled_stop(7727,
            "ghost_exchange: not enough working memory to carry one round of the exported particles",
            __FILE__, __LINE__, __FUNCTION__);
        gizmo_exit_bad_stop_if_requested("ghost_exchange:forward_round");
        return;
    }

    /* Say how the exchange was split, on a change of round count or once the load
       has grown by a quarter since the last word -- the same gating the tree's own
       transport report uses, so a run that never needs a second round says this once
       and a run growing toward one shows it coming. */
    if(ThisTask == 0)
    {
        static int last_rounds_reported = -1;
        static long long last_heaviest_reported = 0;
        if(rounds != last_rounds_reported ||
           heaviest_rank_bytes > last_heaviest_reported + last_heaviest_reported / 4)
        {
            last_rounds_reported    = rounds;
            last_heaviest_reported  = heaviest_rank_bytes;
            printf("Ghost exchange: %d round(s); heaviest rank carries %g MB, working memory free on the "
                   "tightest rank %g MB\n",
                   rounds, (double) heaviest_rank_bytes / (1024.0 * 1024.0),
                   (double) tightest_free / (1024.0 * 1024.0));
            fflush(stdout);
        }
    }

    int *round_send_count = (int *) mymalloc("gx_fwd_sc", NTask * sizeof(int));
    int *round_send_disp  = (int *) mymalloc("gx_fwd_sd", NTask * sizeof(int));
    int *round_recv_count = (int *) mymalloc("gx_fwd_rc", NTask * sizeof(int));
    int *round_recv_disp  = (int *) mymalloc("gx_fwd_rd", NTask * sizeof(int));

    for(int r = 0; r < rounds; r++)
    {
        const int taken = r * slots_per_peer;
        int staged = 0;
        for(int t = 0; t < NTask; t++)
        {
            int s = send_count[t] - taken; if(s < 0) {s = 0;} if(s > slots_per_peer) {s = slots_per_peer;}
            int v = recv_count[t] - taken; if(v < 0) {v = 0;} if(v > slots_per_peer) {v = slots_per_peer;}
            round_send_count[t] = s;
            round_send_disp[t]  = staged;
            staged += s;
            round_recv_count[t] = v;
            /* Straight into the slot this peer's run occupies in the delivered pool,
               so nothing has to be moved afterwards and the result does not depend
               on which round carried it. A peer whose run is already finished keeps
               its run's own start, which is always a real offset: an empty slice must
               still name somewhere inside the buffer, because the transport may form
               the address before it notices the count is zero. */
            round_recv_disp[t]  = (v > 0) ? (recv_disp[t] + taken) : recv_disp[t];
        }

        const size_t staging_slots = (size_t) (staged > 0 ? staged : 1);
        struct particle_data *send_P = (struct particle_data *) mymalloc("gx_fwd_sP",
            staging_slots * sizeof(struct particle_data));
        struct gas_cell_data *send_CellP = (struct gas_cell_data *) mymalloc("gx_fwd_sC",
            staging_slots * sizeof(struct gas_cell_data));

        for(int t = 0; t < NTask; t++)
        {
            /* Only step into the list for a peer this round still has slots for:
               a run that finished in an earlier round would address past its end. */
            if(round_send_count[t] <= 0) {continue;}
            const int *home = send_home_idx + send_disp[t] + taken;
            for(int k = 0; k < round_send_count[t]; k++)
            {
                const int j = home[k];
                const int off = round_send_disp[t] + k;
                /* A run that stopped early on its own guard leaves its remaining slots
                   unclaimed. Send them as zeros: the receiver counts them but no home
                   particle stands behind them, and staging is reused between rounds,
                   so carrying them as they lie would put stale bytes on the wire. */
                if(j < 0) {
                    memset(&send_P[off], 0, sizeof(struct particle_data));
                    memset(&send_CellP[off], 0, sizeof(struct gas_cell_data));
                    continue;
                }
                gx_pack_send_slot(P, CellP, j, &send_P[off], &send_CellP[off]);
            }
        }

        gizmo_mpi_alltoallv_typed(send_P, round_send_count, round_send_disp,
                                  dst_P, round_recv_count, round_recv_disp,
                                  sizeof(struct particle_data), MPI_COMM_WORLD);
        /* Only meaningful when the simulation has any gas particles globally. With
           TotN_gas==0 (N-body / DM-only runs) CellP is allocated to size 0, so writing
           to dst_CellP would dereference out of bounds -- and no gas ghost can exist
           if no gas exists anywhere. */
        if(All.TotN_gas > 0) {
            gizmo_mpi_alltoallv_typed(send_CellP, round_send_count, round_send_disp,
                                      dst_CellP, round_recv_count, round_recv_disp,
                                      sizeof(struct gas_cell_data), MPI_COMM_WORLD);
        }

        myfree(send_CellP);
        myfree(send_P);
    }

    myfree(round_recv_disp);
    myfree(round_recv_count);
    myfree(round_send_disp);
    myfree(round_send_count);
}

struct gx_query_t {
    double pos[3];
    double h;
    int    type;       /* for caller diagnostics; predicate uses h_i directly */
    int    _pad;
};

/* gx_export_envelope_t (the sender->supply wire record) is declared in
 * mesh/neighbor_list.h: the device receiver traversal compiles in another
 * translation unit and consumes the same records. */

/* Walk-export routed producer — the discovery path for every spec that passes
 * gx_walk_export_eligible().  Sender: per local query mode_b_walk_and_export -> per-peer
 * NodeList -> fixed-size envelopes -> Alltoallv.  Receiver: mode_b_walk_from_start_nodes
 * (resume from the exported NodeList) -> gx_pair_accept_wrap_and_test (the ghost-exchange SSOT predicate)
 * -> the send set (which local pool slots go to which peer), which Steps 4-6 install from.
 * MODE-GENERIC (search_mode is passed through): this is the install target for both search
 * modes, so ONEWAY and SYMMETRIC discover on ONE substrate rather than two.
 *
 * COLLECTIVE-SAFE (C MPI buffers): every C allocation preceding a collective — the index
 * arrays before the Alltoall, the envelope buffers before the Alltoallv, the send set's growth
 * — is Allreduce-checked, so a NULL on ANY rank makes ALL ranks return the same status and
 * ranks never diverge across a collective.  NOT covered: the C++ containers
 * (ModeBExportSink, the per-peer envelope vectors) throw std::bad_alloc rank-locally rather
 * than returning NULL, so an allocation failure there aborts that rank instead of returning
 * a uniform status.  On GX_WALK_EXPORT_OK the send set is finished and holds the result; on
 * any other status it holds nothing and its marker is clear.
 *
 * THREADING (host): the sender-query loop and the received-envelope walk are both
 * `omp parallel for` — per-thread export sink and per-thread send buffers on the sender,
 * pre-sized per-envelope candidate slots on the receiver, so each thread writes only its own
 * index.  The walker's lazy node drift is the one shared mutation and is serialized under
 * critical(_modebdrift_).  The candidates' drift is threaded over distinct particles
 * (gx_send_set_accept_rows); merge, accept and send-set emission run serially, which keeps
 * the resulting SET order-independent and therefore deterministic across thread counts.
 *
 * MEMORY SHAPE: the threaded receiver materializes one candidate list per received envelope
 * before the serial accept pass.  That is bounded by the exported envelope volume, which is
 * measured sparse; a pathologically clustered geometry with very large tot_r would grow it,
 * so it is a known scaling watch-point rather than a fixed bound. */
/* RECEIVER_FAIL: some rank's receiver could not build its send set -- it could not
 * grow, or it met a state that has already requested its own stop, printed by the
 * rank that met it. */
enum { GX_WALK_EXPORT_OK = 0, GX_WALK_EXPORT_UNAVAILABLE = 1, GX_WALK_EXPORT_ALLOC_FAIL = 2,
       GX_WALK_EXPORT_RECEIVER_FAIL = 3 };
struct gx_walk_export_result {
    int status;   /* a GX_WALK_EXPORT_* value, uniform across ranks */
};

static void gx_walk_export_discover(
    const struct ghost_exchange_spec_t *spec,
    const struct gx_query_t *local_queries, int n_local_queries,
    int num_pool, unsigned int supply_mask, int search_mode,
    struct ghost_send_set *send_set,
    struct gx_walk_export_result *res)
{
    if(res) memset(res, 0, sizeof(*res));
    /* The walk searches with the SAME j-side reach the accept admits with.
     * The accept uses gx_policy_scaled_h = radius(j,policy) * j_radius_scale *
     * safety_factor, so both factors are folded here; the query side already
     * carries safety (queries are built as h*safety). Were the walk to search
     * with a smaller reach than the accept admits, it would silently discover
     * fewer pairs than the caller asked for -- which is why a safety factor
     * above 1 previously had to fall back to broadcast. */
    const double walker_j_reach_scale = spec->j_radius_scale * spec->safety_factor;
    /* (a) tree availability — collective all-or-none (a rank-local skip would deadlock the
     * envelope Alltoallv below / the caller's compare Allreduce). */
    int ok_local = (All.TreeNodeIndexBase > 0 && Nodes != NULL && Nextnode != NULL) ? 1 : 0;
    int ok_all = 0;
    MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if(!ok_all) { if(res) res->status = GX_WALK_EXPORT_UNAVAILABLE; return; }

    /* (b) SENDER: build per-peer envelope lists (export is a byproduct of the walk).
     * THREADED, following the same shape the runner uses: the topleaf map is built once
     * and READ-ONLY during the walk (shared); each thread uses its own ModeBExportSink + its own
     * per-peer send buffers; the walker's node lazy-drift is race-safe (omp
     * critical(_modebdrift_) + release/acquire).  Merge is serial.  The routed SET is
     * unchanged (order-independent bitmap); only per-peer envelope ORDER differs (D6: FP-reorder only). */
    ModeBTopleafMap map; map.build();
    std::vector<std::vector<struct gx_export_envelope_t>> send(NTask);
    {
#ifdef _OPENMP
        int nthr = omp_get_max_threads();
#else
        int nthr = 1;
#endif
        if(nthr < 1) nthr = 1;
        std::vector<ModeBExportSink> tsink(nthr);
        for(int th = 0; th < nthr; th++) tsink[th].ensure_size(NTask);
        std::vector<std::vector<std::vector<struct gx_export_envelope_t>>> tsend(nthr);
        for(int th = 0; th < nthr; th++) tsend[th].assign(NTask, std::vector<struct gx_export_envelope_t>{});
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
        for(int qi = 0; qi < n_local_queries; qi++) {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            ModeBExportSink &sink = tsink[tid];
            sink.clear_all();
            mode_b_walk_and_export(local_queries[qi].pos, local_queries[qi].h,
                                   supply_mask, search_mode, spec->radius_policy,
                                   /*cand_out=*/NULL, map, sink, walker_j_reach_scale);
            for(int t = 0; t < NTask; t++) {
                if(t == ThisTask) continue;
                long nn = (long)sink.nodes_per_peer[t].size();
                if(nn <= 0) continue;
                const std::vector<int> &nl = sink.nodes_per_peer[t];
                for(long base = 0; base < nn; base += NODELISTLENGTH) {
                    struct gx_export_envelope_t e;
                    e.pos[0] = local_queries[qi].pos[0];
                    e.pos[1] = local_queries[qi].pos[1];
                    e.pos[2] = local_queries[qi].pos[2];
                    e.h = local_queries[qi].h;
                    e.n_nodes = (int)((nn - base < NODELISTLENGTH) ? (nn - base) : NODELISTLENGTH);
                    for(int k = 0; k < e.n_nodes; k++) e.nodes[k] = nl[base + k];
                    e._pad = 0;
                    tsend[tid][t].push_back(e);
                }
            }
        }
        /* serial merge: concatenate per-thread envelopes per peer. */
        for(int th = 0; th < nthr; th++) {
            for(int t = 0; t < NTask; t++)
                send[t].insert(send[t].end(), tsend[th][t].begin(), tsend[th][t].end());
        }
    }

    /* (c) EXCHANGE: counts Alltoall + typed Alltoallv of fixed-size envelopes. */
    int *sc = (int *) calloc(NTask, sizeof(int));
    int *sd = (int *) calloc(NTask, sizeof(int));
    int *rc = (int *) calloc(NTask, sizeof(int));
    int *rd = (int *) calloc(NTask, sizeof(int));
    if(!sc || !sd || !rc || !rd) {
        /* index-array alloc failure is per-rank; make it collective before the Alltoall. */
        int a_ok_local = 0, a_ok_all = 0;
        MPI_Allreduce(&a_ok_local, &a_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        free(sc); free(sd); free(rc); free(rd);
        if(res) res->status = GX_WALK_EXPORT_ALLOC_FAIL; return;
    } else {
        int a_ok_local = 1, a_ok_all = 0;
        MPI_Allreduce(&a_ok_local, &a_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(!a_ok_all) { free(sc); free(sd); free(rc); free(rd); if(res) res->status = GX_WALK_EXPORT_ALLOC_FAIL; return; }
    }
    /* Per-peer counts are bounded by this rank's query set, but the prefix sums and
     * totals are not: the typed Alltoallv below consumes int displacements, so a value
     * past INT_MAX cannot be represented and would wrap into a NEGATIVE displacement,
     * making the memcpy that fills sendbuf an out-of-bounds write before MPI is ever
     * reached.  Same failure mode, and the same collective treatment, as the Step-4
     * transport guard: detect before mutating anything, agree across ranks, and report
     * a uniform status rather than truncating. */
    long tot_s = 0, tot_r = 0;
    int env_range_ok = 1;
    for(int t = 0; t < NTask; t++) {
        long n = (long)send[t].size();
        if(n > INT_MAX) { env_range_ok = 0; n = 0; }
        sc[t] = (int)n;
    }
    MPI_Alltoall(sc, 1, MPI_INT, rc, 1, MPI_INT, MPI_COMM_WORLD);
    for(int t = 0; t < NTask; t++) {
        if(tot_s <= INT_MAX) sd[t] = (int)tot_s; else { sd[t] = 0; env_range_ok = 0; }
        if(tot_r <= INT_MAX) rd[t] = (int)tot_r; else { rd[t] = 0; env_range_ok = 0; }
        tot_s += sc[t];
        tot_r += rc[t];
    }
    if(tot_s > INT_MAX || tot_r > INT_MAX) env_range_ok = 0;
    {
        int range_ok_all = 0;
        MPI_Allreduce(&env_range_ok, &range_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(!range_ok_all) {
            if(ThisTask == 0) {
                printf("ERROR: walk-export envelope counts exceed int MPI transport range (caller=%s)\n",
                       spec->caller_name ? spec->caller_name : "?");
                fflush(stdout);
            }
            free(sc); free(sd); free(rc); free(rd);
            if(res) res->status = GX_WALK_EXPORT_UNAVAILABLE;
            return;
        }
    }
    struct gx_export_envelope_t *sendbuf = (struct gx_export_envelope_t *) malloc((size_t)(tot_s > 0 ? tot_s : 1) * sizeof(struct gx_export_envelope_t));
    struct gx_export_envelope_t *recv    = (struct gx_export_envelope_t *) malloc((size_t)(tot_r > 0 ? tot_r : 1) * sizeof(struct gx_export_envelope_t));
    /* buffer alloc — collective before the Alltoallv (a NULL recv would fault the collective). */
    int buf_ok_local = (sendbuf && recv) ? 1 : 0, buf_ok_all = 0;
    MPI_Allreduce(&buf_ok_local, &buf_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if(!buf_ok_all) { free(sendbuf); free(recv); free(sc); free(sd); free(rc); free(rd);
                      if(res) res->status = GX_WALK_EXPORT_ALLOC_FAIL; return; }
    for(int t = 0; t < NTask; t++)
        if(sc[t] > 0) memcpy(sendbuf + sd[t], send[t].data(), (size_t)sc[t] * sizeof(struct gx_export_envelope_t));
    gizmo_mpi_alltoallv_typed(sendbuf, sc, sd, recv, rc, rd,
                              sizeof(struct gx_export_envelope_t), MPI_COMM_WORLD);
    free(sendbuf); free(sc); free(sd);

    /* (d) RECEIVER: bounded resume-walk from the exported NodeList -> SSOT accept -> send set.
     * Two interchangeable backends feed the same send set, which keeps each (peer, slot)
     * once and reads out in send order, so which one ran is not observable downstream.
     *
     * The device traversal is tried first and answers whenever the tree mirror
     * is current.  It declines rank-locally otherwise, and this window holds no
     * collectives, so a rank that declines simply does the work itself and no
     * other rank needs to agree.  Declining is for an unusable device tree state
     * or an allocation failure -- never for a disagreement, which would be a bug
     * to fix rather than to route around.  A decline happens before the device has
     * emitted anything; a failure after it has emitted is not a decline, because the
     * host walk cannot run over a half-filled set, so it fails this rank instead. */
    int receiver_ok = (gx_send_set_begin(send_set, g_glt_cache.mark, num_pool,
                                         g_glt_cache.j_to_pool, g_glt_cache.NumPart_when_built) == 0);
    int receiver_done_on_device = 0;
    if(receiver_ok && tot_r > 0) {
        std::vector<int> envelope_peer((size_t)tot_r, -1);
        for(int t = 0; t < NTask; t++) {
            for(int r = 0; r < rc[t]; r++) {envelope_peer[(size_t)rd[t] + r] = t;}
        }
        const int device_outcome =
            gx_device_receiver_walk(recv, tot_r, envelope_peer.data(),
                                    supply_mask, search_mode,
                                    spec->radius_policy, spec->j_radius_scale, spec->safety_factor,
                                    g_glt_cache.j_to_pool, g_glt_cache.NumPart_when_built,
                                    num_pool, send_set);
        if(device_outcome == GX_RECEIVER_COMPLETED) {
            receiver_done_on_device = 1;
        } else if(device_outcome != GX_RECEIVER_DECLINED || send_set->used != 0 || send_set->peer >= 0) {
            receiver_ok = 0;
        }
    }
    /* Host backend.  THREADED: the WALK (dominant cost) runs per received envelope into a pre-sized
     * per-envelope cand slot — each thread writes ONLY its own index (no shared write), walker
     * race-safe.  The walk records what may be a neighbour once drifted; the shared accept
     * then drifts those candidates (threaded) and accepts and emits SERIALLY (cheap), peer by
     * peer in ascending order, which is the order the send set requires. */
    if(receiver_ok && !receiver_done_on_device) {
        std::vector<std::vector<int>> per_recv_cands((size_t)(tot_r > 0 ? tot_r : 0));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
        for(long k = 0; k < tot_r; k++) {
            const struct gx_export_envelope_t *e = &recv[k];
            std::vector<int> &cvk = per_recv_cands[k];
            cvk.clear();
            mode_b_walk_from_start_nodes(e->pos, e->h, supply_mask, search_mode,
                                         spec->radius_policy, e->nodes, e->n_nodes,
                                         cvk, walker_j_reach_scale);
        }
        /* One row per envelope, pointing into its own candidate list, in ascending
         * peer order.  The shared accept drifts the distinct candidates of every row
         * together, then accepts at current positions. */
        std::vector<struct gx_candidate_row> rows;
        try {rows.reserve((size_t)tot_r);} catch(const std::bad_alloc &) {receiver_ok = 0;}
        if(receiver_ok) {
            for(int t = 0; t < NTask; t++) {
                if(t == ThisTask) continue;
                for(int r = 0; r < rc[t]; r++) {
                    const long k = (long)rd[t] + r;
                    std::vector<int> &cvk = per_recv_cands[k];
                    struct gx_candidate_row row = {&recv[k], t, cvk.data(), (int)cvk.size()};
                    rows.push_back(row);
                }
            }
            if(gx_send_set_accept_rows(send_set, rows.data(), (long)rows.size(), search_mode,
                                       spec->radius_policy, spec->j_radius_scale, spec->safety_factor) != 0) {receiver_ok = 0;}
        }
    }
    free(recv); free(rc); free(rd);
    if(receiver_ok) {
        gx_send_set_finish_peer(send_set);
        /* Every peer is finished, so every mark must be cleared again; one left set
         * would silently drop that slot for some peer on a later call. */
        if(send_set->marks_outstanding != 0) {
            printf("ERROR: ghost send set on task %d finished with %ld marks still set\n",
                   ThisTask, send_set->marks_outstanding);
            fflush(stdout);
            gizmo_request_controlled_stop(7735, "ghost_exchange: send set finished with its marker not clear",
                                          __FILE__, __LINE__, __FUNCTION__);
            receiver_ok = 0;
        }
    }
    if(!receiver_ok) {gx_send_set_abandon(send_set);}

    /* The send set grows as it is filled, so whether it could be built is known only
     * now.  Agree it here, where every rank arrives whichever backend it used, before
     * the count and payload exchanges that follow. */
    int receiver_ok_all = 0;
    MPI_Allreduce(&receiver_ok, &receiver_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if(!receiver_ok_all) {
        if(receiver_ok) {gx_send_set_abandon(send_set);}
        if(res) res->status = GX_WALK_EXPORT_RECEIVER_FAIL;
        return;
    }
    if(res) res->status = GX_WALK_EXPORT_OK;
}

/* Non-finite test that survives fast-math (isnan/isfinite may fold to false under
 * -ffinite-math-only): exponent all-ones => Inf or NaN. A non-finite query position
 * or radius makes every distance test "match" and would import the whole domain, so
 * the request-driven build below fails closed on it. */
static inline int gx_query_nonfinite(double v) {
    unsigned long long b; memcpy(&b, &v, sizeof(b));
    return (((b >> 52) & 0x7ffULL) == 0x7ffULL);
}

static ghost_exchange_result ghost_exchange_request_driven_impl(const struct ghost_exchange_spec_t *spec)
{
    if(NTask <= 1) return GHOST_EXCHANGE_COMPLETED;
    const double safety_factor = spec->safety_factor;
    const unsigned int request_mask = spec->request_type_mask;
    const unsigned int supply_mask  = spec->supply_type_mask;
    const int  search_mode = spec->search_mode;
    double t_ghost_start = my_second();
    NumPart_before_ghost = NumPart;
    N_gas_before_ghost = N_gas;
    NumGhostParticles = 0;
    g_ghost_pool_current_ti = -1;   /* see the tile impl's note */

    static int gx_call_seq_rd = 0;
    gx_call_seq_rd++;
    int this_call = gx_call_seq_rd;

    /* === Step 1: build local query list === */
    /* Two paths into the wire-format queries:
     *   (a) Caller-explicit (mech_fb migration): spec->n_queries >= 0 — caller
     *       supplied the source list (positions+h). We just copy into
     *       gx_query_t with safety_factor applied. No request_type_mask
     *       filtering — caller did that already in their isactive scan.
     *   (b) Legacy back-compat: spec->n_queries < 0 — scan ActiveParticleList,
     *       filter by spec->request_type_mask, build queries from
     *       P[i].Pos/KernelRadius. Used by ghost_exchange / ghost_exchange_hydro
     *       / ghost_exchange_hydro_oneway wrappers. */
    int n_local_queries = 0;
    struct gx_query_t *local_queries = NULL;
    int q_bad = -1;   /* first query index with a non-finite pos/h; -1 = all finite */
    if(spec->n_queries >= 0) {
        n_local_queries = spec->n_queries;
        local_queries = (struct gx_query_t *)
            malloc((size_t)(n_local_queries > 0 ? n_local_queries : 1) * sizeof(struct gx_query_t));
        for(int q = 0; q < n_local_queries; q++) {
            local_queries[q].pos[0] = spec->query_pos[q][0];
            local_queries[q].pos[1] = spec->query_pos[q][1];
            local_queries[q].pos[2] = spec->query_pos[q][2];
            local_queries[q].h    = spec->query_h[q] * safety_factor;
            local_queries[q].type = -1;   /* unspecified for caller-explicit */
            local_queries[q]._pad = 0;
            if(q_bad < 0 &&
               (gx_query_nonfinite(local_queries[q].pos[0]) || gx_query_nonfinite(local_queries[q].pos[1]) ||
                gx_query_nonfinite(local_queries[q].pos[2]) || gx_query_nonfinite(local_queries[q].h))) q_bad = q;
        }
    } else {
        for(size_t kk = 0; kk < ActiveParticleList.size(); kk++) {
            int i = ActiveParticleList[kk];
            if(i < 0 || i >= NumPart) continue;
            if(P[i].Mass <= 0) continue;
            if(!ghost_type_passes((int)P[i].Type, request_mask)) continue;
            n_local_queries++;
        }
        local_queries = (struct gx_query_t *)
            malloc((size_t)(n_local_queries > 0 ? n_local_queries : 1) * sizeof(struct gx_query_t));
        int q = 0;
        for(size_t kk = 0; kk < ActiveParticleList.size(); kk++) {
            int i = ActiveParticleList[kk];
            if(i < 0 || i >= NumPart) continue;
            if(P[i].Mass <= 0) continue;
            if(!ghost_type_passes((int)P[i].Type, request_mask)) continue;
            local_queries[q].pos[0] = P[i].Pos[0];
            local_queries[q].pos[1] = P[i].Pos[1];
            local_queries[q].pos[2] = P[i].Pos[2];
            double h = (double)P[i].KernelRadius;
            local_queries[q].h    = h * safety_factor;
            local_queries[q].type = (int)P[i].Type;
            local_queries[q]._pad = 0;
            if(q_bad < 0 &&
               (gx_query_nonfinite(local_queries[q].pos[0]) || gx_query_nonfinite(local_queries[q].pos[1]) ||
                gx_query_nonfinite(local_queries[q].pos[2]) || gx_query_nonfinite(local_queries[q].h))) q_bad = q;
            q++;
        }
    }

    /* Fail-closed: a non-finite query position or radius makes every distance test
     * "match" and would import the entire domain (a corrupt query must never turn
     * into "ship everything"). Stop loudly instead; drains collectively at the poll. */
    if(q_bad >= 0) {
        printf("ERROR: non-finite request-driven query on task %d (caller=%s q=%d pos=(%g,%g,%g) h=%g)\n",
               ThisTask, (spec->caller_name ? spec->caller_name : "?"), q_bad,
               local_queries[q_bad].pos[0], local_queries[q_bad].pos[1], local_queries[q_bad].pos[2], local_queries[q_bad].h);
        gizmo_request_controlled_stop(7708, "ghost_exchange (request-driven): non-finite query position or radius", __FILE__, __LINE__, __FUNCTION__);
    }
    gizmo_exit_bad_stop_if_requested("ghost_exchange:query_finite");

        /* === Step 3: per-rank, walk local BVH against each remote rank's queries ===
     *
     * Build a host-side SFC tile + BVH index over the supply-mask-filtered pool.
     * For each remote query, walk the BVH (O(log N + matches) per query, vs the
     * O(N) brute-force scan that had been a tiny-N proof-of-concept). Per-particle
     * acceptance at leaves uses the EXACT predicate (ONEWAY: r²<h_q²; SYMMETRIC:
     * r²<max(h_q,h_j)²) — same predicate the kernel applies later. */

    /* SHARED-TREE build: build over GHOST_TYPE_ALL (all types with mass>0), not
     * just the caller's supply_mask. The walker's per-type hmax filter +
     * per-particle leaf Type-vs-supply_mask check delivers the same imports
     * as the old per-supply-mask build.
     *
     * Bucket 3 (SIDX overlay): the local tree (tiles, pool, bvh,
     * compact_xyzh, pool_types) is cached across calls within a step.
     * Membership is keyed on the supply-identity epoch, which every event that
     * changes who is in the pool already bumps, so within a step the pool is
     * stable and 2nd..Nth calls skip this whole stanza. */

    int *h_pool = NULL;
    int num_pool = 0;
    int from_cache = 0;
    /* All particle types are eligible as ghost sources. */
    unsigned int desired_pool_mask = GHOST_TYPE_ALL;
    /* Mask is compared for EXACT equality, not coverage. Reuse across a narrowed
     * mask is unproven here — a narrower request would also make in-place Type
     * changes membership-relevant, which the epoch does not track — so anything
     * other than the all-types pool falls through to a full rebuild. */
    const int mask_reusable = (desired_pool_mask == GHOST_TYPE_ALL);
    const int identity_valid = (g_glt_cache.valid
                       && g_glt_cache.pool && g_glt_cache.j_to_pool && g_glt_cache.mark
                       && mask_reusable
                       && g_glt_cache.eligible_type_mask_when_built == desired_pool_mask
                       && g_glt_cache.NumPart_when_built == NumPart
                       && g_glt_cache.identity_epoch_when_built == g_supply_identity_epoch);
    int cache_match = identity_valid;

    if(cache_match) {
        h_pool         = g_glt_cache.pool;
        num_pool       = g_glt_cache.num_pool;
        from_cache = 1;
        g_glt_cache_hits++;
    } else {
        /* RANK-LOCAL BRANCH — NO MPI CALLS IN HERE.  cache_match keys on rank-local
         * NumPart, so ranks enter this independently; any collective placed here
         * deadlocks as soon as one rank rebuilds and another does not. */
        g_glt_cache_misses++;
        /* Call-local on purpose: a file-static would latch across calls, and every
         * later rebuild would then skip fill+install while still publishing the
         * pool pointers below -- a NULL walk. */
        int cache_alloc_failed = 0;
        /* Free any stale entry before rebuild (cache key changed). */
        if(g_glt_cache.valid) glt_cache_free();

        /* Fresh build via the existing mymalloc path; the result is copied into
         * malloc-backed cache buffers so it can outlive this function frame
         * without violating mymalloc LIFO ordering. */
        int *tmp_pool = NULL;
        int tmp_num_pool = build_sfc_supply_pool(P, NumPart, (int)desired_pool_mask, &tmp_pool);

        /* Allocate persistent cache buffers + copy. */
        size_t sz_pool = (size_t)(tmp_num_pool > 0 ? tmp_num_pool : 1) * sizeof(int);
        int   *c_pool  = (int *) malloc(sz_pool);
        /* Reverse map j -> pool_pos (-1 if j is not in this build's pool). */
        size_t sz_jtop = (size_t)(NumPart > 0 ? NumPart : 1) * sizeof(int);
        int   *c_jtop  = (int *) malloc(sz_jtop);
        /* The send set's dedup marker, one byte per pool slot.  Zeroed here, once per
         * pool, and kept zero between calls by the send set itself. */
        unsigned char *c_mark = (unsigned char *) calloc((size_t)(tmp_num_pool > 0 ? tmp_num_pool : 1), 1);
        /* An allocation failure here would otherwise be a segfault: the buffers are
         * written unconditionally just below, and this producer is now the only
         * supplier, so there is nothing to fall back to.  The request is RANK-LOCAL
         * on purpose -- this branch is entered per-rank (cache_match keys on
         * rank-local NumPart), so a collective here would deadlock whenever ranks
         * disagree about rebuilding.  The matching drain runs just past the branch,
         * where every rank converges and before anything reads the pool. */
        if(!c_pool || !c_jtop || !c_mark) {
            printf("ERROR: supply-cache allocation failed on task %d (num_pool=%d NumPart=%d)\n",
                   ThisTask, tmp_num_pool, NumPart);
            fflush(stdout);
            free(c_pool); free(c_jtop); free(c_mark);
            c_pool = NULL; c_jtop = NULL; c_mark = NULL;
            if(tmp_pool) myfree(tmp_pool);
            gizmo_request_controlled_stop(7724, "ghost_exchange: supply-cache allocation failed",
                                          __FILE__, __LINE__, __FUNCTION__);
            cache_alloc_failed = 1;
        }
        if(!cache_alloc_failed) {
        for(int j = 0; j < NumPart; j++) c_jtop[j] = -1;
        if(tmp_num_pool > 0) memcpy(c_pool, tmp_pool, (size_t)tmp_num_pool * sizeof(int));
        for(int p = 0; p < tmp_num_pool; p++) {
            int j = tmp_pool[p];
            if(j >= 0 && j < NumPart) c_jtop[j] = p;
        }

        /* Free the mymalloc temp. */
        if(tmp_pool) myfree(tmp_pool);

        g_glt_cache.pool = c_pool;
        g_glt_cache.j_to_pool = c_jtop;
        g_glt_cache.mark = c_mark;
        g_glt_cache.num_pool = tmp_num_pool;
        g_glt_cache.NumPart_when_built = NumPart;
        g_glt_cache.identity_epoch_when_built = g_supply_identity_epoch;
        g_glt_cache.eligible_type_mask_when_built = desired_pool_mask;
        g_glt_cache.valid = 1;
        }   /* end fill+install (skipped when the cache allocation failed) */

        h_pool = c_pool; num_pool = tmp_num_pool;
    }
    /* Both branches converge here, so this poll is reached by every rank: it drains a
     * rank-local supply-cache allocation failure into an all-rank controlled stop
     * BEFORE anything below dereferences the pool. */
    gizmo_exit_bad_stop_if_requested("ghost_exchange:supply_cache_alloc");


    /* The routed walk-export producer (sender export + bounded receiver walk,
     * collective-safe) fills the send set, reading the same g_glt_cache snapshot as
     * the rest of this call.  Membership comes from the SSOT accept
     * (gx_pair_accept_wrap_and_test), so the only way the set can differ from a full
     * walk is routing COVERAGE, which the per-type node band establishes for every
     * radius policy.
     *
     * It is the only producer.  If it could not build the set there is no substrate
     * left that is known to be correct: the walks that once served as fallbacks read
     * cached geometry which has been measured producing wrong densities where that
     * geometry went stale, so reviving one would trade a visible failure for a silent
     * one.  Stop instead.  The producer status is rank-uniform, so all ranks stop
     * together.  The dominant failure mode is allocation under memory pressure; the
     * recovery that fits it is a retry at reduced import padding inside this
     * producer, which does not exist yet -- until it does, the honest outcome is this
     * stop. */
    struct ghost_send_set *send_set = &g_send_set;
    struct gx_walk_export_result walk_export_res;
    gx_walk_export_discover(spec, local_queries, n_local_queries,
                            num_pool, supply_mask, search_mode,
                            send_set, &walk_export_res);
    if(walk_export_res.status != GX_WALK_EXPORT_OK) {
        if(ThisTask == 0) {
            printf("[ghost_exchange call=%d caller=%s: walk-export producer status %d, no send set built]\n",
                   this_call, (spec->caller_name ? spec->caller_name : "?"), walk_export_res.status);
            fflush(stdout);
        }
        gizmo_request_controlled_stop(7723,
            "ghost_exchange: walk-export producer unavailable and no correctness-proven fallback exists",
            __FILE__, __LINE__, __FUNCTION__);
        gizmo_exit_bad_stop_if_requested("ghost_exchange:walk_export_unavailable");
    }

    /* === Step 4: per-peer counts === */
    int *send_count = (int *) mymalloc("gx_rd_sc", NTask * sizeof(int));
    int *recv_count = (int *) mymalloc("gx_rd_rc", NTask * sizeof(int));
    int *send_disp  = (int *) mymalloc("gx_rd_sd", NTask * sizeof(int));
    int *recv_disp  = (int *) mymalloc("gx_rd_rd", NTask * sizeof(int));
    for(int t = 0; t < NTask; t++) { send_count[t] = send_set->count[t]; recv_count[t] = 0; }
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);
    /* CHECKED int64 totals + prefix displacements (same rationale as the tile
     * impl). Request-driven is the last-resort Mode-A discovery — there is NO
     * further fallback — so both a representation overflow and a particle-slot
     * overflow fail HONESTLY via the collective controlled-stop poll, never a
     * silent int wrap or OOB append. */
    long long total_send_ll = 0, total_recv_ll = 0;
    int count_range_ok = 1;
    {
        long long sdisp = 0, rdisp = 0;
        for(int t = 0; t < NTask; t++) {
            if(sdisp <= INT_MAX) send_disp[t] = (int)sdisp; else { send_disp[t] = 0; count_range_ok = 0; }
            if(rdisp <= INT_MAX) recv_disp[t] = (int)rdisp; else { recv_disp[t] = 0; count_range_ok = 0; }
            sdisp += send_count[t];
            rdisp += recv_count[t];
        }
        total_send_ll = sdisp; total_recv_ll = rdisp;
        if(total_send_ll > INT_MAX || total_recv_ll > INT_MAX) count_range_ok = 0;
    }
    int total_send = count_range_ok ? (int)total_send_ll : 0;
    int total_recv = count_range_ok ? (int)total_recv_ll : 0;

    /* Ghosts cannot be refused: this is the last-resort Mode-A discovery and there is nothing
     * further to fall back to, so an import that does not fit the current capacity raises it.
     * The capacity is the size of one memory block, and growing it is legal here with the gravity
     * tree standing and no ghost yet written -- that is what makes it movable mid-step. The need is
     * EXACT, since request-driven discovery already counted it, so nothing is added on top: a
     * capacity only ever rises, and a margin would raise the run's footprint permanently on the
     * strength of one crowded step. Whether the memory exists is the allocator's answer rather than
     * a prediction of it; on failure the resize requests a controlled stop and leaves every array at
     * a capacity the advertised one is backed by. Purely local -- no communication, and nothing on
     * the common path but one comparison. Runs BEFORE the pack so the capacity is settled before
     * anything is written at &P[NumPart]; the guard below still decides whether the append happens. */
    {
        const long long required = (long long) NumPart + total_recv_ll;
        if(count_range_ok && required > (long long) All.MaxPart && required <= (long long) INT_MAX) {
            (void) resize_particle_storage((int) required);
        }
    }

    /* Check space (mirrors legacy guard).  Request-driven is the last-resort
     * Mode-A discovery — there is NO further fallback — so a count/displacement
     * overflow of the int MPI transport range, or ghosts that would not fit
     * P[]/CellP[], fail HONESTLY via the collective controlled-stop poll below.
     * This stays the ONE predicate that decides whether the append may happen: if
     * the growth above was refused or failed, it is this guard that stops the run. */
    if(send_set->used != total_send_ll) {
        /* Step 5 hands the send set's slots over as the send list, one per send
         * position; a different length would pack past it or leave positions unset. */
        printf("ERROR: request-driven ghost exchange on task %d: the send set holds %ld slots but the per-peer "
               "counts send %lld.\n", ThisTask, send_set->used, total_send_ll);
        gizmo_request_controlled_stop(7736, "ghost_exchange (request-driven): send set length disagrees with its counts",
                                      __FILE__, __LINE__, __FUNCTION__);
    } else if(!count_range_ok) {
        printf("ERROR: request-driven ghost exchange counts exceed int MPI transport range on task %d.\n", ThisTask);
        gizmo_request_controlled_stop(7703, "ghost_exchange (request-driven): ghost count/displacement exceeds int MPI transport range", __FILE__, __LINE__, __FUNCTION__);
    } else if(!ghost_particle_slots_fit((long long)NumPart + total_recv_ll)) {
        printf("ERROR: request-driven ghost exchange needs %d ghosts on task %d, only %d free.\n",
               total_recv, ThisTask, All.MaxPart - NumPart);
        gizmo_request_controlled_stop(7702, "ghost_exchange (request-driven): ghost append would exceed the particle capacity and the capacity could not be raised to hold it (add ranks/nodes, or reduce ghost-import demand)", __FILE__, __LINE__, __FUNCTION__);
    }
    /* Per-rank capacity check above is asymmetric; drain it at this all-rank poll
     * BEFORE Step 5, so no rank appends ghosts past MaxPart (OOB) or desyncs the
     * collective pack/exchange. Every rank reaches this unconditionally. */
    gizmo_exit_bad_stop_if_requested("ghost_exchange:capacity_rd");

    /* === Step 5: work out which local particle fills each send position ===
     * The send set already lists each peer's pool slots in send order -- destination
     * rank ascending, then pool index ascending -- in one contiguous run per peer,
     * starting where send_disp[t] says, so mapping it through the pool in place IS
     * the send list.  It is trimmed to its length first: it outlives this call as
     * the refresh provenance, and growth by doubling can leave it up to twice that. */
    int *send_home_idx = send_set->slots;
    if(send_set->capacity > (long)total_send) {
        int *fitted = (int *) realloc(send_home_idx, (size_t)(total_send > 0 ? total_send : 1) * sizeof(int));
        if(fitted) {send_home_idx = fitted;}
    }
    send_set->slots = NULL;
    send_set->capacity = send_set->used = 0;
    for(int off = 0; off < total_send; off++) {send_home_idx[off] = h_pool[send_home_idx[off]];}
    /* Advance the exported particles before copying them, so nothing goes on the
     * wire behind the time its receiver will read it at. */
    const int send_list_current = gx_certify_send_list_current(send_home_idx, total_send, All.Ti_Current);

    /* === Step 6: pack and Alltoallv particles + cells, then home_idx === */
    gx_pack_and_forward_particle_exchange(send_home_idx, send_count, send_disp,
                                          &P[NumPart], &CellP[NumPart], recv_count, recv_disp);

    /* Update counts now so home_idx receive can land at &P[NumPart_before_ghost+...] */
    NumGhostParticles = total_recv;
    NumPart += total_recv;

    /* Mark dirty for compact_xyzh refresh (same as legacy). */
    if(NumGhostParticles > 0) {
        gpu_compact_xyzh_mark_h_dirty_range(NumPart_before_ghost, NumPart);
    }
    /* SIDX lifecycle notify: see comment in tile-overlap impl. Unconditional. */
    /* Every rank advanced its outgoing slots to this same All.Ti_Current above,
       and the exchange is collective, so the pool just installed is current at
       that time -- but only if that advance actually happened. */
    if(send_list_current == 0) {g_ghost_pool_current_ti = All.Ti_Current;}
    gpu_sidx_notify_ghost_imported(NumPart_before_ghost, NumGhostParticles);

    /* Home-index exchange + provenance maps. */
    int *recv_home_idx = (int *) malloc((total_recv > 0 ? total_recv : 1) * sizeof(int));
    gizmo_mpi_alltoallv_typed(send_home_idx, send_count, send_disp,
                              recv_home_idx, recv_count, recv_disp,
                              sizeof(int), MPI_COMM_WORLD);
    ghost_home_rank_map = (int *) malloc((total_recv > 0 ? total_recv : 1) * sizeof(int));
    ghost_home_index_map = recv_home_idx;
    for(int t = 0; t < NTask; t++) {
        for(int g = 0; g < recv_count[t]; g++) {
            ghost_home_rank_map[recv_disp[t] + g] = t;
        }
    }
    /* Preserve comm maps for reverse Alltoallv (ghost writeback). */
    ghost_wb_recv_count = (int *) malloc(NTask * sizeof(int));
    ghost_wb_recv_disp  = (int *) malloc(NTask * sizeof(int));
    ghost_wb_send_count = (int *) malloc(NTask * sizeof(int));
    ghost_wb_send_disp  = (int *) malloc(NTask * sizeof(int));
    memcpy(ghost_wb_recv_count, recv_count, NTask * sizeof(int));
    memcpy(ghost_wb_recv_disp,  recv_disp,  NTask * sizeof(int));
    memcpy(ghost_wb_send_count, send_count, NTask * sizeof(int));
    memcpy(ghost_wb_send_disp,  send_disp,  NTask * sizeof(int));

    /* Preserve send-side provenance for ghost_refresh_values() (take ownership
       of send_home_idx; the free() below then no-ops on NULL). */
    if(ghost_send_home_idx) free(ghost_send_home_idx);
    ghost_send_home_idx   = send_home_idx;
    ghost_send_home_count = total_send;
    send_home_idx         = NULL;
    g_ghost_provenance_epoch++;

    double t_ghost_total = timediff(t_ghost_start, my_second());


    if(ThisTask == 0) {
        /* Every count below belongs to rank 0 alone. There is no global total to hand at this
           point and reducing one for a log line would cost a collective on a hot path, so the
           line says whose numbers these are: a quiet line means rank 0 asked for nothing, NOT
           that the exchange did no work. The discovery walk and its per-peer bookkeeping run
           over the supply pool whatever any one rank asked for.
           The all-rank query total exists only on the broadcast path, where the queries are
           gathered; it is reported when it is real and omitted when it is not, rather than
           printed as a placeholder. */
        PRINT_STATUS("Ghost exchange (request-driven, %s, %s): rank 0 holds %d local + %d ghost "
                     "from %d queries; supply pool %d  [%.4f s]",
                     (spec->caller_name ? spec->caller_name : "?"),
                     (search_mode == NGB_SEARCH_ONEWAY ? "ONEWAY" : "SYMMETRIC"),
                     NumPart_before_ghost, NumGhostParticles,
                     n_local_queries, num_pool, t_ghost_total);
    }

    /* Diagnostic: ghost composition + import-waste ratio (should be ~0% for
     * the request-driven path by construction since per-particle accept ran
     * before pack — provides direct A/B vs the legacy tile-overlap waste). */

    /* Cleanup local. mymalloc requires LIFO free order. Tile/BVH/pool/
     * compact_xyzh/pool_types are now owned by g_glt_cache (malloc-backed)
     * and outlive this frame; do NOT free them here. They're freed at
     * cache invalidation (drift / domain_decomp hooks) via glt_cache_free. */
    myfree(recv_disp);
    myfree(send_disp);
    myfree(recv_count);
    myfree(send_count);
    free(send_home_idx);
    free(local_queries);
    (void)from_cache;
    return GHOST_EXCHANGE_COMPLETED;
}

/* Public wrappers — each fills a spec, calls the single _impl. New callers
 * add a wrapper line; do not duplicate logic.
 *
 * radius_policy + j_radius_scale on the spec are part of the SSOT supply-side
 * contract.  Legacy non-runner wrappers explicitly pass
 * MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES + 1.0 to preserve their pre-policy
 * behavior byte-for-byte (raw P[j].KernelRadius * safety_factor as the
 * supply-side reach).  Runner Mode A passes Spec::radius_policy +
 * nlr_spec_symmetric_j_radius_scale<Spec>() via gizmo_request_filtered_ghost_import_fresh
 * — see ghost_symlist_lifecycle.h. */
void ghost_exchange(double safety_factor)
{
    /* Reach is P[j].KernelRadius (the legacy all-types policy) for every type,
     * times safety_factor. The per-type opener band is seeded from the
     * conservative source union (force_hmax_per_type_particle_radius) and
     * exchanged cross-rank on the nodes the export walk descends, so it bounds
     * this reach per type. */
    struct ghost_exchange_spec_t sp = {GHOST_TYPE_ALL, GHOST_TYPE_ALL, NGB_SEARCH_SYMMETRIC, safety_factor, "all_types", -1, NULL, NULL,
                                       MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES, 1.0};
    ghost_exchange_impl(&sp);
}
void ghost_exchange_hydro(double safety_factor)
{
    /* Gas-only supply at the legacy all-types kernel radius, which the per-type
     * node band is built from, so the band bounds this spec's reach. */
    struct ghost_exchange_spec_t sp = {GHOST_TYPE_0, GHOST_TYPE_0, NGB_SEARCH_SYMMETRIC, safety_factor, "hydro_symmetric", -1, NULL, NULL,
                                       MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES, 1.0};
    ghost_exchange_impl(&sp);
}
void ghost_exchange_hydro_oneway(double safety_factor)
{
    /* ONEWAY opens on the query's own radius; the supply side bounds nothing here. */
    struct ghost_exchange_spec_t sp = {GHOST_TYPE_0, GHOST_TYPE_0, NGB_SEARCH_ONEWAY, safety_factor, "hydro_oneway", -1, NULL, NULL,
                                       MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES, 1.0};
    ghost_exchange_impl(&sp);
}


/*!
 * \brief Remove ghost particles after neighbor operations complete.
 *
 * Resets NumPart and N_gas to pre-exchange values. Must be called after
 * all neighbor loops (density, gradients, hydro force) that use ghosts.
 */
void ghost_exchange_cleanup(void)
{
    if(NumPart_before_ghost < 0) return;
    /* Ghost slots are about to leave scope (NumPart shrinks back to local).
     * No dirty-state scrubbing is done here. Marks are per-cache: one landing
     * outside a cache's registered index range is dropped at mark time, and a
     * cache that WAS registered over these ghost slots is freed -- handle
     * unregistered with it -- by the particle-count change on its next build.
     * Either way no stale ghost-slot bit reaches compact_h_refresh. New ghost
     * slots are marked dirty at import time (mark_h_dirty_range above), so
     * symmetric h-reads on ghosts stay fresh. */
    /* SIDX lifecycle notify BEFORE NumPart shrinks. Called whether or not
     * NumGhostParticles>0 — a cleanup from the no-ghost-imported state is
     * a valid signal that bumps the epoch. */
    gpu_sidx_notify_ghost_cleanup();
    if(NumGhostParticles > GhostEpochHighWater) {GhostEpochHighWater = NumGhostParticles;}
    NumPart = NumPart_before_ghost;
    N_gas = N_gas_before_ghost;
    NumGhostParticles = 0;
    NumPart_before_ghost = -1;
    g_ghost_pool_current_ti = -1;   /* the pool the stamp described is gone */
    /* Free ghost provenance map */
    if(ghost_home_rank_map)  { free(ghost_home_rank_map);  ghost_home_rank_map = NULL; }
    if(ghost_home_index_map) { free(ghost_home_index_map); ghost_home_index_map = NULL; }
    if(ghost_wb_recv_count)  { free(ghost_wb_recv_count);  ghost_wb_recv_count = NULL; }
    if(ghost_wb_recv_disp)   { free(ghost_wb_recv_disp);   ghost_wb_recv_disp = NULL; }
    if(ghost_wb_send_count)  { free(ghost_wb_send_count);  ghost_wb_send_count = NULL; }
    if(ghost_wb_send_disp)   { free(ghost_wb_send_disp);   ghost_wb_send_disp = NULL; }
    if(ghost_send_home_idx)  { free(ghost_send_home_idx);  ghost_send_home_idx = NULL; }
    ghost_send_home_count = 0;
    /* g_ghost_provenance_epoch is a monotonic stamp — NOT reset here. */
}

/* Make refreshed host ghost values visible to the device. With unified-memory
   particles (P/CellP in Kokkos SharedSpace) this is a no-op: a host write to a
   ghost slot is coherent to the next kernel, and compact_xyzh caches Pos+h only
   (unchanged by a value refresh). This helper is ALWAYS in the refresh call path
   (never elided) so a backend with explicit host/device particle buffers (no
   unified-memory coherence) has ONE mandatory place to add an explicit
   host->device ghost copy. Parity with import: import marks compact_xyzh h-dirty
   + notifies SIDX; a value refresh changes neither Pos/h nor the slot set, so it
   does neither — but any explicit device copy import gains MUST be mirrored here. */
static inline void ghost_refresh_make_device_visible(int ghost_base, int ghost_count)
{
    (void)ghost_base; (void)ghost_count;
}

int ghost_refresh_values(void)
{
    /* Fail-closed guards (production callers fall back to full cleanup+reimport).
       A non-NULL ghost_send_home_idx already implies an import happened with no
       cleanup since (cleanup NULLs it). */
    if(NTask <= 1)               return GHOST_REFRESH_SKIP_SERIAL;
    if(NumPart_before_ghost < 0) return GHOST_REFRESH_FAIL_NO_POOL;
    if(!ghost_send_home_idx || !ghost_wb_send_count || !ghost_wb_send_disp ||
       !ghost_wb_recv_count || !ghost_wb_recv_disp)
                                 return GHOST_REFRESH_FAIL_NO_PROVENANCE;
    /* Live-pool consistency: current pool counts must match the preserved
       provenance (topology unchanged; no intervening reimport with different
       totals). */
    long long send_tot = 0, recv_tot = 0;
    for(int t = 0; t < NTask; t++) { send_tot += ghost_wb_send_count[t]; recv_tot += ghost_wb_recv_count[t]; }
    if(NumPart != NumPart_before_ghost + NumGhostParticles) return GHOST_REFRESH_FAIL_POOL_MUTATED;
    if(recv_tot != (long long)NumGhostParticles)            return GHOST_REFRESH_FAIL_POOL_MUTATED;
    if(send_tot != (long long)ghost_send_home_count)        return GHOST_REFRESH_FAIL_POOL_MUTATED;

    int ns = ghost_send_home_count;
    /* Re-pack current owner values via the SAME transport import uses, in the SAME
       send order (ghost_send_home_idx) that produced this pool. */
    const int send_list_current = gx_certify_send_list_current(ghost_send_home_idx, ns, All.Ti_Current);
    /* Replay ONLY the forward transport, overwriting the EXISTING ghost slots at
       [NumPart_before_ghost, NumPart). Slots land at identical offsets by
       construction, so ghost_home_rank/index maps + any built CSR stay valid. */
    gx_pack_and_forward_particle_exchange(ghost_send_home_idx,
                                          ghost_wb_send_count, ghost_wb_send_disp,
                                          &P[NumPart_before_ghost], &CellP[NumPart_before_ghost],
                                          ghost_wb_recv_count, ghost_wb_recv_disp);
    /* The refresh re-packed at the current time, so the pool is current to it --
       re-stamp, or the guarantee would lag the pool it describes. Withheld when the
       advance did not complete, for the same reason the import path withholds it. */
    if(send_list_current == 0) {g_ghost_pool_current_ti = All.Ti_Current;}
    ghost_refresh_make_device_visible(NumPart_before_ghost, NumGhostParticles);
    return GHOST_REFRESH_OK;
}

/* Import-epoch accessor (see g_ghost_provenance_epoch): read by the hydro
   corridor, which fast-paths a value-refresh only when the live pool's epoch
   matches the one its published CSR was built from. */
unsigned long long ghost_provenance_epoch(void) { return g_ghost_provenance_epoch; }

/* Accessors for ghost provenance data — used by ghost_writeback.cc */
/* True iff a ghost import is live (pool materialized, between import and cleanup).
   Distinguishes "live pool with zero ghosts" from "no pool" — callers must not
   infer liveness from ghost_get_num_ghosts(), which returns 0 in both states.
   Used by the neighbor-loop runner to enforce the caller-owned-pool contract
   for external-CSR consumers (see neighbor_loop_runner.h). */
int ghost_pool_is_live(void) { return (NumPart_before_ghost >= 0) ? 1 : 0; }
int ghost_get_num_ghosts(void) { return NumGhostParticles; }
int ghost_get_epoch_high_water(void)
{
    return (GhostPreviousEpochHighWater > GhostEpochHighWater) ? GhostPreviousEpochHighWater
                                                               : GhostEpochHighWater;
}
void ghost_reset_epoch_high_water(void)
{
    GhostPreviousEpochHighWater = GhostEpochHighWater;
    GhostEpochHighWater = 0;
}
int ghost_get_num_local(void)  { return (NumPart_before_ghost >= 0) ? NumPart_before_ghost : NumPart; }
int *ghost_get_home_rank(void)  { return ghost_home_rank_map; }
int *ghost_get_home_index(void) { return ghost_home_index_map; }
int *ghost_get_wb_recv_count(void) { return ghost_wb_recv_count; }
int *ghost_get_wb_recv_disp(void)  { return ghost_wb_recv_disp; }
int *ghost_get_wb_send_count(void) { return ghost_wb_send_count; }
int *ghost_get_wb_send_disp(void)  { return ghost_wb_send_disp; }
