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
#include <vector>
#include "../declarations/allvars.h"
#include "../declarations/lifecycle_counters.h"
#include "../core/proto.h"
#include "../system/mpi_alltoallv_typed.h"
#include "../core/step_phases.h"   /* gizmo_verbose_diag() */
#include "gpu_neighbor_list.h" /* gpu_compact_xyzh_mark_h_dirty_range */
#include "sfc_tiles.h"           /* build_sfc_tiles, build_tile_bvh, sfc_tile_t, tile_bvh_node_t */
#include "neighbor_list.h"       /* NGB_SEARCH_ONEWAY, NGB_SEARCH_SYMMETRIC */
#include "ghost_exchange_functions.h" /* gx_pair_accept (shared geometric accept) */
#include "ghost_writeback.h"     /* ghost_get_num_local (bounded fine-tree walk) */
#include "ghost_exchange_spec.h"
#include "topleaf_router.h"      /* gx_supply_pool_view + ghost_exchange_supply_pool_view (band builder) */
#include "gpu_fine_sidecar.h"    /* L4 S2a device fine-tree sidecar (upload/free/valid/readback) */
#include "../gravity/gpu_gravity_tree.h"  /* gpu_gravity_soa_ensure_drifted (S2b-1 drift stamp) */
#include "mode_b_local_walker.h"  /* mode_b_walk_and_export / mode_b_walk_from_start_nodes */
#ifdef _OPENMP
#include <omp.h>                 /* threaded sender export + receiver walk below */
#endif

/* Defined in neighbor_loop_runner.cc. Declared here rather than including the
 * runner header, which is a heavy template TU. Lets this file's diagnostic output
 * share ONE gate with PHASE0_NLR, so Mode-A and Mode-B runs are collectable in the
 * same configuration and remain comparable. */
bool gizmo_nlr_phase0_diag_enabled(void);

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
static int PreviousGhostCount = 0; /* ghost count from the most recent completed exchange, for domain headroom */

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

/* Persistent local-tree cache (SIDX overlay) for the
 * request-driven ghost exchange path. Within a step, the local pool of
 * particles [0..NumPart_local) is stable across multiple ghost_exchange
 * calls (3-5 calls/step typical). Building tiles+BVH+compact_xyzh once per
 * call costs ~0.19s (gas) / ~0.65s (all-types) on the fire_m11i 6.2M/9.5M
 * pool. Caching them across calls saves N-1 of those builds per step.
 *
 * Invalidation is wired to the same hooks as gpu_step_sidx_invalidate_*:
 *   - run.cc post-drift  -> ghost_exchange_local_tree_invalidate_drift()
 *   - run.cc post-decomp -> ghost_exchange_local_tree_invalidate_full()
 * Drift/h updates mark the cache for exact refit from P[] on the next hit;
 * domain decomposition and pool/ordering changes fully free it.
 *
 * Cache key (NumPart, safety_factor, eligible pool mask) is also checked
 * at use time as a defensive cross-check; pool/ordering changes force a
 * rebuild. Drift/h changes mark the cache for an exact refit from P[].
 *
 * Memory footprint: ~165MB (Type-0/cell pool, 6.2M pool) or ~250MB (all-types,
 * 9.5M pool) per rank. Tolerable on Vista host. Allocated via plain
 * malloc/free to avoid mymalloc-stack LIFO violation when the cache
 * outlives the function frame. */
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
    integertime Ti_when_built;
    double safety_factor_when_built;
    unsigned int eligible_type_mask_when_built;
    int needs_refit;
    /* What this cache entry actually holds (GX_POOL_IDENTITY / GX_POOL_GEOMETRY).
     * Geometry is skipped for calls whose producer never walks it, so `valid`
     * alone does not imply tiles/bvh/compact_xyzh are present. */
    unsigned int caps;
    int ntiles;
    int num_pool;
    int bvh_nnodes;
    int bvh_root;
    sfc_tile_t *tiles;            /* [ntiles] malloc */
    int *pool;                     /* [num_pool] malloc */
    tile_bvh_node_t *bvh;          /* [bvh_nnodes] malloc */
    float *compact_xyzh;           /* [num_pool*4] malloc; h field = supply-side policy * j_scale * safety baked in */
    int *pool_types;               /* [num_pool] malloc */
    int *j_to_pool;                /* [NumPart_when_built] malloc, j -> pool_pos or -1 */
    /* SSOT supply-side reach contract. These four fields, together with the
     * (NumPart, safety, eligible_pool_mask) triple above, form the cache-key
     * invariant: any mismatch on any of them forces a full rebuild (NOT a
     * refit).  Ti and the needs_refit dirty bit continue to trigger
     * glt_cache_refit_from_particles() instead of a rebuild — Ti is NOT a
     * rebuild key: rebuilding every step would be exactly the work
     * explosion this design is supposed to prevent. */
    mode_b_radius_policy_t radius_policy_when_built;
    double j_radius_scale_when_built;
};
/* Designated initialisers: the positional form silently mis-assigns whenever a
 * field is added to the struct above. */
static struct ghost_local_tree_cache_t g_glt_cache = {
    .valid = 0,
    .NumPart_when_built = -1,
    .identity_epoch_when_built = -1,
    .Ti_when_built = -1,
    .safety_factor_when_built = 0.0,
    .eligible_type_mask_when_built = GHOST_TYPE_ALL,
    .needs_refit = 0,
    .caps = 0u,
    .ntiles = 0,
    .num_pool = 0,
    .bvh_nnodes = 0,
    .bvh_root = 0,
    .tiles = NULL,
    .pool = NULL,
    .bvh = NULL,
    .compact_xyzh = NULL,
    .pool_types = NULL,
    .j_to_pool = NULL,
    .radius_policy_when_built = MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES,
    .j_radius_scale_when_built = 1.0,
};

/* gx_policy_scaled_h — SSOT supply-side reach inside ghost_exchange.cc.
 *
 * Single helper used at every site that writes a per-particle supply h into
 * the ghost local-tree cache: fresh build (compact + tile bands via
 * build_sfc_tiles' scale_factor pathway), refit (glt_recompute_tile_), and
 * any future site that needs the supply-side reach.  The output is identical
 * to what build_sfc_tiles aggregates internally when called with
 *   (radius_policy, scale_factor = j_radius_scale * safety_factor)
 * so leaf compact h and BVH band hmax see the exact same supply-side reach
 * for every particle — closing the leaf-vs-band scale gap that would
 * otherwise let the BVH prune pairs the leaf would have accepted.
 *
 * Compile-flag gating on AGS_KernelRadius lives in nlr_radius_policy.h's
 * wrapper; never replicate the `#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE`
 * inside this TU. */
static inline double gx_policy_scaled_h(int j,
                                        mode_b_radius_policy_t radius_policy,
                                        double j_radius_scale,
                                        double safety_factor)
{
    return nlr_particle_symmetric_radius(P[j], radius_policy)
           * j_radius_scale * safety_factor;
}

/* Bucket 3 narrow-refit machinery. Mirrors gpu_neighbor_list.cc's g_dirty_list
 * pattern for the GPU compact_xyzh refresh — same call sites populate both, but
 * different consumers (host ghost_exchange cache vs device GPU NL builder), so
 * the two lists have independent lifecycles.
 *
 *  - g_glt_dirty_all = true  : full refit needed (drift, fresh build seed, or
 *    list overflow promotion). Default = true so first refit is full.
 *  - g_glt_dirty_all = false : refresh only the indices in g_glt_dirty_list.
 *    Indices outside the cache pool (j_to_pool[j] == -1) are skipped cleanly.
 *
 * Promote-to-all threshold matches GPU side (1M indices). When the list grows
 * past that, narrow refit costs ~the same as full, so flip to dirty_all. */
static const int G_GLT_DIRTY_PROMOTE_THRESHOLD = 1 << 20; /* 1M indices */
static bool g_glt_dirty_all = true;
static std::vector<int> g_glt_dirty_list;
static inline void g_glt_dirty_clear_(void)
{
    g_glt_dirty_all = false;
    g_glt_dirty_list.clear();
}
static inline void g_glt_dirty_mark_all_(void)
{
    g_glt_dirty_all = true;
    g_glt_dirty_list.clear();
}

/* Diagnostic counters. */
static long g_glt_cache_hits = 0;
static long g_glt_cache_misses = 0;
static long g_glt_cache_refits = 0;
static long g_glt_cache_narrow_refits = 0;

static void glt_cache_free(void)
{
    if(g_glt_cache.tiles)        { free(g_glt_cache.tiles);        g_glt_cache.tiles = NULL; }
    if(g_glt_cache.pool)         { free(g_glt_cache.pool);         g_glt_cache.pool = NULL; }
    if(g_glt_cache.bvh)          { free(g_glt_cache.bvh);          g_glt_cache.bvh = NULL; }
    if(g_glt_cache.compact_xyzh) { free(g_glt_cache.compact_xyzh); g_glt_cache.compact_xyzh = NULL; }
    if(g_glt_cache.pool_types)   { free(g_glt_cache.pool_types);   g_glt_cache.pool_types = NULL; }
    if(g_glt_cache.j_to_pool)    { free(g_glt_cache.j_to_pool);    g_glt_cache.j_to_pool = NULL; }
    g_glt_cache.valid = 0;
    g_glt_cache.NumPart_when_built = -1;
    g_glt_cache.identity_epoch_when_built = -1;
    g_glt_cache.Ti_when_built = -1;
    g_glt_cache.safety_factor_when_built = 0.0;
    g_glt_cache.eligible_type_mask_when_built = GHOST_TYPE_ALL;
    g_glt_cache.radius_policy_when_built = MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES;
    g_glt_cache.j_radius_scale_when_built = 1.0;
    g_glt_cache.needs_refit = 0;
    g_glt_cache.caps = 0u;
    g_glt_cache.ntiles = 0;
    g_glt_cache.num_pool = 0;
    g_glt_cache.bvh_nnodes = 0;
    g_glt_cache.bvh_root = 0;
    /* Cache gone -> no narrow-refit basis remains; force full on next build. */
    g_glt_dirty_mark_all_();
}

/* Drop ONLY the position/radius-dependent half, keeping pool/j_to_pool/num_pool.
 * Used when a caller needs geometry under a different radius policy or scale:
 * membership does not depend on those, so re-deriving it would be pure waste. */
static void glt_cache_free_geometry_(void)
{
    if(g_glt_cache.tiles)        { free(g_glt_cache.tiles);        g_glt_cache.tiles = NULL; }
    if(g_glt_cache.bvh)          { free(g_glt_cache.bvh);          g_glt_cache.bvh = NULL; }
    if(g_glt_cache.compact_xyzh) { free(g_glt_cache.compact_xyzh); g_glt_cache.compact_xyzh = NULL; }
    if(g_glt_cache.pool_types)   { free(g_glt_cache.pool_types);   g_glt_cache.pool_types = NULL; }
    g_glt_cache.ntiles = 0;
    g_glt_cache.bvh_nnodes = 0;
    g_glt_cache.bvh_root = 0;
    g_glt_cache.caps &= ~GX_POOL_GEOMETRY;
}

extern "C" void ghost_exchange_local_tree_invalidate_drift(void)
{
    if(g_glt_cache.valid) g_glt_cache.needs_refit = 1;
    /* Drift is a pool-wide event (every particle's Pos may have changed):
     * the narrow-refit fast path can't represent that, so promote to full. */
    g_glt_dirty_mark_all_();
}
extern "C" void ghost_exchange_local_tree_invalidate_full(void)  { glt_cache_free(); }

/* Read-only view of the owned-local supply pool + its epoch key, for the
 * top-leaf router band builder.  Returns num_pool, or -1 if the cache is not
 * safe to read (absent, or dirtied and awaiting refit -> compact_xyzh stale).
 * EXCLUDES ghosts by construction (the pool is the owned-local set the local
 * BVH walk uses); compact_xyzh[p*4+3] is the baked gx_policy_scaled_h reach. */
extern "C" int ghost_exchange_supply_pool_view(struct gx_supply_pool_view *out,
                                               unsigned int required_caps)
{
    if(!out || !required_caps) return -1;
    if(!g_glt_cache.valid || g_glt_cache.needs_refit) return -1;
    if((g_glt_cache.caps & required_caps) != required_caps) return -1;
    /* Identity is pool membership/order alone. pool_types is geometry-class — it is
     * read by the leaf supply-mask filter, a geometry walker — so it is required
     * only alongside the rest of the geometry. (A full rebuild happens to fill it
     * for identity-only entries too; it is absent only after the geometry half has
     * been dropped for a policy change.) */
    if(!g_glt_cache.pool) return -1;
    if((required_caps & GX_POOL_GEOMETRY)
       && (!g_glt_cache.compact_xyzh || !g_glt_cache.tiles
           || !g_glt_cache.bvh || !g_glt_cache.pool_types)) return -1;
    out->pool                  = g_glt_cache.pool;
    out->pool_types            = g_glt_cache.pool_types;
    out->compact_xyzh          = g_glt_cache.compact_xyzh;
    out->num_pool              = g_glt_cache.num_pool;
    out->tiles                 = g_glt_cache.tiles;
    out->ntiles                = g_glt_cache.ntiles;
    out->numpart_when_built    = g_glt_cache.NumPart_when_built;
    out->ti_when_built         = (long long)g_glt_cache.Ti_when_built;
    out->safety_when_built     = g_glt_cache.safety_factor_when_built;
    out->eligible_mask_when_built = g_glt_cache.eligible_type_mask_when_built;
    out->radius_policy_when_built = (int)g_glt_cache.radius_policy_when_built;
    out->j_scale_when_built    = g_glt_cache.j_radius_scale_when_built;
    return g_glt_cache.num_pool;
}

extern "C" void ghost_exchange_local_tree_mark_h_dirty_indices(const int *indices, int n)
{
    if(n <= 0 || !indices) return;
    if(g_glt_dirty_all) return; /* already covered by full-refit promotion */
    if((int)(g_glt_dirty_list.size() + (size_t)n) > G_GLT_DIRTY_PROMOTE_THRESHOLD) {
        g_glt_dirty_mark_all_();
        return;
    }
    g_glt_dirty_list.reserve(g_glt_dirty_list.size() + (size_t)n);
    for(int k = 0; k < n; k++) {
        int j = indices[k];
        if(j >= 0) g_glt_dirty_list.push_back(j);
    }
    /* Mark cache for refit-on-next-hit even though it isn't fully invalidated.
     * This wakes up the refit branch in the request-driven build path. */
    if(g_glt_cache.valid) g_glt_cache.needs_refit = 1;
}

extern "C" void ghost_exchange_local_tree_mark_h_dirty_range(int start, int end)
{
    if(end <= start) return;
    if(g_glt_dirty_all) return;
    int n = end - start;
    if((int)(g_glt_dirty_list.size() + (size_t)n) > G_GLT_DIRTY_PROMOTE_THRESHOLD) {
        g_glt_dirty_mark_all_();
        return;
    }
    g_glt_dirty_list.reserve(g_glt_dirty_list.size() + (size_t)n);
    for(int j = start; j < end; j++) g_glt_dirty_list.push_back(j);
    if(g_glt_cache.valid) g_glt_cache.needs_refit = 1;
}

/* TEMPORARY top-leaf-router oracle gate (stripped after the router is blessed).
 * GIZMO_GHOST_ROUTE_ORACLE=1 enables the compute-and-compare oracle in the
 * request-driven path: broadcast is the comparison reference; routed discovery is
 * computed only to verify it imports the SAME ghost set.  Env => identical on
 * every rank => collective-safe. */
static int ghost_route_oracle_enabled(void)
{
    static int initialized = 0;
    static int enabled = 0;
    if(initialized) return enabled;
    initialized = 1;
    const char *e = getenv("GIZMO_GHOST_ROUTE_ORACLE");
    enabled = (e && e[0] && e[0] != '0') ? 1 : 0;
    return enabled;
}

/* Stage-3 producer mode (Step 5).  Selects HOW the matched ghost set is produced;
 * Stages 1-2 (route CSR -> Alltoall -> Alltoallv -> received queries) are shared
 * SSOT regardless of mode.  HOST_ONLY = host walk only; HOST_AND_DEVICE_VALIDATE =
 * host walk + device Stage-3 compared (the Step-4 oracle); DEVICE_ONLY_AUTHORITY =
 * production device Stage-3 with NO host walk / NO broadcast / NO compare (reserved
 * for the device-routed install arm; no caller uses it yet). */
enum gx_producer_mode {
    GX_PRODUCER_HOST_ONLY = 0,
    GX_PRODUCER_HOST_AND_DEVICE_VALIDATE,
    GX_PRODUCER_DEVICE_ONLY_AUTHORITY
};

/* Which producer actually INSTALLED the ghost set this call (telemetry only; the
 * installer downstream is producer-agnostic).  broadcast = Allgatherv+global walk;
 * oneway_routed_bvh = the ONEWAY tile/BVH arm, reachable only as the oracle's
 *   comparison reference (it is NOT a fallback: an eligible spec that fails stops);
 * walk_export = the host walk-export producer, used by both search modes. */
enum gx_installed_producer {
    GX_INSTALLED_BROADCAST = 0,
    GX_INSTALLED_ONEWAY_ROUTED_BVH,
    GX_INSTALLED_WALK_EXPORT        /* routed set from the host walk-export producer (either mode) */
};

static const char *gx_installed_producer_name(enum gx_installed_producer p)
{
    switch(p) {
        case GX_INSTALLED_ONEWAY_ROUTED_BVH: return "oneway_routed_bvh";
        case GX_INSTALLED_WALK_EXPORT:      return "walk_export";
        case GX_INSTALLED_BROADCAST:         return "broadcast";
    }
    return "broadcast";
}

/* TEMPORARY H1 flat-vs-hierarchical owner-set oracle gate (stripped after the
 * hierarchical router is blessed).  GIZMO_GHOST_ROUTE_HIER_ORACLE=1: when routed
 * transport is active (ONEWAY), also build the hierarchical owner sets and verify
 * they EQUAL the flat router's; mismatch => collective controlled stop (the new
 * geometry/traversal is wrong).  Flat stays authoritative — the hierarchical
 * result is only compared, never installed. */
static int ghost_route_hier_oracle_enabled(void)
{
    static int initialized = 0;
    static int enabled = 0;
    if(initialized) return enabled;
    initialized = 1;
    const char *e = getenv("GIZMO_GHOST_ROUTE_HIER_ORACLE");
    enabled = (e && e[0] && e[0] != '0') ? 1 : 0;
    return enabled;
}

/* Recompute one tile's lo/hi/hmax/hmax_by_type fully from current P[] over its
 * pool members.  Also rewrites compact_xyzh + pool_types for those members.
 * Used by both full and narrow refit paths. */
static inline void glt_recompute_tile_(int t)
{
    sfc_tile_t *tile = &g_glt_cache.tiles[t];
    tile->hmax = 0;
    for(int tt = 0; tt < TILE_NUM_PTYPES; tt++) tile->hmax_by_type[tt] = 0;
    /* Empty tiles carry an INVERTED box so they stay neutral under the BVH's
     * min/max union and always fail the sphere-overlap test; a zeroed box would
     * stretch every ancestor to the coordinate origin and defeat pruning there. */
    if(tile->count <= 0) {
        for(int k = 0; k < 3; k++) { tile->lo[k] = MAX_REAL_NUMBER; tile->hi[k] = -MAX_REAL_NUMBER; }
        return;
    }
    /* SSOT supply-side reach pulled from the cache's stored policy/scale.
     * For runner-driven imports this carries the
     * Spec's radius_policy + nlr_spec_symmetric_j_radius_scale<Spec>();
     * for legacy ghost_exchange wrappers it carries
     * MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES + 1.0 → byte-equivalent to the
     * pre-policy code that read P[j].KernelRadius * safety_factor.
     * Leaf compact_xyzh[p*4+3] and tile band hmax_by_type[] use the
     * IDENTICAL gx_policy_scaled_h output → no leaf-vs-band scale gap. */
    const mode_b_radius_policy_t policy = g_glt_cache.radius_policy_when_built;
    const double j_scale = g_glt_cache.j_radius_scale_when_built;
    const double safety  = g_glt_cache.safety_factor_when_built;
    /* Eliminated elements (Mass <= 0) are excluded from the band and the bbox:
     * the pool is mass-filtered when built, but a cached entry outlives a
     * Mass->0 marking until rearrange_particle_sequence() compacts it away, and
     * a dead slot's stale reach would otherwise widen what this rank advertises
     * as supply.  Dropping it can only remove pairs WITH the dead element, which
     * must not be discovered anyway; live pairs carry their own reach.  The leaf
     * walk applies the same test, so leaf h and band hmax stay on one reach.
     * The bbox is seeded from the first live member rather than the tile's first
     * member, which may itself be dead. */
    int seeded = 0;
    for(int s = 0; s < tile->count; s++) {
        int p = tile->first + s;
        int j = g_glt_cache.pool[p];
        int pt = (int)P[j].Type;
        double h = gx_policy_scaled_h(j, policy, j_scale, safety);
        g_glt_cache.compact_xyzh[p*4+0] = (float)P[j].Pos[0];
        g_glt_cache.compact_xyzh[p*4+1] = (float)P[j].Pos[1];
        g_glt_cache.compact_xyzh[p*4+2] = (float)P[j].Pos[2];
        g_glt_cache.compact_xyzh[p*4+3] = (float)h;
        g_glt_cache.pool_types[p] = pt;
        if(P[j].Mass <= 0) continue;
        if(!seeded) {
            for(int k = 0; k < 3; k++) tile->lo[k] = tile->hi[k] = P[j].Pos[k];
            seeded = 1;
        } else {
            for(int k = 0; k < 3; k++) {
                if(P[j].Pos[k] < tile->lo[k]) tile->lo[k] = P[j].Pos[k];
                if(P[j].Pos[k] > tile->hi[k]) tile->hi[k] = P[j].Pos[k];
            }
        }
        if(h > tile->hmax) tile->hmax = h;
        if(pt >= 0 && pt < TILE_NUM_PTYPES && h > tile->hmax_by_type[pt])
            tile->hmax_by_type[pt] = h;
    }
    /* Every member dead: same inverted empty-tile box as count <= 0 above. */
    if(!seeded) { for(int k = 0; k < 3; k++) { tile->lo[k] = MAX_REAL_NUMBER; tile->hi[k] = -MAX_REAL_NUMBER; } }
}

/* Pull one BVH node's lo/hi/hmax/hmax_by_type from its children/leaf-tile.
 * BVH is in children-first order so calling this for indices 0..bvh_nnodes-1
 * in order updates internal nodes after their children. */
static inline void glt_recompute_bvh_node_(int n)
{
    tile_bvh_node_t *node = &g_glt_cache.bvh[n];
    if(node->left < 0) {
        int t = -(node->left + 1);
        if(t < 0 || t >= g_glt_cache.ntiles) return;
        sfc_tile_t *tile = &g_glt_cache.tiles[t];
        for(int k = 0; k < 3; k++) { node->lo[k] = tile->lo[k]; node->hi[k] = tile->hi[k]; }
        node->hmax = tile->hmax;
        for(int tt = 0; tt < TILE_NUM_PTYPES; tt++) node->hmax_by_type[tt] = tile->hmax_by_type[tt];
    } else {
        tile_bvh_node_t *left = &g_glt_cache.bvh[node->left];
        tile_bvh_node_t *right = &g_glt_cache.bvh[node->right];
        for(int k = 0; k < 3; k++) {
            node->lo[k] = DMIN(left->lo[k], right->lo[k]);
            node->hi[k] = DMAX(left->hi[k], right->hi[k]);
        }
        node->hmax = DMAX(left->hmax, right->hmax);
        for(int tt = 0; tt < TILE_NUM_PTYPES; tt++)
            node->hmax_by_type[tt] = DMAX(left->hmax_by_type[tt], right->hmax_by_type[tt]);
    }
}

/* Find the tile containing pool position pp via binary search on tile->first.
 * Tiles are stored in ascending pool-position order; tile.first values are
 * monotonically non-decreasing.  Returns -1 on out-of-range. */
static inline int glt_tile_of_pool_pos_(int pp)
{
    int lo = 0, hi = g_glt_cache.ntiles - 1;
    if(pp < 0 || g_glt_cache.ntiles <= 0) return -1;
    if(pp < g_glt_cache.tiles[0].first) return -1;
    while(lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if(g_glt_cache.tiles[mid].first <= pp) lo = mid; else hi = mid - 1;
    }
    /* Sanity: pp must lie within [first, first+count). */
    sfc_tile_t *tile = &g_glt_cache.tiles[lo];
    if(pp >= tile->first && pp < tile->first + tile->count) return lo;
    return -1;
}

static void glt_cache_refit_from_particles(void)
{
    if(!g_glt_cache.valid || !g_glt_cache.tiles || !g_glt_cache.pool ||
       !g_glt_cache.bvh || !g_glt_cache.compact_xyzh || !g_glt_cache.pool_types) return;

    /* Narrow path: refresh only tiles touched by g_glt_dirty_list, then walk
     * BVH bottom-up updating only ancestor nodes of those tiles. */
    if(!g_glt_dirty_all && !g_glt_dirty_list.empty() && g_glt_cache.j_to_pool) {
        int n_dirty = (int)g_glt_dirty_list.size();
        int ntiles = g_glt_cache.ntiles;
        int bvh_nnodes = g_glt_cache.bvh_nnodes;
        int NumPart_b = g_glt_cache.NumPart_when_built;
        unsigned char *tile_dirty = (unsigned char *) calloc((size_t)(ntiles > 0 ? ntiles : 1), 1);
        unsigned char *node_dirty = (unsigned char *) calloc((size_t)(bvh_nnodes > 0 ? bvh_nnodes : 1), 1);

        /* Phase 1: refresh compact_xyzh + pool_types for each dirty j; mark its tile dirty. */
        for(int k = 0; k < n_dirty; k++) {
            int j = g_glt_dirty_list[k];
            if(j < 0 || j >= NumPart_b) continue;
            int pp = g_glt_cache.j_to_pool[j];
            if(pp < 0 || pp >= g_glt_cache.num_pool) continue;
            int t = glt_tile_of_pool_pos_(pp);
            if(t < 0) continue;
            tile_dirty[t] = 1;
        }

        /* Phase 2: re-scan each touched tile fully (cheaper than tracking which
         * member exactly; tile_size ~64, dirty_count typically ~few-hundred). */
        for(int t = 0; t < ntiles; t++) {
            if(tile_dirty[t]) glt_recompute_tile_(t);
        }

        /* Phase 3: BVH bottom-up single pass.  At each leaf, propagate
         * tile_dirty[tile] -> node_dirty[n].  At each internal node, OR its
         * children's flags; if dirty, recompute lo/hi/hmax/hmax_by_type from
         * children.  Children-first ordering is guaranteed by build_bvh_recursive. */
        for(int n = 0; n < bvh_nnodes; n++) {
            tile_bvh_node_t *node = &g_glt_cache.bvh[n];
            if(node->left < 0) {
                int t = -(node->left + 1);
                if(t >= 0 && t < ntiles && tile_dirty[t]) {
                    node_dirty[n] = 1;
                    glt_recompute_bvh_node_(n);
                }
            } else {
                int L = node->left, R = node->right;
                int dL = (L >= 0 && L < bvh_nnodes) ? node_dirty[L] : 0;
                int dR = (R >= 0 && R < bvh_nnodes) ? node_dirty[R] : 0;
                if(dL || dR) {
                    node_dirty[n] = 1;
                    glt_recompute_bvh_node_(n);
                }
            }
        }

        free(tile_dirty);
        free(node_dirty);
        g_glt_dirty_clear_();
        g_glt_cache.Ti_when_built = All.Ti_Current;
        g_glt_cache.needs_refit = 0;
        g_glt_cache_narrow_refits++;
        return;
    }

    /* Full path: scan every tile + every BVH node. */
    for(int t = 0; t < g_glt_cache.ntiles; t++) glt_recompute_tile_(t);
    for(int n = 0; n < g_glt_cache.bvh_nnodes; n++) glt_recompute_bvh_node_(n);
    g_glt_dirty_clear_();
    g_glt_cache.Ti_when_built = All.Ti_Current;
    g_glt_cache.needs_refit = 0;
    g_glt_cache_refits++;
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
static ghost_exchange_result ghost_exchange_tile_overlap_impl(const struct ghost_exchange_spec_t *spec);
static void gx_print_waste(const struct ghost_exchange_spec_t *spec, int this_call, int total_recv);
static double gx_eff_h(int j, const struct ghost_exchange_spec_t *spec);

/* Is this spec eligible for the walk-export routed producer (sender fine-tree
 * export + bounded receiver walk)?  Keyed on the SEARCH MODE and structural spec
 * fields only — never on a caller name, so every loop of a given class routes.
 *
 * ONEWAY: always eligible.  Its reach is the query's own radius, which the sender
 * knows exactly, so routing needs nothing from the supply side.  The traversal AND
 * export opener are both R_open = h_q (mode_b_local_walker.cc, the R_open branch
 * and the topleaf export re-test), and the accept ignores h_j entirely
 * (ghost_exchange_functions.h gx_pair_accept, ONEWAY branch).  Query h already
 * carries the spec safety factor, so a widened query cannot outgrow the opener.
 * KEEP THOSE THREE IN STEP: if ONEWAY accept ever gains an h_j term, or the opener
 * stops using h_q, this eligibility no longer holds and must be re-derived.
 *
 * SYMMETRIC: eligible only when the supply-side reach is bounded by the per-type
 * node band the sender opener walks against — otherwise a reach beyond the band
 * silently under-imports.  Two structural conditions:
 *   supply_band_dominated  the spec's reach is proven bounded by that band
 *   safety_factor <= 1     a widened query (TURB_DIFF_DYNAMIC) outgrows the
 *                          band until the opener is scaled to match
 * Fails closed: anything unproven keeps the broadcast path it uses today.
 * Rank-uniform — search_mode and both fields are spec constants, identical on
 * every rank, so this never splits ranks across a collective. */
static inline int gx_walk_export_eligible(const struct ghost_exchange_spec_t *spec)
{
    if(!spec) return 0;
    if(spec->search_mode == NGB_SEARCH_ONEWAY) return 1;
    return spec->search_mode == NGB_SEARCH_SYMMETRIC
           && spec->supply_band_dominated
           && spec->safety_factor <= 1.0;
}

/* Announce, ONCE per caller per run, that a SYMMETRIC caller is on broadcast.
 * Once-only so a timing arm is never perturbed by per-call stdout; loud enough
 * that an unpromoted caller cannot hide in a log. */
static void gx_report_symm_broadcast(const struct ghost_exchange_spec_t *spec)
{
    enum { GX_SYMM_REPORT_MAX = 64 };   /* > the number of SYMMETRIC specs in the tree */
    static const char *seen[GX_SYMM_REPORT_MAX];
    static int n_seen = 0;
    static int table_full_reported = 0;
    if(ThisTask != 0 || !spec) return;
    const char *name = spec->caller_name ? spec->caller_name : "?";
    /* compare by CONTENT: callers pass distinct string objects (spec literals,
     * Spec::loop_name), so pointer identity would let one caller report twice. */
    for(int k = 0; k < n_seen; k++) { if(strcmp(seen[k], name) == 0) return; }
    if(n_seen >= GX_SYMM_REPORT_MAX) {
        /* Never degrade to per-call printing: that would flood a log and perturb
         * the very timing this report exists to keep honest. Say so once, then stop. */
        if(!table_full_reported) {
            table_full_reported = 1;
            printf("GHOST_SYMM_BCAST caller=<table full at %d> reason=further_callers_unreported\n",
                   GX_SYMM_REPORT_MAX);
            fflush(stdout);
        }
        return;
    }
    seen[n_seen++] = name;
    printf("GHOST_SYMM_BCAST caller=%s reason=%s\n", name,
           spec->supply_band_dominated ? "safety_gt_1" : "band_unproven");
    fflush(stdout);
}

static void ghost_exchange_impl(const struct ghost_exchange_spec_t *spec)
{
    /* Tiny-N corridor counter: increments on API entry, before any
     * dispatch. Mode B paths in run_neighbor_loop must NOT enter this
     * function. Counts even single-rank early-out cases by design. See
     * declarations/lifecycle_counters.h. */
    g_ghost_import_counter++;

    /* Phase 0 instrumentation: env-gated, all-ranks. Brackets dispatch so
     * both impls are captured without duplication. Off ⇒ no work beyond
     * one static int read. */
    static const char *g_phase0_env_raw = getenv("GIZMO_PHASE0_DIAG");
    static const int phase0_on = (g_phase0_env_raw && g_phase0_env_raw[0] == '1') ? 1 : 0;
    static long long g_phase0_ghost_call_id = 0;
    long long this_phase0_call = 0;
    double t_phase0_start = 0;
    int nlocal_pre = 0;
    if(phase0_on) {
        this_phase0_call = ++g_phase0_ghost_call_id;
        t_phase0_start = my_second();
        nlocal_pre = NumPart;
    }
    /* Dispatch policy: explicit-query callers (runner-issued specs,
     * n_queries >= 0) use the request-driven path — tile-overlap cannot consume
     * an explicit query list (it scans ActiveParticleList filtered by
     * request_type_mask, so a spec with request_type_mask=0u + an explicit list
     * would import zero ghosts).  Every ONEWAY request also uses request-driven,
     * regardless of caller: routed ONEWAY discovery is a property of the search
     * mode, not of any specific loop.  SYMMETRIC requests use request-driven when
     * the spec proves supply-band domination (gx_walk_export_eligible); all other
     * SYMMETRIC callers stay on tile-overlap/broadcast.  Which producer then
     * supplies the matched set inside the request-driven path is a separate,
     * likewise structural decision (same predicate). */
    const int explicit_queries = (spec && spec->n_queries >= 0);
    const int want_request_driven = explicit_queries
                                    || (spec && spec->search_mode == NGB_SEARCH_ONEWAY)
                                    || gx_walk_export_eligible(spec);
    if(spec && spec->search_mode == NGB_SEARCH_SYMMETRIC && !gx_walk_export_eligible(spec))
        gx_report_symm_broadcast(spec);
    const char *selected_impl;
    if(want_request_driven) {
        ghost_exchange_request_driven_impl(spec);
        selected_impl = "request_driven";
    } else {
        ghost_exchange_result result = ghost_exchange_tile_overlap_impl(spec);
        selected_impl = "tile_overlap";
        if(result != GHOST_EXCHANGE_COMPLETED) {
            /* Tile admission failed: it could not fit particle slots (or its counts exceeded
             * the int transport range) COLLECTIVELY — every rank returns the same
             * result (Allreduce in the tile impl), so all ranks take this branch
             * together. Drain any unrelated pending controlled-stop, then fall
             * back to exact request-driven discovery (changes discovery only, not
             * the downstream kernel). */
            const char *reason = (result == GHOST_EXCHANGE_COUNT_RANGE_EXCEEDED)
                                 ? "count_range" : "particle_capacity";
            if(ThisTask == 0) {
                printf("GHOST_ADMIT caller=%s attempted=tile selected=request_driven reason=%s\n",
                       spec->caller_name ? spec->caller_name : "?", reason);
                fflush(stdout);
            }
            gizmo_exit_bad_stop_if_requested("ghost_exchange:tile_fallback");
            ghost_exchange_request_driven_impl(spec);
            selected_impl = "request_driven(fallback)";
        }
    }
    if(phase0_on) {
        double dt_ghost_import = timediff(t_phase0_start, my_second());
        int ghost_added = NumPart - nlocal_pre;
        printf("PHASE0_GHOST rank=%d call=%lld caller=%s impl=%s "
               "nlocal_pre=%d ghost_added=%d ntotal_post=%d dt_ghost_import=%.6f\n",
               ThisTask, this_phase0_call,
               spec->caller_name ? spec->caller_name : "?",
               selected_impl,
               nlocal_pre, ghost_added, NumPart, dt_ghost_import);
        fflush(stdout);
    }
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

/* Tile-build SSOT for the tile ghost-exchange pool: particles are chunked
 * GHOST_TILE_TARGET at a time in pool order, giving the pool + tiling + per-tile
 * supply hmax. Any second consumer of this tiling must use the same value or the
 * two disagree about which particles share a tile. */
static constexpr int GHOST_TILE_TARGET = 64;

/* Does particle i join the tile supply pool? Mass-positive + type in supply_mask. */
static inline int ghost_tile_pool_includes(int i, unsigned int supply_mask)
{
    if(P[i].Mass <= 0) return 0;
    if(!ghost_type_passes((int)P[i].Type, supply_mask)) return 0;
    return 1;
}

/* Effective ghost search radius this rank must SUPPLY for particle j under the
 * tile-build rule. Typed callers use the bare KernelRadius (matches the legacy
 * ngb search; avoids conflating AGS / DM / wind kernels into the
 * hydro/sink/feedback radius). The all-types pool fans in the widest enabled
 * per-type kernel so ghost bboxes cover any active kernel. */
static inline double ghost_tile_effective_radius(int j, unsigned int supply_mask)
{
    if(supply_mask != GHOST_TYPE_ALL) {
        if(!ghost_type_passes((int)P[j].Type, supply_mask)) return 0.0;
        return (double)P[j].KernelRadius;
    }
    double h = P[j].KernelRadius;
#ifdef DM_DISPERSION_LOOP_ACTIVE
    if(P[j].Type == 0 && j < N_gas) {
        if((double)CellP[j].KernelRadiusDM > h) h = CellP[j].KernelRadiusDM;
    }
#endif
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if((double)P[j].AGS_KernelRadius > h) h = P[j].AGS_KernelRadius;
#endif
    return h;
}

/* Pure particle-slot fit predicate: do `required` total (local+ghost) particles
 * fit P[]/CellP[] (All.MaxPart)? Policy (what to do on a miss) lives in the
 * caller, NOT here — keep this free of multi-space/budget logic. */
static inline int ghost_particle_slots_fit(long long required)
{
    return (required <= (long long)All.MaxPart) ? 1 : 0;
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

/* SSOT for the forward particle+cell transport: the two typed Alltoallv calls
   (P always; CellP only when gas exists globally) with element-unit counts.
   Used by BOTH import impls (materialising ghosts at &P[NumPart]) and
   ghost_refresh_values() (overwriting existing ghost slots at
   &P[NumPart_before_ghost]). Verbose per-impl diagnostics stay at the call
   sites; only the transport is factored here. */
static void gx_forward_particle_exchange(const struct particle_data *send_P,
                                         const struct gas_cell_data *send_CellP,
                                         const int *send_count, const int *send_disp,
                                         struct particle_data *dst_P,
                                         struct gas_cell_data *dst_CellP,
                                         const int *recv_count, const int *recv_disp)
{
    gizmo_mpi_alltoallv_typed((void *)send_P, (int *)send_count, (int *)send_disp,
                              dst_P, (int *)recv_count, (int *)recv_disp,
                              sizeof(struct particle_data), MPI_COMM_WORLD);
    if(All.TotN_gas > 0) {
        gizmo_mpi_alltoallv_typed((void *)send_CellP, (int *)send_count, (int *)send_disp,
                                  dst_CellP, (int *)recv_count, (int *)recv_disp,
                                  sizeof(struct gas_cell_data), MPI_COMM_WORLD);
    }
}

static ghost_exchange_result ghost_exchange_tile_overlap_impl(const struct ghost_exchange_spec_t *spec)
{
    const double safety_factor = spec->safety_factor;
    const unsigned int request_mask = spec->request_type_mask;
    const unsigned int supply_mask  = spec->supply_type_mask;
    const int  search_mode = spec->search_mode;
    if(NTask <= 1) return GHOST_EXCHANGE_COMPLETED;
    double t_ghost_start = my_second(), t_ghost_phase;

    /* save current state for cleanup */
    NumPart_before_ghost = NumPart;
    N_gas_before_ghost = N_gas;
    NumGhostParticles = 0;

    int i, k, task;
    int tile_target = GHOST_TILE_TARGET;

    /* ================================================================
       Step 1: Build SFC tiles from local particles.
       Particles are already Peano-Hilbert sorted. Group into tiles of
       ~tile_target particles, computing bbox and hmax per tile.
       Uses malloc (not mymalloc) for tile metadata to avoid stack issues.
       ================================================================ */
    int local_ntiles = 0, num_pool = 0;

    /* Count pool particles. supply_mask gates which Types may join the pool;
     * removes tile-hmax pollution from non-supply types (e.g. DM init-time
     * KernelRadius poisoning hydro tile bboxes). */
    for(i = 0; i < NumPart; i++) {
        if(!ghost_tile_pool_includes(i, supply_mask)) continue;
        num_pool++;
    }

    /* Build pool index array */
    int *pool = (int *) malloc((num_pool > 0 ? num_pool : 1) * sizeof(int));
    int p = 0;
    for(i = 0; i < NumPart; i++) {
        if(!ghost_tile_pool_includes(i, supply_mask)) continue;
        pool[p++] = i;
    }

    local_ntiles = (num_pool + tile_target - 1) / tile_target;
    if(local_ntiles < 1) local_ntiles = 1;

    /* Compact tile metadata for exchange. Two parallel bbox/hmax sets:
     *   (lo, hi, hmax, count)         — over ALL particles in tile
     *   (active_lo, active_hi, active_hmax, active_count) — over ACTIVE only
     *
     * The all-particle set governs what this rank can SUPPLY as ghosts (j may
     * be inactive on its home rank but still a neighbor of an active i on a
     * peer rank — WAKEUP-style semantics demand this).
     *
     * The active-only set governs what this rank actually NEEDS: only tiles
     * containing at least one active particle drive remote-tile imports.
     *
     * Old code (pre-tile, tree-based ngb_treefind_variable_threads) iterated
     * FirstActiveParticle and built per-rank exports off active i's search
     * radius — the new tile path silently dropped that gating and was
     * shipping ~all of the global pool on tiny-N steps. */
    struct tile_meta_t {
        double lo[3], hi[3], hmax;
        int count;
        double active_lo[3], active_hi[3], active_hmax;
        int active_count;
    };

    tile_meta_t *local_meta = (tile_meta_t *) malloc(local_ntiles * sizeof(tile_meta_t));
    int *tile_first = (int *) malloc(local_ntiles * sizeof(int)); /* index into pool[] */

    /* Per-particle active flag — ActiveParticleList lookup at O(1). The active
     * stats below restrict the tile's "what do I need from peers?" criterion
     * to particles that will actually walk neighbors this step. Allocate as
     * char rather than bool to keep the address-stable contract. */
    char *is_active = (char *) calloc(NumPart > 0 ? NumPart : 1, sizeof(char));
    for(size_t kk = 0; kk < ActiveParticleList.size(); kk++) {
        int i_act = ActiveParticleList[kk];
        if(i_act < 0 || i_act >= NumPart) continue;
        /* request_mask gates which types' actives drive imports. */
        if(!ghost_type_passes((int)P[i_act].Type, request_mask)) continue;
        is_active[i_act] = 1;
    }

    for(int t = 0; t < local_ntiles; t++)
    {
        int start = t * tile_target;
        int count = tile_target;
        if(start + count > num_pool) count = num_pool - start;
        if(count < 0) count = 0;
        tile_first[t] = start;
        local_meta[t].count = count;
        local_meta[t].hmax = 0;
        local_meta[t].active_count = 0;
        local_meta[t].active_hmax = 0;
        for(k = 0; k < 3; k++) {
            local_meta[t].active_lo[k] = 0;  /* will be overwritten on first active */
            local_meta[t].active_hi[k] = 0;
        }
        if(count <= 0) {
            /* Rank owns no supply-pool particles for this mask -> empty sentinel
             * tile. Zero bbox/hmax, NO pool[] access (pool may be a 1-elt stub).
             * active_count=0 + hmax=0 means it neither needs nor supplies ghosts
             * downstream. */
            for(k = 0; k < 3; k++) { local_meta[t].lo[k] = 0; local_meta[t].hi[k] = 0; }
            continue;
        }

        int j0 = pool[start];
        for(k = 0; k < 3; k++) local_meta[t].lo[k] = local_meta[t].hi[k] = P[j0].Pos[k];
        {double h0 = ghost_tile_effective_radius(j0, supply_mask); if(h0 > local_meta[t].hmax) local_meta[t].hmax = h0;}
        if(is_active[j0]) {
            local_meta[t].active_count = 1;
            for(k = 0; k < 3; k++) local_meta[t].active_lo[k] = local_meta[t].active_hi[k] = P[j0].Pos[k];
            double h0 = ghost_tile_effective_radius(j0, supply_mask);
            local_meta[t].active_hmax = h0;
        }

        for(int s = 1; s < count; s++) {
            int j = pool[start + s];
            for(k = 0; k < 3; k++) {
                if(P[j].Pos[k] < local_meta[t].lo[k]) local_meta[t].lo[k] = P[j].Pos[k];
                if(P[j].Pos[k] > local_meta[t].hi[k]) local_meta[t].hi[k] = P[j].Pos[k];
            }
            double hj = ghost_tile_effective_radius(j, supply_mask);
            if(hj > local_meta[t].hmax) local_meta[t].hmax = hj;
            if(is_active[j]) {
                if(local_meta[t].active_count == 0) {
                    for(k = 0; k < 3; k++) local_meta[t].active_lo[k] = local_meta[t].active_hi[k] = P[j].Pos[k];
                    local_meta[t].active_hmax = hj;
                } else {
                    for(k = 0; k < 3; k++) {
                        if(P[j].Pos[k] < local_meta[t].active_lo[k]) local_meta[t].active_lo[k] = P[j].Pos[k];
                        if(P[j].Pos[k] > local_meta[t].active_hi[k]) local_meta[t].active_hi[k] = P[j].Pos[k];
                    }
                    if(hj > local_meta[t].active_hmax) local_meta[t].active_hmax = hj;
                }
                local_meta[t].active_count++;
            }
        }
    }
    free(is_active);

    /* Diagnostic: identify top-N tile-hmax outliers and what's setting them.
     * Gated on GIZMO_VERBOSE_DIAG=1 + first 3 calls to keep output manageable.
     * For each top tile we walk its particles to find the one(s) with the
     * largest effective_ghost_radius and dump Type/Mass/h/Pos/ID/TimeBin.
     * Goal: settle whether 582-scale tile hmax comes from DM init values,
     * sink kernel inflation, stale gas, or a real outlier. */
    {
        static int gx_topn_calls = 0;
        if(gizmo_verbose_diag() && gx_topn_calls < 3 && local_ntiles > 0) {
            const int TOPN = 8;
            int top_idx[TOPN]; double top_h[TOPN];
            for(int q = 0; q < TOPN; q++) { top_idx[q] = -1; top_h[q] = -1.0; }
            for(int t = 0; t < local_ntiles; t++) {
                double h = local_meta[t].hmax;
                int slot = -1;
                for(int q = 0; q < TOPN; q++) { if(h > top_h[q]) { slot = q; break; } }
                if(slot >= 0) {
                    for(int q = TOPN - 1; q > slot; q--) { top_h[q] = top_h[q-1]; top_idx[q] = top_idx[q-1]; }
                    top_h[slot] = h; top_idx[slot] = t;
                }
            }
            for(int q = 0; q < TOPN && top_idx[q] >= 0; q++) {
                int t = top_idx[q];
                int start = tile_first[t];
                int n = local_meta[t].count;
                int worst_j = -1; double worst_h = -1.0;
                for(int s = 0; s < n; s++) {
                    int j = pool[start + s];
                    double hj = ghost_tile_effective_radius(j, supply_mask);
                    if(hj > worst_h) { worst_h = hj; worst_j = j; }
                }
                if(worst_j >= 0) {
                    printf("[GX_TOPHMAX rank=%d call=%d t=%d tile_hmax=%.6g particle: idx=%d ID=%llu Type=%d Mass=%.4g KernelRadius=%.6g Pos=(%.4g,%.4g,%.4g) TimeBin=%d active_count=%d tile_count=%d]\n",
                           ThisTask, gx_topn_calls, t, top_h[q],
                           worst_j, (unsigned long long)P[worst_j].ID, (int)P[worst_j].Type,
                           (double)P[worst_j].Mass, (double)P[worst_j].KernelRadius,
                           (double)P[worst_j].Pos[0], (double)P[worst_j].Pos[1], (double)P[worst_j].Pos[2],
                           (int)P[worst_j].TimeBin, local_meta[t].active_count, local_meta[t].count);
                }
            }
            fflush(stdout);
            gx_topn_calls++;
        }
    }

    double t_ghost_tiles = timediff(t_ghost_start, my_second());

    /* ================================================================
       Step 2: Gather tile metadata from all ranks.
       Each rank sends its tile count and metadata to all ranks.
       ================================================================ */
    t_ghost_phase = my_second();
    int *all_ntiles = (int *) malloc(NTask * sizeof(int));
    MPI_Allgather(&local_ntiles, 1, MPI_INT, all_ntiles, 1, MPI_INT, MPI_COMM_WORLD);

    int *tile_disp = (int *) malloc(NTask * sizeof(int));
    int total_tiles = 0;
    for(task = 0; task < NTask; task++) {
        tile_disp[task] = total_tiles;
        total_tiles += all_ntiles[task];
    }

    /* Exchange tile metadata using MPI_Allgatherv with MPI_BYTE */
    tile_meta_t *all_meta = (tile_meta_t *) malloc(total_tiles * sizeof(tile_meta_t));
    int *meta_counts = (int *) malloc(NTask * sizeof(int));
    int *meta_disps = (int *) malloc(NTask * sizeof(int));
    for(task = 0; task < NTask; task++) {
        meta_counts[task] = all_ntiles[task] * sizeof(tile_meta_t);
        meta_disps[task] = tile_disp[task] * sizeof(tile_meta_t);
    }
    MPI_Allgatherv(local_meta, local_ntiles * sizeof(tile_meta_t), MPI_BYTE,
                   all_meta, meta_counts, meta_disps, MPI_BYTE, MPI_COMM_WORLD);
    free(meta_counts); free(meta_disps);

    if(ThisTask == 0) {
        double hmax_min = 1e30, hmax_max = 0, hmax_sum = 0;
        for(int t = 0; t < local_ntiles; t++) {
            if(local_meta[t].hmax < hmax_min) hmax_min = local_meta[t].hmax;
            if(local_meta[t].hmax > hmax_max) hmax_max = local_meta[t].hmax;
            hmax_sum += local_meta[t].hmax;
        }
        PRINT_STATUS("Ghost exchange (tile-based): %d local tiles, %d total across %d ranks, tile hmax=[%.4g, %.4g] avg=%.4g",
                     local_ntiles, total_tiles, NTask, hmax_min, hmax_max, hmax_sum / local_ntiles);
    }

    double t_ghost_meta = timediff(t_ghost_phase, my_second());
    double t_phase_overlap_start = my_second();

    /* Diagnostic: number this ghost_exchange call to track progress in multi-call steps */
    static int ghost_call_seq = 0;
    ghost_call_seq++;
    int this_call = ghost_call_seq;
    if(gizmo_verbose_diag()) {
        printf("[GX rank=%d call=%d] after_allgatherv: local_ntiles=%d total_tiles=%d NumPart=%d\n",
               ThisTask, this_call, local_ntiles, total_tiles, NumPart);
        fflush(stdout);
    }

    /* ================================================================
       Step 3: Per-task tile overlap check.
       For each remote task, check which of its tiles overlap with any
       of our tiles AND which of our tiles overlap with any of its tiles.
       No global Allreduce — each rank independently computes per-task
       recv and send lists using the symmetric overlap criterion.
       ================================================================ */

    /* Per-tile need flags: need_from[total_tiles] = do WE need this remote tile?
       send_to[local_ntiles * NTask] = does task t need our tile lt?
     *
     * Asymmetric criterion:
     *   need_from[rt]            = our active tiles' bbox vs rt's all-bbox
     *   send_to[lt][task]        = task's active tiles' bbox vs our lt's all-bbox
     * Active-aware on each rank's REQUEST side, all-particle on the SUPPLY
     * side (a remote inactive j may still be a neighbor of an active i — see
     * WAKEUP). Tile pairs where the request side has zero actives early-exit
     * — that's where the tiny-N wins come from. */
    int *need_from = (int *) calloc(total_tiles, sizeof(int));
    int *send_to = (int *) calloc(local_ntiles * NTask, sizeof(int));
    int my_tile_start = tile_disp[ThisTask];

    /* Per-axis min-AABB-AABB squared distance under periodic wrap. Returns
     * negative gap on this axis if the AABBs overlap. Inlined for hot loop. */
    auto axis_gap = [&](double c_a, double hw_a, double c_b, double hw_b, int kk) -> double {
#if defined(BOX_PERIODIC)
        int is_periodic = 1;
        double bsize = (kk==0) ? boxSize_X : ((kk==1) ? boxSize_Y : boxSize_Z);
#if defined(BOX_REFLECT_X)
        if(kk==0) is_periodic = 0;
#endif
#if defined(BOX_REFLECT_Y)
        if(kk==1) is_periodic = 0;
#endif
#if defined(BOX_REFLECT_Z)
        if(kk==2) is_periodic = 0;
#endif
#if defined(BOX_OUTFLOW_X)
        if(kk==0) is_periodic = 0;
#endif
#if defined(BOX_OUTFLOW_Y)
        if(kk==1) is_periodic = 0;
#endif
#if defined(BOX_OUTFLOW_Z)
        if(kk==2) is_periodic = 0;
#endif
#else
        int is_periodic = 0;
        double bsize = 0;
#endif
        double dx = fabs(c_a - c_b);
        if(is_periodic && dx > 0.5 * bsize) dx = bsize - dx;
        return dx - hw_a - hw_b;
    };

    /* Pass 1: need_from[rt] — driven by OUR active tiles. Outer loop over
     * local tiles with active_count > 0 only (tiny-N: just a handful). */
    for(int lt_idx = 0; lt_idx < local_ntiles; lt_idx++)
    {
        int lt = my_tile_start + lt_idx;
        tile_meta_t *lm = &all_meta[lt];
        if(lm->active_count == 0) continue;
        double c_lo[3], c_hw[3];
        for(k = 0; k < 3; k++) {
            c_lo[k] = 0.5 * (lm->active_lo[k] + lm->active_hi[k]);
            c_hw[k] = 0.5 * (lm->active_hi[k] - lm->active_lo[k]);
        }
        for(task = 0; task < NTask; task++)
        {
            if(task == ThisTask) continue;
            int t_start = tile_disp[task];
            int t_count = all_ntiles[task];
            for(int rt_idx = 0; rt_idx < t_count; rt_idx++)
            {
                int rt = t_start + rt_idx;
                if(need_from[rt]) continue;            /* already flagged by another lt */
                tile_meta_t *rm = &all_meta[rt];
                /* Search radius depends on caller mode:
                 *   ONEWAY (density): r_ij < h_i — only the LOCAL active's h
                 *     matters. Importing a remote tile because of its own
                 *     particles' large h_j is incorrect for density.
                 *   SYMMETRIC (gradients/sinks/etc.): r_ij < max(h_i, h_j) —
                 *     remote h_j matters because j's kernel may reach back
                 *     to i. */
                double search_r = (search_mode == NGB_SEARCH_ONEWAY)
                                    ? lm->active_hmax * safety_factor
                                    : DMAX(lm->active_hmax, rm->hmax) * safety_factor;
                if(search_r <= 0) continue;
                double search_r2 = search_r * search_r;
                double dist2 = 0;
                int overlaps = 1;
                for(k = 0; k < 3; k++) {
                    double c_r = 0.5 * (rm->lo[k] + rm->hi[k]);
                    double hw_r = 0.5 * (rm->hi[k] - rm->lo[k]);
                    double gap = axis_gap(c_lo[k], c_hw[k], c_r, hw_r, k);
                    if(gap <= 0) continue;
                    if(gap > search_r) { overlaps = 0; break; }
                    dist2 += gap * gap;
                }
                if(overlaps && dist2 < search_r2) need_from[rt] = 1;
            }
        }
    }

    /* Pass 2: derive send_to[lt][task] from peers' need_from directly.
     *
     * Earlier versions recomputed an "active-on-the-other-side" overlap test
     * here. That's algebraically the same predicate as peer's pass 1 — but
     * the two ranks running floating-point min-distance arithmetic on the
     * same all_meta blob can disagree on a tile pair that lies right at the
     * search-radius threshold. When they do, A says send_count[B]=N and B
     * says recv_count[A]=N±1, and the downstream Alltoallv truncates with
     * MPI_ERR_TRUNCATE. Killing the recompute entirely makes the two sides
     * bit-identical by construction.
     *
     * Cost of the Allgather: total_tiles*NTask ints (~1.5 MB at 2 ranks ×
     * 200k tiles, scales as O(NTask^2) — switch to Alltoall of per-rank
     * slices if NTask grows past ~100). The pass-1 cost is unchanged: outer
     * loop active-gated, ~handful of tiles on tiny-N. */
    int *all_need_from = (int *) malloc((size_t)NTask * total_tiles * sizeof(int));
    if(gizmo_verbose_diag()) {
        printf("[GX rank=%d call=%d] BEFORE Allgather(need_from): total_tiles=%d NTask=%d\n",
               ThisTask, this_call, total_tiles, NTask);
        fflush(stdout);
    }
    MPI_Allgather(need_from, total_tiles, MPI_INT,
                  all_need_from, total_tiles, MPI_INT, MPI_COMM_WORLD);
    if(gizmo_verbose_diag()) {
        printf("[GX rank=%d call=%d] AFTER  Allgather(need_from) OK\n", ThisTask, this_call);
        fflush(stdout);
    }
    for(task = 0; task < NTask; task++)
    {
        if(task == ThisTask) continue;
        const int *peer_need = all_need_from + (size_t)task * total_tiles;
        for(int lt_idx = 0; lt_idx < local_ntiles; lt_idx++) {
            send_to[lt_idx * NTask + task] = peer_need[my_tile_start + lt_idx];
        }
    }
    free(all_need_from);

    /* ================================================================
       Step 4: Compute per-task send/recv counts from overlap results.
       ================================================================ */
    int *recv_count = (int *) mymalloc("ghost_rc", NTask * sizeof(int));
    int *send_count = (int *) mymalloc("ghost_sc", NTask * sizeof(int));
    memset(recv_count, 0, NTask * sizeof(int));
    memset(send_count, 0, NTask * sizeof(int));

    for(task = 0; task < NTask; task++) {
        if(task == ThisTask) continue;
        /* Recv: remote tiles we need from this task */
        for(int rt = tile_disp[task]; rt < tile_disp[task] + all_ntiles[task]; rt++) {
            if(need_from[rt]) recv_count[task] += all_meta[rt].count;
        }
        /* Send: our tiles this task needs */
        for(int lt = 0; lt < local_ntiles; lt++) {
            if(send_to[lt * NTask + task]) send_count[task] += local_meta[lt].count;
        }
    }

    /* Compute totals + displacements with CHECKED int64 accumulation. The MPI
     * pack consumes int counts/displacements; a value exceeding INT_MAX cannot
     * be represented and must NOT silently wrap into a false "fits" decision.
     * Such a representation overflow is reported distinctly from a particle-slot
     * overflow (the collective slot admission below). Per-peer counts are bounded by one
     * rank's pool (<= INT_MAX); the prefix sums + totals are the overflow risk. */
    int *recv_disp = (int *) mymalloc("ghost_rd", NTask * sizeof(int));
    int *send_disp = (int *) mymalloc("ghost_sd", NTask * sizeof(int));
    long long total_recv_ll = 0, total_send_ll = 0;
    int count_range_ok = 1;
    {
        long long rdisp = 0, sdisp = 0;
        for(task = 0; task < NTask; task++) {
            if(rdisp <= INT_MAX) recv_disp[task] = (int)rdisp; else { recv_disp[task] = 0; count_range_ok = 0; }
            if(sdisp <= INT_MAX) send_disp[task] = (int)sdisp; else { send_disp[task] = 0; count_range_ok = 0; }
            rdisp += recv_count[task];
            sdisp += send_count[task];
        }
        total_recv_ll = rdisp; total_send_ll = sdisp;
        if(total_recv_ll > INT_MAX || total_send_ll > INT_MAX) count_range_ok = 0;
    }
    int total_recv = count_range_ok ? (int)total_recv_ll : 0;
    int total_send = count_range_ok ? (int)total_send_ll : 0;

    int tiles_needed = 0, tiles_sent = 0;
    for(int rt = 0; rt < total_tiles; rt++) tiles_needed += need_from[rt];
    for(int lt = 0; lt < local_ntiles; lt++) {
        for(task = 0; task < NTask; task++) { if(send_to[lt * NTask + task]) { tiles_sent++; break; } }
    }

    double t_ghost_overlap = timediff(t_phase_overlap_start, my_second()); /* steps 3+4 (overlap + schedule) */
    double t_phase_mpi_start = my_second();

    /* Frees the pre-materialisation tile scratch ONLY (mymalloc LIFO, then
     * malloc). ONE free list, shared by the Stage-0A fallback bail and the
     * normal-completion cleanup; captures only scratch allocated up to here (the
     * later packing allocs are freed explicitly by the normal path). It does NOT
     * touch NumGhostParticles: on normal completion that field already holds the
     * materialised ghost count and MUST survive (ghost-writeback / cleanup /
     * PreviousGhostCount read it); the fallback bail resets it explicitly. */
    auto tile_preflight_cleanup = [&]() {
        myfree(send_disp); myfree(recv_disp);
        myfree(send_count); myfree(recv_count);
        free(send_to); free(need_from);
        free(all_meta); free(tile_disp); free(all_ntiles);
        free(tile_first); free(local_meta); free(pool);
    };

    /* === Admission: collective particle-slot fit ===
     * Decide tile-vs-fallback BEFORE materialising any ghosts (NumPart unchanged
     * to here). Two DISTINCT failure modes — a representation overflow must never
     * masquerade as a capacity decision:
     *   2 = COUNT_RANGE — int MPI count/displacement representation overflowed
     *   1 = CAPACITY    — ghosts would not fit P[]/CellP[] (All.MaxPart)
     * Either, on ANY rank, makes ALL ranks abandon tile (clean rollback) and
     * fall back to exact request-driven discovery in the dispatcher. The Allreduce
     * is a continuing-branch decision, distinct from the terminal controlled-stop
     * poll (which still drains unrelated stops on the fit path below). */
    int local_fail = (!count_range_ok) ? 2
                   : (!ghost_particle_slots_fit((long long)NumPart + total_recv_ll)) ? 1 : 0;
    int global_fail = 0;
    MPI_Allreduce(&local_fail, &global_fail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(global_fail != 0) {
        tile_preflight_cleanup();
        NumGhostParticles = 0;   /* rollback: no ghosts materialised this call */
        return (global_fail == 2) ? GHOST_EXCHANGE_COUNT_RANGE_EXCEEDED
                                  : GHOST_EXCHANGE_PARTICLE_CAPACITY_EXCEEDED;
    }
    /* Tile fits on all ranks. Drain any UNRELATED pending controlled-stop
     * (preserves the original all-rank poll) before materialising ghosts. */
    gizmo_exit_bad_stop_if_requested("ghost_exchange:capacity");

    /* ================================================================
       Step 5: Pack particles from tiles needed by each task.
       ================================================================ */
    struct particle_data *send_P = (struct particle_data *) mymalloc("ghost_sP",
        (total_send > 0 ? total_send : 1) * sizeof(struct particle_data));
    struct gas_cell_data *send_CellP = (struct gas_cell_data *) mymalloc("ghost_sC",
        (total_send > 0 ? total_send : 1) * sizeof(struct gas_cell_data));
    /* Record home index of each sent particle for ghost writeback */
    int *send_home_idx = (int *) malloc((total_send > 0 ? total_send : 1) * sizeof(int));

    int *task_offset = (int *) mymalloc("ghost_toff", NTask * sizeof(int));
    memcpy(task_offset, send_disp, NTask * sizeof(int));

    for(task = 0; task < NTask; task++)
    {
        if(send_count[task] <= 0) continue;
        for(int t = 0; t < local_ntiles; t++)
        {
            if(!send_to[t * NTask + task]) continue;
            for(int s = 0; s < local_meta[t].count; s++)
            {
                if(task_offset[task] >= send_disp[task] + send_count[task]) break;
                int j = pool[tile_first[t] + s];
                int off = task_offset[task]++;
                gx_pack_send_slot(P, CellP, j, &send_P[off], &send_CellP[off]);
                send_home_idx[off] = j; /* record home index for writeback + refresh provenance */
            }
        }
    }

    /* ================================================================
       Step 6: Exchange via MPI_Alltoallv.
       ================================================================ */
    /* Per-particle exchange: send_count/recv_count/send_disp/recv_disp are
     * already in element units. gizmo_mpi_alltoallv_typed builds a contiguous
     * MPI_Datatype per call so element-count int*'s drive the wire — dodging
     * the 2.1 GB per-peer int-overflow that bites fire_m11i at >~6M parts/rank. */
    if(gizmo_verbose_diag()) {
        int ts=0, tr=0;
        for(int tt=0; tt<NTask; tt++) { ts += send_count[tt]; tr += recv_count[tt]; }
        printf("[GX rank=%d call=%d] BEFORE Alltoallv(P): total_send=%d total_recv=%d send[0]=%d recv[0]=%d\n",
               ThisTask, this_call, ts, tr, send_count[0], recv_count[0]);
        fflush(stdout);
    }
    gizmo_mpi_alltoallv_typed(send_P, send_count, send_disp,
                              &P[NumPart], recv_count, recv_disp,
                              sizeof(struct particle_data), MPI_COMM_WORLD);
    if(gizmo_verbose_diag()) {
        printf("[GX rank=%d call=%d] AFTER  Alltoallv(P) OK\n", ThisTask, this_call);
        fflush(stdout);
    }

    /* CellP exchange: only meaningful when the simulation has any gas
       particles globally. With TotN_gas==0 (N-body / DM-only runs), CellP
       is allocated to size 0, so writing to &CellP[NumPart] would dereference
       an out-of-bounds pointer. Skip the CellP alltoallv in that case —
       no gas ghosts can exist if no gas exists anywhere. */
    if(All.TotN_gas > 0) {
        if(gizmo_verbose_diag()) {
            printf("[GX rank=%d call=%d] BEFORE Alltoallv(CellP)\n", ThisTask, this_call);
            fflush(stdout);
        }
        gizmo_mpi_alltoallv_typed(send_CellP, send_count, send_disp,
                                  &CellP[NumPart], recv_count, recv_disp,
                                  sizeof(struct gas_cell_data), MPI_COMM_WORLD);
        if(gizmo_verbose_diag()) {
            printf("[GX rank=%d call=%d] AFTER  Alltoallv(CellP) OK\n", ThisTask, this_call);
            fflush(stdout);
        }
    }

    /* Update counts */
    NumGhostParticles = total_recv;
    NumPart += total_recv;

    /* Diagnostic: ghost composition by Type. If a hydro-context exchange is
     * pulling non-supply Type ghosts back, the type-mask refactor needs to
     * gate them out. Gated on GIZMO_VERBOSE_DIAG=1. */
    if(gizmo_verbose_diag() && total_recv > 0) {
        int by_type[6] = {0,0,0,0,0,0};
        for(int g = 0; g < total_recv; g++) {
            int gi = NumPart_before_ghost + g;
            int tt = (int)P[gi].Type;
            if(tt >= 0 && tt < 6) by_type[tt]++;
        }
        printf("[GX_GHOSTTYPE rank=%d call=%d caller=%s total_recv=%d  T0(cells)=%d T1=%d T2=%d T3=%d T4=%d T5=%d]\n",
               ThisTask, this_call, (spec->caller_name ? spec->caller_name : "?"),
               total_recv,
               by_type[0], by_type[1], by_type[2], by_type[3], by_type[4], by_type[5]);
        fflush(stdout);
    }

    /* Phase-0 import-waste diagnostic — see gx_print_waste(). */
    gx_print_waste(spec, this_call, total_recv);

    /* Multi-rank correctness: ghost slots [NumPart_before_ghost, NumPart) just
     * received fresh particle_data from remote ranks via MPI_Alltoallv. Their
     * KernelRadius values overwrote whatever was in those local P[] slots
     * (uninitialized or stale from prior step's ghost import). Any subsequent
     * cached gpu_ngb_list_build reading compact_xyzh[j*4+3] for h_j on those
     * slots would see stale values without this dirty-mark, silently missing
     * symmetric neighbor pairs where r<h_ghost.  Single-rank runs hit this
     * branch only when NumGhostParticles>0 (rare for true 1-rank), so this
     * fix is functionally a multi-rank correctness guarantee. */
    if(NumGhostParticles > 0) {
        gpu_compact_xyzh_mark_h_dirty_range(NumPart_before_ghost, NumPart);
    }
    /* SIDX lifecycle notify: fires on every rank for every exchange completion,
     * INCLUDING the count==0 / no-receive path. count==0 invalidates any
     * cached ghost segment from a prior import (no stale-ghost survival). */
    gpu_sidx_notify_ghost_imported(NumPart_before_ghost, NumGhostParticles);

    /* ================================================================
       Step 7: Build ghost provenance map for writeback.
       Exchange home indices so each ghost knows its home rank + index.
       ================================================================ */
    {
        /* Exchange home indices: send_home_idx[total_send] → recv_home_idx[total_recv] */
        int *recv_home_idx = (int *) malloc((total_recv > 0 ? total_recv : 1) * sizeof(int));
        if(gizmo_verbose_diag()) {
            printf("[GX rank=%d call=%d] BEFORE Alltoallv(home_idx)\n", ThisTask, this_call);
            fflush(stdout);
        }
        gizmo_mpi_alltoallv_typed(send_home_idx, send_count, send_disp,
                                  recv_home_idx, recv_count, recv_disp,
                                  sizeof(int), MPI_COMM_WORLD);
        if(gizmo_verbose_diag()) {
            printf("[GX rank=%d call=%d] AFTER  Alltoallv(home_idx) OK\n", ThisTask, this_call);
            fflush(stdout);
        }

        /* Build per-ghost provenance: home_rank and home_index */
        ghost_home_rank_map = (int *) malloc((total_recv > 0 ? total_recv : 1) * sizeof(int));
        ghost_home_index_map = recv_home_idx; /* take ownership — freed in cleanup */
        for(task = 0; task < NTask; task++) {
            for(int g = 0; g < recv_count[task]; g++) {
                ghost_home_rank_map[recv_disp[task] + g] = task;
            }
        }

        /* Preserve comm maps for reverse Alltoallv (malloc copies of mymalloc'd arrays) */
        ghost_wb_recv_count = (int *) malloc(NTask * sizeof(int));
        ghost_wb_recv_disp  = (int *) malloc(NTask * sizeof(int));
        ghost_wb_send_count = (int *) malloc(NTask * sizeof(int));
        ghost_wb_send_disp  = (int *) malloc(NTask * sizeof(int));
        memcpy(ghost_wb_recv_count, recv_count, NTask * sizeof(int));
        memcpy(ghost_wb_recv_disp,  recv_disp,  NTask * sizeof(int));
        memcpy(ghost_wb_send_count, send_count, NTask * sizeof(int));
        memcpy(ghost_wb_send_disp,  send_disp,  NTask * sizeof(int));

        /* Preserve send-side provenance for ghost_refresh_values() (take
           ownership of send_home_idx; the free() below then no-ops on NULL). */
        if(ghost_send_home_idx) free(ghost_send_home_idx);
        ghost_send_home_idx   = send_home_idx;
        ghost_send_home_count = total_send;
        send_home_idx         = NULL;
        g_ghost_provenance_epoch++;
    }

    double t_ghost_mpi = timediff(t_phase_mpi_start, my_second()); /* steps 5+6 (pack + MPI + unpack) */
    double t_ghost_end = my_second();
    double t_ghost_total = timediff(t_ghost_start, t_ghost_end);

    /* Active-count diagnostic (gated on GIZMO_VERBOSE_DIAG; collective so all
     * ranks must call). Lets us correlate ghost-exchange wall with how many
     * particles are actually active on this step. */
    int n_active_global = 0;
    if(gizmo_verbose_diag()) {
        int n_active = (int)ActiveParticleList.size();
        printf("[GX rank=%d call=%d] BEFORE MPI_Reduce(n_active=%d)\n", ThisTask, this_call, n_active);
        fflush(stdout);
        MPI_Reduce(&n_active, &n_active_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        printf("[GX rank=%d call=%d] AFTER  MPI_Reduce OK active_global=%d\n", ThisTask, this_call, n_active_global);
        fflush(stdout);
    }
    if(ThisTask == 0) {
        PRINT_STATUS("Ghost exchange: %d local + %d ghost = %d total (recv %d tiles, sent %d/%d) [%.4f s]",
                     NumPart_before_ghost, NumGhostParticles, NumPart,
                     tiles_needed, tiles_sent, local_ntiles, t_ghost_total);
        if(gizmo_verbose_diag()) {
            printf("  ghost_exchange phases: tiles_build=%.4f meta_allgather=%.4f overlap+sched=%.4f pack+mpi=%.4f total=%.4f  active_global=%d\n",
                   t_ghost_tiles, t_ghost_meta, t_ghost_overlap, t_ghost_mpi, t_ghost_total,
                   n_active_global);
            fflush(stdout);
        }
    }
    /* Warn if ghost particles used >80% of available headroom */
    if(NumPart > 0.8 * All.MaxPart) {
        double usage_frac = (double)NumPart / (double)All.MaxPart;
        PRINT_WARNING("Ghost exchange: particle arrays %.0f%% full (%d/%d). "
                      "Consider increasing PartAllocFactor (currently %.2f) to avoid running out of space.",
                      100.0 * usage_frac, NumPart, All.MaxPart, All.PartAllocFactor);
    }

    /* Cleanup: later packing allocs (mymalloc LIFO + malloc), then the shared
     * preflight cleanup frees the discovery scratch (same free list the Stage-0A
     * fallback bail uses). */
    myfree(task_offset);
    myfree(send_CellP); myfree(send_P);
    free(send_home_idx);
    tile_preflight_cleanup();
    return GHOST_EXCHANGE_COMPLETED;
}


/* ============================================================================
 * Request-driven ghost exchange.
 *
 * Replaces the tile-overlap candidate generator + whole-tile MPI payload of
 * ghost_exchange_tile_overlap_impl with:
 *   1. Each rank builds a list of compact query records {pos[3], h, _pad}
 *      for each LOCAL active particle matching spec->request_type_mask.
 *   2. Per-rank query counts via Allgather; queries themselves via Allgatherv.
 *   3. For each remote rank's queries, walk the LOCAL pool and apply the
 *      EXACT predicate per particle (r_ij < h_i for ONEWAY,
 *      r_ij < max(h_i, h_j) for SYMMETRIC). Mark matched local pool indices
 *      in a bitmask per peer rank to dedupe (same j matched by N queries
 *      ships once).
 *   4. Per-peer match list is the basis for Alltoallv pack.
 *   5. Alltoallv particle_data + gas_cell_data + home_idx (existing typed
 *      Alltoallv used for the legacy path).
 *   6. Install as ghosts (NumPart += total_recv, mark dirty for compact_xyzh,
 *      build ghost_home_*_map, save ghost_wb_* arrays).
 *
 * Killshot diagnostic on Vista (job 694703) showed the legacy path imports
 * with 99.97-100% waste — millions of particles imported per call, of which
 * <1000 are actual neighbors under the criterion that imported them. This
 * path's match step performs the EXACT predicate before pack, so per-call
 * waste is zero by construction (modulo a small over-supply factor from
 * ghosts that pass for one query but turn out unused after kernel).
 *
 * For now the local-pool walk is a flat O(N_local_tiles * N_remote_queries)
 * scan with bbox-vs-sphere prune at the tile level; per-particle accept at
 * the leaf. A BVH speedup is straightforward later (the all-particle BVH
 * already used by the GPU neighbor list can be reused) but not a
 * blocker here — even the flat scan is dominated by per-particle work for the
 * tiny-N case that motivates this restructure.
 * ============================================================================ */
struct gx_query_t {
    double pos[3];
    double h;
    int    type;       /* for caller diagnostics; predicate uses h_i directly */
    int    _pad;
};

/* Fixed-size export envelope carried from sender to supply rank: one query plus
 * the start nodes its sender walk reached on that peer.  A (query,peer) needing
 * more than NODELISTLENGTH nodes SPLITS into multiple envelopes (the Mode-B
 * export convention); the receiver walks each independently and dedups through
 * the matched bitmap, so a split costs an extra envelope and nothing else. */
struct gx_export_envelope_t {
    double pos[3];
    double h;
    int    n_nodes;
    int    nodes[NODELISTLENGTH];
    int    _pad;
};

/* Per-particle supply-side reach for the [GX_WASTE] diagnostic.  Policy-aware:
 * under the SSOT supply contract, the diagnostic must report the SAME reach
 * the request-driven impl actually used — otherwise it silently lies about
 * whether ghost imports are oversized.  Returns 0 when
 * j's type is outside the spec's supply_type_mask. */
static double gx_eff_h(int j, const struct ghost_exchange_spec_t *spec)
{
    if(!ghost_type_passes((int)P[j].Type, spec->supply_type_mask)) return 0.0;
    return gx_policy_scaled_h(j, spec->radius_policy,
                              spec->j_radius_scale, spec->safety_factor);
}

/* Print [GX_WASTE] for any ghost_exchange path. Walks (sampled) local actives
 * × imported ghosts, applies ONEWAY and SYMMETRIC predicates, prints the
 * per-ghost OR-aggregate waste ratio. PAIRS_BUDGET caps cost on global steps. */
static void gx_print_waste(const struct ghost_exchange_spec_t *spec, int this_call, int total_recv)
{
    if(!gizmo_verbose_diag()) return;
    if(total_recv <= 0 || NumPart_before_ghost <= 0) return;
    const unsigned int request_mask = spec->request_type_mask;
    const int  search_mode = spec->search_mode;
    const long long PAIRS_BUDGET = 10000000LL;
    int sample_cap = (int)(PAIRS_BUDGET / (long long)total_recv);
    if(sample_cap < 4) sample_cap = 4;
    if(sample_cap > 1024) sample_cap = 1024;
    int n_active_sample = 0;
    int *active_indices = (int *) malloc((size_t)sample_cap * sizeof(int));
    for(size_t kk = 0; kk < ActiveParticleList.size() && n_active_sample < sample_cap; kk++) {
        int i_act = ActiveParticleList[kk];
        if(i_act < 0 || i_act >= NumPart_before_ghost) continue;
        if(!ghost_type_passes((int)P[i_act].Type, request_mask)) continue;
        active_indices[n_active_sample++] = i_act;
    }
    int total_active_full = 0;
    for(size_t kk = 0; kk < ActiveParticleList.size(); kk++) {
        int i_act = ActiveParticleList[kk];
        if(i_act < 0 || i_act >= NumPart_before_ghost) continue;
        if(!ghost_type_passes((int)P[i_act].Type, request_mask)) continue;
        total_active_full++;
    }
    long long pairs_tested = 0;
    long long g_oneway_used = 0, g_symm_used = 0;
    char *used_oneway = (char *) calloc(total_recv, sizeof(char));
    char *used_symm   = (char *) calloc(total_recv, sizeof(char));
    for(int aa = 0; aa < n_active_sample; aa++) {
        int i = active_indices[aa];
        /* Request-side h_i: under SSOT contract the QUERY radius comes from
         * Spec::search_radius (= active_radii[a] * safety) on the runner path.
         * For the diagnostic we don't have direct per-active access to that
         * vector, so we approximate with the same supply-side gx_eff_h reach
         * — accurate when search_radius matches the policy reach (the typical
         * case for the Specs that use this diagnostic). Bounded error for the
         * waste percentage; not used for any correctness gate. */
        double h_i = gx_eff_h(i, spec);
        double h2_i = h_i * h_i;
        double px = P[i].Pos[0], py = P[i].Pos[1], pz = P[i].Pos[2];
        for(int g = 0; g < total_recv; g++) {
            int gi = NumPart_before_ghost + g;
            double dx_raw = px - P[gi].Pos[0];
            double dy_raw = py - P[gi].Pos[1];
            double dz_raw = pz - P[gi].Pos[2];
            MyDouble xtmp = 0; (void)xtmp;
            double adx = NGB_PERIODIC_BOX_LONG_X(dx_raw, dy_raw, dz_raw, 1);
            double ady = NGB_PERIODIC_BOX_LONG_Y(dx_raw, dy_raw, dz_raw, 1);
            double adz = NGB_PERIODIC_BOX_LONG_Z(dx_raw, dy_raw, dz_raw, 1);
            double r2 = adx*adx + ady*ady + adz*adz;
            pairs_tested++;
            if(!used_oneway[g] && r2 < h2_i) used_oneway[g] = 1;
            if(!used_symm[g]) {
                double h_j = gx_eff_h(gi, spec);
                double h2_max = (h2_i > h_j*h_j) ? h2_i : h_j*h_j;
                if(r2 < h2_max) used_symm[g] = 1;
            }
        }
    }
    for(int g = 0; g < total_recv; g++) {
        g_oneway_used += used_oneway[g];
        g_symm_used   += used_symm[g];
    }
    free(used_oneway); free(used_symm); free(active_indices);
    double waste_o = 100.0 * (1.0 - (double)g_oneway_used / (double)total_recv);
    double waste_s = 100.0 * (1.0 - (double)g_symm_used   / (double)total_recv);
    printf("[GX_WASTE rank=%d call=%d caller=%s mode=%s imported=%d n_active=%d (sampled=%d) used_oneway=%lld used_symm=%lld waste_oneway=%.2f%% waste_symm=%.2f%% pairs_tested=%lld]\n",
           ThisTask, this_call,
           (spec->caller_name ? spec->caller_name : "?"),
           (search_mode == NGB_SEARCH_ONEWAY ? "ONEWAY" : "SYMMETRIC"),
           total_recv, total_active_full, n_active_sample,
           g_oneway_used, g_symm_used, waste_o, waste_s, pairs_tested);
    fflush(stdout);
}

/* Host-only BVH bbox-vs-sphere overlap test. Mirrors bbox_overlaps_sphere_gpu
 * (mesh/sfc_tiles_functions.h). For each axis, take the periodic-shortened
 * gap between sphere center and bbox; if any gap exceeds search_r, prune.
 * Otherwise sum-of-squares vs search_r2 for the final accept. */
static inline int gx_bbox_overlaps_sphere(const double bbox_lo[3], const double bbox_hi[3],
                                          const double pos[3], double search_r, double search_r2,
                                          const int periodic_flags[3], const double box_sizes[3])
{
    double dist2 = 0;
    for(int k = 0; k < 3; k++) {
        /* Direct point-to-interval gap from lo/hi (a true lower bound on the
         * point-to-contained-particle distance).  This avoids the center/half-width
         * decomposition (0.5*(lo+hi), 0.5*(hi-lo)) whose rounding could make the gap
         * marginally NON-conservative at the search boundary and prune a reachable
         * particle.  Periodic axes: test the nearest images (tile width < box) and
         * take the minimum, so the gap is the true min-image distance — at least as
         * conservative as the prior center/min-image form. */
        double lo = bbox_lo[k], hi = bbox_hi[k], p = pos[k];
        double gap = (p < lo) ? (lo - p) : ((p > hi) ? (p - hi) : 0.0);
        if(periodic_flags[k] && gap > 0.0) {
            double L = box_sizes[k];
            double pm = p - L, pp = p + L;
            double gm = (pm < lo) ? (lo - pm) : ((pm > hi) ? (pm - hi) : 0.0);
            double gp = (pp < lo) ? (lo - pp) : ((pp > hi) ? (pp - hi) : 0.0);
            if(gm < gap) gap = gm;
            if(gp < gap) gap = gp;
        }
        if(gap <= 0.0) continue;          /* this axis overlaps */
        if(gap > search_r) return 0;      /* prune fast */
        dist2 += gap * gap;
    }
    return (dist2 < search_r2) ? 1 : 0;
}

/* Host-only BVH walk for the request-driven path. Mirrors
 * search_neighbors_sfc_gpu's iterative stack-based walk. For each query,
 * traverses the local BVH, accepts at leaves with per-particle r_ij vs the
 * caller-specified predicate (ONEWAY or SYMMETRIC), marks matched pool
 * indices in match_bitmask. */
/* Compute hmax_eff over the supply types only. Avoids the scalar-hmax
 * contamination where a particle of a type the caller didn't ask for poisons
 * the opener. Inlined; on a hot path. */
static inline double gx_node_hmax_supply(const tile_bvh_node_t *node, unsigned int supply_mask)
{
    double m = 0;
    for(int t = 0; t < TILE_NUM_PTYPES; t++) {
        if((supply_mask & (1u << t)) == 0u) continue;
        if(node->hmax_by_type[t] > m) m = node->hmax_by_type[t];
    }
    return m;
}

static void gx_walk_local_bvh(const float *compact_xyzh,
                              const sfc_tile_t *tiles, int ntiles,
                              const int *pool, int num_pool,
                              const int *pool_types, /* [num_pool] — P[].Type per pool slot, for leaf supply-mask filter */
                              unsigned int supply_mask,
                              const tile_bvh_node_t *bvh, int bvh_root,
                              const double pos_q[3], double h_q, int search_mode,
                              const int periodic_flags[3], const double box_sizes[3],
                              char *match_bitmask /* size num_pool */,
                              long *n_exact_hits /* optional (NULL in production): count EXACT matches
                                                  * this call would produce, computing r2 even for
                                                  * already-set slots so the per-query count is
                                                  * unbiased by the dedup skip (over-route diagnostic) */)
{
    (void)ntiles;
    int stack[TILE_BVH_STACK_SIZE];
    int sp = 0;
    stack[sp++] = bvh_root;
    while(sp > 0) {
        int node_idx = stack[--sp];
        const tile_bvh_node_t *node = &bvh[node_idx];
        /* OPENER criterion. ONEWAY: r_ij<h_q so search_r is just h_q (subtree
         * h's irrelevant). SYMMETRIC: r_ij<max(h_q,h_j); search_r = max(h_q,
         * subtree-max-h-of-supply-types). Per-type hmax filter avoids node
         * hmax being dominated by a type the caller didn't ask for. */
        double node_hmax_eff = gx_node_hmax_supply(node, supply_mask);
        double search_r = (search_mode == NGB_SEARCH_ONEWAY)
                            ? h_q
                            : ((h_q > node_hmax_eff) ? h_q : node_hmax_eff);
        if(search_r <= 0) continue;
        double search_r2 = search_r * search_r;
        if(!gx_bbox_overlaps_sphere(node->lo, node->hi, pos_q, search_r, search_r2,
                                    periodic_flags, box_sizes)) continue;
        if(node->left < 0) {
            /* Leaf: per-particle accept against EXACT predicate. */
            int tile_idx = -(node->left + 1);
            const sfc_tile_t *tile = &tiles[tile_idx];
            for(int s = 0; s < tile->count; s++) {
                int pool_pos = tile->first + s;
                if(pool_pos < 0 || pool_pos >= num_pool) continue;
                int already = match_bitmask[pool_pos];
                if(already && !n_exact_hits) continue;   /* already matched: skip unless the caller wants unbiased hit counts */
                /* Supply-mask filter at leaf: skip particles of a type the
                 * caller didn't ask for (no-op when tree was built with the
                 * same mask, but required when tree is shared across callers). */
                if(pool_types) {
                    int pt = pool_types[pool_pos];
                    if(pt < 0 || pt >= TILE_NUM_PTYPES) continue;
                    if((supply_mask & (1u << (unsigned)pt)) == 0u) continue;
                }
                /* Geometric accept is the shared SSOT predicate; supply-mask
                 * filter + dedup bookkeeping stay caller-side (above/here).
                 * BOTH position AND reach are DOUBLE: P[j].Pos (float absolute
                 * coordinates are invalid for GIZMO's dynamic range) and the double
                 * gx_policy_scaled_h reach (the node/tile hmax the opener prunes with
                 * is built from this same double reach, so the opener conservatively
                 * dominates the leaf — a float leaf h would let the opener under-prune
                 * boundary pairs).  Recompute here for the host path; the production
                 * GPU path needs a double h cache to avoid the per-candidate recompute.
                 * The pool walked here is always g_glt_cache's, so its build-time
                 * policy/scale/safety reproduce the cached double reach exactly. */
                int j_neighbor = pool[pool_pos];
                /* The pool is mass-filtered when it is built, but this pool is
                 * cached across calls and an entry outlives a Mass->0 marking
                 * until rearrange_particle_sequence() compacts it away.  Zero mass
                 * is the code's marker for an eliminated element, which must never
                 * be offered as a ghost source, so re-check it here instead of
                 * trusting build time.  Ahead of the reach recompute so a dead
                 * slot costs one compare. */
                if(P[j_neighbor].Mass <= 0) continue;
                double hj_dbl = gx_policy_scaled_h(j_neighbor, g_glt_cache.radius_policy_when_built,
                                                   g_glt_cache.j_radius_scale_when_built,
                                                   g_glt_cache.safety_factor_when_built);
                if(gx_pair_accept(pos_q, h_q,
                                  P[j_neighbor].Pos[0], P[j_neighbor].Pos[1], P[j_neighbor].Pos[2],
                                  hj_dbl, search_mode, periodic_flags, box_sizes)) {
                    if(n_exact_hits) (*n_exact_hits)++;
                    if(!already) match_bitmask[pool_pos] = 1;
                }
            }
        } else {
            if(sp + 2 > TILE_BVH_STACK_SIZE) break;  /* defensive — shouldn't happen */
            stack[sp++] = node->left;
            stack[sp++] = node->right;
        }
    }
}

/* Conservative sphere-vs-node-cube overlap — trivial pass-through to the shared
 * scalar predicate gx_node_sphere_overlap_center_len (ghost_exchange_functions.h),
 * the SSOT the device fine-tree walk also uses.  Extracts Nodes[] center/len here so
 * the shared helper stays geometry-only (no NODE/globals). */
static inline int gx_node_sphere_overlap(const double pos_q[3], const struct NODE *nop, double R,
                                         const int periodic_flags[3], const double box_sizes[3])
{
    return gx_node_sphere_overlap_center_len(pos_q, (double)nop->center[0], (double)nop->center[1],
                                             (double)nop->center[2], (double)nop->len, R,
                                             periodic_flags, box_sizes);
}



/* ============================================================================
 * TEMPORARY top-leaf-router oracle (GIZMO_GHOST_ROUTE_ORACLE; stripped after the
 * router is blessed).  Recomputes ghost discovery via top-leaf-targeted routing
 * and verifies the resulting imported ghost set is IDENTICAL to the broadcast
 * one.  Broadcast is the comparison reference; routed results are compared, never used.
 * ONEWAY callers only in this stage (band-free).  Collective-safe: every rank
 * runs the same collectives (gated on env + ONEWAY, both uniform); per-rank
 * geometry/route availability is folded into all-or-none Allreduce gates before
 * any routed collective.  Comparison is on sorted/deduped (home_rank,home_index)
 * SETS, not counts; under-route => collective controlled stop.
 * ========================================================================== */
struct gx_oracle_id { int rank; int idx; };
static int gx_oracle_id_cmp(const void *a, const void *b)
{
    const struct gx_oracle_id *x = (const struct gx_oracle_id *)a;
    const struct gx_oracle_id *y = (const struct gx_oracle_id *)b;
    if(x->rank != y->rank) return (x->rank < y->rank) ? -1 : 1;
    if(x->idx  != y->idx)  return (x->idx  < y->idx)  ? -1 : 1;
    return 0;
}
/* sort + unique in place; returns the deduped length. */
static int gx_oracle_sort_unique(struct gx_oracle_id *v, int n)
{
    if(n <= 1) return n;
    qsort(v, (size_t)n, sizeof(struct gx_oracle_id), gx_oracle_id_cmp);
    int w = 1;
    for(int i = 1; i < n; i++) {
        if(gx_oracle_id_cmp(&v[i], &v[w-1]) != 0) v[w++] = v[i];
    }
    return w;
}

/* Checked exclusive-prefix: disp[t]=sum(count[0..t)); returns total, or -1 if any
 * displacement/total exceeds INT_MAX (oracle buffers can be large at FIRE/224-rank
 * scale — fail closed, same no-silent-int-wrap rule as the real RD path). */
static long long gx_oracle_checked_prefix(const int *count, int *disp, int ntask)
{
    long long tot = 0;
    for(int t = 0; t < ntask; t++) {
        if(tot > (long long)INT_MAX) return -1;
        disp[t] = (int)tot;
        tot += count[t];
    }
    return (tot > (long long)INT_MAX) ? -1 : tot;
}

static void ghost_route_oracle_compare(
    const struct ghost_exchange_spec_t *spec, int this_call,
    struct gx_query_t *local_queries, int n_local_queries,
    sfc_tile_t *h_tiles, int ntiles, int *h_pool, int num_pool,
    int *h_pool_types, float *h_compact_xyzh, tile_bvh_node_t *h_bvh, int bvh_root,
    unsigned int supply_mask, int search_mode,
    const int periodic_flags[3], const double box_sizes[3],
    const int *ghost_home_rank_map, const int *ghost_home_index_map, int total_recv)
{
    /* All buffers NULL-initialised so the single cleanup at `done:` is safe on any
     * fail-closed early exit (free(NULL) is a no-op).  Every exit past a collective
     * barrier is reached by ALL ranks together (the barrier is an Allreduce), so the
     * goto-done skips never desync the routed collectives. */
    const char *cname = spec->caller_name ? spec->caller_name : "?";
    int    *route_off = NULL, *route_owners = NULL;
    double *q_pos = NULL, *q_h = NULL;
    int    *rq_sc = NULL, *rq_sd = NULL, *rq_rc = NULL, *rq_rd = NULL;
    struct gx_query_t   *rq_send = NULL, *rq_recv = NULL;
    char   *matched_r = NULL;
    int    *gid_sc = NULL, *gid_sd = NULL, *gid_rc = NULL, *gid_rd = NULL;
    struct gx_oracle_id *gid_send = NULL, *gid_recv = NULL, *bcast = NULL;
    int  oneway = (search_mode == NGB_SEARCH_ONEWAY) ? 1 : 0;
    int  nq = n_local_queries;
    long owners_cap = 0;
    int  aborted = 0, bad_any = 0, skipped = 1;
    long long rq_ts_ll = 0, rq_tr_ll = 0, gid_ts_ll = 0, gid_tr_ll = 0;
    int  rq_total_send = 0, rq_total_recv = 0, gid_total_send = 0, gid_total_recv = 0;
    int  n_routed = 0, n_bcast = 0, local_mismatch = 0, mismatch_any = 0;
    struct gx_oracle_id first_missing; first_missing.rank = -1; first_missing.idx = -1;
    double t_oracle_start = my_second();

    /* Barrier 0 — geometry availability (fresh re-acquire covers all rebuild paths). */
    {
        int local_ok = (topleaf_router_geometry_acquire() == 0) ? 1 : 0;
        int gate = 0;
        MPI_Allreduce(&local_ok, &gate, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(!gate) {
            if(ThisTask == 0) printf("[GX_ROUTE_ORACLE call=%d caller=%s SKIP: geometry unavailable on >=1 rank; broadcast reference used]\n", this_call, cname);
            goto done;
        }
    }

    /* Stage 1 (local) — routed-query CSR + per-dest send list, all guarded. */
    if((long long)nq * (long long)NTask > (long long)INT_MAX) aborted = 1;
    owners_cap   = aborted ? 1 : ((long)(nq > 0 ? nq : 1) * (long)NTask);
    route_off    = (int *)    malloc((size_t)(nq + 1) * sizeof(int));
    route_owners = (int *)    malloc((size_t)(owners_cap > 0 ? owners_cap : 1) * sizeof(int));
    q_pos        = (double *) malloc((size_t)(nq > 0 ? nq : 1) * 3 * sizeof(double));
    q_h          = (double *) malloc((size_t)(nq > 0 ? nq : 1) * sizeof(double));
    rq_sc        = (int *)    calloc((size_t)(NTask > 0 ? NTask : 1), sizeof(int));
    rq_sd        = (int *)    malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    rq_rc        = (int *)    malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    rq_rd        = (int *)    malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    if(!route_off || !route_owners || !q_pos || !q_h || !rq_sc || !rq_sd || !rq_rc || !rq_rd) aborted = 1;
    if(!aborted) {
        for(int i = 0; i < nq; i++) {
            q_pos[i*3+0] = local_queries[i].pos[0];
            q_pos[i*3+1] = local_queries[i].pos[1];
            q_pos[i*3+2] = local_queries[i].pos[2];
            q_h[i]       = local_queries[i].h;
        }
        int rc = topleaf_router_route_queries(q_pos, q_h, nq, supply_mask, oneway,
                                              periodic_flags, box_sizes,
                                              route_off, route_owners, owners_cap, ThisTask);
        if(rc < 0) aborted = 1;
    }
    if(!aborted) {
        for(int i = 0; i < nq; i++)
            for(int k = route_off[i]; k < route_off[i+1]; k++) rq_sc[route_owners[k]]++;
        rq_ts_ll = gx_oracle_checked_prefix(rq_sc, rq_sd, NTask);
        if(rq_ts_ll < 0) aborted = 1; else rq_total_send = (int)rq_ts_ll;
    }
    if(!aborted) {
        rq_send = (struct gx_query_t *) malloc((size_t)(rq_total_send > 0 ? rq_total_send : 1) * sizeof(struct gx_query_t));
        int *toff = (int *) malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
        if(!rq_send || !toff) { aborted = 1; free(toff); }
        else {
            memcpy(toff, rq_sd, (size_t)NTask * sizeof(int));
            for(int i = 0; i < nq; i++)
                for(int k = route_off[i]; k < route_off[i+1]; k++)
                    rq_send[toff[route_owners[k]]++] = local_queries[i];
            free(toff);
        }
    }
    /* Barrier 1 */
    MPI_Allreduce(&aborted, &bad_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(bad_any) {
        if(ThisTask == 0) printf("[GX_ROUTE_ORACLE call=%d caller=%s SKIP: route/alloc/int-range fail on >=1 rank; broadcast reference used]\n", this_call, cname);
        goto done;
    }

    /* Stage 2 — exchange routed queries. */
    MPI_Alltoall(rq_sc, 1, MPI_INT, rq_rc, 1, MPI_INT, MPI_COMM_WORLD);
    rq_tr_ll = gx_oracle_checked_prefix(rq_rc, rq_rd, NTask);
    if(rq_tr_ll < 0) aborted = 1; else rq_total_recv = (int)rq_tr_ll;
    if(!aborted) {
        rq_recv = (struct gx_query_t *) malloc((size_t)(rq_total_recv > 0 ? rq_total_recv : 1) * sizeof(struct gx_query_t));
        if(!rq_recv) aborted = 1;
    }
    /* Barrier 2 */
    MPI_Allreduce(&aborted, &bad_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(bad_any) {
        if(ThisTask == 0) printf("[GX_ROUTE_ORACLE call=%d caller=%s SKIP: recv alloc/int-range fail on >=1 rank; broadcast reference used]\n", this_call, cname);
        goto done;
    }
    gizmo_mpi_alltoallv_typed(rq_send, rq_sc, rq_sd, rq_recv, rq_rc, rq_rd,
                              sizeof(struct gx_query_t), MPI_COMM_WORLD);

    /* Stage 3 — walk routed queries per source segment -> matched_routed[src][p];
     * then build the per-dest routed ghost-id (home_rank,home_index) send list. */
    matched_r = (char *) calloc((size_t)(NTask > 0 ? NTask : 1) * (size_t)(num_pool > 0 ? num_pool : 1), sizeof(char));
    gid_sc    = (int *)  calloc((size_t)(NTask > 0 ? NTask : 1), sizeof(int));
    gid_sd    = (int *)  malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    gid_rc    = (int *)  malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    gid_rd    = (int *)  malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    if(!matched_r || !gid_sc || !gid_sd || !gid_rc || !gid_rd) aborted = 1;
    if(!aborted) {
        for(int s = 0; s < NTask; s++) {
            if(s == ThisTask) continue;
            char *mf = matched_r + (size_t)s * num_pool;
            for(int qi = 0; qi < rq_rc[s]; qi++) {
                const struct gx_query_t *q = &rq_recv[rq_rd[s] + qi];
                gx_walk_local_bvh(h_compact_xyzh, h_tiles, ntiles, h_pool, num_pool,
                                  h_pool_types, supply_mask, h_bvh, bvh_root,
                                  q->pos, q->h, search_mode, periodic_flags, box_sizes, mf, NULL);
            }
            int c = 0; for(int p = 0; p < num_pool; p++) if(mf[p]) c++;
            gid_sc[s] = c;
        }
        gid_ts_ll = gx_oracle_checked_prefix(gid_sc, gid_sd, NTask);
        if(gid_ts_ll < 0) aborted = 1; else gid_total_send = (int)gid_ts_ll;
    }
    if(!aborted) {
        gid_send = (struct gx_oracle_id *) malloc((size_t)(gid_total_send > 0 ? gid_total_send : 1) * sizeof(struct gx_oracle_id));
        int *toff = (int *) malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
        if(!gid_send || !toff) { aborted = 1; free(toff); }
        else {
            memcpy(toff, gid_sd, (size_t)NTask * sizeof(int));
            for(int s = 0; s < NTask; s++) {
                if(s == ThisTask) continue;
                char *mf = matched_r + (size_t)s * num_pool;
                for(int p = 0; p < num_pool; p++) {
                    if(!mf[p]) continue;
                    int off = toff[s]++;
                    gid_send[off].rank = ThisTask;
                    gid_send[off].idx  = h_pool[p];
                }
            }
            free(toff);
        }
    }
    /* Barrier 3 */
    MPI_Allreduce(&aborted, &bad_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(bad_any) {
        if(ThisTask == 0) printf("[GX_ROUTE_ORACLE call=%d caller=%s SKIP: walk/pack alloc/int-range fail on >=1 rank; broadcast reference used]\n", this_call, cname);
        goto done;
    }

    /* Stage 4 — exchange routed ghost-id pairs. */
    MPI_Alltoall(gid_sc, 1, MPI_INT, gid_rc, 1, MPI_INT, MPI_COMM_WORLD);
    gid_tr_ll = gx_oracle_checked_prefix(gid_rc, gid_rd, NTask);
    if(gid_tr_ll < 0) aborted = 1; else gid_total_recv = (int)gid_tr_ll;
    if(!aborted) {
        gid_recv = (struct gx_oracle_id *) malloc((size_t)(gid_total_recv > 0 ? gid_total_recv : 1) * sizeof(struct gx_oracle_id));
        bcast    = (struct gx_oracle_id *) malloc((size_t)(total_recv > 0 ? total_recv : 1) * sizeof(struct gx_oracle_id));
        if(!gid_recv || !bcast) aborted = 1;
    }
    /* Barrier 4 */
    MPI_Allreduce(&aborted, &bad_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(bad_any) {
        if(ThisTask == 0) printf("[GX_ROUTE_ORACLE call=%d caller=%s SKIP: pair recv alloc/int-range fail on >=1 rank; broadcast reference used]\n", this_call, cname);
        goto done;
    }
    gizmo_mpi_alltoallv_typed(gid_send, gid_sc, gid_sd, gid_recv, gid_rc, gid_rd,
                              sizeof(struct gx_oracle_id), MPI_COMM_WORLD);

    /* Stage 5 — sorted/deduped set compare: routed ⊆ broadcast always, so any
     * broadcast pair missing from routed is an UNDER-ROUTE. */
    skipped  = 0;
    n_routed = gx_oracle_sort_unique(gid_recv, gid_total_recv);
    for(int g = 0; g < total_recv; g++) { bcast[g].rank = ghost_home_rank_map[g]; bcast[g].idx = ghost_home_index_map[g]; }
    n_bcast  = gx_oracle_sort_unique(bcast, total_recv);
    local_mismatch = (n_routed != n_bcast) ? 1 : 0;
    {
        int ir = 0;
        for(int ib = 0; ib < n_bcast; ib++) {
            while(ir < n_routed && gx_oracle_id_cmp(&gid_recv[ir], &bcast[ib]) < 0) ir++;
            if(ir >= n_routed || gx_oracle_id_cmp(&gid_recv[ir], &bcast[ib]) != 0) {
                local_mismatch = 1;
                if(first_missing.rank < 0) first_missing = bcast[ib];
                break;
            }
        }
    }
    if(local_mismatch) {
        printf("[GX_ROUTE_ORACLE call=%d caller=%s rank=%d MISMATCH: routed=%d bcast=%d first_missing=(home_rank=%d,home_idx=%d)]\n",
               this_call, cname, ThisTask, n_routed, n_bcast, first_missing.rank, first_missing.idx);
        fflush(stdout);
    }
    MPI_Allreduce(&local_mismatch, &mismatch_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(mismatch_any) {
        gizmo_request_controlled_stop(7704, "ghost route oracle: routed ghost set != broadcast (under-route)",
                                      __FILE__, __LINE__, __FUNCTION__);
    } else if(ThisTask == 0) {
        printf("[GX_ROUTE_ORACLE call=%d caller=%s OK: routed==broadcast ghost set]\n", this_call, cname);
        fflush(stdout);
    }

done:
    free(bcast);
    free(gid_recv); free(gid_rd); free(gid_rc); free(gid_send); free(gid_sd); free(gid_sc);
    free(matched_r);
    free(rq_recv); free(rq_rd); free(rq_rc); free(rq_send); free(rq_sd); free(rq_sc);
    free(q_pos); free(q_h); free(route_off); free(route_owners);

    if(ThisTask == 0) {
        printf("[GX_ROUTE_ORACLE_TIME call=%d caller=%s oracle=%.4f s %s (diagnostic; GX_RD_TIME is inflated while the oracle is on)]\n",
               this_call, cname, timediff(t_oracle_start, my_second()), skipped ? "SKIPPED" : "RAN");
        fflush(stdout);
    }
    /* Drain the controlled-stop request collectively (all ranks reach here). */
    gizmo_exit_bad_stop_if_requested("ghost_exchange:route_oracle");
}

/* Broadcast discovery walk. Walks every remote rank's
 * queries (from the Allgatherv'd all_queries) against the pre-built local supply
 * snapshot and returns the per-peer match bitmask matched[t*num_pool+p] that the
 * shared pack/install steps consume.  Pure local walk over already-exchanged
 * queries -- no routing, no collectives.  CALLER OWNS the returned buffer (free
 * it).  This is the broadcast `matched` producer; the routed producer below adds
 * a routed producer with the identical return layout so the install path stays shared. */
static char *compute_matched_broadcast(
    const struct gx_query_t *all_queries, const int *q_disps, const int *all_q_counts,
    const float *h_compact_xyzh, const sfc_tile_t *h_tiles, int ntiles,
    const int *h_pool, int num_pool, const int *h_pool_types, unsigned int supply_mask,
    const tile_bvh_node_t *h_bvh, int bvh_root, int search_mode,
    const int periodic_flags[3], const double box_sizes[3])
{
    /* Per-peer match bitmask over pool indices (dedup multiple queries → one ghost). */
    char *matched = (char *) calloc((size_t)NTask * (size_t)(num_pool > 0 ? num_pool : 1), sizeof(char));
    for(int t = 0; t < NTask; t++) {
        if(t == ThisTask) continue;
        int q_start = q_disps[t];
        int q_count = all_q_counts[t];
        char *match_for_t = matched + (size_t)t * (size_t)num_pool;
        for(int qi = 0; qi < q_count; qi++) {
            const struct gx_query_t *q = &all_queries[q_start + qi];
            /* q->h already includes safety_factor (set at query-build time).
             * compact_xyzh[*4+3] also includes safety_factor. So leaf r² check
             * uses inflated radii on both sides of max(h_q, h_j). */
            gx_walk_local_bvh(h_compact_xyzh, h_tiles, ntiles, h_pool, num_pool,
                              h_pool_types, supply_mask,
                              h_bvh, bvh_root,
                              q->pos, q->h, search_mode,
                              periodic_flags, box_sizes,
                              match_for_t, NULL);
        }
    }
    return matched;
}

/* Lazy, idempotent collective broadcast of the local query lists to all ranks
 * (the legacy Step-2 Allgather/Allgatherv).  Safe to call multiple times: no-op
 * once *available.  MUST be called collectively (all ranks together).  Fills the
 * three malloc'd arrays (caller frees) + total_queries.  Enables the late
 * fallback: if routed discovery fails AFTER Step 2 was skipped, all ranks return
 * to the same point and run this collectively before the broadcast walk. */
static void ensure_broadcast_queries(int *available,
                                     const struct gx_query_t *local_queries, int n_local_queries,
                                     int **all_q_counts_io, int **q_disps_io,
                                     struct gx_query_t **all_queries_io, int *total_queries_io)
{
    if(*available) return;
    int *all_q_counts = (int *) malloc(NTask * sizeof(int));
    MPI_Allgather(&n_local_queries, 1, MPI_INT, all_q_counts, 1, MPI_INT, MPI_COMM_WORLD);
    int *q_disps = (int *) malloc(NTask * sizeof(int));
    int total_queries = 0;
    for(int t = 0; t < NTask; t++) { q_disps[t] = total_queries; total_queries += all_q_counts[t]; }
    struct gx_query_t *all_queries = (struct gx_query_t *)
        malloc((size_t)(total_queries > 0 ? total_queries : 1) * sizeof(struct gx_query_t));
    int *q_byte_counts = (int *) malloc(NTask * sizeof(int));
    int *q_byte_disps  = (int *) malloc(NTask * sizeof(int));
    for(int t = 0; t < NTask; t++) {
        q_byte_counts[t] = all_q_counts[t] * (int)sizeof(struct gx_query_t);
        q_byte_disps[t]  = q_disps[t]      * (int)sizeof(struct gx_query_t);
    }
    MPI_Allgatherv(local_queries, n_local_queries * (int)sizeof(struct gx_query_t), MPI_BYTE,
                   all_queries, q_byte_counts, q_byte_disps, MPI_BYTE, MPI_COMM_WORLD);
    free(q_byte_counts); free(q_byte_disps);
    *all_q_counts_io = all_q_counts; *q_disps_io = q_disps;
    *all_queries_io = all_queries;   *total_queries_io = total_queries;
    *available = 1;
}

/* Routed `matched` producer: routes each local query to its overlapping
 * top-leaf owners, Alltoallv's the routed queries, and walks ONLY the received
 * queries per source rank -> matched[t*num_pool+p] (IDENTICAL layout to
 * compute_matched_broadcast).  Geometry must already be acquired+validated by the
 * caller's collective gate.  Fully fail-closed: on geometry/INT-overflow/alloc
 * failure on ANY rank, all ranks return NULL together (via Allreduce barriers) so
 * the caller falls back to broadcast collectively.  ONEWAY only (band-free).
 * CALLER OWNS the returned buffer (free it). */
static char *compute_matched_routed(
    const struct gx_query_t *local_queries, int n_local_queries,
    const float *h_compact_xyzh, const sfc_tile_t *h_tiles, int ntiles,
    const int *h_pool, int num_pool, const int *h_pool_types, unsigned int supply_mask,
    const tile_bvh_node_t *h_bvh, int bvh_root, int search_mode,
    const int periodic_flags[3], const double box_sizes[3],
    int use_hier,
    double *t_route_construct, double *t_route_alltoallv, double *t_route_walk,
    long *fanout_owner_sum, long *fanout_owner_max, int *total_recv_out,
    long *diag_pairs, long *diag_pairs_nonzero, long *diag_hit_sum, long *diag_hit_max)
{
    if(fanout_owner_sum) *fanout_owner_sum = 0;
    if(fanout_owner_max) *fanout_owner_max = 0;
    if(total_recv_out)   *total_recv_out   = 0;
    if(diag_pairs)         *diag_pairs         = 0;
    if(diag_pairs_nonzero) *diag_pairs_nonzero = 0;
    if(diag_hit_sum)       *diag_hit_sum       = 0;
    if(diag_hit_max)       *diag_hit_max       = 0;
    int *route_off = NULL, *route_owners = NULL;
    double *q_pos = NULL, *q_h = NULL;
    int *rq_sc = NULL, *rq_sd = NULL, *rq_rc = NULL, *rq_rd = NULL;
    struct gx_query_t *rq_send = NULL, *rq_recv = NULL;
    char *matched = NULL;
    int  oneway = (search_mode == NGB_SEARCH_ONEWAY) ? 1 : 0;
    int  nq = n_local_queries;
    int  aborted = 0, bad_any = 0;
    long owners_cap = 0;
    long long rq_ts = 0, rq_tr = 0;
    int  rq_total_send = 0, rq_total_recv = 0;
    /* Split route timing: construct / query-Alltoallv / walk, so the FIRE
     * perf gate proves the flat O(nq x NTopleaves) cost actually leaves the
     * construct bucket (and doesn't reappear elsewhere).  Decls at top so the
     * goto-fail never jumps over an initialised automatic. */
    double tc0 = my_second(), ta0 = 0.0, tw0 = 0.0;
    if(t_route_construct) *t_route_construct = 0.0;
    if(t_route_alltoallv) *t_route_alltoallv = 0.0;
    if(t_route_walk)      *t_route_walk      = 0.0;

    /* Stage 1 (local): routed-query CSR + per-dest send list (guarded). */
    if((long long)nq * (long long)NTask > (long long)INT_MAX) aborted = 1;
    owners_cap   = aborted ? 1 : ((long)(nq > 0 ? nq : 1) * (long)NTask);
    route_off    = (int *)    malloc((size_t)(nq + 1) * sizeof(int));
    route_owners = (int *)    malloc((size_t)(owners_cap > 0 ? owners_cap : 1) * sizeof(int));
    q_pos        = (double *) malloc((size_t)(nq > 0 ? nq : 1) * 3 * sizeof(double));
    q_h          = (double *) malloc((size_t)(nq > 0 ? nq : 1) * sizeof(double));
    rq_sc        = (int *)    calloc((size_t)(NTask > 0 ? NTask : 1), sizeof(int));
    rq_sd        = (int *)    malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    rq_rc        = (int *)    malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    rq_rd        = (int *)    malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
    if(!route_off || !route_owners || !q_pos || !q_h || !rq_sc || !rq_sd || !rq_rc || !rq_rd) aborted = 1;
    if(!aborted) {
        for(int i = 0; i < nq; i++) {
            q_pos[i*3+0] = local_queries[i].pos[0];
            q_pos[i*3+1] = local_queries[i].pos[1];
            q_pos[i*3+2] = local_queries[i].pos[2];
            q_h[i]       = local_queries[i].h;
        }
        /* Hierarchical TopNodes descent (production); the flat constructor is the
         * reference for the flat-vs-hier oracle. */
        int rc = use_hier
            ? topleaf_router_route_queries_hier(q_pos, q_h, nq, supply_mask, oneway,
                                                periodic_flags, box_sizes,
                                                route_off, route_owners, owners_cap, ThisTask)
            : topleaf_router_route_queries     (q_pos, q_h, nq, supply_mask, oneway,
                                                periodic_flags, box_sizes,
                                                route_off, route_owners, owners_cap, ThisTask);
        if(rc < 0) aborted = 1;
    }
    if(!aborted) {
        for(int i = 0; i < nq; i++)
            for(int k = route_off[i]; k < route_off[i+1]; k++) rq_sc[route_owners[k]]++;
        rq_ts = gx_oracle_checked_prefix(rq_sc, rq_sd, NTask);
        if(rq_ts < 0) aborted = 1; else rq_total_send = (int)rq_ts;
        /* query->owner fanout (H4c perf metric): owners per query = route_off deltas. */
        if(fanout_owner_sum || fanout_owner_max) {
            long fsum = 0, fmax = 0;
            for(int i = 0; i < nq; i++) { long f = route_off[i+1] - route_off[i]; fsum += f; if(f > fmax) fmax = f; }
            if(fanout_owner_sum) *fanout_owner_sum = fsum;
            if(fanout_owner_max) *fanout_owner_max = fmax;
        }
    }
    if(!aborted) {
        rq_send = (struct gx_query_t *) malloc((size_t)(rq_total_send > 0 ? rq_total_send : 1) * sizeof(struct gx_query_t));
        int *toff = (int *) malloc((size_t)(NTask > 0 ? NTask : 1) * sizeof(int));
        if(!rq_send || !toff) { aborted = 1; free(toff); }
        else {
            memcpy(toff, rq_sd, (size_t)NTask * sizeof(int));
            for(int i = 0; i < nq; i++)
                for(int k = route_off[i]; k < route_off[i+1]; k++)
                    rq_send[toff[route_owners[k]]++] = local_queries[i];
            free(toff);
        }
    }
    MPI_Allreduce(&aborted, &bad_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(bad_any) goto fail;
    if(t_route_construct) *t_route_construct = timediff(tc0, my_second());

    /* Stage 2: exchange routed queries. */
    ta0 = my_second();
    MPI_Alltoall(rq_sc, 1, MPI_INT, rq_rc, 1, MPI_INT, MPI_COMM_WORLD);
    rq_tr = gx_oracle_checked_prefix(rq_rc, rq_rd, NTask);
    if(rq_tr < 0) aborted = 1; else { rq_total_recv = (int)rq_tr; if(total_recv_out) *total_recv_out = rq_total_recv; }
    if(!aborted) {
        rq_recv = (struct gx_query_t *) malloc((size_t)(rq_total_recv > 0 ? rq_total_recv : 1) * sizeof(struct gx_query_t));
        if(!rq_recv) aborted = 1;
    }
    MPI_Allreduce(&aborted, &bad_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(bad_any) goto fail;
    gizmo_mpi_alltoallv_typed(rq_send, rq_sc, rq_sd, rq_recv, rq_rc, rq_rd,
                              sizeof(struct gx_query_t), MPI_COMM_WORLD);
    if(t_route_alltoallv) *t_route_alltoallv = timediff(ta0, my_second());

    /* Stage 3: walk routed queries per source segment -> matched[src][p]. */
    tw0 = my_second();
    matched = (char *) calloc((size_t)(NTask > 0 ? NTask : 1) * (size_t)(num_pool > 0 ? num_pool : 1), sizeof(char));
    if(!matched) aborted = 1;
    if(!aborted) {
        /* Over-route diagnostic (H4c): per received (query,owner=this-rank) pair, count
         * EXACT matches so we can report what fraction of routed pairs return zero ghosts
         * = the conservative band's wasted routing.  Only when the caller asks (oracle). */
        int  want_diag = (diag_pairs_nonzero != NULL);
        long d_pairs = 0, d_pairs_nz = 0, d_hit_sum = 0, d_hit_max = 0;
        for(int s = 0; s < NTask; s++) {
            if(s == ThisTask) continue;
            char *mf = matched + (size_t)s * num_pool;
            for(int qi = 0; qi < rq_rc[s]; qi++) {
                const struct gx_query_t *q = &rq_recv[rq_rd[s] + qi];
                if(want_diag) {
                    long hits = 0;
                    gx_walk_local_bvh(h_compact_xyzh, h_tiles, ntiles, h_pool, num_pool,
                                      h_pool_types, supply_mask, h_bvh, bvh_root,
                                      q->pos, q->h, search_mode, periodic_flags, box_sizes, mf, &hits);
                    d_pairs++;
                    if(hits > 0) { d_pairs_nz++; d_hit_sum += hits; if(hits > d_hit_max) d_hit_max = hits; }
                } else {
                    gx_walk_local_bvh(h_compact_xyzh, h_tiles, ntiles, h_pool, num_pool,
                                      h_pool_types, supply_mask, h_bvh, bvh_root,
                                      q->pos, q->h, search_mode, periodic_flags, box_sizes, mf, NULL);
                }
            }
        }
        if(diag_pairs)         *diag_pairs         = d_pairs;
        if(diag_pairs_nonzero) *diag_pairs_nonzero = d_pairs_nz;
        if(diag_hit_sum)       *diag_hit_sum       = d_hit_sum;
        if(diag_hit_max)       *diag_hit_max       = d_hit_max;
    }
    MPI_Allreduce(&aborted, &bad_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(bad_any) goto fail;
    if(t_route_walk) *t_route_walk = timediff(tw0, my_second());

    free(rq_recv); free(rq_rd); free(rq_rc); free(rq_send); free(rq_sd); free(rq_sc);
    free(q_pos); free(q_h); free(route_off); free(route_owners);
    return matched;

fail:
    free(matched);
    free(rq_recv); free(rq_rd); free(rq_rc); free(rq_send); free(rq_sd); free(rq_sc);
    free(q_pos); free(q_h); free(route_off); free(route_owners);
    return NULL;
}

/* Walk-export routed producer — the discovery path for every spec that passes
 * gx_walk_export_eligible().  Sender: per local query mode_b_walk_and_export -> per-peer
 * NodeList -> fixed-size envelopes -> Alltoallv.  Receiver: mode_b_walk_from_start_nodes
 * (resume from the exported NodeList) -> gx_pair_accept (the ghost-exchange SSOT predicate)
 * -> matched[t*num_pool+p] bitmap, the same layout the shared Steps 4-6 install consume.
 * MODE-GENERIC (search_mode is passed through): this is the install target for both search
 * modes, so ONEWAY and SYMMETRIC discover on ONE substrate rather than two.
 *
 * COLLECTIVE-SAFE (C MPI buffers): every C allocation preceding a collective — the index
 * arrays before the Alltoall, the envelope buffers before the Alltoallv, the matched bitmap
 * — is Allreduce-checked, so a NULL on ANY rank makes ALL ranks return the same status and
 * ranks never diverge across a collective.  NOT covered: the C++ containers
 * (ModeBExportSink, the per-peer envelope vectors) throw std::bad_alloc rank-locally rather
 * than returning NULL, so an allocation failure there aborts that rank instead of returning
 * a uniform status.  Returns a caller-owned matched bitmap on GX_WALK_EXPORT_OK, else NULL.
 *
 * THREADING (host): the sender-query loop and the received-envelope walk are both
 * `omp parallel for` — per-thread export sink and per-thread send buffers on the sender,
 * pre-sized per-envelope candidate slots on the receiver, so each thread writes only its own
 * index.  The walker's lazy node drift is the one shared mutation and is serialized under
 * critical(_modebdrift_).  Merge, accept and bitmap set run serially afterwards, which keeps
 * the resulting SET order-independent and therefore deterministic across thread counts.
 *
 * MEMORY SHAPE: the threaded receiver materializes one candidate list per received envelope
 * before the serial accept pass.  That is bounded by the exported envelope volume, which is
 * measured sparse; a pathologically clustered geometry with very large tot_r would grow it,
 * so it is a known scaling watch-point rather than a fixed bound. */
enum { GX_WALK_EXPORT_OK = 0, GX_WALK_EXPORT_UNAVAILABLE = 1, GX_WALK_EXPORT_ALLOC_FAIL = 2 };
struct gx_walk_export_result {
    int status;   /* GX_WALK_EXPORT_OK / _UNAVAILABLE / _ALLOC_FAIL, uniform across ranks */
};

static char *compute_matched_walk_export(
    const struct ghost_exchange_spec_t *spec,
    const struct gx_query_t *local_queries, int n_local_queries,
    int num_pool, unsigned int supply_mask, int search_mode,
    const int periodic_flags[3], const double box_sizes[3],
    struct gx_walk_export_result *res)
{
    if(res) memset(res, 0, sizeof(*res));
    /* (a) tree availability — collective all-or-none (a rank-local skip would deadlock the
     * envelope Alltoallv below / the caller's compare Allreduce). */
    int ok_local = (All.MaxPart > 0 && Nodes != NULL && Nextnode != NULL) ? 1 : 0;
    int ok_all = 0;
    MPI_Allreduce(&ok_local, &ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if(!ok_all) { if(res) res->status = GX_WALK_EXPORT_UNAVAILABLE; return NULL; }

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
                                   /*cand_out=*/NULL, map, sink, spec->j_radius_scale, /*drift_ctr=*/NULL);
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
        if(res) res->status = GX_WALK_EXPORT_ALLOC_FAIL; return NULL;
    } else {
        int a_ok_local = 1, a_ok_all = 0;
        MPI_Allreduce(&a_ok_local, &a_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(!a_ok_all) { free(sc); free(sd); free(rc); free(rd); if(res) res->status = GX_WALK_EXPORT_ALLOC_FAIL; return NULL; }
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
            return NULL;
        }
    }
    struct gx_export_envelope_t *sendbuf = (struct gx_export_envelope_t *) malloc((size_t)(tot_s > 0 ? tot_s : 1) * sizeof(struct gx_export_envelope_t));
    struct gx_export_envelope_t *recv    = (struct gx_export_envelope_t *) malloc((size_t)(tot_r > 0 ? tot_r : 1) * sizeof(struct gx_export_envelope_t));
    /* buffer alloc — collective before the Alltoallv (a NULL recv would fault the collective). */
    int buf_ok_local = (sendbuf && recv) ? 1 : 0, buf_ok_all = 0;
    MPI_Allreduce(&buf_ok_local, &buf_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if(!buf_ok_all) { free(sendbuf); free(recv); free(sc); free(sd); free(rc); free(rd);
                      if(res) res->status = GX_WALK_EXPORT_ALLOC_FAIL; return NULL; }
    for(int t = 0; t < NTask; t++)
        if(sc[t] > 0) memcpy(sendbuf + sd[t], send[t].data(), (size_t)sc[t] * sizeof(struct gx_export_envelope_t));
    gizmo_mpi_alltoallv_typed(sendbuf, sc, sd, recv, rc, rd,
                              sizeof(struct gx_export_envelope_t), MPI_COMM_WORLD);
    free(sendbuf); free(sc); free(sd);

    /* (d) matched bitmap — collective before the RECEIVER walk (the caller reduces over the
     * compare, so a NULL here on one rank would diverge the caller's Allreduce). */
    char *matched_walk_export = (char *) calloc((size_t)NTask * (size_t)(num_pool > 0 ? num_pool : 1), 1);
    int m_ok_local = (matched_walk_export != NULL) ? 1 : 0, m_ok_all = 0;
    MPI_Allreduce(&m_ok_local, &m_ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if(!m_ok_all) { free(matched_walk_export); free(recv); free(rc); free(rd);
                    if(res) res->status = GX_WALK_EXPORT_ALLOC_FAIL; return NULL; }

    /* (e) RECEIVER: bounded resume-walk from the exported NodeList -> SSOT accept -> bitmap.
     * HOST-THREADED: the WALK (dominant cost) runs per received envelope into a pre-sized per-envelope
     * cand slot — each thread writes ONLY its own index (no shared write), walker race-safe.  The
     * ACCEPT + bitmap set run SERIALLY afterward (cheap): this keeps the matched_walk_export bitmap race-free
     * AND the installed set order-independent (a bit is set or not), so the set is deterministic. */
    {
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
                                         cvk, spec->j_radius_scale, NULL);
        }
        for(int t = 0; t < NTask; t++) {
            if(t == ThisTask) continue;
            char *mf = matched_walk_export + (size_t)t * num_pool;
            for(int r = 0; r < rc[t]; r++) {
                const long k = (long)rd[t] + r;
                const struct gx_export_envelope_t *e = &recv[k];
                const std::vector<int> &cvk = per_recv_cands[k];
                for(size_t c = 0; c < cvk.size(); c++) {
                    int j = cvk[c];
                    if(j < 0 || j >= g_glt_cache.NumPart_when_built) continue;
                    int pp = g_glt_cache.j_to_pool ? g_glt_cache.j_to_pool[j] : -1;
                    if(pp < 0 || pp >= num_pool) continue;
                    /* Reach comes from THIS caller's spec, not from whatever policy
                     * the cached pool happened to be built under.  The sender walk
                     * above already uses the spec, so taking it from the cache here
                     * would let a caller inherit another caller's j-side reach now
                     * that pool membership is reused across differing radius
                     * policies.  Only the tile/BVH/compact geometry may use the
                     * build-time policy, because its leaf h was baked with it. */
                    double hj_dbl = gx_policy_scaled_h(j, spec->radius_policy,
                                                       spec->j_radius_scale,
                                                       spec->safety_factor);
                    if(gx_pair_accept(e->pos, e->h, P[j].Pos[0], P[j].Pos[1], P[j].Pos[2],
                                      hj_dbl, search_mode, periodic_flags, box_sizes)) {
                        mf[pp] = 1;   /* idempotent: set semantics, duplicates are a no-op */
                    }
                }
            }
        }
    }
    free(recv); free(rc); free(rd);
    if(res) res->status = GX_WALK_EXPORT_OK;
    return matched_walk_export;
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
    double t_step1_start = my_second();
    int n_local_queries = 0;
    struct gx_query_t *local_queries = NULL;
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
            q++;
        }
    }

    double t_step1 = timediff(t_step1_start, my_second());
    /* === Step 2: LAZY query distribution === The broadcast Allgather/
     * Allgatherv now runs on demand via ensure_broadcast_queries(): SKIPPED in
     * routed-production mode, run up-front for the oracle / when routed is
     * unavailable, and as a collective LATE fallback if routed fails after the
     * skip.  Declarations only here; the gather (if any) happens after the
     * supply-snapshot build below, alongside the routed-vs-broadcast selection. */
    int   *all_q_counts = NULL;
    int   *q_disps      = NULL;
    struct gx_query_t *all_queries = NULL;
    int    total_queries = 0;
    int    bcast_queries_available = 0;
    double t_step2 = 0.0;

    /* === Step 3: per-rank, walk local BVH against each remote rank's queries ===
     *
     * Build a host-side SFC tile + BVH index over the supply-mask-filtered pool.
     * For each remote query, walk the BVH (O(log N + matches) per query, vs the
     * O(N) brute-force scan that had been a tiny-N proof-of-concept). Per-particle
     * acceptance at leaves uses the EXACT predicate (ONEWAY: r²<h_q²; SYMMETRIC:
     * r²<max(h_q,h_j)²) — same predicate the kernel applies later. */
    double t_step3_start = my_second();

    /* SHARED-TREE build: build over GHOST_TYPE_ALL (all types with mass>0), not
     * just the caller's supply_mask. The walker's per-type hmax filter +
     * per-particle leaf Type-vs-supply_mask check delivers the same imports
     * as the old per-supply-mask build.
     *
     * Bucket 3 (SIDX overlay): the local tree (tiles, pool, bvh,
     * compact_xyzh, pool_types) is cached across calls within a step.
     * Cache invalidated via ghost_exchange_local_tree_invalidate_*()
     * hooks at drift / domain_decomp boundaries (run.cc). Within a step
     * the pool is stable, so 2nd..Nth calls skip this whole stanza. */
    sfc_tile_t *h_tiles = NULL;
    int *h_pool = NULL;
    int num_pool = 0;
    int ntiles = 0;
    tile_bvh_node_t *h_bvh = NULL;
    int bvh_nnodes = 0;
    int bvh_root = 0;
    float *h_compact_xyzh = NULL;
    int *h_pool_types = NULL;
    int from_cache = 0;
    /* All particle types are eligible as ghost sources. */
    unsigned int desired_pool_mask = GHOST_TYPE_ALL;
    /* Tile/BVH/compact geometry is built ONLY for a call that actually walks it.
     * The walk-export producer discovers through the tree and reads nothing
     * position-dependent from this cache (it needs pool membership to index the
     * matched set, and takes its supply band from Extnodes), so when it is the
     * authority the geometry would be built and never read.  Structural, keyed on
     * the same predicate as producer selection — no caller names. */
    /* ONE predicate drives BOTH what the cache holds and which producers may run
     * below: geometry is skipped exactly when no consumer of it can execute.  The
     * diag path is part of that condition — it re-enables the broadcast walk, which
     * reads this geometry — so a second, weaker copy of this test here would let a
     * diag run walk NULL tiles.  Keep this as the single definition. */
    const int oracle_on  = ghost_route_oracle_enabled();
    const int walk_export_only    = gx_walk_export_eligible(spec)
                           && !gizmo_nlr_phase0_diag_enabled()
                           && !oracle_on;
    const unsigned int wanted_caps = GX_POOL_IDENTITY | (walk_export_only ? 0u : GX_POOL_GEOMETRY);
    /* TWO validities, because the cache holds two payloads with different
     * dependencies.  IDENTITY (pool, j_to_pool, num_pool) is membership and order:
     * it depends on {NumPart, type mask, epoch} and on nothing positional, so it
     * survives drift and survives a caller arriving with a different radius policy.
     * GEOMETRY (tiles, bvh, compact_xyzh, pool_types) additionally depends on
     * {safety, radius_policy, j_radius_scale} and on live positions/h, so it is
     * rebuilt on a policy change and REFIT (not rebuilt) on drift.
     *
     * Keeping these separate is the point: rebuilding the identity map costs a
     * full pass over the local particles, and callers within a step differ in
     * radius policy far more often than the particle set changes.
     *
     * Mask is compared for EXACT equality, not coverage. Reuse across a narrowed
     * mask is unproven here — a narrower request would also make in-place Type
     * changes membership-relevant, which the epoch does not track — so anything
     * other than the all-types pool falls through to a full rebuild. */
    const int mask_reusable = (desired_pool_mask == GHOST_TYPE_ALL);
    const int identity_valid = (g_glt_cache.valid
                       && (g_glt_cache.caps & GX_POOL_IDENTITY)
                       && g_glt_cache.pool && g_glt_cache.j_to_pool
                       && mask_reusable
                       && g_glt_cache.eligible_type_mask_when_built == desired_pool_mask
                       && g_glt_cache.NumPart_when_built == NumPart
                       && g_glt_cache.identity_epoch_when_built == g_supply_identity_epoch);
    const int geometry_valid = (identity_valid
                       && (g_glt_cache.caps & GX_POOL_GEOMETRY)
                       && g_glt_cache.safety_factor_when_built == safety_factor
                       && g_glt_cache.radius_policy_when_built == spec->radius_policy
                       && g_glt_cache.j_radius_scale_when_built == spec->j_radius_scale);
    const int want_geometry = (wanted_caps & GX_POOL_GEOMETRY) ? 1 : 0;
    int cache_match = identity_valid && (!want_geometry || geometry_valid);
    int geometry_rebuilt = 0;

    /* Identity still good, geometry stale or absent: rebuild geometry ONLY, over
     * the pool we already hold.  This is the case a single fused key could not
     * express, and it is the common one — successive callers in a step share the
     * particle set and differ only in radius policy. */
    if(identity_valid && want_geometry && !geometry_valid) {
        glt_cache_free_geometry_();
        sfc_tile_t *g_tiles = NULL;
        int g_ntiles = build_sfc_tiles_from_pool(P, g_glt_cache.pool, g_glt_cache.num_pool,
                                                 TILE_TARGET_SIZE, &g_tiles,
                                                 spec->radius_policy,
                                                 spec->j_radius_scale * safety_factor);
        tile_bvh_node_t *g_bvh = NULL;
        int g_bvh_nnodes = build_tile_bvh(g_tiles, g_ntiles, &g_bvh);
        size_t gz_tiles   = (size_t)(g_ntiles > 0 ? g_ntiles : 1) * sizeof(sfc_tile_t);
        size_t gz_bvh     = (size_t)(g_bvh_nnodes > 0 ? g_bvh_nnodes : 1) * sizeof(tile_bvh_node_t);
        size_t gz_compact = (size_t)(g_glt_cache.num_pool > 0 ? g_glt_cache.num_pool : 1) * 4 * sizeof(float);
        size_t gz_types   = (size_t)(g_glt_cache.num_pool > 0 ? g_glt_cache.num_pool : 1) * sizeof(int);
        sfc_tile_t      *gc_tiles   = (sfc_tile_t *)      malloc(gz_tiles);
        tile_bvh_node_t *gc_bvh     = (tile_bvh_node_t *) malloc(gz_bvh);
        float           *gc_compact = (float *)           malloc(gz_compact);
        int             *gc_types   = (int *)             malloc(gz_types);
        if(!gc_tiles || !gc_bvh || !gc_compact || !gc_types) {
            /* Same fail-closed shape as the full rebuild below: this producer is the
             * only supplier, so drop the whole entry and let the full path re-decide
             * rather than publishing a half-built one. */
            free(gc_tiles); free(gc_bvh); free(gc_compact); free(gc_types);
            if(g_bvh)   myfree(g_bvh);
            if(g_tiles) myfree(g_tiles);
            glt_cache_free();
        } else {
            if(g_ntiles > 0)      memcpy(gc_tiles, g_tiles, (size_t)g_ntiles * sizeof(sfc_tile_t));
            if(g_bvh_nnodes > 0)  memcpy(gc_bvh,   g_bvh,   (size_t)g_bvh_nnodes * sizeof(tile_bvh_node_t));
            for(int p = 0; p < g_glt_cache.num_pool; p++) {
                int j = g_glt_cache.pool[p];
                gc_compact[p*4+0] = (float)P[j].Pos[0];
                gc_compact[p*4+1] = (float)P[j].Pos[1];
                gc_compact[p*4+2] = (float)P[j].Pos[2];
                gc_compact[p*4+3] = (float)gx_policy_scaled_h(j, spec->radius_policy,
                                                              spec->j_radius_scale, safety_factor);
                gc_types[p] = (int)P[j].Type;
            }
            if(g_bvh)   myfree(g_bvh);
            if(g_tiles) myfree(g_tiles);
            g_glt_cache.tiles        = gc_tiles;
            g_glt_cache.bvh          = gc_bvh;
            g_glt_cache.compact_xyzh = gc_compact;
            g_glt_cache.pool_types   = gc_types;
            g_glt_cache.ntiles       = g_ntiles;
            g_glt_cache.bvh_nnodes   = g_bvh_nnodes;
            g_glt_cache.bvh_root     = g_bvh_nnodes - 1;
            g_glt_cache.safety_factor_when_built  = safety_factor;
            g_glt_cache.radius_policy_when_built  = spec->radius_policy;
            g_glt_cache.j_radius_scale_when_built = spec->j_radius_scale;
            g_glt_cache.Ti_when_built = All.Ti_Current;
            g_glt_cache.needs_refit = 0;
            g_glt_cache.caps |= GX_POOL_GEOMETRY;
            /* Geometry was just seeded from current P[]; marks predating it are
             * obsolete, exactly as after a full build. */
            g_glt_dirty_clear_();
            cache_match = 1;
            geometry_rebuilt = 1;
        }
    }
    if(cache_match) {
        h_tiles        = g_glt_cache.tiles;
        h_pool         = g_glt_cache.pool;
        num_pool       = g_glt_cache.num_pool;
        ntiles         = g_glt_cache.ntiles;
        h_bvh          = g_glt_cache.bvh;
        bvh_nnodes     = g_glt_cache.bvh_nnodes;
        bvh_root       = g_glt_cache.bvh_root;
        h_compact_xyzh = g_glt_cache.compact_xyzh;
        h_pool_types   = g_glt_cache.pool_types;
        from_cache = 1;
        /* A geometry-only rebuild reused the pool but did real build work, so it is
         * not a hit. Counting it as one would overstate the cache and hide the very
         * cost this split exists to measure. */
        if(geometry_rebuilt) g_glt_cache_misses++; else g_glt_cache_hits++;
        /* Refit refreshes position-dependent geometry only; an identity-only
         * entry has none and its membership does not move with the particles. */
        if((g_glt_cache.caps & GX_POOL_GEOMETRY)
           && (g_glt_cache.needs_refit || g_glt_cache.Ti_when_built != All.Ti_Current)) {
            glt_cache_refit_from_particles();
        }
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

        /* Fresh build via existing mymalloc path; we copy the result into
         * malloc-backed cache buffers so it can outlive this function frame
         * without violating mymalloc LIFO ordering. */
        sfc_tile_t *tmp_tiles = NULL;
        int *tmp_pool = NULL;
        int tmp_num_pool = 0;
        /* Tile bands and leaf compact h are computed from the IDENTICAL
         * supply-side reach: nlr_particle_symmetric_radius(P[j], spec->radius_policy)
         * * spec->j_radius_scale * safety_factor.  build_sfc_tiles bakes the
         * scale into tile->hmax / tile->hmax_by_type[]; the loop below bakes it
         * into c_compact[p*4+3] via gx_policy_scaled_h.  No leaf-vs-band scale
         * gap, no BVH-prune-misses-pairs failure mode. */
        int tmp_ntiles = 0;
        tile_bvh_node_t *tmp_bvh = NULL;
        int tmp_bvh_nnodes = 0;
        if(wanted_caps & GX_POOL_GEOMETRY) {
            tmp_ntiles = build_sfc_tiles(P, NumPart, (int)desired_pool_mask, TILE_TARGET_SIZE,
                                         &tmp_tiles, &tmp_pool, &tmp_num_pool,
                                         spec->radius_policy,
                                         spec->j_radius_scale * safety_factor);
            tmp_bvh_nnodes = build_tile_bvh(tmp_tiles, tmp_ntiles, &tmp_bvh);
        } else {
            /* Membership only — same selection, none of the position-dependent work. */
            tmp_num_pool = build_sfc_supply_pool(P, NumPart, (int)desired_pool_mask, &tmp_pool);
        }
        int tmp_bvh_root = tmp_bvh_nnodes - 1;

        /* Allocate persistent cache buffers + copy. */
        size_t sz_tiles   = (size_t)(tmp_ntiles > 0 ? tmp_ntiles : 1) * sizeof(sfc_tile_t);
        size_t sz_pool    = (size_t)(tmp_num_pool > 0 ? tmp_num_pool : 1) * sizeof(int);
        size_t sz_bvh     = (size_t)(tmp_bvh_nnodes > 0 ? tmp_bvh_nnodes : 1) * sizeof(tile_bvh_node_t);
        size_t sz_compact = (size_t)(tmp_num_pool > 0 ? tmp_num_pool : 1) * 4 * sizeof(float);
        size_t sz_types   = (size_t)(tmp_num_pool > 0 ? tmp_num_pool : 1) * sizeof(int);
        const int with_geometry = (wanted_caps & GX_POOL_GEOMETRY) ? 1 : 0;
        sfc_tile_t      *c_tiles   = with_geometry ? (sfc_tile_t *)      malloc(sz_tiles)   : NULL;
        int             *c_pool    =                 (int *)             malloc(sz_pool);
        tile_bvh_node_t *c_bvh     = with_geometry ? (tile_bvh_node_t *) malloc(sz_bvh)     : NULL;
        float           *c_compact = with_geometry ? (float *)           malloc(sz_compact) : NULL;
        int             *c_types   =                 (int *)             malloc(sz_types);
        /* Reverse map j -> pool_pos (-1 if j is not in this build's pool). Sized
         * to NumPart_when_built; bounds-checked at narrow-refit lookup time. */
        size_t sz_jtop = (size_t)(NumPart > 0 ? NumPart : 1) * sizeof(int);
        int             *c_jtop    = (int *)             malloc(sz_jtop);
        /* An allocation failure here would otherwise be a segfault: the buffers are
         * written unconditionally just below, and this producer is now the only
         * supplier, so there is nothing to fall back to.  The request is RANK-LOCAL
         * on purpose -- this branch is entered per-rank (cache_match keys on
         * rank-local NumPart), so a collective here would deadlock whenever ranks
         * disagree about rebuilding.  The matching drain runs just past the branch,
         * where every rank converges and before anything reads the pool. */
        if(!c_pool || !c_types || !c_jtop
           || (with_geometry && (!c_tiles || !c_bvh || !c_compact))) {
            printf("ERROR: supply-cache allocation failed on task %d (num_pool=%d NumPart=%d geometry=%d)\n",
                   ThisTask, tmp_num_pool, NumPart, with_geometry);
            fflush(stdout);
            free(c_tiles); free(c_pool); free(c_bvh); free(c_compact); free(c_types); free(c_jtop);
            c_tiles = NULL; c_pool = NULL; c_bvh = NULL; c_compact = NULL; c_types = NULL; c_jtop = NULL;
            if(tmp_bvh)   myfree(tmp_bvh);
            if(tmp_tiles) myfree(tmp_tiles);
            if(tmp_pool)  myfree(tmp_pool);
            gizmo_request_controlled_stop(7724, "ghost_exchange: supply-cache allocation failed",
                                          __FILE__, __LINE__, __FUNCTION__);
            cache_alloc_failed = 1;
        }
        if(!cache_alloc_failed) {
        for(int j = 0; j < NumPart; j++) c_jtop[j] = -1;
        if(tmp_ntiles > 0)     memcpy(c_tiles, tmp_tiles, (size_t)tmp_ntiles * sizeof(sfc_tile_t));
        if(tmp_num_pool > 0)   memcpy(c_pool,  tmp_pool,  (size_t)tmp_num_pool * sizeof(int));
        if(tmp_bvh_nnodes > 0) memcpy(c_bvh,   tmp_bvh,   (size_t)tmp_bvh_nnodes * sizeof(tile_bvh_node_t));
        for(int p = 0; p < tmp_num_pool; p++) {
            int j = tmp_pool[p];
            if(with_geometry) {
                c_compact[p*4+0] = (float)P[j].Pos[0];
                c_compact[p*4+1] = (float)P[j].Pos[1];
                c_compact[p*4+2] = (float)P[j].Pos[2];
                /* SSOT: leaf compact h uses the SAME formula as build_sfc_tiles'
                 * per-particle aggregation above (= gx_policy_scaled_h).  Result:
                 * leaf h_j == tile band band's contribution from this particle. */
                c_compact[p*4+3] = (float)gx_policy_scaled_h(j, spec->radius_policy,
                                                            spec->j_radius_scale,
                                                            safety_factor);
            }
            c_types[p] = (int)P[j].Type;
            if(j >= 0 && j < NumPart) c_jtop[j] = p;
        }

        /* Free mymalloc temps in LIFO order. */
        if(tmp_bvh)   myfree(tmp_bvh);
        if(tmp_tiles) myfree(tmp_tiles);
        if(tmp_pool)  myfree(tmp_pool);

        /* Install in cache. */
        g_glt_cache.tiles = c_tiles;
        g_glt_cache.pool  = c_pool;
        g_glt_cache.bvh   = c_bvh;
        g_glt_cache.compact_xyzh = c_compact;
        g_glt_cache.pool_types   = c_types;
        g_glt_cache.j_to_pool    = c_jtop;
        g_glt_cache.ntiles    = tmp_ntiles;
        g_glt_cache.num_pool  = tmp_num_pool;
        g_glt_cache.bvh_nnodes = tmp_bvh_nnodes;
        g_glt_cache.bvh_root  = tmp_bvh_root;
        g_glt_cache.NumPart_when_built = NumPart;
        g_glt_cache.identity_epoch_when_built = g_supply_identity_epoch;
        g_glt_cache.Ti_when_built = All.Ti_Current;
        g_glt_cache.safety_factor_when_built = safety_factor;
        g_glt_cache.radius_policy_when_built = spec->radius_policy;
        g_glt_cache.j_radius_scale_when_built = spec->j_radius_scale;
        /* Fresh build seeded compact_xyzh from current P[]; any pre-existing
         * dirty marks are obsolete.  Next refresh decides full vs narrow from
         * marks accumulated AFTER this point. */
        g_glt_dirty_clear_();
        g_glt_cache.eligible_type_mask_when_built = desired_pool_mask;
        g_glt_cache.needs_refit = 0;
        g_glt_cache.caps = wanted_caps;
        g_glt_cache.valid = 1;
        }   /* end fill+install (skipped when the cache allocation failed) */

        h_tiles = c_tiles; h_pool = c_pool; num_pool = tmp_num_pool;
        ntiles = tmp_ntiles; h_bvh = c_bvh; bvh_nnodes = tmp_bvh_nnodes;
        bvh_root = tmp_bvh_root;
        h_compact_xyzh = c_compact; h_pool_types = c_types;
    }
    /* Both branches converge here, so this poll is reached by every rank: it drains a
     * rank-local supply-cache allocation failure into an all-rank controlled stop
     * BEFORE anything below dereferences the pool. */
    gizmo_exit_bad_stop_if_requested("ghost_exchange:supply_cache_alloc");
    (void)bvh_nnodes;

    /* Periodic flags / box sizes for the BVH walker. */
    int periodic_flags[3] = { TILE_PERIODIC_X, TILE_PERIODIC_Y, TILE_PERIODIC_Z };
    double box_sizes[3]   = { boxSize_X, boxSize_Y, boxSize_Z };

    double t_step3_build = timediff(t_step3_start, my_second());
    if(ThisTask == 0 && gizmo_verbose_diag()) {
        printf("[GX_RD_CACHE rank=0 call=%d caller=%s %s build=%.4f num_pool=%d ntiles=%d pool_mask=0x%x hits=%ld misses=%ld refits=%ld narrow=%ld ghost_epoch=%llu pool_epoch=%llu ghost_start=%d ghost_count=%d]\n",
               this_call, (spec->caller_name ? spec->caller_name : "?"),
               /* A geometry-only rebuild reused the pool but paid the full tile+BVH+
                * compact build, so it is neither a hit nor a miss and is labelled
                * for what it is; `build` on this line is that build's cost. */
               (geometry_rebuilt ? "GEOM_REBUILD" : (from_cache ? "HIT" : "MISS")),
               t_step3_build, num_pool, ntiles,
               g_glt_cache.eligible_type_mask_when_built,
               g_glt_cache_hits, g_glt_cache_misses, g_glt_cache_refits, g_glt_cache_narrow_refits,
               (unsigned long long)gpu_sidx_ghost_epoch(),
               (unsigned long long)gpu_sidx_pool_epoch(),
               gpu_sidx_last_ghost_start(), gpu_sidx_last_ghost_count());
        fflush(stdout);
    }
    double t_step3_walk_start = my_second();

    /* Matched producer selection: for ONEWAY callers the routed top-leaf
     * discovery installs when top-leaf geometry is collectively available;
     * broadcast is the fail-closed path for any spec that is not routing-eligible.
     * Broadcast query gather is LAZY (ensure_broadcast_queries) — skipped when
     * routed installs. */
    /* Routed set from the walk-export producer, built further below so it reads the same
     * supply-cache snapshot as the rest of this call. */
    char *matched_walk_export = NULL;
    struct gx_walk_export_result walk_export_res;
    memset(&walk_export_res, 0, sizeof(walk_export_res));

    /* The tile/BVH walk is no longer ONEWAY's discovery path: the walk-export producer
     * supplies the same set from a bounded walk, so running both would be duplicate work.
     * It remains reachable only as the comparison reference when the route oracle is on. */
    int want_routed = (search_mode == NGB_SEARCH_ONEWAY)
                      && (!gx_walk_export_eligible(spec) || ghost_route_oracle_enabled());
    int route_pre_ok_local = want_routed ? ((topleaf_router_geometry_acquire() == 0) ? 1 : 0) : 0;
    int route_pre_available = 0;
    if(want_routed)
        MPI_Allreduce(&route_pre_ok_local, &route_pre_available, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    char *matched = NULL;
    int   used_routed = 0;
    enum gx_installed_producer installed_producer = GX_INSTALLED_BROADCAST;  /* telemetry only */
    int   use_hier = 1;   /* H2: hierarchical route constructor (production) */
    double t_route_construct = 0.0, t_route_alltoallv = 0.0, t_route_walk = 0.0;


    /* walk_export_only (defined with the cache capabilities above) means: no consumer of the
     * tile/BVH geometry runs on this call, so the broadcast walk is skipped and the
     * geometry was never built.  On producer failure there is no geometry-based
     * fallback left, which is why that case is a controlled stop rather than a
     * silent switch to a walk whose inputs are absent.  Rank-uniform (spec constants
     * + diag + oracle env), so ranks never split across the collectives below. */
    /* Up-front broadcast gather when broadcast is needed regardless of routed outcome
     * (oracle compare, or routed unavailable). Collective + uniform.  Skipped when the
     * walk-export producer is the sole supplier (walk_export_only). */
    if((oracle_on || !route_pre_available) && !walk_export_only) {
        double tb = my_second();
        ensure_broadcast_queries(&bcast_queries_available, local_queries, n_local_queries,
                                 &all_q_counts, &q_disps, &all_queries, &total_queries);
        t_step2 += timediff(tb, my_second());
    }

    if(route_pre_available) {
        matched = compute_matched_routed(local_queries, n_local_queries,
                                         h_compact_xyzh, h_tiles, ntiles,
                                         h_pool, num_pool, h_pool_types, supply_mask,
                                         h_bvh, bvh_root, search_mode, periodic_flags, box_sizes,
                                         use_hier, &t_route_construct, &t_route_alltoallv, &t_route_walk,
                                         NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        if(matched) { used_routed = 1; installed_producer = GX_INSTALLED_ONEWAY_ROUTED_BVH; }   /* NULL => routed failed collectively => fall back */
    }


    if(!matched && !walk_export_only) {
        /* LATE collective fallback: routed unavailable or failed after the skip.  (walk_export_only SKIPS this
         * — that path installs at the producer below and its failure is caught by the controlled
         * stop before Step 4, so matched is never NULL entering the count/pack loops.) */
        double tb = my_second();
        ensure_broadcast_queries(&bcast_queries_available, local_queries, n_local_queries,
                                 &all_q_counts, &q_disps, &all_queries, &total_queries);
        t_step2 += timediff(tb, my_second());
        matched = compute_matched_broadcast(all_queries, q_disps, all_q_counts,
                                            h_compact_xyzh, h_tiles, ntiles,
                                            h_pool, num_pool, h_pool_types, supply_mask,
                                            h_bvh, bvh_root, search_mode, periodic_flags, box_sizes);
        /* Broadcast is the safety path — its alloc failing is terminal.  Drain
         * COLLECTIVELY here, BEFORE Step 4 dereferences matched: a per-rank NULL
         * must become an all-rank controlled stop, never a NULL walk/segfault. */
        int bcast_fail_local = (matched == NULL) ? 1 : 0;
        int bcast_fail_any   = 0;
        MPI_Allreduce(&bcast_fail_local, &bcast_fail_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(bcast_fail_any) {
            if(bcast_fail_local) {
                printf("ERROR: request-driven broadcast matched alloc failed on task %d.\n", ThisTask);
                gizmo_request_controlled_stop(7705, "ghost_exchange (request-driven): broadcast matched alloc failed",
                                              __FILE__, __LINE__, __FUNCTION__);
            }
            gizmo_exit_bad_stop_if_requested("ghost_exchange:broadcast_matched_alloc");
        }
    }

    /* Compare-before-install oracle (transport+oracle): routed must EQUAL broadcast,
     * BOTH ways.  Sender-side compare of the exact (dest_rank, home_index=h_pool[p])
     * send identities Steps 4-6 will install:
     *   UNDER-ROUTE (mb && !mr) — broadcast has a ghost the routed set lacks; the
     *       safety-critical direction (missing physics) -> hard stop 7704.
     *   EXTRA-MATCH (mr && !mb) — routed has a ghost broadcast lacks; a predicate/
     *       snapshot bug (over-import; physics still correct but the invariant is
     *       broken) -> hard stop 7722.
     * Both are drained IMMEDIATELY below, BEFORE any Step-4 count/pack/send touches
     * `matched`.
     * SCOPE, and it is narrow: this compares the tile/BVH routed set against broadcast.
     * For a walk-export-eligible spec the set installed further down is the WALK-EXPORT
     * set, which this comparison never inspects -- so passing here says nothing about
     * what ships.  The once-per-run notice below states that rather than leaving the
     * reader to infer it. */
    if(oracle_on && gx_walk_export_eligible(spec) && ThisTask == 0) {
        static int reported_oracle_scope = 0;
        if(!reported_oracle_scope) {
            reported_oracle_scope = 1;
            printf("GHOST_ROUTE_ORACLE scope=reference_only: it does not inspect the "
                   "walk-export set installed for eligible specs, so a clean result here "
                   "is not a check on the shipped set\n");
            fflush(stdout);
        }
    }
    if(installed_producer == GX_INSTALLED_ONEWAY_ROUTED_BVH && oracle_on && matched) {
        char *matched_bcast = compute_matched_broadcast(all_queries, q_disps, all_q_counts,
                                                        h_compact_xyzh, h_tiles, ntiles,
                                                        h_pool, num_pool, h_pool_types, supply_mask,
                                                        h_bvh, bvh_root, search_mode, periodic_flags, box_sizes);
        long under_local = 0, extra_local = 0;   /* under-route / extra-match counts */
        int  under_t = -1, under_p = -1, extra_t = -1, extra_p = -1;   /* first offenders */
        if(matched_bcast == NULL) {
            under_local = 1;   /* guard alloc fail -> treat as un-provable under-route */
        } else {
            for(int t = 0; t < NTask; t++) {
                if(t == ThisTask) continue;
                const char *mr = matched       + (size_t)t * num_pool;
                const char *mb = matched_bcast + (size_t)t * num_pool;
                for(int p = 0; p < num_pool; p++) {
                    if(mb[p] && !mr[p]) { under_local++; if(under_t < 0) { under_t = t; under_p = p; } }
                    else if(mr[p] && !mb[p]) { extra_local++; if(extra_t < 0) { extra_t = t; extra_p = p; } }
                }
            }
        }
        if(under_local > 0)
            printf("[GX_ROUTE_TRANSPORT call=%d caller=%s rank=%d UNDER-ROUTE: %ld broadcast ghosts missing from routed (first dest=%d pool=%d home_idx=%d)]\n",
                   this_call, (spec->caller_name ? spec->caller_name : "?"), ThisTask, under_local,
                   under_t, under_p, (under_t >= 0 && under_p >= 0) ? h_pool[under_p] : -1);
        if(extra_local > 0)
            printf("[GX_ROUTE_TRANSPORT call=%d caller=%s rank=%d EXTRA-MATCH: %ld routed ghosts not in broadcast (first dest=%d pool=%d home_idx=%d)]\n",
                   this_call, (spec->caller_name ? spec->caller_name : "?"), ThisTask, extra_local,
                   extra_t, extra_p, (extra_t >= 0 && extra_p >= 0) ? h_pool[extra_p] : -1);
        if(under_local > 0 || extra_local > 0) fflush(stdout);
        long uv_in[2] = { under_local, extra_local }, uv_any[2] = { 0, 0 };
        MPI_Allreduce(uv_in, uv_any, 2, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
        if(uv_any[0] > 0) {
            gizmo_request_controlled_stop(7704, "ghost route transport: broadcast ghosts missing from routed set (UNDER-ROUTE)",
                                          __FILE__, __LINE__, __FUNCTION__);
        } else if(uv_any[1] > 0) {
            gizmo_request_controlled_stop(7722, "ghost route transport: routed set has ghosts broadcast lacks (EXTRA-MATCH; predicate/snapshot bug)",
                                          __FILE__, __LINE__, __FUNCTION__);
        } else if(ThisTask == 0) {
            printf("[GX_ROUTE_TRANSPORT call=%d caller=%s OK routed==broadcast]\n",
                   this_call, (spec->caller_name ? spec->caller_name : "?"));
            fflush(stdout);
        }
        free(matched_bcast);
        /* Drain BEFORE Step 4 dereferences `matched` (do not defer to a later poll). */
        gizmo_exit_bad_stop_if_requested("ghost_exchange:route_transport_compare");
    }

    /* H1 flat-vs-hierarchical owner-set oracle (gated; flat stays authoritative).
     * Only meaningful when routed transport is active (geometry acquired + flat
     * route exercised).  Builds hierarchical owner sets and verifies they EQUAL
     * the flat router's; mismatch => collective controlled-stop (new geometry/
     * traversal wrong).  All ranks enter together (gate = env + used_routed, both
     * collective-uniform), so the Allreduce stays symmetric. */
    if(used_routed && ghost_route_hier_oracle_enabled() && search_mode == NGB_SEARCH_ONEWAY) {
        int hq_n = n_local_queries;
        double *hq_pos = (double *) malloc((size_t)(hq_n > 0 ? hq_n : 1) * 3 * sizeof(double));
        double *hq_h   = (double *) malloc((size_t)(hq_n > 0 ? hq_n : 1) * sizeof(double));
        int first_bad = -1;
        int rc = -2;   /* default: scratch alloc failure (caller side) => UNAVAILABLE */
        if(hq_pos && hq_h) {
            for(int i = 0; i < hq_n; i++) {
                hq_pos[i*3+0] = local_queries[i].pos[0];
                hq_pos[i*3+1] = local_queries[i].pos[1];
                hq_pos[i*3+2] = local_queries[i].pos[2];
                hq_h[i]       = local_queries[i].h;
            }
            rc = topleaf_router_hier_vs_flat_check(hq_pos, hq_h, hq_n, 0u, 0,
                                                   periodic_flags, box_sizes, ThisTask, &first_bad);
        }
        free(hq_pos); free(hq_h);
        /* rc > 0: mismatch count.  rc == 0: all equal.  rc < 0: UNAVAILABLE
         * (-1 geometry, -2 alloc, -3 stack overflow) — for a validation gate this
         * is FATAL (not a quiet skip): the hierarchical router could not be proven
         * correct.  Track mismatch and unavailable SEPARATELY so a "could not
         * validate" can never masquerade as OK. */
        int local_mismatch = (rc > 0) ? 1 : 0;
        int local_unavail  = (rc < 0) ? 1 : 0;
        if(local_mismatch)
            printf("[GX_HIER_ORACLE call=%d caller=%s rank=%d MISMATCH: %d/%d queries flat!=hier (first q=%d)]\n",
                   this_call, (spec->caller_name ? spec->caller_name : "?"), ThisTask, rc, hq_n, first_bad);
        if(local_unavail)
            printf("[GX_HIER_ORACLE call=%d caller=%s rank=%d UNAVAILABLE rc=%d (%s)]\n",
                   this_call, (spec->caller_name ? spec->caller_name : "?"), ThisTask, rc,
                   (rc == -1 ? "geometry invalid" : rc == -2 ? "alloc fail" : rc == -3 ? "hier stack overflow" : "unknown"));
        int red_in[2]  = { local_mismatch, local_unavail };
        int red_out[2] = { 0, 0 };
        MPI_Allreduce(red_in, red_out, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        long long covered_local  = (rc >= 0) ? (long long)hq_n : 0;   /* queries actually compared */
        long long covered_global = 0;
        MPI_Allreduce(&covered_local, &covered_global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if(red_out[0]) {
            gizmo_request_controlled_stop(7706, "ghost route hier oracle: hierarchical owner set != flat (geometry/traversal bug)",
                                          __FILE__, __LINE__, __FUNCTION__);
        } else if(red_out[1]) {
            gizmo_request_controlled_stop(7707, "ghost route hier oracle: validation UNAVAILABLE (geometry/alloc/stack-overflow) — hierarchical router not provable",
                                          __FILE__, __LINE__, __FUNCTION__);
        } else if(ThisTask == 0) {
            /* Never call zero-coverage "OK" — a validation gate must prove it
             * actually compared owner sets, not that nothing ran. */
            if(covered_global > 0)
                printf("[GX_HIER_ORACLE call=%d caller=%s OK: hierarchical owner sets == flat over %lld queries (global)]\n",
                       this_call, (spec->caller_name ? spec->caller_name : "?"), covered_global);
            else
                printf("[GX_HIER_ORACLE call=%d caller=%s SKIP: 0 queries compared (nothing to validate this call)]\n",
                       this_call, (spec->caller_name ? spec->caller_name : "?"));
            fflush(stdout);
        }
        gizmo_exit_bad_stop_if_requested("ghost_exchange:hier_oracle");
    }


    double t_step3_walk = timediff(t_step3_walk_start, my_second());
    /* An eligible SYMM spec installs the routed set BELOW this print, so used_routed is
     * not yet set for it; report the path that will actually be used or the parity grep
     * reads "bcast" on a routed call. A producer failure after this point prints
     * GX_R1_FALLBACK, which already invalidates the run as a routed timing arm. */
    const char *qdist = (used_routed || gx_walk_export_eligible(spec)) ? "routed" : "bcast";
    if(ThisTask == 0 && gizmo_verbose_diag()) {
        printf("[GX_RD rank=0 step3 build_tiles+bvh+compact=%.4f s discovery=%.4f s ntiles=%d num_pool=%d total_queries=%d qdist=%s]\n",
               t_step3_build, t_step3_walk, ntiles, num_pool,
               (qdist[0] == 'r' ? -1 : total_queries), qdist);
        fflush(stdout);
    }

    /* Walk-export discovery: produce the routed set (sender export + bounded receiver
     * walk, collective-safe) and INSTALL it for an eligible spec via the shared
     * ownership-transfer.  Placed here so it reads the SAME g_glt_cache snapshot as the
     * rest of this call.  Membership comes from the SSOT accept (gx_pair_accept), so the
     * only way this set can differ from a full walk is routing COVERAGE, which is what
     * the per-spec supply-band domination proof establishes. */
    const int walk_export_install = gx_walk_export_eligible(spec);
    if(walk_export_install && NTask > 1) {
        matched_walk_export = compute_matched_walk_export(spec, local_queries, n_local_queries,
                                                  num_pool, supply_mask, search_mode,
                                                  periodic_flags, box_sizes, &walk_export_res);
        /* Install via the shared ownership-transfer, so Steps 4-6 are reached by exactly
         * one path whichever producer supplied the set. */
        if(matched_walk_export && walk_export_res.status == GX_WALK_EXPORT_OK) {
            if(matched) free(matched);
            matched = matched_walk_export; matched_walk_export = NULL;
            installed_producer = GX_INSTALLED_WALK_EXPORT;
            used_routed = 1;
        } else if(ThisTask == 0) {
            /* Report what ACTUALLY happens next, which depends on whether a reference
             * set exists.  In production (walk_export_only) none does, so the stop below
             * fires.  Under diag/oracle the broadcast walk ran, so the call continues on
             * that set -- correct physics, but NOT the routed substrate, so the run is
             * not a valid routed arm either way.
             * The GX_R1_FALLBACK token is kept only because the current run-validity
             * checks grep for it; it does not describe the mechanism.  Rename it to match
             * the surrounding names, or fold it into the general import-failure
             * reporting, once nothing greps for the old spelling. */
            printf("[GX_R1_FALLBACK call=%d caller=%s reason=producer_status_%d -> %s; INVALID as a routed arm]\n",
                   this_call, (spec->caller_name ? spec->caller_name : "?"), walk_export_res.status,
                   walk_export_only ? "controlled stop (no correctness-proven fallback)"
                                    : "continuing on the broadcast reference set");
            fflush(stdout);
        }
        free(matched_walk_export); matched_walk_export = NULL;
    }

    /* An eligible caller skipped the broadcast walk, so if the walk-export producer
     * did not install (UNAVAILABLE/ALLOC_FAIL) there is no set to install and no
     * substrate left that is known to be correct.  The tile/BVH and broadcast walks
     * both read the same cached geometry, which has been measured producing wrong
     * densities on a decomposition where the cached geometry went stale, so falling
     * back to either would trade a visible failure for a silent one.  Stop instead.
     * walk_export_only and the producer status are rank-uniform, so all ranks stop together.
     * The dominant failure mode is envelope allocation under memory pressure; the
     * recovery that fits it is a retry at reduced import padding inside this producer,
     * which does not exist yet — until it does, the honest outcome is this stop. */
    if(walk_export_only && matched == NULL) {
        gizmo_request_controlled_stop(7723,
            "ghost_exchange: walk-export producer unavailable and no correctness-proven fallback exists",
            __FILE__, __LINE__, __FUNCTION__);
        gizmo_exit_bad_stop_if_requested("ghost_exchange:walk_export_unavailable");
    }

    /* === Step 4: per-peer counts + index list === */
    double t_step4_start = my_second();
    int *send_count = (int *) mymalloc("gx_rd_sc", NTask * sizeof(int));
    int *recv_count = (int *) mymalloc("gx_rd_rc", NTask * sizeof(int));
    int *send_disp  = (int *) mymalloc("gx_rd_sd", NTask * sizeof(int));
    int *recv_disp  = (int *) mymalloc("gx_rd_rd", NTask * sizeof(int));
    for(int t = 0; t < NTask; t++) { send_count[t] = 0; recv_count[t] = 0; }
    for(int t = 0; t < NTask; t++) {
        if(t == ThisTask) continue;
        char *match_for_t = matched + (size_t)t * (size_t)num_pool;
        int s = 0;
        for(int p = 0; p < num_pool; p++) if(match_for_t[p]) s++;
        send_count[t] = s;
    }
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

    /* Check space (mirrors legacy guard).  Request-driven is the last-resort
     * Mode-A discovery — there is NO further fallback — so a count/displacement
     * overflow of the int MPI transport range, or ghosts that would not fit
     * P[]/CellP[], fail HONESTLY via the collective controlled-stop poll below. */
    if(!count_range_ok) {
        printf("ERROR: request-driven ghost exchange counts exceed int MPI transport range on task %d.\n", ThisTask);
        gizmo_request_controlled_stop(7703, "ghost_exchange (request-driven): ghost count/displacement exceeds int MPI transport range", __FILE__, __LINE__, __FUNCTION__);
    } else if(!ghost_particle_slots_fit((long long)NumPart + total_recv_ll)) {
        printf("ERROR: request-driven ghost exchange needs %d ghosts on task %d, only %d free.\n",
               total_recv, ThisTask, All.MaxPart - NumPart);
        gizmo_request_controlled_stop(7702, "ghost_exchange (request-driven): ghost append would exceed MaxPart (raise PartAllocFactor)", __FILE__, __LINE__, __FUNCTION__);
    }
    /* Per-rank capacity check above is asymmetric; drain it at this all-rank poll
     * BEFORE Step 5, so no rank appends ghosts past MaxPart (OOB) or desyncs the
     * collective pack/exchange. Every rank reaches this unconditionally. */
    gizmo_exit_bad_stop_if_requested("ghost_exchange:capacity_rd");

    double t_step4 = timediff(t_step4_start, my_second());
    /* === Step 5: pack particle data + cell data + home_idx === */
    double t_step5_start = my_second();
    struct particle_data *send_P = (struct particle_data *) mymalloc("gx_rd_sP",
        (total_send > 0 ? total_send : 1) * sizeof(struct particle_data));
    struct gas_cell_data *send_CellP = (struct gas_cell_data *) mymalloc("gx_rd_sC",
        (total_send > 0 ? total_send : 1) * sizeof(struct gas_cell_data));
    int *send_home_idx = (int *) malloc((total_send > 0 ? total_send : 1) * sizeof(int));
    {
        int *task_offset = (int *) mymalloc("gx_rd_toff", NTask * sizeof(int));
        memcpy(task_offset, send_disp, NTask * sizeof(int));
        for(int t = 0; t < NTask; t++) {
            if(t == ThisTask) continue;
            char *match_for_t = matched + (size_t)t * (size_t)num_pool;
            for(int p = 0; p < num_pool; p++) {
                if(!match_for_t[p]) continue;
                int j = h_pool[p];
                int off = task_offset[t]++;
                gx_pack_send_slot(P, CellP, j, &send_P[off], &send_CellP[off]);
                send_home_idx[off] = j;
            }
        }
        myfree(task_offset);
    }

    double t_step5 = timediff(t_step5_start, my_second());
    /* === Step 6: Alltoallv particles + cells + home_idx === */
    double t_step6_start = my_second();
    gx_forward_particle_exchange(send_P, send_CellP, send_count, send_disp,
                                 &P[NumPart], &CellP[NumPart], recv_count, recv_disp);

    /* Update counts now so home_idx receive can land at &P[NumPart_before_ghost+...] */
    NumGhostParticles = total_recv;
    NumPart += total_recv;

    /* Mark dirty for compact_xyzh refresh (same as legacy). */
    if(NumGhostParticles > 0) {
        gpu_compact_xyzh_mark_h_dirty_range(NumPart_before_ghost, NumPart);
    }
    /* SIDX lifecycle notify: see comment in tile-overlap impl. Unconditional. */
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

    /* Post-install ghost-set comparison, and it runs ONLY when broadcast supplied the
     * installed set (!used_routed) -- its whole premise is that broadcast is the
     * reference, which is false once anything else installed.
     * It therefore does NOT inspect the walk-export set, which is what installs for an
     * eligible spec.  That is not an oversight to fix here: no independent reference is
     * left to compare against once the broadcast walk is skipped, which is why the
     * notice earlier in this call says so out loud.
     * used_routed is uniform across ranks (both producers are collective), so this gate
     * stays collective-safe. */
    if(ghost_route_oracle_enabled() && search_mode == NGB_SEARCH_ONEWAY && !used_routed) {
        ghost_route_oracle_compare(spec, this_call, local_queries, n_local_queries,
                                   h_tiles, ntiles, h_pool, num_pool,
                                   h_pool_types, h_compact_xyzh, h_bvh, bvh_root,
                                   supply_mask, search_mode, periodic_flags, box_sizes,
                                   ghost_home_rank_map, ghost_home_index_map, total_recv);
    }

    double t_step6 = timediff(t_step6_start, my_second());
    double t_ghost_total = timediff(t_ghost_start, my_second());


    /* Per-rank, per-call Step1-Step6 wall breakdown. Pure diagnostic; gated on
     * GIZMO_VERBOSE_DIAG=1 since both ranks emit (so the user can correlate). */
    if(gizmo_verbose_diag()) {
        printf("[GX_RD_TIME rank=%d call=%d caller=%s mode=%s qdist=%s route=%s installed=%s bcast_gather=%s s1_qbuild=%.4f s2_allgather=%.4f s3_build=%.4f s3_walk=%.4f rcon=%.4f ralltoallv=%.4f rwalk=%.4f s4_count=%.4f s5_pack=%.4f s6_alltoallv=%.4f total=%.4f n_local_queries=%d total_queries=%d num_pool=%d ntiles=%d total_send=%d total_recv=%d]\n",
               ThisTask, this_call, (spec->caller_name ? spec->caller_name : "?"),
               (search_mode == NGB_SEARCH_ONEWAY ? "ONEWAY" : "SYMMETRIC"),
               (used_routed ? "routed" : "bcast"),
               (installed_producer == GX_INSTALLED_ONEWAY_ROUTED_BVH ? (use_hier ? "hier" : "flat") : "-"),
               gx_installed_producer_name(installed_producer),
               (bcast_queries_available ? "gathered" : "skipped"),
               t_step1, t_step2, t_step3_build, t_step3_walk,
               t_route_construct, t_route_alltoallv, t_route_walk,
               t_step4, t_step5, t_step6, t_ghost_total,
               n_local_queries, (used_routed ? -1 : total_queries), num_pool, ntiles, total_send, total_recv);
        fflush(stdout);
    }

    if(ThisTask == 0) {
        PRINT_STATUS("Ghost exchange (request-driven, %s, %s, qdist=%s): %d local + %d ghost  queries=%d total_queries=%d num_pool=%d  [%.4f s]",
                     (spec->caller_name ? spec->caller_name : "?"),
                     (search_mode == NGB_SEARCH_ONEWAY ? "ONEWAY" : "SYMMETRIC"),
                     (used_routed ? "routed" : "bcast"),
                     NumPart_before_ghost, NumGhostParticles,
                     n_local_queries, (used_routed ? -1 : total_queries), num_pool, t_ghost_total);
    }

    /* Diagnostic: ghost composition + import-waste ratio (should be ~0% for
     * the request-driven path by construction since per-particle accept ran
     * before pack — provides direct A/B vs the legacy tile-overlap waste). */
    if(gizmo_verbose_diag() && total_recv > 0) {
        int by_type[6] = {0,0,0,0,0,0};
        for(int g = 0; g < total_recv; g++) {
            int gi = NumPart_before_ghost + g;
            int tt = (int)P[gi].Type;
            if(tt >= 0 && tt < 6) by_type[tt]++;
        }
        printf("[GX_GHOSTTYPE rank=%d call=%d caller=%s total_recv=%d  T0(cells)=%d T1=%d T2=%d T3=%d T4=%d T5=%d]\n",
               ThisTask, this_call, (spec->caller_name ? spec->caller_name : "?"),
               total_recv, by_type[0], by_type[1], by_type[2], by_type[3], by_type[4], by_type[5]);
        fflush(stdout);
    }
    gx_print_waste(spec, this_call, total_recv);

    /* Cleanup local. mymalloc requires LIFO free order. Tile/BVH/pool/
     * compact_xyzh/pool_types are now owned by g_glt_cache (malloc-backed)
     * and outlive this frame; do NOT free them here. They're freed at
     * cache invalidation (drift / domain_decomp hooks) via glt_cache_free. */
    myfree(send_CellP);
    myfree(send_P);
    myfree(recv_disp);
    myfree(send_disp);
    myfree(recv_count);
    myfree(send_count);
    free(send_home_idx);
    free(matched);
    free(all_queries); free(q_disps); free(all_q_counts); free(local_queries);
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
    /* supply_band_dominated=1: this spec's reach is P[j].KernelRadius (the legacy
     * all-types policy) for every type, times safety_factor (checked <= 1 at
     * dispatch). The per-type opener band is seeded per particle from the
     * conservative source union capped at All.MaxKernelRadius (ForceSoftening
     * uncapped; force_hmax_per_type_particle_radius) and exchanged cross-rank on
     * the nodes the export walk descends — so under the standing invariant that
     * no P[].KernelRadius of any type exceeds All.MaxKernelRadius (drift clamps
     * in predict.cc; the density-convergence maxsoft clamp; the sink accretion-
     * radius cap in ags_rkern.cc), the band dominates this reach per type.
     * Any future physics path that writes a non-gas KernelRadius must preserve
     * that invariant or routed discovery under-imports. */
    struct ghost_exchange_spec_t sp = {GHOST_TYPE_ALL, GHOST_TYPE_ALL, NGB_SEARCH_SYMMETRIC, safety_factor, "all_types", -1, NULL, NULL,
                                       MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES, 1.0, 1};
    ghost_exchange_impl(&sp);
}
void ghost_exchange_hydro(double safety_factor)
{
    /* supply_band_dominated=1: gas-only supply at the legacy all-types kernel
     * radius. Every KernelRadius is held at or below MaxKernelRadius (density
     * sink setup clamps it; ags_return_maxsoft clamps the drift path), which is
     * exactly the quantity the per-type node band is built from — so the band is
     * a valid upper bound on this spec's reach and routed discovery is complete.
     * safety_factor is checked separately at dispatch (a >1 factor widens the
     * query beyond the band until the opener scales with it). */
    struct ghost_exchange_spec_t sp = {GHOST_TYPE_0, GHOST_TYPE_0, NGB_SEARCH_SYMMETRIC, safety_factor, "hydro_symmetric", -1, NULL, NULL,
                                       MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES, 1.0, 1};
    ghost_exchange_impl(&sp);
}
void ghost_exchange_hydro_oneway(double safety_factor)
{
    /* ONEWAY routes on search mode alone; the flag is unread here. */
    struct ghost_exchange_spec_t sp = {GHOST_TYPE_0, GHOST_TYPE_0, NGB_SEARCH_ONEWAY, safety_factor, "hydro_oneway", -1, NULL, NULL,
                                       MODE_B_RADIUS_LEGACY_KERNEL_ALLTYPES, 1.0, 0};
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
     * Any dirty-list entries for indices in [NumPart_before_ghost, NumPart)
     * would index out of bounds in compact_h_refresh on the next cached call
     * if a smaller-num_total build happens before SIDX invalidation. The old
     * fail-safe escalated to mark_h_dirty_all, which forced a 1.5s full-pool
     * refresh on every active-sink step (3 cleanups -> 3 full refreshes,
     * measured as the dominant tiny-N cost).
     *
     * Surgical fix: filter the dirty list, dropping ghost-slot indices but
     * keeping the valid local-index entries from density iter / lazy drift.
     * Next gpu_ngb_list_build narrow-refreshes those local entries (~ms);
     * the next ghost_exchange will mark new ghost slots dirty at import
     * time (mark_h_dirty_range above) so symmetric h-reads on ghosts stay
     * fresh. */
    if(NumGhostParticles > 0) {
        gpu_compact_xyzh_dirty_drop_above(NumPart_before_ghost);
    }
    /* SIDX lifecycle notify BEFORE NumPart shrinks. Frees any cached ghost
     * segment in the segmented SIDX (later commit). Called whether or not
     * NumGhostParticles>0 — a cleanup from the no-ghost-imported state is
     * a valid signal that bumps epoch with empty range. */
    gpu_sidx_notify_ghost_cleanup();
    PreviousGhostCount = NumGhostParticles; /* save for domain decomposition headroom */
    NumPart = NumPart_before_ghost;
    N_gas = N_gas_before_ghost;
    NumGhostParticles = 0;
    NumPart_before_ghost = -1;
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
    struct particle_data *send_P = (struct particle_data *) malloc((ns > 0 ? ns : 1) * sizeof(struct particle_data));
    struct gas_cell_data *send_CellP = (struct gas_cell_data *) malloc((ns > 0 ? ns : 1) * sizeof(struct gas_cell_data));
    /* Re-pack current owner values via the SAME slot helper import uses, in the
       SAME send order (ghost_send_home_idx) that produced this pool. */
    for(int k = 0; k < ns; k++) {
        gx_pack_send_slot(P, CellP, ghost_send_home_idx[k], &send_P[k], &send_CellP[k]);
    }
    /* Replay ONLY the forward transport, overwriting the EXISTING ghost slots at
       [NumPart_before_ghost, NumPart). Slots land at identical offsets by
       construction, so ghost_home_rank/index maps + any built CSR stay valid. */
    gx_forward_particle_exchange(send_P, send_CellP, ghost_wb_send_count, ghost_wb_send_disp,
                                 &P[NumPart_before_ghost], &CellP[NumPart_before_ghost],
                                 ghost_wb_recv_count, ghost_wb_recv_disp);
    free(send_P); free(send_CellP);
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
int ghost_get_previous_count(void) { return PreviousGhostCount; }
int ghost_get_num_local(void)  { return (NumPart_before_ghost >= 0) ? NumPart_before_ghost : NumPart; }
int *ghost_get_home_rank(void)  { return ghost_home_rank_map; }
int *ghost_get_home_index(void) { return ghost_home_index_map; }
int *ghost_get_wb_recv_count(void) { return ghost_wb_recv_count; }
int *ghost_get_wb_recv_disp(void)  { return ghost_wb_recv_disp; }
int *ghost_get_wb_send_count(void) { return ghost_wb_send_count; }
int *ghost_get_wb_send_disp(void)  { return ghost_wb_send_disp; }
