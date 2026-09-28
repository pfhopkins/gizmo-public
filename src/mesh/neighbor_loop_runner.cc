/* mesh/neighbor_loop_runner.cc — generic NeighborLoopSpec runner.
 *
 * Mode A (GPU NGL pipeline): stages caller-supplied radii and per-call
 * CallScalars host-side pre-arena, runs gpu_particles_arena_acquire +
 * gpu_ngb_list_build (or stages a caller-injected external CSR), then a
 * Kokkos parallel_for calling Spec::load_active to fill ActiveData[] in
 * UVM (same device epoch as the pair walk), then the parametric
 * pair-kernel parallel_for calling Spec::load_neighbor + Spec::pair_kernel.
 * Launches go through gizmo_gpu_kernel_launch (parallel_for + fence +
 * check_last_error).
 *
 * Mode B (request-driven walker, local + cross-rank peer-to-peer) and the
 * host-side invocation with the
 * lazy-drift boundary structurally encoded as collect_candidates_pre_drift
 * -> drift_and_filter_candidates -> evaluate_pairs_post_drift; it uses
 * the SAME drift epoch as Mode B.
 *
 * The Spec contract (hard-required members, hooks, invariants) is
 * documented in mesh/neighbor_loop_runner.h.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <cstdarg>                                   /* va_list, vfprintf for nlr_warn_once_rank0 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <type_traits>
#include <utility>                                   /* std::declval for SFINAE hook detection */
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"          /* per-TU AllDeviceMirror + device-pass `All` redirect; include before allvars.h */
#include "../declarations/allvars.h"
#include "../declarations/lifecycle_counters.h"      /* g_global_drift_counter etc. for Mode B corridor */
#include "../core/proto.h"
#include "../system/gpu_particles_arena.h"
#include "../declarations/gpu_dispatch_templates.h"  /* gizmo_gpu_kernel_launch */

#include "neighbor_loop_runner.h"
#include "gpu_neighbor_list.h"
#include "device_tree_walk.h"          /* the one device traversal; Mode D enters it from the root */
#include "kernel.h"  /* MUST precede sink_env1_loop.h (kernel_main, NEAREST_XYZ) */
#include "ghost_writeback.h"             /* ghost_get_num_local */
#include "ghost_symlist_lifecycle.h"     /* gizmo_request_filtered_ghost_import_fresh, ghost_exchange_cleanup */
#include "mode_b_local_walker.h"         /* mode_b_local_neighbor_walk, brute */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <algorithm>   /* nth_element / max_element (coverage percentiles) */
#include <unordered_map>
#include <cmath>

#include "mode_b_p2p_transport.h"  /* ModeBBoundedExchange (query/reply transport) */

/* Spec instantiations. Each #include declares one Spec type whose explicit
 * template instantiation appears at the bottom of this file. */
#include "../sinks/sink_env1_loop.h"
#include "../sinks/sink_feed_loop.h"
#include "../sinks/sink_swk_loop.h"
#if defined(SINK_PARTICLES) && defined(SINK_GRAVACCRETION) && (SINK_GRAVACCRETION == 0)
#include "../sinks/sink_env2_loop.h"
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
#include "../gravity/ags_density_loop.h"
#include "../gravity/ags_force_loop.h"
#endif

#if defined(CBE_INTEGRATOR_WITHGRADIENTS)
#include "../sidm/cbe_integrator_gradients.h"
#endif

#include "../hydro/density_loop.h"

#ifdef GALSF_FB_MECHANICAL
#include "../galaxy_sf/mechfb_loop.h"
#endif

#ifdef GALSF_FB_THERMAL
#include "../galaxy_sf/thermal_fb_loop.h"
#endif

#ifdef HYDRO_VOLUME_CORRECTIONS
#include "../hydro/cellcorrections_loop.h"
#endif

/* GradientsSpec always built (gradients has no master #ifdef gate — every
 * hydro build runs hydro_gradient_calc). */
#include "../hydro/gradients_loop.h"

/* HydroForceSpec always built (every hydro build runs hydro_force). */
#include "../hydro/hydro_force_loop.h"

#ifdef GALSF_FB_FIRE_RT_LOCALRP
#include "../galaxy_sf/radfb_rp_loop.h"
#endif
#ifdef DM_DISPERSION_LOOP_ACTIVE
#include "../galaxy_sf/dm_dispersion_loop.h"
#endif
#ifdef TURB_DIFF_DYNAMIC
#include "../turb/difffilter_loop.h"
#endif

#ifdef DM_FUZZY
#include "../sidm/dm_fuzzy_loop.h"
#endif

#ifdef DO_FLUID_ALTSPECIES_DRAG_CALCULATION
#include "../solids/grain_physics_loop.h"
#endif

#ifdef RT_SOURCE_INJECTION
#include "../radiation/rt_source_injection_loop.h"
#endif

/* ============================================================================
 * Shared NLR utility helpers (used by env-config and threshold blocks below).
 * File-scope static; TU-local linkage. Defined here so the threshold helpers
 * (which are file-scope `extern` for runner.h API) can call them.
 * ========================================================================== */

/* One-shot rank-0 warning helper.
 *
 * Cached set keyed by string CONTENT (strcmp), not pointer identity, so
 * dynamically-constructed keys dedupe correctly per distinct key rather than
 * per category. Costs a 64x192
 * BSS array (~12 KB) and an O(N) lookup per call; both fine for a warning
 * surface.
 *
 * Cap is generous (64 distinct keys); if hit, subsequent warnings are
 * silently dropped — they are diagnostics, not correctness gates. */
static void nlr_warn_once_rank0(const char *key, const char *fmt, ...)
{
    if(ThisTask != 0) return;
    static char seen[64][192];
    static int  seen_n = 0;
    for(int i = 0; i < seen_n; i++) {
        if(strcmp(seen[i], key) == 0) return;
    }
    if(seen_n < 64) {
        size_t n = strlen(key);
        if(n >= sizeof(seen[0])) n = sizeof(seen[0]) - 1;
        memcpy(seen[seen_n], key, n);
        seen[seen_n][n] = '\0';
        seen_n++;
    }
    fprintf(stderr, "[NLR env] ");
    va_list ap; va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    fflush(stderr);
}

/* ============================================================================
 * Caller-side helpers (declared in runner.h).
 *
 * Out-of-line so changes to global plumbing (NumPart, P, CellP fetch site,
 * ghost-safety-factor source, etc.) are caller-invisible.
 * ========================================================================== */

/* Thin wrapper over the canonical host accessor; retained for callers
 * that already route through this name. With the per-TU AllDeviceMirror
 * scheme's device-pass-gated `#define All AllDeviceMirror`, bare `All.*`
 * in host code reads the host extern unconditionally, so this
 * indirection is no longer load-bearing — kept for stability of
 * existing call sites. */
const struct global_data_all_processes * nlr_host_all_ptr(void)
{
    return gizmo_host_all_ptr();
}

NlrCommonScalars nlr_common_scalars_from_all(void)
{
    const struct global_data_all_processes *h = nlr_host_all_ptr();
    NlrCommonScalars s;
    s.cf_atime                = h->cf_atime;
    s.cf_a2inv                = h->cf_a2inv;
    s.cf_a3inv                = h->cf_a3inv;
    s.cf_hubble_a             = h->cf_hubble_a;
    s.newton_G                = h->G;
    s.hubble                  = h->HubbleParam;
    s.comoving_integration_on = h->ComovingIntegrationOn;
    return s;
}

neighbor_loop_args nlr_default_args(void)
{
    neighbor_loop_args args;
    args.P                   = P;
    args.CellP               = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;
    args.num_total           = NumPart;
    args.active_list         = nullptr;       /* caller fills */
    args.num_active          = 0;             /* caller fills */
    args.active_call_slot    = nullptr;       /* the caller's list is the call's list */
    args.aux                 = nullptr;       /* caller fills */
    args.ghost_safety_factor = gizmo_ghost_safety_factor();
    args.neighbor_type_mask_override = 0;      /* 0 => use Spec::neighbor_type_mask */
    args.external_csr        = nullptr;        /* nullptr => runner builds its own CSR */
    args.dispatch_override   = NlrForceMode::None; /* None => adaptive threshold */
    return args;
}

void nlr_free_active_list(int *active_list)
{
    if(active_list) myfree(active_list);
}

/* ============================================================================
 * Mode-A/B dispatch thresholds
 *
 * Production dispatch policy is the constexpr Spec::modeb_threshold_sum and
 * Spec::modeb_threshold_max in each NeighborLoopSpec — those are code-level
 * dispatch policy constants for the loop, decided alongside the physics.
 *
 * Resolution precedence (first wins):
 *   1. parameterfile NeighborLoopModeBThreshold{Sum,Max} (-1 = unset;
 *      any other value overrides, <= 0 disables Mode B)
 *   2. Spec::modeb_threshold_{sum,max} constexpr (code default)
 *
 * The threshold is settable ONLY from the parameterfile, so the value a run
 * used is recorded in its parameter log.
 * ========================================================================== */

int gizmo_nlr_modeb_threshold_sum_for(const char *loop_name, int spec_default)
{
    (void)loop_name;
    /* Parameterfile override (the supported interface): -1 = unset -> use the
       Spec default; any other value overrides it (<= 0 disables Mode B). */
    { int p = nlr_host_all_ptr()->NeighborLoopModeBThresholdSum; if(p != -1) return p; }
    return spec_default;
}
int gizmo_nlr_modeb_threshold_max_for(const char *loop_name, int spec_default)
{
    (void)loop_name;
    /* Parameterfile override (the supported interface): -1 = unset -> use the
       Spec default; any other value overrides it (<= 0 disables Mode B). */
    { int p = nlr_host_all_ptr()->NeighborLoopModeBThresholdMax; if(p != -1) return p; }
    return spec_default;
}


/* ============================================================================
 * NeighborLoopPlan path predicates — single source of truth keyed on path.
 *
 * New paths (future Mode C, dual-tree large-N, etc.) extend the switch
 * statements below; never add fields to NeighborLoopPlan.
 * ========================================================================== */
bool nlr_path_uses_imported_ghosts(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return true;
        case NeighborLoopPlan::Path::ModeB_Local:  return false;
        case NeighborLoopPlan::Path::ModeB_Remote: return false;
        case NeighborLoopPlan::Path::ModeD_DeviceFused: return false;
    }
    return false;
}

/* Caller-owned ghost pool: a caller supplying external_csr owns the live
 * particle+ghost pool the CSR indexes (see the GHOST-POOL OWNERSHIP contract
 * in neighbor_loop_runner.h). The runner must not import (slot renumbering
 * would silently invalidate the CSR) and must not cleanup (the pool outlives
 * this call). Single ownership signal by design — no separate flag that
 * could be set inconsistently with external_csr. */
static inline bool nlr_caller_owns_ghost_pool(const neighbor_loop_args &args)
{
    return args.external_csr != nullptr;
}

bool nlr_path_uses_gpu_arena(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return true;
        case NeighborLoopPlan::Path::ModeB_Local:  return false;
        case NeighborLoopPlan::Path::ModeB_Remote: return false;
        case NeighborLoopPlan::Path::ModeD_DeviceFused: return true;
    }
    return false;
}

bool nlr_path_permits_global_numpart_mutation(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return true;
        case NeighborLoopPlan::Path::ModeB_Local:  return false;
        case NeighborLoopPlan::Path::ModeB_Remote: return false;
        case NeighborLoopPlan::Path::ModeD_DeviceFused: return false;
    }
    return false;
}

bool nlr_path_uses_lazy_drift(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return false;
        case NeighborLoopPlan::Path::ModeB_Local:  return true;
        case NeighborLoopPlan::Path::ModeB_Remote: return true;
        case NeighborLoopPlan::Path::ModeD_DeviceFused: return true;
    }
    return false;
}

const char *nlr_path_label(NeighborLoopPlan::Path path)
{
    switch(path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl: return "gpu_ngl";
        case NeighborLoopPlan::Path::ModeB_Local:  return "mode_b_local";
        case NeighborLoopPlan::Path::ModeB_Remote: return "mode_b_remote";
        case NeighborLoopPlan::Path::ModeD_DeviceFused: return "mode_d_fused";
    }
    return "unknown";
}



/* ============================================================================
 * SIDX cache kind resolver (private to runner; see neighbor_loop_runner.h
 * SidxCacheKind enum doc).
 * ========================================================================== */

static gpu_spatial_index_t* nlr_resolve_sidx_cache(SidxCacheKind k,
                                                    const char *loop_name)
{
    switch(k) {
        case SidxCacheKind::AllTypes:
            return gpu_step_sidx_alltypes_ptr();
        case SidxCacheKind::GasOnly:
            return gpu_step_sidx_ptr();
        case SidxCacheKind::None:
            return nullptr;
    }
    /* Unreachable today — kept exhaustive for compiler warnings on enum
     * additions. Other kinds land alongside their first caller. */
    fprintf(stderr, "neighbor_loop_runner: SidxCacheKind=%d not implemented "
            "for loop '%s'\n",
            (int)k, loop_name ? loop_name : "?");
    fflush(stderr);
    endrun(81030);
    return nullptr;
}

/* ============================================================================
 * Mode B local self-rank helpers
 *
 * These three helpers STRUCTURALLY ENCODE the lazy-drift invariant from the
 * neighbor-loop binding contract:
 *
 *     collect_candidates_pre_drift<Spec>   — the tree walk records every
 *                                            eligible particle under the nodes
 *                                            it opens, without reading P[j]'s
 *                                            possibly-stale position
 *     drift_and_filter_candidates<Spec>    — drift_particle on every j to
 *                                            All.Ti_Current (a j already
 *                                            current returns at once), then
 *                                            keep exactly the neighbours at
 *                                            current positions
 *     evaluate_pairs_post_drift<Spec>      — calls Spec::pair_kernel via the
 *                                            same KOKKOS_INLINE_FUNCTION
 *                                            Spec::load_active /
 *                                            Spec::load_neighbor used by
 *                                            run_mode_a (host invocation here)
 *
 * Candidate sets are collected before any drift. Annotation
 * on each helper makes the ordering enforcement structural.
 *
 * Self-rank only here: candidates are local real P[] indices in
 * [0, num_local). Cross-rank peer-to-peer comes below.
 *
 * Walker buffer sized to num_local (worst-case SYMMETRIC, no h-bound
 * pre-pruning).
 * ========================================================================== */

/* Mode-B discovery-walk OpenMP threading. The self/receiver walks write into
 * disjoint per-item output slots, so they parallelize over the item index with
 * no shared writes. Threading engages only above a structural work threshold
 * (never a caller name); below it the serial code runs verbatim (tiny-N steps
 * pay nothing). Chunk sizes are separate constants because received-query work
 * variance differs from self-active variance. */
static constexpr int MODEB_OMP_MIN_PER_THREAD = 4;    /* work >= max(64, 4*nthreads) to thread */
static constexpr int MODEB_OMP_CHUNK_ACTIVE   = 16;   /* schedule(dynamic) chunk for self-active loops */
static constexpr int MODEB_OMP_CHUNK_RECV     = 16;   /* schedule(dynamic) chunk for received-query loop */

static inline int nlr_modeb_omp_nthreads(void)
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/* Structural gate: thread iff more than one thread AND enough work to amortize.
 * Reads NOTHING tree- or drift-related; safe to call on any step. */
static inline bool nlr_modeb_use_omp(long long n_items, int nthreads)
{
    if(nthreads <= 1) return false;
    const long long floor_work = (long long)MODEB_OMP_MIN_PER_THREAD * nthreads;
    const long long thresh = (floor_work > 64) ? floor_work : 64;
    return n_items >= thresh;
}

/* WHERE one batch of fused-walk sources runs: device, host threads, or one host
 * core. Every fused walk goes through here -- the source loop and the placement
 * decision exist once, and the five walks that use it (record and evaluate, from
 * the root and resumed from a peer's start nodes, plus the single-rank evaluator)
 * keep their own surrounding protocol and share nothing else. Adding a sixth walk
 * means calling this, not writing another loop.
 *
 * Nothing about the traversal, the leaf policy, the drift, the exchange or the
 * reply protocol varies with the choice -- only the execution space. The body is a
 * `KOKKOS_LAMBDA`, which both device compilers expand to `[=] __host__ __device__`,
 * so one body serves all three arms unchanged.
 *
 * The count is the batch's OWN source count, never a step-level or global one: a
 * rank walking ten sources of its own is a different question from the same rank
 * answering ten thousand imported ones, and on a clustered run those differ by
 * three orders of magnitude within a single call.
 *
 * CAPTURE RESIDENCE -- the one way this can go silently wrong. The body captures by
 * value, so a CallScalars, a device context or a tree view it captures lives INSIDE
 * the closure: on the host stack when called here, device-resident when copied to a
 * launch. A leaf built inside the body therefore points into whichever space is
 * executing, which is what makes all three arms correct. Every other pointer a leaf
 * holds -- particles, accumulators, the touched set, the anomaly word -- is
 * SharedSpace and valid in both. Capturing any of the first group BY REFERENCE
 * would leave a leaf pointing at a host stack object on the device arm, which is a
 * fault or a wrong answer rather than a compile error. Those captures must also
 * stay unmodified for the whole batch, and no pointer into the closure may outlive
 * the call.
 *
 * Two preconditions, both true at present call sites: any device work this batch
 * reads from has already completed (each producer fences at its own launch), and no
 * caller is already inside a parallel region. */
template <class F>
static inline void nlr_walk_for_sources(const char *tag, int n, F &&body)
{
    if(n <= 0) {return;}
    if(n >= GPU_MIN_SOURCES_FOR_WALK_OFFLOAD) {
        gizmo_gpu_kernel_launch(tag, n, std::forward<F>(body));
        return;
    }
    if(nlr_modeb_use_omp(n, nlr_modeb_omp_nthreads())) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, MODEB_OMP_CHUNK_ACTIVE)
#endif
        for(int kk = 0; kk < n; kk++) {body(kk);}
        return;
    }
    for(int kk = 0; kk < n; kk++) {body(kk);}
}

/* Eval-threading policy for evaluate_pairs_post_drift. The production tree eval
 * may thread (BitwiseReadonly or EpsilonAtomic specs, i.e. any non-SerialOnly tier).
 * Passed explicitly at every call site (no default) so production vs
 * reference eval is greppable and can never be silently mis-gated. */
enum class EvalOMPPolicy { AllowProduction, ForceSerialReference };


/* WALK-ONLY. Does NOT mutate P[].Pos/Vel — drift_particle must not be
 * called from inside this helper. (Audited 2026-05-08: walker calls only
 * force_drift_node on tree-internal nodes, which is search-side state.)
 *
 * Per-active variant: walks args.P[i].Pos at radii[aa] for each
 * aa in [0,args.num_active). Self-rank queries.
 */
template <typename Spec>
static void collect_candidates_pre_drift(const neighbor_loop_args& args,
                                          const double *radii,
                                          unsigned int neighbor_type_mask,
                                          DispatchPath backend,
                                          std::vector<std::vector<int>>& per_active_cands)
{
    /* neighbor_type_mask is an explicit caller parameter (mask-threading
     * refactor). Non-iter callers pass Spec::neighbor_type_mask (unchanged
     * behavior); iter dispatch passes sg.j_type_bitmask for per-subgroup walks. */
    const int N = args.num_active;
    const int num_local = ghost_get_num_local();
    /* Per-call ownership: each inner vector owns its capacity for the
     * duration of this call. .assign(N, {}) clobbers any stale residue
     * from a previous call. Walker-append contract: appends via
     * push_back; geometric growth handles any-size match set without
     * imposing the previous full-pool .assign(num_local, 0) cost
     * (~24 MB × N_active on fire_m11i). */
    per_active_cands.assign(N, std::vector<int>{});
    if(num_local <= 0) return;
    if(backend != DispatchPath::ModeB_HostWalker) {
        fprintf(stderr, "neighbor_loop_runner: collect_candidates_pre_drift "
                "called with non-Mode-B backend (%d) for loop '%s'\n",
                (int)backend, Spec::loop_name);
        fflush(stderr);
        endrun(81033);
        return;
    }
    const double jscale = nlr_spec_symmetric_j_radius_scale<Spec>();
    /* Thread the tree walk above the work threshold; below it stays serial. */
    const int nthreads = nlr_modeb_omp_nthreads();
    const bool use_omp = (backend == DispatchPath::ModeB_HostWalker) &&
                         nlr_modeb_use_omp(N, nthreads);
    if(use_omp) {
        /* Each thread mutates only its own per_active_cands[aa] (outer vector
         * pre-sized; no shared push) and its own drift counter (diagnostic
         * only — allocated iff a counter sink was requested). */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, MODEB_OMP_CHUNK_ACTIVE)
#endif
        for(int aa = 0; aa < N; aa++) {
            const int i = args.active_list[aa];
            const double h_q = radii[aa];
            if(h_q <= 0) continue;
            std::vector<int>& cands = per_active_cands[aa];
            cands.clear();
            if(cands.capacity() == 0) cands.reserve(64);
            double pos_arr[3] = {(double)args.P[i].Pos[0],
                                  (double)args.P[i].Pos[1],
                                  (double)args.P[i].Pos[2]};
            mode_b_local_neighbor_walk(pos_arr, h_q, neighbor_type_mask,
                                       Spec::search_mode, Spec::radius_policy,
                                       cands, jscale);
        }
        return;
    }
    for(int aa = 0; aa < N; aa++) {
        const int i = args.active_list[aa];
        const double h_q = radii[aa];
        if(h_q <= 0) continue;
        std::vector<int>& cands = per_active_cands[aa];
        cands.clear();
        if(cands.capacity() == 0) cands.reserve(64); /* small initial; grows geometrically */
        double pos_arr[3] = {(double)args.P[i].Pos[0],
                              (double)args.P[i].Pos[1],
                              (double)args.P[i].Pos[2]};
        mode_b_local_neighbor_walk(pos_arr, h_q, neighbor_type_mask,
                                    Spec::search_mode, Spec::radius_policy,
                                    cands, jscale);
    }
}

/* WALK-ONLY. Peer-side variant: walks against the LOCAL pool using each
 * remote query's pos/h_search drawn from peer_actives[k].{pos,h_search}.
 * Used for queries received from other ranks via the peer-to-peer transport.
 *
 * Pulls pos/h_search directly from Spec::ActiveData fields. The current
 * convention is that every Spec exposes `pos` and `h_search` as flat
 * ActiveData fields (matches sink_env1_loop.h template). Generalizing
 * this access pattern (e.g. to a Spec::query_pos/query_h trait pair) is
 * tracked for the runner-template-hardening pass.
 */
template <typename Spec>
static void collect_candidates_for_remote_queries(
    const std::vector<typename Spec::ActiveData>& peer_actives,
    const std::vector<int>& peer_nodelist_flat,   /* K*NODELISTLENGTH; exported start-nodes per query */
    const std::vector<int>& peer_nnodes,          /* K; valid entries per query's NodeList */
    unsigned int neighbor_type_mask,
    DispatchPath backend,
    std::vector<std::vector<int>>& per_query_cands)
{
    /* neighbor_type_mask is an explicit caller parameter (mask-threading
     * refactor). Non-iter callers pass Spec::neighbor_type_mask; iter
     * dispatch passes sg.j_type_bitmask.
     *
     * ModeB_HostWalker resumes the walk from the exported NodeList
     * start-nodes (legacy mode==1). */
    const int K = (int)peer_actives.size();
    const int num_local = ghost_get_num_local();
    per_query_cands.assign(K, std::vector<int>{});
    if(num_local <= 0) return;
    if(backend != DispatchPath::ModeB_HostWalker) {
        fprintf(stderr, "neighbor_loop_runner: collect_candidates_for_remote_queries"
                " bad backend %d for loop '%s'\n", (int)backend, Spec::loop_name);
        fflush(stderr);
        endrun(81033);
        return;
    }
    const double jscale = nlr_spec_symmetric_j_radius_scale<Spec>();
    /* Thread the received-query tree walk above the work threshold; brute
     * below-threshold stays serial (byte-identical). */
    const int nthreads = nlr_modeb_omp_nthreads();
    const bool use_omp = (backend == DispatchPath::ModeB_HostWalker) &&
                         nlr_modeb_use_omp(K, nthreads);
    if(use_omp) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, MODEB_OMP_CHUNK_RECV)
#endif
        for(int k = 0; k < K; k++) {
            const auto& active = peer_actives[k];
            const double h_q = (double)active.h_search;
            if(h_q <= 0) continue;
            std::vector<int>& cands = per_query_cands[k];
            cands.clear();
            if(cands.capacity() == 0) cands.reserve(64);
            double pos_arr[3] = {(double)active.pos[0], (double)active.pos[1], (double)active.pos[2]};
            if(peer_nnodes[k] > 0) {
                mode_b_walk_from_start_nodes(pos_arr, h_q, neighbor_type_mask,
                                             Spec::search_mode, Spec::radius_policy,
                                             &peer_nodelist_flat[(size_t)k * NODELISTLENGTH],
                                             peer_nnodes[k], cands, jscale);
            } else {
                mode_b_local_neighbor_walk(pos_arr, h_q, neighbor_type_mask,
                                            Spec::search_mode, Spec::radius_policy,
                                            cands, jscale);
            }
        }
        return;
    }
    for(int k = 0; k < K; k++) {
        const auto& active = peer_actives[k];
        const double h_q = (double)active.h_search;
        if(h_q <= 0) continue;
        std::vector<int>& cands = per_query_cands[k];
        cands.clear();
        if(cands.capacity() == 0) cands.reserve(64);
        double pos_arr[3] = {(double)active.pos[0], (double)active.pos[1], (double)active.pos[2]};
        if(peer_nnodes[k] > 0) {
            /* Targeted query: bounded resume from the exported start-nodes. */
            mode_b_walk_from_start_nodes(pos_arr, h_q, neighbor_type_mask,
                                         Spec::search_mode, Spec::radius_policy,
                                         &peer_nodelist_flat[(size_t)k * NODELISTLENGTH],
                                         peer_nnodes[k], cands, jscale);
        } else {
            /* Broadcast query (n_nodes==0): the sender took the broadcast
             * path (uncovered radius policy) — full local
             * walk from root (the prior broadcast behavior). */
            mode_b_local_neighbor_walk(pos_arr, h_q, neighbor_type_mask,
                                        Spec::search_mode, Spec::radius_policy,
                                        cands, jscale);
        }
    }
}

/* Drift each query's walk candidates to the current time and keep exactly its neighbours at their
 * current positions (mode_b_drift_and_filter_candidates).  query_at(aa, pos, h_q) must give the query
 * that list was walked with; every query is read before any candidate is drifted, since a query's own
 * particle may be another query's candidate.  Serial, as drift_particle requires; drift_particle
 * returns at once for a particle already current, so one reached by several queries is drifted once. */
template <typename Spec, typename QueryAt>
static void drift_and_filter_candidates(std::vector<std::vector<int>>& per_active_cands,
                                        unsigned int neighbor_type_mask,
                                        QueryAt query_at)
{
    const double jscale = nlr_spec_symmetric_j_radius_scale<Spec>();
    const size_t n_queries = per_active_cands.size();
    std::vector<double> queries(4 * n_queries);   /* pos[3], h_q per query */
    for(size_t aa = 0; aa < n_queries; aa++) {
        if(per_active_cands[aa].empty()) continue;
        query_at((int)aa, &queries[4 * aa], queries[4 * aa + 3]);
    }
    for(size_t aa = 0; aa < n_queries; aa++) {
        std::vector<int>& v = per_active_cands[aa];
        if(v.empty()) continue;
        mode_b_drift_and_filter_candidates(&queries[4 * aa], queries[4 * aa + 3], neighbor_type_mask,
                                           Spec::search_mode, Spec::radius_policy, jscale, v);
    }
}

/* Evaluates pair_kernel POST-DRIFT for a precomputed ActiveData[] paired with
 * candidate lists. Caller is responsible for having frozen actives[] before
 * any drift, and for having drifted the candidate union before calling.
 *
 * Decoupled from args.active_list: works for self-rank (actives built from
 * args.active_list) and peer-rank (actives built from received envelopes).
 * Same KOKKOS_INLINE_FUNCTION Spec::load_neighbor + Spec::pair_kernel either
 * way.
 */
template <typename Spec, typename DeviceCtx>
static void evaluate_pairs_post_drift(const DeviceCtx& ctx,
                                       const typename Spec::ActiveData *actives,
                                       int N,
                                       const std::vector<std::vector<int>>& per_active_cands,
                                       typename Spec::AccumData *accums,
                                       const typename Spec::CallScalars& cs,
                                       EvalOMPPolicy eval_policy)
{
    using NeighborData = typename Spec::NeighborData;
    using ScatterData  = typename Spec::ScatterData;
    const struct GxMotionTargetSet motion_targets = gx_motion_target_view();

    /* Per-active evaluation. Writes ONLY accums[aa] plus call-local scratch, so
     * distinct aa are independent — the invariant the BitwiseReadonly threading
     * below relies on. */
    auto eval_one = [&](int aa) {
        Spec::zero_accum(accums[aa]);
        ScatterData s{};                              /* NoScatter for ActiveReduceOnly */
        const auto& cands = per_active_cands[aa];
        if constexpr (nlr_spec_has_bind_active_to_eval_context_v<Spec>) {
            /* Spec needs per-eval-pass rebinding of rank-local + eval-pass-local
             * fields embedded in ActiveData (P_base / CellP_base / per-rank
             * gas-delta ptrs / per-rank index bounds). Take a local
             * mutable copy, refresh from the eval ctx, then walk. See
             * nlr_spec_has_bind_active_to_eval_context_v in
             * mesh/neighbor_loop_runner.h for the contract. */
            typename Spec::ActiveData a = actives[aa];
            Spec::bind_active_to_eval_context(ctx, a);
            for(size_t kk = 0; kk < cands.size(); kk++) {
                int j = cands[kk];
                if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {gx_motion_target_mark(motion_targets, j);}
                IdentitySidecar id{};                 /* NoIdentity */
                NeighborData nb = Spec::load_neighbor(ctx, j, id, a);
                Spec::pair_kernel(a, nb, accums[aa], s, cs);
            }
        } else {
            /* No bind hook: walk by const-ref, no copy. */
            const auto& a = actives[aa];
            for(size_t kk = 0; kk < cands.size(); kk++) {
                int j = cands[kk];
                if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {gx_motion_target_mark(motion_targets, j);}
                IdentitySidecar id{};                 /* NoIdentity */
                NeighborData nb = Spec::load_neighbor(ctx, j, id, a);
                Spec::pair_kernel(a, nb, accums[aa], s, cs);
            }
        }
    };

    /* Two tiers thread the production eval over aa; each aa writes only its own
     * accums[aa] plus call-local scratch (fixed per-active neighbor order):
     *   BitwiseReadonly — no j-side writes at all -> bit-identical regardless of
     *     thread count/schedule.
     *   EpsilonAtomic   — j-side scatter goes through Kokkos::atomic_* (the same
     *     atomics Mode A's kernel already applies over concurrent actives), so
     *     threading only re-orders those atomic updates -> ulp-class
     *     nondeterminism, NOT a new race. The per-active AccumData stays
     *     order-independent (no j-write is read back into the kernel).
     * Only the production path threads; the reference eval
     * (ForceSerialReference) and SerialOnly specs run the serial loop verbatim.
     * Below the work threshold no OpenMP region is entered — the serial path is
     * byte-identical to the unthreaded code. */
    if constexpr (nlr_spec_modeb_eval_omp<Spec>() == ModeBEvalOMP::BitwiseReadonly) {
        static_assert(!Spec::uses_ghost_writeback,
                      "BitwiseReadonly eval tier requires uses_ghost_writeback==false "
                      "(no j-side ghost scatter)");
    }
    const int  eval_nthreads = nlr_modeb_omp_nthreads();
    const bool eval_use_omp  =
        (eval_policy == EvalOMPPolicy::AllowProduction) &&
        (nlr_spec_modeb_eval_omp<Spec>() != ModeBEvalOMP::SerialOnly) &&
        nlr_modeb_use_omp((long long)N, eval_nthreads);

    if(eval_use_omp) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, MODEB_OMP_CHUNK_ACTIVE)
#endif
        for(int aa = 0; aa < N; aa++) eval_one(aa);
    } else {
        for(int aa = 0; aa < N; aa++) eval_one(aa);
    }
}

/* ============================================================================
 * run_mode_b_local<Spec> — single-rank Mode B path.
 *
 * Three-epoch staging contract (matches Mode A timing exactly except step 3
 * runs host-side per query instead of inside a Kokkos kernel):
 *   (1) host pre-drift: Spec::search_radius  → radii[num_active]
 *   (2) host pre-drift: Spec::populate_call_scalars → CallScalars cs
 *   (3) host post-drift: Spec::load_active per query (KOKKOS_INLINE_FUNCTION
 *                         invoked host-side; same UVM/Pos epoch as Mode A's
 *                         device kernel)
 *
 * No Kokkos kernels in this path → no fence (project rule from contract).
 * No GPU NGL build, no SIDX, no arena_acquire — Mode B's whole point is
 * escaping that overhead for tiny-N.
 * ========================================================================== */

/* Helper: build the host-frozen ActiveData[] for self-rank actives. Called
 * BEFORE any drift so the active snapshot matches the legacy pack_query
 * epoch. (Note: this is NOT the same epoch as Mode A's device-staged
 * load_active post-NGL-build; documented at run_neighbor_loop dispatch site.)
 */
template <typename Spec, typename DeviceCtx>
static void build_self_actives_host_pre_drift(
    const neighbor_loop_args& args,
    const DeviceCtx& ctx,
    const double *radii,
    const typename Spec::CallScalars& cs,
    typename Spec::ActiveData *actives_out)
{
    const int N = args.num_active;
    for(int aa = 0; aa < N; aa++) {
        const int slot = args.active_call_slot ? args.active_call_slot[aa] : aa;
        actives_out[aa] = Spec::load_active(ctx, slot, args.active_list[aa],
                                             radii[aa], cs);
    }
}

template <typename Spec>
static void run_mode_b_local(const neighbor_loop_args& args, const double *radii)
{
    using ActiveData = typename Spec::ActiveData;
    using AccumData  = typename Spec::AccumData;
    using DeviceCtx  = typename Spec::DeviceContext;

    const int N = args.num_active;
    if(N <= 0) {
        /* Local-zero is legitimate (caller may enter unconditionally if
         * global_num_active > 0 — multi-rank Mode B requires that). */
        return;
    }

    /* (1) Radii are runner-staged and passed in; pointer is call-lifetime
     * only (do not store). See run_neighbor_loop contract in the header. */

    /* (2) Host pre-drift: per-call scalar globals into POD. */
    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    DeviceCtx ctx;
    ctx.P         = args.P;
    ctx.CellP     = args.CellP;
    ctx.num_total = args.num_total;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; the runner captures it by value into Kokkos device lambdas");
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);
    /* A hook that could not stage its buffers has already asked for the stop.
     * This path is host-local with no collectives, so returning leaves no peer
     * waiting; the guard above releases whatever the hook did obtain. */
    if(ctx.populate_failed) { return; }

    /* Freeze active snapshots host-side BEFORE drift. */
    std::vector<ActiveData> actives(N);
    build_self_actives_host_pre_drift<Spec>(args, ctx, radii, cs, actives.data());

    /* Helper layout: collect → drift → evaluate. */
    std::vector<std::vector<int>> cand_modeB;
    const unsigned int modeb_type_mask = nlr_effective_neighbor_type_mask(args, Spec::neighbor_type_mask);
    {
        collect_candidates_pre_drift<Spec>(args, radii, modeb_type_mask,
                                            DispatchPath::ModeB_HostWalker, cand_modeB);
    }
    {
        drift_and_filter_candidates<Spec>(cand_modeB, modeb_type_mask, [&](int aa, double *pos, double &h_q) {
            const int i = args.active_list[aa];
            pos[0] = (double)args.P[i].Pos[0]; pos[1] = (double)args.P[i].Pos[1]; pos[2] = (double)args.P[i].Pos[2];
            h_q = radii[aa];
        });
    }

    std::vector<AccumData> accums(N);
    {
        evaluate_pairs_post_drift<Spec>(ctx, actives.data(), N, cand_modeB, accums.data(), cs, EvalOMPPolicy::AllowProduction);
    }

    /* Host writeback — same code path as Mode A's writeback. */
    {
        for(int aa = 0; aa < N; aa++) {
            Spec::apply_active_writeback(args, aa, args.active_list[aa], accums[aa]);
        }
    }
}






/* ============================================================================
 * Mode B remote (multi-rank) helpers
 *
 * STRONG INVARIANT:
 *   collect-all → drift-union → evaluate-all on each rank.
 *
 * On each rank, the local pool is walked TWICE per call: once for this
 * rank's own self-pair queries (its own actives), and once for queries
 * received from peers. Both candidate sets must be collected pre-drift,
 * then the drift covers their UNION (idempotent so duplicates are free),
 * then evaluation runs post-drift.
 *
 * Sequence:
 *   stage 1 (active rank) build self radii, cs, frozen actives[]
 *   stage 2 (active rank) build envelopes, ALL peers in broadcast pattern
 *   stage 3 (this rank)   collect self candidates pre-drift
 *   stage 4 (collective)  exchange queries (peer-to-peer) -> recv envelopes
 *   stage 5 (this rank)   flatten envelopes to peer_actives[] + provenance[]
 *   stage 6 (this rank)   collect peer candidates pre-drift
 *   stage 7 (this rank)   drift the UNION of the self and peer candidate sets
 *   stage 8 (this rank)   evaluate self candidates
 *   stage 9 (this rank)   evaluate peer candidates
 *   stage 10 (collective) exchange replies (tree result is what ships)
 *   stage 11 (active rank) merge replies via Spec::merge_accum, ascending
 *                                  rank for FP-reproducible order
 *   stage 12 (active rank) writeback
 *
 * Note: the host-frozen actives[] match Mode B /
 * legacy pack_query epoch, NOT Mode A's device-staged post-NGL-build
 * epoch — the two active epochs are not bit-equivalent.
 * ========================================================================== */

/* ============================================================================
 * mode_b_remote_evaluate_into_buffer<Spec> — extracted helper.
 *
 * Mechanical refactor of the existing run_mode_b_remote_impl body, extracting
 * stages 1-12 (queries / collect / drift / evaluate / merge / replies / merge)
 * + the env-gated active_dumps diagnostic into a reusable helper. The caller
 * provides the output AccumData buffer; the helper does NOT call
 * Spec::apply_active_writeback. Final-only writeback decision stays with the
 * caller ("let the caller decide whether to call
 * apply_active_writeback").
 *
 * Two callers:
 *   - run_mode_b_remote_impl (non-iterative wrapper): allocates per-call
 *     AccumData buffer, calls helper (which emits dumps internally on its
 *     locally-built actives[] view), then runs apply_active_writeback.
 *   - nlr_iter_dispatch_subgroup_mode_b_remote (iterative dispatch): calls
 *     helper with driver-owned compacted AccumData buffer; runner does
 *     final-only apply_active_writeback after the iter loop. (Per-iter
 *     active_dumps emits are harmless: env-gated, only fire under spike-test
 *     flags.)
 *
 * Epoch order is preserved EXACTLY: build self/peer
 * queries -> collect all candidate sets pre-drift -> exchange -> drift union
 * -> evaluate -> exchange replies -> merge replies
 * in deterministic peer order. No re-derivation; line-for-line move from
 * the old impl.
 * ========================================================================== */
/* ============================================================================
 * Answering another rank's queries against this rank's particles.
 *
 * A rank that receives a batch of queries has to return one accumulator per
 * query, and there is more than one way to arrive at them.  The host way walks
 * the local tree per query into a candidate list, drifts whatever the walk
 * touched, and then runs the pair kernel over the list.  A device way walks and
 * evaluates in a single pass and never builds a list at all.
 *
 * The seam is drawn around ALL THREE steps rather than around the evaluation
 * alone, and that is the whole point: an evaluation-only hook would take the
 * candidate list as an argument, so any backend behind it would still have to
 * build one -- which is precisely the cost a fused walk exists to avoid.  Drawn
 * here, "how the neighbours are found" and "whether they are ever written down"
 * are the backend's business rather than the transport's.
 *
 * The transport around this -- envelope exchange, provenance, reply exchange,
 * the bounded round loop -- does not know which backend ran and must not learn.
 * There is one body of MPI choreography for the peer wire and it stays that way.
 * ========================================================================== */
template <typename Spec>
struct NlrPeerAnswerHostWalk {
    using AccumData = typename Spec::AccumData;

    /* The host path, unchanged: collect against the local pool, drift what the
     * walk touched, evaluate. */
    static void answer(const typename Spec::DeviceContext& ctx,
                       const typename Spec::CallScalars& cs,
                       const std::vector<typename Spec::ActiveData>& peer_actives,
                       const std::vector<int>& peer_nodelist_flat,
                       const std::vector<int>& peer_nnodes,
                       unsigned int neighbor_type_mask,
                       std::vector<AccumData>& peer_replies_out)
    {
        /* Stage 6: collect PEER candidate sets PRE-DRIFT (against MY local pool). */
        std::vector<std::vector<int>> cand_peer_tree;
        {
            collect_candidates_for_remote_queries<Spec>(peer_actives,
                                                         peer_nodelist_flat, peer_nnodes,
                                                         neighbor_type_mask,
                                                         DispatchPath::ModeB_HostWalker,
                                                         cand_peer_tree);
            /* Reported rather than recorded here: the threaded-walk note is the
             * transport's own diagnostic, and it lives in a lambda that closes
             * over the transport's locals.  A backend that reached for it would
             * be reaching out of its scope, which is what the first draft of
             * this extraction did. */
        }

        /* Stage 7 (peer): drift THIS round's peer candidate sets (self candidates
         * were drifted once before the round loop). Idempotent to All.Ti_Current. */
        {
            drift_and_filter_candidates<Spec>(cand_peer_tree, neighbor_type_mask, [&](int k, double *pos, double &h_q) {
                pos[0] = (double)peer_actives[k].pos[0]; pos[1] = (double)peer_actives[k].pos[1]; pos[2] = (double)peer_actives[k].pos[2];
                h_q = (double)peer_actives[k].h_search;
            });
        }

        /* Stage 9: evaluate PEER queries post-drift -> peer_replies, shipped back
         * to the home rank. */
        const int K = (int)peer_actives.size();
        if(K > 0) {
            evaluate_pairs_post_drift<Spec>(ctx, peer_actives.data(), K,
                                              cand_peer_tree, peer_replies_out.data(), cs, EvalOMPPolicy::AllowProduction);
        }
    }
};

/* The local queries a fused round co-schedules with the ones it received.
 *
 * Defined here rather than beside the walk that consumes it because the
 * transport above fills it in, and it needs nothing the transport does not
 * already have -- no allocator, no leaf policy, just the two Spec types.
 *
 * Empty by default, which is the single-rank case and every round after the one
 * that took this call's locals.  `accums_out` is WRITTEN, not merged: a local
 * query is evaluated exactly once per call, and the reply merge that follows
 * accumulates on top of what this wrote. */
template <typename Spec>
struct NlrModeDLocalSlice {
    const typename Spec::ActiveData *actives = nullptr;
    bool                             device_visible = false;
    int                              n = 0;
    typename Spec::AccumData        *accums_out = nullptr;
};

/* The motion-target set for one call of a loop that writes neighbour
 * velocities: opened before any kernel, raised and closed after the writeback.
 * A loop without the trait touches none of this. */
template <typename Spec>
static void nlr_motion_targets_open(void)
{
    if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {
        const int num_local = ghost_get_num_local();
        if(num_local <= 0) {return;}   /* nothing this rank owns can be marked */
        if(gx_motion_target_ensure(num_local) != 0) {
            /* Without the set no bound can be raised, and a bound not raised
             * under-includes silently. */
            if(ThisTask == 0) {fprintf(stderr, "[%s] FATAL: no memory for the motion-target set.\n", Spec::loop_name); fflush(stderr);}
            endrun(90001040);
            return;
        }
        gx_motion_target_begin_call();
    }
}
template <typename Spec>
static void nlr_motion_targets_close(void)
{
    if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {gx_motion_target_consume();}
}
/* Opened where a call's inputs are final, closed when the call leaves scope
 * -- after its writeback, whichever return it takes. */
template <typename Spec>
struct NlrMotionTargetScope {
    NlrMotionTargetScope()  {nlr_motion_targets_open<Spec>();}
    ~NlrMotionTargetScope() {nlr_motion_targets_close<Spec>();}
};

template <typename Spec>
struct NlrPeerAnswerDeviceFused;

/* One rank's queries answered from its own tree, device-resident end to end.
 * Returns false only when the device buffers could not be had, in which case
 * nothing was evaluated and the caller answers on the host instead. */
template <typename Spec>
static bool nlr_mode_d_evaluate_single_rank(const typename Spec::DeviceContext& ctx,
                                            const typename Spec::CallScalars& cs,
                                            const int *active_idx_host,
                                            const int *active_slot_host,
                                            const double *radii_host,
                                            int n,
                                            unsigned int supply_mask,
                                            const GxDeviceTreeView& tree,
                                            typename Spec::AccumData *accums_out);

/* Holds the queries when they are built where the particles already are.
 *
 * Host-visible as well as device-visible, which is what lets one buffer serve
 * all three readers: the export walk wants each query's position and reach, the
 * envelopes ship the objects to peers, and the walk kernel evaluates them. A
 * plain vector could serve the first two and not the third; this serves all of
 * them, so the queries are built once and not copied again afterwards.
 *
 * That is narrower than "nothing is staged", and deliberately so: the indices
 * and radii the build reads are still staged, and the accumulators still are.
 * What this removes is the host walking P and CellP to fill the queries. */
template <typename Spec>
struct NlrDeviceBuiltActives {
    typename Spec::ActiveData *p = nullptr;
    NlrDeviceBuiltActives() = default;
    NlrDeviceBuiltActives(const NlrDeviceBuiltActives&) = delete;
    NlrDeviceBuiltActives& operator=(const NlrDeviceBuiltActives&) = delete;
    ~NlrDeviceBuiltActives() {
        if(p) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(p);}
    }
};

template <typename Spec, typename DeviceCtx>
static typename Spec::ActiveData *
nlr_build_self_actives_on_device(const neighbor_loop_args& args,
                                 const DeviceCtx& ctx,
                                 const double *radii,
                                 const typename Spec::CallScalars& cs,
                                 int n,
                                 NlrDeviceBuiltActives<Spec>& owner);

/* Which backend answers queries on this call.
 *
 * HostWalk is what Mode B has always done and is the default, so a caller that
 * says nothing gets exactly the path F1 proved unchanged.  DeviceFused is Mode D:
 * the same transport, the same wire, the same reply merge -- a different way of
 * turning a query into an accumulator.  The choice is the CALLER's because it is
 * the same choice that selected the dispatch path, and it was made collectively
 * there; the transport must not re-decide it per rank. */
enum class NlrEvalBackend { HostWalk, DeviceFused };

template <typename Spec, NlrEvalBackend Backend = NlrEvalBackend::HostWalk>
static void mode_b_remote_evaluate_into_buffer(
    const neighbor_loop_args& args,
    const double *radii,
    const typename Spec::CallScalars& cs,
    const typename Spec::DeviceContext& ctx,         /* caller-owned */
    unsigned int neighbor_type_mask,                  /* explicit caller param */
    typename Spec::AccumData *accums_out,             /* size = args.num_active; caller-owned */
    /* Only read on the DeviceFused backend, where it is the tree the collective
     * readiness decision was made against.  Unused by HostWalk. */
    const GxDeviceTreeView *fused_tree = nullptr)
{
    using ActiveData    = typename Spec::ActiveData;
    using AccumData     = typename Spec::AccumData;
    using DeviceCtx     = typename Spec::DeviceContext;
    using Envelope      = NlrQueryEnvelope<ActiveData>;
    using ReplyEnvelope = NlrReplyEnvelope<AccumData>;

    /* Mode predicates (single source of truth). Constexpr so the dead
     * branches are elided per instantiation. */
    static_assert(std::is_trivially_copyable<Envelope>::value,
        "NlrQueryEnvelope must be trivially-copyable for byte-level MPI transfer");
    static_assert(std::is_trivially_copyable<ReplyEnvelope>::value,
        "NlrReplyEnvelope must be trivially-copyable for byte-level MPI transfer");

    /* The fused backend reads the tree through this pointer, and a caller that
     * forgot it would not fail here -- it would fail inside a device kernel,
     * where the symptom is an unbounded read rather than a message.  The backend
     * is a template seam, so the next caller is the one this protects. */
    if constexpr (Backend == NlrEvalBackend::DeviceFused) {
        if(fused_tree == nullptr) {
            if(ThisTask == 0) {
                fprintf(stderr,
                    "[%s] FATAL: the fused evaluation backend was selected without the prepared "
                    "tree it walks. The caller that chose this backend owns that tree.\n",
                    Spec::loop_name);
                fflush(stderr);
            }
            endrun(90001030);
            return;
        }
    }

    const int N    = args.num_active;     /* may be 0 on this rank; collective entry */
    const int nt   = NTask;
    const int rank = ThisTask;
    /* CallScalars and DeviceContext passed in by caller:
     *   - Iterative: driver-owned via NlrIterDriver, populated once at iter-0 entry.
     *   - Non-iterative wrapper: locally-owned, populated + RAII-cleaned by wrapper.
     * Helper does NOT call populate_call_scalars or populate_device_context. */

    /* The queries. One snapshot serves both this rank's own pairs and the copies
     * shipped to peers, so neither can see a different epoch from the other.
     *
     * WHERE they are built is the difference between the backends. The host
     * walker builds them here, one Spec::load_active per active, and each of
     * those is a scattered read of P[i] and CellP[i] -- the canonical arrays
     * surfacing on the host, once per active, on the calls that have the most
     * actives. The fused backend has those same particles resident on the device
     * already, so it builds them there instead, into a buffer the host can still
     * read for the export walk and the envelopes.
     *
     * The objects are the same either way: same function, same inputs, same
     * epoch. The epoch holds because the queries are built from the ACTIVES, and
     * an active is already current when a neighbour loop is entered: the sync
     * point drifts every particle in every active time bin to All.Ti_Current
     * (core/run.cc), ActiveParticleList is built from those same bins, and every
     * Spec that can reach this backend draws its active list from that list.
     * Nothing advances All.Ti_Current again inside the step.
     *
     * That is the invariant the discovery pass below also rests on, so it is
     * worth naming rather than leaving implied: the record pass and the
     * evaluation pass must build the same query from the same particle, and they
     * would not if a drift between them could move an active. The touched-set
     * drift only advances particles that are behind, so an active that is
     * already current cannot be moved by it.
     *
     * Falling back to the host build when the device buffer cannot be had is
     * safe where declining would not be: it changes how the queries are filled,
     * never which collectives this rank enters, so the peers cannot tell. */
    std::vector<ActiveData> actives_host;
    NlrDeviceBuiltActives<Spec> actives_device;
    ActiveData *actives = nullptr;
    bool actives_are_device_visible = false;
    if(N > 0) {
        if constexpr (Backend == NlrEvalBackend::DeviceFused) {
            actives = nlr_build_self_actives_on_device<Spec>(args, ctx, radii, cs, N,
                                                             actives_device);
            actives_are_device_visible = (actives != nullptr);
        }
        if(actives == nullptr) {
            actives_host.resize(N);
            build_self_actives_host_pre_drift<Spec>(args, ctx, radii, cs,
                                                      actives_host.data());
            actives = actives_host.data();
        }
    }

    /* ---- ON THE HOST-WALK BACKEND the self stages run ONCE, BEFORE the peer
     * round loop: self candidate collection + drift + evaluation precede peer
     * evaluation, and only the PEER work streams in bounded rounds below.
     *
     * ⛔ THE FUSED BACKEND DOES NOT DO THIS, AND THE "self-before-peer WRITE
     * order" THIS COMMENT USED TO ASSERT IS NOT A CONTRACT ANYWHERE.  It walks
     * local and received queries in one batch, so a neighbour j can be reached
     * by both in the same launch.  Nothing is lost by that: whether a given
     * local particle is reached by one of this rank's queries or by one a peer
     * sent is an accident of where the domain boundary fell and carries no
     * physical meaning, Mode A has always evaluated local and imported
     * candidates together in a single fused loop, and the only pair kernels
     * this backend admits are those whose neighbour-side writes already go
     * through atomics.  The ordering below is therefore a description of the
     * host path, not a rule the fused path breaks.
     *
     * ORDERING NOTE (host path) — the reorder this introduces vs the prior single-pass form:
     * peer candidate COLLECTION now runs AFTER self EVALUATION (previously it ran
     * before). The peer-candidate walk keys membership on P[j].{Type,Pos,Mass}:
     *   - Type/Pos are NEVER mutated by a pair kernel (only drift writes Pos, and
     *     it is applied identically per round) → position/type membership is
     *     reorder-invariant for ALL specs.
     *   - Mass is the only membership field a pair kernel can change (mass-flux
     *     specs: MFV hydro; feedback add; sink swallow remove). For mass-
     *     PRESERVING specs (density, gradients, MFM hydro_force) the reorder is
     *     EXACTLY neutral. For mass-MUTATING specs a j crossing the Mass>0 gate
     *     during self eval could shift its peer-candidate membership — but those
     *     specs are ALREADY order-dependent, and this stays within that regime,
     *     not a new exactness contract.
     * NO generic j-write-exactness claim is made. Any mass-mutating spec routed
     * through streamed Mode-B at multi-round MUST be re-audited (evrard, the
     * validated case, is mass-preserving MFM). */

    /* Every Mode-B loop routes via targeted export: the per-type node band is
     * cross-rank-fresh and dominates every radius_policy (mode_b_local_walker.h). */
    const double jscale = nlr_spec_symmetric_j_radius_scale<Spec>();

    /* Targeted-export reverse map: topnode indices are stable between builds →
     * build ONCE, reuse for the fused walk. */
    ModeBTopleafMap topleaf_map;   /* shared read-only map, built once */
    ModeBExportSink export_sink;   /* per-query export sink (write-only during the walk) */
    if(N > 0 && nt > 1) { export_sink.ensure_size(nt); topleaf_map.build(); }

    /* Fused-walk export CSR (targeted specs): per active, its per-peer export
     * node-lists, staged ONCE by the fused self walk and marshalled (no second
     * walk) by the peer round loop below. Active-ordered (csr_rec_off), peer-ascending
     * within an active, node order = walk append order. Sized O(total targeted
     * exports), NOT O(NTask*N). */
    struct FusedExportRec { int peer; int node_off; int n_nodes; };
    std::vector<int> csr_rec_off;                    /* size N+1: active aa -> [off[aa],off[aa+1]) recs */
    std::vector<FusedExportRec> csr_recs;
    std::vector<int> csr_nodes;

    const int modeb_nthreads = nlr_modeb_omp_nthreads();

    /* Stage 3: collect SELF candidates PRE-DRIFT. For targeted specs this is the
     * FUSED legacy-mode==0 walk — candidates + export CSR in ONE traversal, keyed
     * on the frozen actives[] snapshot (== the query the receiver walks) so the
     * candidate / export / receiver walks share one query SSOT. Broadcast specs
     * keep the plain candidate walk (they have no export walk to fuse). */
    std::vector<std::vector<int>> cand_self_tree;
    if(N > 0) {
                if(nt > 1) {
            /* want_cands: the tree walk needs candidates; the export CSR
             * for the round loop is built regardless,
             * so the fused walk runs with cand_out=nullptr there. */
            /* The export CSR is built either way; the candidate list is
             * only wanted by a backend that will later walk it.  A fused
             * backend answers its own actives from the tree directly, so
             * collecting them here would be building a list to throw away. */
            const bool want_cands = (Backend == NlrEvalBackend::HostWalk);
            if(want_cands) cand_self_tree.assign(N, std::vector<int>{});
            csr_rec_off.assign(N + 1, 0);
            /* Thread the fused self walk above the work threshold. Each thread
             * walks its actives into its OWN export sink + its OWN CSR segment
             * (no shared push, no lock); a serial prefix-sum then assembles the
             * active-ordered CSR BYTE-IDENTICALLY to the serial build
             * (per-active walk order fixed; peers ascending within an active;
             * node order = walk append order).
             *
             * The threshold is the only thing that decides this. Whether the
             * walk also collects candidates does not: it changes what the walk
             * records, not how its actives divide between threads, and the
             * per-active work is the traversal either way. This test once also
             * required candidates, back when a walk without them was a
             * validation pass whose speed was irrelevant; a walk without them
             * is now the production path that exports to peers, and it has the
             * most actives of any of them. */
            const bool use_omp_self = nlr_modeb_use_omp(N, modeb_nthreads);
            if(use_omp_self) {
                struct AaMeta { int tid; int rec_off; int n_recs; int node_off; int n_nodes; };
                std::vector<AaMeta> meta(N);
                std::vector<ModeBExportSink> tsink(modeb_nthreads);
                for(auto& s : tsink) s.ensure_size(nt);
                std::vector<std::vector<FusedExportRec>> trecs(modeb_nthreads);
                std::vector<std::vector<int>> tnodes(modeb_nthreads);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, MODEB_OMP_CHUNK_ACTIVE)
#endif
                for(int aa = 0; aa < N; aa++) {
#ifdef _OPENMP
                    const int tid = omp_get_thread_num();
#else
                    const int tid = 0;
#endif
                    ModeBExportSink& sink = tsink[tid];
                    std::vector<FusedExportRec>& lrecs = trecs[tid];
                    std::vector<int>& lnodes = tnodes[tid];
                    AaMeta& m = meta[aa];
                    m.tid = tid; m.rec_off = (int)lrecs.size(); m.node_off = (int)lnodes.size();
                    m.n_recs = 0; m.n_nodes = 0;
                    const double h_q = (double)actives[aa].h_search;
                    if(h_q <= 0) continue;
                    double pos_arr[3] = {(double)actives[aa].pos[0],
                                         (double)actives[aa].pos[1],
                                         (double)actives[aa].pos[2]};
                    sink.clear_all();
                    /* No candidate sink when the backend will not walk one:
                     * cand_self_tree is left empty in that case, so taking its
                     * element address would be out of bounds. Both branches
                     * guard it, and this one is now the branch that meets the
                     * case, since a walk that collects nothing threads too. */
                    std::vector<int>* cand_ptr = nullptr;
                    if(want_cands) {
                        cand_ptr = &cand_self_tree[aa];
                        if(cand_ptr->capacity() == 0) cand_ptr->reserve(64);
                    }
                    mode_b_walk_and_export(pos_arr, h_q, neighbor_type_mask,
                                            Spec::search_mode, Spec::radius_policy,
                                            cand_ptr, topleaf_map, sink, jscale);
                    for(int p = 0; p < nt; p++) {
                        if(p == rank) continue;
                        const std::vector<int>& nodes = sink.nodes_per_peer[p];
                        const int nn = (int)nodes.size();
                        if(nn == 0) continue;
                        FusedExportRec rec;
                        rec.peer = p; rec.node_off = (int)lnodes.size(); rec.n_nodes = nn;
                        lrecs.push_back(rec);
                        lnodes.insert(lnodes.end(), nodes.begin(), nodes.end());
                        m.n_recs++;
                        m.n_nodes += nn;
                    }
                }
                /* Deterministic active-ordered merge. */
                size_t total_recs = 0, total_nodes = 0;
                for(int aa = 0; aa < N; aa++) { total_recs += meta[aa].n_recs; total_nodes += meta[aa].n_nodes; }
                csr_recs.resize(total_recs);
                csr_nodes.resize(total_nodes);
                int rec_cursor = 0, node_cursor = 0;
                for(int aa = 0; aa < N; aa++) {
                    csr_rec_off[aa] = rec_cursor;
                    const AaMeta& m = meta[aa];
                    const std::vector<FusedExportRec>& lrecs = trecs[m.tid];
                    const std::vector<int>& lnodes = tnodes[m.tid];
                    for(int rr = 0; rr < m.n_recs; rr++) {
                        FusedExportRec rec = lrecs[m.rec_off + rr];
                        const int local_node_off = rec.node_off;   /* thread-segment-relative */
                        rec.node_off = node_cursor;
                        for(int q = 0; q < rec.n_nodes; q++)
                            csr_nodes[node_cursor++] = lnodes[local_node_off + q];
                        csr_recs[rec_cursor++] = rec;
                    }
                }
                csr_rec_off[N] = rec_cursor;
            } else {
                for(int aa = 0; aa < N; aa++) {
                    csr_rec_off[aa] = (int)csr_recs.size();
                    const double h_q = (double)actives[aa].h_search;
                    if(h_q <= 0) continue;
                    double pos_arr[3] = {(double)actives[aa].pos[0],
                                         (double)actives[aa].pos[1],
                                         (double)actives[aa].pos[2]};
                    export_sink.clear_all();
                    std::vector<int>* cand_ptr = nullptr;
                    if(want_cands) {
                        cand_ptr = &cand_self_tree[aa];
                        if(cand_ptr->capacity() == 0) cand_ptr->reserve(64);
                    }
                    mode_b_walk_and_export(pos_arr, h_q, neighbor_type_mask,
                                            Spec::search_mode, Spec::radius_policy,
                                            cand_ptr, topleaf_map, export_sink, jscale);
                    /* stage this active's per-peer exports into the CSR */
                    for(int p = 0; p < nt; p++) {
                        if(p == rank) continue;
                        const std::vector<int>& nodes = export_sink.nodes_per_peer[p];
                        const int nn = (int)nodes.size();
                        if(nn == 0) continue;
                        FusedExportRec rec;
                        rec.peer = p; rec.node_off = (int)csr_nodes.size(); rec.n_nodes = nn;
                        csr_recs.push_back(rec);
                        csr_nodes.insert(csr_nodes.end(), nodes.begin(), nodes.end());
                    }
                }
                csr_rec_off[N] = (int)csr_recs.size();
            }
        } else {
            /* single rank: no peers to export to → plain candidate walk. */
            collect_candidates_pre_drift<Spec>(args, radii,
                                                neighbor_type_mask,
                                                DispatchPath::ModeB_HostWalker,
                                                cand_self_tree);
        }
    }

    /* Self drift (split from the former self+peer union drift; peer candidates
     * drift per-round below). drift_particle is idempotent to All.Ti_Current
     * (constant across the helper), so a j that is both a self- and peer-
     * candidate drifts once — identical to the old combined union drift. */
    if constexpr (Backend == NlrEvalBackend::HostWalk) {
        /* The same query the walk above used: the frozen actives[] on the multi-rank fused walk, the
         * particle and its radius on the single-rank walk. */
        if (N > 0) drift_and_filter_candidates<Spec>(cand_self_tree, neighbor_type_mask, [&](int aa, double *pos, double &h_q) {
            if(nt > 1) {
                pos[0] = (double)actives[aa].pos[0]; pos[1] = (double)actives[aa].pos[1]; pos[2] = (double)actives[aa].pos[2];
                h_q = (double)actives[aa].h_search;
            } else {
                const int i = args.active_list[aa];
                pos[0] = (double)args.P[i].Pos[0]; pos[1] = (double)args.P[i].Pos[1]; pos[2] = (double)args.P[i].Pos[2];
                h_q = radii[aa];
            }
        });
    }

    /* Stage 8: answer THIS rank's own queries -> accums_out.
     *
     * The host backend evaluates the list it collected above.  The fused backend
     * does NOT answer them here: it carries them into the round loop below and
     * walks them in the same launch as the queries it received, because the two
     * ask the same question of the same tree and a half-empty launch apiece
     * leaves the device waiting twice.  Which of this rank's particles a given
     * query reaches is an accident of where the domain boundary fell, so there
     * is nothing to separate.  `local_cursor` tracks how many have been taken.
     *
     * The fused backend also skips the collected-candidate drift, because it
     * collected no candidates -- not because everything is already current. */
    [[maybe_unused]] int local_cursor = 0;   /* read only on the fused path */
    if constexpr (Backend == NlrEvalBackend::HostWalk) {
        if(N > 0) {
            evaluate_pairs_post_drift<Spec>(ctx, actives, N,
                                              cand_self_tree, accums_out, cs, EvalOMPPolicy::AllowProduction);
        }
        local_cursor = N;   /* answered here; the round loop carries nothing */
    }

    /* ---- PEER round loop (streaming). Build a CommChunkSize-bounded batch of
     * query envelopes from actives[cursor..N), exchange+evaluate+reply, advance
     * cursor, iterate until every rank is drained (legacy do/while +
     * Allreduce(ndone), code_block_xchange_perform_ops.h:9-191). Bounds total
     * in-flight export envelopes so forced-Mode-B is memory-safe at large
     * N_active. This caps the SENDER side only; the full large-N fix needs
     * receiver group staging.
     *
     * ALL Mode-B loops get TARGETED export (gate retired): walk the
     * local tree per active and export ONLY to peers whose remote subtree the query
     * reaches, carrying the exported start-nodes so the receiver resumes a bounded
     * walk. The per-type node band prunes the SYMMETRIC reach and is cross-rank-fresh
     * (via force_update_hmax), so the sender bounds every loop's reach on remote
     * peers. Self-pair handled above; self entry stays empty. */
    /* jscale and the exporter are hoisted above Stage 3 for the fused walk. */
    /* How much this round may carry: one communication chunk, divided by the size of a
     * query and its reply. NlrQueryEnvelope fuses what used to be sent as separate index,
     * node-list and active records, so counting envelopes here counts what the older code
     * counted as export entries. All.CommChunkSize is the parameterfile chunk size. */
    constexpr size_t kReplyBytes = sizeof(ReplyEnvelope);
    const long long kEnvPairBytes = (long long)sizeof(Envelope) + (long long)kReplyBytes;
    long long bunch = ((long long)All.CommChunkSize * 1024 * 1024) /
                      (kEnvPairBytes > 0 ? kEnvPairBytes : 1);
    if(bunch < 1) bunch = 1;

    long long diag_export_qr = 0, diag_node_appends = 0;   /* scalar export volume (NLR diag) */
    int cursor = 0;
    int ndone  = 0;

    do {
        long long round_env_count = 0;
        std::vector<std::vector<Envelope>> queries_per_peer(nt);

        /* Stage 2 (bounded): fill this round's export batch from actives[cursor..).
         * MARSHAL from the fused walk's export CSR — NO walk here (the fused
         * legacy mode==0 walk already ran in Stage 3). Per active, MEASURE its
         * envelope count from the CSR, then commit-or-stop atomically —
         * measure-then-commit gives legacy's all-or-nothing-per-particle
         * rollback (an active never lands half its chunks in one round). An
         * active whose OWN set exceeds `bunch` ships in a solo oversized round
         * (graceful; loud diag) instead of aborting. */
        if(N > 0 && nt > 1) {
                        int aa = cursor;
            for(; aa < N; aa++) {
                const int r0 = csr_rec_off[aa], r1 = csr_rec_off[aa + 1];
                /* envelopes this active would add across all peers */
                long long add = 0;
                for(int r = r0; r < r1; r++)
                    add += (csr_recs[r].n_nodes + NODELISTLENGTH - 1) / NODELISTLENGTH;
                if(round_env_count > 0 && round_env_count + add > bunch) break; /* defer to next round */
                if(round_env_count == 0 && add > bunch) {
                    nlr_warn_once_rank0("modeb_oversize_active",
                        "[mode_b B2a caller=%s] single active's export set (%lld envelopes, "
                        "~%lld bytes) exceeds CommChunkSize bunch (%lld envelopes); shipping a solo "
                        "oversized round — cap ineffective for this call (raise CommChunkSize).",
                        Spec::loop_name, add, add * kEnvPairBytes, bunch);
                }
                /* commit: chunked envelopes per exported peer (CSR records are peer-ascending). */
                int rr = r0;
                for(int p = 0; p < nt; p++) {
                    if(p == rank) continue;
                    if(rr < r1 && csr_recs[rr].peer == p) {
                        const FusedExportRec& rec = csr_recs[rr];
                        const int* nd = &csr_nodes[rec.node_off];
                        const int nn = rec.n_nodes;
                        diag_export_qr++; diag_node_appends += nn;
                        /* Chunk into NODELISTLENGTH-sized records (legacy opens a
                         * fresh export slot when a NodeList fills). Chunks cover
                         * disjoint subtrees → the slot-keyed reply merge sums their
                         * partial results without double counting. All chunks of a
                         * (query,peer) group land in THIS round (all-or-nothing
                         * above), so each group stays contiguous. */
                        for(int c = 0; c < nn; c += NODELISTLENGTH) {
                            Envelope env;
                            env.origin_slot = aa;
                            env.origin_rank = rank;
                            int cnt = 0;
                            for(; cnt < NODELISTLENGTH && (c + cnt) < nn; cnt++) {
                                env.NodeList[cnt] = nd[c + cnt];
                            }
                            env.n_nodes = cnt;
                            env.reserved_wire_padding = 0;
                            for(int t = cnt; t < NODELISTLENGTH; t++) env.NodeList[t] = -1;
                            env.active = actives[aa];
                            queries_per_peer[p].push_back(env);
                        }
                        rr++;
                        continue;
                    }
                }
                round_env_count += add;
            }
            cursor = aa;
        } else {
            cursor = N;   /* nothing to export (N==0 or single rank) */
        }

        /* Stage 4: exchange queries. Every rank participates even if it queued
         * 0 this round (peers may target this rank's pool); the Allreduce(ndone)
         * at the round's end keeps every rank's round count equal, so the
         * exchange stays balanced. begin() posts ALL query-payload Isends + ALL
         * reply Irecvs up front; the receiver then stages, evaluates, and
         * answers incoming queries in memory-bounded WHOLE-PEER groups (the
         * legacy import sub-chunk loop, code_block_xchange_perform_ops.h:96-174)
         * instead of materializing every peer's payload at once. Group budget =
         * the same All.CommChunkSize that sizes the sender bunch. The bound covers
         * TRANSPORT payloads only (envelopes + replies) — candidate vectors and
         * pair-kernel scratch scale with the group's query content, not with
         * NTask. Peers are consumed in ascending rank order, so the
         * concatenated per-group evaluation sequence — and the post-loop reply
         * merge — keep the exact pre-group order (matters for j-writing specs). */
        using XReply = ReplyEnvelope;
        ModeBBoundedExchange<Envelope, XReply> xch;
        {
            xch.begin(queries_per_peer);
        }
        const size_t group_budget_bytes =
            (size_t)((double)All.CommChunkSize * 1024.0 * 1024.0);

        std::vector<int> group_peers;
        std::vector<std::vector<Envelope>> group_queries;
        while(true) {
            bool have_group;
            {
                have_group = xch.next_group(group_budget_bytes, group_peers, group_queries);
            }
            if(!have_group) break;

    /* Stage 5: flatten THIS GROUP's envelopes and build the provenance map.
     * provenance[k] carries:
     *   - source_gidx / source_qi: where in group_queries[][] this k came
     *     from (used for unflattening replies back into per-peer arrays);
     *     source_peer = the peer's rank, for diagnostics.
     *   - origin_slot / origin_rank: copied from the received envelope; ride
     *     into the REPLY envelope so the active rank can merge by slot
     *     without relying on transport ordering (symmetric
     *     query and reply envelopes). */
    std::vector<ActiveData> peer_actives;
    struct Provenance {
        int source_gidx; int source_peer; int source_qi;
        int origin_slot; int origin_rank;
    };
    std::vector<Provenance> peer_provenance;
    /* Parallel to peer_actives[k]: the exported start-node list carried in the
     * received envelope, flattened K*NODELISTLENGTH, so the receiver walk
     * (collect_candidates_for_remote_queries) resumes from those nodes. */
    std::vector<int> peer_nodelist_flat;
    std::vector<int> peer_nnodes;
    size_t total_recv = 0, total_recv_bytes = 0;
    for(size_t gi = 0; gi < group_peers.size(); gi++) {
        total_recv += group_queries[gi].size();
        total_recv_bytes += group_queries[gi].size() * (sizeof(Envelope) + sizeof(XReply));
    }
    peer_actives.reserve(total_recv);
    peer_provenance.reserve(total_recv);
    peer_nodelist_flat.reserve(total_recv * NODELISTLENGTH);
    peer_nnodes.reserve(total_recv);
    for(size_t gi = 0; gi < group_peers.size(); gi++) {
        const int p = group_peers[gi];
        for(int qi = 0; qi < (int)group_queries[gi].size(); qi++) {
            const Envelope& env = group_queries[gi][qi];
            /* Sanity: envelope's origin_rank should equal sender p. */
            if(env.origin_rank != p) {
                fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                        "envelope origin_rank=%d but received from peer=%d "
                        "(qi=%d). Transport corruption?\n",
                        rank, Spec::loop_name, env.origin_rank, p, qi);
                fflush(stderr);
                endrun(81220);
            }
            /* Soft bad-stop on corruption but DO NOT skip: the reply array is sized by
             * recv_counts[p], so omitting this qi would leave a default-initialized reply
             * (bogus origin_rank=0 could be merged into slot 0 on the peer). Keep the entry
             * with the transport-consistent origin rank p (== env.origin_rank when clean);
             * the run drains at the next phase poll with reply choreography intact. */
            peer_actives.push_back(env.active);
            peer_provenance.push_back({(int)gi, p, qi, env.origin_slot, p});
            peer_nnodes.push_back(env.n_nodes);
            for(int t = 0; t < NODELISTLENGTH; t++) peer_nodelist_flat.push_back(env.NodeList[t]);
        }
    }
    /* Note: receiver-side binding of peer_actives (and per-eval-pass rebinding
     * of both self actives and peer actives) is performed inside
     * evaluate_pairs_post_drift, gated on
     * nlr_spec_has_bind_active_to_eval_context_v<Spec>. The bind runs once
     * per active per eval pass with the EXACT eval ctx. No standalone
     * post-flatten rebind needed here. */

    /* Stages 6, 7 and 9 -- find this round's peer neighbours, bring them current,
     * and evaluate -- are one backend operation.  See NlrPeerAnswerHostWalk. */
    const int K = (int)peer_actives.size();
    std::vector<AccumData> peer_replies(K);
    if constexpr (Backend == NlrEvalBackend::HostWalk) {
        NlrPeerAnswerHostWalk<Spec>::answer(ctx, cs, peer_actives,
                                            peer_nodelist_flat, peer_nnodes,
                                            neighbor_type_mask,
                                            peer_replies);
    } else {
        /* Carry whatever of this rank's own queries are still unplaced into the
         * same two launches as this group's received ones -- but only as many as
         * fit under the accumulator high-water the unfused shape already
         * required.
         *
         * WHY THERE IS A BOUND AT ALL.  Answering the locals separately meant
         * their N accumulators were allocated and freed BEFORE this group's K
         * were, so the call's peak was max(N, K).  One batch makes both live at
         * once, which would be N + K.  Holding L + K <= max(N, K) -- i.e.
         * L <= N - K -- keeps the peak exactly what it was, so fusing cannot
         * make a call that used to fit stop fitting.
         *
         * WHAT IT COSTS.  When a group brings MORE queries than this rank has
         * actives (K >= N) the headroom is zero and the locals wait, which for
         * that group is the unfused launch count rather than a regression.  In
         * the ordinary case K is bounded by the group budget and far below N, so
         * the headroom covers essentially every local and one group takes them
         * all.  Whether spreading them across groups instead would fill the
         * device better is a measurement, not something to assume here. */
        NlrModeDLocalSlice<Spec> local_slice;
        const int local_headroom = (N > K) ? (N - K) : 0;
        const int local_remaining = N - local_cursor;
        const int local_take = (local_remaining < local_headroom) ? local_remaining : local_headroom;
        if(local_take > 0) {
            local_slice.actives        = actives + local_cursor;
            local_slice.device_visible = actives_are_device_visible;
            local_slice.n              = local_take;
            local_slice.accums_out     = accums_out + local_cursor;
        }
        NlrPeerAnswerDeviceFused<Spec>::answer(ctx, cs, *fused_tree, peer_actives,
                                               peer_nodelist_flat, peer_nnodes,
                                               neighbor_type_mask,
                                               peer_replies, local_slice);
        local_cursor += local_slice.n;
    }
    /* Stage 10 (per group): build reply envelopes (origin_slot/rank copied from
     * each received query envelope), unflatten into per-peer arrays via the
     * provenance map, then send THIS group's replies — after which the group's
     * buffers are released (the whole point of the group staging). Production
     * replies share one payload type (XReply), so there is one
     * reply exchange per group, never two with identical MPI tags. */
    {
        std::vector<std::vector<XReply>> replies_for_group(group_peers.size());
        for(size_t gi = 0; gi < group_peers.size(); gi++) {
            replies_for_group[gi].assign(group_queries[gi].size(), XReply{});
        }
        for(int k = 0; k < K; k++) {
            const Provenance& pv = peer_provenance[k];
            XReply& re = replies_for_group[pv.source_gidx][pv.source_qi];
            re.origin_slot = pv.origin_slot;
            re.origin_rank = pv.origin_rank;
            re.accum = peer_replies[k];
        }
        /* send_group_replies posts the reply Isends and waits them in finish();
         * move the group's reply buffers into the exchange so they outlive the Isend. */
        xch.send_group_replies(group_peers, std::move(replies_for_group));
    }
        }   /* end whole-peer group loop */

        /* ⛔ ANY OF THIS RANK'S OWN QUERIES STILL UNPLACED ARE ANSWERED HERE,
         * BEFORE THE REPLY DRAIN BELOW -- NEVER AFTER THE ROUND LOOP.
         *
         * accums_out[slot] is ASSIGNED by a local evaluation and MERGED INTO by
         * an incoming reply.  A rank that received nothing this round got no
         * group above, so without this its locals would still be unanswered when
         * the drain merges its replies, and the assignment would then land on top
         * of them and discard them.  Putting the flush here rather than after the
         * do/while is what makes that impossible, and it is also the whole
         * no-received-queries case: a rank that never receives anything from
         * anyone answers all of its own work here, on the first round.
         *
         * Normally a no-op: the group loop above already took them. */
        if constexpr (Backend == NlrEvalBackend::DeviceFused) {
            if(local_cursor < N) {
                std::vector<AccumData> no_replies;
                NlrModeDLocalSlice<Spec> local_slice;
                local_slice.actives        = actives + local_cursor;
                local_slice.device_visible = actives_are_device_visible;
                local_slice.n              = N - local_cursor;
                local_slice.accums_out     = accums_out + local_cursor;
                NlrPeerAnswerDeviceFused<Spec>::answer(ctx, cs, *fused_tree,
                                                       std::vector<ActiveData>{},
                                                       std::vector<int>{}, std::vector<int>{},
                                                       neighbor_type_mask,
                                                       no_replies, local_slice);
                local_cursor = N;
            }
        }

        /* Stage 11: drain the exchange (remaining query-payload Isends + all
         * pre-posted reply Irecvs, byte-count-asserted), then merge replies
         * into accums_out by envelope.origin_slot. Pinned deterministic order:
         * ascending peer rank, ascending qi — identical to the pre-group
         * single-shot merge (the local contribution is already in accums_out:
         * either the group loop above placed it, or the flush just did).
         * Asserts each reply envelope's origin_rank == ThisTask. */
        {
            auto recv_replies = [&]{
                return xch.finish();
            }();
            for (int p = 0; p < nt; p++) {
                if (p == rank) continue;
                const int q_to_p = xch.sent_counts[p];
                for (int qi = 0; qi < q_to_p; qi++) {
                    const XReply& re = recv_replies[p][qi];
                    if(re.origin_rank != rank) {
                        fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                                "reply envelope origin_rank=%d != ThisTask=%d from "
                                "peer=%d qi=%d. Transport/peer-side corruption?\n",
                                rank, Spec::loop_name, re.origin_rank, rank, p, qi);
                        fflush(stderr);
                        /* reply exchange already completed; soft bad-stop + skip the corrupt reply. */
                        endrun(81223); continue;
                    }
                    const int slot = re.origin_slot;
                    if(slot < 0 || slot >= N) {
                        fprintf(stderr, "[neighbor_loop_runner ABORT rank=%d caller=%s] "
                                "reply envelope slot %d out of range [0,%d) from peer %d.\n",
                                rank, Spec::loop_name, slot, N, p);
                        fflush(stderr);
                        /* skip: continuing would merge into accums_out[slot] with slot OOB. */
                        endrun(81224); continue;
                    }
                    Spec::merge_accum(accums_out[slot], re.accum);
                }
            }
        }

        /* Termination (legacy 186-191): done when my cursor drained; the SUM
         * Allreduce makes every rank run the SAME number of rounds so the
         * per-round query/reply exchanges stay collective-balanced. Ranks that
         * finished sending keep entering as 0-send receivers for peers still
         * draining. */
        int ndone_flag = (cursor >= N) ? 1 : 0;
        {
            MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
    } while(ndone < NTask);



    /* End of helper. Caller decides whether to call apply_active_writeback
     * (final-only for iterative; per-call for non-iterative) and whether
     * to emit active_dumps (after writeback for non-iter; deferred for
     * iter until proper diagnostic plumbing lands). */
}

/* ============================================================================
 * run_mode_b_remote_impl<Spec> — non-iterative thin wrapper.
 *
 * Allocates per-call AccumData buffer, calls helper, runs final
 * apply_active_writeback + active_dumps emit. Same epoch order and same
 * writeback per active as the earlier monolithic impl.
 * ========================================================================== */
template <typename Spec>
static void run_mode_b_remote_impl(const neighbor_loop_args& args, const double *radii)
{
    using AccumData    = typename Spec::AccumData;
    using DeviceCtx    = typename Spec::DeviceContext;

    const int N = args.num_active;

    /* Capture CallScalars per call (non-iter wrapper). */
    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    /* Build per-call DeviceContext + RAII cleanup guard (helper
     * no longer builds its own ctx; caller does). Behavior byte-equivalent
     * to the earlier monolith — guard runs Spec::cleanup_device_context at
     * scope exit, after the writeback loop completes. */
    DeviceCtx ctx;
    ctx.P         = args.P;
    ctx.CellP     = args.CellP;
    ctx.num_total = args.num_total;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; "
                  "captured by value into Kokkos device lambdas");
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);
    /* ctx.populate_failed is deliberately NOT honoured on this multi-rank leg:
     * returning here would desync the query/reply transport the peers enter
     * below, and answering their queries needs the very buffers that failed.
     * What makes that safe is that the readers tolerate an absent buffer --
     * a Spec whose load_neighbor indexes a hook-owned array must test it
     * first (sink_feed's binary_merge_eligible is the worked example, and it
     * is sized by num_total rather than by the active count, so it is not
     * small). The stop is already requested with the buffer named, and the
     * loop contributes nothing further. */

    /* Caller-owned per-call AccumData buffer; helper writes into it.
     * Explicit Spec::zero_accum per slot: defensive
     * against future evaluate_pairs variants that don't zero internally.
     * Today evaluate_pairs_post_drift calls Spec::zero_accum at line 812
     * before walking candidates, so this outer zero is redundant — but
     * makes the AccumData contract explicit at the caller level (sentinel-
     * bearing AccumData like sink_feed's Sink_PotentialMinimumOfNeighbors
     * needs Spec::zero_accum, NOT default-construction). */
    std::vector<AccumData> accums_self(N);
    for (int aa = 0; aa < N; aa++) {
        Spec::zero_accum(accums_self[aa]);
    }

    /* Helper runs Stages 1-12; writeback is the wrapper's responsibility
     * (preserves the earlier timing). */
    mode_b_remote_evaluate_into_buffer<Spec>(args, radii, cs, ctx,
                                                       nlr_effective_neighbor_type_mask(args, Spec::neighbor_type_mask),
                                                       (N > 0) ? accums_self.data() : nullptr);

    /* Stage 12 final: writeback per active. */
    {
        for(int aa = 0; aa < N; aa++) {
            Spec::apply_active_writeback(args, aa, args.active_list[aa], accums_self[aa]);
        }
    }

}

template <typename Spec>
static void run_mode_b_remote(const neighbor_loop_args& args, const double *radii) {
    run_mode_b_remote_impl<Spec>(args, radii);
}

#ifdef NEIGHBOR_LOOP_MODE_D
/* Whether a loop can be answered by the fused device walk at all: a one-way
 * search whose pair kernel may run concurrently, and whose active rebinding
 * the device evaluation does not perform.  Compile-time, so the machinery is
 * instantiated only for loops that can reach it. */
template <typename Spec>
constexpr bool nlr_spec_mode_d_eligible_v =
    (Spec::search_mode == MODE_B_SEARCH_ONEWAY) &&
    (nlr_spec_modeb_eval_omp<Spec>() != ModeBEvalOMP::SerialOnly) &&
    !nlr_spec_has_bind_active_to_eval_context_v<Spec>;

/* Mode D for a single-pass loop: the same transport and the same device walk
 * the iterative driver uses, without the per-iteration compaction.  The
 * particles are read from the device-resident arena, as on Mode A, because
 * that is where the walk evaluates them; the tree was prepared and voted on at
 * the dispatch site, once, before this is entered. */
template <typename Spec>
static void run_mode_d(const neighbor_loop_args& args, const double *radii,
                       const GxDeviceTreeView& tree)
{
    using AccumData = typename Spec::AccumData;
    using DeviceCtx = typename Spec::DeviceContext;

    const int N = args.num_active;

    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(args.num_total, args.P, args.CellP);

    DeviceCtx ctx;
    ctx.P         = gpu_particles_arena_P();
    ctx.CellP     = (args.CellP != nullptr) ? gpu_particles_arena_CellP() : nullptr;
    ctx.num_total = args.num_total;
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }

    std::vector<AccumData> accums(N);
    for(int aa = 0; aa < N; aa++) {Spec::zero_accum(accums[aa]);}

    bool evaluated = true;
    {
        /* The guard releases the Spec's context after the writeback below, the
         * order the other wrappers keep. */
        NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);
        /* A hook that could not stage its buffers has already asked for the
         * stop. Alone, this rank enters no collective, so it returns here as
         * Mode A does and the request drains at the caller's next poll. With
         * peers it proceeds, for the reason the host transport states at its
         * own population site: returning would desync the query/reply exchange
         * the peers enter, and every reader of a hook-owned buffer tests for
         * its absence first. */
        if(ctx.populate_failed && NTask == 1) {
            gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);
            return;
        }
        if(NTask > 1) {
            /* Every rank enters the transport, with or without actives of its own. */
            mode_b_remote_evaluate_into_buffer<Spec, NlrEvalBackend::DeviceFused>(
                args, radii, cs, ctx,
                nlr_effective_neighbor_type_mask(args, Spec::neighbor_type_mask),
                (N > 0) ? accums.data() : nullptr, &tree);
        } else if(N > 0) {
            /* Each entry's slot into the call-level staging: carried by a
             * narrowed list, its position otherwise. */
            std::vector<int> slots(N);
            for(int aa = 0; aa < N; aa++) {slots[aa] = args.active_call_slot ? args.active_call_slot[aa] : aa;}
            evaluated = nlr_mode_d_evaluate_single_rank<Spec>(
                ctx, cs, args.active_list, slots.data(), radii, N,
                nlr_effective_neighbor_type_mask(args, Spec::neighbor_type_mask),
                tree, accums.data());
        }
        if(evaluated) {
            for(int aa = 0; aa < N; aa++) {
                Spec::apply_active_writeback(args, aa, args.active_list[aa], accums[aa]);
            }
        }
    }
    gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);

    /* Out of device memory is a reason to answer differently, never to answer
     * less: the host walker takes the call.  Only the single-rank shape can get
     * here -- with peers the transport above has already been entered -- so
     * the change of backend is invisible to every other rank. */
    if(!evaluated) {run_mode_b_local<Spec>(args, radii);}
}
#endif

/* ============================================================================
 * External-CSR staging helpers (hydro corridor support).
 *
 * When args.external_csr is non-null, Mode A skips gpu_ngb_list_build and
 * instead stages the caller's host CSR into Kokkos memory shaped like a
 * gpu_neighbor_list_t — so the rest of run_mode_a is path-agnostic. Only
 * the build site (replaced with this helper) and the free site (the
 * matching helper below) differ between the two paths.
 *
 * The spatial-index fields of gnl (d_tiles / d_bvh / d_pool / ntiles /
 * bvh_root) stay zero/null because the pair_kernel does not read them (it
 * uses nearest_xyz, which reads All.BoxSize_* via the AllDeviceMirror).
 *
 * The runner OWNS the SharedSpace/DeviceSpace allocations made here and
 * frees them in nlr_free_external_csr_gnl(). It does NOT free the caller's
 * host buffers (active_indices / offsets / neighbors). Contract: caller
 * keeps host CSR alive for the duration of every run_neighbor_loop call
 * that injects it; the corridor design owns CSR across multiple consumers
 * by holding it in the gizmo_sym_* globals. */
/* Exhaustion is reported by returning NULL, so a caller NULL-check can
 * controlled-stop with attribution (the label appears in the memory ledger)
 * instead of a hard terminate. */
static void *nlr_shared_alloc_bytes(size_t bytes, const char *label)
{
    if(bytes == 0) { return NULL; }
    return gizmo_gpu_alloc_shared(bytes, label);
}

/* One wording for every runner staging buffer that could not be had. The stop is
 * drained at the next phase boundary; until then the loop simply leaves its actives
 * untouched. Nothing here returns early past a collective — the ghost-writeback
 * pair around the kernel work still fires, so no peer rank is left waiting. */
static void nlr_stop_no_staging_memory(const char *loop_name, const char *what, size_t bytes)
{
    char msg[256];
    snprintf(msg, sizeof(msg),
             "%s: could not allocate %s (%.1f MB); this loop leaves its actives untouched",
             loop_name ? loop_name : "neighbour loop", what,
             (double) bytes / (1024.0 * 1024.0));
    gizmo_request_controlled_stop(7713, msg, __FILE__, __LINE__, __FUNCTION__);
}

static inline bool
nlr_stage_external_csr_into_gnl(const nlr_external_csr *ext,
                                gpu_neighbor_list_t *gnl,
                                const char *loop_name)
{
    /* zero-init everything; we touch only what we own */
    memset(gnl, 0, sizeof(*gnl));

    gnl->num_active  = ext->num_active;
    gnl->total_pairs = ext->total_pairs;

    const size_t off_bytes = (size_t)(ext->num_active + 1) * sizeof(int64_t);
    const size_t act_bytes = (size_t)(ext->num_active > 0 ? ext->num_active : 1)
                             * sizeof(int);
    const int64_t pairs = ext->total_pairs;
    const size_t nbr_bytes = (size_t)(pairs > 0 ? pairs : 1) * sizeof(int);

    /* offsets and d_active in SharedSpace (UVM) — host memcpy is fine */
    gnl->offsets  = (int64_t *) nlr_shared_alloc_bytes(off_bytes, "modea_csr_offsets");
    gnl->d_active = (int *)     nlr_shared_alloc_bytes(act_bytes, "modea_active");
    /* neighbors in DeviceSpace (GPU HBM) — must deep_copy from host view */
    gnl->neighbors = (int *) gizmo_gpu_alloc_device(nbr_bytes, "modea_csr_neighbors");
    /* Report the refusal rather than staging into null: the memcpy below writes
     * the caller's whole CSR, and the pair kernel reads all three. The caller
     * releases what did land and stops. */
    if(!gnl->offsets || !gnl->d_active || !gnl->neighbors) {
        nlr_stop_no_staging_memory(loop_name, "the injected neighbour list",
                                   off_bytes + act_bytes + nbr_bytes);
        return false;
    }
    memcpy(gnl->offsets,  ext->offsets,          off_bytes);
    memcpy(gnl->d_active, ext->active_indices,   act_bytes);

    if(pairs > 0) {
        Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            h_n(ext->neighbors, (size_t)pairs);
        Kokkos::View<int*, GIZMO_KOKKOS_DEVICE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            d_n(gnl->neighbors, (size_t)pairs);
        Kokkos::deep_copy(d_n, h_n);
    }
    return true;
}

static inline void
nlr_free_external_csr_gnl(gpu_neighbor_list_t *gnl)
{
    if(gnl->neighbors) Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(gnl->neighbors);
    if(gnl->d_active)  Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->d_active);
    if(gnl->offsets)   Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(gnl->offsets);
    gnl->neighbors = nullptr;
    gnl->d_active  = nullptr;
    gnl->offsets   = nullptr;
    gnl->num_active = 0;
    gnl->total_pairs = 0;
}

/* ============================================================================
 * run_mode_a<Spec> — generic Mode A path through the GPU NGL pipeline.
 *
 * Three-epoch staging contract (see neighbor_loop_runner.h doc):
 *   (1) host pre-arena: Spec::search_radius  → radii_uvm[num_active]
 *   (2) host pre-arena: Spec::populate_call_scalars → CallScalars cs
 *   (3) device post-NGL-build: Spec::load_active → d_actives[num_active]
 *
 * Both Kokkos launches go through gizmo_gpu_kernel_launch which wraps
 * parallel_for + fence + check_last_error (declarations/gpu_dispatch_templates.h).
 * No additional explicit fence is needed before host-side d_accums readback
 * — the launch helper already fenced.
 *
 * External CSR injection (args.external_csr != nullptr): the caller has
 * already built a symmetric gas CSR (e.g. corridor's gizmo_sym_*) and we
 * stage it into the gnl shape via nlr_stage_external_csr_into_gnl() instead
 * of calling gpu_ngb_list_build. SidxCacheKind::GasOnly only; other Specs
 * MUST leave args.external_csr null. Existing Specs unaffected.
 * ========================================================================== */

/* Per-active Mode-A staging can be CHUNKED so the PER-ACTIVE arrays
 * (d_actives + d_accums) stay bounded regardless of the rank-active count. This
 * cap bounds ONLY those two arrays: radii_uvm and the gnl CSR (offsets/neighbors)
 * are built ONCE over all N and stay full-N, so a loop whose CSR dominates its
 * per-active PODs is not helped by this cap (for the gradient loop the fat
 * accumulator dominates the CSR, which is why it is the useful customer). The
 * bound is an internal per-rank byte target (not a user knob, not a live
 * free-memory query), sized conservatively so ranks_per_node * cap stays under a
 * node's transient headroom for the target rank counts; it is NOT derived from
 * the live rank count, so a very dense packing may need a lower value -- a
 * rank-count-aware node budget is the eventual fix.
 *
 * SAFETY: the Mode-A arena ALIASES the host P/CellP (gpu_particles_arena_acquire
 * sets arena_P = P_host; there is no snapshot), so a per-chunk writeback IS
 * visible to later chunks' stage/pair reads. The unchunked path ran every pair
 * kernel BEFORE any writeback; chunking interleaves them. It is therefore
 * bitwise-safe ONLY when a Spec's apply_active_writeback writes no field that
 * any pair_kernel reads. That is a per-Spec AND per-config property, NOT implied
 * by "i-side" or "no ghost writeback" -- so it is an explicit opt-in trait
 * (Spec::mode_a_chunked_active_staging) the Spec author sets after auditing it.
 * Absent -> the Spec always stages the full set (K == N). */
static constexpr size_t NLR_MODE_A_STAGING_BYTES_CAP = (size_t)128 * 1024 * 1024;

/* Opt-in Spec trait gating the staging chunker above. Absent -> false. */
template <typename Spec, typename = void>
struct nlr_has_mode_a_chunked_active_staging : std::false_type {};
template <typename Spec>
struct nlr_has_mode_a_chunked_active_staging<Spec, decltype((void)Spec::mode_a_chunked_active_staging)>
    : std::true_type {};
template <typename Spec> static constexpr bool nlr_mode_a_chunked_active_staging_v() {
    if constexpr (nlr_has_mode_a_chunked_active_staging<Spec>::value) {
        return Spec::mode_a_chunked_active_staging;
    } else { return false; }
}

/* Defined below with the other lifecycle-trait helpers; forward-declared here as
 * a defensive backstop so a mis-set opt-in on a ghost-writeback Spec cannot chunk. */
template <typename Spec> static constexpr bool nlr_uses_ghost_writeback_v();

/* The staged i-side ActiveData array is written by the staging kernel and read
 * by the pair kernel, both on the device; no host code reads it. Keeping it in
 * device memory therefore costs nothing in reachability and keeps the staging
 * kernel's writes off the migratable shared path, which on a discrete-memory
 * device is where that kernel's cost lives.
 *
 * The accumulator array is NOT eligible and must stay shared: the host reads it
 * directly, in apply_active_writeback on the single-pass path and in the scatter
 * back into the driver's per-slot accumulators on the iterative path.
 *
 * On a host-only backend the two spaces are the same type, so this is a no-op
 * there by construction -- a local run can show that nothing regressed, it
 * cannot exercise the separation. Same non-throwing contract as above. */
static void *nlr_active_stage_alloc_bytes(size_t bytes, const char *label)
{
    if(bytes == 0) { return NULL; }
    return gizmo_gpu_alloc_device(bytes, label);
}

static void nlr_active_stage_free(void *p)
{
    if(p == NULL) { return; }
    Kokkos::kokkos_free<GIZMO_KOKKOS_DEVICE_SPACE>(p);
}

/* ============================================================================
 * TeamRowReduce support: combining lanes' partial accumulators.
 *
 * A Kokkos reducer over Spec::AccumData. Its join IS Spec::merge_accum and its
 * init IS Spec::zero_accum, so the lane combination and the cross-rank Mode-B
 * reply merge are the same algebra by construction rather than by agreement --
 * there is no second copy of the merge rules to drift.
 *
 * This rests on zero_accum being the identity element of merge_accum, which is
 * already required and already holds: Mode-B starts each peer's accumulator from
 * zero_accum before merging, so the pair is a monoid tree-wide. The two cases
 * that could have failed both check out -- DensitySpec::zero_accum sets its
 * sink fields to explicit MIN-reduction sentinels rather than zero, and
 * GradientsSpec byte-zeroes where 0 is a true identity because its Maxima and
 * Minima range over signed deltas, a set that contains the self-delta.
 *
 * join takes the whole accumulator. Decomposing it field-wise would be wrong,
 * not merely slower: SinkFeedSpec merges a coupled minimum-with-payload, taking
 * the peer's position only when the peer's potential is lower, and independent
 * per-field reduction would pair a winning value with a losing payload.
 * ========================================================================== */
template <typename Spec>
struct NlrAccumReducer {
    using reducer          = NlrAccumReducer<Spec>;
    using value_type       = typename Spec::AccumData;
    using result_view_type = Kokkos::View<value_type, Kokkos::AnonymousSpace,
                                          Kokkos::MemoryUnmanaged>;

    /* The lane combination moves an accumulator between threads with a shuffle,
     * and a shuffle of anything larger than one machine word is performed as a
     * whole number of INT-SIZED PIECES: the backends compute that count as
     * sizeof(T) / sizeof(int), truncating. An accumulator whose size is not a
     * whole number of those pieces therefore loses its trailing bytes in the
     * reduction -- silently, with no diagnostic, and only on the device, where
     * whichever field happens to sit last would come back holding whatever the
     * receiving lane had there before.
     *
     * That is a bad failure to leave discoverable only by a wrong answer, so it
     * is a build error instead. If a new accumulator trips this, pad it to a
     * multiple of sizeof(int) rather than reordering its fields to hide the
     * tail -- the next field added would put it back. An accumulator with no
     * members at all is exempt because it carries nothing to lose. */
    static_assert(std::is_empty<value_type>::value ||
                  (sizeof(value_type) % sizeof(int)) == 0,
                  "Spec::AccumData must be a whole number of int-sized pieces: the device lane "
                  "reduction shuffles it as sizeof(AccumData)/sizeof(int) words and would drop "
                  "the remainder. Pad the struct to a multiple of sizeof(int).");

    KOKKOS_INLINE_FUNCTION explicit NlrAccumReducer(value_type& v) : m_value(v) {}

    KOKKOS_INLINE_FUNCTION void join(value_type& dst, const value_type& src) const {
        Spec::merge_accum(dst, src);
    }
    KOKKOS_INLINE_FUNCTION void init(value_type& v) const { Spec::zero_accum(v); }
    KOKKOS_INLINE_FUNCTION value_type&       reference()         const { return m_value; }
    KOKKOS_INLINE_FUNCTION result_view_type  view()              const { return result_view_type(&m_value); }
    KOKKOS_INLINE_FUNCTION bool              references_scalar() const { return true; }

private:
    value_type& m_value;
};

/* The TeamRowReduce pair kernel, shared by both Mode-A dispatch sites.
 *
 * A named functor rather than a lambda for two reasons. It can be instantiated
 * once up front to size the team and again per chunk to run, which keeps the
 * occupancy query out of the chunk loop; and nvcc forbids defining an extended
 * device lambda inside another lambda, which a per-chunk lambda factory would
 * require.
 *
 * The two sites differ only in how a work item reaches its CSR row. The
 * single-pass site walks the active list directly, so row = chunk_base + i. The
 * iterative site walks a compacted active set into a build-time row index, so
 * row = csr_lookup[active_set[i]]. Passing active_set == nullptr selects the
 * former. Both index the staged actives and accumulators by i.
 */
template <typename Spec, typename DeviceCtx>
struct NlrModeATeamPairKernel {
    using ActiveData   = typename Spec::ActiveData;
    using AccumData    = typename Spec::AccumData;
    using ScatterData  = typename Spec::ScatterData;
    using NeighborData = typename Spec::NeighborData;
    using TeamMember   = typename Kokkos::TeamPolicy<>::member_type;

    DeviceCtx      ctx;
    ActiveData    *d_actives;
    AccumData     *d_accums;
    const int64_t *offsets;
    const int     *neighbors;
    const int     *active_set;   /* nullptr on the single-pass site */
    const int     *csr_lookup;   /* used only when active_set != nullptr */
    int            chunk_base;   /* used only when active_set == nullptr */
    /* One snapshot for the whole call, held by value so the device functor
     * carries it without reaching for a global. */
    typename Spec::CallScalars cs;
    struct GxMotionTargetSet   motion_targets;   /* for a loop that writes neighbour motion */

    KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const {
        const int i   = team.league_rank();
        const int row = (active_set != nullptr) ? csr_lookup[active_set[i]]
                                                : chunk_base + i;
        const ActiveData& a = d_actives[i];
        const int64_t start = offsets[row], end = offsets[row + 1];

        /* An empty or malformed row must produce the zero accumulator and
         * nothing else. The serial walk got this for free -- `for(nn = start;
         * nn < end; ...)` simply runs zero times when end <= start -- but a
         * team range takes a COUNT, and a non-positive one is not something to
         * hand it. Restoring the property explicitly keeps the two assignments
         * equivalent on degenerate rows as well as ordinary ones. */
        const int row_len = (end > start) ? (int)(end - start) : 0;
        if(row_len == 0) {
            Kokkos::single(Kokkos::PerTeam(team), [&]() { Spec::zero_accum(d_accums[i]); });
            return;
        }

        AccumData row_accum;
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, row_len),
            [&](int nn, AccumData& lane_accum) {
                ScatterData     s{};
                IdentitySidecar id{};
                if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {gx_motion_target_mark(motion_targets, neighbors[start + nn]);}
                NeighborData    nb = Spec::load_neighbor(ctx, neighbors[start + nn], id, a);
                Spec::pair_kernel(a, nb, lane_accum, s, cs);
            },
            NlrAccumReducer<Spec>(row_accum));

        /* Every lane leaves the reduction holding the combined value; one
         * publishes it. */
        Kokkos::single(Kokkos::PerTeam(team), [&]() { d_accums[i] = row_accum; });
    }
};

/* Lanes per row for a TeamRowReduce Spec, or 1 to select the flat kernel.
 *
 * Keyed only on structural properties -- execution space, the Spec's assignment
 * policy, its search_mode, and sizeof(AccumData). No caller name appears here,
 * so every ONEWAY loop gets the one-way width and every SYMMETRIC loop the
 * symmetric one with no per-loop work, which is the whole point of resolving it
 * in one place.
 *
 * The Kokkos bound is applied last and is a LEGALITY clamp: it reports the
 * largest team the backend can launch for this functor, and does not shrink as
 * the reduction value grows. The accumulator size is handled separately, and
 * deliberately bluntly, by the fat-accumulator cap. */
template <typename Spec, typename Functor>
static int nlr_mode_a_team_width(const Functor& f)
{
    if constexpr (gizmo_gpu_default_space_is_host()) {
        return 1;
    } else if constexpr (nlr_mode_a_pair_policy<Spec>() == ModeAPairAssignment::RowSerial) {
        return 1;
    } else {
        /* Every input is fixed by the Spec and the kernel type, so this resolves
         * once per (Spec, kernel) rather than once per dispatch: the occupancy
         * query behind it is not worth repeating, and the answer cannot change
         * within a run. */
        static const int resolved = [&]() {
            int target;
            if(NUMDIMS < 3) {
                target = NLR_TEAM_WIDTH_LOWDIM;
            } else {
                target = (Spec::search_mode == MODE_B_SEARCH_ONEWAY)
                         ? NLR_TEAM_WIDTH_ONEWAY : NLR_TEAM_WIDTH_SYMMETRIC;
                if(sizeof(typename Spec::AccumData) > NLR_TEAM_FAT_ACCUM_BYTES) {
                    target = NLR_TEAM_WIDTH_FAT_ACCUM;
                }
            }
            /* A non-positive answer means the backend reports NO launchable
             * team size for this functor, not "no limit" -- fall back to the
             * flat kernel rather than launching at the full target. */
            const int hw = gizmo_gpu_team_size_max(f);
            if(hw <= 0)       { return 1; }
            if(hw < target)   { target = hw; }
            int w = 1;
            while((w << 1) <= target) { w <<= 1; }
            return w;
        }();
        return resolved;
    }
}

template <typename Spec>
static void run_mode_a(const neighbor_loop_args& args, const double *radii)
{
    using ActiveData   = typename Spec::ActiveData;
    using AccumData    = typename Spec::AccumData;
    using ScatterData  = typename Spec::ScatterData;
    using NeighborData = typename Spec::NeighborData;
    using CallScalars  = typename Spec::CallScalars;
    using DeviceCtx    = typename Spec::DeviceContext;

    const int N = args.num_active;
    if(N <= 0) {
        /* Defensive: caller already early-outs in collective dispatch when
         * global_num_active == 0, but local-zero with global-positive is
         * legitimate (peer rank with no actives still must enter the
         * collective in Mode B; Mode A has no collective work, just no-op). */
        return;
    }

    /* (1) Host, pre-arena: stage caller-supplied radii into UVM for the
     * device-visible NGL build. Source `radii` is runner-staged on host
     * (call-lifetime only); we copy into shared/UVM so gpu_ngb_list_build
     * and the device pair kernel can read the per-active values. */
    const size_t radii_uvm_bytes = (size_t) N * sizeof(double);
    double *radii_uvm = (double *) nlr_shared_alloc_bytes(radii_uvm_bytes, "modea_radii");
    if(radii_uvm == NULL) {
        /* Nothing has been acquired yet, so there is nothing to release. */
        nlr_stop_no_staging_memory(Spec::loop_name, "the per-active search radii", radii_uvm_bytes);
        return;
    }
    for(int aa = 0; aa < N; aa++) {
        radii_uvm[aa] = radii[aa];
    }

    /* (2) Host, pre-arena: capture per-call scalar globals into a POD. */
    CallScalars cs = Spec::populate_call_scalars(args);

    /* Arena + freshness (matches sinks/sink_environment_gpu.cc:76,86). The
     * caller is responsible for the args.CellP=NULL-when-no-gas decision;
     * runner does not read All.TotN_gas. */
    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(args.num_total, args.P, args.CellP);
    struct particle_data *P_gpu    = gpu_particles_arena_P();
    struct gas_cell_data *CellP_gpu = (args.CellP != nullptr)
                                        ? gpu_particles_arena_CellP() : nullptr;

    /* NGL build using pre-arena radii. SIDX cache resolved from spec.
     * Hydro-corridor external-CSR path: when args.external_csr is non-null
     * the caller has already built a symmetric gas CSR — stage it into gnl
     * shape and skip the build. Contract: GasOnly Specs only. */
    gpu_neighbor_list_t gnl;
    gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                       Spec::loop_name);
    if(args.external_csr != nullptr) {
        /* Runtime checks — cannot be static_assert because that would fire
         * at template instantiation for every NotIterative Spec, including
         * non-GasOnly ones (sink_env1, dm_fuzzy, etc.) whose callers never
         * set external_csr. Compile-time enforcement is impossible since
         * external_csr is a runtime args field, not a Spec constexpr.
         *
         * All checks are UNCONDITIONAL (never diagnostic-gated):
         * external CSR injection is a sharp tool, contract violations
         * cause silent wrong-particle writeback (kernel stages for
         * external_csr->active_indices[aa] but writeback applies
         * d_accums[aa] to args.active_list[aa]). Detect loud always.
         *
         * On violation: soft bad-stop + release the arena and free radii_uvm
         * (both acquired/staged above) + return, skipping the corrupt-CSR
         * staging and the device walk. run_mode_a issues no MPI, so this
         * return cannot desync a peer; the caller's next phase poll drains.
         * else-if chain so a null pointer is never deref'd by a later check
         * (offsets[0] read only once ec->offsets is confirmed non-null). */
        const nlr_external_csr *ec = args.external_csr;
        const char *csr_err = nullptr;
        int         csr_code = 0;
        if(Spec::sidx_cache_kind != SidxCacheKind::GasOnly) { csr_err = "requires GasOnly cache"; csr_code = 7300; }
        else if(ec->num_active != N)                        { csr_err = "num_active mismatch"; csr_code = 7301; }
        else if(!ec->active_indices)                        { csr_err = "null active_indices"; csr_code = 7302; }
        else if(!ec->offsets)                               { csr_err = "null offsets"; csr_code = 7303; }
        else if(ec->total_pairs < 0)                        { csr_err = "negative total_pairs"; csr_code = 7304; }
        else if(ec->total_pairs > 0 && !ec->neighbors)      { csr_err = "null neighbors with total_pairs>0"; csr_code = 7305; }
        else if(N > 0 && ec->offsets[0] != 0)               { csr_err = "offsets[0]!=0"; csr_code = 7306; }
        else if(N > 0 && ec->offsets[N] != ec->total_pairs) { csr_err = "offsets[N]!=total_pairs"; csr_code = 7307; }
        else {
            /* Row order MUST match args.active_list elementwise — otherwise
             * the kernel accumulates for ec->active_indices[aa] but the host
             * writeback re-applies d_accums[aa] to args.active_list[aa]. */
            for(int aa = 0; aa < N; aa++) {
                if(ec->active_indices[aa] != args.active_list[aa]) { csr_err = "row order != active_list"; csr_code = 7308; break; }
                if(ec->offsets[aa+1] < ec->offsets[aa])            { csr_err = "non-monotonic offsets"; csr_code = 7309; break; }
            }
        }
        if(csr_err != nullptr) {
            fprintf(stderr, "[neighbor_loop_runner rank=%d caller=%s] external-CSR contract violation: %s\n",
                    ThisTask, Spec::loop_name, csr_err);
            fflush(stderr);
            endrun(csr_code);
            gpu_particles_arena_release();
            Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm);
            return;
        }
        if(!nlr_stage_external_csr_into_gnl(ec, &gnl, Spec::loop_name)) {
            nlr_free_external_csr_gnl(&gnl);
            Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm);
            gpu_particles_arena_release();
            return;
        }
    } else {
        /* Active-source-in-pool contract: stage explicit P[active_i].Pos for specs
         * whose active sources may be non-pool (else nullptr keeps the compact
         * fast-path). See neighbor_loop_runner.h. Radii are already explicit. */
        std::vector<double> _nlr_srcpos_storage;
        const double* _nlr_srcpos = nlr_stage_explicit_source_positions<Spec>(
            args.P, args.active_list, N, _nlr_srcpos_storage);
        gpu_ngb_list_build(P_gpu, args.num_total,
                           args.active_list, N,
                           Spec::search_mode,
                           (int)nlr_effective_neighbor_type_mask(args, Spec::neighbor_type_mask),
                           &gnl, sidx,
                           1.0, radii_uvm, _nlr_srcpos, Spec::loop_name,
                           nlr_spec_symmetric_j_radius_scale<Spec>(),
                           Spec::radius_policy);
    }

    /* A build that ran out of memory hands back the empty list, without the row
     * index the kernel reads first, and has already asked for the stop. Release
     * this call's own buffers and return; no MPI has been issued here. */
    if(gnl.d_active == nullptr) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm);
        if(args.external_csr != nullptr) { nlr_free_external_csr_gnl(&gnl); }
        else { gpu_ngb_list_free(&gnl, sidx); }
        gpu_particles_arena_release();
        return;
    }

    /* Chunk the fat per-active staging so its transient footprint is bounded.
     * K stays at the full N unless (a) the Spec opts into chunked staging (its
     * writeback is disjoint from every pair read -- see the trait contract) AND
     * (b) staging all N rows would exceed the internal per-rank byte cap. The
     * NoScatter / !uses_ghost_writeback conjuncts are defensive backstops (every
     * current Spec is NoScatter; the opt-in trait is the real contract). K
     * degrades to 1 if a single record exceeds the cap. K == N reproduces the
     * unchunked path exactly. */
    const bool chunk_ok =
        nlr_mode_a_chunked_active_staging_v<Spec>()
        && std::is_same<typename Spec::ScatterData, NoScatter>::value
        && !nlr_uses_ghost_writeback_v<Spec>();
    const size_t rec_bytes = sizeof(ActiveData) + sizeof(AccumData);
    int K = N;
    if(chunk_ok && (size_t)N * rec_bytes > NLR_MODE_A_STAGING_BYTES_CAP) {
        size_t k = NLR_MODE_A_STAGING_BYTES_CAP / rec_bytes;
        if(k < 1) { k = 1; }
        if(k < (size_t)N) { K = (int)k; }
    }

    /* Allocate chunk-sized ActiveData[] and AccumData[] arrays, both through
     * non-throwing allocators. They live in DIFFERENT spaces: the staged
     * actives are device-only (see nlr_active_stage_alloc_bytes), the
     * accumulators must stay host-readable for the writeback below. These are
     * the largest per-active transients (the demonstrated FIF OOM site); an
     * allocation failure here means the corresponding pool is genuinely full
     * (K is byte-capped), so controlled-stop with the buffer named in the
     * ledger rather than a hard terminate. run_mode_a issues no MPI, so the
     * request drains collectively at the caller's next phase poll (same as the
     * external-CSR contract-violation path above). */
    ActiveData *d_actives = (ActiveData *) nlr_active_stage_alloc_bytes((size_t)K * sizeof(ActiveData), "modea_active_data");
    AccumData  *d_accums  = (AccumData  *) nlr_shared_alloc_bytes((size_t)K * sizeof(AccumData),  "modea_accum_data");
    if(d_actives == NULL || d_accums == NULL) {
        if(d_accums)  { Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_accums); }
        nlr_active_stage_free(d_actives);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm);
        if(args.external_csr != nullptr) { nlr_free_external_csr_gnl(&gnl); }
        else { gpu_ngb_list_free(&gnl, sidx); }
        gpu_particles_arena_release();
        gizmo_request_controlled_stop(7710,
            "run_mode_a: Mode-A per-active staging out of memory (modea_active_data is device-resident, "
            "modea_accum_data is host-visible) -- reduce the active set; note that adding ranks per node "
            "does NOT relieve device memory, since the ranks on a node share it",
            __FILE__, __LINE__, __FUNCTION__);
        return;
    }

    /* Build DeviceContext. Specs that extend Spec::DeviceContext beyond
     * NeighborLoopDeviceContextBase get populate_device_context invoked
     * here; base-only Specs skip the call (trait check). */
    DeviceCtx ctx;
    ctx.P         = P_gpu;
    ctx.CellP     = CellP_gpu;
    ctx.num_total = args.num_total;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; the runner captures it by value into Kokkos device lambdas");
    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    NlrDeviceContextCleanupGuard<Spec> _nlr_dctx_cleanup_guard(args, ctx);
    /* A hook that could not stage its buffers has already asked for the stop.
     * Same exit as the staging failure above: release everything this call
     * obtained and return — this function issues no MPI, so the request drains
     * at the caller's next phase poll. The guard releases the hook's buffers. */
    if(ctx.populate_failed) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_accums);
        nlr_active_stage_free(d_actives);   /* device space, not shared -- see nlr_active_stage_alloc_bytes */
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm);
        if(args.external_csr != nullptr) { nlr_free_external_csr_gnl(&gnl); }
        else { gpu_ngb_list_free(&gnl, sidx); }
        gpu_particles_arena_release();
        return;
    }

    /* (3)-(4) Chunked stage -> pair-kernel -> writeback. Each chunk [c0, c0+n)
     * stages into chunk-local slots [0,n): CSR rows / radii are read at the
     * absolute active index aa = c0 + kk, results land in d_actives[kk] /
     * d_accums[kk], and the writeback re-applies them to args.active_list[aa]
     * before the next chunk reuses the arrays. d_active/offsets/neighbors are
     * gnl-resident (built once over all N). K == N is a single pass identical
     * to the unchunked path. */
    int     *d_active_idx = gnl.d_active;
    int64_t *offsets      = gnl.offsets;
    int     *neighbors    = gnl.neighbors;

    /* Team width is a property of the Spec and the kernel type, not of a chunk,
     * so it is resolved once here rather than per chunk -- the occupancy query
     * behind it is not something to repeat inside the loop. */
    using TeamKernel = NlrModeATeamPairKernel<Spec, DeviceCtx>;
    int team_width = 1;
    if constexpr (nlr_mode_a_pair_policy<Spec>() == ModeAPairAssignment::TeamRowReduce) {
        TeamKernel probe{ctx, d_actives, d_accums, offsets, neighbors, nullptr, nullptr, 0, cs};
        team_width = nlr_mode_a_team_width<Spec>(probe);
    }

    for(int c0 = 0; c0 < N; c0 += K) {
        const int n = (N - c0 < K) ? (N - c0) : K;

        /* stage ActiveData for [c0, c0+n) into chunk-local [0,n). */
        gizmo_gpu_kernel_launch("nlr_stage_active", n, KOKKOS_LAMBDA(int kk) {
            const int aa = c0 + kk;
            d_actives[kk] = Spec::load_active(ctx, aa, d_active_idx[aa],
                                              radii_uvm[aa], cs);
        });

        /* pair-kernel over [c0, c0+n) — generic over Spec.
         *
         * Two assignments of the same physics. The flat form gives one work item
         * per active particle, which then walks its whole CSR row in sequence.
         * The team form gives one team per active particle whose lanes stride
         * the row together and combine through Spec::merge_accum. Which one runs
         * is decided by nlr_mode_a_team_width from structural properties only;
         * width 1 selects the flat form, and does so through if constexpr, so a
         * width-1 Spec compiles to exactly the kernel it compiled to before
         * teams existed rather than to a one-lane imitation of a team. */
        {
            const double t_pair_kernel_start = my_second();
            const struct GxMotionTargetSet motion_targets = gx_motion_target_view();

            auto flat_kernel = KOKKOS_LAMBDA(int kk) {
                /* Named here, unconditionally, so that the capture happens outside the
                 * `if constexpr` below: a device lambda may not first-capture a variable
                 * inside one, and a Spec that does not write neighbour motion would
                 * otherwise reach the capture only through the discarded branch. */
                const struct GxMotionTargetSet &targets = motion_targets;
                const int aa = c0 + kk;
                Spec::zero_accum(d_accums[kk]);
                const ActiveData& a = d_actives[kk];
                ScatterData s{};                     /* NoScatter for ActiveReduceOnly */
                int64_t start = offsets[aa], end = offsets[aa + 1];
                for(int64_t nn = start; nn < end; nn++) {
                    int j = neighbors[nn];
                    if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {gx_motion_target_mark(targets, j);}
                    IdentitySidecar id{};            /* NoIdentity */
                    NeighborData nb = Spec::load_neighbor(ctx, j, id, a);
                    Spec::pair_kernel(a, nb, d_accums[kk], s, cs);
                }
            };

            if constexpr (nlr_mode_a_pair_policy<Spec>() == ModeAPairAssignment::TeamRowReduce) {
                if(team_width > 1) {
                    TeamKernel fn{ctx, d_actives, d_accums, offsets, neighbors, nullptr, nullptr, c0, cs, motion_targets};
                    gizmo_gpu_team_kernel_launch(Spec::loop_name, n, team_width, fn);
                } else {
                    gizmo_gpu_kernel_launch(Spec::loop_name, n, flat_kernel);
                }
            } else {
                gizmo_gpu_kernel_launch(Spec::loop_name, n, flat_kernel);
            }
            cpu_charge_child(CPU_PAIR_KERNEL, timediff(t_pair_kernel_start, my_second()));
        }
        /* Launches fenced internally by gizmo_gpu_kernel_launch. UVM coherent ->
         * host reads d_accums[0,n) directly. */

        /* (4) host writeback for [c0, c0+n) before the next chunk reuses arrays. */
        {
            for(int kk = 0; kk < n; kk++) {
                Spec::apply_active_writeback(args, c0 + kk, args.active_list[c0 + kk], d_accums[kk]);
            }
        }
    }


    /* Cleanup. SIDX cache pointer passed so the free leaves cached storage
     * intact for sink_feed/sink_swk reuse (matches existing
     * sink_environment_gpu.cc:261 idiom). External-CSR path frees only what
     * we staged (gnl offsets/neighbors/d_active); caller owns host CSR. */
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_accums);
    nlr_active_stage_free(d_actives);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm);
    if(args.external_csr != nullptr) {
        nlr_free_external_csr_gnl(&gnl);
    } else {
        gpu_ngb_list_free(&gnl, sidx);
    }
    gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);
}

/* ============================================================================
 * SFINAE: detect optional Spec::modeb_threshold_sum / _max constexpr.
 * If absent, fall back to caller-supplied default.
 * ========================================================================== */

template <typename Spec, typename = void>
struct nlr_has_threshold_sum : std::false_type {};
template <typename Spec>
struct nlr_has_threshold_sum<Spec, decltype((void)Spec::modeb_threshold_sum)>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_threshold_max : std::false_type {};
template <typename Spec>
struct nlr_has_threshold_max<Spec, decltype((void)Spec::modeb_threshold_max)>
    : std::true_type {};

template <typename Spec>
static int nlr_spec_threshold_sum(int fallback) {
    if constexpr (nlr_has_threshold_sum<Spec>::value) {
        return (int)Spec::modeb_threshold_sum;
    } else {
        return fallback;
    }
}
template <typename Spec>
static int nlr_spec_threshold_max(int fallback) {
    if constexpr (nlr_has_threshold_max<Spec>::value) {
        return (int)Spec::modeb_threshold_max;
    } else {
        return fallback;
    }
}

/* ============================================================================
 * SFINAE: detect optional Spec::uses_* lifecycle traits (default false) and
 * dispatch to the corresponding hook methods (default no-op). Two channels:
 *
 *   detector  : audit/debug — catches illegal kernel writes to imported
 *               ghosts. Trait `uses_ghost_write_detector`.
 *   writeback : physics state propagation — pre-kernel ghost-state snapshot
 *               + post-kernel reverse-comm of any j-side writes. Trait
 *               `uses_ghost_writeback`. Spec body is the per-flag #ifdef
 *               union of all enabled physics-flag writebacks.
 *
 * Each hook is gated on (a) the trait being true AND (b) whether the chosen
 * path imports ghosts (per-Spec audit decided imported-ghost-only for these
 * hooks; future Specs may differ).
 * ========================================================================== */

/* Trait detection: Spec::uses_ghost_write_detector / _ghost_writeback.
 * Absent ⇒ false. */
template <typename Spec, typename = void>
struct nlr_has_uses_ghost_write_detector : std::false_type {};
template <typename Spec>
struct nlr_has_uses_ghost_write_detector<Spec, decltype((void)Spec::uses_ghost_write_detector)>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_uses_ghost_writeback : std::false_type {};
template <typename Spec>
struct nlr_has_uses_ghost_writeback<Spec, decltype((void)Spec::uses_ghost_writeback)>
    : std::true_type {};

template <typename Spec> static constexpr bool nlr_uses_ghost_write_detector_v() {
    if constexpr (nlr_has_uses_ghost_write_detector<Spec>::value) {
        return Spec::uses_ghost_write_detector;
    } else { return false; }
}
template <typename Spec> static constexpr bool nlr_uses_ghost_writeback_v() {
    if constexpr (nlr_has_uses_ghost_writeback<Spec>::value) {
        return Spec::uses_ghost_writeback;
    } else { return false; }
}

/* Hook-method detection: Spec::ghost_write_detector_begin etc. Absent ⇒
 * runner skips the call (compile-time short-circuit; no runtime cost). */
template <typename Spec, typename = void>
struct nlr_has_hook_gwd_begin : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwd_begin<Spec,
    decltype(Spec::ghost_write_detector_begin(std::declval<const neighbor_loop_args&>(),
                                              std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_hook_gwd_end : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwd_end<Spec,
    decltype(Spec::ghost_write_detector_end(std::declval<const neighbor_loop_args&>(),
                                             std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_hook_gwb_begin : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwb_begin<Spec,
    decltype(Spec::ghost_writeback_begin(std::declval<const neighbor_loop_args&>(),
                                         std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

template <typename Spec, typename = void>
struct nlr_has_hook_gwb_end : std::false_type {};
template <typename Spec>
struct nlr_has_hook_gwb_end<Spec,
    decltype(Spec::ghost_writeback_end(std::declval<const neighbor_loop_args&>(),
                                       std::declval<const NeighborLoopPlan&>()))>
    : std::true_type {};

/* Optional label override: Spec::ghost_write_detector_name. Absent ⇒ runner
 * default labels the detector with Spec::loop_name. Used for the two
 * sink_env Specs whose detector labels predate the runner-template loop_name
 * convention and need to be preserved. */
template <typename Spec, typename = void>
struct nlr_has_ghost_write_detector_name : std::false_type {};
template <typename Spec>
struct nlr_has_ghost_write_detector_name<Spec, decltype((void)Spec::ghost_write_detector_name)>
    : std::true_type {};

template <typename Spec>
static constexpr const char *nlr_ghost_write_detector_name()
{
    if constexpr (nlr_has_ghost_write_detector_name<Spec>::value) {
        return Spec::ghost_write_detector_name;
    } else {
        return Spec::loop_name;
    }
}

/* Dispatch wrappers. Each gates on:
 *   (a) `uses_*` trait true (Spec opted in)
 *   (b) `nlr_path_uses_imported_ghosts(plan.path)` — for SinkEnv1Spec these
 *        hooks are imported-ghost-only. The path gate
 *        IS the policy; the hook trait is the Spec opt-in.
 *   (c) hook method exists (SFINAE) — present ⇒ Spec custom hook fires;
 *        absent ⇒ runner default fires (::ghost_write_detector_begin(name)
 *        / ::ghost_write_detector_end()). The dispatcher enforces that
 *        begin/end hook presence is symmetric: defining begin without end
 *        (or vice versa) is a compile error, not a half-default. */

template <typename Spec>
static void nlr_dispatch_ghost_write_detector_begin(const neighbor_loop_args& args,
                                                    const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_write_detector_v<Spec>()) {
        static_assert(nlr_has_hook_gwd_begin<Spec>::value == nlr_has_hook_gwd_end<Spec>::value,
                      "Spec must define both ghost_write_detector_begin/end or neither");
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwd_begin<Spec>::value) {
                Spec::ghost_write_detector_begin(args, plan);
            } else {
                ::ghost_write_detector_begin(nlr_ghost_write_detector_name<Spec>());
            }
        }
    }
}
template <typename Spec>
static void nlr_dispatch_ghost_write_detector_end(const neighbor_loop_args& args,
                                                  const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_write_detector_v<Spec>()) {
        static_assert(nlr_has_hook_gwd_begin<Spec>::value == nlr_has_hook_gwd_end<Spec>::value,
                      "Spec must define both ghost_write_detector_begin/end or neither");
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwd_end<Spec>::value) {
                Spec::ghost_write_detector_end(args, plan);
            } else {
                ::ghost_write_detector_end();
            }
        }
    }
}
template <typename Spec>
static void nlr_dispatch_ghost_writeback_begin(const neighbor_loop_args& args,
                                               const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_writeback_v<Spec>()) {
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwb_begin<Spec>::value) {
                Spec::ghost_writeback_begin(args, plan);
            }
        }
    }
}
template <typename Spec>
static void nlr_dispatch_ghost_writeback_end(const neighbor_loop_args& args,
                                             const NeighborLoopPlan& plan)
{
    if constexpr (nlr_uses_ghost_writeback_v<Spec>()) {
        if(nlr_path_uses_imported_ghosts(plan.path)) {
            if constexpr (nlr_has_hook_gwb_end<Spec>::value) {
                gx_motion_target_set_armed(nlr_spec_writes_neighbour_motion_v<Spec> ? 1 : 0);
                Spec::ghost_writeback_end(args, plan);
                gx_motion_target_set_armed(0);
            }
        }
    }
}
/* ============================================================================
 * Public entry: run_neighbor_loop<Spec>
 * ========================================================================== */

template <typename Spec>
void run_neighbor_loop(const neighbor_loop_args& args_in)
{
    /* The one live argument view for this call.  A ghost import can raise the particle
     * capacity and move P[]/CellP[], so the cached base pointers are refreshed in this
     * object below; every later reader sees the current values.  There is deliberately no
     * second, unrefreshed copy in scope -- reading a stale P[] reads freed storage. */
    neighbor_loop_args args = args_in;

    /* ---- Compile-time spec consistency ---- */

    /* WritePattern ↔ AccumData/ScatterData consistency. */
    static_assert(Spec::write_pattern != WritePattern::ActiveReduceOnly ||
                  std::is_same<typename Spec::ScatterData, NoScatter>::value,
        "ActiveReduceOnly requires ScatterData == NoScatter");
    static_assert(Spec::write_pattern != WritePattern::NeighborScatter ||
                  std::is_same<typename Spec::AccumData, NoAccum>::value,
        "NeighborScatter requires AccumData == NoAccum");

    /* Active-source-in-pool contract (see neighbor_loop_runner.h). */
    static_assert(nlr_spec_satisfies_source_pool_contract_v<Spec>,
        "Cached-SIDX Spec must declare 'static constexpr bool mode_a_active_sources_in_sidx_pool' "
        "(true = active sources are SIDX-pool members; false = runner stages explicit P[].Pos). "
        "Prevents the stale gas-only-compact source-position bug for non-pool actives.");

    /* POD/device-copy contract — captured by value into Kokkos lambdas
     * and/or staged into UVM arrays. Trivially-copyable is the binding
     * requirement. */
    static_assert(std::is_trivially_copyable<typename Spec::CallScalars>::value,
        "Spec::CallScalars must be trivially-copyable (lambda capture by value)");
    static_assert(std::is_trivially_copyable<typename Spec::ActiveData>::value,
        "Spec::ActiveData must be trivially-copyable (UVM-staged)");
    static_assert(std::is_trivially_copyable<typename Spec::NeighborData>::value,
        "Spec::NeighborData must be trivially-copyable (built per-pair on device)");
    static_assert(std::is_trivially_copyable<typename Spec::AccumData>::value,
        "Spec::AccumData must be trivially-copyable (UVM-staged)");

    /* ---- Dispatch ----
     *
     * Selection precedence (highest first):
     *   1. args.dispatch_override = A -> Mode A unconditionally.
     *   2. args.dispatch_override = B -> Mode B (local if NTask==1, remote else).
     *   3. Threshold dispatch: if (sum_active>0 && sum_active<=TS &&
     *      max_active<=TM) → Mode B; else Mode A. Hierarchy of TS/TM:
     *      parameterfile NeighborLoopModeBThreshold{Sum,Max} >
     *      Spec::modeb_threshold_{sum,max} constexpr defaults (64/64 today).
     *      Setting the parameterfile pair above every active count selects
     *      Mode B for the whole run; setting it to 0 selects Mode A.
     *
     * Note (active-epoch caveat): Mode B host-frozen actives[] are
     * NOT bit-equivalent to Mode A's device-staged post-neighbor-list-build
     * actives. The two modes therefore need not agree bit-for-bit on which
     * particles a call treats as active; consistency between them is a
     * property of the dispatch policy, not of this helper.
     */
    /* Dispatch priority: args.dispatch_override > adaptive threshold.
     * The args field is the corridor mode-decision hook (hydro_corridor.cc): when
     * a corridor consumer sets this to force coherent Mode A or Mode B across the
     * whole hydro corridor (cellcorrections/gradients/hydro_force), the per-call
     * override wins so corridor coherence is enforced, not advisory. */
    const NlrForceMode force_mode = args.dispatch_override;
    const bool force_a   = (force_mode == NlrForceMode::A);
    const bool force_b   = (force_mode == NlrForceMode::B);

    /* Threshold dispatch. Allreduce sum + max of args.num_active.
     * Skipped when the caller supplied a dispatch override (cheap path), which
     * is why the global count is not available on every path: it is computed
     * only where something needs it — the threshold decision itself, and the
     * forced-Mode-B size guard. It is never computed merely to report it. */
    bool select_mode_b = force_b;
    int global_num_active = -1;       /* -1 = not computed on this path */
    if(!force_a && !force_b) {
        int local_act = args.num_active;
        int sum_act = 0, max_act = 0;
        MPI_Allreduce(&local_act, &sum_act, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_act, &max_act, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        global_num_active = sum_act;
        /* Spec::modeb_threshold_{sum,max} via SFINAE; default 64/64. */
        const int spec_default_sum = nlr_spec_threshold_sum<Spec>(64);
        const int spec_default_max = nlr_spec_threshold_max<Spec>(64);
        const int TS = gizmo_nlr_modeb_threshold_sum_for(Spec::loop_name, spec_default_sum);
        const int TM = gizmo_nlr_modeb_threshold_max_for(Spec::loop_name, spec_default_max);
        select_mode_b = (sum_act > 0) && (sum_act <= TS) && (max_act <= TM);
    } else if(force_b) {
        /* Forced Mode B skipped the dispatch Allreduce, but its size guard
         * needs the global count, so it is done here. */
        int local_act = args.num_active;
        int sum_act = 0;
        MPI_Allreduce(&local_act, &sum_act, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        global_num_active = sum_act;
    }
    /* Globally-zero-active call: do NO neighbor work. global_num_active is the
     * dispatch Allreduce of active particles (set on the threshold + force-B
     * paths), so it is identical on every rank -> this return is collective-
     * symmetric (all ranks return together, skipping the Mode-A ghost import /
     * writeback / cleanup as a matched set). NOT the banned local num_active==0
     * early return: the condition is GLOBAL. Without it, a zero-active call
     * falls to Mode A and fires a spurious ghost import with nothing to compute.
     * The forced-Mode-A path does not compute the count, so it does not reach
     * this and still pays that import; giving it the saving means giving it a
     * collective it does not otherwise need, which is a separate decision. */
    if(global_num_active == 0) {
        return;
    }

    /* Compute the execution plan from the dispatch decision. Path is the
     * single-source-of-truth; predicates derive from it. */
    NeighborLoopPlan plan;
    if(force_a) {
        plan.path = NeighborLoopPlan::Path::ModeA_GpuNgl;
    } else if(force_b || select_mode_b) {
        plan.path = (NTask > 1) ? NeighborLoopPlan::Path::ModeB_Remote
                                : NeighborLoopPlan::Path::ModeB_Local;
    } else {
        plan.path = NeighborLoopPlan::Path::ModeA_GpuNgl;
    }
    plan.num_active_global = global_num_active;   /* -1 on dispatch-override paths */

    GxDeviceTreeView mode_d_tree{};
    /* A loop whose kernel reads neighbour state it is itself changing needs the
     * live particle, which a ghost copy is not (nlr_spec_needs_live_neighbours).
     * It is answered where its neighbours live: on the device below when the
     * tree can be described there, by the host walker otherwise. */
    if constexpr (nlr_spec_needs_live_neighbours_v<Spec>) {
        if(plan.path == NeighborLoopPlan::Path::ModeA_GpuNgl) {
            plan.path = (NTask > 1) ? NeighborLoopPlan::Path::ModeB_Remote
                                    : NeighborLoopPlan::Path::ModeB_Local;
        }
    }
#ifdef NEIGHBOR_LOOP_MODE_D
    /* Mode D takes the calls Mode A would have taken, for every loop it can
     * serve, when every rank can describe its tree to the device -- the same
     * decision the iterative driver makes, and for the same reasons (stated
     * there).  A caller-owned neighbour list is Mode A by construction, and an
     * explicit override is honoured as given.  The vote is collective; the
     * label was chosen from a collective too, so every rank reaches it. */
    if constexpr (nlr_spec_mode_d_eligible_v<Spec>) {
        if((plan.path == NeighborLoopPlan::Path::ModeA_GpuNgl ||
            (nlr_spec_needs_live_neighbours_v<Spec> && !select_mode_b)) && !force_a &&
           args.external_csr == nullptr) {
            const int ready_local = (gx_device_fused_walk_prepare(&mode_d_tree, Spec::loop_name) == 0) ? 1 : 0;
            int ready = 0;
            MPI_Allreduce(&ready_local, &ready, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            if(ready) {plan.path = NeighborLoopPlan::Path::ModeD_DeviceFused;}
        }
    }
#endif


    /* ---- Stage radii once ---- */
    /* Computed via Spec::search_radius. Used for any path-conditional
     * prep (Mode A) AND the chosen walker. Pointer is call-lifetime only;
     * see contract in mesh/neighbor_loop_runner.h. */
    std::vector<double> radii(args.num_active);
    for(int aa = 0; aa < args.num_active; ++aa) {
        radii[aa] = Spec::search_radius(args, aa, args.active_list[aa]);
    }
    /* Neighbours whose motion this loop changes are raised when the call ends. */
    NlrMotionTargetScope<Spec> motion_target_scope;

    /* ---- Hard-corridor counter snapshot (always-on, every build) ---- */
    /* Mode B paths must NOT enter move_particles, ghost_exchange_impl, or
     * gpu_particles_arena_acquire, and must NOT mutate NumPart. Counters
     * live in declarations/lifecycle_counters.h, incremented at API entry
     * by the owner TUs. */
    const uint64_t s_drift0 = g_global_drift_counter;
    const uint64_t s_ghost0 = g_ghost_import_counter;
    const uint64_t s_arena0 = g_gpu_arena_acquire_counter;
    const int      s_np0    = NumPart;


    /* ---- Path-conditional prep on imported-ghost paths ---- */
    /* Mode A's substrate today satisfies its freshness + neighbor-pool
     * requirements via gizmo_request_filtered_ghost_import_fresh (full
     * global drift + ghost import). Mode B paths skip this entirely;
     * peer-local pool + lazy candidate drift is sufficient. The dt_prep_import
     * timer is 0 for Mode B paths (genuine 0 — the API isn't called).
     *
     * Caller-owned ghost pool (external_csr non-null): the caller imported
     * the pool the CSR indexes and owns its lifetime — skip the import (a
     * fresh import could renumber ghost slots under the CSR) and, below, the
     * cleanup. Contract enforced loudly here; see neighbor_loop_runner.h. */
    if(args.external_csr != nullptr && !nlr_path_uses_imported_ghosts(plan.path)) {
        /* external_csr is a Mode-A-only input; a Mode B dispatch with a CSR
         * supplied means the caller's dispatch_override and CSR provisioning
         * disagree — fail loudly rather than silently ignoring the CSR. */
        fprintf(stderr, "[neighbor_loop_runner rank=%d caller=%s] external_csr supplied but "
                "dispatch selected a Mode B path — caller contract violation.\n",
                ThisTask, Spec::loop_name);
        fflush(stderr);
        endrun(7312);
    }
    if(nlr_path_uses_imported_ghosts(plan.path)) {
        if(nlr_caller_owns_ghost_pool(args)) {
            /* At NTask>1 the CSR references imported ghost slots: require the
             * caller's pool + provenance to be LIVE, and args to have been
             * built AFTER the caller's import (num_total spans the ghosts).
             * At NTask==1 no pool exists (imports early-out) — nothing to
             * verify. ghost_get_num_ghosts()==0 is NOT a liveness signal
             * (also 0 for a live zero-ghost pool); ghost_pool_is_live() is. */
            if(NTask > 1) {
                const char *own_err = nullptr;
                if(!ghost_pool_is_live())                 own_err = "ghost pool not live";
                else if(ghost_get_home_rank()  == nullptr ||
                        ghost_get_home_index() == nullptr) own_err = "ghost provenance maps absent";
                else if(args.num_total != NumPart)         own_err = "args.num_total != NumPart (args built before caller's import?)";
                if(own_err != nullptr) {
                    fprintf(stderr, "[neighbor_loop_runner rank=%d caller=%s] caller-owned ghost-pool "
                            "contract violation: %s\n", ThisTask, Spec::loop_name, own_err);
                    fflush(stderr);
                    endrun(7313);
                }
            }
        } else {
            gizmo_request_filtered_ghost_import_fresh(Spec::loop_name,
                                                       Spec::search_mode,
                                                       nlr_effective_neighbor_type_mask(args, Spec::neighbor_type_mask),
                                                       args.active_list,
                                                       args.num_active,
                                                       radii.data(),
                                                       args.ghost_safety_factor,
                                                       Spec::radius_policy,
                                                       nlr_spec_symmetric_j_radius_scale<Spec>());
            /* Ghost import grew NumPart and may have realloc'd P/CellP. Refresh
             * the runner's data view; only paths that imported ghosts read this
             * extended view (Mode B paths use the original args via copy). */
            args.num_total = NumPart;
            args.P         = P;
            args.CellP     = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;
        }
    }

    /* ---- Spec lifecycle hooks (begin) ---- */
    /* Mode A imported-ghost ordering invariant (two channels: detector +
     * writeback). The Spec's writeback_begin/_end body is the per-flag
     * #ifdef union of all enabled physics-flag writebacks for that loop:
     *
     *   request_filtered_ghost_import_fresh  (above)
     *   ghost_write_detector_begin           (this hook)
     *   ghost_writeback_begin                (this hook; per-flag union)
     *   <run_mode_a kernel>
     *   ghost_writeback_end                  (reverse order, below)
     *   ghost_write_detector_end
     *   ghost_exchange_cleanup               (below)
     *
     * Each dispatch helper checks both the Spec's `uses_*` trait AND the
     * path-imports-ghosts predicate; on Mode B paths the predicate is
     * false and all four hook calls compile to no-ops. */
    nlr_dispatch_ghost_write_detector_begin<Spec>(args, plan);
    nlr_dispatch_ghost_writeback_begin<Spec>(args, plan);

    /* ---- Path dispatch ---- */
    switch(plan.path) {
        case NeighborLoopPlan::Path::ModeA_GpuNgl:
            run_mode_a<Spec>(args, radii.data());
            break;
        case NeighborLoopPlan::Path::ModeB_Local:
            run_mode_b_local<Spec>(args, radii.data());
            break;
        case NeighborLoopPlan::Path::ModeB_Remote:
            run_mode_b_remote<Spec>(args, radii.data());
            break;
        case NeighborLoopPlan::Path::ModeD_DeviceFused:
#ifdef NEIGHBOR_LOOP_MODE_D
            if constexpr (nlr_spec_mode_d_eligible_v<Spec>) {
                run_mode_d<Spec>(args, radii.data(), mode_d_tree);
                break;
            }
#endif
            /* Never chosen without the build flag and the eligibility above. */
            endrun(90001032);
            break;
    }

    /* ---- Spec lifecycle hooks (end, reverse order) ---- */
    nlr_dispatch_ghost_writeback_end<Spec>(args, plan);
    nlr_dispatch_ghost_write_detector_end<Spec>(args, plan);

    /* ---- Imported-ghost cleanup ---- */
    /* Caller-owned pools (external_csr) outlive this call — the caller tears
     * them down at the end of its span; skip cleanup here. See
     * neighbor_loop_runner.h. */
    if(nlr_path_uses_imported_ghosts(plan.path) && NTask > 1) {
        if(!nlr_caller_owns_ghost_pool(args)) {
            ghost_exchange_cleanup();
        }
    }

    /* ---- Hard-corridor enforcement (Mode B paths) ---- */
    /* HARD ABORT on any counter advance or NumPart change across the path
     * body. Always-on; cheap (4 uint64 compares + one int compare). */
    if(plan.path == NeighborLoopPlan::Path::ModeB_Local ||
       plan.path == NeighborLoopPlan::Path::ModeB_Remote) {
        const bool drift_violation = (g_global_drift_counter      != s_drift0);
        const bool ghost_violation = (g_ghost_import_counter      != s_ghost0);
        const bool arena_violation = (g_gpu_arena_acquire_counter != s_arena0);
        const bool np_violation    = (NumPart != s_np0);
        if(drift_violation || ghost_violation || arena_violation || np_violation) {
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                    "[NLR CORRIDOR ABORT rank=%d caller=%s path=%s] Mode B path "
                    "violated tiny-N corridor invariant during run_neighbor_loop. "
                    "Counter deltas: drift=%llu ghost=%llu arena=%llu NumPart_pre=%d NumPart_post=%d\n",
                    rank, Spec::loop_name, nlr_path_label(plan.path),
                    (unsigned long long)(g_global_drift_counter - s_drift0),
                    (unsigned long long)(g_ghost_import_counter - s_ghost0),
                    (unsigned long long)(g_gpu_arena_acquire_counter - s_arena0),
                    s_np0, NumPart);
            fflush(stderr);
            /* one-shot self-check at runner end (no loop); local PHASE0 emit + return
             * follow, no intervening collective -- soft bad-stop + fall through, drains
             * at the next phase-boundary poll. */
            endrun(81036);
        }
    }

}

/* ============================================================================
 * run_neighbor_loop_iterative<Spec> — iterative entry point.
 *
 * Compile-time spec consistency:
 *   - IterControl == Iterative
 *   - IterScratch trivially copyable
 *   - max_iters >= 1; mode_a_csr_buffer_factor > 1.0
 *   - Spec::after_iter declared and callable (clean diagnostic naming the
 *     missing member if forgot)
 *
 * Runtime checks (fire BEFORE any arena/ghost/session touch):
 *   - num_subgroups >= 1 (caller short-circuits if globally empty)
 *   - num_subgroups > 1 requires Spec::SupportsSubgroups::value
 * ========================================================================== */

/* Helper: compile-time member detector for SupportsSubgroups (default false). */
template <typename Spec, typename = void>
struct nlr_supports_subgroups : std::false_type {};
template <typename Spec>
struct nlr_supports_subgroups<Spec, std::void_t<typename Spec::SupportsSubgroups>>
    : Spec::SupportsSubgroups {};

/* ============================================================================
 * NlrIterDriver<Spec> — constructor + destructor.
 * ========================================================================== */
template <typename Spec>
NlrIterDriver<Spec>::NlrIterDriver(neighbor_loop_args_iterative& a,
                                   const typename Spec::CallScalars& s)
    : args(a), cs(s), iter_index(0),
      local_active_total(0), global_active_total(-1),
      local_active_per_sg (a.num_subgroups, 0),
      global_active_per_sg(a.num_subgroups, 0),
      scratch_uvm(a.num_subgroups, nullptr),
      accum_uvm  (a.num_subgroups, nullptr),
      radii_uvm  (a.num_subgroups, nullptr),
      active_set_uvm  (a.num_subgroups, nullptr),
      active_set_size (a.num_subgroups, 0),
      call_slot_base  (a.num_subgroups, 0),
      /* Mode A iterative cached-CSR state: zero-init per subgroup;
       * UVM allocations land lazily on first Mode A iter dispatch. */
      mode_a_cached_gnl       (a.num_subgroups, gpu_neighbor_list_t{}),
      mode_a_csr_offset_lookup(a.num_subgroups, nullptr),
      mode_a_csr_buffered_h   (a.num_subgroups, nullptr),
      mode_a_csr_valid        (a.num_subgroups, false)
{
    using IterScratch = typename Spec::IterScratch;
    using AccumData   = typename Spec::AccumData;

    /* Build a base neighbor_loop_args view per subgroup so we can call
     * Spec::search_radius for each (subgroup, slot, particle index) tuple.
     * search_radius is host-only and predates any drift / arena work. */
    /* Where each subgroup's entries sit in the call's active_list, which is the
     * layout a Spec's call-level per-active staging is in (NlrSubgroup). A
     * caller that hands in a list which is not the ordered concatenation of its
     * subgroups would have every Spec that stages per-active state read another
     * particle's record, so the identity is checked here rather than assumed. */
    {
        int base = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            const NlrSubgroup& sgr = args.subgroups[sg];
            call_slot_base[sg] = base;
            base += (sgr.num_active_local > 0) ? sgr.num_active_local : 0;
        }
        bool consistent = (base == args.num_active);
        for (int sg = 0; consistent && sg < args.num_subgroups; sg++) {
            const NlrSubgroup& sgr = args.subgroups[sg];
            for (int k = 0; k < sgr.num_active_local; k++) {
                if (args.active_list[call_slot_base[sg] + k] != sgr.active_indices[k]) {consistent = false; break;}
            }
        }
        if (!consistent) {
            fprintf(stderr, "[NlrIterDriver<%s>] FATAL: active_list (%d entries) is not the ordered "
                            "concatenation of the %d subgroups' active_indices; per-active staging "
                            "would be read against the wrong particles.\n",
                    Spec::loop_name, args.num_active, args.num_subgroups);
            fflush(stderr);
            endrun(90001031);
        }
    }

    for (int sg = 0; sg < args.num_subgroups; sg++) {
        const NlrSubgroup& sgr = args.subgroups[sg];
        const int n = sgr.num_active_local;
        active_set_size[sg] = n;
        if (n <= 0) continue;     /* leave pointers as nullptr; pair_kernel is no-op */

        const size_t sg_bytes = (size_t) n * (sizeof(IterScratch) + sizeof(AccumData)
                                              + sizeof(double) + sizeof(int));
        scratch_uvm[sg]    = (IterScratch *) nlr_shared_alloc_bytes((size_t) n * sizeof(IterScratch), "modea_iter_scratch");
        accum_uvm  [sg]    = (AccumData   *) nlr_shared_alloc_bytes((size_t) n * sizeof(AccumData), "modea_accum_data");
        radii_uvm  [sg]    = (double      *) nlr_shared_alloc_bytes((size_t) n * sizeof(double), "modea_radii");
        active_set_uvm[sg] = (int         *) nlr_shared_alloc_bytes((size_t) n * sizeof(int), "modea_active_set");
        /* Without its four arrays this subgroup cannot be walked, so it is left in
         * the same state as one with no actives at all: sized zero, pointers null,
         * which every dispatch and the destructor already handle. The other
         * subgroups keep their state and the run stops at the next phase boundary. */
        if (!scratch_uvm[sg] || !accum_uvm[sg] || !radii_uvm[sg] || !active_set_uvm[sg]) {
            if (scratch_uvm[sg])    { Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_uvm[sg]);    scratch_uvm[sg] = nullptr; }
            if (accum_uvm[sg])      { Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(accum_uvm[sg]);      accum_uvm[sg] = nullptr; }
            if (radii_uvm[sg])      { Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm[sg]);      radii_uvm[sg] = nullptr; }
            if (active_set_uvm[sg]) { Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(active_set_uvm[sg]); active_set_uvm[sg] = nullptr; }
            active_set_size[sg] = 0;
            nlr_stop_no_staging_memory(Spec::loop_name, "the per-subgroup iteration state", sg_bytes);
            continue;
        }

        /* Byte-zero IterScratch ONCE — persists across iters. */
        std::memset(scratch_uvm[sg], 0, n * sizeof(IterScratch));
        /* AccumData NOT zeroed here; runner zeros via Spec::zero_accum at start of each iter. */

        /* Build a sub-args view for Spec::search_radius. Mirrors the per-
         * subgroup args the iterative dispatch will hand to lower-level
         * helpers (so search_radius signature stays a base neighbor_loop_args&). */
        neighbor_loop_args sub = args;          /* base slice; aux/CellP/etc. carry through */
        sub.active_list = sgr.active_indices;
        sub.num_active  = n;

        for (int slot = 0; slot < n; slot++) {
            radii_uvm[sg][slot]      = Spec::search_radius(sub, slot, sgr.active_indices[slot]);
            active_set_uvm[sg][slot] = slot;     /* {0..n-1} initial */
        }
        local_active_total += n;
    }
}

template <typename Spec>
NlrIterDriver<Spec>::~NlrIterDriver()
{
    /* DeviceContext cleanup: only fire if init actually
     * completed. Stubbed/aborted init paths leave ctx_initialized=false →
     * cleanup_device_context is NOT called, preventing free of unallocated
     * resources. Direct Spec call (NOT NlrDeviceContextCleanupGuard) — the
     * guard is RAII-only and wouldn't compose with conditional init.
     *
     * Ordering: cleanup_device_context fires
     * BEFORE the per-subgroup UVM frees below. Spec::cleanup_device_context
     * MUST only free resources it OWNS (e.g., UVM arrays it allocated in
     * populate_device_context). It MUST NOT assume the driver's per-subgroup
     * UVM arrays (scratch/accum/radii/active_set) remain valid AFTER cleanup
     * returns — those frees happen below. In practice cleanup hooks like
     * sink_feed's only free their own UVM, so this ordering is fine; the
     * note exists to flag the contract for future Specs. */
    if (ctx_initialized) {
        if constexpr (nlr_spec_has_cleanup_device_context_v<Spec>) {
            Spec::cleanup_device_context(args, ctx);
        }
    }

    /* Free Mode A cached CSR/lookup state by POINTER STATE (mode_a_csr_valid=false can mean "allocated but invalid,
     * pending rebuild" if the rebuild trigger fired but rebuild itself
     * hadn't completed yet; check pointers, not the flag).
     *
     * gpu_ngb_list_free passes the SIDX pointer so the step-persistent SIDX
     * cache is preserved across iterative calls (matches sink_env1/feed/swk
     * idiom). */
    {
        gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                             Spec::loop_name);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            if (mode_a_cached_gnl[sg].offsets != nullptr ||
                mode_a_cached_gnl[sg].neighbors != nullptr) {
                gpu_ngb_list_free(&mode_a_cached_gnl[sg], sidx);
                mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
            }
            mode_a_csr_valid[sg] = false;
            if (mode_a_csr_offset_lookup[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_offset_lookup[sg]);
                mode_a_csr_offset_lookup[sg] = nullptr;
            }
            if (mode_a_csr_buffered_h[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_buffered_h[sg]);
                mode_a_csr_buffered_h[sg] = nullptr;
            }
        }
    }

    /* Arena cleanup ONCE per call (matches non-iter run_mode_a line 1718). */
    if (arena_acquired) {
        gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);
        arena_acquired = false;
    }

    /* Imported-ghost cleanup ONCE per call (matches non-iter
     * line 2094-2097). Only fires for Mode A multi-rank paths that imported.
     * Runner's Mode A imports an exact-query ghost pool sized to the iter's
     * actives+radii; the caller (e.g. density()) is responsible for any
     * downstream handoff pool it needs AFTER runner-return (post-finalize
     * fresh broad import — not "keep this exact-query pool alive"). */
    if (ghost_import_done && NTask > 1) {
        ghost_exchange_cleanup();
        ghost_import_done = false;
    }

    for (int sg = 0; sg < args.num_subgroups; sg++) {
        if (scratch_uvm[sg])    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(scratch_uvm[sg]);
        if (accum_uvm[sg])      Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(accum_uvm[sg]);
        if (radii_uvm[sg])      Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_uvm[sg]);
        if (active_set_uvm[sg]) Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(active_set_uvm[sg]);
    }
}

/* ============================================================================
 * Path-specific DeviceContext init.
 *
 * Mode B: bind to caller's args.P/CellP — Mode B reads owner-local per-particle
 *         state directly (lazy-drift contract). Static_asserts mirror run_mode_b_*
 *         (DeviceContext base + trivially copyable).
 * Mode A: bind to arena-resident gpu_particles_arena_P()/CellP() — Mode A's
 *         GPU NGL build + kernel walk operate on the arena copy. Body lands
 *         alongside the arena_acquire call site.
 * ========================================================================== */
template <typename Spec>
void NlrIterDriver<Spec>::initialize_device_context_mode_b()
{
    using DeviceCtx = typename Spec::DeviceContext;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; "
                  "captured by value into Kokkos device lambdas");

    /* Double-init guard: accidental call sequences (e.g., path mis-route calling
     * both inits) would silently double-populate / double-cleanup. Catch
     * loudly here. */
    if (ctx_initialized) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::initialize_device_context_mode_b] FATAL: "
                "ctx_initialized=true on entry. Double-init would orphan resources.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        /* Soft bad-stop + return: the valid first-init ctx is left intact
         * (no re-population). This function issues no MPI, so the early
         * return cannot desync a peer; drains at the next phase poll. */
        endrun(81209);
        return;
    }

    ctx.P         = args.P;
    ctx.CellP     = args.CellP;
    ctx.num_total = args.num_total;

    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    ctx_initialized = true;
}

template <typename Spec>
void NlrIterDriver<Spec>::initialize_device_context_mode_a_after_arena()
{
    using DeviceCtx = typename Spec::DeviceContext;
    static_assert(std::is_base_of<NeighborLoopDeviceContextBase, DeviceCtx>::value,
                  "Spec::DeviceContext must publicly derive from NeighborLoopDeviceContextBase");
    static_assert(std::is_trivially_copyable<DeviceCtx>::value,
                  "Spec::DeviceContext must be trivially copyable; "
                  "captured by value into Kokkos device lambdas");

    /* Double-init guard. */
    if (ctx_initialized) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::initialize_device_context_mode_a_after_arena] FATAL: "
                "ctx_initialized=true on entry. Double-init would orphan resources.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        /* Soft bad-stop + return: the valid first-init ctx is left intact
         * (no re-population). This function issues no MPI, so the early
         * return cannot desync a peer; drains at the next phase poll. */
        endrun(81210);
        return;
    }
    /* Arena must already have been acquired. */
    if (!arena_acquired) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::initialize_device_context_mode_a_after_arena] FATAL: "
                "arena not acquired. Call acquire_arena_and_init_ctx_mode_a() instead.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        /* Symmetric lifecycle-contract violation. Soft bad-stop + return WITHOUT binding
         * ctx.P to an unacquired arena; the outer runner's nlr:iter_context_init poll
         * drains all ranks before any device dispatch. */
        endrun(81211); return;
    }

    /* Bind to arena-resident P_gpu / CellP_gpu. Use the driver's
     * args: if ghost import ran above, these
     * point at the refreshed POST-IMPORT global P/CellP/NumPart; otherwise
     * they equal the original base args. CellP=NULL is legitimate (gas-free). */
    if (!nlr_args_view_is_live(args)) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::initialize_device_context_mode_a_after_arena] FATAL: the "
                "argument view no longer describes the live particle storage (P=%p vs %p, "
                "num_total=%d vs %d). Binding a device context to it would read released memory.\n",
                Spec::loop_name, (void *) args.P, (void *) P, args.num_total, NumPart);
            fflush(stderr);
        }
        /* Soft bad-stop + return WITHOUT binding ctx: same shape as the two lifecycle
         * violations above, drained at the runner's next poll before any dispatch. */
        endrun(81215); return;
    }
    ctx.P         = gpu_particles_arena_P();
    ctx.CellP     = (args.CellP != nullptr) ? gpu_particles_arena_CellP() : nullptr;
    ctx.num_total = args.num_total;

    if constexpr (nlr_spec_has_extended_device_context_v<Spec>) {
        Spec::populate_device_context(args, ctx);
    }
    ctx_initialized = true;
}

template <typename Spec>
void NlrIterDriver<Spec>::acquire_arena_and_init_ctx_mode_a()
{
    if (arena_acquired) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[NlrIterDriver<%s>::acquire_arena_and_init_ctx_mode_a] FATAL: "
                "arena already acquired. Single-acquire-per-call contract violated.\n",
                Spec::loop_name);
            fflush(stderr);
        }
        /* Symmetric lifecycle-contract violation. Soft bad-stop + return WITHOUT a second
         * acquire (the first acquire's state stays intact); the outer runner's
         * nlr:iter_context_init poll drains all ranks before any device dispatch. */
        endrun(81212); return;
    }

    /* === (1) Imported-ghost prep ONCE per call ===
     * Mode A requires neighbors that live on peer ranks to be imported as
     * ghosts. Mirrors non-iter run_neighbor_loop lines 2037-2052.
     *
     * Multi-subgroup union semantics: mask = OR of all
     * subgroups' j_type_bitmask; actives = concatenation of all subgroups'
     * full active_indices (initial iter-0 set; later iters' compactions
     * narrow but don't grow); radii = matching concatenation oversized via
     * mode_a_csr_buffer_factor. Single-subgroup case: trivial length-1
     * union = original behavior. Multi-subgroup over-imports in cross-bm
     * cases (deliberate physics/perf tradeoff).
     *
     * On rebuild via rebuild_mode_a_arena_and_ctx_for_current_active_union,
     * the union covers only CURRENT compacted actives (post-Converged-removal)
     * so the rebuild import doesn't over-cover converged slots. */
    if (NTask > 1) {
        std::vector<int>    union_actives;
        std::vector<double> union_radii_oversized;
        unsigned int        mask_union = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            mask_union |= args.subgroups[sg].j_type_bitmask;
            const NlrSubgroup& sgr = args.subgroups[sg];
            /* A subgroup whose per-active radii could not be staged contributes
             * nothing to the import. The import itself is collective, so this
             * rank still enters it -- just asking for no ghosts of its own. */
            if (radii_uvm[sg] == nullptr) { continue; }
            for (int k = 0; k < sgr.num_active_local; k++) {
                union_actives.push_back(sgr.active_indices[k]);
                union_radii_oversized.push_back(radii_uvm[sg][k] * Spec::mode_a_csr_buffer_factor);
            }
        }
        const int union_n = (int)union_actives.size();
        gizmo_request_filtered_ghost_import_fresh(Spec::loop_name,
                                                   Spec::search_mode,
                                                   mask_union,
                                                   (union_n > 0) ? union_actives.data() : nullptr,
                                                   union_n,
                                                   (union_n > 0) ? union_radii_oversized.data() : nullptr,
                                                   args.ghost_safety_factor,
                                                   Spec::radius_policy,
                                                   nlr_spec_symmetric_j_radius_scale<Spec>());
        ghost_import_done = true;

        /* Refresh args from globals (ghost import grew NumPart and
         * may have realloc'd P/CellP). Matches non-iter line 2049-2051. */
        args.num_total = NumPart;
        args.P         = P;
        args.CellP     = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;
    }

    /* === (2) Arena acquire ONCE per call, with refreshed args === */
    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(args.num_total,
                                 args.P, args.CellP);
    arena_acquired = true;

    /* === (3) Bind ctx to arena-resident pointers + populate extended DeviceContext === */
    initialize_device_context_mode_a_after_arena();
}

/* Mode D's context: the same device-resident arrays, WITHOUT the ghost import.
 *
 * Mode A imports because it walks a pool that has to contain the neighbours
 * living on other ranks.  A fused walk never sees a ghost -- the traversal stops
 * at the owned count -- and the queries that need another rank's particles are
 * shipped to that rank instead.  Calling the Mode A entry here would perform the
 * whole collective import and then discard every ghost it fetched, which is the
 * one cost this path exists to avoid.
 *
 * Skipping it is collective-safe because the decision that selected Mode D was
 * itself collective and unanimous: either every rank takes this entry or none
 * does, so no rank is left waiting in an import its peers skipped. */
template <typename Spec>
void NlrIterDriver<Spec>::acquire_arena_and_init_ctx_mode_d()
{
    if (arena_acquired) {return;}

    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(args.num_total, args.P, args.CellP);
    arena_acquired = true;

    initialize_device_context_mode_a_after_arena();
}

/* ============================================================================
 * rebuild_mode_a_arena_and_ctx_for_current_active_union — self-sufficient
 * rebuild method (no parameters).
 *
 * Builds the union from driver's own subgroup state. Used by the Mode A
 * outer-iter pre-dispatch invalidation sweep (step 8). Invalidates ALL
 * subgroup CSR caches as part of the lifecycle (pitfall 2: cross-subgroup
 * neighbor indices become stale after arena re-acquire).
 *
 * Union semantics (pitfall 4): mask = OR of all globally-active subgroups'
 * j_type_bitmask; active list = concatenation of current compacted actives
 * across globally-active subgroups; radii = matching concatenation of
 * radii_uvm[sg][slot] * mode_a_csr_buffer_factor. Over-imports in
 * multi-bm cases (each active gets ghosts of all union types instead of
 * just its own mask). Cosmological FIRE typical = single bm-group = no
 * over-import. Documented as deliberate physics/perf tradeoff.
 *
 * Pre-condition: drv.global_active_per_sg has been computed (set by the
 * outer iter loop's per-iter Allreduce). On iter 0 it may all be 0 (not
 * yet Allreduced); caller should NOT invoke this method at iter 0 — the
 * iter-0 ghost import + arena acquire happens in
 * acquire_arena_and_init_ctx_mode_a per the original path.
 * ========================================================================== */
template <typename Spec>
void NlrIterDriver<Spec>::rebuild_mode_a_arena_and_ctx_for_current_active_union()
{
    if (NTask <= 1) {
        /* Single rank: no ghost pool to manage. Just invalidate all CSR
         * caches; per-subgroup dispatch will rebuild local CSRs on first
         * access against the unchanged arena/pool. */
        gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                             Spec::loop_name);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            if (mode_a_cached_gnl[sg].offsets != nullptr ||
                mode_a_cached_gnl[sg].neighbors != nullptr) {
                gpu_ngb_list_free(&mode_a_cached_gnl[sg], sidx);
                mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
            }
            mode_a_csr_valid[sg] = false;
            if (mode_a_csr_offset_lookup[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_offset_lookup[sg]);
                mode_a_csr_offset_lookup[sg] = nullptr;
            }
            if (mode_a_csr_buffered_h[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_buffered_h[sg]);
                mode_a_csr_buffered_h[sg] = nullptr;
            }
        }
        return;
    }

    /* Build union from current active_set across globally-active subgroups. */
    std::vector<int>    union_actives;
    std::vector<double> union_radii_oversized;
    unsigned int        mask_union = 0;
    for (int sg = 0; sg < args.num_subgroups; sg++) {
        if (global_active_per_sg[sg] <= 0) continue;     /* skip globally-converged */
        mask_union |= args.subgroups[sg].j_type_bitmask;
        const NlrSubgroup& sgr = args.subgroups[sg];
        for (int k = 0; k < active_set_size[sg]; k++) {
            int slot = active_set_uvm[sg][k];
            union_actives.push_back(sgr.active_indices[slot]);
            union_radii_oversized.push_back(radii_uvm[sg][slot] * Spec::mode_a_csr_buffer_factor);
        }
    }

    /* === (0) Invalidate ALL subgroup CSR caches (avoids cross-subgroup
     * staleness after arena teardown). Free by pointer state. */
    {
        gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                             Spec::loop_name);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            if (mode_a_cached_gnl[sg].offsets != nullptr ||
                mode_a_cached_gnl[sg].neighbors != nullptr) {
                gpu_ngb_list_free(&mode_a_cached_gnl[sg], sidx);
                mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
            }
            mode_a_csr_valid[sg] = false;
            if (mode_a_csr_offset_lookup[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_offset_lookup[sg]);
                mode_a_csr_offset_lookup[sg] = nullptr;
            }
            if (mode_a_csr_buffered_h[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(mode_a_csr_buffered_h[sg]);
                mode_a_csr_buffered_h[sg] = nullptr;
            }
        }
    }

    /* === (1) cleanup_device_context for extended Specs === */
    if (ctx_initialized) {
        if constexpr (nlr_spec_has_cleanup_device_context_v<Spec>) {
            Spec::cleanup_device_context(args, ctx);
        }
        ctx_initialized = false;
    }

    /* === (2) END current arena view BEFORE mutating the global ghost pool === */
    if (arena_acquired) {
        gpu_particles_arena_mark_clean_after_scatter(Spec::loop_name);
        arena_acquired = false;
    }

    /* === (3) Destroy current ghost pool === */
    if (ghost_import_done) {
        ghost_exchange_cleanup();
        ghost_import_done = false;
    }

    /* === (4) Reimport with UNION radii / actives / mask (deliberate over-import,
     * an accepted tradeoff) === */
    {
        const int union_n = (int)union_actives.size();
        gizmo_request_filtered_ghost_import_fresh(Spec::loop_name,
                                                   Spec::search_mode,
                                                   mask_union,
                                                   (union_n > 0) ? union_actives.data() : nullptr,
                                                   union_n,
                                                   (union_n > 0) ? union_radii_oversized.data() : nullptr,
                                                   args.ghost_safety_factor,
                                                   Spec::radius_policy,
                                                   nlr_spec_symmetric_j_radius_scale<Spec>());
        ghost_import_done = true;
    }

    /* === (5) Refresh args from post-import globals === */
    args.num_total = NumPart;
    args.P         = P;
    args.CellP     = (gizmo_host_all_ptr()->TotN_gas > 0) ? CellP : nullptr;

    /* === (6) Fresh arena view of new pool === */
    GIZMO_GPU_ENSURE_ALL_FRESH();
    gpu_particles_arena_set_site(Spec::loop_name);
    gpu_particles_arena_acquire(args.num_total,
                                 args.P, args.CellP);
    arena_acquired = true;

    /* === (7) Rebind ctx + (8) populate extended ctx for NEW pool === */
    initialize_device_context_mode_a_after_arena();

    /* Do NOT re-fire reset_per_iter_device_context
     * here. The outer iter loop's step (a-pre) calls this method BEFORE
     * step (a0), which fires the reset hook exactly once per outer iter
     * over the freshly-populated ctx. Re-firing here would cause a
     * double-reset on rebuild iters and violate the "fires once per outer
     * iter" contract used by the synthetic harness validation. */
}

/* ============================================================================
 * Per-iter, per-subgroup dispatch helpers (Mode B local only).
 *
 * Each helper takes the driver, a subgroup index, and runs ONE iteration's
 * pair_kernel walk for that subgroup writing into driver-owned accum_uvm.
 * apply_active_writeback is NOT called here — that's final-only after the
 * outer iter loop completes.
 *
 * The Mode B local helper composes existing lower-level helpers
 * (build_self_actives_host_pre_drift / collect_candidates_pre_drift /
 * drift_and_filter_candidates / evaluate_pairs_post_drift). SSOT preserved with
 * existing run_mode_b_local — same helper chain, just driver-owned output
 * buffer instead of stack vector.
 *
 * Mode A iter lands via the sibling nlr_iter_dispatch_subgroup_mode_a helper.
 *
 * DEVICE CONTEXT LIFETIME:
 * The dispatch helper uses driver-owned drv.ctx (populated once via
 * NlrIterDriver::initialize_device_context_mode_b at iter-0 entry, cleaned
 * once in destructor via cleanup_device_context if extended). Per-iter
 * reset is handled at the outer iter loop via Spec::reset_per_iter_device_context.
 * No per-iter / per-subgroup ctx rebuild here; the older "rebuilds DeviceCtx
 * per-iter + per-subgroup" hack is gone.
 * ========================================================================== */
/* ============================================================================
 * nlr_iter_dispatch_subgroup_mode_b_remote<Spec>.
 *
 * Per-iter Mode B REMOTE dispatch for one subgroup. Composes the same
 * extracted helper (mode_b_remote_evaluate_into_buffer<Spec, false>) the
 * non-iterative wrapper uses — SSOT preserved. The driver-owned compacted
 * AccumData buffer is the only difference. apply_active_writeback NOT
 * called here — that's final-only after the outer iter loop.
 */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_b_remote(NlrIterDriver<Spec>& drv, int sg)
{
    using AccumData = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_size[sg];
    /* Note: even if n_compacted == 0 on this rank, the helper MUST be
     * entered because Mode B remote uses collectives (alltoallv exchange
     * of queries / replies). Empty rank participates with N=0. */

    /* Compact active particle indices + radii for the still-active slots. */
    std::vector<int>    active_particle_indices(n_compacted);
    std::vector<int>    call_slots(n_compacted);
    std::vector<double> radii_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
        call_slots[k]              = drv.call_slot_base[sg] + slot;
        radii_compacted[k]         = drv.radii_uvm[sg][slot];
    }

    neighbor_loop_args sub = drv.args;
    sub.active_list      = (n_compacted > 0) ? active_particle_indices.data() : nullptr;
    sub.active_call_slot = (n_compacted > 0) ? call_slots.data() : nullptr;
    sub.num_active       = n_compacted;

    /* Driver-owned compacted AccumData buffer. Helper writes into this;
     * we scatter back into driver.accum_uvm[sg][slot] for active slots
     * after the helper returns. Explicit Spec::zero_accum per slot:
     * defensive against future evaluate_pairs
     * variants + makes the per-iter AccumData zero contract explicit at
     * the caller level. */
    std::vector<AccumData> accums_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        Spec::zero_accum(accums_compacted[k]);
    }

    mode_b_remote_evaluate_into_buffer<Spec>(
        sub,
        radii_compacted.data(),
        drv.cs,                          /* driver-owned CallScalars */
        drv.ctx,                         /* driver-owned DeviceContext */
        (unsigned int)sgr.j_type_bitmask, /* per-subgroup mask */
        (n_compacted > 0) ? accums_compacted.data() : nullptr);

    /* Scatter compacted accums back into driver-owned per-slot accum_uvm.
     * Slots NOT in active_set_uvm keep their last-evaluated value (the
     * v4.4 invariant: converged slots' final accum persists for
     * apply_active_writeback). */
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        drv.accum_uvm[sg][slot] = accums_compacted[k];
    }
}

/* ============================================================================
 * Mode D — the walk and the pair kernel in one device pass.
 *
 * Mode A brings the neighbours to the seeker: it imports the particles another
 * rank owns, builds a neighbour list over them, and then evaluates.  Mode D
 * sends the seeker to the neighbours instead, and evaluates where they already
 * live, so the list is never built and the particles are never moved.  What that
 * removes is not arithmetic -- the same pair kernel runs on the same pairs --
 * but the import, the list, and the copies that surround them.
 *
 * The seeker cannot tell where it is.  A query that started on this rank and one
 * that arrived from another are the same object to the walk and to the kernel;
 * the only difference is which accumulator the result lands in, and that is the
 * caller's business, not the kernel's.  This file's half is the local one.
 *
 * Availability is a property of the loop, not of the caller's name: the pair
 * kernel has to be one that only reduces into the seeker's own accumulator, or
 * evaluating it on a rank that does not own the seeker would need writes shipped
 * back.  That is exactly what the eval tier already records, so it is what is
 * asserted, and no loop is named here.
 * ========================================================================== */

/* What happens when the walk reaches a particle this rank owns.
 *
 * The traversal decides which indices are reachable and which are ghosts; this
 * decides whether a reachable one is a neighbour, and it applies the same three
 * tests the host walker applies in mesh/mode_b_local_walker.cc -- allowed type,
 * positive mass, and the geometric accept -- through the same shared predicate,
 * so the two cannot answer differently about the same pair.  The reach is the
 * query's alone, which is what makes this the ONEWAY case: no neighbour supply
 * radius is consulted, and none is mirrored to the device to consult. */
template <typename Spec>
struct NlrModeDReduceLeaf {
    const typename Spec::DeviceContext *ctx;
    const typename Spec::ActiveData    *active;
    typename Spec::AccumData           *accum;
    typename Spec::ScatterData         *scatter;
    unsigned int                        supply_mask;
    /* By pointer, not by value: the leaf is built per work item inside the
     * kernel, and a CallScalars is 120-200 B, so copying it into every leaf
     * would put that much in local memory per item for a value every item
     * shares. The pointee is the launching lambda's OWN by-value capture --
     * device-resident -- which is why the leaf must keep being constructed
     * INSIDE the kernel. Hoisting that construction out would leave this
     * pointing at a host stack object and fault on device. */
    const typename Spec::CallScalars   *cs;
    struct GxMotionTargetSet            motion_targets{};   /* for a loop that writes neighbour motion */

    KOKKOS_INLINE_FUNCTION
    void visit(int j, double qx, double qy, double qz, double reach)
    {
        const struct particle_data &Pj = ctx->P[j];
        if(!(supply_mask & (1u << (unsigned int)Pj.Type))) {return;}
        /* Read atomically: a pair kernel that deposits mass into its neighbours
         * may be adding to this one from another lane at this moment. */
        if(Kokkos::atomic_load(&Pj.Mass) <= 0) {return;}
        if(!gx_pair_accept_wrap_and_test(qx - (double)Pj.Pos[0],
                                         qy - (double)Pj.Pos[1],
                                         qz - (double)Pj.Pos[2],
                                         reach, 0.0, NGB_SEARCH_ONEWAY)) {return;}
        if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {gx_motion_target_mark(motion_targets, j);}
        IdentitySidecar id{};
        typename Spec::NeighborData nb = Spec::load_neighbor(*ctx, j, id, *active);
        Spec::pair_kernel(*active, nb, *accum, *scatter, *cs);
    }
};

/* What a walk's anomaly report means, for the three sites that stop the run on
 * one.  The states are distinct and so are their causes, so a single message
 * naming only the tree would send the reader looking in the wrong place. */
static const char *nlr_walk_anomaly_text(int code)
{
    switch(code) {
    case GX_WALK_ANOMALY_MALFORMED_TREE:
        return "a query reached an index in the gap between the particle slots and the node base; the tree is malformed";
    case GX_WALK_ANOMALY_TOUCHED_SET_FULL:
        return "the touched-set list was shorter than the distinct set the recording walk put in it";
    default:
        return "an unrecognised walk anomaly";
    }
}

/* What happens when a RECORDING walk reaches a particle this rank owns.
 *
 * The same traversal as the evaluating leaf above, with the evaluation removed:
 * this one only writes down which local particles the walk touches, so that just
 * those can be brought current before the walk runs again and evaluates them.
 *
 * It records BEFORE the geometric accept, and that is the whole correctness
 * argument rather than an efficiency choice.  The accept test compares positions,
 * and the positions are exactly what has not been brought current yet; a particle
 * that will move into range is rejected here and would then be evaluated on its
 * undrifted position in the second pass, which is worse than not de-fusing at
 * all.  The type and mass tests are safe to apply because a drift changes
 * neither, so they narrow the set without being able to drop anything the second
 * pass can accept.
 *
 * The traversal reads node geometry and the particle-slot links, none of which a
 * drift writes, so the pass that evaluates reaches exactly the set this pass
 * recorded.  That is what makes this complete rather than approximate.
 *
 * It reads no Spec member, so one policy serves every loop the fused backend can
 * take rather than one per Spec. */
struct NlrRecordLeaf {
    const struct particle_data *P;
    struct GxTouchedSet         ts;
    unsigned int                supply_mask;
    int                        *anomaly;

    KOKKOS_INLINE_FUNCTION
    void visit(int j, double, double, double, double) const
    {
        const struct particle_data &Pj = P[j];
        if(!(supply_mask & (1u << (unsigned int)Pj.Type))) {return;}
        if(Pj.Mass <= 0) {return;}
        /* Claim the slot for this generation.  The exchange returns what was
         * there, so exactly one work item sees a value other than the current
         * generation and exactly one item appends -- several actives reaching the
         * same particle is the ordinary case, not the exception. */
        if(Kokkos::atomic_exchange(&ts.seen[j], ts.gen) == ts.gen) {return;}
        const int slot = Kokkos::atomic_fetch_add(ts.counter, 1);
        if(slot < ts.capacity) {
            ts.list[slot] = j;
        } else {
            /* Unreachable: the stamp admits each owned slot once per generation
             * and the list is as long as there are owned slots.  Reported rather
             * than dropped, because a dropped index is a particle silently
             * evaluated at a stale position -- the one failure this design must
             * not be able to have quietly. */
            Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_TOUCHED_SET_FULL);
        }
    }
};

/* Bring current exactly the local particles a fused walk is about to reach.
 *
 * Runs the traversal once with the recording policy, waits for it, advances the
 * distinct set it recorded, and leaves the cursor ready for the next pass.  The
 * evaluation that follows then walks the same tree with the same queries and
 * reaches the same leaves, now current.
 *
 * Called immediately before each evaluation rather than once per call, because
 * the search radius is what the iteration changes: a later pass can reach
 * further than an earlier one.  It costs little to repeat -- a particle advanced
 * for an earlier pass stays current for the rest of the call, and the generation
 * is per call, so the second and later passes record only what is new.
 *
 * One pass covers this rank's own queries and the ones its peers sent, in a
 * single traversal over both.  They ask the same question of the same tree and
 * write the same touched set, and which of the two a given local particle was
 * reached by is an accident of where the domain boundary fell, so separating
 * them bought nothing and cost a launch, a fence and a second drift pass.  The
 * two kinds keep the two ENTRIES the traversal already has -- this rank's root,
 * or the start nodes a peer exported -- and `query` says which this work item
 * takes by yielding a start list or none.
 *
 * `query` is device-callable and yields this work item's position, reach, and
 * either its start-node list or a negative count meaning "walk from the root".
 * Callers keep the two kinds in CONTIGUOUS index ranges so that the entry a
 * work item takes does not alternate between neighbouring lanes. */
template <class QueryFn>
static void nlr_record_and_drift(const struct particle_data *P,
                                 unsigned int supply_mask,
                                 const GxDeviceTreeView &tree,
                                 int n, int *anomaly,
                                 const char *label,
                                 QueryFn query)
{
    if(n <= 0) {return;}
    /* Nothing to discover when the whole rank is already uniform: something else
     * in the step -- a tree build, a decomposition, an output -- has published
     * the full-drift certificate, and every leaf this walk can reach is current
     * by that proof.  Recording them would be a traversal spent to learn that
     * there is no work, which on an all-active call is the largest traversal of
     * the step.  This is the same early-out the neighbour-list hook applies for
     * the same reason.
     *
     * The local certificate alone is the right test here, unlike there: this
     * path imports no ghosts, and the traversal stops at the owned slots, so
     * there is no imported segment to vouch for.
     *
     * Rank-local, and safe to be: neither this test nor the drift it guards
     * enters a collective, so a rank that skips and a rank that does not still
     * meet at the same place. */
    if(gizmo_full_drift_ti() == All.Ti_Current) {return;}
    const struct GxTouchedSet ts = gx_touched_set_view();
    /* The preparation declines collectively when it cannot hold this, so an
     * absent workspace here is not a state to recover from -- but it must not be
     * walked past either, because the drift that would not happen is a particle
     * evaluated at a stale position.  Reported through the channel the caller
     * already treats as fatal, rather than returned, so this adds no path out of
     * an exchange that a peer is waiting on. */
    if(!ts.seen || !ts.list || !ts.counter) {
        Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_TOUCHED_SET_FULL);
        return;
    }
    GIZMO_GPU_ENSURE_ALL_FRESH();
    nlr_walk_for_sources(label, n, KOKKOS_LAMBDA(int kk) {
        double qx = 0, qy = 0, qz = 0, reach = 0;
        const int *start_nodes = nullptr;
        int n_start = -1;                 /* negative: this one walks from the root */
        query(kk, qx, qy, qz, reach, start_nodes, n_start);
        NlrRecordLeaf leaf{P, ts, supply_mask, anomaly};
        if(n_start < 0) {
            gx_device_tree_walk_from_root(qx, qy, qz, reach, tree, leaf, anomaly, supply_mask);
        } else if(n_start > 0) {
            /* Zero start nodes means nothing on this rank was exported to this
             * query, so there is nothing of ours for it to reach -- the same
             * skip the evaluating walk makes at the same data. */
            gx_device_tree_walk(qx, qy, qz, reach, start_nodes, n_start, tree, leaf, anomaly, supply_mask);
        }
    });
    gx_touched_set_drift_and_mark(All.Ti_Current);
}

/* One work item per active: build its query, walk this rank's tree from the
 * root, reduce into its own accumulator.
 *
 * The query is built HERE, on the device, from the resident particle arrays.
 * The host-walker path builds the equivalent on the host, one active at a time,
 * which is right where it runs -- a few actives -- and wrong here, where it
 * would be a serial host pass over the canonical arrays proportional to the
 * active count, on the steps that have the most actives.
 *
 * `radii` rather than the query's own stored radius: the search radius is the
 * runner's, it is what the iteration mutates between passes, and it is what the
 * host walk is given.  Reading it from anywhere else would let the two drift.
 *
 * Ownership is one accumulator per work item, so nothing is shared and nothing
 * needs an atomic.  A traversal that gave several work items a share of one
 * query -- one team per node pair, say -- would not have that property, which is
 * why it is stated rather than assumed. */
template <typename Spec>
static void nlr_mode_d_local_reduce(const typename Spec::DeviceContext &ctx,
                                    const typename Spec::CallScalars &cs,
                                    const int *active_idx,
                                    const int *active_slot,
                                    const double *radii,
                                    int n,
                                    unsigned int supply_mask,
                                    const GxDeviceTreeView &tree,
                                    typename Spec::AccumData *accums_out,
                                    int *anomaly)
{
    using ActiveData  = typename Spec::ActiveData;
    using ScatterData = typename Spec::ScatterData;

    static_assert(nlr_spec_modeb_eval_omp<Spec>() != ModeBEvalOMP::SerialOnly,
                  "Mode D evaluates a seeker's pairs on the rank that owns the neighbours, so a "
                  "pair kernel's neighbour-side writes land on the owner's own particles and need "
                  "no writeback -- but they land from many device lanes at once, so a kernel that "
                  "must run serially (a read-then-write of live neighbour state) cannot be served.");
    static_assert(Spec::search_mode == MODE_B_SEARCH_ONEWAY,
                  "Mode D prunes on the query's reach alone. A symmetric search also needs the "
                  "per-type supply bands, which are not mirrored to the device.");
    static_assert(!nlr_spec_has_bind_active_to_eval_context_v<Spec>,
                  "This Spec rebinds each active to the evaluating context before its pair kernel "
                  "runs, which the host backend does inside evaluate_pairs_post_drift. The fused "
                  "backend evaluates in a device kernel and performs no such rebind, so a received "
                  "active would still carry the sending rank's pointers. Port the rebind into the "
                  "device path before serving this Spec.");

    /* ADMISSION INVARIANT, not checkable here: EVERY ACTIVE THIS BACKEND IS GIVEN
     * IS ALREADY CURRENT AT All.Ti_Current.
     *
     * The three tests above are compile-time properties of the Spec.  This fourth
     * requirement is a property of the LIST the Spec hands in, so it cannot be
     * asserted -- but it is load-bearing in the same way, and this is where a new
     * Spec is admitted, so it is stated here rather than only where the queries
     * get built.
     *
     * It holds today by construction: core/run.cc's sync point drifts every
     * particle in every active time bin to All.Ti_Current, ActiveParticleList is
     * built from those same bins, nlr_build_active_list filters that list, and
     * nothing advances All.Ti_Current again inside the step.
     *
     * What breaks if a future Spec supplies an independent active list: this
     * backend snapshots its queries BEFORE the recording pass -- frozen ActiveData
     * on the transport path, packed envelopes on the peer path -- so a stale
     * active is baked into a snapshot no later drift can repair.  Drifting at the
     * record site would patch only the one site that re-derives its query, and
     * would leave the other two silently wrong while looking defended.  The only
     * correct remedy is to establish currentness BEFORE the first query snapshot
     * or transport packing, as part of admitting that Spec -- a change to this
     * backend's admission, not a defensive line inside the record loop. */


    if(n <= 0) {return;}

    /* The traversal reads All.BoxSize through the wrap macros, which resolve to
     * this unit's mirror; without this the box reads as zero on device and
     * periodic wrapping stops, silently and only on HIP. */
    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* Record what this walk will reach and bring just those current, before it
     * runs for real.  Same tree, same queries, so the same leaves. */
    nlr_record_and_drift(
        ctx.P, supply_mask, tree, n, anomaly, "nlr_mode_d_self_record",
        KOKKOS_LAMBDA(int kk, double &qx, double &qy, double &qz, double &reach,
                      const int *&start_nodes, int &n_start) {
            const ActiveData a_rec = Spec::load_active(ctx, active_slot[kk], active_idx[kk], radii[kk], cs);
            qx = (double)a_rec.pos[0]; qy = (double)a_rec.pos[1]; qz = (double)a_rec.pos[2];
            reach = radii[kk];
            /* Single-rank: every query is this rank's own and walks from the root. */
            start_nodes = nullptr; n_start = -1;
        });

    const struct GxMotionTargetSet motion_targets = gx_motion_target_view();
    nlr_walk_for_sources(Spec::loop_name, n, KOKKOS_LAMBDA(int kk) {
        Spec::zero_accum(accums_out[kk]);
        const int i = active_idx[kk];
        ActiveData  a = Spec::load_active(ctx, active_slot[kk], i, radii[kk], cs);
        ScatterData s{};
        NlrModeDReduceLeaf<Spec> leaf{&ctx, &a, &accums_out[kk], &s, supply_mask, &cs};
        leaf.motion_targets = motion_targets;
        gx_device_tree_walk_from_root((double)a.pos[0], (double)a.pos[1], (double)a.pos[2],
                                      radii[kk], tree, leaf, anomaly, supply_mask);
    });
}

/* Build this rank's queries where its particles already are.
 *
 * The host form of this walks P[i] and CellP[i] once per active, scattered, to
 * fill one ActiveData each. That is correct, and on the transport path it is
 * also the canonical arrays surfacing on the host on exactly the calls with the
 * most actives. The particles are resident on the device, so the same
 * Spec::load_active runs there instead and writes into a buffer the host can
 * still read -- the export walk and the envelope pack both need to.
 *
 * The index and radius arrays it reads from ARE staged host-side, and that is
 * deliberate rather than overlooked: they are runner-owned, contiguous, and one
 * int and one double per active, against the many scattered particle-field reads
 * they replace. The remaining host touch is small and named, not absent.
 *
 * Returns the buffer, or nullptr if the device memory could not be had -- in
 * which case the caller builds on the host and gets identical objects.
 */
template <typename Spec, typename DeviceCtx>
static typename Spec::ActiveData *
nlr_build_self_actives_on_device(const neighbor_loop_args& args,
                                 const DeviceCtx& ctx,
                                 const double *radii,
                                 const typename Spec::CallScalars& cs,
                                 int n,
                                 NlrDeviceBuiltActives<Spec>& owner)
{
    using ActiveData = typename Spec::ActiveData;

    if(n <= 0) {return nullptr;}

    ActiveData *act_d  = (ActiveData *) nlr_shared_alloc_bytes((size_t)n * sizeof(ActiveData), "moded_self_actives");
    int        *idx_d  = (int *)        nlr_shared_alloc_bytes((size_t)n * sizeof(int),        "moded_self_idx");
    double     *rad_d  = (double *)     nlr_shared_alloc_bytes((size_t)n * sizeof(double),     "moded_self_radii");
    int        *slot_d = (int *)        nlr_shared_alloc_bytes((size_t)n * sizeof(int),        "moded_self_slot");

    if(!act_d || !idx_d || !rad_d || !slot_d) {
        if(act_d)  {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(act_d);}
        if(idx_d)  {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(idx_d);}
        if(rad_d)  {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(rad_d);}
        if(slot_d) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(slot_d);}
        return nullptr;
    }

    for(int k = 0; k < n; k++) {
        idx_d[k]  = args.active_list[k];
        rad_d[k]  = radii[k];
        slot_d[k] = args.active_call_slot ? args.active_call_slot[k] : k;
    }

    /* load_active reads All.* through this unit's mirror; without the belt those
     * read as zero on device, silently, and only on HIP. */
    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* Named apart from the walk kernels this loop also launches: they run under
     * the loop's own name, so sharing it would leave a device fault in the query
     * build indistinguishable from one in the traversal. */
    gizmo_gpu_kernel_launch("nlr_mode_d_build_queries", n, KOKKOS_LAMBDA(int k) {
        act_d[k] = Spec::load_active(ctx, slot_d[k], idx_d[k], rad_d[k], cs);
    });

    /* The indices, slots and radii were only ever the kernel's input. The
     * queries themselves outlive this call and belong to the owner. */
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(idx_d);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(rad_d);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(slot_d);

    owner.p = act_d;
    return act_d;
}

/* Answering received queries and this rank's own in ONE device pass, no
 * candidate list.
 *
 * Same contract as NlrPeerAnswerHostWalk for the received half -- queries in,
 * one accumulator each out -- reached by walking each query's exported subtrees
 * on the device and evaluating the pair kernel where the walk lands.  Nothing is
 * COLLECTED, so there is no candidate list to drift -- but the locals this walk
 * lands on are this rank's, and they are brought current by a recording pass
 * over the same queries followed by a drift of exactly what it recorded.  A
 * round can reach locals an earlier one never did, which is why it discovers
 * rather than relying on what an earlier pass found.
 *
 * `local` carries this rank's OWN queries when the round loop has some left to
 * place, and they are walked in the same two launches: one record, one
 * evaluate, over `n_local + K` work items.  Two half-filled launches apiece is
 * what this replaces, and on the steps that matter neither of them came close
 * to filling the device.  The two kinds occupy CONTIGUOUS index ranges --
 * locals first -- so the entry a work item takes does not alternate between
 * neighbouring lanes.
 *
 * A received query arrives already built.  The rank that owns it constructed its
 * ActiveData before shipping it, so the receiver re-uses that rather than
 * rebuilding one from particles it does not own -- and a query is the same
 * object to this kernel whoever sent it, which is what lets one leaf policy
 * serve every half.
 * ========================================================================== */
template <typename Spec>
struct NlrPeerAnswerDeviceFused {
    using AccumData = typename Spec::AccumData;

    static void answer(const typename Spec::DeviceContext& ctx,
                       const typename Spec::CallScalars& cs,
                       const GxDeviceTreeView& tree,
                       const std::vector<typename Spec::ActiveData>& peer_actives,
                       const std::vector<int>& peer_nodelist_flat,
                       const std::vector<int>& peer_nnodes,
                       unsigned int neighbor_type_mask,
                       std::vector<AccumData>& peer_replies_out,
                       const NlrModeDLocalSlice<Spec>& local = {})
    {
        using ActiveData  = typename Spec::ActiveData;
        using ScatterData = typename Spec::ScatterData;

        /* The Mode-D admission tests.  They used to sit on the self walk, which
         * ran on this same Spec ahead of the round loop; that walk is now this
         * one, so they live here or nowhere.  The iterative dispatchers check
         * only the search mode and the eval tier at runtime, so the bind-hook
         * test in particular has no other compile-time home. */
        static_assert(nlr_spec_modeb_eval_omp<Spec>() != ModeBEvalOMP::SerialOnly,
                      "Mode D evaluates a seeker's pairs on the rank that owns the neighbours, so a "
                      "pair kernel's neighbour-side writes land on the owner's own particles and need "
                      "no writeback -- but they land from many device lanes at once, so a kernel that "
                      "must run serially (a read-then-write of live neighbour state) cannot be served.");
        static_assert(Spec::search_mode == MODE_B_SEARCH_ONEWAY,
                      "Mode D prunes on the query's reach alone. A symmetric search also needs the "
                      "per-type supply bands, which are not mirrored to the device.");
        static_assert(!nlr_spec_has_bind_active_to_eval_context_v<Spec>,
                      "This Spec rebinds each active to the evaluating context before its pair kernel "
                      "runs, which the host backend does inside evaluate_pairs_post_drift. The fused "
                      "backend evaluates in a device kernel and performs no such rebind, so a received "
                      "active would still carry the sending rank's pointers. Port the rebind into the "
                      "device path before serving this Spec.");

        const int K = (int)peer_actives.size();
        const int n_local = (local.actives && local.accums_out) ? local.n : 0;
        const int M = n_local + K;
        if(M <= 0) {return;}

        /* The queries and their start-node lists have to be where the kernel can
         * read them.  This is a copy of the QUERIES, which is what Mode D ships
         * anyway -- not of the particles, which is what it exists to avoid.
         *
         * The LOCAL queries are not copied here: the transport normally builds
         * them on the device already, and then both launches read that buffer in
         * place.  They are staged only when it could not be had, which is the
         * same condition the transport reports rather than this guessing from a
         * pointer it cannot interrogate.  The start-node array covers the
         * RECEIVED slice only -- a local query walks from the root and has no
         * list -- so the local half costs no NODELISTLENGTH stride. */
        const bool stage_local = (n_local > 0) && !local.device_visible;
        ActiveData *ql_staged= stage_local
                                 ? (ActiveData *) nlr_shared_alloc_bytes((size_t)n_local * sizeof(ActiveData), "moded_local_q")
                                 : nullptr;
        ActiveData *q_d      = (K > 0) ? (ActiveData *) nlr_shared_alloc_bytes((size_t)K * sizeof(ActiveData), "moded_peer_q") : nullptr;
        int        *nodes_d  = (K > 0) ? (int *)        nlr_shared_alloc_bytes((peer_nodelist_flat.empty() ? 1 : peer_nodelist_flat.size()) * sizeof(int), "moded_peer_nodes") : nullptr;
        int        *nn_d     = (K > 0) ? (int *)        nlr_shared_alloc_bytes((size_t)K * sizeof(int), "moded_peer_nnodes") : nullptr;
        AccumData  *acc_d    = (AccumData *)  nlr_shared_alloc_bytes((size_t)M * sizeof(AccumData), "moded_fused_accum");
        int        *anomaly_d= (int *)        nlr_shared_alloc_bytes(sizeof(int), "moded_peer_anomaly");

        if((K > 0 && (!q_d || !nodes_d || !nn_d)) || (stage_local && !ql_staged) ||
           !acc_d || !anomaly_d) {
            /* Answering short is not an option and neither is answering on the
             * host from inside a path every rank agreed to take, so this is the
             * controlled stop.  The agreement that selected this path is what
             * removes the option of quietly doing something else here. */
            if(ql_staged) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(ql_staged);}
            if(q_d)       {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(q_d);}
            if(nodes_d)   {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(nodes_d);}
            if(nn_d)      {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(nn_d);}
            if(acc_d)     {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(acc_d);}
            if(anomaly_d) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(anomaly_d);}
            if(ThisTask == 0) {
                fprintf(stderr, "[%s] FATAL: no device memory for %d local + %d received queries "
                        "on the fused path.\n", Spec::loop_name, n_local, K);
                fflush(stderr);
            }
            endrun(90001026);
            return;
        }

        /* ⛔ THE NODE LIST IS STRIDED, NOT PACKED.  peer_nodelist_flat is
         * K * NODELISTLENGTH with each query's start nodes at its own fixed
         * offset, and peer_nnodes says how many of that query's slots are valid.
         * Treating it as packed -- cumulative offsets from the counts -- gives
         * every query after the first somebody else's start nodes, so it walks
         * the wrong subtrees and answers short.  The host walker at the same
         * data indexes [k * NODELISTLENGTH]; this has to agree with it. */
        /* Both guards trust the same thing: that the flat list really is
         * K * NODELISTLENGTH.  Clamp the per-query count to what the buffer can
         * actually hold, so a short or ragged list truncates a query instead of
         * walking off the end of nodes_d. */
        const int nodes_stride = NODELISTLENGTH;
        const size_t nodes_have = peer_nodelist_flat.size();
        for(int k = 0; k < K; k++) {
            q_d[k]  = peer_actives[k];
            int nn  = (k < (int)peer_nnodes.size()) ? peer_nnodes[k] : 0;
            if(nn > nodes_stride) {nn = nodes_stride;}
            /* Floor as well as cap.  A count off the wire is not trusted, and a
             * NEGATIVE one would otherwise reach the recording walk as this
             * batch's "no start list, walk from the root" sentinel -- so a
             * corrupt received query would have the whole local tree recorded
             * for it while the evaluating walk skipped it.  Both passes have to
             * make the same decision about the same query. */
            if(nn < 0) {nn = 0;}
            const size_t need = (size_t)(k + 1) * (size_t)nodes_stride;
            if(need > nodes_have) {nn = 0;}
            nn_d[k] = nn;
        }
        for(size_t i = 0; i < peer_nodelist_flat.size(); i++) {nodes_d[i] = peer_nodelist_flat[i];}
        if(stage_local) {
            for(int k = 0; k < n_local; k++) {ql_staged[k] = local.actives[k];}
        }
        const ActiveData *ql = ql_staged ? ql_staged : local.actives;
        *anomaly_d = 0;

        GIZMO_GPU_ENSURE_ALL_FRESH();

        /* Every query in this batch reaches this rank's particles, and they are
         * no more current for a local query than for a received one.  One
         * record-then-drift over the whole batch, each work item entering the
         * traversal the way its own kind does. */
        nlr_record_and_drift(
            ctx.P, neighbor_type_mask, tree, M, anomaly_d, "nlr_mode_d_record",
            KOKKOS_LAMBDA(int kk, double &qx, double &qy, double &qz, double &reach,
                          const int *&start_nodes, int &n_start) {
                const ActiveData& a = (kk < n_local) ? ql[kk] : q_d[kk - n_local];
                qx = (double)a.pos[0]; qy = (double)a.pos[1]; qz = (double)a.pos[2];
                reach = (double)a.h_search;
                if(kk < n_local) {
                    start_nodes = nullptr; n_start = -1;       /* ours: from the root */
                } else {
                    const int kr = kk - n_local;
                    start_nodes = nodes_d + (size_t)kr * NODELISTLENGTH;
                    n_start = nn_d[kr];
                }
            });

        const struct GxMotionTargetSet motion_targets = gx_motion_target_view();
        nlr_walk_for_sources(Spec::loop_name, M, KOKKOS_LAMBDA(int kk) {
            Spec::zero_accum(acc_d[kk]);
            const ActiveData& a = (kk < n_local) ? ql[kk] : q_d[kk - n_local];
            ScatterData s{};
            NlrModeDReduceLeaf<Spec> leaf{&ctx, &a, &acc_d[kk], &s, neighbor_type_mask, &cs};
            leaf.motion_targets = motion_targets;
            if(kk < n_local) {
                gx_device_tree_walk_from_root((double)a.pos[0], (double)a.pos[1], (double)a.pos[2],
                                              (double)a.h_search, tree, leaf, anomaly_d, neighbor_type_mask);
                return;
            }
            const int kr = kk - n_local;
            /* No start nodes: nothing on this rank was exported to this query.
             * ⚠ This is only correct because every Mode-B-wire query is
             * TARGETED: a query carrying no start nodes reached nothing on this
             * rank.  There is no broadcast query shape on this wire.  The host
             * walker skips it rather than walking from anywhere, and so does
             * this -- note it leaves the ZEROED accumulator above in place,
             * which is the right reply. */
            if(nn_d[kr] <= 0) {return;}
            gx_device_tree_walk((double)a.pos[0], (double)a.pos[1], (double)a.pos[2],
                                (double)a.h_search,
                                nodes_d + (size_t)kr * NODELISTLENGTH, nn_d[kr],
                                tree, leaf, anomaly_d, neighbor_type_mask);
        });

        const int anomaly_seen = *anomaly_d;
        /* The local half is ASSIGNED: a local query is evaluated exactly once per
         * call, and the reply merge accumulates on top of what this writes. */
        for(int k = 0; k < n_local; k++) {local.accums_out[k] = acc_d[k];}
        for(int k = 0; k < K; k++) {peer_replies_out[k] = acc_d[n_local + k];}

        if(ql_staged) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(ql_staged);}
        if(q_d)       {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(q_d);}
        if(nodes_d)   {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(nodes_d);}
        if(nn_d)      {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(nn_d);}
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(acc_d);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(anomaly_d);

        if(anomaly_seen != 0) {
            if(ThisTask == 0) {
                fprintf(stderr,
                    "[%s] FATAL: fused query walk: %s.\n", Spec::loop_name,
                    nlr_walk_anomaly_text(anomaly_seen));
                fflush(stderr);
            }
            endrun(90001027);
        }
    }
};

/* ============================================================================
 * nlr_iter_dispatch_subgroup_mode_a<Spec>(drv, sg).
 *
 * Per-iter Mode A dispatch for one subgroup. Cached oversized CSR + rebuild
 * trigger preserve legacy density_gpu.cc:140-220 session semantics.
 *
 * Lifecycle assumptions (set up by run_neighbor_loop_iterative + driver):
 *   - drv.acquire_arena_and_init_ctx_mode_a() ran ONCE before any subgroup
 *     dispatch. drv.ctx.P/CellP point at arena-resident P_gpu/CellP_gpu.
 *   - drv.arena_acquired = true; drv.ctx_initialized = true.
 *
 * Per-iter sequence:
 *   1. If !drv.mode_a_csr_valid[sg]: rebuild CSR with current compacted-
 *      active radii × Spec::mode_a_csr_buffer_factor. Frees old gnl/lookup
 *      (passing SIDX). Populates csr_offset_lookup + csr_buffered_h keyed
 *      on subgroup-slot-at-build-time.
 *   2. ghost_write_detector_begin + ghost_writeback_begin (gated on Spec
 *      traits + path → fires for Mode A iter + uses_ghost_writeback=true,
 *      e.g. ags_density wakeup; no-op for non-j-write Specs / harness).
 *   3. Stage d_actives via Spec::load_active (compacted indices).
 *   4. Pair-kernel launch over compacted active set, writing into a
 *      compacted UVM accums buffer.
 *   5. ghost_writeback_end + ghost_write_detector_end.
 *   6. Scatter compacted accums into driver-owned accum_uvm[sg][slot].
 * NO apply_active_writeback — final-only at driver level (post-iter-loop).
 *
 * CSR row-key invariant: rows in
 * cached_gnl are in COMPACTED-AT-BUILD-TIME order. csr_offset_lookup[sg]
 * maps subgroup-slot → compacted-build-time row index. Active-set
 * compaction in later iters does NOT re-key rows; converged slots simply
 * aren't walked. Rebuild creates a fresh lookup and discards the old one.
 *
 * Buffer-exceedance rebuild trigger: handled at outer-iter level AFTER
 * after_iter (post-radii mutation on AdjustRadius). Helper just consumes
 * the csr_valid[sg] flag.
 * ========================================================================== */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_a(NlrIterDriver<Spec>& drv, int sg)
{
    using ActiveData   = typename Spec::ActiveData;
    using AccumData    = typename Spec::AccumData;
    using ScatterData  = typename Spec::ScatterData;
    using NeighborData = typename Spec::NeighborData;


    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_size[sg];

    /* Collective-symmetry: ghost_writeback_begin/end MUST fire on every rank
     * for uses_ghost_writeback Specs — even
     * empty-actives ranks, otherwise reverse-comm deadlocks. Hoist the
     * dispatch hooks outside the n_compacted<=0 short-circuit. The
     * kernel-side work between begin/end is gated on n_compacted > 0. */
    std::vector<int> active_particle_indices_iter(n_compacted > 0 ? n_compacted : 0);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        active_particle_indices_iter[k] = sgr.active_indices[slot];
    }
    NeighborLoopPlan plan;
    plan.path              = NeighborLoopPlan::Path::ModeA_GpuNgl;
    plan.num_active_global = drv.global_active_total;
    neighbor_loop_args sub = drv.args;
    sub.active_list = (n_compacted > 0) ? active_particle_indices_iter.data() : nullptr;
    sub.num_active  = n_compacted;

    /* ===== (1) Build / rebuild CSR if invalid =====
     * The reimport/re-arena lifecycle lives
     * EXCLUSIVELY in the outer pre-dispatch union rebuild (run_neighbor_loop_iterative
     * step a-pre), which uses the union mask across all globally-active
     * subgroups. Inside this helper, when !mode_a_csr_valid[sg], the ghost
     * pool / arena / ctx are already correct (handled by the outer sweep);
     * we just need to free + rebuild the local CSR for this subgroup. */
    gpu_spatial_index_t *sidx = nlr_resolve_sidx_cache(Spec::sidx_cache_kind,
                                                         Spec::loop_name);
    if (!drv.mode_a_csr_valid[sg] && n_compacted > 0) {
        /* Free old cached CSR / lookup state by POINTER STATE (mode_a_csr_valid=false can mean "allocated
         * but invalid, pending rebuild"; check pointers, not flag). */
        if (drv.mode_a_cached_gnl[sg].offsets != nullptr ||
            drv.mode_a_cached_gnl[sg].neighbors != nullptr) {
            gpu_ngb_list_free(&drv.mode_a_cached_gnl[sg], sidx);
            drv.mode_a_cached_gnl[sg] = gpu_neighbor_list_t{};
        }
        if (drv.mode_a_csr_offset_lookup[sg]) {
            Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(drv.mode_a_csr_offset_lookup[sg]);
            drv.mode_a_csr_offset_lookup[sg] = nullptr;
        }
        if (drv.mode_a_csr_buffered_h[sg]) {
            Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(drv.mode_a_csr_buffered_h[sg]);
            drv.mode_a_csr_buffered_h[sg] = nullptr;
        }

        /* Build oversized CSR with CURRENT active set (compacted) and
         * oversized radii × buffer_factor. */
        std::vector<int>    active_particle_indices_host_build(n_compacted);
        std::vector<double> radii_buffered_host_build(n_compacted);
        for (int k = 0; k < n_compacted; k++) {
            int slot = drv.active_set_uvm[sg][k];
            active_particle_indices_host_build[k] = sgr.active_indices[slot];
            radii_buffered_host_build[k]          = drv.radii_uvm[sg][slot] * Spec::mode_a_csr_buffer_factor;
        }
        const size_t radii_oversized_bytes = (size_t) n_compacted * sizeof(double);
        double *radii_oversized_uvm = (double *) nlr_shared_alloc_bytes(radii_oversized_bytes, "modea_radii");
        /* A refused buffer leaves the CSR unbuilt. The rebuild is skipped rather than
         * returned from: the ghost-writeback pair below fires on every rank, so a
         * return here would leave the peers waiting on a reverse-comm that never
         * comes -- a hang, which is the one outcome worse than the stop itself. */
        if (radii_oversized_uvm != NULL) {
        for (int k = 0; k < n_compacted; k++) radii_oversized_uvm[k] = radii_buffered_host_build[k];

        /* Active-source-in-pool contract: stage explicit P[active_i].Pos for specs
         * whose active sources may be non-pool (else nullptr keeps the compact
         * fast-path). See neighbor_loop_runner.h. Radii are already explicit. */
        std::vector<double> _nlr_srcpos_storage;
        const double* _nlr_srcpos = nlr_stage_explicit_source_positions<Spec>(
            drv.ctx.P, active_particle_indices_host_build.data(), n_compacted, _nlr_srcpos_storage);
        gpu_ngb_list_build(drv.ctx.P, drv.ctx.num_total,
                           active_particle_indices_host_build.data(), n_compacted,
                           Spec::search_mode,
                           (int)sgr.j_type_bitmask,
                           &drv.mode_a_cached_gnl[sg], sidx,
                           1.0, radii_oversized_uvm, _nlr_srcpos, Spec::loop_name,
                           nlr_spec_symmetric_j_radius_scale<Spec>(),
                           Spec::radius_policy);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(radii_oversized_uvm);
        /* A build that ran out of memory hands back the empty list and has already
         * asked for the stop. Its row index is absent, which is exactly what the
         * kernel below reads first, so the CSR must not be marked valid. */
        if (drv.mode_a_cached_gnl[sg].d_active != nullptr) {

        /* Allocate csr_offset_lookup + csr_buffered_h sized to subgroup max.
         * Lookup keyed on subgroup-slot (invariant):
         *   csr_offset_lookup[slot] = build-time row index (0..n_compacted-1),
         *                              or -1 if slot wasn't in build-time set
         *                              (shouldn't happen — rebuild covers all
         *                              current actives, and converged slots
         *                              never re-enter active_set). */
        const int n_max = sgr.num_active_local;
        const size_t csr_lookup_bytes = (size_t) n_max * sizeof(int);
        const size_t csr_h_bytes      = (size_t) n_max * sizeof(double);
        drv.mode_a_csr_offset_lookup[sg] = (int *) nlr_shared_alloc_bytes(csr_lookup_bytes, "modea_csr_lookup");
        drv.mode_a_csr_buffered_h[sg]    = (double *) nlr_shared_alloc_bytes(csr_h_bytes, "modea_csr_h");
        if (!drv.mode_a_csr_offset_lookup[sg] || !drv.mode_a_csr_buffered_h[sg]) {
            if (drv.mode_a_csr_offset_lookup[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(drv.mode_a_csr_offset_lookup[sg]);
                drv.mode_a_csr_offset_lookup[sg] = nullptr;
            }
            if (drv.mode_a_csr_buffered_h[sg]) {
                Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(drv.mode_a_csr_buffered_h[sg]);
                drv.mode_a_csr_buffered_h[sg] = nullptr;
            }
            nlr_stop_no_staging_memory(Spec::loop_name, "the CSR row lookup",
                                       csr_lookup_bytes + csr_h_bytes);
        } else {
        for (int s = 0; s < n_max; s++) {
            drv.mode_a_csr_offset_lookup[sg][s] = -1;
            drv.mode_a_csr_buffered_h[sg][s]    = 0.0;
        }
        for (int k = 0; k < n_compacted; k++) {
            int slot = drv.active_set_uvm[sg][k];
            drv.mode_a_csr_offset_lookup[sg][slot] = k;
            drv.mode_a_csr_buffered_h[sg][slot]    = radii_buffered_host_build[k];
        }
        drv.mode_a_csr_valid[sg] = true;
        }   /* CSR row lookup obtained */
        }   /* neighbour list built */
        } else {
            nlr_stop_no_staging_memory(Spec::loop_name, "the oversized per-active search radii",
                                       radii_oversized_bytes);
        }
    }

    /* ===== (2) Ghost write detector + writeback begin =====
     * UNCONDITIONAL: fires on every rank for
     * uses_ghost_writeback Specs, INCLUDING empty-actives ranks, to keep
     * MPI reverse-comm collectives in lockstep with non-empty ranks.
     * Plan.path == ModeA_GpuNgl gates the dispatch hooks via
     * nlr_path_uses_imported_ghosts; uses args (post-import). */
    nlr_dispatch_ghost_write_detector_begin<Spec>(sub, plan);
    nlr_dispatch_ghost_writeback_begin<Spec>(sub, plan);

    /* The CSR is what the kernel walks, so an unbuilt one means no work this
     * iteration -- but the writeback pair above and below still runs. */
    if (n_compacted > 0 && drv.mode_a_csr_valid[sg]) {
        /* ===== (3) Stage d_actives + (4) pair-kernel launch =====
         * CSR ROW LOOKUP USAGE:
         *   slot = active_set[k]
         *   row  = csr_offset_lookup[slot]    (build-time row index)
         *   i    = d_active[row]              (particle index at build-time)
         *   h    = radii_uvm[slot]            (current radius — kernel filter)
         * CSR build uses drv.ctx.num_total (= post-import effective num_total).
         * Walk gnl.offsets[row]..offsets[row+1]. */
        /* Non-throwing on both, matching run_mode_a. The staged actives are
         * device-resident, and device memory is exhausted PER RANK rather than
         * node-wide, so a throw here would escape a single rank while its peers
         * sat in the reverse-comm collectives below -- a hang instead of a
         * reported stop. Failure is handled like a subgroup with nothing left to
         * do: the work below is skipped, the stop is requested, and control falls
         * through to the UNCONDITIONAL writeback_end, so every rank still enters
         * the collectives the same number of times. */
        ActiveData *d_actives = (ActiveData *) nlr_active_stage_alloc_bytes(
            (size_t)n_compacted * sizeof(ActiveData), "modea_active_data");
        AccumData *d_accums = (AccumData *) nlr_shared_alloc_bytes(
            (size_t)n_compacted * sizeof(AccumData), "modea_accum_data");
        if (d_actives == NULL || d_accums == NULL) {
            if (d_accums) { Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_accums); }
            nlr_active_stage_free(d_actives);
            gizmo_request_controlled_stop(7711,
                "nlr_iter_dispatch_subgroup_mode_a: per-active staging out of memory "
                "(modea_active_data is device-resident, modea_accum_data is host-visible) "
                "-- reduce the active set; note that adding ranks per node does NOT relieve "
                "device memory, since the ranks on a node share it",
                __FILE__, __LINE__, __FUNCTION__);
        } else {
            auto cs_ref = drv.cs;
            const typename Spec::DeviceContext dctx_local = drv.ctx;
            int    *active_set_arr = drv.active_set_uvm[sg];
            int    *csr_lookup     = drv.mode_a_csr_offset_lookup[sg];
            double *radii_arr      = drv.radii_uvm[sg];
            int     *d_active_arr  = drv.mode_a_cached_gnl[sg].d_active;
            int64_t *offsets       = drv.mode_a_cached_gnl[sg].offsets;
            int     *neighbors     = drv.mode_a_cached_gnl[sg].neighbors;

            const int call_slot_base = drv.call_slot_base[sg];
            gizmo_gpu_kernel_launch("nlr_iter_stage_active", n_compacted, KOKKOS_LAMBDA(int k) {
                int slot = active_set_arr[k];
                int row  = csr_lookup[slot];
                int i    = d_active_arr[row];
                double h = radii_arr[slot];
                d_actives[k] = Spec::load_active(dctx_local, call_slot_base + slot, i, h, cs_ref);
            });

            const double t_pair_kernel_start = my_second();
            const struct GxMotionTargetSet motion_targets = gx_motion_target_view();

            /* Same two assignments as the single-pass site; see the commentary
             * there. The only difference is the extra indirection from the
             * compacted active set to the build-time CSR row. */
            auto flat_kernel = KOKKOS_LAMBDA(int k) {
                /* Captured outside the `if constexpr` below; see the single-pass site. */
                const struct GxMotionTargetSet &targets = motion_targets;
                int slot = active_set_arr[k];
                int row  = csr_lookup[slot];
                Spec::zero_accum(d_accums[k]);
                const ActiveData& a = d_actives[k];
                ScatterData s{};
                int64_t start = offsets[row], end = offsets[row + 1];
                for (int64_t nn = start; nn < end; nn++) {
                    int j = neighbors[nn];
                    if constexpr (nlr_spec_writes_neighbour_motion_v<Spec>) {gx_motion_target_mark(targets, j);}
                    IdentitySidecar id{};
                    NeighborData nb = Spec::load_neighbor(dctx_local, j, id, a);
                    Spec::pair_kernel(a, nb, d_accums[k], s, cs_ref);
                }
            };

            if constexpr (nlr_mode_a_pair_policy<Spec>() == ModeAPairAssignment::TeamRowReduce) {
                using TeamKernel = NlrModeATeamPairKernel<Spec, typename Spec::DeviceContext>;
                TeamKernel fn{dctx_local, d_actives, d_accums, offsets, neighbors,
                              active_set_arr, csr_lookup, 0, cs_ref, motion_targets};
                const int team_width = nlr_mode_a_team_width<Spec>(fn);
                if (team_width > 1) {
                    gizmo_gpu_team_kernel_launch(Spec::loop_name, n_compacted, team_width, fn);
                } else {
                    gizmo_gpu_kernel_launch(Spec::loop_name, n_compacted, flat_kernel);
                }
            } else {
                gizmo_gpu_kernel_launch(Spec::loop_name, n_compacted, flat_kernel);
            }
            cpu_charge_child(CPU_PAIR_KERNEL, timediff(t_pair_kernel_start, my_second()));

            /* ===== (6) Scatter compacted accums into driver accum_uvm ===== */
            for (int k = 0; k < n_compacted; k++) {
                int slot = drv.active_set_uvm[sg][k];
                drv.accum_uvm[sg][slot] = d_accums[k];
            }

            Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_accums);
            nlr_active_stage_free(d_actives);
        }
    }

    /* ===== (5) Ghost writeback end + detector end (unconditional) ===== */
    nlr_dispatch_ghost_writeback_end<Spec>(sub, plan);
    nlr_dispatch_ghost_write_detector_end<Spec>(sub, plan);
}

template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_b_local(NlrIterDriver<Spec>& drv, int sg)
{
    using ActiveData  = typename Spec::ActiveData;
    using AccumData   = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_size[sg];
    if (n_compacted <= 0) return;     /* All actives in this subgroup converged. */

    /* Build a base sub-args restricted to the still-active slots in this subgroup.
     * This is per-iter — the active_list/num_active reflect the COMPACTED set
     * after Converged compaction from prior iters. */
    std::vector<int> active_particle_indices(n_compacted);
    std::vector<int> call_slots(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
        call_slots[k]              = drv.call_slot_base[sg] + slot;
    }
    neighbor_loop_args sub = drv.args;
    sub.active_list      = active_particle_indices.data();
    sub.active_call_slot = call_slots.data();
    sub.num_active       = n_compacted;

    /* Per-active radii in compacted order (helpers expect contiguous radii array). */
    std::vector<double> radii_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        radii_compacted[k] = drv.radii_uvm[sg][drv.active_set_uvm[sg][k]];
    }

    /* DeviceContext: driver-owned. Populated once at iter-0 by
     * NlrIterDriver::initialize_device_context_mode_b(). Per-iter reset hook
     * (if Spec declares it) fires at outer iter level, before this dispatch. */

    /* Stage actives + zero accums host-side (compacted). */
    std::vector<ActiveData> actives_compacted(n_compacted);
    build_self_actives_host_pre_drift<Spec>(sub, drv.ctx, radii_compacted.data(),
                                              drv.cs, actives_compacted.data());

    std::vector<AccumData> accums_compacted(n_compacted);
    for (int k = 0; k < n_compacted; k++) {
        Spec::zero_accum(accums_compacted[k]);     /* per-iter zero */
    }

    /* Collect → drift → evaluate (same helper chain as run_mode_b_local).
     * Per-subgroup mask threaded via sgr.j_type_bitmask
     * for multi-subgroup walker support. */
    std::vector<std::vector<int>> cand_modeB;
    collect_candidates_pre_drift<Spec>(sub, radii_compacted.data(),
                                         (unsigned int)sgr.j_type_bitmask,
                                         DispatchPath::ModeB_HostWalker, cand_modeB);
    drift_and_filter_candidates<Spec>(cand_modeB, (unsigned int)sgr.j_type_bitmask, [&](int aa, double *pos, double &h_q) {
        const int i = sub.active_list[aa];
        pos[0] = (double)sub.P[i].Pos[0]; pos[1] = (double)sub.P[i].Pos[1]; pos[2] = (double)sub.P[i].Pos[2];
        h_q = radii_compacted[aa];
    });
    evaluate_pairs_post_drift<Spec>(drv.ctx, actives_compacted.data(), n_compacted,
                                      cand_modeB, accums_compacted.data(), drv.cs, EvalOMPPolicy::AllowProduction);

    /* Scatter compacted accums back into driver-owned per-slot accum_uvm.
     * Slots NOT in active_set_uvm keep their stale values (will not be
     * read this iter — after_iter only fires on still-active slots). */
    for (int k = 0; k < n_compacted; k++) {
        int slot = drv.active_set_uvm[sg][k];
        drv.accum_uvm[sg][slot] = accums_compacted[k];
    }
}


/* ============================================================================
 * Default on_max_iter_exceeded — runner-supplied bad-stop on max iteration.
 * endrun(1155) is now a soft controlled-stop (post Stage-1 macro flip): the
 * run flags + proceeds with un-converged radii to the next phase-boundary
 * poll, then finalizes cleanly (no MPI_Abort). Matches legacy
 * hydro/density.cc:602 and gravity/ags_rkern.cc:414 stop policy.
 * Specs override via `static void on_max_iter_exceeded(const NlrIterDriver<Spec>&);`
 * only when port's legacy policy differs and Phil approved (rare).
 * ========================================================================== */
template <typename Spec>
static void nlr_default_on_max_iter_exceeded(const NlrIterDriver<Spec>& drv)
{
    if (ThisTask == 0) {
        fprintf(stderr,
            "[%s] FATAL: failed to converge in %d iterations (max_iters=%d). "
            "global_active_total=%d remained.\n",
            Spec::loop_name, drv.iter_index, Spec::max_iters, drv.global_active_total);
        fflush(stderr);
    }
    endrun(1155);
}

/* ============================================================================
 * nlr_iter_dispatch_subgroup_mode_d<Spec>(drv, sg) — one iteration, one subgroup.
 *
 * Two shapes, chosen by whether this rank has peers.
 *
 * With peers, the queries travel: this rank's own go out, other ranks' come in,
 * and both halves are answered by walking the tree of whichever rank owns the
 * particles.  That runs through the same transport Mode B uses -- same wire,
 * same bounded rounds, same reply merge -- with the fused backend supplying the
 * answers instead of a host walk over a collected list.
 *
 * Alone, there is nothing to exchange, so the transport is skipped entirely and
 * the queries never leave: they are built in the kernel from resident particles
 * rather than up front on the host.
 *
 * There is no ghost pool and no neighbour list on either shape, and their
 * absence is the whole of the difference from Mode A: the same Spec, the same
 * pair kernel, the same accumulator, reached by walking rather than by
 * importing. Staging is reduced rather than gone -- the queries are built where
 * the particles live instead of on the host, but their indices and radii, and
 * the accumulators coming back, are still staged.
 *
 * WHERE DECLINING IS DECIDED.  Whether this rank can answer on the device is
 * settled COLLECTIVELY before dispatch, when the path is chosen, because the
 * three paths enter different collectives and a rank deciding alone would hang
 * rather than answer differently.  The one rank-local fallback left in this
 * function is for a device allocation that fails, and it is reachable only on
 * the single-rank shape, which issues nothing collective -- with peers the
 * function has already returned through the transport before that point.
 * ========================================================================== */
template <typename Spec>
static void nlr_iter_dispatch_subgroup_mode_d(NlrIterDriver<Spec>& drv, int sg)
{
    using AccumData = typename Spec::AccumData;

    const NlrSubgroup& sgr = drv.args.subgroups[sg];
    const int n_compacted  = drv.active_set_size[sg];

    /* ⛔ NO EARLY RETURN ON AN EMPTY LOCAL SET WHEN THERE ARE PEERS.  The peer
     * path below is collective: a rank with nothing of its own still has to
     * enter it, because its peers are waiting to send it queries and to receive
     * their replies.  The local-only path further down may return early, and
     * does, because it issues nothing collective at all. */

    /* The tree was described, and the rank made current, once for the whole
     * call before any discovery round -- not here.  Re-preparing per subgroup
     * per iteration would repeat a decision the path selection has already acted
     * on, and there is nothing to re-establish: nothing this path touches can go
     * stale within the call. */
    const GxDeviceTreeView &tree = drv.mode_d_tree;

    /* Compact the still-active slots, as every other per-iter dispatch does. */
    std::vector<int>    active_particle_indices(n_compacted);
    std::vector<int>    call_slots(n_compacted);
    std::vector<double> radii_compacted(n_compacted);
    for(int k = 0; k < n_compacted; k++) {
        const int slot = drv.active_set_uvm[sg][k];
        active_particle_indices[k] = sgr.active_indices[slot];
        call_slots[k]              = drv.call_slot_base[sg] + slot;
        radii_compacted[k]         = drv.radii_uvm[sg][slot];
    }
    neighbor_loop_args sub = drv.args;
    sub.active_list      = active_particle_indices.data();
    sub.active_call_slot = call_slots.data();
    sub.num_active       = n_compacted;

    /* With peers, the queries travel and the owners answer them: same wire, same
     * reply merge, same bounded rounds as the host walker uses -- only the way a
     * query becomes an accumulator differs, on both halves.  Entered on every
     * rank, including one with no actives of its own. */
    if(NTask > 1) {
        std::vector<AccumData> accums_compacted(n_compacted);
        for(int k = 0; k < n_compacted; k++) {Spec::zero_accum(accums_compacted[k]);}
        mode_b_remote_evaluate_into_buffer<Spec, NlrEvalBackend::DeviceFused>(
            sub, radii_compacted.data(), drv.cs, drv.ctx,
            (unsigned int)sgr.j_type_bitmask,
            (n_compacted > 0) ? accums_compacted.data() : nullptr,
            &drv.mode_d_tree);
        for(int k = 0; k < n_compacted; k++) {
            drv.accum_uvm[sg][drv.active_set_uvm[sg][k]] = accums_compacted[k];
        }
        return;
    }

    /* Single rank: nothing to exchange, so the queries never leave and are built
     * in the kernel from resident particles rather than up front on the host. */
    if(n_compacted <= 0) {return;}

    std::vector<AccumData> accums_compacted(n_compacted);
    const bool evaluated = nlr_mode_d_evaluate_single_rank<Spec>(
        drv.ctx, drv.cs, active_particle_indices.data(), call_slots.data(),
        radii_compacted.data(), n_compacted,
        /* The SUBGROUP's mask, as the peer branch and every other per-subgroup
         * dispatcher use.  A multi-subgroup Spec would otherwise see different
         * neighbours on one rank than on many. */
        (unsigned int)sgr.j_type_bitmask,
        tree, accums_compacted.data());
    if(!evaluated) {
        /* Out of device memory is a reason to answer differently, never a reason
         * to answer less, so the host walker takes this iteration.
         *
         * Reachable on the single-rank shape only -- the peer shape returned
         * above -- which is what keeps this legal: swapping backends mid-call is
         * invisible to peers, and with peers it would be a hang rather than a
         * different answer.  If this branch is ever given a path that has peers,
         * the decision has to move to the collective site where readiness is
         * already decided. */
        nlr_iter_dispatch_subgroup_mode_b_local<Spec>(drv, sg);
        return;
    }
    for(int k = 0; k < n_compacted; k++) {
        drv.accum_uvm[sg][drv.active_set_uvm[sg][k]] = accums_compacted[k];
    }
}

template <typename Spec>
static bool nlr_mode_d_evaluate_single_rank(const typename Spec::DeviceContext& ctx,
                                            const typename Spec::CallScalars& cs,
                                            const int *active_idx_host,
                                            const int *active_slot_host,
                                            const double *radii_host,
                                            int n,
                                            unsigned int supply_mask,
                                            const GxDeviceTreeView& tree,
                                            typename Spec::AccumData *accums_out)
{
    using AccumData = typename Spec::AccumData;
    if(n <= 0) {return true;}

    int       *idx_d     = (int *)    nlr_shared_alloc_bytes((size_t)n * sizeof(int),    "moded_idx");
    int       *slot_d    = (int *)    nlr_shared_alloc_bytes((size_t)n * sizeof(int),    "moded_slot");
    double    *rad_d     = (double *) nlr_shared_alloc_bytes((size_t)n * sizeof(double), "moded_radii");
    AccumData *acc_d     = (AccumData *) nlr_shared_alloc_bytes((size_t)n * sizeof(AccumData), "moded_accum");
    int       *anomaly_d = (int *)    nlr_shared_alloc_bytes(sizeof(int), "moded_anomaly");

    if(!idx_d || !slot_d || !rad_d || !acc_d || !anomaly_d) {
        if(idx_d)     {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(idx_d);}
        if(slot_d)    {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(slot_d);}
        if(rad_d)     {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(rad_d);}
        if(acc_d)     {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(acc_d);}
        if(anomaly_d) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(anomaly_d);}
        return false;
    }

    for(int k = 0; k < n; k++) {
        idx_d[k]  = active_idx_host[k];
        slot_d[k] = active_slot_host[k];
        rad_d[k]  = radii_host[k];
    }
    *anomaly_d = 0;

    nlr_mode_d_local_reduce<Spec>(ctx, cs, idx_d, slot_d, rad_d, n, supply_mask, tree, acc_d, anomaly_d);

    const int anomaly_seen = *anomaly_d;
    for(int k = 0; k < n; k++) {accums_out[k] = acc_d[k];}

    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(idx_d);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(slot_d);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(rad_d);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(acc_d);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(anomaly_d);

    /* The traversal reports the one state the host walk treats as fatal: an
     * index belonging to none of the three classes, which means the tree is
     * malformed.  The host stops the run there and so does this, rather than
     * keep an answer short by an unknown amount.  Requested after the buffers
     * are released, because the stop drains at the next phase boundary rather
     * than here. */
    if(anomaly_seen != 0) {
        if(ThisTask == 0) {
            fprintf(stderr,
                "[%s] FATAL: device walk: %s.\n", Spec::loop_name,
                nlr_walk_anomaly_text(anomaly_seen));
            fflush(stderr);
        }
        endrun(90001025);
    }
    return true;
}

/* ============================================================================
 * run_neighbor_loop_iterative<Spec> — body (step 2b).
 * ========================================================================== */
template <typename Spec>
void run_neighbor_loop_iterative(const neighbor_loop_args_iterative& args_in)
{
    /* The one live argument view for this call.  A ghost import can raise the particle
     * capacity and move P[]/CellP[], so the cached base pointers are refreshed in this
     * object as they change; the driver holds a reference to it, so the runner and every
     * hook read the same, current values.  The caller's own copy is left untouched. */
    neighbor_loop_args_iterative args = args_in;

    /* ===== Compile-time spec consistency ===== */
    static_assert(std::is_same_v<typename Spec::IterControl, Iterative>,
                  "run_neighbor_loop_iterative requires Spec::IterControl = Iterative. "
                  "NotIterative Specs must call run_neighbor_loop instead.");
    static_assert(nlr_spec_has_after_iter_v<Spec>,
                  "Iterative Spec is missing `static IterResult after_iter(const AfterIterContext<Spec>&, const AccumData&);`. "
                  "Required when IterControl = Iterative.");
    static_assert(std::is_trivially_copyable_v<typename Spec::IterScratch>,
                  "Spec::IterScratch must be std::is_trivially_copyable_v. "
                  "Use `using IterScratch = NoIterScratch;` for the empty form.");
    static_assert(Spec::max_iters >= 1,
                  "Spec::max_iters must be >= 1.");
    /* Active-source-in-pool contract (see neighbor_loop_runner.h). */
    static_assert(nlr_spec_satisfies_source_pool_contract_v<Spec>,
        "Cached-SIDX Spec must declare 'static constexpr bool mode_a_active_sources_in_sidx_pool' "
        "(true = active sources are SIDX-pool members; false = runner stages explicit P[].Pos). "
        "Prevents the stale gas-only-compact source-position bug for non-pool actives.");
    static_assert(Spec::mode_a_csr_buffer_factor >= 1.0,
                  "Spec::mode_a_csr_buffer_factor must be >= 1.0. It oversizes the CSR so an "
                  "iteration that ENLARGES the search radius can reuse the list instead of "
                  "rebuilding it (legacy DENSITY_H_BUFFER_FACTOR = 1.3). A Spec whose radius is "
                  "fixed for the whole call -- no AdjustRadius -- declares 1.0: a larger value "
                  "there is not a buffer, it is a wider import nothing will ever read.");
    /* TRAP-5 carry-forward: same trivially-copyable
     * checks as run_neighbor_loop. Don't let iterative Specs bypass TRAP 5. */
    static_assert(std::is_trivially_copyable_v<typename Spec::CallScalars>,
                  "Spec::CallScalars must be trivially copyable (lambda-captured by value).");
    static_assert(std::is_trivially_copyable_v<typename Spec::ActiveData>,
                  "Spec::ActiveData must be trivially copyable (UVM-staged).");
    static_assert(std::is_trivially_copyable_v<typename Spec::NeighborData>,
                  "Spec::NeighborData must be trivially copyable (built per-pair on device).");
    static_assert(std::is_trivially_copyable_v<typename Spec::AccumData>,
                  "Spec::AccumData must be trivially copyable (UVM-staged).");

    /* ===== Runtime checks BEFORE any state touch ===== */
    if (args.num_subgroups < 1) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[run_neighbor_loop_iterative<%s>] FATAL: num_subgroups=%d < 1. "
                "Caller must short-circuit when global active total is 0.\n",
                Spec::loop_name, args.num_subgroups);
            fflush(stderr);
        }
        endrun(81200);
        return;   /* before any state touch; symmetric caller-contract failure -> graceful return, drains at next poll */
    }
    if (args.num_subgroups > 1 && !nlr_supports_subgroups<Spec>::value) {
        if (ThisTask == 0) {
            fprintf(stderr,
                "[run_neighbor_loop_iterative<%s>] FATAL: num_subgroups=%d > 1 "
                "but Spec did not declare `using SupportsSubgroups = std::true_type;` (TRAP 9).\n",
                Spec::loop_name, args.num_subgroups);
            fflush(stderr);
        }
        endrun(81201);
        return;   /* before any state touch; symmetric caller-contract failure -> graceful return, drains at next poll */
    }
    /* Multi-subgroup walker:
     *   - sg.j_type_bitmask threaded through Mode B local + remote collectors
     *     and Mode A NGL build per-subgroup.
     *   - Mode A multi-subgroup ghost import uses union semantics (mask
     *     OR + concatenation of actives/radii across subgroups). Documented
     *     as deliberate over-import for cross-bm cases.
     *   - Per-subgroup global activity tracking via global_active_per_sg[].
     *   - Pre-dispatch invalidation sweep rebuilds CSR caches once per iter
     *     when any subgroup invalidated.
     *   - The caller must fill subgroups[] from the global_bm_presence union
     *     with identical ordering on every rank, and must place each particle
     *     in at most one subgroup.
     * Per-Spec opt-in still REQUIRED via `using SupportsSubgroups = std::true_type;`
     * (runtime abort 81201 above catches missing trait). */


    /* ===== Path selection at iter 0 (FIXED for whole call) =====
     * Integrates with the canonical override / threshold dispatch
     * (mirrors the selection in run_neighbor_loop). Path is held fixed
     * across all iterations of one iterative call.
     *
     * num_active for threshold uses the UNION across all subgroups (base
     * args.num_active per the doc convention). This selection is single-subgroup-only;
     * multi-subgroup walker mask-threading is handled separately. */
    /* Dispatch priority: args.dispatch_override > adaptive threshold
     * (mirrors the non-iterative site). Iterative Specs (density, ags_density,
     * mechfb, etc.) generally won't set dispatch_override since the
     * corridor mode flows through cellcorrections/gradients/hydro_force; preserved
     * here for completeness and future use. */
    const NlrForceMode force_mode = args.dispatch_override;
    DispatchPath path;
    /* Filled by the Mode D preparation below if that path is taken; handed to
     * the driver so every iteration walks the tree the decision was made on. */
    GxDeviceTreeView mode_d_tree{};
    int forced_modeb_global_active = -1;
    /* Global active-particle sum across all ranks (from the dispatch Allreduce);
     * -1 = not computed (force-A cheap path). Used for the globally-zero-active
     * no-op below. */
    int global_active_sum = -1;
    if (force_mode == NlrForceMode::A) {
        path = DispatchPath::ModeA_GPU_NGL;
    } else if (force_mode == NlrForceMode::B) {
        int local_act = args.num_active;
        MPI_Allreduce(&local_act, &forced_modeb_global_active, 1,
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        global_active_sum = forced_modeb_global_active;
        path = DispatchPath::ModeB_HostWalker;
    } else {
        /* Threshold dispatch: Allreduce sum + max of base args.num_active
         * (= union across subgroups). */
        int local_act = args.num_active;
        int sum_act = 0, max_act = 0;
        MPI_Allreduce(&local_act, &sum_act, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_act, &max_act, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        global_active_sum = sum_act;
        const int spec_default_sum = nlr_spec_threshold_sum<Spec>(64);
        const int spec_default_max = nlr_spec_threshold_max<Spec>(64);
        const int TS = gizmo_nlr_modeb_threshold_sum_for(Spec::loop_name, spec_default_sum);
        const int TM = gizmo_nlr_modeb_threshold_max_for(Spec::loop_name, spec_default_max);
        bool select_mode_b = (sum_act > 0) && (sum_act <= TS) && (max_act <= TM);
        path = select_mode_b ? DispatchPath::ModeB_HostWalker : DispatchPath::ModeA_GPU_NGL;
        /* Answered where its neighbours live, never on ghost copies (stated at
         * the single-pass entry, nlr_spec_needs_live_neighbours). */
        if constexpr (nlr_spec_needs_live_neighbours_v<Spec>) {
            path = DispatchPath::ModeB_HostWalker;
        }
#ifdef NEIGHBOR_LOOP_MODE_D
        /* Mode D takes the calls Mode A would have taken, when the loop is one
         * it can serve and every rank can actually answer on the device.  It
         * serves every one-way search whose pair kernel may run concurrently:
         * the neighbours are always the evaluating rank's own particles, so a
         * kernel's neighbour-side writes are local atomics, exactly as they are
         * on the host walker's remote protocol, which this is the device form of.
         *
         * The verdict has to be UNANIMOUS, not merely identical-looking.  Two of
         * the terms are compile-time Spec traits and one is a path already chosen
         * from a collective, so those cannot differ -- but whether a rank's tree
         * can be described to the device is a local fact that genuinely can
         * differ between ranks, and the three paths enter different collectives.
         * A rank deciding alone would not answer differently, it would hang.
         *
         * The preparation runs HERE, once, before any discovery round, and what
         * it establishes holds for every iteration of the call: no particles are
         * imported on this path and the writeback runs after the last iteration,
         * so nothing it certifies can go stale in between.
         *
         * With the peer half in place this now runs at any task count, and the
         * reduction above stops being a formality: a rank that cannot describe
         * its tree to the device pulls the WHOLE call back to Mode A, because the
         * alternative is one rank walking while its peers wait for queries that
         * never arrive. */
        /* sum_act > 0 FIRST, and it is not a shortcut: preparing mirrors and
         * sweeps the whole tree, and a call with no actives anywhere has no work
         * to justify that.  Without this test such a call reaches here because
         * select_mode_b is false when sum_act == 0, takes the Mode-A label, and
         * pays that sweep on the way to the globally-zero-active no-op that was
         * going to return without computing anything.  That is O(Nnodes) on a
         * call that should cost nothing, appearing only in builds carrying the
         * flag.  It used to also pay a full-rank particle drift here; that is
         * what the touched-set discovery removed, and the node half is what
         * still justifies the gate. */
        if(sum_act > 0 &&
           (path == DispatchPath::ModeA_GPU_NGL ||
            (nlr_spec_needs_live_neighbours_v<Spec> && !select_mode_b)) &&
           nlr_spec_modeb_eval_omp<Spec>() != ModeBEvalOMP::SerialOnly &&
           Spec::search_mode == MODE_B_SEARCH_ONEWAY) {
            const int ready_local = (gx_device_fused_walk_prepare(&mode_d_tree, Spec::loop_name) == 0) ? 1 : 0;
            int ready = 0;
            MPI_Allreduce(&ready_local, &ready, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            if(ready) {path = DispatchPath::ModeD_DeviceFused;}
        }
#endif
    }

    /* Globally-zero-active call: do NO neighbor work. global_active_sum comes
     * from the dispatch Allreduce, so it is identical on every rank -> this
     * return is collective-symmetric (all ranks return together, skipping the
     * ghost import / writeback / cleanup as a matched set). This is NOT the
     * banned local num_active==0 early return: the condition is GLOBAL. Without
     * it, a zero-active call falls to Mode A (sum_act>0 gate fails) and fires a
     * spurious request-driven ghost import with nothing to compute. Not fired on
     * the force-A path (global_active_sum stays -1 there). The normal
     * PHASE0/dispatch summary is intentionally skipped for such calls; a distinct
     * rank-0 marker (diag-gated) keeps them observable. */
    if (global_active_sum == 0) {
        return;
    }


    /* ===== CallScalars captured ONCE for whole call ===== */
    typename Spec::CallScalars cs = Spec::populate_call_scalars(args);

    /* ===== Mode B hard-corridor counter snapshot =====
     * Mode B paths MUST NOT enter move_particles / ghost_exchange_impl /
     * gpu_particles_arena_acquire, and MUST NOT mutate NumPart. Always-on
     * enforcement around the iterative call; mirrors non-iter run_neighbor_loop
     * lines 2018-2026 + 2099-2122. Snapshot here; check post-iter-loop.
     * Mode A paths legitimately advance these counters (one ghost_import +
     * one arena_acquire per call); enforcement skipped. */
    const uint64_t s_drift0 = g_global_drift_counter;
    const uint64_t s_ghost0 = g_ghost_import_counter;
    const uint64_t s_arena0 = g_gpu_arena_acquire_counter;
    const int      s_np0    = NumPart;

    /* ===== Driver lifetime wrapped in inner scope =====
     * Driver destructor runs at scope exit; Mode B hard-corridor check fires
     * AFTER scope exit so the check covers driver cleanup hooks too (defensive
     * — current cleanup_device_context contractually MUST NOT touch
     * drift/ghost/arena globals, but inner-scope wrapping makes the invariant
     * maximally airtight). */
    {
    /* Declared before the driver so it closes after the driver has finished
     * (destructors run in reverse): the raise sees the final velocities. */
    NlrMotionTargetScope<Spec> motion_target_scope;
    NlrIterDriver<Spec> drv(args, cs);

    /* ===== Path-specific DeviceContext init =====
     * Mode B: bind ctx.P=args.P/CellP directly (no arena).
     * Mode A: acquire_arena_and_init_ctx_mode_a — arena_acquire ONCE per call
     *         then bind ctx to arena-resident P_gpu/CellP_gpu. */
    if (path == DispatchPath::ModeB_HostWalker) {
        drv.initialize_device_context_mode_b();
    } else if (path == DispatchPath::ModeD_DeviceFused) {
        drv.mode_d_tree = mode_d_tree;
        drv.acquire_arena_and_init_ctx_mode_d();
    } else if (path == DispatchPath::ModeA_GPU_NGL) {
        /* Both read P/CellP from the device-resident arena; Mode D differs in
         * how it finds the neighbours, not in where the particles live. */
        drv.acquire_arena_and_init_ctx_mode_a();
    } else {
        /* Future path — not currently reachable. */
        if (ThisTask == 0) {
            fprintf(stderr,
                "[run_neighbor_loop_iterative<%s>] FATAL: unhandled path %d.\n",
                Spec::loop_name, (int)path);
            fflush(stderr);
        }
        /* path is Allreduce-symmetric (see dispatch above) -> all ranks hit this
         * together; soft bad-stop + immediate poll (81203 precedent). */
        endrun(81208); gizmo_exit_bad_stop_if_requested("nlr:unhandled_dispatch");
    }
    /* Drain a soft bad-stop raised inside the path-specific init (Mode-A arena
     * lifecycle 81211/81212 return early without binding ctx) BEFORE any device
     * dispatch. All-rank: Mode-A path is symmetric. */
    gizmo_exit_bad_stop_if_requested("nlr:iter_context_init");

    /* ===== Outer iter loop =====
     *
     * INVARIANT: `accum_uvm[sg][slot]` always
     * holds the most recent evaluated AccumData for that slot. Converged
     * slots are NEVER re-touched after they leave active_set; their final
     * iter's result persists in accum_uvm until the post-loop apply_active_writeback
     * reads it.
     *
     * Per-iter zero of AccumData is therefore the dispatch helper's
     * responsibility (it zeroes via Spec::zero_accum on ONLY the compacted
     * active subset before pair_kernel evaluation; see
     * nlr_iter_dispatch_subgroup_mode_b_local).
     * The runner does NOT zero converged slots — doing so would destroy
     * exactly the result apply_active_writeback needs. The earlier "defensive
     * cross-iter zero of all slots" in this position was a bug. */
    /* Partition-by-subgroup contract: record the
     * iter-0 per-subgroup j_type_bitmask values so the per-iter partition
     * assertion below can verify Spec::active_subgroup_key returns the same
     * key for every active in that subgroup on every iter. Only populated
     * when the Spec opts into actives_partition_by_subgroup. Compile-time
     * no-op for Specs that don't. */
    std::vector<unsigned int> partition_expected_bm_key;
    if constexpr (nlr_spec_actives_partition_by_subgroup_v<Spec>) {
        static_assert(nlr_spec_has_active_subgroup_key_v<Spec>,
                      "Spec::actives_partition_by_subgroup=true requires "
                      "Spec::active_subgroup_key(const DeviceContext&, int, "
                      "const CallScalars&) returning int (= bm key).");
        partition_expected_bm_key.resize(args.num_subgroups);
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            partition_expected_bm_key[sg] = args.subgroups[sg].j_type_bitmask;
        }
    }

    for (drv.iter_index = 0; drv.iter_index < Spec::max_iters; drv.iter_index++) {
        if(!nlr_args_view_is_live(args)) {
            if(ThisTask == 0) {
                fprintf(stderr, "[run_neighbor_loop_iterative<%s>] FATAL: the argument view no longer "
                        "describes the live particle storage at iter=%d (P=%p vs %p, num_total=%d vs %d). "
                        "A capacity change during the call was not reflected here.\n",
                        Spec::loop_name, drv.iter_index, (void *) args.P, (void *) P,
                        args.num_total, NumPart);
                fflush(stderr);
            }
            /* Request the stop and stay in lockstep: this loop body holds MPI collectives,
             * so breaking out on one rank alone would hang its peers.  The stop drains at
             * the next phase boundary. */
            endrun(81214);
        }

        /* mode_a_rebuild_csr_every_iter correctness
         * fallback. When the Spec sets
         * this trait true, force-invalidate every subgroup's Mode A CSR
         * cache at the start of every iter > 0, bypassing the buffer-
         * exceedance trigger entirely. The static_assert on
         * mode_a_csr_buffer_factor > 1.0 is unchanged (factor stays a valid
         * number even when unused). Iter 0 is exempt — initial CSR build
         * happens in the per-subgroup dispatch on first invalidation pass. */
        if constexpr (nlr_spec_mode_a_rebuild_csr_every_iter_v<Spec>) {
            if (path == DispatchPath::ModeA_GPU_NGL && drv.iter_index > 0) {
                for (int sg = 0; sg < args.num_subgroups; sg++) {
                    drv.mode_a_csr_valid[sg] = false;
                }
            }
        }

        /* Partition-by-subgroup debug assertion (DEBUG / GIZMO_NLR_ASSERT_PARTITION
         * only; zero production overhead). Checks each active's active_subgroup_key
         * against the iter-0 j_type_bitmask; a mismatch routes to the graceful
         * controlled stop below. Host-side, before per-iter dispatch. */
#if defined(DEBUG) || defined(GIZMO_NLR_ASSERT_PARTITION)
        if constexpr (nlr_spec_actives_partition_by_subgroup_v<Spec>) {
            int local_partition_bad = 0;
            for (int sg = 0; sg < args.num_subgroups && !local_partition_bad; sg++) {
                /* Skip globally-converged subgroups (iter > 0 only). */
                if (drv.iter_index > 0 && drv.global_active_per_sg[sg] <= 0) continue;
                const int n_compacted = drv.active_set_size[sg];
                const NlrSubgroup& sgr = args.subgroups[sg];
                const unsigned int expect = partition_expected_bm_key[sg];
                for (int k = 0; k < n_compacted; k++) {
                    const int slot = drv.active_set_uvm[sg][k];
                    const int i    = sgr.active_indices[slot];
                    const int key  = Spec::active_subgroup_key(drv.ctx, i, drv.cs);
                    if (static_cast<unsigned int>(key) != expect) {
                        printf("[run_neighbor_loop_iterative<%s>] FATAL: "
                               "actives_partition_by_subgroup key drift. "
                               "sg=%d iter=%d slot=%d i=%d expected_bm=%u got_key=%d. "
                               "active_subgroup_key MUST be a pure function of "
                               "state that does not change across iterations.\n",
                               Spec::loop_name, sg, drv.iter_index, slot, i,
                               expect, key); fflush(stdout);
                        local_partition_bad = 1;
                        break;
                    }
                }
            }
            /* Debug-only: collectivize the per-rank assertion result so all ranks
             * reach the controlled stop together BEFORE the per-iter dispatch
             * collectives below (a lone asserting rank would otherwise desync). */
            int global_partition_bad = local_partition_bad;
            if (NTask > 1) {
                MPI_Allreduce(&local_partition_bad, &global_partition_bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (global_partition_bad) {
                endrun(90001020);
                gizmo_exit_bad_stop_if_requested("nlr:partition_key_drift");
            }
        }
#endif

        /* (a-pre) Mode A pre-dispatch invalidation sweep.
         * If any subgroup's CSR is invalid (set by the buffer-exceedance
         * trigger in last iter's step b.5), call the no-arg union-rebuild
         * method ONCE before any per-subgroup dispatch this iter. Prevents
         * mixing old/new arena pools across subgroup dispatches in the same
         * iter. Skips iter 0 (initial arena
         * acquire was via acquire_arena_and_init_ctx_mode_a). */
        if (path == DispatchPath::ModeA_GPU_NGL && drv.iter_index > 0) {
            /* The sweep must be COLLECTIVE.
             * mode_a_csr_valid is rank-local; a rank with no actives for a
             * globally-active subgroup never builds CSR locally, so its
             * mode_a_csr_valid[sg] stays false while another rank's is true.
             * Without an Allreduce, one rank enters ghost_exchange_cleanup +
             * reimport collectives while the other skips them => deadlock.
             * Also skip globally-converged subgroups: a sg with
             * global_active_per_sg[sg]==0 may have invalid CSR from an
             * earlier buffer-exceedance trigger but doesn't need a rebuild. */
            int local_needs_rebuild = 0;
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                if (drv.global_active_per_sg[sg] <= 0) continue;
                /* Note: a rank
                 * with zero local actives in this globally-active subgroup
                 * never builds CSR locally, so its mode_a_csr_valid[sg]
                 * stays false. Without this guard, such a rank would vote
                 * "rebuild" every iteration forever — Allreduce would
                 * force a global rebuild that doesn't actually fix the
                 * vote, infinite loop of useless global rebuilds. The rank
                 * still participates collectively in the Allreduce below
                 * (voting 0); it just doesn't request a rebuild for a CSR
                 * it doesn't need. */
                if (drv.active_set_size[sg] <= 0) continue;
                if (!drv.mode_a_csr_valid[sg]) {
                    local_needs_rebuild = 1;
                    break;
                }
            }
            int global_needs_rebuild = local_needs_rebuild;
            if (NTask > 1) {
                MPI_Allreduce(&local_needs_rebuild, &global_needs_rebuild, 1,
                              MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            }
            if (global_needs_rebuild) {
                drv.rebuild_mode_a_arena_and_ctx_for_current_active_union();
            }
        }

        /* (a0) Per-iter Spec hook: reset_per_iter_device_context.
         * Optional. Runs HOST-side on every rank before subgroup dispatch.
         * Use case: ags_density's per_iter_wakeup_detected counter zero.
         * Hook MUST NOT do MPI. */
        if constexpr (nlr_spec_has_reset_per_iter_device_context_v<Spec>) {
            Spec::reset_per_iter_device_context(args, drv.ctx, drv.iter_index);
        }

        /* (a) Per-subgroup dispatch — fixed path for the call. This
         * implements Mode B local (np=1) and Mode B remote (np>1) via the
         * same DispatchPath::ModeB_HostWalker label; the per-subgroup helper
         * picks local vs remote based on NTask. Mode A iter is still
         * hard-stubbed at the path-selection block. */
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            /* Skip globally-converged subgroups. Saves
             * per-iter no-op collectives. On iter 0 global_active_per_sg
             * is still 0 (not yet Allreduced); use local count as proxy +
             * global presence-of-actives (subgroups[] is filled from
             * global_bm_presence by the caller, so any subgroup present
             * has SOME rank with non-zero actives — entering the dispatch
             * is the correct lockstep collective participation). */
            const bool sg_globally_active = (drv.iter_index == 0)
                ? true   /* iter 0: trust caller's subgroups[] (global union) */
                : (drv.global_active_per_sg[sg] > 0);
            if (!sg_globally_active) continue;

            /* (a) Mode B only (Mode A hard-stubbed at outer entry).
             * Helpers are collective on remote — must be
             * entered on every rank regardless of local n_compacted. */
            switch (path) {
                case DispatchPath::ModeB_HostWalker:
                    if (NTask == 1) {
                        nlr_iter_dispatch_subgroup_mode_b_local<Spec>(drv, sg);
                    } else {
                        nlr_iter_dispatch_subgroup_mode_b_remote<Spec>(drv, sg);
                    }
                    break;
                case DispatchPath::ModeA_GPU_NGL:
                    nlr_iter_dispatch_subgroup_mode_a<Spec>(drv, sg);
                    break;
                case DispatchPath::ModeD_DeviceFused:
                    if constexpr (nlr_spec_modeb_eval_omp<Spec>() != ModeBEvalOMP::SerialOnly &&
                                  Spec::search_mode == MODE_B_SEARCH_ONEWAY) {
                        nlr_iter_dispatch_subgroup_mode_d<Spec>(drv, sg);
                    }
                    break;
            }
        }

        /* (b) Per-active Spec::after_iter — collect IterStatus, compact
         * active_set per subgroup, mutate radii on AdjustRadius. */
        drv.local_active_total = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            int n_compacted   = drv.active_set_size[sg];
            int write_idx     = 0;
            for (int k = 0; k < n_compacted; k++) {
                int  slot = drv.active_set_uvm[sg][k];
                int  i    = args.subgroups[sg].active_indices[slot];
                AfterIterContext<Spec> ctx{
                    args, sg, slot, i, drv.iter_index,
                    /* h_search_current */ drv.radii_uvm[sg][slot],
                    /* scalars (ref to driver-owned) */ drv.cs,
                    /* scratch (mutable ref) */          drv.scratch_uvm[sg][slot]
                };
                IterResult r = Spec::after_iter(ctx, drv.accum_uvm[sg][slot]);
                switch (r.status) {
                    case IterStatus::Converged:
                        /* Drop slot from compacted active_set; its accum_uvm
                         * stays so apply_active_writeback can read it post-loop. */
                        break;
                    case IterStatus::NeedsMore:
                        drv.active_set_uvm[sg][write_idx++] = slot;
                        break;
                    case IterStatus::AdjustRadius:
                        drv.radii_uvm[sg][slot] = r.new_h_search;
                        drv.active_set_uvm[sg][write_idx++] = slot;
                        break;
                    default:
                        /* Unknown enum = Spec bug. Soft bad-stop + drop the
                         * slot as converged (do not re-queue): the run still
                         * stops at the next poll, but stays lockstep with
                         * peers through this iteration's collectives. The
                         * stderr above surfaces the bug — not masked. */
                        if (ThisTask == 0) {
                            fprintf(stderr,
                                "[run_neighbor_loop_iterative<%s>] FATAL: Spec::after_iter "
                                "returned unknown IterStatus=%d at iter=%d sg=%d slot=%d "
                                "(particle index i=%d). Valid values: 0=Converged, "
                                "1=NeedsMore, 2=AdjustRadius.\n",
                                Spec::loop_name, (int)r.status,
                                drv.iter_index, sg, slot, i);
                            fflush(stderr);
                        }
                        endrun(81206);
                        break;
                }
            }
            drv.active_set_size[sg]     = write_idx;
            drv.local_active_per_sg[sg] = write_idx;
            drv.local_active_total     += write_idx;
        }


        /* (b.5) Mode A buffer-exceedance rebuild trigger.
         * After radii mutation on AdjustRadius, check if any still-active
         * slot's radius now exceeds its buffered radius (built oversized
         * at last build). If so, invalidate that subgroup's cached CSR;
         * next iter's dispatch will rebuild. Mirrors legacy density_gpu.cc:
         * 152-188 rebuild trigger. */
        if (path == DispatchPath::ModeA_GPU_NGL) {
            for (int sg = 0; sg < args.num_subgroups; sg++) {
                if (!drv.mode_a_csr_valid[sg]) continue;
                int n_compacted = drv.active_set_size[sg];
                for (int k = 0; k < n_compacted; k++) {
                    int slot = drv.active_set_uvm[sg][k];
                    if (drv.radii_uvm[sg][slot] > drv.mode_a_csr_buffered_h[sg][slot]) {
                        drv.mode_a_csr_valid[sg] = false;
                        break;
                    }
                }
            }
        }

        /* (c) Spec::after_iter_global — host-only, no-MPI (TRAP 7). */
        if constexpr (nlr_spec_has_after_iter_global_v<Spec>) {
            Spec::after_iter_global(args, drv);
        }

        /* (d) Per-iter MPI Allreduce on PER-SUBGROUP local_active counts.
         * Single array Allreduce-SUM gives per-sg
         * global activity; sum is global_active_total for the break check.
         * Globally-converged subgroups (global_active_per_sg[sg]==0) skip
         * the per-iter dispatch + collective on the next iteration. */
        if (NTask > 1) {
            MPI_Allreduce(drv.local_active_per_sg.data(),
                          drv.global_active_per_sg.data(),
                          args.num_subgroups, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        } else {
            drv.global_active_per_sg = drv.local_active_per_sg;
        }
        drv.global_active_total = 0;
        for (int sg = 0; sg < args.num_subgroups; sg++) {
            drv.global_active_total += drv.global_active_per_sg[sg];
        }
        if (drv.global_active_total == 0) break;
    }

    /* ===== Max-iter exceeded check (matches legacy density / ags_density endrun(1155)) ===== */
    if (drv.iter_index >= Spec::max_iters && drv.global_active_total > 0) {
        if constexpr (nlr_spec_has_on_max_iter_exceeded_v<Spec>) {
            Spec::on_max_iter_exceeded(drv);
        } else {
            nlr_default_on_max_iter_exceeded<Spec>(drv);
        }
    }

    /* ===== Final apply_active_writeback (final-only) ===== */
    /* Fires once per active across all subgroups, on the final iteration's
     * accum (whether that iteration was Converged for that slot, or
     * max_iters terminated for everyone). The active_list semantics for
     * apply_active_writeback are the SUBGROUP'S full active_indices —
     * every active particle gets its converged result written back. */
    for (int sg = 0; sg < args.num_subgroups; sg++) {
        const NlrSubgroup& sgr = args.subgroups[sg];
        const int n_total = sgr.num_active_local;
        /* This loop is bounded by the subgroup's ORIGINAL active count, not by
         * how many are still iterating, so it also has to skip a subgroup whose
         * per-active state could never be staged: there are no results to write
         * back, and the constructor has already asked for the stop. */
        if (n_total <= 0 || drv.accum_uvm[sg] == nullptr) continue;
        /* Use args: Mode A iterative refreshed
         * P/CellP/num_total after ghost import; apply_active_writeback hooks
         * may read these for correctness. Mode B leaves args ==
         * base args. */
        neighbor_loop_args sub = drv.args;
        sub.active_list = sgr.active_indices;
        sub.num_active  = n_total;
        for (int slot = 0; slot < n_total; slot++) {
            int i = sgr.active_indices[slot];
            /* If Spec opts into the iterative-variant hook, route the
             * converged radius + IterScratch through it. Specs that don't
             * declare the iterative hook fall through to the original. */
            if constexpr (nlr_spec_has_apply_active_writeback_iterative_v<Spec>) {
                Spec::apply_active_writeback_iterative(
                    sub, slot, i,
                    drv.accum_uvm[sg][slot],
                    drv.radii_uvm[sg][slot],
                    drv.scratch_uvm[sg][slot]);
            } else {
                Spec::apply_active_writeback(sub, slot, i, drv.accum_uvm[sg][slot]);
            }
        }
    }


    }  /* end inner scope: driver destructs HERE */

    /* ===== Mode B hard-corridor enforcement =====
     * Check AFTER final apply_active_writeback AND after driver destruction
     * so the invariant covers the FULL Mode B path including writeback hooks
     * and driver cleanup. Mode A paths legitimately advanced counters; skip
     * enforcement. */
    if (path == DispatchPath::ModeB_HostWalker) {
        const bool drift_violation = (g_global_drift_counter      != s_drift0);
        const bool ghost_violation = (g_ghost_import_counter      != s_ghost0);
        const bool arena_violation = (g_gpu_arena_acquire_counter != s_arena0);
        const bool np_violation    = (NumPart != s_np0);
        if (drift_violation || ghost_violation || arena_violation || np_violation) {
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            fprintf(stderr,
                "[NLR_ITER CORRIDOR ABORT rank=%d caller=%s path=mode_b] Mode B "
                "iterative path violated tiny-N corridor invariant during "
                "run_neighbor_loop_iterative. Counter deltas: drift=%llu "
                "ghost=%llu arena=%llu NumPart_pre=%d NumPart_post=%d\n",
                rank, Spec::loop_name,
                (unsigned long long)(g_global_drift_counter      - s_drift0),
                (unsigned long long)(g_ghost_import_counter      - s_ghost0),
                (unsigned long long)(g_gpu_arena_acquire_counter - s_arena0),
                s_np0, NumPart);
            fflush(stderr);
            /* one-shot self-check at iter-runner end (no loop); function returns next,
             * no intervening collective -- soft bad-stop + fall through, drains at the
             * next phase-boundary poll. */
            endrun(81213);
        }
    }
    /* Driver destructor frees per-subgroup UVM allocations on scope exit. */
}

/* ============================================================================
 * Explicit template instantiations — one per migrated caller.
 *
 * Forgetting an instantiation = clean linker error at the call site. Each
 * Spec must opt in here and is implicitly acknowledging its declared
 * Spec::sidx_cache_kind; cf. nlr_resolve_sidx_cache.
 *
 * ========================================================================== */

#ifdef SINK_PARTICLES
template void run_neighbor_loop<SinkEnv1Spec>(const neighbor_loop_args&);
template void run_neighbor_loop<SinkFeedSpec>(const neighbor_loop_args&);
template void run_neighbor_loop<SinkSwkSpec>(const neighbor_loop_args&);
#if defined(SINK_GRAVACCRETION) && (SINK_GRAVACCRETION == 0)
template void run_neighbor_loop<SinkEnv2Spec>(const neighbor_loop_args&);
#endif
#endif


/* AgsDensitySpec — first production iterative Spec consumer of the
 * runner + partition-by-subgroup + sticky-call-scope wakeup
 * traits. See gravity/ags_density_loop.{h,cc}. */
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
template void run_neighbor_loop_iterative<AgsDensitySpec>(const neighbor_loop_args_iterative&);
/* AgsForceSpec — single-pass-iterative-shaped Spec for
 * non-gas AGS force loop (max_iters=1, after_iter always Converged); uses
 * the iterative path so the multi-subgroup contract is available. See
 * gravity/ags_force_loop.{h,cc}. */
template void run_neighbor_loop_iterative<AgsForceSpec>(const neighbor_loop_args_iterative&);
#endif

/* DensitySpec — hydro density runner port.
 * Single gas-only subgroup, no after_iter P/CellP writes,
 * uses apply_active_writeback_iterative for single-valued radius
 * channeling. See hydro/density_loop.{h,cc}. */
template void run_neighbor_loop_iterative<DensitySpec>(const neighbor_loop_args_iterative&);

/* MechFBSpec — mechanical-feedback runner port
 * (physics-complete pair kernel + full state-machine Spec
 * contract). 6-mode iterative state machine over loop_iteration
 * {-2,-1,0,1,2,3}; mode_a_rebuild_csr_every_iter=false preserves the legacy
 * 1-CSR-shared-across-6-modes optimization. Mode A multi-rank ghost-side
 * writes hit a Kokkos::abort; supporting them needs a custom MechFBGasDelta
 * ghost-writeback callback + lazy d_gas_iter, which are not implemented. See
 * galaxy_sf/mechfb_loop.{h,cc}. */
#ifdef GALSF_FB_MECHANICAL
template void run_neighbor_loop_iterative<MechFBSpec>(const neighbor_loop_args_iterative&);
#endif

/* ThermalFBSpec — thermal-feedback runner port.
 * Non-iterative scatter (Type-4 stars → gas neighbors); ActiveReduceOnly +
 * manifest-bundle ghost_writeback; sink_feed pattern. See
 * galaxy_sf/thermal_fb_loop.{h,cc}. The Spec definition is gated on
 * GALSF_FB_THERMAL in thermal_fb_loop.h, so the explicit instantiation must
 * sit inside the same #ifdef (non-thermal Configs would
 * otherwise hit an undefined type at this template instantiation site). */
#ifdef GALSF_FB_THERMAL
template void run_neighbor_loop<ThermalFBSpec>(const neighbor_loop_args&);
#endif

/* Cellcorrections corridor: CellcorrectionsSpec — first-pass
 * volume corrections (Volume_1 = sum_j Volume_0[j]^2 wk(r, h_j)).
 * NotIterative GasOnly Spec, no j-side writes, no ghost-writeback;
 * first corridor consumer in the chain. See
 * hydro/cellcorrections_loop.{h,cc}. */
#ifdef HYDRO_VOLUME_CORRECTIONS
template void run_neighbor_loop<CellcorrectionsSpec>(const neighbor_loop_args&);
#endif

/* GradientsSpec — runner port of the legacy `gradient_evaluate_gpu`
 * walker. Broad active list (Type==0 && Mass>0) matching the legacy GPU
 * walker; narrow GasGrad_isactive filter stays at the neighbor side inside
 * gradient_accumulate_neighbor. Symmetric gas-gas topology — hydro
 * corridor consumer (after CellcorrectionsSpec, before HydroForceSpec).
 * See hydro/gradients_loop.{h,cc}. */
template void run_neighbor_loop<GradientsSpec>(const neighbor_loop_args&);

#ifdef MHD_CONSTRAINED_GRADIENT
/* GradientsIterSpec — slim variant for the constrained-gradient iterations
 * (grad_iter>0). Accumulates the slim GasGraddata_out_iter_ (FaceDotB +
 * MIDPOINT PhiGrad) instead of the full GasGraddata_out_. Gated with the same
 * #ifdef that gates the Spec definition (guard the instantiation with the
 * Spec's own gate to avoid an undefined type in non-MHD builds). See
 * hydro/gradients_loop.{h,cc}. */
template void run_neighbor_loop<GradientsIterSpec>(const neighbor_loop_args&);
#endif

/* HydroForceSpec — runner port of the legacy `hydro_evaluate_gpu` walker.
 * Final hydro-corridor consumer (after CellcorrectionsSpec and
 * GradientsSpec). uses_ghost_writeback=true with a snapshot-diff bundle
 * (PARTICLE_MAX(wakeup) + MFV GAS_ADD(dMass)) for Mode A imported-ghost
 * lifecycle; Mode B direct-owner-rank j-writes via request-driven P2P.
 * See hydro/hydro_force_loop.{h,cc}. */
template void run_neighbor_loop<HydroForceSpec>(const neighbor_loop_args&);

/* RadFBRPSpec — local radiation-pressure winds.
 * Iterative 2-pass (iter 0 wt_sum aggregation; iter 1 kick application).
 * Ghost-writeback bundle with PARTICLE_ADD_VEC3 + new GAS_ADD_VEC3 ops.
 * iter-gating via Aux::iter_index (set in reset_per_iter_device_context);
 * inter-iter wt_sum staged through IterScratch by after_iter_global. See
 * galaxy_sf/radfb_rp_loop.{h,cc}. */
#ifdef GALSF_FB_FIRE_RT_LOCALRP
template void run_neighbor_loop_iterative<RadFBRPSpec>(const neighbor_loop_args_iterative&);
#endif

/* DMDispersionSpec — DM velocity dispersion runner port.
 * Gas actives (Type==0) iterate bisection on KernelRadiusDM to enclose
 * 64±48 DM (Type==1) neighbors; accumulates unweighted Vel sums for dispersion.
 * ActiveReduceOnly + SidxCacheKind::None (DM tbm matches neither GasOnly nor
 * AllTypes cache). apply_active_writeback_iterative + Aux finalize pattern
 * mirrors DensitySpec exactly. See galaxy_sf/dm_dispersion_loop.{h,cc}. */
#ifdef DM_DISPERSION_LOOP_ACTIVE
template void run_neighbor_loop_iterative<DMDispersionSpec>(const neighbor_loop_args_iterative&);
#endif

/* CBE gradients corrective architecture pivot
 * (sidm/cbe_integrator_gradients.{h,cc}). CBEGradSpec is non-iterative
 * (paralleling DMGradSpec); the two passes (raw LSQ then pairwise BJ-style
 * conservative limiter) are orchestrated by CBEGrad_gradient_calc() at the
 * toplevel via Aux::loop_iteration. Persistent storage on
 * P[i].Gradients_CBE_basis_moments; standard P[] ghost transport carries
 * gradients across ranks (no scratch, no custom Alltoallv). */
#if defined(CBE_INTEGRATOR_WITHGRADIENTS)
template void run_neighbor_loop<CBEGradSpec>(const neighbor_loop_args&);
#endif

/* RtSrcInjectionSpec — radiation source
 * injection runner port. Non-iterative scatter (non-gas sources → gas); the
 * toplevel builds the active list directly (Aux::host_locals) so
 * nlr_build_active_list is not used. SYMMETRIC search matches the legacy GPU
 * evaluator (correctness-required under RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION).
 * Ghost-writeback bundle uses three new generic ops (GAS_ADD_ARRAY,
 * GAS_ADD_2D, GAS_ADD_VEC3_ARRAY); GRAIN_RDI_TESTPROBLEM boundary condition
 * is imposed by an owner-local post-runner fixup in
 * radiation/rt_source_injection.cc (not inside the pair kernel). See
 * radiation/rt_source_injection_loop.{h,cc}. */
#ifdef RT_SOURCE_INJECTION
template void run_neighbor_loop<RtSrcInjectionSpec>(const neighbor_loop_args&);
#endif

/* difffilter: DiffFilterSpec + DynDiffSpec — TURB_DIFF_DYNAMIC
 * velocity-smoothing + dynamic-Smagorinsky loops. Non-iterative,
 * scaled-symmetric gas-gas (symmetric_neighbor_radius_scale = TurbDynamicDiffFac),
 * pure i-side reduce. See turb/difffilter_loop.{h,cc}. */
#ifdef TURB_DIFF_DYNAMIC
template void run_neighbor_loop<DiffFilterSpec>(const neighbor_loop_args&);
template void run_neighbor_loop<DynDiffSpec>(const neighbor_loop_args&);
#endif

/* dm_fuzzy: DMGradSpec — DM_FUZZY higher-order density-gradient
 * estimator. Non-iterative, one-way DM->DM, fixed radius P[i].AGS_KernelRadius,
 * pure i-side reduce. The toplevel DMGrad_gradient_calc issues one call per
 * (gradient-pass, bm group) with args.neighbor_type_mask_override = bm. See
 * sidm/dm_fuzzy_loop.{h,cc}. */
#ifdef DM_FUZZY
template void run_neighbor_loop<DMGradSpec>(const neighbor_loop_args&);
#endif

/* grain_physics: GrainBackrxSpec + GrainRTGasSpec + GrainRTGrainSpec.
 * Three independent non-iterative Specs for the grain_physics neighbor loops
 * (grain→gas backreaction with ghost writeback; gas→grain
 * RT opacity, pure i-side; grain→gas radiation acceleration,
 * pure i-side). Replaces grain_backrx_evaluate_gpu +
 * interpolate_fluxes_opacities_gasgrains_evaluate_gpu in
 * solids/grain_physics_gpu.cc (retired in cleanup commit). See
 * solids/grain_physics_loop.{h,cc}. */
#if defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION) && defined(GRAIN_BACKREACTION)
template void run_neighbor_loop<GrainBackrxSpec>(const neighbor_loop_args&);
#endif
#if defined(DO_FLUID_ALTSPECIES_DRAG_CALCULATION) && defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
template void run_neighbor_loop<GrainRTGasSpec>(const neighbor_loop_args&);
template void run_neighbor_loop<GrainRTGrainSpec>(const neighbor_loop_args&);
#endif

