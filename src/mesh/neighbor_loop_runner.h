/* mesh/neighbor_loop_runner.h — generic neighbor-loop runner for GIZMO.
 *
 * Single source of truth for the neighbor-loop interface every GIZMO physics
 * loop is migrated onto. The interface is the binding contract; new patterns
 * (write modes, iterative loops, identity sidecars) are DECLARED in the
 * interface so callers don't reshape when their needs are added; runtime
 * policy on unsupported features is documented below (TL;DR: Mode A fallback
 * in production; abort if force-Mode-B; static_assert only for internally
 * inconsistent specs).
 *
 * The user-facing copy-pasteable Spec skeleton is below ("NeighborLoopSpec —
 * copy-pasteable skeleton for a new physics loop"). For a worked example,
 * read sinks/sink_env1_loop.h.
 *
 * ============================================================================
 * Three-epoch host/device staging contract
 * ============================================================================
 *
 *   The runner separates per-active staging into three timing points so a
 *   spec can match the byte-exact state epoch the legacy device kernel saw.
 *   For sink_env1 these epochs reproduce the prior sink_environment.cc /
 *   sink_environment_gpu.cc timing exactly.
 *
 *   (1) HOST, pre-arena:
 *         search_radius(args, slot, i) -> double
 *       Reads external state (e.g. P[i].KernelRadius) BEFORE arena_acquire
 *       / drift / freshness. Equivalent to the radii staging that legacy
 *       callers performed in their host loop just after building the
 *       active-particle list. Result is staged into a UVM radii[num_active]
 *       array; Mode A's gpu_ngb_list_build receives it directly, Mode B's
 *       walker reads it per-query.
 *
 *   (2) HOST, pre-arena:
 *         populate_call_scalars(args) -> CallScalars
 *       Captures per-call scalar globals (cosmology factors, gravity
 *       constant, kernel radii, etc.) into a typed POD struct. Captured
 *       once per call; replicated by-value into the device-side lambda
 *       capture so the kernel never reaches into globals. CallScalars
 *       must be std::is_trivially_copyable_v.
 *
 *   (3) DEVICE (or host for Mode B walker), post-NGL-build:
 *         KOKKOS_INLINE_FUNCTION
 *         load_active(ctx, slot, i, h_search, cs) -> ActiveData
 *       Reads ctx.P[i] / ctx.CellP[i] (UVM-resident, post-arena/drift)
 *       and combines with the host-staged h_search and CallScalars.
 *       Runner launches a tiny Kokkos parallel_for that fills
 *       ActiveData[num_active] in UVM. Mode A's pair-kernel launch reads
 *       this array; Mode B/Brute walker calls the SAME function host-side
 *       per query. Same KOKKOS_INLINE_FUNCTION on both paths.
 *
 *   NeighborData is DEVICE-CALLABLE (and host-callable in Mode B):
 *     KOKKOS_INLINE_FUNCTION
 *     load_neighbor(ctx, j, id, active) -> NeighborData
 *     Receives a DeviceContext containing UVM-resident pointers (P, CellP,
 *     plus spec-specific extension if Spec::DeviceContext extends the
 *     base). Same function called from device (Mode A) and host (Mode B,
 *     Brute) — KOKKOS_INLINE_FUNCTION compiles to inline on both.
 *     Implementation must use only device-safe constructs (no std:: I/O,
 *     no host-only helpers); same constraint as the pair_kernel.
 *
 *   Accumulator zeroing:
 *     KOKKOS_INLINE_FUNCTION zero_accum(AccumData& out)
 *     Spec-owned. Generic runner does NOT bake in memset on AccumData;
 *     the spec decides (e.g. memset for POD, structured zeroing for
 *     non-trivial members).
 *
 *   Writebacks are HOST-ONLY:
 *     apply_active_writeback(args, slot, i, accum)
 *     apply_neighbor_writeback(args, j, scatter)   // scatter modes only
 *     After dispatch returns, the runner copies AccumData back to host
 *     (no-op for Mode B/Brute; UVM-coherent for Mode A) and applies
 *     writebacks on the host. Spec is responsible for any further state
 *     scatter (e.g. SinkTempInfo population).
 *
 *   pair_kernel signature (single-source-of-truth lever; contract §"Pair kernel signature"):
 *     KOKKOS_INLINE_FUNCTION
 *     static void pair_kernel(const ActiveData& a,
 *                             const NeighborData& nb,
 *                             AccumData&  active_out,
 *                             ScatterData& nb_out);
 *   Both outputs ALWAYS in the signature. ScatterData=NoScatter for
 *   ActiveReduceOnly compiles to zero-cost. Kernel writes both sides of
 *   any pair interaction in ONE place — never split into "active-side
 *   kernel" + "neighbor-side kernel" (the duplicate-logic trap this
 *   design explicitly avoids).
 *
 *   No `int j` parameter. The kernel cannot tell whether nb came from local
 *   P[j], a Mode A ghost, a Mode B remote reply, or brute. Identity arrives
 *   via the spec's NeighborData (ordinary fields) or via the runner-
 *   populated IdentitySidecar passed to load_neighbor (only fields the
 *   spec's IdentityFields opts in to are valid; others are zero).
 *
 *   Self-skip on identity: physics-meaningful fields like P[j].ID for "is
 *   this the same particle I'm asking about?" SHOULD live as a normal field
 *   on NeighborData, loaded by load_neighbor — they're physics, not
 *   identity in the contract sense. IdentityFields is reserved for
 *   transport-aware metadata (owner_rank, local_index for ownership/
 *   tie-breaking when the same i can be reached from multiple paths).
 *
 * ============================================================================
 * Dispatch policy (the three behaviors)
 * ============================================================================
 *
 *   1. Internally inconsistent spec (compile-time): static_assert at the
 *      top of run_neighbor_loop<Spec> rejects (e.g. ActiveReduceOnly with
 *      non-empty ScatterData; iterative loop without after_iter; etc.).
 *      Caller must fix the spec.
 *
 *   2. Force-Mode-B + unimplemented Mode B feature (runtime abort): when
 *      GIZMO_NLR_FORCE_MODEB=1 OR Spec::force_modeb_required = true, and
 *      the dispatch chooses Mode B but the runner doesn't yet implement
 *      that combination, the runner aborts loudly. This is the gate for
 *      unit-test suites forcing the tiny-N path.
 *
 *   3. Default production with unimplemented Mode B feature (runtime
 *      fallback): dispatch falls through to Mode A and logs a one-time
 *      warning identifying the (loop_name, missing-feature) pair. Lets
 *      callers adopt the runner immediately and pick up Mode B as the
 *      runner gains coverage.
 *
 *   GIZMO_NLR_FORCE_MODEA=1 forces Mode A globally regardless of spec or
 *   threshold (regression baseline / sanity).
 *
 * ============================================================================
 * Per-loop spec carries:
 *
 *     - identifier:       loop_name (drives env gates / PHASE0 labels)
 *     - typed payloads:   ActiveData, NeighborData, AccumData, ScatterData
 *     - device context:   DeviceContext (defaults to base; spec extends if needed)
 *     - search spec:      search_mode, neighbor_type_mask, radius_policy,
 *                         search_radius(active)
 *     - writeback spec:   write_pattern + reduction ops via apply_*_writeback
 *     - identity opt-in:  IdentityFields (defaults NoIdentity)
 *     - iter driver:      IterControl after_iter() (defaults NotIterative)
 *     - oracle compare:   compare_accum, accum_tolerance (compare_scatter +
 *                         scatter_tolerance for scatter specs)
 *     - hooks (host):     search_radius, populate_call_scalars,
 *                         apply_active_writeback,
 *                         apply_neighbor_writeback (scatter only),
 *                         populate_device_context (optional extension)
 *     - hooks (device):   load_neighbor, pair_kernel
 *
 * Generic runner does:
 *     - dispatch:  Mode A vs Mode B vs Brute — runtime choice from
 *                  Allreduce'd Nactive vs threshold + force-mode env
 *     - search:    coarse predicate, type mask + per-type hmax pruning
 *     - transport: ghost (Mode A) or peer-to-peer (Mode B)
 *     - drift:     lazy drift_particle on j candidates pre-kernel
 *     - oracle:    GIZMO_<LOOP>_ORACLE=1 → run Mode B + Brute, compare via
 *                  Spec::compare_accum within Spec::accum_tolerance
 *     - timing:    PHASE0 instrumentation per loop_name
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef GIZMO_NEIGHBOR_LOOP_RUNNER_H
#define GIZMO_NEIGHBOR_LOOP_RUNNER_H

#include <cstdint>
#include <vector>
#include <type_traits>
#include "../declarations/constants.h"  /* NODELISTLENGTH (NlrQueryEnvelope start-node list) */
#include "mode_b_local_walker.h"  /* SearchMode, RadiusPolicy already declared */

/* NLR_INLINE_FUNCTION — private host+device-callable macro for runner-owned
 * accessors (currently just NeighborLoopDeviceContextBase::particle_type).
 *
 * Aliases to KOKKOS_INLINE_FUNCTION when the consuming TU has already
 * included <Kokkos_Core.hpp> (true of all GPU TUs in this codebase: they
 * include Kokkos_Core.hpp before any header chain that pulls runner.h);
 * otherwise falls back to plain `inline` (host-only).
 *
 * Why NOT a fallback `#define KOKKOS_INLINE_FUNCTION inline` (footgun): defining the global Kokkos macro to plain `inline` from this
 * header poisons future includes of <Kokkos_Core.hpp> — Kokkos guards its
 * own definition with `#ifndef`, so a misordered include (runner.h first,
 * then Kokkos_Core.hpp) leaves KOKKOS_INLINE_FUNCTION permanently
 * mis-defined as plain inline for the rest of compilation. Every Kokkos-
 * style accessor in that TU would silently become host-only. The
 * NLR_INLINE_FUNCTION indirection avoids the global poisoning: only
 * runner-tagged accessors are affected by the fallback, and only in TUs
 * that include runner.h before Kokkos (i.e., currently never).
 *
 * Audit invariant: every GPU TU including this header MUST include
 * <Kokkos_Core.hpp> first. Verified true for all 8 current GPU consumers
 * (gravity/ags_density_gpu.cc, gravity/ags_force_gpu.cc, hydro/density_gpu.cc,
 * sinks/sink_feed_loop.cc, sinks/sink_env2_loop.cc, sinks/sink_swk_loop.cc,
 * mesh/neighbor_loop_runner.cc, mesh/test_iter_harness_loop.cc).
 */
#ifdef KOKKOS_INLINE_FUNCTION
#define NLR_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#else
#define NLR_INLINE_FUNCTION inline
#endif
#include "gpu_neighbor_list.h"    /* gpu_neighbor_list_t (NlrIterDriver Mode A fields) */

/* ============================================================================
 * Reused enums:
 *   SearchMode:    MODE_B_SEARCH_ONEWAY, MODE_B_SEARCH_SYMMETRIC
 *   RadiusPolicy:  mode_b_radius_policy_t (MODE_B_RADIUS_GAS_KERNEL, …)
 * ========================================================================== */

/* ============================================================================
 * Writeback specification
 * ========================================================================== */

/* WritePattern describes the runner-managed accum/scatter buffers ONLY.
 * It is ORTHOGONAL to whether the pair_kernel does direct j-side writes
 * (which is governed by `uses_ghost_writeback` and the bundle manifest in
 * the Spec's .cc). A Spec may declare `ActiveReduceOnly` and still write
 * to neighbor_particle.field / neighbor_cell->field from the kernel, as
 * long as it sets `uses_ghost_writeback = true` and supplies a bundle.
 * Mode A's bundle handles reverse-comm of those writes to home ranks;
 * Mode B local + Mode B remote skip the bundle (j-side writes are local
 * on the rank that owns j). sink_feed is the first port to
 * combine `ActiveReduceOnly` with a non-empty ghost-writeback bundle. */
enum class WritePattern : int {
    ActiveReduceOnly       = 0,  /* AccumData filled; ScatterData == NoScatter.
                                    Direct j-side writes still allowed via
                                    uses_ghost_writeback bundle (orthogonal). */
    SymmetricPairScatter   = 1,  /* both filled */
    NeighborScatter        = 2,  /* AccumData == NoAccum; ScatterData filled */
    GhostWritebackRequired = 3,  /* j may live on remote rank; writeback channel */
};

/* ============================================================================
 * Per-field reduction ops
 * ========================================================================== */

enum class ReductionOp : int {
    AssignAdd = 0, Assign = 1, Max = 2, Min = 3,
    AtomicAdd = 4, AtomicMax = 5, AtomicMin = 6, AtomicCAS = 7,
};

/* ============================================================================
 * RemoteEvalMode (Mode B comm strategy)
 * ========================================================================== */

enum class RemoteEvalMode : int {
    RemoteComputesAccum             = 0,
    ReturnNeighborData              = 1,
    ReturnNeighborDataWithWriteback = 2,
};

/* ============================================================================
 * Identity sidecar
 * ========================================================================== */

struct NoIdentity {};
struct IdField_ID {};
struct IdField_LocalIndex {};
struct IdField_OwnerRank {};

template <typename... Tags>
struct WithFields {};

struct IdentitySidecar {
    unsigned long long id;        /* valid iff IdField_ID requested */
    int                local_idx;  /* valid iff IdField_LocalIndex requested */
    int                owner_rank;/* valid iff IdField_OwnerRank requested.
                                   * Mode A ghost-import metadata MUST carry
                                   * owner_rank if any spec declares this. */
};

/* ============================================================================
 * SidxCacheKind — Mode A SIDX (spatial-index cache) selection
 *
 * Each spec declares which step-persistent SIDX cache its NGL build should
 * reuse. Caches are step-scoped and shared across cooperating callers
 * (e.g. sink_env1, sink_feed, sink_swk all share the all-types cache to
 * amortize the SFC sort + BVH build). The runner resolves the kind to a
 * gpu_spatial_index_t* internally; specs never name the global accessor.
 *
 * AllTypes implemented; GasOnly added for the density port. Specs
 * with per-subgroup variable masks use None so each CSR build gets a local
 * SIDX with exactly that subgroup's mask.
 *
 * Type-mask invariant: the chosen cache's recorded
 * type_bitmask MUST match the Spec's neighbor_type_mask. A mismatch
 * (e.g. gas-only Spec routing to the all-types cache) lets the walker
 * return wrong-type neighbors and triggers downstream drift/lazy-drift
 * aborts. gpu_ngb_list_build hard-aborts on mismatch via cache_tbm.
 * ========================================================================== */

enum class SidxCacheKind : int {
    AllTypes = 0,   /* gpu_step_sidx_alltypes_ptr — sink_env1/feed/swk shared (tbm=0x3f) */
    GasOnly  = 1,   /* gpu_step_sidx_ptr — gas-only callers (tbm=1), density */
    None     = 2,   /* no step-persistent cache; required for variable subgroup masks */
    /* future: PerSpecMask — implement if variable-mask callers need cache reuse. */
};

/* ============================================================================
 * ModeBEvalOMP — per-Spec eval-threading tier for Mode-B pair-kernel evaluation.
 *
 * The Mode-B discovery/collection walks are threaded by the runner for every
 * Spec above a work threshold (a runner property, not a per-spec knob). The
 * pair-kernel EVAL is classified per-spec here because only it can write
 * neighbor (j-side) state, so its thread-safety depends on the physics:
 *
 *   BitwiseReadonly — pair_kernel writes ONLY the i-side AccumData; no j-side
 *                     scatter. Threaded eval is bit-identical (per-active
 *                     candidate order is fixed). Requires uses_ghost_writeback
 *                     == false.
 *   EpsilonAtomic   — pair_kernel scatters to j via atomic ops. Threaded eval
 *                     reorders the same atomic updates -> ulp-class
 *                     nondeterminism (the class Mode-A Kokkos kernels already
 *                     exhibit). Validated by ID-sort epsilon + oracle
 *                     membership + accum_tolerance, not bitwise.
 *   SerialOnly      — evaluated on one thread. A justified final SerialOnly
 *                     carries a stated STRUCTURAL reason (order-sensitive
 *                     j-writes, RNG-order dependence, non-atomic scatter). A
 *                     spec whose tier is not yet audited is conservatively
 *                     SerialOnly with a note until its pair-kernel is verified.
 *
 * A Spec that omits the trait resolves to SerialOnly (compile-safe fallback);
 * nlr_modeb_eval_omp_label() reports an explicit SerialOnly distinctly from a
 * missing trait so an unaudited spec is never mistaken for a justified one.
 * ========================================================================== */

enum class ModeBEvalOMP : int {
    BitwiseReadonly = 0,
    EpsilonAtomic   = 1,
    SerialOnly      = 2,
};

/* ============================================================================
 * DeviceContext base
 *
 * Specs may extend by typedef-ing `DeviceContext` to a struct that publicly
 * inherits NeighborLoopDeviceContextBase and adds spec-specific UVM pointers.
 * If extension is used, the spec MUST also provide:
 *
 *   static void populate_device_context(const neighbor_loop_args& args,
 *                                       DeviceContext& ctx);
 *
 * The runner detects extension via type-trait check on Spec::DeviceContext
 * and conditionally invokes populate_device_context.
 * ========================================================================== */

struct NeighborLoopDeviceContextBase {
    struct particle_data *P;
    struct gas_cell_data *CellP;
    int                   num_total;

    /* Path-correct accessor for particle Type. Reads through the ctx's own
     * P pointer (Mode A: arena-resident P_gpu; Mode B: request-driven local
     * slab — both assigned to ctx.P by the path-specific init paths).
     * NEVER reads global P (preserves guideline #3: no globals on Mode B
     * tiny-N path).
     *
     * host+device callable via NLR_INLINE_FUNCTION, which
     * aliases to KOKKOS_INLINE_FUNCTION when Kokkos is available. The
     * partition-key assertion in run_neighbor_loop_iterative calls this
     * from host driver code, and future Specs may call it from device-side
     * helpers via Spec::active_subgroup_key. */
    NLR_INLINE_FUNCTION
    short int particle_type(int i) const { return P[i].Type; }
};

/* ============================================================================
 * NlrCommonScalars — cosmology + gravity scalars almost every Spec needs.
 *
 * Specs include this by composition in their CallScalars struct:
 *
 *     struct MyLoopCallScalars {
 *         NlrCommonScalars common;
 *         double           my_loop_specific_scalar;
 *     };
 *
 *     static MyLoopCallScalars populate_call_scalars(...) {
 *         MyLoopCallScalars s;
 *         s.common = nlr_common_scalars_from_all();
 *         s.my_loop_specific_scalar = ...;
 *         return s;
 *     }
 *
 * Trivially copyable — rides the device lambda capture and any cross-rank
 * transport without special handling. Not opt-out: missing-when-needed is
 * silent floating-point corruption, so the cost of always populating these
 * (a few doubles per call) is the safer default.
 *
 * Slightly overcomplete by design: includes fields not strictly needed by
 * every loop, so a new Spec author copying this template gets the obvious
 * cosmology factors without remembering to re-derive them.
 *
 * Unit conversion factors (All.UnitMass_in_g, UnitLength_in_cm, etc.) are
 * NOT in this struct. They are invariant for the run, so a parallel
 * NlrUnitScalars helper can be added when the first loop needs them in a
 * device-side gas-physics helper. Until then, host-side code reads
 * All.Unit* directly.
 * ========================================================================== */
struct NlrCommonScalars {
    double cf_atime;
    double cf_a2inv;
    double cf_a3inv;
    double cf_hubble_a;
    double newton_G;                /* All.G — Newton's gravitational constant in code units */
    double hubble;                  /* All.HubbleParam */
    int    comoving_integration_on;
};

NlrCommonScalars nlr_common_scalars_from_all(void);

/* Canonical host-side accessor for the global `All` struct.
 *
 * Under the per-TU AllDeviceMirror scheme the device-pass-only macro
 * makes bare `All.*` in host code of a GPU TU safe (reads the host
 * extern), but Spec::populate_call_scalars() and similar host hooks
 * are still encouraged to route through this accessor for
 * explicit-intent and as documentation of the host-snapshot semantics.
 * (The "must use this — NEVER bare All" requirement from the old All_dev
 * trap era is retired; the recommendation now is stylistic.) */
struct global_data_all_processes; /* forward decl — full struct in allvars.h */
const struct global_data_all_processes * nlr_host_all_ptr(void);

/* Forward declaration — full struct neighbor_loop_args is defined later in
 * this header. The DeviceContext trait helpers below need it as a parameter
 * type only (no member access), so the forward decl is sufficient. */
struct neighbor_loop_args;

/* ============================================================================
 * DeviceContext extension trait
 *
 * Detects whether a Spec extends DeviceContext beyond NeighborLoopDeviceContextBase.
 * When true, the runner invokes Spec::populate_device_context(args, ctx) after
 * the base members (P, CellP, num_total) are initialised, letting the Spec
 * stash UVM-resident pointers (e.g., per-active host-staged input buffers)
 * into the extended ctx for device-side load_active / load_neighbor / pair_kernel.
 * After the runner's evaluate phase completes, the runner invokes
 * Spec::cleanup_device_context(args, ctx) on every path so the Spec can free
 * any UVM allocations made during populate. cleanup is also detected
 * conditionally — if the Spec doesn't define it (e.g., a derived ctx that only
 * caches an externally-owned pointer), the runner emits no call.
 *
 * When DeviceContext == base (sink_env1 and any Spec that uses base ctx
 * unmodified), the runner skips both populate and cleanup — no behavior change
 * vs. pre-4.A.0.
 *
 * Spec contract: every Spec MUST declare `using DeviceContext = ...;` (either
 * the base or its own derived type). Specs whose DeviceContext derives MUST
 * also provide `static void populate_device_context(const neighbor_loop_args&,
 * DeviceContext&);` and SHOULD provide
 * `static void cleanup_device_context(const neighbor_loop_args&, DeviceContext&);`
 * if populate_device_context allocated owning resources. populate without
 * cleanup is allowed when populate only caches pointers to memory whose
 * lifetime is externally guaranteed (e.g., args.aux-rooted host buffers).
 *
 * Compile-time DeviceContext invariants:
 *   - Spec::DeviceContext MUST publicly derive from NeighborLoopDeviceContextBase
 *     (or BE the base) so the helper templates can pass it as `const Base&`
 *     when convenient.
 *   - Spec::DeviceContext MUST be trivially copyable: the runner captures it
 *     by value into Kokkos device lambdas. No std::vector, no owning smart
 *     pointers, no objects with non-trivial copy / dtor.
 * Both invariants are checked by static_assert at the runner call sites.
 * ========================================================================== */
template <typename Spec>
struct nlr_spec_has_extended_device_context {
    static constexpr bool value =
        !std::is_same<typename Spec::DeviceContext,
                      NeighborLoopDeviceContextBase>::value;
};

template <typename Spec>
constexpr bool nlr_spec_has_extended_device_context_v =
    nlr_spec_has_extended_device_context<Spec>::value;

/* SFINAE detection of optional Spec::cleanup_device_context. Returns false
 * for Specs that don't define it; the runner uses if constexpr to skip the
 * call when absent. */
template <typename Spec, typename = void>
struct nlr_spec_has_cleanup_device_context : std::false_type {};

template <typename Spec>
struct nlr_spec_has_cleanup_device_context<
    Spec,
    std::void_t<decltype(Spec::cleanup_device_context(
        std::declval<const neighbor_loop_args&>(),
        std::declval<typename Spec::DeviceContext&>()))>>
    : std::true_type {};

template <typename Spec>
constexpr bool nlr_spec_has_cleanup_device_context_v =
    nlr_spec_has_cleanup_device_context<Spec>::value;

/* SFINAE detection of optional Spec::set_oracle_brute_pass. j-side-write
 * Specs use this to suppress side-effects on the oracle's BRUTE pass —
 * without it, the runner's oracle path (run_mode_b_local_with_oracle,
 * run_mode_b_remote_impl<true>) would call the pair body twice (tree + brute)
 * and apply j-side writes (atomic_exchange / atomic_add into P[j] / CellP[j])
 * twice, corrupting additive fields like Injected_Sink_Energy. Oracle is a
 * VALIDATION-ONLY path so this is debug-relevant only; production paths
 * (default, force-A, force-B-without-oracle) never call this hook. */
template <typename Spec, typename = void>
struct nlr_spec_has_set_oracle_brute_pass : std::false_type {};

template <typename Spec>
struct nlr_spec_has_set_oracle_brute_pass<
    Spec,
    std::void_t<decltype(Spec::set_oracle_brute_pass(
        std::declval<typename Spec::DeviceContext&>(), bool{}))>>
    : std::true_type {};

template <typename Spec>
constexpr bool nlr_spec_has_set_oracle_brute_pass_v =
    nlr_spec_has_set_oracle_brute_pass<Spec>::value;

/* SFINAE detection + resolution of the optional per-Spec eval-threading tier
 * Spec::modeb_eval_omp (ModeBEvalOMP). A Spec that omits it resolves to
 * SerialOnly; the *_is_explicit_v trait lets diagnostics report a missing trait
 * separately from a declared SerialOnly (a missing trait is an unaudited spec,
 * not a justified serial one — the Mode-B eval-threading audit closes only when
 * every spec carries an explicit, structurally-justified tier). */
template <typename Spec, typename = void>
struct nlr_spec_has_modeb_eval_omp : std::false_type {};

template <typename Spec>
struct nlr_spec_has_modeb_eval_omp<
    Spec, std::void_t<decltype(Spec::modeb_eval_omp)>> : std::true_type {};

template <typename Spec>
constexpr bool nlr_spec_modeb_eval_omp_is_explicit_v =
    nlr_spec_has_modeb_eval_omp<Spec>::value;

template <typename Spec>
constexpr ModeBEvalOMP nlr_spec_modeb_eval_omp() {
    if constexpr (nlr_spec_has_modeb_eval_omp<Spec>::value) {
        return Spec::modeb_eval_omp;
    } else {
        return ModeBEvalOMP::SerialOnly;
    }
}

/* Human-readable tier label for the GX_MODEB_EXPORT eval-threading audit field.
 * A resolved SerialOnly prints "(explicit)" vs "(missing_trait)" so justified
 * serial rows are distinguishable from unaudited specs that forgot the trait. */
inline const char *nlr_modeb_eval_omp_label(ModeBEvalOMP tier, bool is_explicit) {
    switch(tier) {
        case ModeBEvalOMP::BitwiseReadonly: return "BitwiseReadonly";
        case ModeBEvalOMP::EpsilonAtomic:   return "EpsilonAtomic";
        case ModeBEvalOMP::SerialOnly:
            return is_explicit ? "SerialOnly(explicit)" : "SerialOnly(missing_trait)";
    }
    return "SerialOnly(missing_trait)";
}

/* SFINAE detection of optional Spec::bind_active_to_eval_context.
 *
 * A Spec whose ActiveData carries rank-local context fields (rank-local
 * pointers like P_base / CellP_base, rank-local index bounds, gas-delta
 * write targets, the oracle_dry_run flag — anything that depends on which
 * eval pass / which rank is doing the walk) MUST define this hook. The
 * runner invokes it for every active right before pair_kernel inside
 * evaluate_pairs_post_drift, so the active's rank-local + eval-pass-local
 * snapshots are refreshed to the eval ctx in use.
 *
 * Signature: static void Spec::bind_active_to_eval_context(
 *                                  const DeviceContext& eval_ctx,
 *                                  ActiveData& active);
 *
 * Why per-eval-pass (not just per-flatten on the receiver):
 *   - PEER-eval correctness on a remote rank: ActiveData arrived via MPI
 *     envelope carrying the SENDER's rank-local snapshots; receiver must
 *     overwrite with its own values.
 *   - Oracle correctness on ANY rank: the runner's oracle paths (inline
 *     compare / iterative dual) build a LOCAL ctx_oracle_* by copying ctx
 *     and flipping oracle_dry_run via set_oracle_brute_pass — but the
 *     ActiveData snapshot taken at load_active still carries the OLD
 *     oracle_dry_run. Without per-eval-pass binding, pair_kernel sees
 *     oracle_dry_run=false on the brute pass too, j-side writes fire
 *     twice, and the oracle silently reports no mismatch (atomic-min is
 *     idempotent; atomic-add doubles → 2× coupling).
 *
 * Specs whose ActiveData is purely physical (pos/vel/mass/scalars only) do
 * NOT need this hook; the trait returns false and the runner skips the call.
 *
 * Incident: MechFBActiveState (galaxy_sf/mechfb_loop.h:54) embeds P_base,
 * CellP_base, LocalGasMechFBInfoTemp, d_gas_iter, num_local_gas,
 * num_local_particles, oracle_dry_run — all rank-local or eval-pass-local.
 * Ship-then-deref on the peer caused np=2 wind_singlestar SIGSEGV in
 * mechanical_fb_pair_kernel during PEER eval. */
template <typename Spec, typename = void>
struct nlr_spec_has_bind_active_to_eval_context : std::false_type {};

template <typename Spec>
struct nlr_spec_has_bind_active_to_eval_context<
    Spec,
    std::void_t<decltype(Spec::bind_active_to_eval_context(
        std::declval<const typename Spec::DeviceContext&>(),
        std::declval<typename Spec::ActiveData&>()))>>
    : std::true_type {};

template <typename Spec>
constexpr bool nlr_spec_has_bind_active_to_eval_context_v =
    nlr_spec_has_bind_active_to_eval_context<Spec>::value;


/* nlr_spec_symmetric_j_radius_scale — optional Spec hook.
 *
 * A Spec whose SYMMETRIC search must reach scaled-j neighbors (effective
 * pair support max(h_i, scale*h_j)) declares:
 *
 *   static double symmetric_neighbor_radius_scale();
 *
 * Absent ⇒ 1.0 (plain symmetric search — every pre-existing Spec). The runner
 * threads the returned value into the Mode A NGL build (j_kernel_radius_scale)
 * AND the Mode B local/brute walks (j_radius_scale), so Mode A and Mode B
 * agree on neighbor inclusion. Used by the TURB_DIFF_DYNAMIC wide-filter
 * Specs, which return All.TurbDynamicDiffFac. */
template <typename Spec, typename = void>
struct nlr_spec_has_symmetric_radius_scale : std::false_type {};

template <typename Spec>
struct nlr_spec_has_symmetric_radius_scale<
    Spec,
    std::void_t<decltype(Spec::symmetric_neighbor_radius_scale())>>
    : std::true_type {};

template <typename Spec>
double nlr_spec_symmetric_j_radius_scale() {
    if constexpr (nlr_spec_has_symmetric_radius_scale<Spec>::value) {
        return Spec::symmetric_neighbor_radius_scale();
    } else {
        return 1.0;
    }
}

/* nlr_spec_supply_band_dominated — optional Spec hook.
 *
 * A Spec whose supply-side reach is provably bounded by the per-type node band
 * the sender opener walks against declares:
 *
 *   static constexpr bool supply_band_dominated = true;
 *
 * Absent ⇒ false, which keeps that Spec's SYMMETRIC ghost import on the
 * broadcast path — today's behavior, and the safe answer when the bound is not
 * established. Only a proven Spec sets it, and the proof belongs in a comment
 * at the declaration. The runner threads this into the ghost_exchange spec, so
 * routed SYMMETRIC discovery is selected by this structural property alone and
 * never by which loop is calling. */
template <typename Spec, typename = void>
struct nlr_spec_has_supply_band_dominated : std::false_type {};

template <typename Spec>
struct nlr_spec_has_supply_band_dominated<
    Spec,
    std::void_t<decltype(Spec::supply_band_dominated)>>
    : std::true_type {};

template <typename Spec>
constexpr bool nlr_spec_supply_band_dominated() {
    if constexpr (nlr_spec_has_supply_band_dominated<Spec>::value) {
        return Spec::supply_band_dominated;
    } else {
        return false;
    }
}

/* RAII cleanup guard for the runner's DeviceContext. Construct one right
 * after Spec::populate_device_context returns; destruction at function exit
 * (any path) conditionally invokes Spec::cleanup_device_context if the Spec
 * defined one. No-op for base-ctx Specs. */
template <typename Spec>
struct NlrDeviceContextCleanupGuard {
    const neighbor_loop_args& args;
    typename Spec::DeviceContext& ctx;
    NlrDeviceContextCleanupGuard(const neighbor_loop_args& a,
                                  typename Spec::DeviceContext& c)
        : args(a), ctx(c) {}
    NlrDeviceContextCleanupGuard(const NlrDeviceContextCleanupGuard&) = delete;
    NlrDeviceContextCleanupGuard& operator=(const NlrDeviceContextCleanupGuard&) = delete;
    ~NlrDeviceContextCleanupGuard() {
        if constexpr (nlr_spec_has_cleanup_device_context_v<Spec>) {
            Spec::cleanup_device_context(args, ctx);
        }
    }
};

/* RAII brute-pass guard. Construct at entry of an oracle
 * dispatch helper; destructor flips set_oracle_brute_pass(false) at scope
 * exit even on mid-walk abort. SFINAE-no-op for Specs without
 * Spec::set_oracle_brute_pass. Bound to the OWNING ctx_oracle reference
 * — production ctx is NEVER touched. */
template <typename Spec>
struct NlrOracleBrutePassGuard {
    typename Spec::DeviceContext& ctx_oracle;
    bool toggled = false;

    explicit NlrOracleBrutePassGuard(typename Spec::DeviceContext& c)
        : ctx_oracle(c) {
        if constexpr (nlr_spec_has_set_oracle_brute_pass_v<Spec>) {
            Spec::set_oracle_brute_pass(ctx_oracle, true);
            toggled = true;
        }
    }
    ~NlrOracleBrutePassGuard() {
        if (toggled) {
            Spec::set_oracle_brute_pass(ctx_oracle, false);
        }
    }
    NlrOracleBrutePassGuard(const NlrOracleBrutePassGuard&) = delete;
    NlrOracleBrutePassGuard& operator=(const NlrOracleBrutePassGuard&) = delete;
};

/* SFINAE detection of optional Spec::radius_tolerance (radius convergence tolerance is not the same
 * semantic object as the accumulator comparison tolerance). Defaults to
 * Spec::accum_tolerance when not declared — iterative Specs SHOULD declare
 * their own value. */
template <typename Spec, typename = void>
struct nlr_spec_has_radius_tolerance : std::false_type {};
template <typename Spec>
struct nlr_spec_has_radius_tolerance<
    Spec,
    std::void_t<decltype(Spec::radius_tolerance)>> : std::true_type {};

template <typename Spec>
constexpr double nlr_spec_radius_tolerance_v =
    nlr_spec_has_radius_tolerance<Spec>::value
        ? static_cast<double>(Spec::radius_tolerance)
        : static_cast<double>(Spec::accum_tolerance);

/* ============================================================================
 * Partition-by-subgroup contract
 *
 * Two value-traits + one SFINAE method-detector that together implement the
 * "actives_partition_by_subgroup" contract.
 *
 * Value-trait `actives_partition_by_subgroup` (default false): when true, the
 * runner promises that each active belongs to exactly ONE subgroup (the
 * subgroup it was placed into at iter-0). The Spec MUST also provide
 * `active_subgroup_key(ctx, i, scalars)` returning the bm key (= the
 * subgroup's `j_type_bitmask`) so the runner can assert no drift. Used by
 * ags_density (Type-fixed → bm-fixed); the harness Specs default to false
 * (existing union-merge semantics unchanged).
 *
 * Method-detector `active_subgroup_key`: detected via SFINAE on the 3-arg
 * signature (const DeviceContext&, int, const CallScalars&) returning int.
 * Required when actives_partition_by_subgroup=true; ignored otherwise.
 *
 * Value-trait `mode_a_rebuild_csr_every_iter` (default false): correctness
 * fallback for Specs whose CSR-buffer-reuse semantics can't be validated
 * against legacy/oracle. When true, the runner unconditionally invalidates
 * the cached Mode A CSR at the start of every iter > 0 (bypassing the
 * buffer-exceedance trigger entirely). `Spec::mode_a_csr_buffer_factor`
 * remains required (>1.0 static_assert) but is not consulted on this path.
 * ========================================================================== */

template <typename Spec, typename = void>
struct nlr_spec_actives_partition_by_subgroup : std::false_type {};

template <typename Spec>
struct nlr_spec_actives_partition_by_subgroup<
    Spec, std::void_t<decltype(Spec::actives_partition_by_subgroup)>>
    : std::integral_constant<bool, Spec::actives_partition_by_subgroup> {};

template <typename Spec>
constexpr bool nlr_spec_actives_partition_by_subgroup_v =
    nlr_spec_actives_partition_by_subgroup<Spec>::value;

template <typename Spec, typename = void>
struct nlr_spec_has_active_subgroup_key : std::false_type {};

template <typename Spec>
struct nlr_spec_has_active_subgroup_key<
    Spec,
    std::void_t<decltype(Spec::active_subgroup_key(
        std::declval<const typename Spec::DeviceContext&>(),
        int{},
        std::declval<const typename Spec::CallScalars&>()))>>
    : std::true_type {};

template <typename Spec>
constexpr bool nlr_spec_has_active_subgroup_key_v =
    nlr_spec_has_active_subgroup_key<Spec>::value;

template <typename Spec, typename = void>
struct nlr_spec_mode_a_rebuild_csr_every_iter : std::false_type {};

template <typename Spec>
struct nlr_spec_mode_a_rebuild_csr_every_iter<
    Spec, std::void_t<decltype(Spec::mode_a_rebuild_csr_every_iter)>>
    : std::integral_constant<bool, Spec::mode_a_rebuild_csr_every_iter> {};

template <typename Spec>
constexpr bool nlr_spec_mode_a_rebuild_csr_every_iter_v =
    nlr_spec_mode_a_rebuild_csr_every_iter<Spec>::value;

/* ------------------------------------------------------------------------- *
 * Active-source-in-pool contract (NGL source position/radius correctness).
 *
 * INVARIANT: the cached SIDX `compact_xyzh` holds neighbor-POOL state (positions
 * in [0..2], h in [3]). It is authoritative SOURCE position/radius for an active
 * particle ONLY if that particle is a member of the cached pool (its type is in
 * the pool's type_bitmask). The Mode-A NULL source_positions/radii fast-path reads
 * compact_xyzh[active_index]; for a NON-pool active (e.g. a Type-5 sink doing a
 * gas-neighbor search in a GasOnly density loop) that slot is stale/unrefreshed on
 * a reused cache (incremental drift-refresh only touches pool members), giving a
 * wrong/non-deterministic source. (Fresh-built or AllTypes caches do not have this
 * problem; only cached GasOnly + non-gas active does.)
 *
 * CONTRACT: every cached-SIDX Spec (sidx_cache_kind != None) MUST declare
 *   static constexpr bool mode_a_active_sources_in_sidx_pool = <bool>;
 *     true  -> every active source is a pool member (gas-gas, or AllTypes pool):
 *              keep the compact fast-path (no per-active position copy).
 *     false -> active sources may be non-pool (e.g. sink/star in a GasOnly loop):
 *              the runner stages explicit P[active_i].Pos (radii are already passed
 *              explicitly by the runner). compact_xyzh stays an acceleration
 *              structure only; source coords come from current particle state.
 * The static_assert in the runner bodies makes omission a COMPILE ERROR
 * (safe-by-default; a future non-pool-active GasOnly loop cannot silently inherit
 * the stale-source bug). A complementary runtime guard in gpu_ngb_list_build
 * catches direct (non-runner) callers that pass NULL with non-pool actives.
 * ------------------------------------------------------------------------- */

/* "declared": does the Spec declare mode_a_active_sources_in_sidx_pool at all? */
template <typename Spec, typename = void>
struct nlr_spec_active_sources_in_sidx_pool_declared : std::false_type {};
template <typename Spec>
struct nlr_spec_active_sources_in_sidx_pool_declared<
    Spec, std::void_t<decltype(Spec::mode_a_active_sources_in_sidx_pool)>>
    : std::true_type {};
template <typename Spec>
constexpr bool nlr_spec_active_sources_in_sidx_pool_declared_v =
    nlr_spec_active_sources_in_sidx_pool_declared<Spec>::value;

/* "value": the declared value, or a safe default (true) when absent. The default
 * is only ever reached for sidx_cache_kind==None specs (cached specs are required
 * to declare, enforced by the static_assert), and it is gated out by the cache-kind
 * check below, so it is never actually consumed for an undeclared spec. Reading
 * this trait (not Spec::member directly) avoids hard-instantiating a missing member. */
template <typename Spec, typename = void>
struct nlr_spec_active_sources_in_sidx_pool_value : std::true_type {};
template <typename Spec>
struct nlr_spec_active_sources_in_sidx_pool_value<
    Spec, std::void_t<decltype(Spec::mode_a_active_sources_in_sidx_pool)>>
    : std::integral_constant<bool, Spec::mode_a_active_sources_in_sidx_pool> {};
template <typename Spec>
constexpr bool nlr_spec_active_sources_in_sidx_pool_v =
    nlr_spec_active_sources_in_sidx_pool_value<Spec>::value;

/* Compile-time contract: a cached-SIDX spec must declare the active-in-pool trait. */
template <typename Spec>
constexpr bool nlr_spec_satisfies_source_pool_contract_v =
    (Spec::sidx_cache_kind == SidxCacheKind::None) ||
    nlr_spec_active_sources_in_sidx_pool_declared_v<Spec>;

/* Does this spec need the runner to stage explicit source positions for Mode A? */
template <typename Spec>
constexpr bool nlr_spec_needs_explicit_source_positions_v =
    (Spec::sidx_cache_kind != SidxCacheKind::None) &&
    !nlr_spec_active_sources_in_sidx_pool_v<Spec>;

/* Stage current P[active_i].Pos into `storage` and return it as a source_positions
 * array (layout pos[k*3+axis]) for specs that need explicit positions; returns
 * nullptr (keep the compact fast-path) otherwise. ~3 doubles/active, host-side. */
template <typename Spec>
static inline const double* nlr_stage_explicit_source_positions(
    const struct particle_data* P_host, const int* active_indices, int num_active,
    std::vector<double>& storage)
{
    if (nlr_spec_needs_explicit_source_positions_v<Spec> && num_active > 0) {
        storage.resize((size_t)num_active * 3);
        for (int k = 0; k < num_active; k++) {
            const int i = active_indices[k];
            storage[(size_t)k * 3 + 0] = (double)P_host[i].Pos[0];
            storage[(size_t)k * 3 + 1] = (double)P_host[i].Pos[1];
            storage[(size_t)k * 3 + 2] = (double)P_host[i].Pos[2];
        }
        return storage.data();
    }
    (void)P_host; (void)active_indices;
    return nullptr;
}

/* Remote-helper evaluation mode (SSOT guardrail). Replaces the
 * boolean ORACLE template parameter on mode_b_remote_evaluate_into_buffer.
 * (Distinct from RemoteEvalMode above, which is the Mode B comm strategy
 * enum for callers; this one is the *helper*-internal mode selector.)
 *
 *   Production       — tree walk only; writes tree result into accums_out.
 *                      Production iter dispatch + non-iter non-oracle.
 *   OracleCompare    — tree + brute + internal compare (emit per-slot
 *                      mismatches inline). Writes tree result into
 *                      accums_out. Non-iter oracle wrapper only.
 *   OracleBrutePass  — brute walk only against caller-owned brute-pass-
 *                      guarded ctx; writes brute result into accums_out.
 *                      Iterative oracle dispatch — caller-side per-iter
 *                      4-thing compare handles divergence detection. */
enum class RemoteHelperMode : int {
    Production       = 0,
    OracleCompare    = 1,
    OracleBrutePass  = 2,
    /* OracleIterative: collect tree + brute candidate sets in ONE shared
     * epoch (pre-drift both, drift union, evaluate brute-first then tree).
     * Outputs tree -> accums_out, brute -> accums_oracle_out.  Used by the
     * iterative remote oracle combined dispatch to guarantee epoch
     * equivalence between production and oracle trajectories. */
    OracleIterative  = 3,
};

/* ============================================================================
 * Iteration driver contract
 *
 * Spec opts in by declaring `using IterControl = Iterative;`. NotIterative
 * Specs (every existing sink loop today) declare nothing extra and the runner's
 * iteration machinery is compile-time elided.
 *
 * ALL EXISTING SPECS DECLARE `using IterControl = NotIterative;` AND ARE
 * UNCHANGED. The Iterative path is greenfield for iterative ports.
 * ========================================================================== */

struct NotIterative {};                /* default marker; existing Specs */
struct Iterative   {};                 /* opt-in marker for h-iteration loops */

/* Per-active iteration outcome returned by Spec::after_iter (host-side, post
 * AccumData reduction, before apply_active_writeback). */
enum class IterStatus : int {
    Converged    = 0,                  /* this active is done; runner removes from active_set */
    NeedsMore    = 1,                  /* re-iterate; do NOT change radius */
    AdjustRadius = 2,                  /* re-iterate WITH the new radius below */
};

struct IterResult {
    IterStatus status;
    double     new_h_search;           /* honored only when status == AdjustRadius */
};

/* Empty marker for Specs that don't carry per-active state across iterations.
 * Iterative Specs MUST declare `using IterScratch = ...;` explicitly — there
 * is no magic default — but the empty form is `using IterScratch = NoIterScratch;`. */
struct NoIterScratch {};

/* SupportsSubgroups (EXPLICIT trait):
 *
 * Multi-subgroup-using Specs (ags_density partitions actives by j_type_bitmask
 * and runs one kernel per bm group within each iter; preserves legacy nesting
 * `for iter` outer / `for bm` inner) declare `using SupportsSubgroups = std::true_type;`
 * in the Spec body. Single-subgroup Specs (density, harness) declare nothing —
 * default is `std::false_type`. The CALLER passes a 1-element `subgroups[]`
 * list whose `j_type_bitmask` matches the loop's natural neighbor type mask.
 * (If a future helper synthesizes the 1-element list from base
 * `neighbor_loop_args` for callers who don't want to hand-fill it, that's a
 * separate API.)
 *
 * Runner runtime check (in the iterative driver):
 * if `args.num_subgroups > 1 && !Spec::SupportsSubgroups::value`, hard-abort
 * with a loud message naming the loop. This prevents silent contract drift
 * (a Spec accidentally entering with multi-subgroup args without intending
 * to support the legacy bm-loop semantics). */

/* NlrSubgroup — one entry per j_type_bitmask group within an iterative call.
 *
 * Caller-filled (mirrors legacy gravity/ags_rkern.cc:115-145 partition + line
 * 128 Allreduce-BOR on global_bm_presence). For single-subgroup loops
 * (density, synthetic harness), caller passes a 1-element list with the
 * loop's natural neighbor_type_mask.
 *
 * Collective-symmetry invariant: subgroups[] length and
 * ordering — meaning each subgroups[i].j_type_bitmask value at each index i —
 * MUST be identical across all ranks. A rank may have num_active_local = 0
 * for some entries (its empty bm keys), but the entry MUST be present so the
 * per-iter `for sg in subgroups[]` loop fires the same number of
 * ghost_writeback_* collectives on every rank. Caller-side empty short-
 * circuit: if the global active total across the union is 0, caller MUST
 * skip run_neighbor_loop_iterative entirely; entering with num_subgroups = 0
 * is a caller bug (runner asserts num_subgroups >= 1 at entry). */
struct NlrSubgroup {
    unsigned int j_type_bitmask;       /* per-subgroup neighbor type mask */
    int         *active_indices;       /* host-staged; sub_actives for this bm */
    int          num_active_local;     /* may be 0 on this rank */
    /* Runner-derived per call (filled by NlrIterDriver, NOT by the caller):
     *   per-sub_active radii UVM, active_set UVM, per-subgroup CSR cache
     *   on Mode A, oracle scratch arrays. Out-of-band from this struct. */
};

/* ============================================================================
 * Empty marker types — required to be DISTINCT named types (not just any
 * empty struct) so static_asserts can verify exact intent.
 * ========================================================================== */

struct NoScatter {};
struct NoAccum   {};

/* ============================================================================
 * DispatchPath — runtime choice
 * ========================================================================== */

enum class DispatchPath : int {
    ModeA_GPU_NGL    = 0,
    ModeB_HostWalker = 1,
    Brute_Oracle     = 2,
};

/* ============================================================================
 * NeighborLoopSpec — copy-pasteable skeleton for a new physics loop
 *
 * Copy this block verbatim into your loop's <loop_name>_loop.h, replace the
 * physics, and you're done. The runner contract is what's below; everything
 * else is engine. See sinks/sink_env1_loop.h for a worked example.
 *
 * Per-spec layout convention (mirrors sink_env1_loop.h):
 *   PHYSICS BLOCK   — edit when changing the loop's physics
 *   ENGINE APPARATUS — touch only when changing the runner contract itself
 *   DIAGNOSTICS     — env-gated (see GIZMO_*_DIAG / _ORACLE / _SPIKE_*)
 *
 *   struct MyLoopSpec {
 *     // ============ PHYSICS BLOCK ============
 *
 *     // (1) Loop identity (drives env-var prefixes + diagnostic labels)
 *     static constexpr const char *loop_name = "my_loop";
 *
 *     // (2) Search policy
 *     static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
 *     static constexpr unsigned int            neighbor_type_mask = (1u<<0);
 *     static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;
 *
 *     // (3) Writeback policy
 *     static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
 *     static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::AllTypes;
 *     // REQUIRED for cached-SIDX specs (sidx_cache_kind != None): are all active
 *     // sources pool members? true = keep compact fast-path; false = runner stages
 *     // explicit P[].Pos (non-pool actives, e.g. a sink in a GasOnly loop).
 *     static constexpr bool mode_a_active_sources_in_sidx_pool = true;
 *
 *     // (4) Tolerances
 *     static constexpr double accum_tolerance = 1e-9;
 *
 *     // (5) Active-particle predicate. Caller passes this to nlr_build_active_list.
 *     static bool is_active(int particle_index);
 *
 *     // (6) Per-pair physics types (all trivially copyable; see TRAP 5)
 *     struct CallScalars   { NlrCommonScalars common; ... };
 *     struct ActiveData    { ...; CallScalars scalars; };
 *     struct NeighborData  { const struct particle_data* neighbor_particle; ... };
 *     using  AccumData     = my_loop_accum_t;
 *     using  ScatterData   = NoScatter;        // or your scatter type
 *     using  IdentityFields = NoIdentity;      // see runner header for non-default
 *     using  IterControl   = NotIterative;
 *
 *     // (6b) DeviceContext — REQUIRED. Most Specs declare base unchanged:
 *     using  DeviceContext = NeighborLoopDeviceContextBase;
 *     // To extend (e.g., Specs that need host-staged per-active UVM input
 *     // arrays visible inside load_active / load_neighbor on device):
 *     //   struct MyLoopDeviceContext : NeighborLoopDeviceContextBase {
 *     //       const my_per_active_t *per_active_local;  // UVM
 *     //   };
 *     //   using DeviceContext = MyLoopDeviceContext;
 *     // Then provide:
 *     //   static void populate_device_context(const neighbor_loop_args& args,
 *     //                                       DeviceContext& ctx);
 *     //   static void cleanup_device_context (const neighbor_loop_args& args,
 *     //                                       DeviceContext& ctx);   // optional;
 *     //                                       only required if populate
 *     //                                       allocated owning resources.
 *     // Compile-time invariants enforced by static_assert in the runner:
 *     //   - DeviceContext publicly derives from (or IS) NeighborLoopDeviceContextBase
 *     //   - DeviceContext is trivially copyable (captured by value into Kokkos lambdas)
 *
 *     // (7) Per-active aux passed by caller through neighbor_loop_args::aux
 *     struct Aux { my_loop_accum_t *out_buffer; };
 *
 *     // (8) Per-pair physics body. Forward to a header-inline single-source-of-truth helper if
 *     //     you want the same body callable from device (Mode A) and host
 *     //     (Mode B walker, oracle). See TRAP 3.
 *     KOKKOS_INLINE_FUNCTION
 *     static void pair_kernel(const ActiveData& active, const NeighborData& neighbor,
 *                             AccumData& accum, ScatterData& scatter);
 *
 *     // (9) Per-active and per-call hooks
 *     static double      search_radius(const neighbor_loop_args& args,
 *                                      int active_slot, int i);
 *     static CallScalars populate_call_scalars(const neighbor_loop_args& args);
 *     KOKKOS_INLINE_FUNCTION
 *     static ActiveData  load_active(const NeighborLoopDeviceContextBase& ctx,
 *                                    int active_slot, int i,
 *                                    double h_search, const CallScalars& scalars);
 *     KOKKOS_INLINE_FUNCTION
 *     static NeighborData load_neighbor(const NeighborLoopDeviceContextBase& ctx,
 *                                       int j, const IdentitySidecar& id,
 *                                       const ActiveData& active);
 *     KOKKOS_INLINE_FUNCTION static void zero_accum(AccumData& accum);
 *     static void apply_active_writeback(const neighbor_loop_args& args,
 *                                        int active_slot, int i,
 *                                        const AccumData& accum);
 *     static void merge_accum(AccumData& local_accum,
 *                             const AccumData& peer_accum);
 *
 *     // (10) Ghost-side writeback hooks (often empty). Trait-gated.
 *     static constexpr bool uses_ghost_write_detector = true;
 *     static constexpr bool uses_ghost_writeback      = false;
 *     static void ghost_write_detector_begin(const neighbor_loop_args&,
 *                                             const NeighborLoopPlan&);
 *     static void ghost_write_detector_end  (const neighbor_loop_args&,
 *                                             const NeighborLoopPlan&);
 *     // ...writeback variants as needed
 *
 *     // ============ ITERATIVE BLOCK (only when IterControl = Iterative) ============
 *     //
 *     // ALL of the following REQUIRED for `using IterControl = Iterative;` Specs;
 *     // ALL of the following ABSENT for `using IterControl = NotIterative;` Specs.
 *     //
 *     // using IterScratch — per-active POD that PERSISTS across iterations of a
 *     //   single iterative call. Empty form: `using IterScratch = NoIterScratch;`
 *     //   (no magic default; explicit declaration required). See TRAP 8.
 *     using  IterScratch = NoIterScratch;       // or struct IterScratch { ...; };
 *
 *     // max_iters — hard cap. Runner default behavior on exceedance is
 *     //   endrun(1155) + diagnostic. Override via on_max_iter_exceeded if a
 *     //   port's legacy policy differs (rare; density / ags_density already
 *     //   endrun in legacy).
 *     static constexpr int max_iters = 64;
 *
 *     // mode_a_csr_buffer_factor — Mode A oversized-CSR buffer factor for
 *     //   h-iteration reuse. REQUIRED (no default magic); each port states
 *     //   explicitly. Density legacy = 1.3 (DENSITY_H_BUFFER_FACTOR). Synthetic
 *     //   harness uses 1.05 to force rebuild paths. Static_assert > 1.0.
 *     static constexpr double mode_a_csr_buffer_factor = 1.3;
 *
 *     // SupportsSubgroups (OPTIONAL, default std::false_type) — declare
 *     //   `std::true_type` if the loop legitimately partitions actives by
 *     //   j_type_bitmask (ags_density). Single-subgroup loops (density,
 *     //   harness) omit this. Runner runtime check: num_subgroups > 1
 *     //   without this trait = hard abort.
 *     // using SupportsSubgroups = std::true_type;          // ags_density only
 *
 *     // after_iter — REQUIRED. Runs HOST-side per active after AccumData
 *     //   reduction, before apply_active_writeback. Returns IterResult to
 *     //   tell the runner what to do next for this active.
 *     static IterResult after_iter(const AfterIterContext<MyLoopSpec>& ctx,
 *                                   const AccumData& accum);
 *
 *     // after_iter_global — OPTIONAL. Per-rank post-iter side-effect hook.
 *     //   ===== HOST-ONLY, NO-MPI, LOCAL-RANK-ONLY (TRAP 7) =====
 *     //   Use case: ags_density's NeedToWakeupParticles_local OR-from-counter.
 *     //   Default if absent: runner-supplied no-op.
 *     // static void after_iter_global(const neighbor_loop_args& args,
 *     //                               const NlrIterDriver<MyLoopSpec>& drv);
 *
 *     // on_max_iter_exceeded — OPTIONAL. Default = runner-supplied
 *     //   endrun(1155) + diagnostic. Override only when port's legacy policy
 *     //   differs and Phil has approved.
 *     // static void on_max_iter_exceeded(const NlrIterDriver<MyLoopSpec>& drv);
 *
 *     // ============ ENGINE APPARATUS ============
 *     // (DeviceContext declaration is part of the REQUIRED block above —
 *     //  see (6b) earlier in this skeleton.)
 *
 *     // Optional: per-loop dispatch-threshold constants. These are
 *     // code-level dispatch policy for the loop (chosen alongside the
 *     // physics), NOT user-tunable runtime parameters. Default 64/64
 *     // when omitted. Tester env overrides exist but are not production
 *     // interface; see runner.cc TESTERS' KNOBS block.
 *     // static constexpr int modeb_threshold_sum = 64;
 *     // static constexpr int modeb_threshold_max = 64;
 *
 *     // ============ DIAGNOSTICS (env-gated) ============
 *     static double compare_accum(const AccumData& local, const AccumData& oracle);
 *     static void   diagnostic_dump_active(const ActiveDumpView<MyLoopSpec>& v);
 *     static void   diagnostic_dump_neighbor_list(const NeighborListDumpView<MyLoopSpec>& v);
 *   };
 *
 * ============================================================================
 * TRAPS — invariants a new Spec author will silently violate without these
 * ============================================================================
 *
 *   TRAP 1: populate_call_scalars is host->device VALUE CAPTURE, not
 *           inter-rank sync. Globals enter the kernel via this function only
 *           — the pair_kernel must NEVER read All.* directly. Threading
 *           CallScalars through ActiveData is what makes Mode A and Mode B
 *           bit-identical (single source of truth for per-call state).
 *
 *   TRAP 2: load_active runs on the device. KOKKOS_INLINE_FUNCTION is the
 *           contract: no std::, no host-only helpers. Same for load_neighbor,
 *           pair_kernel, zero_accum.
 *
 *   TRAP 3: pair_kernel body MUST be header-inlined (in <loop>_loop.h). A
 *           definition in the .cc breaks Kokkos device compilation with
 *           obscure linkage errors.
 *
 *   TRAP 4: merge_accum per-field op MUST match the per-field op pair_kernel
 *           writes. Sum for additive fields, MAX for max-reduced fields.
 *           Mismatch = silent multi-rank corruption (no compile error; only
 *           the oracle catches it).
 *
 *   TRAP 5: CallScalars, ActiveData, NeighborData, AccumData must all be
 *           std::is_trivially_copyable_v. Adding a std::vector member compiles
 *           fine and breaks transport silently in Mode B remote. Runner
 *           static_asserts at the explicit-instantiation site.
 *
 *   TRAP 6: diagnostic_dump_neighbor_list is called once per active, in slot
 *           order 0..num_active-1. State machines that rely on slot ordering
 *           (e.g. first-call-only gates) are safe.
 *
 *   TRAP 7: after_iter_global is HOST-ONLY, NO-MPI, LOCAL-RANK-ONLY. NEVER
 *           issue MPI inside it (no MPI_Allreduce / Bcast / Send / Recv /
 *           collective of any kind) — the runner issues the per-iter
 *           Allreduce AFTER the hook returns. Calling MPI inside risks
 *           deadlock and breaks the runner's collective-symmetry contract.
 *           Doc-only enforcement (no static check possible).
 *           Cross-rank iter symmetry IS preserved:
 *           every rank enters the same outer iter count for the call,
 *           INCLUDING empty ranks, until the runner's Allreduce on total
 *           active count breaks the loop. The hook fires once per outer
 *           iter on every rank. What VARIES per rank is the local active
 *           count and which slots are still in the active_set — NOT the
 *           iter count itself. Read those local-only quantities freely;
 *           don't broadcast / Allreduce them inside the hook.
 *
 *   TRAP 8: IterScratch is the ONLY per-active state that PERSISTS across
 *           iterations of a single iterative call. AccumData is byte-zeroed
 *           pre-each-iter via Spec::zero_accum (same channel as today).
 *           radii[slot] is mutated by the runner only (from Spec::search_radius
 *           at iter 0 / from IterResult::new_h_search on AdjustRadius);
 *           Specs MUST NOT touch radii directly. IterScratch is host-only,
 *           visible inside after_iter and after_iter_global; NOT inside
 *           pair_kernel. Conflating AccumData and IterScratch is the classic
 *           density-iteration quiet bug; the type system prevents it here.
 *
 *   TRAP 9: SupportsSubgroups is OPT-IN EXPLICIT (default false_type).
 *           Multi-subgroup Specs (ags_density partitions actives by
 *           j_type_bitmask) MUST declare `using SupportsSubgroups = std::true_type;`
 *           — the runner runtime-aborts loudly if it sees num_subgroups > 1
 *           without the trait. Single-subgroup Specs (density, harness)
 *           omit the trait; the CALLER passes a 1-element subgroups[] list
 *           with the loop's natural neighbor type mask.
 *           Caller-side: subgroups[] length AND ordering MUST match across
 *           ranks; per-rank num_active_local = 0 for some entries is OK
 *           (collective-symmetry filler), but the entry MUST be present.
 *           Caller-side empty short-circuit: if global active total = 0,
 *           caller MUST skip the runner entirely.
 * ========================================================================== */

/* ============================================================================
 * Public runner entry point
 *
 * Linking contract: function template; implementation in
 * neighbor_loop_runner.cc with EXPLICIT INSTANTIATION per migrated caller:
 *   template void run_neighbor_loop<MyLoopSpec>(const neighbor_loop_args&);
 * Forgetting the instantiation = clean linker error (not silent template bloat).
 * ========================================================================== */

/* ============================================================================
 * neighbor_loop_args — physics inputs to run_neighbor_loop<Spec>.
 *
 * The caller fills these with PHYSICS facts ONLY. The caller does not pick
 * the dispatch path, run any prep, or know which substrate (Mode A GPU NGL,
 * Mode B local, Mode B remote, future paths) will execute. The runner owns
 * all of that.
 *
 * Field semantics:
 *   P, CellP        — particle data view at call time. The runner takes a
 *                     local working copy; on paths that mutate the global
 *                     particle arrays (today: only Mode A's ghost-import
 *                     prep), the runner refreshes its working copy after
 *                     mutation. Callers MUST NOT pass pointers that go
 *                     stale across the call.
 *   num_total       — = NumPart at call time. The runner refreshes this
 *                     internally on paths that grow NumPart via ghost
 *                     import (Mode A).
 *   active_list     — caller-owned active-particle indices, valid for the
 *                     duration of the call. args.active_list[slot] =
 *                     particle index.
 *   num_active      — length of active_list.
 *   aux             — Spec-defined POD pointer for per-active side arrays
 *                     (e.g. SinkEnv1Spec::Aux::per_active_accum). Recover
 *                     the typed pointer with nlr_aux<Spec>(args).
 *   ghost_safety_factor — caller fills unconditionally from
 *                     gizmo_ghost_safety_factor(); the runner uses it only
 *                     on paths that import ghosts (Mode A today). Mode B
 *                     paths ignore it.
 *
 * The caller does NOT compute or stage per-active radii. The runner stages
 * radii once via Spec::search_radius and reuses them for any prep
 * (Mode A) and the chosen walker.
 * ========================================================================== */
/* NlrForceMode — A/B dispatch override values. Defined here (early in the
 * header) rather than near the env-var doc block below because struct
 * neighbor_loop_args' in-class default initializer for `dispatch_override`
 * needs the enum to be visible at that point. The companion env-var
 * resolution functions (gizmo_nlr_force_mode, gizmo_nlr_force_mode_for)
 * are declared below near their documentation. */
enum class NlrForceMode { None = 0, A = 1, B = 2 };

/* External symmetric gas-gas CSR injection (hydro corridor support:
 * cellcorrections/gradients/hydro_force share one symmetric gas CSR per
 * step). When the caller has already built a symmetric gas CSR (e.g. the
 * gizmo_gradients_prep_symlist) that's intended to be SHARED across multiple
 * runner-Spec consumers in one step, it can be injected here instead of the
 * runner re-building it per call. The runner stages the host CSR into
 * Kokkos-managed memory for device dispatch but does NOT free the caller's
 * host buffers; the caller owns CSR lifetime across the whole corridor.
 *
 * Contract requirements (UNCONDITIONALLY enforced at runtime via endrun in
 * run_mode_a — external CSR is a sharp tool, violations cause silent
 * wrong-particle writeback because d_accums[aa] is later applied to
 * args.active_list[aa], not to ec->active_indices[aa]):
 *   - active_indices[] equals args.active_list[] elementwise (endrun 7308)
 *   - num_active equals args.num_active (endrun 7301)
 *   - offsets[0] == 0 and offsets[num_active] == total_pairs (7306/7307)
 *   - offsets monotonic non-decreasing per row (7309)
 *   - Symmetric pair structure: for every (a -> j), there exists (b -> i)
 *     where b is the home slot of j and i is in j's neighbor list
 *     (not runtime-checkable — caller-asserted invariant)
 *   - SidxCacheKind::GasOnly Specs only — other kinds rejected (endrun 7300)
 *   - active_indices / offsets non-null; neighbors non-null iff total_pairs>0
 *     (endrun 7302/7303/7304/7305)
 *
 * GHOST-POOL OWNERSHIP (load-bearing): a non-null external_csr means the CSR's
 * neighbor indices refer to the CALLER'S live particle+ghost pool — the caller
 * imported those ghosts and owns their lifetime. The runner therefore does NOT
 * run its own ghost import (a fresh import could renumber ghost slots, leaving
 * the CSR pointing at wrong particles) and does NOT run ghost_exchange_cleanup
 * (the pool outlives this call; the caller tears it down at the end of its
 * span). Ghost writeback + write-detector hooks still run normally — they read
 * the caller's live provenance maps. Enforced at runtime (endrun) at NTask>1:
 * ghost_pool_is_live() + non-null home_rank/home_index provenance +
 * args.num_total == NumPart (args must be built AFTER the caller's import).
 *
 * Mode B ignores external_csr entirely (request-driven walker does not need
 * a CSR). Only consumed by Mode A; dispatch to a Mode B path with a non-null
 * external_csr is a caller error (rejected at runtime). */
struct nlr_external_csr {
    int        *active_indices;  /* len num_active; MUST equal args.active_list */
    int         num_active;
    int64_t    *offsets;         /* len (num_active+1); 64-bit row pointers */
    int        *neighbors;       /* len total_pairs; 32-bit local indices */
    int64_t     total_pairs;
    const char *owner_name;      /* diagnostic, e.g. "corridor:gradients" */
};

struct neighbor_loop_args {
    struct particle_data *P;
    struct gas_cell_data *CellP;
    int    num_total;
    int   *active_list;          /* args.active_list[slot] = particle index */
    int    num_active;
    void  *aux;                  /* spec-defined POD for per-active side arrays;
                                  * recover the typed pointer with nlr_aux<Spec>(args). */
    double ghost_safety_factor;  /* caller fills via gizmo_ghost_safety_factor()
                                  * unconditionally; runner uses ONLY on paths
                                  * that call gizmo_request_filtered_ghost_import. */
    unsigned int neighbor_type_mask_override;
                                 /* 0 (default — set by nlr_default_args) => the
                                  * non-iterative runner uses Spec::neighbor_type_mask.
                                  * Non-zero => that mask is used for THIS call instead.
                                  * Lets a non-iterative caller supply a runtime
                                  * neighbor-type mask without the iterative runner's
                                  * subgroup machinery (dm_fuzzy DMGrad: per-bm calls
                                  * with bm = ags_gravity_kernel_shared_BITFLAG). The
                                  * iterative runner is unaffected — it uses
                                  * NlrSubgroup::j_type_bitmask. */
    const nlr_external_csr *external_csr = nullptr;
                                 /* nullptr (default) => runner builds its own CSR
                                  * via gpu_ngb_list_build AND owns its own ghost
                                  * import + cleanup (legacy behavior).
                                  * non-null => Mode A skips the build, stages the
                                  * provided host CSR into Kokkos memory, and the
                                  * CALLER OWNS THE GHOST POOL: the runner skips
                                  * its ghost import and cleanup (see the
                                  * GHOST-POOL OWNERSHIP contract on struct
                                  * nlr_external_csr above). */
    NlrForceMode dispatch_override = NlrForceMode::None;
                                 /* None (default) => normal priority: env vars
                                  * GIZMO_<LOOPNAME>_FORCE_MODE > GIZMO_NLR_FORCE_MODE
                                  * > adaptive threshold. Non-None => per-call
                                  * override winning ABOVE env vars (corridor mode
                                  * decision in hydro/hydro_corridor.cc uses this
                                  * to enforce coherent Mode A/B across all
                                  * corridor consumers — cellcorrections/gradients/
                                  * hydro_force as they are runner-ported). */
};

/* Effective neighbor-type mask for a non-iterative run_neighbor_loop call:
 * args.neighbor_type_mask_override when set, else the Spec constexpr. MUST be
 * used at EVERY non-iterative mask site (ghost import, Mode A CSR build, Mode B
 * candidate collection, remote eval, oracle/brute) so production and oracle
 * agree on the neighbor set. */
static inline unsigned int
nlr_effective_neighbor_type_mask(const neighbor_loop_args& args,
                                 unsigned int spec_mask)
{
    return args.neighbor_type_mask_override ? args.neighbor_type_mask_override
                                            : spec_mask;
}

/* Iterative-call args — extends neighbor_loop_args with subgroup batch.
 *
 * Used only for `using IterControl = Iterative;` Specs. NotIterative loops
 * pass plain `neighbor_loop_args` and the runner never sees a subgroups list.
 *
 * The caller (e.g., the gravity/ags_rkern.cc port) builds subgroups[]
 * from the global_bm_presence union via Allreduce-BOR, mirroring legacy
 * gravity/ags_rkern.cc:115-128 exactly. For loops without bm partitioning
 * (the hydro/density.cc port; the synthetic harness), the caller
 * passes a 1-element subgroup list whose j_type_bitmask matches the loop's
 * natural neighbor_type_mask.
 *
 * active_list / num_active in the base struct are the UNION across all
 * subgroups — used for caller-side helpers (nlr_build_active_list etc.).
 * Per-subgroup actives live inside subgroups[i].active_indices /
 * num_active_local. The runner does NOT use the union list for dispatch;
 * it walks the subgroups list. */
struct neighbor_loop_args_iterative : neighbor_loop_args {
    NlrSubgroup *subgroups;          /* host array, length num_subgroups */
    int          num_subgroups;      /* >= 1; caller short-circuits if 0 globally */
};

/* ============================================================================
 * Caller-side helpers (collapse boilerplate; physics stays visible)
 *
 * Pattern for the caller of run_neighbor_loop<Spec>:
 *
 *   int *active_list, num_active, num_global_active;
 *   if (!nlr_build_active_list(MyLoopSpec::is_active, &active_list,
 *                              &num_active, &num_global_active,
 *                              "myloop_active")) {
 *       return;                                       // no active particles anywhere
 *   }
 *   MyLoopAccum *out_buffer = (MyLoopAccum*)mymalloc("myloop_out", ...);
 *   MyLoopSpec::Aux aux{out_buffer};
 *
 *   neighbor_loop_args args = nlr_default_args();
 *   args.active_list = active_list;
 *   args.num_active  = num_active;
 *   args.aux         = &aux;
 *   run_neighbor_loop<MyLoopSpec>(args);
 *
 *   // ...physics scatter from out_buffer into your per-loop temp struct...
 *
 *   myfree(out_buffer);
 *   nlr_free_active_list(active_list);
 * ========================================================================== */

/* Recover the typed Spec::Aux* from neighbor_loop_args.aux, replacing
 * hand-written reinterpret_cast at every hook site. */
template <typename Spec>
typename Spec::Aux *nlr_aux(const neighbor_loop_args& args)
{
    return reinterpret_cast<typename Spec::Aux*>(args.aux);
}

/* Count active particles via predicate, MPI_Allreduce, malloc + fill.
 * Returns false if global count is zero (caller should early-return);
 * otherwise *list_out is mymalloc'd and must be released via
 * nlr_free_active_list. malloc_label is forwarded to mymalloc. */
template <typename ActivePredicate>
bool nlr_build_active_list(ActivePredicate is_active,
                           int **list_out,
                           int  *num_local_active_out,
                           int  *num_global_active_out,
                           const char *malloc_label);

void nlr_free_active_list(int *active_list);

/* Fills the engine fields of neighbor_loop_args from current globals
 * (P, CellP, NumPart, ghost_safety_factor). Caller fills active_list,
 * num_active, aux. Out-of-line so changes to global plumbing are
 * caller-invisible. */
neighbor_loop_args nlr_default_args(void);

/* ============================================================================
 * NeighborLoopPlan — runner-internal execution policy descriptor.
 *
 * `path` is the single-source-of-truth. Properties (needs imported ghosts, may acquire GPU
 * arena, may mutate NumPart, etc.) are derived via the nlr_path_*() helpers
 * below — never stored as separate boolean fields on the plan, so a new
 * path adds cases to the helpers and never new fields to this struct.
 * ========================================================================== */
struct NeighborLoopPlan {
    enum class Path {
        ModeA_GpuNgl,    /* GPU NGL pipeline. Current substrate uses imported
                          * ghosts; freshness is satisfied by today's import
                          * helper (gizmo_request_filtered_ghost_import). The
                          * substrate is what the path REQUIRES; the helper is
                          * one (current) way of satisfying that requirement
                          * — future work may swap helpers without changing
                          * the path. Do not equate Mode A with "global drift". */
        ModeB_Local,     /* Host walker on local pool; no global mutation. */
        ModeB_Remote     /* peer-to-peer request/reply; no global mutation. */
        /* Future paths: add cases here AND to nlr_path_*() predicates below. */
    };
    Path path;

    /* Global active-particle count after MPI_Allreduce. Sentinel: -1 means
     * "not populated" — happens on force-mode paths
     * (GIZMO_NLR_FORCE_MODE{A,B}=1) when GIZMO_PHASE0_DIAG is OFF, since the
     * runner skips the Allreduce in that case to avoid an extra collective.
     * Hooks/consumers MUST check for -1 before assuming this field reflects
     * a real global count. Threshold-dispatch paths and any path with
     * PHASE0_DIAG on always populate it. */
    int  num_active_global;
};

/* Path-derived predicate helpers. Single source of truth: keyed on path. */
bool nlr_path_uses_imported_ghosts(NeighborLoopPlan::Path path);
bool nlr_path_uses_gpu_arena(NeighborLoopPlan::Path path);
bool nlr_path_permits_global_numpart_mutation(NeighborLoopPlan::Path path);
bool nlr_path_uses_lazy_drift(NeighborLoopPlan::Path path);

/* Stable string label for diagnostics (PHASE0_NLR `path=...` field).
 * Stable across the runner's lifetime; future path cases extend this. */
const char *nlr_path_label(NeighborLoopPlan::Path path);

/* ============================================================================
 * run_neighbor_loop<Spec>(args) — the runner.
 *
 * Numerical execution policy single-source-of-truth for one neighbor loop. The caller hands
 * over physics; the runner picks the substrate and runs it.
 *
 *   CALLER OWNS                       RUNNER OWNS
 *   -----------                       -----------
 *   Active set                        Selected execution path
 *   Search predicate (Spec traits)    Whether ghost import runs
 *   Per-loop physics data (Spec)      Drift strategy
 *   Write pattern (Spec trait)        Transport strategy
 *   Ghost safety factor (args field)  Per-active radii staging
 *                                     Lifecycle hooks (detector, writeback)
 *                                     Hard-corridor invariant enforcement
 *
 * Execution path is chosen by the runner from a small enum (Mode A GPU NGL,
 * Mode B local, Mode B remote, future paths). Selection precedence: env
 * force-mode > threshold dispatch on
 * num_active_global vs Spec::modeb_threshold_{sum,max} > Spec defaults.
 * Future paths add a case to the dispatch + the path-predicate helpers
 * (nlr_path_uses_imported_ghosts, nlr_path_uses_gpu_arena, etc.) and never
 * add fields to the plan struct itself — path is the single-source-of-truth.
 *
 * Tiny-N hard-corridor invariants (HARD ABORT on violation, always-on,
 * every build, enforced via global lifecycle counters):
 *
 *   On Mode B (local OR remote) paths, from runner entry to runner exit:
 *     - NO move_particles call (g_global_drift_counter unchanged)
 *     - NO ghost_exchange_impl call (g_ghost_import_counter unchanged)
 *     - NO gpu_particles_arena_acquire call
 *       (g_gpu_arena_acquire_counter unchanged)
 *     - NO change to NumPart
 *
 *   These invariants protect the lazy-drift architecture: Mode B reads
 *   owner-local args.P[j] / args.CellP[j] freely during evaluation and
 *   drifts only touched candidates via mode_b_lazy_drift_candidates. The
 *   invariant is NO GLOBAL MUTATION, not no read. A full global drift in
 *   the prep layer would defeat lazy drift; the corridor enforcement
 *   prevents that regression.
 *
 *   Mode A path retains ghost import + global drift (today's substrate).
 *   That cost is appropriate for large-N where the import + drift amortizes
 *   over millions of pair evaluations. Forced Mode A on tiny-N is
 *   tester/baseline territory: expensive but correct.
 *
 * Path -> required-data matrix (today; future paths extend by adding cases
 * to nlr_path_*() predicates):
 *
 *   ModeA_GpuNgl       : may move_particles, may ghost_exchange_impl,
 *                        may grow NumPart, may acquire GPU arena
 *   ModeB_Local        : none of the above
 *   ModeB_Remote       : none of the above (uses peer-to-peer request/reply
 *                        substrate; peer-local walks; lazy drift on
 *                        touched candidates)
 *
 * Lazy-drift contract inside Mode B:
 *   collect candidates pre-drift -> mode_b_lazy_drift_candidates(touched)
 *   -> evaluate pairs post-drift. drift_particle is per-particle and
 *   short-circuits on time1 == time0. Repeat candidates between queries
 *   in the same call are a fast no-op.
 *
 * radii lifetime: the runner stages a std::vector<double> radii of size
 * args.num_active and passes radii.data() to path-specific functions and
 * (on Mode A) to gizmo_request_filtered_ghost_import_fresh. The pointer
 * is RUNNER-CALL-LIFETIME ONLY. Callees MUST NOT store the pointer; if a
 * callee needs to retain values, it copies. The vector goes out of scope
 * at the end of run_neighbor_loop.
 *
 * Spec lifecycle hooks (two channels, optional, SFINAE-detected;
 * default no-op; per-Spec gating decided per-channel by audit):
 *   detector  : Spec::ghost_write_detector_begin/end<Spec> — audit/debug.
 *   writeback : Spec::ghost_writeback_begin/end<Spec> — physics state
 *               propagation. Spec body is the per-flag #ifdef union of
 *               all enabled physics-flag writebacks for that loop.
 * Whether each channel fires on all paths or only on imported-ghost paths
 * is a per-Spec decision settled by reading the underlying call's actual
 * semantics; no blanket rule.
 *
 * For the full Spec contract (load_active, load_neighbor, pair_kernel,
 * accum/scatter writeback, optional iterative loop semantics), see the
 * Spec-shape comment block above this declaration.
 * ========================================================================== */
template <typename Spec>
void run_neighbor_loop(const neighbor_loop_args& args);

/* ============================================================================
 * Iterative entry point (only for `using IterControl = Iterative;`)
 *
 * Caller fills `neighbor_loop_args_iterative` (extends neighbor_loop_args with
 * subgroups[] + num_subgroups). Caller is responsible for the per-iter
 * collective symmetry of subgroups[] across ranks (length + ordering match;
 * see TRAP 9). NotIterative Specs MUST NOT call this entry;
 * they call run_neighbor_loop above.
 *
 * Definition in neighbor_loop_runner.cc.
 * ========================================================================== */
template <typename Spec>
void run_neighbor_loop_iterative(const neighbor_loop_args_iterative& args);

/* ============================================================================
 * NlrIterDriver<Spec> — runner-internal iterative state owner.
 *
 * Owns ALL per-call state for one `run_neighbor_loop_iterative<Spec>` call:
 *
 *   ============ LIFETIME / ZEROING CONTRACT ============
 *
 *   Per-active arrays are keyed by (subgroup_index, subgroup_slot). Every
 *   active particle belongs to exactly ONE subgroup; that subgroup's
 *   `active_set` and per-active arrays are the only place its state lives.
 *
 *   IterScratch (per-active POD)
 *       - Allocated ONCE per call, at iter 0, before any subgroup dispatch.
 *       - Byte-zeroed ONCE before iter 0.
 *       - PERSISTS across iterations — the runner NEVER zeros it between
 *         iters. Spec::after_iter mutates it freely; that's the carry-
 *         forward channel for bisection state, history, etc.
 *       - Freed at end of call.
 *
 *   AccumData (per-active POD)
 *       - Same key shape as IterScratch (subgroup_index, subgroup_slot).
 *       - Zeroed via `Spec::zero_accum` BEFORE EACH iteration's pair_kernel
 *         dispatch (NOT once per call). Resets per iter so each iter's
 *         pair_kernel sees a fresh accumulator. This is the classic
 *         density-iteration contract — stale AccumData across iters is
 *         the quiet bug the type system prevents.
 *       - Lives until end of call (final iter's contents are what
 *         apply_active_writeback reads).
 *
 *   radii (per-active double, UVM)
 *       - Initialized at iter 0 from Spec::search_radius.
 *       - Mutated by RUNNER ONLY, ONLY in response to a Spec::after_iter
 *         returning IterStatus::AdjustRadius with new_h_search. Specs
 *         MUST NOT touch radii directly. h_search_current in
 *         AfterIterContext is a value-snapshot read FROM radii at hook
 *         entry — NOT a reference, NOT mutable through the context.
 *
 *   active_set (per-subgroup int array, UVM)
 *       - Initialized at iter 0 to {0..num_active_local-1} for each subgroup.
 *       - Compacted (Converged slots removed) AFTER each iter's
 *         after_iter pass. Slot identity preserved: the int values are
 *         the subgroup-local slot indices, never re-keyed.
 *
 *   CallScalars (per-call POD)
 *       - Captured ONCE via `Spec::populate_call_scalars(args)` at call
 *         entry. The DRIVER OWNS the CallScalars value for the whole
 *         iterative call (NOT a stack temporary inside any subgroup
 *         loop). All hooks (pair_kernel via load_active, after_iter via
 *         AfterIterContext, after_iter_global via NlrIterDriver const&)
 *         see the SAME object lifetime.
 *
 *   ============ Mode-A-only fields (NOT touched on Mode B paths) ============
 *
 *   csr_cache       — per-subgroup cached oversized CSR + offset_lookup +
 *                     buffered_h. Built at iter 0; rebuilt on h-exceeds-buffer
 *                     trigger (mirrors hydro/density_gpu.cc:152-188 legacy).
 *
 *   ============ Oracle-only fields (only when GIZMO_NLR_ORACLE=1) ============
 *
 *   Separate radii_oracle / scratch_oracle / accum_oracle arrays so brute-
 *   force trajectory is independent of production. Brute dry-run runs FIRST
 *   per iter, production SECOND.
 * ========================================================================== */
template <typename Spec>
struct NlrIterDriver {
    /* Driver-owned for whole call. */
    const neighbor_loop_args_iterative& args;
    typename Spec::CallScalars          cs;          /* value-owned, NOT a reference */

    int iter_index;                                  /* 0-based; shared across all subgroups */
    int local_active_total;                          /* sum of local_active_per_sg; refreshed per iter */
    int global_active_total;                         /* sum of global_active_per_sg; -1 = "not populated yet" */

    /* Per-subgroup activity tracking:
     * legacy ags re-partitions bm_groups each outer iter; with fixed
     * subgroups[] we instead Allreduce-SUM per-subgroup local counts to get
     * global activity per subgroup. Globally-converged subgroups skip the
     * per-iter dispatch + collective. Single-subgroup case: length-1
     * vectors; behavior unchanged. */
    std::vector<int>                  local_active_per_sg;        /* [num_subgroups]; per-rank, post-compaction */
    std::vector<int>                  global_active_per_sg;       /* [num_subgroups]; Allreduce-SUM each iter */

    /* Driver-owned DeviceContext.
     *
     * `ctx` is a host-side trivially-copyable VALUE. Members may point to
     * UVM/device-visible storage (e.g., per-active host-staged UVM arrays
     * from Spec::populate_device_context, or arena-resident P/CellP set by
     * Mode A). The struct ITSELF is captured by VALUE into Kokkos device
     * lambdas, exactly like the per-call DeviceCtx in run_mode_a /
     * run_mode_b_*. `ctx` is NOT UVM-resident.
     *
     * Population is path-dependent and happens AFTER path selection:
     *   Mode B: initialize_device_context_mode_b() at iter-0 entry binds
     *           ctx.P = args.P, ctx.CellP = args.CellP, num_total = args.num_total.
     *   Mode A: initialize_device_context_mode_a_after_arena() AFTER
     *           gpu_particles_arena_acquire binds ctx.P = arena_P(),
     *           ctx.CellP = arena_CellP(). Body lands in step 2c.2.
     * Both methods call Spec::populate_device_context(args, ctx) if extended,
     * then set ctx_initialized = true.
     *
     * `ctx_initialized` guards destructor cleanup:
     * if init was stubbed/hard-aborted before populate_device_context ran,
     * destructor must NOT call cleanup_device_context. Only fire cleanup
     * when init actually completed. */
    typename Spec::DeviceContext ctx;
    bool                         ctx_initialized = false;

    /* Per-subgroup state (DYNAMICALLY SIZED to args.num_subgroups).
     *
     * Each std::vector<...> is host-side, sized at construction to
     * args.num_subgroups. The pointers it holds reference per-subgroup
     * UVM-allocated arrays; pair_kernel / Mode B walker / oracle paths see
     * the right addresses through these pointers.
     *
     * Lifetime: allocated by NlrIterDriver constructor, freed by destructor.
     * No fixed cap — global_bm_presence can span up to 64 bitmask values.
     *
     * Mode A iterative-specific fields (cached CSR + arena/session) and
     * Oracle-specific fields (separate scratch + radii arrays for brute
     * trajectory) are declared further below. */
    std::vector<typename Spec::IterScratch *> scratch_uvm;     /* [num_subgroups][num_active_local] */
    std::vector<typename Spec::AccumData   *> accum_uvm;       /* [num_subgroups][num_active_local] */
    std::vector<double *>                     radii_uvm;       /* [num_subgroups][num_active_local] */
    std::vector<int    *>                     active_set_uvm;  /* [num_subgroups][num_active_local], compacted */
    std::vector<int>                          active_set_size; /* [num_subgroups], shrinks on Converged compaction */

    /* Mode A iterative cached CSR/session state.
     * Allocated lazily on first Mode A iter dispatch; left empty on Mode B paths.
     * Lifecycle: arena_acquire ONCE per call (via acquire_arena_and_init_ctx_mode_a),
     * CSR built once per subgroup at iter 0, rebuilt on h-exceeds-buffer trigger,
     * all freed in driver destructor (passing SIDX to gpu_ngb_list_free so the
     * step-persistent SIDX cache survives — matches sink_env1/feed/swk idiom).
     *
     * CSR row-key invariant: mode_a_csr_offset_lookup[sg][slot]
     * is keyed on subgroup-slot-AT-BUILD-TIME, NOT on the (possibly compacted)
     * current active_set position. Compaction (Converged-slot removal) NEVER
     * rewrites old row identity; rebuild creates a new lookup and discards
     * the old one. */
    bool                              arena_acquired      = false;
    bool                              ghost_import_done   = false;     /* Mode A only; cleaned in destructor */
    std::vector<gpu_neighbor_list_t>  mode_a_cached_gnl;             /* [num_subgroups] */
    std::vector<int    *>             mode_a_csr_offset_lookup;      /* [num_subgroups][num_active_local]; UVM */
    std::vector<double *>             mode_a_csr_buffered_h;         /* [num_subgroups][num_active_local]; UVM */
    std::vector<bool>                 mode_a_csr_valid;              /* [num_subgroups] */

    /* ========================================================================
     * Iterative-oracle state.
     *
     * Allocated iff `oracle_enabled` is true at iter-0 entry AND the production
     * path is ModeB_HostWalker. Mode A iterative oracle is hard-stubbed at the
     * outer body (run_neighbor_loop_iterative) with an explicit
     * endrun + scaffolding-rationale message — oracles are temporary port-
     * validation scaffolding, the Mode A *production* path is validated by
     * the synthetic harness + Mode B oracle on the same Spec.
     *
     * Independent brute trajectory: own IterScratch + radii + accum +
     * active_set per subgroup. Same slot identity as production (driver
     * always references args.subgroups[sg].active_indices) so per-slot
     * compare is unambiguous. Per-iter compare emits four classes of
     * mismatch line on divergence between production and oracle:
     * AccumData / IterStatus / new_h_search / active-set membership.
     * All four share kMismatchPrintCap (=1024) per-call.
     * ====================================================================== */
    bool                                      oracle_enabled = false;
    std::vector<typename Spec::IterScratch *> scratch_oracle_uvm;       /* [num_subgroups][num_active_local] */
    std::vector<typename Spec::AccumData   *> accum_oracle_uvm;         /* same shape */
    std::vector<double *>                     radii_oracle_uvm;         /* same shape */
    std::vector<int    *>                     active_set_oracle_uvm;    /* same shape, compacted */
    std::vector<int>                          active_set_oracle_size;   /* [num_subgroups] */

    /* Per-iter audit instrumentation (4-thing compare).
     * Stored as int (IterStatus) UVM arrays via static_cast for trivially-
     * copyable transit. new_h_*_uvm holds AdjustRadius candidate; for
     * Converged/NeedsMore it holds the current radius (for symmetric init). */
    std::vector<int    *>                     status_prod_uvm;          /* [num_subgroups][num_active_local] */
    std::vector<int    *>                     status_oracle_uvm;        /* same shape */
    std::vector<double *>                     new_h_prod_uvm;           /* same shape */
    std::vector<double *>                     new_h_oracle_uvm;         /* same shape */

    long long                                 oracle_mismatch_count = 0;  /* per-call, shared across origin tags */
    /* TEMPORARY: set true by with_oracle dispatch helpers after comparing from
     * host vectors; b.compare uses this to skip the broken UVM-based accum
     * comparison. Remove once the oracle validation path is removed. */
    bool                                      oracle_accum_compared_in_dispatch = false;

    /* Independent ctx_oracle lifecycle. Own
     * populate_device_context + cleanup_device_context calls; NEVER aliased
     * to production ctx after init. Brute-pass toggle gated by RAII
     * NlrOracleBrutePassGuard<Spec> at every dispatch helper entry. */
    typename Spec::DeviceContext              ctx_oracle;
    bool                                      ctx_oracle_initialized = false;

    /* ========================================================================
     * Mode A diagnostic counters (step 3 harness instrumentation).
     *
     * Public observability for harness validation; production code should
     * never read these (zero overhead under default builds). Counts are
     * per-driver-instance (per iterative call), reset to 0 at construction.
     *
     *   csr_rebuild_count       — # times rebuild_mode_a_arena_and_ctx_for_current_active_union
     *                              completed (Mode A union rebuild fires).
     *   arena_acquire_count     — # times gpu_particles_arena_acquire fired
     *                              (iter-0 init + each union rebuild).
     *   csr_local_rebuild_count — # times a subgroup's CSR was locally
     *                              freed+rebuilt inside
     *                              nlr_iter_dispatch_subgroup_mode_a (per sg, summed).
     * ====================================================================== */
    int                                       csr_rebuild_count       = 0;
    int                                       arena_acquire_count     = 0;
    int                                       csr_local_rebuild_count = 0;

    /* effective_args: writable copy of the base
     * neighbor_loop_args slice. Mode A iter mutates num_total/P/CellP after
     * ghost import to reflect the imported-ghost pool (mirrors non-iter
     * effective_args refresh at runner.cc:2046-2051). Mode B leaves it
     * equal to the original base args (lazy-drift contract — Mode B reads
     * owner-local args.P[j] directly). All per-iter dispatch helpers read
     * effective_args.num_total / P / CellP, NOT args.num_total / P / CellP. */
    neighbor_loop_args                effective_args;

    /* Constructor (runner-private; only invoked inside run_neighbor_loop_iterative<Spec>):
     *   - Resizes per-subgroup vectors to args.num_subgroups.
     *   - Allocates per-subgroup UVM arrays sized to subgroups[sg].num_active_local.
     *   - Initializes radii from Spec::search_radius for every active in every subgroup.
     *   - Byte-zeros IterScratch ONCE at construction (persists across iters).
     *   - Fills active_set with {0..num_active_local-1} per subgroup; active_set_size = num_active_local.
     *   - AccumData NOT zeroed here — runner zeros via Spec::zero_accum at the start of each iter.
     *   - DOES NOT touch ctx (path-dependent init via initialize_device_context_*).
     * Definition in mesh/neighbor_loop_runner.cc. */
    explicit NlrIterDriver(const neighbor_loop_args_iterative& a,
                           const typename Spec::CallScalars& s);

    /* Destructor frees all UVM allocations + cleans DeviceContext IF
     * ctx_initialized. */
    ~NlrIterDriver();

    /* Path-specific DeviceContext init. Called AFTER path
     * selection; binds ctx.P/CellP/num_total + runs Spec::populate_device_context
     * (if extended) + sets ctx_initialized = true.
     *
     * Mode B: caller P/CellP pointers (lazy-drift contract — Mode B reads
     *         args.P[j]/CellP[j] directly).
     * Mode A: arena-resident P_gpu/CellP_gpu (must run
     *         AFTER gpu_particles_arena_acquire).
     *
     * Both methods are runner-private — only called by run_neighbor_loop_iterative
     * after path selection. */
    void initialize_device_context_mode_b();
    void initialize_device_context_mode_a_after_arena();

    /* Single combined Mode A iter-0 init.
     *   1. gpu_particles_arena_acquire ONCE per call (arena_acquired = true).
     *   2. initialize_device_context_mode_a_after_arena() — binds ctx.P/CellP
     *      to arena-resident pointers + populates extended DeviceContext.
     * Driver destructor matches: if arena_acquired, mark_clean once on exit. */
    void acquire_arena_and_init_ctx_mode_a();

    /* Self-sufficient rebuild
     * method. NO PARAMETERS. Driver builds the union from its own subgroup
     * state — current active_set per subgroup, oversized radii via
     * mode_a_csr_buffer_factor, mask = OR of all globally-active subgroups'
     * j_type_bitmask. As part of the lifecycle, invalidates ALL subgroup
     * CSR caches (pitfall 2: any rebuild makes cross-subgroup neighbor
     * indices stale). */
    void rebuild_mode_a_arena_and_ctx_for_current_active_union();

    /* Driver owns Kokkos shared allocations — no copy/move. */
    NlrIterDriver(const NlrIterDriver&) = delete;
    NlrIterDriver& operator=(const NlrIterDriver&) = delete;
};

/* ============================================================================
 * AfterIterContext<Spec> — passed to Spec::after_iter (host-side per active,
 * post-AccumData reduction, pre-apply_active_writeback).
 *
 * Slot identity: subgroup-local slots ARE the runner execution
 * slots. IterScratch / radii / AccumData / oracle arrays / final writeback
 * all keyed by (subgroup_index, subgroup_slot). Particle index `i` always
 * available. Every active belongs to exactly ONE subgroup (legacy invariant
 * — actives partitioned by `j_type_bitmask`). If a future port needs a
 * union-slot key across subgroups, that's an explicit `union_slot_indices`
 * map at port time, NOT inferred.
 *
 * Field semantics (lifetime constraints):
 *   - h_search_current is a VALUE SNAPSHOT of radii_uvm[subgroup_index]
 *     [subgroup_slot] read at hook-entry time. NOT a reference, NOT mutable.
 *     Radius mutation happens ONLY through the IterResult the hook returns;
 *     the runner stages new_h_search into radii_uvm before the next iter.
 *   - scratch is a MUTABLE REFERENCE — the hook is the place it's mutated.
 *   - scalars is a CONST REFERENCE to the driver-owned CallScalars (whole-
 *     call lifetime; never a stack temporary).
 *   - args is the driver's same-call neighbor_loop_args_iterative reference.
 * ========================================================================== */
template <typename Spec>
struct AfterIterContext {
    const neighbor_loop_args_iterative& args;
    int    subgroup_index;                           /* 0 .. args.num_subgroups - 1 */
    int    subgroup_slot;                            /* 0 .. subgroups[subgroup_index].num_active_local - 1 */
    int    i;                                        /* particle_data index */
    int    iter_index;                               /* shared across all subgroups in this iter */
    double h_search_current;                         /* value snapshot; mutate via IterResult only */
    const typename Spec::CallScalars& scalars;       /* driver-owned; whole-call lifetime */
    typename Spec::IterScratch&       scratch;       /* mutable; PERSISTS across iters (allocated/zeroed once at iter 0) */
};

/* ============================================================================
 * SFINAE detectors for iterative Spec hooks.
 *
 * Each detector returns a `_v` constexpr bool used by the runner via
 * `if constexpr (...) { ... }` to fire optional Spec hooks only when declared,
 * and by `static_assert` at run_neighbor_loop_iterative entry to require the
 * mandatory ones for Iterative Specs.
 *
 * `nlr_spec_has_after_iter_v<Spec>`:
 *   REQUIRED for IterControl = Iterative. Detector lets the runner emit a
 *   clean diagnostic naming the missing member (instead of the deeper
 *   "no matching call" error from inside the driver body) when an Iterative
 *   Spec forgets to declare `static IterResult after_iter(...)`.
 *
 * `nlr_spec_has_after_iter_global_v<Spec>`:
 *   OPTIONAL hook. Runner uses if constexpr; default
 *   absent = no-op.
 *
 * `nlr_spec_has_on_max_iter_exceeded_v<Spec>`:
 *   OPTIONAL hook. Runner default = endrun(1155) +
 *   diagnostic; Spec override only when port's legacy policy differs and
 *   Phil approved.
 * ========================================================================== */

template <typename Spec, typename = void>
struct nlr_spec_has_after_iter : std::false_type {};
template <typename Spec>
struct nlr_spec_has_after_iter<
    Spec,
    std::void_t<decltype(Spec::after_iter(
        std::declval<const AfterIterContext<Spec>&>(),
        std::declval<const typename Spec::AccumData&>()))>>
    : std::true_type {};
template <typename Spec>
constexpr bool nlr_spec_has_after_iter_v = nlr_spec_has_after_iter<Spec>::value;

template <typename Spec, typename = void>
struct nlr_spec_has_after_iter_global : std::false_type {};
template <typename Spec>
struct nlr_spec_has_after_iter_global<
    Spec,
    std::void_t<decltype(Spec::after_iter_global(
        std::declval<const neighbor_loop_args&>(),
        std::declval<const NlrIterDriver<Spec>&>()))>>
    : std::true_type {};
template <typename Spec>
constexpr bool nlr_spec_has_after_iter_global_v = nlr_spec_has_after_iter_global<Spec>::value;

/* `nlr_spec_has_apply_active_writeback_iterative_v<Spec>`:
 *   OPTIONAL hook (oracle-contract fix). When
 *   declared, the runner calls this hook INSTEAD OF apply_active_writeback
 *   in the iterative production-only post-iter-loop writeback (around
 *   line 4121 of neighbor_loop_runner.cc), passing the converged radius
 *   and per-active IterScratch alongside the converged AccumData.
 *
 *   Rationale: density needs the converged radius for its post-runner
 *   finalize, but cannot have after_iter write it to Aux (oracle runs
 *   after_iter twice — production + oracle — sharing args.aux, so the
 *   oracle pass would overwrite the production radius). The runner's
 *   final writeback loop IS production-only (oracle has separate
 *   accum_oracle_uvm + radii_oracle_uvm + no writeback call); routing
 *   the radius through here is the only clean channel.
 *
 *   Specs that don't declare this hook continue to use the original
 *   apply_active_writeback path (AGS, harness, all sink Specs). */
template <typename Spec, typename = void>
struct nlr_spec_has_apply_active_writeback_iterative : std::false_type {};
template <typename Spec>
struct nlr_spec_has_apply_active_writeback_iterative<
    Spec,
    std::void_t<decltype(Spec::apply_active_writeback_iterative(
        std::declval<const neighbor_loop_args&>(),
        std::declval<int>(),
        std::declval<int>(),
        std::declval<const typename Spec::AccumData&>(),
        std::declval<double>(),
        std::declval<const typename Spec::IterScratch&>()))>>
    : std::true_type {};
template <typename Spec>
constexpr bool nlr_spec_has_apply_active_writeback_iterative_v =
    nlr_spec_has_apply_active_writeback_iterative<Spec>::value;

template <typename Spec, typename = void>
struct nlr_spec_has_on_max_iter_exceeded : std::false_type {};
template <typename Spec>
struct nlr_spec_has_on_max_iter_exceeded<
    Spec,
    std::void_t<decltype(Spec::on_max_iter_exceeded(
        std::declval<const NlrIterDriver<Spec>&>()))>>
    : std::true_type {};
template <typename Spec>
constexpr bool nlr_spec_has_on_max_iter_exceeded_v = nlr_spec_has_on_max_iter_exceeded<Spec>::value;

/* `nlr_spec_has_reset_per_iter_device_context_v<Spec>`:
 *   OPTIONAL host-side hook called once per outer iter BEFORE per-subgroup
 *   dispatch. Use case: ags_density's `per_iter_wakeup_detected` counter
 *   (zeros it in DeviceContext at start of each iter). For Specs that
 *   don't extend DeviceContext (synthetic harness without per-iter UVM
 *   counters; sink loops which are NotIterative anyway), the hook is
 *   undeclared → runner skips the call. Signature includes iter_index
 *   for diagnostic / harness-validation use. */
template <typename Spec, typename = void>
struct nlr_spec_has_reset_per_iter_device_context : std::false_type {};
template <typename Spec>
struct nlr_spec_has_reset_per_iter_device_context<
    Spec,
    std::void_t<decltype(Spec::reset_per_iter_device_context(
        std::declval<const neighbor_loop_args_iterative&>(),
        std::declval<typename Spec::DeviceContext&>(),
        std::declval<int>()))>>
    : std::true_type {};
template <typename Spec>
constexpr bool nlr_spec_has_reset_per_iter_device_context_v =
    nlr_spec_has_reset_per_iter_device_context<Spec>::value;

/* ============================================================================
 * Compile-time spec consistency checks (used at the top of run_neighbor_loop
 * and run_neighbor_loop_iterative):
 *
 *   ActiveReduceOnly       requires std::is_same_v<ScatterData, NoScatter>
 *   NeighborScatter        requires std::is_same_v<AccumData,   NoAccum>
 *   SymmetricPairScatter   requires both non-empty (NOT NoScatter / NoAccum)
 *   GhostWritebackRequired requires both non-empty
 *
 *   IterControl == NotIterative          requires Spec NOT call run_neighbor_loop_iterative
 *                                          (verified by SFINAE on after_iter availability —
 *                                           NotIterative Specs don't declare it; calling
 *                                           run_neighbor_loop_iterative<NotIterative-Spec>
 *                                           is a clean compile error).
 *   IterControl == Iterative             requires:
 *                                          - `using IterScratch = ...;` declared
 *                                            (NoIterScratch is the explicit empty form;
 *                                             no magic default)
 *                                          - `static constexpr int max_iters >= 1;`
 *                                          - `static constexpr double mode_a_csr_buffer_factor > 1.0;`
 *                                          - `static IterResult after_iter(...)` callable
 *                                          - std::is_trivially_copyable_v<IterScratch>
 *   IterControl == Iterative AND
 *     args.num_subgroups > 1              requires `using SupportsSubgroups = std::true_type;`
 *                                          (RUNTIME check inside run_neighbor_loop_iterative;
 *                                           hard-aborts loudly on mismatch)
 *
 *   ScatterData != NoScatter              requires apply_neighbor_writeback,
 *                                                  compare_scatter
 *   DeviceContext != NeighborLoopDeviceContextBase
 *                                          requires populate_device_context
 *
 *   POD/device-copy contract (enforced at the top of run_neighbor_loop<Spec>):
 *     std::is_trivially_copyable_v<Spec::CallScalars>  (lambda-captured by value)
 *     std::is_trivially_copyable_v<Spec::ActiveData>   (UVM-staged)
 *     std::is_trivially_copyable_v<Spec::NeighborData> (built per-pair on device)
 *     std::is_trivially_copyable_v<Spec::AccumData>    (UVM-staged)
 * ========================================================================== */

/* ============================================================================
 * Dispatch threshold accessors
 *
 * Production dispatch policy lives in the constexpr Spec::modeb_threshold_sum
 * and Spec::modeb_threshold_max in each NeighborLoopSpec — these are
 * code-level dispatch policy constants for the loop, decided alongside the
 * physics. They are NOT user-tunable runtime parameters.
 *
 * Threshold env vars (GIZMO_NLR_MODEB_THRESHOLD_*, GIZMO_<LOOP>_MODEB_THRESHOLD_*)
 * are TESTERS' KNOBS only — temporary overrides for runner-internal testing.
 * They are NOT part of the production interface, NOT promoted to params.txt
 * or Config.sh, and emit a rank-0 one-shot "tester only" warning when
 * actually consumed. See the TESTERS' KNOBS block in
 * mesh/neighbor_loop_runner.cc for the full env-var listing, resolution
 * precedence, and warning semantics.
 *
 * Resolution precedence (per-loop env > global env > Spec constexpr) is
 * preserved by the _for() lookups; invalid env values silently fall through
 * to the next level (no policy change).
 * ========================================================================== */

int  gizmo_nlr_default_modeb_threshold_sum(void);
int  gizmo_nlr_default_modeb_threshold_max(void);
int  gizmo_nlr_modeb_threshold_sum_for(const char *loop_name, int spec_default);
int  gizmo_nlr_modeb_threshold_max_for(const char *loop_name, int spec_default);

/* ============================================================================
 * NLR env config — unified surface for diagnostic, control, and spike vars.
 *
 * Canonical:
 *   GIZMO_NLR_DIAG=<0|1|2|3>          0=off, 1=PHASE0 timing line per call,
 *                                     2=+dispatch trace, 3=reserved (today
 *                                     equivalent to level 2; rank-0 note)
 *   GIZMO_NLR_FORCE_MODE=A|B          tester force-mode override
 *   GIZMO_<LOOP>_FORCE_MODE=A|B       per-loop override; wins over global
 *   GIZMO_NLR_FORCE_MODEB_MAX_ACTIVE=N
 *                                     forced Mode B guardrail cap (default
 *                                     100000); prevents accidental full-N
 *                                     host-walker/oracle runs
 *   GIZMO_NLR_SPIKE_ACCUM_DUMP=1      cross-validation per-active accumulator
 *                                     dump (temporary diagnostic)
 *   GIZMO_NLR_SPIKE_NB_DUMP=1         first-call Mode A neighbor-list dump
 *                                     (temporary diagnostic)
 *   GIZMO_NLR_ORACLE=1                correctness gate (separate concern)
 *   GIZMO_NLR_ORACLE_DUMP=1           field-by-field oracle mismatch dump
 *
 * Aliases (DEPRECATED; accepted for one cycle with rank-0 warning; explicit
 * removal queued for a future cleanup pass):
 *   GIZMO_PHASE0_DIAG=1            -> GIZMO_NLR_DIAG=1
 *   GIZMO_NLR_DISPATCH_TRACE=1     -> GIZMO_NLR_DIAG=2
 *   GIZMO_NLR_FORCE_MODEA=1        -> GIZMO_NLR_FORCE_MODE=A
 *   GIZMO_NLR_FORCE_MODEB=1        -> GIZMO_NLR_FORCE_MODE=B
 *   GIZMO_MODE_B_XVAL_DUMP=1       -> GIZMO_NLR_SPIKE_ACCUM_DUMP=1
 *   GIZMO_MODE_B_XVAL_NB_DUMP=1    -> GIZMO_NLR_SPIKE_NB_DUMP=1
 *
 * Conflict policy (collective; all ranks endrun):
 *   GIZMO_NLR_FORCE_MODE set + any old _FORCE_MODE{A,B} set     -> endrun
 *   _FORCE_MODEA=1 AND _FORCE_MODEB=1                           -> endrun
 *   GIZMO_NLR_FORCE_MODE invalid value (not A or B)             -> endrun
 *   GIZMO_NLR_DIAG invalid (non-integer, <0, or >3)             -> endrun
 *
 * Diagnostic alias precedence (NEW WINS, old-set produces ignore note):
 *   GIZMO_NLR_DIAG set + diagnostic alias set -> new wins, alias ignored
 *
 * Warnings: rank-0 only, cached one-shot per env-var key.
 * ========================================================================== */

/* NlrForceMode is defined above struct neighbor_loop_args because its
 * in-class default initializer for the dispatch_override field needs
 * the enum to be fully visible. See declaration above near line 1030. */

int          gizmo_nlr_diag_level(void);              /* 0..3 (3 == 2 today) */
NlrForceMode gizmo_nlr_force_mode(void);              /* None / A / B */
NlrForceMode gizmo_nlr_force_mode_for(const char *loop_name); /* per-loop override, else global */
bool         gizmo_nlr_spike_accum_dump_enabled(void);
bool         gizmo_nlr_spike_nb_dump_enabled(void);

/* Convenience adapters — preserved for existing callers, all delegating
 * to the unified API above. */
bool gizmo_nlr_force_mode_b_global(void);
bool gizmo_nlr_force_mode_a_global(void);
bool gizmo_nlr_dispatch_trace_enabled(void);

bool gizmo_nlr_oracle_enabled_global(void);
bool gizmo_nlr_oracle_enabled_for(const char *loop_name);

/* ============================================================================
 * NlrQueryEnvelope — runner-owned transport wrapper for cross-rank queries.
 *
 * Wraps Spec::ActiveData with origin metadata so reply merge does NOT depend
 * on transport ordering (symmetric query/reply envelopes). Invariant:
 * carrying origin_slot prevents a nasty future silent bug even
 * though current transport (mode_b_p2p_transport) preserves order — the
 * runner's contract is robust to any future transport substitution.
 *
 * origin_rank kept for diagnostics + sanity asserts (reply receiver can
 * verify envelope.origin_rank == ThisTask).
 *
 * Trivially-copyable iff Spec::ActiveData is — already enforced at the top
 * of run_neighbor_loop<Spec>.
 * ========================================================================== */
template <typename ActiveData>
struct NlrQueryEnvelope {
    int        origin_slot;   /* active_slot on the origin rank (== aa) */
    int        origin_rank;   /* ThisTask of origin rank; for diagnostics */
    /* Targeted-export start-node list (legacy DataNodeList[...].NodeList).
     * n_nodes valid entries in [1, NODELISTLENGTH]; the receiver resumes its walk
     * from these exported nodes and stops at the top-level boundary instead of
     * full-walking from root. A (query,peer) needing more than NODELISTLENGTH
     * nodes is split across multiple envelopes (disjoint subtrees; the slot-keyed
     * reply merge sums them). n_nodes==0 = broadcast query (full local walk on
     * the receiver: uncovered-policy loops, or an oracle probe). */
    int        n_nodes;
    /* Oracle under-route probe (validation scaffolding, oracle mode only): 1 =
     * the sender's targeted routing did NOT select this peer for this query and
     * expects ZERO matches here; the receiver full-walks it and raises a loud
     * SENDER-UNDER-ROUTE alarm if it finds any. 0 on all production envelopes. */
    int        oracle_untargeted_probe;
    int        NodeList[NODELISTLENGTH];
    ActiveData active;
};

/* Symmetric reply envelope: peer copies {origin_slot, origin_rank} back from
 * the received query envelope into its reply, so the active rank can merge
 * by slot WITHOUT relying on transport ordering. Receive-side asserts
 * origin_rank == ThisTask as a sanity check. Trivially-copyable iff
 * AccumData is — already enforced at the top of run_neighbor_loop<Spec>. */
template <typename AccumData>
struct NlrReplyEnvelope {
    int       origin_slot;
    int       origin_rank;
    AccumData accum;
};

/* Temporary iterative-oracle transport payload.
 * This exists only to keep the port-validation oracle from doing two reply
 * exchanges with identical tags. The oracle path is validation scaffolding,
 * not permanent runner architecture; remove this once the oracle path is removed. */
template <typename AccumData>
struct NlrDualReplyEnvelope {
    int       origin_slot;
    int       origin_rank;
    AccumData accum_prod;
    AccumData accum_oracle;
};

/* ============================================================================
 * Diagnostic views and SFINAE-detected optional Spec hooks
 *
 * Optional, never part of the physics API. Specs opt in by defining the
 * static method; the runner SFINAE-detects and only calls when present.
 *
 * View structs carry full identity fields (rank, origin_rank, origin_slot,
 * active_slot, path, call_id) even where some are redundant for self/local
 * today. Locking the shape now avoids a future change when peer-side
 * dumps land.
 *
 * `path` strings preserved byte-identical to legacy:
 *   "gpu_ngl" — Mode A GPU NGL path
 *   "mode_b"  — Mode B local OR Mode B remote (active rank's own actives)
 *
 * Accum dump fires for both Mode A and runner-driven Mode B paths. NB dump
 * fires Mode A only; the Mode B local NB + UVM-readback diagnostic and the
 * Mode B remote peer-side NB dump are not implemented.
 *
 * Cached env-gate accessors below are first-use cached; mid-run env
 * changes do NOT take effect, matching this code's cached-env convention.
 * ========================================================================== */

template <typename Spec>
struct ActiveDumpView {
    int rank;
    int origin_rank;
    int origin_slot;
    int active_slot;
    const char *path;
    int call_id;
    const struct neighbor_loop_args *args;
    const typename Spec::ActiveData *active;     /* may be nullptr if spec/path doesn't supply */
    const typename Spec::AccumData  *accum;
};

template <typename Spec>
struct NeighborListDumpView {
    int rank;
    int origin_rank;
    int origin_slot;
    int active_slot;
    const char *path;
    int call_id;
    const struct neighbor_loop_args *args;
    const typename Spec::ActiveData *active;
    const int *candidate_ids;                    /* host-resident; runner-owned, valid only during call */
    int n_candidates;
};

/* Detection traits. SFINAE on the static call expression. */
template <typename Spec, typename = void>
struct spec_has_dump_active : std::false_type {};
template <typename Spec>
struct spec_has_dump_active<Spec, std::void_t<decltype(
    Spec::diagnostic_dump_active(std::declval<const ActiveDumpView<Spec>&>()))>>
    : std::true_type {};

template <typename Spec, typename = void>
struct spec_has_dump_neighbor_list : std::false_type {};
template <typename Spec>
struct spec_has_dump_neighbor_list<Spec, std::void_t<decltype(
    Spec::diagnostic_dump_neighbor_list(std::declval<const NeighborListDumpView<Spec>&>()))>>
    : std::true_type {};

/* Cached env-gate adapters. Read-once-per-process; mid-run env changes
 * do not take effect. All delegate to the unified gizmo_nlr_diag_level /
 * gizmo_nlr_spike_*_enabled accessors documented at the top of this header.
 *
 * Mapping:
 *   gizmo_nlr_xval_dump_enabled()    -> gizmo_nlr_spike_accum_dump_enabled()
 *   gizmo_nlr_xval_nb_dump_enabled() -> gizmo_nlr_spike_nb_dump_enabled()
 *   gizmo_nlr_phase0_diag_enabled()  -> (gizmo_nlr_diag_level() >= 1)
 *
 * The xval_* names are kept on existing call sites for now and will be
 * renamed alongside the env-name retire in a follow-up cleanup pass. */
bool gizmo_nlr_xval_dump_enabled(void);
bool gizmo_nlr_xval_nb_dump_enabled(void);
bool gizmo_nlr_phase0_diag_enabled(void);

/* ============================================================================
 * RunnerStageTimer — timing accumulator populated when GIZMO_PHASE0_DIAG=1
 *
 * Each path function (run_mode_a / run_mode_b_local / run_mode_b_remote)
 * receives a pointer; nullptr means "phase0 off, skip all timing." Path
 * fills only the fields that apply to its dispatch:
 *
 *   path=gpu_ngl         dt_collect (NGL build) + dt_walk_self (pair kernel)
 *                          + dt_writeback. Others = 0.
 *   path=mode_b_local    dt_collect + dt_drift + dt_walk_self + dt_writeback.
 *                          Others = 0.
 *   path=mode_b_remote   all 8 fields. dt_collect = self+peer pre-drift
 *                          collection; dt_walk_self = self_tree evaluate;
 *                          dt_walk_peer = peer_tree evaluate; dt_exchange_q,
 *                          dt_exchange_r = peer-to-peer comm; dt_reduce = reply merge;
 *                          dt_writeback = writeback loop; dt_drift = lazy
 *                          drift on union.
 *
 * dt_total measured wall-clock from runner entry to runner exit. Now
 * includes the runner-owned prep step (dt_prep_import); the caller has
 * been simplified and no longer runs prep before the runner.
 *
 * dt_prep_import is the wall around gizmo_request_filtered_ghost_import_fresh
 * on Mode A paths; 0 on Mode B paths (genuine 0 — the API isn't called).
 * ========================================================================== */

struct RunnerStageTimer {
    double dt_prep_import;   /* Mode A only; 0 on Mode B paths. */
    double dt_collect;
    double dt_drift;
    double dt_walk_self;
    double dt_walk_peer;
    double dt_exchange_q;
    double dt_exchange_r;
    double dt_reduce;
    double dt_writeback;
    double dt_total;
};

/* ============================================================================
 * nlr_build_active_list inline definition.
 *
 * Header-defined because the predicate is templated. Instantiated in the
 * caller's TU, which already pulls in declarations/allvars.h (NumPart,
 * ActiveParticleList) and core/proto.h (mymalloc) via the loop's spec
 * header.
 * ========================================================================== */
template <typename ActivePredicate>
bool nlr_build_active_list(ActivePredicate is_active,
                           int **list_out,
                           int  *num_local_active_out,
                           int  *num_global_active_out,
                           const char *malloc_label)
{
    int num_local = 0;
    for(int i : ActiveParticleList) { if(is_active(i)) num_local++; }
    int num_global = 0;
    MPI_Allreduce(&num_local, &num_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    *num_local_active_out  = num_local;
    *num_global_active_out = num_global;
    if(num_global == 0) { *list_out = nullptr; return false; }

    int alloc_n = (num_local > 0) ? num_local : 1;
    *list_out = (int *) mymalloc(malloc_label, alloc_n * sizeof(int));
    int slot = 0;
    for(int i : ActiveParticleList) {
        if(is_active(i)) { (*list_out)[slot++] = i; }
    }
    return true;
}

#endif /* GIZMO_NEIGHBOR_LOOP_RUNNER_H */
