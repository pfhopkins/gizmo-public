/* Mode B local neighbor walker — host-side range-walk over the existing
 * gravity tree (Nodes[]/Nextnode[]/sibling). Designed to answer remote
 * Mode B queries from a peer rank's local particle set without touching
 * GPU SIDX state.
 *
 * Two implementations behind one API:
 *   - mode_b_local_neighbor_walk()  — TREE-WALK path. Fast.
 *   - mode_b_local_brute_walk()     — BRUTE-FORCE path. The oracle.
 *
 * Both return the same set of LOCAL real P[] indices intersecting the
 * spherical query (pos, h_q) with type/mode filter. Sorted ascending.
 *
 * Oracle: cross-rank Mode B oracle is owned by the runner via
 * GIZMO_NLR_ORACLE (mesh/neighbor_loop_runner.cc). The runner calls
 * both walks for each query and diffs results.
 *
 * Design constraints (codex 2026-05-07):
 *   - Walk uses Nodes[]/Nextnode[] for pruning. Node bounds come from
 *     force_drift_node(); they are conservative for current Ti_Current.
 *   - Returned candidates are LOCAL real P[] indices in [0, num_local)
 *     where num_local = ghost_get_num_local(). Never LET pseudo nodes,
 *     never ghost imports.
 *   - Particle drift is the CALLER's responsibility. The walker does
 *     NOT drift P[j] before evaluating the predicate.
 *   - Not thread-safe with concurrent particle drift / tree mutation.
 *   - Symmetric mode without per-node max-h tracking degenerates to
 *     "always open"; ONEWAY is the fast path.
 */

#ifndef MODE_B_LOCAL_WALKER_H
#define MODE_B_LOCAL_WALKER_H

#include <vector>

/* search_mode constants — match NGB_SEARCH_ONEWAY/SYMMETRIC in
 * mesh/ghost_exchange_spec.h to keep the existing NGL contract. */
#ifndef MODE_B_SEARCH_ONEWAY
#define MODE_B_SEARCH_ONEWAY    0
#define MODE_B_SEARCH_SYMMETRIC 1
#endif

/* Type-aware "symmetric radius" policy (Phil + codex 2026-05-07).
 *
 * For SYMMETRIC search, the test is r < max(h_query, h_j). But h_j is only
 * physically meaningful for some particle types — and only via specific
 * fields per type. Old GIZMO's tree Hmax was GAS-ONLY because gas cells
 * have finite volume that may overlap a sink even if the cell's mesh-point
 * is outside the sink's KernelRadius. Non-gas types (DM, stars, sinks) are
 * point-like in the relevant physics; their `KernelRadius` field is either
 * stale-from-IC or set by density() purely for THEIR gas-neighbor search,
 * not as a physical extent.
 *
 * Exception: physics with type-mutual interactions (AGSForce, SIDM,
 * grain collisions, dm_fuzzy, cbe) uses `AGS_KernelRadius` as the
 * symmetric radius for non-gas particles.
 *
 * Most callers want the default (gas via KernelRadius, non-gas → 0 →
 * SYMMETRIC degenerates to ONEWAY for those types). AGS-physics callers
 * opt in via the bitflag.
 */
typedef unsigned int mode_b_radius_policy_t;
#define MODE_B_RADIUS_GAS_KERNEL          (1u << 0)  /* Type=0 SYM uses P[j].KernelRadius (always recommended) */
#define MODE_B_RADIUS_AGS_FOR_NONGAS      (1u << 1)  /* non-gas SYM uses P[j].AGS_KernelRadius (AGSForce/SIDM/grain) */
#define MODE_B_RADIUS_DEFAULT             (MODE_B_RADIUS_GAS_KERNEL)

/* Returns the per-j symmetric radius to use under the given policy.
 * Returns 0 when SYMMETRIC should degenerate to ONEWAY for this j's type. */
double mode_b_neighbor_symmetric_radius(int j, mode_b_radius_policy_t policy);

/* Tree-walk path. Fast for spatially-localized queries. APPENDS local-
 * real-particle candidates (P[] indices) onto `out` via push_back; the
 * caller passes a (typically `clear()`'d) vector and the walker grows it
 * geometrically as needed. Does NOT clear `out` itself (so callers may
 * accumulate across queries if desired). Does NOT sort.
 *
 * radius_policy controls how SYMMETRIC search uses h_j for non-gas types.
 * Most callers want MODE_B_RADIUS_DEFAULT (gas KR only).
 *
 * Stage 4 contract change: previous int*+capacity signature replaced by
 * std::vector<int>& to remove the per-query num_local-sized allocation
 * the runner used to make ahead of each call (~24 MB × N_active on
 * fire_m11i — the dominant tiny-N cost post-Stage 2). Geometric growth
 * via push_back handles correctness for any-size match set without
 * imposing full-pool memory traffic on tiny-N. */
/* j_radius_scale: SYMMETRIC-mode multiplier on the j-side kernel radius
 * (1.0 = legacy). TURB_DIFF_DYNAMIC wide-filter loops pass
 * All.TurbDynamicDiffFac so the Mode B reach matches the Mode A scaled-
 * symmetric NGL. See OPEN_3d_difffilter_design.md §3. */
void mode_b_local_neighbor_walk(const double pos[3],
                                double h_q,
                                unsigned int type_mask,
                                int search_mode,
                                mode_b_radius_policy_t radius_policy,
                                std::vector<int>& out,
                                double j_radius_scale = 1.0);

/* Brute-force path. Iterates 0..num_local. Slow but obviously correct.
 * Used as the runner-owned oracle. Same append-oriented contract as the
 * tree walker: appends matches to `out` via push_back; does NOT clear
 * or sort. */
void mode_b_local_brute_walk(const double pos[3],
                             double h_q,
                             unsigned int type_mask,
                             int search_mode,
                             mode_b_radius_policy_t radius_policy,
                             std::vector<int>& out,
                             double j_radius_scale = 1.0);

/* Lazy-drift contract for Mode B (mirrors gpu_ngb_list_build:1542-1580).
 *
 * GPU NGL contract: candidate walk runs on whatever P[j].Pos state exists,
 * THEN drift_particle(j, All.Ti_Current) is called for each j in the CSR,
 * THEN the pair kernel reads drifted P[j]. The walk may see slightly stale
 * positions (h-slack absorbs that); the kernel always reads fresh.
 *
 * Mode B must honor the same contract. Walker collects candidates first;
 * call this helper to drift them; then run the pair kernel.
 *
 * Drifts each j in indices[] to current Ti via drift_particle(). Marks
 * kernel-radius dirty for the GPU SIDX tracker so the next gpu_ngb_list_build
 * call sees fresh compact_h. drift_particle's early-return on time1==time0
 * dedupes naturally.
 *
 * NOT thread-safe with concurrent drift / tree mutation. Caller serializes.
 */
void mode_b_lazy_drift_candidates(const int *indices, int n);

#endif /* MODE_B_LOCAL_WALKER_H */
