/* Mode B local neighbor walker — host-side range-walk + brute-force oracle.
 *
 * See header for design constraints. Both paths return LOCAL real P[]
 * indices in [0, ghost_get_num_local()) intersecting (pos, h_q).
 *
 * Status: WORKING for ONEWAY mode (the density-iter target). SYMMETRIC
 * tree walk currently always opens (no per-node max-h tracking yet);
 * brute walk handles SYMMETRIC correctly. Oracle catches any divergence.
 *
 * Periodic boundaries: handled via NEAREST_XYZ on the displacement.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../gravity/forcetree.h"
#include "ghost_writeback.h"      /* ghost_get_num_local */
#include "gpu_neighbor_list.h"    /* gizmo_mark_kernel_radius_dirty_indices */
#include "mode_b_local_walker.h"

/* Public helper: per-j symmetric radius under policy. See header. */
double mode_b_neighbor_symmetric_radius(int j, mode_b_radius_policy_t policy)
{
    /* Gas: KernelRadius represents finite-volume cell extent. Always
     * physically meaningful for gas. */
    if(P[j].Type == 0) {
        if(policy & MODE_B_RADIUS_GAS_KERNEL) return (double)P[j].KernelRadius;
        return 0.0;
    }
    /* Non-gas: P[j].KernelRadius is NOT physical extent (it's stale-IC or
     * the density-search radius for THIS particle's gas neighbors, not
     * its own physical kernel). For SYMMETRIC pair searches the correct
     * source is AGS_KernelRadius, used only by AGS/SIDM/grain physics. */
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
    if(policy & MODE_B_RADIUS_AGS_FOR_NONGAS) {
        return (double)P[j].AGS_KernelRadius;
    }
#endif
    /* Default for non-gas: 0 → SYMMETRIC degenerates to ONEWAY for this j. */
    return 0.0;
}

/* Predicate: does P[j] satisfy the query (pos, h_q) under search_mode? */
/* j_radius_scale: SYMMETRIC-mode multiplier on the j-side kernel radius
 * (1.0 = legacy). TURB_DIFF_DYNAMIC wide-filter loops pass
 * All.TurbDynamicDiffFac so the Mode B reach matches the Mode A scaled-
 * symmetric NGL. See OPEN_3d_difffilter_design.md §3. */
static inline int particle_passes(int j,
                                  const double pos[3],
                                  double h_q,
                                  unsigned int type_mask,
                                  int search_mode,
                                  mode_b_radius_policy_t radius_policy,
                                  double j_radius_scale)
{
    if(!(type_mask & (1u << P[j].Type))) return 0;
    if(P[j].Mass <= 0) return 0;
    double dx = (double)P[j].Pos[0] - pos[0];
    double dy = (double)P[j].Pos[1] - pos[1];
    double dz = (double)P[j].Pos[2] - pos[2];
    NEAREST_XYZ(dx, dy, dz, 1);
    double r2 = dx*dx + dy*dy + dz*dz;
    double cutoff = h_q;
    if(search_mode == MODE_B_SEARCH_SYMMETRIC) {
        /* Type-aware: only use h_j when it's physically meaningful for
         * j's type under the active radius_policy. For non-gas types
         * without AGS opt-in, h_j contribution is 0 and the test
         * collapses to ONEWAY r < h_q for that j. (Phil + codex 2026-05-07.)
         */
        double hj = mode_b_neighbor_symmetric_radius(j, radius_policy) * j_radius_scale;
        if(hj > cutoff) cutoff = hj;
    }
    return r2 < cutoff * cutoff;
}

void mode_b_local_brute_walk(const double pos[3],
                             double h_q,
                             unsigned int type_mask,
                             int search_mode,
                             mode_b_radius_policy_t radius_policy,
                             std::vector<int>& out,
                             double j_radius_scale)
{
    const int num_local = ghost_get_num_local();
    for(int j = 0; j < num_local; j++) {
        if(!particle_passes(j, pos, h_q, type_mask, search_mode, radius_policy, j_radius_scale)) continue;
        out.push_back(j);
    }
}

/* Sphere-vs-AABB pruning test. Returns 1 if the sphere of radius R
 * centered at pos overlaps the AABB defined by node center + half-len.
 * Uses NEAREST_XYZ for periodic minimum image. */
static inline int sphere_aabb_overlap(const double pos[3],
                                      const struct NODE *nop,
                                      double R)
{
    double dx = (double)nop->center[0] - pos[0];
    double dy = (double)nop->center[1] - pos[1];
    double dz = (double)nop->center[2] - pos[2];
    NEAREST_XYZ(dx, dy, dz, 1);
    /* Conservative bound: cube half-diagonal sqrt(3)/2 * len */
    const double SQRT3_OVER_2 = 0.86602540378443864676;
    double r_max = R + (double)nop->len * SQRT3_OVER_2;
    return (dx*dx + dy*dy + dz*dz) < r_max * r_max;
}

/* Mode B SYMMETRIC effective radius for an internal node, given a type-mask
 * (which types the caller is searching for) and a radius policy (which
 * particle radius source is meaningful per type — see mode_b_neighbor_
 * symmetric_radius). Returns max over types-in-mask of the relevant
 * Extnodes[no].hmax_per_type[t] band, multiplied by a slack factor.
 * Returns 0 if no type contributes (collapses SYMMETRIC pruning to ONEWAY).
 *
 * Slack rationale (matches mesh/gpu_neighbor_list.cc:SIDX_H_SLACK = 0.5):
 * between force_update_hmax() refreshes, per-particle KernelRadius can
 * grow under drift up to factor exp(divv_fac_max/NUMDIMS) ≈ 1.105 (gas)
 * or exp(4/3) ≈ 3.79 (AGS-active). Node hmax decays under a different
 * (looser) clamp and can fall BELOW the live max-particle radius. The
 * 50% inflation absorbs this asymmetry conservatively — over-search is
 * safe (extra candidates → physics kernel filters), under-search is a
 * correctness bug. Same convention as GPU NGL's BVH tile-overlap test. */
static inline double mode_b_node_symmetric_radius(int no,
                                                  unsigned int type_mask,
                                                  mode_b_radius_policy_t policy,
                                                  double j_radius_scale)
{
    static constexpr double H_SLACK = 0.5;  /* matches SIDX_H_SLACK */
    double rmax = 0.0;
    /* Type 0 (gas): always contributes if requested AND policy includes
     * the gas-kernel source (default for sink_env1 / density-style callers). */
    if((type_mask & (1u << 0)) && (policy & MODE_B_RADIUS_GAS_KERNEL)) {
        double v = (double)Extnodes[no].hmax_per_type[0];
        if(v > rmax) rmax = v;
    }
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
    /* Non-gas types (1..5): contribute only if the caller requests them in
     * the type_mask AND the policy permits AGS as the radius source AND
     * the type is in ADAPTIVE_GRAVSOFT_FORALL bitmask (so the band was
     * actually populated at tree-build time). */
    if(policy & MODE_B_RADIUS_AGS_FOR_NONGAS) {
        for(int t = 1; t < 6; t++) {
            if(!(type_mask & (1u << t))) continue;
            if(!((1 << t) & ADAPTIVE_GRAVSOFT_FORALL)) continue;
            double v = (double)Extnodes[no].hmax_per_type[t];
            if(v > rmax) rmax = v;
        }
    }
#endif
    return rmax * (1.0 + H_SLACK) * j_radius_scale;
}

void mode_b_local_neighbor_walk(const double pos[3],
                                double h_q,
                                unsigned int type_mask,
                                int search_mode,
                                mode_b_radius_policy_t radius_policy,
                                std::vector<int>& out,
                                double j_radius_scale)
{
    if(All.MaxPart <= 0 || Nodes == NULL || Nextnode == NULL) return;
    const int num_local = ghost_get_num_local();
    const int max_part  = All.MaxPart;
    const int pseudo_start = max_part + MaxNodes + MaxForeignNodes;

    /* Stage 4 v1: SYMMETRIC pruning uses per-type hmax bands maintained in
     * extNODE (see allvars.h). For each internal node we compute
     *   R_eff = max(h_query, max_{t in type_mask, populated under policy}
     *                          Extnodes[no].hmax_per_type[t])
     * then sphere-vs-AABB overlap with R_eff. Matches old GIZMO's
     * scalar-hmax pattern at the per-type level (legacy_hmax_archaeology.md).
     * ONEWAY ignores hmax entirely and uses h_query alone.
     *
     * NOTE (Stage 3 deferred): cross-rank DomainMoment exchange of per-type
     * bands is not yet wired. At single-rank, host-only, all values are
     * locally consistent. Multi-rank Mode B SYMMETRIC walks reading
     * pseudo-particle nodes with stale hmax_per_type[] could under-prune
     * (return extra candidates → still correct, but slower) or over-prune
     * (return missing candidates → INCORRECT). The walker currently never
     * descends into foreign-pseudo nodes (they're skipped at the bottom of
     * the loop), so this is only a perf concern at np>1 today. */
    const int oneway = (search_mode == MODE_B_SEARCH_ONEWAY);

    int no = max_part; /* root node */

    while(no >= 0) {
        if(no < max_part) {
            /* Particle leaf. Only return if it's a domain-owned local
             * particle (not a ghost import). */
            if(no < num_local && particle_passes(no, pos, h_q, type_mask, search_mode, radius_policy, j_radius_scale)) {
                out.push_back(no);
            }
            no = Nextnode[no];
        } else if(no < pseudo_start) {
            /* Internal node — drift if stale, then prune. */
            struct NODE *nop = &Nodes[no];
            if(nop->Ti_current != All.Ti_Current) {
#ifdef _OPENMP
#pragma omp critical(_modebdrift_)
#endif
                {
                    if(nop->Ti_current != All.Ti_Current) {
                        force_drift_node(no, All.Ti_Current);
                    }
                }
            }
            double R_eff;
            if(oneway) {
                R_eff = h_q;
            } else {
                /* SYMMETRIC: R_eff = max(h_q, per-type hmax under policy/mask) */
                double node_h = mode_b_node_symmetric_radius(no, type_mask, radius_policy, j_radius_scale);
                R_eff = (node_h > h_q) ? node_h : h_q;
            }
            int do_open = sphere_aabb_overlap(pos, nop, R_eff);
            no = do_open ? nop->u.d.nextnode : nop->u.d.sibling;
        } else {
            /* Pseudo-particle node (LET / cross-rank). These are not
             * addressable NODE slots in the Phase-9 layout. Match the
             * legacy force walkers: skip through the corresponding
             * Nextnode[] entry, whose pseudo index is shifted by the
             * local-node and foreign-node reservation. */
            no = Nextnode[no - MaxNodes - MaxForeignNodes];
        }
    }
}

/* Lazy-drift Mode B candidates to current Ti before the pair kernel reads
 * P[j] / CellP[j]. Mirrors gpu_ngb_list_build:1542-1580 contract. */
void mode_b_lazy_drift_candidates(const int *indices, int n)
{
    if(!indices || n <= 0) return;
    const int num_local = ghost_get_num_local();
    integertime time1 = All.Ti_Current;
    for(int k = 0; k < n; k++) {
        int j = indices[k];
        if(j >= 0 && j < num_local) {
            /* drift_particle's early-return on time1==time0 makes repeat
             * calls a fast no-op (e.g. j touched by an earlier query in
             * this same evaluator call). */
            drift_particle(j, time1);
        }
    }
    /* drift_particle mutates Pos and KernelRadius (the *= exp(divv_fac/N)
     * factor in predict.cc:160,229). The next gpu_ngb_list_build call (for
     * non-Mode-B callers) needs to refresh compact_h from these freshly
     * drifted KernelRadius values. */
    gizmo_mark_kernel_radius_dirty_indices(indices, n);
}

/* Per-call same-rank tree-vs-brute oracle is owned by the runner via
 * GIZMO_NLR_ORACLE in mesh/neighbor_loop_runner.cc::run_mode_b_*_with_oracle.
 * The legacy public wrapper mode_b_local_walk_with_oracle() and its
 * GIZMO_MODE_B_ORACLE env helper were retired in this commit; both had
 * zero callers post-3c.3. */
