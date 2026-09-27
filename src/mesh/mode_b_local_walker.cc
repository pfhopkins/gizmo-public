/* Mode B local neighbor walker — host-side range-walk.
 *
 * See header for design constraints. Returns LOCAL real P[]
 * indices in [0, ghost_get_num_local()) intersecting (pos, h_q).
 *
 * SYMMETRIC tree walk prunes internal nodes by the per-type hmax bands
 * (Extnodes[no].hmax_per_type via mode_b_node_symmetric_radius); ONEWAY
 * prunes by h_q alone.
 *
 * The three public tree walks (mode_b_local_neighbor_walk / _walk_and_export /
 * _walk_from_start_nodes) are thin wrappers over the one shared traversal body
 * mode_b_walk_impl below — the same SSOT-body/many-wrappers shape legacy uses
 * (ngb.cc's 8 ngb_treefind_* over one codeblock).
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
#include "../gravity/force_node_drift_sync.h"  /* modeb_node_ti_current_acquire */
#include "ghost_writeback.h"      /* ghost_get_num_local */
#include "ghost_exchange_functions.h" /* gx_extended_overlap_wrap_and_test: canonical-wrap SSOT */
#include "gpu_neighbor_list.h"    /* gizmo_mark_kernel_radius_dirty_indices */
#include "mode_b_local_walker.h"

/* Public helper: per-j symmetric radius under policy. Delegates to the SSOT
 * per-particle wrapper in nlr_radius_policy.h, which owns the
 * AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE compile-flag gate. */
double mode_b_neighbor_symmetric_radius(int j, mode_b_radius_policy_t policy)
{
    return nlr_particle_symmetric_radius(P[j], policy);
}

/* Predicate: does P[j] satisfy the query (pos, h_q) under search_mode? */
/* j_reach_scale: SYMMETRIC-mode multiplier on the j-side kernel radius
 * (1.0 = legacy). TURB_DIFF_DYNAMIC wide-filter loops pass
 * All.TurbDynamicDiffFac so the Mode B reach matches the Mode A scaled-
 * symmetric NGL. It is the TOTAL j-side scale, not one named factor: a caller
 * whose acceptance test applies further multipliers (ghost exchange folds its
 * safety factor in here) must include them, or the walk searches a smaller
 * neighbourhood than the caller then accepts from. */
static inline int particle_passes(int j,
                                  const double pos[3],
                                  double h_q,
                                  unsigned int type_mask,
                                  int search_mode,
                                  mode_b_radius_policy_t radius_policy,
                                  double j_reach_scale)
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
         * collapses to ONEWAY r < h_q for that j.
         */
        double hj = mode_b_neighbor_symmetric_radius(j, radius_policy) * j_reach_scale;
        if(hj > cutoff) cutoff = hj;
    }
    return r2 < cutoff * cutoff;
}

/* Sphere-vs-AABB pruning test. Returns 1 if the sphere of radius R
 * centered at pos overlaps the AABB defined by node center + half-len.
 * Uses NEAREST_XYZ for periodic minimum image. */
static inline int sphere_aabb_overlap(const double pos[3],
                                      const struct NODE *nop,
                                      double R)
{
    /* Node open test.  Both the box wrap (including the shearing-box straddle
     * case) and the legacy acceptance geometry live in the shared predicate, so
     * this walk and the BVH/tile walks cannot drift apart.  A cube's
     * circumradius is 0.866*len, so the shared per-axis-half-width form
     * reproduces the legacy node bound exactly. */
    const double hw = 0.5 * (double)nop->len;
    return gx_extended_overlap_wrap_and_test((double)nop->center[0] - pos[0],
                                             (double)nop->center[1] - pos[1],
                                             (double)nop->center[2] - pos[2],
                                             hw, hw, hw, R);
}

/* Mode B SYMMETRIC effective radius for an internal node. Returns
 *   max over types-in-mask of Extnodes[no].hmax_per_type[t]
 * inflated by a slack factor. Per the invariant in allvars.h, each band is
 * already a conservative upper bound across every leaf-policy-selectable
 * radius source for that type — so node-prune does not need to know the
 * caller's radius_policy. Leaf-level predicate (mode_b_neighbor_symmetric_radius)
 * applies the exact policy. Over-opening here is safe (extra candidates
 * filter at the leaf); under-opening would be a correctness bug.
 *
 * Returns 0 if no requested type has a populated band (degenerates to ONEWAY
 * pruning, which is still correct given the leaf-level filter).
 *
 * Node-open slack: none is carried. A drift grows a particle's radius by at
 * most exp(KERNEL_RADIUS_DRIFT_MAX_LOG_CHANGE/NUMDIMS), and force_drift_node
 * grows the node's hmax and bands by the same capped rule applied to divVmax.
 * divVmax is gathered from gas members only, so a band holding adaptive-
 * softening radii of non-gas members is not grown by its own members'
 * divergence between force_update_hmax refreshes. Over-search is safe (extra
 * candidates filter at the leaf); under-search is a correctness bug. */
static constexpr double MODE_B_NODE_H_SLACK = 0.0;  /* no node-open drift slack; legacy has none (relies on force_update_hmax cadence + 0.866*len node term) */

static inline double mode_b_node_symmetric_radius(int no,
                                                  unsigned int type_mask,
                                                  double j_reach_scale)
{
    double rmax = 0.0;
    for(int t = 0; t < 6; t++) {
        if(!(type_mask & (1u << t))) continue;
        double v = (double)Extnodes[no].hmax_per_type[t];
        if(v > rmax) rmax = v;
    }
    return rmax * (1.0 + MODE_B_NODE_H_SLACK) * j_reach_scale;
}

/* SYMMETRIC node-open reach: a SINGLE per-type band for BOTH candidate traversal
 * and targeted export (mode_b_node_symmetric_radius = max_{t in mask}
 * Extnodes[no].hmax_per_type[t]).  ONEWAY uses h_q (no band).
 *
 * Post-substrate (sink-radius cap + DomainNODE per-type exchange +
 * force_update_hmax post-density per-type exchange), this band is:
 *   - a conservative UPPER BOUND over every leaf-policy-selectable source per type
 *     (allvars.h invariant; the seed is capped at MaxKernelRadius and every kernel/
 *     AGS radius is <= MaxKernelRadius by construction, ForceSoftening seeded
 *     uncapped), so it DOMINATES every radius_policy -- the exact policy filter
 *     lives at the leaf (mode_b_neighbor_symmetric_radius), never the node band;
 *   - cross-rank-correct AND fresh on exactly the nodes the export walk descends:
 *     remote topleaves (DomainNODE pack/apply + force_update_hmax) and their
 *     INTERNAL_TOPLEVEL ancestors (gpu_topnode_moment_resum + post-density
 *     up-propagation) -- so it bounds every loop's j-side reach on remote peers
 *     and the sender never under-routes.
 * Box nesting (box_A superset box_child) + monotone band (band_A = max children)
 * => opening a remote topleaf T opens every ancestor of T, so no export is missed
 * (same strict-refinement guarantee as the LET cover tree).  Hence traversal
 * reach == export reach: ONE open predicate, no scalar/per-type dual path.  The
 * cand_out / export_out sinks only choose WHERE a reached node is recorded (local
 * candidate vs remote-topleaf export), never the reach.
 *
 * (Historically the export decision used the gas-biased cross-rank SCALAR hmax and
 * a fused walk re-tested exports against it; the non-gas SYMMETRIC loops could not
 * be bounded by scalar hmax and stayed on broadcast.  The per-type band above now
 * covers them.) */

/* Build the topleaf reverse map from the DomainNodeIndex SSOT (never from a
 * slot-layout assumption). Sized to the max observed offset; entries outside
 * any topleaf stay -1. O(NTopleaves), rebuilt per export call (topnode indices
 * are stable between tree builds; per-call rebuild avoids any staleness). */
void ModeBTopleafMap::build(void)
{
    const int tree_base = All.TreeNodeIndexBase;
    int max_off = -1;
    for(int i = 0; i < NTopleaves; i++) {
        const int off = DomainNodeIndex[i] - tree_base;
        if(off < 0) { topnode_map_size = 0; return; }   /* malformed map: disable topleaf detection (walk still exports via pseudo branch) */
        if(off > max_off) max_off = off;
    }
    topnode_map_size = max_off + 1;
    leaf_of_topnode.assign(topnode_map_size, -1);
    for(int i = 0; i < NTopleaves; i++) {
        leaf_of_topnode[DomainNodeIndex[i] - tree_base] = i;
    }
}

/* Shared traversal body (SSOT for all three public tree walks).
 *
 * Walks from `start_no`. Local real-particle matches are appended to `cand_out`
 * (nullptr = do not collect). At every remote pseudo-node reached, if
 * `export_out != nullptr` the owner peer + node's DomainNodeIndex is recorded
 * (legacy mode==0 targeted export, ngb_codeblock_after_condition_unthreaded.h:
 * 19-67; modern symbol form matches gravity/forcetree.cc:2173-2210). If
 * `stop_at_toplevel` (legacy mode==1 receiver), the walk returns when it
 * re-enters the top-level tree — the exported subtree is exhausted
 * (after_condition_unthreaded.h:74-81).
 *
 * SYMMETRIC internal-node pruning uses the per-type hmax bands
 * (mode_b_node_symmetric_radius); ONEWAY prunes by h_q alone. Bands are
 * rank-local and re-seeded every build/refresh; a query against another rank's
 * pool is shipped there and answered with that rank's own fresh bands, so no
 * cross-rank band exchange is needed. */
static void mode_b_walk_impl(const double pos[3],
                             double h_q,
                             unsigned int type_mask,
                             int search_mode,
                             mode_b_radius_policy_t radius_policy,
                             double j_reach_scale,
                             int start_no,
                             bool stop_at_toplevel,
                             std::vector<int>* cand_out,
                             const ModeBTopleafMap* topleaf_map,
                             ModeBExportSink* export_out)
{
    if(All.TreeNodeIndexBase <= 0 || Nodes == NULL || Nextnode == NULL) return;
    const int num_local = ghost_get_num_local();
    const int tree_base  = All.TreeNodeIndexBase;
    const int tree_slots = All.TreeParticleSlots;
    const int pseudo_start = tree_base + MaxNodes + MaxForeignNodes;
    const int oneway = (search_mode == MODE_B_SEARCH_ONEWAY);

    int no = start_no;

    while(no >= 0) {
        if(no >= tree_slots && no < tree_base) {/* An index between the particle slots and the node base belongs to neither, so the tree is
             * malformed; stop rather than read a side array or Nodes[] out of bounds. */
            endrun(90001024); no = -1; continue;}
        if(no < tree_slots) {
            /* Particle leaf. Only return domain-owned local particles (not a
             * ghost import). cand_out==nullptr on a pure export-discovery walk. */
            if(cand_out && no < num_local &&
               particle_passes(no, pos, h_q, type_mask, search_mode, radius_policy, j_reach_scale)) {
                cand_out->push_back(no);
            }
            no = Nextnode[no];
        } else if(no < pseudo_start) {
            /* Internal node. */
            struct NODE *nop = &Nodes[no];
            /* Receiver (legacy mode==1): re-entering the top-level tree means
             * the exported branch is done. */
            if(stop_at_toplevel && (nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))) return;
            /* Drift if stale, then prune. Acquire-load Ti_current so a threaded
             * walk that sees it fresh also sees the drifter's fresh geometry
             * (paired with the release store in force_drift_node). */
            if(modeb_node_ti_current_acquire(no) != All.Ti_Current) {
#ifdef _OPENMP
#pragma omp critical(_modebdrift_)
#endif
                {
                    /* Re-check inside the lock: another thread may have drifted
                     * this node between the fast-path load and here. */
                    if(modeb_node_ti_current_acquire(no) != All.Ti_Current) {
                        force_drift_node(no, All.Ti_Current);
                    }
                }
            }
            /* Single open predicate for traversal AND export: the per-type band
             * for SYMMETRIC (dominant + cross-rank-fresh -> bounds every loop's
             * remote reach), h_q for ONEWAY.  No scalar/per-type dual path. */
            double R_open;
            if(oneway) {
                R_open = h_q;
            } else {
                const double node_h = mode_b_node_symmetric_radius(no, type_mask, j_reach_scale);
                R_open = (node_h > h_q) ? node_h : h_q;
            }
            int do_open = sphere_aabb_overlap(pos, nop, R_open);
            /* SENDER export event (legacy pseudo-hit equivalence): the walk OPENED
             * a remote-owned TOPLEAF. Legacy would descend to its pseudo child
             * (export + skip); post-LET that child may be the imported foreign
             * subtree, which holds NO owned-local candidates and must NOT be
             * descended in source-export semantics. Export on the SAME predicate
             * that opened it (traversal reach == export reach) and skip to the
             * sibling REGARDLESS (no owned-local candidates below a remote topleaf). */
            if(do_open && export_out) {
                const int leaf = topleaf_map->topleaf_of(no, tree_base);
                if(leaf >= 0 && DomainTask[leaf] != ThisTask) {
                    export_out->add(DomainTask[leaf], DomainNodeIndex[leaf]);
                    no = nop->u.d.sibling;
                    continue;
                }
            }
            if(do_open) {
                const int child = nop->u.d.nextnode;
                /* Legacy foreign-subtree skip (ALL walk modes): an imported
                 * foreign subtree (post-LET rewire; nodes in [TreeNodeIndexBase+MaxNodes,
                 * pseudo_start)) holds NO owned-local P[] candidates and NO
                 * pseudo-nodes (probe-verified), so no walk mode can gain
                 * anything by descending it — legacy never descends remote
                 * subtrees. Take the sibling. Normally dead code for export-
                 * capable walks (the topleaf branch above skips first); it
                 * covers the malformed-topleaf-map fallback and the plain
                 * local-candidate walk. */
                if(child >= tree_base + MaxNodes && child < pseudo_start) {
                    no = nop->u.d.sibling;
                } else {
                    no = child;
                }
            } else {
                no = nop->u.d.sibling;
            }
        } else {
            /* Pseudo-particle node (cross-rank subtree root; reached only for
             * remote topleaves the LET did not ship/redirect — with the
             * topleaf-boundary export above this branch is normally never hit,
             * but it keeps non-LET / unshipped-leaf / malformed-topleaf-map
             * configs correct). Gate the export on the open predicate against the
             * topleaf node's own geometry (per-type reach for SYMMETRIC), then
             * skip forward exactly as the legacy
             * force walkers: the pseudo index is shifted by the local-node +
             * foreign-node reservation (Phase-9 layout). */
            if(export_out) {
                const int leaf = no - pseudo_start;
                if(leaf >= 0 && leaf < NTopleaves && DomainTask[leaf] != ThisTask) {
                    int do_export = 1;
                    if(!oneway) {   /* re-test at the topleaf's OWN box + per-type reach */
                        const int tl_node = DomainNodeIndex[leaf];
                        if(tl_node >= tree_base && tl_node < pseudo_start) {
                            const double node_h = mode_b_node_symmetric_radius(tl_node, type_mask, j_reach_scale);
                            const double R_exp_tl = (node_h > h_q) ? node_h : h_q;
                            do_export = sphere_aabb_overlap(pos, &Nodes[tl_node], R_exp_tl);
                        }
                    }
                    if(do_export) export_out->add(DomainTask[leaf], DomainNodeIndex[leaf]);
                }
            }
            no = Nextnode[tree_slots + (no - tree_base - MaxNodes - MaxForeignNodes)];
        }
    }
}

/* LOCAL-candidates walk (legacy-mode==0 without export): from root, collect
 * local matches, skip pseudo-nodes. Behavior byte-identical to the prior broadcast-only
 * walker (export_out=nullptr, stop_at_toplevel=false). */
void mode_b_local_neighbor_walk(const double pos[3],
                                double h_q,
                                unsigned int type_mask,
                                int search_mode,
                                mode_b_radius_policy_t radius_policy,
                                std::vector<int>& out,
                                double j_reach_scale)
{
    mode_b_walk_impl(pos, h_q, type_mask, search_mode, radius_policy, j_reach_scale,
                     /*start_no=*/All.TreeNodeIndexBase, /*stop_at_toplevel=*/false,
                     &out, /*topleaf_map=*/nullptr, /*export_out=*/nullptr);
}

/* SENDER walk: from root, record targeted exports at reached remote topleaves
 * / pseudo-nodes; optionally ALSO collect local candidates in the same
 * traversal (cand_out != nullptr = the FUSED legacy-mode==0 walk: one
 * traversal, two sinks, each gated by its own predicate — see the
 * ModeBWalkReach table). Export decisions always use the cross-rank scalar
 * hmax reach so targeting covers remote h_j. */
void mode_b_walk_and_export(const double pos[3],
                            double h_q,
                            unsigned int type_mask,
                            int search_mode,
                            mode_b_radius_policy_t radius_policy,
                            std::vector<int>* cand_out,
                            const ModeBTopleafMap& topleaf_map,
                            ModeBExportSink& sink,
                            double j_reach_scale)
{
    mode_b_walk_impl(pos, h_q, type_mask, search_mode, radius_policy, j_reach_scale,
                     /*start_no=*/All.TreeNodeIndexBase, /*stop_at_toplevel=*/false,
                     cand_out, &topleaf_map, &sink);
}

/* RECEIVER walk: resume from each exported start-node (open its children like
 * legacy density.cc:272), stop at the top-level boundary. Multiple entries
 * (legacy NodeList) cover disjoint exported subtrees. */
void mode_b_walk_from_start_nodes(const double pos[3],
                                  double h_q,
                                  unsigned int type_mask,
                                  int search_mode,
                                  mode_b_radius_policy_t radius_policy,
                                  const int *node_list,
                                  int n_nodes,
                                  std::vector<int>& out,
                                  double j_reach_scale)
{
    if(All.TreeNodeIndexBase <= 0 || Nodes == NULL || Nextnode == NULL) return;
    const int tree_base      = All.TreeNodeIndexBase;
    const int pseudo_start  = tree_base + MaxNodes + MaxForeignNodes;
    for(int k = 0; k < n_nodes; k++) {
        const int nl = node_list[k];
        if(nl < 0) break;   /* -1 terminator (legacy NodeList convention) */
        /* Defensive: a start-node arrives over MPI; it must be an internal
         * node index. Skip a corrupt entry rather than dereference OOB. */
        if(nl < tree_base || nl >= pseudo_start) continue;
        const int start = Nodes[nl].u.d.nextnode;   /* open the exported node */
        mode_b_walk_impl(pos, h_q, type_mask, search_mode, radius_policy, j_reach_scale,
                         start, /*stop_at_toplevel=*/true,
                         &out, /*topleaf_map=*/nullptr, /*export_out=*/nullptr);
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
    /* drift_particle mutates Pos and KernelRadius (by kernel_radius_drift_factor).
     * The next gpu_ngb_list_build call (for
     * non-Mode-B callers) needs to refresh compact_h from these freshly
     * drifted KernelRadius values. */
    gizmo_mark_kernel_radius_dirty_indices(indices, n);
}
