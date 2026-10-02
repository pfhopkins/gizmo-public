/* mesh/device_tree_walk.h -- the one device tree traversal.
 *
 * A query walks a rank's local tree on the device, and the caller decides what
 * happens when the walk reaches a locally-owned particle.  Everything above
 * that decision -- opening nodes, the periodic overlap test, stepping past
 * imported and pseudo subtrees, recognising a malformed tree -- lives here and
 * is written once.
 *
 * There is exactly one of these.  A second device traversal would mean two
 * copies of the periodic-wrap convention and two copies of the rule for which
 * index classes a walk may follow, in the code where a divergence is hardest to
 * see and most expensive to be wrong about.
 *
 * The caller supplies a leaf policy, an object with
 *
 *     KOKKOS_INLINE_FUNCTION
 *     void visit(int local_index, double qx, double qy, double qz, double reach)
 *
 * called once for every locally-owned particle the walk reaches.  The policy
 * applies its own acceptance test and keeps its own result: the ghost-discovery
 * receiver records the accepted supply slots, while a fused loop evaluates its
 * pair kernel and accumulates.  The traversal itself never looks at particle
 * fields and never decides what a neighbour is.
 *
 * Ownership is the traversal's business, not the policy's, because it is a
 * question about index classes and those are decided here.  The particle slots
 * run past the owned locals into the imported ghosts appended behind them, so
 * the walk stops at `local_particle_slots` and a policy never sees a ghost.
 * The host walker draws the same line in the same place.  A policy that had to
 * remember this would eventually be written by someone who did not, and the
 * result -- every ghost counted a second time -- is silent.
 *
 * Where the line falls is the caller's to choose, because it is a property of
 * the tree it hands over, not of the traversal: a walk that answers for its own
 * rank alone sets `local_particle_slots` to the owned count, while one that is
 * meant to see imported copies too sets it to `particle_slots`.  A walk that
 * searches locally and lets other ranks answer for their own particles wants
 * the first, or it counts the far side twice.  Leaving the field unset is not a
 * third option and does not quietly answer short: it is reported as the same
 * fatal state as a malformed tree.
 *
 * Node geometry comes from the gravity tree's device mirror rather than the
 * managed node arrays, because streaming those from a kernel is memory-bound
 * enough to erase the win.
 *
 * The walk cannot drift a stale node -- that needs a lock -- so it is only
 * legal once the node sweep has certified the tree current and no host lazy
 * drift has happened since.  Callers check that immediately before launching,
 * because a discovery walk earlier in the same exchange can withdraw device
 * legality.
 *
 * Only translation units that compile for the device include this file: the
 * traversal is a Kokkos device function, and the periodic-wrap macros it calls
 * read All.BoxSize, which resolves per translation unit to that unit's All
 * mirror.  Any unit that instantiates this walk must therefore call
 * GIZMO_GPU_ENSURE_ALL_FRESH() at its dispatch boundary, or the box size reads
 * as zero on device and periodic wrapping silently stops working.  Kokkos is
 * included directly rather than through the usual annotation fallback, because
 * this walk calls Kokkos itself: a unit that got the fallback would compile the
 * traversal as host code and then fail on the atomic anyway.
 *
 * HOW A WALK ENTERS THE TREE.  There are two entry points and they differ in
 * two ways, so the caller names which one it is.
 *
 * A walk that resumes from start nodes another rank's walk reached begins at
 * each exported node's children, and stops when it re-enters the top-level
 * tree: the querying rank owns everything above that and has covered it
 * already.  A walk that starts from the root begins at the root node itself,
 * and must descend through the top-level tree, because those regions are its
 * own to search.  Stopping there would end such a walk on its first node.
 *
 * Both forms mirror the host walker in mesh/mode_b_local_walker.cc, which
 * takes the same pair of choices as a start node and a stop_at_toplevel flag.
 * Keep them in step: that walker is the reference this one is checked against.
 *
 * Written by Philip F. Hopkins (phopkins@caltech.edu) for GIZMO. */

#ifndef DEVICE_TREE_WALK_H
#define DEVICE_TREE_WALK_H

#include <Kokkos_Core.hpp>

#include "neighbor_list.h"              /* gx_export_envelope_t, GxDeviceTreeView */
#include "ghost_exchange_functions.h"   /* the canonical-wrap overlap predicate */
#include "../gravity/forcetree.h"       /* BITFLAG_TOPLEVEL, NODE_TYPE_PRESENCE */
#include "../core/timestep_functions.h" /* get_drift_factor_impl, DriftKickTableView:
                                       * the SAME interpolator the node sweep uses
                                       * (gpu_force_drift.cc:184), not a second one */

/* The anomaly codes a walk reports through `anomaly` are declared with the tree view in
 * neighbor_list.h, because the tile walk and the recorder claims report the same states.
 * Every one of them is fatal to the walk's caller, which stops the run. */

/* Which entry point a walk is using.  See the entry discussion at the top of
 * this file; the two forms correspond to the host walker's start node and
 * stop_at_toplevel pair. */
enum class GxWalkEntry {
    SubtreeResume,   /* from another rank's start nodes; stop re-entering the top level */
    LocalRoot        /* from this rank's own root; descend the top level */
};

/* Walk one query against the local tree.
 *
 * `anomaly` reports the single state the host walk treats as fatal, an index in
 * the gap that belongs to no class.  The host stops the run there, so this walk
 * cannot simply stop stepping: that would truncate the query silently and
 * return a short answer that looks complete.  It records the state instead and
 * the caller reproduces the host's stop.
 *
 * Callers use the two wrappers below rather than this directly.  `Entry` is a
 * template argument so the entry tests resolve when the walk is compiled and
 * cost nothing per node.  For a root walk `start_nodes` is unused and may be
 * null: the root comes from the tree view. */
template <GxWalkEntry Entry, class LeafPolicy>
KOKKOS_INLINE_FUNCTION
void gx_device_tree_walk_impl(double qx, double qy, double qz, double reach,
                              const int *start_nodes, int n_start,
                              const GxDeviceTreeView &tree,
                              LeafPolicy &leaf_policy,
                              int *anomaly,
                              unsigned int prune_type_mask)
{
    /* An unfilled view would otherwise answer short in silence, which is the one
     * way this walk can be wrong without anything looking wrong. */
    if(tree.local_particle_slots < 0) {Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_MALFORMED_TREE); return;}

    const int n_entries = (Entry == GxWalkEntry::LocalRoot) ? 1 : n_start;

    for(int k = 0; k < n_entries; k++) {
        int no;
        if(Entry == GxWalkEntry::LocalRoot) {
            /* Enter at the root node itself, so its own overlap test runs and
             * the walk descends the top-level tree.  The host local walk enters
             * the same way. */
            no = tree.node_base;
        } else {
            const int start = start_nodes[k];
            if(start < 0) {break;}                   /* -1 terminates the list */
            /* The start list arrived over MPI, so it is validated rather than trusted. */
            if(start < tree.node_base || start >= tree.pseudo_start) {continue;}
            if(start - tree.node_base >= tree.node_capacity) {   /* precondition leaves this unreachable */
                Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_MALFORMED_TREE);
                break;
            }
            no = tree.node_nextnode[start - tree.node_base];   /* open the exported node */
        }

        while(no >= 0) {
            if(no >= tree.particle_slots && no < tree.node_base) {
                Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_MALFORMED_TREE);   /* caller stops the run */
                break;
            }
            if(no < tree.particle_slots) {
                /* Imported ghosts sit above the owned locals in the same slot
                 * range; step over them rather than reporting them twice. */
                if(no < tree.local_particle_slots) {
                    leaf_policy.visit(no, qx, qy, qz, reach);
                }
                no = tree.nextnode_aux[no];
            } else if(no < tree.pseudo_start) {
                const int kn = no - tree.node_base;
                if(kn < 0 || kn >= tree.node_capacity) {
                    Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_MALFORMED_TREE);
                    break;
                }
                /* Re-entering the top-level tree means this exported branch is
                 * exhausted (the querying rank owns everything above it).  A
                 * walk from the root owns those regions itself and descends. */
                if(Entry == GxWalkEntry::SubtreeResume) {
                    if(tree.node_bitflags[kn] & (1u << BITFLAG_TOPLEVEL)) {break;}
                }
                /* WIDEN-ON-OPEN: the sweep's own expression, evaluated lazily.
                 *
                 * `len` is what the mirror last recorded and `node_ti` is WHEN it
                 * recorded it, so a node that has moved since is re-bounded here
                 * rather than eagerly advanced by a whole-tree sweep. Over-widening
                 * only over-includes and the pair kernel re-gates; under-widening
                 * would drop neighbours silently, so the widening is never allowed
                 * to be negative and a non-finite term is refused outright. */
                double len_eff = (double)tree.node_len[kn];
                if(tree.node_vmax && tree.node_ti && tree.drift_tables_ok) {
                    /* The one widening rule (core/timestep_functions.h): the box grows
                     * by how far its fastest member can have moved on the undilated
                     * clock since the mirror was written.  An invalid value is a
                     * defect and is reported through the channel the caller already
                     * treats as fatal, never narrowed. */
                    const double dl = motion_bound_widening((double)tree.node_vmax[kn], tree.node_ti[kn],
                                                            tree.ti_now, &tree.drift_tables);
                    if(!motion_bound_widening_is_valid(dl)) {
                        Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_MALFORMED_TREE);
                    } else {
                        len_eff += dl;
                    }
                }
                const double hw = 0.5 * len_eff;
                int do_open =
                    gx_extended_overlap_wrap_and_test((double)tree.node_center[kn][0] - qx,
                                                      (double)tree.node_center[kn][1] - qy,
                                                      (double)tree.node_center[kn][2] - qz,
                                                      hw, hw, hw, reach);
                /* TYPE PRUNE.  A node holding none of the types this loop consumes has nothing for
                 * it however close it is, so it is skipped whole rather than descended to have
                 * every leaf rejected one at a time.  The same test the leaf already applies,
                 * moved up to where it can save the descent.
                 *
                 * It reads nothing that drifts -- not a position, not a time, not the touched set
                 * -- so a recording walk and the evaluating walk that follows it make the identical
                 * decision at every node.  That is what lets it be added to a pair of walks whose
                 * agreement about which leaves they reach is a correctness requirement. */
                if(do_open && prune_type_mask && tree.type_mask_trusted) {
                    if(!(NODE_TYPE_PRESENCE(tree.node_bitflags[kn]) & prune_type_mask)) {do_open = 0;}
                }
                if(do_open) {
                    const int child = tree.node_nextnode[kn];
                    /* An imported foreign subtree holds no locally-owned
                     * particles, so there is nothing below it to find. */
                    no = (child >= tree.foreign_base && child < tree.pseudo_start)
                             ? tree.node_sibling[kn] : child;
                } else {
                    no = tree.node_sibling[kn];
                }
            } else {
                /* Pseudo-particle: another rank's subtree root, nothing local
                 * below it.  Step past exactly as the host walk does. */
                no = tree.nextnode_aux[tree.particle_slots + (no - tree.pseudo_start)];
            }
        }
    }
}

/* Resume from the start nodes another rank's walk reached, carried in an
 * export envelope. */
template <class LeafPolicy>
KOKKOS_INLINE_FUNCTION
void gx_device_tree_walk(const struct gx_export_envelope_t &env,
                         const GxDeviceTreeView &tree,
                         LeafPolicy &leaf_policy,
                         int *anomaly,
                         unsigned int prune_type_mask = 0)
{
    gx_device_tree_walk_impl<GxWalkEntry::SubtreeResume>(
        env.pos[0], env.pos[1], env.pos[2], env.h,
        env.nodes, env.n_nodes, tree, leaf_policy, anomaly, prune_type_mask);
}

/* The same resumed walk, for a caller that holds the query and its start nodes
 * as plain values rather than as an envelope.
 *
 * A receiver that has already unpacked its incoming batch has the position, the
 * reach and the node list in hand; rebuilding an envelope around them so this
 * function can take it apart again would be a copy per query to satisfy a
 * signature.  Same entry form, same traversal, same rules -- only the argument
 * shape differs, which is why it shares the name. */
template <class LeafPolicy>
KOKKOS_INLINE_FUNCTION
void gx_device_tree_walk(double qx, double qy, double qz, double reach,
                         const int *start_nodes, int n_start,
                         const GxDeviceTreeView &tree,
                         LeafPolicy &leaf_policy,
                         int *anomaly,
                         unsigned int prune_type_mask = 0)
{
    gx_device_tree_walk_impl<GxWalkEntry::SubtreeResume>(
        qx, qy, qz, reach, start_nodes, n_start, tree, leaf_policy, anomaly, prune_type_mask);
}

/* Search this rank's whole tree for a query of its own.  There is no envelope
 * and no start list: the query carries its own position and reach, and the
 * walk begins at the root.
 *
 * A root walk carries the same preconditions as a resumed one and one more.
 * The tree must have been certified current with no host lazy drift since, and
 * the view must be populated and describe that tree.  In addition the root
 * itself must be a real node -- `node_base` below `pseudo_start`, and a node
 * capacity that covers it -- which holds whenever a tree exists at all, and
 * fails only for an empty or half-built one.  Nothing is checked here: this
 * runs per query on the device, and a caller that cannot honour the first
 * precondition cannot honour this one either. */
template <class LeafPolicy>
KOKKOS_INLINE_FUNCTION
void gx_device_tree_walk_from_root(double qx, double qy, double qz, double reach,
                                   const GxDeviceTreeView &tree,
                                   LeafPolicy &leaf_policy,
                                   int *anomaly,
                                   unsigned int prune_type_mask = 0)
{
    gx_device_tree_walk_impl<GxWalkEntry::LocalRoot>(
        qx, qy, qz, reach, nullptr, 0, tree, leaf_policy, anomaly, prune_type_mask);
}

#endif /* DEVICE_TREE_WALK_H */
