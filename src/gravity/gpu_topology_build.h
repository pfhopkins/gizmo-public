/* gpu_topology_build.h
 *
 * GPU tree-build orchestration: assigns local particles to topleaves via
 * device-side Peano walk, computes 128-bit Morton keys, sorts particles
 * within each topleaf range, and emits internal-node topology directly
 * into the SoA `Nodes_dev` mirror.
 *
 * Internal scratch lives in static SharedSpace buffers, reused across
 * tree builds.  Lifecycle: data path allocates / grows on first call;
 * gpu_topology_build_release() frees at shutdown.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef GIZMO_GPU_TOPOLOGY_BUILD_H
#define GIZMO_GPU_TOPOLOGY_BUILD_H


#ifdef __cplusplus
extern "C" {
#endif

/* 6.5c2 data path: for each particle in [0..npart) (= NumPart for the
 * caller), compute its (Peano, Morton) keys, walk TopNodes to find its
 * topleaf, bucket the particle into that topleaf's range in sorted_idx,
 * and Morton-sort within each range.  After this returns, the scratch
 * accessors below provide:
 *
 *   sorted_idx[topleaf_start[t] .. topleaf_start[t] + topleaf_count[t])
 *       -- particle indices in topleaf t, sorted by Morton key
 *   topleaf_start[NTopleaves] = npart   (end-of-buckets sentinel)
 *
 * Pre-conditions:
 *   - gpu_particles_arena_acquire() has been called (P_dev populated).
 *   - TopNodes / DomainNodeIndex are populated on host (i.e.
 *     force_create_empty_nodes has run).
 *
 * mp (optional): when non-NULL, build the tree over the subset P[mp[i].index]
 * for i in 0..npart-1 (e.g. SUBFIND per-species / collective unbinding). When
 * NULL, build over P[0..npart-1] (full/identity). The pipeline is slot-indexed;
 * the slot->particle map is owned by this TU's scratch and consumed through
 * gpu_topology_emit_bfs (freed by gpu_topology_build_release).
 *
 * Returns 0 on success. */
struct unbind_data;
int gpu_topology_build_data_path(int npart, const struct unbind_data *mp);

/* Read-only accessors.  Return NULL / 0 before data_path has run.
 * Buffers live in SharedSpace -- safe for both host and device reads. */
const int *gpu_topology_build_sorted_idx(void);
const int *gpu_topology_build_topleaf_start(void);   /* [NTopleaves + 1] */
const int *gpu_topology_build_topleaf_count(void);   /* [NTopleaves]     */
/* [npart] -- the top-leaf each particle was bucketed under.  That is the leaf its position falls in,
 * except for a particle kept under the leaf it hung from in the standing tree (see below). */
const int *gpu_topology_build_particle_topleaf(void);

/* BFS topology emission: take the Morton-sorted-per-topleaf
 * data laid out by gpu_topology_build_data_path and emit internal-node
 * topology (center, len, father, suns_backup) directly into the
 * gpu_gravity_tree SoA mirror.
 *
 * Pre-conditions:
 *   - gpu_topology_build_data_path has populated sorted_idx / topleaf_*.
 *   - gpu_gravity_tree_acquire has seeded SoA from CPU AoS for the topnode
 *     range (topleaves' geometry already there from force_create_empty_nodes
 *     + seed_from_aos).
 *   - Numnodestree (host) holds the post-topnode-skeleton end of the SoA;
 *     new GPU-built nodes start at this offset.
 *
 * On success: SoA holds the full inside-topleaf topology.  *new_node_count
 * is set to the new value of Numnodestree (= old value + #nodes emitted).
 *
 * Return codes:
 *    0   success
 *    1   MaxNodes overflow during emission (retry path not yet implemented)
 *    2   collocation detected: a sub-range of >1 particles share full LCP
 *        (= 126 bits).  RNG-fallback branch not yet implemented.
 *    >=3 other failure (allocation, missing dependencies, etc.).
 *
 * `start_node_index` -- the SoA index at which the BFS may begin allocating
 * new internal nodes (= Numnodestree at call time). */
int gpu_topology_emit_bfs(int start_node_index, int *new_node_count_out);

/* Helper: copy GPU-built topology (suns_backup, center, len)
 * back into the AoS Nodes_base[] array for the SoA index range
 * [first_soa_idx, last_soa_idx).  Required while force_update_node_recursive
 * still walks AoS u.suns to set sibling/father for the whole tree; retiring
 * force_update_node_recursive on the GPU compile path will eliminate this
 * writeback.
 *
 * Runs as an OMP host loop -- Nodes_base is host malloc, not SharedSpace.
 * The SoA fields live in SharedSpace UVM so reads incur first-touch page
 * faults but no explicit copy.  Cost is bounded by the number of GPU-built
 * inside-topleaf internal nodes and is one-time per tree-build.
 *
 * Returns 0 on success. */
int gpu_topology_writeback_to_aos(int first_soa_idx, int last_soa_idx);

/* Retained attachment, for a whole-tree rebuild that happens without a domain decomposition on the
 * same step.  Particles drift across top-leaf boundaries between decompositions, so by the time such a
 * build runs some of this rank's particles fall geometrically in top-leaves another rank owns.  They
 * cannot be bucketed there -- the pseudo-particle exchange overwrites that node afterwards and the
 * subtree holding them is left unreachable from the root, absent from every rank's multipole moments.
 *
 * gpu_topology_prepare_retained_attachment() runs BEFORE the build, while the standing tree is still
 * intact, and keeps each such particle under the top-leaf it hung from there.  It also computes the
 * keys and leaves the following gpu_topology_build_data_path() would otherwise compute, so the
 * classification costs no extra pass over the particles.
 *
 * topology_valid: whether the standing tree's Father[] links still describe these particles.  When it
 * is false they are not consulted; the caller is told how many crossers there are so it can restore
 * geometric ownership (a repartition) before building.
 *
 * *n_crossed_out     -- particles whose current top-leaf belongs to another rank.
 * *n_unrecovered_out -- of those, how many no owned top-leaf could be recovered for.  Nonzero with a
 *                       valid standing tree means the tree and the particles disagree; it is a stop,
 *                       not a recoverable state.
 * *n_outside_extent_out -- particles outside the extent the domain was built on.  Their keys would be
 *                       wrong, so none is computed for them; nonzero asks for a full decomposition.
 *
 * Returns 0 on success. */
int gpu_topology_prepare_retained_attachment(int npart, int topology_valid,
                                             long *n_crossed_out, long *n_unrecovered_out,
                                             long *n_outside_extent_out);

/* Drop a prepared plan, so a build the stage above does not cover cannot consume keys and attachments
 * computed for a different set of particles. */
void gpu_topology_forget_prepared(void);

/* Grow the nodes holding a retained particle so their cubes cover where it actually is.  Call after
 * gpu_topology_finalize_father (which establishes the particle Father[] links) and before the moments
 * and the pseudo-particle exchange (which carries the grown top-leaf length to the other ranks).
 * Returns 0 on success; a no-op when nothing was retained. */
int gpu_topology_grow_retained_paths(void);

/* Free internal SharedSpace scratch.  Idempotent. */
/* Finish the leaves that hold more than one particle.  Both are no-ops at a leaf size of one.
 * _assign_fathers runs after Father[] is cleared and before the moments; _materialize runs after the
 * successor threading and before the tree is written back or exported. */
int gpu_leaf_chain_assign_fathers(void);
int gpu_leaf_chain_materialize(void);

void gpu_topology_build_release(void);

#ifdef __cplusplus
}
#endif


#endif /* GIZMO_GPU_TOPOLOGY_BUILD_H */
