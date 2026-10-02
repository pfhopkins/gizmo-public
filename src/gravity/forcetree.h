#ifndef FORCETREE_H
#define FORCETREE_H

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially (condensed, new feedback routines added, many different
 * types of walk and calculations added, structures in memory changed,
 * switched options for nodes, optimizations, new physics modules and
 * calcutions, and new variable/memory conventions added)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 * Mike Grudic has also made major revisions to code the Hermitian calculations and binary timestepping.
 */


#define BITFLAG_TOPLEVEL                   0
#define BITFLAG_DEPENDS_ON_LOCAL_ELEMENT   1
#define BITFLAG_INTERNAL_TOPLEVEL          6
#define BITFLAG_MULTIPLEPARTICLES          7
#define BITFLAG_NODEHASBEENKICKED          8
#define BITFLAG_INSIDE_LINKINGLENGTH       9

/* Which particle TYPES are present below a node, one bit per type, packed into the spare high
 * bits of the same bitflags word the topology flags use.  Stored shifted so it cannot collide
 * with the single-bit flags above; NODE_TYPE_PRESENCE() hands it back in the (1u << Type)
 * convention that ghost_exchange_spec.h and every loop's neighbor_type_mask already use, so a
 * consumer compares the two directly.
 *
 * MEANING: a type is marked when a particle of that type OWNED BY THIS RANK sits below the node.
 * NOT the global union -- the device neighbour walks never descend a foreign subtree and only
 * visit owned leaves, so this is the question they ask, and it is tighter than the union would
 * be.  Two places keep it that way: the LET wire copy clears these bits (so foreign nodes carry
 * none), and the top-node re-sum clears them before re-ORing its local children.  A consumer
 * that genuinely needs remote contents needs a different field, not a widening of this one.
 *
 * It may over-claim (a type that has since gone) -- that only costs a missed prune.  It must
 * never under-claim, which is why every mutation is a monotone OR and only a build or moment
 * refresh recomputes it exactly. */
#define BITFLAG_TYPEPRESENT_SHIFT          16
#define BITFLAG_TYPEPRESENT_MASK           (0x3Fu << BITFLAG_TYPEPRESENT_SHIFT)

/* Read on both the host and inside the device tree walk.  Deliberately a macro rather than an
 * annotated inline: the device annotations are themselves defined by whichever header happens to be
 * reached first (vec3.h defines GIZMO_GPU_FUNCTION to nothing behind a bare #ifndef), so an
 * annotated function here would quietly compile host-only in exactly the translation units that
 * need it on the device.  A bit extract has no reason to carry that risk. */
#define NODE_TYPE_PRESENCE(bitflags) \
    ((((unsigned int) (bitflags)) >> BITFLAG_TYPEPRESENT_SHIFT) & 0x3Fu)

void force_update_tree(void);
void force_refresh_node_moments(void);

void force_flag_localnodes(void);

void *gravity_primary_loop(void *p);

/* Walk the tree once for a packet of up to cap targets (see forcetree.cc). workspace is the
 * calling thread's block of force_treewalk_workspace_bytes_per_thread(cap) bytes. Returns 1
 * when every member's result is written (ninter_out[m] = its interaction count), 0 when a
 * packet of several met a pseudo-particle and wrote nothing (walk each member alone), -1 when
 * a packet of one found the import-detector table full. */
int force_treeevaluate(const int *targets, int n_targets, int cap, int *ninter_out, void *workspace,
                       int *exportflag, int *exportnodecount, int *exportindex);
size_t force_treewalk_workspace_bytes_per_thread(int cap);
int force_treeevaluate_ewald_correction(int target, int *exportflag, int *exportnodecount, int *exportindex);
void force_drift_node(int no, integertime time1);

/*! Host-vs-device routing of the gravity walk: yes when this rank has fewer than
 *  All.GravityHostWalkBelowActive active candidates. A count policy only; the dynamic tree
 *  update routes on its own count and guard in force_update_tree. */
int gravity_walk_route_to_host(long long n_local_active);

/*! Time at which a host lazy node drift was last actually performed on this rank
 *  (-1 = never). The device node-drift sweep skips nodes already at its target time, so
 *  it must never run at a time the host has already drifted to: it would leave those
 *  nodes' SoA mirror holding pre-drift geometry. Read by gpu_force_drift_nodes as an
 *  invariant check, not as a control input. */
integertime force_host_lazy_drift_ti(void);
void force_tree_discardpartials(void);
void force_treeupdate_pseudos(int);
void force_update_pseudoparticles(void);
void force_kick_node(int i, Vec3<MyDouble>& dv);
void force_dynamic_update(void);
void force_dynamic_update_node(int no, int mode, MyFloat *minbound, MyFloat *maxbound);
void force_update_hmax(void);
void force_update_hmax_of_node(int no, int mode);
void force_finish_kick_nodes(void);

/* A particle's velocity was changed outside the kick (a feedback receiver, a
 * swallowed cell's neighbour, a scattered dark matter particle): raise the
 * motion bound of every node above each of idx[0..n), and remember the
 * top-level nodes reached so their new bound reaches the other ranks at the
 * next tree-update phase.  The kick itself does not come through here: it has
 * its own route and exchange (force_kick_node / force_finish_kick_nodes). */
void gravity_note_motion_bound(const int *idx, int n);
/* All ranks, once per reused-tree step, after force_update_tree: every rank
 * learns the raised bound of every top-level node any rank changed since the
 * last flush, and applies it up its own copy of the chain.  Nothing else in
 * the node moves: no momentum, no kick flag, no timestamp. */
void gravity_flush_pending_motion_bounds(void);
/* A rebuild sets every bound afresh, so whatever was pending is void. */
void gravity_clear_pending_motion_bounds(void);
int force_create_empty_nodes(int no, int topnode, int bits, peano1D x, peano1D y, peano1D z, int *nodecount, int *nextfree);
int  force_exchange_pseudodata(void);          /* returns complete() status: nonzero = unmatched (caller skips dependent pseudo-update) */
void force_exchange_pseudodata_issue(void);    /* split for non-blocking overlap with LET */
int  force_exchange_pseudodata_complete(void); /* pair to _issue; nonzero = unmatched (pending==NULL) */
void force_insert_pseudo_particles(void);
void force_add_element_to_tree(int igas, int istar);

/* Record that a particle of this particle's CURRENT type sits below every node from its father to
 * the root, so the per-node type-presence bits stay true when a type changes or a particle is
 * inserted while the tree is standing.  Call it wherever a particle's type becomes final and
 * wherever one is added to a live tree; it is monotone and idempotent, so calling it more often
 * than strictly needed can only over-claim, which is safe. */
void force_tree_note_type_presence(int particle);

void   force_costevaluate(void);
int    force_getcost_single(void);
int    force_getcost_quadru(void);
void   force_resetcost(void);
void   force_setupnonrecursive(int no);
/* foreign_node_slots_exact: LET foreign-node CAPACITY to allocate verbatim; negative derives it,
 * which is what every normal caller wants.  Only the restart read passes a value: the node pointers
 * it is about to deserialize encode the writer's capacity (pseudo-particles start at
 * TreeNodeIndexBase + MaxNodes + MaxForeignNodes), so the reader has to reproduce it exactly. */
void   force_treeallocate(int maxnodes, int tree_particle_slots, int foreign_node_slots_exact = -1);
/* Conservative per-particle radius used to seed Extnodes[no].hmax_per_type[Type]
 * bands. Mode B SYMMETRIC tree-prune reads these bands as an upper bound on
 * any leaf-policy-selectable reach for that type, so the band must dominate
 * every per-particle radius source a Spec's radius_policy can pick at the
 * leaf (P[i].KernelRadius, P[i].AGS_KernelRadius when defined, P[i].ForceSoftening).
 * Leaf-level Mode B predicate / Mode A compact_xyzh still apply the exact policy;
 * over-opening here is correct (extra candidates filter at leaf). Capped at
 * All.MaxKernelRadius to match the legacy band semantics. */
double force_hmax_per_type_particle_radius(int i);

/* Monotonic gravity-tree freshness generations.  treebuild_generation bumps on
 * every successful force_treebuild (topology + Father[] + node structure changed);
 * hmax_refresh_generation bumps at the end of force_update_hmax (ancestor node
 * boxes re-drifted + per-type bands re-seeded after density).  Consumers that
 * cache anything derived from the tree geometry (e.g. the ghost-route fine band)
 * key on BOTH plus All.Ti_Current; a mismatch means rebuild / fail-closed.  These
 * are NOT a substitute for the per-data epoch keys, only the tree-side half. */
long   force_treebuild_generation(void);
long   force_hmax_refresh_generation(void);
void   force_bump_hmax_refresh_generation(void);   /* called by force_update_hmax (separate TU) */

/* Returned by force_treebuild() instead of a node count when the build cannot go ahead on the
   decomposition in place: particles have drifted into top-leaves other ranks own, and the standing
   tree cannot say where they were attached -- either because there is none, or because it could not
   place them.  Nothing has been freed or rebuilt at that point, so the caller restores geometric
   ownership (domain_Decomposition_light with the refinement pass suppressed) and asks again, which
   leaves every particle in a top-leaf its own rank owns and so needs nothing retained.  Raised only
   for whole-tree builds. */
#define FORCE_TREE_NEEDS_OWNERSHIP_RESTORE (-2)
/* Returned instead of a node count when a particle lies outside the extent the domain was built on, so
   its Peano key would name another cell.  Checked before any key is formed and before anything is freed;
   the caller does a full decomposition, which re-measures the extent, and asks again.  Raised only for
   whole-tree builds; takes precedence over FORCE_TREE_NEEDS_OWNERSHIP_RESTORE. */
#define FORCE_TREE_NEEDS_DOMAIN_REBUILD (-4)

int    force_treebuild(int npart, struct unbind_data *mp);

/* Whether the standing tree's Father[] links still describe the particles they were built for, which
   is what a whole-tree rebuild without a decomposition needs in order to keep a drifted particle under
   a top-leaf this rank owns.  Dropped when the tree is freed, and by any reordering of the particle
   array the tree does not follow. */
int    force_tree_global_topology_valid(void);
void   force_tree_invalidate_global_topology(void);

/* Let a particle keep its parent when re-sequencing moves it between slots, so the record still
   describes the particles at the next rebuild.  Swaps only Father[]; drops the record if either slot
   is outside the tree's particle slots, and does nothing when there is no record to keep.  Not called
   where MAINTAIN_TREE_IN_REARRANGE is doing the full pointer repair, which already swaps Father. */
void   force_tree_swap_attachment_slots(int i, int j);
int    force_treebuild_single(int npart, struct unbind_data *mp);
int    force_treeevaluate_direct(int target, int mode);
void   force_treefree(void);
int    force_tree_is_allocated(void);   /* nonzero while tree storage is held */

/*! Give the foreign-node range storage for the `foreign_needed` nodes this rank is about to
 *  receive, once the LET exchange has counted them.  force_treeallocate leaves that storage
 *  at zero because the count is not knowable when the tree is allocated.  Call once per tree
 *  build, before the first foreign node is installed; a rank importing nothing passes 0 and
 *  does nothing.  Leaves the tree untouched and returns nonzero if the memory is not there. */
int    force_tree_grow_foreign_storage(long long foreign_needed);

/*! Can a particle created at this index be carried by the tree that is standing right now?
 *  Father[] and Nextnode[]'s particle segment span All.TreeParticleSlots entries and are NOT resized
 *  between rebuilds, so a slot being available in P[] does not by itself mean the live tree can index
 *  it: inserting past those arrays would write out of bounds, and every later pass that reads a
 *  particle's parent would read out of bounds too.  Creation sites ask this in addition to their
 *  existing capacity checks, and decline the way they already decline when storage is short.
 *  A tree built by force_treeallocate sizes these arrays from All.MaxPartExpandable, the run's
 *  ceiling on the particle capacity, so the answer is currently yes for every index P[] can hold and
 *  this costs one comparison.  It is kept rather than deleted because it states the requirement at
 *  the sites that depend on it: if the ceiling ever became something a capacity could reach, these
 *  are exactly the places that must decline.  With no tree standing the answer is also yes, because
 *  there is nothing to outgrow -- not a licence to insert into a tree that is not there, which
 *  force_add_element_to_tree refuses on its own. */
static inline int gizmo_particle_index_fits_live_tree(int index)
{
    return (!force_tree_is_allocated()) || (index < All.TreeParticleSlots);
}
void   force_update_node(int no, int flag);
void   force_update_size_of_parent_node(int no);

void   dump_particles(void);

/* mesh/ngb.cc retired: ngb_treebuild/ngb_treefind_* all dead on the Kokkos path.
   ngb_treebuild() callers replaced with force_treebuild(NumPart, NULL) directly. */

#ifdef BOX_PERIODIC
/* Ewald octant table size. Declared under BOX_PERIODIC to match `EN` in forcetree.cc: the tables
 * + the CPU Ewald functions + the shared interp helper (gravtree_ewald.h) all compile under
 * BOX_PERIODIC (even under GRAVITY_NOT_PERIODIC, where the correction is compiled-but-dead). */
#define GIZMO_EWALD_EN 64
#endif
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
/* Ewald correction table accessor. Returns flat pointers of length (EN+1)^3
 * to the four static look-up tables (fcorrx/y/z/potcorr) inside forcetree.cc,
 * plus fac_intp (= 2*EN/All.BoxSize). Used by gpu_gravtree.cc to mirror the
 * tables into SharedSpace once after ewald_init(). */
void gizmo_get_ewald_tables(const MyFloat **fcorrx_out,
                            const MyFloat **fcorry_out,
                            const MyFloat **fcorrz_out,
                            const MyFloat **potcorr_out,
                            double *fac_intp_out);
#endif

#endif



