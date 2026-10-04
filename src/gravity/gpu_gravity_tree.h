/* gpu_gravity_tree.h
 *
 * GPU-resident SoA mirror of the gravity tree (NODE / extNODE arrays).
 *
 * Build still happens on CPU (force_treebuild → Nodes_base AoS); this layer
 * provides the SoA mirror so a GPU walk kernel can do coalesced
 * reads. Retired once the tree builds directly on the GPU.
 *
 * Lifetime: acquire() takes Nodes_host + Extnodes_host pointers and a
 * capacity. If the SoA already mirrors the latest CPU tree, returns the
 * cached pointers without re-copying. Otherwise reseeds. invalidate() marks
 * the mirror stale (call after force_treebuild, force_update_node_recursive,
 * domain decomp). release() frees the SharedSpace storage.
 *
 * Fields mirrored: the subset the walk reads. center/len for opening, s/mass
 * for force, sibling/nextnode for traversal, bitflags for opening type,
 * maxsoft for adaptive softening. The optional payload families
 * (RT_USE_GRAVTREE, SINK_*, tidal tensor, DM_SCALARFIELD_SCREENING) are
 * mirrored too, gated by the same #ifdefs as the AoS NODE definition; the
 * GPU gravity walk consumes them directly from the SoA.
 *
 * Memory model: GIZMO_KOKKOS_SHARED_SPACE (UVM) so host writes during
 * acquire are visible to device kernels without explicit deep_copy. Per-field
 * SoA arrays (Vec3 fields kept as Vec3 arrays for now; can split into x/y/z
 * if profiling shows coalescing benefit).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef GIZMO_GPU_GRAVITY_TREE_H
#define GIZMO_GPU_GRAVITY_TREE_H

struct NODE;
struct extNODE;


#ifdef __cplusplus
extern "C" {
#endif

/* SoA view exposed to GPU kernels. All pointers live in SharedSpace; indices
 * match the AoS Nodes[] convention (callers index by [no - All.TreeNodeIndexBase] when
 * using the base array, or by absolute Nodes[] index after applying the
 * NTopnodes offset — the exact indexing convention is not yet finalized). */
struct gpu_gravity_tree_soa_t {
    /* Geometric / opening criterion — KEEP MyFloat (double) always. Node
     * center and sidelength drive the opening criterion directly; reducing
     * precision here would lose geometric accuracy independent of flag state. */
    Vec3<MyFloat>  *center;     /* geometric center of node */
    MyFloat        *len;        /* sidelength */
    /* Multipole (currently monopole): MyGravFloat (= float when
     * GIZMO_MIXED_PRECISION_GRAVITY is set, double otherwise). Force-kernel
     * accumulators and moment storage. Narrowing cast happens at seed time
     * from the MyFloat-typed NODE. */
    Vec3<MyGravFloat> *s;       /* center of mass */
    MyGravFloat       *mass;    /* total mass */
    /* Walk traversal — integer bookkeeping, untouched by precision flag. */
    int            *sibling;
    int            *nextnode;
    int            *father;     /* parent node index (needed by GPU moment-refresh dependency-counter walk) */
    unsigned int   *bitflags;
    /* Force kernel */
    MyGravFloat    *maxsoft;
    long           *N_part;
    /* Foreign-leaf identity sidecar (GPU mirror of the host ForeignLeaf* arrays).  Sized
     * AllocatedForeignNodes and indexed by foreign_slot = no - (TreeNodeIndexBase+MaxNodes) == (node SoA idx) -
     * MaxNodes -- a DIFFERENT index than every other array above (which use no - TreeNodeIndexBase).  The GPU
     * walk computes foreign_slot explicitly and bounds-checks it so the two conventions can never be
     * confused.  Populated for the installed foreign range by gpu_scatter_foreign_to_soa. */
    int            *foreign_leaf_tag;   /* LET_LEAF_TAG_* : 0 = descendable node, 1 = real
                                         * single-particle leaf, 2 = truncated aggregate */
    int            *foreign_leaf_type;  /* source particle Type   -> ptype_sec */
    MyFloat        *foreign_leaf_zeta;  /* source particle AGS_zeta -> zeta_sec (double; preserves reciprocity) */
    MyFloat        *foreign_leaf_soft;  /* source particle ForceSoftening -> h_p (pure, NOT node maxsoft) */
    int             foreign_leaf_cap;   /* allocated length of the foreign_leaf_* arrays (== AllocatedForeignNodes) */
    int             nnodes;     /* number of valid entries */
    /* Nextnode[] mirror — used for particle-level traversal: when the walk
     * lands on `no < TreeParticleSlots`, advance via nextnode_aux[no]. Sized for
     * TreeParticleSlots + NTopnodes so the pseudo-particle region is addressable
     * (though Tier-1 GPU walk aborts on pseudo-particle rather than
     * following it). */
    int            *nextnode_aux;
    int             nextnode_aux_size;
    /* Build-time suns[8] per internal node, snapshotted before
     * force_update_node_recursive overwrites the union with the d struct.
     * Sized [capacity * 8].  The GPU nextnode-threading kernel reads from
     * here to reconstruct DFS-order links in parallel. */
    int            *suns_backup;

    /* --- Optional payloads: gated by the same flags as the AoS
     *     NODE definition in allvars.h. Each block is present iff the host
     *     NODE carries the field. Moment/force storage → MyGravFloat.
     *     CHIMES arrays stay double — chemistry tolerances force it. */
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    MyGravFloat    *gasmass;          /* [nnodes] */
#endif
#ifdef RT_USE_GRAVTREE
    MyGravFloat    *stellar_lum;      /* flat [nnodes * N_RT_FREQ_BINS] */
#ifdef CHIMES_STELLAR_FLUXES
    double         *chimes_stellar_lum_G0;  /* flat [nnodes * CHIMES_LOCAL_UV_NBINS] — keep double */
    double         *chimes_stellar_lum_ion; /* flat [nnodes * CHIMES_LOCAL_UV_NBINS] — keep double */
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<MyGravFloat> *rt_source_lum_s;  /* from NODE */
    Vec3<MyGravFloat> *rt_source_lum_vs; /* from extNODE */
#endif
#ifdef SINK_PHOTONMOMENTUM
    MyGravFloat       *sink_lum;
    Vec3<MyGravFloat> *sink_lum_grad;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyGravFloat    *cr_injection;
#endif
#ifdef SINK_CALC_DISTANCES
    MyGravFloat       *sink_mass;
    Vec3<MyGravFloat> *sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
    Vec3<MyGravFloat> *sink_vel;
#endif
#if defined(SPECIAL_POINT_MOTION)
    Vec3<MyGravFloat> *sink_acc;
#endif
#if defined(SINK_NODE_MOTION_TRACKED)
    int            *N_SINK;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    MyGravFloat    *MaxFeedbackVel;
#endif
#endif
    /* Unconditional Extnodes mirrors. Needed by the GPU moment
     * kernel which replaces force_update_node_recursive; these are
     * always computed during moment accumulation regardless of which physics
     * flags are on. The prior SINK_DYNFRICTION_FROMTREE/COMPUTE_JERK_IN_GRAVTREE
     * guard on node_vs was removed — vs is now always present. */
    Vec3<MyGravFloat> *node_vs;       /* mirror of Extnodes[].vs */
    MyGravFloat       *hmax;          /* Extnodes[].hmax — gas kernel extent */
    MyGravFloat       *vmax;          /* Extnodes[].vmax */
    MyGravFloat       *divVmax;       /* Extnodes[].divVmax */
    /* Per-node Ti_current, mirrored for the device.
     *
     * A ONEWAY device walk cannot take a lock, so it cannot drift a node it
     * reaches; instead it widens the node's own bound by how far that node could
     * have moved since it was last advanced -- and that needs the node's time,
     * which no other mirrored field carries.  integertime, not a float: it is
     * compared for equality against All.Ti_Current and fed to get_drift_factor,
     * and a rounded copy of either would be a different node's answer. */
    integertime       *node_ti;       /* mirror of Nodes[].Ti_current */
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    /* 6-component symmetric tensor moment: [xx,yy,zz,xy,xz,yz]. Flat layout
     * [nnodes * 6]. Accumulated by GPU moment kernel; walk consumption not
     * yet wired in (walk guard unchanged). */
    MyGravFloat       *tidal_tensorps;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    MyGravFloat       *mass_dm;
    Vec3<MyGravFloat> *s_dm;
    Vec3<MyGravFloat> *vs_dm;         /* from Extnodes */
#endif
};

/* Acquire SoA mirror sized for at least min_nodes.  Simplified to
 * a pure capacity-grow / pointer-grab.  No AoS->SoA seeding ever happens
 * here — the build pipeline (gpu_nextnode_backup_suns -> emit_bfs ->
 * finalize_father -> finalize_sibling -> moment_refresh -> nextnode_thread)
 * populates the SoA end-to-end, and the GPU drift kernel
 * (gpu_force_drift_nodes) writes both UVM AoS and the SoA mirror in one
 * pass.  Nodes_host / Extnodes_host are unused (kept in the signature for
 * source compatibility; pass Nodes_base / Extnodes_base or NULL). */
void gpu_gravity_tree_acquire(int min_nodes,
                              struct NODE    *Nodes_host,
                              struct extNODE *Extnodes_host);

/* Grow the mirror to min_nodes slots while PRESERVING what it holds -- the counterpart to
 * acquire() for the one point that runs after the build kernels: the LET exchange, which learns
 * how many foreign nodes this rank receives only once the counts have been exchanged, and must
 * add room for them without disturbing the local tree already in the mirror.  Transactional:
 * on failure the previous arrays are still installed and valid.  Returns 1 on success, 0 on
 * failure.  Callers: force_tree_grow_foreign_storage only. */
int gpu_gravity_tree_grow_foreign(int min_nodes);

/* Nextnode[] is allocated in SharedSpace by force_treeallocate
 * (forcetree.cc).  This setter aliases soa->nextnode_aux to that pointer; no
 * separate buffer, no per-walk memcpy.  Pass NULL/0 from force_treefree to
 * clear the alias before the underlying buffer is freed. */
void gpu_gravity_tree_alias_nextnode(int *Nextnode_host, int n);

/* Free SharedSpace storage, plus the drift and moment-refresh pools keyed to it.
 * Called from force_treefree(), so the mirror lives exactly as long as the AoS tree
 * it derives from: it is reallocated once per tree epoch, not held for the run.
 * Nothing may read the mirror between a release and the next build — acquire() marks
 * the arrays valid without seeding them, so a read there sees uninitialized topology
 * rather than the NULL a stale-pointer check would catch. */
void gpu_gravity_tree_release(void);

/* Accessors. Returns NULL pointers / 0 capacity when not held or stale. */
struct gpu_gravity_tree_soa_t *gpu_gravity_tree_soa(void);
int gpu_gravity_tree_capacity(void);
int gpu_gravity_tree_valid(void);
/* The mirror slot of tree node `no`, or -1 when the mirror has none.  The node index range
 * (MaxNodes + AllocatedForeignNodes) can be larger than the mirror, so every host write into
 * the mirror is bounded by the mirror's own capacity through this one test. */
int gpu_gravity_tree_mirror_slot(int no);

/* GPU pre-walk drift kernel — replaces the host loop in
 * gpu_gravtree_walk_primary that called force_drift_node + mark_dirty per
 * stale-Ti_current node.  Mutates Nodes[]/Extnodes[] (UVM) and SoA mirrors
 * in a single kernel; no AoS->SoA reseed afterwards.  Returns 0 on success,
 * nonzero on internal error (SoA not ready). */
/* Sweep node geometry current, also rewriting the SoA mirror of nodes that
   are ALREADY at the target time. A host lazy drift advances a node without
   touching its mirror, so those two can disagree; the ordinary sweep skips
   such nodes and preserves the disagreement. Pass nonzero to pay for a full
   mirror rewrite instead of declining. */
int gpu_force_drift_nodes_ex(integertime time1, int refresh_mirrors_already_current);
int gpu_force_drift_nodes(integertime time1);
void gpu_force_drift_release(void);
void gpu_gravtree_tables_release(void);   /* frees the drift/gravkick table mirror the walk TU keeps */

/* Record that the SoA+AoS node geometry is drifted to `ti` (snapshots the
 * current treebuild generation).  Called by the drift sweep on success so every
 * sweep caller records certification.  The stamp is read by
 * gpu_gravity_tree_nodes_current_at below, and is invalidated on SoA
 * realloc/free/rebuild. */
void gpu_gravity_soa_mark_drift_certified(integertime time1);

/* Record that the tree was BUILT current at `ti` (snapshots the treebuild
 * generation).  Called by the MAIN-STEP tree build only, after force_treebuild
 * returns and only when a full drift has proven the particles reached `ti`.
 * ⛔ Not from force_treebuild itself: that routine also builds group-local and
 * subset trees (subfind), whose geometry must never certify the step's device
 * consumers, and it cannot tell which kind of build it is being asked for. */
void gpu_gravity_tree_mark_born_current(integertime ti);

/* Retire both currency records.  Called by tree teardown: node arrays going
 * away must not leave a record claiming their geometry is current. */
void gpu_gravity_tree_invalidate_currency(void);

/* THE question a device consumer of node geometry should ask: is that geometry
 * current at `ti`?  True if the drift sweep certified it OR the tree was built
 * current at it, AND the mirror those records describe still exists.  Asking
 * only whether the drift sweep ran would miss a tree freshly built at `ti`, which
 * is fully current with no sweep at all.  The converse does not hold: a host tree
 * update (TreeUpdateHostBelowActive) on a reused tree advances only the nodes its
 * active elements reach, so after one only those nodes are current and this
 * answers false. */
int gpu_gravity_tree_nodes_current_at(integertime ti);

/* ============================================================================
 * THE NODE DIRTY SET — class-(a) mirror repair, O(Ndirty).
 *
 * `force_drift_node` advances a node's AoS geometry WITHOUT writing its device
 * mirror, so every such call leaves one node whose mirror is behind its own AoS.
 * Repairing that by sweeping the whole tree costs O(Nnodes) to fix a set measured
 * at ~9,200 nodes per rank per span -- 0.46% of the ~1.99M mirrors a sweep
 * rewrites. This records exactly which nodes are behind, so only those are redone.
 *
 * SAME ALGORITHM AS THE PARTICLE TOUCHED SET (mesh/neighbor_list.h), for the same
 * reason: a generation STAMP that is never cleared, a COMPACTED list that makes
 * the set enumerable in O(Ndirty), and an ATOMIC append cursor. A plain dirty BYTE
 * array cannot be enumerated without an O(Nnodes) scan -- which reintroduces the
 * cost being removed -- and an unsynchronised append list is a data race.
 * It is a separate INSTANCE, not shared storage: this one is owned by gravity
 * (force_drift_node has seven callers, six of them outside Mode-D) and its epoch
 * runs from a host advance until that node's mirror repair completes, a different
 * lifetime from the particle set's one-call span.
 *
 * CAPACITY is MaxNodes + AllocatedForeignNodes -- this rank's LIVE installed
 * range, the same bound the device walk tests, and NOT the run-wide ceiling
 * MaxForeignNodes (the worst rank's import, deliberately generous). Foreign nodes
 * ARE dirtied in practice: 274 per rank per span, present in 91.5% of spans.
 * ========================================================================== */
/* HOST ONLY, storage and both sides of it: the host lazy drift claims, and the claims are
 * answered on the host before each device gravity walk.  No kernel may read the set; its
 * pages then never migrate, every step the host walk runs.
 *
 * What the answer reads, snapshotted once the claim phase is over.  `usable` is 0 when the
 * fail-safe fired this epoch or the storage is missing, and the caller must then sweep. */
struct gpu_node_dirty_view_t {
    const int *list;    /* [count] claimed node indices */
    int        count;
    int        cap;     /* the stamp range: index - base must lie in [0, cap) */
    int        base;    /* All.TreeNodeIndexBase */
    int        usable;
};
struct gpu_node_dirty_view_t gpu_node_dirty_view(void);
void gpu_node_dirty_begin_epoch(void);              /* the claims are ANSWERED: fresh epoch */
void gpu_node_dirty_claim(int no);                  /* host claim, from force_drift_node */
int  gpu_node_dirty_count(void);                    /* claims outstanding in this epoch */
/* Bring every listed node current at `ti` -- drifting the ones behind it and publishing every
 * mirror field the gravity walk reads -- instead of sweeping the whole tree.  Defined beside
 * the sweep (gpu_force_drift.cc) because it runs the sweep's own per-node units.
 * 0 = the listed set stands at `ti`; 1 = the caller must take the full sweep. */
int  gpu_node_dirty_bring_gravity_current(integertime time1);
/* The per-node work of the sweep for an explicit list of node indices (mirror slot = index - base,
 * slot < cap): drift each listed node behind `time1` and publish its mirror.  0 = done, 1 = the
 * caller must sweep instead.  Called only by the device tree update, with the claim list it wrote on the
 * device in its own scratch block -- never with the host's node dirty set, which the host answers above. */
int  gpu_device_node_list_bring_current(const int *list, int n, int base, int cap, integertime time1);
void gpu_node_dirty_grow_to(int cap);   /* keep the set as large as the mirror when foreign storage grows */
void gpu_node_dirty_release(void);
long long gpu_node_dirty_unsafe_events(void);   /* fail-safe firings, run-total; a silent
                                                   permanent revert to sweeping must be visible */

/* Is the device-visible geometry safe for a ONEWAY walk that WIDENS ON OPEN?
 *
 * Deliberately NOT gpu_gravity_tree_nodes_current_at(): that certifies "every
 * mirror is current", which widen-on-open does not need and cannot achieve --
 * reusing it would either keep forcing the sweep or silently weaken its meaning
 * for gravity and the other sweep callers, none of which widen. This is an
 * ADDITIONAL predicate, never a redefinition. */
int  gpu_gravity_tree_oneway_safe_at(integertime ti);


/* GPU moment-refresh kernel. Computes local-tree node moments
 * (mass, COM, vs, hmax, vmax, divVmax, maxsoft, bitflags + all conditional
 * payloads) directly on the device, using dependency-counter atomics on
 * device-local scratch and bulk seed of the SharedSpace SoA.
 *
 * After the kernel returns, the SoA is fully populated for nodes
 * [TreeNodeIndexBase, TreeNodeIndexBase+Numnodestree) AND the Nodes[]/Extnodes[] AoS arrays
 * are written back so that the CPU pseudo-particle path
 * (force_exchange_pseudodata + force_treeupdate_pseudos) sees identical
 * values to what it would have produced.
 *
 * `active_root_node` reserved for a future subtree-hint optimization. Initial
 * callers pass -1 (= whole tree from TreeNodeIndexBase..TreeNodeIndexBase+Numnodestree).
 *
 * Returns 0 on success, nonzero on failure (allocation, bad state, etc). */
int gpu_moment_refresh(int active_root_node);

/* Free the persistent gpu_moment_refresh scratch pools (source-input buffers +
 * Father mirror).  Reached through gpu_gravity_tree_release(), so these pools
 * follow the tree epoch: dropped when the tree is freed, regrown by the next
 * refresh. */
void gpu_moment_refresh_release(void);

/* Bulk write SoA[k=0..n) back into Nodes[TreeNodeIndexBase+k] / Extnodes[TreeNodeIndexBase+k].
 * Invokes from gpu_moment_refresh(); declared here so unit-test scaffolding
 * could call it directly. Caller must have a valid SoA acquired. */
void gpu_moment_writeback_to_aos(int n);

/* Snapshot Nodes_base[k].u.suns[0..7] for k in [0..n) into the
 * SoA's suns_backup buffer.  MUST be called BEFORE force_update_node_recursive
 * overwrites the union with the d struct.  Idempotent; safe to call multiple
 * times during the force_treebuild retry loop (each call refreshes from AoS).
 * Allocates the suns_backup buffer on first use if needed. */
void gpu_nextnode_backup_suns(int n);

/* GPU nextnode-threading kernel.  Recomputes the DFS pre-order
 * `nextnode` link for each internal node (in soa->nextnode + AoS Nodes[].u.d.nextnode)
 * and `Nextnode` for each particle / pseudo-particle (in soa->nextnode_aux + AoS Nextnode[]).
 * Inputs: SoA's suns_backup (populated via gpu_nextnode_backup_suns), sibling[].
 * Returns 0 on success, nonzero on failure. */
int gpu_nextnode_thread(void);

#ifdef __cplusplus
}
#endif


#endif /* GIZMO_GPU_GRAVITY_TREE_H */
