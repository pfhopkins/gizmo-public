/* gpu_force_drift.cc
 *
 * GPU pre-walk drift kernel: replaces the CPU host loop in
 * gpu_gravtree_walk_primary that called force_drift_node + mark_dirty for
 * every node with stale Ti_current.  This kernel:
 *   - mutates Nodes[]/Extnodes[] (UVM AoS) directly inside a Kokkos kernel,
 *   - mirrors the same fields into the SoA used by the GPU walk,
 *   - has zero CPU-side AoS->SoA reseed (the entire dirty_[]/seed_dirty_/
 *     seed_node_/mark_dirty machinery is retired by this commit).
 *
 * Drift factor: the same interpolator the host uses (core/timestep_functions.h).
 * Non-cosmological is a trivial closed form; cosmological reads a SharedSpace mirror
 * of DriftTable and GravKickTable, refreshed from the host on each top-level call.
 * The mirror is 2 x DRIFT_TABLE_LENGTH doubles so the refresh cost is in the noise.
 *
 * Timestep dilation: under USE_TIMESTEP_DILATION_FOR_ZOOMS the dispatcher below
 * fills a per-node dilation_dev cache on the host before the launch; all other
 * configs are provably dilation==1 and skip the cache entirely. Note the node
 * factor itself has a standing gap under SPECIAL_POINT_WEIGHTED_MOTION, where
 * nodes drift undilated -- see return_node_timestep_dilation_factor_P.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../core/timestep_functions.h"
#include "../declarations/gpu_error_check.h"
#include "../system/gpu_particles_arena.h"
#include "gpu_gravity_tree.h"
#include "../declarations/gpu_recorder_claim.h"   /* the node dirty set's view, shared with its device claim */
#include "forcetree.h"


/* USE_TIMESTEP_DILATION_FOR_ZOOMS: node-indexed dilation is supported via a
 * per-call host-pre-compute cache. The dispatcher below allocates a SharedSpace
 * `dilation_dev` array of length n_local_nodes + n_foreign_nodes, populates it
 * on the host using return_node_timestep_dilation_factor(no), then the GPU
 * drift kernel reads the cached value. */

/* --- Drift/GravKick table mirror (cosmological only) ---------------------
 * The storage is per-TU because a device symbol cannot be shared between
 * translation units without relocatable device code; the allocate-and-fill policy
 * lives in the arena so this path and the batched particle drift build the view the
 * same way. See drift_kick_table_mirror_refresh. */
static double *drift_kick_table_dev_ = NULL;   /* SharedSpace, 2 * DRIFT_TABLE_LENGTH doubles */

/* --- dispatcher ---------------------------------------------------------- */

/* ============================================================================
 * THE TWO SHARED PER-NODE UNITS.  The full sweep below and the SUBSET consumer both run
 * exactly these, so a field added to one is added to both by construction; a hand-copied
 * field list in a second consumer is how a subset path silently goes stale.
 * ========================================================================== */

/* Every device-visible mirror field whose canonical AoS value can change BETWEEN TREE
 * BUILDS: drift moves the positions and lengths, a kick raises the velocities and vmax.
 * Everything else the walk reads (mass, maxsoft, the topology, the payload moments) is
 * build- or moment-refresh-owned, and is deliberately absent. */
struct gpu_node_mirror_ptrs_t {
    MyFloat           *len;
    integertime       *node_ti;
    Vec3<MyGravFloat> *s;
    Vec3<MyGravFloat> *vs;
    MyGravFloat       *hmax;
    MyGravFloat       *vmax;
    unsigned int      *bitflags;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<MyGravFloat> *rt_s;
    Vec3<MyGravFloat> *rt_vs;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Vec3<MyGravFloat> *s_dm;
    Vec3<MyGravFloat> *vs_dm;
#endif
#ifdef SINK_NODE_MOTION_TRACKED
    Vec3<MyGravFloat> *sink_pos;
    Vec3<MyGravFloat> *sink_vel;
#endif
};

static inline struct gpu_node_mirror_ptrs_t gpu_node_mirror_ptrs(struct gpu_gravity_tree_soa_t *soa)
{
    struct gpu_node_mirror_ptrs_t m;
    m.len = soa->len; m.node_ti = soa->node_ti; m.s = soa->s; m.vs = soa->node_vs;
    m.hmax = soa->hmax; m.vmax = soa->vmax; m.bitflags = soa->bitflags;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    m.rt_s = soa->rt_source_lum_s; m.rt_vs = soa->rt_source_lum_vs;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    m.s_dm = soa->s_dm; m.vs_dm = soa->vs_dm;
#endif
#ifdef SINK_NODE_MOTION_TRACKED
    m.sink_pos = soa->sink_pos; m.sink_vel = soa->sink_vel;
#endif
    return m;
}

/* Advance ONE node's canonical AoS state to `ti_target`.
 * `fold_kick` is explicit rather than inferred: the host (force_drift_node) does NOT fold a
 * pending kick on a node that already stands at the target time -- it returns early -- while
 * this sweep's refresh pass falls through at dt = 0 and does. Both reach the same state once
 * the next real drift runs, but the node velocity a walk reads in between differs, so the
 * caller states which behaviour it wants instead of inheriting it. */
static KOKKOS_INLINE_FUNCTION void
gpu_node_drift_apply(struct NODE *Nodes_uvm, struct extNODE *Extnodes_uvm, int no,
                     integertime ti_target, double dt_drift, double dt_drift_hmax, double dt_widen,
                     int fold_kick)
{
    /* Fold a pending kick into vs and clear dp -- only when the caller asks (see above). */
    if(fold_kick && (Nodes_uvm[no].u.d.bitflags & (1u << BITFLAG_NODEHASBEENKICKED))) {
        double mass = (double) Nodes_uvm[no].u.d.mass;
        double fac  = (mass > 0) ? (1.0 / mass) : 0.0;

#ifdef RT_SEPARATELY_TRACK_LUMPOS
        double l_tot = 0.0;
        for(int b = 0; b < N_RT_FREQ_BINS; b++) {l_tot += (double)Nodes_uvm[no].stellar_lum[b];}
        double fac_lum = (l_tot > 0) ? (1.0 / l_tot) : 0.0;
#endif
#ifdef DM_SCALARFIELD_SCREENING
        double mass_dm = (double) Nodes_uvm[no].mass_dm;
        double fac_dm  = (mass_dm > 0) ? (1.0 / mass_dm) : 0.0;
#endif

        for(int j = 0; j < 3; j++) {
            Extnodes_uvm[no].vs[j] = (MyFloat)((double)Extnodes_uvm[no].vs[j] + fac * (double)Extnodes_uvm[no].dp[j]);
            Extnodes_uvm[no].dp[j] = 0;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            Extnodes_uvm[no].rt_source_lum_vs[j] = (MyFloat)((double)Extnodes_uvm[no].rt_source_lum_vs[j]
                                                   + fac_lum * (double)Extnodes_uvm[no].rt_source_lum_dp[j]);
            Extnodes_uvm[no].rt_source_lum_dp[j] = 0;
#endif
#ifdef DM_SCALARFIELD_SCREENING
            Extnodes_uvm[no].vs_dm[j] = (MyFloat)((double)Extnodes_uvm[no].vs_dm[j] + fac_dm * (double)Extnodes_uvm[no].dp_dm[j]);
            Extnodes_uvm[no].dp_dm[j] = 0;
#endif
        }
#ifdef SINK_NODE_MOTION_TRACKED
        /* Mirrors the host fold (forcetree_update.cc): normalised by sink_mass, not mass, and
           consumed before the kicked bitflag is cleared below. */
        {
            double sink_mass = (double) Nodes_uvm[no].sink_mass;
            double fac_sink  = (sink_mass > 0) ? (1.0 / sink_mass) : 0.0;
            for(int j = 0; j < 3; j++) {
                Nodes_uvm[no].sink_vel[j] = (MyFloat)((double)Nodes_uvm[no].sink_vel[j] + fac_sink * (double)Extnodes_uvm[no].sink_dp[j]);
                Extnodes_uvm[no].sink_dp[j] = 0;
            }
        }
#endif
        Nodes_uvm[no].u.d.bitflags &= (~(1u << BITFLAG_NODEHASBEENKICKED));
    }

    /* Apply drift to s, len, hmax. */
    for(int j = 0; j < 3; j++) {
        Nodes_uvm[no].u.d.s[j] = (MyFloat)((double)Nodes_uvm[no].u.d.s[j] + (double)Extnodes_uvm[no].vs[j] * dt_drift);
#ifdef SINK_NODE_MOTION_TRACKED
        Nodes_uvm[no].sink_pos[j] = (MyFloat)((double)Nodes_uvm[no].sink_pos[j] + (double)Nodes_uvm[no].sink_vel[j] * dt_drift);
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Nodes_uvm[no].s_dm[j]  = (MyFloat)((double)Nodes_uvm[no].s_dm[j]  + (double)Extnodes_uvm[no].vs_dm[j] * dt_drift);
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Nodes_uvm[no].rt_source_lum_s[j] = (MyFloat)((double)Nodes_uvm[no].rt_source_lum_s[j]
                                            + (double)Extnodes_uvm[no].rt_source_lum_vs[j] * dt_drift);
#endif
    }
    Nodes_uvm[no].len = (MyFloat)((double)Nodes_uvm[no].len
                                  + TREE_DRIFT_VELOCITY_PREFAC * (double)Extnodes_uvm[no].vmax * dt_widen);

    {
        double exp_arg = (double)Extnodes_uvm[no].divVmax * dt_drift_hmax / (double)NUMDIMS;
        if(exp_arg < -1.0) {exp_arg = -1.0;}
        if(exp_arg >  1.0) {exp_arg =  1.0;}
        double decay_fac = exp(exp_arg);
        if(Extnodes_uvm[no].hmax > 0) {
            Extnodes_uvm[no].hmax = (MyFloat)((double)Extnodes_uvm[no].hmax * decay_fac);
        }
        /* Mode B per-type bands: upward-only inflate.
         * Bands include static-ish sources (P[j].ForceSoftening) that
         * don't shrink under drift, so decaying below the actual FS value
         * would under-bound the node-prune. force_update_hmax() re-grows
         * bands per-particle each call; we just must not shrink them
         * here. Scalar `hmax` keeps its legacy bidirectional decay above.
         * Without this guard, expansion regions (positive divVmax) would
         * fail to track via the upward branch. */
        if(decay_fac > 1.0) {
            for(int t = 0; t < 6; t++) {
                if(Extnodes_uvm[no].hmax_per_type[t] > 0) {
                    Extnodes_uvm[no].hmax_per_type[t] = (MyFloat)((double)Extnodes_uvm[no].hmax_per_type[t] * decay_fac);
                }
            }
        }
    }

    Nodes_uvm[no].Ti_current = ti_target;
}

/* Publish one node's mirror.  Plain copies from the AoS, never monotone clamps: the moment
 * refresh copies several of these fields BACK into the AoS, so a mirror value that could ever
 * exceed its AoS source would be copied back and inflate it. */
static KOKKOS_INLINE_FUNCTION void
gpu_node_mirror_publish(const struct gpu_node_mirror_ptrs_t &m, int k, int no,
                        struct NODE *Nodes_uvm, struct extNODE *Extnodes_uvm)
{
    /* SoA mirror update: only the fields the walk reads.  Vec3 narrowing
     * cast for mixed-precision builds (MyGravFloat=float, MyFloat=double). */
    m.len[k]  = Nodes_uvm[no].len;
    /* The time this length was written at: widen-on-open reads the pair. */
    if(m.node_ti) {m.node_ti[k] = Nodes_uvm[no].Ti_current;}
    m.s[k]    = { (MyGravFloat)Nodes_uvm[no].u.d.s[0],
                    (MyGravFloat)Nodes_uvm[no].u.d.s[1],
                    (MyGravFloat)Nodes_uvm[no].u.d.s[2] };
    m.vs[k]   = { (MyGravFloat)Extnodes_uvm[no].vs[0],
                    (MyGravFloat)Extnodes_uvm[no].vs[1],
                    (MyGravFloat)Extnodes_uvm[no].vs[2] };
    m.hmax[k] = (MyGravFloat)Extnodes_uvm[no].hmax;
    m.bitflags[k] = Nodes_uvm[no].u.d.bitflags;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    m.rt_s[k]  = { (MyGravFloat)Nodes_uvm[no].rt_source_lum_s[0],
                     (MyGravFloat)Nodes_uvm[no].rt_source_lum_s[1],
                     (MyGravFloat)Nodes_uvm[no].rt_source_lum_s[2] };
    m.rt_vs[k] = { (MyGravFloat)Extnodes_uvm[no].rt_source_lum_vs[0],
                     (MyGravFloat)Extnodes_uvm[no].rt_source_lum_vs[1],
                     (MyGravFloat)Extnodes_uvm[no].rt_source_lum_vs[2] };
#endif
#ifdef DM_SCALARFIELD_SCREENING
    m.s_dm[k]  = { (MyGravFloat)Nodes_uvm[no].s_dm[0],
                     (MyGravFloat)Nodes_uvm[no].s_dm[1],
                     (MyGravFloat)Nodes_uvm[no].s_dm[2] };
    m.vs_dm[k] = { (MyGravFloat)Extnodes_uvm[no].vs_dm[0],
                     (MyGravFloat)Extnodes_uvm[no].vs_dm[1],
                     (MyGravFloat)Extnodes_uvm[no].vs_dm[2] };
#endif
#ifdef SINK_NODE_MOTION_TRACKED
    /* the device walk reads sink_pos/sink_vel from the SoA, not the AoS, so the fold and drift
       above are invisible to it unless they are mirrored here with the other drifted moments */
    m.sink_pos[k] = { (MyGravFloat)Nodes_uvm[no].sink_pos[0],
                        (MyGravFloat)Nodes_uvm[no].sink_pos[1],
                        (MyGravFloat)Nodes_uvm[no].sink_pos[2] };
    m.sink_vel[k] = { (MyGravFloat)Nodes_uvm[no].sink_vel[0],
                        (MyGravFloat)Nodes_uvm[no].sink_vel[1],
                        (MyGravFloat)Nodes_uvm[no].sink_vel[2] };
#endif
    if(m.vmax) {m.vmax[k] = (MyGravFloat) Extnodes_uvm[no].vmax;}
}

extern "C" int gpu_force_drift_nodes_ex(integertime time1, int refresh_mirrors_already_current)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    if(Numnodestree <= 0) {return 0;}

    /* This sweep skips any node already at time1 and therefore also skips writing that
     * node's SoA mirror. A host lazy drift advances Nodes[].Ti_current without touching
     * the mirror, so if the host has already drifted to time1 the sweep would leave the
     * walk reading pre-drift geometry. The routing keeps the two apart -- a step whose
     * tree update runs on the host also walks on the host -- and this is the check that
     * the two never disagree. The caller treats a nonzero return as a controlled stop
     * taken by every rank at the next poll, so no rank exits a collective alone. */
    if(force_host_lazy_drift_ti() == time1 && !refresh_mirrors_already_current) {
        printf("gpu_force_drift_nodes: task %d already drifted nodes on the host at this time; the device node mirror cannot be brought up to date by a sweep that skips them\n", ThisTask);
        return 1;
    }

    struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa || !soa->len || !soa->s || !soa->node_vs || !soa->bitflags
            || !soa->hmax || !soa->vmax) {
        printf("gpu_force_drift_nodes: SoA mirrors not ready\n");
        return 1;
    }

    /* Refresh the table mirror unconditionally (cheap; 16 KB) and take the log-time
     * bounds from current All.*, so a restart that moved TimeMax cannot leave the
     * mirror describing the old span. */
    struct DriftKickTableView table_view;
    if(drift_kick_table_mirror_refresh(&drift_kick_table_dev_, &table_view) != 0) {return 1;}   /* soft bad-stop propagated: no kernel launch without the tables */

    int      tree_base        = All.TreeNodeIndexBase;
    int      n_local_nodes  = Numnodestree;
    int      n_foreign_nodes = Numforeignnodes;
    int      maxNodes_snap  = MaxNodes;               /* foreign-range base */
    int      n_nodes        = n_local_nodes + n_foreign_nodes;

#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
    /* Pre-compute per-node dilation factors on the host (cheap O(n_nodes) loop
     * with one All.SpecialParticle_Position_ForRefinement comparison per node
     * for SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM). The GPU kernel then just reads
     * dilation_dev[kk] -- no GPU-side host-only field access required. */
    double *dilation_dev = (double *) gizmo_gpu_alloc_shared((size_t)n_nodes * sizeof(double), NULL);
    if(!dilation_dev) {printf("gpu_force_drift_nodes: dilation_dev alloc failed\n"); endrun(929702); return 1;}   /* soft bad-stop: skip the dilation omp loop on NULL; drains at next poll */
#pragma omp parallel for schedule(static)
    for(int kk = 0; kk < n_nodes; kk++) {
        int no_kk;
        if(kk < n_local_nodes) no_kk = tree_base + kk;
        else                   no_kk = tree_base + maxNodes_snap + (kk - n_local_nodes);
        dilation_dev[kk] = return_node_timestep_dilation_factor(no_kk);
    }
#endif


    integertime ti_target   = time1;
    const int   refresh_all = refresh_mirrors_already_current;   /* captured by value */

    /* Use the SHIFTED pointers (Nodes = Nodes_base - tree_base, Extnodes = Extnodes_base - tree_base)
     * so that Nodes_uvm[tree_base + k] == Nodes_base[k] for all k in [0, Numnodestree). */
    struct NODE     *Nodes_uvm    = Nodes;
    struct extNODE  *Extnodes_uvm = Extnodes;

    /* The mirror pointers, packed: the kernel captures ONE value, and the subset consumer
     * fills the same pack from the same SoA so both publish through one unit. */
    const struct gpu_node_mirror_ptrs_t mirror = gpu_node_mirror_ptrs(soa);

    Kokkos::parallel_for("gpu_force_drift_nodes", n_nodes, KOKKOS_LAMBDA(int kk) {
        /* kk in [0, n_local_nodes) drives local nodes; kk in
         * [n_local_nodes, n_local_nodes + n_foreign_nodes) drives foreign
         * nodes installed by LET unpack (slot index in [0, Numforeignnodes)).
         * SoA index `k` differs from iteration index for foreign nodes:
         *   local:    k_soa = kk          ; no = tree_base + kk
         *   foreign:  k_soa = maxNodes_snap + (kk - n_local_nodes)
         *             no    = tree_base + maxNodes_snap + (kk - n_local_nodes) */
        int k, no;
        if(kk < n_local_nodes) {
            k  = kk;
            no = tree_base + kk;
        } else {
            int slot = kk - n_local_nodes;
            k  = maxNodes_snap + slot;
            no = tree_base + maxNodes_snap + slot;
        }
        /* Already at the target time.  Normally there is nothing to do -- but a
         * host lazy drift advances the node WITHOUT writing its mirror, so the
         * mirror can be stale while the node is current, and skipping here is
         * what leaves it that way.  When the caller asks, fall through with a
         * zero drift: s, len and hmax are all advanced by dt, so a zero dt
         * changes none of them, the kick fold is dt-independent and idempotent
         * on an already-folded node, and the tail rewrites the mirror from the
         * node's current values.  Cost is a full mirror rewrite; there is no
         * other effect. */
        const bool node_already_current = (Nodes_uvm[no].Ti_current == ti_target);
        if(node_already_current && !refresh_all) {return;}

        /* Per-node dilation factor (host-pre-computed, see dispatcher above). */
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
        double dilation = dilation_dev[kk];
#else
        double dilation = 1.0;
#endif

        /* Same value as the host get_drift_factor(.., .., no, 1): one interpolator, one view. */
        double dt_drift = node_already_current
                            ? 0.0
                            : get_drift_factor_impl(Nodes_uvm[no].Ti_current, ti_target,
                                                    dilation, &table_view);
        double dt_drift_hmax = dt_drift;
        /* The widening runs on the undilated clock (see force_drift_node): vmax carries each
           member's own dilation.  One interpolation when the two clocks coincide. */
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
        double dt_widen = node_already_current
                            ? 0.0
                            : get_drift_factor_impl(Nodes_uvm[no].Ti_current, ti_target,
                                                    1.0, &table_view);
#else
        double dt_widen = dt_drift;
#endif

        gpu_node_drift_apply(Nodes_uvm, Extnodes_uvm, no, ti_target,
                             dt_drift, dt_drift_hmax, dt_widen, /*fold_kick=*/1);

        gpu_node_mirror_publish(mirror, k, no, Nodes_uvm, Extnodes_uvm);
    });
    Kokkos::fence();
    gizmo_gpu_check_last_error("gpu_force_drift_nodes", n_nodes);
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(dilation_dev);
#endif
    /* The sweep drifted every node's AoS + SoA geometry to time1 → record
     * certification so a consumer can confirm the tree is current at time1
     * without re-sweeping. */
    gpu_gravity_soa_mark_drift_certified(time1);
    /* ⛔ ONLY the full-refresh variant answers the claims. The ordinary sweep
       (refresh_mirrors_already_current == 0) SKIPS precisely the advanced-but-
       unmirrored nodes the claims describe, so retiring the epoch there would
       discard them without the repair they were recorded for -- and the mirror
       would stay behind with nothing left to say so. */
    if(refresh_mirrors_already_current) {gpu_node_dirty_begin_epoch();}
    return 0;
}

/* Bring just the nodes a recorder listed current at `ti`, instead of sweeping every node.
 *
 * The two claimers are the host lazy drift, which advances a node's canonical state and leaves
 * its mirror behind, and the gravity discovery pre-walk, which finds the nodes that walk will
 * reach while they still stand behind it.  Both mean the same thing here -- this node needs
 * attention before the next device gravity walk at `ti` -- and the node's OWN clock says which
 * kind of attention:
 *
 *   already at `ti`  -> the arithmetic is done; publish the mirror.  It does NOT fold a pending
 *                       kick, which is what the host does on a node it finds already current.
 *   behind `ti`      -> the same per-node drift the full sweep runs, folding the kick as the
 *                       sweep does, and then the same publication.
 *
 * Both branches go through the two shared units above, so the field set cannot drift apart from
 * the sweep's: a field added to the publisher is published here in the same edit.
 *
 * Returns 0 when every listed node stands at `ti` with its mirror published, and 1 when the
 * caller must fall back to the full sweep -- the recorder's fail-safe fired, it overflowed, or
 * an index is outside the mirror.  There is no partial success: a nonzero return means nothing
 * here may be relied on, and the caller sweeps before any walk reads the tree. */
extern "C" int gpu_node_dirty_bring_gravity_current(integertime time1)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* ACQUIRE the claim phase before reading what it published, pairing with the fence inside
     * the claim itself. */
    Kokkos::memory_fence();
    const struct gpu_node_dirty_view_t v = gpu_node_dirty_view();
    if(!v.seen || !v.list || !v.ctl) {return 1;}
    if(Kokkos::atomic_load(&v.ctl->unsafe)) {return 1;}

    struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa || !soa->len || !soa->node_ti || !soa->s || !soa->node_vs
            || !soa->bitflags || !soa->hmax || !soa->vmax) {return 1;}

    const int n = Kokkos::atomic_load(&v.ctl->count);
    if(n < 0 || n > v.cap) {return 1;}   /* an overflowed list names fewer nodes than were claimed */

    int       *list = v.list;
    const int  base = v.base;
    const int  cap  = v.cap;
    for(int i = 0; i < n; i++) {
        const int k = list[i] - base;
        if(k < 0 || k >= cap) {return 1;}
    }

    if(n > 0) {
        /* The same interpolator and the same table view the sweep uses. */
        struct DriftKickTableView table_view;
        if(drift_kick_table_mirror_refresh(&drift_kick_table_dev_, &table_view) != 0) {return 1;}

#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
        /* Per-node dilation, host-computed as the sweep does it -- but for the LISTED nodes
         * only, which is the whole point of this path. */
        double *dilation_dev = (double *) gizmo_gpu_alloc_shared((size_t) n * sizeof(double), NULL);
        if(!dilation_dev) {return 1;}
#pragma omp parallel for schedule(static)
        for(int i = 0; i < n; i++) {dilation_dev[i] = return_node_timestep_dilation_factor(list[i]);}
#endif

        struct NODE    *Nodes_uvm    = Nodes;
        struct extNODE *Extnodes_uvm = Extnodes;
        const struct gpu_node_mirror_ptrs_t mirror = gpu_node_mirror_ptrs(soa);
        const integertime ti_target = time1;

        /* One lane per listed node: independent O(1) work, no serial chain, so the second level
         * of parallelism this loop is asked for is over the list itself. */
        Kokkos::parallel_for("gpu_node_subset_bring_current", n, KOKKOS_LAMBDA(int i) {
            const int no = list[i];
            const int k  = no - base;
            if(Nodes_uvm[no].Ti_current != ti_target) {
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
                const double dilation = dilation_dev[i];
#else
                const double dilation = 1.0;
#endif
                const double dt_drift = get_drift_factor_impl(Nodes_uvm[no].Ti_current, ti_target,
                                                              dilation, &table_view);
                /* the widening runs on the undilated clock, exactly as in the sweep */
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
                const double dt_widen = get_drift_factor_impl(Nodes_uvm[no].Ti_current, ti_target,
                                                              1.0, &table_view);
#else
                const double dt_widen = dt_drift;
#endif
                gpu_node_drift_apply(Nodes_uvm, Extnodes_uvm, no, ti_target,
                                     dt_drift, dt_drift, dt_widen, /*fold_kick=*/1);
            }
            gpu_node_mirror_publish(mirror, k, no, Nodes_uvm, Extnodes_uvm);
        });
        Kokkos::fence();
        gizmo_gpu_check_last_error("gpu_node_subset_bring_current", n);
#ifdef USE_TIMESTEP_DILATION_FOR_ZOOMS
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(dilation_dev);
#endif
    }

    /* The claims are answered, so close the epoch: the generation bump invalidates every stamp
     * in O(1) and hands the recorder back to the host claimers for the rest of the step.  A node
     * drifted again afterwards must be able to re-claim, which is what resetting the cursor
     * alone would not allow.
     * This publishes NO whole-tree certificate: what has been brought current is the listed set,
     * and gpu_gravity_tree_nodes_current_at() keeps its global meaning. */
    gpu_node_dirty_begin_epoch();
    return 0;
}

extern "C" void gpu_force_drift_release(void)
{
    if(drift_kick_table_dev_) {
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(drift_kick_table_dev_);
        drift_kick_table_dev_ = NULL;
    }
}

/* The ordinary sweep: advance whatever is behind, leave the rest alone. */
extern "C" int gpu_force_drift_nodes(integertime time1)
{
    return gpu_force_drift_nodes_ex(time1, 0);
}
