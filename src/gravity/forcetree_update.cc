#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../core/timestep_functions.h"   /* dilation, for the motion bound */
#include "force_node_drift_sync.h"
#include "gpu_gravity_tree.h"        /* SoA mirror: vmax coherence for widen-on-open */
#include "gravtree_moment_kernel.h"   /* the shared node-motion arithmetic */

/* GPU replacement for force_update_tree. */
extern "C" void gpu_force_update_tree(void);

/* Atomic max for doubles using integer CAS (clang doesn't support __atomic on floats) */
static inline void atomic_max_double(double* addr, double val) {
    uint64_t val_bits; memcpy(&val_bits, &val, sizeof(double));
    uint64_t old_bits; memcpy(&old_bits, addr, sizeof(double));
    double old_val; memcpy(&old_val, &old_bits, sizeof(double));
    while(val > old_val) {
        if(__atomic_compare_exchange_n((uint64_t*)addr, &old_bits, val_bits, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {break;}
        memcpy(&old_val, &old_bits, sizeof(double));
    }
}
static_assert(sizeof(double) == sizeof(uint64_t), "double must be 64-bit for atomic CAS");



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


/* Time of the last host lazy node drift on this rank; see force_host_lazy_drift_ti(). */
static integertime host_lazy_drift_ti = -1;

integertime force_host_lazy_drift_ti(void)
{
    return __atomic_load_n(&host_lazy_drift_ti, __ATOMIC_RELAXED);
}


/*! Propagate the momentum kicks accumulated in P[i].dp since the last tree update up
 *  each active element's chain of parent nodes, so the existing tree can be reused for
 *  another step instead of rebuilt. Two paths compute the same thing:
 *
 *  host: walk the chain and drift each parent node at the moment it is touched, so the
 *        cost is set by how many nodes the active set actually reaches.
 *  device: bring current, in parallel, exactly the nodes on those chains (and the top-level
 *        tree), then walk the chains in a kernel; when the update is a large fraction of the
 *        rank (TreeUpdateFullSweep_ActiveFraction) it instead drifts every node in one sweep,
 *        which also certifies the whole tree current for the walks that follow. Bringing the
 *        nodes current first is what makes the parallel walk race-free (force_drift_node is
 *        serial by construction). Each device call still carries fixed launch costs, which is
 *        why the smallest steps take the host path.
 *
 *  This runs before the active list is rebuilt for the coming step, so the list it walks
 *  is the one whose kicks have just been closed out — the correct set for this point in
 *  the KDK sequence, and the quantity that sizes the work here.
 *
 *  The route therefore keys on the COMING step's count (All.TreeUpdateHostBelowActive), not
 *  on how many kicks are being propagated. It is independent of the gravity walk's own
 *  routing: a host update drifts nodes lazily and claims them, and a device gravity walk that
 *  follows answers those claims before it reads the mirror, so the two routes need not agree.
 *  When this says host, the device sweep would only cost drifting every node in the tree to
 *  serve however few elements are active. */
void force_update_tree(void)
{
    /* Two tests.  The count: the device path must drift every node in the tree before its
     * parallel kick walk is race-free, and that sweep is sized by the tree, not by how much work
     * this call has to do.  And once any node has been drifted lazily at this time the update
     * stays on the host whatever the count says: the device path's full-tree sweep skips nodes
     * already at its target time, so it cannot bring their mirror up to date and refuses to run
     * (gpu_force_drift_nodes). A tree built after that drift is exempt: the build rewrote every
     * node and every mirror. */
    const int host_drifted_nodes_now = !gpu_gravity_tree_nodes_current_at(All.Ti_Current)
                                       && (force_host_lazy_drift_ti() == All.Ti_Current);
    const int few_updates = All.TreeUpdateHostBelowActive > 0 && NumForceUpdateAtSyncPoint < All.TreeUpdateHostBelowActive;
    const int to_host = host_drifted_nodes_now || few_updates;

    if(!to_host) {gpu_force_update_tree();}
    else
    {
        PRINT_STATUS("Kick-subroutine will prepare for dynamic update of tree");
        GlobFlag++;
        DomainNumChanged = 0;
        DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));
        for(int i : ActiveParticleList)
        {
            force_kick_node(i, P[i].dp);
            P[i].dp = {};
        }
        force_finish_kick_nodes();
        myfree(DomainList);
        DomainList = NULL;
        PRINT_STATUS(" ..Tree has been updated dynamically");
    }

}


/*! Add element i's momentum kick to every parent node above it, drifting each node as it
 *  is reached. Serial by construction: force_drift_node is not thread-safe. Mirrors the
 *  device kick kernel in gpu_force_update.cc field for field; the difference is that the
 *  device version relies on the preceding all-node drift sweep and therefore uses atomics
 *  where this accumulates directly. */

/* Raise the SoA mirror of a node's vmax to match the AoS.
 *
 * vmax is a RUNNING MAX, and the ONEWAY device walk widens its opening bound by
 * TREE_NODE_WIDENING_DELTA(vmax, dt).  A mirror left behind the AoS is therefore SMALLER, the
 * bound is UNDER-widened, and the walk silently under-includes neighbours -- the
 * one failure this contract exists to prevent.  Raising (never lowering) keeps the
 * mirror conservative even if a writer is missed: too large only over-widens, and
 * the pair kernel re-gates.
 *
 * Indexing follows forcetree.cc:1441: slot k = no - All.TreeNodeIndexBase, valid
 * for local nodes (k < MaxNodes) and installed foreign ones, bounded by the mirror
 * that exists rather than by the index range it sits in. */
static inline int force_soa_raise_vmax(int no, MyFloat vmax_aos)   /* returns the mirror slot, or -1 when there is none */
{
    struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa || !soa->vmax) {return -1;}
    /* ⛔ Bound by the mirror that EXISTS. MaxNodes + AllocatedForeignNodes is the
       INDEX range; the allocation can be smaller (gpu_neighbor_list.cc declines a
       walk precisely when capacity < that sum), and writing past it corrupts the
       neighbouring SoA arrays -- which surfaces as the LET walk resolving
       structure the import does not carry, nowhere near this line. */
    const int k = gpu_gravity_tree_mirror_slot(no);
    if(k < 0) {return -1;}
    const MyGravFloat v = (MyGravFloat) vmax_aos;
    if(soa->vmax[k] < v) {soa->vmax[k] = v;}
    return k;
}

/* The same raise for a node that has just been kicked, which also marks the mirror as holding a
 * pending kick, so a walk predicting the node knows to read the impulse from the canonical node.
 * A raise that is not a kick (gravity_note_motion_bound) leaves the mark alone: there is no
 * impulse to read, and the node's kick time was not stamped. */
static inline void force_soa_mark_kick(int no, MyFloat vmax_aos)
{
    const int k = force_soa_raise_vmax(no, vmax_aos);
    struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(k < 0 || !soa->bitflags) {return;}
    soa->bitflags[k] |= (1u << BITFLAG_NODEHASBEENKICKED);
}

void force_kick_node(int i, Vec3<MyDouble>& dp)
{
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    double lum[N_RT_FREQ_BINS];
    Vec3<MyDouble> rt_source_lum_dp = rt_get_source_luminosity(i, -1, lum, P, CellP) ? dp : Vec3<MyDouble>{};
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Vec3<MyDouble> dp_dm = (P[i].Type != 0) ? dp : Vec3<MyDouble>{};
#endif
#ifdef SINK_NODE_MOTION_TRACKED
    /* Same particle type the moment builders sum sink_mass/sink_pos over, so the momentum and the
       mass it is divided by describe the same set. The type is configurable, hence the macro. */
    Vec3<MyDouble> sink_dp = (P[i].Type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) ? dp : Vec3<MyDouble>{};
#endif

    const MyFloat vmax = (MyFloat) particle_motion_speed_bound(i, P, CellP);

    int no = Father[i];
    while(no >= 0)
    {
        force_drift_node(no, All.Ti_Current);

        Extnodes[no].dp += dp;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Extnodes[no].rt_source_lum_dp += rt_source_lum_dp;
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Extnodes[no].dp_dm += dp_dm;
#endif
#ifdef SINK_NODE_MOTION_TRACKED
        Extnodes[no].sink_dp += sink_dp;
#endif
        if(Extnodes[no].vmax < vmax) {Extnodes[no].vmax = vmax;}
        force_soa_mark_kick(no, Extnodes[no].vmax);   /* keep the walk's mirror conservative */
        Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);
        Extnodes[no].Ti_lastkicked = All.Ti_Current;

        if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))    /* top-level node reached: the rest of the chain is shared with other ranks and is handled by the exchange in force_finish_kick_nodes */
        {
            if(Extnodes[no].Flag != GlobFlag)
            {
                Extnodes[no].Flag = GlobFlag;
                DomainList[DomainNumChanged++] = no;
            }
            break;
        }

        no = Nodes[no].u.d.father;
    }
}








/* The kick a top-level node receives: summed momentum, maximum speed. Neither
 * depends on the order its contributions arrive in, which is what lets the apply
 * below total a node's whole subtree once instead of once per contributing node.
 * The sum is therefore reassociated relative to a per-record walk, so it rounds
 * differently in the last bits while describing the same momentum. TU-local. */
namespace {
struct TopNodeKick
{
  MyDouble dp[3];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
  MyDouble rt_dp[3];
#endif
#ifdef DM_SCALARFIELD_SCREENING
  MyDouble dp_dm[3];
#endif
#ifdef SINK_NODE_MOTION_TRACKED
  MyDouble sink_dp[3];
#endif
  MyFloat vmax;
};

inline void kick_accumulate(struct TopNodeKick& into, const struct TopNodeKick& from)
{
  for(int k = 0; k < 3; k++) {into.dp[k] += from.dp[k];}
#ifdef RT_SEPARATELY_TRACK_LUMPOS
  for(int k = 0; k < 3; k++) {into.rt_dp[k] += from.rt_dp[k];}
#endif
#ifdef DM_SCALARFIELD_SCREENING
  for(int k = 0; k < 3; k++) {into.dp_dm[k] += from.dp_dm[k];}
#endif
#ifdef SINK_NODE_MOTION_TRACKED
  for(int k = 0; k < 3; k++) {into.sink_dp[k] += from.sink_dp[k];}
#endif
  if(into.vmax < from.vmax) {into.vmax = from.vmax;}
}

/* One changed top-level node as it crosses the fused Allgatherv. Fixed layout is
 * identical on every rank, so an MPI_BYTE exchange is field-equivalent (any
 * padding bytes are transmitted but ignored). */
struct DomainKickPacked
{
  int node;
  struct TopNodeKick kick;
};

/* The apply below keeps per-node scratch indexed by a node's offset within the
 * replicated top-level tree, so every node it handles must lie in that range:
 * the records arrive from other ranks, and the chain above each one is followed
 * through father links. Both of the two questions asked about a node -- is it in
 * range, and has this call already seen it -- are asked here, because the mark
 * alone does not bound the offset. */
inline int top_level_node_in_range(int no)
{
  return (no >= All.TreeNodeIndexBase) && (no - All.TreeNodeIndexBase < NTopnodes);
}

inline void report_node_out_of_top_level_range(int no)
{
  printf("Task=%d force_finish_kick_nodes: node %d outside the top-level range [%d,%d)\n",
         ThisTask, no, All.TreeNodeIndexBase, All.TreeNodeIndexBase + NTopnodes);
  fflush(stdout);
  endrun(91562);
}

/* Scratch slot of an already-seen top-level node, or -1 if it is neither. */
inline int marked_top_level_slot(const int *node_slot, int no)
{
  if(!top_level_node_in_range(no)) {return -1;}
  if(Extnodes[no].Flag != GlobFlag) {return -1;}
  return node_slot[no - All.TreeNodeIndexBase];
}
}  /* anonymous namespace */

/* ---- Motion bounds raised outside the kick --------------------------------
 * The pending set is the top-level nodes reached by such raises since the last
 * flush, deduplicated by a per-node stamp.  Both live for the tree's lifetime
 * and are re-sized when the top-level tree changes size. */
static int  *g_pending_topnode_list  = NULL;   /* [NTopnodes] */
static int  *g_pending_topnode_stamp = NULL;   /* [NTopnodes]: == g_pending_stamp when listed */
static int   g_pending_topnode_n     = 0;
static int   g_pending_topnode_cap   = 0;
static int   g_pending_stamp         = 1;

static void pending_topnodes_ensure(void)
{
    if(g_pending_topnode_cap == NTopnodes && g_pending_topnode_list) {return;}
    free(g_pending_topnode_list); free(g_pending_topnode_stamp);
    g_pending_topnode_cap   = NTopnodes;
    g_pending_topnode_list  = (int *) calloc((size_t)(NTopnodes > 0 ? NTopnodes : 1), sizeof(int));
    g_pending_topnode_stamp = (int *) calloc((size_t)(NTopnodes > 0 ? NTopnodes : 1), sizeof(int));
    g_pending_topnode_n = 0;
    g_pending_stamp = 1;
}

void gravity_clear_pending_motion_bounds(void)
{
    g_pending_topnode_n = 0;
    g_pending_stamp++;
    if(g_pending_stamp == 0) {g_pending_stamp = 1; if(g_pending_topnode_stamp) {memset(g_pending_topnode_stamp, 0, (size_t)g_pending_topnode_cap * sizeof(int));}}
}

/* Raise one node's bound, in the node and in the walk's mirror. */
static inline void raise_node_motion_bound(int no, MyFloat vmax)
{
    if(Extnodes[no].vmax < vmax) {Extnodes[no].vmax = vmax;}
    (void) force_soa_raise_vmax(no, Extnodes[no].vmax);
}

void gravity_note_motion_bound(const int *idx, int n)
{
    if(n <= 0 || !idx || !Father || !Nodes || !Extnodes) {return;}
    pending_topnodes_ensure();
    for(int k = 0; k < n; k++)
    {
        const int i = idx[k];
        /* Father[] covers the slots the standing tree was built for; a particle
         * created since then has no node above it until the next build. */
        if(i < 0 || i >= NumPart || i >= All.TreeParticleSlots) {continue;}
        const MyFloat vmax = (MyFloat) particle_motion_speed_bound(i, P, CellP);
        int no = Father[i];
        while(no >= 0)
        {
            /* A node already holding the bound has ancestors that do too. */
            if(Extnodes[no].vmax >= vmax) {break;}
            raise_node_motion_bound(no, vmax);
            if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
            {
                /* Shared with every rank: listed once for the flush.  The rest of
                 * the chain is raised here too, so this rank's own walks see it
                 * before the exchange. */
                const int t = no - All.TreeNodeIndexBase;
                if(t >= 0 && t < g_pending_topnode_cap && g_pending_topnode_stamp[t] != g_pending_stamp)
                {
                    g_pending_topnode_stamp[t] = g_pending_stamp;
                    g_pending_topnode_list[g_pending_topnode_n++] = no;
                }
            }
            no = Nodes[no].u.d.father;
        }
    }
}

void gravity_flush_pending_motion_bounds(void)
{
    /* Every rank enters, with or without pending nodes of its own. */
    pending_topnodes_ensure();
    int *counts = (int *) mymalloc("mb_counts", sizeof(int) * NTask);
    int n_local = g_pending_topnode_n;
    MPI_Allgather(&n_local, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
    int total = 0;
    for(int ta = 0; ta < NTask; ta++) {total += counts[ta];}
    if(total > 0)
    {
        struct NodeBound {int node; MyFloat vmax;};
        int *counts_b = (int *) mymalloc("mb_counts_b", sizeof(int) * NTask);
        int *offset_b = (int *) mymalloc("mb_offset_b", sizeof(int) * NTask);
        struct NodeBound *loc = (struct NodeBound *) mymalloc("mb_loc", (size_t)(n_local > 0 ? n_local : 1) * sizeof(struct NodeBound));
        struct NodeBound *all = (struct NodeBound *) mymalloc("mb_all", (size_t)total * sizeof(struct NodeBound));
        for(int k = 0; k < n_local; k++) {loc[k].node = g_pending_topnode_list[k]; loc[k].vmax = Extnodes[g_pending_topnode_list[k]].vmax;}
        for(int ta = 0; ta < NTask; ta++)
        {
            counts_b[ta] = counts[ta] * (int) sizeof(struct NodeBound);
            offset_b[ta] = (ta == 0) ? 0 : offset_b[ta - 1] + counts_b[ta - 1];
        }
        MPI_Allgatherv(loc, n_local * (int) sizeof(struct NodeBound), MPI_BYTE, all, counts_b, offset_b, MPI_BYTE, MPI_COMM_WORLD);
        for(int r = 0; r < total; r++)
        {
            int no = all[r].node;
            /* The record names a top-level node, which every rank holds at the
             * same index; validated rather than trusted, as the kick exchange does. */
            if(!(no >= All.TreeNodeIndexBase && no - All.TreeNodeIndexBase < NTopnodes)) {continue;}
            const MyFloat vmax = all[r].vmax;
            while(no >= 0)
            {
                if(Extnodes[no].vmax >= vmax) {break;}
                raise_node_motion_bound(no, vmax);
                no = Nodes[no].u.d.father;
            }
        }
        myfree(all); myfree(loc); myfree(offset_b); myfree(counts_b);
    }
    myfree(counts);
    gravity_clear_pending_motion_bounds();
}

void force_finish_kick_nodes(void)
{
  int i, no, ta, totDomainNumChanged;
  int *counts, *counts_dp, *offset_dp;

  /* share the momentum-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  counts_dp = (int *) mymalloc("counts_dp", sizeof(int) * NTask);
  offset_dp = (int *) mymalloc("offset_dp", sizeof(int) * NTask);

  /* Exchange per-rank changed-node counts FIRST (all ranks participate in this
   * collective), so a globally-empty update can early-out BEFORE any local
   * payload malloc/pack, the packed Allgatherv, and the apply loop. When
   * totDomainNumChanged==0 every rank's payload is empty, so all that work is
   * mathematically a no-op; the early-out is exact, and symmetric across ranks
   * (totDomainNumChanged is the same global sum on every rank). */
  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0; ta < NTask; ta++)
    totDomainNumChanged += counts[ta];


  if(totDomainNumChanged == 0)
    {   /* no rank changed a top-level node: nothing to pack/exchange/apply */
      myfree(offset_dp);
      myfree(counts_dp);
      myfree(counts);
      return;
    }

  /* Packed single-Allgatherv path: pack each changed node's fields into one
   * contiguous record and exchange in ONE collective instead of 2-5 separate
   * field Allgathervs (fewer collective latencies / skew-absorption points).
   * Field-equivalent apply. */
  {
    struct DomainKickPacked *rec_loc = (struct DomainKickPacked *)
        mymalloc("fut_rec_loc", DomainNumChanged * sizeof(struct DomainKickPacked));
    for(i = 0; i < DomainNumChanged; i++)
      {
        no = DomainList[i];
        rec_loc[i].node = no;
        rec_loc[i].kick.dp[0] = Extnodes[no].dp[0];
        rec_loc[i].kick.dp[1] = Extnodes[no].dp[1];
        rec_loc[i].kick.dp[2] = Extnodes[no].dp[2];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        rec_loc[i].kick.rt_dp[0] = Extnodes[no].rt_source_lum_dp[0];
        rec_loc[i].kick.rt_dp[1] = Extnodes[no].rt_source_lum_dp[1];
        rec_loc[i].kick.rt_dp[2] = Extnodes[no].rt_source_lum_dp[2];
#endif
#ifdef DM_SCALARFIELD_SCREENING
        rec_loc[i].kick.dp_dm[0] = Extnodes[no].dp_dm[0];
        rec_loc[i].kick.dp_dm[1] = Extnodes[no].dp_dm[1];
        rec_loc[i].kick.dp_dm[2] = Extnodes[no].dp_dm[2];
#endif
#ifdef SINK_NODE_MOTION_TRACKED
        for(int k = 0; k < 3; k++) {rec_loc[i].kick.sink_dp[k] = Extnodes[no].sink_dp[k];}
#endif
        rec_loc[i].kick.vmax = Extnodes[no].vmax;
      }
    /* byte counts/offsets for the fixed-size records (reuse counts_dp/offset_dp;
     * counts[] is still the raw per-rank node count here) */
    for(ta = 0; ta < NTask; ta++)
      {
        counts_dp[ta] = counts[ta] * (int) sizeof(struct DomainKickPacked);
        offset_dp[ta] = (ta == 0) ? 0 : offset_dp[ta - 1] + counts[ta - 1] * (int) sizeof(struct DomainKickPacked);
      }
    PRINT_STATUS(" ..exchanged kick momenta for %d top-level nodes out of %d", totDomainNumChanged, NTopleaves);
    struct DomainKickPacked *rec_all = (struct DomainKickPacked *)
        mymalloc("fut_rec_all", totDomainNumChanged * sizeof(struct DomainKickPacked));
    MPI_Allgatherv(rec_loc, DomainNumChanged * (int) sizeof(struct DomainKickPacked), MPI_BYTE,
                   rec_all, counts_dp, offset_dp, MPI_BYTE, MPI_COMM_WORLD);
    /* Apply every rank's records to this rank's copy of the top-level tree. A
     * node's total is its own records plus everything its children received, so
     * each record is added once at its own node and the sums are carried upward
     * in a single sweep -- rather than each record walking its whole chain to
     * the root, which repeats the shared upper part of the chain once per
     * record. Every node reached here is top-level (force_flag_localnodes marks
     * the entire ancestor chain of each top leaf), and the top-level tree is
     * built parent-first from the tree's node index base, so descending node
     * index orders children before their parents. */
    int *uniq = (int *) mymalloc("fut_uniq", (NTopnodes > 0 ? NTopnodes : 1) * sizeof(int));
    /* Slot of a node in uniq[], indexed by the node's offset within the
     * top-level tree. Entries for nodes this call did not reach are never read
     * -- marked_top_level_slot() is the only reader and it range-checks first --
     * so the array needs no initialization. */
    int *node_slot = (int *) mymalloc("fut_node_slot", (NTopnodes > 0 ? NTopnodes : 1) * sizeof(int));
    int nuniq = 0;

    GlobFlag++;
    for(i = 0; i < totDomainNumChanged; i++)
      {
        /* every node is range-checked BEFORE it is used to read anything: the
         * record comes from another rank, and the chain above it is followed
         * through father links */
        no = rec_all[i].node;
        if(!top_level_node_in_range(no)) {report_node_out_of_top_level_range(no); rec_all[i].node = -1; continue;}
        if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_ELEMENT))
          no = Nodes[no].u.d.father;   /* already applied to this node by the local kick */
        rec_all[i].node = no;          /* resolved once here; the scatter below reuses it */
        while(no >= 0)
          {
            if(!top_level_node_in_range(no))
              {   /* the chain left the top-level tree: stop walking it, and stop the run */
                report_node_out_of_top_level_range(no);
                break;
              }
            if(Extnodes[no].Flag == GlobFlag) {break;}   /* this call already took the rest of the chain */
            Extnodes[no].Flag = GlobFlag;
            uniq[nuniq++] = no;
            no = Nodes[no].u.d.father;
          }
      }

    std::sort(uniq, uniq + nuniq, [](int a, int b) {return a > b;});

    struct TopNodeKick *acc = (struct TopNodeKick *)
        mymalloc("fut_acc", (nuniq > 0 ? nuniq : 1) * sizeof(struct TopNodeKick));
    for(int k = 0; k < nuniq; k++) {node_slot[uniq[k] - All.TreeNodeIndexBase] = k; acc[k] = {};}

    for(i = 0; i < totDomainNumChanged; i++)
      {
        const int slot = marked_top_level_slot(node_slot, rec_all[i].node);
        if(slot >= 0) {kick_accumulate(acc[slot], rec_all[i].kick);}
      }

    for(int k = 0; k < nuniq; k++)
      {
        no = uniq[k];
        force_drift_node(no, All.Ti_Current);
        Extnodes[no].dp[0] += acc[k].dp[0];
        Extnodes[no].dp[1] += acc[k].dp[1];
        Extnodes[no].dp[2] += acc[k].dp[2];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Extnodes[no].rt_source_lum_dp[0] += acc[k].rt_dp[0];
        Extnodes[no].rt_source_lum_dp[1] += acc[k].rt_dp[1];
        Extnodes[no].rt_source_lum_dp[2] += acc[k].rt_dp[2];
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Extnodes[no].dp_dm[0] += acc[k].dp_dm[0];
        Extnodes[no].dp_dm[1] += acc[k].dp_dm[1];
        Extnodes[no].dp_dm[2] += acc[k].dp_dm[2];
#endif
#ifdef SINK_NODE_MOTION_TRACKED
        for(int d = 0; d < 3; d++) {Extnodes[no].sink_dp[d] += acc[k].sink_dp[d];}
#endif
        if(Extnodes[no].vmax < acc[k].vmax)
          Extnodes[no].vmax = acc[k].vmax;
        /* The MERGED cross-rank value, not this rank's contribution: this site runs on
           BOTH kick routes (gpu_force_update.cc:279 calls it too), so it is where the
           top-level set gets its final answer. */
        force_soa_mark_kick(no, Extnodes[no].vmax);
        Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);
        Extnodes[no].Ti_lastkicked = All.Ti_Current;

        const int father_slot = marked_top_level_slot(node_slot, Nodes[no].u.d.father);
        if(father_slot >= 0) {kick_accumulate(acc[father_slot], acc[k]);}
      }

    myfree(acc);
    myfree(node_slot);
    myfree(uniq);
    myfree(rec_all);
    myfree(rec_loc);
  }
  myfree(offset_dp);
  myfree(counts_dp);
  myfree(counts);
}



void force_drift_node(int no, integertime time1)
{
  /* Acquire-load: if another thread already drifted this node to time1, we both
   * skip AND observe its published geometry (paired with the release store below). */
  if(time1 == modeb_node_ti_current_acquire(no))
    return;

  /* A kicked node was brought current when its kick was added, so its pending momentum
     describes motion from exactly its own time onward. */
  const int kicked = (Nodes[no].u.d.bitflags & (1 << BITFLAG_NODEHASBEENKICKED)) ? 1 : 0;
  if(kicked && Extnodes[no].Ti_lastkicked != Nodes[no].Ti_current)
    {
      printf("Task=%d Extnodes[no].Ti_lastkicked=%lld  Nodes[no].Ti_current=%lld\n",ThisTask, (long long)Extnodes[no].Ti_lastkicked, (long long)Nodes[no].Ti_current);
      printf("inconsistency in drift node\n"); fflush(stdout); endrun(90001007); return;   /* graceful: skip node drift; bad-stop drains at the next gravity-walk poll */
    }

  /* The arithmetic is the shared node-motion unit (gravtree_moment_kernel.h), the same one the
     device sweep runs. */
  const node_motion_in_arrays node = {Nodes, Extnodes, no};
  if(kicked)
    {
      node_motion_fold_kick(node);
      Nodes[no].u.d.bitflags &= (~(1 << BITFLAG_NODEHASBEENKICKED));
    }

    const double dt_drift = get_drift_factor(Nodes[no].Ti_current, time1, no, 1);
    /* The widening and the kernel-length growth run on the undilated clock: vmax bounds each
       member's motion per unit undilated interval, carrying that member's own dilation, and each
       member's kernel length grows over its own dilated interval, so the node's dilated clock
       (right for its centre of mass) would under-grow both for a member less dilated than the
       node.  The same interval when no dilation is active. */
    const double dt_widen = get_drift_factor_undilated(Nodes[no].Ti_current, time1);

    node_motion_advance(node, dt_drift, dt_widen);
    node_hmax_drift(Extnodes[no], dt_widen);

    /* Record that this rank has now drifted at least one node to time1 without
     * updating that node's device SoA mirror. Relaxed: every caller passes
     * All.Ti_Current, so concurrent writers store the same value, and the only reader
     * is the device sweep's invariant check, which runs between phases. */
    __atomic_store_n(&host_lazy_drift_ti, time1, __ATOMIC_RELAXED);

    /* Release store: publishes Ti_current after all geometry/Extnodes writes so a
     * threaded walk's acquire-load fast path sees fresh Ti => fresh geometry. */
    force_drift_node_publish_current(no, time1);

    /* This call ADVANCED this node without writing its device mirror -- reaching
     * here means exactly that, because the acquire-load at the top returns early
     * for a node already at time1. Record it so the mirror can be repaired for the
     * O(Ndirty) set instead of by sweeping the whole tree.
     *
     * Claimed AFTER the release publish above, so the geometry this claim refers
     * to is already published when a repair observes it. One site, not seven: all
     * seven callers funnel through here, and a future eighth is covered without
     * anyone remembering. */
    gpu_node_dirty_claim(no);

}





/*! This function updates the hmax-values in tree nodes that hold gas cells. These values are needed to find all neighbors in the
 *  hydro-force computation.  Since the KernelRadius-values are potentially changed in the fluid-density computation, force_update_hmax() should be carried
 *  out just before the hydrodynamical forces are computed, i.e. after density(). */
void force_update_hmax(void)
{
  int i, no, ta, totDomainNumChanged;
  int *domainList_all;
  int *counts, *offset_list, *offset_hmax;
  MyFloat *domainHmax_loc, *domainHmax_all;
  /* Per-changed-topleaf exchange record: scalar hmax + divVmax + the 6 Mode-B
   * per-type bands.  The per-type slots ride this SAME post-density exchange so
   * remote topleaf/ancestor per-type bands are as fresh as the scalar hmax (they
   * were locally grown above but, without this, were only cross-rank-fresh at the
   * last full tree build/refresh).  Required for the Mode-B SYMMETRIC targeted
   * export band, which prunes remote topleaves by these per-type bands. */
  enum { HMAX_EXCH_HMAX = 0, HMAX_EXCH_DIVVMAX = 1, HMAX_EXCH_PTYPE0 = 2, HMAX_EXCH_SIZE = 8 };
  int OffsetSIZE = HMAX_EXCH_SIZE;
  double divVel;

  GlobFlag++;

  DomainNumChanged = 0;
  DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));

  /* Phase 1: drift all ancestor nodes (serial — force_drift_node is not thread-safe).
   * Mode B per-type bands now cover every type's leaf-policy-selectable radius
   * (gas via KernelRadius+ForceSoftening; non-gas adds AGS_KernelRadius when defined),
   * so non-gas particles must reach this update path in every build, not only
   * ADAPTIVE_GRAVSOFT_FORALL ones. Scalar Extnodes[no].hmax retains its legacy
   * semantics; it is only grown for the AGS-tracked path below. */
  for (int i : ActiveParticleList)
  {
    if(P[i].Mass > 0)
      {
        no = Father[i];
        while(no >= 0)
        {
            if(Nodes[no].Ti_current == All.Ti_Current) {break;}
            force_drift_node(no, All.Ti_Current);
            no = Nodes[no].u.d.father;
        }
      }
  }
  /* Phase 2: update hmax/divVmax/per-type bands with atomics (parallel). */
#pragma omp parallel for schedule(dynamic)
  for (int idx = 0; idx < (int)ActiveParticleList.size(); idx++)
  {
    int i = ActiveParticleList[idx];
    if(P[i].Mass > 0)
      {
        int no = Father[i];
        double divVel = P[i].Particle_DivVel;

        /* Mode B per-type band: conservative across every leaf-policy-selectable
         * source. Helper covers KernelRadius / ForceSoftening / AGS_KernelRadius
         * (when defined) uniformly per type. */
        int per_type_band = (int)P[i].Type;
        double per_type_htmp = force_hmax_per_type_particle_radius(i);

        /* Scalar `hmax`/`divVmax` eligibility (legacy semantics):
         * non-AGS-FORALL builds → gas only; AGS-FORALL builds → any Mass>0 type.
         * Non-eligible particles still update per-type bands above, but MUST NOT
         * leak their divVel / KernelRadius into the scalar band (which feeds
         * downstream cross-rank exchange and legacy walkers). */
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
        const int scalar_eligible = 1;
#else
        const int scalar_eligible = (P[i].Type == 0);
#endif

        while(no >= 0)
        {
            /* Scalar `hmax` source: legacy gas-KR or AGS-KR per ADAPTIVE_GRAVSOFT_FORALL. */
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
            double kernrad_temp = P[i].AGS_KernelRadius;
            if(P[i].Type == 0) {kernrad_temp = P[i].KernelRadius;}
            double htmp = DMIN(kernrad_temp, All.MaxKernelRadius);
#else
            double htmp = (P[i].Type == 0) ? DMIN(P[i].KernelRadius, All.MaxKernelRadius) : 0.0;
#endif
            int per_type_grew = 0;
            if(per_type_htmp > Extnodes[no].hmax_per_type[per_type_band]) {
                atomic_max_double(&Extnodes[no].hmax_per_type[per_type_band], per_type_htmp);
                per_type_grew = 1;
            }
            int scalar_grew = 0;
            if(scalar_eligible && (htmp > Extnodes[no].hmax || divVel > Extnodes[no].divVmax))
            {
                atomic_max_double(&Extnodes[no].hmax, htmp);
                atomic_max_double(&Extnodes[no].divVmax, divVel);
                scalar_grew = 1;
            }
            if(scalar_grew || per_type_grew)
            {
                if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
                {
                    #pragma omp critical(DomainListAppendHmax)
                    {
                        if(Extnodes[no].Flag != GlobFlag)
                        {
                            Extnodes[no].Flag = GlobFlag;
                            DomainList[DomainNumChanged++] = no;
                        }
                    }
                    break;
                }
            }
            else
                break;

            no = Nodes[no].u.d.father;
        }
      }
  }

  /* share the hmax-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
  offset_hmax = (int *) mymalloc("offset_hmax", sizeof(int) * NTask);

  domainHmax_loc = (MyFloat *) mymalloc("domainHmax_loc", DomainNumChanged * OffsetSIZE * sizeof(MyFloat));

  for(i = 0; i < DomainNumChanged; i++)
    {
      domainHmax_loc[OffsetSIZE * i + HMAX_EXCH_HMAX]    = Extnodes[DomainList[i]].hmax;
      domainHmax_loc[OffsetSIZE * i + HMAX_EXCH_DIVVMAX] = Extnodes[DomainList[i]].divVmax;
      for(int t = 0; t < 6; t++)
          domainHmax_loc[OffsetSIZE * i + HMAX_EXCH_PTYPE0 + t] = Extnodes[DomainList[i]].hmax_per_type[t];
    }


  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_hmax[0] = 0; ta < NTask; ta++)
    {
      totDomainNumChanged += counts[ta];
      if(ta > 0)
	{
	  offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
	  offset_hmax[ta] = offset_hmax[ta - 1] + counts[ta - 1] * OffsetSIZE * sizeof(MyFloat);
	}
    }

  PRINT_STATUS(" ..Hmax exchange: %d topleaves out of %d", totDomainNumChanged, NTopleaves);
  domainHmax_all = (MyFloat *) mymalloc("domainHmax_all", totDomainNumChanged * OffsetSIZE * sizeof(MyFloat));
  domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

  MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
		 domainList_all, counts, offset_list, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    {counts[ta] *= OffsetSIZE * sizeof(MyFloat);}

  MPI_Allgatherv(domainHmax_loc, OffsetSIZE * DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
		 domainHmax_all, counts, offset_hmax, MPI_BYTE, MPI_COMM_WORLD);


  for(i = 0; i < totDomainNumChanged; i++)
    {
        no = domainList_all[i];
        if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_ELEMENT))    {no = Nodes[no].u.d.father;} /* to avoid that the hmax is updated twice */
        
        while(no >= 0)
        {
            force_drift_node(no, All.Ti_Current);

            /* Grow scalar hmax/divVmax + every per-type band; keep walking ancestors
             * while ANY of the 8 fields grew (a remote update that only grows one
             * per-type band must still propagate up the tree, just like scalar hmax). */
            int any_grew = 0;
            if(domainHmax_all[OffsetSIZE * i + HMAX_EXCH_HMAX] > Extnodes[no].hmax)
                {Extnodes[no].hmax = domainHmax_all[OffsetSIZE * i + HMAX_EXCH_HMAX]; any_grew = 1;}
            if(domainHmax_all[OffsetSIZE * i + HMAX_EXCH_DIVVMAX] > Extnodes[no].divVmax)
                {Extnodes[no].divVmax = domainHmax_all[OffsetSIZE * i + HMAX_EXCH_DIVVMAX]; any_grew = 1;}
            for(int t = 0; t < 6; t++)
            {
                MyFloat v = domainHmax_all[OffsetSIZE * i + HMAX_EXCH_PTYPE0 + t];
                if(v > Extnodes[no].hmax_per_type[t]) {Extnodes[no].hmax_per_type[t] = v; any_grew = 1;}
            }
            if(!any_grew) {break;}

            no = Nodes[no].u.d.father;
        }
    }


  myfree(domainList_all);
  myfree(domainHmax_all);
  myfree(domainHmax_loc);
  myfree(offset_hmax);
  myfree(offset_list);
  myfree(counts);
  myfree(DomainList);

  force_bump_hmax_refresh_generation();   /* ancestor boxes re-drifted + per-type bands re-seeded */
  CPU_Step[CPU_TREEHMAXUPDATE] += measure_time();
}
