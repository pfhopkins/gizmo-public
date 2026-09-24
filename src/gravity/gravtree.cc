#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"
#include "../core/proto.h"
#include "gpu_gravtree.h"
#include "gpu_gravity_tree.h"   /* gpu_gravity_tree_mark_born_current */
#include "../mesh/gpu_neighbor_list.h"   /* the touched set's lifecycle counters, reported below */
#include "../system/gpu_particles_arena.h"
#include "../core/timestep_functions.h"   /* Hermite pass state, refreshed below */
#include "../mesh/kernel.h"
#include "./analytic_gravity.h"

/*! Host-vs-device routing for the gravity walk and the dynamic tree update, keyed on the
 *  RANK-LOCAL count of active gravity candidates. The device path must drift every node
 *  in the tree before its parallel walk can be race-free, so its floor is set by the tree
 *  size rather than by the active set; the host walk drifts each node only when it opens
 *  it. Below the threshold the sweep costs more than the walk it enables.
 *
 *  The threshold is conservative against a crossover measured near 6e4 rank-local
 *  candidates on 16-rank FIRE, where routing the whole tree walk to the host cut the
 *  cost of steps with fewer than 1e4 global active elements by a third. Above it the
 *  device path wins and the host walk's serial node drift becomes the bottleneck. */
int gravity_walk_route_to_host(long long n_local_active)
{
    /* Once any node has been drifted lazily at this time, the host owns the rest of the
     * time step: the device sweep skips nodes already at its target time, so it can no
     * longer bring their mirror up to date, and a second gravity evaluation at the same
     * time (a Hermite correction pass, a repeated walk for the opening criterion) would
     * otherwise read that stale geometry.  A tree built after that drift is
     * exempt: the build rewrote every node and every mirror, so the record of an
     * earlier lazy drift no longer describes anything. */
    if(!gpu_gravity_tree_nodes_current_at(All.Ti_Current)
            && force_host_lazy_drift_ti() == All.Ti_Current) {return 1;}

    return (All.GravityHostWalkBelowActive > 0 && n_local_active < (long long)All.GravityHostWalkBelowActive) ? 1 : 0;
}

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local elements, and elements are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */

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

double Ewaldcount, Costtotal;
long long N_nodesinlist;
int Ewald_iter;			/* global in file scope, for simplicity */

/* Per-thread workspace of the host walk loop, allocated once per walk pass: the batch of
 * active-list entries a thread takes at a time, the interaction counts it gets back, and -- only
 * when the primary host walk has candidates -- the walk's own packet workspace. Sized from how
 * many candidates the device pre-pass left to the host, never from the configured packet size
 * alone, so a rank with one host target pays for one and a rank with none pays for the batch. */
#ifndef GRAVITY_PRIMARY_LOOP_BATCH_SIZE
#define GRAVITY_PRIMARY_LOOP_BATCH_SIZE 8
#endif
static char *GravWalkWorkspace = NULL;
static int GravWalkPacketCap = 1, GravWalkBatchCap = GRAVITY_PRIMARY_LOOP_BATCH_SIZE;
static size_t GravWalkThreadBytes = 0, GravWalkBatchBytes = 0;
static void gravity_walk_workspace_allocate(int host_candidates)
{
    GravWalkPacketCap = (host_candidates < TREE_QUERY_PACKET_SIZE) ? host_candidates : TREE_QUERY_PACKET_SIZE;
    GravWalkBatchCap = (GravWalkPacketCap > GRAVITY_PRIMARY_LOOP_BATCH_SIZE) ? GravWalkPacketCap : GRAVITY_PRIMARY_LOOP_BATCH_SIZE;
    GravWalkBatchBytes = ((3 * GravWalkBatchCap * sizeof(int)) + 63) & ~((size_t) 63);   /* indices, list positions, interaction counts */
    GravWalkThreadBytes = GravWalkBatchBytes + ((GravWalkPacketCap > 0) ? force_treewalk_workspace_bytes_per_thread(GravWalkPacketCap) : 0);
    GravWalkWorkspace = (char *) mymalloc("GravWalkWorkspace", (size_t) maxThreads * GravWalkThreadBytes);
}
static void gravity_walk_workspace_free(void) {myfree(GravWalkWorkspace); GravWalkWorkspace = NULL;}


/*! This function computes the gravitational forces for all active elements. If needed, a new tree is constructed, otherwise the dynamically updated
 *  tree is used.  Elements are only exported to other processors when needed. */
/*! Promote the staged tree-opening acceleration scale into the value the criterion reads.
 *
 *  Runs ONCE per tree build, immediately before the tree (and with it the imported ghost tree) is
 *  constructed, and nowhere else.  The import is pruned against this value: a sender ships a node
 *  as a childless multipole precisely when no target on the receiving rank would open it, and
 *  OldAcc enters that test.  Were it to change while the tree still stands, a later walk could
 *  apply a stricter criterion than the sender did, ask to descend a node whose children were never
 *  sent, and either stop on the completeness guard or -- under a topleaf that has a sibling --
 *  silently skip that node's mass.  Promoting here ties the value to the lifetime of the structure
 *  pruned against it, so every walk on a given tree opens exactly what its import covers.
 *
 *  All local particles, not just the active ones: a particle may become active while this tree
 *  still stands, and its opening decisions must be covered by the same import.  For a particle
 *  whose walk measured nothing new this rewrites the identical value.
 */
#ifdef ADAPTIVE_TREEFORCE_UPDATE
/*! The one needs_new_treeforce() answer for this gravity_tree() call, indexed by position in
 *  ActiveParticleList.  Frozen before any walk runs, because the walk writes the fields that
 *  predicate reads; see the call site for what disagrees if it is asked twice.  Sized to the
 *  active set, not to NumPart, so a step with few actives pays for few actives. */
static std::vector<unsigned char> TreeforceCandidateFrozen;

void gravity_freeze_treeforce_candidates(void)
{
    TreeforceCandidateFrozen.assign(ActiveParticleList.size(), 0);
    for(int ii = 0; ii < (int)ActiveParticleList.size(); ii++)
    {
        if(needs_new_treeforce(ActiveParticleList[ii])) {TreeforceCandidateFrozen[ii] = 1;}
    }
}

int gravity_treeforce_candidate_frozen(int ii)
{
    if(ii < 0 || ii >= (int)TreeforceCandidateFrozen.size()) {return 1;}  /* not a member of this call's active set */
    return TreeforceCandidateFrozen[ii];
}
#endif

void refresh_old_acceleration_for_tree_opening(void)
{
    /* OldAcc is overloaded on a 2lpt start: the IC reader leaves the particle masses in it, and
     * init.cc deliberately does not zero it for that reason.  Nothing has been staged yet either,
     * so promoting here would replace those masses with zeros before the first evaluation is done
     * with them.  The same condition suppresses the staging site after the walk; the walk that
     * follows this build does not read OldAcc anyway, because the criterion stays Barnes-Hut for
     * exactly this case. */
    if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0) {return;}
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int i = 0; i < NumPart; i++) {P[i].OldAcc = P[i].OldAcc_LatestWalk;}
}

void gravity_tree(void)
{
    /* initialize variables */
    long long n_exported = 0; int i, j, maxnumnodes, iter; i = 0; j = 0; iter = 0; maxnumnodes=0;
    double t0, t1, timeall = 0, timetree1 = 0, timetree2 = 0, timetree, timewait, timecomm;
    double timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0, sum_costtotal, ewaldtot;
    double maxt, sumt, maxt1, sumt1, maxt2, sumt2, sumcommall, sumwaitall, plb, plb_max;
    CPU_Step[CPU_MISC] += measure_time();

    /* set new softening lengths */
    if(All.ComovingIntegrationOn) {set_softenings();}

    /* Refresh the per-particle ForceSoftening cache for active particles.
     * Single source of truth shared by CPU walk, GPU walk, and tree-build
     * split-scale computation in force_treebuild().  Inputs (KernelRadius,
     * AGS_KernelRadius, tidal_tensor_mag_prev, StarParticleEffectiveSize)
     * are guaranteed valid here: hydro/AGS density already ran earlier in
     * the timestep, set_softenings() above just updated All.ForceSoftening[].
     * Inactive particles retain their cached value from when they were last
     * active -- inputs only mutate during active processing, so the cached
     * value is still correct. */
    compute_all_force_softening(0);

#ifdef HERMITE_INTEGRATION
    /* The pass state both walks read.  Refreshed once here, not per target: the active-bin mask
       is a property of the pass, and rebuilding it per target would put a fixed 60-iteration scan
       on every gravity interaction target.  The drift/kick view is only assembled on a Hermite
       pass, since the predictor returns immediately when HermiteOnlyFlag is 0. */
    HermiteWalk = hermite_walk_state_snapshot();
    if(HermiteOnlyFlag) {HermiteWalkTables = drift_kick_table_view_host();}
#endif

    /* construct tree if needed */
#ifdef HERMITE_INTEGRATION
    if(!HermiteOnlyFlag)
#endif
    if(TreeReconstructFlag)
    {
        PRINT_STATUS("Tree construction initiated (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
        CPU_Step[CPU_MISC] += measure_time();
        move_particles(All.Ti_Current);
        rearrange_particle_sequence();
        refresh_old_acceleration_for_tree_opening();
        gizmo_exit_bad_stop_if_requested("gravtree:before_treebuild"); CPU_Step[CPU_DRIFT] += measure_time(); /* sync before we do the treebuild */
        if(force_treebuild(NumPart, NULL) == FORCE_TREE_NEEDS_OWNERSHIP_RESTORE)
        {
            /* Particles have drifted into top-leaves other ranks own, and the standing tree cannot
             * say where they were attached -- either because there is none, or because it could not
             * place them -- so the decomposition in place cannot carry this build.  A repartition
             * hands each particle back to the rank whose top-leaf it now sits in, after which the
             * ordinary build has nothing to retain.  The refinement pass is suppressed: this is a
             * recovery, not the step's merge/split, and running it here would refine twice. */
            if(ThisTask == 0) {printf("Tree build: restoring geometric particle ownership before rebuilding.\n"); fflush(stdout);}
            domain_Decomposition_light(0, 0);
            if(force_treebuild(NumPart, NULL) == FORCE_TREE_NEEDS_OWNERSHIP_RESTORE) {endrun(91566);}
        }
        /* The tree just built is current by construction: the build set every
         * node's Ti_current to All.Ti_Current and refilled the whole SoA mirror,
         * local and foreign, from those nodes.  Record that here, at the call
         * site that KNOWS this is the main step tree, rather than inside
         * force_treebuild -- which is also used to build group-local and subset
         * trees whose geometry must never be certified for the step's device
         * consumers.  The full-drift test is the remaining precondition: it is
         * what makes the node geometry describe this time rather than merely
         * being freshly written, and move_particles above drifts only the active
         * set.  Absent that proof nothing is recorded and consumers fall back to
         * the host, which is correct but slower -- never wrong.
         * Recorded AFTER the bad-stop drain below: force_treebuild can request a
         * controlled stop during GPU finalize / LET / pseudo handling and still
         * return, and a tree whose build asked to stop must never be recorded as
         * current -- not even for the few statements before the poll exits. */
        gizmo_exit_bad_stop_if_requested("gravtree:after_treebuild"); CPU_Step[CPU_TREEBUILD] += measure_time(); /* and sync after treebuild as well */
        if(gizmo_full_drift_ti() == All.Ti_Current) {gpu_gravity_tree_mark_born_current(All.Ti_Current);}
        report_memory_ledger_on_growth("post-treebuild");  /* after force_treebuild (LET exchange ran); rebuild-only all-rank boundary */
        TreeReconstructFlag = 0;
        TreeMomentsStaleFlag = 0;
        All.NumForcesSinceLastTreeBuild = 0;   /* the counter this build answers */
        PRINT_STATUS(" ..Tree construction done.");
    }

    /* refresh tree moments if stale (e.g. after star formation or sink mass change).
       This must run before ANY gravity evaluation including Hermite calls, since
       stale moments produce wrong forces. Much cheaper than a full treebuild. */
    {
        int TreeMomentsStaleFlag_global;
        MPI_Allreduce(&TreeMomentsStaleFlag, &TreeMomentsStaleFlag_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(TreeMomentsStaleFlag_global)
        {
            CPU_Step[CPU_MISC] += measure_time();
            force_refresh_node_moments();
            gizmo_exit_bad_stop_if_requested("gravtree:after_refresh_moments"); /* drain refresh bad-stop before any gravity walk */
            CPU_Step[CPU_TREEBUILD] += measure_time();
            TreeMomentsStaleFlag = 0;
        }
    }

    CPU_Step[CPU_TREEMISC] += measure_time(); t0 = my_second(); double child0_span = CPU_ChildCharged;
#ifndef SELFGRAVITY_OFF
    /* allocate buffers to arrange communication */
    PRINT_STATUS(" ..Begin tree force. (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
    /* These two tables only RECORD targets whose gravity the locally-built tree cannot supply:
     * the MPI export round-trip is retired (gravity runs on the target-owning rank via the GPU
     * pre-pass and the locally-essential tree), so nothing here is ever sent, and a single entry
     * is enough to stop the run below. They are therefore held to a fixed small capacity rather
     * than to the communication chunk size, which used to reserve ~100 MB per rank on every
     * gravity call for something a healthy run never writes to at all. */
    All.BunchSize = GRAVITY_LET_DETECTOR_ENTRIES;
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    int k, ewald_max, diff, ndone, ndone_flag, place, recvTask; double tstart, tend, ax, ay, az; MPI_Status status;
    Ewaldcount = 0; Costtotal = 0; N_nodesinlist = 0; ewald_max=0;
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
    ewald_max = 1; /* the tree-code will need to iterate to perform the periodic boundary condition corrections */
#endif

    if(GlobNumForceUpdate > All.TreeRebuild_ActiveFraction * All.TotNumPart)
    { /* we have a fresh tree and would like to measure gravity cost */
        /* find the closest level */
        for(i = 1, TakeLevel = 0, diff = abs(All.LevelToTimeBin[0] - All.HighestActiveTimeBin); i < GRAVCOSTLEVELS; i++)
        {
            if(diff > abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin))
                {TakeLevel = i; diff = abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin);}
        }
        if(diff != 0) /* we have not found a matching slot */
        {
            if(All.HighestOccupiedTimeBin - All.HighestActiveTimeBin < GRAVCOSTLEVELS)	/* we should have space */
            {
                /* clear levels that are out of range */
                for(i = 0; i < GRAVCOSTLEVELS; i++)
                {
                    if(All.LevelToTimeBin[i] > All.HighestOccupiedTimeBin) {All.LevelToTimeBin[i] = 0;}
                    if(All.LevelToTimeBin[i] < All.HighestOccupiedTimeBin - (GRAVCOSTLEVELS - 1)) {All.LevelToTimeBin[i] = 0;}
                }
            }
            for(i = 0, TakeLevel = -1; i < GRAVCOSTLEVELS; i++)
            {
                if(All.LevelToTimeBin[i] == 0)
                {
                    All.LevelToTimeBin[i] = All.HighestActiveTimeBin;
                    TakeLevel = i;
                    break;
                }
            }
            if(TakeLevel < 0 && All.HighestOccupiedTimeBin - All.HighestActiveTimeBin < GRAVCOSTLEVELS)	/* we should have space */
                {
                    if(ThisTask == 0) {printf("TakeLevel < 0, even though we should have a slot\n"); fflush(stdout);}
                    endrun(90001008);
                    gizmo_exit_bad_stop_if_requested("gravtree:takelevel_no_slot");  /* symmetric (global LevelToTimeBin + bins): all ranks poll together */
                }
        }
    }
    else
    { /* in this case we do not measure gravity cost. Check whether this time-level
         has previously mean measured. If yes, then delete it so to make sure that it is not out of time */
        for(i = 0; i < GRAVCOSTLEVELS; i++) {if(All.LevelToTimeBin[i] == All.HighestActiveTimeBin) {All.LevelToTimeBin[i] = 0;}}
        TakeLevel = -1;
    }
    if(TakeLevel >= 0) {
        /* Under UVM-canonical particles, arena_P aliases host P[] — an arena
         * mirror write to P_arena_zero[i] is a self-assignment, so the single
         * host write is sufficient for both views. */
        for(i = 0; i < NumPart; i++) { P[i].GravCost[TakeLevel] = 0; }
    } /* re-zero the cost [will be re-summed] */

    /* Decide which particles need a new tree force BEFORE any walk runs, ONCE for the whole call,
       and let every consumer read that one answer. The walk MUTATES the inputs -- it writes
       Min_Sink_FeedbackTime, which needs_new_treeforce() compares against -- so asking again after a
       walk can give a different answer for the same particle. Three things ask: the walk dispatchers
       (via gravity_treewalk_candidate_prewalk), and the finalization loop below, which must take the
       jerk-skip path for exactly the particles the walk skipped. If they disagree, a particle the
       walk computed raw is finalized as though it carried a previous step's G-multiplied values, or
       one the walk never touched is multiplied by G a second time.
       The walk asks more than once per call in two ways that are both real: the Ewald-correction
       pass re-selects after the primary pass has already written those fields, and an import repair
       redoes both passes. */
#ifdef ADAPTIVE_TREEFORCE_UPDATE
    gravity_freeze_treeforce_candidates();
#endif

    /* Import-completeness repair.  The import is pruned when the tree is built, against where the
     * particles were then; the target positions and node centres it was pruned against keep moving
     * while the tree is reused, so a walk can come to resolve structure the import no longer
     * carries.  That cannot be frozen the way the opening scale can -- it is the physical support
     * of the force law, not an estimator -- so it is detected and repaired instead: rebuild the
     * tree and its import against the current positions, and redo this evaluation against them.
     *
     * The whole Ewald_iter loop is redone, in order, because the primary walk ASSIGNS its outputs
     * (so the redo is its own reset) while the Ewald correction ADDS to them, and both run before
     * the finalization below mutates GravAccel in place.
     *
     * ONE repair.  A shortfall that survives a rebuild against current positions is not the tree
     * falling behind the particles, and rebuilding again would not address it. */
    const int gravity_let_repair_max = 1;
    int gravity_let_repair_attempts = 0;
gravity_walk_attempt:

    /* begin main communication and tree-walk loop. note the ewald-iter terms here allow for multiple iterations for periodic-tree corrections if needed */
    for(Ewald_iter = 0; Ewald_iter <= ewald_max; Ewald_iter++)
    {
        NextParticle = 0;	/* begin with this index */
        memset(ProcessedFlag, 0, All.MaxPart * sizeof(unsigned char));
        BufferCollisionFlag = 0; /* set to zero before operations begin */

        /* Speculative GPU pre-pass: walks the local tree on GPU for each
         * active particle; on success, writes GravAccel and marks
         * ProcessedFlag so the CPU primary loop below skips it.
         * On pseudo-particle hit, leaves the particle untouched for the
         * CPU loop + MPI export machinery to handle unchanged. Ewald_iter
         * splits primary (==0) vs Ewald-correction (==1) walks; both are
         * active on all Kokkos builds. */
        /* Time the device walk.  It is the primary work of this routine and it sat outside every
         * timer: the published work-load balance reduces timetree1, which brackets only the CPU
         * leftover loop below, so on an accelerated build it reported milliseconds for a walk
         * costing seconds -- and, far worse, a reassuring balance figure for the most imbalanced
         * phase in the run.  The two accumulators already mean "primary walk" and "Ewald-correction
         * walk"; this gives them their real contents, so work-load balance and rel1to2 describe the
         * walk that actually happened on both accelerated and host-only builds. */
        tstart = my_second();
        int host_candidates = 0;   /* the Ewald-correction walk takes every target alone, so its host loop needs no packets */
        if(Ewald_iter == 0) {gpu_gravtree_walk_primary(&host_candidates);}
        else                {gpu_ewald_walk_primary();}
        tend = my_second();
        if(Ewald_iter == 0) {timetree1 += timediff(tstart, tend);}
        else                {timetree2 += timediff(tstart, tend);}
        gravity_walk_workspace_allocate(host_candidates);

        do /* primary point-element loop */
        {
            iter++;
            BufferFullFlag = 0; Nexport = 0; tstart = my_second();

#ifdef _OPENMP
#pragma omp parallel
#endif
            {
#ifdef _OPENMP
                int mainthreadid = omp_get_thread_num();
#else
                int mainthreadid = 0;
#endif
                gravity_primary_loop(&mainthreadid);	/* do local particles and prepare export list */
            }
            tend = my_second(); timetree1 += timediff(tstart, tend);

            /* ============================================================
             * CPU gravity export round-trip is RETIRED.
             * ------------------------------------------------------------
             * The GPU pre-pass + Locally Essential Tree supply all
             * foreign-rank gravity on the target-owning rank, so the legacy
             * MPI export round-trip (compaction, MPI_Alltoall/Sendrecv,
             * imported-particle walk, scatter-back) and its gravdata_in/out
             * buffers are fully removed.  What survives is a DETECTOR:
             * DataIndexTable/DataNodeList are no longer shipped -- the tree
             * walk only records LET-incompleteness into Nexport.
             *
             * If Nexport > 0 a particle's gravity is not covered by LET, and
             * we request a graceful controlled-stop here (drained at the
             * all-rank poll below, no retry, no new collective).  It is not a
             * capacity shortfall: an import too large for the foreign-node
             * index range is reported by the exchange, which raises the range
             * and rebuilds before this walk runs.
             * ============================================================ */
            /* Import completeness is NOT resolved here.  Every walk of this pass records into one
             * ledger; whether the shortfall is repairable is a collective question, so it is
             * reduced across ranks and answered once after the pass, below the Ewald_iter loop. */

            if(Nexport > 0) {
                printf("The locally essential tree did not cover the gravity of %ld particles on rank %d. Stopping.\n", Nexport, ThisTask);
                fflush(stdout);
                /* Graceful soft-stop: the export round-trip is retired, so we cannot service
                 * these particles -- but the same loop iteration reaches the all-rank
                 * MPI_Allreduce + gizmo_exit_bad_stop_if_requested poll below, which drains
                 * this flag cleanly (no retry, no new collective). */
                gizmo_request_controlled_stop(914040, "gravtree: the locally essential tree did not cover some targets' gravity", __FILE__, __LINE__, __FUNCTION__);
            }

            /* Export-back loop is retired under GPU offload, so the arena
             * is not invalidated by host-side P[] writes here.
             * There is no gpu_particles_arena_invalidate() call at this point:
             * the arena is coherent here (gpu_gravtree_walk_primary invalidates
             * internally if its host scatter happens). A redundant
             * double-invalidate would cost nothing but blocks fast-path
             * acquires after the arena refresh, so it is intentionally absent.
             *
             * If a pre-acquire host mutation ever makes the arena stale at this
             * point, the GIZMO_GPU_ARENA_DEBUG=1 byte-compare guard will abort
             * with the site name -- rely on that trip wire instead of a
             * cargo-cult invalidate. */
            if(NextParticle >= (int)ActiveParticleList.size()) {ndone_flag = 1;} else {ndone_flag = 0;} /* figure out if we are done with the particular active set here */
            tstart = my_second();
            MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /* call an allreduce to figure out if all tasks are also done here, otherwise we need to iterate */
            tend = my_second(); timewait2 += timediff(tstart, tend);
            gizmo_exit_bad_stop_if_requested("gravtree:tree_export_loop"); /* drain a buffer-too-small bad-stop here instead of retrying the export with zero progress */
        }
        while(ndone < NTask);
        gravity_walk_workspace_free();
    } /* Ewald_iter */

    /* Resolve this pass's import completeness.  Every rank left the loop above only once all of
     * them were done, so they arrive here together and reduce the same question; the decision is
     * taken from the reduced value alone, so every rank takes the same branch and the rebuild
     * below stays collective. */
    {
        long long shortfall_local = gravity_incomplete_import_count();
        long long unshippable_local = gravity_unshippable_import_count();
        long long tally[3] = {shortfall_local, (shortfall_local > 0) ? 1 : 0, unshippable_local};
        long long total[3] = {0, 0, 0};
        MPI_Allreduce(tally, total, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if(total[0] > 0)
        {
            /* A shortfall is worth a rebuild only if a rebuild could answer it.  Where the walk
             * opened a subtree the sender never owned, the import is short for a topological
             * reason: the same pack against current positions reaches the same conclusion, so the
             * repair would cost a full tree and exchange and then stop here anyway.  Say so at the
             * first pass instead of at the second. */
            int all_unrepairable = (total[2] >= total[0]);
            int repairing = (gravity_let_repair_attempts < gravity_let_repair_max) && !all_unrepairable;
            /* One line for the whole repair, not one per rank: the decision is collective and every
             * rank that saw a shortfall would otherwise say the same thing about it.  Pick the
             * lowest-numbered rank that actually has an example -- a rank whose shortfall came only
             * from a device walk has a count but nothing to show -- and let it supply the example. */
            char example[256]; example[0] = '\0';
            int have = gravity_incomplete_import_example(example, (int) sizeof(example));
            int pick[2] = {have ? ThisTask : NTask, ThisTask}, chosen[2] = {NTask, 0};
            MPI_Allreduce(pick, chosen, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
            if(chosen[0] < NTask) {MPI_Bcast(example, (int) sizeof(example), MPI_CHAR, chosen[0], MPI_COMM_WORLD);}
            else {example[0] = '\0';}
            if(ThisTask == 0)
            {
                if(repairing) {
                    printf("The gravity walk needed to descend %lld imported node(s) that arrived without children, on "
                           "%lld of %d ranks. The import is pruned when the tree is built, against where the particles "
                           "were then, and they have since moved far enough that the walk resolves structure it no "
                           "longer carries.%s%s Rebuilding the tree and redoing this evaluation against it. Frequent "
                           "repairs mean the tree is being reused too long -- lower TreeRebuild_ActiveFraction.\n",
                           total[0], total[1], NTask, example[0] ? " First case: " : "", example);
                } else if(all_unrepairable) {
                    printf("The gravity walk needed to descend %lld imported node(s) that arrived without children, "
                           "on %lld of %d ranks, and %lld of them are subtrees the sending rank never owned: every "
                           "child was a pseudo-particle or another rank's foreign node, so no rank could have "
                           "shipped them whole.%s%s Rebuilding reproduces that exactly. "
                           "Reaching this means a target opens a node whose contents live on a third rank, which "
                           "the one-shot import does not express. Stopping.\n",
                           total[0], total[1], NTask, total[2], example[0] ? " First case: " : "", example);
                } else {
                    printf("The gravity walk still needed to descend %lld imported node(s) that arrived without "
                           "children, on %lld of %d ranks, after the tree and its import were rebuilt against the "
                           "current particle positions.%s%s That is not the tree falling behind the particles, so "
                           "rebuilding again would not fix it -- the receiver cover, the wire format or the import "
                           "install is not shipping what this walk opens. Stopping.\n",
                           total[0], total[1], NTask, example[0] ? " First case: " : "", example);
                }
                fflush(stdout);
            }
            gravity_clear_incomplete_import();
            if(!repairing)
            {
                /* Drain here rather than letting this fall through.  Every rank reaches this on the
                 * same reduced value, so the poll is symmetric, and the finalization below would
                 * otherwise run its in-place mutation of GravAccel, the potential and the RT fields
                 * over a pass that has just been declared not computable. */
                gizmo_request_controlled_stop(90000087, "gravtree: the imported tree did not carry the structure the walk resolved",
                                              __FILE__, __LINE__, __FUNCTION__);
                gizmo_exit_bad_stop_if_requested("gravtree:let_repair_exhausted");
            }
            else
            {
                const double t_repair_start = my_second();
                const double child0_repair = CPU_ChildCharged;
                gravity_let_repair_attempts++;

                /* The detector tables were allocated after the tree, so they sit above it in the
                 * arena.  force_treebuild frees and reallocates the tree when it has to grow the
                 * node arena or the foreign-node range, and the arena is LIFO -- so hand those two
                 * back first and take them again afterwards. */
                myfree(DataNodeList); myfree(DataIndexTable);

                /* Same build the start of this call would have done, minus the drift and
                 * re-sequencing: the particles are already at All.Ti_Current, and re-sequencing
                 * here would move indices under the active list this walk is iterating.  The
                 * rebuild flags are deliberately NOT cleared -- this repair does not satisfy
                 * whatever else asked for a rebuild, and the next step is entitled to see it.
                 * It also stands on the domain frame already in force, so under RANDOMIZE_GRAVTREE
                 * it does not draw a new one: a repair is not a scheduled rebuild, and re-keying
                 * the particles here would move them out from under the walk in progress. */
                refresh_old_acceleration_for_tree_opening();
                gizmo_exit_bad_stop_if_requested("gravtree:before_repair_treebuild");
                /* This rebuild stands on the tree just built here, whose attachments are intact, so it
                 * cannot ask for ownership to be restored -- and must not, with a walk in progress. */
                if(force_treebuild(NumPart, NULL) < 0) {endrun(91567);}
                gizmo_exit_bad_stop_if_requested("gravtree:after_repair_treebuild");
                if(gizmo_full_drift_ti() == All.Ti_Current) {gpu_gravity_tree_mark_born_current(All.Ti_Current);}
                TreeMomentsStaleFlag = 0;   /* the build just refreshed every moment */

                All.BunchSize = GRAVITY_LET_DETECTOR_ENTRIES;
                DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
                DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

                /* The redone pass re-counts its own work, so drop what the abandoned one counted.
                 * GravCost is not just a diagnostic -- it is the per-particle weight the next
                 * domain decomposition balances on.  Only the active set can have been written, so
                 * only the active set is cleared: a step with few actives must not pay for NumPart. */
                Costtotal = 0; Ewaldcount = 0; N_nodesinlist = 0;
                if(TakeLevel >= 0) {for(int ii = 0; ii < (int)ActiveParticleList.size(); ii++) {P[ActiveParticleList[ii]].GravCost[TakeLevel] = 0;}}

                /* The repair goes through the same build and the same walks as any other, so its
                 * cost belongs in the same rows; charge the build here so the walk row it sits
                 * inside does not absorb it. */
                cpu_charge_child(CPU_TREEBUILD, cpu_minus_children(timediff(t_repair_start, my_second()), child0_repair));
                goto gravity_walk_attempt;
            }
        }
    }

    myfree(DataNodeList); myfree(DataIndexTable);

    /* assign node cost to particles */
    if(TakeLevel >= 0) {
        /* Modern GPU/LET gravity executes work on the target-owning rank, so
         * gpu_gravtree_walk_primary() records target-side interaction counts
         * directly in P[target].GravCost[TakeLevel]. */
    }


    /* now perform final operations on results [communication loop is done] */
#ifndef GRAVITY_HYBRID_OPENING_CRIT  // in collisional systems we don't want to rely on the relative opening criterion alone, because aold can be dominated by a binary companion but we still want accurate contributions from distant nodes. Thus we combine BH and relative criteria. - MYG
    /* Switch to the relative opening criterion for the following force computations.
     * (Second-order ICs keep Barnes-Hut on the very first step.) */
    double errtol_before = All.ErrTolTheta;
    int enable_relative_opening =
        (All.TypeOfOpeningCriterion == 1) &&
        !(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0);
    if(enable_relative_opening) { All.ErrTolTheta = 0; }
    /* The opening criterion just changed; the installed LET was built/exported under the previous
     * criterion and is not valid for the next walk. Rebuild the tree+LET before it is reused. */
    if(errtol_before != 0 && All.ErrTolTheta == 0) { TreeReconstructFlag = 1; }
#endif

#ifdef SINGLE_STAR_DIRECT_GRAVITY
    /* the exact star-star sum, which the tree walk above deliberately left out. Here and not later:
       the tree's own contributions are complete (imports included) but not yet multiplied by All.G
       in the loop below, and the direct sum is in those same G-free units. */
    CPU_Step[CPU_TREEMISC] += measure_time();
    star_direct_gravity_build_table();
    star_direct_gravity_compute();
    star_direct_gravity_free_table();
    CPU_Step[CPU_TREEWALK1] += measure_time();
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int ii = 0; ii < (int)ActiveParticleList.size(); ii++)
    {
        int i = ActiveParticleList[ii];
#ifdef HERMITE_INTEGRATION
        if(HermiteOnlyFlag) {if(!eligible_for_hermite(i, P)) continue;} /* if we are completing an extra loop required for the Hermite integration, all of the below would be double-calculated, so skip it */
#endif      
#ifdef ADAPTIVE_TREEFORCE_UPDATE
        double dt = get_particle_timestep_in_physical(i, P);
        if(!gravity_treeforce_candidate_frozen(ii)) { // the one decision this call froze before any walk ran, the same one the walk dispatchers used
            P[i].GravAccel += P[i].GravJerk * (dt * All.cf_a2inv); // a^-1 from converting velocity term in the jerk to physical; a^-3 from the 1/r^3; a^2 from converting the physical dt * j increment to GravAccel back to the units for GravAccel; result is a^-2; note that Ewald and PMGRID terms are neglected from the jerk at present
            P[i].time_since_last_treeforce += dt;
            continue;
        } else {
            P[i].time_since_last_treeforce = dt;
        }
#endif
        /* before anything: multiply by G for correct units [be sure operations above/below are aware of this!] */
        P[i].GravAccel *= All.G;
#if (SINGLE_STAR_TIMESTEPPING > 0)
        P[i].COM_GravAccel *= All.G;
#endif

#ifdef EVALPOTENTIAL
        P[i].Potential *= All.G;
#ifdef BOX_PERIODIC
        if(All.ComovingIntegrationOn) {P[i].Potential -= All.G * 2.8372975 * pow(P[i].Mass, 2.0 / 3) * pow(All.OmegaMatter * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G), 1.0 / 3);} else {if(All.OmegaLambda>0) {P[i].Potential -= 0.5*All.OmegaLambda*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits * (P[i].Pos.norm_sq());}}
#endif
#ifdef PMGRID
        P[i].Potential += P[i].PM_Potential; /* add in long-range potential */
#endif
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
        P[i].TreeMass += P[i].Mass;
        if(P[i].Type == 5) printf("Particle %d sees mass %g in the gravity tree\n", P[i].ID, P[i].TreeMass);
#endif

#ifdef SPECIAL_POINT_WEIGHTED_MOTION
        if(P[i].Type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
        {
            P[i].vel_of_nearest_special /= P[i].weight_sum_for_special_point_smoothing;
            P[i].acc_of_nearest_special /= P[i].weight_sum_for_special_point_smoothing;
            /* now reset the local values for this to actually match these, recalling the special particle in this module is just a tracer element */
            double dtime_phys = (All.Time - P[i].Time_Of_Last_SmoothedVelUpdate) / All.cf_hubble_a; /* want to convert to physical units */
            if(dtime_phys > 0) {
                P[i].Acc_Total_PrevStep = (P[i].vel_of_nearest_special - P[i].Vel) / (All.cf_atime * dtime_phys * All.cf_a2inv); /* converting to cosmological units here */
                P[i].Vel = P[i].vel_of_nearest_special;
            }
        }
#endif

        /* Measure the acceleration scale the relative tree-opening criterion will use, HERE: after
         * the G multiplication above, and before the companion subtraction, radiation pressure and
         * analytic-gravity terms below.  That keeps it a property of the force the gravity tree
         * computes -- the error the opening criterion exists to control -- rather than of every
         * force acting on the particle, which would let rapidly varying radiation or an external
         * field decide how finely the tree is opened.  It is only staged here; the tree build
         * promotes it into OldAcc, so the value cannot shift under a tree whose import was pruned
         * against it.  (Particles that skipped the walk above never reach this and keep the scale
         * from their last real walk, as before.) */
        if(!(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0)) /* 2lpt ICs keep masses in OldAcc until the first evaluation */
        {
            auto accel_for_opening = P[i].GravAccel;
#ifdef PMGRID
            accel_for_opening += P[i].GravPM;
#endif
            P[i].OldAcc_LatestWalk = accel_for_opening.norm() / All.G;   /* back to non-G units, matching the predicate */
        }

#if (SINGLE_STAR_TIMESTEPPING > 0) /* Subtract component of force from companion if in binary, because we will operator-split this */
        if((P[i].Type == 5) && (P[i].is_in_a_binary == 1)) {subtract_companion_gravity(i);}
#endif

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* final operations to compute the tidal tensor and related quantities */
        P[i].tidal_tensorps *= All.G; /* give this the proper units */
#ifdef COMPUTE_JERK_IN_GRAVTREE
        P[i].GravJerk *= All.G; /* units */
#endif
#if defined(PMGRID) && !defined(ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)
        P[i].tidal_tensorps += P[i].tidal_tensorpsPM; /* add the long-range (pm-grid) contribution; but make sure to do this after the unit multiplication by G above, since the PM term already has G built into it */
#endif
#endif /* COMPUTE_TIDAL_TENSOR_IN_GRAVTREE */

#if defined(RT_OTVET) /* normalize the Eddington tensors we just calculated by walking the tree (normalize to trace=1) */
        if(P[i].Type == 0) {
            int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
            {double trace = CellP[i].ET[k_freq].trace();
                if(!isnan(trace) && (trace>0)) {CellP[i].ET[k_freq] /= trace;} else {CellP[i].ET[k_freq].set_isotropic(1./3.);}}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) /* normalize to energy density with C, and multiply by volume to use standard 'finite volume-like' quantity as elsewhere in-code */
        if(P[i].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[i].Rad_E_gamma[kf] *= P[i].Mass/(CellP[i].Density*All.cf_a3inv * C_LIGHT_CODE_REDUCED);}}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        if(P[i].Type==0) {CellP[i].SubGrid_CosmicRayEnergyDensity *= cr_get_source_shieldfac(i, P, CellP);}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX) /* multiply by volume to use standard 'finite volume-like' quantity as elsewhere in-code */
        if(P[i].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[i].Rad_Flux[kf] *= P[i].Mass/(CellP[i].Density*All.cf_a3inv);}} // convert to standard finite-volume-like units //
#if !defined(RT_DISABLE_RAD_PRESSURE) // if we save the fluxes, we didnt apply forces on-the-spot, which means we appky them here //
        if((P[i].Type==0) && (P[i].Mass>0)
#ifdef HYDRO_MULTIFLUID_DM
           && (P[i].FluidType != FLUID_DM)   /* dark fluid does not feel baryonic RT radiation pressure */
#endif
          )
        {
            int kfreq; double vol_inv=CellP[i].Density*All.cf_a3inv/P[i].Mass, h_i=P[i].Get_Particle_Size()*All.cf_atime, sigma_eff_i=P[i].Mass/(h_i*h_i);
            Vec3<double> radacc={};
            for(kfreq=0; kfreq<N_RT_FREQ_BINS; kfreq++)
            {
                double f_slab=1, erad_i=0, kappa_rad=rt_kappa(i,kfreq, P, CellP), tau_eff=kappa_rad*sigma_eff_i; if(tau_eff > 1.e-4) {f_slab = (1.-exp(-tau_eff)) / tau_eff;} // account for optically thick local 'slabs' self-shielding themselves
                double acc_norm = kappa_rad * f_slab / C_LIGHT_CODE_REDUCED; // pre-factor for radiation pressure acceleration
#if defined(RT_LEBRON)
                acc_norm *= All.PhotonMomentum_Coupled_Fraction; // allow user to arbitrarily increase/decrease strength of RP forces for testing
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                erad_i = CellP[i].Rad_E_gamma_Pred[kfreq]*vol_inv; // if can, include the O[v/c] terms
#endif
                Vec3<double> flux_i = CellP[i].Rad_Flux_Pred[kfreq] * vol_inv;
                Vec3<double> vel_i = CellP[i].VelPred * (1.0/All.cf_atime);
                double flux_mag2 = flux_i.norm_sq() + MIN_REAL_NUMBER, vdotflux = dot(vel_i, flux_i); // initialize a bunch of variables we will need
                Vec3<double> vdot_h = (vel_i + flux_i * (vdotflux/flux_mag2)) * erad_i; // calculate volume integral of scattering coefficient t_inv * (gas_vel . [e_rad*I + P_rad_tensor]), which gives an additional time-derivative term. this is the P term //
                radacc += (flux_i - vdot_h) * acc_norm; // note these 'vdoth' terms shouldn't be included in FLD, since its really assuming the entire right-hand-side of the flux equation reaches equilibrium with the pressure tensor, which gives the expression in rt_utilities
            }
#if defined(RT_RAD_PRESSURE_OUTPUT)
            CellP[i].Rad_Accel = radacc; // here units are the same as hydroaccel, so no extra comoving units
#else
            P[i].GravAccel += radacc * (1.0/All.cf_a2inv); // convert into our code units for GravAccel, which are comoving gm/r^2 units //
#endif
        }
#endif
#endif

#ifdef RT_USE_TREECOL_FOR_NH  /* compute the effective column density that gives equivalent attenuation of a uniform background: -log(avg(exp(-tau)))/kappa */
        double attenuation=0, minimum_column=MAX_REAL_NUMBER; int kbin;
        double kappa_photoelectric = 500. * DMAX(1e-4, (P[i].Metallicity[0]/All.SolarAbundances[0])*return_dust_to_metals_ratio_vs_solar(i,0, P, CellP)); // dust opacity in cgs
        for(kbin=0; kbin<RT_USE_TREECOL_FOR_NH; kbin++) {
	      attenuation += exp(DMAX(-P[i].ColumnDensityBins[kbin] * UNIT_SURFDEN_IN_CGS * kappa_photoelectric,-100));
	      minimum_column = DMIN(minimum_column,P[i].ColumnDensityBins[kbin]);
	    } // we put a floor here to avoid underflow errors where exp(-large) = 0 - will just return a very high surface density that will be in the highly optically thick regime where both the ISRF and cooling radiation escape will be negligible
        P[i].SigmaEff = -log(attenuation/RT_USE_TREECOL_FOR_NH) / (kappa_photoelectric * UNIT_SURFDEN_IN_CGS);
	    if(P[i].SigmaEff < minimum_column) {P[i].SigmaEff = minimum_column;} // if in the overflowing regime just take the minimum column density to extrapolate better to the IR-thick regime
#ifdef GIZMO_TREECOL_DIAG
        if(ThisTask==0 && P[i].Type==0 && i<10 && All.NumCurrentTiStep<4) {
            double csum=0; for(int kb=0;kb<RT_USE_TREECOL_FOR_NH;kb++) csum+=P[i].ColumnDensityBins[kb];
            printf("TREECOL_DIAG step=%d i=%d ID=%llu SigmaEff=%g binsum=%g bins=[%g,%g,%g,%g,%g,%g] LET=%d\n",
                   All.NumCurrentTiStep,(int)i,(unsigned long long)P[i].ID,P[i].SigmaEff,csum,
                   P[i].ColumnDensityBins[0],P[i].ColumnDensityBins[1],P[i].ColumnDensityBins[2],
                   P[i].ColumnDensityBins[3],P[i].ColumnDensityBins[4],P[i].ColumnDensityBins[5],
                   (MaxForeignNodes>0)?1:0); fflush(stdout);
        }
#endif
#endif

#if !defined(BOX_PERIODIC) && !defined(PMGRID) /* some factors here in case we are trying to do comoving simulations in a non-periodic box (special use cases) */
        if(All.ComovingIntegrationOn) {P[i].GravAccel += P[i].Pos * (0.5*All.OmegaMatter*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits);}
        if(All.ComovingIntegrationOn==0) {P[i].GravAccel += P[i].Pos * (All.OmegaLambda*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits);}
#ifdef EVALPOTENTIAL
        if(All.ComovingIntegrationOn) {P[i].Potential -= 0.5*All.OmegaMatter*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits * P[i].Pos.norm_sq();}
#endif
#endif

    } /* end of loop over active particles*/

    /* Arena mirror-update is a no-op under UVM-canonical (arena_P
     * aliases host P[]); the post-loop above already wrote canonical state. */

#endif /* end SELFGRAVITY operations (check if SELFGRAVITY_OFF not enabled) */


    add_analytic_gravitational_forces(); /* add analytic terms, which -CAN- be enabled even if self-gravity is not */


    /* Now the force computation is finished: gather timing and diagnostic information */
    t1 = my_second(); cpu_chain_sync(t1); timeall = cpu_minus_children(timediff(t0, t1), child0_span);
    timetree = timetree1 + timetree2; timewait = timewait1 + timewait2; timecomm = timecommsumm1 + timecommsumm2;
    MPI_Reduce(&timetree, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timetree, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timetree1, &sumt1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timetree1, &maxt1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timetree2, &sumt2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timetree2, &maxt2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timewait, &sumwaitall, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timecomm, &sumcommall, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Costtotal, &sum_costtotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Ewaldcount, &ewaldtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    sumup_longs(1, &n_exported, &n_exported);
    sumup_longs(1, &N_nodesinlist, &N_nodesinlist);
    All.TotNumOfForces += GlobNumForceUpdate;
    plb = (NumPart / ((double) All.TotNumPart)) * NTask;
    MPI_Reduce(&plb, &plb_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    /* Packets the device engine gave up on this call, by reason, summed over ranks and also taken
     * at its worst rank: a budget cliff is a property of the rank with the divergent subtree, and
     * a sum over ranks alone would average it away.  Collective, so it sits with the other
     * reductions rather than inside the rank-0 report below. */
    /* The two stamped recorders' fail-safe totals ride in the same two reductions rather than
     * adding their own.  Each exists so that a permanent silent revert -- to sweeping every node,
     * or to drifting every particle -- is visible rather than indistinguishable from the
     * optimisation working, and nothing read either of them; but observability may not charge the
     * path it is watching, and this call already carries collectives.  Three extra words on an
     * existing reduction is free; two more reductions per gravity call would not be. */
    /* Scoped, not #define: a macro here would be translation-unit-wide despite sitting in a
       function. */
    constexpr int recorder_slots = 5;
    constexpr int report_slots   = GRAV_PACKET_FAIL_REASON_SLOTS + recorder_slots;
    long long packet_fail[report_slots] = {0};
    long long packet_fail_sum[report_slots] = {0};
    long long packet_fail_max[report_slots] = {0};
    gpu_gravtree_packet_failures(packet_fail, GRAV_PACKET_FAIL_REASON_SLOTS);
    packet_fail[GRAV_PACKET_FAIL_REASON_SLOTS + 0] = gpu_node_dirty_unsafe_events();
    packet_fail[GRAV_PACKET_FAIL_REASON_SLOTS + 1] = gx_touched_set_refused_epochs();
    packet_fail[GRAV_PACKET_FAIL_REASON_SLOTS + 2] = gx_touched_set_retire_faults();
    gpu_gravtree_subset_drift_counts(&packet_fail[GRAV_PACKET_FAIL_REASON_SLOTS + 3],
                                     &packet_fail[GRAV_PACKET_FAIL_REASON_SLOTS + 4]);
    MPI_Reduce(packet_fail, packet_fail_sum, report_slots, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(packet_fail, packet_fail_max, report_slots, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Numnodestree, &maxnumnodes, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    /* The span from the end of the build to here is the force walk.  Only the
     * wait is separately measured inside it, so the walk row is the rest of the
     * span; charging the remainder to `misc` instead reported the whole walk as
     * unattributed.  The send/recv and second-walk timers this routine used to
     * split out have no writer any more and are not charged. */
    CPU_Step[CPU_TREEWALK1] += timeall - timewait2;
    CPU_Step[CPU_TREEWAIT2] += timewait2;
#ifdef OUTPUT_ADDITIONAL_RUNINFO
    if(ThisTask == 0)
    {
        fprintf(FdTimings, "Step= %lld  t= %.16g  dt= %.16g \n",(long long) All.NumCurrentTiStep, All.Time, All.TimeStep);
        fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g (%g) iter= %d\n", (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000), n_exported / ((double) GlobNumForceUpdate), N_nodesinlist / ((double) n_exported + 1.0e-10), iter); /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */
        fprintf(FdTimings, "work-load balance: %g (%g %g) rel1to2=%g   max=%g avg=%g\n", maxt / (1.0e-6 + sumt / NTask), maxt1 / (1.0e-6 + sumt1 / NTask), maxt2 / (1.0e-6 + sumt2 / NTask), sumt1 / (1.0e-6 + sumt1 + sumt2), maxt, sumt / NTask);
        fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
        fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes, maxnumnodes / ((double) MaxNodes));
        fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", GlobNumForceUpdate / (sumt + 1.0e-20), GlobNumForceUpdate / (1.0e-6 + maxt * NTask), ((double) (sum_costtotal)) / (1.0e-20 + GlobNumForceUpdate), ((double) ewaldtot) / (1.0e-20 + GlobNumForceUpdate)); {int packet_team, packet_q_dev; gpu_gravtree_packet_shape(&packet_team, &packet_q_dev); fprintf(FdTimings, "packet: Q=%d T=%d Qdev=%d\n", TREE_QUERY_PACKET_SIZE, packet_team, packet_q_dev);
            fprintf(FdTimings, "packet-gaveup: malformed=%lld stale=%lld pseudo=%lld nocont=%lld record=%lld (worst rank: %lld %lld %lld %lld %lld)\n",
                    packet_fail_sum[1], packet_fail_sum[2], packet_fail_sum[3], packet_fail_sum[4], packet_fail_sum[5],
                    packet_fail_max[1], packet_fail_max[2], packet_fail_max[3], packet_fail_max[4], packet_fail_max[5]);
            /* Only when there is something to say.  These are fail-safe EVENTS, not telemetry:
               a zero line on every call would be noise in the artifact the track reads, while a
               nonzero one is the whole point -- it says the run has quietly stopped using the
               mechanism being priced. */
            if(packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 0] ||
               packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 1] ||
               packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 2]) {
                fprintf(FdTimings, "recorder-failsafe: node-unsafe=%lld touched-refused=%lld touched-retire-faults=%lld (worst rank: %lld %lld %lld)\n",
                        packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 0],
                        packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 1],
                        packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 2],
                        packet_fail_max[GRAV_PACKET_FAIL_REASON_SLOTS + 0],
                        packet_fail_max[GRAV_PACKET_FAIL_REASON_SLOTS + 1],
                        packet_fail_max[GRAV_PACKET_FAIL_REASON_SLOTS + 2]);
            }
            /* How the sources were brought current, on the same cadence as the packet shape
               beside it: this is the row's own quantity, not a fail-safe, so it prints whenever
               the device route ran rather than only when something went wrong. */
            if(packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 3] ||
               packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 4]) {
                fprintf(FdTimings, "subset-drift: taken=%lld declined=%lld (worst rank: %lld %lld)\n",
                        packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 3],
                        packet_fail_sum[GRAV_PACKET_FAIL_REASON_SLOTS + 4],
                        packet_fail_max[GRAV_PACKET_FAIL_REASON_SLOTS + 3],
                        packet_fail_max[GRAV_PACKET_FAIL_REASON_SLOTS + 4]);
            }} fprintf(FdTimings, "\n");
        fflush(FdTimings);
    }
    double costtotal_new = 0, sum_costtotal_new;
    if(TakeLevel >= 0)
    {
        for(i = 0; i < NumPart; i++) {costtotal_new += P[i].GravCost[TakeLevel];}
        MPI_Reduce(&costtotal_new, &sum_costtotal_new, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        /* Both walks accumulate the same per-target interaction count into GravCost and
         * into Costtotal, so the two totals describe the same quantity and this should be
         * at round-off whichever path each target took. A non-negligible value means the
         * two are no longer measuring the same thing. */
        if(sum_costtotal>0) {PRINT_STATUS(" ..relative error in the total number of tree-gravity interactions = %g", (sum_costtotal - sum_costtotal_new) / sum_costtotal);}
    }
#endif
    CPU_Step[CPU_TREEMISC] += measure_time();
}




void *gravity_primary_loop(void *p)
{
    int i, j, thread_id = *(int *) p, *exportflag, *exportnodecount, *exportindex;
    exportflag = Exportflag + thread_id * NTask; exportnodecount = Exportnodecount + thread_id * NTask; exportindex = Exportindex + thread_id * NTask;
    for(j = 0; j < NTask; j++) {exportflag[j] = -1;} /* Note: exportflag is local to each thread */
#ifdef _OPENMP
    if(BufferCollisionFlag && thread_id) {return NULL;} /* force to serial for this subloop if threads simultaneously cross the Nexport bunchsize threshold */
#endif
    /* this thread's workspace (see gravity_walk_workspace_allocate) */
    char *thread_ws = GravWalkWorkspace + (size_t) thread_id * GravWalkThreadBytes;
    int *batch = (int *) thread_ws, *batch_pos = batch + GravWalkBatchCap, *ninter = batch + 2 * GravWalkBatchCap;
    void *walk_ws = thread_ws + GravWalkBatchBytes;   /* present only when packet_cap > 0, i.e. when the primary walk has host candidates */
    const int packet_cap = GravWalkPacketCap;
    /* what a completed walk hands to the rest of the step: the interaction count is the work
     * weight for the next domain decomposition (the device walk records the same quantity,
     * gpu_gravtree.cc, so a step whose walks are split between the two paths feeds one
     * consistent measure to domain_particle_costfactor()); each thread writes only its own
     * targets, so no synchronization is needed beyond the shared total */
    auto commit_target = [&](int target, int n_interactions)
    {
        if(TakeLevel >= 0) {P[target].GravCost[TakeLevel] = n_interactions;}
#ifdef _OPENMP
#pragma omp atomic
#endif
        Costtotal += n_interactions;
        ProcessedFlag[target] = 1;
    };
    while(1)
    {
        /* The active-list position travels with the particle: the frozen candidacy this call
         * decided is keyed by it, and the batch would otherwise carry only the particle index. */
        int batch_count = 0;
#ifdef _OPENMP
#pragma omp critical(_nextlistgravprim_)
#endif
        {
            while(batch_count < GravWalkBatchCap && BufferFullFlag == 0 && NextParticle < (int)ActiveParticleList.size())
            {
                int pos = NextParticle, idx = ActiveParticleList[NextParticle]; NextParticle++;
                if(!ProcessedFlag[idx]) {batch_pos[batch_count] = pos; batch[batch_count++] = idx;}
            }
        }
        if(batch_count == 0) {break;}
        int buffer_full = 0;
        /* SSOT pre-walk candidacy (Mass>0 + Hermite eligibility + needs_new_treeforce);
         * non-candidates are marked done so the finalization loop skips them, and the
         * candidates close up in place so consecutive ones form a packet below. */
        int n_candidates = 0;
        for(int b = 0; b < batch_count; b++)
        {
            i = batch[b];
            if(!gravity_treewalk_candidate_prewalk(i, batch_pos[b])) {ProcessedFlag[i]=1; continue;}
            batch[n_candidates++] = i;
        }

#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
        if(Ewald_iter)
        {
            for(int b = 0; b < n_candidates; b++)
            {
                i = batch[b];
                int ret = force_treeevaluate_ewald_correction(i, exportflag, exportnodecount, exportindex);
                if(ret >= 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    Ewaldcount += ret;
                } else {buffer_full = 1; break;}
                ProcessedFlag[i] = 1;
            }
        }
        else
#endif
        {
            /* consecutive candidates walk the tree together, up to the packet capacity at a time; a
             * candidate the pre-pass did not count cannot exist, so a zero capacity means none */
            if(n_candidates > 0 && packet_cap <= 0) {endrun(90001055); break;}
            for(int start = 0; start < n_candidates && !buffer_full; start += packet_cap)
            {
                const int *packet = batch + start;
                int n_packet = n_candidates - start; if(n_packet > packet_cap) {n_packet = packet_cap;}
                int ret = force_treeevaluate(packet, n_packet, packet_cap, ninter, walk_ws, exportflag, exportnodecount, exportindex);
                if(ret > 0) {for(int m = 0; m < n_packet; m++) {commit_target(packet[m], ninter[m]);} continue;}
                if(ret < 0) {buffer_full = 1; break;}
                /* the packet met a pseudo-particle and wrote nothing: each member walks alone,
                 * which records the pseudo-particle for that member as a single walk always has */
                for(int m = 0; m < n_packet; m++)
                {
                    ret = force_treeevaluate(packet + m, 1, packet_cap, ninter, walk_ws, exportflag, exportnodecount, exportindex);
                    if(ret < 0) {buffer_full = 1; break;}
                    commit_target(packet[m], ninter[0]);
                }
            }
        }
        if(buffer_full) {break;}
    } // while loop
    return NULL;
}






/*! This function sets the (comoving) softening length of all particle types in the table All.ForceSoftening[...].
 We check that the physical softening length is bounded by the Softening-MaxPhys values */
void set_softenings(void)
{
    int i; double soft[6];
    soft[0] = All.SofteningGas;
    soft[1] = All.SofteningHalo;
    soft[2] = All.SofteningDisk;
    soft[3] = All.SofteningBulge;
    soft[4] = All.SofteningStars;
    soft[5] = All.SofteningBndry;
    if(All.ComovingIntegrationOn)
    {
        double soft_temp[6], cf_atime = 1./All.Time;
        soft_temp[0] = All.SofteningGasMaxPhys * cf_atime;
        soft_temp[1] = All.SofteningHaloMaxPhys * cf_atime;
        soft_temp[2] = All.SofteningDiskMaxPhys * cf_atime;
        soft_temp[3] = All.SofteningBulgeMaxPhys * cf_atime;
        soft_temp[4] = All.SofteningStarsMaxPhys * cf_atime;
        soft_temp[5] = All.SofteningBndryMaxPhys * cf_atime;
        for(i=0; i<6; i++) {if(soft_temp[i]<soft[i]) {soft[i]=soft_temp[i];}}
    }
    for(i=0; i<6; i++) {All.ForceSoftening[i] = soft[i] / KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER;}
    All.MinKernelRadius = All.MinGasKernelRadiusFractional * All.ForceSoftening[0]; /* set the minimum gas kernel length to be used this timestep */
#ifndef SELFGRAVITY_OFF
    if(All.MinKernelRadius <= 5.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0]) {All.MinKernelRadius = 5.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0];}
#endif
}


/* The DataIndexTable sorter (data_index_compare / mysort_dataindex) is retired with
 * the gravity export round-trip: the detector never ships or sorts its records. */


#if (SINGLE_STAR_TIMESTEPPING > 0)
void subtract_companion_gravity(int i)
{
    /* Remove contribution to gravitational field and tidal tensor from the stars in the binary to the center of mass */
    double u, dr, fac, fac2, h, h_inv, h3_inv, u2; SymmetricTensor2<MyFloat> tidal_tensorps; int i1, i2;
    dr = P[i].comp_dx.norm();
    h = SinkParticle_GravityKernelRadius;  h_inv = 1.0 / h; h3_inv = h_inv*h_inv*h_inv; u = dr*h_inv; u2=u*u;
    fac = P[i].comp_Mass / (dr*dr*dr); fac2 = 3.0 * P[i].comp_Mass / (dr*dr*dr*dr*dr); /* no softening nonsense */
    if(dr < h) /* second derivatives needed -> calculate them from softened potential */
    {
	    fac = P[i].comp_Mass * kernel_gravity(u, h_inv, h3_inv, 1);
        fac2 = P[i].comp_Mass * kernel_gravity(u, h_inv, h3_inv, 2);
    }
    P[i].COM_GravAccel = P[i].GravAccel - P[i].comp_dx * (fac * All.G); /* this assumes the 'G' has been put into the units for the grav accel */

    /* Adjusting tidal tensor according to terms above */
    tidal_tensorps = P[i].tidal_tensorps - fac2 * outer_product(P[i].comp_dx);
    tidal_tensorps[0][0] += fac; tidal_tensorps[1][1] += fac; tidal_tensorps[2][2] += fac;

#ifdef SINK_OUTPUT_MOREINFO
    printf("Corrected center of mass acceleration %g %g %g tidal tensor diagonal elements %g %g %g \n", P[i].COM_GravAccel[0], P[i].COM_GravAccel[1], P[i].COM_GravAccel[2], tidal_tensorps[0][0],tidal_tensorps[1][1],tidal_tensorps[2][2]);
#endif
    P[i].COM_dt_tidal = sqrt(1.0 / (All.G * tidal_tensorps.frobenius_norm()));
}
#endif

#ifdef ADAPTIVE_TREEFORCE_UPDATE
int needs_new_treeforce(int n){
    if(P[n].Type > 0){ // in this implementation we only do the lazy updating for gas cells whose timesteps are otherwise constrained by multiphysics (e.g. radiation, feedback)
        return 1;
    } else {
        if(P[n].time_since_last_treeforce >= P[n].tdyn_step_for_treeforce * ADAPTIVE_TREEFORCE_UPDATE) {return 1;}
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        else if(P[n].time_since_last_treeforce >= P[n].Min_Sink_FeedbackTime) {return 1;} // we want ejecta to re-calculate their feedback time so they don't get stuck on a short timestep
#endif        
        else {return 0;}
    }
}
#endif

/* SSOT pre-walk gravity tree-walk candidacy: true iff active particle i will
 * receive a real tree-force walk this step (and thus consume the installed LET).
 * Shared by the CPU primary, GPU primary, and GPU Ewald walk filters (and, later,
 * the LET-freshness check) so the freshness basis matches the actual LET consumer.
 * Mass>0 enforces the scheduler contract (Mass<=0 = scheduled-for-deletion, never a
 * valid gravity target) uniformly -- a defensive parity guard (the Ewald walk
 * already filtered it; the primary walks relied on the active-list builder and the
 * device early-return). ProcessedFlag is NOT part of candidacy -- it is per-walk
 * done-bookkeeping each caller keeps separately. */
int gravity_treewalk_candidate_prewalk(int i, int ii)
{
    if(P[i].Mass <= 0) {return 0;}
#ifdef HERMITE_INTEGRATION
    if(HermiteOnlyFlag && !eligible_for_hermite(i, P)) {return 0;}
#endif
#ifdef ADAPTIVE_TREEFORCE_UPDATE
    if(!gravity_treeforce_candidate_frozen(ii)) {return 0;}
#else
    (void) ii;
#endif
    return 1;
}
