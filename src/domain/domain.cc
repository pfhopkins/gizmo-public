#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../core/wakeup_sidecar.h"
#include "../system/gpu_particles_arena.h"
#include "../mesh/gpu_neighbor_list.h" /* gizmo_mark_kernel_radius_dirty_range */
#include "../mesh/topleaf_router.h"    /* topleaf_router_geometry_invalidate */


/*! Announce that the local particle array was re-indexed: particles were
 *  exchanged between ranks and the remaining ones compacted, so a given slot
 *  in P[] generally holds a different particle than it did before.
 *
 *  Everything cached against local slot identity or ordering is therefore
 *  invalid: the supply pools keyed on membership/order, the ghost-exchange
 *  local tree, the GPU spatial index (its tiles and pool reference slots by
 *  number), and the top-leaf router's node geometry. None of them can detect
 *  this on their own -- a redistribution can leave a rank's particle count
 *  unchanged, and it changes no particle's own properties, so counts and
 *  per-particle state both look untouched.
 *
 *  This lives here, at the point where the re-indexing happens, because
 *  callers cannot be relied on to follow every decomposition with the same
 *  bookkeeping: the initial decomposition, the restart path, FOF and SUBFIND
 *  all call it from outside the main loop, and a stale index there produced
 *  neighbour lists that silently omitted real neighbours.
 *
 *  Not needed for the GPU particle arena: it aliases P[] directly rather than
 *  holding a copy, so there is nothing in it to go stale. */
static void domain_particle_layout_changed(const char *reason)
{
    /* Record the event first, so a consumer that only compares the epoch sees
     * it regardless of what the cache-freeing calls below do. */
    ghost_exchange_supply_identity_changed(reason);
    ghost_exchange_local_tree_invalidate_full();
    gpu_step_sidx_invalidate_full();
    topleaf_router_geometry_invalidate();
}


/*! \file domain.c
 *  \brief code for domain decomposition
 *
 *  This file contains the code for the domain decomposition of the
 *  simulation volume.  The domains are constructed from disjoint subsets
 *  of the leaves of a fiducial top-level tree that covers the full
 *  simulation volume. Domain boundaries hence run along tree-node
 *  divisions of a fiducial global BH tree. As a result of this method, the
 *  tree force are in principle strictly independent of the way the domains
 *  are cut. The domain decomposition can be carried out for an arbitrary
 *  number of CPUs. Individual domains are not cubical, but spatially
 *  coherent since the leaves are traversed in a Peano-Hilbert order and
 *  individual domains form segments along this order.  This also ensures
 *  that each domain has a small surface to volume ratio, which minimizes
 *  communication.
 */


/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * significantly by Phil Hopkins (phopkins@caltech.edu) for GIZMO; these
 * modifications do not change the most basic algorithmic choices in the domain
 * decomposion, but have optimized it in
 * some places, changed relative weighting factors for different levels in the
 * domain decomposition, and similar details. Also how some memory issues are
 * handled has been updated to reflect the newer more general parallelization
 * structures in GIZMO. Changed buffer structures. Many variable name changes,
 * and changes in how cells vs non-cell elements are treated and weighted
 * and assigned to memory. Removed extensive non-functional code not relevant
 * for GIZMO.
 */


/*! toGo[task*NTask + partner] gives the number of particles in task 'task'
 *  that have to go to task 'partner'
 */
static int *toGo, *toGoGas;
static int *toGet, *toGetGas;
static int *list_NumPart;
static int *list_N_gas;
static int *list_load;
static int *list_loadgas;
static double *list_work;
static double *list_workgas;
extern int old_MaxPart, new_MaxPart;
#define N_DOMAINDECOMP_QUEUES 4

/* Adaptive domain balance weights: these are updated after each decomposition based
   on measured imbalance, following the GADGET-4 approach of dynamically adjusting
   the relative importance of work vs memory balance. */
static double domain_fac_work = 0.5;    /* weight for gravity work */
static double domain_fac_workgas = 0.0; /* weight for hydro work (set >0 when gas is present) */
static double domain_fac_load = 0.5;    /* weight for particle count (memory) */

static struct local_topnode_data
{
  peanokey Size;		/*!< number of Peano-Hilbert mesh-cells represented by top-level node */
  peanokey StartKey;		/*!< first Peano-Hilbert key in top-level node */
  long long Count;		/*!< counts the number of particles in this top-level node */
  double Cost;
  double GasCost;
  int Daughter;			/*!< index of first daughter cell (out of 8) of top-level node */
  int Leaf;			/*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  int Parent;
  int PIndex;			/*!< first particle in node */
}
 *topNodes;			/*!< points to the root node of the top-level tree */

static struct peano_hilbert_data
{
  peanokey key;
  int index;
}
 *mp;

static int domain_insertnode(struct local_topnode_data *treeA, struct local_topnode_data *treeB, int noA, int noB);
static void domain_add_cost(struct local_topnode_data *treeA, int noA, long long count, double cost, double gascost);

/*! Walk the top tree to find the leaf node for a given Peano-Hilbert key, using bitwise operations instead of 128-bit division */
static inline int domain_toptree_leaf(peanokey key, struct local_topnode_data *tNodes)
{
    int no = 0;
    peanokey mask = ((peanokey)7) << (3 * (BITS_PER_DIMENSION - 1));
    int shift = 3 * (BITS_PER_DIMENSION - 1);
    while(tNodes[no].Daughter >= 0)
    {
        no = tNodes[no].Daughter + (int)((key & mask) >> shift);
        mask >>= 3;
        shift -= 3;
    }
    return tNodes[no].Leaf;
}

static float *particle_total_cost;  /*!< cached per-particle total work cost: (1+multiplier)*costfactor */
static float *particle_costfactor;  /*!< cached per-particle base cost factor (for gas work accounting) */
static float *domainWork;	/*!< a table that gives the total "work" due to the particles stored by each processor */
static float *domainWorkGas;	/*!< a table that gives the total "work" due to the particles stored by each processor */
static int *domainCount;	/*!< a table that gives the total number of particles held by each processor */
static int *domainCountGas;	/*!< a table that gives the total number of gas cells held by each processor */
static int domain_allocated_flag = 0;
static int maxLoad, maxLoadgas;
static double totgravcost, gravcost, totgascost, gascost;
static long long totpartcount;
static int UseAllParticles;
static peanokey *PersistentKey = NULL; /*!< persistent Peano-Hilbert keys surviving between domain decompositions, used by lightweight repartition */
static int PersistentKeySize = 0;     /*!< allocated size of PersistentKey array */
static int LightRepartitionCount = 0; /*!< number of consecutive lightweight repartitions since last full decomposition */
#define MAX_LIGHT_REPARTITIONS 20     /*!< force a full domain decomposition after this many consecutive lightweight ones, to adapt top tree to changed particle distribution */

#if (DOMAIN_TIMEBINS == 1)
/* Per-timebin cost tracking for Gadget-4-style domain decomposition.
   Instead of a single composite cost, we track gravity and hydro costs separately
   for each occupied timebin, then balance all timebins simultaneously during
   domain assignment. */
static int NumTimeBinsToBeBalanced;
static int ListOfTimeBinsToBeBalanced[TIMEBINS];
static double GravCostPerListedTimeBin[TIMEBINS];
static double GravCostNormFactors[TIMEBINS];
static double HydroCostPerListedTimeBin[TIMEBINS];
static double HydroCostNormFactors[TIMEBINS];
static double NormFactorLoad, NormFactorLoadGas;
static float *domainBinGravCost = NULL;  /* [NumTimeBinsToBeBalanced * NTopleaves] per-timebin gravity cost per leaf */
static float *domainBinHydroCost = NULL; /* [NumTimeBinsToBeBalanced * NTopleaves] per-timebin hydro cost per leaf */

/*! Determine which timebins need to be individually balanced and compute their
 *  normalization factors. Follows Gadget-4's domain_init_sum_cost() and
 *  domain_find_total_cost(). */
void domain_init_timebin_costs(void)
{
    /* Determine which timebins to balance: all occupied timebins from
       HighestOccupiedTimeBin down to the lowest with particles */
    NumTimeBinsToBeBalanced = 0;
    long long tot_count[TIMEBINS], tot_count_gas[TIMEBINS];
    sumup_large_ints(TIMEBINS, TimeBinCount, tot_count);
    sumup_large_ints(TIMEBINS, TimeBinCountGas, tot_count_gas);

    /* Always include the highest active timebin */
    ListOfTimeBinsToBeBalanced[0] = All.HighestActiveTimeBin;
    NumTimeBinsToBeBalanced = 1;

    /* Add all lower timebins that have particles, with exponentially increasing weight */
    for(int i = All.HighestActiveTimeBin - 1; i >= 0; i--)
    {
        if(tot_count[i] > 0 || tot_count_gas[i] > 0)
        {
            ListOfTimeBinsToBeBalanced[NumTimeBinsToBeBalanced] = i;
            NumTimeBinsToBeBalanced++;
        }
    }

    /* Compute global per-timebin cost totals and normalization factors.
       A particle on timebin b contributes to all listed timebins with bin >= b. */
    for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
    {
        GravCostPerListedTimeBin[n] = 0;
        HydroCostPerListedTimeBin[n] = 0;
    }

    for(int i = 0; i < NumPart; i++)
    {
        double gc = (double)particle_total_cost[i];
        if(gc <= 0) gc = 1.0;
        for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
        {
            int bin = ListOfTimeBinsToBeBalanced[n];
            if(bin >= P[i].TimeBin)
                GravCostPerListedTimeBin[n] += gc;
            if(P[i].Type == 0 && bin >= P[i].TimeBin)
                HydroCostPerListedTimeBin[n] += 1.0;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, GravCostPerListedTimeBin, NumTimeBinsToBeBalanced, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, HydroCostPerListedTimeBin, NumTimeBinsToBeBalanced, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* Normalization: each timebin's total cost normalizes to ~1, so all timebins
       contribute equally to the composite cost used for tree splitting */
    for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
    {
        GravCostNormFactors[n] = (GravCostPerListedTimeBin[n] > 0) ? 1.0 / GravCostPerListedTimeBin[n] : 0.0;
        HydroCostNormFactors[n] = (HydroCostPerListedTimeBin[n] > 0) ? 1.0 / HydroCostPerListedTimeBin[n] : 0.0;
    }

    /* Load normalization */
    NormFactorLoad = (totpartcount > 0) ? 1.0 / (double)totpartcount : 0.0;
    long long totgas = 0; for(int i = 0; i < 6; i++) {if(i == 0) totgas = Ntype[0];}
    NormFactorLoadGas = (totgas > 0) ? 1.0 / (double)totgas : 0.0;

    if(ThisTask == 0)
    {
        printf("DOMAIN_TIMEBINS: balancing %d timebins:", NumTimeBinsToBeBalanced);
        for(int n = 0; n < NumTimeBinsToBeBalanced; n++)
            printf(" [bin=%d grav=%.3g hydro=%.3g]", ListOfTimeBinsToBeBalanced[n],
                   GravCostPerListedTimeBin[n], HydroCostPerListedTimeBin[n]);
        printf("\n");
    }
}
#endif /* DOMAIN_TIMEBINS == 1 */

/*! This is the main routine for the domain decomposition.  It acts as a driver routine that allocates various temporary buffers, maps the
 *  particles back onto the periodic box if needed, and then does the domain decomposition, and a final Peano-Hilbert order of all particles as a tuning measure. */
void domain_Decomposition(int UseAllTimeBins, int SaveKeys, int do_particle_mergesplit_key)
{
    int i, ret, retsum, diff, highest_bin_to_include; size_t bytes, all_bytes; double t0, t1;
    
    /* call first -before- a merge-split, to be sure particles are in the correct order in the tree */
    // TO: we don't have to call this before merge_and_split particles()
    // Actually we shouldn't because there are tree-walks in merge_and_split_particles().
    //rearrange_particle_sequence();
    const double child0_drift = CPU_ChildCharged;
    double t_drift_start = my_second(), t_mergesplit=0, t_rearrange=0, t_drift_loop=0, t_treefree=0, t_boxwrap=0, t_barrier=0;
    if((All.Ti_Current > All.TimeBegin)&&(do_particle_mergesplit_key==1))
    {
        merge_and_split_particles(); /* do the particle split/merge operations: only do this on tree-building super-steps */
    }
    t_mergesplit = timediff(t_drift_start, my_second());
    double t_tmp = my_second();
    rearrange_particle_sequence(); /* must be called after merge_and_split_particles, and should always be called before new domains are built */
    t_rearrange = timediff(t_tmp, my_second());

    UseAllParticles = UseAllTimeBins;

    t_tmp = my_second();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for(i = 0; i < NumPart; i++) {if(P[i].Ti_current != All.Ti_Current) {drift_particle(i, All.Ti_Current);}}
    gizmo_mark_kernel_radius_dirty_range(0, NumPart);
    t_drift_loop = timediff(t_tmp, my_second());

    t_tmp = my_second();
    force_treefree();
    domain_free();
    t_treefree = timediff(t_tmp, my_second());

    if(old_MaxPart) {All.MaxPart = new_MaxPart; old_MaxPart = 0;}

#ifdef BOX_PERIODIC
    t_tmp = my_second();
    do_box_wrapping();		/* map the particles back onto the box */
    t_boxwrap = timediff(t_tmp, my_second());
#endif

    t_tmp = my_second();
    MPI_Barrier(MPI_COMM_WORLD);
    t_barrier = timediff(t_tmp, my_second());
    const double t_drift_end = my_second();
    double t_drift_total = cpu_minus_children(timediff(t_drift_start, t_drift_end), child0_drift);
    CPU_Step[CPU_DRIFT] += t_drift_total;
    cpu_chain_sync(t_drift_end);
    if(ThisTask == 0) {
        printf("  domain_Decomp drift breakdown: mergesplit=%.4f rearrange=%.4f drift_loop=%.4f treefree=%.4f boxwrap=%.4f barrier=%.4f total=%.4f\n",
               t_mergesplit, t_rearrange, t_drift_loop, t_treefree, t_boxwrap, t_barrier, t_drift_total);
    }
    
    TreeReconstructFlag = 1;	/* ensures that new tree will be constructed */
#ifdef SINGLE_STAR_SINK_DYNAMICS
    All.NumForcesSinceLastDomainDecomp = 0;
#endif
    
    /* we take the closest cost factor */
    if(UseAllParticles) {highest_bin_to_include = All.HighestOccupiedTimeBin;} else {highest_bin_to_include = All.HighestActiveTimeBin;}
    
    for(i = 1, TakeLevel = 0, diff = abs(All.LevelToTimeBin[0] - highest_bin_to_include); i < GRAVCOSTLEVELS; i++)
        {if(diff > abs(All.LevelToTimeBin[i] - highest_bin_to_include)) {TakeLevel = i; diff = abs(All.LevelToTimeBin[i] - highest_bin_to_include);}}
    
    PRINT_STATUS("Domain decomposition building... LevelToTimeBin[TakeLevel=%d]=%d  (presently allocated=%g MB)", TakeLevel, All.LevelToTimeBin[TakeLevel], AllocatedBytes / (1024.0 * 1024.0));
    t0 = my_second();

    do
    {
      domain_allocate();

      all_bytes = 0;

      Key = (peanokey *) mymalloc("domain_key", bytes = (sizeof(peanokey) * All.MaxPart));
      all_bytes += bytes;

      toGo = (int *) mymalloc("toGo", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      toGoGas = (int *) mymalloc("toGoGas", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      toGet = (int *) mymalloc("toGet", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      toGetGas = (int *) mymalloc("toGetGas", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_NumPart = (int *) mymalloc("list_NumPart", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_N_gas = (int *) mymalloc("list_N_gas", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_load = (int *) mymalloc("list_load", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_loadgas = (int *) mymalloc("list_loadgas", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_work = (double *) mymalloc("list_work", bytes = (sizeof(double) * NTask));
      all_bytes += bytes;
      list_workgas = (double *) mymalloc("list_workgas", bytes = (sizeof(double) * NTask));
      all_bytes += bytes;
      domainWork = (float *) mymalloc("domainWork", bytes = (MaxTopNodes * sizeof(float)));
      all_bytes += bytes;
      domainWorkGas = (float *) mymalloc("domainWorkGas", bytes = (MaxTopNodes * sizeof(float)));
      all_bytes += bytes;
      domainCount = (int *) mymalloc("domainCount", bytes = (MaxTopNodes * sizeof(int)));
      all_bytes += bytes;
      domainCountGas = (int *) mymalloc("domainCountGas", bytes = (MaxTopNodes * sizeof(int)));
      all_bytes += bytes;

      topNodes = (struct local_topnode_data *) mymalloc("topNodes", bytes =
							(MaxTopNodes * sizeof(struct local_topnode_data)));
      all_bytes += bytes;

	  PRINT_STATUS(" ..using %g MB of temporary storage for domain decomposition... (presently allocated=%g MB)",all_bytes / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));

      maxLoad = (int) (All.MaxPart * REDUC_FAC_FOR_MEMORY_IN_DOMAIN);
      maxLoadgas = (int) (All.MaxPartGas * REDUC_FAC_FOR_MEMORY_IN_DOMAIN);
      /* Adaptive ghost headroom: reduce maxLoad by the previous ghost count (with safety margin)
         so domain decomposition leaves room for ghost particles on boundary ranks.
         Uses MPI_MAX of per-rank ghost counts as a conservative global estimate. */
      {
          int prev_ghost_local = ghost_get_previous_count();
          int prev_ghost_max;
          MPI_Allreduce(&prev_ghost_local, &prev_ghost_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
          if(prev_ghost_max > 0) {
              int ghost_headroom = (int)(prev_ghost_max * 1.3); /* 30% safety margin for h-growth */
              int maxLoad_new = All.MaxPart - ghost_headroom;
              if(maxLoad_new < All.MaxPart / 2) {maxLoad_new = All.MaxPart / 2;} /* never reserve more than half */
              if(maxLoad_new < maxLoad) {
                  if(ThisTask == 0) {PRINT_STATUS("Domain: ghost headroom %d (prev_max=%d * 1.3), maxLoad %d -> %d",
                                                   ghost_headroom, prev_ghost_max, maxLoad, maxLoad_new);}
                  maxLoad = maxLoad_new;
              }
          }
      }

      report_memory_usage(&HighMark_domain, "DOMAIN");

      ret = domain_decompose();
        
      /* copy what we need for the topnodes */
      for(i = 0; i < NTopnodes; i++)
      {
          TopNodes[i].StartKey = topNodes[i].StartKey;
          TopNodes[i].Size = topNodes[i].Size;
          TopNodes[i].Daughter = topNodes[i].Daughter;
          TopNodes[i].Leaf = topNodes[i].Leaf;
      }

      myfree(topNodes);

      myfree(domainCountGas);
      myfree(domainCount);
      myfree(domainWorkGas);
      myfree(domainWork);
      myfree(list_workgas);
      myfree(list_work);
      myfree(list_loadgas);
      myfree(list_load);
      myfree(list_N_gas);
      myfree(list_NumPart);
      myfree(toGetGas);
      myfree(toGet);
      myfree(toGoGas);
      myfree(toGo);

      MPI_Allreduce(&ret, &retsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(retsum)
      {
        myfree(Key);
        domain_free();

        if(ThisTask == 0) {printf("Increasing TopNodeAllocFactor=%g  ", All.TopNodeAllocFactor);}

        All.TopNodeAllocFactor *= 1.3;

        PRINT_STATUS("..new value=%g", All.TopNodeAllocFactor);
        double topnode_limit = 1000;
#if (NUMDIMS < 3)
        topnode_limit = 100000; /* 1D/2D problems need many more top nodes because the 3D Peano-Hilbert decomposition wastes refinement on degenerate dimensions */
#endif
        /* Symmetric (follows MPI_Allreduce(&ret,&retsum); All.TopNodeAllocFactor
         * grows identically on every rank). Key + domain are already freed above,
         * so a break would use-after-free below -- soft bad-stop + an immediate
         * all-rank poll drains here before any continuation on freed domain state. */
        if(All.TopNodeAllocFactor > topnode_limit) {printf("something seems to be going seriously wrong here. Stopping.\n"); fflush(stdout); endrun(90000010); gizmo_exit_bad_stop_if_requested("domain:topnode_alloc_factor");}
      }
    }
    while(retsum);

    t1 = my_second();

    PRINT_STATUS(" ..domain decomposition done. (took %g sec)", timediff(t0, t1));
    CPU_Step[CPU_DOMAIN] += measure_time();

    /* Per-rank particle-type validation. A bad type must not continue into
     * peano_hilbert_order / topnode compaction / tree alloc below, so soft
     * bad-stop on the first offender and poll AFTER the loop -- all ranks call
     * the poll (one domain-phase Allreduce on the clean path). */
    for(i = 0; i < NumPart; i++) {
        if(P[i].Type > 5 || P[i].Type < 0) {
            printf("task=%d:  P[i=%d].Type=%d\n", ThisTask, i, P[i].Type);
            endrun(90000011);
            break;
        }
    }
    gizmo_exit_bad_stop_if_requested("domain:particle_type_check");

#ifdef SUBFIND
    if(GrNr < 0)			/* we don't do it when SUBFIND is executed for a certain group */
#endif
    {peano_hilbert_order();}
    CPU_Step[CPU_PEANO] += measure_time();

  LightRepartitionCount = 0; /* reset counter: top tree is fresh */

  /* save keys persistently for potential lightweight repartition */
  if(PersistentKeySize < All.MaxPart) {
      if(PersistentKey) {free(PersistentKey);}
      PersistentKey = (peanokey *) malloc(All.MaxPart * sizeof(peanokey));
      PersistentKeySize = All.MaxPart;
  }
  memcpy(PersistentKey, Key, NumPart * sizeof(peanokey));

  myfree(Key);
  memmove(TopNodes + NTopnodes, DomainTask, NTopnodes * sizeof(int));
  TopNodes = (struct topnode_data *) myrealloc(TopNodes, bytes = (NTopnodes * sizeof(struct topnode_data) + NTopnodes * sizeof(int)));
  PRINT_STATUS(" ..freed %g MByte in top-level domain structure", (MaxTopNodes - NTopnodes) * sizeof(struct topnode_data) / (1024.0 * 1024.0));
  DomainTask = (int *) (TopNodes + NTopnodes);
  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
  gizmo_exit_bad_stop_if_requested("domain:treeallocate"); /* drain a tree-alloc UVM OOM (all-rank) before any tree use */
  reconstruct_timebins();
  gpu_particles_arena_invalidate(); /* P[] reordered across ranks; arena stale */
  wakeup_sidecar_invalidate();      /* P[] reindexed across ranks → rebuild WakeupDirty from P[] next scan */
  domain_particle_layout_changed("domain_Decomposition");
}


/*! Lightweight domain repartition: reuses the existing top-tree structure and Peano-Hilbert keys,
 *  only recomputing particle costs, re-splitting, and exchanging particles that changed domain.
 *  Skips: key computation, sorting, top-tree building/combining, and PH reorder.
 *  This is O(N) instead of O(N log N) and avoids expensive MPI tree combination. */
void domain_Decomposition_light(int UseAllTimeBins)
{
    int i, no; size_t bytes; double t0, t1;

    /* fall back to full decomposition if persistent state is not available, or if
       too many consecutive lightweight repartitions have occurred (top tree may be stale) */
    if(!PersistentKey || !domain_allocated_flag || LightRepartitionCount >= MAX_LIGHT_REPARTITIONS) {domain_Decomposition(UseAllTimeBins, 0, 1); return;}
    LightRepartitionCount++;

    const double child0_light = CPU_ChildCharged;
    double t_light_start = my_second(), t_light_rearrange=0, t_light_drift=0, t_light_boxwrap=0, t_light_barrier=0;
    rearrange_particle_sequence();
    t_light_rearrange = timediff(t_light_start, my_second());
    UseAllParticles = UseAllTimeBins;

    double t_tmp2 = my_second();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
    for(i = 0; i < NumPart; i++) {if(P[i].Ti_current != All.Ti_Current) {drift_particle(i, All.Ti_Current);}}
    gizmo_mark_kernel_radius_dirty_range(0, NumPart);
    t_light_drift = timediff(t_tmp2, my_second());

#ifdef BOX_PERIODIC
    t_tmp2 = my_second();
    do_box_wrapping();
    t_light_boxwrap = timediff(t_tmp2, my_second());
#endif

    t_tmp2 = my_second();
    MPI_Barrier(MPI_COMM_WORLD);
    t_light_barrier = timediff(t_tmp2, my_second());
    const double t_light_end = my_second();
    double t_light_total = cpu_minus_children(timediff(t_light_start, t_light_end), child0_light);
    CPU_Step[CPU_DRIFT] += t_light_total;
    cpu_chain_sync(t_light_end);
    if(ThisTask == 0) {
        printf("  domain_light drift breakdown: rearrange=%.4f drift_loop=%.4f boxwrap=%.4f barrier=%.4f total=%.4f\n",
               t_light_rearrange, t_light_drift, t_light_boxwrap, t_light_barrier, t_light_total);
    }

    /* we take the closest cost factor */
    int diff, highest_bin_to_include;
    if(UseAllParticles) {highest_bin_to_include = All.HighestOccupiedTimeBin;} else {highest_bin_to_include = All.HighestActiveTimeBin;}
    for(i = 1, TakeLevel = 0, diff = abs(All.LevelToTimeBin[0] - highest_bin_to_include); i < GRAVCOSTLEVELS; i++)
        {if(diff > abs(All.LevelToTimeBin[i] - highest_bin_to_include)) {TakeLevel = i; diff = abs(All.LevelToTimeBin[i] - highest_bin_to_include);}}

    PRINT_STATUS("Domain decomposition (lightweight)... LevelToTimeBin[TakeLevel=%d]=%d", TakeLevel, All.LevelToTimeBin[TakeLevel]);
    t0 = my_second();

    /* free force tree but keep domain structures (TopNodes, DomainTask) */
    force_treefree();

    TreeReconstructFlag = 1;
#ifdef SINGLE_STAR_SINK_DYNAMICS
    All.NumForcesSinceLastDomainDecomp = 0;
#endif

    int multipledomains = MULTIPLEDOMAINS;

    /* recompute keys for particles that may have drifted across cell boundaries */
    Key = PersistentKey; /* reuse persistent key storage directly */
    for(i = 0; i < NumPart; i++)
    {
        peano1D xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) / DomainLen) + 1.0);
        peano1D yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) / DomainLen) + 1.0);
        peano1D zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) / DomainLen) + 1.0);
        Key[i] = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
    }

    /* temporarily reconstruct the local topNodes from the persistent TopNodes for domain_sumCost.
       We need the local_topnode_data format with Daughter/Leaf/StartKey/Size fields. */
    topNodes = (struct local_topnode_data *) mymalloc("topNodes_light", NTopnodes * sizeof(struct local_topnode_data));
    for(i = 0; i < NTopnodes; i++)
    {
        topNodes[i].StartKey = TopNodes[i].StartKey;
        topNodes[i].Size = TopNodes[i].Size;
        topNodes[i].Daughter = TopNodes[i].Daughter;
        topNodes[i].Leaf = TopNodes[i].Leaf;
    }

    /* recompute per-particle costs */
    for(i = 0; i < 6; i++) {NtypeLocal[i] = 0;}
    particle_total_cost = (float *) mymalloc("particle_total_cost", NumPart * sizeof(float));
    particle_costfactor = (float *) mymalloc("particle_costfactor", NumPart * sizeof(float));
    for(i = 0, gravcost = gascost = 0; i < NumPart; i++)
    {
        NtypeLocal[P[i].Type]++;
        double wt_0 = domain_particle_costfactor(i);
        double wt_mult = domain_particle_cost_multiplier(i);
        particle_costfactor[i] = (float)wt_0;
        particle_total_cost[i] = (float)((1 + wt_mult) * wt_0);
        gravcost += particle_total_cost[i];
        if(P[i].Type == 0) {if(TimeBinActive[P[i].TimeBin] || UseAllParticles) {gascost += wt_0;}}
    }
    sumup_large_ints(6, NtypeLocal, Ntype);
    for(i = 0, totpartcount = 0; i < 6; i++) {totpartcount += Ntype[i];}
    MPI_Allreduce(&gravcost, &totgravcost, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gascost, &totgascost, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#if (DOMAIN_TIMEBINS == 1)
    domain_init_timebin_costs();
#endif

    /* allocate work/count arrays */
    domainWork = (float *) mymalloc("domainWork", NTopnodes * sizeof(float));
    domainWorkGas = (float *) mymalloc("domainWorkGas", NTopnodes * sizeof(float));
    domainCount = (int *) mymalloc("domainCount", NTopnodes * sizeof(int));
    domainCountGas = (int *) mymalloc("domainCountGas", NTopnodes * sizeof(int));

    /* recompute costs per top-tree leaf */
    domain_sumCost();

    /* free cost caches (LIFO: domainCountGas is on top, but we need them for the split below,
       so free cost arrays first since they sit below domainWork on the stack) */
    /* Actually cost arrays are above domainWork... let's just keep them all alive until the end. */

    /* allocate list arrays for the split/assignment */
    toGo = (int *) mymalloc("toGo", NTask * sizeof(int));
    toGoGas = (int *) mymalloc("toGoGas", NTask * sizeof(int));
    toGet = (int *) mymalloc("toGet", NTask * sizeof(int));
    toGetGas = (int *) mymalloc("toGetGas", NTask * sizeof(int));
    list_NumPart = (int *) mymalloc("list_NumPart", NTask * sizeof(int));
    list_N_gas = (int *) mymalloc("list_N_gas", NTask * sizeof(int));
    list_load = (int *) mymalloc("list_load", NTask * sizeof(int));
    list_loadgas = (int *) mymalloc("list_loadgas", NTask * sizeof(int));
    list_work = (double *) mymalloc("list_work", NTask * sizeof(double));
    list_workgas = (double *) mymalloc("list_workgas", NTask * sizeof(double));

    maxLoad = (int) (All.MaxPart * REDUC_FAC_FOR_MEMORY_IN_DOMAIN);
    maxLoadgas = (int) (All.MaxPartGas * REDUC_FAC_FOR_MEMORY_IN_DOMAIN);
    {   /* adaptive ghost headroom (same logic as full decomposition above) */
        int prev_ghost_local = ghost_get_previous_count();
        int prev_ghost_max;
        MPI_Allreduce(&prev_ghost_local, &prev_ghost_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(prev_ghost_max > 0) {
            int ghost_headroom = (int)(prev_ghost_max * 1.3);
            int maxLoad_new = All.MaxPart - ghost_headroom;
            if(maxLoad_new < All.MaxPart / 2) {maxLoad_new = All.MaxPart / 2;}
            if(maxLoad_new < maxLoad) {maxLoad = maxLoad_new;}
        }
    }

    /* re-split and re-assign */
    domain_findSplit_work_balanced(multipledomains * NTask, NTopleaves);
    domain_assign_load_or_work_balanced(1, multipledomains);

#if (DOMAIN_TIMEBINS == 1)
    if(domainBinGravCost) {free(domainBinGravCost); domainBinGravCost = NULL;}
    if(domainBinHydroCost) {free(domainBinHydroCost); domainBinHydroCost = NULL;}
#endif

    int status = domain_check_memory_bound(multipledomains);
    if(status != 0)
    {
        domain_findSplit_load_balanced(multipledomains * NTask, NTopleaves);
        domain_assign_load_or_work_balanced(0, multipledomains);
        status = domain_check_memory_bound(multipledomains);
        if(status != 0) {if(ThisTask == 0) {printf("Lightweight repartition: memory bound violated.\n");}}
    }

    /* flag particles that need to move */
    for(i = 0; i < NumPart; i++)
    {
        no = domain_toptree_leaf(Key[i], topNodes);
        int task = DomainTask[no];
        if(task != ThisTask) {P[i].Type |= 32;}
    }

    /* exchange particles */
    int ret; size_t exchange_limit;
    const size_t exchange_overhead = NTask * (24 * sizeof(int) + 16 * sizeof(MPI_Request));
    const size_t exchange_min_package = sizeof(struct particle_data) + sizeof(struct gas_cell_data) + sizeof(peanokey);
    do
    {
        /* arena-capacity guard. FreeBytes is size_t, so the old `FreeBytes - overhead`
         * underflowed (never <=0) when the arena was too small -- a dead check that let
         * exchange-OOM fall through to the rank-divergent allocator. Guard explicitly and
         * request a graceful symmetric stop (drained at the all-rank poll below). */
        if(FreeBytes <= exchange_overhead + exchange_min_package) {
            printf("task=%d: domain exchange arena exhausted (FreeBytes=%g MB)\n", ThisTask, FreeBytes / (1024.0 * 1024.0));
            endrun(90000012);
            exchange_limit = exchange_min_package;   /* avoid size_t underflow; poll below exits before use */
        } else {
            exchange_limit = FreeBytes - exchange_overhead;
        }
        ret = domain_countToGo(exchange_limit);
        /* predicted post-exchange occupancy from the toGo/toGet just computed -- guard
         * BEFORE domain_exchange() receives into P[]/CellP[] so an overflow never overruns
         * the arrays (per-rank condition, drained at the all-rank poll). */
        {
            int n; long long ng_out = 0, ng_in = 0, gas_out = 0, gas_in = 0;
            for(n = 0; n < NTask; n++) {ng_out += toGo[n]; ng_in += toGet[n]; gas_out += toGoGas[n]; gas_in += toGetGas[n];}
            if((long long)NumPart - ng_out + ng_in > All.MaxPart)    {endrun(90000015);}
            if((long long)N_gas   - gas_out + gas_in > All.MaxPartGas) {endrun(90000016);}
        }
        gizmo_exit_bad_stop_if_requested("domain:exchange_capacity");
        domain_exchange();
        gizmo_exit_bad_stop_if_requested("domain:exchange");
    }
    while(ret > 0);

    /* free everything in LIFO order */
    myfree(list_workgas); myfree(list_work); myfree(list_loadgas); myfree(list_load);
    myfree(list_N_gas); myfree(list_NumPart);
    myfree(toGetGas); myfree(toGet); myfree(toGoGas); myfree(toGo);
    myfree(domainCountGas); myfree(domainCount); myfree(domainWorkGas); myfree(domainWork);
    myfree(particle_costfactor); myfree(particle_total_cost);
    myfree(topNodes);

    t1 = my_second();
    PRINT_STATUS(" ..lightweight domain repartition done. (took %g sec)", timediff(t0, t1));
    CPU_Step[CPU_DOMAIN] += measure_time();

    /* update persistent keys for moved particles */
    for(i = 0; i < NumPart; i++)
    {
        peano1D xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) / DomainLen) + 1.0);
        peano1D yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) / DomainLen) + 1.0);
        peano1D zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) / DomainLen) + 1.0);
        PersistentKey[i] = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
    }
    Key = NULL; /* no longer valid as a mymalloc pointer */

    force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
    gizmo_exit_bad_stop_if_requested("domain:treeallocate_light"); /* drain a tree-alloc UVM OOM (all-rank) before any tree use */
    reconstruct_timebins();
    wakeup_sidecar_invalidate();   /* light repartition rearranged + exchanged particles → rebuild WakeupDirty next scan */
    domain_particle_layout_changed("domain_Decomposition_light");
}


/*! This function allocates all the stuff that will be required for the tree-construction/walk later on */
void domain_allocate(void)
{
  size_t bytes, all_bytes = 0;

  MaxTopNodes = (int) (All.TopNodeAllocFactor * All.MaxPart + 1);

  DomainStartList = (int *) mymalloc("DomainStartList", bytes = (NTask * MULTIPLEDOMAINS * sizeof(int)));
  all_bytes += bytes;

  DomainEndList = (int *) mymalloc("DomainEndList", bytes = (NTask * MULTIPLEDOMAINS * sizeof(int)));
  all_bytes += bytes;

  TopNodes = (struct topnode_data *) mymalloc("TopNodes", bytes = (MaxTopNodes * sizeof(struct topnode_data) + MaxTopNodes * sizeof(int)));
  all_bytes += bytes;

  DomainTask = (int *) (TopNodes + MaxTopNodes);

  PRINT_STATUS(" ..allocated %g MByte for top-level domain structure", all_bytes / (1024.0 * 1024.0));

  domain_allocated_flag = 1;
}

void domain_free(void)
{
  if(domain_allocated_flag)
    {
      myfree(TopNodes);
      myfree(DomainEndList);
      myfree(DomainStartList);
      domain_allocated_flag = 0;
    }
}

static struct topnode_data *save_TopNodes;
static int *save_DomainStartList, *save_DomainEndList;

void domain_free_trick(void)
{
  if(domain_allocated_flag)
    {
      save_TopNodes = TopNodes;
      save_DomainEndList = DomainEndList;
      save_DomainStartList = DomainStartList;
      domain_allocated_flag = 0;
    }
  else
    {endrun(90000013); return;} /* not allocated: soft bad-stop + return (void fn, MPI-free, symmetric invariant); drains at the next domain poll */
}

void domain_allocate_trick(void)
{
  domain_allocated_flag = 1;
  TopNodes = save_TopNodes;
  DomainEndList = save_DomainEndList;
  DomainStartList = save_DomainStartList;
}


/* this function determines how particle work-costs are 'weighted' for load-balancing. if you 
    have additional, expensive physics which only apply to a subset of particles, it may be worth 
    up-weighting those particles here, so the code knows to try and spread them around. otherwise, 
    they may end up all bunched onto the same processor */
double domain_particle_cost_multiplier(int i)
{
    double multiplier = 0;
    
    if(P[i].Type == 0) /* for gas, weight particles with large neighbor number more, since they require more work */
    {
        double nngb_reduced = P[i].NumNgb; /* remember, in density.c we reduce this by pow(1/NUMDIMS), for use in other routines: need to correct here */
#if (NUMDIMS==3)
        multiplier = nngb_reduced*nngb_reduced*nngb_reduced / All.DesNumNgb;
#elif (NUMDIMS==2)
        multiplier = nngb_reduced*nngb_reduced / All.DesNumNgb;
#else
        multiplier = nngb_reduced / All.DesNumNgb;
#endif
        if(multiplier < 0.5) {multiplier = 0.5;} // floor //
    } // end gas check

#if defined(GALSF) /* with star formation active, we will up-weight star particles which are active feedback sources */
#ifndef CHIMES /* With CHIMES, the chemistry dominates the cost, so we boost (dense) gas but not stars. */
    if(((P[i].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[i].Type == 2)||(P[i].Type==3))))&&(P[i].Mass>0))
    {
        double star_age = evaluate_stellar_age_Gyr(i);
        if(star_age>0.1) {multiplier = 3.125;} else {if(star_age>0.035) {multiplier = 5.;} else {multiplier = 10.;}}
    }
#endif 
#endif

#ifdef CHIMES 
    /* With CHIMES, cost is dominated by the chemistry, particularly in dense gas. We therefore boost the cost factor of gas particles with nH >~ 1 cm^-3. */
    if(P[i].Type == 0) {double nH_cgs = CellP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_NHCGS; if(nH_cgs > 1) {multiplier = 10.0;}}
#endif
    
#ifdef CRFLUID_EVOLVE_SPECTRUM // again, cost totally dominated by dense gas here, this helps significantly
    if(P[i].Type == 0) {double nH_cgs = CellP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_NHCGS; if(nH_cgs > 1) {multiplier *= 100.;} else {multiplier *= 10.;}}
#endif
    
    return multiplier;
}


/* simple function to return costfactor for pure gravity calculation: based just on gravcost calculation, with constant for safety */
double domain_particle_costfactor(int i)
{
    return 0.1 + P[i].GravCost[TakeLevel];
}



/*! Update the adaptive domain balance weights based on measured imbalance from the
 *  current decomposition. Weights are shifted toward the metric with the worst imbalance,
 *  following the GADGET-4 approach. Uses exponential moving average with smoothing factor alpha. */
void domain_update_adaptive_weights(void)
{
    double alpha = 0.3; /* smoothing factor: 0 = keep old weights, 1 = fully replace */
    double imbal_work = 0, imbal_load = 0, imbal_workgas = 0;
    double sumwork = 0, sumload = 0, sumworkgas = 0;
    double maxwork = 0, maxloadd = 0, maxworkgas = 0;

    for(int i = 0; i < NTask; i++) {
        sumwork += list_work[i];
        sumload += list_load[i];
        sumworkgas += list_workgas[i];
        if(list_work[i] > maxwork) maxwork = list_work[i];
        if(list_load[i] > maxloadd) maxloadd = list_load[i];
        if(list_workgas[i] > maxworkgas) maxworkgas = list_workgas[i];
    }

    /* imbalance = max/avg - 1: 0 = perfect, higher = worse */
    if(sumwork > 0) imbal_work = maxwork / (sumwork / NTask) - 1.0;
    if(sumload > 0) imbal_load = maxloadd / (sumload / NTask) - 1.0;
    if(sumworkgas > 0) imbal_workgas = maxworkgas / (sumworkgas / NTask) - 1.0;

    /* target weights proportional to measured imbalance (shift attention to worst-balanced metric) */
    double eps = 0.05; /* minimum weight floor to prevent any metric from being totally ignored */
    double target_work = DMAX(imbal_work, eps);
    double target_load = DMAX(imbal_load, eps);
    double target_workgas = (sumworkgas > 0) ? DMAX(imbal_workgas, eps) : 0.0;

    double target_total = target_work + target_load + target_workgas;
    target_work /= target_total;
    target_load /= target_total;
    target_workgas /= target_total;

    /* exponential moving average update */
    domain_fac_work = (1.0 - alpha) * domain_fac_work + alpha * target_work;
    domain_fac_load = (1.0 - alpha) * domain_fac_load + alpha * target_load;
    domain_fac_workgas = (1.0 - alpha) * domain_fac_workgas + alpha * target_workgas;

    /* ensure minimum floor */
    if(domain_fac_work < eps) domain_fac_work = eps;
    if(domain_fac_load < eps) domain_fac_load = eps;

    if(ThisTask == 0) {
        printf("Domain adaptive weights: work=%.3f  load=%.3f  workgas=%.3f  (imbalance: work=%.3f load=%.3f gas=%.3f)\n",
               domain_fac_work, domain_fac_load, domain_fac_workgas, imbal_work, imbal_load, imbal_workgas);
    }
}


/*! This function carries out the actual domain decomposition for all
 *  particle types. It will try to balance the work-load for each domain,
 *  as estimated based on the P[i]-GravCost values.  The decomposition will
 *  respect the maximum allowed memory-imbalance given by the value of
 *  PartAllocFactor.
 */
int domain_decompose(void)
{
    int i, no, status;
    long long sumtogo, sumload, sumloadgas;
    int maxload, maxloadgas, multipledomains = MULTIPLEDOMAINS;
    double sumwork, maxwork, sumworkgas, maxworkgas;

    for(i = 0; i < 6; i++) {NtypeLocal[i] = 0;}

    /* compute and cache per-particle costs once, reused in domain_check_for_local_refine and domain_sumCost.
       particle_total_cost: unweighted cost used for top-tree refinement (preserves spatial structure).
       For domain_sumCost, we apply timestep frequency weighting there to balance sub-step work
       without distorting the tree refinement. */
    particle_total_cost = (float *) mymalloc("particle_total_cost", NumPart * sizeof(float));
    particle_costfactor = (float *) mymalloc("particle_costfactor", NumPart * sizeof(float));
    for(i = 0, gravcost = gascost = 0; i < NumPart; i++)
    {
#ifdef SUBFIND
        if(GrNr >= 0 && P[i].GrNr != GrNr) {continue;}
#endif
        NtypeLocal[P[i].Type]++;
        double wt_0 = domain_particle_costfactor(i);
        double wt_mult = domain_particle_cost_multiplier(i);
        particle_costfactor[i] = (float)wt_0;
        particle_total_cost[i] = (float)((1 + wt_mult) * wt_0);
        gravcost += particle_total_cost[i];
        if(P[i].Type == 0) {if(TimeBinActive[P[i].TimeBin] || UseAllParticles) {gascost += wt_0;}}
    }
    /* because Ntype[] is of type `long long', we cannot do a simple MPI_Allreduce() to sum the total particle numbers */
    sumup_large_ints(6, NtypeLocal, Ntype);

    for(i = 0, totpartcount = 0; i < 6; i++) {totpartcount += Ntype[i];}

    MPI_Allreduce(&gravcost, &totgravcost, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gascost, &totgascost, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#if (DOMAIN_TIMEBINS == 1)
    domain_init_timebin_costs(); /* compute per-timebin normalization before tree building / cost accumulation */
#endif

    /* determine global dimensions of domain grid */
    domain_findExtent();
    if(domain_determineTopTree()) {myfree(particle_costfactor); myfree(particle_total_cost); return 1;}
    myfree(particle_costfactor); myfree(particle_total_cost);

    /* find the split of the domain grid */
    domain_findSplit_work_balanced(multipledomains * NTask, NTopleaves);
    domain_assign_load_or_work_balanced(1,multipledomains);

#if (DOMAIN_TIMEBINS == 1)
    free(domainBinHydroCost); free(domainBinGravCost);
    domainBinHydroCost = domainBinGravCost = NULL;
#endif

    status = domain_check_memory_bound(multipledomains);

    if(status != 0)		/* the optimum balanced solution violates memory constraint, let's try something different */
    {
      if(ThisTask == 0) {printf("Note: the domain decomposition is suboptimum because the ceiling for memory-imbalance is reached\n");}

      domain_findSplit_load_balanced(multipledomains * NTask, NTopleaves);
      domain_assign_load_or_work_balanced(0,multipledomains);
      status = domain_check_memory_bound(multipledomains);

      if(status != 0)
      {
          if(ThisTask == 0) {printf("No domain decomposition that stays within memory bounds is possible.\n");}
          endrun(0);
      }
    }

    if(ThisTask == 0)
    {
        sumload = maxload = sumloadgas = maxloadgas = 0;
        sumwork = sumworkgas = maxwork = maxworkgas = 0;

        for(i = 0; i < NTask; i++)
        {
            sumload += list_load[i];
            sumloadgas += list_loadgas[i];
            sumwork += list_work[i];
            sumworkgas += list_workgas[i];

            if(list_load[i] > maxload) {maxload = list_load[i];}
            if(list_loadgas[i] > maxloadgas) {maxloadgas = list_loadgas[i];}
            if(list_work[i] > maxwork) {maxwork = list_work[i];}
            if(list_workgas[i] > maxworkgas) {maxworkgas = list_workgas[i];}
        }

        printf("Balance: gravity work-load balance=%g   memory-balance=%g   hydro work-load balance=%g\n",
	     maxwork / (sumwork / NTask), maxload / (((double) sumload) / NTask), maxworkgas / ((sumworkgas + 1.0e-30) / NTask));
    }

    /* update adaptive weights for next decomposition based on current imbalance */
    if(ThisTask == 0) {domain_update_adaptive_weights();}
    MPI_Bcast(&domain_fac_work, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&domain_fac_load, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&domain_fac_workgas, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* flag the particles that need to be exported */

    for(i = 0; i < NumPart; i++)
    {
#ifdef SUBFIND
      if(GrNr >= 0 && P[i].GrNr != GrNr) {continue;}
#endif

      no = domain_toptree_leaf(Key[i], topNodes);
      int task = DomainTask[no];
      if(task != ThisTask) {P[i].Type |= 32;}
    }

    int iter = 0, ret;
    size_t exchange_limit;
    const size_t exchange_overhead = NTask * (24 * sizeof(int) + 16 * sizeof(MPI_Request));
    const size_t exchange_min_package = sizeof(struct particle_data) + sizeof(struct gas_cell_data) + sizeof(peanokey);

    do
    {
        /* arena-capacity guard (see domain_Decomposition_light): the old size_t
         * `FreeBytes - overhead <= 0` underflowed and never fired. Guard explicitly and
         * request a graceful symmetric stop, drained at the all-rank poll below. */
        if(FreeBytes <= exchange_overhead + exchange_min_package) {
            printf("task=%d: domain exchange arena exhausted (FreeBytes=%g MB)\n", ThisTask, FreeBytes / (1024.0 * 1024.0));
            endrun(90000014);
            exchange_limit = exchange_min_package;   /* avoid size_t underflow; poll below exits before use */
        } else {
            exchange_limit = FreeBytes - exchange_overhead;
        }

        /* determine for each cpu how many particles have to be shifted to other cpus */
        ret = domain_countToGo(exchange_limit);

        for(i = 0, sumtogo = 0; i < NTask; i++) {sumtogo += toGo[i];}

        sumup_longs(1, &sumtogo, &sumtogo);

        PRINT_STATUS(" ..iter=%d exchange of %d%09d particles (ret=%d)", iter, (int) (sumtogo / 1000000000), (int) (sumtogo % 1000000000), ret);

        /* predicted post-exchange occupancy -- guard BEFORE domain_exchange() receives
         * into P[]/CellP[] so an overflow never overruns the arrays (drained at the poll). */
        {
            int n; long long ng_out = 0, ng_in = 0, gas_out = 0, gas_in = 0;
            for(n = 0; n < NTask; n++) {ng_out += toGo[n]; ng_in += toGet[n]; gas_out += toGoGas[n]; gas_in += toGetGas[n];}
            if((long long)NumPart - ng_out + ng_in > All.MaxPart)    {endrun(90000015);}
            if((long long)N_gas   - gas_out + gas_in > All.MaxPartGas) {endrun(90000016);}
        }
        gizmo_exit_bad_stop_if_requested("domain:exchange_capacity");

        domain_exchange();
        gizmo_exit_bad_stop_if_requested("domain:exchange");

        iter++;
    }
    while(ret > 0);

    return 0;
}






int domain_check_memory_bound(int multipledomains)
{
  int ta, m, i;
  int load, gasload, max_load, max_gasload;
  double work, workgas;

  max_load = max_gasload = 0;

  for(ta = 0; ta < NTask; ta++)
    {
      load = gasload = 0;
      work = workgas = 0;

      for(m = 0; m < multipledomains; m++)
      for(i = DomainStartList[ta * multipledomains + m]; i <= DomainEndList[ta * multipledomains + m]; i++)
	  {
	    load += domainCount[i];
	    gasload += domainCountGas[i];
	    work += domainWork[i];
	    workgas += domainWorkGas[i];
	  }

      list_load[ta] = load;
      list_loadgas[ta] = gasload;
      list_work[ta] = work;
      list_workgas[ta] = workgas;

      if(load > max_load) {max_load = load;}
      if(gasload > max_gasload) {max_gasload = gasload;}
    }

#ifdef SUBFIND
  if(GrNr >= 0)
    {
      load = max_load;
      gasload = max_gasload;

      for(i = 0; i < NumPart; i++)
	  {
	    if(P[i].GrNr != GrNr)
	    {
	      load++;
	      if(P[i].Type == 0) {gasload++;}
	    }
	  }
      MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&gasload, &max_gasload, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
#endif

    if(max_load > maxLoad)
    {
      if(ThisTask == 0) {printf("desired memory imbalance=%g  (limit=%d, needed=%d)\n", (max_load * All.PartAllocFactor) / maxLoad, maxLoad, max_load);}
      return 1;
    }

    if(max_gasload > maxLoadgas)
    {
      if(ThisTask == 0) {printf("desired memory imbalance=%g  (GAS/FLUID) (limit=%d, needed=%d)\n", (max_gasload * All.PartAllocFactor) / maxLoadgas, maxLoadgas, max_gasload);}
      return 1;
    }

  return 0;
}


void domain_exchange(void)
{
  long count_togo = 0, count_togo_gas = 0, count_get = 0, count_get_gas = 0;
  long *count, *count_gas, *offset, *offset_gas;
  long *count_recv, *count_recv_gas, *offset_recv, *offset_recv_gas;
  long i, n, ngrp, no, target;
  struct particle_data *partBuf;
  struct gas_cell_data *gasBuf;
  peanokey *keyBuf;

  count = (long *) mymalloc("count", NTask * sizeof(long));
  count_gas = (long *) mymalloc("count_gas", NTask * sizeof(long));
  offset = (long *) mymalloc("offset", NTask * sizeof(long));
  offset_gas = (long *) mymalloc("offset_gas", NTask * sizeof(long));

  count_recv = (long *) mymalloc("count_recv", NTask * sizeof(long));
  count_recv_gas = (long *) mymalloc("count_recv_gas", NTask * sizeof(long));
  offset_recv = (long *) mymalloc("offset_recv", NTask * sizeof(long));
  offset_recv_gas = (long *) mymalloc("offset_recv_gas", NTask * sizeof(long));


  long prec_offset, prec_count;
  long *decrease;

  decrease = (long *) mymalloc("decrease", NTask * sizeof(long));

  for(i = 1, offset_gas[0] = 0, decrease[0] = 0; i < NTask; i++)
    {
      offset_gas[i] = offset_gas[i - 1] + toGoGas[i - 1];
      decrease[i] = toGoGas[i - 1];
    }

  prec_offset = offset_gas[NTask - 1] + toGoGas[NTask - 1];


  offset[0] = prec_offset;
  for(i = 1; i < NTask; i++) {offset[i] = offset[i - 1] + (toGo[i - 1] - decrease[i]);}

  myfree(decrease);

  for(i = 0; i < NTask; i++)
    {
      count_togo += toGo[i];
      count_togo_gas += toGoGas[i];

      count_get += toGet[i];
      count_get_gas += toGetGas[i];

    }

  partBuf = (struct particle_data *) mymalloc("partBuf", count_togo * sizeof(struct particle_data));
  gasBuf = (struct gas_cell_data *) mymalloc("gasBuf", count_togo_gas * sizeof(struct gas_cell_data));
#ifdef CHIMES 
  struct gasVariables *gasChimesBuf;
  ChimesFloat *gasAbundancesBuf, *gasAbundancesRecvBuf, *tempAbundanceArray;
  int abunIndex; 
  gasChimesBuf = (struct gasVariables *) mymalloc("chiBuf", count_togo_gas * sizeof(struct gasVariables));
  gasAbundancesBuf = (ChimesFloat *) mymalloc("abunBuf", count_togo_gas * ChimesGlobalVars.totalNumberOfSpecies * sizeof(ChimesFloat));
  gasAbundancesRecvBuf = (ChimesFloat *) mymalloc("xRecBuf", count_get_gas * ChimesGlobalVars.totalNumberOfSpecies * sizeof(ChimesFloat));
  tempAbundanceArray = (ChimesFloat *) malloc(ChimesGlobalVars.totalNumberOfSpecies * sizeof(ChimesFloat));
#endif
  keyBuf = (peanokey *) mymalloc("keyBuf", count_togo * sizeof(peanokey));

  for(i = 0; i < NTask; i++) {count[i] = count_gas[i] = 0;}


  for(n = 0; n < NumPart; n++)
    {
      if((P[n].Type & (32 + 16)) == (32 + 16)) /* flagged with both 16 and 32 */
	{
	  P[n].Type &= 15; /* clear 16 and 32 */

	  no = domain_toptree_leaf(Key[n], topNodes);

	  target = DomainTask[no];

	  if(P[n].Type == 0)
	    {
	      partBuf[offset_gas[target] + count_gas[target]] = P[n];
	      keyBuf[offset_gas[target] + count_gas[target]] = Key[n];
#ifdef CHIMES
	      for(i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) {gasAbundancesBuf[((offset_gas[target] + count_gas[target]) * ChimesGlobalVars.totalNumberOfSpecies) + i] = ChimesGasVars[n].abundances[i];}
	      free_gas_abundances_memory(&(ChimesGasVars[n]), &ChimesGlobalVars);
	      ChimesGasVars[n].abundances = NULL;
	      ChimesGasVars[n].isotropic_photon_density = NULL;
	      ChimesGasVars[n].G0_parameter = NULL;
	      ChimesGasVars[n].H2_dissocJ = NULL;
	      gasChimesBuf[offset_gas[target] + count_gas[target]] = ChimesGasVars[n];
#endif
	      gasBuf[offset_gas[target] + count_gas[target]] = CellP[n];
	      count_gas[target]++;
	    }
	  else
	    {
	      partBuf[offset[target] + count[target]] = P[n];
	      keyBuf[offset[target] + count[target]] = Key[n];
	      count[target]++;
	    }


	  if(P[n].Type == 0)
	    {
	      P[n] = P[N_gas - 1];
	      CellP[n] = CellP[N_gas - 1];
	      Key[n] = Key[N_gas - 1];

#ifdef CHIMES
	      if (n < N_gas - 1)
		{
		  for(abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
		    {tempAbundanceArray[abunIndex] = ChimesGasVars[N_gas - 1].abundances[abunIndex];}
		  free_gas_abundances_memory(&(ChimesGasVars[N_gas - 1]), &ChimesGlobalVars);
		  ChimesGasVars[N_gas - 1].abundances = NULL;
		  ChimesGasVars[N_gas - 1].isotropic_photon_density = NULL;
		  ChimesGasVars[N_gas - 1].G0_parameter = NULL;
		  ChimesGasVars[N_gas - 1].H2_dissocJ = NULL;
		  ChimesGasVars[n] = ChimesGasVars[N_gas - 1];
		  allocate_gas_abundances_memory(&(ChimesGasVars[n]), &ChimesGlobalVars);
		  for (abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
		    {ChimesGasVars[n].abundances[abunIndex] = tempAbundanceArray[abunIndex];}
		}
#endif

	      P[N_gas - 1] = P[NumPart - 1];
	      Key[N_gas - 1] = Key[NumPart - 1];

	      NumPart--;
	      N_gas--;
	      n--;
	    }
	  else
	    {
	      P[n] = P[NumPart - 1];
	      Key[n] = Key[NumPart - 1];
	      NumPart--;
	      n--;
	    }
	}
    }

#ifdef CHIMES 
  free(tempAbundanceArray); 
#endif 

  long count_totget;

  count_totget = count_get_gas;

  if(count_totget)
    {
      memmove(P + N_gas + count_totget, P + N_gas, (NumPart - N_gas) * sizeof(struct particle_data));
      memmove(Key + N_gas + count_totget, Key + N_gas, (NumPart - N_gas) * sizeof(peanokey));
    }


  for(i = 0; i < NTask; i++)
    {
      count_recv_gas[i] = toGetGas[i];
      count_recv[i] = toGet[i] - toGetGas[i];
    }


  for(i = 1, offset_recv_gas[0] = N_gas; i < NTask; i++)
    {offset_recv_gas[i] = offset_recv_gas[i - 1] + count_recv_gas[i - 1];}
  prec_count = N_gas + count_get_gas;


  offset_recv[0] = NumPart - N_gas + prec_count;

  for(i = 1; i < NTask; i++)
    {offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];}


  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_gas[target] > 0 || count_recv_gas[target] > 0)
	    {
            MPI_Sizelimited_Sendrecv(partBuf + offset_gas[target], count_gas[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA_GAS,
			   P + offset_recv_gas[target], count_recv_gas[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA_GAS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Sizelimited_Sendrecv(gasBuf + offset_gas[target], count_gas[target] * sizeof(struct gas_cell_data),
			   MPI_BYTE, target, TAG_GASDATA,
			   CellP + offset_recv_gas[target],
			   count_recv_gas[target] * sizeof(struct gas_cell_data), MPI_BYTE, target,
			   TAG_GASDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef CHIMES 
            MPI_Sizelimited_Sendrecv(gasChimesBuf + offset_gas[target], count_gas[target] * sizeof(struct gasVariables),
			   MPI_BYTE, target, TAG_CHIMESDATA, ChimesGasVars + offset_recv_gas[target],
			   count_recv_gas[target] * sizeof(struct gasVariables), MPI_BYTE, target,
			   TAG_CHIMESDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#ifdef CHIMES_USE_DOUBLE_PRECISION
            MPI_Sizelimited_Sendrecv(gasAbundancesBuf + (offset_gas[target] * ChimesGlobalVars.totalNumberOfSpecies),
			   count_gas[target] * ChimesGlobalVars.totalNumberOfSpecies, MPI_DOUBLE, target, TAG_ABUNDATA, 
			   gasAbundancesRecvBuf + ((offset_recv_gas[target] - offset_recv_gas[0]) * ChimesGlobalVars.totalNumberOfSpecies),
			   count_recv_gas[target] * ChimesGlobalVars.totalNumberOfSpecies, MPI_DOUBLE, target, 
			   TAG_ABUNDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else 
            MPI_Sizelimited_Sendrecv(gasAbundancesBuf + (offset_gas[target] * ChimesGlobalVars.totalNumberOfSpecies),
			   count_gas[target] * ChimesGlobalVars.totalNumberOfSpecies, MPI_FLOAT, target, TAG_ABUNDATA, 
			   gasAbundancesRecvBuf + ((offset_recv_gas[target] - offset_recv_gas[0]) * ChimesGlobalVars.totalNumberOfSpecies),
			   count_recv_gas[target] * ChimesGlobalVars.totalNumberOfSpecies, MPI_FLOAT, target, 
			   TAG_ABUNDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif 
#endif 

            MPI_Sizelimited_Sendrecv(keyBuf + offset_gas[target], count_gas[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY_GAS,
			   Key + offset_recv_gas[target], count_recv_gas[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY_GAS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }


	  if(count[target] > 0 || count_recv[target] > 0)
	    {
	      /* Use the chunking wrapper — fire_m11i with ~6M non-gas parts/rank
	       * × sizeof(struct particle_data)~400 B exceeds int_max for the
	       * MPI_Sendrecv byte count. The gas branch above already uses this
	       * wrapper for the same reason. */
	      MPI_Sizelimited_Sendrecv(partBuf + offset[target], count[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA,
			   P + offset_recv[target], count_recv[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      MPI_Sizelimited_Sendrecv(keyBuf + offset[target], count[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY,
			   Key + offset_recv[target], count_recv[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

#ifdef CHIMES 
  /* Loop through received gas cells and read in abundances from the buffer. */ 
  for (target = 0; target < NTask; target++)
    {
      if(count_recv_gas[target] > 0)
	{
	  for (i = 0; i < count_recv_gas[target]; i++)
	    {
	      allocate_gas_abundances_memory(&(ChimesGasVars[offset_recv_gas[target] + i]), &ChimesGlobalVars); 
	      for (abunIndex = 0; abunIndex < ChimesGlobalVars.totalNumberOfSpecies; abunIndex++)
          {ChimesGasVars[offset_recv_gas[target] + i].abundances[abunIndex] = gasAbundancesRecvBuf[((offset_recv_gas[target] - offset_recv_gas[0] + i) * ChimesGlobalVars.totalNumberOfSpecies) + abunIndex];}
	    }
	}
    }
#endif 

  NumPart += count_get;
  N_gas += count_get_gas;


  /* Defensive post-exchange occupancy checks. The predictive caller guard
   * (domain:exchange_capacity) already prevents reaching domain_exchange() with an
   * overflow, so these should be unreachable; kept as soft belt-and-suspenders, drained
   * at the caller's post-exchange poll. */
  if(NumPart > All.MaxPart)
    {
      printf("Task=%d NumPart=%d All.MaxPart=%d\n", ThisTask, NumPart, All.MaxPart);
      endrun(90000090);
    }

  if(N_gas > All.MaxPartGas) {endrun(90000091);}


  myfree(keyBuf);
#ifdef CHIMES 
  myfree(gasAbundancesRecvBuf);
  myfree(gasAbundancesBuf);
  myfree(gasChimesBuf);
#endif 
  myfree(gasBuf);
  myfree(partBuf);


  myfree(offset_recv_gas);
  myfree(offset_recv);
  myfree(count_recv_gas);
  myfree(count_recv);

  myfree(offset_gas);
  myfree(offset);
  myfree(count_gas);
  myfree(count);
}



void domain_findSplit_work_balanced(int ncpu, int ndomain)
{
  int i, start, end;
  double work, workgas, workavg, work_before, workavg_before;
  double load, fac_work, fac_load, fac_workgas;

  for(i = 0, work = load = workgas = 0; i < ndomain; i++)
    {
      work += domainWork[i];
      load += domainCount[i];
      workgas += domainWorkGas[i];
    }
 
      /* Use adaptive weights that are updated after each decomposition based on measured
         imbalance (GADGET-4 approach). The static variables domain_fac_work/workgas/load
         are adjusted in domain_update_adaptive_weights() after each decomposition. */
      if(workgas <= 0) {
          /* no gas work: redistribute gas weight to gravity and load */
          double total = domain_fac_work + domain_fac_load + domain_fac_workgas;
          fac_work = (domain_fac_work + 0.5 * domain_fac_workgas) / (total * work);
          fac_load = (domain_fac_load + 0.5 * domain_fac_workgas) / (total * load);
          fac_workgas = 0.0;
      } else {
          double total = domain_fac_work + domain_fac_load + domain_fac_workgas;
          fac_work = domain_fac_work / (total * work);
          fac_load = domain_fac_load / (total * load);
          fac_workgas = domain_fac_workgas / (total * workgas);
      }

  workavg = 1.0 / ncpu;

  work_before = workavg_before = 0;

  start = 0;

  for(i = 0; i < ncpu; i++)
    {
      work = 0;
      end = start;

      work += fac_work * domainWork[end] + fac_load * domainCount[end] + fac_workgas * domainWorkGas[end];
      while((work + work_before < workavg + workavg_before) || (i == ncpu - 1 && end < ndomain - 1))
	{
	  if((ndomain - end) > (ncpu - i))
	    {end++;}
	  else
	    {break;}

	  work += fac_work * domainWork[end] + fac_load * domainCount[end] + fac_workgas * domainWorkGas[end];
	}

      DomainStartList[i] = start;
      DomainEndList[i] = end;

      work_before += work;
      workavg_before += workavg;
      start = end + 1;
    }
}

static struct domain_segments_data
{
  int task, start, end;
  double work;
  double load;
  double load_activegas;
  double load_gas;         /* gas particle count for memory constraint */
  double normalized_load;
#if (DOMAIN_TIMEBINS == 1)
  double bin_GravCost[TIMEBINS];   /* per-timebin gravity cost for this segment */
  double bin_HydroCost[TIMEBINS];  /* per-timebin hydro cost for this segment */
#endif
}
 *domainAssign;


struct queue_data
{
  int first, last;
  int *next;
  int *previous;
  double *value;
}
queues[N_DOMAINDECOMP_QUEUES];

struct tasklist_data
{
  double work;
  double load;
  double load_activegas;
  double load_gas;         /* gas particle count for memory constraint */
  int count;
#if (DOMAIN_TIMEBINS == 1)
  double bin_GravCost[TIMEBINS];   /* per-timebin gravity cost accumulated on this task */
  double bin_HydroCost[TIMEBINS];  /* per-timebin hydro cost accumulated on this task */
#endif
}
 *tasklist;

int domain_sort_task(const void *a, const void *b)
{
  if(((struct domain_segments_data *) a)->task < (((struct domain_segments_data *) b)->task)) {return -1;}
  if(((struct domain_segments_data *) a)->task > (((struct domain_segments_data *) b)->task)) {return +1;}
  return 0;
}

int domain_sort_load(const void *a, const void *b)
{
  if(((struct domain_segments_data *) a)->normalized_load >
     (((struct domain_segments_data *) b)->normalized_load)) {return -1;}
  if(((struct domain_segments_data *) a)->normalized_load <
     (((struct domain_segments_data *) b)->normalized_load)) {return +1;}
  return 0;
}

void domain_assign_load_or_work_balanced(int mode, int multipledomains)
{
  double target_work_balance, target_load_balance, target_load_activegas_balance, target_load_gas_balance;
  double value, target_max_balance, best_balance;
  double tot_work, tot_load, tot_loadactivegas, tot_loadgas;


  int best_queue, target, next, prev;
  int i, n, q, ta;

  /* memory limits: reject assignments that would push a rank above this fraction of max capacity */
  double memory_safety_frac = 0.9;
  double max_load_per_task = memory_safety_frac * maxLoad;
  double max_gasload_per_task = memory_safety_frac * maxLoadgas;

  domainAssign = (struct domain_segments_data *) mymalloc("domainAssign",
							  multipledomains * NTask *
							  sizeof(struct domain_segments_data));

  tasklist = (struct tasklist_data *) mymalloc("tasklist", NTask * sizeof(struct tasklist_data));

  for(ta = 0; ta < NTask; ta++)
    {
      tasklist[ta].work = 0;
      tasklist[ta].load = 0;
      tasklist[ta].load_activegas = 0;
      tasklist[ta].load_gas = 0;
      tasklist[ta].count = 0;
#if (DOMAIN_TIMEBINS == 1)
      for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {tasklist[ta].bin_GravCost[k] = 0; tasklist[ta].bin_HydroCost[k] = 0;}
#endif
    }

  tot_work = 0;
  tot_load = 0;
  tot_loadactivegas = 0;
  tot_loadgas = 0;

#if (DOMAIN_TIMEBINS == 1)
  /* Per-timebin cost totals for normalization during imbalance evaluation */
  double tot_binGravCost[TIMEBINS], tot_binHydroCost[TIMEBINS];
  for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {tot_binGravCost[k] = 0; tot_binHydroCost[k] = 0;}
#endif

  for(n = 0; n < multipledomains * NTask; n++)
    {
      domainAssign[n].start = DomainStartList[n];
      domainAssign[n].end = DomainEndList[n];
      domainAssign[n].work = 0;
      domainAssign[n].load = 0;
      domainAssign[n].load_activegas = 0;
      domainAssign[n].load_gas = 0;
#if (DOMAIN_TIMEBINS == 1)
      for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {domainAssign[n].bin_GravCost[k] = 0; domainAssign[n].bin_HydroCost[k] = 0;}
#endif

      for(i = DomainStartList[n]; i <= DomainEndList[n]; i++)
	{
	  domainAssign[n].work += domainWork[i];
	  domainAssign[n].load += domainCount[i];
	  domainAssign[n].load_activegas += domainWorkGas[i];
	  domainAssign[n].load_gas += domainCountGas[i];
#if (DOMAIN_TIMEBINS == 1)
	  for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {
	      domainAssign[n].bin_GravCost[k] += domainBinGravCost[k * NTopleaves + i];
	      domainAssign[n].bin_HydroCost[k] += domainBinHydroCost[k * NTopleaves + i];
	  }
#endif
	}

      tot_work += domainAssign[n].work;
      tot_load += domainAssign[n].load;
      tot_loadactivegas += domainAssign[n].load_activegas;
      tot_loadgas += domainAssign[n].load_gas;
#if (DOMAIN_TIMEBINS == 1)
      for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {
          tot_binGravCost[k] += domainAssign[n].bin_GravCost[k];
          tot_binHydroCost[k] += domainAssign[n].bin_HydroCost[k];
      }
#endif
    }

  for(n = 0; n < multipledomains * NTask; n++)
    {
        if(mode==1) {
            domainAssign[n].normalized_load = domainAssign[n].work / (tot_work + 1.0e-30) + domainAssign[n].load_activegas / (tot_loadactivegas + 1.0e-30);
        } else {
            domainAssign[n].normalized_load = domainAssign[n].load / (tot_load + 1.0e-30);
        }
    }

  qsort(domainAssign, multipledomains * NTask, sizeof(struct domain_segments_data), domain_sort_load);

  /* initialize queues (4 queues: work, load, active gas, gas memory) */
  for(q = 0; q < N_DOMAINDECOMP_QUEUES; q++)
    {
      queues[q].next = (int *) mymalloc("queues[q].next", NTask * sizeof(int));
      queues[q].previous = (int *) mymalloc("queues[q].previous", NTask * sizeof(int));
      queues[q].value = (double *) mymalloc("queues[q].value", NTask * sizeof(double));

      for(ta = 0; ta < NTask; ta++)
	{
	  queues[q].next[ta] = ta + 1;
	  queues[q].previous[ta] = ta - 1;
	  queues[q].value[ta] = 0;
	}
      queues[q].previous[0] = -1;
      queues[q].next[NTask - 1] = -1;
      queues[q].first = 0;
      queues[q].last = NTask - 1;
    }

  for(n = 0; n < multipledomains * NTask; n++)
    {
      /* need to decide, which of the tasks that has the lowest load in one of the queues is best */
      for(q = 0, best_balance = 1.0e30, best_queue = 0; q < N_DOMAINDECOMP_QUEUES; q++)
	{
	  target = queues[q].first;

	  while(tasklist[target].count == multipledomains) {target = queues[q].next[target];}

	  /* check memory hard constraint: skip targets that would exceed memory limits */
	  if(mode == 1) {
	      double candidate_load = tasklist[target].load + domainAssign[n].load;
	      double candidate_gasload = tasklist[target].load_gas + domainAssign[n].load_gas;
	      if(candidate_load > max_load_per_task || candidate_gasload > max_gasload_per_task) {
	          /* try next targets in this queue until we find one that fits */
	          int orig_target = target;
	          target = queues[q].next[target];
	          while(target >= 0) {
	              if(tasklist[target].count < multipledomains) {
	                  candidate_load = tasklist[target].load + domainAssign[n].load;
	                  candidate_gasload = tasklist[target].load_gas + domainAssign[n].load_gas;
	                  if(candidate_load <= max_load_per_task && candidate_gasload <= max_gasload_per_task) break;
	              }
	              target = queues[q].next[target];
	          }
	          if(target < 0) target = orig_target; /* fallback: accept overload rather than crash */
	      }
	  }

	  target_work_balance = (domainAssign[n].work + tasklist[target].work) / (tot_work + 1.0e-30);
	  target_load_balance = (domainAssign[n].load + tasklist[target].load) / (tot_load + 1.0e-30);
	  target_load_activegas_balance = (domainAssign[n].load_activegas + tasklist[target].load_activegas) / (tot_loadactivegas + 1.0e-30);
	  target_load_gas_balance = (domainAssign[n].load_gas + tasklist[target].load_gas) / (tot_loadgas + 1.0e-30);

        if(mode==1) {
            target_max_balance = target_work_balance;
            if(target_max_balance < target_load_balance) {target_max_balance = target_load_balance;}
            if(target_max_balance < target_load_activegas_balance) {target_max_balance = target_load_activegas_balance;}
            if(target_max_balance < target_load_gas_balance) {target_max_balance = target_load_gas_balance;}
#if (DOMAIN_TIMEBINS == 1)
            /* Also check per-timebin cost imbalance: the assignment must balance
               each timebin individually, not just the composite cost */
            for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {
                double bg = (domainAssign[n].bin_GravCost[k] + tasklist[target].bin_GravCost[k]) / (tot_binGravCost[k] + 1.0e-30);
                double bh = (domainAssign[n].bin_HydroCost[k] + tasklist[target].bin_HydroCost[k]) / (tot_binHydroCost[k] + 1.0e-30);
                if(bg > target_max_balance) {target_max_balance = bg;}
                if(bh > target_max_balance) {target_max_balance = bh;}
            }
#endif
        } else {
            target_max_balance = target_load_balance;
        }

	  if(target_max_balance < best_balance)
	    {
	      best_balance = target_max_balance;
	      best_queue = q;
	    }
	}

      /* Now we know the best queue, and hence the best target task. Assign this piece to this task */
      target = queues[best_queue].first;

      while(tasklist[target].count == multipledomains) {target = queues[best_queue].next[target];}

      /* re-check memory constraint for the selected target */
      if(mode == 1) {
          double candidate_load = tasklist[target].load + domainAssign[n].load;
          double candidate_gasload = tasklist[target].load_gas + domainAssign[n].load_gas;
          if(candidate_load > max_load_per_task || candidate_gasload > max_gasload_per_task) {
              int next_t = queues[best_queue].next[target];
              while(next_t >= 0) {
                  if(tasklist[next_t].count < multipledomains &&
                     tasklist[next_t].load + domainAssign[n].load <= max_load_per_task &&
                     tasklist[next_t].load_gas + domainAssign[n].load_gas <= max_gasload_per_task) {
                      target = next_t; break;
                  }
                  next_t = queues[best_queue].next[next_t];
              }
          }
      }

      domainAssign[n].task = target;
      tasklist[target].work += domainAssign[n].work;
      tasklist[target].load += domainAssign[n].load;
      tasklist[target].load_activegas += domainAssign[n].load_activegas;
      tasklist[target].load_gas += domainAssign[n].load_gas;
#if (DOMAIN_TIMEBINS == 1)
      for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {
          tasklist[target].bin_GravCost[k] += domainAssign[n].bin_GravCost[k];
          tasklist[target].bin_HydroCost[k] += domainAssign[n].bin_HydroCost[k];
      }
#endif
      tasklist[target].count++;

      /* now we need to remove the element 'target' from the queues and reinsert it */
    for(q = 0; q < N_DOMAINDECOMP_QUEUES; q++)
    {
	  switch (q)
	    {
	    case 0:
	      value = tasklist[target].work;
	      break;
	    case 1:
	      value = tasklist[target].load;
	      break;
	    case 2:
	      value = tasklist[target].load_activegas;
	      break;
	    case 3:
	      value = tasklist[target].load_gas;
	      break;
	    default:
	      value = 0;
	      break;
	    }

	  /* now remove the element target */
	  prev = queues[q].previous[target];
	  next = queues[q].next[target];

	  if(prev >= 0)		/* previous exists */
	    {queues[q].next[prev] = next;}
	  else
	    {queues[q].first = next;}	/* we remove the head of the queue */


	  if(next >= 0)		/* next exists */
	    {queues[q].previous[next] = prev;}
	  else
	    {queues[q].last = prev;}	/* we remove the end of the queue */

	  /* now we insert the element again, in an ordered fashion, starting from the end of the queue */
	  if(queues[q].last >= 0)
	    {
	      ta = queues[q].last;

	      while(value < queues[q].value[ta])
		{
		  ta = queues[q].previous[ta];
		  if(ta < 0) {break;}
		}

	      if(ta < 0)	/* we insert the element as the first element */
		{
		  queues[q].next[target] = queues[q].first;
		  queues[q].previous[queues[q].first] = target;
		  queues[q].first = target;
		}
	      else
		{
		  /* insert behind ta */
		  queues[q].next[target] = queues[q].next[ta];
		  if(queues[q].next[ta] >= 0)
		    {queues[q].previous[queues[q].next[ta]] = target;}
		  else
		    {queues[q].last = target;}	/* we insert a new last element */
		  queues[q].previous[target] = ta;
		  queues[q].next[ta] = target;
		}
	    }
	  else
	    {
	      /* queue was empty */
	      queues[q].previous[target] = queues[q].next[target] = -1;
	      queues[q].first = queues[q].last = target;
	    }

	  queues[q].value[target] = value;
	}
    }

  /* Iterative refinement (GADGET-4 approach): try random swaps of domain segments
     between tasks and keep swaps that reduce the maximum imbalance across all metrics.
     This improves the greedy solution, especially for inhomogeneous particle distributions. */
  if(mode == 1 && multipledomains * NTask > 1)
  {
      int nswaps_accepted = 0;
      int n_segments = multipledomains * NTask;
      int max_iterations = 200 * n_segments; /* scale attempts with problem size */
      unsigned int seed = 42 + NTask; /* deterministic seed for reproducibility */

      /* compute current max imbalance */
      auto compute_max_imbalance = [&]() -> double {
          double max_frac_work = 0, max_frac_load = 0, max_frac_gas = 0, max_frac_gasload = 0;
          for(int t = 0; t < NTask; t++) {
              double fw = tasklist[t].work / (tot_work + 1.0e-30);
              double fl = tasklist[t].load / (tot_load + 1.0e-30);
              double fg = tasklist[t].load_activegas / (tot_loadactivegas + 1.0e-30);
              double fgl = tasklist[t].load_gas / (tot_loadgas + 1.0e-30);
              if(fw > max_frac_work) max_frac_work = fw;
              if(fl > max_frac_load) max_frac_load = fl;
              if(fg > max_frac_gas) max_frac_gas = fg;
              if(fgl > max_frac_gasload) max_frac_gasload = fgl;
          }
          double result = max_frac_work;
          if(max_frac_load > result) result = max_frac_load;
          if(max_frac_gas > result) result = max_frac_gas;
          if(max_frac_gasload > result) result = max_frac_gasload;
#if (DOMAIN_TIMEBINS == 1)
          for(int t = 0; t < NTask; t++) {
              for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {
                  double bg = tasklist[t].bin_GravCost[k] / (tot_binGravCost[k] + 1.0e-30);
                  double bh = tasklist[t].bin_HydroCost[k] / (tot_binHydroCost[k] + 1.0e-30);
                  if(bg > result) result = bg;
                  if(bh > result) result = bh;
              }
          }
#endif
          return result;
      };

      double current_imbalance = compute_max_imbalance();

      for(int iter = 0; iter < max_iterations; iter++)
      {
          /* pick two random segments assigned to different tasks */
          seed = seed * 1103515245 + 12345; int s1 = (seed >> 16) % n_segments;
          seed = seed * 1103515245 + 12345; int s2 = (seed >> 16) % n_segments;
          if(s1 == s2) continue;
          int t1 = domainAssign[s1].task, t2 = domainAssign[s2].task;
          if(t1 == t2) continue;

          /* trial swap: move s1 to t2 and s2 to t1 */
          tasklist[t1].work += domainAssign[s2].work - domainAssign[s1].work;
          tasklist[t1].load += domainAssign[s2].load - domainAssign[s1].load;
          tasklist[t1].load_activegas += domainAssign[s2].load_activegas - domainAssign[s1].load_activegas;
          tasklist[t1].load_gas += domainAssign[s2].load_gas - domainAssign[s1].load_gas;
          tasklist[t2].work += domainAssign[s1].work - domainAssign[s2].work;
          tasklist[t2].load += domainAssign[s1].load - domainAssign[s2].load;
          tasklist[t2].load_activegas += domainAssign[s1].load_activegas - domainAssign[s2].load_activegas;
          tasklist[t2].load_gas += domainAssign[s1].load_gas - domainAssign[s2].load_gas;
#if (DOMAIN_TIMEBINS == 1)
          for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {
              tasklist[t1].bin_GravCost[k] += domainAssign[s2].bin_GravCost[k] - domainAssign[s1].bin_GravCost[k];
              tasklist[t1].bin_HydroCost[k] += domainAssign[s2].bin_HydroCost[k] - domainAssign[s1].bin_HydroCost[k];
              tasklist[t2].bin_GravCost[k] += domainAssign[s1].bin_GravCost[k] - domainAssign[s2].bin_GravCost[k];
              tasklist[t2].bin_HydroCost[k] += domainAssign[s1].bin_HydroCost[k] - domainAssign[s2].bin_HydroCost[k];
          }
#endif

          double new_imbalance = compute_max_imbalance();

          if(new_imbalance < current_imbalance) {
              /* accept swap */
              domainAssign[s1].task = t2;
              domainAssign[s2].task = t1;
              current_imbalance = new_imbalance;
              nswaps_accepted++;
          } else {
              /* revert */
              tasklist[t1].work -= domainAssign[s2].work - domainAssign[s1].work;
              tasklist[t1].load -= domainAssign[s2].load - domainAssign[s1].load;
              tasklist[t1].load_activegas -= domainAssign[s2].load_activegas - domainAssign[s1].load_activegas;
              tasklist[t1].load_gas -= domainAssign[s2].load_gas - domainAssign[s1].load_gas;
              tasklist[t2].work -= domainAssign[s1].work - domainAssign[s2].work;
              tasklist[t2].load -= domainAssign[s1].load - domainAssign[s2].load;
              tasklist[t2].load_activegas -= domainAssign[s1].load_activegas - domainAssign[s2].load_activegas;
              tasklist[t2].load_gas -= domainAssign[s1].load_gas - domainAssign[s2].load_gas;
#if (DOMAIN_TIMEBINS == 1)
              for(int k = 0; k < NumTimeBinsToBeBalanced; k++) {
                  tasklist[t1].bin_GravCost[k] -= domainAssign[s2].bin_GravCost[k] - domainAssign[s1].bin_GravCost[k];
                  tasklist[t1].bin_HydroCost[k] -= domainAssign[s2].bin_HydroCost[k] - domainAssign[s1].bin_HydroCost[k];
                  tasklist[t2].bin_GravCost[k] -= domainAssign[s1].bin_GravCost[k] - domainAssign[s2].bin_GravCost[k];
                  tasklist[t2].bin_HydroCost[k] -= domainAssign[s1].bin_HydroCost[k] - domainAssign[s2].bin_HydroCost[k];
              }
#endif
          }
      }
      if(ThisTask == 0 && nswaps_accepted > 0) {
          PRINT_STATUS(" ..domain iterative refinement: %d swaps accepted out of %d attempts", nswaps_accepted, max_iterations);
      }
  }

  qsort(domainAssign, multipledomains * NTask, sizeof(struct domain_segments_data), domain_sort_task);

  for(n = 0; n < multipledomains * NTask; n++)
    {
      DomainStartList[n] = domainAssign[n].start;
      DomainEndList[n] = domainAssign[n].end;

      for(i = DomainStartList[n]; i <= DomainEndList[n]; i++) {DomainTask[i] = domainAssign[n].task;}
    }

  /* free the queues */
  for(q = N_DOMAINDECOMP_QUEUES-1; q >= 0; q--)
    {
      myfree(queues[q].value);
      myfree(queues[q].previous);
      myfree(queues[q].next);
    }

  myfree(tasklist);

  myfree(domainAssign);
}


void domain_findSplit_load_balanced(int ncpu, int ndomain)
{
  int i, start, end;
  double load, loadavg, load_before, loadavg_before, fac_load, fac;
  for(i = 0, load = 0; i < ndomain; i++)
    {
      load += domainCount[i];
    }
  fac = 1;
  fac_load = fac / load;
  loadavg = 1.0 / ncpu;
  load_before = loadavg_before = 0;
  start = 0;
  for(i = 0; i < ncpu; i++)
    {
      load = 0;
      end = start;
      load += fac_load * domainCount[end];
      while((load + load_before < loadavg + loadavg_before) || (i == ncpu - 1 && end < ndomain - 1))
	{
	  if((ndomain - end) > (ncpu - i)) {end++;} else {break;}
	  load += fac_load * domainCount[end];
	}
      DomainStartList[i] = start;
      DomainEndList[i] = end;
      load_before += load;
      loadavg_before += loadavg;
      start = end + 1;
    }
}







/*! This function determines how many particles that are currently stored
 *  on the local CPU have to be moved off according to the domain
 *  decomposition.
 */
int domain_countToGo(size_t nlimit)
{
    int n, no, ret, retsum; size_t package;

    for(n = 0; n < NTask; n++)
    {
        toGo[n] = 0;
        toGoGas[n] = 0;
    }

    package = (sizeof(struct particle_data) + sizeof(struct gas_cell_data) + sizeof(peanokey));
    /* cannot pack even one particle into the exchange budget. Soft bad-stop: the
     * packing loop below no-ops (package>=nlimit => 0 iterations, toGo stays zero so the
     * Alltoalls move nothing and stay matched) and the caller's exchange-capacity poll
     * drains it. The caller's arena guard normally catches this first. */
    if(package >= nlimit) {endrun(90000017);}

    for(n = 0; n < NumPart && package < nlimit; n++)
    {
#ifdef SUBFIND
        if(GrNr >= 0 && P[n].GrNr != GrNr) {continue;}
#endif
        if(P[n].Type & 32)
        {
            no = domain_toptree_leaf(Key[n], topNodes);

            if(DomainTask[no] != ThisTask)
            {
                toGo[DomainTask[no]] += 1;
                nlimit -= sizeof(struct particle_data) + sizeof(peanokey);
                
                if((P[n].Type & 15) == 0)
                {
                    toGoGas[DomainTask[no]] += 1;
                    nlimit -= sizeof(struct gas_cell_data);
                }
                P[n].Type |= 16;    /* flag this particle for export */
            }
        }
    }

    MPI_Alltoall(toGo, 1, MPI_INT, toGet, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Alltoall(toGoGas, 1, MPI_INT, toGetGas, 1, MPI_INT, MPI_COMM_WORLD);

    if(package >= nlimit) {ret = 1;} else {ret = 0;}
    MPI_Allreduce(&ret, &retsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(retsum)
    {
        /* in this case, we are not guaranteed that the temporary state after the partial exchange will actually observe the particle limits on all
         processors... we need to test this explicitly and rework the exchange such that this is guaranteed. This is actually a rather non-trivial constraint. */
        
        MPI_Allgather(&NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&N_gas, 1, MPI_INT, list_N_gas, 1, MPI_INT, MPI_COMM_WORLD);
        int flag, flagsum, ntoomany, ta, i, target;
        int count_togo, count_toget, count_togo_gas, count_toget_gas;
        
        do
        {
            flagsum = 0;
            do
            {
                flag = 0;
                for(ta = 0; ta < NTask; ta++)
                {
                    if(ta == ThisTask)
                    {
                        count_togo = count_toget = 0;
                        count_togo_gas = count_toget_gas = 0;
                        for(i = 0; i < NTask; i++)
                        {
                            count_togo += toGo[i];
                            count_toget += toGet[i];
                            count_togo_gas += toGoGas[i];
                            count_toget_gas += toGetGas[i];
                        }
                    }
                    MPI_Bcast(&count_togo, 1, MPI_INT, ta, MPI_COMM_WORLD);
                    MPI_Bcast(&count_toget, 1, MPI_INT, ta, MPI_COMM_WORLD);
                    MPI_Bcast(&count_togo_gas, 1, MPI_INT, ta, MPI_COMM_WORLD);
                    MPI_Bcast(&count_toget_gas, 1, MPI_INT, ta, MPI_COMM_WORLD);
                    
                    int ifntoomany;
                    ntoomany = list_N_gas[ta] + count_toget_gas - count_togo_gas - All.MaxPartGas;
                    ifntoomany = (ntoomany > 0);
                    if(ifntoomany)
                    {
                        if(ntoomany > 0 && ThisTask==0)
                        {
                            if(flagsum < 25) {PRINT_STATUS(" ..domain exchange must be modified - cannot receive %d gas elements on task=%d (iter=%d)", ntoomany, ta, flagsum+1);}
                            else {
                                printf(" ..domain exchange must be modified - cannot receive %d gas elements on task=%d (iter=%d)\n", ntoomany, ta, flagsum+1);
                                printf(" ..list_N_gas[ta=%d]=%d  count_toget_gas=%d count_togo_gas=%d MaxPartGas=%d NTask=%d flagsum=%d\n", ta, list_N_gas[ta], count_toget_gas, count_togo_gas,All.MaxPartGas,NTask,flagsum); fflush(stdout);
                            }
                        }
                        flag = 1;
                        i = flagsum % NTask;
                        while(ifntoomany)
                        {
                            if(i == ThisTask)
                            {
                                if(toGoGas[ta] > 0)
                                    if(ntoomany > 0)
                                    {
                                        toGoGas[ta]--;
                                        count_toget_gas--;
                                        count_toget--;
                                        ntoomany--;
                                    }
                            }
                            
                            MPI_Bcast(&ntoomany, 1, MPI_INT, i, MPI_COMM_WORLD);
                            MPI_Bcast(&count_toget, 1, MPI_INT, i, MPI_COMM_WORLD);
                            MPI_Bcast(&count_toget_gas, 1, MPI_INT, i, MPI_COMM_WORLD);
                            i++;
                            if(i >= NTask) {i = 0;}
                            
                            ifntoomany = (ntoomany > 0);
                        }
                    }
                    
                    ntoomany = list_NumPart[ta] + count_toget - count_togo - All.MaxPart;
                    ifntoomany = (ntoomany > 0);
                    if(ifntoomany) {
                        if(ntoomany > 0 && ThisTask==0) {
                            if(flagsum < 25) {PRINT_STATUS(" ..domain exchange must be modified - cannot receive %d elements on task=%d (iter=%d)", ntoomany, ta, flagsum+1);}
                            else {
                                printf(" ..domain exchange must be modified - cannot receive %d elements on task=%d (iter=%d)\n", ntoomany, ta, flagsum+1);
                                printf(" ..list_NumPart[ta=%d]=%d count_toget=%d count_togo=%d MaxPart=%d NTask=%d flagsum=%d \n", ta, list_NumPart[ta], count_toget, count_togo, All.MaxPart, NTask, flagsum); fflush(stdout);
                            }
                        }
                    
                        flag = 1;
                        i = flagsum % NTask;
                        while(ntoomany)
                        {
                            if(i == ThisTask)
                            {
                                if(toGo[ta] > 0)
                                {
                                    toGo[ta]--;
                                    count_toget--;
                                    ntoomany--;
                                }
                            }
                            
                            MPI_Bcast(&ntoomany, 1, MPI_INT, i, MPI_COMM_WORLD);
                            MPI_Bcast(&count_toget, 1, MPI_INT, i, MPI_COMM_WORLD);
                            
                            i++;
                            if(i >= NTask) {i = 0;}
                        }
                    }
                }
                flagsum += flag;
                
                /* Symmetric: flagsum is driven only by MPI_Bcast/MPI_Allgather-
                 * synchronized counters, so it is identical on every rank and this
                 * branch fires on all ranks together. Soft bad-stop + an immediate
                 * all-rank poll (error-branch only, no normal-path cost) drains here
                 * instead of re-entering the inner do/while forever. */
                if(flagsum > 100) {if(ThisTask==0) {printf("Failed to converge in domain.c, flagsum=%d",flagsum); fflush(stdout);} endrun(90000018); gizmo_exit_bad_stop_if_requested("domain:countToGo_converge");}
                MPI_Alltoall(toGo, 1, MPI_INT, toGet, 1, MPI_INT, MPI_COMM_WORLD);
                MPI_Alltoall(toGoGas, 1, MPI_INT, toGetGas, 1, MPI_INT, MPI_COMM_WORLD);
            }
            while(flag);
            
            if(flagsum)
            {
                int *local_toGo, *local_toGoGas;
                local_toGo = (int *) mymalloc("          local_toGo", NTask * sizeof(int));
                local_toGoGas = (int *) mymalloc("          local_toGoGas", NTask * sizeof(int));
                
                for(n = 0; n < NTask; n++)
                {
                    local_toGo[n] = 0;
                    local_toGoGas[n] = 0;
                }
                
                for(n = 0; n < NumPart; n++)
                {
                    if(P[n].Type & 32)
                    {
                        P[n].Type &= (15 + 32);    /* clear 16 */
                        
                        no = domain_toptree_leaf(Key[n], topNodes);
                        target = DomainTask[no];
                        
                        if((P[n].Type & 15) == 0)
                        {
                            if(local_toGoGas[target] < toGoGas[target] && local_toGo[target] < toGo[target])
                            {
                                local_toGo[target] += 1;
                                local_toGoGas[target] += 1;
                                P[n].Type |= 16;
                            }
                        }
                        else
                        {
                            if(local_toGo[target] < toGo[target])
                            {
                                local_toGo[target] += 1;
                                P[n].Type |= 16;
                            }
                        }
                    }
                }
                
                for(n = 0; n < NTask; n++)
                {
                    toGo[n] = local_toGo[n];
                    toGoGas[n] = local_toGoGas[n];
                }
                
                MPI_Alltoall(toGo, 1, MPI_INT, toGet, 1, MPI_INT, MPI_COMM_WORLD);
                MPI_Alltoall(toGoGas, 1, MPI_INT, toGetGas, 1, MPI_INT, MPI_COMM_WORLD);
                myfree(local_toGoGas);
                myfree(local_toGo);
            }
        }
        while(flagsum);
        
        return 1;
    }
    else
        {return 0;}
}






/*! This function walks the global top tree in order to establish the
 *  number of leaves it has. These leaves are distributed to different
 *  processors.
 */
void domain_walktoptree(int no)
{
  int i;
  if(topNodes[no].Daughter == -1)
    {
      topNodes[no].Leaf = NTopleaves;
      NTopleaves++;
    }
  else
    {
      for(i = 0; i < 8; i++) {domain_walktoptree(topNodes[no].Daughter + i);}
    }
}


int domain_compare_key(const void *a, const void *b)
{
  if(((struct peano_hilbert_data *) a)->key < (((struct peano_hilbert_data *) b)->key)) {return -1;}
  if(((struct peano_hilbert_data *) a)->key > (((struct peano_hilbert_data *) b)->key)) {return +1;}
  return 0;
}


int domain_check_for_local_refine(int i, double countlimit, double costlimit)
{
  int j, p, sub, flag = 0;
  
  if(topNodes[i].Parent >= 0)
    {
      if(topNodes[i].Count > 0.8 * topNodes[topNodes[i].Parent].Count || topNodes[i].Cost > 0.8 * topNodes[topNodes[i].Parent].Cost) {flag = 1;}
    }

  if((topNodes[i].Count > countlimit || topNodes[i].Cost > costlimit || flag == 1) && topNodes[i].Size >= 8)
    {
      if(topNodes[i].Size >= 8)
	{
	  if((NTopnodes + 8) <= MaxTopNodes)
	    {
	      topNodes[i].Daughter = NTopnodes;

	      for(j = 0; j < 8; j++)
		{
		  sub = topNodes[i].Daughter + j;
		  topNodes[sub].Daughter = -1;
		  topNodes[sub].Parent = i;
		  topNodes[sub].Size = (topNodes[i].Size >> 3);
		  topNodes[sub].StartKey = topNodes[i].StartKey + j * topNodes[sub].Size;
		  topNodes[sub].PIndex = topNodes[i].PIndex;
		  topNodes[sub].Count = 0;
		  topNodes[sub].Cost = 0;

		}

	      NTopnodes += 8;

	      sub = topNodes[i].Daughter;

	      for(p = topNodes[i].PIndex, j = 0; p < topNodes[i].PIndex + topNodes[i].Count; p++)
		{
		  if(j < 7)
		    while(mp[p].key >= topNodes[sub + 1].StartKey)
		      {
			j++;
			sub++;
			topNodes[sub].PIndex = p;
			if(j >= 7) {break;}
		      }

		  topNodes[sub].Cost += particle_total_cost[mp[p].index];
		  topNodes[sub].Count++;
		}

	      for(j = 0; j < 8; j++)
		{
		  sub = topNodes[i].Daughter + j;

		  if(domain_check_for_local_refine(sub, countlimit, costlimit)) {return 1;}
		}
	    }
	  else
	    {return 1;}
	}
    }

  return 0;
}


int domain_recursively_combine_topTree(int start, int ncpu)
{
  int i, nleft, nright, errflag = 0;
  int recvTask, ntopnodes_import;
  int domainkey_top_left, domainkey_top_right;
  struct local_topnode_data *topNodes_import = 0, *topNodes_temp;

  nleft = ncpu / 2;
  nright = ncpu - nleft;

  if(ncpu > 2)
    {
      errflag += domain_recursively_combine_topTree(start, nleft);
      errflag += domain_recursively_combine_topTree(start + nleft, nright);
    }

  if(ncpu >= 2)
    {
      domainkey_top_left = start;
      domainkey_top_right = start + nleft;
      /* Symmetric (start/ncpu identical on every rank); impossible under correct
       * logic (nleft=ncpu/2>=1 when ncpu>=2). Soft bad-stop + status-return BEFORE
       * the guarded Sendrecv: all ranks return together (no deadlock) and errflag
       * propagates to the MPI_Allreduce(&errflag,&errsum) -> TopNodeAllocFactor retry. */
      if(domainkey_top_left == domainkey_top_right) {endrun(90000019); return errflag + 1;}

      if(ThisTask == domainkey_top_left || ThisTask == domainkey_top_right)
	{
	  if(ThisTask == domainkey_top_left)
	    {recvTask = domainkey_top_right;}
	  else
	    {recvTask = domainkey_top_left;}

	  /* inform each other about the length of the trees */
	  MPI_Sendrecv(&NTopnodes, 1, MPI_INT, recvTask, TAG_GRAV_A,
		       &ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MPI_COMM_WORLD,
		       MPI_STATUS_IGNORE);


	  topNodes_import = (struct local_topnode_data *) mymalloc("topNodes_import",
						   IMAX(ntopnodes_import, NTopnodes) * sizeof(struct local_topnode_data));

	  /* exchange the trees */
	  MPI_Sendrecv(topNodes,
		       NTopnodes * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B,
		       topNodes_import,
		       ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      if(ThisTask == domainkey_top_left)
	{
	  for(recvTask = domainkey_top_left + 1; recvTask < domainkey_top_left + nleft; recvTask++)
	    {
	      MPI_Send(&ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MPI_COMM_WORLD);
	      MPI_Send(topNodes_import,
		       ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B, MPI_COMM_WORLD);
	    }
	}

      if(ThisTask == domainkey_top_right)
	{
	  for(recvTask = domainkey_top_right + 1; recvTask < domainkey_top_right + nright; recvTask++)
	    {
	      MPI_Send(&ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MPI_COMM_WORLD);
	      MPI_Send(topNodes_import,
		       ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B, MPI_COMM_WORLD);
	    }
	}

      if(ThisTask > domainkey_top_left && ThisTask < domainkey_top_left + nleft)
	{
	  MPI_Recv(&ntopnodes_import, 1, MPI_INT, domainkey_top_left, TAG_GRAV_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  topNodes_import =
	    (struct local_topnode_data *) mymalloc("topNodes_import",
						   IMAX(ntopnodes_import,
							NTopnodes) * sizeof(struct local_topnode_data));

	  MPI_Recv(topNodes_import,
		   ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		   domainkey_top_left, TAG_GRAV_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}


      if(ThisTask > domainkey_top_right && ThisTask < domainkey_top_right + nright)
	{
	  MPI_Recv(&ntopnodes_import, 1, MPI_INT, domainkey_top_right, TAG_GRAV_A, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

	  topNodes_import =
	    (struct local_topnode_data *) mymalloc("topNodes_import",
						   IMAX(ntopnodes_import,
							NTopnodes) * sizeof(struct local_topnode_data));

	  MPI_Recv(topNodes_import,
		   ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		   domainkey_top_right, TAG_GRAV_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      if(ThisTask >= domainkey_top_left && ThisTask < domainkey_top_left + nleft)
	{
	  /* swap the two trees so that result will be equal on all cpus */

	  topNodes_temp =
	    (struct local_topnode_data *) mymalloc("topNodes_temp",
						   NTopnodes * sizeof(struct local_topnode_data));
	  memcpy(topNodes_temp, topNodes, NTopnodes * sizeof(struct local_topnode_data));
	  memcpy(topNodes, topNodes_import, ntopnodes_import * sizeof(struct local_topnode_data));
	  memcpy(topNodes_import, topNodes_temp, NTopnodes * sizeof(struct local_topnode_data));
	  myfree(topNodes_temp);
	  i = NTopnodes;
	  NTopnodes = ntopnodes_import;
	  ntopnodes_import = i;
	}

      if(ThisTask >= start && ThisTask < start + ncpu)
	{
	  if(errflag == 0)
	    {
	      if((NTopnodes + ntopnodes_import) <= MaxTopNodes)
		{
		  errflag += domain_insertnode(topNodes, topNodes_import, 0, 0);
		}
	      else
		{
		  errflag += 1;
		}
	    }

	  myfree(topNodes_import);
	}
    }

  return errflag;
}


#ifdef ALT_QSORT
#define KEY_TYPE struct peano_hilbert_data
#define KEY_BASE_TYPE peanokey
#define KEY_GETVAL(pk) ((pk)->key)
#define KEY_COPY(pk1,pk2)       \
  {                               \
    (pk2)->key = (pk1)->key;      \
    (pk2)->index = (pk1)->index;  \
  }
#define QSORT qsort_domain
#include "system/myqsort.h"
#endif

/*! This function constructs the global top-level tree node that is used
 *  for the domain decomposition. This is done by considering the string of
 *  Peano-Hilbert keys for all particles, which is recursively chopped off
 *  in pieces of eight segments until each segment holds at most a certain
 *  number of particles.
 */
int domain_determineTopTree(void)
{
  int i, count, j, sub, ngrp;
  int recvTask, sendTask, ntopnodes_import, errflag, errsum;
  struct local_topnode_data *topNodes_import, *topNodes_temp;
  double costlimit, countlimit;
  MPI_Status status;
  int multipledomains = MULTIPLEDOMAINS;

  mp = (struct peano_hilbert_data *) mymalloc("mp", sizeof(struct peano_hilbert_data) * NumPart);

  /* Compute Peano-Hilbert keys for all particles */
#ifdef SUBFIND
  /* With SUBFIND, some particles may be skipped, so compute keys in parallel then compact */
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(i = 0; i < NumPart; i++)
    {
      peano1D xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) / DomainLen) + 1.0);
      peano1D yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) / DomainLen) + 1.0);
      peano1D zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) / DomainLen) + 1.0);
      Key[i] = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
    }
  for(i = 0, count = 0; i < NumPart; i++)
    {
      if(GrNr >= 0 && P[i].GrNr != GrNr) {continue;}
      mp[count].key = Key[i];
      mp[count].index = i;
      count++;
    }
#else
  /* Without SUBFIND, count == i always, so the loop is embarrassingly parallel */
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(i = 0; i < NumPart; i++)
    {
      peano1D xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) / DomainLen) + 1.0);
      peano1D yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) / DomainLen) + 1.0);
      peano1D zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) / DomainLen) + 1.0);
      mp[i].key = Key[i] = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
      mp[i].index = i;
    }
  count = NumPart;
#endif

#ifdef SUBFIND
  if(GrNr >= 0 && count != NumPartGroup)
    endrun(1222);   /* SUBFIND group-count invariant; soft stop drains at the caller's poll after domain_Decomposition (continuation here is bounded: mysort_domain uses count<=NumPart) */
#endif

  mysort_domain(mp, count, sizeof(struct peano_hilbert_data));

  NTopnodes = 1;
  topNodes[0].Daughter = -1;
  topNodes[0].Parent = -1;
  topNodes[0].Size = PEANOCELLS;
  topNodes[0].StartKey = 0;
  topNodes[0].PIndex = 0;
  topNodes[0].Count = count;
  topNodes[0].Cost = gravcost;

  costlimit = totgravcost / (TOPNODEFACTOR * multipledomains * NTask);
  countlimit = totpartcount / (TOPNODEFACTOR * multipledomains * NTask);

  errflag = domain_check_for_local_refine(0, countlimit, costlimit);

  myfree(mp);

  MPI_Allreduce(&errflag, &errsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(errsum)
    {
      if(ThisTask == 0) printf("We are out of Topnodes. We'll try to repeat with a higher value than All.TopNodeAllocFactor=%g\n", All.TopNodeAllocFactor);
      return errsum;
    }


  /* we now need to exchange tree parts and combine them as needed */

  if(NTask == (1 << PTask))	/* the following algoritm only works for power of 2 */
    {
      for(ngrp = 1, errflag = 0; ngrp < (1 << PTask); ngrp <<= 1)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      /* inform each other about the length of the trees */
	      MPI_Sendrecv(&NTopnodes, 1, MPI_INT, recvTask, TAG_GRAV_A,
			   &ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);


	      topNodes_import =
		(struct local_topnode_data *) mymalloc("topNodes_import",
						       IMAX(ntopnodes_import,
							    NTopnodes) * sizeof(struct local_topnode_data));

	      /* exchange the trees */
	      MPI_Sendrecv(topNodes,
			   NTopnodes * sizeof(struct local_topnode_data), MPI_BYTE,
			   recvTask, TAG_GRAV_B,
			   topNodes_import,
			   ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
			   recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

	      if(sendTask > recvTask)	/* swap the two trees so that result will be equal on all cpus */
		{
		  topNodes_temp =
		    (struct local_topnode_data *) mymalloc("topNodes_temp",
							   NTopnodes * sizeof(struct local_topnode_data));
		  memcpy(topNodes_temp, topNodes, NTopnodes * sizeof(struct local_topnode_data));
		  memcpy(topNodes, topNodes_import, ntopnodes_import * sizeof(struct local_topnode_data));
		  memcpy(topNodes_import, topNodes_temp, NTopnodes * sizeof(struct local_topnode_data));
		  myfree(topNodes_temp);
		  i = NTopnodes;
		  NTopnodes = ntopnodes_import;
		  ntopnodes_import = i;
		}


	      if(errflag == 0)
		{
		  if((NTopnodes + ntopnodes_import) <= MaxTopNodes)
		    {
		      errflag += domain_insertnode(topNodes, topNodes_import, 0, 0);
		    }
		  else
		    {
		      errflag = 1;
		    }
		}

	      myfree(topNodes_import);
	    }
	}
    }
  else
    {
      errflag = domain_recursively_combine_topTree(0, NTask);
    }

  MPI_Allreduce(&errflag, &errsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(errsum) {if(ThisTask == 0) printf("Can't combine trees due to lack of storage. Will try again.\n"); return errsum;}

  /* now let's see whether we should still append more nodes, based on the estimated cumulative cost/count in each cell */

  PRINT_STATUS(" ..NTopNodes before=%d", NTopnodes);
    
  for(i = 0, errflag = 0; i < NTopnodes; i++)
    {
      if(topNodes[i].Daughter < 0)
	if(topNodes[i].Count > countlimit || topNodes[i].Cost > costlimit)	/* ok, let's add nodes if we can */
	  if(topNodes[i].Size > 1)
	    {
	      if((NTopnodes + 8) <= MaxTopNodes)
		{
		  topNodes[i].Daughter = NTopnodes;

		  for(j = 0; j < 8; j++)
		    {
		      sub = topNodes[i].Daughter + j;
		      topNodes[sub].Size = (topNodes[i].Size >> 3);
		      topNodes[sub].Count = topNodes[i].Count / 8;
		      topNodes[sub].Cost = topNodes[i].Cost / 8;
		      topNodes[sub].Daughter = -1;
		      topNodes[sub].Parent = i;
		      topNodes[sub].StartKey = topNodes[i].StartKey + j * topNodes[sub].Size;
		    }

		  NTopnodes += 8;
		}
	      else
		{
		  errflag = 1;
		  break;
		}
	    }
    }

  MPI_Allreduce(&errflag, &errsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(errsum) {return errsum;}

  PRINT_STATUS(" ..NTopnodes after=%d", NTopnodes);
  /* count toplevel leaves */
  domain_sumCost();

  /* Symmetric (NTopleaves is identical on every rank after the replicated
   * tree-combine + domain_sumCost). Soft bad-stop + status-return feeds the
   * existing nonzero-return -> TopNodeAllocFactor retry path in domain_decompose. */
  if(NTopleaves < multipledomains * NTask) {endrun(90000021); return 1;}

  return 0;
}



void domain_sumCost(void)
{
  int i, n, no;
  float *local_domainWork;
  float *local_domainWorkGas;
  int *local_domainCount;
  int *local_domainCountGas;

  local_domainWork = (float *) mymalloc("local_domainWork", NTopnodes * sizeof(float));
  local_domainWorkGas = (float *) mymalloc("local_domainWorkGas", NTopnodes * sizeof(float));
  local_domainCount = (int *) mymalloc("local_domainCount", NTopnodes * sizeof(int));
  local_domainCountGas = (int *) mymalloc("local_domainCountGas", NTopnodes * sizeof(int));

#if (DOMAIN_TIMEBINS == 1)
  /* Per-timebin cost arrays: use malloc (not mymalloc) because they must persist across
     mymalloc stack frames until after domain_assign_load_or_work_balanced finishes */
  int ntb = NumTimeBinsToBeBalanced;
  if(domainBinGravCost) {free(domainBinGravCost); domainBinGravCost = NULL;}
  if(domainBinHydroCost) {free(domainBinHydroCost); domainBinHydroCost = NULL;}
  domainBinGravCost = (float *) calloc(ntb * NTopnodes, sizeof(float));
  domainBinHydroCost = (float *) calloc(ntb * NTopnodes, sizeof(float));
#endif

  NTopleaves = 0;
  domain_walktoptree(0);

  for(i = 0; i < NTopleaves; i++)
    {
      local_domainWork[i] = 0;
      local_domainWorkGas[i] = 0;
      local_domainCount[i] = 0;
      local_domainCountGas[i] = 0;
    }
#if (DOMAIN_TIMEBINS == 1)
  for(i = 0; i < ntb * NTopleaves; i++) {domainBinGravCost[i] = 0; domainBinHydroCost[i] = 0;}
#endif

    PRINT_STATUS(" ..NTopleaves= %d  NTopnodes=%d (space for %d)", NTopleaves, NTopnodes, MaxTopNodes);

  /* Cost accumulation modes controlled by DOMAIN_TIMEBINS:
     - undefined: unweighted (original GIZMO scheme)
     - DOMAIN_TIMEBINS=0: frequency-weighted by 2^(HighestBin - particleBin)
     - DOMAIN_TIMEBINS=1: per-timebin cost accumulation (Gadget-4 scheme),
       with composite cost for tree splitting and separate per-timebin arrays for assignment */

  /* Macro for the per-particle cost accumulation body, shared between OMP and serial paths */
#if defined(DOMAIN_TIMEBINS) && (DOMAIN_TIMEBINS == 0)
  #define DOMAIN_SUMCOST_PARTICLE_BODY(n, no, wk, wkg, cnt, cntg) \
    { float freq_weight = 1.0f; \
      if(All.HighestOccupiedTimeBin > P[n].TimeBin) { \
          int dbin = All.HighestOccupiedTimeBin - P[n].TimeBin; \
          if(dbin > 20) {dbin = 20;} \
          freq_weight = (float)(1 << dbin); } \
      wk[no] += particle_total_cost[n] * freq_weight; \
      cnt[no] += 1; \
      if(P[n].Type == 0) { \
          if(TimeBinActive[P[n].TimeBin] || UseAllParticles) {wkg[no] += particle_costfactor[n] * freq_weight;} \
          cntg[no] += 1;} }
#elif defined(DOMAIN_TIMEBINS) && (DOMAIN_TIMEBINS == 1)
  #define DOMAIN_SUMCOST_PARTICLE_BODY(n, no, wk, wkg, cnt, cntg) \
    { double gc = (double)particle_total_cost[n]; \
      if(gc <= 0) gc = 1.0; \
      float composite_cost = 0; \
      for(int k_ = 0; k_ < ntb; k_++) { \
          int bin_ = ListOfTimeBinsToBeBalanced[k_]; \
          if(bin_ >= P[n].TimeBin) { \
              float contrib = (float)(GravCostNormFactors[k_] * gc); \
              composite_cost += contrib; \
              my_binGravCost[k_ * NTopleaves + no] += contrib; } \
          if(P[n].Type == 0 && bin_ >= P[n].TimeBin) { \
              float hcontrib = (float)HydroCostNormFactors[k_]; \
              composite_cost += hcontrib; \
              my_binHydroCost[k_ * NTopleaves + no] += hcontrib; } } \
      wk[no] += composite_cost; \
      cnt[no] += 1; \
      if(P[n].Type == 0) { \
          if(TimeBinActive[P[n].TimeBin] || UseAllParticles) {wkg[no] += particle_costfactor[n];} \
          cntg[no] += 1;} }
#else
  #define DOMAIN_SUMCOST_PARTICLE_BODY(n, no, wk, wkg, cnt, cntg) \
    { wk[no] += particle_total_cost[n]; \
      cnt[no] += 1; \
      if(P[n].Type == 0) { \
          if(TimeBinActive[P[n].TimeBin] || UseAllParticles) {wkg[no] += particle_costfactor[n];} \
          cntg[no] += 1;} }
#endif

#ifdef _OPENMP
  #pragma omp parallel
  {
    float *my_domainWork = (float *) calloc(NTopleaves, sizeof(float));
    float *my_domainWorkGas = (float *) calloc(NTopleaves, sizeof(float));
    int *my_domainCount = (int *) calloc(NTopleaves, sizeof(int));
    int *my_domainCountGas = (int *) calloc(NTopleaves, sizeof(int));
#if (DOMAIN_TIMEBINS == 1)
    float *my_binGravCost = (float *) calloc(ntb * NTopleaves, sizeof(float));
    float *my_binHydroCost = (float *) calloc(ntb * NTopleaves, sizeof(float));
#endif

    #pragma omp for schedule(static)
    for(n = 0; n < NumPart; n++)
      {
#ifdef SUBFIND
        if(GrNr >= 0 && P[n].GrNr != GrNr) {continue;}
#endif
        no = domain_toptree_leaf(Key[n], topNodes);
        DOMAIN_SUMCOST_PARTICLE_BODY(n, no, my_domainWork, my_domainWorkGas, my_domainCount, my_domainCountGas)
      }

    #pragma omp critical
    {
      for(i = 0; i < NTopleaves; i++) {
        local_domainWork[i] += my_domainWork[i];
        local_domainWorkGas[i] += my_domainWorkGas[i];
        local_domainCount[i] += my_domainCount[i];
        local_domainCountGas[i] += my_domainCountGas[i];
      }
#if (DOMAIN_TIMEBINS == 1)
      for(i = 0; i < ntb * NTopleaves; i++) {
        domainBinGravCost[i] += my_binGravCost[i];
        domainBinHydroCost[i] += my_binHydroCost[i];
      }
#endif
    }
    free(my_domainWork); free(my_domainWorkGas);
    free(my_domainCount); free(my_domainCountGas);
#if (DOMAIN_TIMEBINS == 1)
    free(my_binGravCost); free(my_binHydroCost);
#endif
  }
#else
  for(n = 0; n < NumPart; n++)
    {
#ifdef SUBFIND
      if(GrNr >= 0 && P[n].GrNr != GrNr) {continue;}
#endif
      no = domain_toptree_leaf(Key[n], topNodes);
#if (DOMAIN_TIMEBINS == 1)
      /* In the serial path, my_binGravCost/my_binHydroCost are just the global arrays */
      float *my_binGravCost = domainBinGravCost;
      float *my_binHydroCost = domainBinHydroCost;
#endif
      DOMAIN_SUMCOST_PARTICLE_BODY(n, no, local_domainWork, local_domainWorkGas, local_domainCount, local_domainCountGas)
    }
#endif
#undef DOMAIN_SUMCOST_PARTICLE_BODY

  MPI_Allreduce(local_domainWork, domainWork, NTopleaves, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_domainWorkGas, domainWorkGas, NTopleaves, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_domainCount, domainCount, NTopleaves, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_domainCountGas, domainCountGas, NTopleaves, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#if (DOMAIN_TIMEBINS == 1)
  /* Allreduce per-timebin cost arrays in-place (they use malloc, not mymalloc) */
  MPI_Allreduce(MPI_IN_PLACE, domainBinGravCost, ntb * NTopleaves, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, domainBinHydroCost, ntb * NTopleaves, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif

  myfree(local_domainCountGas);
  myfree(local_domainCount);
  myfree(local_domainWorkGas);
  myfree(local_domainWork);
}


/*! Coordinate conversion to integer. d is coordinate in double precision. returns coordinate in integer of type peano1D. written as part of arepo code dev, from public arepo code by V Springel. */
peano1D domain_double_to_int(double d)
{
    union
    {
        double d;
        unsigned long long ull;
    } u;
    u.d = d;
    return (peano1D)((u.ull & 0xFFFFFFFFFFFFFllu) >> (52 - BITS_PER_DIMENSION));
}


/*! This routine finds the extent of the global domain grid.
 */
void domain_findExtent(void)
{
  int i, j;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];

  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
      xmin[j] = MAX_REAL_NUMBER;
      xmax[j] = -MAX_REAL_NUMBER;
    }

  for(i = 0; i < NumPart; i++)
    {
#ifdef SUBFIND
      if(GrNr >= 0 && P[i].GrNr != GrNr) {continue;}
#endif

      for(j = 0; j < 3; j++)
	{
	  if(xmin[j] > P[i].Pos[j]) {xmin[j] = P[i].Pos[j];}
	  if(xmax[j] < P[i].Pos[j]) {xmax[j] = P[i].Pos[j];}
	}
    }

  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  len = 0;
  for(j = 0; j < 3; j++) {if(xmax_glob[j] - xmin_glob[j] > len) {len = xmax_glob[j] - xmin_glob[j];}}

  len *= 1.001;

  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }
#ifdef RANDOMIZE_GRAVTREE // double the size of the root node and pick a random offset for its center, so that forcetree errors get decorrelated each time the tree is rebuilt
  double dx[3]; 
  if(ThisTask == 0) { for(j = 0; j < 3; j++) {dx[j] = len * (get_random_number((MyIDType) (All.NumCurrentTiStep) + j) - 0.5);}}
  MPI_Bcast(dx, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for(j=0; j<3; j++) {
      DomainCenter[j] += dx[j];
      DomainCorner[j] = DomainCenter[j] - len;
  }
  len *= 2;
#endif  
  DomainLen = len;
  DomainFac = 1.0 / len * (((peanokey) 1) << (BITS_PER_DIMENSION));
}




void domain_add_cost(struct local_topnode_data *treeA, int noA, long long count, double cost, double gascost)
{
  int i, sub;
  long long countA, countB;

  countB = count / 8;
  countA = count - 7 * countB;

  cost = cost / 8;
  gascost = gascost / 8;

  for(i = 0; i < 8; i++)
    {
      sub = treeA[noA].Daughter + i;
      if(i == 0) {count = countA;} else {count = countB;}

      treeA[sub].Count += count;
      treeA[sub].Cost += cost;
      treeA[sub].GasCost += gascost;

      if(treeA[sub].Daughter >= 0) {domain_add_cost(treeA, sub, count, cost, gascost);}
    }
}


/* Pure local recursive combine of two replicated top-trees -- no MPI here, so it
 * runs identically on every rank. Returns nonzero (status) on a capacity/invariant
 * failure so the caller can fold it into errflag -> MPI_Allreduce(errsum) ->
 * TopNodeAllocFactor retry instead of a hard abort. */
int domain_insertnode(struct local_topnode_data *treeA, struct local_topnode_data *treeB, int noA, int noB)
{
  int j, sub;
  long long count, countA, countB;
  double cost, costA, costB;

  if(treeB[noB].Size < treeA[noA].Size)
    {
      if(treeA[noA].Daughter < 0)
	{
	  if((NTopnodes + 8) <= MaxTopNodes)
	    {
	      count = treeA[noA].Count - treeB[treeB[noB].Parent].Count;
	      countB = count / 8;
	      countA = count - 7 * countB;

	      cost = treeA[noA].Cost - treeB[treeB[noB].Parent].Cost;
	      costB = cost / 8;
	      costA = cost - 7 * costB;

	      treeA[noA].Daughter = NTopnodes;
	      for(j = 0; j < 8; j++)
		{
		  if(j == 0)
		    {
		      count = countA;
		      cost = costA;
		    }
		  else
		    {
		      count = countB;
		      cost = costB;
		    }

		  sub = treeA[noA].Daughter + j;
		  topNodes[sub].Size = (treeA[noA].Size >> 3);
		  topNodes[sub].Count = count;
		  topNodes[sub].Cost = cost;
		  topNodes[sub].Daughter = -1;
		  topNodes[sub].Parent = noA;
		  topNodes[sub].StartKey = treeA[noA].StartKey + j * treeA[sub].Size;
		}
	      NTopnodes += 8;
	    }
	  else
	    {endrun(90000022); return 1;} /* out of TopNodes: soft bad-stop + status-return (skips the Daughter<0 deref below); folds into caller errflag -> TopNodeAllocFactor retry */
	}

      sub = treeA[noA].Daughter + (treeB[noB].StartKey - treeA[noA].StartKey) / (treeA[noA].Size >> 3);
      return domain_insertnode(treeA, treeB, sub, noB);
    }
  else if(treeB[noB].Size == treeA[noA].Size)
    {
      treeA[noA].Count += treeB[noB].Count;
      treeA[noA].Cost += treeB[noB].Cost;

      if(treeB[noB].Daughter >= 0)
	{
	  for(j = 0; j < 8; j++)
	    {
	      sub = treeB[noB].Daughter + j;
	      if(domain_insertnode(treeA, treeB, noA, sub)) {return 1;}
	    }
	}
      else
	{
	  if(treeA[noA].Daughter >= 0)
	    domain_add_cost(treeA, noA, treeB[noB].Count, treeB[noB].Cost, treeB[noB].GasCost);
	}
    }
  else
    {endrun(90000023); return 1;} /* unexpected node size (symmetric invariant): soft bad-stop + status-return -> caller errflag -> TopNodeAllocFactor retry */

  return 0;
}



static void parallel_sort_phdata_domain(struct peano_hilbert_data *arr, size_t n)
{
  auto cmp = [](const peano_hilbert_data &a, const peano_hilbert_data &b) { return a.key < b.key; };
  const size_t CUTOFF = 10000;
  if(n <= CUTOFF) { std::sort(arr, arr + n, cmp); return; }
  size_t mid = n / 2;
#ifdef _OPENMP
  #pragma omp task shared(arr) if(n > CUTOFF)
#endif
  parallel_sort_phdata_domain(arr, mid);
#ifdef _OPENMP
  #pragma omp task shared(arr) if(n > CUTOFF)
#endif
  parallel_sort_phdata_domain(arr + mid, n - mid);
#ifdef _OPENMP
  #pragma omp taskwait
#endif
  std::inplace_merge(arr, arr + mid, arr + n, cmp);
}

void mysort_domain(void *b, size_t n, size_t s)
{
  struct peano_hilbert_data *arr = (struct peano_hilbert_data *)b;
#ifdef _OPENMP
  #pragma omp parallel
  {
    #pragma omp single
    parallel_sort_phdata_domain(arr, n);
  }
#else
  parallel_sort_phdata_domain(arr, n);
#endif
}
