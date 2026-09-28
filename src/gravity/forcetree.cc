#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "forcetree.h"               /* GIZMO_EWALD_EN + Ewald-table accessor decls */
#include "gravtree_force_kernel.h"   /* shared CPU/GPU accepted-source contribution physics (SSOT) */
#include "gravtree_moment_kernel.h"  /* shared node moment/payload construction physics (SSOT); plain primitives only here */
#include "gravtree_moment_sources.h" /* shared per-particle RT/sink/CR source-input gates (SSOT) */
#include "gravtree_ewald.h"          /* shared CPU/GPU Ewald image-correction trilinear interp (SSOT) */
#include "pm_highres_region.h"       /* pmforce_is_particle_high_res SSOT (device-callable) */
#include "let_data.h"   /* LET wire format + per-rank payload structs */
#include "../mesh/gpu_neighbor_list.h" /* gizmo_mark_kernel_radius_dirty_indices */
#include "../mesh/nlr_radius_policy.h" /* SSOT helper for force_hmax_per_type_particle_radius */
#include "../core/timestep_functions.h" /* Hermite source eligibility + prediction, inline (a cross-TU call here costs ~16% of the walk) */
#ifdef SUBFIND
#include "../structure/subfind/subfind.h"
#endif
#include "gpu_gravity_tree.h"
#include "gpu_peano_walk.h"
#include "gpu_topology_build.h"
#include "gpu_topology_finalize.h"
#include "gpu_pseudo_update.h"

/*! \file forcetree.c
 *  \brief gravitational tree and code for Ewald correction
 *
 *  This file contains the computation of the gravitational force by means
 *  of a tree. The type of tree implemented is a geometrical oct-tree,
 *  starting from a cube encompassing all particles. This cube is
 *  automatically found in the domain decomposition, which also splits up
 *  the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. Tree nodes can
 *  be dynamically updated in drift/kick operations to avoid having to
 *  reconstruct the tree every timestep.
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


/* Compute the force-softening kernel radius for one particle.
 *
 * Single source of truth for the per-particle softening logic.  All call sites
 * (CPU walk, GPU walk, tree-build split-scale, etc.) read the cached result via
 * ForceSoftening_KernelRadius() / gpu_force_softening_kernel_radius(); this
 * routine is the only place that performs the actual computation.
 *
 * The cache is refreshed once per gravity_tree() call by compute_all_force_softening()
 * (active-particle loop) and seeded over all particles at startup in init.c.  Inputs
 * (KernelRadius, AGS_KernelRadius, tidal_tensor_mag_prev, StarParticleEffectiveSize)
 * only mutate when the particle is active, so cached values for inactive particles
 * remain correct between active steps. */

double compute_force_softening_kernel_radius(int p)
{
    double soft = All.ForceSoftening[P[p].Type];   /* default if no branch below is active */

#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
    if(P[p].Type == 4) {soft = P[p].StarParticleEffectiveSize;} // this variable is defined in force softening terms
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORALL)
    if((1 << P[p].Type) & (ADAPTIVE_GRAVSOFT_FORALL)) {soft = P[p].AGS_KernelRadius;}
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(SELFGRAVITY_OFF) /* softening scale still appears in timestep criterion for problems without self-gravity, so set it adaptively */
#ifdef ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT
    if(P[p].Type == 0) {soft = DMIN(P[p].KernelRadius, ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT/All.cf_atime);}
#else
    if(P[p].Type == 0) {soft = P[p].KernelRadius;}
#endif
#endif

#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    if(P[p].Type == 4) {soft = All.ForceSoftening[P[p].Type] * DMIN(100., DMAX(1., pow(P[p].Mass*UNIT_MASS_IN_SOLAR/100. , 0.33)));}
#endif

#if defined(ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION) /* still playing with criterion below, highly experimental for now */
    if((1 << P[p].Type) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)) {if((P[p].tidal_tensor_mag_prev>0) && (All.Time>All.TimeBegin)) {soft = DMIN(1.e2*All.ForceSoftening[P[p].Type] , DMAX(All.ForceSoftening[P[p].Type] , All.ForceSoftening[P[p].Type] + 1.25 * pow( (All.DesNumNgb * All.G * P[p].Mass / P[p].tidal_tensor_mag_prev) , 1./3. )));} else {soft = 100.*All.ForceSoftening[P[p].Type];}}
#endif

#ifdef PMGRID
    if(All.Asmth[0] > 0) {soft = DMIN(soft, 0.5*All.Asmth[0]);} /* maximum softening before PMGRID cutoff invalidates softening */
#ifdef PM_PLACEHIGHRESREGION
    if(pmforce_is_particle_high_res(P[p].Type, P[p].Pos)) {soft = DMIN(soft, 0.5*All.Asmth[1]);}
#endif
#endif

    return soft;
}

/* Public accessor: returns the cached value populated by compute_all_force_softening().
 * Both CPU and GPU walks read this same value, so the softening logic above lives in
 * exactly one place. The body is ForceSoftening_KernelRadius_P in core/proto.h, which
 * device code calls directly. */
double ForceSoftening_KernelRadius(int p)
{
    return ForceSoftening_KernelRadius_P(p, P);
}

/* Refresh the per-particle ForceSoftening cache.  Called from gravity_tree() at
 * the start of every walk dispatch (active particles only) and from init.c during
 * startup (all particles).  Inputs to compute_force_softening_kernel_radius() only
 * change for active particles within a timestep, so an active-particle pass is
 * sufficient for steady-state operation; the init pass seeds inactive particles
 * loaded from the IC file or spawned mid-run. */
void compute_all_force_softening(int mode)
{
    /* mode = 0 : active particles only (FirstActiveParticle list)
     * mode = 1 : every particle in [0, NumPart) -- used at startup / after restart  */
    if(mode == 1)
    {
        int i;
#pragma omp parallel for schedule(static)
        for(i = 0; i < NumPart; i++) {P[i].ForceSoftening = compute_force_softening_kernel_radius(i);}
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int ii = 0; ii < (int)ActiveParticleList.size(); ii++)
        {
            int i = ActiveParticleList[ii];
            P[i].ForceSoftening = compute_force_softening_kernel_radius(i);
        }
    }
}



/*! auxiliary variable used to set-up non-recursive walk */
static int last;

/* NEIGHBORS_MUST_BE_COMPUTED_EXPLICITLY_IN_FORCETREE is defined globally in
 * precompiler_logic.h (via allvars.h) so the CPU and GPU tree walks share the
 * same leaf-opening criterion. */

/*! length of look-up table for short-range force kernel in TreePM algorithm */
#define NTAB GRAVTREE_SHORTRANGE_NTAB   /* table length owned by gravtree_force_kernel.h (shared with the GPU walk) */
/*! variables for short-range lookup table.  Non-static so the GPU gravity
 *  walk in gpu_gravtree.cc can read them via extern declarations.  Sized at
 *  NTAB floats = 4 KB each — fine to leave in host memory on Kokkos OMP;
 *  for true device offload they will need mirroring. */
float shortrange_table[NTAB], shortrange_table_potential[NTAB];
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
/* Non-static so the GPU walk in gpu_gravtree.cc can read this via extern (mirrors
 * the shortrange_table / shortrange_table_potential pattern above). */
float shortrange_table_tidal[NTAB];
#endif
/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;
static int tree_allocated_flag = 0;

/*! Whether tree storage is currently held. The flag is the only record of that, so callers who
 *  must know (the particle-storage resize, which may not release capacity under a standing tree)
 *  ask here rather than inferring it from MaxNodes or a live pointer. */
int force_tree_is_allocated(void) {return tree_allocated_flag;}

#ifdef BOX_PERIODIC
/*! Size of 3D look-up table for Ewald correction force */
#define EN  64
/*! 3D look-up table for Ewald correction to force and potential. Only one octant is stored, the rest constructed by using the symmetry of the problem */
static MyFloat fcorrx[EN + 1][EN + 1][EN + 1];
static MyFloat fcorry[EN + 1][EN + 1][EN + 1];
static MyFloat fcorrz[EN + 1][EN + 1][EN + 1];
static MyFloat potcorr[EN + 1][EN + 1][EN + 1];
static double fac_intp;
#if !defined(GRAVITY_NOT_PERIODIC)
/* Accessor for gpu_gravtree.cc — see forcetree.h. */
void gizmo_get_ewald_tables(const MyFloat **fcorrx_out, const MyFloat **fcorry_out,
                            const MyFloat **fcorrz_out, const MyFloat **potcorr_out,
                            double *fac_intp_out)
{
    *fcorrx_out  = &fcorrx[0][0][0];
    *fcorry_out  = &fcorry[0][0][0];
    *fcorrz_out  = &fcorrz[0][0][0];
    *potcorr_out = &potcorr[0][0][0];
    *fac_intp_out = fac_intp;
}
#endif
#endif


/* Gravity box-wrap policy is centralized in the shared helpers (single source used by both
 * the CPU tree walk here and the GPU walks). The gravity gate (wrap only for periodic box
 * AND periodic gravity, else no-op/abs) lives inside gravity_box_distance.h. The GRAVITY_*
 * macros below are thin wrappers so the existing call sites are unchanged. */
#include "gravity_box_distance.h"
#include "gravtree_opening.h"
#define GRAVITY_NEAREST_XYZ(x,y,z,sign)             gravity_box_nearest_image(x,y,z,sign)
#define GRAVITY_NGB_PERIODIC_BOX_LONG_X(x,y,z,sign) gravity_box_long_abs_x(x,y,z,sign)
#define GRAVITY_NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) gravity_box_long_abs_y(x,y,z,sign)
#define GRAVITY_NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign) gravity_box_long_abs_z(x,y,z,sign)

/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
/*! Mode B per-type hmax host-side re-seed.
 *
 *  Why this exists: gpu_moment_refresh() writes the scalar Extnodes[no].hmax
 *  to AoS but not the per-type bands Extnodes[no].hmax_per_type[].  The GPU
 *  SoA intentionally does not carry per-type bands -- their only consumer is
 *  this host-side Mode B walker (mesh/mode_b_local_walker.cc) -- so the GPU
 *  moment path bypasses the host moment loop that would otherwise seed them.
 *  Without this pass every full force_treebuild and every
 *  force_refresh_node_moments would leave the bands at zero, and Mode B's
 *  SYMMETRIC walker reading zero bands would over-prune (collapse to ONEWAY)
 *  -- the neighbour-set loss observed on fire_m11i.
 *
 *  Behavior: zero all internal-node bands, leaf-seed each particle's
 *  conservative radius into Father[i]'s band, then bottom-up max-over-children
 *  (cheap host loop ~O(NumPart + Nnodes)).  Caller-restriction: must run AFTER
 *  gpu_moment_refresh has populated Father[] / Nodes[].u.d.father, since the
 *  bottom-up step walks the father chain.
 */
/* Conservative per-particle radius (node-prune upper bound) for
 * Extnodes[no].hmax_per_type[Type] band seeding. See forcetree.h docstring.
 * Used by every site that grows or seeds a per-type band:
 * force_refresh_hmax_per_type_host, force_update_node_recursive,
 * force_refresh_node_moments, force_update_hmax (forcetree_update.cc),
 * force_add_element_to_tree.
 *
 * SSOT: routes through nlr_particle_symmetric_radius_capped — no
 * AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE branch here.  MODE_B_RADIUS_ALL_SOURCES
 * is the conservative union (KernelRadius / AGS_KernelRadius / ForceSoftening
 * across all types).  Capping at All.MaxKernelRadius applies to kernel radii
 * only; ForceSoftening is uncapped so the band dominates the leaf-policy reach
 * even when FS > All.MaxKernelRadius. */
double force_hmax_per_type_particle_radius(int i)
{
    return nlr_particle_symmetric_radius_capped(P[i],
                                                MODE_B_RADIUS_ALL_SOURCES,
                                                (double)All.MaxKernelRadius);
}

static void force_refresh_hmax_per_type_host(int Numnodestree)
{
    /* Step 1: zero per-type bands in all internal nodes. */
    for(int no = All.TreeNodeIndexBase; no < All.TreeNodeIndexBase + Numnodestree; no++) {
        for(int t = 0; t < 6; t++) Extnodes[no].hmax_per_type[t] = 0;
    }
    /* Step 2: leaf seed — Father[i] only (per-particle to its immediate
     * parent), conservative across every leaf-policy-selectable source.
     * Stops at the tree's own particle slots rather than at NumPart: this also runs as a refresh on
     * a tree that is already standing, and by then P[] can reach further than Father[] does, since
     * the particle capacity may be raised mid-step.  Imported ghosts, which also sit above the local
     * particles, are skipped by the test just below instead: every slot the tree does not contain
     * carries -1, which is what a ghost's slot holds. */
    const int nseed = (NumPart < All.TreeParticleSlots) ? NumPart : All.TreeParticleSlots;
    for(int i = 0; i < nseed; i++) {
        int no = Father[i];
        if(no < 0) continue;
        struct particle_data *pa = &P[i];
        if(pa->Mass <= 0) continue;
        double htmp = force_hmax_per_type_particle_radius(i);
        int t = (int)pa->Type;
        if(htmp > Extnodes[no].hmax_per_type[t]) Extnodes[no].hmax_per_type[t] = (MyFloat)htmp;
    }
    /* Step 3: bottom-up max-over-children via father chain. Children always
     * allocated at higher indices than parents, so reverse iteration gives
     * children-before-parent order without an explicit DAG sort. */
    for(int no = All.TreeNodeIndexBase + Numnodestree - 1; no >= All.TreeNodeIndexBase; no--) {
        int father = Nodes[no].u.d.father;
        if(father < All.TreeNodeIndexBase || father >= All.TreeNodeIndexBase + Numnodestree) continue;
        for(int t = 0; t < 6; t++) {
            if(Extnodes[no].hmax_per_type[t] > Extnodes[father].hmax_per_type[t]) {
                Extnodes[father].hmax_per_type[t] = Extnodes[no].hmax_per_type[t];
            }
        }
    }
}



/* Gravity-tree freshness generations (see forcetree.h). Plain host counters,
 * SSOT in this TU; force_update_hmax (forcetree_update.cc) bumps the hmax one
 * via force_bump_hmax_refresh_generation(). */
static long g_force_treebuild_generation = 0;
static long g_force_hmax_refresh_generation = 0;
long force_treebuild_generation(void)        { return g_force_treebuild_generation; }
long force_hmax_refresh_generation(void)      { return g_force_hmax_refresh_generation; }
void force_bump_hmax_refresh_generation(void) { g_force_hmax_refresh_generation++; }

/* Whether the standing tree's Father[] links still describe the particles they were built for.  A
 * whole-tree rebuild that happens without a domain decomposition on the same step needs them: a
 * particle that has drifted into a top-leaf owned by another rank is kept under the leaf it hung from
 * here, and these links are the only record of that.  Set by a successful whole-tree build, dropped
 * when the tree goes away or when the particle array is reordered underneath it. */
static int g_force_global_topology_valid = 0;
int  force_tree_global_topology_valid(void) {return g_force_global_topology_valid;}
void force_tree_invalidate_global_topology(void) {g_force_global_topology_valid = 0;}

/*! Let a particle keep its parent when it changes slots.  Re-sequencing moves particles between
 *  slots without moving what the tree says about them, so without this the record stops describing
 *  the particles and the next rebuild has to fall back to restoring ownership instead -- which is a
 *  domain repartition, the very thing the tree's own cadence exists to avoid.  Only the parent moves:
 *  the walk's own links are not repaired, because repairing them means a full tree traversal per
 *  swap.  What makes that safe is that a re-sequencing which moved slots asks for a rebuild, and
 *  until that rebuild happens nothing may combine the carried parents with the old walk links --
 *  force_add_element_to_tree declines to, which is the one operation that reads both together.
 *
 *  The standing contract for creation: a particle that is CREATED must either be inserted into the
 *  live tree by force_add_element_to_tree, or the record must be dropped.  Every creation path
 *  satisfies one or the other -- the wind spawn and star formation insert, and the refinement pass
 *  runs only inside a decomposition, which frees the tree.  A future creation path that does neither
 *  would leave a parent belonging to some other particle, and the rebuild would attach this one
 *  wherever that parent happens to be. */
void force_tree_swap_attachment_slots(int i, int j)
{
    if(!g_force_global_topology_valid) {return;}
    if(i < 0 || j < 0 || i >= All.TreeParticleSlots || j >= All.TreeParticleSlots || !Father)
    {
        /* The attachment could not be moved, so the record no longer describes the particles.
           Saying so is the whole job: returning quietly would leave it advertised as valid and the
           next rebuild would trust it. */
        force_tree_invalidate_global_topology();
        return;
    }
    const int no = Father[i]; Father[i] = Father[j]; Father[j] = no;
}

/* How much the foreign-node index ceiling is padded above measured demand when it grows.  It is
 * a large fraction because the ceiling is cheap: it buys index range and its Nextnode ints, not
 * the nodes, which are allocated to the exact import.  The pad is what stops a build that grows
 * slightly from forcing another full rebuild. */
static const double LET_FOREIGN_PAD_FRACTION = 0.5;

/* Why a tree build failed, as seen by the retry in force_treebuild.  A non-negative return is
 * the node count.  Only OUT_OF_NODES may be answered by growing the node arena: every other
 * cause is indifferent to how many nodes are reserved, so retrying it larger cannot succeed and
 * instead ratchets MaxNodes -- and MaxForeignNodes with it, since that is derived from it --
 * until the arena itself cannot be allocated.  That is how one build failure became a 409x
 * growth and a dead run.  The outer loop reduces this with MPI_Allreduce(MIN), so HARD_FAILURE
 * on any rank outranks OUT_OF_NODES everywhere and all ranks stop together. */
#define FORCE_TREEBUILD_OUT_OF_NODES  (-1)
/* NOT -2: that is FORCE_TREE_NEEDS_OWNERSHIP_RESTORE, which force_treebuild returns to its own
 * callers and which asks them to restore ownership and build again.  The two are consumed in
 * different places today, but sharing a value would let a later refactor answer a fatal
 * allocation failure with a decomposition and a retry. */
#define FORCE_TREEBUILD_HARD_FAILURE  (-3)


int force_treebuild(int npart, struct unbind_data *mp)
{
    int flag;
    /* Adaptive LET foreign-arena retry: force_treebuild owns the rebuild loop (same
     * idiom as the TreeAllocFactor-overflow retry below). On a retryable LET overflow
     * we grow the adaptive floor and rebuild the whole tree from scratch -- the foreign
     * arena is contiguous inside Nodes_base, so growing it means reallocating + rebuilding.
     * One rebuild per ratchet event suffices (the in-call domain is fixed, so the rebuilt
     * LET needs the same capacity and the pad covers it); the bound is a backstop. */
    int let_retry = 0;
    const int LET_MAX_RETRY = 3;

    /* Retained attachment.  Particles drift across top-leaf boundaries between decompositions, so a
     * whole-tree rebuild that happens without one finds some of this rank's particles falling
     * geometrically in top-leaves another rank owns.  Bucketing them there detaches them: the
     * pseudo-particle exchange below overwrites those nodes, and the subtree holding them is left
     * unreachable from the root, so their mass enters no rank's multipole moments and the forces on
     * everything else are wrong by it.  They are kept under the top-leaf they hung from in the
     * standing tree instead.
     *
     * This runs here, before the loop, because it is the last point at which the standing tree still
     * exists: a retry frees the tree, taking Father[] and DomainNodeIndex with it.  Its result is
     * reused by every attempt -- positions, the top tree and the retained attachment do not change
     * while the build retries for a larger arena. */
    long crossed_local = 0, unrecovered_local = 0;
    /* Whether THIS rank is building its whole local tree.  The test is rank-local -- a group tree the
       halo finder builds can happen to hold as many members as a rank has particles -- so the
       reduction below carries it too, and a build that is not the whole tree everywhere keeps nothing
       from this stage. */
    const int whole_tree_local = (mp == NULL && npart == NumPart) ? 1 : 0;
    if(whole_tree_local)
    {
        if(gpu_topology_prepare_retained_attachment(npart, force_tree_global_topology_valid(),
                                                    &crossed_local, &unrecovered_local) != 0)
        {
            printf("force_treebuild: task %d could not prepare retained top-leaf attachments\n", ThisTask);
            endrun(91564);
        }
    }
    else {gpu_topology_forget_prepared();}

    /* Three globally-reduced answers, in one collective every rank reaches.  "Restore ownership" is a
       request the caller can honour with a repartition, and is raised only when there is something to
       restore: a build with no crossers needs nothing from the standing tree, which is why an ordinary
       build after a decomposition is unaffected.  "Unrecovered" is a valid standing tree that
       nevertheless cannot place a particle -- the tree and the particles disagree, and there is nothing
       to fall back to. */
    long counts_local[3], counts_any[3];
    counts_local[0] = (crossed_local > 0 && !force_tree_global_topology_valid()) ? 1 : 0;
    counts_local[1] = unrecovered_local;
    counts_local[2] = whole_tree_local ? 0 : 1;
    MPI_Allreduce(counts_local, counts_any, 3, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
    if(counts_any[2] > 0)
    {
        gpu_topology_forget_prepared();   /* not the whole tree somewhere: this build keeps nothing */
    }
    else
    {
        if(counts_any[1] > 0 && unrecovered_local > 0)
        {
            /* Said only by the ranks that actually hold one; the others would print zeros and bury it.
               This one is worth reporting on its own: the standing tree was declared usable and then
               could not place a particle, so something about that declaration is wrong even though the
               state itself is repairable below. */
            printf("force_treebuild: task %d holds %ld particles whose standing-tree attachment cannot be "
                   "used -- it names no top-leaf this rank owns, or a top-leaf whose node lies outside the "
                   "tree -- from %ld that have crossed a top-leaf boundary. The tree and the particles "
                   "disagree.\n", ThisTask, unrecovered_local, crossed_local);
        }
        if(counts_any[0] > 0 || counts_any[1] > 0)
        {
            /* Both states have the same repair and it already exists: recomputing ownership from the
               particles' actual positions leaves every particle in a top-leaf its own rank owns, after
               which the build needs nothing from the standing tree at all.  Nothing has been freed or
               rebuilt yet, so the caller can do that and ask again; it stops only if the repair itself
               does not take. */
            gpu_topology_forget_prepared();
            return FORCE_TREE_NEEDS_OWNERSHIP_RESTORE;
        }
    }

let_build_attempt:
    /* reset force_add_element insertion counter at each full build. */
    ForceAddElementToTree_CallsSinceBuild = 0;
    /* How many times the node arena may be GROWN for one build; the build that fails after the
     * last of those stops instead of growing again.  Each growth multiplies the arena by 1.15, so
     * this allows 5.4x, far more than a genuine shortfall has needed, while keeping the loop
     * finite.  Without a bound, a cause the arena cannot fix ratchets until the arena itself
     * cannot be allocated: 43 growths once took MaxNodes from 2.1e6 to 8.6e8, a 2.1 TB request. */
    const int TREEBUILD_MAX_NODE_RETRY = 12;
    int treebuild_node_retry = 0;
    do
    {
        Numnodestree = force_treebuild_single(npart, mp);
        MPI_Allreduce(&Numnodestree, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(flag < 0 && flag != FORCE_TREEBUILD_OUT_OF_NODES)
        {
            /* A cause that reserving more nodes cannot address.  The rank that hit it has already
             * said which one; stop here rather than growing the arena against it.  Symmetric:
             * flag comes from the Allreduce, so every rank takes this branch together. */
            if(ThisTask == 0)
            {
                printf("The tree could not be built, for a reason that more tree nodes would not fix (see the failure\n"
                       "reported above). Not growing the node arena. Stopping.\n");
                fflush(stdout);
            }
            endrun(90000060);
            /* endrun only REQUESTS a stop; drain it here, because the finalize stages below index
             * the node array with Numnodestree, which is negative on a failed build.  Symmetric:
             * flag is the reduced value, so every rank reaches this branch with the same verdict. */
            gizmo_exit_bad_stop_if_requested("gravtree:treebuild");
            break;
        }
        if(flag == FORCE_TREEBUILD_OUT_OF_NODES && treebuild_node_retry >= TREEBUILD_MAX_NODE_RETRY)
        {
            if(ThisTask == 0)
            {
                printf("The tree still ran out of nodes after %d successive growths of the node arena, each 1.15x the\n"
                       "last (MaxNodes=%d now). The demand is outrunning the growth, so this is not a shortfall that a\n"
                       "further retry settles; running on more ranks or nodes reduces what each one has to hold. Stopping.\n",
                       treebuild_node_retry, MaxNodes);
                fflush(stdout);
            }
            endrun(90000061);
            gizmo_exit_bad_stop_if_requested("gravtree:treebuild");   /* as above: stop before the finalize stages read a negative node count */
            break;
        }
        if(flag == FORCE_TREEBUILD_OUT_OF_NODES)
        {
            treebuild_node_retry++;
            /* Grow the node arena by the same factor the ratchet applies, measured from the
               current allocation rather than re-derived from All.MaxPart: that keeps the retry
               on whatever basis this tree was sized from (the domain's local-particle cap for
               the gravity tree, All.MaxPart for the trees subfind builds over its own particle
               sets), so the ratcheted factor converges on the size the next build will ask for. */
            int maxnodes_grown = (int) (1.15 * (double) (MaxNodes - NTopnodes)) + NTopnodes;
            if(maxnodes_grown <= MaxNodes) {maxnodes_grown = MaxNodes + 1;}   /* integer truncation makes the 1.15x a no-op below 7 nodes; always make progress */
            force_treefree();
            if(ThisTask == 0) {printf("Increasing TreeAllocFactor=%g", All.TreeAllocFactor);}
            All.TreeAllocFactor *= 1.15;
            if(ThisTask == 0) {printf(" new value=%g\n", All.TreeAllocFactor);}
            force_treeallocate(maxnodes_grown, All.MaxPartExpandable);
            /* drain a tree-alloc UVM OOM before force_treebuild_single re-runs on a NULL-based
             * tree. Symmetric: all ranks enter this block together (flag from Allreduce(MIN)). */
            gizmo_exit_bad_stop_if_requested("gravtree:treeallocate");
        }
    }
    while(flag == FORCE_TREEBUILD_OUT_OF_NODES);
    /* GPU finalize stage replaces force_update_node_recursive's
     * sibling/father/Father[] outputs.  Order matters:
     *   1. finalize_father: writes soa->father for all internal nodes
     *      (covers topnodes — emit_bfs only set inside-topleaf), and
     *      writes Father[i] for every particle child.  Must run before
     *      moment_refresh (which reads soa->father in its dependency walk).
     *   2. finalize_sibling: writes soa->sibling for all internal nodes.
     *      Must run before nextnode_thread (which reads soa->sibling).
     *   3. moment_refresh: writes moments + Extnodes/N_part/maxsoft/bitflags.
     *   4. nextnode_thread: writes nextnode + Nextnode[] from suns_backup.
     *   5. writeback_d_to_aos: pushes soa->sibling/father into AoS u.d for
     *      legacy CPU walks.  Clobbers u.suns via union, but suns_backup
     *      in SoA is the truth.  Runs last so prior steps reading SoA see
     *      consistent state.  topnode-range center/len was already
     *      pulled into SoA by gpu_nextnode_backup_suns inside
     *      force_treebuild_single. */
    /* The GPU tree-finalize steps below are rank-local; on failure they set
     * a soft bad-stop and fall through the matched, topology-driven pseudo/
     * LET collectives, which drain at the gravtree:after_treebuild poll
     * before the GPU gravity walk reads any moments. */
    if(gpu_topology_finalize_father(Numnodestree)  != 0) {endrun(90000065);}
    /* Every further member of a multi-particle leaf takes the same father as its head. Must precede
     * the moment pass, which accumulates through Father[] one particle at a time. */
    if(gpu_leaf_chain_assign_fathers() != 0) {endrun(90000081);}
    if(gpu_topology_finalize_sibling(Numnodestree) != 0) {endrun(90000066);}
    /* Cover the retained particles' true positions.  gpu_topology_finalize_father above is what
     * establishes the particle Father[] links this walks, and the pseudo-particle exchange further
     * down is what carries a grown top-leaf length to the other ranks, so it belongs between them. */
    if(gpu_topology_grow_retained_paths() != 0) {endrun(90000090);}
    /* GPU kernel resets GravCost + ephemeral fields for all
     * nodes.  On the CPU path FUNR does this work inline; on the GPU path
     * FUNR is retired so the kernel takes its place.  Replaces a
     * host loop over Numnodestree -- the worst sparse-active scaling. */
    if(gpu_node_reset_ephemeral(Numnodestree) != 0) {endrun(90000067);}
    if(gpu_moment_refresh(-1) != 0) {endrun(90000068);}
    if(gpu_nextnode_thread() != 0) {endrun(90000069);}
    /* Link each multi-particle leaf's members now that its successor is known. ⛔ Must run BEFORE the
     * writeback and the tree export below: the export enumerates a leaf by following these links, so
     * exporting an unlinked leaf would ship its first particle and silently omit the rest. */
    if(gpu_leaf_chain_materialize() != 0) {endrun(90000079);}
    if(gpu_topology_writeback_d_to_aos(Numnodestree) != 0) {endrun(90000070);}
    /* Mode B: GPU moment refresh writes scalar hmax but not per-type bands;
     * re-seed those host-side now. MUST run AFTER gpu_topology_writeback_d_to_aos
     * because force_refresh_hmax_per_type_host's Step 3 propagation walks via
     * Nodes[no].u.d.father, which is only valid post-writeback (the SoA→AoS
     * writeback overwrites the union slot from the build-time u.suns layout). */
    force_refresh_hmax_per_type_host(Numnodestree);
    /* set TOPLEVEL/INTERNAL_TOPLEVEL/DEPENDS bitflags in SoA
     * (and mirror to AoS for force_exchange_pseudodata / force_treeupdate_pseudos
     * which still run on CPU). */
    if(gpu_force_flag_localnodes() != 0) {endrun(90000071);}
    /* Non-blocking overlap: post the pseudo-data Iallgathervs first, then run the
     * LET MPI round concurrently, then wait/unpack pseudo-data and resum.
     * LET pack reads only LOCAL Nodes/Extnodes (which are already valid from
     * gpu_moment_refresh above) and does not depend on foreign topleaves; so
     * the two MPI exchanges can overlap.  Latency drops from sum to max of
     * the two collectives' wall-times. */
    force_exchange_pseudodata_issue();
    /* LET exchange returns a typed status + (on overflow) the foreign-node capacity the
     * receiver needed. The pseudodata Iallgatherv posted just above must still be
     * completed regardless, so no early return here. */
    long long foreign_needed = 0;
    let_exchange_status_t let_status = let_run_exchange(&foreign_needed);
    int pseudo_status = force_exchange_pseudodata_complete();

    /* Classify globally (worst status wins). A send-buffer malloc failure
     * (LET_PACK_OOM), a malformed exchange (LET_UNPACK_INTERNAL), or a pseudodata
     * failure are HARD -- graceful stop, no retry (a bigger arena cannot fix them).
     * A LET_OVERFLOW_RETRYABLE on any rank triggers an arena grow + full rebuild.
     * Deciding globally keeps every rank on the same path (downstream + retry are
     * collective). */
    int flags_local[3], flags_any[3] = {0, 0, 0};
    flags_local[0] = (let_status == LET_OVERFLOW_RETRYABLE);
    flags_local[1] = (let_status == LET_PACK_OOM) || (let_status == LET_ARENA_SHORT) || (let_status == LET_FOREIGN_STORAGE_SHORT) || (let_status == LET_UNPACK_INTERNAL) || (pseudo_status != 0);
    /* Whether this is a build over the whole particle set, for the integrity check at the end.
       It rides here rather than in its own reduction, and must be reduced UNCONDITIONALLY: the
       predicate is rank-local, and SUBFIND calls force_treebuild(NumPartGroup, NULL) after an
       exchange that leaves NumPart unequal across ranks, so testing it before deciding to
       communicate would have some ranks reduce while their peers did not, offsetting every later
       collective. Any rank answering "not the whole set" vetoes, hence MAX of the negation. */
    flags_local[2] = (mp == NULL && npart == NumPart) ? 0 : 1;
    long long need_max = 0;
    MPI_Allreduce(flags_local, flags_any, 3, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&foreign_needed, &need_max, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    int overflow_any = flags_any[0], hardfail_any = flags_any[1], not_whole_tree_any = flags_any[2];

    if(hardfail_any)
    {
        /* Non-retryable. The failing rank records the cause; pseudodata failures
         * already soft-stopped inside force_exchange_pseudodata_complete(). Skip the
         * foreign-moment scatter/finalize/resum; drains at gravtree:after_treebuild. */
        if(let_status == LET_PACK_OOM)
        {
            printf("LET wire OOM on rank %d: ran out of memory reserving ghost-import (LET) transport buffers "
                   "(the memory ledger's LET-wire 'failed' bytes shows how much). This configuration's ghost import "
                   "does not fit in node memory; feasible: fewer ranks/node, or lower resolution. Stopping cleanly.\n",
                   ThisTask);
            fflush(stdout);
        }
        /* LET_ARENA_SHORT and LET_FOREIGN_STORAGE_SHORT already reported themselves, from the rank
         * that is actually short and with the sizes involved, at the point the shortfall was found.
         * Every rank carries those statuses because the decision is collective, so anything printed
         * here would print NTask times and say nothing the short rank has not said more precisely. */
        if(let_status == LET_PACK_OOM || let_status == LET_ARENA_SHORT
           || let_status == LET_FOREIGN_STORAGE_SHORT || let_status == LET_UNPACK_INTERNAL) {endrun(90000072);}
    }
    else if(overflow_any)
    {
        if(let_retry >= LET_MAX_RETRY)
        {
            /* Should not happen: one rebuild per ratchet suffices (fixed in-call domain
             * + 1.5x headroom). Each overflowing rank reports its own need, then stop.
             * Which rank is short is asked of the capacity itself: the exchange answers
             * whether ANY rank overflowed, so its status is the same on every rank and
             * cannot say who. */
            if(foreign_needed > (long long) MaxForeignNodes)
                printf("The foreign-node index range could not be made large enough in %d rebuilds on rank=%d: needed %lld nodes against a range of %d (local nodes %d). The range grows itself when an import does not fit, so reaching this means the demand outran it every time; running on more ranks or nodes reduces what each one imports. Stopping.\n",
                       LET_MAX_RETRY, ThisTask, foreign_needed, MaxForeignNodes, MaxNodes);
            fflush(stdout);
            endrun(90000062);   /* graceful drain; skip the foreign-moment steps */
        }
        else
        {
            /* Raise the index ceiling to the demand plus a pad, so one rebuild settles the retry
             * instead of ratcheting again on the next slightly larger step.  The pad is generous
             * on purpose: what it buys is four bytes per slot in Nextnode[], since the nodes
             * themselves are no longer sized from this ceiling but from the exact import
             * (force_tree_grow_foreign_storage).  A rebuild is far more expensive than the ints,
             * so the ceiling is the one place here worth being loose about.  The floor is monotone
             * and restart-persisted, which is now equally cheap for the same reason. */
            long long want = need_max + (long long) ceil(LET_FOREIGN_PAD_FRACTION * (double) need_max);
            if(want > RuntimeMinLETForeignNodes) {RuntimeMinLETForeignNodes = want;}
            if(ThisTask == 0)
                printf("LET foreign-node index range too small (needed up to %lld nodes); raising the ceiling to %lld and rebuilding the tree (retry %d/%d).\n",
                       need_max, RuntimeMinLETForeignNodes, let_retry + 1, LET_MAX_RETRY);
            fflush(stdout);
            let_retry++;
            int maxnodes_same = MaxNodes;   /* only the foreign floor grew; keep the local node sizing */
            force_treefree();
            force_treeallocate(maxnodes_same, All.MaxPartExpandable);
            /* drain a tree-alloc UVM OOM before rebuilding on a NULL-based tree (symmetric) */
            gizmo_exit_bad_stop_if_requested("gravtree:treeallocate");
            goto let_build_attempt;   /* rebuild the whole tree with the larger arena */
        }
    }
    else
    {
        /* Success: foreign moments complete. Scatter AoS→SoA, finalize LET
         * completeness, then re-sum ancestor topnode moments. */
        if(gpu_scatter_pseudo_to_soa() != 0)    {endrun(90000073);}
        let_finalize_unredirected_foreign_topleaves();
        if(gpu_topnode_moment_resum() != 0)     {endrun(90000074);}
    }
    /* Tree-integrity invariant. The root node's count is re-accumulated across foreign domains
       above, so after the resum it must equal the global particle number. A particle inserted
       under a top node owned by another task is detached and then appears in NO rank's tree:
       absent from every multipole moment, contributing to no force, and otherwise silent. The
       test is on the count rather than the mass because node masses are MyFloat, so one particle
       lost in a million falls below the accumulated summation error. Warn rather than abort, so a
       false positive cannot kill a long run. Group finding builds over subsets, whose root
       legitimately holds fewer, hence the veto reduced above. */
    if(!not_whole_tree_any)
    {
        long long in_tree = (long long) Nodes[All.TreeNodeIndexBase].N_part;
        if(in_tree != (long long) All.TotNumPart && ThisTask == 0)
        {
            printf("WARNING: tree integrity check failed: root node holds %lld particles, expected %lld. "
                   "Particles missing from the tree contribute to no rank's multipole moments, so the "
                   "forces on everything else are wrong by their mass.\n",
                   in_tree, (long long) All.TotNumPart);
            fflush(stdout);
        }
    }

    TimeOfLastTreeConstruction = All.Time;
    g_force_treebuild_generation++;   /* topology + Father[] + node structure changed */
    /* Only a whole-tree build leaves Father[] describing every local particle; a subset build writes
       links for its own members and leaves the rest at -1, which is no use as an attachment record.
       A whole-tree build whose retained stage was vetoed because some other rank was building a subset
       does not qualify either: its crossers went under foreign top-leaves, so its Father[] would send
       the next build's recovery to a leaf this rank does not own. */
    g_force_global_topology_valid = (whole_tree_local && counts_any[2] == 0) ? 1 : 0;
    return Numnodestree;
}



/*! Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.TreeNodeIndexBase.... All.TreeNodeIndexBase+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.TreeNodeIndexBase] gives the first tree node. Finally, node indices
 *  with values 'All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the hydro
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
/* What gpu_topology_emit_bfs's nonzero codes mean, so a failure says which one happened instead
 * of printing a bare number the reader has to go and look up. */
static const char *force_emit_bfs_reason(int rc)
{
    switch(rc)
    {
        case 1:  return "the node arena overflowed";
        case 3: return "a build allocation had already failed, or the leaves needing a member list did not fit the room reserved for them";
        case 4:  return "the breadth-first emit hit its depth guard, so some node would not subdivide";
        default: return "an unrecognised failure";
    }
}

int force_treebuild_single(int npart, struct unbind_data *mp)
{
    int i, j, k, subnode = 0, shift, parent, numnodes, rep, nfree, th, nn, no;
    struct NODE *nfreep;
    MyFloat lenhalf;
    peanokey key, morton, th_key, *morton_list;
    
    /* create an empty root node  */
    nfree = All.TreeNodeIndexBase;        /* index of first free node */
    nfreep = &Nodes[nfree];    /* select first node */
    nfreep->len = DomainLen;
    nfreep->center = {(MyFloat)DomainCenter[0], (MyFloat)DomainCenter[1], (MyFloat)DomainCenter[2]};
    for(j = 0; j < 8; j++) {nfreep->u.suns[j] = -1;}
    numnodes = 1;
    nfreep++;
    nfree++;
    
    /* create a set of empty nodes corresponding to the top-level domain grid. We need to generate these nodes first to make sure that we have a
     * complete top-level tree which allows the easy insertion of the pseudo-particles at the right place */
    
    /* Root topnode 0 maps to the root tree node All.TreeNodeIndexBase; children are mapped
     * inside the recursion (top-leaf router geometry SSOT, H0). */
    TopNodeNodeIndex[0] = All.TreeNodeIndexBase;
    if(TopNodes[0].Daughter < 0) {
        /* Degenerate root-is-leaf (single top-cell): the recursion below sets no
         * children, so set the leaf's DomainNodeIndex explicitly to keep both maps
         * complete + consistent. */
        DomainNodeIndex[TopNodes[0].Leaf] = All.TreeNodeIndexBase;
    }
    if(force_create_empty_nodes(All.TreeNodeIndexBase, 0, 1, 0, 0, 0, &numnodes, &nfree) < 0) {return FORCE_TREEBUILD_OUT_OF_NODES;}
    /* H0 post-build validation: every topnode must map to a valid Nodes[] slot. */
    {
        const int node_lo = All.TreeNodeIndexBase, node_hi = All.TreeNodeIndexBase + MaxNodes;
        for(int tnchk = 0; tnchk < NTopnodes; tnchk++) {
            if(TopNodeNodeIndex[tnchk] < node_lo || TopNodeNodeIndex[tnchk] >= node_hi) {
                printf("force_treebuild: TopNodeNodeIndex[%d]=%d out of range [%d,%d) (NTopnodes=%d) — top-leaf router map incomplete\n",
                       tnchk, TopNodeNodeIndex[tnchk], node_lo, node_hi, NTopnodes);
                endrun(91561);
            }
        }
    }
    /* if a high-resolution region in a global tree is used, we need to generate an additional set empty nodes to make sure that we have a complete top-level tree for the high-resolution inset */

    /* GPU tree-build replaces the per-particle CPU
     * insertion loop for inside-topleaf topology.  Order on GPU compile:
     *   1. force_insert_pseudo_particles (modifies foreign-topleaf u.suns).
     *   2. Acquire SoA + Peano-walk mirrors.
     *   3. gpu_topology_build_data_path: per-particle Peano walk + Morton
     *      sort within each topleaf.
     *   4. gpu_topology_emit_bfs: BFS from each topleaf root, emits
     *      inside-topleaf internal-node topology into SoA suns_backup,
     *      center, len, father.
     *   5. Writeback inside-topleaf range AoS u.suns/center/len.  This
     *      writeback exists (with force_update_node_recursive retired) so
     *      that any non-GPU CPU consumer of AoS u.suns sees complete topology
     *      between force_treebuild_single and the final
     *      gpu_topology_writeback_d_to_aos in force_treebuild.
     *
     * On overflow (rc=1), return -1 so force_treebuild's outer loop grows
     * TreeAllocFactor and retries -- same contract as the CPU path. */
    {
        force_insert_pseudo_particles();

        /* the old mark_all_dirty + acquire pair triggered seed_full_
         * to copy AoS topnode center/len into SoA before BFS.  That seeding now
         * happens inside gpu_nextnode_backup_suns below (single GPU kernel reads
         * UVM AoS, writes SoA suns_backup + center + len for [0..numnodes)). */
        if(gpu_peano_walk_acquire() != 0)
        {
            printf("force_treebuild_single: rank %d could not acquire the peano-walk mirrors. That is an allocation\n"
                   "failure, not a shortage of tree nodes, so the build is not retried with a larger node arena.\n", ThisTask);
            fflush(stdout);
            return FORCE_TREEBUILD_HARD_FAILURE;
        }
        /* Snapshot topnode u.suns -> SoA suns_backup.  At this point u.suns
         * for intermediate topnodes is populated by force_create_empty_nodes;
         * force_insert_pseudo_particles set u.suns[0] for foreign topleafs;
         * local topleaf u.suns are uninitialized (BFS will overwrite their
         * suns_backup entries with the local particle subtree topology). */
        gpu_nextnode_backup_suns(numnodes);
        /* It reports an allocation failure by invalidating the SoA and requesting a stop, and it returns
         * void.  Reading that here names the real site: otherwise the first symptom is the topology emit
         * finding a null SoA, which reads as a topology fault rather than as running out of memory. */
        if(!gpu_gravity_tree_valid())
        {
            printf("force_treebuild_single: rank %d could not allocate the tree SoA for %d nodes. That is an\n"
                   "allocation failure, not a shortage of tree nodes, so the build is not retried larger.\n", ThisTask, numnodes);
            fflush(stdout);
            return FORCE_TREEBUILD_HARD_FAILURE;
        }

        if(gpu_topology_build_data_path(npart, mp) != 0)
        {
            printf("force_treebuild_single: rank %d could not build the topology data path -- an allocation failure,\n"
                   "a key/top-tree disagreement, or a missing DomainTask mirror; the site above says which. None of\n"
                   "them is a shortage of tree nodes, so the build is not retried with a larger node arena.\n", ThisTask);
            fflush(stdout);
            return FORCE_TREEBUILD_HARD_FAILURE;
        }
        int new_numnodes = numnodes;
        int rc = gpu_topology_emit_bfs(numnodes, &new_numnodes);
        if(rc == 1) {return FORCE_TREEBUILD_OUT_OF_NODES;}   /* the one cause a larger node arena fixes */
        if(rc != 0) {
            printf("force_treebuild_single: the topology emit failed on rank %d (rc=%d): %s. A larger node arena\n"
                   "does not address that, so the build is not retried.\n", ThisTask, rc, force_emit_bfs_reason(rc));
            fflush(stdout);
            return FORCE_TREEBUILD_HARD_FAILURE;
        }
        int topnode_end = numnodes;
        numnodes = new_numnodes;

        /* Writeback GPU-built suns / center / len to AoS for the FULL range
         * (0..numnodes), not just the new BFS nodes (topnode_end..numnodes).
         * The topleaf nodes (0..topnode_end) have their soa->suns_backup entries
         * updated by BFS to hold the local-particle subtree root indices.
         * Without flushing these back to AoS, force_update_node_recursive reads
         * the stale (uninitialized) topleaf suns, never reaches local particles,
         * and Father[i] is never set -- causing an infinite loop in
         * setup_smoothinglengths which walks Nodes[Father[i]].u.d.father. */
        if(gpu_topology_writeback_to_aos(0, numnodes) != 0)
        {
            printf("force_treebuild_single: rank %d could not write the built topology back to the node array.\n"
                   "That is not a shortage of tree nodes, so the build is not retried larger.\n", ThisTask);
            fflush(stdout);
            return FORCE_TREEBUILD_HARD_FAILURE;
        }
    }

    /* now compute the multipole moments recursively */
    last = -1;
    /* force_update_node_recursive retired on GPU build.  The GPU
     * finalize stage in force_treebuild (gpu_topology_finalize_father,
     * gpu_topology_finalize_sibling, gpu_moment_refresh, gpu_nextnode_thread)
     * now produces all of FUNR's outputs (sibling, father, Father[], moments,
     * nextnode).  The Nextnode[last]=-1 tail fixup is redundant: sibling
     * for the root is -1, which propagates through the DFS chain in
     * gpu_nextnode_thread to give the last DFS particle Nextnode[] = -1
     * automatically.  The second gpu_nextnode_backup_suns is removed too:
     * snapshot #1 (in the GPU build block above) plus emit_bfs's direct
     * SoA writes give complete suns_backup coverage, and nothing clobbers
     * AoS u.suns until gpu_topology_writeback_d_to_aos at the very end of
     * force_treebuild (by which point all SoA readers are done). */

    return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 */
int force_create_empty_nodes(int no, int topnode, int bits, peano1D x, peano1D y, peano1D z, int *nodecount,
                              int *nextfree)
{
    int i, j, k, n, sub, count;
    MyFloat lenhalf;

    if(TopNodes[topnode].Daughter >= 0)
    {
        for(i = 0; i < 2; i++)
            for(j = 0; j < 2; j++)
                for(k = 0; k < 2; k++)
                {
                    sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

                    count = i + 2 * j + 4 * k;

                    Nodes[no].u.suns[count] = *nextfree;

                    /* H0 (top-leaf router geometry SSOT): map this child topnode
                     * (PH offset Daughter+sub) to the physical Nodes[] slot
                     * (*nextfree) whose exact center/len is set just below.  NOTE
                     * the topnode child offset `sub` is Peano-Hilbert while the
                     * Nodes suns index `count` is Morton i+2j+4k — they differ, so
                     * the router MUST read geometry via this map, never derive it. */
                    TopNodeNodeIndex[TopNodes[topnode].Daughter + sub] = *nextfree;

                    lenhalf = 0.25 * Nodes[no].len;
                    Nodes[*nextfree].len = 0.5 * Nodes[no].len;
                    Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
                    Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
                    Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;

                    for(n = 0; n < 8; n++)
                        Nodes[*nextfree].u.suns[n] = -1;

                    if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
                    {
                        DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;
                        /* H0 SSOT consistency: a leaf topnode's router slot must equal
                         * its DomainNodeIndex entry (both are *nextfree here). Guards
                         * against future desync of the two maps. */
                        if(TopNodeNodeIndex[TopNodes[topnode].Daughter + sub] !=
                           DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf])
                        {
                            printf("force_create_empty_nodes: TopNodeNodeIndex/DomainNodeIndex mismatch "
                                   "(child topnode %d, leaf %d)\n",
                                   TopNodes[topnode].Daughter + sub,
                                   TopNodes[TopNodes[topnode].Daughter + sub].Leaf);
                            endrun(91560);
                        }
                    }

                    *nextfree = *nextfree + 1;
                    *nodecount = *nodecount + 1;

                    if((*nodecount) >= MaxNodes)
                    {
                        printf("task %d: maximum number MaxNodes=%d of tree-nodes reached."
                               "MaxTopNodes=%d NTopnodes=%d NTopleaves=%d nodecount=%d\n",
                               ThisTask, MaxNodes, MaxTopNodes, NTopnodes, NTopleaves, *nodecount);
                        printf("in create empty nodes\n");
                        if(All.TreeAllocFactor > 5.0)
                        {
                            dump_particles();
                            endrun(11);
                        }
                        return -1; /* signal to caller to retry with larger TreeAllocFactor */
                    }

                    if(force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
                                             bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree) < 0)
                        return -1;
                }
    }
    return 0;
}



/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 */
void force_insert_pseudo_particles(void)
{
    int i, index;
    
    for(i = 0; i < NTopleaves; i++)
    {
        index = DomainNodeIndex[i];
        
        if(DomainTask[i] != ThisTask)
            Nodes[index].u.suns[0] = All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes + i;    /* pseudo-particles live above the foreign-node range */
    }
}






/*! Pseudo-particle exchange wire format.  Lifted to file scope so
 *  the issue/complete halves of force_exchange_pseudodata can share the type
 *  across the LET overlap window. */
struct DomainNODE
    {
        Vec3<MyFloat> s;
        Vec3<MyFloat> vs;
        MyFloat mass;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        MyFloat gasmass;
#endif
        MyFloat hmax;
        MyFloat hmax_per_type[6];   /* cross-rank per-type h band (mirror of Extnodes.hmax_per_type) */
        /* The owner's sidelength for this top-leaf.  Every rank builds the top tree to the same
           nominal geometry, but the owner may have grown this leaf to cover a particle it kept there
           after the particle drifted out of the leaf's nominal cube.  That growth is not a function of
           anything else on the wire -- unlike the drift growth, which every rank reproduces from vmax
           by the same rule -- so without it a remote walk can accept the nominal box and be wrong. */
        MyFloat len;
        MyFloat vmax;
        MyFloat divVmax;
        long N_part;
        MyFloat maxsoft;
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        MyFloat cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
        MyFloat stellar_lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
        double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS];
        double chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Vec3<MyFloat> rt_source_lum_s;
        Vec3<MyFloat> rt_source_lum_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
        MyFloat sink_lum; Vec3<MyFloat> sink_lum_grad;
#endif
#ifdef SINK_CALC_DISTANCES
        MyFloat sink_mass;
        Vec3<MyFloat> sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
        int N_SINK;
        Vec3<MyFloat> sink_vel;
#ifdef SPECIAL_POINT_MOTION
        Vec3<MyFloat> sink_acc;
#endif
#ifdef  SINGLE_STAR_FB_TIMESTEPLIMIT
        MyFloat MaxFeedbackVel;
#endif
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        SymmetricTensor2<MyFloat> tidal_tensorps_prevstep;
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Vec3<MyFloat> s_dm;
        Vec3<MyFloat> vs_dm;
        MyFloat mass_dm;
#endif
        unsigned int bitflags;
#ifdef PAD_STRUCTURES
        int pad[3];
#endif
    };

/*! State shared between force_exchange_pseudodata_issue() and
 *  ..._complete() for the non-blocking-overlap pattern.  let_run_exchange
 *  runs concurrently with the pseudodata Iallgathervs in the GPU build path. */
static struct DomainNODE *DomainMoment_pending = NULL;
static MPI_Request *pseudo_requests_pending = NULL;
static int *pseudo_recvcounts_pending = NULL;
static int *pseudo_recvoffset_pending = NULL;
static int  pseudo_n_requests_pending = 0;

void force_exchange_pseudodata_issue(void)
{
    int i, no, m;
    /* Re-entrant issue() (a prior issue had no matching complete()): soft
     * bad-stop + return BEFORE allocating/posting a second Iallgatherv, so the
     * already-pending exchange is left intact for its complete(). No poll here
     * (a nonblocking collective may be outstanding); drains at a later poll. */
    if(DomainMoment_pending != NULL) {endrun(90000075); return;}

    DomainMoment_pending = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));
    struct DomainNODE *DomainMoment = DomainMoment_pending;

    for(m = 0; m < All.DomainSegmentsPerRank; m++)
        for(i = DomainStartList[ThisTask * All.DomainSegmentsPerRank + m];
            i <= DomainEndList[ThisTask * All.DomainSegmentsPerRank + m]; i++)
        {
            no = DomainNodeIndex[i];
            
            /* read out the multipole moments from the local base cells */
            DomainMoment[i].s = Nodes[no].u.d.s;
            DomainMoment[i].vs = Extnodes[no].vs;
            DomainMoment[i].mass = Nodes[no].u.d.mass;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
            DomainMoment[i].gasmass = Nodes[no].gasmass;
#endif
            DomainMoment[i].hmax = Extnodes[no].hmax;
            for(int t = 0; t < 6; t++) DomainMoment[i].hmax_per_type[t] = Extnodes[no].hmax_per_type[t];
            DomainMoment[i].len = Nodes[no].len;
            DomainMoment[i].vmax = Extnodes[no].vmax;
            DomainMoment[i].divVmax = Extnodes[no].divVmax;
            DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
            DomainMoment[i].N_part = Nodes[no].N_part;
            DomainMoment[i].maxsoft = Nodes[no].maxsoft;
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            DomainMoment[i].cr_injection = Nodes[no].cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
            int k; for(k=0;k<N_RT_FREQ_BINS;k++) {DomainMoment[i].stellar_lum[k] = Nodes[no].stellar_lum[k];}
#ifdef CHIMES_STELLAR_FLUXES
            for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++)
            {
                DomainMoment[i].chimes_stellar_lum_G0[k] = Nodes[no].chimes_stellar_lum_G0[k];
                DomainMoment[i].chimes_stellar_lum_ion[k] = Nodes[no].chimes_stellar_lum_ion[k];
            }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            DomainMoment[i].rt_source_lum_s = Nodes[no].rt_source_lum_s;
            DomainMoment[i].rt_source_lum_vs = Extnodes[no].rt_source_lum_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
            DomainMoment[i].sink_lum = Nodes[no].sink_lum;
            DomainMoment[i].sink_lum_grad = Nodes[no].sink_lum_grad;
#endif
#ifdef SINK_CALC_DISTANCES
            DomainMoment[i].sink_mass = Nodes[no].sink_mass;
            DomainMoment[i].sink_pos = Nodes[no].sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
            DomainMoment[i].sink_vel = Nodes[no].sink_vel;
            DomainMoment[i].N_SINK = Nodes[no].N_SINK;
#ifdef SPECIAL_POINT_MOTION
            DomainMoment[i].sink_acc = Nodes[no].sink_acc;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
            DomainMoment[i].MaxFeedbackVel = Nodes[no].MaxFeedbackVel;
#endif
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
            DomainMoment[i].tidal_tensorps_prevstep = Nodes[no].tidal_tensorps_prevstep;
#endif
#ifdef DM_SCALARFIELD_SCREENING
            DomainMoment[i].s_dm = Nodes[no].s_dm;
            DomainMoment[i].mass_dm = Nodes[no].mass_dm;
            DomainMoment[i].vs_dm = Extnodes[no].vs_dm;
#endif
        }

    /* Post one MPI_Iallgatherv per All.DomainSegmentsPerRank slice; the requests
     * are stored in static pseudo_requests_pending and waited on in _complete().
     * Per-slice recvcounts/recvoffset arrays must remain valid until Wait, so
     * we allocate one set per slice and free them all in _complete(). */
    pseudo_n_requests_pending = All.DomainSegmentsPerRank;
    pseudo_requests_pending = (MPI_Request *) mymalloc("pseudo_requests",
                                  All.DomainSegmentsPerRank * sizeof(MPI_Request));
    pseudo_recvcounts_pending = (int *) mymalloc("pseudo_recvcounts",
                                  All.DomainSegmentsPerRank * NTask * sizeof(int));
    pseudo_recvoffset_pending = (int *) mymalloc("pseudo_recvoffset",
                                  All.DomainSegmentsPerRank * NTask * sizeof(int));
    for(m = 0; m < All.DomainSegmentsPerRank; m++)
    {
        int *rc = pseudo_recvcounts_pending + m * NTask;
        int *ro = pseudo_recvoffset_pending + m * NTask;
        for(int recvTask = 0; recvTask < NTask; recvTask++)
        {
            rc[recvTask] =
                (DomainEndList[recvTask * All.DomainSegmentsPerRank + m] -
                 DomainStartList[recvTask * All.DomainSegmentsPerRank + m] + 1)
                * sizeof(struct DomainNODE);
            /* MPI_Iallgatherv takes int byte counts and displacements, so the whole pseudodata
             * block has to stay under 2 GB.  That ceiling is a property of this exchange, not of
             * the caller, and silently wrapping it would hand MPI a negative displacement -- so
             * check it here, where the number is formed. */
            const long long offset_bytes =
                (long long) DomainStartList[recvTask * All.DomainSegmentsPerRank + m]
                * (long long) sizeof(struct DomainNODE);
            if(offset_bytes > (long long) INT_MAX)
              {
                if(ThisTask == 0)
                  {
                    printf("Pseudo-particle exchange needs a %lld byte offset, beyond what MPI's int displacements can carry.\n", offset_bytes);
                    printf("There are %d top-tree leaves; lower DOMAIN_SEGMENTS_SCALE or run on fewer ranks.\n", NTopleaves);
                    fflush(stdout);
                  }
                endrun(90000025);
                return;
              }
            ro[recvTask] = (int) offset_bytes;
        }
        MPI_Iallgatherv(MPI_IN_PLACE, rc[ThisTask], MPI_BYTE,
                        &DomainMoment[0], rc, ro, MPI_BYTE, MPI_COMM_WORLD,
                        &pseudo_requests_pending[m]);
    }
}

/*! Waits on the Iallgathervs posted by force_exchange_pseudodata_issue() and
 *  unpacks the received topleaf moments into the AoS Nodes_base / Extnodes_base.
 *  GPU build path calls let_run_exchange() in between to overlap MPI; CPU and
 *  refresh paths call the sync wrapper force_exchange_pseudodata() below. */
int force_exchange_pseudodata_complete(void)
{
    /* Unmatched complete (pending==NULL = complete without a matching issue, or a
     * double-complete): symmetric control-flow invariant. Soft bad-stop + status-return
     * (1) so the caller skips the foreign-moment scatter/finalize/resum -- which would run
     * on un-exchanged moments and could itself fatal hard-exit in let_finalize -- and
     * drains at the next poll. */
    if(DomainMoment_pending == NULL) {endrun(90000076); return 1;}
    struct DomainNODE *DomainMoment = DomainMoment_pending;

    const int n_requests_issued = pseudo_n_requests_pending;
    MPI_Waitall(pseudo_n_requests_pending, pseudo_requests_pending, MPI_STATUSES_IGNORE);

    /* Free request/count buffers (LIFO order: ro, rc, requests). */
    myfree(pseudo_recvoffset_pending);
    myfree(pseudo_recvcounts_pending);
    myfree(pseudo_requests_pending);
    pseudo_requests_pending = NULL;
    pseudo_recvcounts_pending = NULL;
    pseudo_recvoffset_pending = NULL;
    pseudo_n_requests_pending = 0;

    int i, no, m, ta;
    /* The segment count this exchange was posted with, not whatever it is now: the domain lists
     * being walked here are the ones that were current at issue time. */
    const int segments = n_requests_issued;
    for(ta = 0; ta < NTask; ta++)
        if(ta != ThisTask)
            for(m = 0; m < segments; m++)
                for(i = DomainStartList[ta * segments + m]; i <= DomainEndList[ta * segments + m]; i++)
                {
                    no = DomainNodeIndex[i];

                    Nodes[no].u.d.s = DomainMoment[i].s;
                    Extnodes[no].vs = DomainMoment[i].vs;
                    Nodes[no].u.d.mass = DomainMoment[i].mass;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                    Nodes[no].gasmass = DomainMoment[i].gasmass;
#endif
                    Extnodes[no].hmax = DomainMoment[i].hmax;
                    for(int t = 0; t < 6; t++) Extnodes[no].hmax_per_type[t] = DomainMoment[i].hmax_per_type[t];
                    /* The receiver's own copy is the nominal cube, so the maximum is the owner's value;
                       taken as a maximum rather than an assignment so the node can only ever bound more
                       of its contents, never less. */
                    if(DomainMoment[i].len > Nodes[no].len) {Nodes[no].len = DomainMoment[i].len;}
                    Extnodes[no].vmax = DomainMoment[i].vmax;
                    Extnodes[no].divVmax = DomainMoment[i].divVmax;
                    Nodes[no].N_part = DomainMoment[i].N_part;
                    Nodes[no].u.d.bitflags = (Nodes[no].u.d.bitflags & (~((1 << BITFLAG_MULTIPLEPARTICLES)))) | (DomainMoment[i].bitflags & ((1 << BITFLAG_MULTIPLEPARTICLES)));
                    Nodes[no].maxsoft = DomainMoment[i].maxsoft;
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                    Nodes[no].cr_injection = DomainMoment[i].cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
                    int k; for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k] = DomainMoment[i].stellar_lum[k];}
#ifdef CHIMES_STELLAR_FLUXES
                    for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++)
                    {
                        Nodes[no].chimes_stellar_lum_G0[k] = DomainMoment[i].chimes_stellar_lum_G0[k];
                        Nodes[no].chimes_stellar_lum_ion[k] = DomainMoment[i].chimes_stellar_lum_ion[k];
                    }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
                    Nodes[no].rt_source_lum_s = DomainMoment[i].rt_source_lum_s;
                    Extnodes[no].rt_source_lum_vs = DomainMoment[i].rt_source_lum_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
                    Nodes[no].sink_lum = DomainMoment[i].sink_lum;
                    Nodes[no].sink_lum_grad = DomainMoment[i].sink_lum_grad;
#endif
#ifdef SINK_CALC_DISTANCES
                    Nodes[no].sink_mass = DomainMoment[i].sink_mass;
                    Nodes[no].sink_pos = DomainMoment[i].sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
                    Nodes[no].sink_vel = DomainMoment[i].sink_vel;
                    Nodes[no].N_SINK = DomainMoment[i].N_SINK;
#ifdef SPECIAL_POINT_MOTION
                    Nodes[no].sink_acc = DomainMoment[i].sink_acc;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                    Nodes[no].MaxFeedbackVel = DomainMoment[i].MaxFeedbackVel;
#endif
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                    Nodes[no].tidal_tensorps_prevstep = DomainMoment[i].tidal_tensorps_prevstep;
#endif
#ifdef DM_SCALARFIELD_SCREENING
                    Nodes[no].s_dm = DomainMoment[i].s_dm;
                    Nodes[no].mass_dm = DomainMoment[i].mass_dm;
                    Extnodes[no].vs_dm = DomainMoment[i].vs_dm;
#endif
                }

    myfree(DomainMoment);
    DomainMoment_pending = NULL;
    return 0;
}

/*! Synchronous wrapper (begin + complete in one call) for the CPU and
 *  refresh code paths (which do not have a LET round to overlap with).
 *  Returns the complete() status (nonzero = unmatched; caller skips dependent
 *  pseudo-update work and drains at its poll). */
int force_exchange_pseudodata(void)
{
    force_exchange_pseudodata_issue();
    return force_exchange_pseudodata_complete();
}



/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_pseudos(int no)
{
    int j, p, count_particles, multiple_flag;
    MyFloat hmax, vmax;
    MyFloat hmax_per_type[6];
    MyFloat divVmax;
    Vec3<MyFloat> s, vs; MyFloat mass;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    MyFloat gasmass = 0;
#endif

#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double cr_injection = 0;
#endif
#ifdef RT_USE_GRAVTREE
    MyFloat stellar_lum[N_RT_FREQ_BINS]={0};
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS]={0};
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<MyFloat> rt_source_lum_s, rt_source_lum_vs;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Vec3<MyFloat> s_dm, vs_dm; MyFloat mass_dm;
#endif
    
    MyFloat maxsoft;
    
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    rt_source_lum_s = rt_source_lum_vs = {};
#endif
#ifdef SINK_PHOTONMOMENTUM
    MyFloat sink_lum = 0; Vec3<MyFloat> sink_lum_grad = {};
#endif
#ifdef SINK_CALC_DISTANCES
    MyFloat sink_mass=0;
    Vec3<MyFloat> sink_pos_times_mass = {};
#if defined(SINK_NODE_MOTION_TRACKED)
    Vec3<MyFloat> sink_mom = {};
    int N_SINK = 0;
#ifdef SPECIAL_POINT_MOTION
    Vec3<MyFloat> sink_force = {};
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    MyFloat max_feedback_vel=0;
#endif
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    SymmetricTensor2<MyFloat> tidal_tensorps_prevstep = {};
#endif
#ifdef DM_SCALARFIELD_SCREENING
    mass_dm = 0;
    s_dm = vs_dm = {};
#endif
    mass = 0;
    s = vs = {};
    hmax = 0;
    for(int t = 0; t < 6; t++) hmax_per_type[t] = 0;
    vmax = 0;
    divVmax = 0;
    count_particles = 0;
    maxsoft = 0;

    p = Nodes[no].u.d.nextnode;
    
    for(j = 0; j < 8; j++)    /* since we are dealing with top-level nodes, we now that there are 8 consecutive daughter nodes */
    {
        if(p >= All.TreeNodeIndexBase && p < All.TreeNodeIndexBase + MaxNodes)    /* internal node */
        {
            if(Nodes[p].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL)) {force_treeupdate_pseudos(p);}
            
            mass += (Nodes[p].u.d.mass);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
            gasmass += Nodes[p].gasmass;
#endif
            s += Nodes[p].u.d.mass * Nodes[p].u.d.s;
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            cr_injection += Nodes[p].cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
            int k; for(k=0;k<N_RT_FREQ_BINS;k++) {stellar_lum[k] += (Nodes[p].stellar_lum[k]);}
#ifdef CHIMES_STELLAR_FLUXES
            for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++)
            {
                chimes_stellar_lum_G0[k] += Nodes[p].chimes_stellar_lum_G0[k];
                chimes_stellar_lum_ion[k] += Nodes[p].chimes_stellar_lum_ion[k];
            }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            double l_tot=0; for(k=0;k<N_RT_FREQ_BINS;k++) {l_tot += (Nodes[p].stellar_lum[k]);}
            rt_source_lum_s += l_tot * Nodes[p].rt_source_lum_s;
            rt_source_lum_vs += l_tot * Extnodes[p].rt_source_lum_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
            sink_lum += Nodes[p].sink_lum;
            sink_lum_grad += Nodes[p].sink_lum * Nodes[p].sink_lum_grad;
#endif
#ifdef SINK_CALC_DISTANCES
            sink_mass += Nodes[p].sink_mass;
            sink_pos_times_mass += Nodes[p].sink_mass * Nodes[p].sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
            N_SINK += Nodes[p].N_SINK;
            sink_mom += Nodes[p].sink_mass * Nodes[p].sink_vel;
#ifdef SPECIAL_POINT_MOTION
            sink_force += Nodes[p].sink_mass * Nodes[p].sink_acc;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
            if(Nodes[p].sink_mass > 0) {max_feedback_vel = DMAX(max_feedback_vel, Nodes[p].MaxFeedbackVel);}
#endif
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
            {int k; for(k=0;k<6;k++) {tidal_tensorps_prevstep.data[k] += Nodes[p].u.d.mass * Nodes[p].tidal_tensorps_prevstep.data[k];}}
#endif
#ifdef DM_SCALARFIELD_SCREENING
            mass_dm += (Nodes[p].mass_dm);
            s_dm += Nodes[p].mass_dm * Nodes[p].s_dm;
            vs_dm += Nodes[p].mass_dm * Extnodes[p].vs_dm;
#endif
            vs += Nodes[p].u.d.mass * Extnodes[p].vs;
            
            if(Extnodes[p].hmax > hmax) {hmax = Extnodes[p].hmax;}
            for(int t = 0; t < 6; t++) {
                if(Extnodes[p].hmax_per_type[t] > hmax_per_type[t]) {hmax_per_type[t] = Extnodes[p].hmax_per_type[t];}
            }
            if(Extnodes[p].vmax > vmax) {vmax = Extnodes[p].vmax;}
            if(Extnodes[p].divVmax > divVmax) {divVmax = Extnodes[p].divVmax;}
            /* Count every particle, including those whose mass is zero (swallowed sinks awaiting
               cleanup). The device moment kernel counts them unconditionally, so gating here made
               the same tree report two different root counts depending on which path built it, and
               with it two different BITFLAG_MULTIPLEPARTICLES: a zero-mass node holding two
               particles read as single-particle, which the LET packer resolves by shipping one
               particle's leaf identity for the pair. N_part also divides node mass for the
               dynamical-friction effective mass, which diverged the same way. */
            count_particles += Nodes[p].N_part;
            if(Nodes[p].maxsoft > maxsoft) {maxsoft = Nodes[p].maxsoft;}
        }
        else
        {   /* invalid node type: soft bad-stop + break before the corrupt
             * Nodes[p].sibling deref below; drains at a gravtree poll */
            endrun(90000077);
            break;
        }

        p = Nodes[p].u.d.sibling;
    }
    
    if(mass)
    {
        s /= mass;
        vs /= mass;
    }
    else
    {
        s = Nodes[no].center;
        vs = {};
    }
    
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    double l_tot=0; int kfreq; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {l_tot += stellar_lum[kfreq];}
    if(l_tot)
    {
        rt_source_lum_s /= l_tot;
        rt_source_lum_vs /= l_tot;
    }
    else
    {
        rt_source_lum_s = Nodes[no].center;
        rt_source_lum_vs = {};
    }
#endif
#ifdef SINK_PHOTONMOMENTUM
    if(sink_lum)
    {
        sink_lum_grad /= sink_lum;
    }
    else
    {
        sink_lum_grad = {0, 0, 1};
    }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    if(mass_dm)
    {
        s_dm /= mass_dm;
        vs_dm /= mass_dm;
    }
    else
    {
        s_dm = Nodes[no].center;
        vs_dm = {};
    }
#endif


    Nodes[no].u.d.s = s;
    Extnodes[no].vs = vs;
    Nodes[no].u.d.mass = mass;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    Nodes[no].gasmass = gasmass;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    Nodes[no].cr_injection = cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
    int k; for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k] = stellar_lum[k];}
#ifdef CHIMES_STELLAR_FLUXES
    for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++)
    {
        Nodes[no].chimes_stellar_lum_G0[k] = chimes_stellar_lum_G0[k];
        Nodes[no].chimes_stellar_lum_ion[k] = chimes_stellar_lum_ion[k];
    }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Nodes[no].rt_source_lum_s = rt_source_lum_s;
    Extnodes[no].rt_source_lum_vs = rt_source_lum_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
    Nodes[no].sink_lum = sink_lum;
    Nodes[no].sink_lum_grad = sink_lum_grad;
#endif
#ifdef SINK_CALC_DISTANCES
    Nodes[no].sink_mass = sink_mass;
#if defined(SINK_NODE_MOTION_TRACKED)
    Nodes[no].N_SINK = N_SINK;
#endif
    if(sink_mass > 0)
    {
        Nodes[no].sink_pos = sink_pos_times_mass / sink_mass;
#if defined(SINK_NODE_MOTION_TRACKED)
        Nodes[no].sink_vel = sink_mom / sink_mass;
#if defined(SPECIAL_POINT_MOTION)
        Nodes[no].sink_acc = sink_force / sink_mass;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        Nodes[no].MaxFeedbackVel = max_feedback_vel;
#endif
#endif
    }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    {MyFloat inv_mass = 1.0/(mass+MIN_REAL_NUMBER); int k; for(k=0;k<6;k++) {Nodes[no].tidal_tensorps_prevstep.data[k] = tidal_tensorps_prevstep.data[k] * inv_mass;}}
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Nodes[no].s_dm = s_dm;
    Nodes[no].mass_dm = mass_dm;
    Extnodes[no].vs_dm = vs_dm;
#endif
    
    Extnodes[no].hmax = hmax;
    for(int t = 0; t < 6; t++) Extnodes[no].hmax_per_type[t] = hmax_per_type[t];
    Extnodes[no].vmax = vmax;
    Extnodes[no].divVmax = divVmax;
    Extnodes[no].Flag = GlobFlag;
    Nodes[no].N_part = count_particles; // record
    if(count_particles > 1) {multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);} else {multiple_flag = 0;}
    Nodes[no].u.d.bitflags &= (~((1 << BITFLAG_MULTIPLEPARTICLES)));    /* this clears the bits */
    Nodes[no].u.d.bitflags |= multiple_flag;
    Nodes[no].maxsoft = maxsoft;
    /* No GPU SoA hook here — force_treeupdate_pseudos is the CPU
     * pseudo-path; on GPU builds it is replaced by gpu_topnode_moment_resum
     * (gpu_pseudo_update.cc) which writes the SoA directly.  This function
     * compiles on both, but the GPU build never calls it. */
}



/*! This function flags nodes in the top-level tree that are dependent on
 *  local particle data.
 */
void force_flag_localnodes(void)
{
    int no, i, m;
    
    /* mark all top-level nodes */
    
    for(i = 0; i < NTopleaves; i++)
    {
        no = DomainNodeIndex[i];
        
        while(no >= 0)
        {
            if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL)) {break;}
            
            Nodes[no].u.d.bitflags |= (1 << BITFLAG_TOPLEVEL);
            
            no = Nodes[no].u.d.father;
        }
        
        /* mark also internal top level nodes */
        
        no = DomainNodeIndex[i];
        no = Nodes[no].u.d.father;
        
        while(no >= 0)
        {
            if(Nodes[no].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL)) {break;}
            
            Nodes[no].u.d.bitflags |= (1 << BITFLAG_INTERNAL_TOPLEVEL);
            
            no = Nodes[no].u.d.father;
        }
    }
    
    /* mark top-level nodes that contain local particles */
    
    for(m = 0; m < All.DomainSegmentsPerRank; m++)
        for(i = DomainStartList[ThisTask * All.DomainSegmentsPerRank + m];
            i <= DomainEndList[ThisTask * All.DomainSegmentsPerRank + m]; i++)
        {
            no = DomainNodeIndex[i];
            
            if(DomainTask[i] != ThisTask) {endrun(90000078); continue;} /* soft bad-stop + skip the foreign entry instead of mismarking its DEPENDS bitflags; drains at a gravtree poll */
            
            while(no >= 0)
            {
                if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_ELEMENT)) {break;}
                
                Nodes[no].u.d.bitflags |= (1 << BITFLAG_DEPENDS_ON_LOCAL_ELEMENT);
                
                no = Nodes[no].u.d.father;
            }
        }
}


/*! When a new additional resolution element is created, we can put it into the
 *  tree at the position of the spawning element. This is possible
 *  because the Nextnode[] array essentially describes the full tree walk as a
 *  link list. Multipole moments of tree nodes need not be changed (the new
 *  particle inherits the parent's position so mass+CoM are preserved at the
 *  insertion site; the 9.6 ForceAddElementToTree_CallsSinceBuild guardrail
 *  bounds drift to ancestor nodes between full rebuilds).
 *
 *  Father/Nextnode/Extnodes are
 *  UVM (SharedSpace), so the CPU mutations below
 *  are GPU-visible without a copy.  We additionally mirror the parent node's
 *  hmax/vmax/len into the SoA walk view, because the next gpu_force_drift_nodes
 *  early-outs when the parent's Ti_current already matches All.Ti_Current
 *  (insertions between drifts at the same Ti_current would otherwise leave
 *  the SoA stale until the next full rebuild).
 */
void force_add_element_to_tree(int iparent, int ichild)
{
    /* Both indices are written into the tree's particle-side arrays below, so both must be inside
     * them.  Creation sites already decline to make a particle the standing tree cannot carry, so
     * reaching this is a caller that skipped that test, not a tight fit -- stop rather than write
     * past Father[]/Nextnode[]. */
    if(!force_tree_is_allocated() || ichild >= All.TreeParticleSlots || iparent >= All.TreeParticleSlots)
    {
        printf("force_add_element_to_tree: task=%d asked to insert particle %d under %d, but the tree %s (particle slots=%d).\n",
               ThisTask, ichild, iparent, force_tree_is_allocated() ? "does not reach that far" : "is not standing", All.TreeParticleSlots);
        fflush(stdout);
        gizmo_request_controlled_stop(90000101, "force_add_element_to_tree: particle index outside the live tree's particle slots", __FILE__, __LINE__, __FUNCTION__);
        return;
    }
#ifndef MAINTAIN_TREE_IN_REARRANGE
    /* A rebuild is already required, so this tree will not be walked again.  Keep only what that
     * rebuild needs -- the new element belongs with its parent -- and leave the walk links alone: a
     * re-sequencing may have moved Father[] with the particles without moving them, and combining
     * the two would link this element into one node's traversal while growing another node's
     * bounds, which a neighbour search pruning on the first bound would then miss.
     *
     * This is a guard for a configuration rather than a repair of an observed fault.  Of the two
     * callers, the wind spawn inserts before its own routine asks for the rebuild, and star
     * formation reaches here only when GALSF_GENERATIONS exceeds one, since a single generation
     * converts the gas element in place instead of creating one.  So the case below arises for a
     * multi-generation star formation model running alongside spawned winds, and costs nothing to
     * carry until then.  The element is left for the rebuild to place, so when there is no record
     * to inherit the parent is set to the value a slot outside the tree carries, which every reader
     * of the whole array already skips. */
    if(TreeReconstructFlag)
    {
        Father[ichild] = g_force_global_topology_valid ? Father[iparent] : -1;
        return;
    }
#endif
    int father = Father[iparent];
    int no = Nextnode[iparent];
    Nextnode[iparent] = ichild; // insert new particle into linked list
    Nextnode[ichild] = no; // order correctly
    Father[ichild] = father; // set parent node to be the same
    // update parent node properties [maximum softening, speed] for opening criteria
    MyFloat new_hmax = DMAX(Extnodes[father].hmax, (MyFloat) moment_gas_hmax_from_kernelradius(P[iparent].KernelRadius, All.MaxKernelRadius));
    Extnodes[father].hmax = new_hmax;
    /* Mode B per-type incremental update — only the band corresponding to
     * iparent's type. Other bands unchanged (correct: this insertion adds
     * one particle, only its type's band can grow). Conservative across every
     * leaf-policy-selectable source via the shared helper. */
    {
        int ptype = (int)P[iparent].Type;
        double htmp = force_hmax_per_type_particle_radius(iparent);
        if(htmp > Extnodes[father].hmax_per_type[ptype]) {
            Extnodes[father].hmax_per_type[ptype] = (MyFloat)htmp;
        }
    }
    double new_vmax = DMAX((double) Extnodes[father].vmax, particle_motion_speed_bound(ichild, P, CellP));
    Extnodes[father].vmax = (MyFloat) new_vmax;

    /* Keep SoA walk-mirror coherent with the AoS Extnodes
     * change above.  hmax and vmax are read by the walk's opening criteria
     * (and vmax drives bbox expansion in subsequent drifts via Nodes[].len).
     * Indexing matches the local-or-foreign convention from gpu_force_drift_nodes:
     * SoA slot k = father - All.TreeNodeIndexBase for both local nodes (k < MaxNodes) and
     * foreign nodes (MaxNodes <= k < MaxNodes + AllocatedForeignNodes) -- bounded by the mirror
     * that exists, not by the index range it sits in. */
    {
        struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
        int k_soa = father - All.TreeNodeIndexBase;
        if(soa && k_soa >= 0 && k_soa < MaxNodes + AllocatedForeignNodes) {
            if(soa->hmax) {soa->hmax[k_soa] = (MyGravFloat) new_hmax;}
            if(soa->vmax) {soa->vmax[k_soa] = (MyGravFloat) new_vmax;}
            /* ⛔ node_ti is deliberately NOT written here. It must pair with the
               LENGTH it describes, and this site does not write soa->len: stamping a
               newer time onto an older length makes the walk widen by too little --
               silent under-inclusion. A fresher vmax is safe (running max, so it can
               only over-widen); a fresher time is not. */
        }
    }

    /* Each insertion stales the LET / pseudo-particle
     * moments shipped on the last full build.  Mass+CoM remain conserved at
     * the insertion site, but ancestor topnodes (and any rank's foreign view
     * of them) carry the pre-insertion moments until the next rebuild. */
    ForceAddElementToTree_CallsSinceBuild++;
}



/*! This routine computes the gravitational force for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
/*! The modern version of this routine handles both the PM-grid and non-PM
 *  cases, unlike the previous version (which used two, redundant, algorithms)
 */
/*! In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= PM_RCUT*PM_ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access panelty (which reduces cache performance) incurred by the
 *  table.
 */
/* Import-completeness record for the gravity walks.
 *
 * A walk that opens a foreign node the sender shipped multipole-only has outrun its import: the
 * force for that target is not computable from the nodes in hand.  The walk cannot report this
 * where it happens -- it runs threaded over targets, one incompleteness can involve many nodes,
 * and endrun() only REQUESTS a stop and returns, so the walk keeps going.  So every walk records
 * here instead, host and device, and the caller speaks once for the whole pass.
 *
 * The count is read before the pass is resolved, because whether this is repairable is a
 * collective question: gravity_tree() reduces it across ranks and only then decides between
 * rebuilding the tree and redoing the evaluation, or stopping.  The host walks also keep the first
 * case as an example; the device walks can only count, so an example is not guaranteed. */
static long long IncompleteImportCount = 0;
/* Counted separately because it answers a different question: not "is the import short" but "can a
 * rebuild fix it".  An essential subtree whose children the sender never owned is short for a
 * topological reason, so rebuilding reproduces it and the repair would only burn a build. */
static long long UnshippableImportCount = 0;
static int    IncompleteImportNode = -1, IncompleteImportType = -1;
static int    IncompleteImportHaveExample = 0;
static unsigned long long IncompleteImportID = 0;
static double IncompleteImportLen = 0, IncompleteImportMass = 0;

void gravity_note_incomplete_import(int node, unsigned long long id, int ptype, double len, double mass)
{
#ifdef _OPENMP
#pragma omp critical(_incomplete_import_)
#endif
    {
        if(!IncompleteImportHaveExample)
        {
            IncompleteImportNode = node; IncompleteImportID = id; IncompleteImportType = ptype;
            IncompleteImportLen = len;   IncompleteImportMass = mass;
            IncompleteImportHaveExample = 1;
        }
        IncompleteImportCount++;
    }
}

/* For the device walks, which count inside the kernel and have no per-node detail to carry out. */
void gravity_note_incomplete_import_count(long long n)
{
    if(n <= 0) {return;}
#ifdef _OPENMP
#pragma omp critical(_incomplete_import_)
#endif
    {
        IncompleteImportCount += n;
    }
}

/* For both host and device walks: how many of the shortfalls were of the unrepairable kind. */
void gravity_note_unshippable_import(long long n)
{
    if(n <= 0) {return;}
#ifdef _OPENMP
#pragma omp critical(_incomplete_import_)
#endif
    {
        UnshippableImportCount += n;
    }
}

long long gravity_unshippable_import_count(void) {return UnshippableImportCount;}

long long gravity_incomplete_import_count(void) {return IncompleteImportCount;}

/* One example of what this rank could not descend, for the single line the caller prints.  Returns
 * whether there is one: the device walks can only count, so a rank may have a shortfall and no
 * example, and the caller picks a rank that has one. */
int gravity_incomplete_import_example(char *buf, int buflen)
{
    if(!IncompleteImportHaveExample || buflen <= 0) {return 0;}
    snprintf(buf, (size_t) buflen, "node %d, target ID=%llu type=%d, node size %g mass %g.",
             IncompleteImportNode, IncompleteImportID, IncompleteImportType,
             IncompleteImportLen, IncompleteImportMass);
    return 1;
}

void gravity_clear_incomplete_import(void)
{
    UnshippableImportCount = 0;
    IncompleteImportCount = 0;
    IncompleteImportHaveExample = 0;
}

/* How many accepted elements a walk holds between traversal and evaluation. The run lives
 * in the walking thread's workspace and is evaluated whenever it fills, so its size bounds
 * memory, never the number of interactions; a small run only evaluates sooner. */
#define GRAVTREE_WALK_RECORD_CHUNK 256

/* One target of a packet walk: its inputs to the opening decision and to the pair evaluation,
 * fixed for the walk, and everything the walk accumulates for it. The packet's members share
 * one traversal of the tree; each member keeps its own opening decisions, so its accepted
 * elements are exactly those its own walk would accept, in the same order. */
struct grav_walk_open_inputs_t {
    Vec3<double> pos; int ptype; double soft, aold;
#ifdef PMGRID
    double rcut, rcut2;
#endif
};
struct grav_walk_member_t {
    int target;                 /* particle index */
    int resume_at;              /* GRAV_WALK_MEMBER_IN_PLAY; the index at which this member re-joins the traversal; or GRAV_WALK_MEMBER_NEVER for a massless target that takes part in nothing and writes nothing */
    struct grav_walk_open_inputs_t open;   /* what the opening decision reads at every node */
    double pmass, zeta;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    Vec3<double> vel;
#endif
    grav_pair_tgt_t tgt;
#ifdef RT_USE_GRAVTREE
    int valid_gas_particle_for_rt;
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    double fac_stellum[N_RT_FREQ_BINS];
#endif
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    double r_for_total_menclosed, m_enc_in_rcrit;
#endif
    /* accumulators */
    grav_pair_acc_t out;
#ifdef RT_USE_TREECOL_FOR_NH
    double treecol_angular_bins[RT_USE_TREECOL_FOR_NH];
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    MyFloat tree_mass;
#endif
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_flux_G0[CHIMES_LOCAL_UV_NBINS], chimes_flux_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#ifdef RT_OTVET
    SymmetricTensor2<double> RT_ET[N_RT_FREQ_BINS];
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    double incident_flux_uv, incident_flux_euv;
#endif
#ifdef SINK_COMPTON_HEATING
    double incident_flux_agn;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double SubGrid_CosmicRayEnergyDensity;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    double Rad_E_gamma[N_RT_FREQ_BINS];
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    Vec3<double> Rad_Flux[N_RT_FREQ_BINS];
#endif
#ifdef SINK_CALC_DISTANCES
    grav_sink_prox_accum_t sink_prox;
#endif
};
#define GRAV_WALK_MEMBER_IN_PLAY (-2)   /* resume_at value of a member taking part in the traversal (-1 is a valid tree index meaning "end of walk") */
#define GRAV_WALK_MEMBER_NEVER   (-3)   /* resume_at value of a member that never takes part */

/* Per-thread workspace of a packet walk: the members, the run of accepted records, and one
 * mask word set per record naming the members that accepted it. Members and masks are sized
 * from the packet capacity the caller chose, so a large configured packet size costs nothing
 * on a rank with one host target. */
#define GRAV_WALK_MASK_BITS 64
static inline int grav_walk_mask_words(int cap) {return (cap + GRAV_WALK_MASK_BITS - 1) / GRAV_WALK_MASK_BITS;}   /* per record: one bit per member the workspace can hold */
static inline size_t grav_walk_round64(size_t n) {return (n + 63) & ~((size_t) 63);}
size_t force_treewalk_workspace_bytes_per_thread(int cap)
{
    if(cap < 1) {cap = 1;}
    return grav_walk_round64(cap * sizeof(struct grav_walk_member_t))
         + grav_walk_round64(GRAVTREE_WALK_RECORD_CHUNK * sizeof(grav_walk_record_t))
         + grav_walk_round64((size_t) GRAVTREE_WALK_RECORD_CHUNK * grav_walk_mask_words(cap) * sizeof(unsigned long long))
         + grav_walk_round64(2 * grav_walk_mask_words(cap) * sizeof(unsigned long long))
#if defined(RT_USE_GRAVTREE) || defined(COSMIC_RAY_SUBGRID_LEBRON)
         + grav_walk_round64(GRAVTREE_WALK_RECORD_CHUNK * sizeof(struct gravtree_source_inputs_t))
         + grav_walk_round64(GRAVTREE_WALK_RECORD_CHUNK)
#endif
         ;
}

/* Walk the tree once for a packet of targets and accumulate every target's tree force.
 *
 * The traversal is shared: a node is descended when any member needs it opened. The
 * decisions are not shared: at every node each member still taking part applies its own
 * opening criterion (its own separation, softening, acceleration and type), and a member
 * that accepts or passes a node the packet then descends sits out until the traversal
 * reaches that node's sibling, which in this tree's depth-first order is the first index
 * visited after the node's subtree. Each member therefore judges exactly the nodes its
 * own walk would visit, in the same order, and its accepted elements are the same sequence
 * its own walk would accumulate. Terminal elements (particles, imported leaves, imported
 * aggregates that arrived without children) never descend and are judged per member as
 * before.
 *
 * Every source is drifted during the traversal, once, before any record naming it is
 * appended, and drifting is idempotent within a step, so the state an element is evaluated
 * at is the state it had when it was accepted, whichever member accepted it.
 *
 * Nothing a packet of several targets produces reaches the particles, the cost counters or
 * the import diagnostics until the whole packet has completed. A packet that meets a
 * pseudo-particle returns 0 having written nothing, and the caller walks each of its
 * members alone: a packet of one records the pseudo-particle as before and returns 1, or -1
 * when the detector table is full (then, as before, nothing is written for that target but
 * the diagnostics it noted stand). On 1, ninter_out[m] holds each member's interaction count.
 * workspace is this thread's block of force_treewalk_workspace_bytes_per_thread(cap) bytes;
 * n_targets <= cap. */
int force_treeevaluate(const int *targets, int n_targets, int cap, int *ninter_out, void *workspace,
                       int *exportflag, int *exportnodecount, int *exportindex)
{
    struct NODE *nop = 0;
    int no, nexp, task, treeBase = All.TreeNodeIndexBase, treeSlots = All.TreeParticleSlots;
    long bunchSize = All.BunchSize; int maxNodes = MaxNodes; int maxForeignNodes = MaxForeignNodes; integertime ti_Current = All.Ti_Current;    /* maxForeignNodes shifts pseudo-particle range above the foreign-node range */
    const int target = targets[0];   /* the target the import detector records for a packet of one (n_targets == 1 is the only case that reaches it) */
#ifdef RT_USE_TREECOL_FOR_NH
    double angular_bin_size = 4*M_PI / RT_USE_TREECOL_FOR_NH;
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    double gasmass;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    double center[3]={0};
#ifdef BOX_PERIODIC
    center[0] = 0.5 * boxSize_X; center[1] = 0.5 * boxSize_Y; center[2] = 0.5 * boxSize_Z;
#endif
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double cr_injection = 0;
    /* per-call CR gate + max stream time, hoisted from the accumulation (the time helper is host-only; value is interaction-independent) */
    int cr_active_gate = (All.Time > All.TimeBegin) ? 1 : 0; double cr_t_max = 0;
    if(cr_active_gate) {cr_t_max = DMIN(1., evaluate_time_since_t_initial_in_Gyr(All.TimeBegin))/UNIT_TIME_IN_GYR;}
#endif
#ifdef RT_USE_GRAVTREE
    double mass_stellarlum[N_RT_FREQ_BINS]; int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {mass_stellarlum[k_freq]=0;}
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_mass_stellarlum_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_mass_stellarlum_ion[CHIMES_LOCAL_UV_NBINS]={0};
#endif
#endif

    /* the workspace: the members, the run of records and their masks, and where compiled the
       per-run payload cache */
    if(cap < 1) {cap = 1;}
    if(n_targets < 1 || n_targets > cap || cap > TREE_QUERY_PACKET_SIZE) {endrun(90001054); return -1;}
    /* the members taking part; with a packet size of one this is the constant 1, so every loop
       over members below collapses to its single iteration, and the per-record bookkeeping that
       only several members need (which members accepted a record; whether a leaf's payload has
       already been derived for another member) is left out at compile time */
    constexpr int packet_of_one = (TREE_QUERY_PACKET_SIZE == 1);
    const int n_members = packet_of_one ? 1 : n_targets;
    const int mask_words = packet_of_one ? 1 : grav_walk_mask_words(cap);   /* sized by what this call's workspace holds, not by the configured packet size */
    char *ws = (char *) workspace;
    struct grav_walk_member_t *members = (struct grav_walk_member_t *) ws; ws += grav_walk_round64(cap * sizeof(struct grav_walk_member_t));
    grav_walk_record_t *records = (grav_walk_record_t *) ws; ws += grav_walk_round64(GRAVTREE_WALK_RECORD_CHUNK * sizeof(grav_walk_record_t));
    unsigned long long *masks = (unsigned long long *) ws; ws += grav_walk_round64((size_t) GRAVTREE_WALK_RECORD_CHUNK * mask_words * sizeof(unsigned long long));
    /* the node's scratch masks: one word each for a packet of one (a register), otherwise a word per
       64 members the workspace holds, kept in the workspace */
    unsigned long long single_accept_mask[1], single_opened_mask[1];
    unsigned long long *accept_mask = packet_of_one ? single_accept_mask : (unsigned long long *) ws;
    unsigned long long *opened_mask = packet_of_one ? single_opened_mask : (unsigned long long *) ws + mask_words;
    ws += grav_walk_round64(2 * mask_words * sizeof(unsigned long long));
#if defined(RT_USE_GRAVTREE) || defined(COSMIC_RAY_SUBGRID_LEBRON)
    /* a leaf's RT / sink / CR payload, derived once per run of records for the members that consume it */
    struct gravtree_source_inputs_t *payloads = (struct gravtree_source_inputs_t *) ws; ws += grav_walk_round64(GRAVTREE_WALK_RECORD_CHUNK * sizeof(struct gravtree_source_inputs_t));
    unsigned char *payload_loaded = (unsigned char *) ws;
#endif
    int n_records = 0;
#define GRAV_WALK_MASK_SET(mask, m)  ((mask)[(m) / GRAV_WALK_MASK_BITS] |= (1ULL << ((m) % GRAV_WALK_MASK_BITS)))
#define GRAV_WALK_MASK_TEST(mask, m) (((mask)[(m) / GRAV_WALK_MASK_BITS] >> ((m) % GRAV_WALK_MASK_BITS)) & 1ULL)

    /* Import-completeness notes are held here until the packet completes (a packet of several
     * that is abandoned for a pseudo-particle must leave no trace of the walk it did not
     * finish), then handed to the shared counters through the same calls as before. */
    struct { long long count, unshippable; int have_example, node, ptype; unsigned long long id; double len, mass; } notes = {0, 0, 0, 0, 0, 0ULL, 0.0, 0.0};
    auto commit_notes = [&]()
    {
        if(notes.have_example) {gravity_note_incomplete_import(notes.node, notes.id, notes.ptype, notes.len, notes.mass);}
        gravity_note_incomplete_import_count(notes.count - (notes.have_example ? 1 : 0));
        gravity_note_unshippable_import(notes.unshippable);
    };

    /* the members: each target's inputs to the walk, and zeroed accumulators */
    for(int m = 0; m < n_members; m++)
    {
        struct grav_walk_member_t &mem = members[m];
        const int tm = targets[m];
        mem.target = tm; mem.resume_at = GRAV_WALK_MEMBER_IN_PLAY;
        mem.zeta = 0;
        mem.open.pos = P[tm].Pos;
        mem.open.ptype = P[tm].Type;
        mem.open.soft = ForceSoftening_KernelRadius(tm);
        mem.open.aold = All.ErrTolForceAcc * P[tm].OldAcc;
        mem.pmass = P[tm].Mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
        mem.vel = P[tm].Vel;
#endif
        double sink_mass = 0; (void) sink_mass;
#if defined(SINK_DYNFRICTION_FROMTREE)
        if(mem.open.ptype == 5) {sink_mass = P[tm].Sink_Mass;}
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        grav_target_select_soft_and_zeta(mem.open.ptype, P[tm].AGS_zeta, mem.open.soft, mem.zeta);
#endif
#ifdef PMGRID
        double rcut = All.Rcut[0], asmth = All.Asmth[0];
#if defined(PM_PLACEHIGHRESREGION)
        if(pmforce_is_particle_high_res(mem.open.ptype, P[tm].Pos)) {rcut = All.Rcut[1]; asmth = All.Asmth[1];}
#endif
        mem.open.rcut = rcut; mem.open.rcut2 = rcut * rcut;
#endif
        if(mem.pmass <= 0) {mem.resume_at = GRAV_WALK_MEMBER_NEVER;} /* a massless target: nothing to compute, nothing to write */
        int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(mem.open.ptype); // determine allowed particle types for correction terms for adaptive gravitational softening terms
        /* read-only PM short-range config for the shared force helpers (empty when !PMGRID;
         * built once per target after the PM_PLACEHIGHRESREGION rcut/asmth override above). */
        grav_pm_shortrange_t pm{};
#ifdef PMGRID
        pm.rcut = rcut; pm.rcut2 = mem.open.rcut2; pm.asmthfac = grav_pm_asmthfac(asmth); pm.shortrange_tab = shortrange_table;
#ifdef EVALPOTENTIAL
        pm.shortrange_pot_tab = shortrange_table_potential;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        pm.shortrange_tidal_tab = shortrange_table_tidal;
#endif
#endif
        /* the target's inputs to the shared pair evaluation, fixed for this walk */
        grav_pair_tgt_t tgt{}; tgt.ptype = mem.open.ptype; tgt.pmass = mem.pmass; tgt.h = mem.open.soft; tgt.zeta = mem.zeta; tgt.ags_bitflag = AGS_kernel_shared_BITFLAG; tgt.pm = pm;
#ifdef SINK_DYNFRICTION_FROMTREE
        tgt.sink_mass = sink_mass;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        tgt.i_zeta_tidal_tensorps_prevstep = P[tm].tidal_tensorps_prevstep;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
        tgt.pos = mem.open.pos; tgt.center[0] = center[0]; tgt.center[1] = center[1]; tgt.center[2] = center[2];
#endif
        mem.tgt = tgt;
#ifdef RT_USE_GRAVTREE
        mem.valid_gas_particle_for_rt = grav_target_valid_gas_for_rt(mem.open.ptype, mem.open.soft, mem.pmass);
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        if(mem.valid_gas_particle_for_rt)
        {
            double kappa_eff[N_RT_FREQ_BINS]; int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {kappa_eff[kf] = rt_kappa(-1,kf, P, CellP);} // rt_kappa is in physical code units (needs the walk's particle pointers, so evaluated here)
            grav_target_rt_fac_stellum(mem.open.soft, mem.pmass, kappa_eff, mem.fac_stellum);
        }
#endif
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
        mem.m_enc_in_rcrit = 0; mem.r_for_total_menclosed = grav_target_menc_radius(mem.open.soft); /* baseline Rcrit_min applied in the helper, otherwise we get statistics that are very noisy */
#endif
        /* the accumulators the shared pair evaluation writes (acceleration, potential, interaction
           count, and the tidal / jerk / tidal-zeta terms where compiled), and the walker's own */
        grav_pair_acc_init(mem.out);
#ifdef RT_USE_TREECOL_FOR_NH
        {int k; for(k=0; k<RT_USE_TREECOL_FOR_NH; k++) {mem.treecol_angular_bins[k] = 0;}}
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
        mem.tree_mass = 0;
#endif
#ifdef CHIMES_STELLAR_FLUXES
        {int kc; for(kc=0; kc<CHIMES_LOCAL_UV_NBINS; kc++) {mem.chimes_flux_G0[kc] = 0; mem.chimes_flux_ion[kc] = 0;}}
#endif
#ifdef RT_OTVET
        {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {mem.RT_ET[kf] = {};}}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        mem.incident_flux_uv = 0; mem.incident_flux_euv = 0;
#endif
#ifdef SINK_COMPTON_HEATING
        mem.incident_flux_agn = 0;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        mem.SubGrid_CosmicRayEnergyDensity = 0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {mem.Rad_E_gamma[kf] = 0;}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        {int kf; for(kf=0; kf<N_RT_FREQ_BINS; kf++) {mem.Rad_Flux[kf] = {};}}
#endif
#ifdef SINK_CALC_DISTANCES
        grav_sink_prox_accum_init(mem.sink_prox); /* nearest-sink + single-star timestep/binary accumulators (gravtree_force_kernel.h) */
#endif
    }


    /* Evaluate a run of accepted elements for every member that accepted each, member by member:
       a member's inputs and accumulators are held in locals for the whole run (the compiler keeps
       them in registers, which is what the pair evaluation's cost depends on), and each element
       it accepted is loaded again from the tree, from its index alone, in acceptance order: a
       particle leaf from P[], a node from Nodes[]/Extnodes[] and the foreign-leaf sidecars. A
       leaf's luminosity / cosmic-ray payload is derived once per run and shared by the members
       that consume it. Every source drift happened during the traversal that accepted the
       element, and drifting is idempotent within a step, so the state read here is the state
       the element had when it was accepted; a Hermite-owned source is predicted from that state,
       as it was before the traversal and the evaluation were separated. A source found not
       drifted, or a node whose classification no longer matches its record, is a broken
       invariant and stops the run. */
    auto evaluate_records = [&](const grav_walk_record_t *rec, const unsigned long long *rec_masks, int n_rec)
    {
#if defined(RT_USE_GRAVTREE) || defined(COSMIC_RAY_SUBGRID_LEBRON)
        if(!packet_of_one) {for(int irec = 0; irec < n_rec; irec++) {payload_loaded[irec] = 0;}}
#endif
        for(int m = 0; m < n_members; m++)
        {
            struct grav_walk_member_t &mem = members[m];
            if(mem.resume_at == GRAV_WALK_MEMBER_NEVER) {continue;}
            if(!packet_of_one)
            {   /* nothing in this run for a member that accepted none of it */
                unsigned long long any = 0;
                for(int irec = 0; irec < n_rec; irec++) {any |= rec_masks[(size_t) irec * mask_words + m / GRAV_WALK_MASK_BITS] & (1ULL << (m % GRAV_WALK_MASK_BITS));}
                if(!any) {continue;}
            }
            /* the member's inputs, fixed for the walk */
            const Vec3<double> pos = mem.open.pos; const int ptype = mem.open.ptype; const double soft = mem.open.soft, pmass = mem.pmass; (void) soft; (void) pmass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
            const Vec3<double> vel = mem.vel;
#endif
            const grav_pair_tgt_t tgt = mem.tgt;
#ifdef RT_USE_GRAVTREE
            const int valid_gas_particle_for_rt = mem.valid_gas_particle_for_rt;
            Vec3<double> d_stellarlum = {};
#ifdef SINK_PHOTONMOMENTUM
            double mass_sinklumwt_forradfb=0; // convert bh luminosity to our tree units
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
            Vec3<double> d_dm = {}; double mass_dm = 0;
#endif
            /* the member's accumulators, taken up for this run and put back after it */
            grav_pair_acc_t out = mem.out;
#ifdef COUNT_MASS_IN_GRAVTREE
            MyFloat tree_mass = mem.tree_mass;
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
            double m_enc_in_rcrit = mem.m_enc_in_rcrit; const double r_for_total_menclosed = mem.r_for_total_menclosed;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
            double incident_flux_uv = mem.incident_flux_uv, incident_flux_euv = mem.incident_flux_euv;
#endif
#ifdef SINK_COMPTON_HEATING
            double incident_flux_agn = mem.incident_flux_agn;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            double SubGrid_CosmicRayEnergyDensity = mem.SubGrid_CosmicRayEnergyDensity;
#endif
#ifdef SINK_CALC_DISTANCES
            grav_sink_prox_accum_t sink_prox = mem.sink_prox;
#endif

            for(int irec = 0; irec < n_rec; irec++)
            {
                if(!packet_of_one && !GRAV_WALK_MASK_TEST(rec_masks + (size_t) irec * mask_words, m)) {continue;}
                const int no = rec[irec].no;
                grav_pair_src_t src;   /* dr, r2, mass are assigned on every path that reaches the evaluation */
                Vec3<double> &dr = src.dr; double &r2 = src.r2, &mass = src.mass, &h_p = src.h_p, &zeta_sec = src.zeta_sec; int &ptype_sec = src.ptype_sec;
                h_p = -1; ptype_sec = -1; zeta_sec = 0;
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
                Vec3<double> &dv = src.dv;
#endif
#if defined(SINK_DYNFRICTION_FROMTREE)
                double &m_j_eff_for_df = src.m_j_eff_for_df;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
                double &r_source = src.r_source;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                SymmetricTensor2<double> &j_zeta_tidal_tensorps_prevstep = src.j_zeta_tidal_tensorps_prevstep;
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                gasmass=0; /* reset per interaction: non-gas leaf sources carry NO gas mass. Without this, a
                            * non-gas leaf inherits the stale gasmass of an earlier source (or garbage before the
                            * first assignment) and the TREECOL column estimate fabricates contributions through
                            * dark-matter/star particles. Nodes assign unconditionally below; only gas leaves
                            * (and the sink alpha-disk reservoir, where enabled) carry gas mass. */
#endif
                if(no < treeSlots) /* this is a particle, we will use it */
                {
                    if(P[no].Ti_current != ti_Current) {endrun(90001051);}   /* accepted during traversal, which drifted it */
                    /* the source state this interaction is evaluated at, which is the drifted state
                       except where the Hermite predictor below replaces it */
                    Vec3<double> src_pos = P[no].Pos;
                    Vec3<double> src_vel = P[no].Vel;   /* unconditional: the sink-proximity block below reads it under SINK_CALC_DISTANCES, which several flags reach without the jerk or dynamical-friction terms */
#ifdef HERMITE_INTEGRATION
                    /* On a Hermite pass a source the Hermite integrator owns but is not advancing this
                       step is second-order wrong where it stands; evaluate it from its own start-of-step
                       state instead. Single sources only: one absorbed into a node multipole still
                       contributes from the node's drifted centre of mass. Under
                       SINGLE_STAR_DIRECT_GRAVITY_RADIUS the close star pairs this matters most for are
                       force-opened to singles and so do take this branch. Nothing is written back. */
                    if(hermite_source_needs_prediction(no, P, HermiteWalk)) {
                        hermite_predict_source_state(no, P, HermiteWalk, &HermiteWalkTables, src_pos, src_vel);
                    }
#endif
                    dr = src_pos - pos;
                    GRAVITY_NEAREST_XYZ(dr[0],dr[1],dr[2],-1);
                    r2 = dr.norm_sq();
                    mass = P[no].Mass;

#ifdef GRAVITY_SPHERICAL_SYMMETRY
                    r_source = grav_spherical_symmetry_r_from_center(src_pos[0],src_pos[1],src_pos[2],center[0],center[1],center[2]);
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
                    dv = src_vel - vel;
#endif
#if defined(SINK_DYNFRICTION_FROMTREE)
                    m_j_eff_for_df = mass;
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                    if(P[no].Type == 0) {gasmass = P[no].Mass;}
#if defined(SINK_ALPHADISK_ACCRETION) && defined(RT_USE_TREECOL_FOR_NH)
                    if(P[no].Type == 5) {gasmass = P[no].Sink_Mass_Reservoir;} // gas at the inner edge of a disk should not see a hole due to the sink
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                    j_zeta_tidal_tensorps_prevstep=P[no].tidal_tensorps_prevstep;
#endif

                    /* only proceed if the mass is positive and there is separation! */
                    if((r2 > 0) && (mass > 0))
                    {

#ifdef SINK_CALC_DISTANCES
                        /* nearest-sink + single-star timestep/binary tracking via the shared helper (gravtree_force_kernel.h) */
                        grav_sink_prox_target_t prox_target = {}; prox_target.ptype = ptype; prox_target.pmass = pmass; prox_target.soft = soft;
#if defined(SINGLE_STAR_TIMESTEPPING)
                        prox_target.vel = vel;
#endif
                        grav_sink_prox_leaf_src_t prox_src = {}; prox_src.src_type = P[no].Type; prox_src.src_mass = P[no].Mass; prox_src.motion.vel = src_vel;   /* the state this interaction was evaluated at, so the pair (dr, vel) feeding Min_Sink_Approach_Time stays mutually consistent on a Hermite pass */
#if defined(SPECIAL_POINT_MOTION) || defined(SPECIAL_POINT_WEIGHTED_MOTION)
                        prox_src.motion.acc = P[no].Acc_Total_PrevStep;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
                        prox_src.motion.max_feedback_vel = P[no].MaxFeedbackVel;
#endif
                        grav_sink_prox_leaf_accumulate(r2, dr, prox_target, prox_src, sink_prox);
#endif // SINK_CALC_DISTANCES

#if defined(RT_USE_GRAVTREE) || defined(COSMIC_RAY_SUBGRID_LEBRON)
                        /* the source's RT / sink / CR payload through the shared gates (gravtree_moment_sources.h),
                           from the drifted state; derived once per run for the members that consume it, so a
                           target that is not a valid RT receiver does no luminosity work it never did */
                        int need_source_payload = 0;
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                        need_source_payload = 1;
#else
                        need_source_payload = valid_gas_particle_for_rt;
#endif
                        const struct gravtree_source_inputs_t &source_payload = payloads[irec];
                        if(need_source_payload && (packet_of_one || !payload_loaded[irec])) {gravtree_fill_particle_source_inputs(no, P, CellP, &payloads[irec]); if(!packet_of_one) {payload_loaded[irec] = 1;}}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                        cr_injection = source_payload.cr_inject;
#endif

#ifdef RT_USE_GRAVTREE
                        if(valid_gas_particle_for_rt)    /* we have a (valid) gas particle as target */
                        {
                            d_stellarlum=dr;
                            int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(source_payload.rt_active) {mass_stellarlum[kf]=source_payload.src_lum[kf];} else {mass_stellarlum[kf]=0;}}
#ifdef CHIMES_STELLAR_FLUXES
                            for(kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++)
                            {
                                if(source_payload.rt_active) {chimes_mass_stellarlum_G0[kf] = source_payload.src_lum_G0[kf]; chimes_mass_stellarlum_ion[kf] = source_payload.src_lum_ion[kf];} else {chimes_mass_stellarlum_G0[kf] = 0; chimes_mass_stellarlum_ion[kf] = 0;}
                            }
#endif
#ifdef SINK_PHOTONMOMENTUM
                            mass_sinklumwt_forradfb=0;
                            if(P[no].Type == 5)
                            {
                                double bhlum_t = source_payload.bh_active ? (double) source_payload.bh_lum : 0.0;
                                Vec3<double> bh_angle = source_payload.bh_active ? Vec3<double>{(double) source_payload.bh_angle[0], (double) source_payload.bh_angle[1], (double) source_payload.bh_angle[2]} : Vec3<double>{0,0,0};
                                mass_sinklumwt_forradfb = sink_fb_angleweight(bhlum_t, bh_angle, dr[0],dr[1],dr[2]);
                            }
#endif
                        }
#endif // RT_USE_GRAVTREE

#ifdef DM_SCALARFIELD_SCREENING
                        if(ptype != 0) {if(P[no].Type == 1) {d_dm = dr; mass_dm = mass;} else {d_dm = {}; mass_dm = 0;}} /* we have a dark matter particle as target */
#endif

                        h_p = ForceSoftening_KernelRadius(no);
                        ptype_sec=P[no].Type; zeta_sec=0; /* set secondary softening and zeta term */
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                        if(ptype_sec==0) {zeta_sec=P[no].AGS_zeta;}
#elif defined(ADAPTIVE_GRAVSOFT_FORALL)
                        zeta_sec=P[no].AGS_zeta;
#endif
                    } // closes (if((r2 > 0) && (mass > 0))) check

                }
                else /* we have an internal node the traversal accepted */
                {
                    struct NODE *nop = &Nodes[no];
                    if(nop->Ti_current != ti_Current) {endrun(90001052);}   /* accepted during traversal, which drifted it */
                    int in_foreign = (no >= treeBase + maxNodes && no < treeBase + maxNodes + maxForeignNodes);
                    /* Foreign-leaf identity lookup (host sidecar; foreign_slot = no-(treeBase+maxNodes),
                     * EXPLICIT and bounds-checked -- not the node index no-treeBase). */
                    int    fl_tag = 0, fl_type = -1;
                    double fl_zeta = 0.0, fl_soft = 0.0;
                    if(in_foreign && ForeignLeafTag) {
                        int fs = no - (treeBase + maxNodes);
                        if(fs >= 0 && fs < AllocatedForeignNodes) {
                            fl_tag  = ForeignLeafTag[fs];
                            fl_type = ForeignLeafType[fs];
                            fl_zeta = (double) ForeignLeafZeta[fs];
                            fl_soft = (double) ForeignLeafSoft[fs];
                        }
                    }
                    if(fl_tag != rec[irec].leaf_tag || grav_classify_node(in_foreign, fl_tag, nop->u.d.nextnode) != rec[irec].kind) {endrun(90001053);}
                    mass = nop->u.d.mass;

#ifdef SINGLE_STAR_DIRECT_GRAVITY
                    /* Remove the sinks from this node for a star target: star-star pairs come exactly from
                       star_direct_gravity_compute(), so taking them here too would double them. This tree
                       carries monopoles only (u.d.mass at u.d.s -- struct NODE has no quadrupole moments),
                       so the subtraction is exact rather than approximate: drop the sink mass and move the
                       center of mass to that of what remains. Both terms are on the same clock, since
                       SINK_NODE_MOTION_TRACKED drifts sink_pos with sink_vel exactly as u.d.s is drifted
                       with vs. The traversal judged the opening criteria on this reduced mass too. */
                    if((ptype == 5) && (nop->sink_mass > 0))
                    {
                        double mass_nosink = mass - nop->sink_mass;
                        dr = (nop->u.d.s * mass - nop->sink_pos * nop->sink_mass) / mass_nosink - pos;
                        mass = mass_nosink;
                    }
                    else {dr = nop->u.d.s - pos;}
#else
                    dr = nop->u.d.s - pos;
#endif
                    GRAVITY_NEAREST_XYZ(dr[0],dr[1],dr[2],-1);
                    r2 = dr.norm_sq();

                    /* ok we will be using this node, can now set variables that depend on it */
                    h_p = nop->maxsoft;
                    zeta_sec = 0; ptype_sec = -1; /* set secondary softening and zeta terms */
                    /* A tagged real foreign single-particle leaf is consumed with particle-leaf
                     * secondary semantics -- restore the Type + AGS_zeta the node moment cannot carry,
                     * via the shared seam (identical to the GPU walk). */
                    if(fl_tag == 1) { grav_apply_foreign_leaf_identity(fl_tag, fl_type, fl_zeta, fl_soft, &ptype_sec, &zeta_sec, &h_p); }
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                    gasmass = nop->gasmass;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
                    r_source = grav_spherical_symmetry_r_from_center(nop->u.d.s[0],nop->u.d.s[1],nop->u.d.s[2],center[0],center[1],center[2]);
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
                    dv = Extnodes[no].vs - vel;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                    cr_injection = nop->cr_injection;
#endif

#ifdef RT_USE_GRAVTREE
                    if(valid_gas_particle_for_rt)    /* we have a (valid) gas particle as target */
                    {
                        int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {mass_stellarlum[kf] = nop->stellar_lum[kf];}
#ifdef CHIMES_STELLAR_FLUXES
                        for(kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++)
                        {
                            chimes_mass_stellarlum_G0[kf] = nop->chimes_stellar_lum_G0[kf];
                            chimes_mass_stellarlum_ion[kf] = nop->chimes_stellar_lum_ion[kf];
                        }
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
                        d_stellarlum = nop->rt_source_lum_s - pos;
                        GRAVITY_NEAREST_XYZ(d_stellarlum[0],d_stellarlum[1],d_stellarlum[2],-1);
#else
                        d_stellarlum = dr;
#endif
#ifdef SINK_PHOTONMOMENTUM
                        mass_sinklumwt_forradfb = sink_fb_angleweight(nop->sink_lum, nop->sink_lum_grad, d_stellarlum[0],d_stellarlum[1],d_stellarlum[2]);
#endif
                    }
#endif // RT_USE_GRAVTREE

#ifdef DM_SCALARFIELD_SCREENING
                    if(ptype != 0) {d_dm = nop->s_dm - pos; mass_dm = nop->mass_dm;} else {d_dm = {}; mass_dm = 0;} /* we have a dark matter particle as target */
#endif
#if defined(SINK_DYNFRICTION_FROMTREE)
                    m_j_eff_for_df = (nop->u.d.mass) / (nop->N_part);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                    j_zeta_tidal_tensorps_prevstep=nop->tidal_tensorps_prevstep;
#endif

#ifdef SINK_CALC_DISTANCES // NOTE: moved this to AFTER the checks for node opening, because we only want to record BH positions from the nodes that actually get used for the force calculation - MYG
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
                    grav_sink_prox_node_specialweighted(r2, Extnodes[no].vs, ptype, sink_prox);
#endif
                    if(nop->sink_mass > 0)        /* found a node with non-zero BH mass */
                    {
                        Vec3<double> sink_dr = nop->sink_pos - pos;  /* SHEA:  now using sink_pos instead of center */
                        GRAVITY_NEAREST_XYZ(sink_dr[0],sink_dr[1],sink_dr[2],-1);
                        grav_sink_prox_target_t prox_target = {}; prox_target.ptype = ptype; prox_target.pmass = pmass; prox_target.soft = soft;
#if defined(SINGLE_STAR_TIMESTEPPING)
                        prox_target.vel = vel;
#endif
                        grav_sink_prox_node_src_t prox_src = {}; prox_src.sink_mass = nop->sink_mass;
#if defined(SINGLE_STAR_FIND_BINARIES)
                        prox_src.n_sink = (int)nop->N_SINK;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SPECIAL_POINT_MOTION)
                        prox_src.motion.vel = nop->sink_vel;
#endif
#if defined(SPECIAL_POINT_MOTION)
                        prox_src.motion.acc = nop->sink_acc;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
                        prox_src.motion.max_feedback_vel = nop->MaxFeedbackVel;
#endif
                        grav_sink_prox_node_accumulate(r2, sink_dr, prox_src, prox_target, sink_prox);
                    }
#endif // SINK_CALC_DISTANCES

                } /* the node's inputs are loaded */


                if((r2 > 0) && (mass > 0)) // only go forward if mass positive and there is separation -- this is check for the whole block below, which should no include 'self' terms
                {
#if defined(EVALPOTENTIAL) && defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
                    /* periodic-image potential correction, from the separation before the shared
                       evaluation's spherical-symmetry override; added right after the pair potential */
                    double pot_periodic_image = mass * ewald_pot_corr(dr[0], dr[1], dr[2]);
#endif
                    /* pair-wise gravity, PM truncation, and the accumulations inside the PM short-range
                     * gate, via the shared evaluation (gravtree_force_kernel.h), the single home for the
                     * pair physics on both walks */
                    grav_pair_result_t res = grav_pair_evaluate_core(tgt, src, out);
                    const double r = res.r, fac_accel = res.fac_accel; (void) r; (void) fac_accel;
#if defined(EVALPOTENTIAL) && defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
                    out.pot += pot_periodic_image;
#endif

#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
                    if(r < r_for_total_menclosed) {m_enc_in_rcrit += mass;}
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
                    tree_mass += mass;
#endif
#ifdef RT_USE_TREECOL_FOR_NH
                    grav_treecol_accumulate(dr, r, fac_accel, gasmass, mass, angular_bin_size, mem.treecol_angular_bins);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                    grav_cr_lebron_accumulate(ptype, r, soft, cr_injection, cr_active_gate, cr_t_max, tgt.pm, SubGrid_CosmicRayEnergyDensity);
#endif
#ifdef RT_USE_GRAVTREE
                    if(valid_gas_particle_for_rt)    /* we have a (valid) gas particle as target; payload formulas in the shared helper */
                    {
                        grav_rt_src_t rt_src = {}; rt_src.d_stellarlum = d_stellarlum; rt_src.soft = soft; rt_src.mass_stellarlum = mass_stellarlum;
#ifdef CHIMES_STELLAR_FLUXES
                        rt_src.chimes_mass_stellarlum_G0 = chimes_mass_stellarlum_G0; rt_src.chimes_mass_stellarlum_ion = chimes_mass_stellarlum_ion;
#endif
#ifdef SINK_PHOTONMOMENTUM
                        rt_src.mass_sinklumwt_forradfb = mass_sinklumwt_forradfb;
#endif
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
                        rt_src.fac_stellum = mem.fac_stellum;
#endif
                        grav_rt_accum_t rt_accum = {};
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                        rt_accum.Rad_E_gamma = mem.Rad_E_gamma;
#endif
#ifdef CHIMES_STELLAR_FLUXES
                        rt_accum.chimes_flux_G0 = mem.chimes_flux_G0; rt_accum.chimes_flux_ion = mem.chimes_flux_ion;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
                        rt_accum.incident_flux_uv = &incident_flux_uv; rt_accum.incident_flux_euv = &incident_flux_euv;
#endif
#ifdef SINK_COMPTON_HEATING
                        rt_accum.incident_flux_agn = &incident_flux_agn;
#endif
#ifdef RT_OTVET
                        rt_accum.RT_ET = mem.RT_ET;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
                        rt_accum.Rad_Flux = mem.Rad_Flux;
#endif
                        grav_rt_payload_accumulate(rt_src, rt_accum, out.acc);
                    } // closes if(valid_gas_particle_for_rt)

#endif // RT_USE_GRAVTREE


#ifdef DM_SCALARFIELD_SCREENING
                    if(ptype != 0)    /* we have a dark matter particle as target */
                    {
                        grav_dm_scalarfield_accumulate(d_dm, mass_dm, tgt.h, tgt.pm, out.acc);
                    } // closes if(ptype != 0)
#endif // DM_SCALARFIELD_SCREENING //

                } // closes (if((r2 > 0) && (mass > 0))) check
            } // closes the record loop

            /* put the member's accumulators back */
            mem.out = out;
#ifdef COUNT_MASS_IN_GRAVTREE
            mem.tree_mass = tree_mass;
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
            mem.m_enc_in_rcrit = m_enc_in_rcrit;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
            mem.incident_flux_uv = incident_flux_uv; mem.incident_flux_euv = incident_flux_euv;
#endif
#ifdef SINK_COMPTON_HEATING
            mem.incident_flux_agn = incident_flux_agn;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            mem.SubGrid_CosmicRayEnergyDensity = SubGrid_CosmicRayEnergyDensity;
#endif
#ifdef SINK_CALC_DISTANCES
            mem.sink_prox = sink_prox;
#endif
        } // closes the member loop
    };


    /* Traverse the tree once for the packet, recording each accepted element and the members
       that accepted it into a bounded run that is evaluated in acceptance order whenever it
       fills and once more when the traversal ends. The run is private to this call; nothing
       about it outlives the walk of one packet. */
    no = treeBase;        /* root node */
    /* A packet of one holds its member's opening inputs in a local for the whole traversal, so
       they stay in registers across the loop's calls; a packet of several reads each member's at
       the node. One loop body serves both: open_of(m) names whichever applies. */
    const struct grav_walk_open_inputs_t single_open = members[0].open;
    auto open_of = [&](int m) -> const struct grav_walk_open_inputs_t & {return (TREE_QUERY_PACKET_SIZE == 1) ? single_open : members[m].open;};
    int single_resume_at = members[0].resume_at;   /* likewise the member's standing in the traversal */
    auto resume_of = [&](int m) -> int & {return (TREE_QUERY_PACKET_SIZE == 1) ? single_resume_at : members[m].resume_at;};
#ifdef GRAVITY_HYBRID_OPENING_CRIT
    const int pred_is_first_step = (All.Ti_Current == 0 && RestartFlag != 1);
#else
    const int pred_is_first_step = 0;
#endif

    while(1)
    {
        /* one place evaluates the run: when it is full, and once more when the traversal ends */
        if(no < 0 || n_records == GRAVTREE_WALK_RECORD_CHUNK) {evaluate_records(records, masks, n_records); n_records = 0; if(no < 0) {break;}}
        /* a member that sat out a subtree re-joins the traversal at the index it named, before
           anything is decided about that index */
        for(int m = 0; m < n_members; m++) {int &resume_at = resume_of(m); if(resume_at == no) {resume_at = GRAV_WALK_MEMBER_IN_PLAY;}}

        if(no >= treeSlots && no < treeBase) {/* An index between the particle slots and the node base belongs to neither, so the tree is
         * malformed; stop rather than read a side array or Nodes[] out of bounds. */
            endrun(90001024); no = -1; continue;}
        if(no < treeSlots) /* this is a particle, we will use it */
        {
            /* the index of the node is the index of the particle */
            if(P[no].Ti_current != ti_Current)
            {
#ifdef _OPENMP
#pragma omp critical(_particledriftforce_)
#endif
                {
                    if(P[no].Ti_current != ti_Current) {
                        drift_particle(no, ti_Current);
                        gizmo_mark_kernel_radius_dirty_indices(&no, 1);
                    }
                }
            }
            for(int w = 0; w < mask_words; w++) {accept_mask[w] = 0;}
            int any_accept = 0;
            for(int m = 0; m < n_members; m++)
            {
                if(resume_of(m) != GRAV_WALK_MEMBER_IN_PLAY) {continue;}
#ifdef SINGLE_STAR_DIRECT_GRAVITY
                /* star-star pairs are summed exactly in star_direct_gravity_compute(); taking them
                   here as well would double every such force */
                if((open_of(m).ptype == 5) && (P[no].Type == 5)) {continue;}
#endif
                GRAV_WALK_MASK_SET(accept_mask, m); any_accept = 1;
            }
            if(any_accept)
            {
                records[n_records].no = no; records[n_records].kind = GRAV_NODE_LOCAL; records[n_records].leaf_tag = LET_LEAF_TAG_NODE;
                if(!packet_of_one) {for(int w = 0; w < mask_words; w++) {masks[(size_t) n_records * mask_words + w] = accept_mask[w];}}
                n_records++;
            }
        }
        else /* we have an  internal node */
        {
            if(no >= treeBase + maxNodes + maxForeignNodes) /* pseudo particle (foreign-node range below pseudos) -- this will not be used for calculations below, but needs to be parsed here */
            {
                /* A packet of several targets cannot record this per target; it stops here with
                 * nothing written, and the caller walks each member alone. */
                if(n_targets > 1) {return 0;}
                /* LET-incompleteness DETECTOR (not an export system: the MPI round-trip is
                 * retired). Reaching a non-empty pseudo means this target's gravity is not
                 * covered by the local LET; record it so Nexport>0 and gravity_tree() can
                 * raise a graceful controlled-stop. The DataIndexTable/DataNodeList entries
                 * are never shipped -- they only count. */
                if(exportflag[task = DomainTask[no - (treeBase + maxNodes + maxForeignNodes)]] != target)
                {
                    exportflag[task] = target;
                    exportnodecount[task] = NODELISTLENGTH;
                }
                if(exportnodecount[task] == NODELISTLENGTH)
                {
                    int exitFlag = 0;
#ifdef _OPENMP
#pragma omp critical(_nexportforce_)
#endif
                    {
                        if(Nexport >= bunchSize)
                        {
                            /* The table is full, so this target cannot even be recorded. That is the
                             * same failure as recording it -- its gravity is not covered by the local
                             * tree -- so ask for the same controlled stop here. Without this the walk
                             * would return below having silently dropped the targets it had already
                             * taken from the active list, since there is no longer a retry pass to
                             * pick them up again. */
                            BufferFullFlag = 1;
                            exitFlag = 1;
                            gizmo_request_controlled_stop(914040, "gravtree: the locally essential tree did not cover these targets' gravity, and there were too many of them to record", __FILE__, __LINE__, __FUNCTION__);
                        }
                        else
                        {
                            nexp = Nexport;
                            Nexport++;
                        }
                    }
                    if(exitFlag) {commit_notes(); return -1;} /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine. The walk is abandoned: nothing accumulated for this target is written, whether already evaluated or still recorded. */
                    exportnodecount[task] = 0;
                    exportindex[task] = nexp;
                    DataIndexTable[nexp].Task = task;
                    DataIndexTable[nexp].Index = target;
                    DataIndexTable[nexp].IndexGet = nexp;
                }
                DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
                DomainNodeIndex[no - (treeBase + maxNodes + maxForeignNodes)];
                if(exportnodecount[task] < NODELISTLENGTH) {DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;}
                no = Nextnode[treeSlots + (no - treeBase - maxNodes - maxForeignNodes)];
                continue;
            }
            /* ok we have an internal node on the local processor, need to decide if we open it and go further or keep it */
            nop = &Nodes[no];
            int in_foreign = (no >= treeBase + maxNodes && no < treeBase + maxNodes + maxForeignNodes);
            /* Foreign-leaf identity lookup (host sidecar; foreign_slot = no-(treeBase+maxNodes),
             * EXPLICIT and bounds-checked -- not the node index no-treeBase). Only the tag is
             * needed here; the evaluation reads the identity fields again. */
            int fl_tag = 0;
            if(in_foreign && ForeignLeafTag) {
                int fs = no - (treeBase + maxNodes);
                if(fs >= 0 && fs < AllocatedForeignNodes) {fl_tag = ForeignLeafTag[fs];}
            }

            const double node_mass = nop->u.d.mass;
            if(node_mass <= 0) /* nothing in the node */
            {
                no = nop->u.d.sibling;
                continue;
            }
            /* Classify BEFORE anything that could descend.  The wire tag is the authority on
             * whether this node's children were shipped; both the single-particle branch just
             * below and the acceptance predicate further down consult this one answer, so no
             * path can follow nextnode on a node whose children the sender never sent. */
            grav_node_kind_t node_kind = grav_classify_node(in_foreign, fl_tag, nop->u.d.nextnode);
            const int node_is_terminal = grav_node_is_terminal(node_kind);
            const int node_sibling = nop->u.d.sibling;
            for(int w = 0; w < mask_words; w++) {accept_mask[w] = 0; opened_mask[w] = 0;}
            int any_accept = 0, any_open = 0;
#ifdef SINGLE_STAR_DIRECT_GRAVITY
            /* A star target must take no star mass from the tree, since star_direct_gravity_compute()
               supplies every star-star pair exactly. Nodes made entirely of stars therefore have
               nothing left for such a target; it is done with this node before the drift and the
               opening criteria below, whatever the packet's other members do with it. */
            const int node_is_pure_star = ((nop->sink_mass > 0) && (node_mass - nop->sink_mass <= 0));
#else
            const int node_is_pure_star = 0;
#endif
            //if(nop->N_part <= 1)
            if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
            {
                if(node_mass) /* open cell: descend to the particle this node holds */
                {
                    /* Only a local node, or a foreign node shipped WITH its children, has
                     * anything below it here.  A foreign node that reaches this branch with the
                     * multi-particle bit clear was shipped multipole-only (the packer sets that
                     * bit on every leaf it does ship), so its nextnode is the continuation past
                     * the subtree, not a child -- following it would skip the node's mass. */
                    if(node_kind == GRAV_NODE_LOCAL || node_kind == GRAV_NODE_FOREIGN_OPENABLE)
                    {
                        /* descend only for a member that is not done with this node; a star member
                           looking at a pure-star node sits out the descent, and if no one else is
                           taking part the node is passed over as a single walk passes it */
                        int n_descending = 0;
                        for(int m = 0; m < n_members; m++)
                        {
                            int &resume_at = resume_of(m);
                            if(resume_at != GRAV_WALK_MEMBER_IN_PLAY) {continue;}
                            if(node_is_pure_star && open_of(m).ptype == 5) {resume_at = node_sibling;} else {n_descending++;}
                        }
                        if(n_descending == 0) {no = node_sibling; continue;}   /* the members sat out re-join there, at the loop head */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
            }
            if(nop->Ti_current != ti_Current) // add this so that threads arriving here after the the node has been drifted do not have to enter critical at all!
            {
#ifdef _OPENMP
#pragma omp critical(_nodedriftforce_)
#endif
                {
                    if(nop->Ti_current != ti_Current) {force_drift_node(no, ti_Current);}
                }
            }

            for(int m = 0; m < n_members; m++)
            {
                if(resume_of(m) != GRAV_WALK_MEMBER_IN_PLAY) {continue;}
                const struct grav_walk_open_inputs_t &open = open_of(m);
                if(node_is_pure_star && open.ptype == 5) {continue;} /* done with this node (see above) */
                const int ptype = open.ptype; const Vec3<double> &pos = open.pos; const double soft = open.soft, h = open.soft, aold = open.aold;
                double mass = node_mass;
                /* the node's centre of mass and mass as the opening criteria judge them; the
                   evaluation loads both again the same way */
                Vec3<double> dr;
#ifdef SINGLE_STAR_DIRECT_GRAVITY
                /* Remove the sinks from this node for a star target (see the evaluation above for why
                   the subtraction is exact); this sits after force_drift_node so both terms are on the
                   same clock, and mass is reduced before the opening criteria below so they judge the
                   node on the mass actually being used. Pure-star nodes were already skipped above. */
                if((ptype == 5) && (nop->sink_mass > 0))
                {
                    double mass_nosink = mass - nop->sink_mass;
                    dr = (nop->u.d.s * mass - nop->sink_pos * nop->sink_mass) / mass_nosink - pos;
                    mass = mass_nosink;
                }
                else {dr = nop->u.d.s - pos;}
#else
                dr = nop->u.d.s - pos;
#endif
                GRAVITY_NEAREST_XYZ(dr[0],dr[1],dr[2],-1);
                double r2 = dr.norm_sq();
                /* Acceptance geometry via the shared predicate (gravtree_opening.h), the single home
                 * for the node opening decision. The caller owns the wrapped dr/r2 and the
                 * foreign-multipole policy; the predicate is foreign-blind geometry. PM short-range
                 * cull, neighbour sphere-box / softening-open, the angular and relative opening
                 * criteria, and the sink-direct gate all live in the predicate. */
                double cen0 = nop->center[0] - pos[0];
                double cen1 = nop->center[1] - pos[1];
                double cen2 = nop->center[2] - pos[2];
#ifdef PMGRID
                double pred_rcut = open.rcut, pred_rcut2 = open.rcut2;
#else
                double pred_rcut = 0.0, pred_rcut2 = 0.0;
#endif
#if (defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)) && defined(SINGLE_STAR_DIRECT_GRAVITY_RADIUS)
                int pred_n_sink = (int)nop->N_SINK;
#else
                int pred_n_sink = 0;
#endif
                gravtree_open_t pred = gravtree_open_decision_from_distances(
                    r2, cen0, cen1, cen2, soft, h, aold, ptype,
                    nop->len, mass, nop->maxsoft,
                    pred_rcut, pred_rcut2, pred_n_sink, pred_is_first_step);
                /* Foreign LET policy (mirrors the GPU walk exactly).  A terminal foreign node --
                 * a tagged single-particle leaf, or an aggregate the sender shipped multipole-only --
                 * has no children here: its nextnode is the continuation PAST the subtree, not a
                 * child.  So a predicate OPEN on one cannot mean "descend"; for a leaf it means
                 * "accept this already-leaf source with leaf semantics" (restored at evaluation), and for a
                 * truncated aggregate it means the import no longer covers what this walk asks of
                 * it.  Only a node shipped WITH its children takes nextnode. */
                if(pred == GRAV_SKIP_NODE) {continue;}
                if(pred == GRAV_OPEN_NODE && !node_is_terminal) {GRAV_WALK_MASK_SET(opened_mask, m); any_open = 1; continue;}
                /* Import-completeness guard.  Accepting this multipole would drop the sub-node
                 * structure the target resolves, silently and asymmetrically between ranks, so the
                 * walk notes it instead.  Counted rather than printed: this can fire per opened
                 * node per target from inside the threaded walk, and endrun() only REQUESTS a stop
                 * and returns, so printing here would bury the run in interleaved lines.  The count
                 * is reported once and drained after the walk (gravtree.cc). */
                if(pred == GRAV_OPEN_NODE && (node_kind == GRAV_NODE_FOREIGN_TRUNCATED
                                              || node_kind == GRAV_NODE_FOREIGN_UNSHIPPABLE))
                {
                    notes.count++;
                    if(!notes.have_example) {notes.have_example = 1; notes.node = no; notes.id = (unsigned long long) P[members[m].target].ID; notes.ptype = ptype; notes.len = (double) nop->len; notes.mass = (double) nop->u.d.mass;}
                    if(node_kind == GRAV_NODE_FOREIGN_UNSHIPPABLE) {notes.unshippable++;}
                }
                /* ok this member will be using this node */
                GRAV_WALK_MASK_SET(accept_mask, m); any_accept = 1;
            }
            if(any_accept)
            {
                records[n_records].no = no; records[n_records].kind = node_kind; records[n_records].leaf_tag = fl_tag;
                if(!packet_of_one) {for(int w = 0; w < mask_words; w++) {masks[(size_t) n_records * mask_words + w] = accept_mask[w];}}
                n_records++;
            }
            if(any_open)
            {
                /* the packet descends for the members that opened; everyone else still taking part
                   is done with this node and sits out its subtree */
                for(int m = 0; m < n_members; m++)
                {
                    int &resume_at = resume_of(m);
                    if(resume_at != GRAV_WALK_MEMBER_IN_PLAY || GRAV_WALK_MASK_TEST(opened_mask, m)) {continue;}
                    resume_at = node_sibling;
                }
                no = nop->u.d.nextnode;
                continue;
            }
        } /* ok we've completed all the opening criteria -- we will keep this node or particle as-is */

        /* advance for used nodes: note this used to be above, now handled down here so we can use the 'no/nop' structures above */
        if(no < treeSlots) {
            no = Nextnode[no];
        } else {
            no = nop->u.d.sibling;
        }

    } // closes the traversal
#undef GRAV_WALK_MASK_SET
#undef GRAV_WALK_MASK_TEST

    /* the packet completed: its import notes and every member's results are committed now */
    if(TREE_QUERY_PACKET_SIZE == 1) {members[0].resume_at = single_resume_at;}
    commit_notes();
    for(int m = 0; m < n_members; m++)
    {
        const struct grav_walk_member_t &mem = members[m];
        const int target_m = mem.target;
        ninter_out[m] = 0;
        if(mem.resume_at == GRAV_WALK_MEMBER_NEVER) {continue;} /* nothing was computed, nothing is written */
        const grav_pair_acc_t &out = mem.out;
#ifdef SINK_CALC_DISTANCES
        const grav_sink_prox_accum_t &sink_prox = mem.sink_prox;
#endif
#ifdef RT_USE_GRAVTREE
        const int valid_gas_particle_for_rt = mem.valid_gas_particle_for_rt;
#endif
        /* store result at the proper place (local target only; the imported-particle export path is retired) */
        P[target_m].GravAccel = out.acc;
#ifdef RT_USE_TREECOL_FOR_NH
        int k; for(k=0; k < RT_USE_TREECOL_FOR_NH; k++) P[target_m].ColumnDensityBins[k] = mem.treecol_angular_bins[k];
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
        P[target_m].TreeMass = mem.tree_mass;
#endif
#ifdef RT_OTVET
        if(valid_gas_particle_for_rt) {int k; for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[target_m].ET[k] = mem.RT_ET[k];}} else {if(P[target_m].Type==0) {int k; for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[target_m].ET[k] = {};}}}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        if(valid_gas_particle_for_rt) {CellP[target_m].Rad_Flux_UV = mem.incident_flux_uv;}
        if(valid_gas_particle_for_rt) {CellP[target_m].Rad_Flux_EUV = mem.incident_flux_euv;}
#endif
#ifdef CHIMES_STELLAR_FLUXES
        if(valid_gas_particle_for_rt)
        {
            int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {CellP[target_m].Chimes_G0[kc] = mem.chimes_flux_G0[kc]; CellP[target_m].Chimes_fluxPhotIon[kc] = mem.chimes_flux_ion[kc];}
        }
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
        P[target_m].MencInRcrit = mem.m_enc_in_rcrit;
#endif
#ifdef SINK_COMPTON_HEATING
        if(valid_gas_particle_for_rt) {CellP[target_m].Rad_Flux_AGN = mem.incident_flux_agn;}
#endif
#if defined(COSMIC_RAY_SUBGRID_LEBRON)
        if(P[target_m].Type==0) {CellP[target_m].SubGrid_CosmicRayEnergyDensity = mem.SubGrid_CosmicRayEnergyDensity;}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        if(valid_gas_particle_for_rt) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[target_m].Rad_E_gamma[kf] = mem.Rad_E_gamma[kf];}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        if(valid_gas_particle_for_rt) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[target_m].Rad_Flux[kf] = mem.Rad_Flux[kf];}}
#endif
#ifdef EVALPOTENTIAL
        P[target_m].Potential = out.pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        P[target_m].tidal_tensorps = out.tidal_tensorps;
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        P[target_m].tidal_zeta = out.tidal_zeta;
#endif
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        P[target_m].GravJerk = out.jerk;
#endif
#ifdef SINK_CALC_DISTANCES
        P[target_m].Min_Distance_to_Sink = sqrt( sink_prox.Min_Distance_to_Sink2 );
        P[target_m].Min_xyz_to_Sink = sink_prox.Min_xyz_to_Sink;   /* remember, dr = x_SINK - myx */
#ifdef SPECIAL_POINT_MOTION
        {
            P[target_m].vel_of_nearest_special = sink_prox.vel_of_nearest_special;
            P[target_m].acc_of_nearest_special = sink_prox.acc_of_nearest_special;
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
            P[target_m].weight_sum_for_special_point_smoothing = sink_prox.weight_sum_for_special_point_smoothing; /* weighted sum needed */
#endif
        }
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
        P[target_m].is_in_a_binary=0; P[target_m].Min_Sink_OrbitalTime=sink_prox.Min_Sink_OrbitalTime; //orbital time for binary
        if (sink_prox.Min_Sink_OrbitalTime<MAX_REAL_NUMBER)
        {
            P[target_m].is_in_a_binary=1; P[target_m].comp_Mass=sink_prox.comp_Mass; //mass of binary companion
            P[target_m].comp_dx = sink_prox.comp_dx; P[target_m].comp_dv = sink_prox.comp_dv;
        }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
        P[target_m].Min_Sink_Approach_Time = sqrt(sink_prox.Min_Sink_Approach_Time);
        P[target_m].Min_Sink_Freefall_time = sqrt(sqrt(sink_prox.Min_Sink_Freefall_time)/All.G);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        P[target_m].Min_Sink_FeedbackTime = sqrt(sink_prox.Min_Sink_FeedbackTime);
#endif
#endif
#endif // SINK_CALC_DISTANCES
        ninter_out[m] = out.ninter;
    }

    return 1;
}





#ifdef BOX_PERIODIC
/*! This function computes the Ewald correction, and is needed if periodic
 *  boundary conditions together with a pure tree algorithm are used. Note
 *  that the ordinary tree walk does not carry out this correction directly
 *  as it was done in Gadget-1.1. Instead, the tree is walked a second
 *  time. This is actually faster because the "Ewald-Treewalk" can use a
 *  different opening criterion than the normal tree walk. In particular,
 *  the Ewald correction is negligible for particles that are very close,
 *  but it is large for particles that are far away (this is quite
 *  different for the normal direct force). So we can here use a different
 *  opening criterion. Sufficient accuracy is usually obtained if the node
 *  length has dropped to a certain fraction ~< 0.25 of the
 *  BoxLength. However, we may only short-cut the interaction list of the
 *  normal full Ewald tree walk if we are sure that the whole node and all
 *  daughter nodes "lie on the same side" of the periodic boundary,
 *  i.e. that the real tree walk would not find a daughter node or particle
 *  that was mapped to a different nearest neighbour position when the tree
 *  walk would be further refined.
 */
int force_treeevaluate_ewald_correction(int target, int *exportflag, int *exportnodecount, int *exportindex)
{
    struct NODE *nop = 0;
    int signx, signy, signz, nexp, openflag, task, no, cost;
    double mass, r2, u;   /* u: scratch in the periodic-boundary shortcut; interp locals now in the SSOT helper */
    double boxsize, boxhalf, aold, xtmp; xtmp=0;
    Vec3<double> pos, dr;
    Vec3<MyDouble> acc = {};

    boxsize = All.BoxSize;
    boxhalf = 0.5 * All.BoxSize;

    cost = 0;
    pos = P[target].Pos;
    aold = All.ErrTolForceAcc * P[target].OldAcc;

    no = All.TreeNodeIndexBase;        /* root node */

    while(no >= 0)   /* outer loop runs once: the mode-1 imported-NodeList iteration is retired */
    {
        while(no >= 0)
        {
            /* What the imported tree allows this walk to do with the node.  Set once per node, above
             * every branch that could descend it, exactly as the primary walk does: this walk runs on
             * the same installed tree and reaches the same imported nodes, so it owes the same
             * terminal-node contract (see let_data.h).  A particle leaf never reaches the branches
             * that read it, so a local node is the right default. */
            grav_node_kind_t node_kind = GRAV_NODE_LOCAL;
            if(no >= All.TreeParticleSlots && no < All.TreeNodeIndexBase) {/* An index between the particle slots and the node base belongs to neither, so the tree is
             * malformed; stop rather than read a side array or Nodes[] out of bounds. */
                endrun(90001024); no = -1; continue;}
            if(no < All.TreeParticleSlots)    /* single particle */
            {
                /* the index of the node is the index of the particle */
                /* observe the sign */
                if(P[no].Ti_current != All.Ti_Current)
                {
#ifdef _OPENMP
#pragma omp critical(_particledriftewald_)
#endif
                    {
                        if(P[no].Ti_current != All.Ti_Current) {
                            drift_particle(no, All.Ti_Current);
                            gizmo_mark_kernel_radius_dirty_indices(&no, 1);
                        }
                    }
                }

                dr = P[no].Pos - pos;
                mass = P[no].Mass;
            }
            else            /* we have an  internal node */
            {
                if(no >= All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes)    /* pseudo particle (foreign-node range below pseudos) */
                {
                    /* LET-incompleteness DETECTOR (the MPI export round-trip is retired): reaching a
                     * non-empty pseudo means this target's Ewald correction is not covered by the local
                     * LET; record it so Nexport>0 triggers gravity_tree()'s graceful controlled-stop. */
                    if(exportflag[task = DomainTask[no - (All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes)]] != target)
                    {
                        exportflag[task] = target;
                        exportnodecount[task] = NODELISTLENGTH;
                    }

                    if(exportnodecount[task] == NODELISTLENGTH)
                    {
                        int exitFlag = 0;
#ifdef _OPENMP
#pragma omp critical(_nexportewald_)
#endif
                        {
                            if(Nexport >= All.BunchSize)
                            {
                                /* Table full: same failure as recording the entry (see the primary walk),
                                 * and the dropped targets have no retry pass, so stop the run here too. */
                                BufferFullFlag = 1;
                                exitFlag = 1;
                                gizmo_request_controlled_stop(914040, "gravtree: the locally essential tree did not cover these targets' gravity, and there were too many of them to record", __FILE__, __LINE__, __FUNCTION__);
                            }
                            else
                            {
                                nexp = Nexport;
                                Nexport++;
                            }
                        }
                        if(exitFlag) {return -1;} /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */

                        exportnodecount[task] = 0;
                        exportindex[task] = nexp;
                        DataIndexTable[nexp].Task = task;
                        DataIndexTable[nexp].Index = target;
                        DataIndexTable[nexp].IndexGet = nexp;
                    }

                    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] = DomainNodeIndex[no - (All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes)];

                    if(exportnodecount[task] < NODELISTLENGTH) {DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;}
                    no = Nextnode[All.TreeParticleSlots + (no - All.TreeNodeIndexBase - MaxNodes - MaxForeignNodes)];
                    continue;
                }

                nop = &Nodes[no];
                /* The wire tag is the authority on whether this node's children were shipped.  Read
                 * it here, before either descent below, so neither can follow a continuation pointer
                 * into children that were never sent.  Only the tag is needed: the Ewald correction
                 * depends on mass and separation alone, so a foreign leaf carries no identity this
                 * walk has to restore. */
                {
                    int in_foreign = (no >= All.TreeNodeIndexBase + MaxNodes &&
                                      no <  All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes);
                    int fl_tag = 0;
                    if(in_foreign && ForeignLeafTag)
                    {
                        int fs = no - (All.TreeNodeIndexBase + MaxNodes);
                        if(fs >= 0 && fs < AllocatedForeignNodes) {fl_tag = ForeignLeafTag[fs];}
                    }
                    node_kind = grav_classify_node(in_foreign, fl_tag, nop->u.d.nextnode);
                }

                //if(nop->N_part <= 1) /* open cell */
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* Only a local node, or a foreign one shipped WITH its children, has anything
                     * below it.  A terminal foreign node holds a single particle whose multipole is
                     * that particle exactly, so fall through and use it where it stands. */
                    if(node_kind == GRAV_NODE_LOCAL || node_kind == GRAV_NODE_FOREIGN_OPENABLE)
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
                if(nop->Ti_current != All.Ti_Current)
                {
#ifdef _OPENMP
#pragma omp critical(_nodedriftewald_)
#endif
                    {
                        if(nop->Ti_current != All.Ti_Current) {force_drift_node(no, All.Ti_Current);}
                    }
                }

                mass = nop->u.d.mass;
                dr = nop->u.d.s - pos;
            }
            GRAVITY_NEAREST_XYZ(dr[0],dr[1],dr[2],-1);

            if(no < All.TreeParticleSlots)
            {no = Nextnode[no];}
            else            /* we have an internal node. Need to check opening criterion */
            {
                openflag = 0;
                r2 = dr.norm_sq();
                if(r2 <= 0) {r2=MIN_REAL_NUMBER;}
                if(All.ErrTolTheta)    /* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                        openflag = 1;
                    }
                }
#ifndef GRAVITY_HYBRID_OPENING_CRIT
                else        /* check relative opening criterion */
#else
                    if(!(All.Ti_Current == 0 && RestartFlag != 1))
#endif
                    {
                        if(mass * nop->len * nop->len > r2 * r2 * aold)
                        {
                            openflag = 1;
                        }
                        else
                        {
                            if(GRAVITY_NGB_PERIODIC_BOX_LONG_X(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1) < 0.60 * nop->len)
                            {
                                if(GRAVITY_NGB_PERIODIC_BOX_LONG_Y(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1) < 0.60 * nop->len)
                                {
                                    if(GRAVITY_NGB_PERIODIC_BOX_LONG_Z(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1) < 0.60 * nop->len)
                                    {
                                        openflag = 1;
                                    }
                                }
                            }
                        }
                    }

                if(openflag)
                {
                    /* The multipole may still stand in for the whole node, but only if the node lies
                     * wholly on one side of every periodic boundary -- otherwise its daughters map to
                     * different nearest images and the correction has to be resolved below it -- and
                     * only if it is small compared with the box. */
                    int must_refine = 0;
                    for(int kdim = 0; kdim < 3 && !must_refine; kdim++)
                    {
                        u = nop->center[kdim] - pos[kdim];
                        if(u > boxhalf) {u -= boxsize;}
                        if(u < -boxhalf) {u += boxsize;}
                        if(fabs(u) > 0.5 * (boxsize - nop->len)) {must_refine = 1;}
                    }
                    if(!must_refine && nop->len > 0.20 * boxsize) {must_refine = 1;} /* cell is too large */

                    if(must_refine)
                    {
                        /* Only a node shipped with its children can be resolved below.  A foreign
                         * leaf is a single particle, so refining it would change nothing and using
                         * it is exact.  A truncated aggregate is neither: the import no longer
                         * carries the structure this walk resolves, so record it -- the count is
                         * reported once and drained after the walk (gravtree.cc), for the same
                         * reason as the primary walk's, this one being threaded over targets too. */
                        if(!grav_node_is_terminal(node_kind)) {no = nop->u.d.nextnode; continue;}
                        if(node_kind == GRAV_NODE_FOREIGN_TRUNCATED || node_kind == GRAV_NODE_FOREIGN_UNSHIPPABLE)
                        {
                            gravity_note_incomplete_import(no, (unsigned long long) P[target].ID, P[target].Type,
                                                           (double) nop->len, (double) nop->u.d.mass);
                            if(node_kind == GRAV_NODE_FOREIGN_UNSHIPPABLE) {gravity_note_unshippable_import(1);}
                        }
                    }
                }

                no = nop->u.d.sibling;    /* ok, node can be used */
            }
            
            /* compute the Ewald correction force */

            if(dr[0] < 0)
            {
                dr[0] = -dr[0];
                signx = +1;
            }
            else
            {signx = -1;}
            if(dr[1] < 0)
            {
                dr[1] = -dr[1];
                signy = +1;
            }
            else
            {signy = -1;}
            if(dr[2] < 0)
            {
                dr[2] = -dr[2];
                signz = +1;
            }
            else
            {signz = -1;}
            /* trilinear interp of the Ewald force octant tables via the shared SSOT helper
             * (gravtree_ewald.h): the index + 8 weights are built once from |dr| and applied to
             * all three force tables; the odd-force per-component signs stay here. */
            grav_ewald_interp_weights ew = grav_ewald_interp_setup(dr[0], dr[1], dr[2], fac_intp);
            acc[0] += mass * signx * grav_ewald_interp_apply(&fcorrx[0][0][0], ew);
            acc[1] += mass * signy * grav_ewald_interp_apply(&fcorry[0][0][0], ew);
            acc[2] += mass * signz * grav_ewald_interp_apply(&fcorrz[0][0][0], ew);
            cost++;
        }
        
    }

    /* add the result at the proper place (local target only; the imported-particle export path is retired) */
    P[target].GravAccel += acc;

    return cost;
}
#endif // #ifdef BOX_PERIODIC //









#ifdef SUBFIND
/* Local-only Barnes-Hut potential walk for SUBFIND unbinding. Foreign-rank mass is
 * supplied by the locally-installed LET nodes (LET is mandatory and force_treebuild
 * hard-stops on an incomplete LET), so no export/import round-trip is needed: each
 * rank fully computes its own particles' DM_Potential. Reaching a pseudo-particle
 * would mean the LET is incomplete -- a hard correctness failure, since the export
 * fallback is retired -- so we request a controlled stop and abandon the walk
 * (return nonzero; the caller drains at its next collective poll before any
 * potential-derived logic). */
int subfind_force_treeevaluate_potential(int target)
{
    struct NODE *nop = 0;
    MyDouble pot = 0;
    int no = All.TreeNodeIndexBase;    /* root node */
    double r2, dx, dy, dz, mass, r, u, h, h_inv, pos_x, pos_y, pos_z;

    pos_x = P[target].Pos[0];
    pos_y = P[target].Pos[1];
    pos_z = P[target].Pos[2];
    h = ForceSoftening_KernelRadius(target); h_inv = 1.0 / h;

    while(no >= 0)
    {
        /* What the imported tree allows this walk to do with the node; set once per node, above every
         * branch that could descend it (see let_data.h).  Particle leaves never reach the branches
         * that read it, so a local node is the right default. */
        grav_node_kind_t node_kind = GRAV_NODE_LOCAL;
        if(no >= All.TreeParticleSlots && no < All.TreeNodeIndexBase) {/* An index between the particle slots and the node base belongs to neither, so the tree is
             * malformed; stop rather than read a side array or Nodes[] out of bounds. */
            endrun(90001024); no = -1; continue;}
        if(no < All.TreeParticleSlots)    /* single particle: node index is the particle index */
        {
            dx = P[no].Pos[0] - pos_x;
            dy = P[no].Pos[1] - pos_y;
            dz = P[no].Pos[2] - pos_z;
            mass = P[no].Mass;
        }
        else
        {
            if(no >= All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes)    /* pseudo particle */
            {
                /* LET supplies complete foreign coverage in all allowed builds; the
                 * export fallback is retired. A pseudo here means the LET is incomplete. */
                endrun(90000080);    /* graceful stop request; does NOT return -- abandon this walk */
                return 1;
            }

            nop = &Nodes[no];
            mass = nop->u.d.mass;
            /* Same terminal-node contract as the gravity walks: this walk runs on the same installed
             * tree, so the wire tag decides whether a node's children are here to descend into. */
            {
                int in_foreign = (no >= All.TreeNodeIndexBase + MaxNodes &&
                                  no <  All.TreeNodeIndexBase + MaxNodes + MaxForeignNodes);
                int fl_tag = 0;
                if(in_foreign && ForeignLeafTag)
                {
                    int fs = no - (All.TreeNodeIndexBase + MaxNodes);
                    if(fs >= 0 && fs < AllocatedForeignNodes) {fl_tag = ForeignLeafTag[fs];}
                }
                node_kind = grav_classify_node(in_foreign, fl_tag, nop->u.d.nextnode);
            }
            if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
            {
                /* A terminal foreign node holds one particle; its multipole is that particle, so use
                 * it where it stands rather than following a pointer past the subtree. */
                if(mass && (node_kind == GRAV_NODE_LOCAL || node_kind == GRAV_NODE_FOREIGN_OPENABLE))
                {no = nop->u.d.nextnode; continue;}    /* open cell */
            }

            dx = nop->u.d.s[0] - pos_x;
            dy = nop->u.d.s[1] - pos_y;
            dz = nop->u.d.s[2] - pos_z;
        }
        GRAVITY_NEAREST_XYZ(dx,dy,dz,-1);
        r2 = dx * dx + dy * dy + dz * dz;
        if(no < All.TreeParticleSlots)
        {
            no = Nextnode[no];
        }
        else            /* internal node: check Barnes-Hut opening criterion */
        {
            double ErrTolThetaSubfind = All.ErrTolTheta;
            if(nop->len * nop->len > r2 * ErrTolThetaSubfind * ErrTolThetaSubfind)
            {
                /* A terminal foreign node has nothing below it to resolve, so it is consumed where
                 * it stands: a leaf IS the single particle, and an aggregate the sender shipped
                 * multipole-only is used as that multipole.  Unlike the gravity walks, this does not
                 * stop the run.  This potential is an unbinding estimate, not a force the run
                 * integrates, and its opening test is not the one the import was pruned against, so
                 * the aggregate would be refused far more often than the import is genuinely short.
                 * Its multipole is bounded by the accuracy the sender pruned to, which is better
                 * than what following the continuation pointer would do -- leave the node's mass out
                 * of the potential altogether. */
                if(mass && !grav_node_is_terminal(node_kind)) {no = nop->u.d.nextnode; continue;}    /* open cell */
            }
            no = nop->u.d.sibling;    /* node can be used */
        }

        r = sqrt(r2);
        if(r >= h)
        {pot += (-mass / r);}
        else
        {
            u = r * h_inv;
            pot += ( mass * kernel_gravity(u, h_inv, 1, -1) );
        }
    }

    P[target].u.DM_Potential = pot;
    return 0;
}
#endif // SUBFIND //




/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually,
 *  maxnodes approximately equal to 0.7*tree_particle_slots is sufficient to
 *  store the tree for up to tree_particle_slots particles.
 */
void force_treeallocate(int maxnodes, int tree_particle_slots, int foreign_node_slots_exact)
{
    int i;
    size_t bytes;
    double allbytes = 0, allbytes_topleaves = 0;
    double u;
    
    tree_allocated_flag = 1;
    /* Reset the tree-node byte total for a fresh tree. Every caller frees the prior
       tree before rebuilding (domain and subfind free first; the LET/tree retry frees
       inside its loop) or starts from a fresh state (restart), so the counter is
       already 0 here in normal flow; the reset additionally clears any stale total left
       by a prior partial-allocation controlled-stop. Bytes accumulate per successful
       UVM array below and zero again on free. DomainNodeIndex / TopNodeNodeIndex go
       through mymalloc (Base arena) and are not counted here, to avoid double-counting
       the Base total. */
    gizmo_mem_account_set(GIZMO_MEM_TREE_NODES, 0);
    DomainNodeIndex = (int *) mymalloc("DomainNodeIndex", bytes = NTopleaves * sizeof(int));
    allbytes_topleaves += bytes;
    /* Top-leaf-router geometry SSOT (allocated immediately after DomainNodeIndex
     * so the mymalloc LIFO free order in force_treefree is simply the reverse).
     * Sized for ALL topnodes (internal + leaf); populated in
     * force_create_empty_nodes. */
    TopNodeNodeIndex = (int *) mymalloc("TopNodeNodeIndex", bytes = (NTopnodes > 0 ? NTopnodes : 1) * sizeof(int));
    allbytes_topleaves += bytes;
    for(i = 0; i < NTopnodes; i++) TopNodeNodeIndex[i] = -1;  /* sentinel: post-build validation requires all populated */
    /* The tree's storage contract, enforced at the single point every caller passes through:
     *     NumPart  <=  tree_particle_slots  <=  All.TreeNodeIndexBase.
     * The upper leg keeps particle and node indices from ever overlapping (the base is fixed for
     * the run); the lower leg keeps every local particle index addressable in Father[] and
     * Nextnode[], whatever capacity a caller chooses. */
    if(tree_particle_slots > All.TreeNodeIndexBase)
    {
        printf("force_treeallocate: %d particle slots exceed the tree node index base %d, so particle and node indices would overlap.\n",
               tree_particle_slots, All.TreeNodeIndexBase);
        fflush(stdout);
        gizmo_request_controlled_stop(90001023, "force_treeallocate: particle slots exceed the tree node index base (restarting with a much larger PartAllocFactor is the usual cause; restart with the original value, or start a fresh run at the larger one)", __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    if(NumPart > tree_particle_slots)
    {
        printf("force_treeallocate: %d local particles exceed the %d particle slots this tree is being allocated for.\n",
               NumPart, tree_particle_slots);
        fflush(stdout);
        gizmo_request_controlled_stop(90001025, "force_treeallocate: local particle count exceeds the tree's particle slots", __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    All.TreeParticleSlots = tree_particle_slots;
    MaxNodes = maxnodes;
    /* LET: foreign-node headroom in Nodes_base/Extnodes_base/Nextnode.
     * Index map (single source of truth).  Write Base for All.TreeNodeIndexBase:
     *   [0,                              Base)                                      -> particles
     *   [Base,                           Base+MaxNodes)                             -> local tree nodes
     *   [Base+MaxNodes,                  Base+MaxNodes+MaxForeignNodes)             -> foreign tree nodes (LET unpack)
     *   [Base+MaxNodes+MaxForeignNodes,  Base+MaxNodes+MaxForeignNodes+NTopleaves)  -> pseudo-particles
     * SoA slot for any node index `no` (local OR foreign): idx = no - Base (same formula).
     * These are INDICES.  Storage is separate and follows real counts: Nextnode[]/Father[] hold
     * All.TreeParticleSlots particle slots (with Nextnode[] carrying the pseudo segment after them).
     * Nothing is ever sized from Base.
     * The foreign range is where index and storage part company, and the distinction is the whole
     * reason the arena is affordable.  MaxForeignNodes is a CEILING on the index range: it fixes
     * where the pseudo-particles start, so it must be identical on every rank, stable for the whole
     * tree epoch, and reproducible from a restart file.  What it costs to raise is arithmetic plus
     * four bytes per slot in Nextnode[] -- so it is set generously and left alone.  The NODES
     * themselves are the expensive part, and how many of them a rank receives is known only after
     * the import counts have been exchanged, long after this function runs.  So Nodes_base and
     * Extnodes_base are allocated here for the local tree ALONE, and the foreign extent (with the
     * foreign-leaf sidecars and the GPU node mirror) is added by force_tree_grow_foreign_storage
     * once the exact count is in hand -- see AllocatedForeignNodes.  Sizing storage from the ceiling
     * instead would reserve every rank the worst rank's import, for the whole run.
     * MaxForeignNodes starts at MaxNodes + synth_overhead, where synth_overhead accounts
     * for synthesized particle leaves (one entry per particle that is a direct child of an essential
     * multi-particle node).  Synthesis overhead ≤ NumPart_per_rank per received rank; using
     * 2 × (All.MaxPartAssignable / PartAllocFactor) = 2 × TotNumPart / NTask as headroom covers NTask=2
     * worst case (both ranks overlap entirely) while staying modest for large NTask.  The
     * factor is fixed within a run, which is required: this term sets MaxForeignNodes and so
     * the pseudo-particle index base, which restart-serialized node pointers are written
     * against.  That is why the ASSIGNMENT cap appears here and not the storage capacity:
     * storage may be raised to hold an unusually large ghost import, on one rank and not
     * another, so deriving from it would make the index base differ between ranks and move
     * during a run.  The assignment cap is the same number everywhere and holds still for as long
     * as any tree does -- it is worked out again only at a restart, before that run builds a tree
     * of its own -- and it is brought up to date there rather than left at whatever the initial
     * conditions held.  The ratio is the balanced particle count on any run whose cap still tracks
     * that count.  It is only an approximation to it on a run that has grown past what it can
     * index, where the cap holds still instead of following: the term is then whatever the cap and
     * the factor in force imply, which errs HIGH if the factor was also lowered there, and high is
     * the safe direction -- this bounds an index range, not a storage size.  That run is at
     * the edge of its index space, and the runtime foreign floor (RuntimeMinLETForeignNodes), which
     * ratchets up when an import does not fit and travels in the restart file, is what covers it.
     * Do NOT re-base this on All.TotNumPart directly,
     * which splits/spawns/eliminations mutate mid-epoch and would move the pseudo range under
     * pointers already written against it.
     * The restart read does not derive at all: it passes the capacity stored in the file
     * (foreign_node_slots_exact), because inside the read window All.MaxPart is the WRITER's
     * (restored for the serialized Nextnode layout) while PartAllocFactor may already be a
     * reader-edited value, so a derivation there would mix provenance and move the pseudo range
     * out from under pointers that are already written.  Every later tree is built from scratch
     * and serializes nothing, so it derives normally.
     * Numforeignnodes (current count, <= AllocatedForeignNodes <= MaxForeignNodes) is reset on each
     * LET exchange.  The range is empty (MaxForeignNodes = 0) only before the first tree is
     * allocated, or when a restart file carries none. */
    if(foreign_node_slots_exact >= 0)
    {
        /* Restart read: the index ceiling the file's node pointers were written against, not an
         * occupancy -- Numforeignnodes is the current count and is not serialized here, and neither
         * are the foreign nodes themselves, which the first tree build after the restart re-imports. */
        MaxForeignNodes = foreign_node_slots_exact;
    }
    else
    {
        double synth_overhead = 2.0 * (double)All.MaxPartAssignable / (double)All.PartAllocFactor;
        /* A starting guess only.  What a rank actually imports is not predictable from the size of
         * its own tree, so this is deliberately not tuned: the runtime floor below is raised to the
         * demand the first time an import does not fit, and that is what sizes the range from then
         * on.  Being wrong here costs one extra tree rebuild, and the range itself costs four bytes
         * per slot. */
        long long base = (long long) ceil((double)MaxNodes + synth_overhead);
        /* Take the larger of the starting guess and the runtime floor (raised by force_treebuild
         * when an import does not fit the index range). This is the SOLE place MaxForeignNodes is
         * derived. */
        long long want = (base > RuntimeMinLETForeignNodes) ? base : RuntimeMinLETForeignNodes;
        if(want < 0) want = 0;
        if(want > (long long)INT_MAX)
        {
            printf("force_treeallocate: LET foreign-node capacity %lld exceeds INT_MAX "
                   "(MaxNodes=%d, RuntimeMinLETForeignNodes=%lld). Stopping.\n",
                   want, MaxNodes, RuntimeMinLETForeignNodes);
            fflush(stdout);
            endrun(90000082);   /* graceful drain; same family as the foreign-arena UVM OOM below */
            want = (long long)INT_MAX;
        }
        MaxForeignNodes = (int) want;
    }
    if(MaxForeignNodes < 0) {MaxForeignNodes = 0;}
    {
        long long index_space_top = (long long) All.TreeNodeIndexBase + (long long) MaxNodes
                                  + (long long) MaxForeignNodes + (long long) NTopleaves;
        if(index_space_top >= (long long) INT_MAX)
        {
            printf("force_treeallocate: tree index space %lld exceeds the int range (base=%d, MaxNodes=%d, MaxForeignNodes=%d, NTopleaves=%d).\n",
                   index_space_top, All.TreeNodeIndexBase, MaxNodes, MaxForeignNodes, NTopleaves);
            fflush(stdout);
            /* Return before allocating anything, the same shape the UVM-OOM paths below use: the
             * node pointers stay NULL from force_treefree, so every caller either polls the stop
             * or NULL-checks and skips its payload before touching the tree.  The base is fixed
             * for the run, so unlike the LET foreign arena this cannot be retried smaller. */
            gizmo_request_controlled_stop(90001022, "force_treeallocate: tree node index space exceeds the int range (add ranks/nodes, or lower PartAllocFactor so the tree index base is smaller)", __FILE__, __LINE__, __FUNCTION__);
            return;
        }
    }
    Numforeignnodes = 0;
    AllocatedForeignNodes = 0;   /* no foreign storage until the import count is known */
    /* FOF/SUBFIND/twopoint pseudo-particle threshold uses
     * `TreeNodeIndexBase+MaxNodes+MaxForeignNodes` (matches the
     * forcetree.cc/let_pack.cc convention).  Halo finding uses two different trees:
     * subfind_loctree builds a private one with no pseudo-particles and no LET, while
     * the group-subset walks go through force_treebuild, which does run
     * let_run_exchange().  The threshold is an index ceiling, not a count, so it is right in
     * both regimes whatever Numforeignnodes happens to be. */
    long long total_node_slots = (long long) MaxNodes + 1LL;   /* local tree only; the foreign extent is added later */
    /* Nodes_base / Extnodes_base live in SharedSpace (UVM) so GPU
     * kernels can read/write them directly.  Same pattern as Father[] and
     * Nextnode[] below.  Skip mymalloc accounting; kokkos_malloc has its
     * own. */
    bytes = (size_t) total_node_slots * sizeof(struct NODE);
    Nodes_base = (struct NODE *) gpu_tree_alloc_bytes(bytes, "tree_nodes");
    if(!Nodes_base)
    {
        printf("failed to allocate %d tree-nodes (%g MB) in SharedSpace.\n",
               MaxNodes, bytes / (1024.0 * 1024.0));
        /* UVM OOM (per-rank). Controlled-stop request + return BEFORE `Nodes = Nodes_base - TreeNodeIndexBase`
         * so the wild alias is never formed; Nodes_base stays NULL for the caller's check.
         * DomainNodeIndex + tree_allocated_flag are already set (partial allocation) -- safe
         * not because nothing was allocated, but because the caller immediately polls (all-rank
         * sites) or NULL-checks + skips payload (restart turn) before any tree use. */
        gizmo_request_controlled_stop(90000082, "force_treeallocate: tree Nodes UVM/SharedSpace OOM (add ranks/nodes or reduce tree/LET-foreign demand)", __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    gizmo_mem_account_add(GIZMO_MEM_TREE_NODES, (long long) bytes);   /* Nodes_base */
    bytes = (size_t) total_node_slots * sizeof(struct extNODE);
    Extnodes_base = (struct extNODE *) gpu_tree_alloc_bytes(bytes, "tree_extnodes");
    if(!Extnodes_base)
    {
        printf("failed to allocate %d tree-extnodes (%g MB) in SharedSpace.\n",
               MaxNodes, bytes / (1024.0 * 1024.0));
        /* UVM OOM (per-rank). Controlled-stop request + return BEFORE `Extnodes = Extnodes_base - TreeNodeIndexBase`
         * so the wild alias is never formed; Extnodes_base stays NULL for the caller's check. */
        gizmo_request_controlled_stop(90000083, "force_treeallocate: tree Extnodes UVM/SharedSpace OOM (add ranks/nodes or reduce tree/LET-foreign demand)", __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    gizmo_mem_account_add(GIZMO_MEM_TREE_NODES, (long long) bytes);   /* Extnodes_base */
    Nodes = Nodes_base - All.TreeNodeIndexBase;
    Extnodes = Extnodes_base - All.TreeNodeIndexBase;
    /* Nextnode also in SharedSpace; soa->nextnode_aux is aliased to
     * this pointer (no separate buffer, no per-walk memcpy).
     * Pseudo-particles live at Nextnode[All.TreeParticleSlots + pseudo ordinal]
     * after the foreign range; foreign nodes carry their next-sibling pointers in NODE.u.d
     * directly so they don't consume Nextnode[] slots, but we extend the buffer so the
     * pseudo-particle range stays in bounds after the index shift. */
    long long nextnode_slots = (long long) tree_particle_slots + (long long) NTopnodes + (long long) MaxForeignNodes;
    bytes = (size_t) nextnode_slots * sizeof(int);
    Nextnode = (int *) gpu_tree_alloc_bytes(bytes, "tree_nextnode");
    if(!Nextnode)
    {
        printf("Failed to allocate %lld 'Nextnode' slots (%g MB) in SharedSpace\n",
               nextnode_slots, bytes / (1024.0 * 1024.0));
        /* UVM OOM (per-rank). Controlled-stop request + return BEFORE gpu_gravity_tree_alias_nextnode()
         * so a NULL Nextnode is never registered; Nextnode stays NULL for the caller's check. */
        gizmo_request_controlled_stop(90000084, "force_treeallocate: tree Nextnode UVM/SharedSpace OOM (add ranks/nodes or reduce tree/LET-foreign demand)", __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    gpu_gravity_tree_alias_nextnode(Nextnode, (int) nextnode_slots);
    gizmo_mem_account_add(GIZMO_MEM_TREE_NODES, (long long) bytes);   /* Nextnode */
    /* Father[] is UVM (SharedSpace) so the GPU father kernel can
     * write into it directly and host readers (setup_smoothinglengths etc.)
     * page-fault on touch.  No per-tree-build deep_copy needed.  Skip the
     * mymalloc accounting; the matching gpu_father_free lives in force_treefree. */
    bytes = (size_t)tree_particle_slots * sizeof(int);
    Father = gpu_father_alloc(tree_particle_slots);
    if(!Father)
    {
        printf("Failed to allocate %d spaces for 'Father' array (%g MB) in SharedSpace\n",
               tree_particle_slots, bytes / (1024.0 * 1024.0));
        /* UVM OOM (per-rank). Controlled-stop request + return BEFORE any Father[i] write/use;
         * Father stays NULL for the caller's check. */
        gizmo_request_controlled_stop(90000085, "force_treeallocate: tree Father UVM/SharedSpace OOM (add ranks/nodes or reduce tree/LET-foreign demand)", __FILE__, __LINE__, __FUNCTION__);
        return;
    }
    gizmo_mem_account_add(GIZMO_MEM_TREE_NODES, (long long) bytes);   /* Father */
    /* The foreign-leaf identity sidecar follows the foreign node storage, so like it there is
     * nothing to allocate yet: force_tree_grow_foreign_storage creates both together, once the
     * import count is known.  Until then the pointers are NULL, which every reader of them
     * already tests for. */
    ForeignLeafTag = ForeignLeafType = NULL; ForeignLeafZeta = ForeignLeafSoft = NULL;

    /* Don't add to allbytes — kokkos_malloc accounting is separate. */
    if(first_flag == 0)
    {
        first_flag = 1;
        if(ThisTask == 0)
        {
            printf
            ("Allocated %g MByte for tree, and %g Mbyte for top-leaves.  (presently allocated %g MB)\n",
             allbytes / (1024.0 * 1024.0), allbytes_topleaves / (1024.0 * 1024.0),
             AllocatedBytes / (1024.0 * 1024.0));
            /* Gravity-node sizing audit. Lists active payload #ifdefs
             * that inflate NODE/extNODE; informs the compact-node variants. */
            printf("Gravity tree node sizes: sizeof(NODE)=%zu B, sizeof(extNODE)=%zu B; "
                   "MyGravFloat=%zu B (mixed-precision gravity %s). Active payload flags:",
                   sizeof(struct NODE), sizeof(struct extNODE), sizeof(MyGravFloat),
#ifdef GIZMO_MIXED_PRECISION_GRAVITY
                   "ON"
#else
                   "OFF"
#endif
                   );
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
            printf(" GRAVTREE_CALCULATE_GAS_MASS_IN_NODE");
#endif
#ifdef RT_USE_GRAVTREE
            printf(" RT_USE_GRAVTREE");
#endif
#ifdef CHIMES_STELLAR_FLUXES
            printf(" CHIMES_STELLAR_FLUXES");
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            printf(" RT_SEPARATELY_TRACK_LUMPOS");
#endif
#ifdef SINK_PHOTONMOMENTUM
            printf(" SINK_PHOTONMOMENTUM");
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
            printf(" COSMIC_RAY_SUBGRID_LEBRON");
#endif
#ifdef SINK_CALC_DISTANCES
            printf(" SINK_CALC_DISTANCES");
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
            printf(" ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION");
#endif
#ifdef DM_SCALARFIELD_SCREENING
            printf(" DM_SCALARFIELD_SCREENING");
#endif
            printf("\n");
        }
        for(i = 0; i < NTAB; i++)
        {
            u = 3.0 / NTAB * (i + 0.5);
            shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
            shortrange_table_potential[i] = erfc(u);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
            shortrange_table_tidal[i] = 4.0 * u * u * u / sqrt(M_PI) * exp(-u * u);
#endif
        }
    }
}


/*! Give the foreign-node range the storage it needs, now that the LET exchange has counted
 *  exactly how many nodes this rank is about to receive.
 *
 *  force_treeallocate cannot do this: the count depends on the tree that has yet to be built,
 *  and on what every other rank decides to send.  Sizing from the index ceiling instead would
 *  hand every rank the worst rank's import for the whole run, which is what this exists to
 *  stop -- the demand is strongly uneven between ranks, and which rank is the busy one changes
 *  from one build to the next, so no fixed or remembered figure fits it.  Waiting costs one
 *  copy of the local tree, which is a small fraction of the arena being sized.
 *
 *  Called once per tree build, before the first foreign node is installed, on every rank
 *  (a rank importing nothing passes 0 and does nothing).  Growth is one-way within a build:
 *  the count is exact, so there is never a second ask.
 *
 *  Each array is replaced only once its successor is allocated and filled, so a failure here
 *  leaves the tree exactly as it was, with a status for the caller to reduce and report.
 *  Returns 0 on success, nonzero on failure. */
int force_tree_grow_foreign_storage(long long foreign_needed)
{
    if(foreign_needed <= (long long) AllocatedForeignNodes) {return 0;}   /* nothing to add */
    if(foreign_needed > (long long) MaxForeignNodes)
    {
        /* The caller checks this against the ceiling first and retries the whole build with a
         * larger one, so reaching here means the two disagree.  Refuse rather than allocate
         * storage for indices the pseudo-particle range already occupies. */
        printf("force_tree_grow_foreign_storage: asked for %lld foreign nodes on rank %d, above the index ceiling %d.\n",
               foreign_needed, ThisTask, MaxForeignNodes);
        fflush(stdout);
        return 1;
    }
    const int n_foreign = (int) foreign_needed;
    const long long local_slots = (long long) MaxNodes + 1LL;   /* the +1 sentinel travels with the local tree */
    const long long new_slots   = local_slots + (long long) n_foreign;
    const int old_foreign = AllocatedForeignNodes;

    /* Every allocation is made, and the built tree carried across, BEFORE anything the rest of
     * the code can see changes.  So a shortfall anywhere below leaves the tree exactly as it
     * was, still complete and still usable, and the only consequence is the status returned. */
    struct NODE    *new_nodes    = NULL;
    struct extNODE *new_extnodes = NULL;
    int     *new_tag  = NULL;
    int     *new_type = NULL;
    MyFloat *new_zeta = NULL;
    MyFloat *new_soft = NULL;

    /* Carry across everything the old arrays held -- the local tree, and any foreign nodes a
     * previous call already made room for.  Today the exchange sizes once per build so the
     * foreign part is empty, but copying it costs nothing and makes this correct on its own
     * terms rather than by the caller's discipline. */
    const long long keep_slots = local_slots + (long long) old_foreign;
    new_nodes = (struct NODE *) gpu_tree_alloc_bytes((size_t) new_slots * sizeof(struct NODE), "tree_nodes");
    if(new_nodes) {memcpy(new_nodes, Nodes_base, (size_t) keep_slots * sizeof(struct NODE));}
    new_extnodes = (struct extNODE *) gpu_tree_alloc_bytes((size_t) new_slots * sizeof(struct extNODE), "tree_extnodes");
    if(new_extnodes) {memcpy(new_extnodes, Extnodes_base, (size_t) keep_slots * sizeof(struct extNODE));}

    /* Foreign-leaf identity sidecar: one entry per foreign slot, describing THIS import, so
     * there is nothing to carry over -- allocate fresh and zero, which leaves every slot at the
     * leaf_tag=0 (plain node) default until the install writes it. */
    new_tag  = (int *)     gpu_tree_alloc_bytes((size_t) n_foreign * sizeof(int),     "tree_foreign_tag");
    new_type = (int *)     gpu_tree_alloc_bytes((size_t) n_foreign * sizeof(int),     "tree_foreign_type");
    new_zeta = (MyFloat *) gpu_tree_alloc_bytes((size_t) n_foreign * sizeof(MyFloat), "tree_foreign_zeta");
    new_soft = (MyFloat *) gpu_tree_alloc_bytes((size_t) n_foreign * sizeof(MyFloat), "tree_foreign_soft");

    /* The mirror sizes its own copy of the sidecar from AllocatedForeignNodes, so that has to
     * read as the new count while it allocates -- and go back if it cannot. */
    int mirror_grown = 0;
    if(new_nodes && new_extnodes && new_tag && new_type && new_zeta && new_soft)
    {
        AllocatedForeignNodes = n_foreign;
        mirror_grown = gpu_gravity_tree_grow_foreign((int) new_slots);
        if(!mirror_grown) {AllocatedForeignNodes = old_foreign;}
        /* The dirty set is bounded by the mirror it repairs; grow it with the mirror, or the
         * claims the walk makes on the new foreign slots fall outside it. */
        if(mirror_grown) {gpu_node_dirty_grow_to((int) new_slots);}
    }
    if(!mirror_grown)
    {
        printf("force_tree_grow_foreign_storage: rank %d could not find room for %d foreign nodes. "
               "Their node arrays alone are %g MB; the foreign-leaf sidecars and this rank's slice of "
               "the GPU node mirror are on top of that, and the mirror is usually the larger of the two. "
               "The ghost import this step does not fit in node memory; feasible: fewer ranks/node, or "
               "lower resolution.\n",
               ThisTask, n_foreign,
               (double)((size_t) n_foreign * (sizeof(struct NODE) + sizeof(struct extNODE))) / (1024.0 * 1024.0));
        fflush(stdout);
        if(new_soft)     {gpu_tree_free_bytes(new_soft);}
        if(new_zeta)     {gpu_tree_free_bytes(new_zeta);}
        if(new_type)     {gpu_tree_free_bytes(new_type);}
        if(new_tag)      {gpu_tree_free_bytes(new_tag);}
        if(new_extnodes) {gpu_tree_free_bytes(new_extnodes);}
        if(new_nodes)    {gpu_tree_free_bytes(new_nodes);}
        return 1;
    }
    memset(new_tag,  0, (size_t) n_foreign * sizeof(int));
    memset(new_type, 0, (size_t) n_foreign * sizeof(int));
    memset(new_zeta, 0, (size_t) n_foreign * sizeof(MyFloat));
    memset(new_soft, 0, (size_t) n_foreign * sizeof(MyFloat));

    /* Everything is in hand: swap, then release what was superseded. */
    gpu_tree_free_bytes(Nodes_base);
    gpu_tree_free_bytes(Extnodes_base);
    if(ForeignLeafTag)  {gpu_tree_free_bytes(ForeignLeafTag);}
    if(ForeignLeafType) {gpu_tree_free_bytes(ForeignLeafType);}
    if(ForeignLeafZeta) {gpu_tree_free_bytes(ForeignLeafZeta);}
    if(ForeignLeafSoft) {gpu_tree_free_bytes(ForeignLeafSoft);}
    Nodes_base    = new_nodes;
    Extnodes_base = new_extnodes;
    Nodes    = Nodes_base    - All.TreeNodeIndexBase;
    Extnodes = Extnodes_base - All.TreeNodeIndexBase;
    ForeignLeafTag = new_tag; ForeignLeafType = new_type;
    ForeignLeafZeta = new_zeta; ForeignLeafSoft = new_soft;

    gizmo_mem_account_add(GIZMO_MEM_TREE_NODES,
                          (long long)((long long)(n_foreign - old_foreign)
                                      * (long long)(sizeof(struct NODE) + sizeof(struct extNODE)
                                                    + 2 * sizeof(int) + 2 * sizeof(MyFloat))));
    return 0;
}


/*! Memory-ledger provider (diagnostic): the CAPACITY / REQUESTED byte breakdown of
 *  the tree Kokkos allocations -- local node arrays vs foreign-LET node arrays (incl.
 *  foreign-leaf sidecars) vs Father/Nextnode aux -- plus foreign capacity, since-start
 *  used high-water, and the adaptive floor.
 *  The foreign figure is the storage that EXISTS (AllocatedForeignNodes), which is what the
 *  node is actually paying and what the ledger's totals have to add up to.  It is reported
 *  beside the index ceiling MaxForeignNodes, which costs only the Nextnode ints counted in the
 *  aux term: without both, a reader cannot tell a rank with a small import from a rank whose
 *  ceiling happens to be small.  The local term is DEMAND from MaxNodes, since
 *  `tree_allocated_flag` is set at the TOP of force_treeallocate, so on a partial-allocation
 *  controlled-stop this reports the capacity that could not be met -- exactly what is wanted
 *  at a stop.  Returns zeros when no tree is allocated. Called at ledger print time (host).
 *  Byte model mirrors force_treeallocate plus force_tree_grow_foreign_storage:
 *  Nodes_base/Extnodes_base span MaxNodes + AllocatedForeignNodes + 1 (the +1 sentinel folded
 *  into local); DomainNodeIndex/TopNodeNodeIndex are mymalloc (Base family), excluded here. */
void gizmo_tree_mem_breakdown(double *local_mb, double *foreign_cap_mb, double *aux_mb,
                              long long *foreign_alloc_nodes, long long *foreign_used_hw_nodes,
                              long long *foreign_floor_nodes, long long *foreign_ceiling_nodes)
{
    *local_mb = *foreign_cap_mb = *aux_mb = 0.0;
    *foreign_alloc_nodes = *foreign_used_hw_nodes = *foreign_floor_nodes = *foreign_ceiling_nodes = 0;
    if(!tree_allocated_flag) {return;}
    const double MB = 1024.0 * 1024.0;
    long long node_pair = (long long) sizeof(struct NODE) + (long long) sizeof(struct extNODE);
    *local_mb = (double)(((long long) MaxNodes + 1LL) * node_pair) / MB;  /* +1 sentinel slot */
    long long sidecar = 2LL * (long long) sizeof(int) + 2LL * (long long) sizeof(MyFloat);  /* Tag+Type+Zeta+Soft */
    *foreign_cap_mb = (double)((long long) AllocatedForeignNodes * (node_pair + sidecar)) / MB;
    long long nextnode_slots = (long long) All.TreeParticleSlots + (long long) NTopnodes + (long long) MaxForeignNodes;
    *aux_mb = (double)(((long long) All.TreeParticleSlots + nextnode_slots) * (long long) sizeof(int)) / MB;  /* Father + Nextnode */
    *foreign_alloc_nodes   = AllocatedForeignNodes;
    *foreign_used_hw_nodes = Numforeignnodes_highwater;
    *foreign_floor_nodes   = RuntimeMinLETForeignNodes;
    *foreign_ceiling_nodes = MaxForeignNodes;
}

/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
    if(tree_allocated_flag)
    {
        /* SharedSpace (UVM) frees for GPU-addressable tree
         * storage.  Order is reverse-of-alloc (LIFO discipline preserved for
         * the residual mymalloc'd DomainNodeIndex). */
        if(Father)        {gpu_father_free(Father); Father = NULL;}
        gpu_gravity_tree_alias_nextnode(NULL, 0);  /* clear SoA alias before free */
        /* Freeing the tree invalidates every GPU representation derived from it: the
         * node mirror the device walk reads, the drift and moment-refresh state keyed
         * to that mirror, and the records saying the node geometry is current. Nothing
         * may still claim currency for nodes that are going away, so both records are
         * retired here rather than left for the next build's generation bump to mask;
         * release() reaches them through free_arrays_, and the explicit call keeps that
         * guarantee at this site even if the release path is restructured. The mirror
         * and its pools are therefore reallocated once per tree epoch rather than held
         * for the run. */
        gpu_gravity_tree_invalidate_currency();
        gpu_gravity_tree_release();
        if(Nextnode)      {gpu_tree_free_bytes(Nextnode);      Nextnode      = NULL;}
        if(Extnodes_base) {gpu_tree_free_bytes(Extnodes_base); Extnodes_base = NULL;}
        if(Nodes_base)    {gpu_tree_free_bytes(Nodes_base);    Nodes_base    = NULL;}
        /* Free the foreign-leaf sidecar (SharedSpace, same allocator as the foreign-node arena). */
        if(ForeignLeafTag)  {gpu_tree_free_bytes(ForeignLeafTag);  ForeignLeafTag  = NULL;}
        if(ForeignLeafType) {gpu_tree_free_bytes(ForeignLeafType); ForeignLeafType = NULL;}
        if(ForeignLeafZeta) {gpu_tree_free_bytes(ForeignLeafZeta); ForeignLeafZeta = NULL;}
        if(ForeignLeafSoft) {gpu_tree_free_bytes(ForeignLeafSoft); ForeignLeafSoft = NULL;}
        myfree(TopNodeNodeIndex);   /* LIFO: allocated right after DomainNodeIndex, so freed right before it */
        myfree(DomainNodeIndex);
        gizmo_mem_account_set(GIZMO_MEM_TREE_NODES, 0);   /* whole-family teardown */
        AllocatedForeignNodes = 0;   /* the foreign storage went with the arrays above */
        tree_allocated_flag = 0;
        /* Father[] went with the tree, so there is no attachment record left to consult. */
        g_force_global_topology_valid = 0;
    }
}





/*! This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 */
void dump_particles(void)
{
    FILE *fd;
    char buffer[DEFAULT_PATH_BUFFERSIZE_TOUSE];
    int i;
    
    snprintf(buffer, DEFAULT_PATH_BUFFERSIZE_TOUSE, "particles%d.dat", ThisTask);
    fd = fopen(buffer, "w");
    my_fwrite(&NumPart, 1, sizeof(int), fd);
    for(i = 0; i < NumPart; i++)
        my_fwrite(&P[i].Pos[0], 3, sizeof(MyFloat), fd);
    for(i = 0; i < NumPart; i++)
        my_fwrite(&P[i].Vel[0], 3, sizeof(MyFloat), fd);
    for(i = 0; i < NumPart; i++)
        my_fwrite(&P[i].ID, 1, sizeof(int), fd);
    fclose(fd);
}



#ifdef BOX_PERIODIC

/*! This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point mass located
 *  at the origin. These corrections are obtained by Ewald summation. (See
 *  Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) The correction fields
 *  are used to obtain the full periodic force if periodic boundaries
 *  combined with the pure tree algorithm are used. For the TreePM
 *  algorithm, the Ewald correction is not used.
 *
 *  The correction fields are stored on disk once they are computed. If a
 *  corresponding file is found, they are loaded from disk to speed up the
 *  initialization.  The Ewald summation is done in parallel, i.e. the
 *  processors share the work to compute the tables if needed.
 */
void ewald_init(void)
{
#ifndef SELFGRAVITY_OFF
    int i, j, k, beg, len, size, n, task, count;
    double x[3], force[3];
    char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
    FILE *fd;
    
    if(ThisTask == 0) {printf("Initializing Ewald correction...\n");}
    
    snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "ewald_spc_table_%d_dbl.dat", EN);
    if((fd = fopen(buf, "r")))
    {
        my_fread(&fcorrx[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        my_fread(&fcorry[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        my_fread(&fcorrz[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        my_fread(&potcorr[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        fclose(fd);
    }
    else
    {
        if(ThisTask == 0) {printf("\nNo Ewald tables in file `%s' found.\nRecomputing them...\n", buf);}
        
        /* ok, let's recompute things. Actually, we do that in parallel. */
        
        size = (EN + 1) * (EN + 1) * (EN + 1) / NTask;
        beg = ThisTask * size;
        len = size;
        if(ThisTask == (NTask - 1))
            len = (EN + 1) * (EN + 1) * (EN + 1) - beg;
        for(i = 0, count = 0; i <= EN; i++)
            for(j = 0; j <= EN; j++)
                for(k = 0; k <= EN; k++)
                {
                    n = (i * (EN + 1) + j) * (EN + 1) + k;
                    if(n >= beg && n < (beg + len))
                    {
                        if((count % (len / 20)) == 0) {PRINT_STATUS("%4.1f percent done", count / (len / 100.0));}
                        x[0] = 0.5 * ((double) i) / EN;
                        x[1] = 0.5 * ((double) j) / EN;
                        x[2] = 0.5 * ((double) k) / EN;
                        ewald_force(i, j, k, x, force);
                        fcorrx[i][j][k] = force[0];
                        fcorry[i][j][k] = force[1];
                        fcorrz[i][j][k] = force[2];
                        if(i + j + k == 0)
                            potcorr[i][j][k] = 2.8372975;
                        else
                            potcorr[i][j][k] = ewald_psi(x);
                        count++;
                    }
                }
        
        for(task = 0; task < NTask; task++)
        {
            beg = task * size;
            len = size;
            if(task == (NTask - 1))
                len = (EN + 1) * (EN + 1) * (EN + 1) - beg;
            MPI_Bcast(&fcorrx[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
            MPI_Bcast(&fcorry[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
            MPI_Bcast(&fcorrz[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
            MPI_Bcast(&potcorr[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
        }
        
        if(ThisTask == 0)
        {
            printf("\nwriting Ewald tables to file `%s'\n", buf);
            if((fd = fopen(buf, "w")))
            {
                my_fwrite(&fcorrx[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                my_fwrite(&fcorry[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                my_fwrite(&fcorrz[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                my_fwrite(&potcorr[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                fclose(fd);
            }
        }
    }
    
    fac_intp = 2 * EN / All.BoxSize;
    for(i = 0; i <= EN; i++)
        for(j = 0; j <= EN; j++)
            for(k = 0; k <= EN; k++)
            {
                potcorr[i][j][k] /= All.BoxSize;
                fcorrx[i][j][k] /= All.BoxSize * All.BoxSize;
                fcorry[i][j][k] /= All.BoxSize * All.BoxSize;
                fcorrz[i][j][k] /= All.BoxSize * All.BoxSize;
            }
    
    if(ThisTask == 0) {printf(" ..initialization of periodic boundaries finished.\n");}
#endif // #ifndef SELFGRAVITY_OFF
}


/*! This function looks up the correction potential due to the infinite
 *  number of periodic particle/node images. We here use tri-linear
 *  interpolation to get it from the precomputed table, which contains
 *  one octant around the target particle at the origin. The other
 *  octants are obtained from it by exploiting symmetry properties.
 */
double ewald_pot_corr(double dx, double dy, double dz)
{
    /* trilinear interp of the potential octant table via the shared SSOT helper (gravtree_ewald.h) */
    grav_ewald_interp_weights w = grav_ewald_interp_setup(dx, dy, dz, fac_intp);
    return grav_ewald_interp_apply(&potcorr[0][0][0], w);
}



/*! This function computes the potential correction term by means of Ewald
 *  summation.
 */
double ewald_psi(double x[3])
{
    double alpha, psi;
    double r, sum1, sum2, hdotx;
    double dx[3];
    int i, n[3], h[3], h2;
    
    alpha = 2.0;
    for(n[0] = -4, sum1 = 0; n[0] <= 4; n[0]++)
        for(n[1] = -4; n[1] <= 4; n[1]++)
            for(n[2] = -4; n[2] <= 4; n[2]++)
            {
                for(i = 0; i < 3; i++)
                    dx[i] = x[i] - n[i];
                r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
                sum1 += erfc(alpha * r) / r;
            }
    
    for(h[0] = -4, sum2 = 0; h[0] <= 4; h[0]++)
        for(h[1] = -4; h[1] <= 4; h[1]++)
            for(h[2] = -4; h[2] <= 4; h[2]++)
            {
                hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
                h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
                if(h2 > 0)
                    sum2 += 1 / (M_PI * h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * cos(2 * M_PI * hdotx);
            }
    
    r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    psi = M_PI / (alpha * alpha) - sum1 - sum2 + 1 / r;
    return psi;
}


/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 */
void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
    double alpha, r2;
    double r, val, hdotx, dx[3];
    int i, h[3], n[3], h2;
    
    alpha = 2.0;
    for(i = 0; i < 3; i++)
        force[i] = 0;
    if(iii == 0 && jjj == 0 && kkk == 0)
        return;
    r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    for(i = 0; i < 3; i++)
        force[i] += x[i] / (r2 * sqrt(r2));
    for(n[0] = -4; n[0] <= 4; n[0]++)
        for(n[1] = -4; n[1] <= 4; n[1]++)
            for(n[2] = -4; n[2] <= 4; n[2]++)
            {
                for(i = 0; i < 3; i++)
                    dx[i] = x[i] - n[i];
                r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
                val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);
                for(i = 0; i < 3; i++)
                    force[i] -= dx[i] / (r * r * r) * val;
            }
    
    for(h[0] = -4; h[0] <= 4; h[0]++)
        for(h[1] = -4; h[1] <= 4; h[1]++)
            for(h[2] = -4; h[2] <= 4; h[2]++)
            {
                hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
                h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
                if(h2 > 0)
                {
                    val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);
                    for(i = 0; i < 3; i++)
                        force[i] -= h[i] * val;
                }
            }
}
#endif // #ifdef BOX_PERIODIC //


/*! Refresh tree node moments without rebuilding the tree structure. Uses bottom-up accumulation
 *  via Father[] pointers instead of u.suns[] (which are destroyed after the initial tree build
 *  since they share a union with u.d). Nodes are processed from high to low index, which gives
 *  bottom-up order since children are always allocated with higher indices than parents.
 *  Use this when particle properties (mass, type, luminosity) have changed but particles haven't
 *  moved, e.g. after star formation or sink SN events. */
void force_refresh_node_moments(void)
{
    int i, k, no;
    PRINT_STATUS("Refreshing tree node moments (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));

    /* GPU moment-refresh kernel computes local-tree node
     * moments + writes back to AoS. After this returns, Nodes[] /
     * Extnodes[] are in the same state CPU steps 1-4 below would
     * produce, so the CPU pseudo-particle path can run unchanged. */
    {
        /* Reset GravCost/Ti_current/Flag/Ti_lastkicked/dp/dp_dm/dp_stellarlum
         * fields that the GPU kernel does not own. These mirror the
         * non-moment lines in CPU step 1 (forcetree.cc:3837..3848).
         *
         * A node that is behind the current time is drifted to it FIRST, by the one
         * routine that owns the lazy node drift: its length widens for the interval
         * since it was last drifted, its pending kick is applied.  Stamping it current
         * without that would erase the only record of that interval, and every later
         * widening would start from a length that no longer encloses the particles that
         * moved in it.  A build is preceded by a full drift, so there this is a no-op;
         * a refresh runs mid-step, where inactive particles and their nodes are behind. */
        for(no = All.TreeNodeIndexBase; no < All.TreeNodeIndexBase + Numnodestree; no++) {
            if(Nodes[no].Ti_current != All.Ti_Current) {force_drift_node(no, All.Ti_Current);}
            Nodes[no].GravCost = 0;
            Nodes[no].Ti_current = All.Ti_Current;
            Extnodes[no].dp = {};
            Extnodes[no].Ti_lastkicked = All.Ti_Current;
            Extnodes[no].Flag = GlobFlag;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            Extnodes[no].rt_source_lum_dp = {};
#endif
#ifdef DM_SCALARFIELD_SCREENING
            Extnodes[no].dp_dm = {};
#endif
#ifdef SINK_NODE_MOTION_TRACKED
            Extnodes[no].sink_dp = {};   /* sink_pos/sink_vel are set fresh by the moment pass */
#endif
        }
        /* Rank-local GPU refresh steps: on failure set a soft bad-stop and
         * fall through force_exchange_pseudodata (matched, topology-driven);
         * the gravtree:after_refresh_moments poll drains before the walk. */
        if(gpu_moment_refresh(-1) != 0)          {endrun(90000086);}
        /* Mode B: re-seed per-type bands; gpu_moment_refresh wrote scalar
         * hmax to AoS but not per-type. Without this, hmax_per_type[] are
         * left at zero by the GPU bypass and Mode B's SYMMETRIC walker
         * over-prunes (observed as lost neighbours). */
        force_refresh_hmax_per_type_host(Numnodestree);
        if(gpu_force_flag_localnodes() != 0)     {endrun(90000087);}
        int pseudo_status = force_exchange_pseudodata();
        /* skip dependent pseudo-update on an unmatched complete (soft bad-stop set);
         * drains at gravtree:after_refresh_moments. */
        if(!pseudo_status) {
            if(gpu_scatter_pseudo_to_soa() != 0)     {endrun(90000088);}
            if(gpu_topnode_moment_resum() != 0)      {endrun(90000089);}
        }
        PRINT_STATUS(" ..tree node moments refreshed (GPU).");
        return;
    }
}
