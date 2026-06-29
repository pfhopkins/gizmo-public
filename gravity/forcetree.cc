#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#ifdef SUBFIND
#include "../structure/subfind/subfind.h"
#endif

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


/* function to return the value of the force softening for a given cell/particle, depending on the physics and numerical options */
double ForceSoftening_KernelRadius(int p)
{
#ifdef GALSF_MERGER_STARCLUSTER_PARTICLES
    if(P[p].Type == 4) {return P[p].StarParticleEffectiveSize;} // this variable is defined in force softening terms
    //if(P[p].Type == 4) {return All.ForceSoftening[4] * pow(P[p].Mass / (0.5*(All.MaxMassForParticleSplit/3.01+All.MinMassForParticleMerger/0.49)),0.333);} // alternative 'adaptive' version for constant-resolution runs
    //if(P[p].Type == 4) {return All.ForceSoftening[4] * pow(P[p].Mass*UNIT_MASS_IN_SOLAR / (GALSF_MERGER_STARCLUSTER_PARTICLES),0.333);}
#endif
    
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
    if((1 << P[p].Type) & (ADAPTIVE_GRAVSOFT_FORALL)) {return P[p].AGS_KernelRadius;}
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(SELFGRAVITY_OFF) /* softening scale still appears in timestep criterion for problems without self-gravity, so set it adaptively */
#ifdef ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT
    if(P[p].Type == 0) {return DMIN(P[p].KernelRadius, ADAPTIVE_GRAVSOFT_MAX_SOFT_HARD_LIMIT/All.cf_atime);}
#else
    if(P[p].Type == 0) {return P[p].KernelRadius;}
#endif
#endif

#if defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
    if(P[p].Type == 4) {return All.ForceSoftening[P[p].Type] * DMIN(100., DMAX(1., pow(P[p].Mass*UNIT_MASS_IN_SOLAR/100. , 0.33)));}
#endif

#if defined(ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION) /* still playing with criterion below, highly experimental for now */
    if((1 << P[p].Type) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)) {if((P[p].tidal_tensor_mag_prev>0) && (All.Time>All.TimeBegin)) {return DMIN(1.e2*All.ForceSoftening[P[p].Type] , DMAX(All.ForceSoftening[P[p].Type] , All.ForceSoftening[P[p].Type] + 1.25 * pow( (All.DesNumNgb * All.G * P[p].Mass / P[p].tidal_tensor_mag_prev) , 1./3. )));} else {return 100.*All.ForceSoftening[P[p].Type];}}
#endif

    return All.ForceSoftening[P[p].Type]; // this is the default if nothing was active above
}



/*! auxiliary variable used to set-up non-recursive walk */
static int last;

/* some modules compute neighbor fluxes explicitly within the force-tree: in these cases, we need to
 take extra care about opening leaves to ensure possible neighbors are not missed, so defined a flag below for it */
#if (defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(SINGLE_STAR_SINK_DYNAMICS) || defined(GRAVITY_ACCURATE_FEWBODY_INTEGRATION) || defined(ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION))
#define NEIGHBORS_MUST_BE_COMPUTED_EXPLICITLY_IN_FORCETREE
#endif

/*! length of look-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000
/*! variables for short-range lookup table */
static float shortrange_table[NTAB], shortrange_table_potential[NTAB];
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
static float shortrange_table_tidal[NTAB];
#endif
/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;
static int tree_allocated_flag = 0;

#ifdef BOX_PERIODIC
/*! Size of 3D look-up table for Ewald correction force */
#define EN  64
/*! 3D look-up table for Ewald correction to force and potential. Only one octant is stored, the rest constructed by using the symmetry of the problem */
static MyFloat fcorrx[EN + 1][EN + 1][EN + 1];
static MyFloat fcorry[EN + 1][EN + 1][EN + 1];
static MyFloat fcorrz[EN + 1][EN + 1][EN + 1];
static MyFloat potcorr[EN + 1][EN + 1][EN + 1];
static double fac_intp;
#endif


#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) /* need to do box-wrapping, just refer to our standard box-wrapping macros */
#define GRAVITY_NEAREST_XYZ(x,y,z,sign) NEAREST_XYZ(x,y,z,sign)
#define GRAVITY_NGB_PERIODIC_BOX_LONG_X(x,y,z,sign) NGB_PERIODIC_BOX_LONG_X(x,y,z,sign)
#define GRAVITY_NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign)
#define GRAVITY_NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign) NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign)
#else /* either the box is not periodic, OR gravity is not, in either case no box-wrapping is needed */
#define GRAVITY_NEAREST_XYZ(x,y,z,sign) /* this is an empty macro: nothing will happen to the variables input here */
#define GRAVITY_NGB_PERIODIC_BOX_LONG_X(x,y,z,sign) (fabs(x)) /* just return absolute values */
#define GRAVITY_NGB_PERIODIC_BOX_LONG_Y(x,y,z,sign) (fabs(y))
#define GRAVITY_NGB_PERIODIC_BOX_LONG_Z(x,y,z,sign) (fabs(z))
#endif

/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
int force_treebuild(int npart, struct unbind_data *mp)
{
    int flag, treebuild_retries = 0;
    do
    {
        Numnodestree = force_treebuild_single(npart, mp);
        /* -1 from the build is a genuine node-arena capacity failure (both return sites are
           out-of-nodes); values -2 and -3 are RESERVED for a future needs-ownership-restore
           signal and a hard (non-capacity) failure respectively -- cf. gizmo-cpp 98949019 for
           the taxonomy; the MIN reduce below already gives a more-negative code precedence. */
        MPI_Allreduce(&Numnodestree, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(flag == -1)
        {
            /* Bound and abort decided HERE, from the reduced flag and the shared TreeAllocFactor,
               so every rank takes the same branch. Growing past 5.0 (or ~18 growths from the 0.45
               default) means something other than capacity is wrong -- a particle pileup or
               corrupted positions -- and more nodes cannot fix it. */
            if(All.TreeAllocFactor > 5.0 || ++treebuild_retries > 24)
            {
                printf("task %d: force_treebuild: node arena still insufficient at TreeAllocFactor=%g after %d growths -- not a capacity problem. Dumping particles and aborting.\n", ThisTask, All.TreeAllocFactor, treebuild_retries);
                dump_particles();
                endrun(1);
            }
            force_treefree();
            if(ThisTask == 0) {printf("Increasing TreeAllocFactor=%g", All.TreeAllocFactor);}
            All.TreeAllocFactor *= 1.15;
            if(ThisTask == 0) {printf(" new value=%g\n", All.TreeAllocFactor);}
            force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
        }
    }
    while(flag == -1);
    force_flag_localnodes();
    force_exchange_pseudodata();
    force_treeupdate_pseudos(All.MaxPart);
    TreeWalkValidatePending = 0; /* fresh tree: nothing pending from pre-build rearranges */
    if(mp == NULL && npart == NumPart) {TreeOpsCount[TREEOPS_BUILD]++; ForceAddElementToTree_CallsSinceBuild = 0; All.NumForcesSinceLastDomainDecomp = 0;} /* whole-tree builds only; subset builds are fof-internal. Both cadence counters measure tree age, so they reset HERE, not at decomposition (every decomposition is followed by a build anyway) */
#ifdef TREE_INTEGRITY_AUDITS
    if(mp == NULL && npart == NumPart) {force_tree_full_audit(1, "build");} /* whole-tree builds only: subset builds legitimately leave particles unreached */
#endif

    /* Tree-integrity invariant. N_part is exchanged with the pseudo-particle data and re-accumulated
     * across foreign domains just above, so the root node now counts every particle globally and must
     * agree with All.TotNumPart. A particle inserted under a top node owned by another task gets
     * detached by force_insert_pseudo_particles() and then appears in NO rank's tree: its mass is
     * absent from every multipole moment and the error is otherwise completely silent. This is an
     * integer comparison, which matters -- node masses are MyFloat, so one lost particle in ~1e6
     * is below the accumulated summation error and a mass-based check would miss it.
     * Only meaningful for whole-tree builds; FOF/SUBFIND build over subsets (mp != NULL).
     * Zero-mass particles (e.g. swallowed sinks awaiting cleanup) are legitimately uncounted when
     * they occupy a subtree by themselves, since internal nodes only propagate N_part when mass > 0,
     * so they are allowed for before flagging. Warn rather than abort: a false positive should not
     * be able to kill a long production run, and a warning is enough to stop this being silent. */
    if(mp == NULL)
    {
        int i; long long red_loc[3] = {0, 0, 0}, red_tot[3] = {0, 0, 0};
        for(i = 0; i < NumPart; i++) {if(P[i].Mass <= 0) {red_loc[0]++;}}
        /* whole-set veto, reduced with (not gated on) the zero-mass count: collective SUBFIND builds
           pass mp == NULL with npart = NumPartGroup != NumPart, whose root legitimately holds fewer
           than TotNumPart -- comparing there is a false alarm, not a detached particle. The predicate
           is rank-local, so it rides this reduction (SUM of 0/1 flags is an OR) rather than gating it,
           which would let ranks disagree about entering the collective. */
        red_loc[1] = (npart != NumPart) ? 1 : 0;
        red_loc[2] = NumPart; /* bookkeeping cross-check: sum vs All.TotNumPart, which is hand-maintained at ~6 mutation sites */
        MPI_Allreduce(red_loc, red_tot, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        long long n_zero_tot = red_tot[0], not_whole_any = red_tot[1];
        if(!not_whole_any && red_tot[2] != (long long) All.TotNumPart && ThisTask == 0)
        {
            printf("WARNING: particle bookkeeping drift: sum(NumPart)=%lld but All.TotNumPart=%lld. "
                   "One of the creation/destruction sites mis-updated the global count.\n", red_tot[2], (long long) All.TotNumPart);
            fflush(stdout);
        }
        /* The deficit is measured against the reduced count of particles that actually exist
           (red_tot[2]), not against All.TotNumPart: that global is hand-maintained, and once it has
           drifted -- the case the warning above reports -- a comparison against it would see a
           permanent deficit and buy a decomposition on every step without repairing anything. */
        long long in_tree = (long long) Nodes[All.MaxPart].N_part;
        long long deficit = red_tot[2] - in_tree;
        if(not_whole_any) {deficit = 0;} /* subset build: the count comparison is meaningless */
        if(deficit > n_zero_tot || deficit < 0)
        {
            DomainReconstructFlag = 1; /* self-heal: a deficit at a whole-tree build means ownership has
                drifted from geometry (escaped particles' subtrees are destroyed by the pseudo exchange).
                Rank-uniform: computed from the reduced totals. The next ladder decomposes, and the build
                after it re-covers the missing mass. Audited builds die at the audit instead. */
            if(ThisTask == 0)
            {
                printf("WARNING: tree integrity check failed: root node holds %lld particles, expected %lld "
                       "(deficit %lld, of which %lld zero-mass are allowed; All.TotNumPart=%lld). Particles missing "
                       "from the tree contribute to no rank's multipole moments, so forces are wrong by their mass. "
                       "Requesting a domain decomposition to restore ownership.\n",
                       in_tree, red_tot[2], deficit, n_zero_tot, (long long) All.TotNumPart);
                fflush(stdout);
            }
        }
    }

    TimeOfLastTreeConstruction = All.Time;
    return Numnodestree;
}



/*! Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.MaxPart.... All.MaxPart+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.MaxPart] gives the first tree node. Finally, node indices
 *  with values 'All.MaxPart + MaxNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the hydro
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
int force_treebuild_single(int npart, struct unbind_data *mp)
{
    int i, j, k, subnode = 0, shift, parent, numnodes, rep, nfree, th, nn, no;
    struct NODE *nfreep;
    MyFloat lenhalf;
    peanokey key, morton, th_key, *morton_list;
    
    /* create an empty root node  */
    nfree = All.MaxPart;        /* index of first free node */
    nfreep = &Nodes[nfree];    /* select first node */
    nfreep->len = DomainLen;
    nfreep->center = {(MyFloat)DomainCenter[0], (MyFloat)DomainCenter[1], (MyFloat)DomainCenter[2]};
    for(j = 0; j < 8; j++) {nfreep->u.suns[j] = -1;}
    numnodes = 1;
    nfreep++;
    nfree++;
    
    /* create a set of empty nodes corresponding to the top-level domain grid. We need to generate these nodes first to make sure that we have a
     * complete top-level tree which allows the easy insertion of the pseudo-particles at the right place */
    
    if(force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree) < 0) {return -1;}
    /* if a high-resolution region in a global tree is used, we need to generate an additional set empty nodes to make sure that we have a complete top-level tree for the high-resolution inset */
    nfreep = &Nodes[nfree];
    parent = -1;            /* note: will not be used below before it is changed */
    morton_list = (peanokey *) mymalloc("morton_list", NumPart * sizeof(peanokey));
    
    /* now we insert all particles */
    for(k = 0; k < npart; k++)
    {
        if(mp) {i = mp[k].index;} else {i = k;}
        rep = 0;
        /* must be the same mapping the domain decomposition assigned ownership with, hence the
         * shared definition rather than a local copy of the expression -- see domain_peano_key() */
        key = domain_peano_key(i, &morton);
        morton_list[i] = morton;
        shift = 3 * (BITS_PER_DIMENSION - 1);
        no = 0;
        while(TopNodes[no].Daughter >= 0)
        {
            no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);
            shift -= 3;
            rep++;
        }
        no = TopNodes[no].Leaf;
        th = DomainNodeIndex[no];
        
        while(1)
        {
            if(th >= All.MaxPart)    /* we are dealing with an internal node */
            {
                if(shift >= 0) {subnode = ((morton >> shift) & 7);}
                else
                {
                    subnode = 0;
                    if(P[i].Pos[0] > Nodes[th].center[0]) {subnode += 1;}
                    if(P[i].Pos[1] > Nodes[th].center[1]) {subnode += 2;}
                    if(P[i].Pos[2] > Nodes[th].center[2]) {subnode += 4;}
                }
                
                if(Nodes[th].len < EPSILON_FOR_TREERND_SUBNODE_SPLITTING * ForceSoftening_KernelRadius(i))
                {
                    /* seems like we're dealing with particles at identical (or extremely close) locations. Randomize subnode index to allow tree construction. Note: Multipole moments
                     * of tree are still correct, but this will only happen well below gravitational softening length-scale anyway. */
                    subnode = (int) (8.0 * get_random_number(P[i].ID));
                    if(subnode >= 8) {subnode = 7;}
                }
                
                nn = Nodes[th].u.suns[subnode];
                shift -= 3;
                
                if(nn >= 0)    /* ok, something is in the daughter slot already, need to continue */
                {
                    parent = th;
                    th = nn;
                    rep++;
                }
                else
                {
                    /* here we have found an empty slot where we can attach the new particle as a leaf. */
                    Nodes[th].u.suns[subnode] = i;
                    break;    /* done for this particle */
                }
            }
            else
            {
                /* We try to insert into a leaf with a single particle.  Need to generate a new internal node at this point. */
                Nodes[parent].u.suns[subnode] = nfree;
                nfreep->len = 0.5 * Nodes[parent].len;
                lenhalf = 0.25 * Nodes[parent].len;
                
                if(subnode & 1) {nfreep->center[0] = Nodes[parent].center[0] + lenhalf;}
                else {nfreep->center[0] = Nodes[parent].center[0] - lenhalf;}
                
                if(subnode & 2) {nfreep->center[1] = Nodes[parent].center[1] + lenhalf;}
                else {nfreep->center[1] = Nodes[parent].center[1] - lenhalf;}
                
                if(subnode & 4) {nfreep->center[2] = Nodes[parent].center[2] + lenhalf;}
                else {nfreep->center[2] = Nodes[parent].center[2] - lenhalf;}
                
                nfreep->u.suns[0] = -1;
                nfreep->u.suns[1] = -1;
                nfreep->u.suns[2] = -1;
                nfreep->u.suns[3] = -1;
                nfreep->u.suns[4] = -1;
                nfreep->u.suns[5] = -1;
                nfreep->u.suns[6] = -1;
                nfreep->u.suns[7] = -1;
                
                if(shift >= 0)
                {
                    th_key = morton_list[th];
                    subnode = ((th_key >> shift) & 7);
                }
                else
                {
                    subnode = 0;
                    if(P[th].Pos[0] > nfreep->center[0]) {subnode += 1;}
                    if(P[th].Pos[1] > nfreep->center[1]) {subnode += 2;}
                    if(P[th].Pos[2] > nfreep->center[2]) {subnode += 4;}
                }
                
                if(nfreep->len < EPSILON_FOR_TREERND_SUBNODE_SPLITTING * ForceSoftening_KernelRadius(th))
                {
                    /* seems like we're dealing with particles at identical (or extremely close) locations. Randomize subnode index to allow tree construction. Note: Multipole moments
                     * of tree are still correct, but this will only happen well below gravitational softening length-scale anyway. */
                    subnode = (int) (8.0 * get_random_number(P[th].ID));
                    if(subnode >= 8) {subnode = 7;}
                }
                
                nfreep->u.suns[subnode] = th;
                th = nfree;    /* resume trying to insert the new particle at the newly created internal node */
                numnodes++;
                nfree++;
                nfreep++;
                
                if((numnodes) >= MaxNodes)
                {
                    printf("task %d: maximum number %d of tree-nodes reached for particle %d.\n", ThisTask, MaxNodes, i);
                    /* always report capacity to the caller: the abort decision lives in the
                       force_treebuild driver, taken from the REDUCED flag so every rank agrees --
                       the old rank-local >5.0 endrun here could leave ranks divergent on whether
                       to abort */
                    myfree(morton_list);
                    return -1;
                }
            }
        }
    }
    
    myfree(morton_list);
    
    /* insert the pseudo particles that represent the mass distribution of other domains */
    force_insert_pseudo_particles();
    
    /* Seed Father[] for the recursion below, which is its only writer and only reaches particles
     * the u.suns[] walk visits. A local particle whose position keys into a top-node region owned
     * by another task is inserted under that node, and force_insert_pseudo_particles() above then
     * overwrites u.suns[0] there, detaching the subtree it sits in, so the recursion misses it.
     * Father comes from the LIFO mymalloc arena, so an unwritten slot holds the previous build's
     * value for that particle: usually a valid node index, which consumers then follow, silently
     * attributing the particle to the wrong node (and segfaulting when it is stale enough to fall
     * outside the node array). Consumers already treat a negative entry as "no parent"; seeding
     * -1 makes that convention real, so a missed particle is skipped instead.
     * This bounds the damage, it does not fix the orphaning: the particle is still absent from the
     * tree moments, so forces are still wrong by its mass. */
    for(i = 0; i < All.MaxPart; i++) {Father[i] = -1;}

    /* now compute the multipole moments recursively */
    last = -1;
    force_update_node_recursive(All.MaxPart, -1, -1);
    
    if(last >= All.MaxPart)
    {
        if(last >= All.MaxPart + MaxNodes) {Nextnode[last - MaxNodes] = -1;}    /* a pseudo-particle */
        else {Nodes[last].u.d.nextnode = -1;}
    }
    else {Nextnode[last] = -1;}
    
    return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 */
int force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
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

                    lenhalf = 0.25 * Nodes[no].len;
                    Nodes[*nextfree].len = 0.5 * Nodes[no].len;
                    Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
                    Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
                    Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;

                    for(n = 0; n < 8; n++)
                        Nodes[*nextfree].u.suns[n] = -1;

                    if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
                        DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;

                    *nextfree = *nextfree + 1;
                    *nodecount = *nodecount + 1;

                    if((*nodecount) >= MaxNodes)
                    {
                        printf("task %d: maximum number MaxNodes=%d of tree-nodes reached."
                               "MaxTopNodes=%d NTopnodes=%d NTopleaves=%d nodecount=%d\n",
                               ThisTask, MaxNodes, MaxTopNodes, NTopnodes, NTopleaves, *nodecount);
                        printf("in create empty nodes\n");
                        return -1; /* capacity signal; the collective abort decision lives in the force_treebuild driver */
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
            Nodes[index].u.suns[0] = All.MaxPart + MaxNodes + i;
    }
}


/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  Note that the bitflags-variable for each node is used to store in the
 *  lowest bits some special information: Bit 0 flags whether the node
 *  belongs to the top-level tree corresponding to the domain
 *  decomposition, while Bit 1 signals whether the top-level node is
 *  dependent on local mass/resolution elements.
 */
void force_update_node_recursive(int no, int sib, int father)
{
    int j, jj, k, p, pp, nextsib, suns[8], count_particles, multiple_flag;
    MyFloat hmax, vmax, v, divVmax, divVel, mass;
    Vec3<MyFloat> s, vs;
    struct particle_data *pa;
    
#ifdef DM_SCALARFIELD_SCREENING
    Vec3<MyFloat> s_dm, vs_dm; MyFloat mass_dm;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double cr_injection = 0;
#endif
#ifdef RT_USE_GRAVTREE
    MyFloat stellar_lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS]={0};
#endif
    for(j=0;j<N_RT_FREQ_BINS;j++) {stellar_lum[j]=0;}
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<MyFloat> rt_source_lum_s, rt_source_lum_vs;
#endif

    MyFloat maxsoft;

    if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)    /* internal node */
    {
        for(j = 0; j < 8; j++)
            suns[j] = Nodes[no].u.suns[j];    /* this "backup" is necessary because the nextnode entry will overwrite one element (union!) */
        if(last >= 0)
        {
            if(last >= All.MaxPart)
            {
                if(last >= All.MaxPart + MaxNodes)    /* a pseudo-particle */
                    Nextnode[last - MaxNodes] = no;
                else
                    Nodes[last].u.d.nextnode = no;
            }
            else
            {Nextnode[last] = no;}
        }
        
        last = no;
        
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        MyFloat gasmass = 0;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        cr_injection = 0;
#endif
#ifdef RT_USE_GRAVTREE
        for(j=0;j<N_RT_FREQ_BINS;j++) {stellar_lum[j]=0;}
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        rt_source_lum_s = rt_source_lum_vs = {};
#endif
#ifdef SINK_PHOTONMOMENTUM
        MyFloat sink_lum = 0; Vec3<MyFloat> sink_lum_grad = {};
#endif
#ifdef SINK_CALC_DISTANCES
        MyFloat sink_mass=0; Vec3<MyFloat> sink_pos_times_mass = {};   /* position of each sink particle in the node times its mass; divide by total mass at the end to get COM */
#ifdef SINK_NODE_MOTION_TRACKED
        Vec3<MyFloat> sink_mom = {}, sink_force = {}; int N_SINK = 0;
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
        vmax = 0;
        divVmax = 0;
        count_particles = 0;
        maxsoft = 0;
        
        for(j = 0; j < 8; j++)
        {
            if((p = suns[j]) >= 0)
            {
                /* check if we have a sibling on the same level */
                for(jj = j + 1; jj < 8; jj++)
                    if((pp = suns[jj]) >= 0)
                        break;
                
                if(jj < 8)    /* yes, we do */
                    nextsib = pp;
                else
                    nextsib = sib;
                
                force_update_node_recursive(p, nextsib, no);
                
                if(p >= All.MaxPart)    /* an internal node or pseudo particle */
                {
                    if(p >= All.MaxPart + MaxNodes)    /* a pseudo particle */
                    {
                        /* nothing to be done here because the mass of the
                         * pseudo-particle is still zero. This will be changed
                         * later.
                         */
                    }
                    else
                    {
                        mass += (Nodes[p].u.d.mass);
                        s += Nodes[p].u.d.mass * Nodes[p].u.d.s;
                        vs += Nodes[p].u.d.mass * Extnodes[p].vs;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                        gasmass += Nodes[p].gasmass;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                        cr_injection += Nodes[p].cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
                        for(k=0;k<N_RT_FREQ_BINS;k++) {stellar_lum[k] += (Nodes[p].stellar_lum[k]);}
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
#ifdef SINK_NODE_MOTION_TRACKED
                        sink_mom += Nodes[p].sink_mass * Nodes[p].sink_vel;
#ifdef SPECIAL_POINT_MOTION
                        sink_force += Nodes[p].sink_mass * Nodes[p].sink_acc;
#endif
                        N_SINK += Nodes[p].N_SINK;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                        if(Nodes[p].sink_mass > 0) {max_feedback_vel = DMAX(Nodes[p].MaxFeedbackVel, max_feedback_vel);}
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
                        if(Nodes[p].u.d.mass > 0) {count_particles += Nodes[p].N_part;} // we're saving the number of particles in the node, so simply add it
                        if(Extnodes[p].hmax > hmax) {hmax = Extnodes[p].hmax;}
                        if(Extnodes[p].vmax > vmax) {vmax = Extnodes[p].vmax;}
                        if(Extnodes[p].divVmax > divVmax) {divVmax = Extnodes[p].divVmax;}
                        
                        /* update of the maximum gravitational softening in the node */
                        if(Nodes[p].maxsoft > maxsoft) {maxsoft = Nodes[p].maxsoft;}
                        
                    }
                }
                else        /* a particle */
                {
                    count_particles++;
                    
                    pa = &P[p];
                    
                    mass += (pa->Mass);
                    s += pa->Mass * pa->Pos;
                    vs += pa->Mass * pa->Vel;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                    if(pa->Type == 0) gasmass += pa->Mass;
#if defined(SINK_ALPHADISK_ACCRETION) && defined(RT_USE_TREECOL_FOR_NH)
                    if(pa->Type == 5) gasmass += P[p].Sink_Mass_Reservoir; // gas at the inner edge of a disk should not see a hole due to the sink
#endif
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                    cr_injection += cr_get_source_injection_rate(p, P, CellP);
#endif
#ifdef RT_USE_GRAVTREE
                    double lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
                    double chimes_lum_G0[CHIMES_LOCAL_UV_NBINS];
                    double chimes_lum_ion[CHIMES_LOCAL_UV_NBINS];
                    int active_check = rt_get_source_luminosity_chimes(p,1,lum,chimes_lum_G0, chimes_lum_ion, P, CellP);
#else
                    int active_check = rt_get_source_luminosity(p,1,lum, P, CellP);
#endif
                    if(active_check)
                    {
                        double l_sum = 0; for(k=0;k<N_RT_FREQ_BINS;k++) {stellar_lum[k] += lum[k]; l_sum += lum[k];}
#ifdef CHIMES_STELLAR_FLUXES
                        for (k = 0; k < CHIMES_LOCAL_UV_NBINS; k++)
                        {
                            chimes_stellar_lum_G0[k] += chimes_lum_G0[k];
                            chimes_stellar_lum_ion[k] += chimes_lum_ion[k];
                        }
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
                        rt_source_lum_s += l_sum * pa->Pos;
                        rt_source_lum_vs += l_sum * pa->Vel;
#endif
                    }
#endif
                    
                    
                    
#ifdef SINK_PHOTONMOMENTUM
                    if(pa->Type == 5)
                    {
                        if((pa->Mass>0)&&(pa->DensityAroundParticle>0)&&(pa->Sink_Mdot>0))
                        {
                            double BHLum = sink_lum_bol(pa->Sink_Mdot, pa->Sink_Mass, p);
                            sink_lum += BHLum;
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
                            sink_lum_grad += pa->Sink_Specific_AngMom * BHLum;
#else
                            sink_lum_grad += pa->GradRho * BHLum;
#endif
                        }
                    }
#endif
#ifdef SINK_CALC_DISTANCES
                    if(pa->Type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
                    {
                        sink_mass += pa->Mass;    /* actual value is not used for distances */
                        sink_pos_times_mass += pa->Mass * pa->Pos;  /* positition times mass; divide by total mass later */
#ifdef SINK_NODE_MOTION_TRACKED
                        N_SINK += 1;
#endif
#ifdef SINK_NODE_MOTION_TRACKED
                        sink_mom += pa->Mass * pa->Vel;
#ifdef SPECIAL_POINT_MOTION
                        sink_force += pa->Mass * pa->Acc_Total_PrevStep;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                        max_feedback_vel = DMAX(pa->MaxFeedbackVel, max_feedback_vel);
#endif
#endif
                    }
#endif
                    
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                    {int k; for(k=0;k<6;k++) {tidal_tensorps_prevstep.data[k] += pa->Mass * pa->tidal_tensorps_prevstep.data[k];}}
#endif
                    
#ifdef DM_SCALARFIELD_SCREENING
                    if(pa->Type != 0)
                    {
                        mass_dm += (pa->Mass);
                        s_dm += pa->Mass * pa->Pos;
                        vs_dm += pa->Mass * pa->Vel;
                    }
#endif
                    if(pa->Type == 0)
                    {
                        double htmp = DMIN(All.MaxKernelRadius, P[p].KernelRadius);
                        if(htmp > hmax) {hmax = htmp;}
                        divVel = P[p].Particle_DivVel;
                        if(divVel > divVmax) {divVmax = divVel;}
                    }
                    
                    for(k = 0; k < 3; k++) {if((v = fabs(pa->Vel[k])) > vmax) {vmax = v;}}
                    
                    /* update of the maximum gravitational softening  */
                    double soft_p = ForceSoftening_KernelRadius(p);
                    if(soft_p > maxsoft) {maxsoft = soft_p;}
#ifdef SINGLE_STAR_SINK_DYNAMICS
                    if(pa->Type == 5) if(P[p].KernelRadius > maxsoft) {maxsoft = P[p].KernelRadius;}
#endif
                }
            }
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
        double l_tot=0; for(k=0;k<N_RT_FREQ_BINS;k++) {l_tot += stellar_lum[k];}
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
        } else {
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
        
        
        Nodes[no].Ti_current = All.Ti_Current;
        Nodes[no].u.d.mass = mass;
        Nodes[no].u.d.s = s;
        Nodes[no].GravCost = 0;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        Nodes[no].gasmass = gasmass;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        Nodes[no].cr_injection = cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
        for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k] = stellar_lum[k];}
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
        Extnodes[no].rt_source_lum_dp = {};
#endif
#ifdef SINK_PHOTONMOMENTUM
        Nodes[no].sink_lum = sink_lum;
        Nodes[no].sink_lum_grad = sink_lum_grad;
#endif
#ifdef SINK_CALC_DISTANCES
        Nodes[no].sink_mass = sink_mass;
#ifdef SINK_NODE_MOTION_TRACKED
        Nodes[no].N_SINK = N_SINK;
#endif
        if(sink_mass > 0)
        {
            Nodes[no].sink_pos = sink_pos_times_mass / sink_mass;  /* weighted position is sum(pos*mass)/sum(mass) */
#ifdef SINK_NODE_MOTION_TRACKED
            Nodes[no].sink_vel = sink_mom / sink_mass;
#endif
#if defined(SPECIAL_POINT_MOTION)
            Nodes[no].sink_acc = sink_force / sink_mass;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
            Nodes[no].MaxFeedbackVel = max_feedback_vel;
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
        Extnodes[no].dp_dm = {};
#endif
        
        Extnodes[no].Ti_lastkicked = All.Ti_Current;
        Extnodes[no].Flag = GlobFlag;
        Extnodes[no].vs = vs;
        Extnodes[no].hmax = hmax;
        Extnodes[no].vmax = vmax;
        Extnodes[no].divVmax = divVmax;
        Extnodes[no].dp = {};
#ifdef SINK_NODE_MOTION_TRACKED
        Extnodes[no].sink_dp = {}; /* sink_pos/sink_vel are set fresh here, so there is no pending kick to carry */
#endif

        Nodes[no].N_part = count_particles; /* save this value */
        if(count_particles > 1) {multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);} else {multiple_flag = 0;} /* this flags that the node represents more than one particle */
        Nodes[no].u.d.bitflags = multiple_flag;
        Nodes[no].maxsoft = maxsoft;
        Nodes[no].u.d.sibling = sib;
        Nodes[no].u.d.father = father;
    }
    else                /* single particle or pseudo particle */
    {
        if(last >= 0)
        {
            if(last >= All.MaxPart)
            {
                if(last >= All.MaxPart + MaxNodes) {Nextnode[last - MaxNodes] = no;} else {Nodes[last].u.d.nextnode = no;}    /* a pseudo-particle */
            } else {Nextnode[last] = no;}
        }
        last = no;
        if(no < All.MaxPart) {Father[no] = father;}    /* only set it for single particles */
    }
}




/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void force_exchange_pseudodata(void)
{
    int i, no, m, ta, recvTask;
    int *recvcounts, *recvoffset;
    struct DomainNODE
    {
        Vec3<MyFloat> s;
        Vec3<MyFloat> vs;
        MyFloat mass;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        MyFloat gasmass;
#endif
        MyFloat hmax;
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
#ifdef SINK_NODE_MOTION_TRACKED
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
    }
    *DomainMoment;
    
    
    DomainMoment = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));
    
    for(m = 0; m < MULTIPLEDOMAINS; m++)
        for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
            i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
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
#ifdef SINK_NODE_MOTION_TRACKED
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

    /* share the pseudo-particle data accross CPUs */
    recvcounts = (int *) mymalloc("recvcounts", sizeof(int) * NTask);
    recvoffset = (int *) mymalloc("recvoffset", sizeof(int) * NTask);
    
    for(m = 0; m < MULTIPLEDOMAINS; m++)
    {
        for(recvTask = 0; recvTask < NTask; recvTask++)
        {
            recvcounts[recvTask] =
            (DomainEndList[recvTask * MULTIPLEDOMAINS + m] - DomainStartList[recvTask * MULTIPLEDOMAINS + m] +
             1) * sizeof(struct DomainNODE);
            recvoffset[recvTask] = DomainStartList[recvTask * MULTIPLEDOMAINS + m] * sizeof(struct DomainNODE);
        }
        MPI_Allgatherv(MPI_IN_PLACE, recvcounts[ThisTask], MPI_BYTE, &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);
    }
    
    myfree(recvoffset);
    myfree(recvcounts);
    
    
    for(ta = 0; ta < NTask; ta++)
        if(ta != ThisTask)
            for(m = 0; m < MULTIPLEDOMAINS; m++)
                for(i = DomainStartList[ta * MULTIPLEDOMAINS + m]; i <= DomainEndList[ta * MULTIPLEDOMAINS + m]; i++)
                {
                    no = DomainNodeIndex[i];
                    
                    Nodes[no].u.d.s = DomainMoment[i].s;
                    Extnodes[no].vs = DomainMoment[i].vs;
                    Nodes[no].u.d.mass = DomainMoment[i].mass;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                    Nodes[no].gasmass = DomainMoment[i].gasmass;
#endif
                    Extnodes[no].hmax = DomainMoment[i].hmax;
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
#ifdef SINK_NODE_MOTION_TRACKED
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
}



/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_pseudos(int no)
{
    int j, p, count_particles, multiple_flag;
    MyFloat hmax, vmax;
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
#ifdef SINK_NODE_MOTION_TRACKED
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
    vmax = 0;
    divVmax = 0;
    count_particles = 0;
    maxsoft = 0;
    
    p = Nodes[no].u.d.nextnode;
    
    for(j = 0; j < 8; j++)    /* since we are dealing with top-level nodes, we now that there are 8 consecutive daughter nodes */
    {
        if(p >= All.MaxPart && p < All.MaxPart + MaxNodes)    /* internal node */
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
#ifdef SINK_NODE_MOTION_TRACKED
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
            if(Extnodes[p].vmax > vmax) {vmax = Extnodes[p].vmax;}
            if(Extnodes[p].divVmax > divVmax) {divVmax = Extnodes[p].divVmax;}
            if(Nodes[p].u.d.mass > 0) {count_particles += Nodes[p].N_part;} // saved, so directly add
            if(Nodes[p].maxsoft > maxsoft) {maxsoft = Nodes[p].maxsoft;}
        }
        else
            endrun(6767);        /* may not happen */
        
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
#ifdef SINK_NODE_MOTION_TRACKED
    Nodes[no].N_SINK = N_SINK;
#endif
    if(sink_mass > 0)
    {
        Nodes[no].sink_pos = sink_pos_times_mass / sink_mass;
#ifdef SINK_NODE_MOTION_TRACKED
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
    Extnodes[no].vmax = vmax;
    Extnodes[no].divVmax = divVmax;
    Extnodes[no].Flag = GlobFlag;
    Nodes[no].N_part = count_particles; // record
    if(count_particles > 1) {multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);} else {multiple_flag = 0;}
    Nodes[no].u.d.bitflags &= (~((1 << BITFLAG_MULTIPLEPARTICLES)));    /* this clears the bits */
    Nodes[no].u.d.bitflags |= multiple_flag;
    Nodes[no].maxsoft = maxsoft;
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
    
    for(m = 0; m < MULTIPLEDOMAINS; m++)
        for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
            i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
        {
            no = DomainNodeIndex[i];
            
            if(DomainTask[i] != ThisTask) {endrun(131231231);}
            
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
 *  link list. Multipole moments of tree nodes need not be changed -- but ONLY
 *  because every caller debits the child's mass from a co-located parent, so
 *  node total mass and COM are conserved. Type-conditioned moments (gasmass,
 *  stellar_lum, sink_*) are NOT updated here and go stale until the next
 *  refresh/rebuild, and hmax/vmax are bumped on the immediate parent only.
 */
void force_add_element_to_tree(int iparent, int ichild)
{
    ForceAddElementToTree_CallsSinceBuild++;
#ifndef MAINTAIN_TREE_IN_REARRANGE
    /* A pending rebuild means the walk links are (or are about to be) condemned. Splicing into
       them would let a neighbour search prune on one node's bounds while the traversal delivers
       the child through another. Keep only what the rebuild needs -- the child belongs with its
       parent -- and touch nothing else. (cf. gizmo-cpp a06e0073) */
    if(TreeReconstructFlag) {Father[ichild] = Father[iparent]; return;}
#endif
    int no;
    no = Nextnode[iparent];
    Nextnode[iparent] = ichild; // insert new particle into linked list
    Nextnode[ichild] = no; // order correctly
    Father[ichild] = Father[iparent]; // set parent node to be the same
    /* Father[] is seeded to -1 per treebuild, so a negative entry means the parent was never
     * reached by force_update_node_recursive (see the note there) and has no node to update.
     * Indexing Extnodes with it would be a wild write, so leave the opening-criterion fields
     * alone; the child inherits the same -1 and is skipped by the other consumers too. */
    if(Father[iparent] < 0) {return;}
    /* update parent node properties [maximum softening, speed] for opening criteria.
       The hmax bound must cover the CHILD's kernel: for a sink-spawn call the parent is the sink,
       whose kernel radius is the wrong particle's -- it both misses a larger gas child and
       pollutes the gas-hmax bound with sink kernels (the unsafe direction for neighbour
       completeness is the miss). */
    if(P[ichild].Type == 0) {Extnodes[Father[iparent]].hmax = DMAX(Extnodes[Father[iparent]].hmax, DMIN(P[ichild].KernelRadius, All.MaxKernelRadius));}
    double vmax = Extnodes[Father[iparent]].vmax;
    int k; for(k=0; k<3; k++) {if(fabs(P[ichild].Vel[k]) > vmax) {vmax = fabs(P[ichild].Vel[k]);}}
    Extnodes[Father[iparent]].vmax = vmax;
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
int force_treeevaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex)
{
    struct NODE *nop = 0;
    int no, nodesinlist=0, ptype, ninteractions=0, nexp, task, listindex = 0, maxPart = All.MaxPart;
    long bunchSize = All.BunchSize; int maxNodes = MaxNodes; integertime ti_Current = All.Ti_Current;
    double soft, r2, mass, r, fac_accel, u, h=0, h_p=0, h_inv, h3_inv, h_p_inv, h_p3_inv, u_p, xtmp, aold; xtmp=0; soft=0;
    Vec3<double> pos, dr; Vec3<MyDouble> acc = {};
    double pmass;
    double zeta=0, zeta_sec=0; int ptype_sec=-1;
#ifdef RT_USE_TREECOL_FOR_NH
    double angular_bin_size = 4*M_PI / RT_USE_TREECOL_FOR_NH, treecol_angular_bins[RT_USE_TREECOL_FOR_NH] = {0};
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    Vec3<double> dv;
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    double gasmass;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    Vec3<double> jerk = {};
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    Vec3<double> vel;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    double r_source, r_target, center[3]={0};
#ifdef BOX_PERIODIC
    center[0] = 0.5 * boxSize_X; center[1] = 0.5 * boxSize_Y; center[2] = 0.5 * boxSize_Z;
#endif
#endif
#ifdef PMGRID
    int tabindex; double eff_dist, rcut, asmth, asmthfac, rcut2, dist; dist = 0; rcut = All.Rcut[0]; asmth = All.Asmth[0];
    if(mode != 0 && mode != 1) {printf("%d %d %d %d %d\n", target, mode, *exportflag, *exportnodecount, *exportindex); endrun(444);}
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    MyFloat tree_mass = 0;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    int i1, i2; double fac2_tidal, fac_tidal; SymmetricTensor2<MyDouble> tidal_tensorps;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double cr_injection = 0;
#endif
#ifdef RT_USE_GRAVTREE
    double mass_stellarlum[N_RT_FREQ_BINS]; int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {mass_stellarlum[k_freq]=0;}
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_mass_stellarlum_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_mass_stellarlum_ion[CHIMES_LOCAL_UV_NBINS]={0}, chimes_flux_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_flux_ion[CHIMES_LOCAL_UV_NBINS]={0};
#endif
    Vec3<double> d_stellarlum = {}; int valid_gas_particle_for_rt = 0;
#ifdef RT_OTVET
    SymmetricTensor2<double> RT_ET[N_RT_FREQ_BINS]={};
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    double mass_sinklumwt_forradfb=0; // convert bh luminosity to our tree units
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
    double incident_flux_uv=0, incident_flux_euv=0;
#endif
#ifdef SINK_COMPTON_HEATING
    double incident_flux_agn=0;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    double SubGrid_CosmicRayEnergyDensity = 0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    double Rad_E_gamma[N_RT_FREQ_BINS]={0};
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    Vec3<double> Rad_Flux[N_RT_FREQ_BINS]; {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {Rad_Flux[kf] = {};}}
#endif
#ifdef SINK_CALC_DISTANCES
    double Min_Distance_to_Sink2=MAX_REAL_NUMBER; Vec3<double> Min_xyz_to_Sink = {MAX_REAL_NUMBER,MAX_REAL_NUMBER,MAX_REAL_NUMBER};
#ifdef SPECIAL_POINT_MOTION
    Vec3<double> vel_of_nearest_special = {}, acc_of_nearest_special = {};
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
    double weight_sum_for_special_point_smoothing = 0;
#endif
#endif
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
    double Min_Sink_OrbitalTime=MAX_REAL_NUMBER, comp_Mass; Vec3<double> comp_dx, comp_dv;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    double Min_Sink_Approach_Time = MAX_REAL_NUMBER, Min_Sink_Freefall_time = MAX_REAL_NUMBER;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    double Min_Sink_FeedbackTime = MAX_REAL_NUMBER;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Vec3<double> d_dm = {}; double mass_dm = 0;
#endif
#if defined(SINK_DYNFRICTION_FROMTREE)
    double sink_mass = 0, m_j_eff_for_df = 0;
#endif
#ifdef EVALPOTENTIAL
    double fac_pot; MyDouble pot; pot = 0;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    tidal_tensorps = {};
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    double tidal_zeta=0; SymmetricTensor2<MyFloat> i_zeta_tidal_tensorps_prevstep, j_zeta_tidal_tensorps_prevstep;
    if(mode==0) {i_zeta_tidal_tensorps_prevstep=P[target].tidal_tensorps_prevstep;} else {i_zeta_tidal_tensorps_prevstep=GravDataGet[target].tidal_tensorps_prevstep;}
#endif
#endif
    
#ifdef KETJU_REGULARIZATION
    MyIDType source_KetjuChainID = 0;
#endif
    if(mode == 0)
    {
        pos = P[target].Pos;
        ptype = P[target].Type;
        soft = ForceSoftening_KernelRadius(target);
        aold = All.ErrTolForceAcc * P[target].OldAcc;
        pmass = P[target].Mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
        vel = P[target].Vel;
#endif
#if defined(SINK_DYNFRICTION_FROMTREE)
        if(ptype == 5) {sink_mass = P[target].Sink_Mass;}
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
        if(ptype == 0) {if(soft > All.ForceSoftening[P[target].Type]) {zeta = P[target].AGS_zeta;} else {soft=All.ForceSoftening[P[target].Type]; zeta=0;}}
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
        if(soft > All.ForceSoftening[P[target].Type]) {zeta = P[target].AGS_zeta;} else {soft=All.ForceSoftening[P[target].Type]; zeta=0;}
#endif
#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
        if(pmforce_is_particle_high_res(ptype, P[target].Pos)) {rcut = All.Rcut[1]; asmth = All.Asmth[1];}
#endif
#ifdef KETJU_REGULARIZATION
        source_KetjuChainID = P[target].KetjuChainID;
#endif
    }
    else
    {
        pos = GravDataGet[target].Pos;
        ptype = GravDataGet[target].Type;
        soft = GravDataGet[target].Soft;
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
        pmass = GravDataGet[target].Mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
        vel = GravDataGet[target].Vel;
#endif
#if defined(SINK_DYNFRICTION_FROMTREE)
        if(ptype == 5) {sink_mass = GravDataGet[target].Sink_Mass;}
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        zeta = GravDataGet[target].AGS_zeta;
#endif
#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
        if(pmforce_is_particle_high_res(ptype, GravDataGet[target].Pos)) {rcut = All.Rcut[1]; asmth = All.Asmth[1];}
#endif
#ifdef KETJU_REGULARIZATION
        source_KetjuChainID = GravDataGet[target].KetjuChainID;
#endif
    }
    
    
    if(pmass<=0) {return 0;} /* quick check if particle has mass: if not, we won't deal with it */
    int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(ptype); // determine allowed particle types for correction terms for adaptive gravitational softening terms
#ifdef PMGRID
    rcut2 = rcut * rcut; asmthfac = 0.5 / asmth * (NTAB / 3.0);
#endif
#ifdef RT_USE_GRAVTREE
    if(ptype==0) {if((soft>0)&&(pmass>0)) {valid_gas_particle_for_rt = 1;}}
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    double fac_stellum[N_RT_FREQ_BINS];
    if(valid_gas_particle_for_rt)
    {
        double h_eff_phys = soft * pow(VOLUME_NORM_COEFF_FOR_NDIMS/All.DesNumNgb,1./NUMDIMS) * All.cf_atime; // convert from softening kernel extent to effective size, assuming 3D here, and convert to physical code units
        double sigma_particle =  pmass / (h_eff_phys*h_eff_phys); // quick estimate of effective surface density of the target, in physical code units
        double fac_stellum_0 = -All.PhotonMomentum_Coupled_Fraction / (4.*M_PI * C_LIGHT_CODE_REDUCED * sigma_particle * All.G); // this will be multiplied by L/r^2 below, giving acceleration, then extra G because code thinks this is gravity, so put extra G here. everything is in -physical- code units here //
        int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {fac_stellum[kf] = fac_stellum_0*(1 - exp(-rt_kappa(-1,kf, P, CellP)*sigma_particle));} // rt_kappa is in physical code units, so sigma_eff_abs should be also -- approximate surface-density through particle (for checking if we enter optically-thick limit)
    }
#endif
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    double m_enc_in_rcrit = 0, r_for_total_menclosed = soft; r_for_total_menclosed = DMAX( r_for_total_menclosed , 0.1/(UNIT_LENGTH_IN_KPC*All.cf_atime) ); /* set a baseline Rcrit_min, otherwise we get statistics that are very noisy */
#endif
    
    
    if(mode == 0)
    {
        no = maxPart;        /* root node */
    }
    else
    {
        nodesinlist++;
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;    /* open it */
    }
    
    while(no >= 0)
    {
        while(no >= 0)
        {
            h=soft; h_p=-1; /* initialize h and h_p, for use below: make sure to do so at the top of each iteration */
            
            if(no < maxPart) /* this is a particle, we will use it */
            {
#ifdef SINGLE_STAR_DIRECT_GRAVITY
                /* star-star pairs are summed exactly in star_direct_gravity_compute(); taking them
                   here as well would double every such force */
                if((ptype == 5) && (P[no].Type == 5)) {no = Nextnode[no]; continue;}
#endif
                /* the index of the node is the index of the particle */
                if(P[no].Ti_current != ti_Current)
                {
#ifdef _OPENMP
#pragma omp critical(_particledriftforce_)
#endif
                    {
                        if(P[no].Ti_current != ti_Current) {drift_particle(no, ti_Current);}
                    }
                }
                dr = P[no].Pos - pos;
                GRAVITY_NEAREST_XYZ(dr[0],dr[1],dr[2],-1);
                r2 = dr.norm_sq();
                mass = P[no].Mass;
                
#ifdef GRAVITY_SPHERICAL_SYMMETRY
                r_source = sqrt(pow(P[no].Pos[0] - center[0],2) + pow(P[no].Pos[1] - center[1],2) + pow(P[no].Pos[2] - center[2],2));
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
                dv = P[no].Vel - vel;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
                /* the source velocity this interaction is evaluated at; the Hermite predictor below
                   replaces it, so the (dr, dv) pair feeding Min_Sink_Approach_Time and the binary
                   pairing stays mutually consistent on a Hermite pass instead of mixing a predicted
                   dr with the whole-step-kicked stored Vel */
                Vec3<double> src_vel_asused = P[no].Vel;
#endif
#ifdef HERMITE_INTEGRATION
                /* Hermite-only passes: an inactive Hermite source sits at its KDK-drifted position
                 * with a whole-step-kicked Vel (see do_the_kick), O(dt^2) wrong mid-step, so
                 * re-predict it from its own Old* state as do_hermite_prediction does for the
                 * target. Single-particle branch only: a source inside a node multipole keeps its
                 * drifted position, but SINGLE_STAR_DIRECT_GRAVITY_RADIUS opens star-bearing nodes
                 * nearby, so the close pairs where this matters take this branch. Nothing is
                 * written back. Active sources are skipped: at HermiteOnlyFlag==1 their Old* are
                 * stale and their live Pos is already correct. The KDK pass is untouched, since the
                 * leapfrog's cancellation assumes plain drifted positions.
                 *
                 * The Old* polynomial only interpolates over the source's own step. Past that it
                 * extrapolates outside its convergence radius and, for a close pair, biases the
                 * source along the orbit, which rectifies into COM drift -- so outside the span the
                 * drifted state is kept: lower order, but unbiased. The integer span checks come
                 * first because eligible_for_hermite() is the costly part: it depends only on the
                 * source but runs per pair, and being cross-TU the compiler cannot hoist it. */
                if(HermiteOnlyFlag && !TimeBinActive[P[no].TimeBin]
                   && (ti_Current >= P[no].Ti_begstep)
                   && ((ti_Current - P[no].Ti_begstep) <= P[no].integertime_step())
                   && eligible_for_hermite(no))
                {
                    double hD = get_gravkick_factor(P[no].Ti_begstep, ti_Current, no, 0);
                    dr = P[no].OldPos + (P[no].OldVel + (P[no].Hermite_OldAcc + P[no].OldJerk * (hD/3)) * (hD/2)) * hD - pos;
                    GRAVITY_NEAREST_XYZ(dr[0],dr[1],dr[2],-1);
                    r2 = dr.norm_sq();
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
                    dv = P[no].OldVel + (P[no].Hermite_OldAcc + P[no].OldJerk * (hD/2)) * hD - vel;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
                    src_vel_asused = P[no].OldVel + (P[no].Hermite_OldAcc + P[no].OldJerk * (hD/2)) * hD;
#endif
                }
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
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
                    if(ptype == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
                    {
                        int kx; double wt_special = weight_function_for_weighted_motion_smoothing(sqrt(r2),1);
                        weight_sum_for_special_point_smoothing += wt_special;
                        vel_of_nearest_special += wt_special * P[no].Vel;
                        acc_of_nearest_special += wt_special * P[no].Acc_Total_PrevStep;
                    }
#endif
                    if(P[no].Type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) /* found a BH particle in grav calc */
                    {
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
                        if(ptype != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
#endif
                            if(r2 < Min_Distance_to_Sink2)    /* is this the closest BH part I've found yet? */
                            {
                                Min_Distance_to_Sink2 = r2;   /* if yes: adjust min bh dist */
                                Min_xyz_to_Sink = dr; /* remember, dr = x_SINK - myx */
#ifdef SPECIAL_POINT_MOTION
                                int kx;
                                vel_of_nearest_special = P[no].Vel;
                                acc_of_nearest_special = P[no].Acc_Total_PrevStep;
#endif
                            }
#ifdef SINGLE_STAR_TIMESTEPPING
                        Vec3<double> sink_dv=src_vel_asused-vel; double vSqr=sink_dv.norm_sq(), M_total=P[no].Mass+pmass, r2soft=SinkParticle_GravityKernelRadius;
                        r2soft = DMAX(r2soft, soft);
                        r2soft *= KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER;
                        r2soft = r2 + r2soft*r2soft;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                        if(ptype == 0) {
                            double tSqr_fb = r2soft /(P[no].MaxFeedbackVel * P[no].MaxFeedbackVel + MIN_REAL_NUMBER);
                            if(tSqr_fb < Min_Sink_FeedbackTime) {Min_Sink_FeedbackTime = tSqr_fb;}
                        } // for gas, add the signal velocity of feedback from the star
#endif
                        double tSqr = r2soft/(vSqr + MIN_REAL_NUMBER), tff4 = r2soft*r2soft*r2soft/(M_total*M_total);

                        /* NOTE: chain-ID-based skip for KETJU is implemented
                         * (source_KetjuChainID vs P[no].KetjuChainID) but currently
                         * disabled — see timestep-optimization note in timestep.cc.
                         * When enabled, same-chain Type-5 neighbors are excluded so
                         * dt_2body reflects only non-chain encounters. */
                        if(tSqr < Min_Sink_Approach_Time) {Min_Sink_Approach_Time = tSqr;}
                        if(tff4 < Min_Sink_Freefall_time) {Min_Sink_Freefall_time = tff4;}
#ifdef SINGLE_STAR_FIND_BINARIES
                        if(ptype == 5) // only for BH particles and for non center of mass calculation
                        {
                            double r_p5=sqrt(r2), specific_energy = 0.5*vSqr - All.G*M_total/r_p5;
                            if(r2 < SinkParticle_GravityKernelRadius*SinkParticle_GravityKernelRadius)
                            {
                                double hinv_p5 = 1. / SinkParticle_GravityKernelRadius;
                                specific_energy = 0.5*vSqr + All.G*M_total*kernel_gravity(r_p5*hinv_p5, hinv_p5, hinv_p5*hinv_p5*hinv_p5, -1);
                            }
                            if (specific_energy < 0)
                            {
                                double semimajor_axis= -All.G*M_total/(2.*specific_energy);
                                double t_orbital = 2.*M_PI*sqrt( semimajor_axis*semimajor_axis*semimajor_axis / (All.G*M_total) );
                                if(t_orbital < Min_Sink_OrbitalTime) /* Save parameters of companion */
                                {
                                    Min_Sink_OrbitalTime=t_orbital; comp_Mass=P[no].Mass;
                                    comp_dx=dr; comp_dv=sink_dv;
                                }
                            } /* specific_energy < 0 */
                        } /* ptype == 5 */
#endif //#ifdef SINGLE_STAR_FIND_BINARIES
#endif //#ifdef SINGLE_STAR_TIMESTEPPING
                    }
#endif // SINK_CALC_DISTANCES

#ifdef COSMIC_RAY_SUBGRID_LEBRON
                    cr_injection = cr_get_source_injection_rate(no, P, CellP);
#endif

#ifdef RT_USE_GRAVTREE
                    if(valid_gas_particle_for_rt)    /* we have a (valid) gas particle as target */
                    {
                        d_stellarlum=dr;
                        double lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
                        double chimes_lum_G0[CHIMES_LOCAL_UV_NBINS], chimes_lum_ion[CHIMES_LOCAL_UV_NBINS];
                        int active_check = rt_get_source_luminosity_chimes(no,1,lum, chimes_lum_G0, chimes_lum_ion, P, CellP);
#else
                        int active_check = rt_get_source_luminosity(no,1,lum, P, CellP);
#endif
                        int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(active_check) {mass_stellarlum[kf]=lum[kf];} else {mass_stellarlum[kf]=0;}}
#ifdef CHIMES_STELLAR_FLUXES
                        for(kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++)
                        {
                            if(active_check) {chimes_mass_stellarlum_G0[kf] = chimes_lum_G0[kf]; chimes_mass_stellarlum_ion[kf] = chimes_lum_ion[kf];} else {chimes_mass_stellarlum_G0[kf] = 0; chimes_mass_stellarlum_ion[kf] = 0;}
                        }
#endif
#ifdef SINK_PHOTONMOMENTUM
                        mass_sinklumwt_forradfb=0;
                        if(P[no].Type == 5)
                        {
                            double bhlum_t = sink_lum_bol(P[no].Sink_Mdot, P[no].Sink_Mass, no);
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
                            mass_sinklumwt_forradfb = sink_fb_angleweight(bhlum_t, P[no].Sink_Specific_AngMom, dr[0],dr[1],dr[2]);
#else
                            mass_sinklumwt_forradfb = sink_fb_angleweight(bhlum_t, P[no].GradRho, dr[0],dr[1],dr[2]);
#endif
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
            else /* we have an  internal node */
            {
                if(no >= maxPart + maxNodes) /* pseudo particle -- this will not be used for calculations below, but needs to be parsed here */
                {
                    if(mode == 0)
                    {
                        if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
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
                                    /* out of buffer space. Need to discard work for this particle and interrupt */
                                    BufferFullFlag = 1;
                                    exitFlag = 1;
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
                        DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
                        DomainNodeIndex[no - (maxPart + maxNodes)];
                        if(exportnodecount[task] < NODELISTLENGTH) {DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;}
                    }
                    no = Nextnode[no - maxNodes];
                    continue;
                }
                /* ok we have an internal node on the local processor, need to decide if we open it and go further or keep it */
                nop = &Nodes[no];
                
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }
                mass = nop->u.d.mass;
                if(mass <= 0) /* nothing in the node */
                {
                    no = nop->u.d.sibling;
                    continue;
                }
#ifdef SINGLE_STAR_DIRECT_GRAVITY
                /* A star target must take no star mass from the tree, since star_direct_gravity_compute()
                   supplies every star-star pair exactly. This tree carries monopoles only (u.d.mass at
                   u.d.s; there are no quadrupole moments in struct NODE), so removing the sinks is exact
                   rather than an approximation: drop their mass, and shift the node's center of mass to
                   the center of mass of what is left. Both terms are on the same clock -- SINK_NODE_MOTION_TRACKED
                   drifts sink_pos with sink_vel exactly as u.d.s is drifted with vs, and kicks sink_vel
                   from the sinks' own momentum -- so the difference below really is the non-sink COM.
                   Nodes made entirely of stars vanish and are skipped; everything else stays usable at
                   monopole level, so a star target keeps the tree's O(log N) walk instead of having to
                   open every node that happens to contain a star. */
                if((ptype == 5) && (nop->sink_mass > 0) && (mass - nop->sink_mass <= 0)) {no = nop->u.d.sibling; continue;} /* pure-star node: nothing here for us. Checked before the drift below so we never bother drifting a node we discard; the COM shift itself has to wait until after it. */
#endif
                //if(nop->N_part <= 1)
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    if(mass) /* open cell */
                    {
                        no = nop->u.d.nextnode;
                        continue;
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

#ifdef SINGLE_STAR_DIRECT_GRAVITY
                /* Remove the sinks from this node for a star target: star-star pairs come exactly from
                   star_direct_gravity_compute(), so taking them here too would double them. The tree
                   carries monopoles only (u.d.mass at u.d.s -- struct NODE has no quadrupole moments),
                   so this subtraction is exact rather than approximate: drop the sink mass and move the
                   center of mass to that of what remains. Both terms are on the same clock, since
                   SINK_NODE_MOTION_TRACKED drifts sink_pos with sink_vel exactly as u.d.s is drifted with vs
                   -- which is also why this sits after force_drift_node above. Doing it this way means a
                   star target keeps the ordinary O(log N) walk; the alternative, opening every node that
                   contains a star, costs an extra O(N_star log N) node visits per star target.
                   Note mass is reduced here, before the opening criteria below, so they judge the node on
                   the mass actually being used. */
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
#ifdef PMGRID
                if(r2 > rcut2) /* check whether we can stop walking along this branch */
                {
                    eff_dist = rcut + 0.5 * nop->len;
                    dist = GRAVITY_NGB_PERIODIC_BOX_LONG_X(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1);
                    if(dist > eff_dist) {no = nop->u.d.sibling; continue;}
                    dist = GRAVITY_NGB_PERIODIC_BOX_LONG_Y(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1);
                    if(dist > eff_dist) {no = nop->u.d.sibling; continue;}
                    dist = GRAVITY_NGB_PERIODIC_BOX_LONG_Z(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1);
                    if(dist > eff_dist) {no = nop->u.d.sibling; continue;}
                }
#endif // PMGRID //
#ifdef NEIGHBORS_MUST_BE_COMPUTED_EXPLICITLY_IN_FORCETREE
                Vec3<double> d_nc = nop->center - pos;
                GRAVITY_NEAREST_XYZ(d_nc[0],d_nc[1],d_nc[2],-1); /* find the closest image in the given box size  */
                double dist_to_center2 = d_nc.norm_sq();
                double dist_to_open = DMAX(soft , nop->maxsoft) + nop->len*1.73205/2.0;
                if(dist_to_center2  < dist_to_open*dist_to_open) /* check if any portion the cell lies within the interaction range, then open cell */
                {
                    no = nop->u.d.nextnode;
                    continue;
                }
#else
                if(h < nop->maxsoft) // compare primary softening to node maximum
                {
                    if(r2 < nop->maxsoft * nop->maxsoft) {no = nop->u.d.nextnode; continue;} // inside node maxsoft, continue down tree
                }
#endif
                if(All.ErrTolTheta)    /* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta) /* open cell */
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#ifndef GRAVITY_HYBRID_OPENING_CRIT
                else        /* check relative opening criterion */
#else
                    if(!(All.Ti_Current == 0 && RestartFlag != 1))
#endif
                    {
                        /* force node to open if we are within the gravitational softening length */
                        if((r2 < (soft+0.6*nop->len)*(soft+0.6*nop->len)) || (r2 < (nop->maxsoft+0.6*nop->len)*(nop->maxsoft+0.6*nop->len)))
                        {
                            no = nop->u.d.nextnode;
                            continue;
                        }
                        if(mass * nop->len * nop->len > r2 * r2 * aold) /* open cell */
                        {
                            no = nop->u.d.nextnode;
                            continue;
                        }
                        /* check in addition whether we lie inside the cell */
                        if(GRAVITY_NGB_PERIODIC_BOX_LONG_X(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1) < 0.60 * nop->len)
                        {
                            if(GRAVITY_NGB_PERIODIC_BOX_LONG_Y(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1) < 0.60 * nop->len)
                            {
                                if(GRAVITY_NGB_PERIODIC_BOX_LONG_Z(nop->center[0] - pos[0], nop->center[1] - pos[1], nop->center[2] - pos[2], -1) < 0.60 * nop->len)
                                {
                                    no = nop->u.d.nextnode;
                                    continue;
                                }
                            }
                        }
#if (defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)) && defined(SINGLE_STAR_DIRECT_GRAVITY_RADIUS)
                        if(ptype == 5) {
                            if((nop->N_SINK > 0) && (r2 < pow(SINGLE_STAR_DIRECT_GRAVITY_RADIUS/UNIT_LENGTH_IN_AU + 0.6*nop->len,2))) // we are a star looking at another star within the specified radius, open cell to get direct force summation
                            {
                                no = nop->u.d.nextnode;
                                continue;
                            }
                        }
#endif
                    }
                
                /* ok we will be using this node, can now set variables that depend on it */
                h_p = nop->maxsoft;
                zeta_sec = 0; ptype_sec = -1; /* set secondary softening and zeta terms */
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                gasmass = nop->gasmass;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
                r_source = sqrt(pow(nop->u.d.s[0] - center[0],2) + pow(nop->u.d.s[1] - center[1],2) + pow(nop->u.d.s[2] - center[2],2));
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
                if(ptype == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
                {
                    int kx; double wt_special = weight_function_for_weighted_motion_smoothing(sqrt(r2),1);
                    weight_sum_for_special_point_smoothing += wt_special;
                    vel_of_nearest_special = Extnodes[no].vs;
                    acc_of_nearest_special = {}; /* no accel for now, that will be computed later, but keep this if needed */
                }
#endif
                if(nop->sink_mass > 0)        /* found a node with non-zero BH mass */
                {
                    Vec3<double> sink_dr = nop->sink_pos - pos;  /* SHEA:  now using sink_pos instead of center */
                    GRAVITY_NEAREST_XYZ(sink_dr[0],sink_dr[1],sink_dr[2],-1);
                    double sink_r2 = sink_dr.norm_sq(); // + (nop->len)*(nop->len);
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
                    if(ptype != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
#endif
                        if(sink_r2 < Min_Distance_to_Sink2)
                        {
                            Min_Distance_to_Sink2 = sink_r2; Min_xyz_to_Sink = sink_dr; /* remember, dr = x_SINK - myx */
#ifdef SPECIAL_POINT_MOTION
                            int kx;
                            vel_of_nearest_special = nop->sink_vel;
                            acc_of_nearest_special = nop->sink_acc;
#endif
                        }
#ifdef SINGLE_STAR_TIMESTEPPING
                    Vec3<double> sink_dv = nop->sink_vel - vel; double vSqr=sink_dv.norm_sq(), M_total=nop->sink_mass+pmass, r2soft;
                    r2soft = DMAX(SinkParticle_GravityKernelRadius, soft) * KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER; r2soft = sink_r2 + r2soft*r2soft;
                    double tSqr = r2soft/(vSqr + MIN_REAL_NUMBER), tff4 = r2soft*r2soft*r2soft/(M_total*M_total);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                    if(ptype == 0) {
                        double tSqr_fb = r2soft /(nop->MaxFeedbackVel * nop->MaxFeedbackVel + MIN_REAL_NUMBER);
                        if(tSqr_fb < Min_Sink_FeedbackTime) {Min_Sink_FeedbackTime = tSqr_fb;}
                    } // for gas, add the signal velocity of feedback from the star
#endif
                    if(tSqr < Min_Sink_Approach_Time) {Min_Sink_Approach_Time = tSqr;}
                    if(tff4 < Min_Sink_Freefall_time) {Min_Sink_Freefall_time = tff4;}
#ifdef SINGLE_STAR_FIND_BINARIES
                    if(ptype == 5 && nop->N_SINK == 1) // only do it if we're looking at a single star in the node
                    {
                        double specific_energy = 0.5*vSqr - All.G*M_total/sqrt(r2);
                        if (specific_energy<0)
                        {
                            double semimajor_axis= -All.G*M_total/(2.*specific_energy);
                            double t_orbital = 2.*M_PI*sqrt( semimajor_axis*semimajor_axis*semimajor_axis / (All.G*M_total) );
                            if(t_orbital < Min_Sink_OrbitalTime) /* Save parameters of companion */
                            {
                                Min_Sink_OrbitalTime=t_orbital; comp_Mass=nop->sink_mass;
                                comp_dx = sink_dr; comp_dv = sink_dv;
                            }
                        } /* specific_energy < 0 */
                    } /* ptype == 5 */
#endif //#ifdef SINGLE_STAR_FIND_BINARIES
#endif //#ifdef SINGLE_STAR_TIMESTEPPING
                }
#endif // SINK_CALC_DISTANCES
                
            } /* ok we've completed all the opening criteria -- we will keep this node or particle as-is */
            
            
            
            
            if((r2 > 0) && (mass > 0)) // only go forward if mass positive and there is separation -- this is check for the whole block below, which should no include 'self' terms
            {
                r = sqrt(r2);
                
                /* now we compute the actual pair-wise gravity terms */
                if((r >= h) && (r >= h_p)) // can safely do a purely-Newtonian force (can be done by kernel-gravity as well, but no need to enter all the conditional statements below so just do it here for simplicity //
                {
                    fac_accel = mass / (r2 * r);
#ifdef EVALPOTENTIAL
                    fac_pot = -mass / r;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                    fac_tidal = fac_accel; fac2_tidal = 3.0 * mass / (r2 * r2 * r); /* second derivative of potential needs this factor */
#endif
                }
                else
                {
                    double h_grav = h;
#if !defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
                    if(h_p > h_grav) {h_grav = h_p;} // in this case, symmetrize by taking the maximum here always
#endif
                    h_inv=1./h_grav; h3_inv=h_inv*h_inv*h_inv; u=r*h_inv; // set here to ensure this is using the correct values //
                    fac_accel = mass * kernel_gravity(u, h_inv, h3_inv, 1);
#ifdef EVALPOTENTIAL
                    fac_pot = mass * kernel_gravity(u, h_inv, h3_inv, -1);
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* second derivatives needed -> calculate them from softened potential */
                    fac_tidal = fac_accel; fac2_tidal = mass * kernel_gravity(u, h_inv, h3_inv, 2);  /* save original fac_accel without shortrange_table factor or zeta terms (needed for tidal field calculation) */
#endif
                    
#if defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
                    if(h_p > 0) // first, appropriately symmetrize the forces between particles. only do this is secondary is a particle, so has a type and softening! //
                    {
                        int symmetrize_by_averaging = 0; // default here to symmetrize by taking the maximum, but this will vary below //
                        // the 'zeta' terms for conservation with adaptive softening assume kernel-scale forces are averaged to symmetrize, to make them continuous
                        if(ptype_sec>=0) {if((1 << ptype_sec) & (AGS_kernel_shared_BITFLAG)) {symmetrize_by_averaging=1;}} // symmetrize by averaging only for particles which have a shared AGS structure since this is how our correction terms are derived //
#ifdef SINGLE_STAR_SINK_DYNAMICS
                        if((ptype!=0) || (ptype_sec!=0)) {symmetrize_by_averaging=0;} // we don't want to do the symmetrization below for sink interactions because it can create very noisy interactions between tiny sink particles and diffuse gas. However we do want it for gas-gas interactions so we keep the below
#endif
                        double prefac_corr_p=1., prefac_corr_orig=1.; // this will give a symmetrized pair by linear averaging
                        if(symmetrize_by_averaging==0) {prefac_corr_p=2; prefac_corr_orig=0.;} // symmetrize instead with the old method of simply taking the larger of the pair. here only act if the softening of the particle whose force is being summed is greater than the target //
                        if((symmetrize_by_averaging==1) || (h_p>h)) // condition to need to evaluate the alternate particle ('p' side)
                        {
                            h_p_inv=1./h_p; h_p3_inv=h_p_inv*h_p_inv*h_p_inv; u_p=r*h_p_inv;
                            fac_accel = 0.5 * (prefac_corr_orig * fac_accel + prefac_corr_p * mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 1)); // average with neighbor
#ifdef EVALPOTENTIAL
                            fac_pot = 0.5 * (prefac_corr_orig * fac_pot + prefac_corr_p * mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, -1)); // average with neighbor
#endif
#if defined(COMPUTE_TIDAL_TENSOR_IN_GRAVTREE)
                            fac_tidal = fac_accel; fac2_tidal = 0.5 * (prefac_corr_orig * fac2_tidal + prefac_corr_p * mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 2)); // average forces -> average in tidal tensor as well. also save updated fac_tidal
#endif
                        }
                    } // closes if((h_p > 0)) clause
#endif // closes clause to symmetrize by averaging instead of taking the larger softening
                    
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                    double fac_corr = 0; int add_ags_zeta_terms_primary=1, add_ags_zeta_terms_secondary=1; u_p=r/h_p; // define the correction factor but also a clause to see if we should apply any 'correction' term at all
                    if((r<=0) || (pmass<=0) || (mass<=0) || (ptype_sec<0)) {add_ags_zeta_terms_primary=0; add_ags_zeta_terms_secondary=0;} // define conditions to add these terms at all (don't go forward if any of these conditions are met)
                    if((zeta == 0) || (u >= 1) || (h <= 0)) {add_ags_zeta_terms_primary=0;} // other conditions that mean -dont- use the term for the ab side
                    if((zeta_sec == 0) || (u_p >= 1) || (h_p <= 0)) {add_ags_zeta_terms_secondary=0;} // other conditions that mean -dont- use the term for the ba side
                    // correction only applies to 'shared-kernel' particles: so this needs to check if these are the same particles for which the 'shared' kernel lengths are computed
#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
                    if(ptype != 0 || ptype_sec != 0) {add_ags_zeta_terms_primary=0; add_ags_zeta_terms_secondary=0;} // primary and secondary must be gas for ab side or ba side
#else
                    if(!((1 << ptype) & (ADAPTIVE_GRAVSOFT_FORALL)) || !((1 << ptype_sec) & (ags_gravity_kernel_shared_BITFLAG(ptype)))) {add_ags_zeta_terms_primary=0;} // primary must be a valid ags particle and 'see' secondary for ab side
                    if(!((1 << ptype_sec) & (ADAPTIVE_GRAVSOFT_FORALL)) || !((1 << ptype) & (ags_gravity_kernel_shared_BITFLAG(ptype_sec)))) {add_ags_zeta_terms_secondary=0;} // secondary must be a valid ags particle and 'see' primary for ba side
#endif
                    if(add_ags_zeta_terms_primary) // ab side
                    {
                        double dWdr, wp; kernel_main(u, h3_inv, h3_inv*h_inv, &wp, &dWdr, 1);
                        fac_corr += -(zeta/pmass) * dWdr / r; // go ahead and add the term
                    }
                    if(add_ags_zeta_terms_secondary) // ba side
                    {
                        double dWdr, wp; h_p_inv=1./h_p; h_p3_inv=h_p_inv*h_p_inv*h_p_inv; kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wp, &dWdr, 1);
                        fac_corr += -(zeta_sec/pmass) * dWdr / r; // go ahead and add the term
                    }
                    if(!isnan(fac_corr)) {fac_accel += fac_corr;}
#endif
                } // closes r < h (else) clause [where we need to deal with inside-the-softening factors]
                
                
#ifdef PMGRID
                tabindex = (int) (asmthfac * r);
                if(tabindex < NTAB && tabindex >= 0)
#endif // PMGRID //
                {
#ifdef PMGRID
                    fac_accel *= shortrange_table[tabindex];
#endif
#ifdef EVALPOTENTIAL
#ifdef PMGRID
                    fac_pot *= shortrange_table_potential[tabindex];
#endif
                    pot += (fac_pot);
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
                    pot += (mass * ewald_pot_corr(dr[0], dr[1], dr[2]));
#endif
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
                    r_target = sqrt(pow(pos[0] - center[0],2) + pow(pos[1] - center[1],2) + pow(pos[2] - center[2],2)); // distance of target point from box center
                    if(r_source < r_target) {dr[0] = center[0] - pos[0]; dr[1] = center[1] - pos[1]; dr[2] = center[2] - pos[2]; fac_accel = mass/pow(DMAX(GRAVITY_SPHERICAL_SYMMETRY,DMAX(r_target,h)),3);} else {fac_accel = 0;}
#endif
                    
                    /* actually add the accelerations, now that we've corrected for the ewald and other terms */
                    acc += fac_accel * dr;
                    
                    
#if defined(SINK_DYNFRICTION_FROMTREE)
                    if( (fac_accel>MIN_REAL_NUMBER) && (ptype==5) && (mass>MIN_REAL_NUMBER) )
                    {
                        double dv2=dv.norm_sq();
                        if((dv2 > MIN_REAL_NUMBER) && (sink_mass > MIN_REAL_NUMBER))
                        {
                            double dv0=sqrt(dv2); Vec3<double> dv_h=dv/dv0;
                            double rdotvhat=dot(dr,dv_h);
                            Vec3<double> b_im=dr-rdotvhat*dv_h; double b_impact=b_im.norm();
                            double a_im=(b_impact*All.cf_atime)*(dv2*All.cf_a2inv)/(All.G*sink_mass), fac_df=fac_accel*b_impact*a_im/(1.+a_im*a_im); // need to convert to fully-physical units to ensure this has the correct dimensions
                            /* this is where we can insert an ad-hoc renormalization to avoid double-counting if we have a genuinely very massive BH (so DF is well-resolved) */
                            {
                                double m_j=m_j_eff_for_df; /* estimate mean mass of the particles in the node */
                                if(sink_mass > 14.251*m_j) {fac_df *= DMIN(1.,DMAX(0.,(-1.+3./log10(sink_mass/m_j))/1.6));} /* approximate correction factor estimated by linhao */
                            }
                            if((m_j_eff_for_df <= MIN_REAL_NUMBER) || (b_impact <= MIN_REAL_NUMBER) || (dv2 <= MIN_REAL_NUMBER)) {fac_df = 0;}
                            /* parallel deflection component: dv = V[distant particle/node] - V[bh], sign here is set to accelerate towards V[ext], as needed */
                            acc += fac_df * dv_h;
                            /* perpendicular deflection component b_im = P[distant particle/node] - P[bh], so positive = accel -towards- P[ext], but this is the residual term (after subtracting the homogeneous term), which points in the opposite direction */
                            double fac_df_p = -fac_df / (b_impact * a_im + MIN_REAL_NUMBER);
                            if(fabs(fac_df_p)<MAX_REAL_NUMBER && isfinite(fac_df_p)) {acc += fac_df_p * b_im;}
                        }
                    }
#endif
                    
                    
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION /* these are the 'correction' terms for variable smoothing lengths (analogous to the ags-zeta terms above). need to adjust for variable ptypes using these */
                    int primary_uses_tidal_criterion=0, secondary_uses_tidal_criterion=0;
                    if(mass > 0 && r2 > 0)
                    {
                        if((1 << ptype) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)) {primary_uses_tidal_criterion=1;} /* check if the primary particle uses the tidal softening */
                        if((1 << ptype_sec) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION)) {secondary_uses_tidal_criterion=1;} /* check if the secondary particle uses the tidal softening */
                        double prefac_tt=0.5, h_touse=h, u_tt=sqrt(r2)/h_touse; // this corresponds to the result of symmetrizing by averaging
#if !defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
                        if(h >= h_p) {prefac_tt=1;} else {prefac_tt=1; h_touse=h_p; u_tt=sqrt(r2)/h_touse;} // this corresponds to adopting the MAX criterion for the softening
#endif
                        if(u_tt<1 && prefac_tt>0) {tidal_zeta += prefac_tt * mass * kernel_gravity(u_tt,1./h,1./(h*h*h),0);} // simple sum to calculate this contribution, only from particles inside the kernel of the primary -- this is up here instead of below the if below because it needs to include the 'self' contribution here
                    }
                    if(primary_uses_tidal_criterion || secondary_uses_tidal_criterion) // primary or secondary has associated correction terms here
                    { // now this is correct, but always need to carefully ensure correction terms are only applied in the correct 'direction' if we have a mixed-particle-type pair //
                        double h_touse = DMAX(h, h_p), f_b = -r*fac2_tidal, f_a = (6./r)*fac_tidal; // these will give the correct factors for the correction terms below, automatically symmetrized appropriately based on the same symmetry rules we use to define the tidal tensor itself in the first place
                        if(r < h_touse)
                        {
                            double dwk,wk,f_a_corr; u=r/h_touse; kernel_main(u,1.,1.,&wk,&dwk,0); // gather the remaining kernel terms (note these derivatives come from the laplacian and its derivative, so can be reconstructed from our usual wk and dwk terms //
                            f_a_corr = 4.*M_PI*mass * (dwk - (2./u)*wk) / (h_touse*h_touse*h_touse*h_touse); // default to symmetrize by taking the maximum, here
#if defined(ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING)
                            if(h<h_p) {h_touse=h;} else {h_touse=h_p;}
                            u=r/h_touse; kernel_main(u,1.,1.,&wk,&dwk,0);
                            f_a_corr = 0.5 * (f_a_corr + 4.*M_PI*mass * (dwk - (2./u)*wk) / pow(h_touse,4)); // symmetrize by averaging since thats what we did above
#endif
                            f_a += f_a_corr; // add this to the relevant function to use below
                        }
                        int ki,kj,kk; double acc_corr_zeta[3]={0}; Vec3<double> rh = dr / r;
                        for(kk=0;kk<3;kk++)
                        {
                            for(ki=0;ki<3;ki++)
                            {
                                for(kj=0;kj<3;kj++)
                                {
                                    double q0=rh[ki]*rh[kj]*rh[kk], fb_rh_add=0; /* first compute di dj dk [phi_kernel] -- this is generic */
                                    if(ki==kj) {fb_rh_add+=rh[kk];} /* these are the delta_ij terms */
                                    if(ki==kk) {fb_rh_add+=rh[kj];}
                                    if(kj==kk) {fb_rh_add+=rh[ki];}
                                    double qfun = f_a * q0 + f_b * (-3.*q0 + fb_rh_add); /* now double-dot this properly to get the sum we need - note only here need the combination of TT and zeta terms */
                                    acc_corr_zeta[kk] += primary_uses_tidal_criterion * i_zeta_tidal_tensorps_prevstep[ki][kj] * qfun; // only non-zero here if primary is a tidal-softening-active particle
                                    acc_corr_zeta[kk] += secondary_uses_tidal_criterion * j_zeta_tidal_tensorps_prevstep[ki][kj] * qfun; // only non-zero here if secondary is a tidal-softening-active particle
                                }
                            }
                        }
                        acc[0]+=acc_corr_zeta[0]; acc[1]+=acc_corr_zeta[1]; acc[2]+=acc_corr_zeta[2]; // final assignment
                    }
#endif
                    
                    
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                    /* tidal_tensorps[][] = Matrix of second derivatives of grav. potential, symmetric:
                     |Txx Txy Txz|   |tidal_tensorps[0][0] tidal_tensorps[0][1] tidal_tensorps[0][2]|
                     |Tyx Tyy Tyz| = |tidal_tensorps[1][0] tidal_tensorps[1][1] tidal_tensorps[1][2]|
                     |Tzx Tzy Tzz|   |tidal_tensorps[2][0] tidal_tensorps[2][1] tidal_tensorps[2][2]|  */
#ifdef GRAVITY_SPHERICAL_SYMMETRY
                    if(r_source < r_target) {fac2_tidal = 3 * mass / pow(DMAX(GRAVITY_SPHERICAL_SYMMETRY,DMAX(r_target,h)),5);} else {fac2_tidal = 0;}
#endif
#ifdef PMGRID
                    tidal_tensorps[0][0] += ((-fac_tidal + dr[0] * dr[0] * fac2_tidal) * shortrange_table[tabindex]) +
                    dr[0] * dr[0] * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                    tidal_tensorps[0][1] += ((dr[0] * dr[1] * fac2_tidal) * shortrange_table[tabindex]) +
                    dr[0] * dr[1] * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                    tidal_tensorps[0][2] += ((dr[0] * dr[2] * fac2_tidal) * shortrange_table[tabindex]) +
                    dr[0] * dr[2] * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                    tidal_tensorps[1][1] += ((-fac_tidal + dr[1] * dr[1] * fac2_tidal) * shortrange_table[tabindex]) +
                    dr[1] * dr[1] * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                    tidal_tensorps[1][2] += ((dr[1] * dr[2] * fac2_tidal) * shortrange_table[tabindex]) +
                    dr[1] * dr[2] * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                    tidal_tensorps[2][2] += ((-fac_tidal + dr[2] * dr[2] * fac2_tidal) * shortrange_table[tabindex]) +
                    dr[2] * dr[2] * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
#else
                    tidal_tensorps[0][0] += (-fac_tidal + dr[0] * dr[0] * fac2_tidal);
                    tidal_tensorps[0][1] += (dr[0] * dr[1] * fac2_tidal);
                    tidal_tensorps[0][2] += (dr[0] * dr[2] * fac2_tidal);
                    tidal_tensorps[1][1] += (-fac_tidal + dr[1] * dr[1] * fac2_tidal);
                    tidal_tensorps[1][2] += (dr[1] * dr[2] * fac2_tidal);
                    tidal_tensorps[2][2] += (-fac_tidal + dr[2] * dr[2] * fac2_tidal);
#endif
#endif // COMPUTE_TIDAL_TENSOR_IN_GRAVTREE //
#ifdef COMPUTE_JERK_IN_GRAVTREE
#ifndef ADAPTIVE_TREEFORCE_UPDATE // we want the jerk if we're doing lazy force updates
                    if(ptype > 0)
#endif
                    {
                        double dv_dot_dr = dot(dv, dr);
                        jerk += fac_accel * dv - dv_dot_dr * fac2_tidal * dr;
                    }
#endif
                } // closes TABINDEX<NTAB
                
                ninteractions++;
                
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
                if(r < r_for_total_menclosed) {m_enc_in_rcrit += mass;}
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
                tree_mass += mass;
#endif
#ifdef RT_USE_TREECOL_FOR_NH
                if(gasmass>0)
                {
                    int bin; // Here we do a simple six-bin angular binning scheme
                    if((fabs(dr[0]) > fabs(dr[1])) && (fabs(dr[0])>fabs(dr[2]))) {if (dr[0] > 0) {bin = 0;} else {bin=1;}
                    } else if (fabs(dr[1])>fabs(dr[2])){if (dr[1] > 0) {bin = 2;} else {bin=3;}
                    } else {if (dr[2] > 0) {bin = 4;} else {bin = 5;}}
                    treecol_angular_bins[bin] += fac_accel*gasmass*r / (angular_bin_size*mass); // in our binning scheme, we stretch the gas mass over a patch  of the sphere located at radius r subtending solid angle equal to the bin size - thus the area is r^2 * angular_bin_size, so sigma = m/(r^2 * angular bin size) = fac_accel/r / angular bin size. Factor of gasmass / mass corrects the gravitational mass to the gas mass
                }
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                if(ptype==0 && r>0 && cr_injection>0 && All.Time>All.TimeBegin)
                {
                    double kappa_0 = All.CosmicRay_Subgrid_Kappa_0, vst_0 = All.CosmicRay_Subgrid_Vstream_0; // in code units
                    double r_phys = sqrt(r*r + soft*soft/4.) * All.cf_atime, t_max = DMIN(1., evaluate_time_since_t_initial_in_Gyr(All.TimeBegin))/UNIT_TIME_IN_GYR; // make sure we're working in physical code units, and assign max time to formation at begin time, and include very crude 'softening' term here to prevent divergennce as r->0: for our default parameters can't be too large here or we get unphysically large CR halos compared to reality, b/c of large effective streaming terms
                    double r_max = 0.5*t_max*vst_0 * (1. + sqrt(1. + 16.*kappa_0/(vst_0*vst_0*t_max))); // maximum stream distance
#ifdef PMGRID
                    r_max = DMIN(r_max , 0.5*rcut*All.cf_atime); // truncate before reach the boundary of the grid to avoid numerical errors there
#endif
                    double fac_cr_distance = 1./(4.*M_PI*r_phys*(kappa_0 + vst_0*r_phys)) * exp(-DMIN(r_phys*r_phys/(1.e-6*r_phys*r_phys+r_max*r_max),50.));
                    if(fac_cr_distance>0) {SubGrid_CosmicRayEnergyDensity += fac_cr_distance * cr_injection / All.cf_a3inv;} // convert to appropriate code units for an energy density or pressure
                }
#endif
#ifdef RT_USE_GRAVTREE
                if(valid_gas_particle_for_rt)    /* we have a (valid) gas particle as target */
                {
                    r2 = d_stellarlum.norm_sq(); r = sqrt(r2); double fac_rt;
                    if(r >= soft) {fac_rt=1./(r2*r);} else {double h_inv_rt=1./soft, h3_inv_rt=h_inv_rt*h_inv_rt*h_inv_rt, u_rt=r*h_inv_rt; fac_rt=kernel_gravity(u_rt,h_inv_rt,h3_inv_rt,1);}
                    if((soft>r)&&(soft>0)) fac_rt *= (r2/(soft*soft)); // don't allow cross-section > r2
                    double fac_intensity; fac_intensity = fac_rt * r * All.cf_a2inv / (4.*M_PI); // ~L/(4pi*r^2), in -physical- units, since L is physical
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                    {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {Rad_E_gamma[kf] += fac_intensity * mass_stellarlum[kf];}}
#endif
#ifdef CHIMES_STELLAR_FLUXES
                    int chimes_k; double chimes_fac = fac_intensity / (UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS);  // 1/(4 * pi * r^2), in cm^-2
                    for (chimes_k = 0; chimes_k < CHIMES_LOCAL_UV_NBINS; chimes_k++)
                    {
                        chimes_flux_G0[chimes_k] += chimes_fac * chimes_mass_stellarlum_G0[chimes_k];   // Habing flux units
                        chimes_flux_ion[chimes_k] += chimes_fac * chimes_mass_stellarlum_ion[chimes_k]; // cm^-2 s^-1
                    }
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
                    incident_flux_uv += fac_intensity * mass_stellarlum[RT_FREQ_BIN_FIRE_UV];// * shortrange_table[tabindex];
                    if((mass_stellarlum[RT_FREQ_BIN_FIRE_IR]<mass_stellarlum[RT_FREQ_BIN_FIRE_UV])&&(mass_stellarlum[RT_FREQ_BIN_FIRE_IR]>0)) // if this -isn't- satisfied, no chance you are optically thin to EUV //
                    {
                        // here, use ratio and linear scaling of escape with tau to correct to the escape fraction for the correspondingly higher EUV kappa: factor ~2000 here comes from the ratio of (kappa_euv/kappa_uv)
                        incident_flux_euv += fac_intensity * mass_stellarlum[RT_FREQ_BIN_FIRE_UV] * (All.PhotonMomentum_fUV + (1-All.PhotonMomentum_fUV) *
                                                                                                     ((mass_stellarlum[RT_FREQ_BIN_FIRE_UV] + mass_stellarlum[RT_FREQ_BIN_FIRE_IR]) /
                                                                                                      (mass_stellarlum[RT_FREQ_BIN_FIRE_UV] + 2042.6*mass_stellarlum[RT_FREQ_BIN_FIRE_IR])));
                    } else {
                        // here, just enforce a minimum escape fraction //
                        double m_lum_total = 0; int ks_q; for(ks_q=0;ks_q<N_RT_FREQ_BINS;ks_q++) {m_lum_total += mass_stellarlum[ks_q];}
                        incident_flux_euv += All.PhotonMomentum_fUV * fac_intensity * m_lum_total;
                    }
                    // don't multiply by shortrange_table since that is to prevent 2x-counting by PMgrid (which never happens here) //
#endif
#ifdef SINK_PHOTONMOMENTUM
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                    Rad_E_gamma[RT_FREQ_BIN_FIRE_IR] += fac_intensity * mass_sinklumwt_forradfb;
#endif
#ifdef SINK_COMPTON_HEATING
                    incident_flux_agn += fac_intensity * mass_sinklumwt_forradfb; // L/(4pi*r*r) analog
#endif
#endif
                    
#ifdef RT_OTVET
                    /* use the information we have here from the gravity tree (optically thin incident fluxes) to estimate the Eddington tensor */
                    if(r>0)
                    {
                        double fac_otvet_sum=0; int kf_rt;
                        for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++)
                        {
                            fac_otvet_sum = mass_stellarlum[kf_rt];
                            fac_otvet_sum *= fac_rt / (1.e-37 + r); // units are not important, since ET will be dimensionless, but final ET should scale as ~luminosity/r^2
                            RT_ET[kf_rt] += fac_otvet_sum * outer_product(d_stellarlum);
                        }
                    }
                    
#endif
                    
#ifdef RT_LEBRON /* now we couple radiation pressure [single-scattering] terms within this module */
#ifdef GALSF_FB_FIRE_RT_LONGRANGE /* we only allow the momentum to couple over some distance to prevent bad approximations when the distance between points here is enormous */
                    if(r*UNIT_LENGTH_IN_KPC*All.cf_atime > 50.) {fac_rt=0;}
#endif
                    int kf_rt; double lum_force_fac=0;
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX) /* save the fluxes for use below, where we will calculate their RP normally */
                    double fac_flux = -fac_rt * All.cf_a2inv / (4.*M_PI); // ~L/(4pi*r^3), in -physical- units (except for last r, cancelled by dx_stellum), since L is physical
                    for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++) {Rad_Flux[kf_rt] += mass_stellarlum[kf_rt]*fac_flux*d_stellarlum;}
#else /* simply apply an on-the-spot approximation and do the absorption and RP force now */
                    for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++) {lum_force_fac += mass_stellarlum[kf_rt] * fac_stellum[kf_rt];} // add directly to forces. appropriate normalization (and sign) in 'fac_stellum'
#endif
#ifdef SINK_PHOTONMOMENTUM /* divide out PhotoMom_coupled_frac here b/c we have our own SINK_Rad_Mom factor, and don't want to double-count */
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
                    Rad_Flux[RT_FREQ_BIN_FIRE_IR] += mass_sinklumwt_forradfb*fac_flux*d_stellarlum;
#elif !defined(RT_DISABLE_RAD_PRESSURE)
                    lum_force_fac += (All.Sink_Rad_MomentumFactor / (MIN_REAL_NUMBER + All.PhotonMomentum_Coupled_Fraction)) * mass_sinklumwt_forradfb * fac_stellum[N_RT_FREQ_BINS-1];
#endif
#endif
                    if(lum_force_fac>0) {acc += (fac_rt*lum_force_fac) * d_stellarlum;}
#endif
                } // closes if(valid_gas_particle_for_rt)
                
#endif // RT_USE_GRAVTREE
                
                
#ifdef DM_SCALARFIELD_SCREENING
                if(ptype != 0)    /* we have a dark matter particle as target */
                {
                    GRAVITY_NEAREST_XYZ(d_dm[0],d_dm[1],d_dm[2],-1);
                    r2 = d_dm.norm_sq();
                    r = sqrt(r2); double fac_dmsf, h_inv_dmsf, h3inv_dmsf, u_dmsf;
                    if(r >= h) {fac_dmsf = mass_dm / (r2 * r);} else {
                        h_inv_dmsf = 1.0 / h; h3inv_dmsf = h_inv_dmsf * h_inv_dmsf * h_inv_dmsf; u_dmsf = r * h_inv_dmsf;
                        fac_dmsf = mass_dm * kernel_gravity(u_dmsf, h_inv_dmsf, h3inv_dmsf, 1);
                    }
                    /* assemble force with strength, screening length, and target charge.  */
                    fac_dmsf *= All.ScalarBeta * (1 + r / All.ScalarScreeningLength) * exp(-r / All.ScalarScreeningLength);
#ifdef PMGRID
                    tabindex = (int) (asmthfac * r);
                    if(tabindex < NTAB && tabindex >= 0)
#endif
                    {
#ifdef PMGRID
                        fac_dmsf *= shortrange_table[tabindex];
#endif
                        acc += fac_dmsf * d_dm;
                    }
                } // closes if(ptype != 0)
#endif // DM_SCALARFIELD_SCREENING //
                
            } // closes (if((r2 > 0) && (mass > 0))) check
            
            
            /* advance for used nodes: note this used to be above, now handled down here so we can use the 'no/nop' structures above */
            if(no < maxPart) {
                if(TakeLevel >= 0) {P[no].GravCost[TakeLevel] += 1.0;} /* node was used */
                no = Nextnode[no];
            } else {
                if(TakeLevel >= 0) {nop->GravCost += 1.0;}
                no = nop->u.d.sibling;
            }
            
        } // closes inner (while(no>=0)) check
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0)
                {
                    nodesinlist++;
                    no = Nodes[no].u.d.nextnode;    /* open it */
                }
            }
        } // closes (mode == 1) check
    } // closes outer (while(no>=0)) check
    
    
    /* store result at the proper place */
    if(mode == 0)
    {
        P[target].GravAccel = acc;
#ifdef RT_USE_TREECOL_FOR_NH
        int k; for(k=0; k < RT_USE_TREECOL_FOR_NH; k++) P[target].ColumnDensityBins[k] = treecol_angular_bins[k];
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
        P[target].TreeMass = tree_mass;
#endif
#ifdef RT_OTVET
        if(valid_gas_particle_for_rt) {int k; for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[target].ET[k] = RT_ET[k];}} else {if(P[target].Type==0) {int k; for(k=0;k<N_RT_FREQ_BINS;k++) {CellP[target].ET[k] = {};}}}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        if(valid_gas_particle_for_rt) {CellP[target].Rad_Flux_UV = incident_flux_uv;}
        if(valid_gas_particle_for_rt) {CellP[target].Rad_Flux_EUV = incident_flux_euv;}
#endif
#ifdef CHIMES_STELLAR_FLUXES
        if(valid_gas_particle_for_rt)
        {
            int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {CellP[target].Chimes_G0[kc] = chimes_flux_G0[kc]; CellP[target].Chimes_fluxPhotIon[kc] = chimes_flux_ion[kc];}
        }
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
        P[target].MencInRcrit = m_enc_in_rcrit;
#endif
#ifdef SINK_COMPTON_HEATING
        if(valid_gas_particle_for_rt) {CellP[target].Rad_Flux_AGN = incident_flux_agn;}
#endif
#if defined(COSMIC_RAY_SUBGRID_LEBRON)
        if(P[target].Type==0) {CellP[target].SubGrid_CosmicRayEnergyDensity = SubGrid_CosmicRayEnergyDensity;}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        if(valid_gas_particle_for_rt) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[target].Rad_E_gamma[kf] = Rad_E_gamma[kf];}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        if(valid_gas_particle_for_rt) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[target].Rad_Flux[kf] = Rad_Flux[kf];}}
#endif
#ifdef EVALPOTENTIAL
        P[target].Potential = pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        P[target].tidal_tensorps = tidal_tensorps;
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        P[target].tidal_zeta = tidal_zeta;
#endif
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        P[target].GravJerk = jerk;
#endif
#ifdef SINK_CALC_DISTANCES
        P[target].Min_Distance_to_Sink = sqrt( Min_Distance_to_Sink2 );
        P[target].Min_xyz_to_Sink = Min_xyz_to_Sink;   /* remember, dr = x_SINK - myx */
#ifdef SPECIAL_POINT_MOTION
        {
            P[target].vel_of_nearest_special = vel_of_nearest_special;
            P[target].acc_of_nearest_special = acc_of_nearest_special;
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
            P[target].weight_sum_for_special_point_smoothing = weight_sum_for_special_point_smoothing; /* weighted sum needed */
#endif
        }
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
        P[target].is_in_a_binary=0; P[target].Min_Sink_OrbitalTime=Min_Sink_OrbitalTime; //orbital time for binary
        if (Min_Sink_OrbitalTime<MAX_REAL_NUMBER)
        {
            P[target].is_in_a_binary=1; P[target].comp_Mass=comp_Mass; //mass of binary companion
            P[target].comp_dx = comp_dx; P[target].comp_dv = comp_dv;
        }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
        P[target].Min_Sink_Approach_Time = sqrt(Min_Sink_Approach_Time);
        P[target].Min_Sink_Freefall_time = sqrt(sqrt(Min_Sink_Freefall_time)/All.G);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        P[target].Min_Sink_FeedbackTime = sqrt(Min_Sink_FeedbackTime);
#endif
#endif
#endif // SINK_CALC_DISTANCES
    }
    else
    {
        GravDataResult[target].Acc = acc;
#ifdef COUNT_MASS_IN_GRAVTREE
        GravDataResult[target].TreeMass = tree_mass;
#endif
#ifdef RT_USE_TREECOL_FOR_NH
        {int k; for(k=0;k<RT_USE_TREECOL_FOR_NH;k++) GravDataResult[target].ColumnDensityBins[k] = treecol_angular_bins[k];}
#endif
#ifdef RT_OTVET
        {int k; for(k=0;k<N_RT_FREQ_BINS;k++) {GravDataResult[target].ET[k] = RT_ET[k];}}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        GravDataResult[target].Rad_Flux_UV = incident_flux_uv;
        GravDataResult[target].Rad_Flux_EUV = incident_flux_euv;
#endif
#ifdef CHIMES_STELLAR_FLUXES
        int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {GravDataResult[target].Chimes_G0[kc] = chimes_flux_G0[kc]; GravDataResult[target].Chimes_fluxPhotIon[kc] = chimes_flux_ion[kc];}
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
        GravDataResult[target].MencInRcrit = m_enc_in_rcrit;
#endif
#ifdef SINK_COMPTON_HEATING
        GravDataResult[target].Rad_Flux_AGN = incident_flux_agn;
#endif
#if defined(COSMIC_RAY_SUBGRID_LEBRON)
        GravDataResult[target].SubGrid_CosmicRayEnergyDensity = SubGrid_CosmicRayEnergyDensity;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {GravDataResult[target].Rad_E_gamma[kf] = Rad_E_gamma[kf];}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {GravDataResult[target].Rad_Flux[kf] = Rad_Flux[kf];}}
#endif
#ifdef EVALPOTENTIAL
        GravDataResult[target].Potential = pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        GravDataResult[target].tidal_tensorps = tidal_tensorps;
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        GravDataResult[target].tidal_zeta = tidal_zeta;
#endif
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        GravDataResult[target].GravJerk = jerk;
#endif
#ifdef SINK_CALC_DISTANCES
        GravDataResult[target].Min_Distance_to_Sink = sqrt( Min_Distance_to_Sink2 );
        GravDataResult[target].Min_xyz_to_Sink = Min_xyz_to_Sink;   /* remember, dr = x_SINK - myx */
#ifdef SPECIAL_POINT_MOTION
        {
            GravDataResult[target].vel_of_nearest_special = vel_of_nearest_special;
            GravDataResult[target].acc_of_nearest_special = acc_of_nearest_special;
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
            GravDataResult[target].weight_sum_for_special_point_smoothing = weight_sum_for_special_point_smoothing; /* weighted sum needed */
#endif
        }
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
        GravDataResult[target].is_in_a_binary=0; GravDataResult[target].Min_Sink_OrbitalTime=Min_Sink_OrbitalTime; // orbital time for binary
        if (Min_Sink_OrbitalTime<MAX_REAL_NUMBER)
        {
            GravDataResult[target].is_in_a_binary = 1; GravDataResult[target].comp_Mass=comp_Mass; //mass of binary companion
            GravDataResult[target].comp_dx = comp_dx; GravDataResult[target].comp_dv = comp_dv;
        }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
        GravDataResult[target].Min_Sink_Approach_Time = sqrt(Min_Sink_Approach_Time);
        GravDataResult[target].Min_Sink_Freefall_time = sqrt(sqrt(Min_Sink_Freefall_time)/All.G);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        GravDataResult[target].Min_Sink_FeedbackTime = sqrt(Min_Sink_FeedbackTime);
#endif
#endif
#endif // SINK_CALC_DISTANCES
        *exportflag = nodesinlist;
    }
    
    return ninteractions;
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
int force_treeevaluate_ewald_correction(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex)
{
    struct NODE *nop = 0;
    int signx, signy, signz, nexp, i, j, k, openflag, task, no, cost, listindex = 0;
    double mass, r2, u, v, w, f1, f2, f3, f4, f5, f6, f7, f8;
    double boxsize, boxhalf, aold, xtmp; xtmp=0;
    Vec3<double> pos, dr;
    Vec3<MyDouble> acc = {};

    boxsize = All.BoxSize;
    boxhalf = 0.5 * All.BoxSize;

    cost = 0;
    if(mode == 0)
    {
        pos = P[target].Pos;
        aold = All.ErrTolForceAcc * P[target].OldAcc;
    }
    else
    {
        pos = GravDataGet[target].Pos;
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
    }
    
    if(mode == 0)
    {
        no = All.MaxPart;        /* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;    /* open it */
    }
    
    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)    /* single particle */
            {
                /* the index of the node is the index of the particle */
                /* observe the sign */
                if(P[no].Ti_current != All.Ti_Current)
                {
#ifdef _OPENMP
#pragma omp critical(_particledriftewald_)
#endif
                    {
                        if(P[no].Ti_current != All.Ti_Current) {drift_particle(no, All.Ti_Current);}
                    }
                }

                dr = P[no].Pos - pos;
                mass = P[no].Mass;
            }
            else            /* we have an  internal node */
            {
                if(no >= All.MaxPart + MaxNodes)    /* pseudo particle */
                {
                    if(mode == 0)
                    {
                        if(exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
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
                                    /* out if buffer space. Need to discard work for this particle and interrupt */
                                    BufferFullFlag = 1;
                                    exitFlag = 1;
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
                        
                        DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] = DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                        
                        if(exportnodecount[task] < NODELISTLENGTH) {DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;}
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }
                
                nop = &Nodes[no];
                
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }
                //if(nop->N_part <= 1) /* open cell */
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    no = nop->u.d.nextnode;
                    continue;
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

            if(no < All.MaxPart)
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
                    /* now we check if we can avoid opening the cell */

                    u = nop->center[0] - pos[0];
                    if(u > boxhalf) {u -= boxsize;}
                    if(u < -boxhalf) {u += boxsize;}
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

                    u = nop->center[1] - pos[1];
                    if(u > boxhalf) {u -= boxsize;}
                    if(u < -boxhalf) {u += boxsize;}
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

                    u = nop->center[2] - pos[2];
                    if(u > boxhalf) {u -= boxsize;}
                    if(u < -boxhalf) {u += boxsize;}
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                    
                    /* if the cell is too large, we need to refine it further */
                    if(nop->len > 0.20 * boxsize)
                    {
                        /* cell is too large */
                        no = nop->u.d.nextnode;
                        continue;
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
            u = dr[0] * fac_intp;
            i = (int) u;
            if(i >= EN) {i = EN - 1;}
            u -= i;
            v = dr[1] * fac_intp;
            j = (int) v;
            if(j >= EN) {j = EN - 1;}
            v -= j;
            w = dr[2] * fac_intp;
            k = (int) w;
            if(k >= EN) {k = EN - 1;}
            w -= k;
            /* compute factors for trilinear interpolation */
            f1 = (1 - u) * (1 - v) * (1 - w);
            f2 = (1 - u) * (1 - v) * (w);
            f3 = (1 - u) * (v) * (1 - w);
            f4 = (1 - u) * (v) * (w);
            f5 = (u) * (1 - v) * (1 - w);
            f6 = (u) * (1 - v) * (w);
            f7 = (u) * (v) * (1 - w);
            f8 = (u) * (v) * (w);
            acc[0] += (mass * signx * (fcorrx[i][j][k] * f1 +
                                      fcorrx[i][j][k + 1] * f2 +
                                      fcorrx[i][j + 1][k] * f3 +
                                      fcorrx[i][j + 1][k + 1] * f4 +
                                      fcorrx[i + 1][j][k] * f5 +
                                      fcorrx[i + 1][j][k + 1] * f6 +
                                      fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8));
            acc[1] +=
            (mass * signy *
             (fcorry[i][j][k] * f1 + fcorry[i][j][k + 1] * f2 +
              fcorry[i][j + 1][k] * f3 + fcorry[i][j + 1][k + 1] * f4 + fcorry[i + 1]
              [j][k] * f5 + fcorry[i + 1][j][k + 1] * f6 + fcorry[i + 1][j + 1][k] *
              f7 + fcorry[i + 1][j + 1][k + 1] * f8));
            acc[2] +=
            (mass * signz *
             (fcorrz[i][j][k] * f1 + fcorrz[i][j][k + 1] * f2 +
              fcorrz[i][j + 1][k] * f3 + fcorrz[i][j + 1][k + 1] * f4 + fcorrz[i + 1]
              [j][k] * f5 + fcorrz[i + 1][j][k + 1] * f6 + fcorrz[i + 1][j + 1][k] *
              f7 + fcorrz[i + 1][j + 1][k + 1] * f8));
            cost++;
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0) {no = Nodes[no].u.d.nextnode;}    /* open it */
            }
        }
    }
    
    /* add the result at the proper place */
    
    if(mode == 0)
    {
        P[target].GravAccel += acc;
    }
    else
    {
        GravDataResult[target].Acc = acc;
    }

    return cost;
}
#endif // #ifdef BOX_PERIODIC //




/*! This routine computes the gravitational potential by walking the
 *  tree. The same opening criteria is used as for the gravitational force
 *  walk.
 */
/*! This function also computes the short-range potential when the TreePM
 *  algorithm is used. This potential is the Newtonian potential, modified
 *  by a complementary error function.
 */
int force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
{
    struct NODE *nop = 0;
    MyDouble pot;
    int no, ptype, task, nexport_save, listindex = 0;
    double r2, dx, dy, dz, mass, r, u, h, h_inv, pos_x, pos_y, pos_z, aold, fac_pot, dxx, dyy, dzz, soft = 0;
#ifdef PMGRID
    int tabindex;
    double eff_dist, rcut, asmth, asmthfac;
#endif
    
    nexport_save = *nexport;
    pot = 0;
#ifdef PMGRID
    rcut = All.Rcut[0];
    asmth = All.Asmth[0];
#endif
    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        ptype = P[target].Type;
        aold = All.ErrTolForceAcc * P[target].OldAcc;
        soft = ForceSoftening_KernelRadius(target);
#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
        if(pmforce_is_particle_high_res(ptype, P[target].Pos))
        {
            rcut = All.Rcut[1];
            asmth = All.Asmth[1];
        }
#endif
    }
    else
    {
        pos_x = GravDataGet[target].Pos[0];
        pos_y = GravDataGet[target].Pos[1];
        pos_z = GravDataGet[target].Pos[2];
        ptype = GravDataGet[target].Type;
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
        soft = GravDataGet[target].Soft;
#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
        if(pmforce_is_particle_high_res(ptype, GravDataGet[target].Pos))
        {
            rcut = All.Rcut[1];
            asmth = All.Asmth[1];
        }
#endif
    }
    
#ifdef PMGRID
    asmthfac = 0.5 / asmth * (NTAB / 3.0);
#endif
    if(mode == 0)
    {
        no = All.MaxPart;        /* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;    /* open it */
    }
    
    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)    /* single particle */
            {
                /* the index of the node is the index of the particle */
                /* observe the sign  */
                if(P[no].Ti_current != All.Ti_Current) {drift_particle(no, All.Ti_Current);}
                dx = P[no].Pos[0] - pos_x;
                dy = P[no].Pos[1] - pos_y;
                dz = P[no].Pos[2] - pos_z;
                mass = P[no].Mass;
            }
            else
            {
                if(no >= All.MaxPart + MaxNodes)    /* pseudo particle */
                {
                    if(mode == 0)
                    {
                        if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                        {
                            Exportflag[task] = target;
                            Exportnodecount[task] = NODELISTLENGTH;
                        }
                        
                        if(Exportnodecount[task] == NODELISTLENGTH)
                        {
                            if(*nexport >= All.BunchSize)
                            {
                                *nexport = nexport_save;
                                if(nexport_save == 0) {endrun(13002);} /* in this case, the buffer is too small to process even a single particle */
                                for(task = 0; task < NTask; task++) {nsend_local[task] = 0;}
                                for(no = 0; no < nexport_save; no++) {nsend_local[DataIndexTable[no].Task]++;}
                                return -1; /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */
                            }
                            Exportnodecount[task] = 0;
                            Exportindex[task] = *nexport;
                            DataIndexTable[*nexport].Task = task;
                            DataIndexTable[*nexport].Index = target;
                            DataIndexTable[*nexport].IndexGet = *nexport;
                            *nexport = *nexport + 1;
                            nsend_local[task]++;
                        }
                        
                        DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                        if(Exportnodecount[task] < NODELISTLENGTH) {DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;}
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }
                
                nop = &Nodes[no];
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }
                //if(nop->N_part <= 1) /* open cell */
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    no = nop->u.d.nextnode;
                    continue;
                }
                if(nop->Ti_current != All.Ti_Current) {force_drift_node(no, All.Ti_Current);}
                mass = nop->u.d.mass;
                dx = nop->u.d.s[0] - pos_x;
                dy = nop->u.d.s[1] - pos_y;
                dz = nop->u.d.s[2] - pos_z;
            }
            GRAVITY_NEAREST_XYZ(dx,dy,dz,-1);
            r2 = dx * dx + dy * dy + dz * dz;
            if(no < All.MaxPart)
            {
                h = soft; /* set softening */
                no = Nextnode[no];
            }
            else            /* we have an internal node. Need to check opening criterion */
            {
#ifdef PMGRID
                /* check whether we can stop walking along this branch */
                if(no >= All.MaxPart + MaxNodes)    /* pseudo particle */
                {
                    if(mode == 0)
                    {
                        if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                        {
                            Exportflag[task] = target;
                            DataIndexTable[*nexport].Index = target;
                            DataIndexTable[*nexport].Task = task;    /* Destination task */
                            *nexport = *nexport + 1;
                            nsend_local[task]++;
                        }
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }
                
                eff_dist = rcut + 0.5 * nop->len;
                dxx = nop->center[0] - pos_x;    /* observe the sign ! */
                dyy = nop->center[1] - pos_y;    /* this vector is -y in my thesis notation */
                dzz = nop->center[2] - pos_z;
                GRAVITY_NEAREST_XYZ(dxx,dyy,dzz,-1);
                if(dxx < -eff_dist || dxx > eff_dist)
                {
                    no = nop->u.d.sibling;
                    continue;
                }
                
                if(dyy < -eff_dist || dyy > eff_dist)
                {
                    no = nop->u.d.sibling;
                    continue;
                }
                
                if(dzz < -eff_dist || dzz > eff_dist)
                {
                    no = nop->u.d.sibling;
                    continue;
                }
#else // PMGRID
                dxx = nop->center[0] - pos_x;    /* observe the sign ! */
                dyy = nop->center[1] - pos_y;    /* this vector is -y in my thesis notation */
                dzz = nop->center[2] - pos_z;
                GRAVITY_NEAREST_XYZ(dxx,dyy,dzz,-1);
#endif // PMGRID
                
                if(All.ErrTolTheta)    /* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#ifndef GRAVITY_HYBRID_OPENING_CRIT
                else        /* check relative opening criterion */
#else
                    if(!(All.Ti_Current == 0 && RestartFlag != 1))
#endif
                    {
                        
                        /* force node to open if we are within the gravitational softening length */
                        if((r2 < (soft+0.6*nop->len)*(soft+0.6*nop->len)) || (r2 < (nop->maxsoft+0.6*nop->len)*(nop->maxsoft+0.6*nop->len)))
                        {
                            no = nop->u.d.nextnode;
                            continue;
                        }
                        
                        if(mass * nop->len * nop->len > r2 * r2 * aold)
                        {
                            /* open cell */
                            no = nop->u.d.nextnode;
                            continue;
                        }
                        
                        if(fabs(dxx) < 0.60 * nop->len)
                        {
                            if(fabs(dyy) < 0.60 * nop->len)
                            {
                                if(fabs(dzz) < 0.60 * nop->len)
                                {
                                    no = nop->u.d.nextnode;
                                    continue;
                                }
                            }
                        }
                    }
                
                h = soft; // set h if not already set above
                if(h < nop->maxsoft) {if(r2 < nop->maxsoft * nop->maxsoft) {no = nop->u.d.nextnode; continue;}}
                no = nop->u.d.sibling;    /* node can be used */
            }
            
            r = sqrt(r2);
#ifdef PMGRID
            tabindex = (int) (r * asmthfac);
            if(tabindex < NTAB && tabindex >= 0)
#endif
            {
#ifdef PMGRID
                fac_pot = shortrange_table_potential[tabindex];
#else
                fac_pot = 1;
#endif
                if(r >= h) {pot += (-fac_pot * mass / r);} else {
                    h_inv = 1.0 / h; u = r * h_inv;
                    pot += ( fac_pot * mass * kernel_gravity(u, h_inv, 1, -1) );
                }
            }
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
            pot += (mass * ewald_pot_corr(dx, dy, dz));
#endif
        }
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0) {no = Nodes[no].u.d.nextnode;}    /* open it */
            }
        }
    }
    
    /* store result at the proper place */
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL)
    if(mode == 0) {P[target].Potential = pot;} else {PotDataResult[target].Potential = pot;}
#endif
    return 0;
}





#ifdef SUBFIND
int subfind_force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
{
    struct NODE *nop = 0;
    MyDouble pot;
    int no, ptype, task, nexport_save, listindex = 0;
    double r2, dx, dy, dz, mass, r, u, h, h_inv, pos_x, pos_y, pos_z, soft=0;
    
    nexport_save = *nexport;
    pot = 0;
    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        ptype = P[target].Type;
        soft  = ForceSoftening_KernelRadius(target);
    }
    else
    {
        pos_x = GravDataGet[target].Pos[0];
        pos_y = GravDataGet[target].Pos[1];
        pos_z = GravDataGet[target].Pos[2];
        ptype = GravDataGet[target].Type;
        soft  = GravDataGet[target].Soft;
    }
    
    h = soft; h_inv = 1.0 / h;
    
    if(mode == 0)
    {
        no = All.MaxPart;        /* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;    /* open it */
    }
    
    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)    /* single particle */
            {
                /* the index of the node is the index of the particle */
                /* observe the sign */
                
                dx = P[no].Pos[0] - pos_x;
                dy = P[no].Pos[1] - pos_y;
                dz = P[no].Pos[2] - pos_z;
                mass = P[no].Mass;
            }
            else
            {
                if(no >= All.MaxPart + MaxNodes)    /* pseudo particle */
                {
                    if(mode == 0)
                    {
                        if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                        {
                            Exportflag[task] = target;
                            Exportnodecount[task] = NODELISTLENGTH;
                        }
                        
                        if(Exportnodecount[task] == NODELISTLENGTH)
                        {
                            if(*nexport >= All.BunchSize)
                            {
                                *nexport = nexport_save;
                                if(nexport_save == 0) {endrun(13001);} /* in this case, the buffer is too small to process even a single particle */
                                for(task = 0; task < NTask; task++) {nsend_local[task] = 0;}
                                for(no = 0; no < nexport_save; no++) {nsend_local[DataIndexTable[no].Task]++;}
                                return -1; /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */
                            }
                            Exportnodecount[task] = 0;
                            Exportindex[task] = *nexport;
                            DataIndexTable[*nexport].Task = task;
                            DataIndexTable[*nexport].Index = target;
                            DataIndexTable[*nexport].IndexGet = *nexport;
                            *nexport = *nexport + 1;
                            nsend_local[task]++;
                        }
                        
                        DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                        if(Exportnodecount[task] < NODELISTLENGTH) {DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;}
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }
                
                nop = &Nodes[no];
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }
                
                mass = nop->u.d.mass;
                //if(nop->N_part <= 1)
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    if(mass) /* open cell */
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
                
                dx = nop->u.d.s[0] - pos_x;
                dy = nop->u.d.s[1] - pos_y;
                dz = nop->u.d.s[2] - pos_z;
            }
            GRAVITY_NEAREST_XYZ(dx,dy,dz,-1);
            r2 = dx * dx + dy * dy + dz * dz;
            if(no < All.MaxPart)
            {
                no = Nextnode[no];
            }
            else            /* we have an internal node. Need to check opening criterion */
            {
                /* check Barnes-Hut opening criterion */
                double ErrTolThetaSubfind = All.ErrTolTheta;
                if(nop->len * nop->len > r2 * ErrTolThetaSubfind * ErrTolThetaSubfind)
                {
                    /* open cell */
                    if(mass)
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
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
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0) {no = Nodes[no].u.d.nextnode;}    /* open it */
            }
        }
    }
    
    /* store result at the proper place */
    
    if(mode == 0)
        P[target].u.DM_Potential = pot;
    else
        PotDataResult[target].Potential = pot;
    return 0;
}
#endif // SUBFIND //




/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually,
 *  maxnodes approximately equal to 0.7*maxpart is sufficient to store the
 *  tree for up to maxpart particles.
 */
void force_treeallocate(int maxnodes, int maxpart)
{
    int i;
    size_t bytes;
    double allbytes = 0, allbytes_topleaves = 0;
    double u;
    
    tree_allocated_flag = 1;
    DomainNodeIndex = (int *) mymalloc("DomainNodeIndex", bytes = NTopleaves * sizeof(int));
    allbytes_topleaves += bytes;
    MaxNodes = maxnodes;
    if(!(Nodes_base = (struct NODE *) mymalloc("Nodes_base", bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
        printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
        endrun(3);
    }
    allbytes += bytes;
    if(!
       (Extnodes_base =
        (struct extNODE *) mymalloc("Extnodes_base", bytes = (MaxNodes + 1) * sizeof(struct extNODE))))
    {
        printf("failed to allocate memory for %d tree-extnodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
        endrun(3);
    }
    allbytes += bytes;
    Nodes = Nodes_base - All.MaxPart;
    Extnodes = Extnodes_base - All.MaxPart;
    if(!(Nextnode = (int *) mymalloc("Nextnode", bytes = (maxpart + NTopnodes) * sizeof(int))))
    {
        printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n",
               maxpart + NTopnodes, bytes / (1024.0 * 1024.0));
        endrun(8267342);
    }
    allbytes += bytes;
    if(!(Father = (int *) mymalloc("Father", bytes = (maxpart) * sizeof(int))))
    {
        printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
        endrun(438965237);
    }
    allbytes += bytes;
    if(first_flag == 0)
    {
        first_flag = 1;
        if(ThisTask == 0)
            printf
            ("Allocated %g MByte for tree, and %g Mbyte for top-leaves.  (presently allocated %g MB)\n",
             allbytes / (1024.0 * 1024.0), allbytes_topleaves / (1024.0 * 1024.0),
             AllocatedBytes / (1024.0 * 1024.0));
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


/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
    if(tree_allocated_flag)
    {
        myfree(Father);
        myfree(Nextnode);
        myfree(Extnodes_base);
        myfree(Nodes_base);
        myfree(DomainNodeIndex);
        tree_allocated_flag = 0;
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
    int i, j, k;
    double u, v, w;
    double f1, f2, f3, f4, f5, f6, f7, f8;
    
    if(dx < 0)
        dx = -dx;
    if(dy < 0)
        dy = -dy;
    if(dz < 0)
        dz = -dz;
    u = dx * fac_intp;
    i = (int) u;
    if(i >= EN)
        i = EN - 1;
    u -= i;
    v = dy * fac_intp;
    j = (int) v;
    if(j >= EN)
        j = EN - 1;
    v -= j;
    w = dz * fac_intp;
    k = (int) w;
    if(k >= EN)
        k = EN - 1;
    w -= k;
    f1 = (1 - u) * (1 - v) * (1 - w);
    f2 = (1 - u) * (1 - v) * (w);
    f3 = (1 - u) * (v) * (1 - w);
    f4 = (1 - u) * (v) * (w);
    f5 = (u) * (1 - v) * (1 - w);
    f6 = (u) * (1 - v) * (w);
    f7 = (u) * (v) * (1 - w);
    f8 = (u) * (v) * (w);
    return potcorr[i][j][k] * f1 +
    potcorr[i][j][k + 1] * f2 +
    potcorr[i][j + 1][k] * f3 +
    potcorr[i][j + 1][k + 1] * f4 +
    potcorr[i + 1][j][k] * f5 +
    potcorr[i + 1][j][k + 1] * f6 + potcorr[i + 1][j + 1][k] * f7 + potcorr[i + 1][j + 1][k + 1] * f8;
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
    TreeOpsCount[TREEOPS_REFRESH]++;
    int i, k, no;
    PRINT_STATUS("Refreshing tree node moments (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));

    /* Step 1: zero all node moment fields (preserving structural fields: nextnode, sibling, father) */
    for(no = All.MaxPart; no < All.MaxPart + Numnodestree; no++)
    {
        unsigned int saved_bitflags = Nodes[no].u.d.bitflags; /* preserve topology flags (TOPLEVEL etc) */
        Nodes[no].u.d.mass = 0;
        Nodes[no].u.d.s = {};
        Nodes[no].GravCost = 0;
        Nodes[no].Ti_current = All.Ti_Current;
        Nodes[no].N_part = 0;
        Nodes[no].maxsoft = 0;
        Nodes[no].u.d.bitflags = saved_bitflags & ((1 << BITFLAG_TOPLEVEL) | (1 << BITFLAG_DEPENDS_ON_LOCAL_ELEMENT) | (1 << BITFLAG_INTERNAL_TOPLEVEL));
        Extnodes[no].vs = {};
        Extnodes[no].hmax = 0;
        Extnodes[no].vmax = 0;
        Extnodes[no].divVmax = 0;
        Extnodes[no].dp = {};
#ifdef SINK_NODE_MOTION_TRACKED
        Extnodes[no].sink_dp = {};
#endif
        Extnodes[no].Ti_lastkicked = All.Ti_Current;
        Extnodes[no].Flag = GlobFlag;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        Nodes[no].gasmass = 0;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        Nodes[no].cr_injection = 0;
#endif
#ifdef RT_USE_GRAVTREE
        for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k]=0;}
#ifdef CHIMES_STELLAR_FLUXES
        for(k=0;k<CHIMES_LOCAL_UV_NBINS;k++) {Nodes[no].chimes_stellar_lum_G0[k]=0; Nodes[no].chimes_stellar_lum_ion[k]=0;}
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Nodes[no].rt_source_lum_s = {};
        Extnodes[no].rt_source_lum_vs = {};
        Extnodes[no].rt_source_lum_dp = {};
#endif
#ifdef SINK_PHOTONMOMENTUM
        Nodes[no].sink_lum = 0;
        Nodes[no].sink_lum_grad = {};
#endif
#ifdef SINK_CALC_DISTANCES
        Nodes[no].sink_mass = 0;
        Nodes[no].sink_pos = {};
#ifdef SINK_NODE_MOTION_TRACKED
        Nodes[no].sink_vel = {};
        Nodes[no].N_SINK = 0;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        Nodes[no].MaxFeedbackVel = 0;
#endif
#endif
#ifdef SPECIAL_POINT_MOTION
        Nodes[no].sink_acc = {};
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        Nodes[no].tidal_tensorps_prevstep = {};
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Nodes[no].mass_dm = 0;
        Nodes[no].s_dm = {};
        Extnodes[no].vs_dm = {};
        Extnodes[no].dp_dm = {};
#endif
    }

    /* Step 2: accumulate particle contributions into their immediate parent nodes.
       Store s as mass-weighted position sum (NOT center of mass yet) for propagation. */
    for(i = 0; i < NumPart; i++)
    {
        no = Father[i];
        if(no < 0) {continue;}
        struct particle_data *pa = &P[i];

        Nodes[no].u.d.mass += pa->Mass;
        Nodes[no].u.d.s += pa->Mass * pa->Pos;
        Extnodes[no].vs += pa->Mass * pa->Vel;
        Nodes[no].N_part++;

        MyFloat v, vmax_p = 0;
        for(k = 0; k < 3; k++) {if((v = fabs(pa->Vel[k])) > vmax_p) {vmax_p = v;}}
        if(vmax_p > Extnodes[no].vmax) {Extnodes[no].vmax = vmax_p;}

        double soft_p = ForceSoftening_KernelRadius(i);
        if(soft_p > Nodes[no].maxsoft) {Nodes[no].maxsoft = soft_p;}
#ifdef SINGLE_STAR_SINK_DYNAMICS
        if(pa->Type == 5) {if(P[i].KernelRadius > Nodes[no].maxsoft) {Nodes[no].maxsoft = P[i].KernelRadius;}}
#endif

        if(pa->Type == 0)
        {
            double htmp = DMIN(All.MaxKernelRadius, P[i].KernelRadius);
            if(htmp > Extnodes[no].hmax) {Extnodes[no].hmax = htmp;}
            if(P[i].Particle_DivVel > Extnodes[no].divVmax) {Extnodes[no].divVmax = P[i].Particle_DivVel;}
        }

#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        if(pa->Type == 0) {Nodes[no].gasmass += pa->Mass;}
#if defined(SINK_ALPHADISK_ACCRETION) && defined(RT_USE_TREECOL_FOR_NH)
        if(pa->Type == 5) {Nodes[no].gasmass += P[i].Sink_Mass_Reservoir;}
#endif
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        Nodes[no].cr_injection += cr_get_source_injection_rate(i, P, CellP);
#endif
#ifdef RT_USE_GRAVTREE
        {double lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
        double chimes_lum_G0[CHIMES_LOCAL_UV_NBINS], chimes_lum_ion[CHIMES_LOCAL_UV_NBINS];
        int active_check = rt_get_source_luminosity_chimes(i,1,lum,chimes_lum_G0,chimes_lum_ion, P, CellP);
#else
        int active_check = rt_get_source_luminosity(i,1,lum, P, CellP);
#endif
        if(active_check) {
            double l_sum = 0;
            for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k] += lum[k]; l_sum += lum[k];}
#ifdef CHIMES_STELLAR_FLUXES
            for(k=0;k<CHIMES_LOCAL_UV_NBINS;k++) {Nodes[no].chimes_stellar_lum_G0[k] += chimes_lum_G0[k]; Nodes[no].chimes_stellar_lum_ion[k] += chimes_lum_ion[k];}
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            Nodes[no].rt_source_lum_s += l_sum * pa->Pos;
            Extnodes[no].rt_source_lum_vs += l_sum * pa->Vel;
#endif
        }}
#endif
#ifdef SINK_PHOTONMOMENTUM
        if(pa->Type == 5 && pa->Mass > 0 && pa->DensityAroundParticle > 0 && pa->Sink_Mdot > 0) {
            double BHLum = sink_lum_bol(pa->Sink_Mdot, pa->Sink_Mass, i);
            Nodes[no].sink_lum += BHLum;
#if defined(SINK_FOLLOW_ACCRETED_ANGMOM)
            Nodes[no].sink_lum_grad += pa->Sink_Specific_AngMom * BHLum;
#else
            Nodes[no].sink_lum_grad += pa->GradRho * BHLum;
#endif
        }
#endif
#ifdef SINK_CALC_DISTANCES
        if(pa->Type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES) {
            Nodes[no].sink_mass += pa->Mass;
            Nodes[no].sink_pos += pa->Mass * pa->Pos; /* store as mass-weighted sum, normalize later */
#ifdef SINK_NODE_MOTION_TRACKED
            Nodes[no].N_SINK += 1;
#endif
#ifdef SINK_NODE_MOTION_TRACKED
            Nodes[no].sink_vel += pa->Mass * pa->Vel; /* mass-weighted, normalize later */
#endif
#ifdef SPECIAL_POINT_MOTION
            Nodes[no].sink_acc += pa->Mass * pa->Acc_Total_PrevStep;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
            if(pa->MaxFeedbackVel > Nodes[no].MaxFeedbackVel) {Nodes[no].MaxFeedbackVel = pa->MaxFeedbackVel;}
#endif
        }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        {for(k=0;k<6;k++) {Nodes[no].tidal_tensorps_prevstep.data[k] += pa->Mass * pa->tidal_tensorps_prevstep.data[k];}}
#endif
#ifdef DM_SCALARFIELD_SCREENING
        if(pa->Type != 0) {Nodes[no].mass_dm += pa->Mass; Nodes[no].s_dm += pa->Mass * pa->Pos; Extnodes[no].vs_dm += pa->Mass * pa->Vel;}
#endif
    }

    /* Step 3: propagate node moments bottom-up. Children have higher indices than parents,
       so iterating in reverse order processes children before parents. At this stage s, vs,
       sink_pos, sink_vel, etc. are stored as mass-weighted sums (not yet normalized). */
    for(no = All.MaxPart + Numnodestree - 1; no >= All.MaxPart; no--)
    {
        int father = Nodes[no].u.d.father;
        if(father < All.MaxPart || father >= All.MaxPart + Numnodestree) {continue;} /* root or invalid */

        Nodes[father].u.d.mass += Nodes[no].u.d.mass;
        Nodes[father].u.d.s += Nodes[no].u.d.s; /* still mass-weighted position sum */
        Extnodes[father].vs += Extnodes[no].vs;
        Nodes[father].N_part += Nodes[no].N_part;
        if(Extnodes[no].hmax > Extnodes[father].hmax) {Extnodes[father].hmax = Extnodes[no].hmax;}
        if(Extnodes[no].vmax > Extnodes[father].vmax) {Extnodes[father].vmax = Extnodes[no].vmax;}
        if(Extnodes[no].divVmax > Extnodes[father].divVmax) {Extnodes[father].divVmax = Extnodes[no].divVmax;}
        if(Nodes[no].maxsoft > Nodes[father].maxsoft) {Nodes[father].maxsoft = Nodes[no].maxsoft;}
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        Nodes[father].gasmass += Nodes[no].gasmass;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        Nodes[father].cr_injection += Nodes[no].cr_injection;
#endif
#ifdef RT_USE_GRAVTREE
        for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[father].stellar_lum[k] += Nodes[no].stellar_lum[k];}
#ifdef CHIMES_STELLAR_FLUXES
        for(k=0;k<CHIMES_LOCAL_UV_NBINS;k++) {Nodes[father].chimes_stellar_lum_G0[k] += Nodes[no].chimes_stellar_lum_G0[k]; Nodes[father].chimes_stellar_lum_ion[k] += Nodes[no].chimes_stellar_lum_ion[k];}
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Nodes[father].rt_source_lum_s += Nodes[no].rt_source_lum_s;
        Extnodes[father].rt_source_lum_vs += Extnodes[no].rt_source_lum_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
        Nodes[father].sink_lum += Nodes[no].sink_lum;
        Nodes[father].sink_lum_grad += Nodes[no].sink_lum_grad; /* still lum-weighted sum */
#endif
#ifdef SINK_CALC_DISTANCES
        Nodes[father].sink_mass += Nodes[no].sink_mass;
        Nodes[father].sink_pos += Nodes[no].sink_pos; /* propagate mass-weighted sum */
#ifdef SINK_NODE_MOTION_TRACKED
        Nodes[father].N_SINK += Nodes[no].N_SINK;
        Nodes[father].sink_vel += Nodes[no].sink_vel;
#endif
#ifdef SPECIAL_POINT_MOTION
        Nodes[father].sink_acc += Nodes[no].sink_acc;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        if(Nodes[no].sink_mass > 0 && Nodes[no].MaxFeedbackVel > Nodes[father].MaxFeedbackVel) {Nodes[father].MaxFeedbackVel = Nodes[no].MaxFeedbackVel;}
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        {for(k=0;k<6;k++) {Nodes[father].tidal_tensorps_prevstep.data[k] += Nodes[no].tidal_tensorps_prevstep.data[k];}}
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Nodes[father].mass_dm += Nodes[no].mass_dm;
        Nodes[father].s_dm += Nodes[no].s_dm;
        Extnodes[father].vs_dm += Extnodes[no].vs_dm;
#endif
    }

    /* Step 4: normalize mass-weighted sums to get actual COM, velocities, etc. */
    for(no = All.MaxPart; no < All.MaxPart + Numnodestree; no++)
    {
        MyFloat mass = Nodes[no].u.d.mass;
        if(mass > 0) {
            Nodes[no].u.d.s /= mass;
            Extnodes[no].vs /= mass;
        } else {
            Nodes[no].u.d.s = Nodes[no].center;
            Extnodes[no].vs = {};
        }
        if(Nodes[no].N_part > 1) {Nodes[no].u.d.bitflags |= (1 << BITFLAG_MULTIPLEPARTICLES);} else {Nodes[no].u.d.bitflags &= ~(1 << BITFLAG_MULTIPLEPARTICLES);}
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        {double l_tot=0; for(k=0;k<N_RT_FREQ_BINS;k++) {l_tot += Nodes[no].stellar_lum[k];}
        if(l_tot > 0) {Nodes[no].rt_source_lum_s /= l_tot; Extnodes[no].rt_source_lum_vs /= l_tot;}
        else {Nodes[no].rt_source_lum_s = Nodes[no].center; Extnodes[no].rt_source_lum_vs = {};}}
#endif
#ifdef SINK_PHOTONMOMENTUM
        if(Nodes[no].sink_lum > 0) {Nodes[no].sink_lum_grad /= Nodes[no].sink_lum;} else {Nodes[no].sink_lum_grad = {0,0,1};}
#endif
#ifdef SINK_CALC_DISTANCES
        if(Nodes[no].sink_mass > 0) {
            Nodes[no].sink_pos /= Nodes[no].sink_mass;
#ifdef SINK_NODE_MOTION_TRACKED
            Nodes[no].sink_vel /= Nodes[no].sink_mass;
#endif
#ifdef SPECIAL_POINT_MOTION
            Nodes[no].sink_acc /= Nodes[no].sink_mass;
#endif
        }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        if(mass > 0) {MyFloat inv_mass = 1.0/(mass+MIN_REAL_NUMBER); for(k=0;k<6;k++) {Nodes[no].tidal_tensorps_prevstep.data[k] *= inv_mass;}}
#endif
#ifdef DM_SCALARFIELD_SCREENING
        if(Nodes[no].mass_dm > 0) {Nodes[no].s_dm /= Nodes[no].mass_dm; Extnodes[no].vs_dm /= Nodes[no].mass_dm;} else {Nodes[no].s_dm = Nodes[no].center; Extnodes[no].vs_dm = {};}
#endif
    }

    /* Step 5: sync pseudo-particle data across MPI ranks */
    force_exchange_pseudodata();
    force_treeupdate_pseudos(All.MaxPart);

    PRINT_STATUS(" ..tree node moments refreshed.");
#ifdef TREE_INTEGRITY_AUDITS
    /* STARFORGE-type runs refresh near-every step and the moments re-derivation costs about one
       extra refresh, so sample it -- EXCEPT on the first refresh after a particle-list mutation,
       which must always be audited: a corrupted Father[] shows up as a moment mismatch right
       here, and a pure count-stride can let the corruption reach the dynamics first. The cheap
       exactly-once topology walk runs every time. Builds (rarer) always audit fully. */
    {
        static int refresh_audit_counter = 0;
        int audit_moments_now = ((++refresh_audit_counter % 8) == 0) || TreeAuditMomentsPending;
        TreeAuditMomentsPending = 0;
        force_tree_full_audit(audit_moments_now, "refresh");
    }
#endif
}


/*! Tier-0 (always-on) walk-threading validator. Runs only when a maintained rearrange noted a
 *  change (TreeWalkValidatePending); costs nothing otherwise. Walks the Nextnode threading from
 *  the root exactly as the gravity/neighbour walks will, counts the local particles it can reach,
 *  and bounds the iteration so a cycle cannot hang it. On failure it WARNS and condemns the tree
 *  (TreeReconstructFlag -- self-healing at the next ladder, whose reduce uniformizes the local
 *  raise); the current step's walk still proceeds, matching the warn-don't-abort philosophy of
 *  the root-count check. Under TREE_INTEGRITY_AUDITS a failure is fatal instead. */
void force_validate_tree_links(const char *tag)
{
    if(!TreeWalkValidatePending) {return;}
    TreeWalkValidatePending = 0;
    if(Numnodestree <= 0) {return;} /* no standing tree */
    long long nvisit = 0, iter = 0, bound = (long long) NumPart + Numnodestree + NTopleaves + 64;
    /* parent stack: cross-check Father[]/u.d.father against the threading. The walk links and the
       Father pointers are maintained by the same mutation code but are independent structures --
       divergence is the signature of a partial update, and everything downstream that climbs
       Father (moment refresh included) silently follows the wrong pointer. */
    int *nstack = (int *) mymalloc("tree_validate_nstack", 2 * ((size_t) Numnodestree + 2) * sizeof(int));
    int *sstack = nstack + Numnodestree + 2, ntop = 0, diverged = 0;
    int no = All.MaxPart;
    while(no >= 0 && iter < bound)
    {
        iter++;
        while(ntop > 0 && no == sstack[ntop - 1]) {ntop--;} /* walk left the top subtree(s) */
        if(no < All.MaxPart)
        {
            if(ntop > 0 && Father[no] != nstack[ntop - 1]) {diverged = 1; break;}
            nvisit++; no = Nextnode[no];
        }
        else if(no < All.MaxPart + MaxNodes)
        {
            if(Nodes[no].u.d.father != ((ntop > 0) ? nstack[ntop - 1] : -1)) {diverged = 1; break;}
            if(ntop >= Numnodestree + 2) {diverged = 1; break;}
            nstack[ntop] = no; sstack[ntop] = Nodes[no].u.d.sibling; ntop++;
            no = Nodes[no].u.d.nextnode;
        }
        else {no = Nextnode[no - MaxNodes];}
    }
    myfree(nstack);
    if(iter >= bound || diverged || nvisit != (long long) NumPart)
    {
        printf("WARNING task %d [%s]: tree-link validation failed: reached %lld of NumPart=%d local particles%s%s. Condemning the tree for rebuild.\n",
               ThisTask, tag, nvisit, NumPart, (iter >= bound) ? " (iteration bound hit: cycle in the threading?)" : "",
               diverged ? " (Father pointers diverge from the walk threading)" : "");
        fflush(stdout);
        TreeReconstructFlag = 1;
#ifdef TREE_INTEGRITY_AUDITS
        endrun(91570);
#endif
    }
}


#ifdef TREE_INTEGRITY_AUDITS
/*! Tier-1 deep audit (TREE_INTEGRITY_AUDITS builds only; fatal). Verifies:
 *  (a) every local particle is reachable from the root exactly once via the Nextnode threading,
 *      with in-range indices and no cycle;
 *  (b) every Father[] entry is -1 or a valid node index;
 *  (c) with audit_moments (valid only immediately after a build or refresh, when moments are
 *      freshly derived): per-node mass and N_part re-derived from the particle arrays via the
 *      Father chains match the stored values -- restricted to nodes below the top level, whose
 *      stored moments carry no cross-rank pseudo contributions -- and every gas particle's kernel
 *      is covered by its parent's hmax. A wrong Father shows up here as a node-mass mismatch, so
 *      this doubles as the Father-consistency check. Failures are fatal with a tagged report. */
void force_tree_full_audit(int audit_moments, const char *tag)
{
    if(Numnodestree <= 0) {return;}
    int i, fail = 0;
    char *visited = (char *) mymalloc("tree_audit_visited", (size_t) NumPart * sizeof(char));
    memset(visited, 0, (size_t) NumPart * sizeof(char));
    long long iter = 0, bound = (long long) NumPart + Numnodestree + NTopleaves + 64;
    /* parent stack: verify Father[]/u.d.father against the walk threading itself. Without this the
       moments re-derivation below is circular -- it climbs the same Father chain the refresh does,
       so a corrupted Father gives stored and re-derived moments the same wrong answer. */
    int *nstack = (int *) mymalloc("tree_audit_nstack", 2 * ((size_t) Numnodestree + 2) * sizeof(int));
    int *sstack = nstack + Numnodestree + 2, ntop = 0;
    int no = All.MaxPart;
    while(no >= 0 && iter < bound)
    {
        iter++;
        while(ntop > 0 && no == sstack[ntop - 1]) {ntop--;} /* walk left the top subtree(s) */
        if(no < All.MaxPart)
        {
            if(no >= NumPart) {printf("TREE-AUDIT task %d [%s]: walk reached particle index %d >= NumPart=%d\n", ThisTask, tag, no, NumPart); fail = 1; break;}
            if(visited[no]) {printf("TREE-AUDIT task %d [%s]: particle %d visited twice (threading loop)\n", ThisTask, tag, no); fail = 1; break;}
            if(ntop > 0 && Father[no] != nstack[ntop - 1])
            {printf("TREE-AUDIT task %d [%s]: particle %d chained under node %d but Father=%d (threading/Father divergence)\n", ThisTask, tag, no, nstack[ntop - 1], Father[no]); fail = 1; break;}
            visited[no] = 1; no = Nextnode[no];
        }
        else if(no < All.MaxPart + MaxNodes)
        {
            int fwant = (ntop > 0) ? nstack[ntop - 1] : -1;
            if(Nodes[no].u.d.father != fwant)
            {printf("TREE-AUDIT task %d [%s]: node %d threaded under %d but u.d.father=%d\n", ThisTask, tag, no, fwant, Nodes[no].u.d.father); fail = 1; break;}
            if(ntop >= Numnodestree + 2) {printf("TREE-AUDIT task %d [%s]: parent stack overflow at node %d\n", ThisTask, tag, no); fail = 1; break;}
            nstack[ntop] = no; sstack[ntop] = Nodes[no].u.d.sibling; ntop++;
            no = Nodes[no].u.d.nextnode;
        }
        else {no = Nextnode[no - MaxNodes];}
    }
    if(iter >= bound) {printf("TREE-AUDIT task %d [%s]: iteration bound hit -- cycle in the walk threading\n", ThisTask, tag); fail = 1;}
    if(!fail)
    {
        long long nmiss = 0;
        for(i = 0; i < NumPart; i++)
        {
            if(!visited[i]) {if(nmiss < 8) {printf("TREE-AUDIT task %d [%s]: particle %d (Type=%d Mass=%g Father=%d) unreachable from root\n", ThisTask, tag, i, P[i].Type, P[i].Mass, Father[i]);} nmiss++;}
            int f = Father[i];
            if(f != -1 && (f < All.MaxPart || f >= All.MaxPart + Numnodestree)) {printf("TREE-AUDIT task %d [%s]: particle %d has invalid Father=%d\n", ThisTask, tag, i, f); fail = 1;}
        }
        if(nmiss) {printf("TREE-AUDIT task %d [%s]: %lld particles unreachable\n", ThisTask, tag, nmiss); fail = 1;}
    }
    if(!fail && audit_moments)
    {
        double *m_scr = (double *) mymalloc("tree_audit_mscr", (size_t) Numnodestree * sizeof(double));
        long long *np_scr = (long long *) mymalloc("tree_audit_npscr", (size_t) Numnodestree * sizeof(long long));
        memset(m_scr, 0, (size_t) Numnodestree * sizeof(double)); memset(np_scr, 0, (size_t) Numnodestree * sizeof(long long));
        for(i = 0; i < NumPart; i++)
        {
            int f = Father[i];
            while(f >= All.MaxPart && f < All.MaxPart + Numnodestree)
            {
                m_scr[f - All.MaxPart] += P[i].Mass;
                np_scr[f - All.MaxPart]++; /* unconditional, mirroring force_refresh_node_moments -- zero-mass locals are normally gone by audit time (rearrange precedes build/refresh) */
                f = Nodes[f].u.d.father;
            }
            if(P[i].Type == 0 && Father[i] >= All.MaxPart && Father[i] < All.MaxPart + Numnodestree)
            {
                double hneed = DMIN(P[i].KernelRadius, All.MaxKernelRadius);
                if(Extnodes[Father[i]].hmax < 0.999999 * hneed)
                {printf("TREE-AUDIT task %d [%s]: gas particle %d kernel %g exceeds parent node hmax %g\n", ThisTask, tag, i, hneed, Extnodes[Father[i]].hmax); fail = 1;}
            }
        }
        long long nbad = 0;
        for(no = All.MaxPart; no < All.MaxPart + Numnodestree; no++)
        {
            unsigned int bf = Nodes[no].u.d.bitflags;
            if(bf & ((1 << BITFLAG_TOPLEVEL) | (1 << BITFLAG_INTERNAL_TOPLEVEL))) {continue;} /* stored moments include cross-rank pseudo mass here */
            double ms = m_scr[no - All.MaxPart], m0 = Nodes[no].u.d.mass;
            /* mass tolerance 1e-4 relative: stored moments are MyFloat accumulated bottom-up, the
               re-derivation is double accumulated per-particle -- summation-order noise alone
               reaches ~depth*eps_float. N_part is integer-exact. */
            if(np_scr[no - All.MaxPart] != (long long) Nodes[no].N_part || fabs(ms - m0) > 1.e-4 * (fabs(m0) + 1.e-300))
            {
                if(nbad < 8) {printf("TREE-AUDIT task %d [%s]: node %d moments mismatch: N_part stored=%lld rederived=%lld, mass stored=%g rederived=%g\n",
                                     ThisTask, tag, no, (long long) Nodes[no].N_part, np_scr[no - All.MaxPart], m0, ms);}
                nbad++; fail = 1;
            }
        }
        if(nbad) {printf("TREE-AUDIT task %d [%s]: %lld nodes with mismatched moments\n", ThisTask, tag, nbad);}
        myfree(np_scr); myfree(m_scr);
    }
    myfree(nstack);
    myfree(visited);
    if(fail) {endrun(91571);}
}
#endif
