#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef SUBFIND
#include "../structure/subfind/subfind.h"
#endif
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
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
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially (condensed, new feedback routines added,
 * some optimizatins, and new variable/memory conventions added)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

/*! auxiliary variable used to set-up non-recursive walk */
static int last;


/* some modules compute neighbor fluxes explicitly within the force-tree: in these cases, we need to
    take extra care about opening leaves to ensure possible neighbors are not missed, so defined a flag below for it */
#if (defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(SINGLE_STAR_SINK_DYNAMICS) || defined(GRAVITY_ACCURATE_FEWBODY_INTEGRATION))
#define NEIGHBORS_MUST_BE_COMPUTED_EXPLICITLY_IN_FORCETREE
#endif

#define ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING /* comment out to revert to behavior of taking 'greater' softening in pairwise kernel interactions with adaptive softenings enabled */

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


#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport, mutex_partnodedrift;

#define LOCK_NEXPORT         pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT       pthread_mutex_unlock(&mutex_nexport);
#define LOCK_PARTNODEDRIFT   pthread_mutex_lock(&mutex_partnodedrift);
#define UNLOCK_PARTNODEDRIFT pthread_mutex_unlock(&mutex_partnodedrift);
/*! The cost computation for the tree-gravity (required for the domain
 decomposition) is not exactly thread-safe if THREAD_SAFE_COSTS is not defined.
 However using locks for an exactly thread-safe cost computiation results in a
 significant (~25%) performance penalty in the tree-walk while having only an
 extremely small effect on the obtained costs. The domain decomposition should
 thus not be significantly changed if THREAD_SAFE_COSTS is not used.
 (modern code no longer includes this option - need to consult legacy code) */
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#define LOCK_PARTNODEDRIFT
#define UNLOCK_PARTNODEDRIFT
#endif



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

    int flag;
    do
    {
        Numnodestree = force_treebuild_single(npart, mp);

        MPI_Allreduce(&Numnodestree, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

        if(flag == -1)
        {
            force_treefree();

            if(ThisTask == 0)
                printf("Increasing TreeAllocFactor=%g", All.TreeAllocFactor);

            All.TreeAllocFactor *= 1.15;

            if(ThisTask == 0)
                printf("new value=%g\n", All.TreeAllocFactor);

            force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
        }
    }
    while(flag == -1);

    force_flag_localnodes();

    force_exchange_pseudodata();

    force_treeupdate_pseudos(All.MaxPart);

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
    int i, j, k, subnode = 0, shift, parent, numnodes, rep;
    int nfree, th, nn, no;
    struct NODE *nfreep;
    MyFloat lenhalf;
    peanokey key, morton, th_key, *morton_list;


    /* create an empty root node  */
    nfree = All.MaxPart;		/* index of first free node */
    nfreep = &Nodes[nfree];	/* select first node */

    nfreep->len = DomainLen;
    for(j = 0; j < 3; j++)
        nfreep->center[j] = DomainCenter[j];
    for(j = 0; j < 8; j++)
        nfreep->u.suns[j] = -1;


    numnodes = 1;
    nfreep++;
    nfree++;

    /* create a set of empty nodes corresponding to the top-level domain
     * grid. We need to generate these nodes first to make sure that we have a
     * complete top-level tree which allows the easy insertion of the
     * pseudo-particles at the right place
     */

    force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);

    /* if a high-resolution region in a global tree is used, we need to generate
     * an additional set empty nodes to make sure that we have a complete
     * top-level tree for the high-resolution inset
     */

    nfreep = &Nodes[nfree];
    parent = -1;			/* note: will not be used below before it is changed */

    morton_list = (peanokey *) mymalloc("morton_list", NumPart * sizeof(peanokey));

    /* now we insert all particles */
    for(k = 0; k < npart; k++)
    {
        if(mp)
            i = mp[k].index;
        else
            i = k;

        rep = 0;

        key = peano_and_morton_key((int) ((P[i].Pos[0] - DomainCorner[0]) * DomainFac),
                                   (int) ((P[i].Pos[1] - DomainCorner[1]) * DomainFac),
                                   (int) ((P[i].Pos[2] - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION,
                                   &morton);
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
            if(th >= All.MaxPart)	/* we are dealing with an internal node */
            {
                if(shift >= 0)
                {
                    subnode = ((morton >> shift) & 7);
                }
                else
                {
                    subnode = 0;
                    if(P[i].Pos[0] > Nodes[th].center[0])
                        subnode += 1;
                    if(P[i].Pos[1] > Nodes[th].center[1])
                        subnode += 2;
                    if(P[i].Pos[2] > Nodes[th].center[2])
                        subnode += 4;
                }

#ifndef NOTREERND
                if(Nodes[th].len < EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[P[i].Type])
                {
                    /* seems like we're dealing with particles at identical (or extremely close)
                     * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
                     * of tree are still correct, but this will only happen well below gravitational softening
                     * length-scale anyway.
                     */
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
                    subnode = (int) (8.0 * get_random_number((P[i].ID + rep) % (RNDTABLE + (rep & 3))));
#else
                    subnode = (int) (8.0 * get_random_number(P[i].ID));
#endif

                    if(subnode >= 8)
                        subnode = 7;
                }
#endif

                nn = Nodes[th].u.suns[subnode];

                shift -= 3;

                if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
                {
                    parent = th;
                    th = nn;
                    rep++;
                }
                else
                {
                    /* here we have found an empty slot where we can attach
                     * the new particle as a leaf.
                     */
                    Nodes[th].u.suns[subnode] = i;
                    break;	/* done for this particle */
                }
            }
            else
            {
                /* We try to insert into a leaf with a single particle.  Need
                 * to generate a new internal node at this point.
                 */
                Nodes[parent].u.suns[subnode] = nfree;

                nfreep->len = 0.5 * Nodes[parent].len;
                lenhalf = 0.25 * Nodes[parent].len;

                if(subnode & 1)
                    nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
                else
                    nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

                if(subnode & 2)
                    nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
                else
                    nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

                if(subnode & 4)
                    nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
                else
                    nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

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
                    if(P[th].Pos[0] > nfreep->center[0])
                        subnode += 1;
                    if(P[th].Pos[1] > nfreep->center[1])
                        subnode += 2;
                    if(P[th].Pos[2] > nfreep->center[2])
                        subnode += 4;
                }

#ifndef NOTREERND
                if(nfreep->len < EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[P[th].Type])
                {
                    /* seems like we're dealing with particles at identical (or extremely close)
                     * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
                     * of tree are still correct, but this will only happen well below gravitational softening
                     * length-scale anyway.
                     */
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
                    subnode = (int) (8.0 * get_random_number((P[th].ID + rep) % (RNDTABLE + (rep & 3))));
#else
                    subnode = (int) (8.0 * get_random_number(P[th].ID));
#endif

                    if(subnode >= 8)
                        subnode = 7;
                }
#endif
                nfreep->u.suns[subnode] = th;

                th = nfree;	/* resume trying to insert the new particle at
                             * the newly created internal node
                             */

                numnodes++;
                nfree++;
                nfreep++;

                if((numnodes) >= MaxNodes)
                {
                    printf("task %d: maximum number %d of tree-nodes reached for particle %d.\n", ThisTask, MaxNodes, i);

                    if(All.TreeAllocFactor > 5.0)
                    {
                        printf("task %d: looks like a serious problem for particle %d, stopping with particle dump.\n", ThisTask, i);
                        dump_particles();
                        endrun(1);
                    }
                    else
                    {
                        myfree(morton_list);
                        return -1;
                    }
                }
            }
        }
    }

    myfree(morton_list);


    /* insert the pseudo particles that represent the mass distribution of other domains */
    force_insert_pseudo_particles();


    /* now compute the multipole moments recursively */
    last = -1;

    force_update_node_recursive(All.MaxPart, -1, -1);

    if(last >= All.MaxPart)
    {
        if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
            Nextnode[last - MaxNodes] = -1;
        else
            Nodes[last].u.d.nextnode = -1;
    }
    else
        Nextnode[last] = -1;

    return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 */
void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
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

                    if((*nodecount) >= MaxNodes || (*nodecount) >= MaxTopNodes)
                    {
                        printf("task %d: maximum number MaxNodes=%d of tree-nodes reached."
                               "MaxTopNodes=%d NTopnodes=%d NTopleaves=%d nodecount=%d\n",
                               ThisTask, MaxNodes, MaxTopNodes, NTopnodes, NTopleaves, *nodecount);
                        printf("in create empty nodes\n");
                        dump_particles();
                        endrun(11);
                    }

                    force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
                                             bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
                }
    }
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
 *  dependent on local mass.
 */
void force_update_node_recursive(int no, int sib, int father)
{
    int j, jj, k, p, pp, nextsib, suns[8], count_particles, multiple_flag;
    MyFloat hmax, vmax, v, divVmax, divVel, s[3], vs[3], mass;
    struct particle_data *pa;

#ifdef DM_SCALARFIELD_SCREENING
    MyFloat s_dm[3], vs_dm[3], mass_dm;
#endif
#ifdef RT_USE_GRAVTREE
    MyFloat stellar_lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS]={0};
#endif
    for(j=0;j<N_RT_FREQ_BINS;j++) {stellar_lum[j]=0;}
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MyFloat rt_source_lum_s[3];
    MyFloat rt_source_lum_vs[3];
#endif

    MyFloat maxsoft;

    if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
        for(j = 0; j < 8; j++)
            suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will overwrite one element (union!) */
        if(last >= 0)
        {
            if(last >= All.MaxPart)
            {
                if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
                    Nextnode[last - MaxNodes] = no;
                else
                    Nodes[last].u.d.nextnode = no;
            }
            else
                Nextnode[last] = no;
        }

        last = no;

#ifdef RT_USE_TREECOL_FOR_NH
        MyFloat gasmass = 0;
#endif
#ifdef RT_USE_GRAVTREE
        for(j=0;j<N_RT_FREQ_BINS;j++) {stellar_lum[j]=0;}
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        rt_source_lum_s[0] = rt_source_lum_vs[0] = 0;
        rt_source_lum_s[1] = rt_source_lum_vs[1] = 0;
        rt_source_lum_s[2] = rt_source_lum_vs[2] = 0;
#endif
#ifdef BH_PHOTONMOMENTUM
        MyFloat bh_lum,bh_lum_grad[3]; bh_lum=bh_lum_grad[0]=bh_lum_grad[1]=bh_lum_grad[2]=0;
#endif
#ifdef BH_CALC_DISTANCES
        MyFloat bh_mass=0, bh_pos_times_mass[3]={0,0,0};   /* position of each black hole in the node times its mass; divide by total mass at the end to get COM */
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
        MyFloat bh_mom[3] = {0,0,0}; int N_BH = 0;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        MyFloat max_feedback_vel=0;
#endif        
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
        mass_dm = 0;
        s_dm[0] = vs_dm[0] = 0;
        s_dm[1] = vs_dm[1] = 0;
        s_dm[2] = vs_dm[2] = 0;
#endif
        mass = 0;
        s[0] = 0;
        s[1] = 0;
        s[2] = 0;
        vs[0] = 0;
        vs[1] = 0;
        vs[2] = 0;
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

                if(jj < 8)	/* yes, we do */
                    nextsib = pp;
                else
                    nextsib = sib;

                force_update_node_recursive(p, nextsib, no);

                if(p >= All.MaxPart)	/* an internal node or pseudo particle */
                {
                    if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
                    {
                        /* nothing to be done here because the mass of the
                         * pseudo-particle is still zero. This will be changed
                         * later.
                         */
                    }
                    else
                    {
                        mass += (Nodes[p].u.d.mass);
                        s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
                        s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
                        s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
                        vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
                        vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
                        vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);
#ifdef RT_USE_TREECOL_FOR_NH
                        gasmass += Nodes[p].gasmass;
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
                        rt_source_lum_s[0] += (l_tot * Nodes[p].rt_source_lum_s[0]);
                        rt_source_lum_s[1] += (l_tot * Nodes[p].rt_source_lum_s[1]);
                        rt_source_lum_s[2] += (l_tot * Nodes[p].rt_source_lum_s[2]);
                        rt_source_lum_vs[0] += (l_tot * Extnodes[p].rt_source_lum_vs[0]);
                        rt_source_lum_vs[1] += (l_tot * Extnodes[p].rt_source_lum_vs[1]);
                        rt_source_lum_vs[2] += (l_tot * Extnodes[p].rt_source_lum_vs[2]);
#endif
#ifdef BH_PHOTONMOMENTUM
                        bh_lum += Nodes[p].bh_lum;
                        bh_lum_grad[0] += Nodes[p].bh_lum * Nodes[p].bh_lum_grad[0];
                        bh_lum_grad[1] += Nodes[p].bh_lum * Nodes[p].bh_lum_grad[1];
                        bh_lum_grad[2] += Nodes[p].bh_lum * Nodes[p].bh_lum_grad[2];
#endif
#ifdef BH_CALC_DISTANCES
                        bh_mass += Nodes[p].bh_mass;
                        bh_pos_times_mass[0] += Nodes[p].bh_pos[0] * Nodes[p].bh_mass;
                        bh_pos_times_mass[1] += Nodes[p].bh_pos[1] * Nodes[p].bh_mass;
                        bh_pos_times_mass[2] += Nodes[p].bh_pos[2] * Nodes[p].bh_mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
                        bh_mom[0] += Nodes[p].bh_vel[0] * Nodes[p].bh_mass;
                        bh_mom[1] += Nodes[p].bh_vel[1] * Nodes[p].bh_mass;
                        bh_mom[2] += Nodes[p].bh_vel[2] * Nodes[p].bh_mass;
                        N_BH += Nodes[p].N_BH;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                        if(Nodes[p].bh_mass > 0) {max_feedback_vel = DMAX(Nodes[p].MaxFeedbackVel, max_feedback_vel);}
#endif                        
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
                        mass_dm += (Nodes[p].mass_dm);
                        s_dm[0] += (Nodes[p].mass_dm * Nodes[p].s_dm[0]);
                        s_dm[1] += (Nodes[p].mass_dm * Nodes[p].s_dm[1]);
                        s_dm[2] += (Nodes[p].mass_dm * Nodes[p].s_dm[2]);
                        vs_dm[0] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[0]);
                        vs_dm[1] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[1]);
                        vs_dm[2] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[2]);
#endif
                        if(Nodes[p].u.d.mass > 0)
                        {
                            if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
                                count_particles += 2;
                            else
                                count_particles++;
                        }

                        if(Extnodes[p].hmax > hmax)
                            hmax = Extnodes[p].hmax;

                        if(Extnodes[p].vmax > vmax)
                            vmax = Extnodes[p].vmax;
                        if(Extnodes[p].divVmax > divVmax)
                            divVmax = Extnodes[p].divVmax;

                        /* update of the maximum gravitational softening in the node */
                        if(Nodes[p].maxsoft > maxsoft)
                            maxsoft = Nodes[p].maxsoft;

                    }
                }
                else		/* a particle */
                {
                    count_particles++;

                    pa = &P[p];

                    mass += (pa->Mass);
                    s[0] += (pa->Mass * pa->Pos[0]);
                    s[1] += (pa->Mass * pa->Pos[1]);
                    s[2] += (pa->Mass * pa->Pos[2]);
                    vs[0] += (pa->Mass * pa->Vel[0]);
                    vs[1] += (pa->Mass * pa->Vel[1]);
                    vs[2] += (pa->Mass * pa->Vel[2]);
#ifdef RT_USE_TREECOL_FOR_NH
                    if(pa->Type == 0) gasmass += pa->Mass;
#ifdef BH_ALPHADISK_ACCRETION
                    if(pa->Type == 5) gasmass += BPP(p).BH_Mass_AlphaDisk; // gas at the inner edge of a disk should not see a hole due to the sink
#endif
#endif
#ifdef RT_USE_GRAVTREE
                    double lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
                    double chimes_lum_G0[CHIMES_LOCAL_UV_NBINS];
                    double chimes_lum_ion[CHIMES_LOCAL_UV_NBINS];
                    int active_check = rt_get_source_luminosity_chimes(p,1,lum,chimes_lum_G0, chimes_lum_ion);
#else
                    int active_check = rt_get_source_luminosity(p,1,lum);
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
                        rt_source_lum_s[0] += (l_sum * pa->Pos[0]);
                        rt_source_lum_s[1] += (l_sum * pa->Pos[1]);
                        rt_source_lum_s[2] += (l_sum * pa->Pos[2]);
                        rt_source_lum_vs[0] += (l_sum * pa->Vel[0]);
                        rt_source_lum_vs[1] += (l_sum * pa->Vel[1]);
                        rt_source_lum_vs[2] += (l_sum * pa->Vel[2]);
#endif
                    }
#endif



#ifdef BH_PHOTONMOMENTUM
                    if(pa->Type == 5)
                    {
                        if((pa->Mass>0)&&(pa->DensAroundStar>0)&&(pa->BH_Mdot>0))
                        {
                            double BHLum = bh_lum_bol(pa->BH_Mdot, pa->BH_Mass, p);
                            bh_lum += BHLum;
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                            for(k=0;k<3;k++) {bh_lum_grad[k] += BHLum * pa->BH_Specific_AngMom[k];}
#else
                            for(k=0;k<3;k++) {bh_lum_grad[k] += BHLum * pa->GradRho[k];}
#endif
                        }
                    }
#endif
#ifdef BH_CALC_DISTANCES
                    if(pa->Type == 5)
                    {
                        bh_mass += pa->Mass;    /* actual value is not used for distances */
                        bh_pos_times_mass[0] += pa->Pos[0] * pa->Mass;  /* positition times mass; divide by total mass later */
                        bh_pos_times_mass[1] += pa->Pos[1] * pa->Mass;
                        bh_pos_times_mass[2] += pa->Pos[2] * pa->Mass;
#ifdef SINGLE_STAR_TIMESTEPPING
                        bh_mom[0] += pa->Vel[0] * pa->Mass;
                        bh_mom[1] += pa->Vel[1] * pa->Mass;
                        bh_mom[2] += pa->Vel[2] * pa->Mass;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                        max_feedback_vel = DMAX(pa->MaxFeedbackVel, max_feedback_vel);
#endif                        
#endif
                    }
#endif


#ifdef DM_SCALARFIELD_SCREENING
                    if(pa->Type != 0)
                    {
                        mass_dm += (pa->Mass);
                        s_dm[0] += (pa->Mass * pa->Pos[0]);
                        s_dm[1] += (pa->Mass * pa->Pos[1]);
                        s_dm[2] += (pa->Mass * pa->Pos[2]);
                        vs_dm[0] += (pa->Mass * pa->Vel[0]);
                        vs_dm[1] += (pa->Mass * pa->Vel[1]);
                        vs_dm[2] += (pa->Mass * pa->Vel[2]);
                    }
#endif
                    if(pa->Type == 0)
                    {
                        if(PPP[p].Hsml > hmax) {hmax = PPP[p].Hsml;}
                        divVel = P[p].Particle_DivVel;
                        if(divVel > divVmax) {divVmax = divVel;}
                    }

                    for(k = 0; k < 3; k++) {if((v = fabs(pa->Vel[k])) > vmax) {vmax = v;}}

                    /* update of the maximum gravitational softening  */
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                    if(PPP[p].AGS_Hsml > maxsoft) {maxsoft = PPP[p].AGS_Hsml;}
#else
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
                    if(All.ForceSoftening[pa->Type] > maxsoft) {maxsoft = All.ForceSoftening[pa->Type];}
#else
                    if(pa->Type == 0)
                    {
                        if(PPP[p].Hsml > maxsoft) {maxsoft = PPP[p].Hsml;}
                    } else {
                        if(All.ForceSoftening[pa->Type] > maxsoft) {maxsoft = All.ForceSoftening[pa->Type];}
                    }
#endif
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
                    if(pa->Type == 5) if(PPP[p].Hsml > maxsoft) {maxsoft = PPP[p].Hsml;}
#endif
                }
            }
        }


        if(mass)
        {
            s[0] /= mass;
            s[1] /= mass;
            s[2] /= mass;
            vs[0] /= mass;
            vs[1] /= mass;
            vs[2] /= mass;
        }
        else
        {
            s[0] = Nodes[no].center[0];
            s[1] = Nodes[no].center[1];
            s[2] = Nodes[no].center[2];
            vs[0] = 0;
            vs[1] = 0;
            vs[2] = 0;
        }

#ifdef RT_SEPARATELY_TRACK_LUMPOS
        double l_tot=0; for(k=0;k<N_RT_FREQ_BINS;k++) {l_tot += stellar_lum[k];}
        if(l_tot)
        {
            rt_source_lum_s[0] /= l_tot;
            rt_source_lum_s[1] /= l_tot;
            rt_source_lum_s[2] /= l_tot;
            rt_source_lum_vs[0] /= l_tot;
            rt_source_lum_vs[1] /= l_tot;
            rt_source_lum_vs[2] /= l_tot;
        }
        else
        {
            rt_source_lum_s[0] = Nodes[no].center[0];
            rt_source_lum_s[1] = Nodes[no].center[1];
            rt_source_lum_s[2] = Nodes[no].center[2];
            rt_source_lum_vs[0] = 0;
            rt_source_lum_vs[1] = 0;
            rt_source_lum_vs[2] = 0;
        }
#endif
#ifdef BH_PHOTONMOMENTUM
        if(bh_lum)
        {
            bh_lum_grad[0] /= bh_lum; bh_lum_grad[1] /= bh_lum; bh_lum_grad[2] /= bh_lum;
        } else {
            bh_lum_grad[0]=bh_lum_grad[1]=0; bh_lum_grad[2]=1;
        }
#endif
#ifdef DM_SCALARFIELD_SCREENING
        if(mass_dm)
        {
            s_dm[0] /= mass_dm;
            s_dm[1] /= mass_dm;
            s_dm[2] /= mass_dm;
            vs_dm[0] /= mass_dm;
            vs_dm[1] /= mass_dm;
            vs_dm[2] /= mass_dm;
        }
        else
        {
            s_dm[0] = Nodes[no].center[0];
            s_dm[1] = Nodes[no].center[1];
            s_dm[2] = Nodes[no].center[2];
            vs_dm[0] = 0;
            vs_dm[1] = 0;
            vs_dm[2] = 0;
        }
#endif


        Nodes[no].Ti_current = All.Ti_Current;
        Nodes[no].u.d.mass = mass;
        Nodes[no].u.d.s[0] = s[0];
        Nodes[no].u.d.s[1] = s[1];
        Nodes[no].u.d.s[2] = s[2];
        Nodes[no].GravCost = 0;
#ifdef RT_USE_TREECOL_FOR_NH
        Nodes[no].gasmass = gasmass;
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
        Nodes[no].rt_source_lum_s[0] = rt_source_lum_s[0];
        Nodes[no].rt_source_lum_s[1] = rt_source_lum_s[1];
        Nodes[no].rt_source_lum_s[2] = rt_source_lum_s[2];
        Extnodes[no].rt_source_lum_vs[0] = rt_source_lum_vs[0];
        Extnodes[no].rt_source_lum_vs[1] = rt_source_lum_vs[1];
        Extnodes[no].rt_source_lum_vs[2] = rt_source_lum_vs[2];
        Extnodes[no].rt_source_lum_dp[0] = 0;
        Extnodes[no].rt_source_lum_dp[1] = 0;
        Extnodes[no].rt_source_lum_dp[2] = 0;
#endif
#ifdef BH_PHOTONMOMENTUM
        Nodes[no].bh_lum = bh_lum;
        Nodes[no].bh_lum_grad[0] = bh_lum_grad[0];
        Nodes[no].bh_lum_grad[1] = bh_lum_grad[1];
        Nodes[no].bh_lum_grad[2] = bh_lum_grad[2];
#endif
#ifdef BH_CALC_DISTANCES
        Nodes[no].bh_mass = bh_mass;
        if(bh_mass > 0)
            {
                Nodes[no].bh_pos[0] = bh_pos_times_mass[0] / bh_mass;  /* weighted position is sum(pos*mass)/sum(mass) */
                Nodes[no].bh_pos[1] = bh_pos_times_mass[1] / bh_mass;
                Nodes[no].bh_pos[2] = bh_pos_times_mass[2] / bh_mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
                Nodes[no].bh_vel[0] = bh_mom[0] / bh_mass;
                Nodes[no].bh_vel[1] = bh_mom[1] / bh_mass;
                Nodes[no].bh_vel[2] = bh_mom[2] / bh_mass;
                Nodes[no].N_BH = N_BH;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                Nodes[no].MaxFeedbackVel = max_feedback_vel;
#endif                        
#endif
            }
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Nodes[no].s_dm[0] = s_dm[0];
        Nodes[no].s_dm[1] = s_dm[1];
        Nodes[no].s_dm[2] = s_dm[2];
        Nodes[no].mass_dm = mass_dm;
        Extnodes[no].vs_dm[0] = vs_dm[0];
        Extnodes[no].vs_dm[1] = vs_dm[1];
        Extnodes[no].vs_dm[2] = vs_dm[2];
        Extnodes[no].dp_dm[0] = 0;
        Extnodes[no].dp_dm[1] = 0;
        Extnodes[no].dp_dm[2] = 0;
#endif

        Extnodes[no].Ti_lastkicked = All.Ti_Current;
        Extnodes[no].Flag = GlobFlag;
        Extnodes[no].vs[0] = vs[0];
        Extnodes[no].vs[1] = vs[1];
        Extnodes[no].vs[2] = vs[2];
        Extnodes[no].hmax = hmax;
        Extnodes[no].vmax = vmax;
        Extnodes[no].divVmax = divVmax;
        Extnodes[no].dp[0] = 0;
        Extnodes[no].dp[1] = 0;
        Extnodes[no].dp[2] = 0;

        if(count_particles > 1)	/* this flags that the node represents more than one particle */
            multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
        else
            multiple_flag = 0;

        Nodes[no].u.d.bitflags = multiple_flag;
        Nodes[no].maxsoft = maxsoft;
        Nodes[no].u.d.sibling = sib;
        Nodes[no].u.d.father = father;
    }
    else				/* single particle or pseudo particle */
    {
        if(last >= 0)
        {
            if(last >= All.MaxPart)
            {
                if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
                    Nextnode[last - MaxNodes] = no;
                else
                    Nodes[last].u.d.nextnode = no;
            }
            else
                Nextnode[last] = no;
        }

        last = no;

        if(no < All.MaxPart)	/* only set it for single particles */
            Father[no] = father;
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
        MyFloat s[3];
        MyFloat vs[3];
        MyFloat mass;
        MyFloat hmax;
        MyFloat vmax;
        MyFloat divVmax;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        MyFloat maxsoft;
#endif
#ifdef RT_USE_GRAVTREE
        MyFloat stellar_lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
        double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS];
        double chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS];
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        MyFloat rt_source_lum_s[3];
        MyFloat rt_source_lum_vs[3];
#endif
#ifdef BH_PHOTONMOMENTUM
        MyFloat bh_lum,bh_lum_grad[3];
#endif
#ifdef BH_CALC_DISTANCES
        MyFloat bh_mass;
        MyFloat bh_pos[3];
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
        MyFloat bh_vel[3];
        int N_BH;
#ifdef  SINGLE_STAR_FB_TIMESTEPLIMIT
        MyFloat MaxFeedbackVel;
#endif        
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
        MyFloat s_dm[3];
        MyFloat vs_dm[3];
        MyFloat mass_dm;
#endif
        unsigned int bitflags;
#ifdef PAD_STRUCTURES
#ifndef DOUBLEPRECISION
        int pad[5];
#else
#if (DOUBLEPRECISION+0) == 2
        /* mixed precision */
        int pad[5];
#else
        int pad[3];
#endif
#endif				/* DOUBLEPRECISION  */
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
            DomainMoment[i].s[0] = Nodes[no].u.d.s[0];
            DomainMoment[i].s[1] = Nodes[no].u.d.s[1];
            DomainMoment[i].s[2] = Nodes[no].u.d.s[2];
            DomainMoment[i].vs[0] = Extnodes[no].vs[0];
            DomainMoment[i].vs[1] = Extnodes[no].vs[1];
            DomainMoment[i].vs[2] = Extnodes[no].vs[2];
            DomainMoment[i].mass = Nodes[no].u.d.mass;
            DomainMoment[i].hmax = Extnodes[no].hmax;
            DomainMoment[i].vmax = Extnodes[no].vmax;
            DomainMoment[i].divVmax = Extnodes[no].divVmax;
            DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
            DomainMoment[i].maxsoft = Nodes[no].maxsoft;
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
            DomainMoment[i].rt_source_lum_s[0] = Nodes[no].rt_source_lum_s[0];
            DomainMoment[i].rt_source_lum_s[1] = Nodes[no].rt_source_lum_s[1];
            DomainMoment[i].rt_source_lum_s[2] = Nodes[no].rt_source_lum_s[2];
            DomainMoment[i].rt_source_lum_vs[0] = Extnodes[no].rt_source_lum_vs[0];
            DomainMoment[i].rt_source_lum_vs[1] = Extnodes[no].rt_source_lum_vs[1];
            DomainMoment[i].rt_source_lum_vs[2] = Extnodes[no].rt_source_lum_vs[2];
#endif
#ifdef BH_PHOTONMOMENTUM
            DomainMoment[i].bh_lum = Nodes[no].bh_lum;
            DomainMoment[i].bh_lum_grad[0] = Nodes[no].bh_lum_grad[0];
            DomainMoment[i].bh_lum_grad[1] = Nodes[no].bh_lum_grad[1];
            DomainMoment[i].bh_lum_grad[2] = Nodes[no].bh_lum_grad[2];
#endif
#ifdef BH_CALC_DISTANCES
            DomainMoment[i].bh_mass = Nodes[no].bh_mass;
            DomainMoment[i].bh_pos[0] = Nodes[no].bh_pos[0];
            DomainMoment[i].bh_pos[1] = Nodes[no].bh_pos[1];
            DomainMoment[i].bh_pos[2] = Nodes[no].bh_pos[2];
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
            DomainMoment[i].bh_vel[0] = Nodes[no].bh_vel[0];
            DomainMoment[i].bh_vel[1] = Nodes[no].bh_vel[1];
            DomainMoment[i].bh_vel[2] = Nodes[no].bh_vel[2];
            DomainMoment[i].N_BH = Nodes[no].N_BH;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
            DomainMoment[i].MaxFeedbackVel = Nodes[no].MaxFeedbackVel;
#endif            
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
            DomainMoment[i].s_dm[0] = Nodes[no].s_dm[0];
            DomainMoment[i].s_dm[1] = Nodes[no].s_dm[1];
            DomainMoment[i].s_dm[2] = Nodes[no].s_dm[2];
            DomainMoment[i].mass_dm = Nodes[no].mass_dm;
            DomainMoment[i].vs_dm[0] = Extnodes[no].vs_dm[0];
            DomainMoment[i].vs_dm[1] = Extnodes[no].vs_dm[1];
            DomainMoment[i].vs_dm[2] = Extnodes[no].vs_dm[2];
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
#ifdef USE_MPI_IN_PLACE
        MPI_Allgatherv(MPI_IN_PLACE, recvcounts[ThisTask],
                       MPI_BYTE, &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);
#else
        MPI_Allgatherv(&DomainMoment[DomainStartList[ThisTask * MULTIPLEDOMAINS + m]], recvcounts[ThisTask],
                       MPI_BYTE, &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);
#endif
    }

    myfree(recvoffset);
    myfree(recvcounts);


    for(ta = 0; ta < NTask; ta++)
        if(ta != ThisTask)
            for(m = 0; m < MULTIPLEDOMAINS; m++)
                for(i = DomainStartList[ta * MULTIPLEDOMAINS + m]; i <= DomainEndList[ta * MULTIPLEDOMAINS + m]; i++)
                {
                    no = DomainNodeIndex[i];

                    Nodes[no].u.d.s[0] = DomainMoment[i].s[0];
                    Nodes[no].u.d.s[1] = DomainMoment[i].s[1];
                    Nodes[no].u.d.s[2] = DomainMoment[i].s[2];
                    Extnodes[no].vs[0] = DomainMoment[i].vs[0];
                    Extnodes[no].vs[1] = DomainMoment[i].vs[1];
                    Extnodes[no].vs[2] = DomainMoment[i].vs[2];
                    Nodes[no].u.d.mass = DomainMoment[i].mass;
                    Extnodes[no].hmax = DomainMoment[i].hmax;
                    Extnodes[no].vmax = DomainMoment[i].vmax;
                    Extnodes[no].divVmax = DomainMoment[i].divVmax;
                    Nodes[no].u.d.bitflags =
                    (Nodes[no].u.d.bitflags & (~BITFLAG_MASK)) | (DomainMoment[i].bitflags & BITFLAG_MASK);
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                    Nodes[no].maxsoft = DomainMoment[i].maxsoft;
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
                    Nodes[no].rt_source_lum_s[0] = DomainMoment[i].rt_source_lum_s[0];
                    Nodes[no].rt_source_lum_s[1] = DomainMoment[i].rt_source_lum_s[1];
                    Nodes[no].rt_source_lum_s[2] = DomainMoment[i].rt_source_lum_s[2];
                    Extnodes[no].rt_source_lum_vs[0] = DomainMoment[i].rt_source_lum_vs[0];
                    Extnodes[no].rt_source_lum_vs[1] = DomainMoment[i].rt_source_lum_vs[1];
                    Extnodes[no].rt_source_lum_vs[2] = DomainMoment[i].rt_source_lum_vs[2];
#endif
#ifdef BH_PHOTONMOMENTUM
                    Nodes[no].bh_lum = DomainMoment[i].bh_lum;
                    Nodes[no].bh_lum_grad[0] = DomainMoment[i].bh_lum_grad[0];
                    Nodes[no].bh_lum_grad[1] = DomainMoment[i].bh_lum_grad[1];
                    Nodes[no].bh_lum_grad[2] = DomainMoment[i].bh_lum_grad[2];
#endif
#ifdef BH_CALC_DISTANCES
                    Nodes[no].bh_mass = DomainMoment[i].bh_mass;
                    Nodes[no].bh_pos[0] = DomainMoment[i].bh_pos[0];
                    Nodes[no].bh_pos[1] = DomainMoment[i].bh_pos[1];
                    Nodes[no].bh_pos[2] = DomainMoment[i].bh_pos[2];
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
                    Nodes[no].bh_vel[0] = DomainMoment[i].bh_vel[0];
                    Nodes[no].bh_vel[1] = DomainMoment[i].bh_vel[1];
                    Nodes[no].bh_vel[2] = DomainMoment[i].bh_vel[2];
                    Nodes[no].N_BH = DomainMoment[i].N_BH;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                    Nodes[no].MaxFeedbackVel = DomainMoment[i].MaxFeedbackVel;
#endif                        
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
                    Nodes[no].s_dm[0] = DomainMoment[i].s_dm[0];
                    Nodes[no].s_dm[1] = DomainMoment[i].s_dm[1];
                    Nodes[no].s_dm[2] = DomainMoment[i].s_dm[2];
                    Nodes[no].mass_dm = DomainMoment[i].mass_dm;
                    Extnodes[no].vs_dm[0] = DomainMoment[i].vs_dm[0];
                    Extnodes[no].vs_dm[1] = DomainMoment[i].vs_dm[1];
                    Extnodes[no].vs_dm[2] = DomainMoment[i].vs_dm[2];
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
    MyFloat s[3], vs[3], mass;

#ifdef RT_USE_GRAVTREE
    MyFloat stellar_lum[N_RT_FREQ_BINS]={0};
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_stellar_lum_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_stellar_lum_ion[CHIMES_LOCAL_UV_NBINS]={0};
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MyFloat rt_source_lum_s[3];
    MyFloat rt_source_lum_vs[3];
#endif
#ifdef DM_SCALARFIELD_SCREENING
    MyFloat s_dm[3], vs_dm[3], mass_dm;
#endif

    MyFloat maxsoft;

#ifdef RT_SEPARATELY_TRACK_LUMPOS
    rt_source_lum_s[0] = 0;
    rt_source_lum_s[1] = 0;
    rt_source_lum_s[2] = 0;
    rt_source_lum_vs[0] = 0;
    rt_source_lum_vs[1] = 0;
    rt_source_lum_vs[2] = 0;
#endif
#ifdef BH_PHOTONMOMENTUM
    MyFloat bh_lum,bh_lum_grad[3]; bh_lum=bh_lum_grad[0]=bh_lum_grad[1]=bh_lum_grad[2]=0;
#endif
#ifdef BH_CALC_DISTANCES
    MyFloat bh_mass=0;
    MyFloat bh_pos_times_mass[3]={0,0,0};
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
    MyFloat bh_mom[3] = {0,0,0};
    int N_BH = 0;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    MyFloat max_feedback_vel=0;
#endif    
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
    mass_dm = 0;
    s_dm[0] = vs_dm[0] = 0;
    s_dm[1] = vs_dm[1] = 0;
    s_dm[2] = vs_dm[2] = 0;
#endif
    mass = 0;
    s[0] = 0;
    s[1] = 0;
    s[2] = 0;
    vs[0] = 0;
    vs[1] = 0;
    vs[2] = 0;
    hmax = 0;
    vmax = 0;
    divVmax = 0;
    count_particles = 0;
    maxsoft = 0;

    p = Nodes[no].u.d.nextnode;

    for(j = 0; j < 8; j++)	/* since we are dealing with top-level nodes, we now that there are 8 consecutive daughter nodes */
    {
        if(p >= All.MaxPart && p < All.MaxPart + MaxNodes)	/* internal node */
        {
            if(Nodes[p].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
                force_treeupdate_pseudos(p);

            mass += (Nodes[p].u.d.mass);
            s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
            s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
            s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
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
            rt_source_lum_s[0] += (l_tot * Nodes[p].rt_source_lum_s[0]);
            rt_source_lum_s[1] += (l_tot * Nodes[p].rt_source_lum_s[1]);
            rt_source_lum_s[2] += (l_tot * Nodes[p].rt_source_lum_s[2]);
            rt_source_lum_vs[0] += (l_tot * Extnodes[p].rt_source_lum_vs[0]);
            rt_source_lum_vs[1] += (l_tot * Extnodes[p].rt_source_lum_vs[1]);
            rt_source_lum_vs[2] += (l_tot * Extnodes[p].rt_source_lum_vs[2]);
#endif
#ifdef BH_PHOTONMOMENTUM
            bh_lum += Nodes[p].bh_lum;
            bh_lum_grad[0] += Nodes[p].bh_lum * Nodes[p].bh_lum_grad[0];
            bh_lum_grad[1] += Nodes[p].bh_lum * Nodes[p].bh_lum_grad[1];
            bh_lum_grad[2] += Nodes[p].bh_lum * Nodes[p].bh_lum_grad[2];
#endif
#ifdef BH_CALC_DISTANCES
            bh_mass += Nodes[p].bh_mass;
            bh_pos_times_mass[0] += Nodes[p].bh_pos[0] * Nodes[p].bh_mass;
            bh_pos_times_mass[1] += Nodes[p].bh_pos[1] * Nodes[p].bh_mass;
            bh_pos_times_mass[2] += Nodes[p].bh_pos[2] * Nodes[p].bh_mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
            bh_mom[0] += Nodes[p].bh_vel[0] * Nodes[p].bh_mass;
            bh_mom[1] += Nodes[p].bh_vel[1] * Nodes[p].bh_mass;
            bh_mom[2] += Nodes[p].bh_vel[2] * Nodes[p].bh_mass;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
            if(Nodes[p].bh_mass > 0) {max_feedback_vel = DMAX(max_feedback_vel, Nodes[p].MaxFeedbackVel);}
#endif
            N_BH += Nodes[p].N_BH;
#endif
#endif
#ifdef DM_SCALARFIELD_SCREENING
            mass_dm += (Nodes[p].mass_dm);
            s_dm[0] += (Nodes[p].mass_dm * Nodes[p].s_dm[0]);
            s_dm[1] += (Nodes[p].mass_dm * Nodes[p].s_dm[1]);
            s_dm[2] += (Nodes[p].mass_dm * Nodes[p].s_dm[2]);
            vs_dm[0] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[0]);
            vs_dm[1] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[1]);
            vs_dm[2] += (Nodes[p].mass_dm * Extnodes[p].vs_dm[2]);
#endif
            vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
            vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
            vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);

            if(Extnodes[p].hmax > hmax)
                hmax = Extnodes[p].hmax;
            if(Extnodes[p].vmax > vmax)
                vmax = Extnodes[p].vmax;
            if(Extnodes[p].divVmax > divVmax)
                divVmax = Extnodes[p].divVmax;

            if(Nodes[p].u.d.mass > 0)
            {
                if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
                    count_particles += 2;
                else
                    count_particles++;
            }

            if(Nodes[p].maxsoft > maxsoft)
                maxsoft = Nodes[p].maxsoft;
        }
        else
            endrun(6767);		/* may not happen */

        p = Nodes[p].u.d.sibling;
    }

    if(mass)
    {
        s[0] /= mass;
        s[1] /= mass;
        s[2] /= mass;
        vs[0] /= mass;
        vs[1] /= mass;
        vs[2] /= mass;
    }
    else
    {
        s[0] = Nodes[no].center[0];
        s[1] = Nodes[no].center[1];
        s[2] = Nodes[no].center[2];
        vs[0] = 0;
        vs[1] = 0;
        vs[2] = 0;
    }

#ifdef RT_SEPARATELY_TRACK_LUMPOS
    double l_tot=0; int kfreq; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {l_tot += stellar_lum[kfreq];}
    if(l_tot)
    {
        rt_source_lum_s[0] /= l_tot;
        rt_source_lum_s[1] /= l_tot;
        rt_source_lum_s[2] /= l_tot;
        rt_source_lum_vs[0] /= l_tot;
        rt_source_lum_vs[1] /= l_tot;
        rt_source_lum_vs[2] /= l_tot;
    }
    else
    {
        rt_source_lum_s[0] = Nodes[no].center[0];
        rt_source_lum_s[1] = Nodes[no].center[1];
        rt_source_lum_s[2] = Nodes[no].center[2];
        rt_source_lum_vs[0] = 0;
        rt_source_lum_vs[1] = 0;
        rt_source_lum_vs[2] = 0;
    }
#endif
#ifdef BH_PHOTONMOMENTUM
    if(bh_lum)
    {
        bh_lum_grad[0] /= bh_lum; bh_lum_grad[1] /= bh_lum; bh_lum_grad[2] /= bh_lum;
    }
    else
    {
        bh_lum_grad[0]=bh_lum_grad[1]=0; bh_lum_grad[2]=1;
    }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    if(mass_dm)
    {
        s_dm[0] /= mass_dm;
        s_dm[1] /= mass_dm;
        s_dm[2] /= mass_dm;
        vs_dm[0] /= mass_dm;
        vs_dm[1] /= mass_dm;
        vs_dm[2] /= mass_dm;
    }
    else
    {
        s_dm[0] = Nodes[no].center[0];
        s_dm[1] = Nodes[no].center[1];
        s_dm[2] = Nodes[no].center[2];
        vs_dm[0] = 0;
        vs_dm[1] = 0;
        vs_dm[2] = 0;
    }
#endif


    Nodes[no].u.d.s[0] = s[0];
    Nodes[no].u.d.s[1] = s[1];
    Nodes[no].u.d.s[2] = s[2];
    Extnodes[no].vs[0] = vs[0];
    Extnodes[no].vs[1] = vs[1];
    Extnodes[no].vs[2] = vs[2];
    Nodes[no].u.d.mass = mass;
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
    Nodes[no].rt_source_lum_s[0] = rt_source_lum_s[0];
    Nodes[no].rt_source_lum_s[1] = rt_source_lum_s[1];
    Nodes[no].rt_source_lum_s[2] = rt_source_lum_s[2];
    Extnodes[no].rt_source_lum_vs[0] = rt_source_lum_vs[0];
    Extnodes[no].rt_source_lum_vs[1] = rt_source_lum_vs[1];
    Extnodes[no].rt_source_lum_vs[2] = rt_source_lum_vs[2];
#endif
#ifdef BH_PHOTONMOMENTUM
    Nodes[no].bh_lum = bh_lum;
    Nodes[no].bh_lum_grad[0] = bh_lum_grad[0];
    Nodes[no].bh_lum_grad[1] = bh_lum_grad[1];
    Nodes[no].bh_lum_grad[2] = bh_lum_grad[2];
#endif
#ifdef BH_CALC_DISTANCES
    Nodes[no].bh_mass = bh_mass;
    if(bh_mass > 0)
        {
            Nodes[no].bh_pos[0] = bh_pos_times_mass[0] / bh_mass;
            Nodes[no].bh_pos[1] = bh_pos_times_mass[1] / bh_mass;
            Nodes[no].bh_pos[2] = bh_pos_times_mass[2] / bh_mass;
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(SINGLE_STAR_FIND_BINARIES)
            Nodes[no].bh_vel[0] = bh_mom[0] / bh_mass;
            Nodes[no].bh_vel[1] = bh_mom[1] / bh_mass;
            Nodes[no].bh_vel[2] = bh_mom[2] / bh_mass;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
            Nodes[no].MaxFeedbackVel = max_feedback_vel;
#endif            
            Nodes[no].N_BH = N_BH;
#endif
        }
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Nodes[no].s_dm[0] = s_dm[0];
    Nodes[no].s_dm[1] = s_dm[1];
    Nodes[no].s_dm[2] = s_dm[2];
    Nodes[no].mass_dm = mass_dm;
    Extnodes[no].vs_dm[0] = vs_dm[0];
    Extnodes[no].vs_dm[1] = vs_dm[1];
    Extnodes[no].vs_dm[2] = vs_dm[2];
#endif

    Extnodes[no].hmax = hmax;
    Extnodes[no].vmax = vmax;
    Extnodes[no].divVmax = divVmax;
    Extnodes[no].Flag = GlobFlag;


    if(count_particles > 1)
        multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
    else
        multiple_flag = 0;

    Nodes[no].u.d.bitflags &= (~BITFLAG_MASK);	/* this clears the bits */
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
            if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
                break;

            Nodes[no].u.d.bitflags |= (1 << BITFLAG_TOPLEVEL);

            no = Nodes[no].u.d.father;
        }

        /* mark also internal top level nodes */

        no = DomainNodeIndex[i];
        no = Nodes[no].u.d.father;

        while(no >= 0)
        {
            if(Nodes[no].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
                break;

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

            if(DomainTask[i] != ThisTask)
                endrun(131231231);

            while(no >= 0)
            {
                if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))
                    break;

                Nodes[no].u.d.bitflags |= (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS);

                no = Nodes[no].u.d.father;
            }
        }
}


/*! When a new additional star particle is created, we can put it into the
 *  tree at the position of the spawning gas particle. This is possible
 *  because the Nextnode[] array essentially describes the full tree walk as a
 *  link list. Multipole moments of tree nodes need not be changed.
 */
void force_add_star_to_tree(int igas, int istar)
{
    int no;
    no = Nextnode[igas];
    Nextnode[igas] = istar; // insert new particle into linked list
    Nextnode[istar] = no; // order correctly
    Father[istar] = Father[igas]; // set parent node to be the same
    // update parent node properties [maximum softening, speed] for opening criteria
    Extnodes[Father[igas]].hmax = DMAX(Extnodes[Father[igas]].hmax, P[igas].Hsml);
    double vmax = Extnodes[Father[igas]].vmax;
    int k; for(k=0; k<3; k++) {if(fabs(P[istar].Vel[k]) > vmax) {vmax = fabs(P[istar].Vel[k]);}}
    Extnodes[Father[igas]].vmax = vmax;
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
    int no, nodesinlist, ptype, ninteractions, nexp, task, listindex = 0;
    double r2, dx, dy, dz, mass, r, fac, u, h=0, h_inv, h3_inv, xtmp; xtmp=0;
#ifdef RT_USE_TREECOL_FOR_NH
    double gasmass, angular_bin_size = 4*M_PI / RT_USE_TREECOL_FOR_NH, treecol_angular_bins[RT_USE_TREECOL_FOR_NH] = {0};
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
    double dvx, dvy, dvz;
    double jerk[3] = {0,0,0};
#endif
    double pos_x, pos_y, pos_z, aold;

#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(BH_DYNFRICTION_FROMTREE)
    double vel_x, vel_y, vel_z;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    double r_source, r_target, center[3]={0};
#ifdef BOX_PERIODIC
    center[0] = 0.5 * boxSize_X;
    center[1] = 0.5 * boxSize_Y;
    center[2] = 0.5 * boxSize_Z;
#endif    
#endif
#ifdef PMGRID
    int tabindex;
    double eff_dist, rcut, asmth, asmthfac, rcut2, dist;
    dist = 0;
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
    MyFloat tree_mass = 0;
#endif
    MyLongDouble acc_x, acc_y, acc_z;
    // cache some global vars in local vars to help compiler with alias analysis
    int maxPart = All.MaxPart;
    long bunchSize = All.BunchSize;
    int maxNodes = MaxNodes;
    integertime ti_Current = All.Ti_Current;
    double errTol2 = All.ErrTolTheta * All.ErrTolTheta;
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    int i1, i2; double fac2_tidal, fac_tidal; MyDouble tidal_tensorps[3][3];
#endif
#if defined(REDUCE_TREEWALK_BRANCHING) && defined(PMGRID)
    double dxx, dyy, dzz, pdxx, pdyy, pdzz;
#endif
#ifdef RT_USE_GRAVTREE
    double mass_stellarlum[N_RT_FREQ_BINS];
    int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {mass_stellarlum[k_freq]=0;}
#ifdef CHIMES_STELLAR_FLUXES
    double chimes_mass_stellarlum_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_mass_stellarlum_ion[CHIMES_LOCAL_UV_NBINS]={0};
    double chimes_flux_G0[CHIMES_LOCAL_UV_NBINS]={0}, chimes_flux_ion[CHIMES_LOCAL_UV_NBINS]={0};
#endif
    double dx_stellarlum=0, dy_stellarlum=0, dz_stellarlum=0;
    int valid_gas_particle_for_rt = 0;
#ifdef RT_OTVET
    double RT_ET[N_RT_FREQ_BINS][6]={{0}};
#endif
#endif

#ifdef BH_PHOTONMOMENTUM
    double mass_bhlum=0; // convert bh luminosity to our tree units
#endif
#ifdef BH_COMPTON_HEATING
    double incident_flux_agn=0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
    double Rad_E_gamma[N_RT_FREQ_BINS]={0};
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    double Rad_Flux[N_RT_FREQ_BINS][3]; {int kf,k2; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k2=0;k2<3;k2++) {Rad_Flux[kf][k2]=0;}}}
#endif

#ifdef BH_CALC_DISTANCES
    double min_dist_to_bh2=1.e37;
    double min_xyz_to_bh[3]={1.e37,1.e37,1.e37};
#ifdef SINGLE_STAR_FIND_BINARIES
    double min_bh_t_orbital=MAX_REAL_NUMBER, comp_dx[3], comp_dv[3], comp_Mass;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
    double min_bh_approach_time = MAX_REAL_NUMBER;
    double min_bh_freefall_time = MAX_REAL_NUMBER;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
    double min_bh_fb_time = MAX_REAL_NUMBER;
#endif
#endif
#endif


#ifdef DM_SCALARFIELD_SCREENING
    double dx_dm = 0, dy_dm = 0, dz_dm = 0, mass_dm = 0;
#endif
#if defined(BH_DYNFRICTION_FROMTREE)
    double bh_mass = 0;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    double soft=0, pmass;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
    double h_p_inv=0, h_p3_inv=0, u_p=0, zeta, zeta_sec=0;
    int ptype_sec=-1;
#endif
#endif
#ifdef EVALPOTENTIAL
    double facpot;
    MyLongDouble pot;
    pot = 0;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    for(i1 = 0; i1 < 3; i1++) {for(i2 = 0; i2 < 3; i2++) {tidal_tensorps[i1][i2] = 0.0;}}
#endif

    acc_x = 0;
    acc_y = 0;
    acc_z = 0;
    ninteractions = 0;
    nodesinlist = 0;

#ifdef PMGRID
    rcut = All.Rcut[0];
    asmth = All.Asmth[0];
    if(mode != 0 && mode != 1)
    {
        printf("%d %d %d %d %d\n", target, mode, *exportflag, *exportnodecount, *exportindex);
        endrun(444);
    }
#endif

    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        ptype = P[target].Type;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING)
        pmass = P[target].Mass;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(BH_DYNFRICTION_FROMTREE)
        vel_x = P[target].Vel[0];
        vel_y = P[target].Vel[1];
        vel_z = P[target].Vel[2];
#endif
#if defined(BH_DYNFRICTION_FROMTREE)
        if(ptype==5) {bh_mass = P[target].BH_Mass;}
#endif
        aold = All.ErrTolForceAcc * P[target].OldAcc;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        soft = All.ForceSoftening[ptype];
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
        if((ptype == 0) && (PPP[target].Hsml>All.ForceSoftening[ptype]))
        {
            soft = PPP[target].Hsml; zeta = PPPZ[target].AGS_zeta;
        } else {
            soft = All.ForceSoftening[ptype]; zeta = 0;
        }
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
        soft = PPP[target].AGS_Hsml;
        zeta = PPPZ[target].AGS_zeta;
#endif
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
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING)
        pmass = GravDataGet[target].Mass;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(BH_DYNFRICTION_FROMTREE)
        vel_x = GravDataGet[target].Vel[0];
        vel_y = GravDataGet[target].Vel[1];
        vel_z = GravDataGet[target].Vel[2];
#endif
        ptype = GravDataGet[target].Type;
#if defined(BH_DYNFRICTION_FROMTREE)
        if(ptype==5) {bh_mass = GravDataGet[target].BH_Mass;}
#endif
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        soft = GravDataGet[target].Soft;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
        zeta = GravDataGet[target].AGS_zeta;
#endif
#endif
#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
        if(pmforce_is_particle_high_res(ptype, GravDataGet[target].Pos))
        {
            rcut = All.Rcut[1];
            asmth = All.Asmth[1];
        }
#endif
    }
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
    /* quick check if particle has mass: if not, we won't deal with it */
    if(pmass<=0) return 0;
    int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(ptype); // determine allowed particle types for correction terms for adaptive gravitational softening terms
    int j0_sec_for_ags = -1;
#endif
#ifdef PMGRID
    rcut2 = rcut * rcut;
    asmthfac = 0.5 / asmth * (NTAB / 3.0);
#endif


#ifdef NEIGHBORS_MUST_BE_COMPUTED_EXPLICITLY_IN_FORCETREE
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
    double targeth_si = soft;
#else
    double targeth_si = All.ForceSoftening[ptype];
#endif
#endif




#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    h=soft;
#else
    if(ptype==0) {h=soft;} else {h=All.ForceSoftening[ptype];}
#endif
    h_inv = 1.0 / h;
    h3_inv = h_inv * h_inv * h_inv;
#endif


#ifdef RT_USE_GRAVTREE
    if(ptype==0) {if((soft>0)&&(pmass>0)) {valid_gas_particle_for_rt = 1;}}
#if defined(RT_LEBRON) && !defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
    double fac_stellum[N_RT_FREQ_BINS];
    if(valid_gas_particle_for_rt)
    {
        double h_eff_phys = soft * pow(NORM_COEFF/All.DesNumNgb,1./NUMDIMS) * All.cf_atime; // convert from softening kernel extent to effective size, assuming 3D here, and convert to physical code units
        double sigma_particle =  pmass / (h_eff_phys*h_eff_phys); // quick estimate of effective surface density of the target, in physical code units
        double fac_stellum_0 = -All.PhotonMomentum_Coupled_Fraction / (4.*M_PI * C_LIGHT_CODE_REDUCED * sigma_particle * All.G); // this will be multiplied by L/r^2 below, giving acceleration, then extra G because code thinks this is gravity, so put extra G here. everything is in -physical- code units here //
        int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {fac_stellum[kf] = fac_stellum_0*(1 - exp(-rt_kappa(-1,kf)*sigma_particle));} // rt_kappa is in physical code units, so sigma_eff_abs should be also -- approximate surface-density through particle (for checking if we enter optically-thick limit)
    }
#endif
#endif


#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
    double m_enc_in_rcrit = 0, r_for_total_menclosed = h;
    if(r_for_total_menclosed <= 0) {r_for_total_menclosed=All.ForceSoftening[ptype];}
    r_for_total_menclosed = DMAX( r_for_total_menclosed , 0.1/(UNIT_LENGTH_IN_KPC*All.cf_atime) ); /* set a baseline Rcrit_min, otherwise we get statistics that are very noisy */
#endif


    if(mode == 0)
    {
        no = maxPart;		/* root node */
    }
    else
    {
        nodesinlist++;
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;	/* open it */
    }

    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < maxPart)
            {
                /* the index of the node is the index of the particle */
                if(P[no].Ti_current != ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    drift_particle(no, ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }
                dx = P[no].Pos[0] - pos_x;
                dy = P[no].Pos[1] - pos_y;
                dz = P[no].Pos[2] - pos_z;
#ifdef GRAVITY_SPHERICAL_SYMMETRY
		r_source = sqrt(pow(P[no].Pos[0] - center[0],2) + pow(P[no].Pos[1] - center[1],2) + pow(P[no].Pos[2] - center[2],2));
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(BH_DYNFRICTION_FROMTREE)
                dvx = P[no].Vel[0] - vel_x;
                dvy = P[no].Vel[1] - vel_y;
                dvz = P[no].Vel[2] - vel_z;
#endif
                GRAVITY_NEAREST_XYZ(dx,dy,dz,-1);
                r2 = dx * dx + dy * dy + dz * dz;
                mass = P[no].Mass;
#ifdef RT_USE_TREECOL_FOR_NH
                if(P[no].Type == 0) gasmass = P[no].Mass;
#ifdef BH_ALPHADISK_ACCRETION
                if(P[no].Type == 5) gasmass = BPP(no).BH_Mass_AlphaDisk; // gas at the inner edge of a disk should not see a hole due to the sink
#endif
#endif
                /* only proceed if the mass is positive and there is separation! */
                if((r2 > 0) && (mass > 0))
                {

#ifdef BH_CALC_DISTANCES
                if(P[no].Type == 5)             /* found a BH particle in grav calc */
                {
                    if(r2 < min_dist_to_bh2)    /* is this the closest BH part I've found yet? */
                    {
                        min_dist_to_bh2 = r2;   /* if yes: adjust min bh dist */
                        min_xyz_to_bh[0] = dx;  /* remember, dx = x_BH - myx */
                        min_xyz_to_bh[1] = dy;
                        min_xyz_to_bh[2] = dz;
                    }
#ifdef SINGLE_STAR_TIMESTEPPING
                    double bh_dvx=P[no].Vel[0]-vel_x, bh_dvy=P[no].Vel[1]-vel_y, bh_dvz=P[no].Vel[2]-vel_z, vSqr=bh_dvx*bh_dvx+bh_dvy*bh_dvy+bh_dvz*bh_dvz, M_total=P[no].Mass+pmass, r2soft=All.ForceSoftening[5];
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
                    r2soft = DMAX(r2soft, soft);
#endif                    
                    r2soft *= KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER;
                    r2soft = r2 + r2soft*r2soft;
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                    if(ptype == 0) {
                        double tSqr_fb = r2soft /(P[no].MaxFeedbackVel * P[no].MaxFeedbackVel + MIN_REAL_NUMBER);
                        if(tSqr_fb < min_bh_fb_time) {min_bh_fb_time = tSqr_fb;}
                    } // for gas, add the signal velocity of feedback from the star
#endif                    
                    double tSqr = r2soft/(vSqr + MIN_REAL_NUMBER), tff4 = r2soft*r2soft*r2soft/(M_total*M_total);

                    if(tSqr < min_bh_approach_time) {min_bh_approach_time = tSqr;}
                    if(tff4 < min_bh_freefall_time) {min_bh_freefall_time = tff4;}
#ifdef SINGLE_STAR_FIND_BINARIES
                    if(ptype == 5) // only for BH particles and for non center of mass calculation
                    {
                        double r_p5=sqrt(r2), specific_energy = 0.5*vSqr - All.G*M_total/r_p5;
                        if(r2 < All.ForceSoftening[5]*All.ForceSoftening[5])
                        {
                            double hinv_p5 = 1/All.ForceSoftening[5];
                            specific_energy = 0.5*vSqr + All.G*M_total*kernel_gravity(r_p5*hinv_p5, hinv_p5, hinv_p5*hinv_p5*hinv_p5, -1);
                        }
                        if (specific_energy < 0)
                        {
                            double semimajor_axis= -All.G*M_total/(2.*specific_energy);
                            double t_orbital = 2.*M_PI*sqrt( semimajor_axis*semimajor_axis*semimajor_axis / (All.G*M_total) );
                            if(t_orbital < min_bh_t_orbital) /* Save parameters of companion */
                            {
                                min_bh_t_orbital=t_orbital; comp_Mass=P[no].Mass;
                                comp_dx[0]=dx; comp_dx[1]=dy; comp_dx[2]=dz; comp_dv[0]=bh_dvx; comp_dv[1]=bh_dvy; comp_dv[2]=bh_dvz;
                            }
                        } /* specific_energy < 0 */
                    } /* ptype == 5 */
#endif //#ifdef SINGLE_STAR_FIND_BINARIES
#endif //#ifdef SINGLE_STAR_TIMESTEPPING
                }
#endif

                    
#ifdef RT_USE_GRAVTREE
                if(valid_gas_particle_for_rt)	/* we have a (valid) gas particle as target */
                {
                    dx_stellarlum=dx; dy_stellarlum=dy; dz_stellarlum=dz;
                    double lum[N_RT_FREQ_BINS];
#ifdef CHIMES_STELLAR_FLUXES
                    double chimes_lum_G0[CHIMES_LOCAL_UV_NBINS];
                    double chimes_lum_ion[CHIMES_LOCAL_UV_NBINS];
                    int active_check = rt_get_source_luminosity_chimes(no,1,lum, chimes_lum_G0, chimes_lum_ion);
#else
                    int active_check = rt_get_source_luminosity(no,1,lum);
#endif
                    int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(active_check) {mass_stellarlum[kf]=lum[kf];} else {mass_stellarlum[kf]=0;}}
#ifdef CHIMES_STELLAR_FLUXES
                    for (kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++)
                    {
                        if(active_check) {chimes_mass_stellarlum_G0[kf] = chimes_lum_G0[kf]; chimes_mass_stellarlum_ion[kf] = chimes_lum_ion[kf];} else {chimes_mass_stellarlum_G0[kf] = 0; chimes_mass_stellarlum_ion[kf] = 0;}
                    }
#endif
#ifdef BH_PHOTONMOMENTUM
                    mass_bhlum=0;
		            if(P[no].Type == 5)
		            {
			            double bhlum_t = bh_lum_bol(P[no].BH_Mdot, P[no].BH_Mass, no);
#if defined(BH_FOLLOW_ACCRETED_ANGMOM)
                        mass_bhlum = bh_angleweight(bhlum_t, P[no].BH_Specific_AngMom, dx,dy,dz);
#else
			            mass_bhlum = bh_angleweight(bhlum_t, P[no].GradRho, dx,dy,dz);
#endif
		            }
#endif
                }
#endif


#ifdef DM_SCALARFIELD_SCREENING
                if(ptype != 0)	/* we have a dark matter particle as target */
                {
                    if(P[no].Type == 1)
                    {
                        dx_dm = dx;
                        dy_dm = dy;
                        dz_dm = dz;
                        mass_dm = mass;
                    }
                    else
                    {
                        mass_dm = 0;
                        dx_dm = dy_dm = dz_dm = 0;
                    }
                }
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL) /* set secondary softening and zeta term */
                    ptype_sec = P[no].Type; j0_sec_for_ags = no;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                    if(ptype_sec == 0) {h_p_inv=1./PPP[no].Hsml; zeta_sec=PPPZ[no].AGS_zeta;} else {h_p_inv=1./All.ForceSoftening[P[no].Type]; zeta_sec=0;}
#else
                    h_p_inv=1./PPP[no].AGS_Hsml; zeta_sec=PPPZ[no].AGS_zeta;
#endif
#else
                    h = All.ForceSoftening[ptype];
                    if(h < All.ForceSoftening[P[no].Type]) {h = All.ForceSoftening[P[no].Type];}
#endif
                } // closes (if((r2 > 0) && (mass > 0))) check

                if(TakeLevel >= 0) {P[no].GravCost[TakeLevel] += 1.0;}
                no = Nextnode[no];
            }
            else			/* we have an  internal node */
            {
                if(no >= maxPart + maxNodes)	/* pseudo particle */
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
                            LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
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
                            UNLOCK_NEXPORT;
                            if(exitFlag) {return -1;} /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */

                            exportnodecount[task] = 0;
                            exportindex[task] = nexp;
                            DataIndexTable[nexp].Task = task;
                            DataIndexTable[nexp].Index = target;
                            DataIndexTable[nexp].IndexGet = nexp;
                        }

                        DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
                        DomainNodeIndex[no - (maxPart + maxNodes)];

                        if(exportnodecount[task] < NODELISTLENGTH)
                            DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
                    }
                    no = Nextnode[no - maxNodes];
                    continue;
                }

                nop = &Nodes[no];

                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }

                mass = nop->u.d.mass;
#ifdef RT_USE_TREECOL_FOR_NH
                gasmass = nop->gasmass;
#endif
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* open cell */
                    if(mass)
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }

                if(nop->Ti_current != ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    force_drift_node(no, ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }

                dx = nop->u.d.s[0] - pos_x;
                dy = nop->u.d.s[1] - pos_y;
                dz = nop->u.d.s[2] - pos_z;
#ifdef GRAVITY_SPHERICAL_SYMMETRY
		r_source = sqrt(pow(nop->u.d.s[0] - center[0],2) + pow(nop->u.d.s[1] - center[1],2) + pow(nop->u.d.s[2] - center[2],2));
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(BH_DYNFRICTION_FROMTREE)
                dvx = Extnodes[no].vs[0] - vel_x;
                dvy = Extnodes[no].vs[1] - vel_y;
                dvz = Extnodes[no].vs[2] - vel_z;
#endif
                GRAVITY_NEAREST_XYZ(dx,dy,dz,-1);
                r2 = dx * dx + dy * dy + dz * dz;

                
#ifdef RT_USE_GRAVTREE
                if(valid_gas_particle_for_rt)	/* we have a (valid) gas particle as target */
                {
                    int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {mass_stellarlum[kf] = nop->stellar_lum[kf];}
#ifdef CHIMES_STELLAR_FLUXES
                    for (kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++)
                    {
                        chimes_mass_stellarlum_G0[kf] = nop->chimes_stellar_lum_G0[kf];
                        chimes_mass_stellarlum_ion[kf] = nop->chimes_stellar_lum_ion[kf];
                    }
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
                    dx_stellarlum = nop->rt_source_lum_s[0] - pos_x; dy_stellarlum = nop->rt_source_lum_s[1] - pos_y; dz_stellarlum = nop->rt_source_lum_s[2] - pos_z;
                    GRAVITY_NEAREST_XYZ(dx_stellarlum,dy_stellarlum,dz_stellarlum,-1);
#else
                    dx_stellarlum = dx; dy_stellarlum = dy; dz_stellarlum = dz;
#endif
#ifdef BH_PHOTONMOMENTUM
                    mass_bhlum = bh_angleweight(nop->bh_lum, nop->bh_lum_grad, dx_stellarlum,dy_stellarlum,dz_stellarlum);
#endif
                }
#endif


#ifdef DM_SCALARFIELD_SCREENING
                if(ptype != 0)	/* we have a dark matter particle as target */
                {
                    dx_dm = nop->s_dm[0] - pos_x;
                    dy_dm = nop->s_dm[1] - pos_y;
                    dz_dm = nop->s_dm[2] - pos_z;
                    mass_dm = nop->mass_dm;
                }
                else
                {
                    mass_dm = 0;
                    dx_dm = dy_dm = dz_dm = 0;
                }
#endif

#ifdef PMGRID
#ifdef REDUCE_TREEWALK_BRANCHING
                dxx = (nop->center[0] - pos_x);
                dyy = (nop->center[1] - pos_y);
                dzz = (nop->center[2] - pos_z);
                eff_dist = rcut + 0.5 * nop->len;
                pdxx = GRAVITY_NGB_PERIODIC_BOX_LONG_X(dxx,dyy,dzz,-1);
                pdyy = GRAVITY_NGB_PERIODIC_BOX_LONG_Y(dxx,dyy,dzz,-1);
                pdzz = GRAVITY_NGB_PERIODIC_BOX_LONG_Z(dxx,dyy,dzz,-1);
                /* check whether we can stop walking along this branch */
                if((r2 > rcut2) & ((pdxx > eff_dist) | (pdyy > eff_dist) | (pdzz > eff_dist)))
                {
                    no = nop->u.d.sibling;
                    continue;
                }
#else
                /* check whether we can stop walking along this branch */
                if(r2 > rcut2)
                {
                    eff_dist = rcut + 0.5 * nop->len;
                    dist = GRAVITY_NGB_PERIODIC_BOX_LONG_X(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1);
                    if(dist > eff_dist)
                    {
                        no = nop->u.d.sibling;
                        continue;
                    }
                    dist = GRAVITY_NGB_PERIODIC_BOX_LONG_Y(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1);
                    if(dist > eff_dist)
                    {
                        no = nop->u.d.sibling;
                        continue;
                    }
                    dist = GRAVITY_NGB_PERIODIC_BOX_LONG_Z(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1);
                    if(dist > eff_dist)
                    {
                        no = nop->u.d.sibling;
                        continue;
                    }
                }
#endif
#endif // PMGRID //


#ifdef NEIGHBORS_MUST_BE_COMPUTED_EXPLICITLY_IN_FORCETREE
                {
                    double dx_nc = nop->center[0] - pos_x;
                    double dy_nc = nop->center[1] - pos_y;
                    double dz_nc = nop->center[2] - pos_z;
                    GRAVITY_NEAREST_XYZ(dx_nc,dy_nc,dz_nc,-1); /* find the closest image in the given box size  */
                    double dist_to_center2 = dx_nc*dx_nc +  dy_nc*dy_nc + dz_nc*dz_nc;
                    /* check if any portion the cell lies within the interaction range */
		            double dist_to_open = 2.0*targeth_si + nop->len*1.73205/2.0;
                    if(dist_to_center2  < dist_to_open*dist_to_open)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#endif


                if(errTol2)	/* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * errTol2)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#ifndef GRAVITY_HYBRID_OPENING_CRIT
                else		/* check relative opening criterion */
#else
                if(!(All.Ti_Current == 0 && RestartFlag != 1))
#endif
                {
                    /* force node to open if we are within the gravitational softening length */
#if !(defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(FLAG_NOT_IN_PUBLIC_CODE))
                    double soft = All.ForceSoftening[ptype];
#endif
                    if((r2 < (soft+0.6*nop->len)*(soft+0.6*nop->len)) || (r2 < (nop->maxsoft+0.6*nop->len)*(nop->maxsoft+0.6*nop->len)))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

#if defined(REDUCE_TREEWALK_BRANCHING) && defined(PMGRID)
                    if((mass * nop->len * nop->len > r2 * r2 * aold) |
                       ((pdxx < 0.60 * nop->len) & (pdyy < 0.60 * nop->len) & (pdzz < 0.60 * nop->len)))
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
#else
                    if(mass * nop->len * nop->len > r2 * r2 * aold)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }

                    /* check in addition whether we lie inside the cell */

                    if(GRAVITY_NGB_PERIODIC_BOX_LONG_X(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1) < 0.60 * nop->len)
                    {
                        if(GRAVITY_NGB_PERIODIC_BOX_LONG_Y(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1) < 0.60 * nop->len)
                        {
                            if(GRAVITY_NGB_PERIODIC_BOX_LONG_Z(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1) < 0.60 * nop->len)
                            {
                                no = nop->u.d.nextnode;
                                continue;
                            }
                        }
                    }
#endif
                }

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                /* set secondary softening and zeta term */
                if(nop->maxsoft > 0) {h_p_inv = 1.0 / nop->maxsoft;} else {h_p_inv = 0;}
                zeta_sec = 0; ptype_sec = -1; j0_sec_for_ags = -1;
                if(h < nop->maxsoft) // compare primary softening to node maximum
                {
                    if(r2 < nop->maxsoft * nop->maxsoft) // inside node maxsoft! continue down tree
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#else
                h = All.ForceSoftening[ptype];
                if(h < nop->maxsoft)
                {
                    h = nop->maxsoft;
                    if(r2 < h * h)
                    {
                        if(maskout_different_softening_flag(nop->u.d.bitflags))	/* bit-5 signals that there are particles of different softening in the node */
                        {
                            no = nop->u.d.nextnode;
                            continue;
                        }
                    }
                }
#endif

                if(TakeLevel >= 0) {nop->GravCost += 1.0;}
                no = nop->u.d.sibling;	/* ok, node can be used */

#ifdef BH_CALC_DISTANCES // NOTE: moved this to AFTER the checks for node opening, because we only want to record BH positions from the nodes that actually get used for the force calculation - MYG
                if(nop->bh_mass > 0)        /* found a node with non-zero BH mass */
                {
                    double bh_dx = nop->bh_pos[0] - pos_x;      /* SHEA:  now using bh_pos instead of center */
                    double bh_dy = nop->bh_pos[1] - pos_y;
                    double bh_dz = nop->bh_pos[2] - pos_z;
                    GRAVITY_NEAREST_XYZ(bh_dx,bh_dy,bh_dz,-1);
                    double bh_r2 = bh_dx * bh_dx + bh_dy * bh_dy + bh_dz * bh_dz; // + (nop->len)*(nop->len);
                    if(bh_r2 < min_dist_to_bh2)
                    {
                        min_dist_to_bh2 = bh_r2;
                        min_xyz_to_bh[0] = bh_dx;    /* remember, dx = x_BH - myx */
                        min_xyz_to_bh[1] = bh_dy;
                        min_xyz_to_bh[2] = bh_dz;
                    }
#ifdef SINGLE_STAR_TIMESTEPPING
                    double bh_dvx=nop->bh_vel[0]-vel_x, bh_dvy=nop->bh_vel[1]-vel_y, bh_dvz=nop->bh_vel[2]-vel_z, vSqr=bh_dvx*bh_dvx+bh_dvy*bh_dvy+bh_dvz*bh_dvz, M_total=nop->bh_mass+pmass, r2soft;
                    r2soft = DMAX(All.ForceSoftening[5], soft) * KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER;
                    r2soft = r2 + r2soft*r2soft;
                    double tSqr = r2soft/(vSqr + MIN_REAL_NUMBER), tff4 = r2soft*r2soft*r2soft/(M_total*M_total);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                    if(ptype == 0) {
                        double tSqr_fb = r2soft /(nop->MaxFeedbackVel * nop->MaxFeedbackVel + MIN_REAL_NUMBER);
                        if(tSqr_fb < min_bh_fb_time) {min_bh_fb_time = tSqr_fb;}
                    } // for gas, add the signal velocity of feedback from the star
#endif                                                            
                    if(tSqr < min_bh_approach_time) {min_bh_approach_time = tSqr;}
                    if(tff4 < min_bh_freefall_time) {min_bh_freefall_time = tff4;}
#ifdef SINGLE_STAR_FIND_BINARIES
                    if(ptype == 5 && nop->N_BH == 1) // only do it if we're looking at a single star in the node
                    {
                        double specific_energy = 0.5*vSqr - All.G*M_total/sqrt(r2);
                        if (specific_energy<0)
                        {
                            double semimajor_axis= -All.G*M_total/(2.*specific_energy);
                            double t_orbital = 2.*M_PI*sqrt( semimajor_axis*semimajor_axis*semimajor_axis / (All.G*M_total) );
                            if(t_orbital < min_bh_t_orbital) /* Save parameters of companion */
                            {
                                min_bh_t_orbital=t_orbital; comp_Mass=nop->bh_mass;
                                comp_dx[0]=bh_dx; comp_dx[1]=bh_dy; comp_dx[2]=bh_dz; comp_dv[0]=bh_dvx; comp_dv[1]=bh_dvy; comp_dv[2]=bh_dvz;
                            }
                        } /* specific_energy < 0 */
                    } /* ptype == 5 */
#endif //#ifdef SINGLE_STAR_FIND_BINARIES
#endif //#ifdef SINGLE_STAR_TIMESTEPPING
                }
#endif
            }
	    	   
            if((r2 > 0) && (mass > 0)) // only go forward if mass positive and there is separation
            {
            r = sqrt(r2);
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
            if((r >= h) && !((ptype_sec > -1) && (r < 1/h_p_inv))) // can only do the Newtonian force if the field source is outside our own softening, and we are not within the softening of a field source particle
#else
            if(r >= h)
#endif
            {
                fac = mass / (r2 * r);
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                fac2_tidal = 3.0 * mass / (r2 * r2 * r); /* second derivative of potential needs this factor */
#endif
#ifdef EVALPOTENTIAL
                facpot = -mass / r;
#endif
            }
            else
            {
#if !defined(ADAPTIVE_GRAVSOFT_FORALL) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
                h_inv = 1.0 / h;
                h3_inv = h_inv * h_inv * h_inv;
#endif
                u = r * h_inv;
                fac = mass * kernel_gravity(u, h_inv, h3_inv, 1);
#ifdef EVALPOTENTIAL
                facpot = mass * kernel_gravity(u, h_inv, h3_inv, -1);
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* second derivatives needed -> calculate them from softened potential. NOTE this is here -assuming- a cubic spline, will be inconsistent for different kernels used! */
                fac2_tidal = mass * kernel_gravity(u, h_inv, h3_inv, 2);
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                // first, appropriately symmetrize the forces between particles //
                if((h_p_inv > 0) && (ptype_sec > -1))
                {
                    int symmetrize_by_averaging = 0;
#ifdef ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING // the 'zeta' terms for conservation with adaptive softening assume kernel-scale forces are averaged to symmetrize, to make them continuous
                    //symmetrize_by_averaging = 1; // always symmetrize by averaging //
                    if((1 << ptype_sec) & (AGS_kernel_shared_BITFLAG)) {symmetrize_by_averaging=1;} // symmetrize by averaging only for particles which have a shared AGS structure since this is how our correction terms are derived //
#ifdef SINGLE_STAR_SINK_DYNAMICS
                    if((ptype!=0) || (ptype_sec!=0)) {symmetrize_by_averaging=0;} // we don't want to do the symmetrization below for sink interactions because it can create very noisy interactions between tiny sink particles and diffuse gas. However we do want it for gas-gas interactions so we keep the below
#endif
#endif
                    if(symmetrize_by_averaging==1)
                    {
                    h_p3_inv = h_p_inv * h_p_inv * h_p_inv; u_p = r * h_p_inv;
                    fac = 0.5 * (fac + mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 1)); // average with neighbor
#ifdef EVALPOTENTIAL
                    facpot = 0.5 * (facpot + mass * kernel_gravity(u, h_p_inv, h_p3_inv, -1)); // average with neighbor
#endif
#if defined(COMPUTE_TIDAL_TENSOR_IN_GRAVTREE)
                    fac2_tidal = 0.5 * (fac2_tidal + mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 2)); // average forces -> average in tidal tensor as well
#endif
                    // correction only applies to 'shared-kernel' particles: so this needs to check if
                    // these are the same particles for which the 'shared' kernel lengths are computed
                    if((1 << ptype_sec) & (AGS_kernel_shared_BITFLAG))
                    {
                        if((r>0) && (pmass>0)) // checks that these aren't the same particle or test particle
                        {
                            double dWdr, wp, fac_corr = 0;
                            if((zeta != 0) && (u < 1)) // in kernel [zeta non-zero]
                            {
                                kernel_main(u, h3_inv, h3_inv*h_inv, &wp, &dWdr, 1);
                                fac_corr += -(zeta/pmass) * dWdr / r;   // 0.5 * zeta * omega * dWdr / r;
                            }
                            if((zeta_sec != 0) && (u_p < 1)) // in kernel [ zeta non-zero]
                            {
                                kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wp, &dWdr, 1);
                                fac_corr += -(zeta_sec/pmass) * dWdr / r;   // 0.5 * zeta * omega * dWdr / r;
                            }
                            if(!isnan(fac_corr)) {fac += fac_corr;}
                        }
                    } // if(ptype==ptype_sec)
                    } // closes block for symmetrizing forces by averaging //
                    else
                    { // open block to symmetrize instead with the old method of simply taking the larger of the pair //
                    if(h_p_inv < h_inv) // if the softening of the particle whose force is being summed is greater than the target
                    {
                        h_p3_inv = h_p_inv * h_p_inv * h_p_inv; u_p = r * h_p_inv;
                        fac = mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 1);
#ifdef EVALPOTENTIAL
                        facpot = mass * kernel_gravity(u, h_p_inv, h_p3_inv, -1);
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                        fac2_tidal = mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 2);
#endif
                    }
                    // correction only applies to 'shared-kernel' particles: so this needs to check if
                    // these are the same particles for which the 'shared' kernel lengths are computed
                    if((1 << ptype_sec) & (AGS_kernel_shared_BITFLAG))
                    {
                        if((r>0) && (pmass>0)) // checks that these aren't the same particle or test particle
                        {
                            double dWdr, wp, fac_corr=0;
                            if(h_p_inv >= h_inv)
                            {
                                if((zeta != 0) && (u < 1)) // in kernel [zeta non-zero]
                                {
                                    kernel_main(u, h3_inv, h3_inv*h_inv, &wp, &dWdr, 1);
                                    fac_corr += -2. * (zeta/pmass) * dWdr / r;   // 0.5 * zeta * omega * dWdr / r;
                                }
                            } else {
                                if((zeta_sec != 0) && (u_p < 1)) // in kernel [ zeta non-zero]
                                {
                                    kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wp, &dWdr, 1);
                                    fac_corr += -2. * (zeta_sec/pmass) * dWdr / r;
                                }
                            }
                            if(!isnan(fac_corr)) {fac += fac_corr;}
                        }
                    } // if(ptype==ptype_sec)
                    } // closes else{} block for when take the larger of two softenings or symmetrize by averaging
                } // closes (if((h_p_inv > 0) && (ptype_sec > -1)))
#endif // #if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL) //
            } // closes r < h (else) clause


#ifdef PMGRID
            tabindex = (int) (asmthfac * r);
            if(tabindex < NTAB && tabindex >= 0)
#endif // PMGRID //
            {
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                fac_tidal = fac; /* save original fac without shortrange_table factor (needed for tidal field calculation) */
#endif

#ifdef PMGRID
                fac *= shortrange_table[tabindex];
#endif

#ifdef EVALPOTENTIAL
#ifdef PMGRID
                facpot *= shortrange_table_potential[tabindex];
#endif
                pot += FLT(facpot);
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
                pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
		r_target = sqrt(pow(pos_x - center[0],2) + pow(pos_y - center[1],2) + pow(pos_z - center[2],2)); // distance of target point from box center
		if(r_source < r_target){
		    dx = center[0] - pos_x;
   		    dy = center[1] - pos_y;
		    dz = center[2] - pos_z;
		    fac = mass/pow(DMAX(GRAVITY_SPHERICAL_SYMMETRY,DMAX(r_target,h)),3);
		} else {
           	    fac = 0;
		}
#endif
                acc_x += FLT(dx * fac);
                acc_y += FLT(dy * fac);
                acc_z += FLT(dz * fac);


#if defined(BH_DYNFRICTION_FROMTREE)
                if(ptype==5)
                {
                    double dv2=dvx*dvx+dvy*dvy+dvz*dvz;
                    if(dv2 > 0)
                    {
                        double dv0=sqrt(dv2),dvx_h=dvx/dv0,dvy_h=dvy/dv0,dvz_h=dvz/dv0,rdotvhat=dx*dvx_h+dy*dvy_h+dz*dvz_h;
                        double bx_im=dx-rdotvhat*dvx_h,by_im=dy-rdotvhat*dvy_h,bz_im=dz-rdotvhat*dvz_h,b_impact=sqrt(bx_im*bx_im+by_im*by_im+bz_im*bz_im);
                        double a_im=(b_impact*All.cf_atime)*(dv2*All.cf_a2inv)/(All.G*bh_mass), fac_df=fac*b_impact*a_im/(1.+a_im*a_im); // need to convert to fully-physical units to ensure this has the correct dimensions
                        acc_x += fac_df * dvx_h;
                        acc_y += fac_df * dvy_h;
                        acc_z += fac_df * dvz_h;
                    }
                }
#endif


#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                /*
                 tidal_tensorps[][] = Matrix of second derivatives of grav. potential, symmetric:
                 |Txx Txy Txz|   |tidal_tensorps[0][0] tidal_tensorps[0][1] tidal_tensorps[0][2]|
                 |Tyx Tyy Tyz| = |tidal_tensorps[1][0] tidal_tensorps[1][1] tidal_tensorps[1][2]|
                 |Tzx Tzy Tzz|   |tidal_tensorps[2][0] tidal_tensorps[2][1] tidal_tensorps[2][2]|
                 */
#ifdef GRAVITY_SPHERICAL_SYMMETRY
		if(r_source < r_target){
		    fac2_tidal = 3 * mass / pow(DMAX(GRAVITY_SPHERICAL_SYMMETRY,DMAX(r_target,h)),5);
		} else {
   		    fac2_tidal = 0;
		}
#endif
#ifdef PMGRID
                tidal_tensorps[0][0] += ((-fac_tidal + dx * dx * fac2_tidal) * shortrange_table[tabindex]) +
                    dx * dx * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                tidal_tensorps[0][1] += ((dx * dy * fac2_tidal) * shortrange_table[tabindex]) +
                    dx * dy * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                tidal_tensorps[0][2] += ((dx * dz * fac2_tidal) * shortrange_table[tabindex]) +
                    dx * dz * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                tidal_tensorps[1][1] += ((-fac_tidal + dy * dy * fac2_tidal) * shortrange_table[tabindex]) +
                    dy * dy * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                tidal_tensorps[1][2] += ((dy * dz * fac2_tidal) * shortrange_table[tabindex]) +
                    dy * dz * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
                tidal_tensorps[2][2] += ((-fac_tidal + dz * dz * fac2_tidal) * shortrange_table[tabindex]) +
                    dz * dz * fac2_tidal / 3.0 * shortrange_table_tidal[tabindex];
#else
                tidal_tensorps[0][0] += (-fac_tidal + dx * dx * fac2_tidal);
                tidal_tensorps[0][1] += (dx * dy * fac2_tidal);
                tidal_tensorps[0][2] += (dx * dz * fac2_tidal);
                tidal_tensorps[1][1] += (-fac_tidal + dy * dy * fac2_tidal);
                tidal_tensorps[1][2] += (dy * dz * fac2_tidal);
                tidal_tensorps[2][2] += (-fac_tidal + dz * dz * fac2_tidal);
#endif
                tidal_tensorps[1][0] = tidal_tensorps[0][1];
                tidal_tensorps[2][0] = tidal_tensorps[0][2];
                tidal_tensorps[2][1] = tidal_tensorps[1][2];
#endif // COMPUTE_TIDAL_TENSOR_IN_GRAVTREE //
#ifdef COMPUTE_JERK_IN_GRAVTREE
#ifndef ADAPTIVE_TREEFORCE_UPDATE // we want the jerk if we're doing lazy force updates
		if(ptype > 0)
#endif                    
                {
		    double dv_dot_dx = dx*dvx + dy*dvy + dz*dvz;
		    jerk[0] += dvx * fac - dv_dot_dx * fac2_tidal * dx;
		    jerk[1] += dvy * fac - dv_dot_dx * fac2_tidal * dy;
		    jerk[2] += dvz * fac - dv_dot_dx * fac2_tidal * dz;
		}
#endif
            } // closes TABINDEX<NTAB

            ninteractions++;


#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
            if(r < r_for_total_menclosed) {m_enc_in_rcrit += mass;}
#endif

#ifdef COUNT_MASS_IN_GRAVTREE
            tree_mass += mass;
#endif
#ifdef RT_USE_TREECOL_FOR_NH
            if(gasmass>0){
                int bin; // Here we do a simple six-bin angular binning scheme
                if ((fabs(dx) > fabs(dy)) && (fabs(dx)>fabs(dz))){if (dx > 0) {bin = 0;} else {bin=1;}
                } else if (fabs(dy)>fabs(dz)){if (dy > 0) {bin = 2;} else {bin=3;}
                } else {if (dz > 0) {bin = 4;} else {bin = 5;}}
                treecol_angular_bins[bin] += fac*gasmass*r / (angular_bin_size*mass); // in our binning scheme, we stretch the gas mass over a patch */ of the sphere located at radius r subtending solid angle equal to the bin size - thus the area is r^2 * angular_bin_size, so sigma = m/(r^2 * angular bin size) = fac/r / angular bin size. Factor of gasmass / mass corrects the gravitational mass to the gas mass
            }
#endif
#ifdef RT_USE_GRAVTREE
            if(valid_gas_particle_for_rt)	/* we have a (valid) gas particle as target */
            {
                r2 = dx_stellarlum*dx_stellarlum + dy_stellarlum*dy_stellarlum + dz_stellarlum*dz_stellarlum; r = sqrt(r2);
                if(r >= soft) {fac=1./(r2*r);} else {h_inv=1./soft; h3_inv=h_inv*h_inv*h_inv; u=r*h_inv; fac=kernel_gravity(u,h_inv,h3_inv,1);}
                if((soft>r)&&(soft>0)) fac *= (r2/(soft*soft)); // don't allow cross-section > r2
                double fac_intensity; fac_intensity = fac * r * All.cf_a2inv / (4.*M_PI); // ~L/(4pi*r^2), in -physical- units, since L is physical
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {Rad_E_gamma[kf] += fac_intensity * mass_stellarlum[kf];}}
#endif
#ifdef CHIMES_STELLAR_FLUXES
                int chimes_k;
                double chimes_fac = fac_intensity / (UNIT_LENGTH_IN_CGS*UNIT_LENGTH_IN_CGS);  // 1/(4 * pi * r^2), in cm^-2
                for (chimes_k = 0; chimes_k < CHIMES_LOCAL_UV_NBINS; chimes_k++)
                {
                    chimes_flux_G0[chimes_k] += chimes_fac * chimes_mass_stellarlum_G0[chimes_k];   // Habing flux units
                    chimes_flux_ion[chimes_k] += chimes_fac * chimes_mass_stellarlum_ion[chimes_k]; // cm^-2 s^-1
                }
#endif 
#ifdef BH_PHOTONMOMENTUM
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                Rad_E_gamma[RT_FREQ_BIN_FIRE_IR] += fac_intensity * mass_bhlum;
#endif
#ifdef BH_COMPTON_HEATING
                incident_flux_agn += fac_intensity * mass_bhlum; // L/(4pi*r*r) analog
#endif
#endif

#ifdef RT_OTVET
                /* use the information we have here from the gravity tree (optically thin incident fluxes) to estimate the Eddington tensor */
                // for now, just one tensor; so we use the sum of luminosities to determine the weights in the Eddington tensor
                if(r>0)
                {
                    double fac_sum=0;
                    int kf_rt;
                    for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++)
                    {
                        fac_sum = mass_stellarlum[kf_rt];
                        fac_sum *= fac / (1.e-37 + r); // units are not important, since ET will be dimensionless, but final ET should scale as ~luminosity/r^2
                        RT_ET[kf_rt][0] += dx_stellarlum * dx_stellarlum * fac_sum;
                        RT_ET[kf_rt][1] += dy_stellarlum * dy_stellarlum * fac_sum;
                        RT_ET[kf_rt][2] += dz_stellarlum * dz_stellarlum * fac_sum;
                        RT_ET[kf_rt][3] += dx_stellarlum * dy_stellarlum * fac_sum;
                        RT_ET[kf_rt][4] += dy_stellarlum * dz_stellarlum * fac_sum;
                        RT_ET[kf_rt][5] += dz_stellarlum * dx_stellarlum * fac_sum;
                    }
                }

#endif

#ifdef RT_LEBRON /* now we couple radiation pressure [single-scattering] terms within this module */
                int kf_rt; double lum_force_fac=0;
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX) /* save the fluxes for use below, where we will calculate their RP normally */
                double fac_flux = -fac * All.cf_a2inv / (4.*M_PI); // ~L/(4pi*r^3), in -physical- units (except for last r, cancelled by dx_stellum), since L is physical
                for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++) {Rad_Flux[kf_rt][0]+=mass_stellarlum[kf_rt]*fac_flux*dx_stellarlum; Rad_Flux[kf_rt][1]+=mass_stellarlum[kf_rt]*fac_flux*dy_stellarlum; Rad_Flux[kf_rt][2]+=mass_stellarlum[kf_rt]*fac_flux*dz_stellarlum;}
#else /* simply apply an on-the-spot approximation and do the absorption and RP force now */
                for(kf_rt=0;kf_rt<N_RT_FREQ_BINS;kf_rt++) {lum_force_fac += mass_stellarlum[kf_rt] * fac_stellum[kf_rt];} // add directly to forces. appropriate normalization (and sign) in 'fac_stellum'
#endif
#ifdef BH_PHOTONMOMENTUM /* divide out PhotoMom_coupled_frac here b/c we have our own BH_Rad_Mom factor, and don't want to double-count */
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
                Rad_Flux[RT_FREQ_BIN_FIRE_IR][0]+=mass_bhlum*fac_flux*dx_stellarlum; Rad_Flux[RT_FREQ_BIN_FIRE_IR][1]+=mass_bhlum*fac_flux*dy_stellarlum; Rad_Flux[RT_FREQ_BIN_FIRE_IR][2]+=mass_bhlum*fac_flux*dz_stellarlum;
#elif !defined(RT_DISABLE_RAD_PRESSURE)
                lum_force_fac += (All.BH_Rad_MomentumFactor / (MIN_REAL_NUMBER + All.PhotonMomentum_Coupled_Fraction)) * mass_bhlum * fac_stellum[N_RT_FREQ_BINS-1];
#endif
#endif
                if(lum_force_fac>0) {acc_x += FLT(dx_stellarlum * fac*lum_force_fac); acc_y += FLT(dy_stellarlum * fac*lum_force_fac); acc_z += FLT(dz_stellarlum * fac*lum_force_fac);}
#endif
            } // closes if(valid_gas_particle_for_rt)

#endif // RT_USE_GRAVTREE


#ifdef DM_SCALARFIELD_SCREENING
            if(ptype != 0)	/* we have a dark matter particle as target */
            {
                GRAVITY_NEAREST_XYZ(dx_dm,dy_dm,dz_dm,-1);
                r2 = dx_dm * dx_dm + dy_dm * dy_dm + dz_dm * dz_dm;
                r = sqrt(r2);
                if(r >= h)
                    fac = mass_dm / (r2 * r);
                else
                {
                    h_inv = 1.0 / h;
                    h3_inv = h_inv * h_inv * h_inv;
                    u = r * h_inv;
                    fac = mass_dm * kernel_gravity(u, h_inv, h3_inv, 1);
                }
                /* assemble force with strength, screening length, and target charge.  */
                fac *= All.ScalarBeta * (1 + r / All.ScalarScreeningLength) * exp(-r / All.ScalarScreeningLength);
#ifdef PMGRID
                tabindex = (int) (asmthfac * r);
                if(tabindex < NTAB && tabindex >= 0)
#endif
                {
#ifdef PMGRID
                    fac *= shortrange_table[tabindex];
#endif
                    acc_x += FLT(dx_dm * fac);
                    acc_y += FLT(dy_dm * fac);
                    acc_z += FLT(dz_dm * fac);
                }
            } // closes if(ptype != 0)
#endif // DM_SCALARFIELD_SCREENING //

        } // closes (if((r2 > 0) && (mass > 0))) check

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
                    no = Nodes[no].u.d.nextnode;	/* open it */
                }
            }
        } // closes (mode == 1) check
    } // closes outer (while(no>=0)) check


    /* store result at the proper place */
    if(mode == 0)
    {
        P[target].GravAccel[0] = acc_x;
        P[target].GravAccel[1] = acc_y;
        P[target].GravAccel[2] = acc_z;
#ifdef RT_USE_TREECOL_FOR_NH
        int k;
        for(k=0; k < RT_USE_TREECOL_FOR_NH; k++) P[target].ColumnDensityBins[k] = treecol_angular_bins[k];
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
        P[target].TreeMass = tree_mass;
#endif
#ifdef RT_OTVET
        if(valid_gas_particle_for_rt) {int k,k_et; for(k=0;k<N_RT_FREQ_BINS;k++) for(k_et=0;k_et<6;k_et++) {SphP[target].ET[k][k_et] = RT_ET[k][k_et];}} else {if(P[target].Type==0) {int k,k_et; for(k=0;k<N_RT_FREQ_BINS;k++) for(k_et=0;k_et<6;k_et++) {SphP[target].ET[k][k_et]=0;}}}
#endif
#ifdef CHIMES_STELLAR_FLUXES
        if(valid_gas_particle_for_rt)
        {
            int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {
                SphP[target].Chimes_G0[kc] = chimes_flux_G0[kc]; SphP[target].Chimes_fluxPhotIon[kc] = chimes_flux_ion[kc];}
        }
#endif
#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
        P[target].MencInRcrit = m_enc_in_rcrit;
#endif
#ifdef BH_COMPTON_HEATING
        if(valid_gas_particle_for_rt) {SphP[target].Rad_Flux_AGN = incident_flux_agn;}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        if(valid_gas_particle_for_rt) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {SphP[target].Rad_E_gamma[kf] = Rad_E_gamma[kf];}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        if(valid_gas_particle_for_rt) {int kf,k2; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k2=0;k2<3;k2++) {SphP[target].Rad_Flux[kf][k2] = Rad_Flux[kf][k2];}}}
#endif
#ifdef EVALPOTENTIAL
        P[target].Potential = pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        {int i1,i2; for(i1 = 0; i1 < 3; i1++) {for(i2 = 0; i2 < 3; i2++) {P[target].tidal_tensorps[i1][i2] = tidal_tensorps[i1][i2];}}}
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        {int i1; for(i1 = 0; i1 < 3; i1++) {P[target].GravJerk[i1] = jerk[i1];}}
#endif
#ifdef BH_CALC_DISTANCES
        P[target].min_dist_to_bh = sqrt( min_dist_to_bh2 );
        P[target].min_xyz_to_bh[0] = min_xyz_to_bh[0];   /* remember, dx = x_BH - myx */
        P[target].min_xyz_to_bh[1] = min_xyz_to_bh[1];
        P[target].min_xyz_to_bh[2] = min_xyz_to_bh[2];
#ifdef SINGLE_STAR_FIND_BINARIES
        P[target].is_in_a_binary=0; P[target].min_bh_t_orbital=min_bh_t_orbital; //orbital time for binary
        if (min_bh_t_orbital<MAX_REAL_NUMBER)
        {
	        P[target].is_in_a_binary=1; P[target].comp_Mass=comp_Mass; //mass of binary companion
            int i1; for(i1=0;i1<3;i1++) {P[target].comp_dx[i1]=comp_dx[i1]; P[target].comp_dv[i1]=comp_dv[i1];}
        }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
        P[target].min_bh_approach_time = sqrt(min_bh_approach_time);
        P[target].min_bh_freefall_time = sqrt(sqrt(min_bh_freefall_time)/All.G);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        P[target].min_bh_fb_time = sqrt(min_bh_fb_time);
#endif  
#endif
#endif // BH_CALC_DISTANCES        
    }
    else
    {
        GravDataResult[target].Acc[0] = acc_x;
        GravDataResult[target].Acc[1] = acc_y;
        GravDataResult[target].Acc[2] = acc_z;
#ifdef COUNT_MASS_IN_GRAVTREE
        GravDataResult[target].TreeMass = tree_mass;
#endif
#ifdef RT_USE_TREECOL_FOR_NH
        {int k; for(k=0;k<RT_USE_TREECOL_FOR_NH;k++) GravDataResult[target].ColumnDensityBins[k] = treecol_angular_bins[k];}
#endif
#ifdef RT_OTVET
        {int k,k_et; for(k=0;k<N_RT_FREQ_BINS;k++) for(k_et=0;k_et<6;k_et++) {GravDataResult[target].ET[k][k_et] = RT_ET[k][k_et];}}
#endif
#ifdef CHIMES_STELLAR_FLUXES
        int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {GravDataResult[target].Chimes_G0[kc] = chimes_flux_G0[kc]; GravDataResult[target].Chimes_fluxPhotIon[kc] = chimes_flux_ion[kc];}
#endif
#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
        GravDataResult[target].MencInRcrit = m_enc_in_rcrit;
#endif
#ifdef BH_COMPTON_HEATING
        GravDataResult[target].Rad_Flux_AGN = incident_flux_agn;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {GravDataResult[target].Rad_E_gamma[kf] = Rad_E_gamma[kf];}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        {int kf,k2; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k2=0;k2<3;k2++) {GravDataResult[target].Rad_Flux[kf][k2] = Rad_Flux[kf][k2];}}}
#endif
#ifdef EVALPOTENTIAL
        GravDataResult[target].Potential = pot;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        {int i1,i2; for(i1 = 0; i1 < 3; i1++) {for(i2 = 0; i2 < 3; i2++) {GravDataResult[target].tidal_tensorps[i1][i2] = tidal_tensorps[i1][i2];}}}
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        {int i1; for(i1 = 0; i1 < 3; i1++) {GravDataResult[target].GravJerk[i1] = jerk[i1];}}
#endif
#ifdef BH_CALC_DISTANCES
        GravDataResult[target].min_dist_to_bh = sqrt( min_dist_to_bh2 );
        GravDataResult[target].min_xyz_to_bh[0] = min_xyz_to_bh[0];   /* remember, dx = x_BH - myx */
        GravDataResult[target].min_xyz_to_bh[1] = min_xyz_to_bh[1];
        GravDataResult[target].min_xyz_to_bh[2] = min_xyz_to_bh[2];
#ifdef SINGLE_STAR_FIND_BINARIES
        GravDataResult[target].is_in_a_binary=0; GravDataResult[target].min_bh_t_orbital=min_bh_t_orbital; // orbital time for binary
        if (min_bh_t_orbital<MAX_REAL_NUMBER)
        {
		    GravDataResult[target].is_in_a_binary = 1; GravDataResult[target].comp_Mass=comp_Mass; //mass of binary companion
            int i1; for(i1=0;i1<3;i1++) {GravDataResult[target].comp_dx[i1]=comp_dx[i1]; GravDataResult[target].comp_dv[i1]=comp_dv[i1];}
	    }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
        GravDataResult[target].min_bh_approach_time = sqrt(min_bh_approach_time);
        GravDataResult[target].min_bh_freefall_time = sqrt(sqrt(min_bh_freefall_time)/All.G);
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        GravDataResult[target].min_bh_fb_time = sqrt(min_bh_fb_time);
#endif        
#endif
#endif // BH_CALC_DISTANCES        
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
    double dx, dy, dz, mass, r2, u, v, w, f1, f2, f3, f4, f5, f6, f7, f8;
    double boxsize, boxhalf, pos_x, pos_y, pos_z, aold;
    MyLongDouble acc_x, acc_y, acc_z, xtmp; xtmp=0;

    boxsize = All.BoxSize;
    boxhalf = 0.5 * All.BoxSize;

    acc_x = 0;
    acc_y = 0;
    acc_z = 0;
    cost = 0;
    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        aold = All.ErrTolForceAcc * P[target].OldAcc;
    }
    else
    {
        pos_x = GravDataGet[target].Pos[0];
        pos_y = GravDataGet[target].Pos[1];
        pos_z = GravDataGet[target].Pos[2];
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
    }

    if(mode == 0)
    {
        no = All.MaxPart;		/* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;	/* open it */
    }

    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)	/* single particle */
            {
                /* the index of the node is the index of the particle */
                /* observe the sign */
                if(P[no].Ti_current != All.Ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    drift_particle(no, All.Ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }

                dx = P[no].Pos[0] - pos_x;
                dy = P[no].Pos[1] - pos_y;
                dz = P[no].Pos[2] - pos_z;
                mass = P[no].Mass;
            }
            else			/* we have an  internal node */
            {
                if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
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
                            LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
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
                            UNLOCK_NEXPORT;
                            if(exitFlag) {return -1;} /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */

                            exportnodecount[task] = 0;
                            exportindex[task] = nexp;
                            DataIndexTable[nexp].Task = task;
                            DataIndexTable[nexp].Index = target;
                            DataIndexTable[nexp].IndexGet = nexp;
                        }

                        DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
                        DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

                        if(exportnodecount[task] < NODELISTLENGTH)
                            DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }

                nop = &Nodes[no];

                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }

                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* open cell */
                    no = nop->u.d.nextnode;
                    continue;
                }

                if(nop->Ti_current != All.Ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    force_drift_node(no, All.Ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }

                mass = nop->u.d.mass;
                dx = nop->u.d.s[0] - pos_x;
                dy = nop->u.d.s[1] - pos_y;
                dz = nop->u.d.s[2] - pos_z;
            }
            GRAVITY_NEAREST_XYZ(dx,dy,dz,-1);

            if(no < All.MaxPart)
                no = Nextnode[no];
            else			/* we have an  internal node. Need to check opening criterion */
            {
                openflag = 0;
                r2 = dx * dx + dy * dy + dz * dz;
                if(r2 <= 0) {r2=MIN_REAL_NUMBER;}
                if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                        openflag = 1;
                    }
                }
#ifndef GRAVITY_HYBRID_OPENING_CRIT
                else		/* check relative opening criterion */
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
                        if(GRAVITY_NGB_PERIODIC_BOX_LONG_X(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1) < 0.60 * nop->len)
                        {
                            if(GRAVITY_NGB_PERIODIC_BOX_LONG_Y(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1) < 0.60 * nop->len)
                            {
                                if(GRAVITY_NGB_PERIODIC_BOX_LONG_Z(nop->center[0] - pos_x, nop->center[1] - pos_y, nop->center[2] - pos_z, -1) < 0.60 * nop->len)
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

                    u = nop->center[0] - pos_x;
                    if(u > boxhalf)
                        u -= boxsize;
                    if(u < -boxhalf)
                        u += boxsize;
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

                    u = nop->center[1] - pos_y;
                    if(u > boxhalf)
                        u -= boxsize;
                    if(u < -boxhalf)
                        u += boxsize;
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

                    u = nop->center[2] - pos_z;
                    if(u > boxhalf)
                        u -= boxsize;
                    if(u < -boxhalf)
                        u += boxsize;
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

                    /* if the cell is too large, we need to refine
                     * it further
                     */
                    if(nop->len > 0.20 * boxsize)
                    {
                        /* cell is too large */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }

                no = nop->u.d.sibling;	/* ok, node can be used */
            }

            /* compute the Ewald correction force */

            if(dx < 0)
            {
                dx = -dx;
                signx = +1;
            }
            else
                signx = -1;
            if(dy < 0)
            {
                dy = -dy;
                signy = +1;
            }
            else
                signy = -1;
            if(dz < 0)
            {
                dz = -dz;
                signz = +1;
            }
            else
                signz = -1;
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
            /* compute factors for trilinear interpolation */
            f1 = (1 - u) * (1 - v) * (1 - w);
            f2 = (1 - u) * (1 - v) * (w);
            f3 = (1 - u) * (v) * (1 - w);
            f4 = (1 - u) * (v) * (w);
            f5 = (u) * (1 - v) * (1 - w);
            f6 = (u) * (1 - v) * (w);
            f7 = (u) * (v) * (1 - w);
            f8 = (u) * (v) * (w);
            acc_x += FLT(mass * signx * (fcorrx[i][j][k] * f1 +
                                         fcorrx[i][j][k + 1] * f2 +
                                         fcorrx[i][j + 1][k] * f3 +
                                         fcorrx[i][j + 1][k + 1] * f4 +
                                         fcorrx[i + 1][j][k] * f5 +
                                         fcorrx[i + 1][j][k + 1] * f6 +
                                         fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8));
            acc_y +=
            FLT(mass * signy *
                (fcorry[i][j][k] * f1 + fcorry[i][j][k + 1] * f2 +
                 fcorry[i][j + 1][k] * f3 + fcorry[i][j + 1][k + 1] * f4 + fcorry[i +
                                                                                  1]
                 [j][k] * f5 + fcorry[i + 1][j][k + 1] * f6 + fcorry[i + 1][j +
                                                                            1][k] *
                 f7 + fcorry[i + 1][j + 1][k + 1] * f8));
            acc_z +=
            FLT(mass * signz *
                (fcorrz[i][j][k] * f1 + fcorrz[i][j][k + 1] * f2 +
                 fcorrz[i][j + 1][k] * f3 + fcorrz[i][j + 1][k + 1] * f4 + fcorrz[i +
                                                                                  1]
                 [j][k] * f5 + fcorrz[i + 1][j][k + 1] * f6 + fcorrz[i + 1][j +
                                                                            1][k] *
                 f7 + fcorrz[i + 1][j + 1][k + 1] * f8));
            cost++;
        }

        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0)
                    no = Nodes[no].u.d.nextnode;	/* open it */
            }
        }
    }

    /* add the result at the proper place */

    if(mode == 0)
    {
        P[target].GravAccel[0] += acc_x;
        P[target].GravAccel[1] += acc_y;
        P[target].GravAccel[2] += acc_z;
    }
    else
    {
        GravDataResult[target].Acc[0] = acc_x;
        GravDataResult[target].Acc[1] = acc_y;
        GravDataResult[target].Acc[2] = acc_z;
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
    MyLongDouble pot;
    int no, ptype, task, nexport_save, listindex = 0;
    double r2, dx, dy, dz, mass, r, u, h, h_inv;
    double pos_x, pos_y, pos_z, aold;
    double fac, dxx, dyy, dzz;
#ifdef PMGRID
    int tabindex;
    double eff_dist, rcut, asmth, asmthfac;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    double soft = 0;
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
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
        if((ptype == 0) && (PPP[target].Hsml > All.ForceSoftening[ptype]))
        {
            soft = PPP[target].Hsml;
        } else {
            soft = All.ForceSoftening[ptype];
        }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        soft = PPP[target].AGS_Hsml;
#endif
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
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
        if(ptype == 0) {soft = GravDataGet[target].Soft;}
#endif
#if defined(PMGRID) && defined(PM_PLACEHIGHRESREGION)
        if(pmforce_is_particle_high_res(ptype, GravDataGet[target].Pos))
        {
            rcut = All.Rcut[1];
            asmth = All.Asmth[1];
        }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        soft = GravDataGet[target].Soft;
#endif
    }

#ifdef PMGRID
    asmthfac = 0.5 / asmth * (NTAB / 3.0);
#endif
    if(mode == 0)
    {
        no = All.MaxPart;		/* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;	/* open it */
    }

    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)	/* single particle */
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
                if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
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

                        DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                        DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                        if(Exportnodecount[task] < NODELISTLENGTH)
                            DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }

                nop = &Nodes[no];
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }

                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* open cell */
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
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                h = soft; /* set softening */
#else
                if(ptype == 0) {h = soft;} else {h = All.ForceSoftening[ptype];} /* set softening */
#endif
#else
                h = All.ForceSoftening[ptype];
                if(h < All.ForceSoftening[P[no].Type]) {h = All.ForceSoftening[P[no].Type];}
#endif
                no = Nextnode[no];
            }
            else			/* we have an internal node. Need to check opening criterion */
            {
#ifdef PMGRID
                /* check whether we can stop walking along this branch */
                if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
                {
                    if(mode == 0)
                    {
                        if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                        {
                            Exportflag[task] = target;
                            DataIndexTable[*nexport].Index = target;
                            DataIndexTable[*nexport].Task = task;	/* Destination task */
                            *nexport = *nexport + 1;
                            nsend_local[task]++;
                        }
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }

                eff_dist = rcut + 0.5 * nop->len;
                dxx = nop->center[0] - pos_x;	/* observe the sign ! */
                dyy = nop->center[1] - pos_y;	/* this vector is -y in my thesis notation */
                dzz = nop->center[2] - pos_z;
                GRAVITY_NEAREST_XYZ(dxx,dyy,dzz,-1);
#ifdef REDUCE_TREEWALK_BRANCHING
                if((fabs(dxx) > eff_dist) | (fabs(dyy) > eff_dist) | (fabs(dzz) > eff_dist))
                {
                    no = nop->u.d.sibling;
                    continue;
                }
#else
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
#endif // REDUCE_TREEWALK_BRANCHING
#else // PMGRID
                dxx = nop->center[0] - pos_x;	/* observe the sign ! */
                dyy = nop->center[1] - pos_y;	/* this vector is -y in my thesis notation */
                dzz = nop->center[2] - pos_z;
                GRAVITY_NEAREST_XYZ(dxx,dyy,dzz,-1);
#endif // PMGRID

                if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#ifndef GRAVITY_HYBRID_OPENING_CRIT
                else		/* check relative opening criterion */
#else
                if(!(All.Ti_Current == 0 && RestartFlag != 1))
#endif
                {

                    /* force node to open if we are within the gravitational softening length */
#if defined(SELFGRAVITY_OFF) || defined(RT_SELFGRAVITY_OFF) || (!(defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)))
                    double soft = All.ForceSoftening[ptype];
#endif
                    if((r2 < (soft+0.6*nop->len)*(soft+0.6*nop->len)) || (r2 < (nop->maxsoft+0.6*nop->len)*(nop->maxsoft+0.6*nop->len)))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

#ifdef REDUCE_TREEWALK_BRANCHING
                    if((mass * nop->len * nop->len > r2 * r2 * aold) |
                       ((fabs(dxx) < 0.60 * nop->len) & (fabs(dyy) < 0.60 * nop->len) & (fabs(dzz) < 0.60 * nop->len)))
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
#else
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
#endif // REDUCE_TREEWALK_BRANCHING //
                }

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                h = soft;
#else
                if(ptype == 0) {h = soft;} else {h = All.ForceSoftening[ptype];}
#endif

                if(h < nop->maxsoft)
                {
                    //h = nop->maxsoft; // only applies if symmetrizing with MAX(h_i,h_j)
                    if(r2 < nop->maxsoft * nop->maxsoft)
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#else
                h = All.ForceSoftening[ptype];
                if(h < nop->maxsoft)
                {
                    h = nop->maxsoft;
                    if(r2 < h * h)
                    {
                        /* bit-5 signals that there are particles of
                         * different softening in the node
                         */
                        if(maskout_different_softening_flag(nop->u.d.bitflags))
                        {
                            no = nop->u.d.nextnode;
                            continue;
                        }
                    }
                }
#endif // #if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL) //
                no = nop->u.d.sibling;	/* node can be used */
            }

            r = sqrt(r2);
#ifdef PMGRID
            tabindex = (int) (r * asmthfac);
            if(tabindex < NTAB && tabindex >= 0)
#endif
            {
#ifdef PMGRID
                fac = shortrange_table_potential[tabindex];
#else
                fac = 1;
#endif
                if(r >= h)
                {
                    pot += FLT(-fac * mass / r);
                } else {
                    h_inv = 1.0 / h;
                    u = r * h_inv;
                    pot += FLT( fac * mass * kernel_gravity(u, h_inv, 1, -1) );
                }
            }
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
            pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
        }
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0)
                    no = Nodes[no].u.d.nextnode;	/* open it */
            }
        }
    }

    /* store result at the proper place */
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL)
    if(mode == 0)
        P[target].Potential = pot;
    else
        PotDataResult[target].Potential = pot;
#endif
    return 0;
}





#ifdef SUBFIND
int subfind_force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
{
    struct NODE *nop = 0;
    MyLongDouble pot;
    int no, ptype, task, nexport_save, listindex = 0;
    double r2, dx, dy, dz, mass, r, u, h, h_inv;
    double pos_x, pos_y, pos_z;

    nexport_save = *nexport;
    pot = 0;
    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        ptype = P[target].Type;
    }
    else
    {
        pos_x = GravDataGet[target].Pos[0];
        pos_y = GravDataGet[target].Pos[1];
        pos_z = GravDataGet[target].Pos[2];
        ptype = GravDataGet[target].Type;
    }

    h = All.ForceSoftening[ptype];
    h_inv = 1.0 / h;

    if(mode == 0)
    {
        no = All.MaxPart;		/* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;	/* open it */
    }

    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)	/* single particle */
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
                if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
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

                        DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                        DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                        if(Exportnodecount[task] < NODELISTLENGTH)
                            DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }

                nop = &Nodes[no];
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }

                mass = nop->u.d.mass;
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* open cell */
                    if(mass)
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
            else			/* we have an internal node. Need to check opening criterion */
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
                no = nop->u.d.sibling;	/* node can be used */
            }

            r = sqrt(r2);
            if(r >= h)
                pot += FLT(-mass / r);
            else
            {
                u = r * h_inv;
                pot += FLT( mass * kernel_gravity(u, h_inv, 1, -1) );
            }
        }
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0)
                    no = Nodes[no].u.d.nextnode;	/* open it */
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
        printf("failed to allocate memory for %d tree-extnodes (%g MB).\n",
               MaxNodes, bytes / (1024.0 * 1024.0));
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
    char buffer[200];
    int i;

    sprintf(buffer, "particles%d.dat", ThisTask);
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
    char buf[200];
    FILE *fd;

    if(ThisTask == 0) {printf("Initializing Ewald correction...\n");}

#ifdef DOUBLEPRECISION
    sprintf(buf, "ewald_spc_table_%d_dbl.dat", EN);
#else
    sprintf(buf, "ewald_spc_table_%d.dat", EN);
#endif
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
