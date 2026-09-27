#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"

/*! \file star_direct_gravity.cc
 *  \brief exact direct-summation gravity between star (type-5) particles.
 *
 *  The tree approximates the field of a group of distant stars by its monopole, which is fine for
 *  gas but destroys the few-body dynamics stars have with each other: the error is uncorrelated
 *  between steps, so it does not average out, and it is largest exactly where the orbits are
 *  tightest. With SINGLE_STAR_DIRECT_GRAVITY every star-star pair is instead evaluated exactly.
 *
 *  This is only tractable because the stars are few: the table is O(N_star) on EVERY task, and the
 *  sum is O(N_star^2). Do not enable it for a run whose star count is comparable to its cell count.
 *
 *  Correctness here depends on the tree NOT also supplying star->star forces, or every pair would
 *  be counted twice. force_treeevaluate() drops those contributions when the target is type 5; see
 *  the SINGLE_STAR_DIRECT_GRAVITY blocks there.
 */

#ifdef SINGLE_STAR_DIRECT_GRAVITY

/* same guard forcetree.cc uses: wrap separations only when gravity itself is periodic */
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
#define GRAVITY_NEAREST_XYZ(x,y,z,sign) NEAREST_XYZ(x,y,z,sign)
#else
#define GRAVITY_NEAREST_XYZ(x,y,z,sign) /* no box-wrapping needed */
#endif

struct star_direct_data *StarDirect;
int N_StarDirect;

/*! Gather every star in the simulation onto every task.
 *
 *  Positions are drifted to All.Ti_Current first: the tree walk drifts the particles it touches
 *  on demand (drift_particle, under a critical section), but a brute-force pass has no tree to
 *  hang that off, and mixing drifted and undrifted positions in a pair sum would break Newton's
 *  third law -- a self-accelerating cluster, not a small error. */
void star_direct_gravity_build_table(void)
{
    int i, n_local = 0;
    for(i = 0; i < NumPart; i++) {if(P[i].Type == 5 && P[i].Mass > 0) {n_local++;}}

    /* mymalloc is a LIFO stack allocator, so StarDirect -- the one buffer that has to outlive this
       function -- must sit at the bottom, underneath every temporary. That means we need the total
       star count BEFORE any array is allocated, hence an Allreduce here and the per-task counts
       gathered further down rather than one combined Allgather. */
    MPI_Allreduce(&n_local, &N_StarDirect, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* mymalloc refuses a zero-byte request, and a star-free run (before the first sink forms) is
       the normal early state of every STARFORGE simulation, not an edge case */
    size_t nalloc = (N_StarDirect > 0) ? (size_t)N_StarDirect : 1;
    StarDirect = (struct star_direct_data *) mymalloc("StarDirect", nalloc * sizeof(struct star_direct_data));
    if(N_StarDirect == 0) {return;}

    int *counts = (int *) mymalloc("StarDirectCounts", NTask * sizeof(int));
    int *offsets = (int *) mymalloc("StarDirectOffsets", NTask * sizeof(int));
    MPI_Allgather(&n_local, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
    int running = 0;
    for(i = 0; i < NTask; i++) {offsets[i] = running; running += counts[i];}

    struct star_direct_data *sendbuf = (struct star_direct_data *) mymalloc("StarDirectSend", ((n_local > 0) ? n_local : 1) * sizeof(struct star_direct_data));
    int k = 0;
    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Type != 5 || P[i].Mass <= 0) {continue;}
        if(P[i].Ti_current != All.Ti_Current) {drift_particle(i, All.Ti_Current);}
        sendbuf[k].Pos = P[i].Pos;
        sendbuf[k].Vel = P[i].Vel;
#ifdef HERMITE_INTEGRATION
        /* Same correction as the tree walk's single-particle branch (forcetree.cc): during the
         * Hermite-only passes, send an inactive star's Old*-predicted state instead of its
         * KDK-drifted one, which is O(dt^2) off mid-step with a whole-step-kicked Vel. Only the
         * send buffer changes, so the KDK path is unaffected. Active stars keep their live state
         * (their Old* are stale at HermiteOnlyFlag==1), and so does any star past its own step
         * span, where the Old* extrapolation would be biased along the orbit. */
        if(HermiteOnlyFlag && !TimeBinActive[P[i].TimeBin]
           && (All.Ti_Current >= P[i].Ti_begstep)
           && ((All.Ti_Current - P[i].Ti_begstep) <= P[i].integertime_step())
           && eligible_for_hermite(i))
        {
            double hD = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current, i, 0);
            sendbuf[k].Pos = P[i].OldPos + (P[i].OldVel + (P[i].Hermite_OldAcc + P[i].OldJerk * (hD/3)) * (hD/2)) * hD;
            sendbuf[k].Vel = P[i].OldVel + (P[i].Hermite_OldAcc + P[i].OldJerk * (hD/2)) * hD;
        }
#endif
        sendbuf[k].Mass = P[i].Mass;
        sendbuf[k].Soft = ForceSoftening_KernelRadius(i);
        sendbuf[k].ID = P[i].ID;
#ifdef KETJU_REGULARIZATION
        sendbuf[k].KetjuChainID = P[i].KetjuChainID;
#endif
        k++;
    }

    /* byte counts, so the struct travels as an opaque blob and this does not need an MPI datatype
       rebuilt every time a member is added to it */
    int *bcounts = (int *) mymalloc("StarDirectBCounts", NTask * sizeof(int));
    int *boffsets = (int *) mymalloc("StarDirectBOffsets", NTask * sizeof(int));
    for(i = 0; i < NTask; i++) {bcounts[i] = counts[i] * sizeof(struct star_direct_data); boffsets[i] = offsets[i] * sizeof(struct star_direct_data);}
    MPI_Allgatherv(sendbuf, n_local * sizeof(struct star_direct_data), MPI_BYTE, StarDirect, bcounts, boffsets, MPI_BYTE, MPI_COMM_WORLD);

    myfree(boffsets); myfree(bcounts); myfree(sendbuf); myfree(offsets); myfree(counts);
}


void star_direct_gravity_free_table(void)
{
    myfree(StarDirect); StarDirect = NULL; N_StarDirect = 0;
}


/*! Add the exact star-star terms to every active local star.
 *
 *  Must run after the tree walk has finished accumulating (including the imported contributions
 *  from other tasks) and before gravtree.cc multiplies by All.G, so the sums here are in the same
 *  G-free units the tree walk uses.
 *
 *  Pairs are identified by ID, not by index: the table is a global concatenation, so a target's own
 *  entry sits at an offset that depends on which task owns it. */
void star_direct_gravity_compute(void)
{
    if(N_StarDirect <= 0) {return;}
    int n_active = (int) ActiveParticleList.size();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int ii = 0; ii < n_active; ii++)
    {
        int i = ActiveParticleList[ii];
        if(P[i].Type != 5 || P[i].Mass <= 0) {continue;}
#ifdef HERMITE_INTEGRATION
        if(HermiteOnlyFlag) {if(!eligible_for_hermite(i)) {continue;}}
#endif
        Vec3<double> pos = P[i].Pos, acc = {};
        double h_i = ForceSoftening_KernelRadius(i);
        MyIDType id_i = P[i].ID;
#ifdef KETJU_REGULARIZATION
        MyIDType chain_i = P[i].KetjuChainID;
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_CALC_DISTANCES)
        Vec3<double> vel = P[i].Vel;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        Vec3<double> jerk = {};
#endif
#ifdef EVALPOTENTIAL
        double pot = 0;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        SymmetricTensor2<double> tidal = {};
#endif
#ifdef SINK_CALC_DISTANCES
        double Min_Distance_to_Sink2 = MAX_REAL_NUMBER; Vec3<double> Min_xyz_to_Sink = {MAX_REAL_NUMBER,MAX_REAL_NUMBER,MAX_REAL_NUMBER};
        double Min_Sink_Approach_Time = MAX_REAL_NUMBER, Min_Sink_Freefall_time = MAX_REAL_NUMBER;
#endif

        for(int j = 0; j < N_StarDirect; j++)
        {
            if(StarDirect[j].ID == id_i) {continue;} /* self */
            double mass = StarDirect[j].Mass; if(mass <= 0) {continue;}
            Vec3<double> dr = StarDirect[j].Pos - pos;
            GRAVITY_NEAREST_XYZ(dr[0],dr[1],dr[2],-1);
            double r2 = dr.norm_sq(); if(r2 <= 0) {continue;}
            double r = sqrt(r2);

            /* symmetrize by taking the larger softening, matching the non-averaged branch the tree
               uses for sink pairs (ADAPTIVE_GRAVSOFT_SYMMETRIZE_FORCE_BY_AVERAGING is explicitly
               disabled for sink interactions in force_treeevaluate) */
            double h = DMAX(h_i, StarDirect[j].Soft);
            double fac_accel, fac2_tidal;
#ifdef EVALPOTENTIAL
            double fac_pot;
#endif
            if(r >= h)
            {
                fac_accel = mass / (r2 * r);
                fac2_tidal = 3.0 * mass / (r2 * r2 * r);
#ifdef EVALPOTENTIAL
                fac_pot = -mass / r;
#endif
            }
            else
            {
                double h_inv = 1./h, h3_inv = h_inv*h_inv*h_inv, u = r*h_inv;
                fac_accel = mass * kernel_gravity(u, h_inv, h3_inv, 1);
                fac2_tidal = mass * kernel_gravity(u, h_inv, h3_inv, 2);
#ifdef EVALPOTENTIAL
                fac_pot = mass * kernel_gravity(u, h_inv, h3_inv, -1);
#endif
            }
            acc += fac_accel * dr;
#ifdef EVALPOTENTIAL
            pot += fac_pot;
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_CALC_DISTANCES)
            Vec3<double> dv = StarDirect[j].Vel - vel;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
            jerk += fac_accel * dv - dot(dv, dr) * fac2_tidal * dr;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
            double fac_tidal = fac_accel;
            tidal[0][0] += (-fac_tidal + dr[0]*dr[0]*fac2_tidal);
            tidal[0][1] += (dr[0]*dr[1]*fac2_tidal);
            tidal[0][2] += (dr[0]*dr[2]*fac2_tidal);
            tidal[1][1] += (-fac_tidal + dr[1]*dr[1]*fac2_tidal);
            tidal[1][2] += (dr[1]*dr[2]*fac2_tidal);
            tidal[2][2] += (-fac_tidal + dr[2]*dr[2]*fac2_tidal);
#endif
#ifdef SINK_CALC_DISTANCES
            /* the tree normally supplies these for star targets from the sink nodes it visits; those
               are skipped now, so they have to be reproduced here or the sink timestep criteria
               silently see no other stars at all */
            if(r2 < Min_Distance_to_Sink2) {Min_Distance_to_Sink2 = r2; Min_xyz_to_Sink = dr;}
            double r2soft = DMAX(SinkParticle_GravityKernelRadius, h_i) * KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER;
            r2soft = r2 + r2soft*r2soft;
            double vSqr = dv.norm_sq(), M_total = mass + P[i].Mass;
            double tSqr = r2soft/(vSqr + MIN_REAL_NUMBER), tff4 = r2soft*r2soft*r2soft/(M_total*M_total);
#ifdef KETJU_REGULARIZATION
            /* the tree walk's rule (forcetree.cc): MSTAR integrates a same-chain pair closer than the
               region radius, so it does not set the host step */
            if((chain_i == 0) || (StarDirect[j].KetjuChainID != chain_i) || (r2 >= All.KetjuRegionRadius * All.KetjuRegionRadius))
#endif
            {
            if(tSqr < Min_Sink_Approach_Time) {Min_Sink_Approach_Time = tSqr;}
            if(tff4 < Min_Sink_Freefall_time) {Min_Sink_Freefall_time = tff4;}
            }
#endif
        }

        P[i].GravAccel += acc;
#ifdef EVALPOTENTIAL
        P[i].Potential += pot;
#endif
#ifdef COMPUTE_JERK_IN_GRAVTREE
        P[i].GravJerk += jerk;
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
        P[i].tidal_tensorps[0][0] += tidal[0][0]; P[i].tidal_tensorps[0][1] += tidal[0][1]; P[i].tidal_tensorps[0][2] += tidal[0][2];
        P[i].tidal_tensorps[1][1] += tidal[1][1]; P[i].tidal_tensorps[1][2] += tidal[1][2]; P[i].tidal_tensorps[2][2] += tidal[2][2];
#endif
#ifdef SINK_CALC_DISTANCES
        /* min, not sum: take whichever is closer/shorter, the tree's answer (gas, and any star the
           tree still reported) or ours */
        if(sqrt(Min_Distance_to_Sink2) < P[i].Min_Distance_to_Sink) {P[i].Min_Distance_to_Sink = sqrt(Min_Distance_to_Sink2); P[i].Min_xyz_to_Sink = Min_xyz_to_Sink;}
#ifdef SINGLE_STAR_TIMESTEPPING
        if(sqrt(Min_Sink_Approach_Time) < P[i].Min_Sink_Approach_Time) {P[i].Min_Sink_Approach_Time = sqrt(Min_Sink_Approach_Time);}
        {double tff = sqrt(sqrt(Min_Sink_Freefall_time)/All.G); if(tff < P[i].Min_Sink_Freefall_time) {P[i].Min_Sink_Freefall_time = tff;}}
#endif
#endif
    }
}

#endif // SINGLE_STAR_DIRECT_GRAVITY
