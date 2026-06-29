#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <gsl/gsl_eigen.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "./analytic_gravity.h"

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
void sum_top_level_node_costfactors(void);


/*! This function computes the gravitational forces for all active elements. If needed, a new tree is constructed, otherwise the dynamically updated
 *  tree is used.  Elements are only exported to other processors when needed. */
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
        MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_DRIFT] += measure_time(); /* sync before we do the treebuild */
        force_treebuild(NumPart, NULL);
        MPI_Barrier(MPI_COMM_WORLD); CPU_Step[CPU_TREEBUILD] += measure_time(); /* and sync after treebuild as well */
        TreeReconstructFlag = 0;
        TreeMomentsStaleFlag = 0;
        PRINT_STATUS(" ..Tree construction done.");
    }

    force_validate_tree_links("gravtree"); /* tier-0: no-op unless a maintained rearrange changed the threading since the last build */

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
            CPU_Step[CPU_TREEBUILD] += measure_time();
            TreeMomentsStaleFlag = 0;
        }
    }

    CPU_Step[CPU_TREEMISC] += measure_time(); t0 = my_second();
#if !defined(SELFGRAVITY_OFF) || defined(RT_USE_GRAVTREE) || defined(RT_USE_TREECOL_FOR_NH) || defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    /* allocate buffers to arrange communication */
    PRINT_STATUS(" ..Begin tree force. (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (long) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct gravdata_in) + sizeof(struct gravdata_out) +
                                             sizemax(sizeof(struct gravdata_in),sizeof(struct gravdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) {if(ThisTask == 0) printf(" ..All.BunchSize=%ld\n", All.BunchSize);}
    int k, ewald_max, diff, save_NextParticle, ndone, ndone_flag, place, recvTask; double tstart, tend, ax, ay, az; MPI_Status status;
    Ewaldcount = 0; Costtotal = 0; N_nodesinlist = 0; ewald_max=0;
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
    ewald_max = 1; /* the tree-code will need to iterate to perform the periodic boundary condition corrections */
#endif

    if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)
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
                {terminate("TakeLevel < 0, even though we should have a slot");}
        }
    }
    else
    { /* in this case we do not measure gravity cost. Check whether this time-level
         has previously mean measured. If yes, then delete it so to make sure that it is not out of time */
        for(i = 0; i < GRAVCOSTLEVELS; i++) {if(All.LevelToTimeBin[i] == All.HighestActiveTimeBin) {All.LevelToTimeBin[i] = 0;}}
        TakeLevel = -1;
    }
    if(TakeLevel >= 0) {for(i = 0; i < NumPart; i++) {P[i].GravCost[TakeLevel] = 0;}} /* re-zero the cost [will be re-summed] */

    /* cache which particles need a new tree force BEFORE the tree walk runs: the tree walk can modify
       quantities like Min_Sink_FeedbackTime that needs_new_treeforce() depends on, so re-calling it
       in the post-processing loop could give a different answer, causing a particle that was computed
       by the tree walk (raw GravAccel, raw GravJerk) to incorrectly take the jerk-skip path (which
       expects G-multiplied GravAccel/GravJerk from a previous step). */
#ifdef ADAPTIVE_TREEFORCE_UPDATE
    std::vector<int> treeforce_skip_flag(ActiveParticleList.size(), 0);
    for(int ii = 0; ii < (int)ActiveParticleList.size(); ii++) {
        int i = ActiveParticleList[ii];
        if(!needs_new_treeforce(i)) {treeforce_skip_flag[ii] = 1;}
    }
#endif

    /* begin main communication and tree-walk loop. note the ewald-iter terms here allow for multiple iterations for periodic-tree corrections if needed */
    for(Ewald_iter = 0; Ewald_iter <= ewald_max; Ewald_iter++)
    {
        NextParticle = 0;	/* begin with this index */
        memset(ProcessedFlag, 0, All.MaxPart * sizeof(unsigned char));
        BufferCollisionFlag = 0; /* set to zero before operations begin */
        do /* primary point-element loop */
        {
            iter++;
            BufferFullFlag = 0; Nexport = 0; save_NextParticle = NextParticle; tstart = my_second();

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

            if(BufferFullFlag) /* we've filled the buffer or reached the end of the list, prepare for communications */
            {
                int last_nextparticle = NextParticle;
                int processed_particles = 0;
                int first_unprocessedparticle = -1;
                NextParticle = save_NextParticle; /* figure out where we are */
                while(NextParticle < (int)ActiveParticleList.size())
                {
                    if(NextParticle == last_nextparticle) {break;}
                    int pindex = ActiveParticleList[NextParticle];
#ifndef _OPENMP
                    if(ProcessedFlag[pindex] != 1) {break;}
#else
                    if(ProcessedFlag[pindex] == 0 && first_unprocessedparticle < 0) {first_unprocessedparticle = NextParticle;}
                    if(ProcessedFlag[pindex] == 1)
#endif
                    {
                        processed_particles++;
                        ProcessedFlag[pindex] = 2;
                    }
                    NextParticle++;
                }
#ifdef _OPENMP
                if(first_unprocessedparticle >= 0) {NextParticle = first_unprocessedparticle;} /* reset the neighbor list properly for the next group since we can get 'jumps' with openmp active */
                if(processed_particles == 0 && NextParticle == save_NextParticle && NextParticle < (int)ActiveParticleList.size()) {
                    BufferCollisionFlag++; if(BufferCollisionFlag < 2) {continue;}} /* we overflowed without processing a single particle, but this could be because of a collision, try once with the serialized approach, but if it fails then, we're truly stuck */
                else if(processed_particles && BufferCollisionFlag) {BufferCollisionFlag = 0;} /* we had a problem in a previous iteration but things worked, reset to normal operations */
#endif
                if(processed_particles <= 0 && NextParticle == save_NextParticle) {endrun(114408);} /* in this case, the buffer is too small to process even a single particle */

                int new_export = 0; /* actually calculate exports [so we can tell other tasks] */
                for(j = 0, k = 0; j < Nexport; j++)
                {
                    if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                    {
                        if(k < j + 1) {k = j + 1;}
                        for(; k < Nexport; k++)
                            if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                            {
                                int old_index = DataIndexTable[j].Index;
                                DataIndexTable[j] = DataIndexTable[k]; DataNodeList[j] = DataNodeList[k]; DataIndexTable[j].IndexGet = j; new_export++;
                                DataIndexTable[k].Index = old_index; k++;
                                break;
                            }
                    }
                    else {new_export++;}
                }
                Nexport = new_export; /* counting exports... */
            }
            n_exported += Nexport;
            for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
            for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
            mysort_dataindex(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare); /* construct export count tables */
            tstart = my_second();
            MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD); /* broadcast import/export counts */
            tend = my_second(); timewait1 += timediff(tstart, tend);

            for(j = 0, Send_offset[0] = 0; j < NTask; j++) {if(j > 0) {Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];}} /* calculate export table offsets */
            GravDataIn = (struct gravdata_in *) mymalloc("GravDataIn", Nexport * sizeof(struct gravdata_in));
            GravDataOut = (struct gravdata_out *) mymalloc("GravDataOut", Nexport * sizeof(struct gravdata_out));
            for(j = 0; j < Nexport; j++) /* prepare particle data for export [fill in the structures to be passed] */
            {
                place = DataIndexTable[j].Index;

                /* assign values (input-function to pass in memory) */
                GravDataIn[j].Pos = P[place].Pos;
                GravDataIn[j].Type = P[place].Type;
                GravDataIn[j].Soft = ForceSoftening_KernelRadius(place);
                GravDataIn[j].OldAcc = P[place].OldAcc;
                GravDataIn[j].Mass = P[place].Mass;
#if defined(SINK_DYNFRICTION_FROMTREE)
                if(P[place].Type==5) {GravDataIn[j].Sink_Mass = P[place].Sink_Mass;}
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
                GravDataIn[j].Vel = P[place].Vel;
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
                if(P[place].Type == 5)
                {
                    GravDataIn[j].Min_Sink_OrbitalTime = P[place].Min_Sink_OrbitalTime; //orbital time for binary
                    GravDataIn[j].comp_Mass = P[place].comp_Mass; //mass of binary companion
                    GravDataIn[j].is_in_a_binary = P[place].is_in_a_binary; // 1 if we're in a binary, 0 if not
                    GravDataIn[j].comp_dx = P[place].comp_dx; GravDataIn[j].comp_dv = P[place].comp_dv;
                }
                else {P[place].is_in_a_binary=0; /* setting values to zero just to be sure */}
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                if((P[place].Type == 0) && (P[place].KernelRadius > All.ForceSoftening[P[place].Type])) {GravDataIn[j].AGS_zeta = P[place].AGS_zeta;} else {GravDataIn[j].AGS_zeta = 0;}
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                GravDataIn[j].Soft = P[place].AGS_KernelRadius;
                GravDataIn[j].AGS_zeta = P[place].AGS_zeta;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                GravDataIn[j].tidal_tensorps_prevstep=P[place].tidal_tensorps_prevstep;
#endif
#ifdef KETJU_REGULARIZATION
                GravDataIn[j].KetjuChainID = P[place].KetjuChainID;
#endif
                memcpy(GravDataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
            }

            /* ok now we have to figure out if there is enough memory to handle all the tasks sending us their data, and if not, break it into sub-chunks */
            int N_chunks_for_import, ngrp_initial, ngrp;
            for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) /* sub-chunking loop opener */
            {
                int flagall;
                N_chunks_for_import = (1 << PTask) - ngrp_initial;
                do {
                    int flag = 0; Nimport = 0;
                    for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++)
                    {
                        recvTask = ThisTask ^ ngrp;
                        if(recvTask < NTask) {if(Recv_count[recvTask] > 0) {Nimport += Recv_count[recvTask];}}
                    }
                    size_t space_needed = Nimport * sizeof(struct gravdata_in) + Nimport * sizeof(struct gravdata_out) + 16384; /* extra bitflag is a padding, to avoid overflows */
                    if(space_needed > FreeBytes) {flag = 1;}
                    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                    if(flagall) {N_chunks_for_import /= 2;} else {break;}
                } while(N_chunks_for_import > 0);
                if(N_chunks_for_import == 0) {printf("Memory is insufficient for even one import-chunk: N_chunks_for_import=%d  ngrp_initial=%d  Nimport=%ld  FreeBytes=%lld , but we need to allocate=%lld \n",N_chunks_for_import, ngrp_initial, Nimport, (long long)FreeBytes,(long long)(Nimport * sizeof(struct gravdata_in) + Nimport * sizeof(struct gravdata_out) + 16384)); endrun(9966);}
                if(flagall) {if(ThisTask==0) PRINT_WARNING("Splitting import operation into sub-chunks as we are hitting memory limits (check this isn't imposing large communication cost)");}

                /* now allocated the import and results buffers */
                GravDataGet = (struct gravdata_in *) mymalloc("GravDataGet", Nimport * sizeof(struct gravdata_in));
                GravDataResult = (struct gravdata_out *) mymalloc("GravDataResult", Nimport * sizeof(struct gravdata_out));

                tstart = my_second(); Nimport = 0; /* reset because this will be cycled below to calculate the recieve offsets (Recv_offset) */
                for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) /* exchange particle data */
                {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask)
                    {
                        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) /* get the particles */
                        {
                            MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE, recvTask, TAG_GRAV_A,
                                         &GravDataGet[Nimport], Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE, recvTask, TAG_GRAV_A,
                                         MPI_COMM_WORLD, &status);
                            Nimport += Recv_count[recvTask];
                        }
                    }
                }
                tend = my_second(); timecommsumm1 += timediff(tstart, tend);
                report_memory_usage(&HighMark_gravtree, "GRAVTREE");

                /* now do the particles that were sent to us */
                tstart = my_second(); NextJ = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
                {
#ifdef _OPENMP
                    int mainthreadid = omp_get_thread_num();
#else
                    int mainthreadid = 0;
#endif
                    gravity_secondary_loop(&mainthreadid);
                }
                tend = my_second(); timetree2 += timediff(tstart, tend); tstart = my_second();
                MPI_Barrier(MPI_COMM_WORLD); /* insert MPI Barrier here - will be forced by comms below anyways but this allows for clean timing measurements */
                tend = my_second(); timewait2 += timediff(tstart, tend);

                tstart = my_second(); Nimport = 0;
                for(ngrp = ngrp_initial; ngrp < ngrp_initial + N_chunks_for_import; ngrp++) /* send the results for imported elements back to their host tasks */
                {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask)
                    {
                        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                        {
                            MPI_Sendrecv(&GravDataResult[Nimport], Recv_count[recvTask] * sizeof(struct gravdata_out), MPI_BYTE, recvTask, TAG_GRAV_B,
                                         &GravDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct gravdata_out), MPI_BYTE, recvTask, TAG_GRAV_B,
                                         MPI_COMM_WORLD, &status);
                            Nimport += Recv_count[recvTask];
                        }
                    }
                }
                tend = my_second(); timecommsumm2 += timediff(tstart, tend);
                myfree(GravDataResult); myfree(GravDataGet); /* free the structures used to send data back to tasks, its sent */

            } /* close the sub-chunking loop: for(ngrp_initial = 1; ngrp_initial < (1 << PTask); ngrp_initial += N_chunks_for_import) */

            /* we have all our results back from the elements we exported: add the result to the local elements */
            tstart = my_second();
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index;
                P[place].GravAccel += GravDataOut[j].Acc;
                if(Ewald_iter > 0) continue; /* everything below is ONLY evaluated if we are in the first sub-loop, not the periodic correction, or else we will get un-allocated memory or un-physical values */

#ifdef EVALPOTENTIAL
                P[place].Potential += GravDataOut[j].Potential;
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
                P[place].TreeMass += GravDataOut[j].TreeMass;
#endif
#ifdef SINK_CALC_DISTANCES /* GravDataOut[j].Min_Distance_to_Sink contains the min dist to particle "P[place]" on another task.  We now check if it is smaller than the current value */
                if(GravDataOut[j].Min_Distance_to_Sink < P[place].Min_Distance_to_Sink)
                {
                    P[place].Min_Distance_to_Sink = GravDataOut[j].Min_Distance_to_Sink;
                    P[place].Min_xyz_to_Sink = GravDataOut[j].Min_xyz_to_Sink;
#ifdef SPECIAL_POINT_MOTION
                    if(P[place].Type != SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
                    {
                        P[place].vel_of_nearest_special = GravDataOut[j].vel_of_nearest_special;
                        P[place].acc_of_nearest_special = GravDataOut[j].acc_of_nearest_special;
                    }
#endif
                }
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
                if(P[place].Type == SPECIAL_POINT_TYPE_FOR_NODE_DISTANCES)
                {
                    P[place].vel_of_nearest_special += GravDataOut[j].vel_of_nearest_special; /* this is the weighted sum of the velocity around that cell */
                    P[place].acc_of_nearest_special += GravDataOut[j].acc_of_nearest_special; /* this is the weighted sum of the velocity around that cell */
                    P[place].weight_sum_for_special_point_smoothing += GravDataOut[j].weight_sum_for_special_point_smoothing; /* weighted sum needed */
                }
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
                if(GravDataOut[j].Min_Sink_Approach_Time < P[place].Min_Sink_Approach_Time) {P[place].Min_Sink_Approach_Time = GravDataOut[j].Min_Sink_Approach_Time;}
                if(GravDataOut[j].Min_Sink_Freefall_time < P[place].Min_Sink_Freefall_time) {P[place].Min_Sink_Freefall_time = GravDataOut[j].Min_Sink_Freefall_time;}
#ifdef SINGLE_STAR_FIND_BINARIES
                if((P[place].Type == 5) && (GravDataOut[j].Min_Sink_OrbitalTime < P[place].Min_Sink_OrbitalTime))
                {
                    P[place].Min_Sink_OrbitalTime = GravDataOut[j].Min_Sink_OrbitalTime;
                    P[place].comp_Mass = GravDataOut[j].comp_Mass;
                    P[place].is_in_a_binary = GravDataOut[j].is_in_a_binary;
                    P[place].comp_dx=GravDataOut[j].comp_dx; P[place].comp_dv=GravDataOut[j].comp_dv;
                }
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                if(GravDataOut[j].Min_Sink_FeedbackTime < P[place].Min_Sink_FeedbackTime) {P[place].Min_Sink_FeedbackTime = GravDataOut[j].Min_Sink_FeedbackTime;}
#endif                
#endif
#endif // SINK_CALC_DISTANCES

#ifdef RT_USE_TREECOL_FOR_NH
                int kbin=0; for(kbin=0; kbin < RT_USE_TREECOL_FOR_NH; kbin++) {P[place].ColumnDensityBins[kbin] += GravDataOut[j].ColumnDensityBins[kbin];}
#endif
#ifdef SINK_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
                P[place].MencInRcrit += GravDataOut[j].MencInRcrit;
#endif
#ifdef RT_OTVET
                if(P[place].Type==0) {int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {CellP[place].ET[k_freq] += GravDataOut[j].ET[k_freq];}}
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
                if(P[place].Type==0) {CellP[place].Rad_Flux_UV += GravDataOut[j].Rad_Flux_UV;}
                if(P[place].Type==0) {CellP[place].Rad_Flux_EUV += GravDataOut[j].Rad_Flux_EUV;}
#ifdef CHIMES
                if(P[place].Type == 0)
                {
                    int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++)
                    {
                        CellP[place].Chimes_G0[kc] += GravDataOut[j].Chimes_G0[kc];
                        CellP[place].Chimes_fluxPhotIon[kc] += GravDataOut[j].Chimes_fluxPhotIon[kc];
                    }
                }
#endif
#endif
#ifdef SINK_COMPTON_HEATING
                if(P[place].Type==0) CellP[place].Rad_Flux_AGN += GravDataOut[j].Rad_Flux_AGN;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                if(P[place].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[place].Rad_E_gamma[kf] += GravDataOut[j].Rad_E_gamma[kf];}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
                if(P[place].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP[place].Rad_Flux[kf] += GravDataOut[j].Rad_Flux[kf];}}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
                if(P[place].Type==0) {CellP[place].SubGrid_CosmicRayEnergyDensity += GravDataOut[j].SubGrid_CosmicRayEnergyDensity;}
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                P[place].tidal_tensorps += GravDataOut[j].tidal_tensorps;
#ifdef COMPUTE_JERK_IN_GRAVTREE
                {int i1tt; for(i1tt=0; i1tt<3; i1tt++) P[place].GravJerk[i1tt] += GravDataOut[j].GravJerk[i1tt];}
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                P[place].tidal_zeta += GravDataOut[j].tidal_zeta;
#endif
#endif
            }
            tend = my_second(); timetree1 += timediff(tstart, tend);
            myfree(GravDataOut); myfree(GravDataIn);

            if(NextParticle >= (int)ActiveParticleList.size()) {ndone_flag = 1;} else {ndone_flag = 0;} /* figure out if we are done with the particular active set here */
            tstart = my_second();
            MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /* call an allreduce to figure out if all tasks are also done here, otherwise we need to iterate */
            tend = my_second(); timewait2 += timediff(tstart, tend);
        }
        while(ndone < NTask);
    } /* Ewald_iter */
    myfree(DataNodeList); myfree(DataIndexTable);

#ifdef SINGLE_STAR_DIRECT_GRAVITY
    /* the exact star-star sum, which the tree walk above deliberately left out. Here and not later:
       the tree's own contributions are complete (imports included) but not yet multiplied by All.G,
       and the direct sum is in those same G-free units. */
    CPU_Step[CPU_TREEMISC] += measure_time();
    star_direct_gravity_build_table();
    star_direct_gravity_compute();
    star_direct_gravity_free_table();
    CPU_Step[CPU_TREEWALK1] += measure_time();
#endif

    /* assign node cost to particles */
    if(TakeLevel >= 0) {
        sum_top_level_node_costfactors();
        for(i = 0; i < NumPart; i++)
        {
            int no = Father[i];
            while(no >= 0)
            {
                /* Father[] is seeded to -1 per treebuild, so no>=0 means the recursion really did
                 * reach this particle and the chain should be walkable. An index outside the node
                 * array therefore means genuine corruption of u.d.father; say so instead of
                 * dereferencing it and dying in a way that points nowhere. */
                if(no < All.MaxPart || no >= All.MaxPart + MaxNodes)
                {
                    printf("task %d: corrupt tree father chain for particle %d (ID=%llu): node index %d outside [%d, %d)\n",
                           ThisTask, i, (unsigned long long) P[i].ID, no, All.MaxPart, All.MaxPart + MaxNodes);
                    fflush(stdout); endrun(918273);
                }
                if(Nodes[no].u.d.mass > 0) {P[i].GravCost[TakeLevel] += Nodes[no].GravCost * P[i].Mass / Nodes[no].u.d.mass;}
                no = Nodes[no].u.d.father;
            }
        }
    }


    /* now perform final operations on results [communication loop is done] */
#ifndef GRAVITY_HYBRID_OPENING_CRIT  // in collisional systems we don't want to rely on the relative opening criterion alone, because aold can be dominated by a binary companion but we still want accurate contributions from distant nodes. Thus we combine BH and relative criteria. - MYG
    if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS) {if(!(All.Ti_Current == 0 && RestartFlag == 0)) {if(All.TypeOfOpeningCriterion == 1) {All.ErrTolTheta = 0;}}} else {if(All.TypeOfOpeningCriterion == 1) {All.ErrTolTheta = 0;}} /* This will switch to the relative opening criterion for the following force computations */
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int ii = 0; ii < (int)ActiveParticleList.size(); ii++)
    {
        int i = ActiveParticleList[ii];
#ifdef HERMITE_INTEGRATION
        if(HermiteOnlyFlag) {if(!eligible_for_hermite(i)) continue;} /* if we are completing an extra loop required for the Hermite integration, all of the below would be double-calculated, so skip it */
#endif      
#ifdef ADAPTIVE_TREEFORCE_UPDATE
        double dt = get_particle_timestep_in_physical(i);
        if(treeforce_skip_flag[ii]) { // use cached decision from BEFORE tree walk to avoid mismatch if tree walk modified quantities that needs_new_treeforce depends on
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

        /* calculate 'old acceleration' for use in the relative tree-opening criterion */
        if(!(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0)) /* to prevent that we overwrite OldAcc in the first evaluation for 2lpt ICs */
            {
                auto accel = P[i].GravAccel;
#ifdef PMGRID
                accel += P[i].GravPM;
#endif
                P[i].OldAcc = accel.norm() / All.G; /* convert back to the non-G units for convenience to match units in loops assumed */
            }

#if (SINGLE_STAR_TIMESTEPPING > 0) /* Subtract component of force from companion if in binary, because we will operator-split this */
#ifdef KETJU_REGULARIZATION
        if((P[i].Type == 5) && (P[i].is_in_a_binary == 1) && !P[i].KetjuIntegrated) {subtract_companion_gravity(i);}
#else
        if((P[i].Type == 5) && (P[i].is_in_a_binary == 1)) {subtract_companion_gravity(i);}
#endif
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
        if((P[i].Type==0) && (P[i].Mass>0))
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
#endif

#if !defined(BOX_PERIODIC) && !defined(PMGRID) /* some factors here in case we are trying to do comoving simulations in a non-periodic box (special use cases) */
        if(All.ComovingIntegrationOn) {P[i].GravAccel += P[i].Pos * (0.5*All.OmegaMatter*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits);}
        if(All.ComovingIntegrationOn==0) {P[i].GravAccel += P[i].Pos * (All.OmegaLambda*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits);}
#ifdef EVALPOTENTIAL
        if(All.ComovingIntegrationOn) {P[i].Potential -= 0.5*All.OmegaMatter*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits * P[i].Pos.norm_sq();}
#endif
#endif

    } /* end of loop over active particles*/


#else /* the tree walk that would normally zero-then-accumulate these never runs here, so they would
         otherwise hold whatever was in the malloc slot -- and consumers read them regardless of whether
         gravity is compiled in (HERMITE_INTEGRATION's do_the_kick() copies GravAccel/GravJerk into
         Hermite_OldAcc/OldJerk every step). GravAccel is zeroed again shortly after this, inside
         add_analytic_gravitational_forces(); that is deliberate redundancy, so this block does not
         depend on what that routine happens to do. */
    for (int i : ActiveParticleList) {
        P[i].GravAccel = {};
#ifdef COMPUTE_JERK_IN_GRAVTREE
        P[i].GravJerk = {};
#endif
#ifdef EVALPOTENTIAL
        P[i].Potential = 0;
#endif
    }
#endif /* end tree walk (guarded by SELFGRAVITY_OFF unless non-gravity tree features are active) */


    add_analytic_gravitational_forces(); /* add analytic terms, which -CAN- be enabled even if self-gravity is not */


    /* Now the force computation is finished: gather timing and diagnostic information */
    t1 = WallclockTime = my_second(); timeall = timediff(t0, t1);
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
    MPI_Reduce(&Numnodestree, &maxnumnodes, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    CPU_Step[CPU_TREEMISC] += timeall - (timetree + timewait + timecomm);
    CPU_Step[CPU_TREEWALK1] += timetree1; CPU_Step[CPU_TREEWALK2] += timetree2;
    CPU_Step[CPU_TREESEND] += timecommsumm1; CPU_Step[CPU_TREERECV] += timecommsumm2;
    CPU_Step[CPU_TREEWAIT1] += timewait1; CPU_Step[CPU_TREEWAIT2] += timewait2;
#ifdef OUTPUT_ADDITIONAL_RUNINFO
    if(ThisTask == 0)
    {
        fprintf(FdTimings, "Step= %lld  t= %.16g  dt= %.16g \n",(long long) All.NumCurrentTiStep, All.Time, All.TimeStep);
        fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g (%g) iter= %d\n", (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000), n_exported / ((double) GlobNumForceUpdate), N_nodesinlist / ((double) n_exported + 1.0e-10), iter); /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */
        fprintf(FdTimings, "work-load balance: %g (%g %g) rel1to2=%g   max=%g avg=%g\n", maxt / (1.0e-6 + sumt / NTask), maxt1 / (1.0e-6 + sumt1 / NTask), maxt2 / (1.0e-6 + sumt2 / NTask), sumt1 / (1.0e-6 + sumt1 + sumt2), maxt, sumt / NTask);
        fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
        fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes, maxnumnodes / (All.TreeAllocFactor * All.MaxPart + NTopnodes));
        fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", GlobNumForceUpdate / (sumt + 1.0e-20), GlobNumForceUpdate / (1.0e-6 + maxt * NTask), ((double) (sum_costtotal)) / (1.0e-20 + GlobNumForceUpdate), ((double) ewaldtot) / (1.0e-20 + GlobNumForceUpdate)); fprintf(FdTimings, "\n");
        fflush(FdTimings);
    }
    double costtotal_new = 0, sum_costtotal_new;
    if(TakeLevel >= 0)
    {
        for(i = 0; i < NumPart; i++) {costtotal_new += P[i].GravCost[TakeLevel];}
        MPI_Reduce(&costtotal_new, &sum_costtotal_new, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(sum_costtotal>0) {PRINT_STATUS(" ..relative error in the total number of tree-gravity interactions = %g", (sum_costtotal - sum_costtotal_new) / sum_costtotal);} /* can be non-zero if THREAD_SAFE_COSTS is not used (and due to round-off errors). */
    }
#endif
    CPU_Step[CPU_TREEMISC] += measure_time();
}




void *gravity_primary_loop(void *p)
{
    int i, j, ret, thread_id = *(int *) p, *exportflag, *exportnodecount, *exportindex;
    exportflag = Exportflag + thread_id * NTask; exportnodecount = Exportnodecount + thread_id * NTask; exportindex = Exportindex + thread_id * NTask;
    for(j = 0; j < NTask; j++) {exportflag[j] = -1;} /* Note: exportflag is local to each thread */
#ifdef _OPENMP
    if(BufferCollisionFlag && thread_id) {return NULL;} /* force to serial for this subloop if threads simultaneously cross the Nexport bunchsize threshold */
#endif
#ifndef GRAVITY_PRIMARY_LOOP_BATCH_SIZE
#define GRAVITY_PRIMARY_LOOP_BATCH_SIZE 8
#endif
    while(1)
    {
        int batch[GRAVITY_PRIMARY_LOOP_BATCH_SIZE], batch_count = 0;
#ifdef _OPENMP
#pragma omp critical(_nextlistgravprim_)
#endif
        {
            while(batch_count < GRAVITY_PRIMARY_LOOP_BATCH_SIZE && BufferFullFlag == 0 && NextParticle < (int)ActiveParticleList.size())
            {
                int idx = ActiveParticleList[NextParticle]; NextParticle++;
                if(!ProcessedFlag[idx]) {batch[batch_count++] = idx;}
            }
        }
        if(batch_count == 0) {break;}
        int buffer_full = 0;
        for(int b = 0; b < batch_count; b++)
        {
            i = batch[b];
#ifdef HERMITE_INTEGRATION /* if we are in the Hermite extra loops and a particle is not flagged for this, simply mark it done and move on */
            if(HermiteOnlyFlag && !eligible_for_hermite(i)) {ProcessedFlag[i]=1; continue;}
#endif
#ifdef ADAPTIVE_TREEFORCE_UPDATE
            if(!needs_new_treeforce(i)) {ProcessedFlag[i]=1; continue;}
#endif

#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
            if(Ewald_iter)
            {
                ret = force_treeevaluate_ewald_correction(i, 0, exportflag, exportnodecount, exportindex);
                if(ret >= 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    Ewaldcount += ret;
                } else {buffer_full = 1; break;}
            }
            else
#endif
            {
                ret = force_treeevaluate(i, 0, exportflag, exportnodecount, exportindex);
                if(ret < 0) {buffer_full = 1; break;}
#ifdef _OPENMP
#pragma omp atomic
#endif
                Costtotal += ret;
            }
            ProcessedFlag[i] = 1;
        }
        if(buffer_full) {break;}
    } // while loop
    return NULL;
}


void *gravity_secondary_loop(void *p)
{
    int j, nodesinlist, dummy, ret;
#ifndef GRAVITY_SECONDARY_LOOP_BATCH_SIZE
#define GRAVITY_SECONDARY_LOOP_BATCH_SIZE 8
#endif
    while(1)
    {
        int batch[GRAVITY_SECONDARY_LOOP_BATCH_SIZE], batch_count = 0;
#ifdef _OPENMP
#pragma omp critical(_nextlistgravsec_)
#endif
        {
            while(batch_count < GRAVITY_SECONDARY_LOOP_BATCH_SIZE && NextJ < Nimport)
            {
                batch[batch_count++] = NextJ; NextJ++;
            }
        }
        if(batch_count == 0) {break;}
        for(int b = 0; b < batch_count; b++)
        {
            j = batch[b];
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
            if(Ewald_iter)
            {
                int cost = force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy, &dummy);
#ifdef _OPENMP
#pragma omp atomic
#endif
                Ewaldcount += cost;
            }
            else
#endif
            {
                ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy);
#ifdef _OPENMP
#pragma omp atomic
#endif
                N_nodesinlist += nodesinlist;
#ifdef _OPENMP
#pragma omp atomic
#endif
                Costtotal += ret;
            }
        }
    }
    return NULL;
}


void sum_top_level_node_costfactors(void)
{
    double *costlist = (double*)mymalloc("costlist", NTopnodes * sizeof(double));
    double *costlist_all = (double*)mymalloc("costlist_all", NTopnodes * sizeof(double));
    int i; for(i = 0; i < NTopnodes; i++) {costlist[i] = Nodes[All.MaxPart + i].GravCost;}
    MPI_Allreduce(costlist, costlist_all, NTopnodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(i = 0; i < NTopnodes; i++) {Nodes[All.MaxPart + i].GravCost = costlist_all[i];}
    myfree(costlist_all); myfree(costlist);
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


/*! This function is used as a comparison kernel in a sort routine. It is used to group particles in the communication buffer that are going to be sent to the same CPU */
int data_index_compare(const void *a, const void *b)
{
    if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task)) {return -1;}
    if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task)) {return +1;}
    if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index)) {return -1;}
    if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index)) {return +1;}
    if(((struct data_index *) a)->IndexGet < (((struct data_index *) b)->IndexGet)) {return -1;}
    if(((struct data_index *) a)->IndexGet > (((struct data_index *) b)->IndexGet)) {return +1;}
    return 0;
}


static void msort_dataindex_with_tmp(struct data_index *b, size_t n, struct data_index *t)
{
    if(n <= 1) {return;}
    struct data_index *tmp;
    struct data_index *b1, *b2;
    size_t n1, n2;
    n1 = n / 2;
    n2 = n - n1;
    b1 = b;
    b2 = b + n1;
    msort_dataindex_with_tmp(b1, n1, t);
    msort_dataindex_with_tmp(b2, n2, t);
    tmp = t;
    while(n1 > 0 && n2 > 0)
    {
        if(b1->Task < b2->Task || (b1->Task == b2->Task && b1->Index <= b2->Index))
        {
            --n1;
            *tmp++ = *b1++;
        }
        else
        {
            --n2;
            *tmp++ = *b2++;
        }
    }
    if(n1 > 0) {memcpy(tmp, b1, n1 * sizeof(struct data_index));}
    memcpy(b, t, (n - n2) * sizeof(struct data_index));
}


void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
    const size_t size = n * s;
    struct data_index *tmp = (struct data_index *) mymalloc("struct data_index *tmp", size);
    msort_dataindex_with_tmp((struct data_index *) b, n, tmp);
    myfree(tmp);
}


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
#ifdef HERMITE_INTEGRATION
    if((1 << P[n].Type) & HERMITE_INTEGRATION) {return 1;} // Hermite-integrated particles must always get a fresh tree force; the lazy-cache path is incompatible with the predictor-corrector sub-stepping
#endif
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
