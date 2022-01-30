#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <gsl/gsl_eigen.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#include "./analytic_gravity.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
pthread_mutex_t mutex_nexport;
pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

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
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * The code has been modified substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
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
        PRINT_STATUS(" ..Tree construction done.");
    }

    CPU_Step[CPU_TREEMISC] += measure_time(); t0 = my_second();
#ifndef SELFGRAVITY_OFF
    /* allocate buffers to arrange communication */
    PRINT_STATUS(" ..Begin tree force. (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
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

    /* begin main communication and tree-walk loop. note the ewald-iter terms here allow for multiple iterations for periodic-tree corrections if needed */
    for(Ewald_iter = 0; Ewald_iter <= ewald_max; Ewald_iter++)
    {
        NextParticle = FirstActiveParticle;	/* begin with this index */
        do /* primary point-element loop */
        {
            iter++;
            BufferFullFlag = 0; Nexport = 0; save_NextParticle = NextParticle; tstart = my_second();

#ifdef PTHREADS_NUM_THREADS
            pthread_t mythreads[PTHREADS_NUM_THREADS - 1]; int threadid[PTHREADS_NUM_THREADS - 1];
            pthread_attr_t attr; pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            pthread_mutex_init(&mutex_nexport, NULL);
            pthread_mutex_init(&mutex_partnodedrift, NULL); TimerFlag = 0;
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            {
                threadid[j] = j + 1;
                pthread_create(&mythreads[j], &attr, gravity_primary_loop, &threadid[j]);
            }
#endif
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
#ifdef PTHREADS_NUM_THREADS
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) pthread_join(mythreads[j], NULL);
#endif
            tend = my_second(); timetree1 += timediff(tstart, tend);

            if(BufferFullFlag) /* we've filled the buffer or reached the end of the list, prepare for communications */
            {
                int last_nextparticle = NextParticle; NextParticle = save_NextParticle;
                while(NextParticle >= 0)
                {
                    if(NextParticle == last_nextparticle) {break;}
                    if(ProcessedFlag[NextParticle] != 1) {break;}
                    ProcessedFlag[NextParticle] = 2; NextParticle = NextActiveParticle[NextParticle];
                }
                if(NextParticle == save_NextParticle) {endrun(114408);} /* in this case, the buffer is too small to process even a single particle */

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
            MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare); /* construct export count tables */
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
                GravDataIn[j].Type = P[place].Type;
                GravDataIn[j].OldAcc = P[place].OldAcc;
                for(k = 0; k < 3; k++) {GravDataIn[j].Pos[k] = P[place].Pos[k];}
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(SINGLE_STAR_TIMESTEPPING)
                GravDataIn[j].Mass = P[place].Mass;
#endif
#if defined(BH_DYNFRICTION_FROMTREE)
                if(P[place].Type==5) {GravDataIn[j].BH_Mass = P[place].BH_Mass;}
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(BH_DYNFRICTION_FROMTREE)
                for(k = 0; k < 3; k++) {GravDataIn[j].Vel[k] = P[place].Vel[k];}
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
                if(P[place].Type == 5)
                {
                    GravDataIn[j].min_bh_t_orbital = P[place].min_bh_t_orbital; //orbital time for binary
                    GravDataIn[j].comp_Mass = P[place].comp_Mass; //mass of binary companion
                    GravDataIn[j].is_in_a_binary = P[place].is_in_a_binary; // 1 if we're in a binary, 0 if not
                    for(k=0;k<3;k++) {GravDataIn[j].comp_dx[k]=P[place].comp_dx[k]; GravDataIn[j].comp_dv[k]=P[place].comp_dv[k];}
                }
                else {P[place].is_in_a_binary=0; /* setting values to zero just to be sure */}
#endif
#if defined(SINGLE_STAR_TIMESTEPPING)
                GravDataIn[j].Soft = All.ForceSoftening[P[place].Type];
#endif
#if defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(FLAG_NOT_IN_PUBLIC_CODE)
                if( (P[place].Type == 0) && (PPP[place].Hsml > All.ForceSoftening[P[place].Type]) ) {GravDataIn[j].Soft = PPP[place].Hsml;} else {GravDataIn[j].Soft = All.ForceSoftening[P[place].Type];}
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                if((P[place].Type == 0) && (PPP[place].Hsml > All.ForceSoftening[P[place].Type])) {GravDataIn[j].AGS_zeta = PPPZ[place].AGS_zeta;} else {GravDataIn[j].AGS_zeta = 0;}
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                GravDataIn[j].Soft = PPP[place].AGS_Hsml;
                GravDataIn[j].AGS_zeta = PPPZ[place].AGS_zeta;
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
#ifdef PTHREADS_NUM_THREADS
                for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) pthread_create(&mythreads[j], &attr, gravity_secondary_loop, &threadid[j]);
#endif
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
#ifdef PTHREADS_NUM_THREADS
                for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
                pthread_mutex_destroy(&mutex_partnodedrift); pthread_mutex_destroy(&mutex_nexport); pthread_attr_destroy(&attr);
#endif
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
                for(k=0;k<3;k++) {P[place].GravAccel[k] += GravDataOut[j].Acc[k];}
                if(Ewald_iter > 0) continue; /* everything below is ONLY evaluated if we are in the first sub-loop, not the periodic correction, or else we will get un-allocated memory or un-physical values */

#ifdef EVALPOTENTIAL
                P[place].Potential += GravDataOut[j].Potential;
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
                P[place].TreeMass += GravDataOut[j].TreeMass;
#endif
#ifdef BH_CALC_DISTANCES /* GravDataOut[j].min_dist_to_bh contains the min dist to particle "P[place]" on another task.  We now check if it is smaller than the current value */
                if(GravDataOut[j].min_dist_to_bh < P[place].min_dist_to_bh)
                {
                    P[place].min_dist_to_bh = GravDataOut[j].min_dist_to_bh;
                    P[place].min_xyz_to_bh[0] = GravDataOut[j].min_xyz_to_bh[0];
                    P[place].min_xyz_to_bh[1] = GravDataOut[j].min_xyz_to_bh[1];
                    P[place].min_xyz_to_bh[2] = GravDataOut[j].min_xyz_to_bh[2];
                }
#ifdef SINGLE_STAR_TIMESTEPPING
                if(GravDataOut[j].min_bh_approach_time < P[place].min_bh_approach_time) {P[place].min_bh_approach_time = GravDataOut[j].min_bh_approach_time;}
                if(GravDataOut[j].min_bh_freefall_time < P[place].min_bh_freefall_time) {P[place].min_bh_freefall_time = GravDataOut[j].min_bh_freefall_time;}
#ifdef SINGLE_STAR_FIND_BINARIES
                if((P[place].Type == 5) && (GravDataOut[j].min_bh_t_orbital < P[place].min_bh_t_orbital))
                {
                    P[place].min_bh_t_orbital = GravDataOut[j].min_bh_t_orbital;
                    P[place].comp_Mass = GravDataOut[j].comp_Mass;
                    P[place].is_in_a_binary = GravDataOut[j].is_in_a_binary;
                    for(k=0;k<3;k++) {P[place].comp_dx[k]=GravDataOut[j].comp_dx[k]; P[place].comp_dv[k]=GravDataOut[j].comp_dv[k];}
                }
#endif
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
                if(GravDataOut[j].min_bh_fb_time < P[place].min_bh_fb_time) {P[place].min_bh_fb_time = GravDataOut[j].min_bh_fb_time;}
#endif                
#endif
#endif // BH_CALC_DISTANCES

#ifdef RT_USE_TREECOL_FOR_NH
                int kbin=0; for(kbin=0; kbin < RT_USE_TREECOL_FOR_NH; kbin++) {P[place].ColumnDensityBins[kbin] += GravDataOut[j].ColumnDensityBins[kbin];}
#endif
#ifdef BH_SEED_FROM_LOCALGAS_TOTALMENCCRITERIA
                P[place].MencInRcrit += GravDataOut[j].MencInRcrit;
#endif
#ifdef RT_OTVET
                if(P[place].Type==0) {int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) for(k=0;k<6;k++) SphP[place].ET[k_freq][k] += GravDataOut[j].ET[k_freq][k];}
#endif
#ifdef BH_COMPTON_HEATING
                if(P[place].Type==0) SphP[place].Rad_Flux_AGN += GravDataOut[j].Rad_Flux_AGN;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                if(P[place].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {SphP[place].Rad_E_gamma[kf] += GravDataOut[j].Rad_E_gamma[kf];}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
                if(P[place].Type==0) {int kf,k2; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k2=0;k2<3;k2++) {SphP[place].Rad_Flux[kf][k2] += GravDataOut[j].Rad_Flux[kf][k2];}}}
#endif
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                {int i1tt,i2tt; for(i1tt=0;i1tt<3;i1tt++) {for(i2tt=0;i2tt<3;i2tt++) {P[place].tidal_tensorps[i1tt][i2tt] += GravDataOut[j].tidal_tensorps[i1tt][i2tt];}}}
#ifdef COMPUTE_JERK_IN_GRAVTREE
                {int i1tt; for(i1tt=0; i1tt<3; i1tt++) P[place].GravJerk[i1tt] += GravDataOut[j].GravJerk[i1tt];}
#endif
#endif
            }
            tend = my_second(); timetree1 += timediff(tstart, tend);
            myfree(GravDataOut); myfree(GravDataIn);

            if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;} /* figure out if we are done with the particular active set here */
            tstart = my_second();
            MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /* call an allreduce to figure out if all tasks are also done here, otherwise we need to iterate */
            tend = my_second(); timewait2 += timediff(tstart, tend);
        }
        while(ndone < NTask);
    } /* Ewald_iter */
    myfree(DataNodeList); myfree(DataIndexTable);

    /* assign node cost to particles */
    if(TakeLevel >= 0) {
        sum_top_level_node_costfactors();
        for(i = 0; i < NumPart; i++)
        {
            int no = Father[i];
            while(no >= 0)
            {
                if(Nodes[no].u.d.mass > 0) {P[i].GravCost[TakeLevel] += Nodes[no].GravCost * P[i].Mass / Nodes[no].u.d.mass;}
                no = Nodes[no].u.d.father;
            }
        }
    }


    /* now perform final operations on results [communication loop is done] */
#ifndef GRAVITY_HYBRID_OPENING_CRIT  // in collisional systems we don't want to rely on the relative opening criterion alone, because aold can be dominated by a binary companion but we still want accurate contributions from distant nodes. Thus we combine BH and relative criteria. - MYG
    if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS) {if(!(All.Ti_Current == 0 && RestartFlag == 0)) {if(All.TypeOfOpeningCriterion == 1) {All.ErrTolTheta = 0;}}} else {if(All.TypeOfOpeningCriterion == 1) {All.ErrTolTheta = 0;}} /* This will switch to the relative opening criterion for the following force computations */
#endif
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef HERMITE_INTEGRATION
        if(HermiteOnlyFlag) {if(!eligible_for_hermite(i)) continue;} /* if we are completing an extra loop required for the Hermite integration, all of the below would be double-calculated, so skip it */
#endif      
#ifdef ADAPTIVE_TREEFORCE_UPDATE
        double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
        if(!needs_new_treeforce(i)) { // if we don't yet need a new tree pass, just update GravAccel according to the jerk term, increment the counter, and go to the next particle           
            for(j=0; j<3; j++) {P[i].GravAccel[j] += dt * P[i].GravJerk[j] * All.cf_a2inv;} // a^-1 from converting velocity term in the jerk to physical; a^-3 from the 1/r^3; a^2 from converting the physical dt * j increment to GravAccel back to the units for GravAccel; result is a^-2; note that Ewald and PMGRID terms are neglected from the jerk at present
            P[i].time_since_last_treeforce += dt;
            continue;
        } else {
            P[i].time_since_last_treeforce = dt;
        }
#endif
        /* before anything: multiply by G for correct units [be sure operations above/below are aware of this!] */
        for(j=0;j<3;j++) {P[i].GravAccel[j] *= All.G;}        
#if (SINGLE_STAR_TIMESTEPPING > 0)
        for(j=0;j<3;j++) {P[i].COM_GravAccel[j] *= All.G;}
#endif

#ifdef EVALPOTENTIAL
        P[i].Potential *= All.G;
#ifdef BOX_PERIODIC
        if(All.ComovingIntegrationOn) {P[i].Potential -= All.G * 2.8372975 * pow(P[i].Mass, 2.0 / 3) * pow(All.OmegaMatter * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G), 1.0 / 3);} else {if(All.OmegaLambda>0) {P[i].Potential -= 0.5*All.OmegaLambda*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits * (P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2]);}}
#endif
#ifdef PMGRID
        P[i].Potential += P[i].PM_Potential; /* add in long-range potential */
#endif
#endif
#ifdef COUNT_MASS_IN_GRAVTREE
        P[i].TreeMass += P[i].Mass;
        if(P[i].Type == 5) printf("Particle %d sees mass %g in the gravity tree\n", P[i].ID, P[i].TreeMass);
#endif

        /* calculate 'old acceleration' for use in the relative tree-opening criterion */
        if(!(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0)) /* to prevent that we overwrite OldAcc in the first evaluation for 2lpt ICs */
            {
                ax=P[i].GravAccel[0]; ay=P[i].GravAccel[1]; az=P[i].GravAccel[2];
#ifdef PMGRID
                ax+=P[i].GravPM[0]; ay+=P[i].GravPM[1]; az+=P[i].GravPM[2];
#endif
                P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az) / All.G; /* convert back to the non-G units for convenience to match units in loops assumed */
            }

#if (SINGLE_STAR_TIMESTEPPING > 0) /* Subtract component of force from companion if in binary, because we will operator-split this */
        if((P[i].Type == 5) && (P[i].is_in_a_binary == 1)) {subtract_companion_gravity(i);}
#endif

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* final operations to compute the diagonalized tidal tensor and related quantities */
#if (defined(TIDAL_TIMESTEP_CRITERION) || defined(GALSF_SFR_TIDAL_HILL_CRITERION)) // diagonalize the tidal tensor so we can use its invariants, which don't change with rotation
        double tt[9]; for(j=0; j<3; j++) {for (k=0; k<3; k++) tt[3*j+k] = P[i].tidal_tensorps[j][k];}
#ifdef PMGRID
        for(j=0; j<3; j++) {for (k=0; k<3; k++) tt[3*j+k] += P[i].tidal_tensorpsPM[j][k];}
#endif
        gsl_matrix_view m = gsl_matrix_view_array (tt, 3, 3);
        gsl_vector *eval = gsl_vector_alloc (3);
        gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (3);
        gsl_eigen_symm(&m.matrix, eval,  w);
        for(k=0; k<3; k++) P[i].tidal_tensorps[k][k] = gsl_vector_get(eval,k); // set diagonal elements to eigenvalues
        P[i].tidal_tensorps[0][1] = P[i].tidal_tensorps[1][0] = P[i].tidal_tensorps[1][2] = P[i].tidal_tensorps[2][1] = P[i].tidal_tensorps[0][2] = P[i].tidal_tensorps[2][0] = 0; //zero out off-diagonal elements
        gsl_eigen_symm_free(w); gsl_vector_free(eval);
#endif
#ifdef GDE_DISTORTIONTENSOR /* for GDE implementation, want to include particle self-tide contribution -- for timestep or hill criteria, on the other hand, this is not necessary */
        if(All.ComovingIntegrationOn) {P[i].tidal_tensorps[0][0] -= All.TidalCorrection/All.G; P[i].tidal_tensorps[1][1] -= All.TidalCorrection/All.G; P[i].tidal_tensorps[2][2] -= All.TidalCorrection/All.G;} // subtract Hubble flow terms //
        /* Diagonal terms of tidal tensor need correction, because tree is running over all particles -> also over target particle -> extra term -> correct it */
        P[i].tidal_tensorps[0][0] += P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type]) * 10.666666666667;
        P[i].tidal_tensorps[1][1] += P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type]) * 10.666666666667;
        P[i].tidal_tensorps[2][2] += P[i].Mass / (All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type] * All.ForceSoftening[P[i].Type]) * 10.666666666667;
#endif
        for(j=0;j<3;j++) {int i2tt; for(i2tt=0;i2tt<3;i2tt++) {P[i].tidal_tensorps[j][i2tt] *= All.G;}} // units //
#ifdef COMPUTE_JERK_IN_GRAVTREE
        for(j=0;j<3;j++) P[i].GravJerk[j] *= All.G;
#endif
#endif /* COMPUTE_TIDAL_TENSOR_IN_GRAVTREE */

#if defined(RT_OTVET) /* normalize the Eddington tensors we just calculated by walking the tree (normalize to trace=1) */
        if(P[i].Type == 0) {
            int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
            {double trace = SphP[i].ET[k_freq][0] + SphP[i].ET[k_freq][1] + SphP[i].ET[k_freq][2];
                if(!isnan(trace) && (trace>0)) {for(k=0;k<6;k++) {SphP[i].ET[k_freq][k]/=trace;}} else {SphP[i].ET[k_freq][0]=SphP[i].ET[k_freq][1]=SphP[i].ET[k_freq][2]=1./3.; SphP[i].ET[k_freq][3]=SphP[i].ET[k_freq][4]=SphP[i].ET[k_freq][5]=0;}}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) /* normalize to energy density with C, and multiply by volume to use standard 'finite volume-like' quantity as elsewhere in-code */
        if(P[i].Type==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {SphP[i].Rad_E_gamma[kf] *= P[i].Mass/(SphP[i].Density*All.cf_a3inv * C_LIGHT_CODE_REDUCED);}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX) /* multiply by volume to use standard 'finite volume-like' quantity as elsewhere in-code */
        if(P[i].Type==0) {int kf,k2; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k2=0;k2<3;k2++) {SphP[i].Rad_Flux[kf][k2] *= P[i].Mass/(SphP[i].Density*All.cf_a3inv);}}} // convert to standard finite-volume-like units //
#if !defined(RT_DISABLE_RAD_PRESSURE) // if we save the fluxes, we didnt apply forces on-the-spot, which means we appky them here //
        if((P[i].Type==0) && (P[i].Mass>0))
        {
            int k,kfreq; double vol_inv=SphP[i].Density*All.cf_a3inv/P[i].Mass, radacc[3]={0}, h_i=Get_Particle_Size(i)*All.cf_atime, sigma_eff_i=P[i].Mass/(h_i*h_i);
            for(kfreq=0; kfreq<N_RT_FREQ_BINS; kfreq++)
            {
                double f_slab=1, erad_i=0, vel_i[3]={0}, vdot_h[3]={0}, flux_i[3]={0}, flux_mag2=MIN_REAL_NUMBER, vdotflux=0, kappa_rad=rt_kappa(i,kfreq), tau_eff=kappa_rad*sigma_eff_i; if(tau_eff > 1.e-4) {f_slab = (1.-exp(-tau_eff)) / tau_eff;} // account for optically thick local 'slabs' self-shielding themselves
                double acc_norm = kappa_rad * f_slab / C_LIGHT_CODE_REDUCED; // pre-factor for radiation pressure acceleration
#if defined(RT_LEBRON)
                acc_norm *= All.PhotonMomentum_Coupled_Fraction; // allow user to arbitrarily increase/decrease strength of RP forces for testing
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
                erad_i = SphP[i].Rad_E_gamma_Pred[kfreq]*vol_inv; // if can, include the O[v/c] terms
#endif
                for(k=0;k<3;k++) {flux_i[k]=SphP[i].Rad_Flux_Pred[kfreq][k]*vol_inv; flux_mag2+=flux_i[k]*flux_i[k]; vel_i[k]=SphP[i].VelPred[k]/All.cf_atime; vdotflux+=vel_i[k]*flux_i[k];} // initialize a bunch of variables we will need
                for(k=0;k<3;k++) {vdot_h[k] = erad_i * (vel_i[k] + vdotflux*flux_i[k]/flux_mag2);} // calculate volume integral of scattering coefficient t_inv * (gas_vel . [e_rad*I + P_rad_tensor]), which gives an additional time-derivative term. this is the P term //
                for(k=0;k<3;k++) {radacc[k] += acc_norm * (flux_i[k] - vdot_h[k]);} // note these 'vdoth' terms shouldn't be included in FLD, since its really assuming the entire right-hand-side of the flux equation reaches equilibrium with the pressure tensor, which gives the expression in rt_utilities
            }
            for(k=0;k<3;k++) {P[i].GravAccel[k] += radacc[k] / All.cf_a2inv;} // convert into our code units for GravAccel, which are comoving gm/r^2 units //
        }
#endif
#endif

#ifdef RT_USE_TREECOL_FOR_NH  /* compute the effective column density that gives equivalent attenuation of a uniform background: -log(avg(exp(-tau)))/kappa */
        double attenuation=0; int kbin; // first do a sum of the columns and express columns in units of that sum, so that we're plugging O(1) values into exp and avoid overflow when we have unfortunate units. Then we just multiply by the sum at the end.
        double kappa_photoelectric = 500. * DMAX(1e-4, (P[i].Metallicity[0]/All.SolarAbundances[0])*return_dust_to_metals_ratio_vs_solar(i)); // dust opacity in cgs
        for(kbin=0; kbin<RT_USE_TREECOL_FOR_NH; kbin++) {attenuation += exp(-P[i].ColumnDensityBins[kbin] * UNIT_SURFDEN_IN_CGS * kappa_photoelectric);}
        P[i].SigmaEff = -log(attenuation/RT_USE_TREECOL_FOR_NH) / (kappa_photoelectric * UNIT_SURFDEN_IN_CGS);       
#endif

#if !defined(BOX_PERIODIC) && !defined(PMGRID) /* some factors here in case we are trying to do comoving simulations in a non-periodic box (special use cases) */
        if(All.ComovingIntegrationOn) {for(j=0;j<3;j++) {P[i].GravAccel[j] += 0.5*All.OmegaMatter *All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits * P[i].Pos[j];}}
        if(All.ComovingIntegrationOn==0) {for(j=0;j<3;j++) {P[i].GravAccel[j] += All.OmegaLambda*All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits * P[i].Pos[j];}}
#ifdef EVALPOTENTIAL
        if(All.ComovingIntegrationOn) {for(j=0;j<3;j++) {P[i].Potential -= 0.5*All.OmegaMatter *All.Hubble_H0_CodeUnits*All.Hubble_H0_CodeUnits * P[i].Pos[j]*P[i].Pos[j];}}
#endif
#endif

    } /* end of loop over active particles*/


#endif /* end SELFGRAVITY operations (check if SELFGRAVITY_OFF not enabled) */


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
#ifndef IO_REDUCED_MODE
    if(ThisTask == 0)
    {
        fprintf(FdTimings, "Step= %lld  t= %g  dt= %g \n",(long long) All.NumCurrentTiStep, All.Time, All.TimeStep);
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

    while(1)
    {
        int exitFlag = 0;
        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
        {
        if(BufferFullFlag != 0 || NextParticle < 0) {exitFlag=1;}
            else {i=NextParticle; ProcessedFlag[i]=0; NextParticle=NextActiveParticle[NextParticle];}
        }
        UNLOCK_NEXPORT;
        if(exitFlag) {break;}

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
            if(ret >= 0) {Ewaldcount += ret; /* note: ewaldcount may be slightly incorrect for multiple threads if buffer gets filled up */} else {break; /* export buffer has filled up */}
        }
        else
#endif
        {
            ret = force_treeevaluate(i, 0, exportflag, exportnodecount, exportindex);
            if(ret < 0) {break;} /* export buffer has filled up */
            Costtotal += ret;
        }
        ProcessedFlag[i] = 1;	/* particle successfully finished */
    } // while loop
    return NULL;
}


void *gravity_secondary_loop(void *p)
{
    int j, nodesinlist, dummy, ret;
    while(1)
    {
        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
        {
            j = NextJ;
            NextJ++;
        }
        UNLOCK_NEXPORT;
        if(j >= Nimport) {break;}

#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(PMGRID)
        if(Ewald_iter)
        {
            int cost = force_treeevaluate_ewald_correction(j, 1, &dummy, &dummy, &dummy);
            Ewaldcount += cost;
        }
        else
#endif
        {
            ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy);
            N_nodesinlist += nodesinlist; Costtotal += ret;
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


/*! This function sets the (comoving) softening length of all particle types in the table All.SofteningTable[...].
 We check that the physical softening length is bounded by the Softening-MaxPhys values */
void set_softenings(void)
{
    if(All.ComovingIntegrationOn)
    {
        if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys) {All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;} else {All.SofteningTable[0] = All.SofteningGas;}
        if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys) {All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;} else {All.SofteningTable[1] = All.SofteningHalo;}
        if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys) {All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;} else {All.SofteningTable[2] = All.SofteningDisk;}
        if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys) {All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;} else {All.SofteningTable[3] = All.SofteningBulge;}
        if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys) {All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;} else {All.SofteningTable[4] = All.SofteningStars;}
        if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys) {All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;} else {All.SofteningTable[5] = All.SofteningBndry;}
    }
    else
    {
        All.SofteningTable[0] = All.SofteningGas;
        All.SofteningTable[1] = All.SofteningHalo;
        All.SofteningTable[2] = All.SofteningDisk;
        All.SofteningTable[3] = All.SofteningBulge;
        All.SofteningTable[4] = All.SofteningStars;
        All.SofteningTable[5] = All.SofteningBndry;
    }
    int i; for(i = 0; i < 6; i++) {All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];} 
    /* set the minimum gas kernel length to be used this timestep */
    All.MinHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
#ifndef SELFGRAVITY_OFF
    if(All.MinHsml <= 5.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0])
        {All.MinHsml = 5.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0];}
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
    double u, dr, fac, fac2, h, h_inv, h3_inv, u2, tidal_tensorps[3][3]; int i1, i2;
    dr = sqrt(P[i].comp_dx[0]*P[i].comp_dx[0] + P[i].comp_dx[1]*P[i].comp_dx[1] + P[i].comp_dx[2]*P[i].comp_dx[2]);
    h = All.ForceSoftening[5];  h_inv = 1.0 / h; h3_inv = h_inv*h_inv*h_inv; u = dr*h_inv; u2=u*u;
    fac = P[i].comp_Mass / (dr*dr*dr); fac2 = 3.0 * P[i].comp_Mass / (dr*dr*dr*dr*dr); /* no softening nonsense */
    if(dr < h) /* second derivatives needed -> calculate them from softened potential. NOTE this is here -assuming- a cubic spline, will be inconsistent for different kernels used! */
    {
	    fac = P[i].comp_Mass * kernel_gravity(u, h_inv, h3_inv, 1);
        fac2 = P[i].comp_Mass * kernel_gravity(u, h_inv, h3_inv, 2);
    }
    for(i1=0;i1<3;i1++) {P[i].COM_GravAccel[i1] = P[i].GravAccel[i1] - P[i].comp_dx[i1] * fac * All.G;} /* this assumes the 'G' has been put into the units for the grav accel */

    /* Adjusting tidal tensor according to terms above */
    tidal_tensorps[0][0] = P[i].tidal_tensorps[0][0] - (-fac + P[i].comp_dx[0] * P[i].comp_dx[0] * fac2);
    tidal_tensorps[0][1] = P[i].tidal_tensorps[0][1] - (P[i].comp_dx[0] * P[i].comp_dx[1] * fac2);
    tidal_tensorps[0][2] = P[i].tidal_tensorps[0][2] - (P[i].comp_dx[0] * P[i].comp_dx[2] * fac2);
    tidal_tensorps[1][1] = P[i].tidal_tensorps[1][1] - (-fac + P[i].comp_dx[1] * P[i].comp_dx[1] * fac2);
    tidal_tensorps[1][2] = P[i].tidal_tensorps[1][2] - (P[i].comp_dx[1] * P[i].comp_dx[2] * fac2);
    tidal_tensorps[2][2] = P[i].tidal_tensorps[2][2] - (-fac + P[i].comp_dx[2] * P[i].comp_dx[2] * fac2);
    tidal_tensorps[1][0]=tidal_tensorps[0][1]; tidal_tensorps[2][0]=tidal_tensorps[0][2]; tidal_tensorps[2][1]=tidal_tensorps[1][2]; /* symmetric so just set these now */

#ifdef BH_OUTPUT_MOREINFO
    printf("Corrected center of mass acceleration %g %g %g tidal tensor diagonal elements %g %g %g \n", P[i].COM_GravAccel[0], P[i].COM_GravAccel[1], P[i].COM_GravAccel[2], tidal_tensorps[0][0],tidal_tensorps[1][1],tidal_tensorps[2][2]);
#endif
    P[i].COM_dt_tidal = 0; for(i1=0;i1<3;i1++) for(i2=0;i2<3;i2++) {P[i].COM_dt_tidal += tidal_tensorps[i1][i2]*tidal_tensorps[i1][i2];}
    P[i].COM_dt_tidal = sqrt(1.0 / (All.G * sqrt(P[i].COM_dt_tidal)));
}
#endif

#ifdef ADAPTIVE_TREEFORCE_UPDATE
int needs_new_treeforce(int n){
    if(P[n].Type > 0){ // in this implementation we only do the lazy updating for gas cells whose timesteps are otherwise constrained by multiphysics (e.g. radiation, feedback)
        return 1;
    } else {
        if(P[n].time_since_last_treeforce >= P[n].tdyn_step_for_treeforce * ADAPTIVE_TREEFORCE_UPDATE) {return 1;}
#ifdef SINGLE_STAR_FB_TIMESTEPLIMIT
        else if(P[n].time_since_last_treeforce >= P[n].min_bh_fb_time) {return 1;} // we want ejecta to re-calculate their feedback time so they don't get stuck on a short timestep
#endif        
        else {return 0;}
    }
}
#endif
