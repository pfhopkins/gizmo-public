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

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#ifdef PTHREADS_NUM_THREADS
pthread_mutex_t mutex_nexport;
pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif



double Ewaldcount, Costtotal;
long long N_nodesinlist;


int Ewald_iter;			/* global in file scope, for simplicity */



void sum_top_level_node_costfactors(void);



/*! This function computes the gravitational forces for all active particles.
 *  If needed, a new tree is constructed, otherwise the dynamically updated
 *  tree is used.  Particles are only exported to other processors when really
 *  needed, thereby allowing a good use of the communication buffer.
 */
void gravity_tree(void)
{
    long long n_exported = 0; int i, j, maxnumnodes, iter; iter = 0;
    double t0, t1, timeall = 0, timetree1 = 0, timetree2 = 0, timetree, timewait, timecomm;
    double timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0, sum_costtotal, ewaldtot;
    double maxt, sumt, maxt1, sumt1, maxt2, sumt2, sumcommall, sumwaitall, plb, plb_max;
    
#ifdef FIXEDTIMEINFIRSTPHASE
    int counter; double min_time_first_phase, min_time_first_phase_glob;
#endif
#ifndef SELFGRAVITY_OFF
    int k, ewald_max, diff, save_NextParticle, ndone, ndone_flag, ngrp, place, recvTask; double tstart, tend, ax, ay, az; MPI_Status status;
#endif

    CPU_Step[CPU_MISC] += measure_time();
    
    /* set new softening lengths */
    if(All.ComovingIntegrationOn)
        set_softenings();
    
    /* construct tree if needed */
    if(TreeReconstructFlag)
    {
#ifndef IO_REDUCED_MODE
        if(ThisTask == 0) printf("Tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
#endif
        CPU_Step[CPU_MISC] += measure_time();

        move_particles(All.Ti_Current);

        force_treebuild(NumPart, NULL);
        
        CPU_Step[CPU_TREEBUILD] += measure_time();

        TreeReconstructFlag = 0;
        
#ifndef IO_REDUCED_MODE
	if(ThisTask == 0) printf("Tree construction done.\n");
#endif
    }
    
#ifndef SELFGRAVITY_OFF
    
    /* allocate buffers to arrange communication */
#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    if(ThisTask == 0) printf("Begin tree force.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

    
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct gravdata_in) + sizeof(struct gravdata_out) +
                                             sizemax(sizeof(struct gravdata_in),sizeof(struct gravdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    if(ThisTask == 0) printf("All.BunchSize=%d\n", All.BunchSize);
    
    Ewaldcount = 0;
    Costtotal = 0;
    N_nodesinlist = 0;
    
    
    CPU_Step[CPU_TREEMISC] += measure_time();
    t0 = my_second();
    
#if defined(BOX_PERIODIC) && !defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC)
    ewald_max = 1;
#else
    ewald_max = 0;
#endif
    

        if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)
        {
            /* we have a fresh tree and would like to measure gravity cost */
            
            /* find the closest level */
            for(i = 1, TakeLevel = 0, diff = abs(All.LevelToTimeBin[0] - All.HighestActiveTimeBin);
                i < GRAVCOSTLEVELS; i++)
            {
                if(diff > abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin))
                {
                    TakeLevel = i;
                    diff = abs(All.LevelToTimeBin[i] - All.HighestActiveTimeBin);
                }
            }
            
            if(diff != 0)		/* we have not found a matching slot */
            {
                
                if(All.HighestOccupiedTimeBin - All.HighestActiveTimeBin < GRAVCOSTLEVELS)	/* we should have space */
                {
                    /* clear levels that are out of range */
                    for(i = 0; i < GRAVCOSTLEVELS; i++)
                    {
                        if(All.LevelToTimeBin[i] > All.HighestOccupiedTimeBin)
                            All.LevelToTimeBin[i] = 0;
                        if(All.LevelToTimeBin[i] < All.HighestOccupiedTimeBin - (GRAVCOSTLEVELS - 1))
                            All.LevelToTimeBin[i] = 0;
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
                    terminate("TakeLevel < 0, even though we should have a slot");
            }
        }
        else
        {
            /* in this case we do not measure gravity cost. Check whether this time-level
             has previously mean measured. If yes, then delete it so to make sure that it is not out of time */
            
            for(i = 0; i < GRAVCOSTLEVELS; i++)
                if(All.LevelToTimeBin[i] == All.HighestActiveTimeBin)
                    All.LevelToTimeBin[i] = 0;
            
            TakeLevel = -1;
        }
        
        if(TakeLevel >= 0)
            for(i = 0; i < NumPart; i++)
                P[i].GravCost[TakeLevel] = 0;
        
        
        for(Ewald_iter = 0; Ewald_iter <= ewald_max; Ewald_iter++)
        {
            
            NextParticle = FirstActiveParticle;	/* begin with this index */
            
            do
            {
                iter++;
                BufferFullFlag = 0;
                Nexport = 0;
                save_NextParticle = NextParticle;
                
                tstart = my_second();
                
#ifdef PTHREADS_NUM_THREADS
                pthread_t mythreads[PTHREADS_NUM_THREADS - 1];
                int threadid[PTHREADS_NUM_THREADS - 1];
                pthread_attr_t attr;
                
                pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                pthread_mutex_init(&mutex_nexport, NULL);
                pthread_mutex_init(&mutex_partnodedrift, NULL);
                
                TimerFlag = 0;
                
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
                for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                    pthread_join(mythreads[j], NULL);
#endif
                
                tend = my_second();
                timetree1 += timediff(tstart, tend);
                
                
                if(BufferFullFlag)
                {
                    int last_nextparticle = NextParticle;
                    
                    NextParticle = save_NextParticle;
                    
                    while(NextParticle >= 0)
                    {
                        if(NextParticle == last_nextparticle)
                            break;
                        
                        if(ProcessedFlag[NextParticle] != 1)
                            break;
                        
                        ProcessedFlag[NextParticle] = 2;
                        
                        NextParticle = NextActiveParticle[NextParticle];
                    }
                    
                    if(NextParticle == save_NextParticle)
                    {
                        /* in this case, the buffer is too small to process even a single particle */
                        endrun(114408);
                    }
                    
                    
                    int new_export = 0;
                    
                    for(j = 0, k = 0; j < Nexport; j++)
                        if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                        {
                            if(k < j + 1)
                                k = j + 1;
                            
                            for(; k < Nexport; k++)
                                if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                                {
                                    int old_index = DataIndexTable[j].Index;
                                    
                                    DataIndexTable[j] = DataIndexTable[k];
                                    DataNodeList[j] = DataNodeList[k];
                                    DataIndexTable[j].IndexGet = j;
                                    new_export++;
                                    
                                    DataIndexTable[k].Index = old_index;
                                    k++;
                                    break;
                                }
                        }
                        else
                            new_export++;
                    
                    Nexport = new_export;
                    
                }
                
                
                n_exported += Nexport;
                
                for(j = 0; j < NTask; j++)
                    Send_count[j] = 0;
                for(j = 0; j < Nexport; j++)
                    Send_count[DataIndexTable[j].Task]++;
                
                
                MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
                
                tstart = my_second();
                
                MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
                
                tend = my_second();
                timewait1 += timediff(tstart, tend);
                
                
                for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
                {
                    Nimport += Recv_count[j];
                    
                    if(j > 0)
                    {
                        Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                        Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                    }
                }
                
                GravDataGet = (struct gravdata_in *) mymalloc("GravDataGet", Nimport * sizeof(struct gravdata_in));
                GravDataIn = (struct gravdata_in *) mymalloc("GravDataIn", Nexport * sizeof(struct gravdata_in));
                
                /* prepare particle data for export */
                
                for(j = 0; j < Nexport; j++)
                {
                    place = DataIndexTable[j].Index;
                    
                    GravDataIn[j].Type = P[place].Type;
                    GravDataIn[j].OldAcc = P[place].OldAcc;
                    for(k = 0; k < 3; k++) {GravDataIn[j].Pos[k] = P[place].Pos[k];}
#if defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
                    GravDataIn[j].Mass = P[place].Mass;
#endif
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(COMPUTE_JERK_IN_GRAVTREE)
                    for(k = 0; k < 3; k++) {GravDataIn[j].Vel[k] = P[place].Vel[k];}		    
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
                    if (P[place].Type==5)
                    {
                        GravDataIn[j].min_bh_t_orbital = P[place].min_bh_t_orbital; //orbital time for binary
                        GravDataIn[j].comp_Mass = P[place].comp_Mass; //mass of binary companion
                        GravDataIn[j].is_in_a_binary = P[place].is_in_a_binary; // 1 if we're in a binary, 0 if not
                        for(k=0;k<3;k++) {GravDataIn[j].comp_dx[k]=P[place].comp_dx[k]; GravDataIn[j].comp_dv[k]=P[place].comp_dv[k];}
                    }
                    else {P[place].is_in_a_binary=0; /* setting values to zero just to be sure */}
#endif
#if defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
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
                
                
                /* exchange particle data */
                
                tstart = my_second();
                for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
                {
                    recvTask = ThisTask ^ ngrp;
                    
                    if(recvTask < NTask)
                    {
                        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                        {
                            /* get the particles */
                            MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]],
                                         Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
                                         recvTask, TAG_GRAV_A,
                                         &GravDataGet[Recv_offset[recvTask]],
                                         Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
                                         recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
                        }
                    }
                }
                tend = my_second();
                timecommsumm1 += timediff(tstart, tend);
                
                
                myfree(GravDataIn);
                GravDataResult =
                (struct gravdata_out *) mymalloc("GravDataResult", Nimport * sizeof(struct gravdata_out));
                GravDataOut =
                (struct gravdata_out *) mymalloc("GravDataOut", Nexport * sizeof(struct gravdata_out));
                
                report_memory_usage(&HighMark_gravtree, "GRAVTREE");

                /* now do the particles that were sent to us */
                tstart = my_second();
                
                NextJ = 0;
                
#ifdef PTHREADS_NUM_THREADS
                for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                    pthread_create(&mythreads[j], &attr, gravity_secondary_loop, &threadid[j]);
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
                for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                    pthread_join(mythreads[j], NULL);
                
                pthread_mutex_destroy(&mutex_partnodedrift);
                pthread_mutex_destroy(&mutex_nexport);
                pthread_attr_destroy(&attr);
#endif
                
                tend = my_second();
                timetree2 += timediff(tstart, tend);
                
                if(NextParticle < 0)
                    ndone_flag = 1;
                else
                    ndone_flag = 0;
                
                tstart = my_second();
                MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                tend = my_second();
                timewait2 += timediff(tstart, tend);
                
                
                /* get the result */
                tstart = my_second();
                for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
                {
                    recvTask = ThisTask ^ ngrp;
                    if(recvTask < NTask)
                    {
                        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                        {
                            /* send the results */
                            MPI_Sendrecv(&GravDataResult[Recv_offset[recvTask]],
                                         Recv_count[recvTask] * sizeof(struct gravdata_out),
                                         MPI_BYTE, recvTask, TAG_GRAV_B,
                                         &GravDataOut[Send_offset[recvTask]],
                                         Send_count[recvTask] * sizeof(struct gravdata_out),
                                         MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);
                        }
                    }
                    
                }
                tend = my_second();
                timecommsumm2 += timediff(tstart, tend);
                
                
                /* add the results to the local particles */
                tstart = my_second();
                for(j = 0; j < Nexport; j++)
                {
                    place = DataIndexTable[j].Index;
                    for(k = 0; k < 3; k++) {P[place].GravAccel[k] += GravDataOut[j].Acc[k];}

#ifdef BH_CALC_DISTANCES
                    /* GravDataOut[j].min_dist_to_bh contains the min dist to particle "P[place]" on another
                     task.  We now check if it is smaller than the current value */
                    if(Ewald_iter==0)
                    {
                        if(GravDataOut[j].min_dist_to_bh < P[place].min_dist_to_bh)
                        {
                            P[place].min_dist_to_bh = GravDataOut[j].min_dist_to_bh;
                            P[place].min_xyz_to_bh[0] = GravDataOut[j].min_xyz_to_bh[0];
                            P[place].min_xyz_to_bh[1] = GravDataOut[j].min_xyz_to_bh[1];
                            P[place].min_xyz_to_bh[2] = GravDataOut[j].min_xyz_to_bh[2];
                        }
                    }
#endif
                    if(Ewald_iter==0) /* don't allow for an infinite hierarchy of these moments, or you will get nonsense */
                    {
#ifdef RT_OTVET
                        if(P[place].Type==0) {int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) for(k=0;k<6;k++) SphP[place].ET[k_freq][k] += GravDataOut[j].ET[k_freq][k];}
#endif
                    }
                    
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
                    if(Ewald_iter==0) {int i1tt,i2tt; for(i1tt=0;i1tt<3;i1tt++) {for(i2tt=0;i2tt<3;i2tt++) {P[place].tidal_tensorps[i1tt][i2tt] += GravDataOut[j].tidal_tensorps[i1tt][i2tt];}}}
#ifdef COMPUTE_JERK_IN_GRAVTREE
                    if(Ewald_iter==0) {int i1tt; for(i1tt=0; i1tt<3; i1tt++) P[place].GravJerk[i1tt] += GravDataOut[j].GravJerk[i1tt];}
#endif		    		    
#endif
                    
#ifdef EVALPOTENTIAL
                    P[place].Potential += GravDataOut[j].Potential;
#endif
                }
                tend = my_second();
                timetree1 += timediff(tstart, tend);
                
                myfree(GravDataOut);
                myfree(GravDataResult);
                myfree(GravDataGet);
            }
            while(ndone < NTask);
        }			/* Ewald_iter */
        


    myfree(DataNodeList);
    myfree(DataIndexTable);


        
    /* assign node cost to particles */
    if(TakeLevel >= 0)
    {
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
    
    /* now add things for comoving integration */
    
#ifndef BOX_PERIODIC
#ifndef PMGRID
    if(All.ComovingIntegrationOn)
    {
        double fac = 0.5 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits * All.Omega0 / All.G;
        
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        {
            for(j = 0; j < 3; j++)
                P[i].GravAccel[j] += fac * P[i].Pos[j];
        }
    }
#endif
#endif
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef PMGRID
        ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
        ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
        az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
        ax = P[i].GravAccel[0];
        ay = P[i].GravAccel[1];
        az = P[i].GravAccel[2];
#endif
        
        if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS && All.Ti_Current == 0 && RestartFlag == 0)
            continue;		/* to prevent that we overwrite OldAcc in the first evaluation for 2lpt ICs */
        
        P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
    }

#ifndef GRAVITY_HYBRID_OPENING_CRIT  // in collisional systems we don't want to rely on the relative opening criterion alone, because aold can be dominated by a binary companion but we still want accurate contributions from distant nodes. Thus we combine BH and relative criteria. - MYG
    if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
    {
        if(!(All.Ti_Current == 0 && RestartFlag == 0)) {if(All.TypeOfOpeningCriterion == 1) {All.ErrTolTheta = 0;}}	/* This will switch to the relative opening criterion for the following force computations */
    }
    else
    {
        if(All.TypeOfOpeningCriterion == 1) {All.ErrTolTheta = 0;} /* This will switch to the relative opening criterion for the following force computations */
    }
#endif
    
    /*  multiply by G */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
	    		
        for(j=0;j<3;j++) {P[i].GravAccel[j] *= All.G;}

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
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
        gsl_eigen_symm_free(w);
        gsl_vector_free (eval);
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

        
#ifdef EVALPOTENTIAL /* remove self-potential */
	    //	P[i].Potential += P[i].Mass / All.SofteningTable[P[i].Type];
#ifdef BOX_PERIODIC
        if(All.ComovingIntegrationOn) {P[i].Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) * pow(All.Omega0 * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G), 1.0 / 3);}
#endif
        P[i].Potential *= All.G;
#ifdef PMGRID
        P[i].Potential += P[i].PM_Potential;	/* add in long-range potential */
#endif
        
        if(All.ComovingIntegrationOn)
        {
#ifndef BOX_PERIODIC
            double fac, r2; fac = -0.5 * All.Omega0 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits;
            for(k = 0, r2 = 0; k < 3; k++) {r2 += P[i].Pos[k] * P[i].Pos[k];}
            P[i].Potential += fac * r2;
#endif
        } else {
            double fac, r2; fac = -0.5 * All.OmegaLambda * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits;
            if(fac != 0)
            {
                for(k = 0, r2 = 0; k < 3; k++) {r2 += P[i].Pos[k] * P[i].Pos[k];}
                P[i].Potential += fac * r2;
            }
        }
#endif
    }
    
    
#if defined(RT_OTVET)
    /* normalize the Eddington tensors we just calculated by walking the tree */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        if(P[i].Type == 0)
        {
            int k_freq;
            for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++)
            {
                double trace = SphP[i].ET[k_freq][0] + SphP[i].ET[k_freq][1] + SphP[i].ET[k_freq][2];
#ifdef RT_FLUXLIMITEDDIFFUSION
                trace = 0; /* force the code to always use the isotropic Eddington tensor */
#endif
                if(!isnan(trace) && (trace > 0))
                {
                    for(k = 0; k < 6; k++) SphP[i].ET[k_freq][k] /= trace;
                }
                else
                {
                    for(k = 0; k < 6; k++) {SphP[i].ET[k_freq][k] = 0.0;}
                    SphP[i].ET[k_freq][0] = SphP[i].ET[k_freq][1] = SphP[i].ET[k_freq][2] = 1./3.;
                }
            }
        }
#endif
    
    /* Finally, the following factor allows a computation of a cosmological simulation
     with vacuum energy in physical coordinates */
#ifndef BOX_PERIODIC
#ifndef PMGRID
    if(All.ComovingIntegrationOn == 0)
    {
        double fac = All.OmegaLambda * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits;
        
        for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]){
            for(j = 0; j < 3; j++)
                P[i].GravAccel[j] += fac * P[i].Pos[j];
	}
    }
#endif
#endif
    
        
#else /* gravity is switched off */
    t0 = my_second();
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {for(j = 0; j < 3; j++) {P[i].GravAccel[j] = 0;}}
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {P[i].tidal_tensorps[0][0] = P[i].tidal_tensorps[0][1] = P[i].tidal_tensorps[0][2] = P[i].tidal_tensorps[1][0] = P[i].tidal_tensorps[1][1] = P[i].tidal_tensorps[1][2] = P[i].tidal_tensorps[2][0] = P[i].tidal_tensorps[2][1] = P[i].tidal_tensorps[2][2] = 0.0;}
#endif
#endif /* end of SELFGRAVITY_OFF */
    
    
    
    
   
    
#ifdef RT_SELFGRAVITY_OFF
    /* if this is set, we zero out gravity here, just after computing it! */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
	P[i].GravAccel[0]=P[i].GravAccel[1]=P[i].GravAccel[2]=0;
    }
#endif
    
    add_analytic_gravitational_forces();
    
    
    /* Now the force computation is finished */
    
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    
    /*  gather some diagnostic information */
    
    timetree = timetree1 + timetree2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    
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
    CPU_Step[CPU_TREEWALK1] += timetree1;
    CPU_Step[CPU_TREEWALK2] += timetree2;
    CPU_Step[CPU_TREESEND] += timecommsumm1;
    CPU_Step[CPU_TREERECV] += timecommsumm2;
    CPU_Step[CPU_TREEWAIT1] += timewait1;
    CPU_Step[CPU_TREEWAIT2] += timewait2;
    
    
#ifdef FIXEDTIMEINFIRSTPHASE
    MPI_Reduce(&min_time_first_phase, &min_time_first_phase_glob, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        printf("FIXEDTIMEINFIRSTPHASE=%g  min_time_first_phase_glob=%g\n",
               FIXEDTIMEINFIRSTPHASE, min_time_first_phase_glob);
    }
#endif
    
#ifndef IO_REDUCED_MODE
    if(ThisTask == 0)
    {
        fprintf(FdTimings, "Step= %lld  t= %g  dt= %g \n",(long long) All.NumCurrentTiStep, All.Time, All.TimeStep);
        fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g (%g) iter= %d\n",
                (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
                (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
                n_exported / ((double) GlobNumForceUpdate), N_nodesinlist / ((double) n_exported + 1.0e-10),
                iter);
        /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */
        
        fprintf(FdTimings, "work-load balance: %g (%g %g) rel1to2=%g   max=%g avg=%g\n",
                maxt / (1.0e-6 + sumt / NTask), maxt1 / (1.0e-6 + sumt1 / NTask),
                maxt2 / (1.0e-6 + sumt2 / NTask), sumt1 / (1.0e-6 + sumt1 + sumt2), maxt, sumt / NTask);
        fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
        fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
                maxnumnodes / (All.TreeAllocFactor * All.MaxPart + NTopnodes));
        fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", GlobNumForceUpdate / (sumt + 1.0e-20),
                GlobNumForceUpdate / (1.0e-6 + maxt * NTask),
                ((double) (sum_costtotal)) / (1.0e-20 + GlobNumForceUpdate),
                ((double) ewaldtot) / (1.0e-20 + GlobNumForceUpdate));
        fprintf(FdTimings, "\n");
        
        fflush(FdTimings);
    }
#endif
    
    CPU_Step[CPU_TREEMISC] += measure_time();
    
#ifndef IO_REDUCED_MODE
    double costtotal_new = 0, sum_costtotal_new;
    if(TakeLevel >= 0)
    {
        for(i = 0; i < NumPart; i++)
            costtotal_new += P[i].GravCost[TakeLevel];
        MPI_Reduce(&costtotal_new, &sum_costtotal_new, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(ThisTask == 0)
            printf("relative error in the total number of tree-gravity interactions = %g\n",
                   (sum_costtotal - sum_costtotal_new) / sum_costtotal);
        /* can be non-zero if THREAD_SAFE_COSTS is not used (and due to round-off errors). */
    }
#endif
}




void *gravity_primary_loop(void *p)
{
    int i, j, ret;
    int thread_id = *(int *) p;
    
    int *exportflag, *exportnodecount, *exportindex;
    
    exportflag = Exportflag + thread_id * NTask;
    exportnodecount = Exportnodecount + thread_id * NTask;
    exportindex = Exportindex + thread_id * NTask;
    
    /* Note: exportflag is local to each thread */
    for(j = 0; j < NTask; j++)
        exportflag[j] = -1;
    
#ifdef FIXEDTIMEINFIRSTPHASE
    int counter = 0;
    double tstart;
    
    if(thread_id == 0)
    {
        tstart = my_second();
    }
#endif
    
    while(1)
    {
        int exitFlag = 0;
        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
        if(BufferFullFlag != 0 || NextParticle < 0)
        {
            exitFlag = 1;
        }
        else
        {
            i = NextParticle;
            ProcessedFlag[i] = 0;
            NextParticle = NextActiveParticle[NextParticle];
        }
        UNLOCK_NEXPORT;
        if(exitFlag)
            break;
		    
#if !defined(PMGRID)
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
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
#else
        
#ifdef NEUTRINOS
        if(P[i].Type != 2)
#endif
        {
            ret = force_treeevaluate(i, 0, exportflag, exportnodecount, exportindex);
            if(ret < 0) {break;} /* export buffer has filled up */
            Costtotal += ret;
        }
        
#endif
        ProcessedFlag[i] = 1;	/* particle successfully finished */
        
#ifdef FIXEDTIMEINFIRSTPHASE
        if(thread_id == 0)
        {
            counter++;
            if((counter & 255) == 0)
            {
                if(timediff(tstart, my_second()) > FIXEDTIMEINFIRSTPHASE)
                {
                    TimerFlag = 1;
                    break;
                }
            }
        }
        else
        {
            if(TimerFlag)
                break;
        }
#endif
    }
    
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
        
        if(j >= Nimport)
            break;
        
#if !defined(PMGRID)
#if defined(BOX_PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
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
#else
        ret = force_treeevaluate(j, 1, &nodesinlist, &dummy, &dummy);
        N_nodesinlist += nodesinlist; Costtotal += ret;
#endif
    }
    
    return NULL;
}



void sum_top_level_node_costfactors(void)
{
    int i;
    
    double *costlist = (double*)mymalloc("costlist", NTopnodes * sizeof(double));
    double *costlist_all = (double*)mymalloc("costlist_all", NTopnodes * sizeof(double));
    
    for(i = 0; i < NTopnodes; i++)
        costlist[i] = Nodes[All.MaxPart + i].GravCost;
    
    MPI_Allreduce(costlist, costlist_all, NTopnodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for(i = 0; i < NTopnodes; i++)
        Nodes[All.MaxPart + i].GravCost = costlist_all[i];
    
    myfree(costlist_all);
    myfree(costlist);
}





/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
    int i;
    
    if(All.ComovingIntegrationOn)
    {
        if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
            All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
        else
            All.SofteningTable[0] = All.SofteningGas;
        
        if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
            All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
        else
            All.SofteningTable[1] = All.SofteningHalo;
        
        if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
            All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
        else
            All.SofteningTable[2] = All.SofteningDisk;
        
        if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
            All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
        else
            All.SofteningTable[3] = All.SofteningBulge;
        
        if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
            All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
        else
            All.SofteningTable[4] = All.SofteningStars;
        
        if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
            All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
        else
            All.SofteningTable[5] = All.SofteningBndry;
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
    for(i = 0; i < 6; i++)
    {
        All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];
    }
    /* set the minimum gas kernel length to be used this timestep */
    All.MinHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
#ifndef SELFGRAVITY_OFF
    if(All.MinHsml <= 5.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0])
        All.MinHsml = 5.0*EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[0];
#endif
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int data_index_compare(const void *a, const void *b)
{
    if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task))
        return -1;
    
    if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task))
        return +1;
    
    if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index))
        return -1;
    
    if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index))
        return +1;
    
    if(((struct data_index *) a)->IndexGet < (((struct data_index *) b)->IndexGet))
        return -1;
    
    if(((struct data_index *) a)->IndexGet > (((struct data_index *) b)->IndexGet))
        return +1;
    
    return 0;
}

static void msort_dataindex_with_tmp(struct data_index *b, size_t n, struct data_index *t)
{
    struct data_index *tmp;
    struct data_index *b1, *b2;
    size_t n1, n2;
    
    if(n <= 1)
        return;
    
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
    
    if(n1 > 0)
        memcpy(tmp, b1, n1 * sizeof(struct data_index));
    
    memcpy(b, t, (n - n2) * sizeof(struct data_index));
}

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
    const size_t size = n * s;
    
    struct data_index *tmp = (struct data_index *) mymalloc("struct data_index *tmp", size);
    
    msort_dataindex_with_tmp((struct data_index *) b, n, tmp);
    
    myfree(tmp);
}

