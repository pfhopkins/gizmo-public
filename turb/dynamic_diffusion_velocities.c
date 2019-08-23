#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif



/*! \file dynamic_diffusion_velocities.c
 *  \brief need filtered velocity information to calculate filtered gradients
 *
 */
/*
 * This file was rewritten by Doug Rennehan (douglas.rennehan@gmail.com) for GIZMO, and was
 * copied with modifications from gradients.c, which was written by Phil Hopkins 
 * (phopkins@caltech.edu) for GIZMO.
 */

#ifdef TURB_DIFF_DYNAMIC

#define ASSIGN_ADD_PRESET(x,y,mode) (x+=y)
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#define SHOULD_I_USE_SPH_GRADIENTS(condition_number) ((condition_number > CONDITION_NUMBER_DANGER) ? (1):(0))

#ifdef OMP_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

#define NV_MYSIGN(x) (( x > 0 ) - ( x < 0 ))

struct kernel_DiffFilter {
    double dp[3], r, wk_i, wk_j, dwk_i, dwk_j, h_i;
};

struct DiffFilterdata_in {
    MyDouble Pos[3];
    MyFloat Mass;
    MyFloat Hsml;
    MyDouble Density;
    integertime Timestep;
#ifndef DONOTUSENODELIST
    int NodeList[NODELISTLENGTH];
#endif
#ifdef GALSF_SUBGRID_WINDS
    MyFloat DelayTime; 
#endif
    MyDouble VelPred[3];
}
*DiffFilterDataIn, *DiffFilterDataGet;

struct DiffFilterdata_out {
    MyDouble Norm_hat;
    MyDouble Velocity_bar[3];
    MyFloat FilterWidth_bar;
    MyFloat MaxDistance_for_grad;
}
*DiffFilterDataResult, *DiffFilterDataOut;

/* These functions will handle setting the calculated information */
static inline void particle2in_DiffFilter(struct DiffFilterdata_in *in, int i);
static inline void out2particle_DiffFilter(struct DiffFilterdata_out *out, int i, int mode);

static inline void particle2in_DiffFilter(struct DiffFilterdata_in *in, int i) {
    int k, v;
    for (k = 0; k < 3; k++) {
        in->Pos[k] = P[i].Pos[k];
        in->VelPred[k] = SphP[i].VelPred[k];
    }

    in->Density = SphP[i].Density;
    in->Hsml = PPP[i].Hsml;
    in->Mass = P[i].Mass;

#ifdef GALSF_SUBGRID_WINDS
    in->DelayTime = SphP[i].DelayTime;
#endif

    if (in->Mass < 0) {
        in->Mass = 0;
    }

    in->Timestep = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0);
}


#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))

static inline void out2particle_DiffFilter(struct DiffFilterdata_out *out, int i, int mode) {
    int k, v;
    for (k = 0; k < 3; k++) {
        ASSIGN_ADD_PRESET(SphP[i].Velocity_bar[k], out->Velocity_bar[k], mode);
    }

    MAX_ADD(SphP[i].FilterWidth_bar, out->FilterWidth_bar, mode);
    MAX_ADD(SphP[i].MaxDistance_for_grad, out->MaxDistance_for_grad, mode);
    ASSIGN_ADD_PRESET(SphP[i].Norm_hat, out->Norm_hat, mode);
}


/**
 *
 *
 *
 */
void dynamic_diff_vel_calc(void) {
    mpi_printf("start velocity smoothing computation...\n");
    int i, j, k, v, k1, ngrp, ndone, ndone_flag;
    double shear_factor, prev_coeff_inv;
    double smoothInv = 1.0 / All.TurbDynamicDiffSmoothing;
    int recvTask, place;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0, timewait3 = 0;
    double timecomp, timecomm, timewait, tstart, tend, t0, t1;
    int save_NextParticle;
    long long n_exported = 0;
    
    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct DiffFilterdata_in) +
                                                             sizeof(struct DiffFilterdata_out) +
                                                             sizemax(sizeof(struct DiffFilterdata_in),
                                                                     sizeof(struct DiffFilterdata_out))));
    CPU_Step[CPU_IMPROVDIFFMISC] += measure_time();
    t0 = my_second();
    
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    /* Because of the smoothing operation, need to set bar quantity to current SPH value first */
    for (i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if (P[i].Type == 0) {
            SphP[i].Norm_hat = 0;
            SphP[i].h_turb = Get_Particle_Size(i); // All.cf_atime unnecessary, will multiply later
            SphP[i].FilterWidth_bar = 0;
            SphP[i].MaxDistance_for_grad = 0;

            for (k = 0; k < 3; k++) {
                SphP[i].Velocity_bar[k] = SphP[i].VelPred[k] * smoothInv;
            }
        }
    }

    NextParticle = FirstActiveParticle;	/* begin with this index */
    
    do {    
        BufferFullFlag = 0;
        Nexport = 0;
        save_NextParticle = NextParticle;
            
        for (j = 0; j < NTask; j++) {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
            
        /* do local particles and prepare export list */
        tstart = my_second();
            
#ifdef OMP_NUM_THREADS
        pthread_t mythreads[OMP_NUM_THREADS - 1];
        int threadid[OMP_NUM_THREADS - 1];
        pthread_attr_t attr;
            
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        pthread_mutex_init(&mutex_nexport, NULL);
        pthread_mutex_init(&mutex_partnodedrift, NULL);
            
        TimerFlag = 0;
            
        for (j = 0; j < OMP_NUM_THREADS - 1; j++) {
            threadid[j] = j + 1;
            pthread_create(&mythreads[j], &attr, DiffFilter_evaluate_primary, &threadid[j]);
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
            DiffFilter_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
        }
            
#ifdef OMP_NUM_THREADS
        for (j = 0; j < OMP_NUM_THREADS - 1; j++) pthread_join(mythreads[j], NULL);
#endif
            
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
            
        if (BufferFullFlag) {
            int last_nextparticle = NextParticle;
                
            NextParticle = save_NextParticle;
                
            while (NextParticle >= 0) {
                if (NextParticle == last_nextparticle) break;
                if (ProcessedFlag[NextParticle] != 1) break;
                    
                ProcessedFlag[NextParticle] = 2;
                NextParticle = NextActiveParticle[NextParticle];
            }
                
            if (NextParticle == save_NextParticle) {
                /* in this case, the buffer is too small to process even a single particle */
                endrun(113308);
            }
                
            int new_export = 0;
                
            for (j = 0, k = 0; j < Nexport; j++) {
                if (ProcessedFlag[DataIndexTable[j].Index] != 2) {
                    if (k < j + 1) k = j + 1;
                        
                    for (; k < Nexport; k++) {
                        if (ProcessedFlag[DataIndexTable[k].Index] == 2) {
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
                }
                else {
                    new_export++;
                }                
            }

            Nexport = new_export;       
        }
            
        n_exported += Nexport;
            
        for (j = 0; j < NTask; j++) Send_count[j] = 0;
        for (j = 0; j < Nexport; j++) Send_count[DataIndexTable[j].Task]++;
            
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
            
        tstart = my_second();
            
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
            
        tend = my_second();
        timewait1 += timediff(tstart, tend);
            
        for (j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++) {
            Nimport += Recv_count[j];
                
            if (j > 0) {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
            
        DiffFilterDataGet = (struct DiffFilterdata_in *) mymalloc("DiffFilterDataGet", Nimport * sizeof(struct DiffFilterdata_in));
        DiffFilterDataIn = (struct DiffFilterdata_in *) mymalloc("DiffFilterDataIn", Nexport * sizeof(struct DiffFilterdata_in));
            
        /* prepare particle data for export */
            
        for (j = 0; j < Nexport; j++) {
            place = DataIndexTable[j].Index;
            particle2in_DiffFilter(&DiffFilterDataIn[j], place);
#ifndef DONOTUSENODELIST
            memcpy(DiffFilterDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif
        }
            
        /* exchange particle data */
        tstart = my_second();
        for (ngrp = 1; ngrp < (1 << PTask); ngrp++) {
            recvTask = ThisTask ^ ngrp;
                
            if (recvTask < NTask) {
                if (Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {
                    /* get the particles */
                    MPI_Sendrecv(&DiffFilterDataIn[Send_offset[recvTask]],
                                Send_count[recvTask] * sizeof(struct DiffFilterdata_in), MPI_BYTE,
                                recvTask, TAG_GRADLOOP_A,
                                &DiffFilterDataGet[Recv_offset[recvTask]],
                                Recv_count[recvTask] * sizeof(struct DiffFilterdata_in), MPI_BYTE,
                                recvTask, TAG_GRADLOOP_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

        tend = my_second();
        timecommsumm1 += timediff(tstart, tend);
            
        myfree(DiffFilterDataIn);

        DiffFilterDataResult = (struct DiffFilterdata_out *) mymalloc("DiffFilterDataResult", Nimport * sizeof(struct DiffFilterdata_out));
        DiffFilterDataOut = (struct DiffFilterdata_out *) mymalloc("DiffFilterDataOut", Nexport * sizeof(struct DiffFilterdata_out));
            
        /* now do the particles that were sent to us */
        tstart = my_second();
        NextJ = 0;
            
#ifdef OMP_NUM_THREADS
        for (j = 0; j < OMP_NUM_THREADS - 1; j++) pthread_create(&mythreads[j], &attr, DiffFilter_evaluate_secondary, &threadid[j]);
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
            DiffFilter_evaluate_secondary(&mainthreadid);
        }
            
#ifdef OMP_NUM_THREADS
        for (j = 0; j < OMP_NUM_THREADS - 1; j++) pthread_join(mythreads[j], NULL);
            
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
            
        tend = my_second();
        timecomp2 += timediff(tstart, tend);
            
        if (NextParticle < 0) {
            ndone_flag = 1;
        }
        else {
            ndone_flag = 0;
        }
            
        tstart = my_second();
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        tend = my_second();
        timewait2 += timediff(tstart, tend);
            
        /* get the result */
        tstart = my_second();
        for (ngrp = 1; ngrp < (1 << PTask); ngrp++) {
            recvTask = ThisTask ^ ngrp;
            if (recvTask < NTask) {
                if (Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) {
                    /* send the results */
                    MPI_Sendrecv(&DiffFilterDataResult[Recv_offset[recvTask]],
                                Recv_count[recvTask] * sizeof(struct DiffFilterdata_out),
                                MPI_BYTE, recvTask, TAG_GRADLOOP_B,
                                &DiffFilterDataOut[Send_offset[recvTask]],
                                Send_count[recvTask] * sizeof(struct DiffFilterdata_out),
                                MPI_BYTE, recvTask, TAG_GRADLOOP_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

        tend = my_second();
        timecommsumm2 += timediff(tstart, tend);
            
        /* add the result to the local particles */
        tstart = my_second();
        for (j = 0; j < Nexport; j++) {
            place = DataIndexTable[j].Index;
            out2particle_DiffFilter(&DiffFilterDataOut[j], place, 1);
        }

        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        myfree(DiffFilterDataOut);
        myfree(DiffFilterDataResult);

        myfree(DiffFilterDataGet);
    }
    while(ndone < NTask);
 
    tstart = my_second();

    /* Must wait for ALL tasks for finish each iteration in order to converge */
    /* MPI_Barrier(MPI_COMM_WORLD); */

    tend = my_second();
    timewait3 += timediff(tstart, tend); 
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2 + timewait3;
    timecomm = timecommsumm1 + timecommsumm2;
    
    CPU_Step[CPU_IMPROVDIFFCOMPUTE] += timecomp;
    CPU_Step[CPU_IMPROVDIFFWAIT] += timewait;
    CPU_Step[CPU_IMPROVDIFFCOMM] += timecomm;
    CPU_Step[CPU_IMPROVDIFFMISC] += timeall - (timecomp + timewait + timecomm);
    mpi_printf("velocity smoothing done.\n");
}


/**
 *
 *  Computes the smoothed velocity field Velocity_bar according to eq. 2.17 in 
 *  Monaghan 2011 (turbulence for SPH). 
 *  - D. Rennehan
 *
 */
int DiffFilter_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist) {
    int startnode, numngb, listindex = 0;
    int j, k, v, k2, n, swap_to_j;
    double hinv, hinv3, hinv4, r2, u, hinv_j, hinv3_j, hinv4_j;
    double shear_factor;
    struct kernel_DiffFilter kernel;
    struct DiffFilterdata_in local;
    struct DiffFilterdata_out out;
    memset(&out, 0, sizeof(struct DiffFilterdata_out));

    memset(&kernel, 0, sizeof(struct kernel_DiffFilter));
    
    if (mode == 0) {
        particle2in_DiffFilter(&local, target);
    }
    else {
        local = DiffFilterDataGet[target];
    }
  
    /* check if we should bother doing a neighbor loop */
    if (local.Hsml <= 0) return 0;
    if (local.Mass == 0) return 0;
    if (local.Density <= 0) return 0;
    
    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.Hsml;
    double h2_i = kernel.h_i * kernel.h_i;
    int kernel_mode_i = 0;
    int kernel_mode_j = 0;

    kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
    
    /* Now start the actual neighbor computation for this particle */ 
    if (mode == 0) {
        startnode = All.MaxPart;	/* root node */
    }
    else {
        startnode = DiffFilterDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
   
    while (startnode >= 0) {
        while (startnode >= 0) {
          
            numngb = ngb_treefind_pairs_threads(local.Pos, All.TurbDynamicDiffFac * kernel.h_i, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);

            if (numngb < 0) return -1;
            
            for (n = 0; n < numngb; n++) {
                j = ngblist[n];
                integertime TimeStep_J = (P[j].TimeBin ? (((integertime) 1) << P[j].TimeBin) : 0);
#ifndef SHEARING_BOX // (shearing box means the fluxes at the boundaries are not actually symmetric, so can't do this) //
                if (local.Timestep > TimeStep_J) continue; /* compute from particle with smaller timestep */
                /* use relative positions to break degeneracy */
                if (local.Timestep == TimeStep_J) {
                    int n0 = 0; 
                    if(local.Pos[n0] == P[j].Pos[n0]) {n0++; if(local.Pos[n0] == P[j].Pos[n0]) n0++;}
                    if(local.Pos[n0] < P[j].Pos[n0]) continue;
                }

                swap_to_j = TimeBinActive[P[j].TimeBin];
#else
                swap_to_j = 0;
#endif
#ifdef GALSF_SUBGRID_WINDS
                if (local.DelayTime == 0 && SphP[j].DelayTime > 0) continue;
                if (local.DelayTime > 0 && SphP[j].DelayTime == 0) continue;
#endif
                if (P[j].Mass <= 0) continue;
                if (SphP[j].Density <= 0) continue;
                
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC			/*  now find the closest image in the given box size  */
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                double mean_weight = 0.5 * (local.Density + SphP[j].Density) / (local.Density * SphP[j].Density);
                double h_j = PPP[j].Hsml;
                double V_j = P[j].Mass * mean_weight;
                double V_i = local.Mass * mean_weight;
                if (r2 <= 0) continue;

                double h_avg = 0.5 * (kernel.h_i + h_j);
                double hhat = All.TurbDynamicDiffFac * kernel.h_i;
                double hhat_j = All.TurbDynamicDiffFac * h_j;
                double hhat_avg = All.TurbDynamicDiffFac * h_avg;
                double hhat2 = hhat * hhat;
                double hhatj2 = hhat_j * hhat_j;

                if ((r2 >= hhat2) && (r2 >= hhatj2)) continue;

                double hhatinv, hhatinv3, hhatinv4, uhat, wkhat, dwkhat, rhat;
                double hhatinv_j, hhatinv3_j, hhatinv4_j, wkhat_j, dwkhat_j;

                rhat = sqrt(r2);

                if (rhat < hhat) {
                    kernel_hinv(hhat, &hhatinv, &hhatinv3, &hhatinv4);
                    uhat = DMIN(rhat * hhatinv, 1.0);
                    kernel_main(uhat, hhatinv3, hhatinv4, &wkhat, &dwkhat, kernel_mode_i);
                }
                else {
                    wkhat = dwkhat = 0;
                }

                if (rhat < hhat_j && (swap_to_j)) {
                    kernel_hinv(hhat_j, &hhatinv_j, &hhatinv3_j, &hhatinv4_j);
                    uhat = DMIN(rhat * hhatinv_j, 1.0);
                    kernel_main(uhat, hhatinv3_j, hhatinv4_j, &wkhat_j, &dwkhat_j, kernel_mode_j);
                }
                else {
                    wkhat_j = dwkhat_j = 0;
                }

                if (rhat < hhat) {
                    out.Norm_hat += P[j].Mass * wkhat;
                }

                if (rhat < hhat_j && (swap_to_j)) {
                    SphP[j].Norm_hat += local.Mass * wkhat_j;
                }

                if ((r2 >= h2_i) && (r2 >= (h_j * h_j))) continue;

                kernel.r = sqrt(r2);

                if (kernel.r > out.FilterWidth_bar) {
                    out.FilterWidth_bar = kernel.r;
                }

                // This is very, very important for supersonic flows,
                // or any flow with highly varying smoothing lengths.
                // The FilterWidth_bar (h_bar) *must* extend out to the
                // maximum interaction distance. 
                if (kernel.r > SphP[j].FilterWidth_bar && (swap_to_j)) {
                    SphP[j].FilterWidth_bar = kernel.r;
                }

                if (kernel.r > out.MaxDistance_for_grad) {
                    out.MaxDistance_for_grad = kernel.r;
                }

                if (kernel.r > SphP[j].MaxDistance_for_grad && (swap_to_j)) {
                    SphP[j].MaxDistance_for_grad = kernel.r;
                }

                kernel_hinv(h_avg, &hinv, &hinv3, &hinv4);
                u = DMIN(kernel.r * hinv, 1.0);
                kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, kernel_mode_i);

                double weight_i = kernel.wk_i * V_j;
                double weight_j = kernel.wk_i * V_i;
                double VelPred_diff[3];

                /* Because we are using the average h value between i,j: W_ij = W_ji */
                if (kernel.r < h_avg) {
                    for (k = 0; k < 3; k++) {
                        VelPred_diff[k] = SphP[j].VelPred[k] - local.VelPred[k];
                        out.Velocity_bar[k] += VelPred_diff[k] * weight_i;
                    }

                    if (swap_to_j) {
                        for (k = 0; k < 3; k++) {
                            SphP[j].Velocity_bar[k] -= VelPred_diff[k] * weight_j;
                        }
                    }
                }
            } // numngb loop
        } // while(startnode)
        
#ifndef DONOTUSENODELIST
        if (mode == 1) {
            listindex++;
            if (listindex < NODELISTLENGTH) {
                startnode = DiffFilterDataGet[target].NodeList[listindex];
                if (startnode >= 0) startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
#endif
    }
    
    /* ------------------------------------------------------------------------------------------------ */
    /* Now collect the result at the right place */
    if (mode == 0) {
        out2particle_DiffFilter(&out, target, 0);
    }
    else {
        DiffFilterDataResult[target] = out;
    }
    /* ------------------------------------------------------------------------------------------------ */
    
    return 0;
}


void *DiffFilter_evaluate_primary(void *p) {
    int thread_id = *(int *) p;
    int i, j;
    int *exportflag, *exportnodecount, *exportindex, *ngblist;
    ngblist = Ngblist + thread_id * NumPart;
    exportflag = Exportflag + thread_id * NTask;
    exportnodecount = Exportnodecount + thread_id * NTask;
    exportindex = Exportindex + thread_id * NTask;
    
    /* Note: exportflag is local to each thread */
    for (j = 0; j < NTask; j++) exportflag[j] = -1;
    
    while (1) {
        int exitFlag = 0;
        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
        {
            if (BufferFullFlag != 0 || NextParticle < 0) {
                exitFlag = 1;
            }
            else {
                i = NextParticle;
                ProcessedFlag[i] = 0;
                NextParticle = NextActiveParticle[NextParticle];
            }
        }

        UNLOCK_NEXPORT;
        if (exitFlag) break;
        
        if (P[i].Type == 0) {
	    if (DiffFilter_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0) break;		/* export buffer has filled up */
        }

        ProcessedFlag[i] = 1; /* particle successfully finished */
    }

    return NULL;
}



void *DiffFilter_evaluate_secondary(void *p) {
    int thread_id = *(int *) p;
    int j, dummy, *ngblist;
    ngblist = Ngblist + thread_id * NumPart;

    while (1) {
        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
        {
            j = NextJ;
            NextJ++;
        }

        UNLOCK_NEXPORT;
        
        if (j >= Nimport) break;

        DiffFilter_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

    return NULL;
}

#endif /* End TURB_DIFF_DYNAMIC */

