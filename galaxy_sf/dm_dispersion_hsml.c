#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

/*! \file dm_dispersion_hsml
 *  \brief smoothing length and velocity dispersion calculation for dark matter particles around gas particles
 *
 *  This file contains a loop modeled on the gas density computation which
 *    determines kernel lengths for dark matter particles around a given set of gas particles;
 *    this is used by the flag GALSF_SUBGRID_WIND_SCALING==2 to estimate the local dark
 *    matter velocity dispersion around a given gas particle, which (in turn) is used to set the
 *    sub-grid wind velocity and mass loading. The loop here needs to be called for these models (note this
 *    in general will require a different smoothing length from, say, the dm-dm force softening, or the
 *    gas kernel length for hydro, hence it requires a whole additional loop, even though the loop is
 *    functionally identical - modulo which particles are used - to the loop for adaptive gravitational softening
 *
 * This file was written by Qirong Zhu, for GIZMO, based on Phil Hopkins's adaptive gravitational softening
 *    routine. It has been modified by Phil on re-merger into the main branch of GIZMO with various optimizations
 *
 */

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct disp_densdata_in
{
    MyDouble Pos[3];
    MyFloat HsmlDM;
    int NodeList[NODELISTLENGTH];
}
*DISP_DensDataIn, *DISP_DensDataGet;

static struct disp_densdata_out
{
    MyLongDouble Ngb, DM_Vel_Disp, DM_Vx, DM_Vy, DM_Vz;
}
*DISP_DensDataResult, *DISP_DensDataOut;

void disp_particle2in_density(struct disp_densdata_in *in, int i);
void disp_out2particle_density(struct disp_densdata_out *out, int i, int mode);

void disp_particle2in_density(struct disp_densdata_in *in, int i)
{
    int k;
    for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    in->HsmlDM = SphP[i].HsmlDM;
}

void disp_out2particle_density(struct disp_densdata_out *out, int i, int mode)
{
    ASSIGN_ADD(SphP[i].DM_Vx, out->DM_Vx, mode);
    ASSIGN_ADD(SphP[i].DM_Vy, out->DM_Vy, mode);
    ASSIGN_ADD(SphP[i].DM_Vz, out->DM_Vz, mode);
    ASSIGN_ADD(SphP[i].DM_VelDisp, out->DM_Vel_Disp, mode);
    ASSIGN_ADD(SphP[i].NumNgbDM, out->Ngb, mode);
}

void disp_density(void)
{
    MyFloat *Left, *Right;
    int i, j, k, ndone, ndone_flag, npleft, iter = 0;
    int ngrp, recvTask, place;
    long long ntot;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
    double timecomp, timecomm, timewait;
    double tstart, tend, t0, t1;
    double desnumngb, desnumngbdev;
    int save_NextParticle;
    long long n_exported = 0;
    int redo_particle;
    
    CPU_Step[CPU_AGSDENSMISC] += measure_time();
    
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(disp_density_isactive(i))
        {
            SphP[i].NumNgbDM = 0;
            Left[i] = Right[i] = 0;
        }
    }
    
    /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                           sizeof(struct disp_densdata_in) + sizeof(struct disp_densdata_out) +
                                                           sizemax(sizeof(struct disp_densdata_in),sizeof(struct disp_densdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    t0 = my_second();
    
    desnumngb = 64;
    desnumngbdev = 48;
    
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        
        NextParticle = FirstActiveParticle;	/* begin with this index */
        
        do
        {
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
                pthread_create(&mythreads[j], &attr, disp_density_evaluate_primary, &threadid[j]);
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
                disp_density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
            }
            
#ifdef PTHREADS_NUM_THREADS
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                pthread_join(mythreads[j], NULL);
#endif
            
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
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
                    printf("DM disp: Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n",ThisTask,P[NextParticle].Type,
                           P[NextParticle].Pos[0],P[NextParticle].Pos[1],P[NextParticle].Pos[2],P[NextParticle].Mass);
                    
                    endrun(111008);
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
            
            DISP_DensDataGet = (struct disp_densdata_in *) mymalloc("DISP_DensDataGet", Nimport * sizeof(struct disp_densdata_in));
            DISP_DensDataIn = (struct disp_densdata_in *) mymalloc("DISP_DensDataIn", Nexport * sizeof(struct disp_densdata_in));
            
            /* prepare particle data for export */
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index;
                
                disp_particle2in_density(&DISP_DensDataIn[j], place);
                
                memcpy(DISP_DensDataIn[j].NodeList,
                       DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
                        MPI_Sendrecv(&DISP_DensDataIn[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct disp_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DMDENS_A,
                                     &DISP_DensDataGet[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct disp_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DMDENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second();
            timecommsumm1 += timediff(tstart, tend);
            
            myfree(DISP_DensDataIn);
            DISP_DensDataResult = (struct disp_densdata_out *) mymalloc("DISP_DensDataResult", Nimport * sizeof(struct disp_densdata_out));
            DISP_DensDataOut = (struct disp_densdata_out *) mymalloc("DISP_DensDataOut", Nexport * sizeof(struct disp_densdata_out));
            
            /* now do the particles that were sent to us */
            
            tstart = my_second();
            
            NextJ = 0;
            
#ifdef PTHREADS_NUM_THREADS
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                pthread_create(&mythreads[j], &attr, disp_density_evaluate_secondary, &threadid[j]);
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
                disp_density_evaluate_secondary(&mainthreadid);
            }
            
#ifdef PTHREADS_NUM_THREADS
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                pthread_join(mythreads[j], NULL);
            
            pthread_mutex_destroy(&mutex_partnodedrift);
            pthread_mutex_destroy(&mutex_nexport);
            pthread_attr_destroy(&attr);
#endif
            
            tend = my_second();
            timecomp2 += timediff(tstart, tend);
            
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
                        MPI_Sendrecv(&DISP_DensDataResult[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct disp_densdata_out),
                                     MPI_BYTE, recvTask, TAG_DMDENS_B,
                                     &DISP_DensDataOut[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct disp_densdata_out),
                                     MPI_BYTE, recvTask, TAG_DMDENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                
            }
            tend = my_second();
            timecommsumm2 += timediff(tstart, tend);
            
            
            /* add the result to the local particles */
            tstart = my_second();
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index;
                disp_out2particle_density(&DISP_DensDataOut[j], place, 1);
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
            
            myfree(DISP_DensDataOut);
            myfree(DISP_DensDataResult);
            myfree(DISP_DensDataGet);
        }
        while(ndone < NTask);
        
        
        /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(disp_density_isactive(i))
            {
                redo_particle = 0;
                /* now check whether we have enough neighbours, and are below the maximum search radius */
                double maxsoft = DMIN(All.MaxHsml, 10.0*PPP[i].Hsml);
                if(((SphP[i].NumNgbDM < desnumngb - desnumngbdev) || (SphP[i].NumNgbDM > (desnumngb + desnumngbdev)))
                   && (Right[i]-Left[i] > 0.001*Left[i] || Left[i]==0 || Right[i]==0))
                {
                    redo_particle = 1;
                }
                if(SphP[i].HsmlDM >= maxsoft)
                {
                    SphP[i].HsmlDM = maxsoft;
                    redo_particle = 0;
                }

                if(redo_particle)
                {
                    /* need to redo this particle */
                    npleft++;
                    
                    if(SphP[i].NumNgbDM < desnumngb-desnumngbdev) {Left[i]=DMAX(SphP[i].HsmlDM, Left[i]);}
                    if(SphP[i].NumNgbDM > desnumngb+desnumngbdev) {if(Right[i]>0) {Right[i]=DMIN(Right[i],SphP[i].HsmlDM);} else {Right[i]=SphP[i].HsmlDM;}}
                    
                    if(iter >= MAXITER - 10)
                    {
                        printf
                        ("DM disp: i=%d task=%d ID=%llu Type=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                         i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, SphP[i].HsmlDM, Left[i], Right[i],
                         (float) SphP[i].NumNgbDM, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                    }
                    
                    // right/left define upper/lower bounds from previous iterations
                    if(Right[i] > 0 && Left[i] > 0)
                    {
                        // geometric interpolation between right/left //
                        if(SphP[i].NumNgbDM > 1)
                        {
                            SphP[i].HsmlDM *= pow( desnumngb / SphP[i].NumNgbDM , 1./NUMDIMS );
                        } else {
                            SphP[i].HsmlDM *= 2.0;
                        }
                        if((SphP[i].HsmlDM<Right[i])&&(SphP[i].HsmlDM>Left[i]))
                        {
                            SphP[i].HsmlDM = pow(SphP[i].HsmlDM*SphP[i].HsmlDM*SphP[i].HsmlDM*SphP[i].HsmlDM * Left[i]*Right[i] , 1./6.);
                        } else {
                            if(SphP[i].HsmlDM>Right[i]) SphP[i].HsmlDM=Right[i];
                            if(SphP[i].HsmlDM<Left[i]) SphP[i].HsmlDM=Left[i];
                            SphP[i].HsmlDM = pow(SphP[i].HsmlDM * Left[i] * Right[i] , 1.0/3.0);
                        }
                    }
                    else
                    {
                        if(Right[i] == 0 && Left[i] == 0)
                        {
                            char buf[1000];
                            sprintf(buf, "DM disp: Right[i] == 0 && Left[i] == 0 && SphP[i].HsmlDM=%g\n", SphP[i].HsmlDM);
                            terminate(buf);
                        }
                        double fac;
                        if(Right[i] == 0 && Left[i] > 0)
                        {
                            if(SphP[i].NumNgbDM > 1) {fac = log( desnumngb / SphP[i].NumNgbDM ) / NUMDIMS;} else {fac=1.4;}
                            if((SphP[i].NumNgbDM < 2*desnumngb)&&(SphP[i].NumNgbDM > 0.1*desnumngb)) {SphP[i].HsmlDM *= exp(fac);} else {SphP[i].HsmlDM *= 1.26;}
                        }
                        
                        if(Right[i] > 0 && Left[i] == 0)
                        {
                            if(SphP[i].NumNgbDM > 1) {fac = log( desnumngb / SphP[i].NumNgbDM ) / NUMDIMS;} else {fac=1.4;}
                            fac = DMAX(fac,-1.535);
                            if((SphP[i].NumNgbDM < 2*desnumngb)&&(SphP[i].NumNgbDM > 0.1*desnumngb)) {SphP[i].HsmlDM *= exp(fac);} else {SphP[i].HsmlDM /= 1.26;}
                        }
                    }
                }
                else
                    P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
            } //  if(disp_density_isactive(i))
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        
        if(ntot > 0)
        {
            iter++;
            if(iter > 0 && ThisTask == 0)
            {
#ifdef IO_REDUCED_MODE
                if(iter > 10)
#endif
                printf("DM disp: ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
                       (int) (ntot / 1000000000), (int) (ntot % 1000000000));
            }
            if(iter > MAXITER)
            {
                printf("DM disp: failed to converge in neighbour iteration in disp_density()\n");
                fflush(stdout);
                endrun(1155);
            }
        }
    }
    while(ntot > 0);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Right);
    myfree(Left);
    myfree(Ngblist);
    
    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0)
            P[i].TimeBin = -P[i].TimeBin - 1;
    }
    
    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(disp_density_isactive(i))
        {
            if(SphP[i].NumNgbDM > 0)
            {
                SphP[i].DM_Vx /= SphP[i].NumNgbDM;
                SphP[i].DM_Vy /= SphP[i].NumNgbDM;
                SphP[i].DM_Vz /= SphP[i].NumNgbDM;
                SphP[i].DM_VelDisp /= SphP[i].NumNgbDM;
                SphP[i].DM_VelDisp = (1./All.cf_atime) * sqrt(SphP[i].DM_VelDisp -
                                                              SphP[i].DM_Vx * SphP[i].DM_Vx -
                                                              SphP[i].DM_Vy * SphP[i].DM_Vy -
                                                              SphP[i].DM_Vz * SphP[i].DM_Vz) / 1.732;//	   1d velocity dispersion
            } else {
                if((SphP[i].DM_VelDisp <= 0) || isnan(SphP[i].DM_VelDisp))
                {
                    SphP[i].DM_VelDisp = sqrt(P[i].Vel[0]*P[i].Vel[0]+P[i].Vel[1]*P[i].Vel[1]+P[i].Vel[2]*P[i].Vel[2])/All.cf_atime;
                }
            }
        }
    }
    
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp;
    CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm;
    CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}






/*! This function represents the core of the density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int disp_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, n;
    int startnode, numngb_inbox, listindex = 0;
    struct disp_densdata_in local;
    struct disp_densdata_out out;
    memset(&out, 0, sizeof(struct disp_densdata_out));
    
    if(mode == 0)
        disp_particle2in_density(&local, target);
    else
        local = DISP_DensDataGet[target];

    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = DISP_DensDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.HsmlDM, target, &startnode, mode, exportflag,
                                                                 exportnodecount, exportindex, ngblist, 2); // search for high-res DM particles only: 2^1 = 2
            if(numngb_inbox < 0) return -1;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Mass <= 0) continue;
                out.DM_Vx += P[j].Vel[0];
                out.DM_Vy += P[j].Vel[1];
                out.DM_Vz += P[j].Vel[2];
                out.DM_Vel_Disp += (P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);
                out.Ngb++;
            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DISP_DensDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }
    
    if(mode == 0)
        disp_out2particle_density(&out, target, 0);
    else
        DISP_DensDataResult[target] = out;
    
    return 0;
}



void *disp_density_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(disp_density_isactive(i))
#define EVALUATION_CALL disp_density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *disp_density_evaluate_secondary(void *p)
{
#define EVALUATION_CALL disp_density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}

#endif
#endif


/* routine to determine if we need to use disp_density to calculate Hsml */
int disp_density_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0;
    if(P[i].Type > 0) return 0; // only gas particles //
    if(P[i].Mass <= 0) return 0;
    return 1;
}


