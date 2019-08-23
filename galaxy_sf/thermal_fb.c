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



/* Routines for pure thermal/scalar feedback/enrichment models: these are intended to represent
    extremely simplified models for stellar feedback manifest as a "pure thermal energy dump"
    (potentially with some cooling turnoff). This is -not- a model for mechanical feedback
    (which, critically, must include the actual momentum and solve for wind/SNe shock properties
    at the interface with the ISM). Those physics are included in the mechanical_fb.c file and
    algorithms therein. This also uses an extremely simple kernel-weighting, rather than a
    more self-consistent area weighting. */

/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#if defined(GALSF_FB_THERMAL)

/* routine that evaluates whether a FB event occurs in a given particle, in a given timestep */
void determine_where_addthermalFB_events_occur(void)
{
    int i; double check = 0;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 4) continue;
        if(P[i].Mass <= 0) continue;
        check += mechanical_fb_calculate_eventrates(i,1); // this should do the calculation and add to number of SNe as needed //
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //
}



/* define structures to use below */
struct addthermalFBdata_in
{
    MyDouble Pos[3], Hsml, Msne, Esne, wt_sum;
#ifdef METALS
    MyDouble yields[NUM_METAL_SPECIES];
#endif
    int NodeList[NODELISTLENGTH];
}
*addthermalFBDataIn, *addthermalFBDataGet;

           
/* define properties to be injected. these must be scalar-only -- the simple routine below will not conserve vector inputs/ejecta (e.g. momentum) */
void particle2in_addthermalFB(struct addthermalFBdata_in *in, int i);
void particle2in_addthermalFB(struct addthermalFBdata_in *in, int i)
{
    if((P[i].SNe_ThisTimeStep<=0)||(P[i].DensAroundStar<=0)) {in->Msne=0; return;} // trap for no sne
    int k; in->Hsml=PPP[i].Hsml; in->wt_sum=P[i].DensAroundStar; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];} // simple kernel-weighted deposition
    struct addFBdata_in local; particle2in_addFB_fromstars(&local,i,0); // get feedback properties from generic routine //
    in->Msne = local.Msne; in->Esne = 0.5 * local.Msne * local.SNe_v_ejecta*local.SNe_v_ejecta; // assign mass and energy to be used below
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=local.yields[k];} // assign yields //
#endif
}



struct addthermalFBdata_out
{
    MyFloat M_coupled;
}
*addthermalFBDataResult, *addthermalFBDataOut;

void out2particle_addthermalFB(struct addthermalFBdata_out *out, int i, int mode);
void out2particle_addthermalFB(struct addthermalFBdata_out *out, int i, int mode)
{
    P[i].Mass -= out->M_coupled;
    if((P[i].Mass<0)||(isnan(P[i].Mass))) {P[i].Mass=0;}
}

struct kernel_addthermalFB {double dp[3], r, wk, dwk, hinv, hinv3, hinv4;};



int addthermalFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double u,r2,h2,kernel_zero,wk;
    struct kernel_addthermalFB kernel;
    struct addthermalFBdata_in local;
    struct addthermalFBdata_out out;
    memset(&out, 0, sizeof(struct addthermalFBdata_out));
    kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1);
    
    /* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_addthermalFB(&local, target);} else {local = addthermalFBDataGet[target];}
    if(local.Msne<=0) return 0; // no SNe for the master particle! nothing to do here //
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml;
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    
    /* Now start the actual FB computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = addthermalFBDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) return -1;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
#endif
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) continue; // same particle //
                if(r2>=h2) continue; // outside kernel //
                // calculate kernel quantities //
                kernel.r = sqrt(r2);
                if(kernel.r <= 0) continue;
                u = kernel.r * kernel.hinv;
                kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                if((kernel.wk <= 0)||(isnan(kernel.wk))) continue;
                wk = P[j].Mass * kernel.wk / local.wt_sum; // normalized weight function
                
                /* inject mass */
                double dM_ejecta_in = wk * local.Msne;
                if(P[j].Hsml<=0) {if(SphP[j].Density>0){SphP[j].Density*=(1+dM_ejecta_in/P[j].Mass);} else {SphP[j].Density=dM_ejecta_in*kernel.hinv3;}} else {SphP[j].Density+=kernel_zero*dM_ejecta_in/(P[j].Hsml*P[j].Hsml*P[j].Hsml);}
                SphP[j].Density *= 1 + dM_ejecta_in/P[j].Mass; // inject mass at constant particle volume //
                P[j].Mass += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[j].MassTrue += dM_ejecta_in;
#endif
#ifdef METALS
                /* inject metals */
                for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k]=(1-dM_ejecta_in/P[j].Mass)*P[j].Metallicity[k] + dM_ejecta_in/P[j].Mass*local.yields[k];}
#endif
                /* inject energy */
                SphP[j].InternalEnergy += wk * local.Esne / P[j].Mass;
                SphP[j].InternalEnergyPred += wk * local.Esne / P[j].Mass;
#ifdef GALSF_FB_TURNOFF_COOLING
                /* if the sub-grid 'cooling turnoff' model is enabled, turn off cooling for the 'blastwave timescale',
                 which is physically the timescale for the blastwave to be completely stopped by ISM ram-pressure
                 (much longer than the actual cooling time of the blastwave) */
                double Esne51 = local.Esne * (All.UnitEnergy_in_cgs/All.HubbleParam) / 1.e51;
                double density_to_n = All.cf_a3inv*All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / PROTONMASS;
                double pressure_to_p4 = (1/All.cf_afac1)*density_to_n*(All.UnitEnergy_in_cgs/All.UnitMass_in_g*PROTONMASS/BOLTZMANN) / 1.0e4;
                double dt_ram = 7.08 * pow(Esne51*SphP[j].Density*density_to_n,0.34) * pow(SphP[j].Pressure*pressure_to_p4,-0.70) / (All.UnitTime_in_Megayears/All.HubbleParam);
                if(dt_ram > SphP[j].DelayTimeCoolingSNe) SphP[j].DelayTimeCoolingSNe = dt_ram;
#endif
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
{
                startnode = addthermalFBDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}    /* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addthermalFB(&out, target, 0);} else {addthermalFBDataResult[target] = out;}
    return 0;
} // int addthermalFB_evaluate



void thermal_fb_calc(void)
{
    int j, k, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
    long long n_exported = 0;
    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct addthermalFBdata_in) + sizeof(struct addthermalFBdata_out) + sizemax(sizeof(struct addthermalFBdata_in),sizeof(struct addthermalFBdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    NextParticle = FirstActiveParticle;	/* begin with this index */
    do
    {
        BufferFullFlag = 0;
        Nexport = 0;
        save_NextParticle = NextParticle;
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local particles and prepare export list */
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
            pthread_create(&mythreads[j], &attr, addthermalFB_evaluate_primary, &threadid[j]);
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
            addthermalFB_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
        }
        
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) pthread_join(mythreads[j], NULL);
#endif
        if(BufferFullFlag)
        {
            int last_nextparticle = NextParticle;
            NextParticle = save_NextParticle;
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle) break;
                if(ProcessedFlag[NextParticle] != 1) break;
                ProcessedFlag[NextParticle] = 2;
                NextParticle = NextActiveParticle[NextParticle];
            }
            if(NextParticle == save_NextParticle) {endrun(116608);} /* in this case, the buffer is too small to process even a single particle */
            int new_export = 0;
            for(j = 0, k = 0; j < Nexport; j++)
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1) k = j + 1;
                    
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
                else {new_export++;}
            Nexport = new_export;
        }
        n_exported += Nexport;
        for(j = 0; j < NTask; j++) Send_count[j] = 0;
        for(j = 0; j < Nexport; j++) Send_count[DataIndexTable[j].Task]++;
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            Nimport += Recv_count[j];
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        /* prepare particle data for export */
        addthermalFBDataGet = (struct addthermalFBdata_in *) mymalloc("addthermalFBDataGet", Nimport * sizeof(struct addthermalFBdata_in));
        addthermalFBDataIn = (struct addthermalFBdata_in *) mymalloc("addthermalFBDataIn", Nexport * sizeof(struct addthermalFBdata_in));
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_addthermalFB(&addthermalFBDataIn[j], place);
            memcpy(addthermalFBDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        /* exchange particle data */
        int TAG_TO_USE = TAG_FBLOOP_1A;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&addthermalFBDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct addthermalFBdata_in), MPI_BYTE, recvTask, TAG_TO_USE,
                                 &addthermalFBDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct addthermalFBdata_in), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        myfree(addthermalFBDataIn);
        addthermalFBDataResult = (struct addthermalFBdata_out *) mymalloc("addthermalFBDataResult", Nimport * sizeof(struct addthermalFBdata_out));
        addthermalFBDataOut = (struct addthermalFBdata_out *) mymalloc("addthermalFBDataOut", Nexport * sizeof(struct addthermalFBdata_out));
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, addthermalFB_evaluate_secondary, &threadid[j]);
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
            addthermalFB_evaluate_secondary(&mainthreadid);
        }
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) pthread_join(mythreads[j], NULL);
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /* get the result */
        TAG_TO_USE = TAG_FBLOOP_1B;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&addthermalFBDataResult[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct addthermalFBdata_out), MPI_BYTE, recvTask, TAG_TO_USE,
                                 &addthermalFBDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct addthermalFBdata_out), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        /* add the result to the local particles */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_addthermalFB(&addthermalFBDataOut[j], place, 1);
        }
        myfree(addthermalFBDataOut);
        myfree(addthermalFBDataResult);
        myfree(addthermalFBDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
}



int addthermalFB_evaluate_active_check(int i);
int addthermalFB_evaluate_active_check(int i)
{
    if(P[i].Type != 4) return 0;
    if(P[i].Mass <= 0) return 0;
    if(PPP[i].Hsml <= 0) return 0;
    if(PPP[i].NumNgb <= 0) return 0;
    if(P[i].SNe_ThisTimeStep>0) {return 1;}
    return 0;
}


void *addthermalFB_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(addthermalFB_evaluate_active_check(i)==1)
#define EVALUATION_CALL addthermalFB_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *addthermalFB_evaluate_secondary(void *p)
{
#define EVALUATION_CALL addthermalFB_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}






#endif /* GALSF_FB_THERMAL */

