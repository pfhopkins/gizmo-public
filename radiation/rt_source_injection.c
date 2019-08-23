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

/*! \file rt_source_injection.c
 *  \brief inject luminosity from point sources to neighboring gas particles
 *
 *  This file contains a loop modeled on the gas density computation which will
 *    share luminosity from non-gas particles to the surrounding gas particles,
 *    so that it can be treated within e.g. the flux-limited diffusion or other
 *    radiation-hydrodynamic approximations. Basically the same concept as
 *    injecting the radiation 'in a cell' surrounding the particle
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#if defined(RT_RAD_PRESSURE_FORCES)
#define RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION // adds correction for un-resolved extinction which cannot generate photon momentum with M1, FLD, OTVET, etc.
#endif

#if defined(GALSF) && !defined(RT_INJECT_PHOTONS_DISCRETELY)
#define RT_INJECT_PHOTONS_DISCRETELY
#endif

#ifdef RT_SOURCE_INJECTION

/*! Structure for communication during the kernel computation. Holds data that is sent to other processors  */
static struct rt_sourcedata_in
{
    MyDouble Pos[3];
    MyFloat Hsml;
    MyFloat KernelSum_Around_RT_Source;
    MyFloat Luminosity[N_RT_FREQ_BINS];
    int NodeList[NODELISTLENGTH];
}
*RT_SourceDataIn, *RT_SourceDataGet;

/* declare subroutines */
void rt_particle2in_source(struct rt_sourcedata_in *in, int i);
int rt_sourceinjection_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void *rt_sourceinjection_evaluate_primary(void *p);
void *rt_sourceinjection_evaluate_secondary(void *p);


/*! subroutine to insert the data needed to be passed to other processors: here for convenience, match to structure above  */
void rt_particle2in_source(struct rt_sourcedata_in *in, int i)
{
    int k;
    for(k=0; k<3; k++) {in->Pos[k] = P[i].Pos[k];}
    in->Hsml = PPP[i].Hsml;
    //if(P[i].Type==0) {in->KernelSum_Around_RT_Source = SphP[i].Density;} else {in->KernelSum_Around_RT_Source = P[i].DensAroundStar;}
    in->KernelSum_Around_RT_Source = P[i].KernelSum_Around_RT_Source;
    /* luminosity is set to zero here for gas particles because their self-illumination is handled trivially in a single loop, earlier */
    double lum[N_RT_FREQ_BINS];
    int active_check = rt_get_source_luminosity(i,0,lum);
    double dt = 1; // make this do nothing unless flags below are set:
#if defined(RT_INJECT_PHOTONS_DISCRETELY)
#ifndef WAKEUP
    dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
#else
    dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a;
#endif
#endif
    for(k=0; k<N_RT_FREQ_BINS; k++)
    {
        if(P[i].Type==0 || active_check==0) {in->Luminosity[k]=0;} else {in->Luminosity[k] = lum[k] * dt;}
    }
}

/* routine to do the master loop over particles, for the source injection (photons put into surrounding gas) */
void rt_source_injection(void)
{
    int j, k, ngrp, ndone, ndone_flag;
    int recvTask, place;
    int save_NextParticle;
    long long n_exported = 0;
    
    /* first, we do a loop over the gas particles themselves. these are trivial -- they don't need to share any information,
     they just determine their own source functions. so we don't need to do any loops. and we can zero everything before the loop below. */
    for(j=0;j<N_gas;j++)
    {
        if(P[j].Type==0)
        {
            double lum[N_RT_FREQ_BINS];
            for(k=0;k<N_RT_FREQ_BINS;k++) {SphP[j].Je[k]=0;} // need to zero -before- calling injection //
            int active_check = rt_get_source_luminosity(j,0,lum);
            /* here is where we would need to code some source luminosity for the gas */
            for(k=0;k<N_RT_FREQ_BINS;k++) if(active_check) {SphP[j].Je[k]=lum[k];}
        }
    }

    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct rt_sourcedata_in) + sizeof(struct rt_sourcedata_in)));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    NextParticle = FirstActiveParticle;	/* begin with this index */
    do
    {
        BufferFullFlag = 0;
        Nexport = 0;
        save_NextParticle = NextParticle;
        for(j = 0; j < NTask; j++) {Send_count[j] = 0; Exportflag[j] = -1;}
        
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
            pthread_create(&mythreads[j], &attr, rt_sourceinjection_evaluate_primary, &threadid[j]);
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
            rt_sourceinjection_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
        }
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
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
            if(NextParticle == save_NextParticle)
            {
                /* in this case, the buffer is too small to process even a single particle */
                endrun(116610);
            }
            int new_export = 0;
            for(j = 0, k = 0; j < Nexport; j++)
            {
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1) {k = j + 1;}
                    for(; k < Nexport; k++)
                    {
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
                } else {new_export++;}
            }
            Nexport = new_export;
        }
        n_exported += Nexport;
        
        for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
        for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
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
        
        RT_SourceDataGet = (struct rt_sourcedata_in *) mymalloc("RT_SourceDataGet", Nimport * sizeof(struct rt_sourcedata_in));
        RT_SourceDataIn = (struct rt_sourcedata_in *) mymalloc("RT_SourceDataIn", Nexport * sizeof(struct rt_sourcedata_in));
        /* prepare particle data for export */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            rt_particle2in_source(&RT_SourceDataIn[j], place);
#ifndef DONOTUSENODELIST
            memcpy(RT_SourceDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif
        }
        /* exchange particle data */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&RT_SourceDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct rt_sourcedata_in), MPI_BYTE, recvTask, TAG_RT_C,
                                 &RT_SourceDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct rt_sourcedata_in), MPI_BYTE, recvTask, TAG_RT_C, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        myfree(RT_SourceDataIn);
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, rt_sourceinjection_evaluate_secondary, &threadid[j]);
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
            rt_sourceinjection_evaluate_secondary(&mainthreadid);
        }
        
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        myfree(RT_SourceDataGet);
    }
    while(ndone < NTask);
    /* free memory */
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
}


/* subroutine that actually distributes the luminosity as desired to neighbor particles in the kernel */
int rt_sourceinjection_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    /* Load the data for the particle */
    int j, k, n, startnode, numngb_inbox, listindex = 0;
    struct rt_sourcedata_in local;
    if(mode == 0) {rt_particle2in_source(&local, target);} else {local = RT_SourceDataGet[target];}
    /* basic calculations */
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    double hinv, hinv3, hinv4, h2=local.Hsml*local.Hsml;
    kernel_hinv(local.Hsml, &hinv, &hinv3, &hinv4);
    
    /* Now start the actual operations for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = RT_SourceDataGet[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;/* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -1;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                double dp[3]; for(k=0; k<3; k++) {dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef BOX_PERIODIC	/* find the closest image in the given box size  */
                NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif
                double r2=0; for(k=0;k<3;k++) {r2 += dp[k]*dp[k];}
                if(r2<=0) continue; // same particle //
                if(r2>=h2) continue; // outside kernel //
                // calculate kernel quantities //
                //double wk, dwk, u = sqrt(r2) * hinv;
                //kernel_main(u, hinv3, hinv4, &wk, &dwk, -1); // traditional kernel
                //wk *= P[j].Mass / local.KernelSum_Around_RT_Source;
                double wk = (1 - r2*hinv*hinv) / local.KernelSum_Around_RT_Source;
#if defined(RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION)
                double r = sqrt(r2);
                double dv0 = -1. / (RT_SPEEDOFLIGHT_REDUCTION * (C / All.UnitVelocity_in_cm_per_s) * r);
                double lmax_0 = DMAX(local.Hsml, r);
#ifdef RT_EVOLVE_INTENSITIES
                int kx; double angle_wt_Inu_sum=0, angle_wt_Inu[N_RT_INTENSITY_BINS];
                // pre-compute a set of weights based on the projection of the particle position along the radial direction for the radiation direction //
                for(kx=0;kx<N_RT_INTENSITY_BINS;kx++)
                {
                    double cos_t=0; int kq; for(kq=0;kq<3;kq++) {cos_t+=dp[kq]*All.RT_Intensity_Direction[kx][kq];}
                    cos_t *= -1/r;
                    double wt_function = cos_t*cos_t*cos_t*cos_t;
                    if(cos_t < 0) {wt_function=0;}
                    angle_wt_Inu[kx] = wt_function; angle_wt_Inu_sum += angle_wt_Inu[kx];
                }
#endif
#endif
                // now actually apply the kernel distribution
                for(k=0;k<N_RT_FREQ_BINS;k++) 
                {
                    double dE = wk * local.Luminosity[k];
#if defined(RT_INJECT_PHOTONS_DISCRETELY)
                    SphP[j].E_gamma[k] += dE;
#ifdef RT_EVOLVE_NGAMMA
                    SphP[j].E_gamma_Pred[k] += dE; // dump discreetly (noisier, but works smoothly with large timebin hierarchy)
#endif
#if defined(RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION)
                    // add discrete photon momentum from un-resolved absorption //
                    double x_abs = 2. * SphP[j].Kappa_RT[k] * (SphP[j].Density*All.cf_a3inv) * (DMAX(2.*Get_Particle_Size(j),lmax_0)*All.cf_atime); // effective optical depth through particle
                    double slabfac_x = x_abs * slab_averaging_function(x_abs); // 1-exp(-x)
                    if(isnan(slabfac_x)||(slabfac_x<=0)) {slabfac_x=0;}
                    if(slabfac_x>1) {slabfac_x=1;}
                    double dv = slabfac_x * dv0 * dE / P[j].Mass; // total absorbed momentum (needs multiplication by dp[kv] for directionality)
                    int kv; for(kv=0;kv<3;kv++) {P[j].Vel[kv] += dv*dp[kv]; SphP[j].VelPred[kv] += dv*dp[kv];}
#if defined(RT_EVOLVE_FLUX)
                    double dflux = -dE * (RT_SPEEDOFLIGHT_REDUCTION * (C / All.UnitVelocity_in_cm_per_s)) / r;
                    for(kv=0;kv<3;kv++) {SphP[j].Flux[k][kv] += dflux*dp[kv]; SphP[j].Flux_Pred[k][kv] += dflux*dp[kv];}
#endif
#ifdef RT_EVOLVE_INTENSITIES
                    double dflux = dE * (RT_SPEEDOFLIGHT_REDUCTION * (C / All.UnitVelocity_in_cm_per_s)) / angle_wt_Inu_sum;
                    for(kv=0;kv<N_RT_INTENSITY_BINS;kv++) {SphP[j].Intensity[k][kv] += dflux * angle_wt_Inu[N_RT_INTENSITY_BINS]; SphP[j].Intensity_Pred[k][kv] += dflux * angle_wt_Inu[N_RT_INTENSITY_BINS];}
#endif
#endif
#else
                    SphP[j].Je[k] += dE; // treat continuously
#endif
                }
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
#ifndef DONOTUSENODELIST
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = RT_SourceDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        } // if(mode == 1)
#endif
    } // while(startnode >= 0)
    return 0;
} // int rt_sourceinjection_evaluate


int rt_sourceinjection_active_check(int i);
int rt_sourceinjection_active_check(int i)
{
    if(PPP[i].NumNgb <= 0) return 0;
    if(PPP[i].Hsml <= 0) return 0;
    if(PPP[i].Mass <= 0) return 0;
    double lum[N_RT_FREQ_BINS];
    return rt_get_source_luminosity(i,-1,lum);
}

/* routine for initial loop of particles on local processor (and determination of which need passing) */
void *rt_sourceinjection_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(rt_sourceinjection_active_check(i)==1)
#define EVALUATION_CALL rt_sourceinjection_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *rt_sourceinjection_evaluate_secondary(void *p)
{
#define EVALUATION_CALL rt_sourceinjection_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}




#endif
