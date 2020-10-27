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

/*! \file rt_CGmethod.c
 *  \brief solve the implicit sparse matrix inversion problem for RT fluxes via CG iteraiton
 *
 *  This file contains a loop modeled on the gas density computation which is fully parallel
 *    and calculates the matrix and its inversion, using the conjugate-gradient method (coupled to
 *    the method of steepest descent). The core of the module can be used for any global sparse
 *    matrix inversion. Here it is applied to obtain the implicit solution of the diffusion
 *    problem for radiative transfer, following Petkova & Springel 2008
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO, heavily modifying some of the
 *   original routines from GADGET-3 and re-writing the actual parallelization of the module.
 */

#ifdef RT_DIFFUSION_CG


/* declare some of the global variables and functions that are going to be used in this solution */

/*! global variables to be used */
#define MAX_ITER 10000
#define ACCURACY 1.0e-2
static double **ZVec, **XVec, **QVec, **DVec, **Residue, **Diag, **Diag2;

/*! structure for communication. holds data that is sent to other processors  */
static struct rt_cg_data_in
{
    int NodeList[NODELISTLENGTH];
    MyDouble Pos[3];
    MyFloat Mass;
    MyFloat Density;
    MyFloat Hsml;
    MyFloat ET[N_RT_FREQ_BINS][6];
    MyDouble RT_DiffusionCoeff[N_RT_FREQ_BINS];
    //MyDouble Lambda[N_RT_FREQ_BINS];
}
*rt_cg_DataIn, *rt_cg_DataGet;

/*! structure for communication. holds data to be returned to the spawning processor */
struct rt_cg_data_out
{
    MyDouble matrixmult_out[N_RT_FREQ_BINS];
    MyDouble matrixmult_sum[N_RT_FREQ_BINS];
}
*rt_cg_DataResult, *rt_cg_DataOut;

/*! declare functions */
double rt_diffusion_cg_vector_multiply(double *a, double *b);
double rt_diffusion_cg_vector_sum(double *a);
void rt_diffusion_cg_matrix_multiply(double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum);
int rt_diffusion_cg_evaluate(int target, int mode, double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist);
void particle2in_rt_cg(struct rt_cg_data_in *in, int i);
void *rt_diffusion_cg_evaluate_primary(void *p, double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum);
void *rt_diffusion_cg_evaluate_secondary(void *p, double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum);


/*! subroutine to insert the data needed to be passed to other processors: here for convenience, match to structure above  */
void particle2in_rt_cg(struct rt_cg_data_in *in, int i)
{
    int k;
    for(k=0; k<3; k++) {in->Pos[k] = P[i].Pos[k];}
    int kET; for(k=0;k<N_RT_FREQ_BINS;k++) for(kET=0; kET<6; kET++) {in->ET[k][kET] = SphP[i].ET[k][kET];}
    in->Hsml = PPP[i].Hsml;
    in->Mass = P[i].Mass;
    in->Density = SphP[i].Density;
    for(k=0; k<N_RT_FREQ_BINS; k++) in->RT_DiffusionCoeff[k] = rt_diffusion_coefficient(i,k);
}

/* internal product of two vectors (for all gas particles) */
double rt_diffusion_cg_vector_multiply(double *a, double *b)
{
    int i; double sum, sumall;
    for(i = 0, sum = 0; i < N_gas; i++)
        if(P[i].Type == 0)
            sum += a[i] * b[i];
    MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sumall;
}


/* absolute sum of vector elements (for all gas particles) */
double rt_diffusion_cg_vector_sum(double *a)
{
    int i; double sum, sumall;
    for(i = 0, sum = 0; i < N_gas; i++)
        if(P[i].Type == 0)
            sum += fabs(a[i]);
    MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sumall;
}

/* define a convenient macro for allocating the required arrays below */
#define MALLOC_CG(x) {\
x = (double **) malloc(N_RT_FREQ_BINS * sizeof(double *));\
for(k=0;k<N_RT_FREQ_BINS;k++) x[k] = (double *) malloc(N_gas * sizeof(double));\
for(k=0;k<N_RT_FREQ_BINS;k++) memset(x[k], 0, N_gas * sizeof(double));}


/*! routine to do the top-level loop for the CG iteration - this is the actual solver; it calls various subroutines
 to do the weights/matrix calculation on all particles */
void rt_diffusion_cg_solve(void)
{
    PRINT_STATUS("start CG iteration for radiative transfer (diffusion equation)...");
    int k, j; double alpha_cg, beta, sum, rel, res, maxrel, glob_maxrel, DQ;
    double dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * UNIT_INTEGERTIME_IN_PHYSICAL;
    
    /* initialization for the CG method */
    MALLOC_CG(ZVec); MALLOC_CG(XVec); MALLOC_CG(QVec); MALLOC_CG(DVec); MALLOC_CG(Residue); MALLOC_CG(Diag); MALLOC_CG(Diag2); // allocate and zero all the arrays
    for(j = 0; j < N_gas; j++)
        if(P[j].Type == 0)
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                XVec[k][j] = SphP[j].Rad_E_gamma[k] * SphP[j].Density / (1.e-37+P[j].Mass); /* define the coefficients: note we need energy densities for this operation */
                SphP[j].Rad_E_gamma[k] += dt * SphP[j].Rad_Je[k]; /* -then- add the source terms */
            }
 
    /* do a first pass of our 'workhorse' routine, which lets us pre-condition to improve convergence */
    rt_diffusion_cg_matrix_multiply(XVec, Residue, Diag);
    /* take the diagonal matrix elements as a Jacobi preconditioner */
    double delta_new_initial[N_RT_FREQ_BINS];
    for(k = 0; k < N_RT_FREQ_BINS; k++)
    {
        for(j = 0; j < N_gas; j++)
            if(P[j].Type == 0)
            {                
                Residue[k][j] = SphP[j].Rad_E_gamma[k] * SphP[j].Density / (1.e-37+P[j].Mass) - Residue[k][j]; // note: source terms have been added here to Rad_E_gamma //
                /* note: in principle we would have to substract the w_ii term, but this is zero by definition */
                ZVec[k][j] = Residue[k][j] / Diag[k][j];
                DVec[k][j] = ZVec[k][j];
            }
        delta_new_initial[k] = rt_diffusion_cg_vector_multiply(ZVec[k], Residue[k]);
    }
    
    /* begin the CG method iteration */
    int iter=0, ndone=0, done_key[N_RT_FREQ_BINS]; 
    double delta_new[N_RT_FREQ_BINS], delta_old[N_RT_FREQ_BINS];
    for(k = 0; k < N_RT_FREQ_BINS; k++) {done_key[k]=0; delta_new[k]=delta_new_initial[k]; delta_old[k]=delta_new_initial[k];}
    do
    {
        /* this is the 'workhorse' routine with the neighbor communication and actual calculation */
        rt_diffusion_cg_matrix_multiply(DVec, QVec, Diag2);
        ndone = 0; 
        for(k = 0; k < N_RT_FREQ_BINS; k++)
        {
            /* define residues */
            DQ = rt_diffusion_cg_vector_multiply(DVec[k], QVec[k]);
            if(DQ == 0) {alpha_cg = 0;} else {alpha_cg = delta_new[k] / DQ;}
            for(j = 0, maxrel = 0; j < N_gas; j++)
            {
                XVec[k][j] += alpha_cg * DVec[k][j];
                Residue[k][j] -= alpha_cg * QVec[k][j];
                ZVec[k][j] = Residue[k][j] / Diag[k][j];
                rel = fabs(alpha_cg * DVec[k][j]) / (XVec[k][j] + 1.0e-10);
                if(rel > maxrel) {maxrel = rel;}
            }
            delta_old[k] = delta_new[k];
            delta_new[k] = rt_diffusion_cg_vector_multiply(ZVec[k], Residue[k]);
            
            /* sum up residues to define next step */
            sum = rt_diffusion_cg_vector_sum(XVec[k]);
            res = rt_diffusion_cg_vector_sum(Residue[k]);
            if(delta_old[k]) {beta = delta_new[k] / delta_old[k];} else {beta = 0;}
            for(j = 0; j < N_gas; j++) {DVec[k][j] = ZVec[k][j] + beta * DVec[k][j];}
            
            /* broadcast and decide if we need to keep iterating */
            MPI_Allreduce(&maxrel, &glob_maxrel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            PRINT_STATUS("CG iteration: iter=%3d  |res|/|x|=%12.6g  maxrel=%12.6g  |x|=%12.6g |res|=%12.6g\n", iter, res / sum, glob_maxrel, sum, res);
            if(iter >= 1 && (res <= ACCURACY * sum || iter >= MAX_ITER)) {done_key[k]=1; ndone++;}
        }
        iter++;
        if(iter > MAX_ITER) {terminate("failed to converge in CG iteration");}
    }
    while(ndone < N_RT_FREQ_BINS);
    
    /* success! */
    PRINT_STATUS("%d iterations performed\n", iter);
    /* update the intensity */
    for(j = 0; j < N_gas; j++)
        if(P[j].Type == 0)
            for(k = 0; k < N_RT_FREQ_BINS; k++)
                SphP[j].Rad_E_gamma[k] = DMAX(XVec[k][j],0) * P[j].Mass / SphP[j].Density; // convert back to an absolute energy, instead of a density //
    
    /* free memory */
    free(Diag2);
    free(Diag);
    free(Residue);
    free(DVec);
    free(QVec);
    free(XVec);
    free(ZVec);
}







/* this function computes the vector b(matrixmult_out) given the vector x(in) such as Ax = b, where A is a matrix */
void rt_diffusion_cg_matrix_multiply(double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum)
{
    /* allocate buffers to arrange communication */
    int j, k, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
    long long n_exported = 0;
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                             sizeof(struct rt_cg_data_in) + sizeof(struct rt_cg_data_out) + sizemax(sizeof(struct rt_cg_data_in),sizeof(struct rt_cg_data_out))));
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
            pthread_create(&mythreads[j], &attr, rt_diffusion_cg_evaluate_primary, &threadid[j]);
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
            rt_diffusion_cg_evaluate_primary(&mainthreadid, matrixmult_in, matrixmult_out, matrixmult_sum);	/* do local particles and prepare export list */
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
                endrun(116609);
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
        rt_cg_DataGet = (struct rt_cg_data_in *) mymalloc("rt_cg_DataGet", Nimport * sizeof(struct rt_cg_data_in));
        rt_cg_DataIn = (struct rt_cg_data_in *) mymalloc("rt_cg_DataIn", Nexport * sizeof(struct rt_cg_data_in));
        /* prepare particle data for export */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_rt_cg(&rt_cg_DataIn[j], place);
#ifndef DONOTUSENODELIST
            memcpy(rt_cg_DataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
                    MPI_Sendrecv(&rt_cg_DataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct rt_cg_data_in), MPI_BYTE, recvTask, TAG_RT_A,
                                 &rt_cg_DataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct rt_cg_data_in), MPI_BYTE, recvTask, TAG_RT_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        myfree(rt_cg_DataIn);
        rt_cg_DataResult = (struct rt_cg_data_out *) mymalloc("rt_cg_DataResult", Nimport * sizeof(struct rt_cg_data_out));
        rt_cg_DataOut = (struct rt_cg_data_out *) mymalloc("rt_cg_DataOut", Nexport * sizeof(struct rt_cg_data_out));
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, rt_diffusion_cg_evaluate_secondary, &threadid[j]);
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
            rt_diffusion_cg_evaluate_secondary(&mainthreadid, matrixmult_in, matrixmult_out, matrixmult_sum);
        }
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /* get the result */
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&rt_cg_DataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct rt_cg_data_out), MPI_BYTE, recvTask, TAG_RT_B,
                                 &rt_cg_DataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct rt_cg_data_out), MPI_BYTE, recvTask, TAG_RT_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        
        /* add the result to the local particles */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                matrixmult_out[k][place] += rt_cg_DataOut[j].matrixmult_out[k];
                matrixmult_sum[k][place] += rt_cg_DataOut[j].matrixmult_sum[k];
            }
        }
        myfree(rt_cg_DataOut);
        myfree(rt_cg_DataResult);
        myfree(rt_cg_DataGet);
    }
    while(ndone < NTask);
    /* free memory */
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);

    
    /* do final operations on results */
    {double dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * UNIT_INTEGERTIME_IN_PHYSICAL; int i;
    for(i = 0; i < N_gas; i++)
        if(P[i].Type == 0) {
            for(k = 0; k < N_RT_FREQ_BINS; k++) {
                double fac_i = dt * rt_absorption_rate(i,k); 
                if((1 + fac_i + matrixmult_sum[k][i]) < 0) {printf("1 + matrixmult_sum + rate= %g   matrixmult_sum=%g rate=%g i =%d\n", 1 + fac_i + matrixmult_sum[k][i], matrixmult_sum[k][i], fac_i, i); endrun(11111111);}
                /* the "1" here accounts for the fact that we must start from the previous photon number (the matrix includes only the "dt" term); 
                    the fac_i term here accounts for sinks [here, the rate of photon absorption]; the in*sum part below accounts for the re-arrangement 
                    of indices [swapping indices i and j in the relevant equations so we account for both sides of the difference terms */
                matrixmult_sum[k][i] += 1.0 + fac_i; matrixmult_out[k][i] += matrixmult_in[k][i] * matrixmult_sum[k][i];
            }}}
}


/* subroutine that actually does the neighbor calculation: this needs to be customized for your problem! */
/*   -- this subroutine contains no writes to shared memory -- */
int rt_diffusion_cg_evaluate(int target, int mode, double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    /* Load the data for the particle */
    int j, k, n, startnode, numngb_inbox, listindex = 0;
    struct rt_cg_data_in local;
    if(mode == 0) {particle2in_rt_cg(&local, target);} else {local = rt_cg_DataGet[target];}
    struct rt_cg_data_out out;
    memset(&out, 0, sizeof(struct rt_cg_data_out));

    /* basic calculations */
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    double hinv, hinv3, hinv4, h2=local.Hsml*local.Hsml;
    kernel_hinv(local.Hsml, &hinv, &hinv3, &hinv4);
    double dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * UNIT_INTEGERTIME_IN_PHYSICAL;
#ifdef RT_DIFFUSION_CG_MODIFY_EDDINGTON_TENSOR
    /*modify Eddington tensor */
    for(j=0;j<N_RT_FREQ_BINS;j++)
    {
        double ET[6];
        int kET; for(kET = 0; k < 6; k++) {ET[k] = local.ET[j][k];}
        local.ET[j][0] = 2.*ET[0] - 0.5*ET[1] - 0.5*ET[2];
        local.ET[j][1] = 2.*ET[1] - 0.5*ET[2] - 0.5*ET[0];
        local.ET[j][2] = 2.*ET[2] - 0.5*ET[0] - 0.5*ET[1];
        for(k=3;k<6;k++) {local.ET[j][k] = 2.5*ET[k];}
    }
#endif
    
    /* Now start the actual operations for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = rt_cg_DataGet[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;/* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                double dp[3]; for(k=0; k<3; k++) {dp[k] = local.Pos[k] - P[j].Pos[k];}
                NEAREST_XYZ(dp[0],dp[1],dp[2],1); /* find the closest image in the given box size */
                double r2=0; for(k=0;k<3;k++) {r2 += dp[k]*dp[k];}
                if(r2<=0) continue; // same particle //
                if((r2>h2)||(r2>PPP[j].Hsml*PPP[j].Hsml)) continue; // outside kernel //
                // calculate kernel quantities //
                double r = sqrt(r2), wk, dwk_i=0, dwk_j=0;
                if(r<local.Hsml)
                {
                    kernel_main(r*hinv, hinv3, hinv4, &wk, &dwk_i, 1);
                }
                if(r<PPP[j].Hsml)
                {
                    double hinv_j,hinv3_j,hinv4_j; kernel_hinv(PPP[j].Hsml, &hinv_j, &hinv3_j, &hinv4_j);
                    kernel_main(r*hinv_j, hinv3_j, hinv4_j, &wk, &dwk_j, 1);
                }
                
                double tensor_norm = -dt * (dwk_i*local.Mass/local.Density + dwk_j*P[j].Mass/SphP[j].Density) / r;
                if(tensor_norm > 0)
                {
                    for(k=0;k<N_RT_FREQ_BINS;k++)
                    {
                
                        double ET_ij[6];
#ifdef RT_DIFFUSION_CG_MODIFY_EDDINGTON_TENSOR
                        double ET_j[6];
                        ET_j[0] = 2.*SphP[j].ET[k][0] - 0.5*SphP[j].ET[k][1] - 0.5*SphP[j].ET[k][2];
                        ET_j[1] = 2.*SphP[j].ET[k][1] - 0.5*SphP[j].ET[k][2] - 0.5*SphP[j].ET[k][0];
                        ET_j[2] = 2.*SphP[j].ET[k][2] - 0.5*SphP[j].ET[k][0] - 0.5*SphP[j].ET[k][1];
                        int kET;
                        for(kET=3;kET<6;kET++) {ET_j[kET] = 2.5*SphP[j].ET[k][kET];}
                        for(kET=0;kET<6;kET++) {ET_ij[kET] = 0.5 * (local.ET[k][kET] + ET_j[kET]);}
#else
                        int kET; for(kET=0;kET<6;kET++) {ET_ij[kET] = 0.5 * (local.ET[k][kET] + SphP[j].ET[k][kET]);}
#endif
                        double tensor = (ET_ij[0]*dp[0]*dp[0] + ET_ij[1]*dp[1]*dp[1] + ET_ij[2]*dp[2]*dp[2]
                                         + 2.*ET_ij[3]*dp[0]*dp[1] + 2.*ET_ij[4]*dp[1]*dp[2] + 2.*ET_ij[5]*dp[2]*dp[0]) / r2;
                        double kappa_ij = 0.5*(local.RT_DiffusionCoeff[k] + rt_diffusion_coefficient(j,k));
                        double fac = tensor_norm * tensor * kappa_ij;
                        out.matrixmult_out[k] -= fac * matrixmult_in[k][j];
                        out.matrixmult_sum[k] += fac;
                    }
                }
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
#ifndef DONOTUSENODELIST
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = rt_cg_DataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        } // if(mode == 1)
#endif
    } // while(startnode >= 0)
    /* Now collect the result at the right place */
    if(mode == 0)
    {
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            matrixmult_out[k][target] = out.matrixmult_out[k];
            matrixmult_sum[k][target] = out.matrixmult_sum[k];
        }
    }
    else
        rt_cg_DataResult[target] = out;
    return 0;
}



/* routine for initial loop of particles on local processor (and determination of which need passing) */
void *rt_diffusion_cg_evaluate_primary(void *p, double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum)
{
#define CONDITION_FOR_EVALUATION if((P[i].Type==0)&&(PPP[i].NumNgb>0)&&(PPP[i].Hsml>0)&&(P[i].Mass>0))
#define EVALUATION_CALL rt_diffusion_cg_evaluate(i,0,matrixmult_in,matrixmult_out,matrixmult_sum,exportflag,exportnodecount,exportindex,ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *rt_diffusion_cg_evaluate_secondary(void *p, double **matrixmult_in, double **matrixmult_out, double **matrixmult_sum)
{
#define EVALUATION_CALL rt_diffusion_cg_evaluate(j, 1, matrixmult_in, matrixmult_out, matrixmult_sum, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}

#endif
