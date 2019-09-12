#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#include "../mesh_motion.h"
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

/*! \file density.c
 *  \brief SPH density computation and kernel length determination
 *
 *  This file contains the "first hydro loop", where the gas densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the kernel length if needed.
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * substantially (condensed, different criteria for kernel lengths,
 * some optimizatins, and new variable/memory conventions added)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

struct kernel_density
{
  double dp[3],dv[3],r;
  double wk, dwk;
  double hinv, hinv3, hinv4;
  double mj_wk, mj_dwk_r;
};


/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct densdata_in
{
  MyDouble Pos[3];
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
  MyFloat Accel[3];
#endif
  MyFloat Vel[3];
  MyFloat Hsml;
#ifdef GALSF_SUBGRID_WINDS
  MyFloat DelayTime;
#endif
  int NodeList[NODELISTLENGTH];
  int Type;
}
 *DensDataIn, *DensDataGet;

static struct densdata_out
{
    MyLongDouble Ngb;
    MyLongDouble Rho;
    MyLongDouble DhsmlNgb;
    MyLongDouble Particle_DivVel;
    MyFloat NV_T[3][3];
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
    MyDouble ParticleVel[3];
#endif
#ifdef HYDRO_SPH
    MyLongDouble DhsmlHydroSumFactor;
#endif
#ifdef RT_SOURCE_INJECTION
    MyLongDouble KernelSum_Around_RT_Source;
#endif
    
#ifdef HYDRO_PRESSURE_SPH
    MyLongDouble EgyRho;
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
    MyFloat AGS_zeta;
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    MyFloat NV_D[3][3];
    MyFloat NV_A[3][3];
#endif

#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    MyFloat GradRho[3];
#endif
    
#if defined(BLACK_HOLES)
    int BH_TimeBinGasNeighbor;
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    MyDouble BH_dr_to_NearestGasNeighbor;
#endif 
#endif

#if defined(TURB_DRIVING) || defined(GRAIN_FLUID)
    MyDouble GasVel[3];
#endif
#if defined(GRAIN_FLUID)
    MyDouble Gas_InternalEnergy;
#ifdef GRAIN_LORENTZFORCE
    MyDouble Gas_B[3];
#endif
#endif

}
 *DensDataResult, *DensDataOut;

void particle2in_density(struct densdata_in *in, int i);
void out2particle_density(struct densdata_out *out, int i, int mode);
void density_evaluate_extra_physics_gas(struct densdata_in *local, struct densdata_out *out,
					struct kernel_density *kernel, int j);


void particle2in_density(struct densdata_in *in, int i)
{
    int k;
    in->Type = P[i].Type;
    in->Hsml = PPP[i].Hsml;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        if(P[i].Type==0) {in->Vel[k]=SphP[i].VelPred[k];} else {in->Vel[k]=P[i].Vel[k];}
    }
    
    if(P[i].Type == 0)
    {
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        for(k = 0; k < 3; k++)
            in->Accel[k] = All.cf_a2inv*P[i].GravAccel[k] + SphP[i].HydroAccel[k]; // PHYSICAL units //
#endif
        
#ifdef GALSF_SUBGRID_WINDS
        in->DelayTime = SphP[i].DelayTime;
#endif
    }
}


void out2particle_density(struct densdata_out *out, int i, int mode)
{
    int j,k;
    ASSIGN_ADD(PPP[i].NumNgb, out->Ngb, mode);
    ASSIGN_ADD(PPP[i].DhsmlNgbFactor, out->DhsmlNgb, mode);
    ASSIGN_ADD(P[i].Particle_DivVel, out->Particle_DivVel,   mode);
    
    if(P[i].Type == 0)
    {
        ASSIGN_ADD(SphP[i].Density, out->Rho, mode);
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
        for(k=0;k<3;k++) ASSIGN_ADD(SphP[i].ParticleVel[k], out->ParticleVel[k],   mode);
#endif
        for(k = 0; k < 3; k++) {for(j = 0; j < 3; j++) {ASSIGN_ADD(SphP[i].NV_T[k][j], out->NV_T[k][j], mode);}}

#ifdef HYDRO_SPH
        ASSIGN_ADD(SphP[i].DhsmlHydroSumFactor, out->DhsmlHydroSumFactor, mode);
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
        ASSIGN_ADD(PPPZ[i].AGS_zeta, out->AGS_zeta,   mode);
#endif

#ifdef HYDRO_PRESSURE_SPH
        ASSIGN_ADD(SphP[i].EgyWtDensity,   out->EgyRho,   mode);
#endif

#if defined(TURB_DRIVING)
        for(k = 0; k < 3; k++)
            ASSIGN_ADD(SphP[i].SmoothedVel[k], out->GasVel[k], mode);
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        for(k = 0; k < 3; k++)
            for(j = 0; j < 3; j++)
            {
                ASSIGN_ADD(SphP[i].NV_D[k][j], out->NV_D[k][j], mode);
                ASSIGN_ADD(SphP[i].NV_A[k][j], out->NV_A[k][j], mode);
            }
#endif
    } // P[i].Type == 0 //

#if defined(GRAIN_FLUID)
    if(P[i].Type > 0)
    {
        ASSIGN_ADD(P[i].Gas_Density, out->Rho, mode);
        ASSIGN_ADD(P[i].Gas_InternalEnergy, out->Gas_InternalEnergy, mode);
        for(k = 0; k<3; k++) {ASSIGN_ADD(P[i].Gas_Velocity[k], out->GasVel[k], mode);}
#ifdef GRAIN_LORENTZFORCE
        for(k = 0; k<3; k++) {ASSIGN_ADD(P[i].Gas_B[k], out->Gas_B[k], mode);}
#endif
    }
#endif

#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    if(P[i].Type != 0)
    {
        ASSIGN_ADD(P[i].DensAroundStar, out->Rho, mode);
        for(k = 0; k<3; k++) {ASSIGN_ADD(P[i].GradRho[k], out->GradRho[k], mode);}
    }
#endif
    
#if defined(RT_SOURCE_INJECTION)
    if((1 << P[i].Type) & (RT_SOURCES)) {ASSIGN_ADD(P[i].KernelSum_Around_RT_Source, out->KernelSum_Around_RT_Source, mode);}
#endif
    

#ifdef BLACK_HOLES
    if(P[i].Type == 5)
    {
        if(mode == 0) {BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;} else {if(BPP(i).BH_TimeBinGasNeighbor > out->BH_TimeBinGasNeighbor) {BPP(i).BH_TimeBinGasNeighbor = out->BH_TimeBinGasNeighbor;}}
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        if(mode == 0) {BPP(i).BH_dr_to_NearestGasNeighbor = out->BH_dr_to_NearestGasNeighbor;} else {if(BPP(i).BH_dr_to_NearestGasNeighbor > out->BH_dr_to_NearestGasNeighbor) {BPP(i).BH_dr_to_NearestGasNeighbor = out->BH_dr_to_NearestGasNeighbor;}}
#endif
    } /* if(P[i].Type == 5) */
#endif
}



/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current kernel radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its kernel region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the kernel length is adjusted accordingly, and the
 * density() computation is called again.  Note that the kernel length is
 * not allowed to fall below the lower bound set by MinHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void density(void)
{
  MyFloat *Left, *Right;
  int i, j, k, k1, k2, ndone, ndone_flag, npleft, iter = 0;
  int ngrp, recvTask, place;
  long long ntot;
  double fac, fac_lim;
  double Tinv[3][3], detT, CNumHolder=0, ConditionNumber=0;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  double desnumngb, desnumngbdev;
  int save_NextParticle;
  long long n_exported = 0;
  int redo_particle;
  int particle_set_to_minhsml_flag = 0;
  int particle_set_to_maxhsml_flag = 0;

  CPU_Step[CPU_DENSMISC] += measure_time();

  long long NTaskTimesNumPart;
  NTaskTimesNumPart = maxThreads * NumPart;
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(density_isactive(i))
        {
            Left[i] = Right[i] = 0;
#ifdef BLACK_HOLES
            P[i].SwallowID = 0;   
#ifdef SINGLE_STAR_SINK_DYNAMICS
	    P[i].SwallowTime = MAX_REAL_NUMBER;
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)	    
	    P[i].BH_Ngb_Flag = 0;
#endif		
#endif
        }
    } /* done with intial zero-out loop */

  /* allocate buffers to arrange communication */
  size_t MyBufferSize = All.BufferSize;
  All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct densdata_in) + sizeof(struct densdata_out) +
					     sizemax(sizeof(struct densdata_in),sizeof(struct densdata_out))));
  DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  t0 = my_second();

  desnumngb = All.DesNumNgb;
  desnumngbdev = All.MaxNumNgbDeviation;
    /* in the initial timestep and iteration, use a much more strict tolerance for the neighbor number */
    if(All.Time==All.TimeBegin) {if(All.MaxNumNgbDeviation > 0.05) desnumngbdev=0.05;}
    double desnumngbdev_0 = desnumngbdev;

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
	      pthread_create(&mythreads[j], &attr, density_evaluate_primary, &threadid[j]);
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
	    density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
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
		  printf("Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n",ThisTask,P[NextParticle].Type,
			 P[NextParticle].Pos[0],P[NextParticle].Pos[1],P[NextParticle].Pos[2],P[NextParticle].Mass);
		  if(P[NextParticle].Type == 0)
		    printf("   rho=%g hsml=%g\n",SphP[NextParticle].Density,PPP[NextParticle].Hsml);

		  endrun(112208);
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

	  DensDataGet = (struct densdata_in *) mymalloc("DensDataGet", Nimport * sizeof(struct densdata_in));
	  DensDataIn = (struct densdata_in *) mymalloc("DensDataIn", Nexport * sizeof(struct densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      particle2in_density(&DensDataIn[j], place);

	      memcpy(DensDataIn[j].NodeList,
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
		      MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = my_second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(DensDataIn);
	  DensDataResult =
	    (struct densdata_out *) mymalloc("DensDataResult", Nimport * sizeof(struct densdata_out));
	  DensDataOut =
	    (struct densdata_out *) mymalloc("DensDataOut", Nexport * sizeof(struct densdata_out));

	  report_memory_usage(&HighMark_sphdensity, "SPH_DENSITY");

	  /* now do the particles that were sent to us */

	  tstart = my_second();

	  NextJ = 0;

#ifdef PTHREADS_NUM_THREADS
	  for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
	    pthread_create(&mythreads[j], &attr, density_evaluate_secondary, &threadid[j]);
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
	    density_evaluate_secondary(&mainthreadid);
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
		      MPI_Sendrecv(&DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
	      out2particle_density(&DensDataOut[j], place, 1);
	    }
	  tend = my_second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(DensDataOut);
	  myfree(DensDataResult);
	  myfree(DensDataGet);
	}
      while(ndone < NTask);


        /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(density_isactive(i))
            {
                if(PPP[i].NumNgb > 0)
                {
                    PPP[i].DhsmlNgbFactor *= PPP[i].Hsml / (NUMDIMS * PPP[i].NumNgb);
                    P[i].Particle_DivVel /= PPP[i].NumNgb;
                    /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                    PPP[i].NumNgb *= NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS);
                } else {
                    PPP[i].NumNgb = PPP[i].DhsmlNgbFactor = P[i].Particle_DivVel = 0;
                }
#if defined(ADAPTIVE_GRAVSOFT_FORALL) /* if particle is AGS-active and non-gas, set DivVel to zero because it will be reset in ags_hsml routine */
                if(ags_density_isactive(i) && (P[i].Type > 0)) {PPP[i].Particle_DivVel = 0;}
#endif
                
                // inverse of SPH volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(PPP[i].DhsmlNgbFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
                    PPP[i].DhsmlNgbFactor = 1 / (1 + PPP[i].DhsmlNgbFactor);
                else
                    PPP[i].DhsmlNgbFactor = 1;
                P[i].Particle_DivVel *= PPP[i].DhsmlNgbFactor;
            
                if(P[i].Type == 0)
                {
                    /* fill in the missing elements of NV_T (it's symmetric, so we saved time not computing these directly) */
                    SphP[i].NV_T[1][0]=SphP[i].NV_T[0][1]; SphP[i].NV_T[2][0]=SphP[i].NV_T[0][2]; SphP[i].NV_T[2][1]=SphP[i].NV_T[1][2];
                    /* Now invert the NV_T matrix we just measured */
                    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
                        this will tell us how robust our procedure is (and let us know if we need to expand the neighbor number */
                    ConditionNumber=CNumHolder=0;
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {ConditionNumber += SphP[i].NV_T[k1][k2]*SphP[i].NV_T[k1][k2];}}
#if (NUMDIMS==1)
                    /* one-dimensional case */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Tinv[k1][k2]=0;}}
                    detT = SphP[i].NV_T[0][0];
                    if(SphP[i].NV_T[0][0]!=0 && !isnan(SphP[i].NV_T[0][0])) Tinv[0][0] = 1/detT; /* only one non-trivial element in 1D! */
#endif
#if (NUMDIMS==2)
                    /* two-dimensional case */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Tinv[k1][k2]=0;}}
                    detT = SphP[i].NV_T[0][0]*SphP[i].NV_T[1][1] - SphP[i].NV_T[0][1]*SphP[i].NV_T[1][0];
                    if((detT != 0)&&(!isnan(detT)))
                    {
                        Tinv[0][0] = SphP[i].NV_T[1][1] / detT;
                        Tinv[0][1] = -SphP[i].NV_T[0][1] / detT;
                        Tinv[1][0] = -SphP[i].NV_T[1][0] / detT;
                        Tinv[1][1] = SphP[i].NV_T[0][0] / detT;
                    }
#endif
#if (NUMDIMS==3)
                    /* three-dimensional case */
                    detT = SphP[i].NV_T[0][0] * SphP[i].NV_T[1][1] * SphP[i].NV_T[2][2] +
                        SphP[i].NV_T[0][1] * SphP[i].NV_T[1][2] * SphP[i].NV_T[2][0] +
                        SphP[i].NV_T[0][2] * SphP[i].NV_T[1][0] * SphP[i].NV_T[2][1] -
                        SphP[i].NV_T[0][2] * SphP[i].NV_T[1][1] * SphP[i].NV_T[2][0] -
                        SphP[i].NV_T[0][1] * SphP[i].NV_T[1][0] * SphP[i].NV_T[2][2] -
                        SphP[i].NV_T[0][0] * SphP[i].NV_T[1][2] * SphP[i].NV_T[2][1];
                    /* check for zero determinant */
                    if((detT != 0) && !isnan(detT))
                    {
                        Tinv[0][0] = (SphP[i].NV_T[1][1] * SphP[i].NV_T[2][2] - SphP[i].NV_T[1][2] * SphP[i].NV_T[2][1]) / detT;
                        Tinv[0][1] = (SphP[i].NV_T[0][2] * SphP[i].NV_T[2][1] - SphP[i].NV_T[0][1] * SphP[i].NV_T[2][2]) / detT;
                        Tinv[0][2] = (SphP[i].NV_T[0][1] * SphP[i].NV_T[1][2] - SphP[i].NV_T[0][2] * SphP[i].NV_T[1][1]) / detT;
                        Tinv[1][0] = (SphP[i].NV_T[1][2] * SphP[i].NV_T[2][0] - SphP[i].NV_T[1][0] * SphP[i].NV_T[2][2]) / detT;
                        Tinv[1][1] = (SphP[i].NV_T[0][0] * SphP[i].NV_T[2][2] - SphP[i].NV_T[0][2] * SphP[i].NV_T[2][0]) / detT;
                        Tinv[1][2] = (SphP[i].NV_T[0][2] * SphP[i].NV_T[1][0] - SphP[i].NV_T[0][0] * SphP[i].NV_T[1][2]) / detT;
                        Tinv[2][0] = (SphP[i].NV_T[1][0] * SphP[i].NV_T[2][1] - SphP[i].NV_T[1][1] * SphP[i].NV_T[2][0]) / detT;
                        Tinv[2][1] = (SphP[i].NV_T[0][1] * SphP[i].NV_T[2][0] - SphP[i].NV_T[0][0] * SphP[i].NV_T[2][1]) / detT;
                        Tinv[2][2] = (SphP[i].NV_T[0][0] * SphP[i].NV_T[1][1] - SphP[i].NV_T[0][1] * SphP[i].NV_T[1][0]) / detT;
                    } else {
                        for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Tinv[k1][k2]=0;}}
                    }
#endif
                    
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {CNumHolder += Tinv[k1][k2]*Tinv[k1][k2];}}
                    ConditionNumber = sqrt(ConditionNumber*CNumHolder) / NUMDIMS;
                    if(ConditionNumber<1) ConditionNumber=1;
                    /* this = sqrt( ||NV_T^-1||*||NV_T|| ) :: should be ~1 for a well-conditioned matrix */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {SphP[i].NV_T[k1][k2]=Tinv[k1][k2];}}
                    /* now NV_T holds the inverted matrix elements, for use in hydro */
                } // P[i].Type == 0 //
                
                /* now check whether we had enough neighbours */
                double ncorr_ngb = 1.0;
                double cn=1;
                double c0 = 0.1 * (double)CONDITION_NUMBER_DANGER;
                if(P[i].Type==0)
                {
                    /* use the previous timestep condition number to correct how many neighbors we should use for stability */
                    if((iter==0)&&(ConditionNumber>SphP[i].ConditionNumber))
                    {
                        /* if we find ourselves with a sudden increase in condition number - check if we have a reasonable 
                            neighbor number for the previous iteration, and if so, use the new (larger) correction */
                        ncorr_ngb=1; cn=SphP[i].ConditionNumber; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
                        double dn_ngb = fabs(PPP[i].NumNgb-All.DesNumNgb*ncorr_ngb)/(desnumngbdev_0*ncorr_ngb);
                        ncorr_ngb=1; cn=ConditionNumber; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
                        double dn_ngb_alt = fabs(PPP[i].NumNgb-All.DesNumNgb*ncorr_ngb)/(desnumngbdev_0*ncorr_ngb);
                        dn_ngb = DMIN(dn_ngb,dn_ngb_alt);
                        if(dn_ngb < 10.0) SphP[i].ConditionNumber = ConditionNumber;
                    }
                    ncorr_ngb=1; cn=SphP[i].ConditionNumber; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
                }
                desnumngb = All.DesNumNgb * ncorr_ngb;
                desnumngbdev = desnumngbdev_0 * ncorr_ngb;
                /* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}

#ifdef BLACK_HOLES
                if(P[i].Type == 5)
                {
                    desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;
#ifdef SINGLE_STAR_SINK_DYNAMICS		    
                    desnumngbdev = (All.BlackHoleNgbFactor+1);
#else
                    desnumngbdev = 4 * (All.BlackHoleNgbFactor+1);     
#endif		    
                }
#endif

#ifdef GRAIN_FLUID
                /* for the grains, we only need to estimate neighboring gas properties, we don't need to worry about 
                    condition numbers or conserving an exact neighbor number */
                if(P[i].Type>0)
                {
                    desnumngb = All.DesNumNgb;
                    desnumngbdev = All.DesNumNgb / 4;
#ifdef GRAIN_BACKREACTION
                    desnumngbdev = desnumngbdev_0;
#endif
                }
#endif

                double minsoft = All.MinHsml;
                double maxsoft = All.MaxHsml;

#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
                /* use a much looser check for N_neighbors when the central point is a star particle,
                 since the accuracy is limited anyways to the coupling efficiency -- the routines use their
                 own estimators+neighbor loops, anyways, so this is just to get some nearby particles */
                if((P[i].Type!=0)&&(P[i].Type!=5))
                {
                    desnumngb = All.DesNumNgb;
#if defined(RT_SOURCE_INJECTION)
                    if(desnumngb < 64.0) {desnumngb = 64.0;} // we do want a decent number to ensure the area around the particle is 'covered'
#endif
#ifdef GALSF
                    if(desnumngb < 64.0) {desnumngb = 64.0;} // we do want a decent number to ensure the area around the particle is 'covered'
                    // if we're finding this for feedback routines, there isn't any good reason to search beyond a modest physical radius //
                    double unitlength_in_kpc=All.UnitLength_in_cm/All.HubbleParam/3.086e21*All.cf_atime;
                    maxsoft = 2.0 / unitlength_in_kpc;
#endif
                    desnumngbdev = desnumngb / 2; // enforcing exact number not important
                }
#endif
                
#ifdef BLACK_HOLES
                if(P[i].Type == 5) {maxsoft = All.BlackHoleMaxAccretionRadius / All.cf_atime;}  // MaxAccretionRadius is now defined in params.txt in PHYSICAL units
#ifdef SINGLE_STAR_SINK_DYNAMICS
		        if(P[i].Type == 5) {minsoft = All.ForceSoftening[5] / All.cf_atime;} // we should always find all neighbours within the softening kernel/accretion radius, which is a lower bound on the accretion radius
#endif		
#endif

                redo_particle = 0;
                
                /* check if we are in the 'normal' range between the max/min allowed values */
                if((PPP[i].NumNgb < (desnumngb - desnumngbdev) && PPP[i].Hsml < 0.999*maxsoft) ||
                   (PPP[i].NumNgb > (desnumngb + desnumngbdev) && PPP[i].Hsml > 1.001*minsoft))
                    redo_particle = 1;
                
                /* check maximum kernel size allowed */
                particle_set_to_maxhsml_flag = 0;
                if((PPP[i].Hsml >= 0.999*maxsoft) && (PPP[i].NumNgb < (desnumngb - desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].Hsml == maxsoft)
                    {
                        /* iteration at the maximum value is already complete */
                        particle_set_to_maxhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        PPP[i].Hsml = maxsoft;
                        particle_set_to_maxhsml_flag = 1;
                    }
                }
                
                /* check minimum kernel size allowed */
                particle_set_to_minhsml_flag = 0;
                if((PPP[i].Hsml <= 1.001*minsoft) && (PPP[i].NumNgb > (desnumngb + desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].Hsml == minsoft)
                    {
                        /* this means we've already done an iteration with the MinHsml value, so the
                         neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
                        particle_set_to_minhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        PPP[i].Hsml = minsoft;
                        particle_set_to_minhsml_flag = 1;
                    }
                }
                
#ifdef GALSF
                if((All.ComovingIntegrationOn)&&(All.Time>All.TimeBegin))
                {
                    if((P[i].Type==4)&&(iter>1)&&(PPP[i].NumNgb>4)&&(PPP[i].NumNgb<100)&&(redo_particle==1)) {redo_particle=0;}
                }
#endif    
                
                if((redo_particle==0)&&(P[i].Type == 0))
                {
                    /* ok we have reached the desired number of neighbors: save the condition number for next timestep */
                    if(ConditionNumber > 1000.0 * (double)CONDITION_NUMBER_DANGER)
                    {
#ifndef IO_REDUCED_MODE
                        printf("Warning: Condition number=%g CNum_prevtimestep=%g Num_Ngb=%g desnumngb=%g Hsml=%g Hsml_min=%g Hsml_max=%g\n",
                               ConditionNumber,SphP[i].ConditionNumber,PPP[i].NumNgb,desnumngb,PPP[i].Hsml,All.MinHsml,All.MaxHsml);
#endif
                    }
                    SphP[i].ConditionNumber = ConditionNumber;
                }
                
                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
                        printf("i=%d task=%d ID=%llu Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)\n",
                               i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, PPP[i].Hsml, PPP[i].DhsmlNgbFactor, Left[i], Right[i],
                               (float) PPP[i].NumNgb, Right[i] - Left[i], particle_set_to_maxhsml_flag, particle_set_to_minhsml_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                    }
                    
                    /* need to redo this particle */
                    npleft++;
                    
                    if(Left[i] > 0 && Right[i] > 0)
                        if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                            SphP[i].ConditionNumber = ConditionNumber;
                            continue;
                        }
                    
                    if((particle_set_to_maxhsml_flag==0)&&(particle_set_to_minhsml_flag==0))
                    {
                        if(PPP[i].NumNgb < (desnumngb - desnumngbdev))
                            Left[i] = DMAX(PPP[i].Hsml, Left[i]);
                        else
                        {
                            if(Right[i] != 0)
                            {
                                if(PPP[i].Hsml < Right[i])
                                    Right[i] = PPP[i].Hsml;
                            }
                            else
                                Right[i] = PPP[i].Hsml;
                        }
                        
                        // right/left define upper/lower bounds from previous iterations
                        if(Right[i] > 0 && Left[i] > 0)
                        {
                            // geometric interpolation between right/left //
                            double maxjump=0;
                            if(iter>1) {maxjump = 0.2*log(Right[i]/Left[i]);}
                            if(PPP[i].NumNgb > 1)
                            {
                                double jumpvar = PPP[i].DhsmlNgbFactor * log( desnumngb / PPP[i].NumNgb ) / NUMDIMS;
                                if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
                                PPP[i].Hsml *= exp(jumpvar);
                            } else {
                                PPP[i].Hsml *= 2.0;
                            }
                            if((PPP[i].Hsml<Right[i])&&(PPP[i].Hsml>Left[i]))
                            {
                                if(iter > 1)
                                {
                                    double hfac = exp(maxjump);
                                    if(PPP[i].Hsml > Right[i] / hfac) {PPP[i].Hsml = Right[i] / hfac;}
                                    if(PPP[i].Hsml < Left[i] * hfac) {PPP[i].Hsml = Left[i] * hfac;}
                                }
                            } else {
                                if(PPP[i].Hsml>Right[i]) PPP[i].Hsml=Right[i];
                                if(PPP[i].Hsml<Left[i]) PPP[i].Hsml=Left[i];
                                PPP[i].Hsml = pow(PPP[i].Hsml * Left[i] * Right[i] , 1.0/3.0);
                            }
                        }
                        else
                        {
                            if(Right[i] == 0 && Left[i] == 0)
                            {
                                char buf[1000];
                                sprintf(buf, "Right[i] == 0 && Left[i] == 0 && PPP[i].Hsml=%g\n", PPP[i].Hsml);
                                terminate(buf);
                            }
                            
                            if(Right[i] == 0 && Left[i] > 0)
                            {
                                if (PPP[i].NumNgb > 1)
                                    fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (+0.231=2x desnumngb)
                                else
                                    fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                                
                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4)
                                        if(PPP[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                    
                                    if(fac < fac_lim+0.231)
                                    {
                                        PPP[i].Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                    {
                                        PPP[i].Hsml *= exp(fac_lim+0.231);
                                        // fac~0.26 leads to expected doubling of number if density is constant,
                                        //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                    }
                                }
                                else
                                    PPP[i].Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                            
                            if(Right[i] > 0 && Left[i] == 0)
                            {
                                if (PPP[i].NumNgb > 1)
                                    fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                                else
                                    fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                                
                                if (fac_lim < -1.535) fac_lim = -1.535; // decreasing N_ngb by factor ~100
                                
                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4)
                                        if(PPP[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                    
                                    if(fac > fac_lim-0.231)
                                    {
                                        PPP[i].Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                        PPP[i].Hsml *= exp(fac_lim-0.231); // limiter to prevent --too-- far a jump in a single iteration
                                }
                                else
                                    PPP[i].Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                        } // closes if[particle_set_to_max/minhsml_flag]
                    } // closes redo_particle
                    /* resets for max/min values */
                    if(PPP[i].Hsml < minsoft) PPP[i].Hsml = minsoft;
                    if(particle_set_to_minhsml_flag==1) PPP[i].Hsml = minsoft;
                    if(PPP[i].Hsml > maxsoft) PPP[i].Hsml = maxsoft;
                    if(particle_set_to_maxhsml_flag==1) PPP[i].Hsml = maxsoft;
                }
                else
                    P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
            }
        }
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
                printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
                       (int) (ntot / 1000000000), (int) (ntot % 1000000000));
            }
            if(iter > MAXITER)
            {
                printf("failed to converge in neighbour iteration in density()\n");
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
        if(P[i].TimeBin < 0) {P[i].TimeBin = -P[i].TimeBin - 1;}
    }
    
    
    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results
     ( any quantities that only need to be evaluated once, on the final iteration --
     won't save much b/c the real cost is in the neighbor loop for each particle, but it's something )
     -- also, some results (for example, viscosity suppression below) should not be calculated unless
     the quantities are 'stabilized' at their final values -- */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(density_isactive(i))
        {
            if(P[i].Type == 0 && P[i].Mass > 0)
            {
                if(SphP[i].Density > 0)
                {
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
                    /* set motion of the mesh-generating points */
#if (HYDRO_FIX_MESH_MOTION==4)
                    set_mesh_motion(i); // use user-specified analytic function to define mesh motions //
#elif ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                    double eps_pvel = 0.3; // normalization for how much 'weight' to give to neighbors (unstable if >=0.5)
                    for(k=0;k<3;k++) {SphP[i].ParticleVel[k] = (1.-eps_pvel)*SphP[i].VelPred[k] + eps_pvel*SphP[i].ParticleVel[k]/SphP[i].Density;} // assign mixture velocity
#elif (HYDRO_FIX_MESH_MOTION==7)
                    for(k=0;k<3;k++) {SphP[i].ParticleVel[k] = SphP[i].VelPred[k];} // move with fluid
#endif
#endif

#ifdef HYDRO_SPH
#ifdef HYDRO_PRESSURE_SPH
                    if(SphP[i].InternalEnergyPred > 0)
                    {
                        SphP[i].EgyWtDensity /= SphP[i].InternalEnergyPred;
                    } else {
                        SphP[i].EgyWtDensity = 0;
                    }
#endif
                    /* need to divide by the sum of x_tilde=1, i.e. numden_ngb */
                    if((PPP[i].Hsml > 0)&&(PPP[i].NumNgb > 0))
                    {
                        double numden_ngb = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
                        SphP[i].DhsmlHydroSumFactor *= PPP[i].Hsml / (NUMDIMS * numden_ngb);
                        SphP[i].DhsmlHydroSumFactor *= -PPP[i].DhsmlNgbFactor; /* now this is ready to be called in hydro routine */
                    } else {
                        SphP[i].DhsmlHydroSumFactor = 0;
                    }
#endif
                    
              
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
                    for(k1 = 0; k1 < 3; k1++)
                        for(k2 = 0; k2 < 3; k2++)
                        {
                            SphP[i].NV_D[k2][k1] *= All.cf_a2inv; // converts to physical velocity/length
                            SphP[i].NV_A[k2][k1] /= All.cf_atime; // converts to physical accel/length
                        }
                    // all quantities below in this block should now be in proper PHYSICAL units, for subsequent operations //
                    double dtDV[3][3], A[3][3], V[3][3], S[3][3];
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            V[k1][k2] = SphP[i].NV_D[k1][0]*SphP[i].NV_T[0][k2] + SphP[i].NV_D[k1][1]*SphP[i].NV_T[1][k2] + SphP[i].NV_D[k1][2]*SphP[i].NV_T[2][k2];
                            A[k1][k2] = SphP[i].NV_A[k1][0]*SphP[i].NV_T[0][k2] + SphP[i].NV_A[k1][1]*SphP[i].NV_T[1][k2] + SphP[i].NV_A[k1][2]*SphP[i].NV_T[2][k2];
                        }
                    SphP[i].NV_DivVel = V[0][0] + V[1][1] + V[2][2];
                    SphP[i].NV_trSSt = 0;
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            dtDV[k1][k2] = A[k1][k2] - (V[k1][0]*V[0][k2] + V[k1][1]*V[1][k2] + V[k1][2]*V[2][k2]);
                            /* S = 0.5*(V+V_transpose) - delta_ij*div_v/3 */
                            S[k1][k2] = 0.5 * (V[k1][k2] + V[k2][k1]);
                            if(k2==k1) S[k1][k2] -= SphP[i].NV_DivVel / NUMDIMS;
                            /* Trace[S*S_transpose] = SSt[0][0]+SSt[1][1]+SSt[2][2] = |S|^2 = sum(Sij^2) */
                            SphP[i].NV_trSSt += S[k1][k2]*S[k1][k2];
                        }
                    SphP[i].NV_dt_DivVel = dtDV[0][0] + dtDV[1][1] + dtDV[2][2];
#endif
                    
                    
#if defined(TURB_DRIVING)
                    if(SphP[i].Density > 0)
                    {
                        SphP[i].SmoothedVel[0] /= SphP[i].Density;
                        SphP[i].SmoothedVel[1] /= SphP[i].Density;
                        SphP[i].SmoothedVel[2] /= SphP[i].Density;
                    } else {
                        SphP[i].SmoothedVel[0] = SphP[i].SmoothedVel[1] = SphP[i].SmoothedVel[2] = 0;
                    }
#endif
                }
                
#ifndef HYDRO_SPH
                if((PPP[i].Hsml > 0)&&(PPP[i].NumNgb > 0))
                {
                    SphP[i].Density = P[i].Mass * PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) ); // divide mass by volume
                } else {
                    if(PPP[i].Hsml <= 0)
                    {
                        SphP[i].Density = 0; // in this case, give up, no meaningful volume
                    } else {
                        SphP[i].Density = P[i].Mass / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) ); // divide mass (lone particle) by volume
                    }
                }
#endif
                SphP[i].Pressure = get_pressure(i);		// should account for density independent pressure

            } // P[i].Type == 0

            
#if defined(GRAIN_FLUID)
            if(P[i].Type > 0)
            {
                if(P[i].Gas_Density > 0)
                {
                    P[i].Gas_InternalEnergy /= P[i].Gas_Density;
                    for(k = 0; k<3; k++) {P[i].Gas_Velocity[k] /= P[i].Gas_Density;}
                } else {
                    P[i].Gas_InternalEnergy = 0;
                    for(k = 0; k<3; k++) {P[i].Gas_Velocity[k] = 0;}
#ifdef GRAIN_LORENTZFORCE
                    for(k = 0; k<3; k++) {P[i].Gas_B[k] = 0;}
#endif
                }
            }
#endif
            
            
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
            /* non-gas particles are handled separately, in the ags_hsml routine */
            if(P[i].Type==0)
            {
                PPPZ[i].AGS_zeta = 0;
                double zeta_0 = 0; // 2.0 * P[i].Mass*P[i].Mass * PPP[i].Hsml*PPP[i].Hsml; // self-value of zeta if no neighbors are found //
                if((PPP[i].Hsml > 0)&&(PPP[i].NumNgb > 0))
                {
                    /* the zeta terms ONLY control errors if we maintain the 'correct' neighbor number: for boundary
                        particles, it can actually be worse. so we need to check whether we should use it or not */
                    if((PPP[i].Hsml > 1.001*All.MinHsml) && (PPP[i].Hsml < 0.999*All.MaxHsml) &&
                        (fabs(PPP[i].NumNgb-All.DesNumNgb)/All.DesNumNgb < 0.05))
                    {
                        double ndenNGB = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].Hsml,NUMDIMS) );
                        PPPZ[i].AGS_zeta *= 0.5 * P[i].Mass * PPP[i].Hsml / (NUMDIMS * ndenNGB) * PPP[i].DhsmlNgbFactor;
                    } else {
                        PPPZ[i].AGS_zeta = zeta_0;
                    }
                } else {
                    PPPZ[i].AGS_zeta = zeta_0;
                }
            }
#endif
            
#ifdef PM_HIRES_REGION_CLIPPING
#ifdef BLACK_HOLES
            if (P[i].Type != 5)
            {
#endif
                if(P[i].Type == 0) if ((SphP[i].Density <= 0) || (PPP[i].NumNgb <= 0)) P[i].Mass = 0;
                if ((PPP[i].Hsml <= 0) || (PPP[i].Hsml >= PM_HIRES_REGION_CLIPPING)) P[i].Mass = 0;
                double vmag=0; for(k=0;k<3;k++) vmag+=P[i].Vel[k]*P[i].Vel[k]; vmag = sqrt(vmag);
                if(vmag>5.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) P[i].Mass=0;
                if(vmag>1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) for(k=0;k<3;k++) P[i].Vel[k]*=(1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s)/vmag;
#ifdef BLACK_HOLES
            }
#endif // BLACK_HOLES
#endif // ifdef PM_HIRES_REGION_CLIPPING
            
            
         /* finally, convert NGB to the more useful format, NumNgb^(1/NDIMS),
            which we can use to obtain the corrected particle sizes. Because of how this number is used above, we --must-- make 
            sure that this operation is the last in the loop here */
            if(PPP[i].NumNgb > 0) {PPP[i].NumNgb=pow(PPP[i].NumNgb,1./NUMDIMS);} else {PPP[i].NumNgb=0;}
            
        } // density_isactive(i)
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    
    
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_DENSCOMPUTE] += timecomp;
    CPU_Step[CPU_DENSWAIT] += timewait;
    CPU_Step[CPU_DENSCOMM] += timecomm;
    CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}






/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
                     int *ngblist)
{
    int j, n;
    int startnode, numngb_inbox, listindex = 0;
    double r2, h2, u, mass_j, wk;
    struct kernel_density kernel;
    struct densdata_in local;
    struct densdata_out out;
    memset(&out, 0, sizeof(struct densdata_out));
#if defined(BLACK_HOLES)
    out.BH_TimeBinGasNeighbor = TIMEBINS;
#ifdef BH_ACCRETE_NEARESTFIRST
    out.BH_dr_to_NearestGasNeighbor = MAX_REAL_NUMBER;
#endif 
#endif
    
    if(mode == 0)
        particle2in_density(&local, target);
    else
        local = DensDataGet[target];
    h2 = local.Hsml * local.Hsml;
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = DensDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            
            if(numngb_inbox < 0) return -1;
            
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
#ifdef GALSF_SUBGRID_WINDS
                if(SphP[j].DelayTime > 0)	/* partner is a wind particle */
                    if(!(local.DelayTime > 0))	/* if I'm not wind, then ignore the wind particle */
                        continue;
#endif
                if(P[j].Mass <= 0) continue;
                
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                
                if(r2 < h2)
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                    mass_j = P[j].Mass;
                    kernel.mj_wk = FLT(mass_j * kernel.wk);
                    
                    out.Ngb += kernel.wk;
                    out.Rho += kernel.mj_wk;

#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                    if(local.Type == 0 && kernel.r==0) {int kv; for(kv=0;kv<3;kv++) {out.ParticleVel[kv] += kernel.mj_wk * SphP[j].VelPred[kv];}} // just the self-contribution //
#endif
#if defined(RT_SOURCE_INJECTION)
                    if((1 << local.Type) & (RT_SOURCES)) {out.KernelSum_Around_RT_Source += 1.-u*u;}
#endif
                    out.DhsmlNgb += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#ifdef HYDRO_SPH
                    double mass_eff = mass_j;
#ifdef HYDRO_PRESSURE_SPH
                    mass_eff *= SphP[j].InternalEnergyPred;
                    out.EgyRho += kernel.wk * mass_eff;
#endif
                    out.DhsmlHydroSumFactor += -mass_eff * (NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#endif
                    
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || (ADAPTIVE_GRAVSOFT_FORALL & 1)
                    if(local.Type == 0) {out.AGS_zeta += mass_j * kernel_gravity(u, kernel.hinv, kernel.hinv3, 0);}
#endif
                    /* for everything below, we do NOT include the particle self-contribution! */
                    if(kernel.r > 0)
                    {
                        if(local.Type == 0)
                        {
                            wk = kernel.wk; /* MAKE SURE THIS MATCHES CHOICE IN GRADIENTS.c!!! */
                            /* the weights for the MLS tensor used for gradient estimation */
                            out.NV_T[0][0] +=  wk * kernel.dp[0] * kernel.dp[0];
                            out.NV_T[0][1] +=  wk * kernel.dp[0] * kernel.dp[1];
                            out.NV_T[0][2] +=  wk * kernel.dp[0] * kernel.dp[2];
                            out.NV_T[1][1] +=  wk * kernel.dp[1] * kernel.dp[1];
                            out.NV_T[1][2] +=  wk * kernel.dp[1] * kernel.dp[2];
                            out.NV_T[2][2] +=  wk * kernel.dp[2] * kernel.dp[2];
                        }
                        kernel.dv[0] = local.Vel[0] - SphP[j].VelPred[0];
                        kernel.dv[1] = local.Vel[1] - SphP[j].VelPred[1];
                        kernel.dv[2] = local.Vel[2] - SphP[j].VelPred[2];
#ifdef BOX_SHEARING
                        if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
                        if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
#endif
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                        // do neighbor contribution to smoothed particle velocity here, after wrap, so can account for shearing boxes correctly //
                        {int kv; for(kv=0;kv<3;kv++) {out.ParticleVel[kv] += kernel.mj_wk * (local.Vel[kv] - kernel.dv[kv]);}}
#endif
                        out.Particle_DivVel -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how Hsml will evolve (particle drift) */
                        
                        density_evaluate_extra_physics_gas(&local, &out, &kernel, j);
                    } // kernel.r > 0 //
                }
            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DensDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }
    
    if(mode == 0)
        out2particle_density(&out, target, 0);
    else
        DensDataResult[target] = out;
    
    return 0;
}



void *density_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(density_isactive(i))
#define EVALUATION_CALL density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *density_evaluate_secondary(void *p)
{
#define EVALUATION_CALL density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}




int density_isactive(int n)
{
    /* first check our 'marker' for particles which have finished iterating to an Hsml solution (if they have, dont do them again) */
    if(P[n].TimeBin < 0) return 0;
    
#if defined(GRAIN_FLUID)
    /* all particles can potentially interact with the gas in this mode, if drag > 0 */
    if(P[n].Type >= 0) return 1;
#endif
    
#if defined(RT_SOURCE_INJECTION)
    if((1 << P[n].Type) & (RT_SOURCES)) 
    {
#if defined(GALSF)
       if(((P[n].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[n].Type == 2)||(P[n].Type==3))))&&(P[n].Mass>0)) 
        {
            double star_age = evaluate_stellar_age_Gyr(P[n].StellarAge);
	        if((star_age < 0.1)&&(star_age > 0)&&(!isnan(star_age))) return 1;
        }
#else
        if(Flag_FullStep) {return 1;} // only do on full timesteps
#endif
    }
#endif
    
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
    if(((P[n].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[n].Type == 2)||(P[n].Type==3))))&&(P[n].Mass>0))
    {
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
        /* check if there is going to be a SNe this timestep, in which case, we want the density info! */
        if(P[n].SNe_ThisTimeStep>0) return 1;
#endif
#if defined(GALSF)
        if(P[n].DensAroundStar<=0) return 1;
        // only do stellar age evaluation if we have to //
        if(All.ComovingIntegrationOn==0)
        {
            float star_age=0;
            star_age = evaluate_stellar_age_Gyr(P[n].StellarAge);
            if(star_age < 0.035) return 1;
        }
#endif
    }
#endif
    
#ifdef BLACK_HOLES
    if(P[n].Type == 5) return 1;
#endif
    
    if(P[n].Type == 0 && P[n].Mass > 0) return 1;
    
    return 0;
}





void density_evaluate_extra_physics_gas(struct densdata_in *local, struct densdata_out *out,
                                        struct kernel_density *kernel, int j)
{
    kernel->mj_dwk_r = P[j].Mass * kernel->dwk / kernel->r;


    if(local->Type != 0)
    {
        
#if defined(GRAIN_FLUID)
        out->Gas_InternalEnergy += kernel->mj_wk * SphP[j].InternalEnergyPred;
        out->GasVel[0] += kernel->mj_wk * (local->Vel[0]-kernel->dv[0]);
        out->GasVel[1] += kernel->mj_wk * (local->Vel[1]-kernel->dv[1]);
        out->GasVel[2] += kernel->mj_wk * (local->Vel[2]-kernel->dv[2]);
#ifdef GRAIN_LORENTZFORCE
        out->Gas_B[0] += kernel->wk * SphP[j].BPred[0];
        out->Gas_B[1] += kernel->wk * SphP[j].BPred[1];
        out->Gas_B[2] += kernel->wk * SphP[j].BPred[2];
#endif
#endif
        
#if defined(BLACK_HOLES)
        if(local->Type == 5)
        {
            P[j].SwallowID = 0;  // this way we don't have to do a global loop over local particles in blackhole_accretion() to reset these quantities...
            short int TimeBin_j = P[j].TimeBin; if(TimeBin_j < 0) {TimeBin_j = -TimeBin_j - 1;} // need to make sure we correct for the fact that TimeBin is used as a 'switch' here to determine if a particle is active for iteration, otherwise this gives nonsense!
            if(out->BH_TimeBinGasNeighbor > TimeBin_j) {out->BH_TimeBinGasNeighbor = TimeBin_j;}
#if (SINGLE_STAR_SINK_FORMATION & 8)
	    P[j].BH_Ngb_Flag = 1;
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
	    P[j].SwallowTime = MAX_REAL_NUMBER;
#endif	    
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(FLAG_NOT_IN_PUBLIC_CODE)
            double dr_eff_wtd = Get_Particle_Size(j); dr_eff_wtd=sqrt(dr_eff_wtd*dr_eff_wtd + (kernel->r)*(kernel->r)); /* effective distance for Gaussian-type kernel, weighted by density */
            if((dr_eff_wtd < out->BH_dr_to_NearestGasNeighbor) && (P[j].Mass > 0)) {out->BH_dr_to_NearestGasNeighbor = dr_eff_wtd;}
#endif
        }
#endif
        
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
        /* this is here because for the models of BH growth and self-shielding of stars, we
         just need a quick-and-dirty, single-pass approximation for the gradients (the error from
         using this as opposed to the higher-order gradient estimators is small compared to the
         Sobolev approximation): use only for -non-gas- particles */
        out->GradRho[0] += kernel->mj_dwk_r * kernel->dp[0];
        out->GradRho[1] += kernel->mj_dwk_r * kernel->dp[1];
        out->GradRho[2] += kernel->mj_dwk_r * kernel->dp[2];
#endif
        
    } else { /* local.Type == 0 */

#if defined(TURB_DRIVING)
        out->GasVel[0] += kernel->mj_wk * (local->Vel[0]-kernel->dv[0]);
        out->GasVel[1] += kernel->mj_wk * (local->Vel[1]-kernel->dv[1]);
        out->GasVel[2] += kernel->mj_wk * (local->Vel[2]-kernel->dv[2]);
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        double wk = kernel->wk;
        out->NV_A[0][0] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dp[0] * wk;
        out->NV_A[0][1] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dp[1] * wk;
        out->NV_A[0][2] += (local->Accel[0] - All.cf_a2inv*P[j].GravAccel[0] - SphP[j].HydroAccel[0]) * kernel->dp[2] * wk;
        out->NV_A[1][0] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dp[0] * wk;
        out->NV_A[1][1] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dp[1] * wk;
        out->NV_A[1][2] += (local->Accel[1] - All.cf_a2inv*P[j].GravAccel[1] - SphP[j].HydroAccel[1]) * kernel->dp[2] * wk;
        out->NV_A[2][0] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dp[0] * wk;
        out->NV_A[2][1] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dp[1] * wk;
        out->NV_A[2][2] += (local->Accel[2] - All.cf_a2inv*P[j].GravAccel[2] - SphP[j].HydroAccel[2]) * kernel->dp[2] * wk;
        
        out->NV_D[0][0] += kernel->dv[0] * kernel->dp[0] * wk;
        out->NV_D[0][1] += kernel->dv[0] * kernel->dp[1] * wk;
        out->NV_D[0][2] += kernel->dv[0] * kernel->dp[2] * wk;
        out->NV_D[1][0] += kernel->dv[1] * kernel->dp[0] * wk;
        out->NV_D[1][1] += kernel->dv[1] * kernel->dp[1] * wk;
        out->NV_D[1][2] += kernel->dv[1] * kernel->dp[2] * wk;
        out->NV_D[2][0] += kernel->dv[2] * kernel->dp[0] * wk;
        out->NV_D[2][1] += kernel->dv[2] * kernel->dp[1] * wk;
        out->NV_D[2][2] += kernel->dv[2] * kernel->dp[2] * wk;
#endif
    
    } // Type = 0 check
}




