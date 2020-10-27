#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
/*
* This file was originally part of the GADGET3 code developed by Volker Springel.
* It has been updated significantly by PFH for basic compatibility with GIZMO,
* as well as code cleanups, and accommodating new GIZMO functionality for various
* other operations. See notes in subfind.c and GIZMO User Guide for details.
*/

#ifdef SUBFIND

#include "../fof.h"
#include "subfind.h"

static MyOutputFloat *R200, *M200;
struct Subfind_DensityOtherPropsEval_data_out *Subfind_DensityOtherPropsEval_DataResult, *Subfind_DensityOtherPropsEval_DataOut, *Subfind_DensityOtherPropsEval_GlobalPasser;


/*!
 The following functions are designed to flexibly allow adding/removing additional properties computed in groups of varying defintiions/radii/overdensities, etc.
 To add similar computations, follow their template (e.g. the "SUBFIND_ADDIO_..." options). This first subroutine is the core computation of the relevant properties
 in the group. compute them/add them into the "out" structure and the code here and other two scripts below should take care of the rest
*/
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here!!! -- */
/*! first define a short structure needed to pass in the group info here */
static struct Subfind_DensityOtherPropsEval_data_in {MyDouble Pos[3]; MyOutputFloat R200; int NodeList[NODELISTLENGTH]; /* all needed for any version */} *Subfind_DensityOtherPropsEval_DataIn, *Subfind_DensityOtherPropsEval_DataGet;
/*! now the main routine */
int Subfind_DensityOtherProps_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
    int ngb,n,j,k,startnode,listindex=0; double *subhalo_pos, Hsearch; struct Subfind_DensityOtherPropsEval_data_out out; memset(&out, 0, sizeof(struct Subfind_DensityOtherPropsEval_data_out));
    if(mode == 0) {subhalo_pos=Group[target].Pos; Hsearch=R200[target];} else {subhalo_pos=Subfind_DensityOtherPropsEval_DataGet[target].Pos; Hsearch=Subfind_DensityOtherPropsEval_DataGet[target].R200;}
    if(mode == 0) {startnode = All.MaxPart;} else {startnode = Subfind_DensityOtherPropsEval_DataGet[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;}
    while(startnode >= 0) {while(startnode >= 0) {
      ngb = ngb_treefind_variable_targeted(subhalo_pos, Hsearch, target, &startnode, mode, nexport, nsend_local, 63);
      if(ngb < 0) {return -2;}
      for(n = 0; n < ngb; n++)
        { j = Ngblist[n];
            double dp[3]; for(k=0;k<3;k++) {dp[k]=P[j].Pos[k]-subhalo_pos[k];} NEAREST_XYZ(dp[0],dp[1],dp[2],-1); double r2=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2]; if(r2>Hsearch*Hsearch) continue; /* position offset */
            out.M200+=P[j].Mass;
            
#ifdef SUBFIND_ADDIO_VELDISP
            for(k=0;k<3;k++) {double dv=P[j].Vel[k]/All.cf_atime+All.cf_hubble_a*All.cf_atime*dp[k]; out.V200[k]+=P[j].Mass*dv; out.Disp200+=P[j].Mass*dv*dv;}
#endif
#ifdef SUBFIND_ADDIO_BARYONS
            if(P[j].Type==0)
            {
                double temp_keV = 6.14e-16 * (SphP[j].InternalEnergy/UNIT_SPECEGY_IN_CGS); /* temp in keV, for fully-ionized primordial gas */
                out.gas_mass += P[j].Mass; out.temp += P[j].Mass * temp_keV;
                out.xlum += 1.52e-20 * (P[j].Mass*UNIT_MASS_IN_CGS) * (SphP[j].Density*All.cf_a3inv*UNIT_DENSITY_IN_CGS) * sqrt(temp_keV); /* converts to 1e44 erg/s assuming thermal brems for fully-ionized primordial composition */
            } else if (P[j].Type==4) {out.star_mass += P[j].Mass;}
#endif
            
        }
    } if(mode==1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = Subfind_DensityOtherPropsEval_DataGet[target].NodeList[listindex]; if(startnode >= 0) startnode = Nodes[startnode].u.d.nextnode;}}}
    if(mode == 0) {out2particle_Subfind_DensityOtherPropsEval(&out, target, 0);} else {Subfind_DensityOtherPropsEval_DataResult[target] = out;} /* collects the result at the right place */
  return 0;
}

/*! add data from other processors back to shared structure after computation. usually ASSIGN_ADD for normal 'set it equal if initial step, then add', but here in case you need to take maxes or mins or other complicated operations */
static inline void out2particle_Subfind_DensityOtherPropsEval(struct Subfind_DensityOtherPropsEval_data_out *out, int i, int mode)
{
    ASSIGN_ADD(Subfind_DensityOtherPropsEval_GlobalPasser[i].M200,out->M200,mode);
#ifdef SUBFIND_ADDIO_VELDISP
    {int k; for(k=0;k<3;k++) {ASSIGN_ADD(Subfind_DensityOtherPropsEval_GlobalPasser[i].V200[k],out->V200[k],mode);}}
    ASSIGN_ADD(Subfind_DensityOtherPropsEval_GlobalPasser[i].Disp200,out->Disp200,mode);
#endif
#ifdef SUBFIND_ADDIO_BARYONS
    ASSIGN_ADD(Subfind_DensityOtherPropsEval_GlobalPasser[i].gas_mass,out->gas_mass,mode); ASSIGN_ADD(Subfind_DensityOtherPropsEval_GlobalPasser[i].temp,out->temp,mode); ASSIGN_ADD(Subfind_DensityOtherPropsEval_GlobalPasser[i].xlum,out->xlum,mode);
    ASSIGN_ADD(Subfind_DensityOtherPropsEval_GlobalPasser[i].star_mass,out->star_mass,mode);
#endif
}

/*! do any final opertations on the data after it comes back from  Subfind_DensityOtherProps_evaluate  */
void Subfind_DensityOtherProps_finaloperations(struct Subfind_DensityOtherPropsEval_data_out *in)
{
    if(in->M200 <= 0) return; /* no members found */
#ifdef SUBFIND_ADDIO_VELDISP
    in->Disp200/=in->M200; int k; for(k=0;k<3;k++) {in->V200[k]/=in->M200; in->Disp200-=(in->V200[k])*(in->V200[k]);}
    in->Disp200=sqrt(in->Disp200/3.); /* convert to 1D velocity dispersion */
#endif
#ifdef SUBFIND_ADDIO_BARYONS
    if(in->gas_mass>0) {in->temp/=in->gas_mass;}
#endif
}

    
    
    

/*! ok everything below is the 'generic' part of this set of routines */



void Subfind_DensityOtherProps_Loop(void)
{
    long long ntot; int i, j, ndone, ndone_flag, npleft, dummy, rep, iter;
    MyFloat *Left, *Right; char *Todo; int ngrp, recvTask, place, nexport, nimport;
    double t0, t1, t2, t3, rguess, overdensity, Deltas[SUBFIND_ADDIO_NUMOVERDEN], z;

    /* allocate buffers to arrange communication */
    Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
    Left = (MyFloat *) mymalloc("Left", sizeof(MyFloat) * Ngroups);
    Right = (MyFloat *) mymalloc("Right", sizeof(MyFloat) * Ngroups);
    R200 = (MyOutputFloat *) mymalloc("R200", sizeof(MyOutputFloat) * Ngroups);
    M200 = (MyOutputFloat *) mymalloc("M200", sizeof(MyOutputFloat) * Ngroups);
    Subfind_DensityOtherPropsEval_GlobalPasser = (struct Subfind_DensityOtherPropsEval_data_out *) mymalloc("Subfind_DensityOtherPropsEval_GlobalPasser",Ngroups * sizeof(struct Subfind_DensityOtherPropsEval_data_out));
    Todo = (char *)mymalloc("Todo", sizeof(char) * Ngroups);
    size_t MyBufferSize = (size_t)All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct Subfind_DensityOtherPropsEval_data_in) + sizeof(struct Subfind_DensityOtherPropsEval_data_out) +
                     sizemax(sizeof(struct Subfind_DensityOtherPropsEval_data_in), sizeof(struct Subfind_DensityOtherPropsEval_data_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    if(All.ComovingIntegrationOn) {z = 1 / All.Time - 1;} else {z = 0;}
    double rhoback = 3 * All.OmegaMatter * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G), zplusone=1.+z;
    double omegaz = All.OmegaMatter * pow(zplusone,3) / (All.OmegaRadiation * pow(zplusone,4) + All.OmegaMatter * pow(zplusone,3) + (1 - All.OmegaMatter - All.OmegaLambda - All.OmegaRadiation) * pow(zplusone,2) + All.OmegaLambda);
    double x = omegaz - 1, DeltaTopHat = (18 * M_PI * M_PI + 82 * x - 39 * x * x) / omegaz;

    double Delta_ToEvalList[10];
    Delta_ToEvalList[0] = 200;             /* standard fixed overdensity with respect to background mean density */
    Delta_ToEvalList[1] = DeltaTopHat;     /* spherical tophat collapse-model overdensity with respect to background */
    Delta_ToEvalList[2] = 200/omegaz;      /* overdensity of 200 relative to critical, expressed relative to background density */
    Delta_ToEvalList[3] = 500/omegaz;      /* in this case use crit-500 with respect to background */
    Delta_ToEvalList[4] = 1000/omegaz;     /* in this case use crit-1000 with respect to background */
    Delta_ToEvalList[5] = 2500/omegaz;     /* in this case use crit-2500 with respect to background */
    Delta_ToEvalList[6] = 500;             /* in this case use mean-500 with respect to background */
    Delta_ToEvalList[7] = 1000;            /* in this case use mean-1000 with respect to background */
    Delta_ToEvalList[8] = 2500;            /* in this case use mean-2500 with respect to background */
    Delta_ToEvalList[9] = 5000;            /* in this case use mean-5000 with respect to background */

  for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {Deltas[j]=Delta_ToEvalList[j];} /* list we will use */

  for(rep = 0; rep < SUBFIND_ADDIO_NUMOVERDEN; rep++)	/* repeat for all three overdensity values */
    {
      t2 = my_second();
      for(i = 0; i < Ngroups; i++)
    {
      if(Group[i].Nsubs > 0)
        {rguess = pow(All.G * Group[i].Mass / (100 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits), 1.0 / 3); Right[i] = 3 * rguess; Left[i] = 0; Todo[i] = 1;} else {Todo[i] = 0;}
    }
      iter = 0;

      /* we will repeat the whole thing for those groups where we didn't converge to a SO radius yet */
      do
	{
	  t0 = my_second();
	  i = 0;		/* begin with this index */
	  do
	    {
	      for(j = 0; j < NTask; j++)
		{
		  Send_count[j] = 0;
		  Exportflag[j] = -1;
		}

	      /* do local particles and prepare export list */
	      for(nexport = 0; i < Ngroups; i++)
		{
		  if(Todo[i])
		    {
		      R200[i] = 0.5 * (Left[i] + Right[i]);
		      if(Subfind_RvirMvir_evaluate(i, 0, &nexport, Send_count) < 0)
			break;
		    }
		}
	      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
	      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
	      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
		{
		  nimport += Recv_count[j];
		  if(j > 0)
		    {
		      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		    }
		}
	      Subfind_DensityOtherPropsEval_DataGet = (struct Subfind_DensityOtherPropsEval_data_in *) mymalloc("	      Subfind_DensityOtherPropsEval_DataGet", nimport * sizeof(struct Subfind_DensityOtherPropsEval_data_in));
	      Subfind_DensityOtherPropsEval_DataIn = (struct Subfind_DensityOtherPropsEval_data_in *) mymalloc("	      Subfind_DensityOtherPropsEval_DataIn", nexport * sizeof(struct Subfind_DensityOtherPropsEval_data_in));

	      /* prepare particle data for export */
	      for(j = 0; j < nexport; j++)
		{
		  place = DataIndexTable[j].Index;

		  Subfind_DensityOtherPropsEval_DataIn[j].Pos[0] = Group[place].Pos[0];
		  Subfind_DensityOtherPropsEval_DataIn[j].Pos[1] = Group[place].Pos[1];
		  Subfind_DensityOtherPropsEval_DataIn[j].Pos[2] = Group[place].Pos[2];
		  Subfind_DensityOtherPropsEval_DataIn[j].R200 = R200[place];

		  memcpy(Subfind_DensityOtherPropsEval_DataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
		}

	      /* exchange data */
	      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
		{
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
			{
			  /* get the data */
			  MPI_Sendrecv(&Subfind_DensityOtherPropsEval_DataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_in), MPI_BYTE, recvTask, TAG_DENS_A,
				       &Subfind_DensityOtherPropsEval_DataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		    }
		}

	      myfree(Subfind_DensityOtherPropsEval_DataIn);
	      Subfind_DensityOtherPropsEval_DataResult = (struct Subfind_DensityOtherPropsEval_data_out *) mymalloc("	      Subfind_DensityOtherPropsEval_DataResult", nimport * sizeof(struct Subfind_DensityOtherPropsEval_data_out));
	      Subfind_DensityOtherPropsEval_DataOut = (struct Subfind_DensityOtherPropsEval_data_out *) mymalloc("	      Subfind_DensityOtherPropsEval_DataOut", nexport * sizeof(struct Subfind_DensityOtherPropsEval_data_out));

	      /* now do the locations that were sent to us */
	      for(j = 0; j < nimport; j++) Subfind_RvirMvir_evaluate(j, 1, &dummy, &dummy);

	      if(i >= Ngroups)
		ndone_flag = 1;
	      else
		ndone_flag = 0;

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
			  MPI_Sendrecv(&Subfind_DensityOtherPropsEval_DataResult[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_out),MPI_BYTE, recvTask, TAG_DENS_B,
				       &Subfind_DensityOtherPropsEval_DataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_out), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		    }
		}

	      /* add the result to the local particles */
	      for(j = 0; j < nexport; j++)
		{
		  place = DataIndexTable[j].Index;
		  M200[place] += Subfind_DensityOtherPropsEval_DataOut[j].M200;
		}

	      myfree(Subfind_DensityOtherPropsEval_DataOut);
	      myfree(Subfind_DensityOtherPropsEval_DataResult);
	      myfree(Subfind_DensityOtherPropsEval_DataGet);
	    }
	  while(ndone < NTask);


	  /* do final operations on results */
	  for(i = 0, npleft = 0; i < Ngroups; i++)
	    {
	      if(Todo[i])
		{
		  overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

		  if((Right[i] - Left[i]) > 1.0e-4 * Left[i])
		    {
		      /* need to redo this group */
		      npleft++;
                if(overdensity > Deltas[rep]) {Left[i] = R200[i];} else {Right[i] = R200[i];}
		      if(iter >= MAXITER - 10)
			{
			  printf
			    ("gr=%d task=%d  R200=%g Left=%g Right=%g Menclosed=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			     i, ThisTask, R200[i], Left[i], Right[i], M200[i], Right[i] - Left[i], Group[i].Pos[0], Group[i].Pos[1], Group[i].Pos[2]);
			  fflush(stdout);
			}
		    }
		  else
		    Todo[i] = 0;
		}
	    }

	  sumup_large_ints(1, &npleft, &ntot);

	  t1 = my_second();

	  if(ntot > 0)
	    {
	      iter++;

	      if(iter > 0 && ThisTask == 0)
		{
		  printf("SO iteration %d: need to repeat for %d%09d particles. (took %g sec)\n", iter,
			 (int) (ntot / 1000000000), (int) (ntot % 1000000000), timediff(t0, t1));
		  fflush(stdout);
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

      MPI_Barrier(MPI_COMM_WORLD);

      t3 = my_second();
      if(ThisTask == 0)
        {
          printf("SO Radius calculation took %g sec\n", timediff(t2, t3));
          fflush(stdout);
        }

      t0 = my_second();
      i = 0;			/* begin with this index */
      do
	{
	  for(j = 0; j < NTask; j++)
	    {
	      Send_count[j] = 0;
	      Exportflag[j] = -1;
	    }

	  /* do local particles and prepare export list */

	  for(nexport = 0; i < Ngroups; i++)
	    {
	      if(Subfind_DensityOtherProps_evaluate(i, 0, &nexport, Send_count) < 0)
		break;
	    }
	  qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

	  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	    {
	      nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

	  Subfind_DensityOtherPropsEval_DataGet = (struct Subfind_DensityOtherPropsEval_data_in *) mymalloc("	  Subfind_DensityOtherPropsEval_DataGet", nimport * sizeof(struct Subfind_DensityOtherPropsEval_data_in));
	  Subfind_DensityOtherPropsEval_DataIn = (struct Subfind_DensityOtherPropsEval_data_in *) mymalloc("	  Subfind_DensityOtherPropsEval_DataIn", nexport * sizeof(struct Subfind_DensityOtherPropsEval_data_in));

	  /* prepare particle data for export */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      Subfind_DensityOtherPropsEval_DataIn[j].Pos[0] = Group[place].Pos[0];
	      Subfind_DensityOtherPropsEval_DataIn[j].Pos[1] = Group[place].Pos[1];
	      Subfind_DensityOtherPropsEval_DataIn[j].Pos[2] = Group[place].Pos[2];
	      Subfind_DensityOtherPropsEval_DataIn[j].R200 = R200[place];

	      memcpy(Subfind_DensityOtherPropsEval_DataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	    }

	  /* exchange data */
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the data */
		      MPI_Sendrecv(&Subfind_DensityOtherPropsEval_DataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &Subfind_DensityOtherPropsEval_DataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  myfree(Subfind_DensityOtherPropsEval_DataIn);
	  Subfind_DensityOtherPropsEval_DataResult = (struct Subfind_DensityOtherPropsEval_data_out *) mymalloc("	  Subfind_DensityOtherPropsEval_DataResult", nimport * sizeof(struct Subfind_DensityOtherPropsEval_data_out));
	  Subfind_DensityOtherPropsEval_DataOut = (struct Subfind_DensityOtherPropsEval_data_out *) mymalloc("	  Subfind_DensityOtherPropsEval_DataOut", nexport * sizeof(struct Subfind_DensityOtherPropsEval_data_out));

	  /* now do the locations that were sent to us */
	  for(j = 0; j < nimport; j++)
	    Subfind_DensityOtherProps_evaluate(j, 1, &dummy, &dummy);

	  if(i >= Ngroups)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

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
		      MPI_Sendrecv(&Subfind_DensityOtherPropsEval_DataResult[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_out), MPI_BYTE, recvTask, TAG_DENS_B,
				   &Subfind_DensityOtherPropsEval_DataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct Subfind_DensityOtherPropsEval_data_out), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  /* add the result to the local particles */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
            out2particle_Subfind_DensityOtherPropsEval(Subfind_DensityOtherPropsEval_DataOut, place, 1);
	    }

	  myfree(Subfind_DensityOtherPropsEval_DataOut);
	  myfree(Subfind_DensityOtherPropsEval_DataResult);
	  myfree(Subfind_DensityOtherPropsEval_DataGet);
	}
      while(ndone < NTask);

      MPI_Barrier(MPI_COMM_WORLD);

      t1 = my_second();

      if(ThisTask == 0)
	{
	  printf("secondary sunfind loop for additional information took %g sec\n", timediff(t0, t1));
	  fflush(stdout);
	}


      t0 = my_second();
      for(i = 0; i < Ngroups; i++)
	{
	  if(Group[i].Nsubs > 0)
	    {
	      overdensity = M200[i] / (4.0 * M_PI / 3.0 * R200[i] * R200[i] * R200[i]) / rhoback;

	      if((overdensity - Deltas[rep]) > 0.1 * Deltas[rep])
		{
		  R200[i] = M200[i] = 0;
		}
	      else if(M200[i] < 5 * Group[i].Mass / Group[i].Len)
		{
		  R200[i] = M200[i] = 0;
		}
	    }
	  else
	    R200[i] = M200[i] = 0;

        Subfind_DensityOtherPropsEval_GlobalPasser[i].R200 = R200[i];
        Subfind_DensityOtherProps_finaloperations(&Subfind_DensityOtherPropsEval_GlobalPasser[i]); /* any other final operations before saving */
        memcpy(&Group[i].SubHaloProps_vsDelta[rep], &Subfind_DensityOtherPropsEval_GlobalPasser[i], sizeof(struct Subfind_DensityOtherPropsEval_data_out)); /* save it in memory */

    }
    }
  t1 = my_second();
  if(ThisTask == 0) {printf("Saving data took %g sec\n", timediff(t0, t1)); fflush(stdout);}

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Todo);
  myfree(Subfind_DensityOtherPropsEval_GlobalPasser);
  myfree(M200);
  myfree(R200);
  myfree(Right);
  myfree(Left);

  myfree(Ngblist);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int Subfind_RvirMvir_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, listindex = 0;
  double h, mass, massret;
  MyDouble *pos;


  if(mode == 0)
    {
      pos = Group[target].Pos;
      h = R200[target];
    }
  else
    {
      pos = Subfind_DensityOtherPropsEval_DataGet[target].Pos;
      h = Subfind_DensityOtherPropsEval_DataGet[target].R200;
    }

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = Subfind_DensityOtherPropsEval_DataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  mass = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  massret = subfind_ovderdens_treefind(pos, h, target, &startnode, mode, nexport, nsend_local);
	  if(massret < 0) {return -2;}
	  mass += massret;
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = Subfind_DensityOtherPropsEval_DataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    M200[target] = mass;
  else
    Subfind_DensityOtherPropsEval_DataResult[target].M200 = mass;

  return 0;
}


double subfind_ovderdens_treefind(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				  int mode, int *nexport, int *nsend_local)
{
  int no, p, task, nexport_save;
  struct NODE *current;
  double mass;
  MyDouble dx, dy, dz, dist, r2;

#define FACT2 0.86602540
  nexport_save = *nexport;

  mass = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

        dist = hsml; double xtmp; xtmp=0;
      dx = NGB_PERIODIC_BOX_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
	  if(dx > dist)
	    continue;
      dy = NGB_PERIODIC_BOX_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
	  if(dy > dist)
	    continue;
      dz = NGB_PERIODIC_BOX_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2], -1);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  mass += P[p].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(mode == 0)
		{
		  if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= All.BunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0) {endrun(13005);}	/* in this case, the buffer is too small to process even a single particle */
			  for(task = 0; task < NTask; task++) {nsend_local[task] = 0;}
			  for(no = 0; no < nexport_save; no++) {nsend_local[DataIndexTable[no].Task]++;}
			  return -1; /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */
			}
		      Exportnodecount[task] = 0;
		      Exportindex[task] = *nexport;
		      DataIndexTable[*nexport].Task = task;
		      DataIndexTable[*nexport].Index = target;
		      DataIndexTable[*nexport].IndexGet = *nexport;
		      *nexport = *nexport + 1;
		      nsend_local[task]++;
		    }

		  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		    DomainNodeIndex[no - (All.MaxPart + MaxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return mass;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */
        dist = hsml + 0.5 * current->len; double xtmp; xtmp=0;
      dx = NGB_PERIODIC_BOX_LONG_X(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
	  if(dx > dist)
	    continue;
      dy = NGB_PERIODIC_BOX_LONG_Y(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
	  if(dy > dist)
	    continue;
      dz = NGB_PERIODIC_BOX_LONG_Z(current->center[0] - searchcenter[0], current->center[1] - searchcenter[1], current->center[2] - searchcenter[2], -1);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
	    continue;

	  if((current->u.d.bitflags & ((1 << BITFLAG_TOPLEVEL) + (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))) == 0)	/* only use fully local nodes */
	    {
	      /* test whether the node is contained within the sphere */
	      dist = hsml - FACT2 * current->len;
	      if(dist > 0)
		if(r2 < dist * dist)
		  {
		    mass += current->u.d.mass;
		    continue;
		  }
	    }

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return mass;
}

#endif
