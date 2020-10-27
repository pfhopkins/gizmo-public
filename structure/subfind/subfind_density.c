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


/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct densdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *DensDataIn, *DensDataGet;


static struct densdata_out
{
  MyFloat Rho;
  MyFloat VelDisp, Vx, Vy, Vz;
  int Ngb;
}
 *DensDataResult, *DensDataOut;


static MyFloat *DM_Vx, *DM_Vy, *DM_Vz;
static long long Ntotal;

void subfind_density(int j_in)
{
  long long ntot;
  int i, j, ndone, ndone_flag, npleft, dummy, iter = 0;
  MyFloat *Left, *Right;
  char *Todo;
  int ngrp, recvTask, place, nexport, nimport;
  double vel_to_phys, atime, t0, t1;

#ifdef FOF_DENSITY_SPLIT_TYPES
  int j_target = j_in;
  if(j_in < 0) j_target = -j_in - 1;
#endif

  if(ThisTask == 0)
    {
      printf("finding densities, mode %d \n", j_in);
      fflush(stdout);
    }

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));
  Dist2list = (double *) mymalloc("Dist2list", NumPart * sizeof(double));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct densdata_in) + sizeof(struct densdata_out) +
					     sizemax(sizeof(struct densdata_in),
						     sizeof(struct densdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  Left = (MyFloat *)mymalloc("Left", sizeof(MyFloat) * NumPart);
  Right = (MyFloat *)mymalloc("Right", sizeof(MyFloat) * NumPart);
  Todo = (char *)mymalloc("Todo", sizeof(char) * NumPart);

  DM_Vx = (MyFloat *)mymalloc("DM_Vx", sizeof(MyFloat) * NumPart);
  DM_Vy = (MyFloat *)mymalloc("DM_Vy", sizeof(MyFloat) * NumPart);
  DM_Vz = (MyFloat *)mymalloc("DM_Vz", sizeof(MyFloat) * NumPart);

  for(i = 0; i < NumPart; i++)
    {
      Left[i] = Right[i] = 0;
      P[i].DM_NumNgb = 0;
      Todo[i] = 1;
    }

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
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

	  for(nexport = 0; i < NumPart; i++)
	    {
	      if(Todo[i])
#ifdef FOF_DENSITY_SPLIT_TYPES
		if(P[i].Type == j_target)
#else
        if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
#endif
		  {
		    if(subfind_density_evaluate(i, 0, &nexport, Send_count, j_in) < 0)
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

	  DensDataGet =
	    (struct densdata_in *) mymalloc("	  DensDataGet", nimport * sizeof(struct densdata_in));
	  DensDataIn =
	    (struct densdata_in *) mymalloc("	  DensDataIn", nexport * sizeof(struct densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      DensDataIn[j].Pos[0] = P[place].Pos[0];
	      DensDataIn[j].Pos[1] = P[place].Pos[1];
	      DensDataIn[j].Pos[2] = P[place].Pos[2];
	      DensDataIn[j].Hsml = P[place].DM_Hsml;

	      memcpy(DensDataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
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
		      MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  myfree(DensDataIn);
	  DensDataResult =
	    (struct densdata_out *) mymalloc("	  DensDataResult", nimport * sizeof(struct densdata_out));
	  DensDataOut =
	    (struct densdata_out *) mymalloc("	  DensDataOut", nexport * sizeof(struct densdata_out));


	  /* now do the particles that were sent to us */
	  for(j = 0; j < nimport; j++)
	    subfind_density_evaluate(j, 1, &dummy, &dummy, j_in);

        if(i >= NumPart)
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
		      MPI_Sendrecv(&DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  /* add the result to the local particles */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      P[place].DM_NumNgb += DensDataOut[j].Ngb;
	      P[place].u.DM_Density += DensDataOut[j].Rho;

            if(j_in >= 0)	/* sum up VelDisp only for own species */
                P[place].v.DM_VelDisp += DensDataOut[j].VelDisp;
            
	      DM_Vx[place] += DensDataOut[j].Vx;
	      DM_Vy[place] += DensDataOut[j].Vy;
	      DM_Vz[place] += DensDataOut[j].Vz;
	    }

	  myfree(DensDataOut);
	  myfree(DensDataResult);
	  myfree(DensDataGet);
	}
      while(ndone < NTask);

      if(ThisTask == 0) printf("checking smoothing lenghs ... \n");

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
	{
	  /* now check whether we had enough neighbours */
#ifdef FOF_DENSITY_SPLIT_TYPES
	  if(j_in < 0) {Todo[i] = 0;}
	  if((P[i].Type != j_target) || (j_in < 0)) continue;
#else
	  if(!((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES))) continue;
#endif

	  if(Todo[i])
	    {
	      if(P[i].DM_NumNgb != All.DesLinkNgb &&
		 ((Right[i] - Left[i]) > 1.0e-4 * Left[i] || Left[i] == 0 || Right[i] == 0))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(P[i].DM_NumNgb < All.DesLinkNgb)
		    Left[i] = DMAX(P[i].DM_Hsml, Left[i]);
		  else
		    {
		      if(Right[i] != 0)
			{
			  if(P[i].DM_Hsml < Right[i])
			    Right[i] = P[i].DM_Hsml;
			}
		      else
			Right[i] = P[i].DM_Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%llu Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (unsigned long long) P[i].ID, P[i].DM_Hsml, Left[i], Right[i],
			 (double) P[i].DM_NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(Right[i] > 0 && Left[i] > 0)
		    P[i].DM_Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		  else
		    {
		      if(Right[i] == 0 && Left[i] == 0)
			endrun(8187);	/* can't occur */

		      if(Right[i] == 0 && Left[i] > 0)
			P[i].DM_Hsml *= 1.26;

		      if(Right[i] > 0 && Left[i] == 0)
			P[i].DM_Hsml /= 1.26;
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
	      printf("ngb iteration %d: need to repeat for %d%09d particles. (took %g sec)\n", iter,
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

  if(All.ComovingIntegrationOn)
    atime = All.Time;
  else
    atime = 1;

  vel_to_phys = 1.0 / atime;

  if(ThisTask == 0) printf("final operations ... \n");

  for(i = 0; i < NumPart; i++)
#ifdef FOF_DENSITY_SPLIT_TYPES
    if(P[i].Type == j_target && j_in >= 0)
#else
    if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
#endif
      {
	if(P[i].DM_NumNgb > 0)
    {
	DM_Vx[i] /= P[i].DM_NumNgb;
	DM_Vy[i] /= P[i].DM_NumNgb;
	DM_Vz[i] /= P[i].DM_NumNgb;
	P[i].v.DM_VelDisp /= P[i].DM_NumNgb;
	P[i].v.DM_VelDisp = vel_to_phys * sqrt(P[i].v.DM_VelDisp - DM_Vx[i] * DM_Vx[i] - DM_Vy[i] * DM_Vy[i] - DM_Vz[i] * DM_Vz[i]);
    }
      }

  myfree(DM_Vz);
  myfree(DM_Vy);
  myfree(DM_Vx);

  myfree(Todo);
  myfree(Right);
  myfree(Left);

  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Dist2list);
  myfree(Ngblist);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here!!! -- */
int subfind_density_evaluate(int target, int mode, int *nexport, int *nsend_local, int tp)
{
  int j, n;
  int startnode, numngb, ngb, listindex = 0;
  double hmax;
  double h, h2, hinv, hinv3, hinv4;
  double rho, wk, dwk;
  double r, r2, u, mass_j, v2, vx, vy, vz;
  MyDouble *pos;

  rho = 0;
  numngb = 0;
  v2 = vx = vy = vz = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = P[target].DM_Hsml;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      h = DensDataGet[target].Hsml;
    }

  if(h == 0) endrun(5622);

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
#ifdef FOF_DENSITY_SPLIT_TYPES
	  if(tp < 0)
	    {ngb = subfind_ngb_treefind_linkpairs(pos, h, target, &startnode, mode, &hmax, nexport, nsend_local);}
	  else
	    {ngb = subfind_ngb_treefind_linkngb(pos, h, target, &startnode, mode, &hmax, nexport, nsend_local);}
#else
	  ngb = subfind_ngb_treefind_linkngb(pos, h, target, &startnode, mode, &hmax, nexport, nsend_local);
#endif
	  if(ngb < 0) {return -2;}

	  if(tp >= 0)
	    if(mode == 0 && hmax > 0)
	      {
              P[target].DM_Hsml = hmax;
              h = hmax;
              if(ngb != All.DesLinkNgb) {endrun(121);}
	      }

	  numngb += ngb;

	  for(n = 0; n < ngb; n++)
	    {
	      j = Ngblist[n];
	      r2 = Dist2list[n];

#ifdef FOF_DENSITY_SPLIT_TYPES
	      if(tp < 0) {h = P[j].DM_Hsml;}
#endif

	      h2 = h * h;
	      if(r2 < h2)
		{
		  r = sqrt(r2);
		  kernel_hinv(h, &hinv, &hinv3, &hinv4);
		  u = r * hinv;
		  kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
		  mass_j = P[j].Mass;
		  rho += (mass_j * wk);
		}

	      vx += P[j].Vel[0];
	      vy += P[j].Vel[1];
	      vz += P[j].Vel[2];

	      v2 += P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2];
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
    {
      P[target].DM_NumNgb = numngb;
      P[target].u.DM_Density = rho;
      if(tp >= 0) {P[target].v.DM_VelDisp = v2;}		/* sum up VelDisp only for own species */
          
      DM_Vx[target] = vx;
      DM_Vy[target] = vy;
      DM_Vz[target] = vz;
    }
  else
    {
      DensDataResult[target].Ngb = numngb;
      DensDataResult[target].Rho = rho;
      DensDataResult[target].VelDisp = v2;
      DensDataResult[target].Vx = vx;
      DensDataResult[target].Vy = vy;
      DensDataResult[target].Vz = vz;
    }

  return 0;
}


void subfind_setup_smoothinglengths(int j)
{
  int i, no, p;

  for(i = 0; i < NumPart; i++)
    {
#ifdef FOF_DENSITY_SPLIT_TYPES
      if(P[i].Type == j)
#else
      if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)) || ((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES)))
#endif
	{
	  no = Father[i];

	  /* Not a good guess for gas/stars component, need more thought ! */
	  while(10 * All.DesLinkNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }
#ifdef FOF_DENSITY_SPLIT_TYPES
	  if(P[i].Type == 0) {P[i].DM_Hsml = PPP[i].Hsml * pow( 1.*All.DesLinkNgb / All.DesNumNgb, 1./3.);}
	  else {P[i].DM_Hsml = pow(3.0 / (4 * M_PI) * All.DesLinkNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;}
#else
	  P[i].DM_Hsml = pow(3.0 / (4 * M_PI) * All.DesLinkNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#endif

	}
    }
}


static int Nhsml;

static struct hsml_data
{
  float Hsml;
  float Density;
  float VelDisp;
  MyIDType ID;
}
 *Hsml_list;

int subfind_compare_hsml_data(const void *a, const void *b)
{
  if(((struct hsml_data *) a)->ID < ((struct hsml_data *) b)->ID) {return -1;}
  if(((struct hsml_data *) a)->ID > ((struct hsml_data *) b)->ID) {return +1;}
  return 0;
}


void subfind_save_densities(int num)
{
  int i, nprocgroup, primaryTask, groupTask;
  char buf[1000];
  double t0, t1;

  if(ThisTask == 0)
    {
      printf("start saving smoothing lengths and densities\n");
      fflush(stdout);
    }

  for(i = 0, Nhsml = 0; i < NumPart; i++)
#ifdef FOF_DENSITY_SPLIT_TYPES
    if(((1 << P[i].Type) & (FOF_DENSITY_SPLIT_TYPES)))
#else
    if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
#endif
      Nhsml++;

  MPI_Allgather(&Nhsml, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 1, Send_offset[0] = 0; i < NTask; i++)
    Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];

  sumup_large_ints(1, &Nhsml, &Ntotal);

  Hsml_list = (struct hsml_data *)mymalloc("Hsml_list", Nhsml * sizeof(struct hsml_data));

  for(i = 0, Nhsml = 0; i < NumPart; i++)
#ifdef FOF_DENSITY_SPLIT_TYPES
    if(((1 << P[i].Type) & (FOF_DENSITY_SPLIT_TYPES)))
#else
    if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
#endif
      {
	Hsml_list[Nhsml].Hsml = P[i].DM_Hsml;
	Hsml_list[Nhsml].Density = P[i].u.DM_Density;
	Hsml_list[Nhsml].VelDisp = P[i].v.DM_VelDisp;
	Hsml_list[Nhsml].ID = P[i].ID;
	Nhsml++;
      }

  t0 = my_second();
  parallel_sort(Hsml_list, Nhsml, sizeof(struct hsml_data), subfind_compare_hsml_data);
  t1 = my_second();

  if(ThisTask == 0)
    {
      printf("Sorting of densities in ID sequence took = %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/hsmldir_%03d", All.OutputDir, num);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);

  if(NTask < All.NumFilesWrittenInParallel)
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(241931);
    }

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  primaryTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (primaryTask + groupTask))	/* ok, it's this processor's turn */
	subfind_save_local_densities(num);
      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }

  myfree(Hsml_list);

}

void subfind_save_local_densities(int num)
{
  char fname[1000];
  int i;
  float *tmp;
  FILE *fd;


  sprintf(fname, "%s/hsmldir_%03d/%s_%03d.%d", All.OutputDir, num, "hsml", num, ThisTask);
  if(!(fd = fopen(fname, "w")))
    {
      printf("can't open file `%s`\n", fname);
      endrun(1183);
    }

  my_fwrite(&Nhsml, sizeof(int), 1, fd);
  my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);	/* this is the number of IDs in previous files */
  my_fwrite(&Ntotal, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);

  tmp = (float *)mymalloc("tmp", Nhsml * sizeof(float));

  for(i = 0; i < Nhsml; i++)
    tmp[i] = Hsml_list[i].Hsml;
  my_fwrite(tmp, sizeof(float), Nhsml, fd);

  for(i = 0; i < Nhsml; i++)
    tmp[i] = Hsml_list[i].Density;
  my_fwrite(tmp, sizeof(float), Nhsml, fd);

  for(i = 0; i < Nhsml; i++)
    tmp[i] = Hsml_list[i].VelDisp;
  my_fwrite(tmp, sizeof(float), Nhsml, fd);

  myfree(tmp);

  fclose(fd);
}


#endif
