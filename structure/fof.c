#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>
#include "../allvars.h"
#include "../proto.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*! \file fof.c
 *  \brief parallel FoF group finder
 */
/*
* This file was originally part of the GADGET3 code developed by Volker Springel.
* It has been updated significantly by PFH for basic compatibility with GIZMO,
* as well as code cleanups, and accommodating new GIZMO functionality for various
* other operations. 
*/


#ifdef FOF
#include "fof.h"
#ifdef SUBFIND
#include "subfind/subfind.h"
#endif


int Ngroups, TotNgroups;
long long TotNids;

group_properties *Group;



static struct fofdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyIDType MinID;
  MyIDType MinIDTask;
  int NodeList[NODELISTLENGTH];
}
 *FoFDataIn, *FoFDataGet;

static struct fofdata_out
{
  MyFloat Distance;
  MyIDType MinID;
  MyIDType MinIDTask;
}
 *FoFDataResult, *FoFDataOut;


static struct fof_particle_list
{
  MyIDType MinID;
  MyIDType MinIDTask;
  int Pindex;
}
 *FOF_PList;

static fof_group_list *FOF_GList;

static fof_id_list *ID_list;


static double LinkL;
static int NgroupsExt, Nids;

static int MyFOF_PRIMARY_LINK_TYPES;
static int MyFOF_SECONDARY_LINK_TYPES;
static int MyFOF_GROUP_MIN_SIZE;

static MyIDType *Head, *Len, *Next, *Tail, *MinID, *MinIDTask;
static char *NonlocalFlag;


static float *fof_nearest_distance;
static float *fof_nearest_hsml;


void fof_fof(int num)
{
  int i, ndm, start, lenloc, largestgroup, n;
  double mass, masstot, rhodm, t0, t1;
  struct unbind_data *d;
  long long ndmtot;

#ifdef IO_SUBFIND_READFOF_FROMIC
  read_fof(num);
#endif

  MyFOF_PRIMARY_LINK_TYPES = FOF_PRIMARY_LINK_TYPES;
  MyFOF_SECONDARY_LINK_TYPES = FOF_SECONDARY_LINK_TYPES;
  MyFOF_GROUP_MIN_SIZE = FOF_GROUP_MIN_SIZE;

  if(ThisTask == 0)
    {
      printf("\nBegin to compute FoF group catalogues...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
    }

  CPU_Step[CPU_MISC] += measure_time();

  domain_Decomposition(1, 0, 0);

  force_treefree();

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
      {
        ndm++;
        mass += P[i].Mass;
      }
  sumup_large_ints(1, &ndm, &ndmtot);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(All.TotN_gas)
    {rhodm = (All.OmegaMatter - All.OmegaBaryon) * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);}
  else
    {rhodm = All.OmegaMatter * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);}

  LinkL = LINKLENGTH * pow(masstot / ndmtot / rhodm, 1.0 / 3);

    PRINT_STATUS("Comoving linking length: %g : (presently allocated=%g MB) ",LinkL,AllocatedBytes / (1024.0 * 1024.0))

  FOF_PList =
    (struct fof_particle_list *) mymalloc("FOF_PList", NumPart *
					  sizemax(sizeof(struct fof_particle_list), 3 * sizeof(MyIDType)));

  MinID = (MyIDType *) FOF_PList;
  MinIDTask = MinID + NumPart;
  Head = MinIDTask + NumPart;
  Len = (MyIDType *) mymalloc("Len", NumPart * sizeof(MyIDType));
  Next = (MyIDType *) mymalloc("Next", NumPart * sizeof(MyIDType));
  Tail = (MyIDType *) mymalloc("Tail", NumPart * sizeof(MyIDType));

  CPU_Step[CPU_FOF] += measure_time();

  if(ThisTask == 0)
    printf("Tree construction.\n");

  /* build index list of particles of selected primary species */
  d = (struct unbind_data *) mymalloc("d", NumPart * sizeof(struct unbind_data));
  for(i = 0, n = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
      d[n++].index = i;

  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
    
  force_treebuild(n, d);

  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i] = 1;
      Next[i] = -1;
      MinID[i] = P[i].ID;
      MinIDTask[i] = ThisTask;
    }


  t0 = my_second();

  fof_find_groups();

  t1 = my_second();
  if(ThisTask == 0)
    printf("group finding took = %g sec\n", timediff(t0, t1));


  t0 = my_second();

  fof_find_nearest_dmparticle();

  t1 = my_second();
  if(ThisTask == 0)
    printf("attaching gas and star particles to nearest dm particles took = %g sec\n", timediff(t0, t1));


  t0 = my_second();

  for(i = 0; i < NumPart; i++)
    {
      Next[i] = MinID[Head[i]];
      Tail[i] = MinIDTask[Head[i]];

      if(Tail[i] >= (MyIDType)NTask)	/* it appears that the Intel C 9.1 on Itanium2 produces incorrect code if
				   this if-statemet is omitted. Apparently, the compiler then joins the two loops,
				   but this is here not permitted because storage for FOF_PList actually overlaps
				   (on purpose) with MinID/MinIDTask/Head */
	{
	  printf("oh no: ThisTask=%d i=%d Head[i]=%d  NumPart=%d MinIDTask[Head[i]]=%d\n",
		 ThisTask, i, (int) Head[i], NumPart, (int) MinIDTask[Head[i]]);
	  fflush(stdout);
	  endrun(8812);
	}
    }

  for(i = 0; i < NumPart; i++)
    {
      FOF_PList[i].MinID = Next[i];
      FOF_PList[i].MinIDTask = Tail[i];
      FOF_PList[i].Pindex = i;
    }

  force_treefree();

  myfree(d);
  
  myfree(Tail);
  myfree(Next);
  myfree(Len);

  FOF_GList = (fof_group_list *) mymalloc("FOF_GList", sizeof(fof_group_list) * NumPart);

  fof_compile_catalogue();

  t1 = my_second();
  if(ThisTask == 0)
    printf("compiling local group data and catalogue took = %g sec\n", timediff(t0, t1));


  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &Nids, &TotNids);

  if(TotNgroups > 0)
    {
      int largestloc = 0;

      for(i = 0; i < NgroupsExt; i++)
	if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount > largestloc)
	  largestloc = FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
  else
    largestgroup = 0;

  if(ThisTask == 0)
    {
      printf("\nTotal number of groups with at least %d particles: %d\n", MyFOF_GROUP_MIN_SIZE, TotNgroups);
      if(TotNgroups > 0)
	{
	  printf("Largest group has %d particles.\n", largestgroup);
	  printf("Total number of particles in groups: %d%09d\n\n",
		 (int) (TotNids / 1000000000), (int) (TotNids % 1000000000));
	}
    }

  t0 = my_second();

  Group = (group_properties *) mymalloc("Group", sizeof(group_properties) *
					 IMAX(NgroupsExt, TotNgroups / NTask + 1));

  PRINT_STATUS("group properties are now allocated.. (presently allocated=%g MB)",AllocatedBytes / (1024.0 * 1024.0));
    
  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
	{
	  start++;
	  if(start > NumPart)
	    endrun(78);
	}

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
	endrun(123);

      for(lenloc = 0; start + lenloc < NumPart;)
	if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
	  lenloc++;
	else
	  break;

      Group[i].MinID = FOF_GList[i].MinID;
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_group_properties(i, start, lenloc);

      start += lenloc;
    }

  fof_exchange_group_data();

  fof_finish_group_properties();

  t1 = my_second();
  if(ThisTask == 0)
    printf("computation of group properties took = %g sec\n", timediff(t0, t1));

#ifdef BH_SEED_FROM_FOF
  if(num < 0){   // Make BHs in every call to fof_fof (including the group finding for each snapshot)
      if(All.Time < 1.0/(1.0+All.SeedBlackHoleMinRedshift)) { fof_make_black_holes(); } else {  printf("skipping black hole seeding at a = %g \n", All.Time); }
  }
#endif

#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==1)
  if(num < 0)
    fof_assign_HostHaloMass();
#endif
#endif

  CPU_Step[CPU_FOF] += measure_time();

  if(num >= 0)
    {
      fof_save_groups(num);
#ifdef SUBFIND
      if(DumpFlag != 2)
	subfind(num);
#endif
    }

  myfree(Group);

  myfree(FOF_GList);
  myfree(FOF_PList);

  PRINT_STATUS("Finished computing FoF groups.  (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));

  CPU_Step[CPU_FOF] += measure_time();

#ifdef SUBFIND
  domain_Decomposition(1, 0, 0);
#else
  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
#endif

  if(ThisTask == 0)
    printf("Tree construction.\n");
  force_treebuild(NumPart, NULL);

  TreeReconstructFlag = 0;
}



void fof_find_groups(void)
{
  int i, j, ndone_flag, link_count, dummy, nprocessed;
  int ndone, ngrp, recvTask, place, nexport, nimport, link_across;
  int npart, marked;
  long long totmarked, totnpart;
  long long link_across_tot, ntot;
  MyIDType *MinIDOld;
  char *FoFDataOut, *FoFDataResult, *MarkedFlag, *ChangedFlag;
  double t0, t1;

  PRINT_STATUS("Start linking particles (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     2 * sizeof(struct fofdata_in)));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  NonlocalFlag = (char *) mymalloc("NonlocalFlag", NumPart * sizeof(char));
  MarkedFlag = (char *) mymalloc("MarkedFlag", NumPart * sizeof(char));
  ChangedFlag = (char *) mymalloc("ChangedFlag", NumPart * sizeof(char));
  MinIDOld = (MyIDType *) mymalloc("MinIDOld", NumPart * sizeof(MyIDType));

  t0 = my_second();

  /* first, link only among local particles */
  for(i = 0, marked = 0, npart = 0; i < NumPart; i++)
    {
      if(((1 << P[i].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
	{
	  fof_find_dmparticles_evaluate(i, -1, &dummy, &dummy);

	  npart++;

	  if(NonlocalFlag[i])
	    marked++;
	}
    }


  sumup_large_ints(1, &marked, &totmarked);
  sumup_large_ints(1, &npart, &totnpart);

  t1 = my_second();


  PRINT_STATUS("links on local processor done (took %g sec).\nMarked=%d%09d out of the %d%09d primaries which are linked",timediff(t0, t1),(int) (totmarked / 1000000000), (int) (totmarked % 1000000000),(int) (totnpart / 1000000000), (int) (totnpart % 1000000000));
  PRINT_STATUS("\nlinking across processors (presently allocated=%g MB)",AllocatedBytes / (1024.0 * 1024.0));
    
  for(i = 0; i < NumPart; i++)
    {
      MinIDOld[i] = MinID[Head[i]];
      MarkedFlag[i] = 1;
    }

  do
    {
      t0 = my_second();

      for(i = 0; i < NumPart; i++)
	{
	  ChangedFlag[i] = MarkedFlag[i];
	  MarkedFlag[i] = 0;
	}

      i = 0;			/* begin with this index */
      link_across = 0;
      nprocessed = 0;

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
	      if(((1 << P[i].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
		{
		  if(NonlocalFlag[i] && ChangedFlag[i])
		    {
		      if(fof_find_dmparticles_evaluate(i, 0, &nexport, Send_count) < 0)
			break;

		      nprocessed++;
		    }
		}
	    }

	  MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

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

	  FoFDataGet = (struct fofdata_in *) mymalloc("FoFDataGet", nimport * sizeof(struct fofdata_in));
	  FoFDataIn = (struct fofdata_in *) mymalloc("FoFDataIn", nexport * sizeof(struct fofdata_in));


	  /* prepare particle data for export */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      FoFDataIn[j].Pos[0] = P[place].Pos[0];
	      FoFDataIn[j].Pos[1] = P[place].Pos[1];
	      FoFDataIn[j].Pos[2] = P[place].Pos[2];
	      FoFDataIn[j].MinID = MinID[Head[place]];
	      FoFDataIn[j].MinIDTask = MinIDTask[Head[place]];

	      memcpy(FoFDataIn[j].NodeList,
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
		      MPI_Sendrecv(&FoFDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
				   recvTask, TAG_FOF_A,
				   &FoFDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
				   recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  myfree(FoFDataIn);
	  FoFDataResult = (char *) mymalloc("FoFDataResult", nimport * sizeof(char));
	  FoFDataOut = (char *) mymalloc("FoFDataOut", nexport * sizeof(char));

	  /* now do the particles that were sent to us */

	  for(j = 0; j < nimport; j++)
	    {
	      link_count = fof_find_dmparticles_evaluate(j, 1, &dummy, &dummy);
	      link_across += link_count;
	      if(link_count)
		FoFDataResult[j] = 1;
	      else
		FoFDataResult[j] = 0;
	    }

	  /* exchange data */
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&FoFDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(char),
				   MPI_BYTE, recvTask, TAG_FOF_B,
				   &FoFDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(char),
				   MPI_BYTE, recvTask, TAG_FOF_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  /* need to mark the particle if it induced a link */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
	      if(FoFDataOut[j])
		MarkedFlag[place] = 1;
	    }

	  myfree(FoFDataOut);
	  myfree(FoFDataResult);
	  myfree(FoFDataGet);

	  if(i >= NumPart)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
      while(ndone < NTask);


      sumup_large_ints(1, &link_across, &link_across_tot);
      sumup_large_ints(1, &nprocessed, &ntot);

      t1 = my_second();

	  PRINT_STATUS("have done %d%09d cross links (processed %d%09d, took %g sec)",(int) (link_across_tot / 1000000000), (int) (link_across_tot % 1000000000),(int) (ntot / 1000000000), (int) (ntot % 1000000000), timediff(t0, t1));

      /* let's check out which particles have changed their MinID */
      for(i = 0; i < NumPart; i++)
	if(NonlocalFlag[i])
	  {
	    if(MinID[Head[i]] != MinIDOld[i])
	      MarkedFlag[i] = 1;

	    MinIDOld[i] = MinID[Head[i]];
	  }

    }
  while(link_across_tot > 0);

  myfree(MinIDOld);
  myfree(ChangedFlag);
  myfree(MarkedFlag);
  myfree(NonlocalFlag);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

PRINT_STATUS("Local groups found.");
}


/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here!!! -- */
int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, links, p, s, ss, listindex = 0;
  int startnode, numngb_inbox;
  MyDouble *pos;

  links = 0;

  if(mode == 0 || mode == -1)
    pos = P[target].Pos;
  else
    pos = FoFDataGet[target].Pos;

  if(mode == 0 || mode == -1)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = FoFDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  if(mode == -1) {*nexport = 0;}

	  numngb_inbox = ngb_treefind_fof_primary(pos, LinkL, target, &startnode, mode, nexport, nsend_local, MyFOF_PRIMARY_LINK_TYPES);
	  if(numngb_inbox < 0) {return -2;}

	  if(mode == -1) {if(*nexport == 0) {NonlocalFlag[target] = 0;} else {NonlocalFlag[target] = 1;}}

	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = Ngblist[n];

	      if(mode == 0 || mode == -1)
		{
		  if(Head[target] != Head[j])	/* only if not yet linked */
		    {

		      if(mode == 0)
			endrun(87654);

		      if(Len[Head[target]] > Len[Head[j]])	/* p group is longer */
			{
			  p = target;
			  s = j;
			}
		      else
			{
			  p = j;
			  s = target;
			}
		      Next[Tail[Head[p]]] = Head[s];

		      Tail[Head[p]] = Tail[Head[s]];

		      Len[Head[p]] += Len[Head[s]];

		      ss = Head[s];
		      do
			Head[ss] = Head[p];
		      while((ss = Next[ss]) >= 0);

		      if(MinID[Head[s]] < MinID[Head[p]])
			{
			  MinID[Head[p]] = MinID[Head[s]];
			  MinIDTask[Head[p]] = MinIDTask[Head[s]];
			}
		    }
		}
	      else		/* mode is 1 */
		{
		  if(MinID[Head[j]] > FoFDataGet[target].MinID)
		    {
		      MinID[Head[j]] = FoFDataGet[target].MinID;
		      MinIDTask[Head[j]] = FoFDataGet[target].MinIDTask;
		      links++;
		    }
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = FoFDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  return links;
}



void fof_compile_catalogue(void)
{
  int i, j, start, nimport, ngrp, recvTask;
  fof_group_list *get_FOF_GList;

  /* sort according to MinID */
  qsort(FOF_PList, NumPart, sizeof(struct fof_particle_list), fof_compare_FOF_PList_MinID);

  for(i = 0; i < NumPart; i++)
    {
      FOF_GList[i].MinID = FOF_PList[i].MinID;
      FOF_GList[i].MinIDTask = FOF_PList[i].MinIDTask;
      if(FOF_GList[i].MinIDTask == (MyIDType)ThisTask)
	{
	  FOF_GList[i].LocCount = 1;
	  FOF_GList[i].ExtCount = 0;
#ifdef FOF_DENSITY_SPLIT_TYPES
	  if(((1 << P[FOF_PList[i].Pindex].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
	    FOF_GList[i].LocDMCount = 1;
	  else
	    FOF_GList[i].LocDMCount = 0;
	  FOF_GList[i].ExtDMCount = 0;
#endif
	}
      else
	{
	  FOF_GList[i].LocCount = 0;
	  FOF_GList[i].ExtCount = 1;
#ifdef FOF_DENSITY_SPLIT_TYPES
	  FOF_GList[i].LocDMCount = 0;
	  if(((1 << P[FOF_PList[i].Pindex].Type) & (MyFOF_PRIMARY_LINK_TYPES)))
	    FOF_GList[i].ExtDMCount = 1;
	  else
	    FOF_GList[i].ExtDMCount = 0;
#endif
	}
    }

  /* eliminate duplicates in FOF_GList with respect to MinID */

  if(NumPart)
    NgroupsExt = 1;
  else
    NgroupsExt = 0;

  for(i = 1, start = 0; i < NumPart; i++)
    {
      if(FOF_GList[i].MinID == FOF_GList[start].MinID)
	{
	  FOF_GList[start].LocCount += FOF_GList[i].LocCount;
	  FOF_GList[start].ExtCount += FOF_GList[i].ExtCount;
#ifdef FOF_DENSITY_SPLIT_TYPES
	  FOF_GList[start].LocDMCount += FOF_GList[i].LocDMCount;
	  FOF_GList[start].ExtDMCount += FOF_GList[i].ExtDMCount;
#endif
	}
      else
	{
	  start = NgroupsExt;
	  FOF_GList[start] = FOF_GList[i];
	  NgroupsExt++;
	}
    }


  /* sort the remaining ones according to task */
  qsort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_MinIDTask);
  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      if(j == ThisTask)		/* we will not exchange the ones that are local */
	Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  get_FOF_GList =
    (fof_group_list *) mymalloc("get_FOF_GList", nimport * sizeof(fof_group_list));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&FOF_GList[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(fof_group_list), MPI_BYTE,
			   recvTask, TAG_FOF_C,
			   &get_FOF_GList[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(fof_group_list), MPI_BYTE,
			   recvTask, TAG_FOF_C, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = i;


  /* sort the groups according to MinID */
  qsort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_MinID);
  qsort(get_FOF_GList, nimport, sizeof(fof_group_list), fof_compare_FOF_GList_MinID);
  /* merge the imported ones with the local ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(7973);
	}

      if(get_FOF_GList[i].LocCount != 0)
	endrun(123);

      if(FOF_GList[start].MinIDTask != (MyIDType)ThisTask)
	endrun(124);

      FOF_GList[start].ExtCount += get_FOF_GList[i].ExtCount;
#ifdef FOF_DENSITY_SPLIT_TYPES
      FOF_GList[start].ExtDMCount += get_FOF_GList[i].ExtDMCount;
#endif
    }

  /* copy the size information back into the list, to inform the others */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(797831);
	}

      get_FOF_GList[i].ExtCount = FOF_GList[start].ExtCount;
      get_FOF_GList[i].LocCount = FOF_GList[start].LocCount;
#ifdef FOF_DENSITY_SPLIT_TYPES
      get_FOF_GList[i].ExtDMCount = FOF_GList[start].ExtDMCount;
      get_FOF_GList[i].LocDMCount = FOF_GList[start].LocDMCount;
#endif
    }

  /* sort the imported/exported list according to MinIDTask */
  qsort(get_FOF_GList, nimport, sizeof(fof_group_list), fof_compare_FOF_GList_MinIDTask);
  qsort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_MinIDTask);

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = ThisTask;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&get_FOF_GList[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(fof_group_list), MPI_BYTE,
			   recvTask, TAG_FOF_D,
			   &FOF_GList[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(fof_group_list), MPI_BYTE,
			   recvTask, TAG_FOF_D, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  myfree(get_FOF_GList);

  /* eliminate all groups that are too small, and count local groups */
  for(i = 0, Ngroups = 0, Nids = 0; i < NgroupsExt; i++)
    {
#ifdef FOF_DENSITY_SPLIT_TYPES
      if(FOF_GList[i].LocDMCount + FOF_GList[i].ExtDMCount < MyFOF_GROUP_MIN_SIZE)
#else
      if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount < MyFOF_GROUP_MIN_SIZE)
#endif
	{
	  FOF_GList[i] = FOF_GList[NgroupsExt - 1];
	  NgroupsExt--;
	  i--;
	}
      else
	{
	  if(FOF_GList[i].MinIDTask == (MyIDType)ThisTask)
	    {
	      Ngroups++;
	      Nids += FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
	    }
	}
    }

  /* sort the group list according to MinID */
  qsort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_MinID);
}



void fof_compute_group_properties(int gr, int start, int len)
{
  int j, k, index;
  double xyz[3];

  Group[gr].Len = 0;
  Group[gr].Mass = 0;
#ifdef GALSF
  Group[gr].Sfr = 0;
#endif
#ifdef BLACK_HOLES
  Group[gr].BH_Mass = 0;
  Group[gr].BH_Mdot = 0;
#ifdef BH_SEED_FROM_FOF
  Group[gr].MinPot = BHPOTVALUEINIT;
  Group[gr].index_maxdens = Group[gr].task_maxdens = -1;
#endif
#endif

  for(k = 0; k < 3; k++)
    {
      Group[gr].CM[k] = 0;
      Group[gr].Vel[k] = 0;
      Group[gr].FirstPos[k] = P[FOF_PList[start].Pindex].Pos[k];
    }

  for(k = 0; k < 6; k++)
    {
      Group[gr].LenType[k] = 0;
      Group[gr].MassType[k] = 0;
    }

  for(k = 0; k < len; k++)
    {
      index = FOF_PList[start + k].Pindex;

      Group[gr].Len++;
      Group[gr].Mass += P[index].Mass;
      Group[gr].LenType[P[index].Type]++;
      Group[gr].MassType[P[index].Type] += P[index].Mass;


#ifdef GALSF
      if(P[index].Type == 0)
	Group[gr].Sfr += SphP[index].Sfr;
#endif

#ifdef BLACK_HOLES
      if(P[index].Type == 5)
	{
	  Group[gr].BH_Mdot += BPP(index).BH_Mdot;
	  Group[gr].BH_Mass += BPP(index).BH_Mass;
	}

#ifdef BH_SEED_FROM_FOF
#if (BH_SEED_FROM_FOF==0)
      if(P[index].Type==0)
#elif (BH_SEED_FROM_FOF==1)
      if(P[index].Type==4)
#endif
         if(P[index].Potential < Group[gr].MinPot)
        {
          Group[gr].MinPot = P[index].Potential;
          Group[gr].index_maxdens = index;
          Group[gr].task_maxdens = ThisTask;
        }
#endif
#endif // BLACK_HOLES

        for(j = 0; j < 3; j++) {xyz[j] = P[index].Pos[j] - Group[gr].FirstPos[j];}
        NEAREST_XYZ(xyz[0],xyz[1],xyz[2],-1);
        for(j = 0; j < 3; j++)
        {
            Group[gr].CM[j] += P[index].Mass * xyz[j];
            Group[gr].Vel[j] += P[index].Mass * P[index].Vel[j];
        }
    }
}


void fof_exchange_group_data(void)
{
  group_properties *get_Group;
  int i, j, ngrp, recvTask, nimport, start;
  double xyz[3];

  /* sort the groups according to task */
  qsort(Group, NgroupsExt, sizeof(group_properties), fof_compare_Group_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      if(j == ThisTask)		/* we will not exchange the ones that are local */
	Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  get_Group = (group_properties *) mymalloc("get_Group", sizeof(group_properties) * nimport);

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group data */
	      MPI_Sendrecv(&Group[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(group_properties), MPI_BYTE,
			   recvTask, TAG_FOF_E,
			   &get_Group[Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(group_properties), MPI_BYTE,
			   recvTask, TAG_FOF_E, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  /* sort the groups again according to MinID */
  qsort(Group, NgroupsExt, sizeof(group_properties), fof_compare_Group_MinID);
  qsort(get_Group, nimport, sizeof(group_properties), fof_compare_Group_MinID);
  /* now add in the partial imported group data to the main ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(Group[start].MinID < get_Group[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(797890);
	}

      Group[start].Len += get_Group[i].Len;
      Group[start].Mass += get_Group[i].Mass;

      for(j = 0; j < 6; j++)
	{
	  Group[start].LenType[j] += get_Group[i].LenType[j];
	  Group[start].MassType[j] += get_Group[i].MassType[j];
	}

#ifdef GALSF
      Group[start].Sfr += get_Group[i].Sfr;
#endif
#ifdef BLACK_HOLES
      Group[start].BH_Mdot += get_Group[i].BH_Mdot;
      Group[start].BH_Mass += get_Group[i].BH_Mass;
#ifdef BH_SEED_FROM_FOF
      if(get_Group[i].MinPot < Group[start].MinPot)
        {
          Group[start].MinPot = get_Group[i].MinPot;
          Group[start].index_maxdens = get_Group[i].index_maxdens;     // "index" and "task" refer to MinPot here
          Group[start].task_maxdens = get_Group[i].task_maxdens;
        }
#endif
#endif

        for(j = 0; j < 3; j++) {xyz[j] = get_Group[i].CM[j] / get_Group[i].Mass + get_Group[i].FirstPos[j] - Group[start].FirstPos[j];}
        NEAREST_XYZ(xyz[0],xyz[1],xyz[2],1);
        for(j = 0; j < 3; j++)
        {
            Group[start].CM[j] += get_Group[i].Mass * xyz[j];
            Group[start].Vel[j] += get_Group[i].Vel[j];
        }
    }

  myfree(get_Group);
}

void fof_finish_group_properties(void)
{
  double cm[3];
  int i, j, ngr;

  for(i = 0; i < NgroupsExt; i++)
    {
      if(Group[i].MinIDTask == (MyIDType)ThisTask)
	{
	  for(j = 0; j < 3; j++)
	    {
	      Group[i].Vel[j] /= Group[i].Mass;

            cm[j] = Group[i].CM[j] / Group[i].Mass + Group[i].FirstPos[j];
            cm[j] = WRAP_POSITION_UNIFORM_BOX(cm[j]);
	      Group[i].CM[j] = cm[j];
	    }
	}
    }

  /* eliminate the non-local groups */
  for(i = 0, ngr = NgroupsExt; i < ngr; i++)
    {
      if(Group[i].MinIDTask != (MyIDType)ThisTask)
	{
	  Group[i] = Group[ngr - 1];
	  i--;
	  ngr--;
	}
    }

  if(ngr != Ngroups)
    endrun(876889);

  qsort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_MinID);
}



void fof_save_groups(int num)
{
  int i, j, start, lenloc, nprocgroup, primaryTask, groupTask, ngr, totlen;
  long long totNids;
  char buf[500];
  double t0, t1;

  PRINT_STATUS("start global sorting of group catalogues");
    
  t0 = my_second();

  /* assign group numbers (at this point, both Group and FOF_GList are sorted by MinID) */
  for(i = 0; i < NgroupsExt; i++)
    {
      FOF_GList[i].LocCount += FOF_GList[i].ExtCount;	/* total length */
      FOF_GList[i].ExtCount = ThisTask;	/* original task */
#ifdef FOF_DENSITY_SPLIT_TYPES
      FOF_GList[i].LocDMCount += FOF_GList[i].ExtDMCount;	/* total length */
      FOF_GList[i].ExtDMCount = ThisTask;	/* not longer needed/used (hopefully) */
#endif
    }

  parallel_sort(FOF_GList, NgroupsExt, sizeof(fof_group_list),
		fof_compare_FOF_GList_LocCountTaskDiffMinID);

  for(i = 0, ngr = 0; i < NgroupsExt; i++)
    {
      if((MyIDType)FOF_GList[i].ExtCount == FOF_GList[i].MinIDTask)
	ngr++;

      FOF_GList[i].GrNr = ngr;
    }

  MPI_Allgather(&ngr, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  for(i = 0; i < NgroupsExt; i++)
    FOF_GList[i].GrNr += Send_offset[ThisTask];


  MPI_Allreduce(&ngr, &i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(i != TotNgroups)
    {
      printf("i=%d\n", i);
      endrun(123123);
    }

  /* bring the group list back into the original order */
  parallel_sort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_ExtCountMinID);

  /* Assign the group numbers to the group properties array */
  for(i = 0, start = 0; i < Ngroups; i++)
    {
      while(FOF_GList[start].MinID < Group[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    endrun(7297890);
	}
      Group[i].GrNr = FOF_GList[start].GrNr;
    }

  /* sort the groups according to group-number */
  parallel_sort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);

  /* fill in the offset-values */
  for(i = 0, totlen = 0; i < Ngroups; i++)
    {
      if(i > 0)
	Group[i].Offset = Group[i - 1].Offset + Group[i - 1].Len;
      else
	Group[i].Offset = 0;
      totlen += Group[i].Len;
    }

  MPI_Allgather(&totlen, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  unsigned int *uoffset = (unsigned int *)mymalloc("uoffset", NTask * sizeof(unsigned int));

  for(j = 1, uoffset[0] = 0; j < NTask; j++)
    uoffset[j] = uoffset[j - 1] + Send_count[j - 1];

  for(i = 0; i < Ngroups; i++)
    Group[i].Offset += uoffset[ThisTask];

  myfree(uoffset);

  /* prepare list of ids with assigned group numbers */

  ID_list = (fof_id_list *)mymalloc("ID_list", sizeof(fof_id_list) * NumPart);

#ifdef SUBFIND
  for(i = 0; i < NumPart; i++)
    P[i].GrNr = TotNgroups + 1;	/* will mark particles that are not in any group */
#endif

  for(i = 0, start = 0, Nids = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
	{
	  start++;
	  if(start > NumPart)
	    endrun(78);
	}

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
	endrun(1313);

      for(lenloc = 0; start + lenloc < NumPart;)
	if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
	  {
	    ID_list[Nids].GrNr = FOF_GList[i].GrNr;
	    ID_list[Nids].ID = P[FOF_PList[start + lenloc].Pindex].ID;
#ifdef SUBFIND
	    P[FOF_PList[start + lenloc].Pindex].GrNr = FOF_GList[i].GrNr;
#endif
	    Nids++;
	    lenloc++;
	  }
	else
	  break;

      start += lenloc;
    }

  sumup_large_ints(1, &Nids, &totNids);

  MPI_Allgather(&Nids, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];


  if(totNids != TotNids)
    {
      printf("Task=%d Nids=%d totNids=%d TotNids=%d\n", ThisTask, Nids, (int) totNids, (int) TotNids);
      endrun(12);
    }

  /* sort the particle IDs according to group-number */
  parallel_sort(ID_list, Nids, sizeof(fof_id_list), fof_compare_ID_list_GrNrID);

  t1 = my_second();
  PRINT_STATUS("Group catalogues globally sorted. took = %g sec. Started saving of group catalogue", timediff(t0, t1));
  t0 = my_second();

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);


  if(NTask < All.NumFilesWrittenInParallel) {printf("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n"); endrun(241931);}

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  primaryTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (primaryTask + groupTask))	/* ok, it's this processor's turn */
	fof_save_local_catalogue(num);
      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }

  myfree(ID_list);

  t1 = my_second();
  PRINT_STATUS("Group catalogues saved. took = %g sec", timediff(t0, t1));
}



void fof_save_local_catalogue(int num)
{
  FILE *fd;
  float *mass, *cm, *vel;
  char fname[500];
  int i, j, *len;
  MyIDType *ids;

  sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);
  if(!(fd = fopen(fname, "w")))
    {
      printf("can't open file `%s`\n", fname);
      endrun(1183);
    }

  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);

  /* group len */
  len = (int *)mymalloc("len", Ngroups * sizeof(int));
  for(i = 0; i < Ngroups; i++)
    len[i] = Group[i].Len;
  my_fwrite(len, Ngroups, sizeof(int), fd);
  myfree(len);

  /* offset into id-list */
  len = (int *)mymalloc("len", Ngroups * sizeof(int));
  for(i = 0; i < Ngroups; i++)
    len[i] = Group[i].Offset;
  my_fwrite(len, Ngroups, sizeof(int), fd);
  myfree(len);

  /* mass */
  mass = (float *)mymalloc("mass", Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].Mass;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* CM */
  cm = (float *)mymalloc("cm", Ngroups * 3 * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    for(j = 0; j < 3; j++)
      cm[i * 3 + j] = Group[i].CM[j];
  my_fwrite(cm, Ngroups, 3 * sizeof(float), fd);
  myfree(cm);

  /* vel */
  vel = (float *)mymalloc("vel", Ngroups * 3 * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    for(j = 0; j < 3; j++)
      vel[i * 3 + j] = Group[i].Vel[j];
  my_fwrite(vel, Ngroups, 3 * sizeof(float), fd);
  myfree(vel);

  /* group len for each type */
  len = (int *)mymalloc("len", Ngroups * 6 * sizeof(int));
  for(i = 0; i < Ngroups; i++)
    for(j = 0; j < 6; j++)
      len[i * 6 + j] = Group[i].LenType[j];
  my_fwrite(len, Ngroups, 6 * sizeof(int), fd);
  myfree(len);

  /* group mass for each type */
  mass = (float *)mymalloc("mass", Ngroups * 6 * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    for(j = 0; j < 6; j++)
      mass[i * 6 + j] = Group[i].MassType[j];
  my_fwrite(mass, Ngroups, 6 * sizeof(float), fd);
  myfree(mass);

#ifdef GALSF
  /* sfr */
  mass = (float *)mymalloc("mass", Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].Sfr;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);
#endif

#ifdef BLACK_HOLES
  /* BH_Mass */
  mass = (float *)mymalloc("mass", Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].BH_Mass;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* BH_Mdot */
  mass = (float *)mymalloc("mass", Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].BH_Mdot;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);
#endif

  fclose(fd);


  ids = (MyIDType *) ID_list;
  for(i = 0; i < Nids; i++)
    ids[i] = ID_list[i].ID;

  sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, ThisTask);
  if(!(fd = fopen(fname, "w")))
    {
      printf("can't open file `%s`\n", fname);
      endrun(1184);
    }

  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);	/* this is the number of IDs in previous files */
  my_fwrite(ids, sizeof(MyIDType), Nids, fd);
  fclose(fd);
}


void fof_find_nearest_dmparticle(void)
{
  int i, j, n, ntot, dummy;
  int ndone, ndone_flag, ngrp, recvTask, place, nexport, nimport, npleft, iter;

  PRINT_STATUS("Start finding nearest dm-particle (presently allocated=%g MB)",AllocatedBytes / (1024.0 * 1024.0));
  fof_nearest_distance = (float *) mymalloc("fof_nearest_distance", sizeof(float) * NumPart);
  fof_nearest_hsml = (float *) mymalloc("fof_nearest_hsml", sizeof(float) * NumPart);

  for(n = 0; n < NumPart; n++)
    {
      if(((1 << P[n].Type) & (MyFOF_SECONDARY_LINK_TYPES)))
	{
	  fof_nearest_distance[n] = 1.0e30;
	  fof_nearest_hsml[n] = 0.1 * LinkL;
	}
    }

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct fofdata_in) + sizeof(struct fofdata_out) +
					     sizemax(sizeof(struct fofdata_in), sizeof(struct fofdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  iter = 0;
  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
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
	    if(((1 << P[i].Type) & (MyFOF_SECONDARY_LINK_TYPES)))
	      {
		if(fof_nearest_distance[i] > 1.0e29)
		  {
		    if(fof_find_nearest_dmparticle_evaluate(i, 0, &nexport, Send_count) < 0)
		      break;
		  }
	      }

	  MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

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

	  FoFDataGet = (struct fofdata_in *) mymalloc("FoFDataGet", nimport * sizeof(struct fofdata_in));
	  FoFDataIn = (struct fofdata_in *) mymalloc("FoFDataIn", nexport * sizeof(struct fofdata_in));

      PRINT_STATUS("still finding nearest... (presently allocated=%g MB)",AllocatedBytes / (1024.0 * 1024.0));
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      FoFDataIn[j].Pos[0] = P[place].Pos[0];
	      FoFDataIn[j].Pos[1] = P[place].Pos[1];
	      FoFDataIn[j].Pos[2] = P[place].Pos[2];
	      FoFDataIn[j].Hsml = fof_nearest_hsml[place];

	      memcpy(FoFDataIn[j].NodeList,
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
		      MPI_Sendrecv(&FoFDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
				   recvTask, TAG_FOF_F,
				   &FoFDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
				   recvTask, TAG_FOF_F, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  myfree(FoFDataIn);
	  FoFDataResult =
	    (struct fofdata_out *) mymalloc("FoFDataResult", nimport * sizeof(struct fofdata_out));
	  FoFDataOut = (struct fofdata_out *) mymalloc("FoFDataOut", nexport * sizeof(struct fofdata_out));

	  for(j = 0; j < nimport; j++)
	    {
	      fof_find_nearest_dmparticle_evaluate(j, 1, &dummy, &dummy);
	    }

	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&FoFDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct fofdata_out),
				   MPI_BYTE, recvTask, TAG_FOF_G,
				   &FoFDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct fofdata_out),
				   MPI_BYTE, recvTask, TAG_FOF_G, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }

	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      if(FoFDataOut[j].Distance < fof_nearest_distance[place])
		{
		  fof_nearest_distance[place] = FoFDataOut[j].Distance;
		  MinID[place] = FoFDataOut[j].MinID;
		  MinIDTask[place] = FoFDataOut[j].MinIDTask;
		}
	    }

	  if(i >= NumPart)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	  myfree(FoFDataOut);
	  myfree(FoFDataResult);
	  myfree(FoFDataGet);
	}
      while(ndone < NTask);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
	{
	  if(((1 << P[i].Type) & (MyFOF_SECONDARY_LINK_TYPES)))
	    {
	      if(fof_nearest_distance[i] > 1.0e29)
		{
		  if(fof_nearest_hsml[i] < 4 * LinkL)	/* we only search out to a maximum distance */
		    {
		      /* need to redo this particle */
		      npleft++;
		      fof_nearest_hsml[i] *= 2.0;
		      if(iter >= MAXITER - 10)
			{
			  printf("i=%d task=%d ID=%llu Hsml=%g  pos=(%g|%g|%g)\n",
				 i, ThisTask, (unsigned long long) P[i].ID, fof_nearest_hsml[i],
				 P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
			  fflush(stdout);
			}
		    }
		  else
		    {
		      fof_nearest_distance[i] = 0;	/* we not continue to search for this particle */
		    }
		}
	    }
	}

      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(ntot > 0)
	{
	  iter++;
	  if(iter > 10 && ThisTask == 0) {printf("fof-nearest iteration %d: need to repeat for %d particles.\n", iter, ntot); fflush(stdout);}
	  if(iter > MAXITER) {printf("failed to converge in fof-nearest\n"); fflush(stdout); endrun(1159);}
	}
    }
  while(ntot > 0);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  myfree(fof_nearest_hsml);
  myfree(fof_nearest_distance);

  PRINT_STATUS("done finding nearest dm-particle");
}


/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here!!! -- */
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n, index, listindex = 0;
  int startnode, numngb_inbox;
  double h, r2max;
  double dx, dy, dz, r2;
  MyDouble *pos;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = fof_nearest_hsml[target];
    }
  else
    {
      pos = FoFDataGet[target].Pos;
      h = FoFDataGet[target].Hsml;
    }

  index = -1;
  r2max = 1.0e30;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = FoFDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
        numngb_inbox = ngb_treefind_variable_targeted(pos, h, target, &startnode, mode, nexport, nsend_local, MyFOF_PRIMARY_LINK_TYPES); // MyFOF_PRIMARY_LINK_TYPES defines which types of particles we search for
        if(numngb_inbox < 0) {return -2;}

	  for(n = 0; n < numngb_inbox; n++)
	    {
            j = Ngblist[n];
            dx = pos[0] - P[j].Pos[0];
            dy = pos[1] - P[j].Pos[1];
            dz = pos[2] - P[j].Pos[2];
            NEAREST_XYZ(dx,dy,dz,1);
            r2 = dx * dx + dy * dy + dz * dz;
            if(r2 < r2max && r2 < h * h)
		{
		  index = j;
		  r2max = r2;
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = FoFDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  if(mode == 0)
    {
      if(index >= 0)
	{
	  fof_nearest_distance[target] = sqrt(r2max);
	  MinID[target] = MinID[Head[index]];
	  MinIDTask[target] = MinIDTask[Head[index]];
	}
    }
  else
    {
      if(index >= 0)
	{
	  FoFDataResult[target].Distance = sqrt(r2max);
	  FoFDataResult[target].MinID = MinID[Head[index]];
	  FoFDataResult[target].MinIDTask = MinIDTask[Head[index]];
	}
      else
	FoFDataResult[target].Distance = 2.0e30;
    }
  return 0;
}




#ifdef BH_SEED_FROM_FOF

void fof_make_black_holes(void)
{
  int i, j, n, ntot;
  long nexport, nimport;
  int recvTask, level;
  int *import_indices, *export_indices;
  gsl_rng *random_generator_forbh;
  double random_number_forbh=0, unitmass_in_msun;

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
#if (BH_SEED_FROM_FOF==0)
    if(Group[i].MassType[1] >= (All.OmegaMatter - All.OmegaBaryon) / All.OmegaMatter * All.MinFoFMassForNewSeed)
#elif (BH_SEED_FROM_FOF==1)
    if(Group[i].MassType[4] > All.MinFoFMassForNewSeed)
#endif
	if(Group[i].LenType[5] == 0)
	  {
	    if(Group[i].index_maxdens >= 0)
	      Send_count[Group[i].task_maxdens]++;
	  }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
	{
	  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	}
    }

  import_indices = mymalloc("import_indices", nimport * sizeof(int));
  export_indices = mymalloc("export_indices", nexport * sizeof(int));

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
#if (BH_SEED_FROM_FOF==0)
        if(Group[i].MassType[1] >= (All.OmegaMatter - All.OmegaBaryon) / All.OmegaMatter * All.MinFoFMassForNewSeed)
#elif (BH_SEED_FROM_FOF==1)
        if(Group[i].MassType[4] > All.MinFoFMassForNewSeed)
#endif
	if(Group[i].LenType[5] == 0)
	  {
	    if(Group[i].index_maxdens >= 0)
	      export_indices[Send_offset[Group[i].task_maxdens] +
			     Send_count[Group[i].task_maxdens]++] = Group[i].index_maxdens;
	  }
    }

  memcpy(&import_indices[Recv_offset[ThisTask]], &export_indices[Send_offset[ThisTask]],
	 Send_count[ThisTask] * sizeof(int));

  for(level = 1; level < (1 << PTask); level++)
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&export_indices[Send_offset[recvTask]],
		     Send_count[recvTask] * sizeof(int),
		     MPI_BYTE, recvTask, TAG_FOF_I,
		     &import_indices[Recv_offset[recvTask]],
		     Recv_count[recvTask] * sizeof(int),
		     MPI_BYTE, recvTask, TAG_FOF_I, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  MPI_Allreduce(&nimport, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  PRINT_STATUS("Making %d new black hole particles", ntot);
  All.TotBHs += ntot;

  for(n = 0; n < nimport; n++)
    {
#if (BH_SEED_FROM_FOF==0)
        if(P[import_indices[n]].Type != 0)
#elif (BH_SEED_FROM_FOF==1)
            if(P[import_indices[n]].Type != 4)
#endif
                endrun(7772);
        
        P[import_indices[n]].Type = 5;    /* make it a black hole particle */
#ifdef GALSF
        P[import_indices[n]].StellarAge = All.Time; /* reset formation time to match BH formation */
#endif
        /* generate BH mass */
        if(All.SeedBlackHoleMassSigma > 0)
        {
            /* compute gaussian random number: mean=0, sigma=All.SeedBlackHoleMassSigma */
            random_generator_forbh = gsl_rng_alloc(gsl_rng_ranlxd1);
            gsl_rng_set(random_generator_forbh,P[import_indices[n]].ID+17);
            random_number_forbh = gsl_ran_gaussian(random_generator_forbh, All.SeedBlackHoleMassSigma);
            BPP(import_indices[n]).BH_Mass = pow( 10., log10(All.SeedBlackHoleMass) + random_number_forbh );
            unitmass_in_msun = UNIT_MASS_IN_SOLAR;
            if( BPP(import_indices[n]).BH_Mass < 100./unitmass_in_msun )
                BPP(import_indices[n]).BH_Mass = 100./unitmass_in_msun;      // enforce lower limit of Mseed = 100 x Msun
        } else {
            BPP(import_indices[n]).BH_Mass = All.SeedBlackHoleMass;
        }
        BPP(import_indices[n]).BH_Mdot = 0;
        /* set hydro-ish variables */
        if(BPP(import_indices[n]).Type == 0){
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            P[import_indices[n]].Mass = SphP[import_indices[n]].MassTrue + SphP[import_indices[n]].dMass;
#endif
            P[import_indices[n]].DensAroundStar = SphP[import_indices[n]].Density;
        }
        /* set some specific BH variables that are needed below */
#ifdef BH_INCREASE_DYNAMIC_MASS
        P[import_indices[n]].Mass *= BH_INCREASE_DYNAMIC_MASS;
#endif
#ifdef BH_ALPHADISK_ACCRETION
        BPP(import_indices[n]).BH_Mass_AlphaDisk = All.SeedAlphaDiskMass;
#endif
#ifdef BH_WIND_SPAWN
        BPP(import_indices[n]).unspawned_wind_mass = 0;
#endif
#ifdef BH_COUNTPROGS
        BPP(import_indices[n]).BH_CountProgs = 1;
#endif
        /* record that we actually made a BH, count numbers for book-keeping in domains */
#if (BH_SEED_FROM_FOF != 1)
        Stars_converted++;
        TimeBinCountSph[P[import_indices[n]].TimeBin]--;
#endif

    }

  All.TotN_gas -= ntot;

  myfree(export_indices);
  myfree(import_indices);
}

#endif // BH_SEED_FROM_FOF



#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==1)

struct group_mass_MinID
{
  unsigned long long MinID;
  double mass;
};

int compare_group_mass_ID(const void *a, const void *b)
{
  if(((struct group_mass_MinID *) a)->MinID < (((struct group_mass_MinID *) b)->MinID)) {return -1;}
  if(((struct group_mass_MinID *) a)->MinID > (((struct group_mass_MinID *) b)->MinID)) {return +1;}
  return 0;
}

void fof_assign_HostHaloMass(void)	/* assigns mass of host FoF group to SphP[].HostHaloMass for SPH particles */
{
  int i, j, k, start, lenloc, nimport;
  struct group_mass_MinID *required_groups, *groups_to_export;

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)	/* loop over all groups for which at least one particle is on this task */
    Send_count[FOF_GList[i].MinIDTask]++;	/* its FoF group properties are stored on Task = MinIDTask */

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

  qsort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_MinIDTask_MinID);

  required_groups =
    (struct group_mass_MinID *) mymalloc("required_groups", NgroupsExt * sizeof(struct group_mass_MinID));

  MPI_Datatype mpi_groups_mass_MinID;
  MPI_Datatype used_types[2] = { MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE };
  int used_blocklen[2] = { 1, 1 };
  MPI_Aint disp[2], cur_addr, start_addr;
  MPI_Get_address(&required_groups[0], &start_addr);
  MPI_Get_address(&required_groups[0].MinID, &cur_addr);
  disp[0] = cur_addr - start_addr;
  MPI_Get_address(&required_groups[0].mass, &cur_addr);
  disp[1] = cur_addr - start_addr;
  MPI_Type_create_struct(2, used_blocklen, disp, used_types, &mpi_groups_mass_MinID);
  MPI_Type_commit(&mpi_groups_mass_MinID);	/* defines an MPI datatpye containing group mass and group MinID */

  for(i = 0; i < NgroupsExt; i++)
    required_groups[i].MinID = FOF_GList[i].MinID;

  groups_to_export =
    (struct group_mass_MinID *) mymalloc("groups_to_export", nimport * sizeof(struct group_mass_MinID));

  /* send list of groups for which we need the masses */
  MPI_Alltoallv(required_groups, Send_count, Send_offset, mpi_groups_mass_MinID,
		groups_to_export, Recv_count, Recv_offset, mpi_groups_mass_MinID, MPI_COMM_WORLD);

  for(j = 0, start = 0; j < NTask; j++)
    {
      i = 0;
      k = 0;

      while(i < Recv_count[j] && k < Ngroups)
	{
	  if(groups_to_export[start].MinID == Group[k].MinID)
	    {
	      groups_to_export[start].mass = Group[k].Mass;
	      i++;
	      k++;
	      start++;
	    }
	  else
	    k++;
	}
    }
  if(start != nimport)
    terminate("start != nimport");

  /* send group masses to requesting tasks */
  MPI_Alltoallv(groups_to_export, Recv_count, Recv_offset, mpi_groups_mass_MinID,
		required_groups, Send_count, Send_offset, mpi_groups_mass_MinID, MPI_COMM_WORLD);

  myfree(groups_to_export);
  MPI_Type_free(&mpi_groups_mass_MinID);

  qsort(required_groups, NgroupsExt, sizeof(struct group_mass_MinID), compare_group_mass_ID);

    for(i = 0; i < N_gas; i++) {if(P[i].Type==0) {SphP[i].HostHaloMass = 0;}}

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < required_groups[i].MinID)
	{
	  start++;
	  if(start > NumPart)
	    terminate("start > NumPart");
	}

      if(FOF_PList[start].MinID != required_groups[i].MinID)
	terminate("FOF_PList[start].MinID != required_groups[i].MinID");

      for(lenloc = 0; start + lenloc < NumPart;)
	if(FOF_PList[start + lenloc].MinID == required_groups[i].MinID)
	  {
	    if(P[FOF_PList[start + lenloc].Pindex].Type == 0)
	      SphP[FOF_PList[start + lenloc].Pindex].HostHaloMass = required_groups[i].mass;

	    lenloc++;
	  }
	else
	  break;

      start += lenloc;
    }

  myfree(required_groups);

  qsort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_MinID);	/* restore original order */
}

#endif
#endif // defined(GALSF_SUBGRID_WINDS) && defined(GALSF_SUBGRID_WIND_SCALING==1)




int fof_compare_FOF_PList_MinID(const void *a, const void *b)
{
  if(((struct fof_particle_list *) a)->MinID < ((struct fof_particle_list *) b)->MinID) {return -1;}
  if(((struct fof_particle_list *) a)->MinID > ((struct fof_particle_list *) b)->MinID) {return +1;}
  return 0;
}

int fof_compare_FOF_GList_MinID(const void *a, const void *b)
{
  if(((fof_group_list *) a)->MinID < ((fof_group_list *) b)->MinID) {return -1;}
  if(((fof_group_list *) a)->MinID > ((fof_group_list *) b)->MinID) {return +1;}
  return 0;
}

int fof_compare_FOF_GList_MinIDTask(const void *a, const void *b)
{
  if(((fof_group_list *) a)->MinIDTask < ((fof_group_list *) b)->MinIDTask) {return -1;}
  if(((fof_group_list *) a)->MinIDTask > ((fof_group_list *) b)->MinIDTask) {return +1;}
  return 0;
}

int fof_compare_FOF_GList_MinIDTask_MinID(const void *a, const void *b)
{
  if(((fof_group_list *) a)->MinIDTask < ((fof_group_list *) b)->MinIDTask) {return -1;}
  if(((fof_group_list *) a)->MinIDTask > ((fof_group_list *) b)->MinIDTask) {return +1;}
  if(((fof_group_list *) a)->MinID < ((fof_group_list *) b)->MinID) {return -1;}
  if(((fof_group_list *) a)->MinID > ((fof_group_list *) b)->MinID) {return +1;}
  return 0;
}

int fof_compare_FOF_GList_LocCountTaskDiffMinID(const void *a, const void *b)
{
  if(((fof_group_list *) a)->LocCount > ((fof_group_list *) b)->LocCount) {return -1;}
  if(((fof_group_list *) a)->LocCount < ((fof_group_list *) b)->LocCount) {return +1;}
  if(((fof_group_list *) a)->MinID < ((fof_group_list *) b)->MinID) {return -1;}
  if(((fof_group_list *) a)->MinID > ((fof_group_list *) b)->MinID) {return +1;}
  if((((fof_group_list *) a)->ExtCount - ((fof_group_list *) a)->MinIDTask) <
     (((fof_group_list *) b)->ExtCount - ((fof_group_list *) b)->MinIDTask)) {return -1;}
  if((((fof_group_list *) a)->ExtCount - ((fof_group_list *) a)->MinIDTask) >
     (((fof_group_list *) b)->ExtCount - ((fof_group_list *) b)->MinIDTask)) {return +1;}
  return 0;
}

int fof_compare_FOF_GList_ExtCountMinID(const void *a, const void *b)
{
  if(((fof_group_list *) a)->ExtCount < ((fof_group_list *) b)->ExtCount) {return -1;}
  if(((fof_group_list *) a)->ExtCount > ((fof_group_list *) b)->ExtCount) {return +1;}
  if(((fof_group_list *) a)->MinID < ((fof_group_list *) b)->MinID) {return -1;}
  if(((fof_group_list *) a)->MinID > ((fof_group_list *) b)->MinID) {return +1;}
  return 0;
}

int fof_compare_Group_MinID(const void *a, const void *b)
{
  if(((group_properties *) a)->MinID < ((group_properties *) b)->MinID) {return -1;}
  if(((group_properties *) a)->MinID > ((group_properties *) b)->MinID) {return +1;}
  return 0;
}

int fof_compare_Group_GrNr(const void *a, const void *b)
{
  if(((group_properties *) a)->GrNr < ((group_properties *) b)->GrNr) {return -1;}
  if(((group_properties *) a)->GrNr > ((group_properties *) b)->GrNr) {return +1;}
  return 0;
}

int fof_compare_Group_MinIDTask(const void *a, const void *b)
{
  if(((group_properties *) a)->MinIDTask < ((group_properties *) b)->MinIDTask) {return -1;}
  if(((group_properties *) a)->MinIDTask > ((group_properties *) b)->MinIDTask) {return +1;}
  return 0;
}

int fof_compare_Group_MinIDTask_MinID(const void *a, const void *b)
{
  if(((group_properties *) a)->MinIDTask < ((group_properties *) b)->MinIDTask) {return -1;}
  if(((group_properties *) a)->MinIDTask > ((group_properties *) b)->MinIDTask) {return +1;}
  if(((group_properties *) a)->MinID < ((group_properties *) b)->MinID) {return -1;}
  if(((group_properties *) a)->MinID > ((group_properties *) b)->MinID) {return +1;}

  return 0;
}


int fof_compare_Group_Len(const void *a, const void *b)
{
  if(((group_properties *) a)->Len > ((group_properties *) b)->Len) {return -1;}
  if(((group_properties *) a)->Len < ((group_properties *) b)->Len) {return +1;}
  return 0;
}



int fof_compare_ID_list_GrNrID(const void *a, const void *b)
{
  if(((fof_id_list *) a)->GrNr < ((fof_id_list *) b)->GrNr) {return -1;}
  if(((fof_id_list *) a)->GrNr > ((fof_id_list *) b)->GrNr) {return +1;}
  if(((fof_id_list *) a)->ID < ((fof_id_list *) b)->ID) {return -1;}
  if(((fof_id_list *) a)->ID > ((fof_id_list *) b)->ID) {return +1;}

  return 0;
}








#ifdef IO_SUBFIND_READFOF_FROMIC		/* read already existing FOF instead of recomputing it */

void read_fof(int num)
{
  FILE *fd;
  double t0, t1;
  char fname[500];
  float *mass, *cm;
  int i, j, ntask, *len, count;
  MyIDType *ids;
  int *list_of_ngroups, *list_of_nids, *list_of_allgrouplen;
  int *recvoffset;
  int grnr, ngrp, sendTask, recvTask;
  int nprocgroup, primaryTask, groupTask, nid_previous;
  int fof_compare_P_SubNr(const void *a, const void *b);
    PRINT_STATUS("Trying to read preexisting FoF group catalogues...  (presently allocated=%g MB)",AllocatedBytes / (1024.0 * 1024.0));
  domain_Decomposition(1, 0, 0);

  force_treefree();


  /* start reading of group catalogue */

  if(ThisTask == 0)
    {
      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, 0);
      if(!(fd = fopen(fname, "r")))
	{
	  printf("can't read file `%s`\n", fname);
	  endrun(11831);
	}

      my_fread(&Ngroups, sizeof(int), 1, fd);
      my_fread(&TotNgroups, sizeof(int), 1, fd);
      my_fread(&Nids, sizeof(int), 1, fd);
      my_fread(&TotNids, sizeof(long long), 1, fd);
      my_fread(&ntask, sizeof(int), 1, fd);
      fclose(fd);
    }

  MPI_Bcast(&ntask, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotNgroups, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotNids, sizeof(long long), MPI_BYTE, 0, MPI_COMM_WORLD);

  t0 = my_second();

  if(NTask != ntask)
    {
      if(ThisTask == 0)
	printf
	  ("number of files (%d) in group catalogues does not match MPI-Tasks, I'm working around this.\n",
	   ntask);

      Group =
	(group_properties *) mymalloc("Group",
					     sizeof(group_properties) * (TotNgroups / NTask + NTask));

      ID_list = mymalloc("ID_list", (TotNids / NTask + NTask) * sizeof(fof_id_list));


      int filenr, target, ngroups, nids, nsend, stored;

      int *ngroup_to_get = mymalloc("ngroup_to_get", NTask * sizeof(NTask));
      int *nids_to_get = mymalloc("nids_to_get", NTask * sizeof(NTask));
      int *ngroup_obtained = mymalloc("ngroup_obtained", NTask * sizeof(NTask));
      int *nids_obtained = mymalloc("nids_obtained", NTask * sizeof(NTask));

      for(i = 0; i < NTask; i++)
	ngroup_obtained[i] = nids_obtained[i] = 0;

      for(i = 0; i < NTask - 1; i++)
	ngroup_to_get[i] = TotNgroups / NTask;
      ngroup_to_get[NTask - 1] = TotNgroups - (NTask - 1) * (TotNgroups / NTask);

      for(i = 0; i < NTask - 1; i++)
	nids_to_get[i] = (int) (TotNids / NTask);
      nids_to_get[NTask - 1] = (int) (TotNids - (NTask - 1) * (TotNids / NTask));

      Ngroups = ngroup_to_get[ThisTask];
      Nids = nids_to_get[ThisTask];



      if(ThisTask == 0)
	{
	  for(filenr = 0; filenr < ntask; filenr++)
	    {
	      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, filenr);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(11831);
		}

	      printf("reading '%s'\n", fname);
	      my_fread(&ngroups, sizeof(int), 1, fd);
	      my_fread(&TotNgroups, sizeof(int), 1, fd);
	      my_fread(&nids, sizeof(int), 1, fd);
	      my_fread(&TotNids, sizeof(long long), 1, fd);
	      my_fread(&ntask, sizeof(int), 1, fd);

	      group_properties *tmpGroup =
		(group_properties *) mymalloc("tmpGroup", sizeof(group_properties) * ngroups);

	      /* group len */
	      len = mymalloc("len", ngroups * sizeof(int));
	      my_fread(len, ngroups, sizeof(int), fd);
	      for(i = 0; i < ngroups; i++)
		tmpGroup[i].Len = len[i];
	      myfree(len);

	      /* offset into id-list */
	      len = mymalloc("len", ngroups * sizeof(int));
	      my_fread(len, ngroups, sizeof(int), fd);
	      for(i = 0; i < ngroups; i++)
		tmpGroup[i].Offset = len[i];
	      myfree(len);

	      /* mass */
	      mass = mymalloc("mass", ngroups * sizeof(float));
	      my_fread(mass, ngroups, sizeof(float), fd);
	      for(i = 0; i < ngroups; i++)
		tmpGroup[i].Mass = mass[i];
	      myfree(mass);

	      /* CM */
	      cm = mymalloc("cm", ngroups * 3 * sizeof(float));
	      my_fread(cm, ngroups, 3 * sizeof(float), fd);
	      for(i = 0; i < ngroups; i++)
		for(j = 0; j < 3; j++)
		  tmpGroup[i].CM[j] = cm[i * 3 + j];
	      myfree(cm);

	      fclose(fd);

	      target = 0;
	      stored = 0;
	      while(ngroups > 0)
		{
		  while(ngroup_to_get[target] == 0)
		    target++;

		  if(ngroups > ngroup_to_get[target])
		    nsend = ngroup_to_get[target];
		  else
		    nsend = ngroups;

		  if(target == 0)
		    memcpy(&Group[ngroup_obtained[target]], &tmpGroup[stored],
			   nsend * sizeof(group_properties));
		  else
		    {
		      MPI_Send(&nsend, 1, MPI_INT, target, TAG_N, MPI_COMM_WORLD);
		      MPI_Send(&tmpGroup[stored], nsend * sizeof(group_properties), MPI_BYTE,
			       target, TAG_PDATA, MPI_COMM_WORLD);
		    }

		  ngroup_to_get[target] -= nsend;
		  ngroup_obtained[target] += nsend;
		  ngroups -= nsend;
		  stored += nsend;
		}

	      myfree(tmpGroup);
	    }



	      /**** now ids ****/
	  for(filenr = 0; filenr < ntask; filenr++)
	    {
	      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, filenr);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(1184132);
		}
	      my_fread(&ngroups, sizeof(int), 1, fd);
	      my_fread(&TotNgroups, sizeof(int), 1, fd);
	      my_fread(&nids, sizeof(int), 1, fd);
	      my_fread(&TotNids, sizeof(long long), 1, fd);
	      my_fread(&ntask, sizeof(int), 1, fd);
	      my_fread(&nid_previous, sizeof(int), 1, fd);	/* this is the number of IDs in previous files */


	      fof_id_list *tmpID_list = mymalloc("tmpID_list", nids * sizeof(fof_id_list));

	      ids = mymalloc("ids", nids * sizeof(MyIDType));

	      my_fread(ids, sizeof(MyIDType), nids, fd);

	      for(i = 0; i < nids; i++)
		tmpID_list[i].ID = ids[i];

	      myfree(ids);

	      fclose(fd);

	      target = 0;
	      stored = 0;
	      while(nids > 0)
		{
		  while(nids_to_get[target] == 0)
		    target++;

		  if(nids > nids_to_get[target])
		    nsend = nids_to_get[target];
		  else
		    nsend = nids;

		  if(target == 0)
		    memcpy(&ID_list[nids_obtained[target]], &tmpID_list[stored],
			   nsend * sizeof(fof_id_list));
		  else
		    {
		      MPI_Send(&nsend, 1, MPI_INT, target, TAG_HEADER, MPI_COMM_WORLD);
		      MPI_Send(&tmpID_list[stored], nsend * sizeof(fof_id_list), MPI_BYTE,
			       target, TAG_SPHDATA, MPI_COMM_WORLD);
		    }

		  nids_to_get[target] -= nsend;
		  nids_obtained[target] += nsend;
		  nids -= nsend;
		  stored += nsend;
		}

	      myfree(tmpID_list);
	    }
	}
      else
	{
	  while(ngroup_to_get[ThisTask])
	    {
	      MPI_Recv(&nsend, 1, MPI_INT, 0, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      MPI_Recv(&Group[ngroup_obtained[ThisTask]], nsend * sizeof(group_properties), MPI_BYTE,
		       0, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      ngroup_to_get[ThisTask] -= nsend;
	      ngroup_obtained[ThisTask] += nsend;
	    }

	  while(nids_to_get[ThisTask])
	    {
	      MPI_Recv(&nsend, 1, MPI_INT, 0, TAG_HEADER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      MPI_Recv(&ID_list[nids_obtained[ThisTask]], nsend * sizeof(fof_id_list), MPI_BYTE,
		       0, TAG_SPHDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      nids_to_get[ThisTask] -= nsend;
	      nids_obtained[ThisTask] += nsend;
	    }
	}

      myfree(nids_obtained);
      myfree(ngroup_obtained);
      myfree(nids_to_get);
      myfree(ngroup_to_get);
    }

  else
    {
      /* read routine can continue in parallel */

      nprocgroup = NTask / All.NumFilesWrittenInParallel;
      if((NTask % All.NumFilesWrittenInParallel))
	nprocgroup++;
      primaryTask = (ThisTask / nprocgroup) * nprocgroup;
      for(groupTask = 0; groupTask < nprocgroup; groupTask++)
	{
	  if(ThisTask == (primaryTask + groupTask))	/* ok, it's this processor's turn */
	    {

	      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(11831);
		}
	      my_fread(&Ngroups, sizeof(int), 1, fd);
	      my_fread(&TotNgroups, sizeof(int), 1, fd);
	      my_fread(&Nids, sizeof(int), 1, fd);
	      my_fread(&TotNids, sizeof(long long), 1, fd);
	      my_fread(&ntask, sizeof(int), 1, fd);
	      if(NTask != ntask)
		{
		  if(ThisTask == 0)
		    printf("number of files in group catalogues needs to match MPI-Tasks\n");
		  endrun(0);
		}

	      if(ThisTask == 0)
		printf("TotNgroups=%d\n", TotNgroups);

	      Group = (group_properties *) mymalloc("Group", sizeof(group_properties) *
							   IMAX(Ngroups, TotNgroups / NTask + 1));

	      /* group len */
	      len = mymalloc("len", Ngroups * sizeof(int));
	      my_fread(len, Ngroups, sizeof(int), fd);
	      for(i = 0; i < Ngroups; i++)
		Group[i].Len = len[i];
	      myfree(len);

	      /* offset into id-list */
	      len = mymalloc("len", Ngroups * sizeof(int));
	      my_fread(len, Ngroups, sizeof(int), fd);
	      for(i = 0; i < Ngroups; i++)
		Group[i].Offset = len[i];
	      myfree(len);

	      /* mass */
	      mass = mymalloc("mass", Ngroups * sizeof(float));
	      my_fread(mass, Ngroups, sizeof(float), fd);
	      for(i = 0; i < Ngroups; i++)
		Group[i].Mass = mass[i];
	      myfree(mass);

	      /* CM */
	      cm = mymalloc("cm", Ngroups * 3 * sizeof(float));
	      my_fread(cm, Ngroups, 3 * sizeof(float), fd);
	      for(i = 0; i < Ngroups; i++)
		for(j = 0; j < 3; j++)
		  Group[i].CM[j] = cm[i * 3 + j];
	      myfree(cm);

	      fclose(fd);

	      sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, ThisTask);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(1184132);
		}
	      my_fread(&Ngroups, sizeof(int), 1, fd);
	      my_fread(&TotNgroups, sizeof(int), 1, fd);
	      my_fread(&Nids, sizeof(int), 1, fd);
	      my_fread(&TotNids, sizeof(long long), 1, fd);
	      my_fread(&ntask, sizeof(int), 1, fd);
	      my_fread(&nid_previous, sizeof(int), 1, fd);	/* this is the number of IDs in previous files */

	      ID_list = mymalloc("ID_list", Nids * sizeof(fof_id_list));
	      ids = mymalloc("ids", Nids * sizeof(MyIDType));

	      my_fread(ids, sizeof(MyIDType), Nids, fd);

	      for(i = 0; i < Nids; i++)
		ID_list[i].ID = ids[i];

	      myfree(ids);

	      fclose(fd);
	    }

	  MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
	}
    }

  t1 = my_second();
  if(ThisTask == 0)
    printf("reading  took %g sec\n", timediff(t0, t1));



  t0 = my_second();

  /* now need to assign group numbers */


  list_of_ngroups = mymalloc("list_of_ngroups", NTask * sizeof(int));
  list_of_nids = mymalloc("list_of_nids", NTask * sizeof(int));

  MPI_Allgather(&Ngroups, 1, MPI_INT, list_of_ngroups, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&Nids, 1, MPI_INT, list_of_nids, 1, MPI_INT, MPI_COMM_WORLD);

  list_of_allgrouplen = mymalloc("list_of_allgrouplen", TotNgroups * sizeof(int));

  recvoffset = mymalloc("recvoffset", NTask * sizeof(int));
  for(i = 1, recvoffset[0] = 0; i < NTask; i++)
    recvoffset[i] = recvoffset[i - 1] + list_of_ngroups[i - 1];
  len = mymalloc("len", Ngroups * sizeof(int));
  for(i = 0; i < Ngroups; i++)
    len[i] = Group[i].Len;
  MPI_Allgatherv(len, Ngroups, MPI_INT, list_of_allgrouplen, list_of_ngroups, recvoffset, MPI_INT,
		 MPI_COMM_WORLD);
  myfree(len);
  myfree(recvoffset);

  /* do a check */
  long long totlen;

  for(i = 0, totlen = 0; i < TotNgroups; i++)
    totlen += list_of_allgrouplen[i];
  if(totlen != TotNids)
    endrun(8881);


  for(i = 0, count = 0; i < ThisTask; i++)
    count += list_of_ngroups[i];

  for(i = 0; i < Ngroups; i++)
    Group[i].GrNr = 1 + count + i;

  /* fix Group.Offset (may have overflown) */
  if(Ngroups > 0)
    for(i = 0, count = 0, Group[0].Offset = 0; i < ThisTask; i++)
      for(j = 0; j < list_of_ngroups[i]; j++)
	Group[0].Offset += list_of_allgrouplen[count++];

  for(i = 1; i < Ngroups; i++)
    Group[i].Offset = Group[i - 1].Offset + Group[i - 1].Len;


  long long *idoffset = mymalloc("idoffset", NTask * sizeof(long long));

  for(i = 1, idoffset[0] = 0; i < NTask; i++)
    idoffset[i] = idoffset[i - 1] + list_of_nids[i - 1];

  count = 0;

  for(i = 0, grnr = 1, totlen = 0; i < TotNgroups; totlen += list_of_allgrouplen[i++], grnr++)
    {
      if(totlen + list_of_allgrouplen[i] - 1 >= idoffset[ThisTask] && totlen < idoffset[ThisTask] + Nids)
	{
	  for(j = 0; j < list_of_allgrouplen[i]; j++)
	    {
	      if((totlen + j) >= idoffset[ThisTask] && (totlen + j) < (idoffset[ThisTask] + Nids))
		{
		  ID_list[(totlen + j) - idoffset[ThisTask]].GrNr = grnr;
		  count++;
		}
	    }
	}
    }
  if(count != Nids)
    endrun(1231);

  myfree(idoffset);

  t1 = my_second();
  if(ThisTask == 0)
    printf("assigning took %g sec\n", timediff(t0, t1));



  t0 = my_second();

  for(i = 0; i < NumPart; i++)
    P[i].SubNr = i;

  qsort(P, NumPart, sizeof(struct particle_data), io_compare_P_ID);
  qsort(ID_list, Nids, sizeof(fof_id_list), fof_compare_ID_list_ID);

  for(i = 0; i < NumPart; i++)
    P[i].GrNr = TotNgroups + 1;	/* will mark particles that are not in any group */

  t1 = my_second();
  if(ThisTask == 0)
    printf("sorting took %g sec\n", timediff(t0, t1));


  static fof_id_list *recv_ID_list;

  t0 = my_second();

  int matches = 0;

  /* exchange  data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(list_of_nids[sendTask] > 0 || list_of_nids[recvTask] > 0)
	    {
	      if(ngrp == 0)
		{
		  recv_ID_list = ID_list;
		}
	      else
		{
		  recv_ID_list = mymalloc("recv_ID_list", list_of_nids[recvTask] * sizeof(fof_id_list));

		  /* get the particles */
		  MPI_Sendrecv(ID_list, Nids * sizeof(fof_id_list), MPI_BYTE, recvTask, TAG_FOF_H,
			       recv_ID_list, list_of_nids[recvTask] * sizeof(fof_id_list), MPI_BYTE,
			       recvTask, TAG_FOF_H, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

	      for(i = 0, j = 0; i < list_of_nids[recvTask]; i++)
		{
		  while(j < NumPart - 1 && P[j].ID < recv_ID_list[i].ID)
		    j++;

		  if(recv_ID_list[i].ID == P[j].ID)
		    {
		      P[j].GrNr = recv_ID_list[i].GrNr;
		      matches++;
		    }
		}

	      if(ngrp != 0)
		myfree(recv_ID_list);
	    }
	}
    }

  sumup_large_ints(1, &matches, &totlen);
  if(totlen != TotNids)
    endrun(543);

  t1 = my_second();
  if(ThisTask == 0)
    printf("assigning GrNr to P[] took %g sec\n", timediff(t0, t1));

  MPI_Barrier(MPI_COMM_WORLD);

  myfree(list_of_allgrouplen);
  myfree(list_of_nids);
  myfree(list_of_ngroups);

  /* restore peano-hilbert order */
  qsort(P, NumPart, sizeof(struct particle_data), fof_compare_P_SubNr);
  subfind(num);
  endrun(0);

}



int fof_compare_ID_list_ID(const void *a, const void *b)
{
  if(((fof_id_list *) a)->ID < ((fof_id_list *) b)->ID) {return -1;}
  if(((fof_id_list *) a)->ID > ((fof_id_list *) b)->ID) {return +1;}
  return 0;
}


int fof_compare_P_SubNr(const void *a, const void *b)
{
  if(((struct particle_data *) a)->SubNr < (((struct particle_data *) b)->SubNr)) {return -1;}
  if(((struct particle_data *) a)->SubNr > (((struct particle_data *) b)->SubNr)) {return +1;}
  return 0;
}

#endif // IO_SUBFIND_READFOF_FROMIC

#endif /* of FOF */
