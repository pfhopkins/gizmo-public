#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

/*
* This file was originally part of the GADGET3 code developed by Volker Springel.
* It has been updated significantly by PFH for basic compatibility with GIZMO,
* as well as code cleanups, and accommodating new GIZMO functionality for various
* other operations. See GIZMO User Guide for discussion of the state of the code.
* Some parts have been heavily re-written by PFH, in particular the ADDIO flags, which
* are designed to more flexibly use the modern GIZMO structures and minimize the user
* having to add large amounts of redundant code for new operations. But most of the module
* remains in the GADGET3 format and parallelization style, with only basic modifications
* by PFH to remove bugs, dead code, correct incorrect function calls, update variable names
* to the modern GIZMO conventions, replace variables no longer used, rewrite macro syntax to
* GIZMO conventions, and ensure compilation without errors on modern compilers.
*/


/* by default, always use snapshot format for outputs, otherwise these get unwieldy very fast */
#if !defined(IO_SUBFIND_IN_OLD_ASCII_FORMAT)
#define SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
#endif
#if defined(HAVE_HDF5)&& defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
#include <hdf5.h>
#endif
#include "../../allvars.h"
#include "../../proto.h"
#include "../../domain.h"


#ifdef SUBFIND

#include "../fof.h"
#include "subfind.h"

static struct subfind_id_list
{
  MyIDType ID;
  int GrNr;
  int SubNr;
  float BindingEgy;
#ifdef SUBFIND_SAVE_PARTICLEDATA
  float Pos[3];
  float Vel[3];
  int Type;
    float Mass;
#endif
}
 *ID_list;

static int Nids;


void subfind(int num)
{
  double t0, t1, tstart, tend;
  int i, gr, nlocid, offset, limit, ncount, ntotingrouplocal, nminingrouplocal, nmaxingrouplocal;

#ifdef FOF_DENSITY_SPLIT_TYPES
  struct unbind_data *d;
  int j, n, count[6], countall[6];
#endif

  if(ThisTask == 0)
    printf("\nWe now execute a parallel version of SUBFIND.\n");

  tstart = my_second();

#ifdef FOF_DENSITY_SPLIT_TYPES
  for(j = 0; j < 6; j++)
    count[j] = 0;

  /* let's count number of particles of selected species */
  for(i = 0; i < NumPart; i++)
    count[P[i].Type]++;

  MPI_Allreduce(count, countall, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* do first loop: basically just defining the hsml for different species */
  for(j = 0; j < 6; j++)
    {
      if((1 << j) & (FOF_DENSITY_SPLIT_TYPES))
	{
	  //if(j == 5) {countall[j] = 0;}	/* this will prevent that the black holes are treated separately */

	  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

	  if(countall[j] > All.DesLinkNgb)
	    {
	      /* build index list of particles of selected species */
	      d = (struct unbind_data *) mymalloc("	      d", count[j] * sizeof(struct unbind_data));
	      for(i = 0, n = 0; i < NumPart; i++)
		if(P[i].Type == j)
		  d[n++].index = i;

	      t0 = my_second();
	      if(ThisTask == 0)
		printf("Tree construction for species %d (%d).\n", j, countall[j]);

	      CPU_Step[CPU_FOF] += measure_time();

	      force_treebuild(count[j], d);

	      myfree(d);

	      t1 = my_second();
	      if(ThisTask == 0)
		printf("tree build for species %d took %g sec\n", j, timediff(t0, t1));
	    }
	  else
	    {
	      t0 = my_second();
	      if(ThisTask == 0)
		printf("Tree construction.\n");

	      CPU_Step[CPU_FOF] += measure_time();

	      force_treebuild(NumPart, NULL);

	      t1 = my_second();
	      if(ThisTask == 0)
		printf("tree build took %g sec\n", timediff(t0, t1));
	    }


	  /* let's determine the local densities */
	  t0 = my_second();
	  subfind_setup_smoothinglengths(j);
	  subfind_density(j);
	  t1 = my_second();
	  if(ThisTask == 0)
	    printf("density and smoothing length for species %d took %g sec\n", j, timediff(t0, t1));

	  force_treefree();

	  /* let's save density contribution of own species */
	  for(i = 0; i < NumPart; i++)
	    if(P[i].Type == j)
	      P[i].w.density_sum = P[i].u.DM_Density;

	}
    }

  /* do second loop: now calculate all density contributions */
  for(j = 0; j < 6; j++)
    {
      if((1 << j) & (FOF_DENSITY_SPLIT_TYPES))
	{
	  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

	  /* build index list of particles of selectes species */
	  d = (struct unbind_data *) mymalloc("	  d", count[j] * sizeof(struct unbind_data));
	  for(i = 0, n = 0; i < NumPart; i++)
	    if(P[i].Type == j)
	      d[n++].index = i;

	  t0 = my_second();
	  if(ThisTask == 0)
	    printf("Tree construction for species %d (%d).\n", j, countall[j]);

	  CPU_Step[CPU_FOF] += measure_time();

	  force_treebuild(count[j], d);

	  myfree(d);

	  t1 = my_second();
	  if(ThisTask == 0)
	    printf("tree build for species %d took %g sec\n", j, timediff(t0, t1));

	  /* let's determine the local densities */
	  t0 = my_second();
	  for(i = 0; i < 6; i++)
	    if((1 << i) & (FOF_DENSITY_SPLIT_TYPES))
	      if(j != i)
		{
		  if(countall[i] > All.DesLinkNgb)
		    {
		      if(ThisTask == 0)
			printf("calculating density contribution of species %d to species %d\n", j, i);
		      subfind_density(-(i + 1));
		    }
		}
	  t1 = my_second();
	  if(ThisTask == 0)
	    printf("density() of species %d took %g sec\n", j, timediff(t0, t1));

	  force_treefree();

	  /* let's sum up density contribution */
	  for(i = 0; i < NumPart; i++)
	    if((1 << P[i].Type) & (FOF_DENSITY_SPLIT_TYPES))
	      if(j != P[i].Type)
		if(countall[P[i].Type] > All.DesLinkNgb)
		  P[i].w.density_sum += P[i].u.DM_Density;
	}
    }


  for(i = 0; i < NumPart; i++)
    {
      P[i].u.DM_Density = P[i].w.density_sum;

      if(P[i].Type == 0)
          P[i].w.int_energy = SphP[i].InternalEnergy;
      else
	P[i].w.int_energy = 0;

    }
#else
  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  t0 = my_second();
  if(ThisTask == 0)
    printf("Tree construction.\n");

  CPU_Step[CPU_FOF] += measure_time();

  force_treebuild(NumPart, NULL);

  t1 = my_second();
  if(ThisTask == 0)
    printf("tree build took %g sec\n", timediff(t0, t1));


  /* let's determine the local dark matter densities */
  t0 = my_second();
  subfind_setup_smoothinglengths(0);
  subfind_density(0);
  t1 = my_second();
  if(ThisTask == 0)
    printf("dark matter density() took %g sec\n", timediff(t0, t1));

  force_treefree();
#endif /* FOF_DENSITY_SPLIT_TYPES */

    {
      /* let's save the densities to a file (for making images) */
      t0 = my_second();
      subfind_save_densities(num);
      t1 = my_second();
      if(ThisTask == 0)
	printf("saving densities took %g sec\n", timediff(t0, t1));
    }

  /* count how many groups we have that should be done collectively */
  limit = (int)(0.6 * All.TotNumPart / NTask);

  for(i = 0, ncount = 0, ntotingrouplocal = 0; i < Ngroups; i++)
    if(Group[i].Len >= limit)
      ncount++;
    else
      ntotingrouplocal += Group[i].Len;

  MPI_Allreduce(&ncount, &Ncollective, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ntotingrouplocal, &nminingrouplocal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&ntotingrouplocal, &nmaxingrouplocal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("\nNumber of FOF halos treated with collective SubFind code = %d\n", Ncollective);
      printf("(the adopted size-limit for the collective algorithm was %d particles.)\n", limit);
      printf("the other %d FOF halos are treated in parallel with serial code\n\n", TotNgroups - Ncollective);
      printf("Unbalance in total number of particles in FOF halos is %d - %d \n\n", nminingrouplocal,
	     nmaxingrouplocal);
    }

  /*  to decide on which task a group should be:
   *  if   GrNr <= Ncollective:  collective groupfinding.
   *  the task where the group info is put is TaskNr = (GrNr - 1) % NTask
   */

  /* now we distribute the particles such that small groups are assigned in
   *  total to certain CPUs, and big groups are left where they are 
   */

  t0 = my_second();

  for(i = 0; i < NumPart; i++)
    {
      P[i].origintask2 = ThisTask;

      if(P[i].GrNr > Ncollective && P[i].GrNr <= TotNgroups)	/* particle is in small group */
	P[i].targettask = (P[i].GrNr - 1) % NTask;
      else
	P[i].targettask = ThisTask;
    }

  subfind_exchange();		/* distributes gas particles as well if needed */

  t1 = my_second();
  if(ThisTask == 0)
    printf("subfind_exchange()() took %g sec\n", timediff(t0, t1));

  subfind_distribute_groups();

  qsort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);

  for(i = 0; i < NumPart; i++)
    if(P[i].GrNr > Ncollective && P[i].GrNr <= TotNgroups)
      if(((P[i].GrNr - 1) % NTask) != ThisTask)
	{
	  printf("i=%d %d task=%d type=%d\n", i, P[i].GrNr, ThisTask, P[i].Type);
	  endrun(87);
	}

  /* lets estimate the maximum number of substructures we need to store on the local CPU */
  for(i = 0, nlocid = 0; i < Ngroups; i++)
    nlocid += Group[i].Len;

  MaxNsubgroups = nlocid / All.DesLinkNgb + NTask;	/* this is a quite conservative upper limit */
  Nsubgroups = 0;
  SubGroup =
    (struct subgroup_properties *) mymalloc("SubGroup", MaxNsubgroups * sizeof(struct subgroup_properties));

  for(i = 0; i < NumPart; i++)
    P[i].SubNr = (1 << 30);	/* default */

  /* we begin by applying the collective version of subfind to distributed groups */
  t0 = my_second();
  for(GrNr = 1; GrNr <= Ncollective; GrNr++)
    subfind_process_group_collectively(num);
  t1 = my_second();
  if(ThisTask == 0)
    printf("processing of collective halos took %g sec\n", timediff(t0, t1));

  for(i = 0; i < NumPart; i++)
    {
      P[i].origindex = i;
      P[i].origintask = ThisTask;
    }

  t0 = my_second();
  qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_GrNr_DM_Density);
  t1 = my_second();
  if(ThisTask == 0)
    printf("sort of local particles()() took %g sec\n", timediff(t0, t1));


  /* now we have the particles of groups consecutively, but SPH particles are
     not aligned. They can however be accessed via SphP[P[i].originindex] */


  /* let's count how many local particles we have in small groups */
  for(i = 0, nlocid = 0; i < NumPart; i++)
    if(P[i].GrNr > Ncollective && P[i].GrNr <= Ngroups)	/* particle is in small group */
      nlocid++;

  if(ThisTask == 0)
    printf("contructing tree for serial subfind of local groups\n");

  subfind_loctree_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  if(ThisTask == 0)
    printf("Start to do local groups with serial subfind algorithm\n");

  t0 = my_second();

  /* we now apply a serial version of subfind to the local groups */
  for(gr = 0, offset = 0; gr < Ngroups; gr++)
    {
      if(Group[gr].GrNr > Ncollective)
	{
	  if(((Group[gr].GrNr - 1) % NTask) == ThisTask)
	    offset = subfind_process_group_serial(gr, offset);
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);

  t1 = my_second();
  if(ThisTask == 0)
    printf("\nprocessing of local groups took took %g sec\n\n", timediff(t0, t1));


  subfind_loctree_treefree();


  /* bringing back particles in original positions, such that gas particles are aligned */
  t0 = my_second();
  qsort(P, NumPart, sizeof(struct particle_data), subfind_compare_P_origindex);
  t1 = my_second();
  if(ThisTask == 0)
    printf("unsorting of local particles()() took %g sec\n", timediff(t0, t1));


  GrNr = -1;			/* to ensure that domain decomposition acts normally again */

  /* now determine the remaining spherical overdensity values for the non-local groups */

  domain_free_trick();

  CPU_Step[CPU_FOF] += measure_time();


#ifdef FOF_DENSITY_SPLIT_TYPES
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].origintask != ThisTask)
	printf("Task %d: Holding particle of task %d !\n", ThisTask, P[i].origintask);
      if(P[i].origindex != i)
	printf("Task %d: Particles is in wrong position (is=%d, was=%d) !\n", ThisTask, i, P[i].origindex);
    }
#endif



  t0 = my_second();

  for(i = 0; i < NumPart; i++)
    P[i].targettask = P[i].origintask2;

  subfind_exchange();		/* distributes gas particles as well if needed */

  t1 = my_second();
  if(ThisTask == 0)
    printf("subfind_exchange() (for return to original CPU)  took %g sec\n", timediff(t0, t1));

  domain_Decomposition(1, 0, 0);

  force_treebuild(NumPart, NULL);


  /* compute spherical overdensities for FOF groups */
  t0 = my_second();

  Subfind_DensityOtherProps_Loop();

  t1 = my_second();
  if(ThisTask == 0)
    printf("determining spherical overdensity masses took %g sec\n", timediff(t0, t1));


  /* determine which halos are contaminated by boundary particles */
  t0 = my_second();

  subfind_contamination();

  t1 = my_second();
  if(ThisTask == 0)
    printf("determining contamination of halos took %g sec\n", timediff(t0, t1));


  force_treefree();
  domain_free();

  domain_allocate_trick();

  /* now assemble final output */
  subfind_save_final(num);

  tend = my_second();

  if(ThisTask == 0)
    printf("\nFinished with SUBFIND.  (total time=%g sec)\n\n", timediff(tstart, tend));

  myfree(SubGroup);

  CPU_Step[CPU_FOF] += measure_time();
}




void subfind_save_final(int num)
{
  int i, j, totsubs, primaryTask, groupTask, nprocgroup;
  char buf[1000];
  double t0, t1;

  /* prepare list of ids with assigned group numbers */

  t0 = my_second();

  t0 = my_second();
  parallel_sort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_GrNr);
  t1 = my_second();
  if(ThisTask == 0)
    {
      printf("Global sort of Groups took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  t0 = my_second();
  parallel_sort(SubGroup, Nsubgroups, sizeof(struct subgroup_properties),subfind_compare_SubGroup_GrNr_SubNr);
  t1 = my_second();
  if(ThisTask == 0)
    {
      printf("Global sort of SubGroups took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  ID_list = (struct subfind_id_list *)mymalloc("ID_list", sizeof(struct subfind_id_list) * (int)(All.MaxPart * 1.5));

  for(i = 0, Nids = 0; i < NumPart; i++)
    {
      if(P[i].GrNr <= TotNgroups)
	{
	  ID_list[Nids].GrNr = P[i].GrNr;
	  ID_list[Nids].SubNr = P[i].SubNr;
	  ID_list[Nids].BindingEgy = P[i].v.DM_BindingEnergy;
	  ID_list[Nids].ID = P[i].ID;
#ifdef SUBFIND_SAVE_PARTICLEDATA
	  for(j = 0; j < 3; j++)
	    {
	      ID_list[Nids].Pos[j] = P[i].Pos[j];
	      ID_list[Nids].Vel[j] = P[i].Vel[j];
	    }
	  ID_list[Nids].Type = P[i].Type;
	  ID_list[Nids].Mass = P[i].Mass;
#endif
	  Nids++;
	}
    }

  t0 = my_second();

  parallel_sort(ID_list, Nids, sizeof(struct subfind_id_list), subfind_compare_ID_list);
  t1 = my_second();
  if(ThisTask == 0)
    {
      printf("Global sort of IDs took %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  MPI_Allreduce(&Nsubgroups, &TotNsubgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  /* fill in the FirstSub-values */
  for(i = 0, totsubs = 0; i < Ngroups; i++)
    {
      if(i > 0)
	Group[i].FirstSub = Group[i - 1].FirstSub + Group[i - 1].Nsubs;
      else
	Group[i].FirstSub = 0;
      totsubs += Group[i].Nsubs;
    }

  MPI_Allgather(&totsubs, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  for(i = 0; i < Ngroups; i++)
    Group[i].FirstSub += Send_offset[ThisTask];


  MPI_Allgather(&Nids, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);
  for(j = 1, Send_offset[0] = 0; j < NTask; j++)
    Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);


  if(NTask < All.NumFilesWrittenInParallel)
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(241931);
    }

  t0 = my_second();

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  primaryTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (primaryTask + groupTask))	/* ok, it's this processor's turn */
	subfind_save_local_catalogue(num);
      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }

  t1 = my_second();

  if(ThisTask == 0)
    {
      printf("Subgroup catalogues saved. took = %g sec\n", timediff(t0, t1));
      fflush(stdout);
    }

  myfree(ID_list);
}


int get_sub_entrytype_of_block(enum siofields blocknr)
{
  int type = 0;

  switch (blocknr)
    {
    case SIO_SLEN:
    case SIO_SOFF:
    case SIO_PFOF:
    case SIO_MSUB:
    case SIO_SPOS:
    case SIO_SVEL:
    case SIO_SCM:
    case SIO_SPIN:
    case SIO_DSUB:
    case SIO_VMAX:
    case SIO_RVMAX:
    case SIO_RHMS:
    case SIO_MBID:
    case SIO_GRNR:
    case SIO_SMST:
    case SIO_SLUM:
    case SIO_SLATT:
    case SIO_SLOBS:
    case SIO_DUST:
    case SIO_SAGE:
    case SIO_SZ:
    case SIO_SSFR:
      type = 1;
      break;
    case SIO_PPOS:
    case SIO_PVEL:
    case SIO_PTYP:
    case SIO_PMAS:
    case SIO_PID:
      type = 2;
      break;
    case SIO_BGPOS:
    case SIO_BGMTOP:
    case SIO_BGRTOP:
      type = 3;
      break;
    default:
      type = 0;
      break;
    }

  return type;
}

int get_values_per_sub(enum siofields blocknr)
{
  int n = 1;

  switch (blocknr)
    {
    case SIO_GPOS:
    case SIO_SPOS:
    case SIO_SVEL:
    case SIO_SCM:
    case SIO_SPIN:
    case SIO_PPOS:
    case SIO_PVEL:
    case SIO_BGPOS:
      n = 3;
      break;
    case SIO_SMST:
      n = 6;
      break;
    case SIO_SLUM:
    case SIO_SLATT:
    case SIO_SLOBS:
      n = 1;
      break;
    case SIO_DUST:
      n = 1;
      break;
    case SIO_DELTA_MSUB:
    case SIO_DELTA_RSUB:
    case SIO_DELTA_DISPSUB:
    case SIO_DELTA_MGASSUB:
    case SIO_DELTA_MSTSUB:
    case SIO_DELTA_TEMPSUB:
    case SIO_DELTA_LXSUB:
      n = SUBFIND_ADDIO_NUMOVERDEN;
      break;
    default:
      n = 1;
      break;
    }

  return n;
}

#ifdef HAVE_HDF5
void get_IO_Label_HDF5_sub(enum siofields blocknr, char *label)
{
  switch (blocknr)
    {
    case SIO_GLEN:
      strcpy(label, "GroupLen");
      break;
    case SIO_GOFF:
      strcpy(label, "GroupOffset");
      break;
    case SIO_MTOT:
      strcpy(label, "GroupMass");
      break;
    case SIO_GPOS:
      strcpy(label, "GroupPos");
      break;
    case SIO_DELTA_MSUB:
      strcpy(label, "Group_vsDelta_Mass");
      break;
    case SIO_DELTA_RSUB:
      strcpy(label, "Group_vsDelta_Radius");
      break;
    case SIO_DELTA_DISPSUB:
      strcpy(label, "Group_vsDelta_1DVelDispersion");
      break;
    case SIO_DELTA_MGASSUB:
      strcpy(label, "Group_vsDelta_GasMass");
      break;
    case SIO_DELTA_MSTSUB:
      strcpy(label, "Group_vsDelta_StellarMass");
      break;
    case SIO_DELTA_TEMPSUB:
      strcpy(label, "Group_vsDelta_MeanGasTempInKeV");
      break;
    case SIO_DELTA_LXSUB:
      strcpy(label, "Group_vsDelta_LXrayGasIn1e44");
      break;
    case SIO_NCON:
      strcpy(label, "GroupContaminationCoun");
      break;
    case SIO_MCON:
      strcpy(label, "GroupContaminationMass");
      break;
    case SIO_BGPOS:
      strcpy(label, "BGPO");
      break;
    case SIO_BGMTOP:
      strcpy(label, "BGMA");
      break;
    case SIO_BGRTOP:
      strcpy(label, "BGRA");
      break;
    case SIO_NSUB:
      strcpy(label, "GroupNsubs");
      break;
    case SIO_FSUB:
      strcpy(label, "GroupFirstSub");
      break;
    case SIO_SLEN:
      strcpy(label, "SubhaloLen");
      break;
    case SIO_SOFF:
      strcpy(label, "SubhaloOffset");
      break;
    case SIO_PFOF:
      strcpy(label, "SubhaloParent");
      break;
    case SIO_MSUB:
      strcpy(label, "SubhaloMass");
      break;
    case SIO_SPOS:
      strcpy(label, "SubhaloPos");
      break;
    case SIO_SVEL:
      strcpy(label, "SubhaloVel");
      break;
    case SIO_SCM:
      strcpy(label, "SubhaloCM");
      break;
    case SIO_SPIN:
      strcpy(label, "SubhaloSpin");
      break;
    case SIO_DSUB:
      strcpy(label, "SubhaloVelDisp");
      break;
    case SIO_VMAX:
      strcpy(label, "SubhaloVmax");
      break;
    case SIO_RVMAX:
      strcpy(label, "SubhaloVmaxRad");
      break;
    case SIO_RHMS:
      strcpy(label, "SubhaloHalfmassRad");
      break;
    case SIO_MBID:
      strcpy(label, "SubhaloIDMostbound");
      break;
    case SIO_GRNR:
      strcpy(label, "SubhaloGrNr");
      break;
    case SIO_SMST:
      strcpy(label, "SMST");
      break;
    case SIO_SLUM:
      strcpy(label, "SLUM");
      break;
    case SIO_SLATT:
      strcpy(label, "SLAT");
      break;
    case SIO_SLOBS:
      strcpy(label, "SLOB");
      break;
    case SIO_DUST:
      strcpy(label, "DUST");
      break;
    case SIO_SAGE:
      strcpy(label, "SAGE");
      break;
    case SIO_SZ:
      strcpy(label, "SZ  ");
      break;
    case SIO_SSFR:
      strcpy(label, "SSFR");
      break;
    case SIO_PPOS:
      strcpy(label, "PPOS");
      break;
    case SIO_PVEL:
      strcpy(label, "PVEL");
      break;
    case SIO_PTYP:
      strcpy(label, "PTYP");
      break;
    case SIO_PMAS:
      strcpy(label, "PMAS");
      break;
    case SIO_PID:
      strcpy(label, "PID ");
      break;
    default:
      endrun(987453);
      break;
    }
}
#endif

void get_IO_Label_sub(enum siofields blocknr, char *label)
{
  switch (blocknr)
    {
    case SIO_GLEN:
      strncpy(label, "GLEN", 4);
      break;
    case SIO_GOFF:
      strncpy(label, "GOFF", 4);
      break;
    case SIO_MTOT:
      strncpy(label, "MTOT", 4);
      break;
    case SIO_GPOS:
      strncpy(label, "GPOS", 4);
      break;
    case SIO_DELTA_MSUB:
      strncpy(label, "DMSB", 4);
      break;
    case SIO_DELTA_RSUB:
      strncpy(label, "DRSB", 4);
      break;
    case SIO_DELTA_DISPSUB:
      strncpy(label, "DDSB", 4);
      break;
    case SIO_DELTA_MGASSUB:
      strncpy(label, "DMGS", 4);
      break;
    case SIO_DELTA_MSTSUB:
      strncpy(label, "DMST", 4);
      break;
    case SIO_DELTA_TEMPSUB:
      strncpy(label, "DTMP", 4);
      break;
    case SIO_DELTA_LXSUB:
      strncpy(label, "DLXS", 4);
      break;
    case SIO_NCON:
      strncpy(label, "NCON", 4);
      break;
    case SIO_MCON:
      strncpy(label, "MCON", 4);
      break;
    case SIO_BGPOS:
      strncpy(label, "BGPO", 4);
      break;
    case SIO_BGMTOP:
      strncpy(label, "BGMA", 4);
      break;
    case SIO_BGRTOP:
      strncpy(label, "BGRA", 4);
      break;
    case SIO_NSUB:
      strncpy(label, "NSUB", 4);
      break;
    case SIO_FSUB:
      strncpy(label, "FSUB", 4);
      break;
    case SIO_SLEN:
      strncpy(label, "SLEN", 4);
      break;
    case SIO_SOFF:
      strncpy(label, "SOFF", 4);
      break;
    case SIO_PFOF:
      strncpy(label, "SSUB", 4);
      break;
    case SIO_MSUB:
      strncpy(label, "MSUB", 4);
      break;
    case SIO_SPOS:
      strncpy(label, "SPOS", 4);
      break;
    case SIO_SVEL:
      strncpy(label, "SVEL", 4);
      break;
    case SIO_SCM:
      strncpy(label, "SCM ", 4);
      break;
    case SIO_SPIN:
      strncpy(label, "SPIN", 4);
      break;
    case SIO_DSUB:
      strncpy(label, "DSUB", 4);
      break;
    case SIO_VMAX:
      strncpy(label, "VMAX", 4);
      break;
    case SIO_RVMAX:
      strncpy(label, "RMAX", 4);
      break;
    case SIO_RHMS:
      strncpy(label, "RHMS", 4);
      break;
    case SIO_MBID:
      strncpy(label, "MBID", 4);
      break;
    case SIO_GRNR:
      strncpy(label, "GRNR", 4);
      break;
    case SIO_SMST:
      strncpy(label, "SMST", 4);
      break;
    case SIO_SLUM:
      strncpy(label, "SLUM", 4);
      break;
    case SIO_SLATT:
      strncpy(label, "SLAT", 4);
      break;
    case SIO_SLOBS:
      strncpy(label, "SLOB", 4);
      break;
    case SIO_DUST:
      strncpy(label, "DUST", 4);
      break;
    case SIO_SAGE:
      strncpy(label, "SAGE", 4);
      break;
    case SIO_SZ:
      strncpy(label, "SZ  ", 4);
      break;
    case SIO_SSFR:
      strncpy(label, "SSFR", 4);
      break;
    case SIO_PPOS:
      strncpy(label, "PPOS", 4);
      break;
    case SIO_PVEL:
      strncpy(label, "PVEL", 4);
      break;
    case SIO_PTYP:
      strncpy(label, "PTYP", 4);
      break;
    case SIO_PMAS:
      strncpy(label, "PMAS", 4);
      break;
    case SIO_PID:
      strncpy(label, "PID ", 4);
      break;
    default:
      endrun(987453);
      break;
    }
}

int get_datatype_in_sub(enum siofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case SIO_MBID:
    case SIO_PID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;

    case SIO_GLEN:
    case SIO_GOFF:
    case SIO_NCON:
    case SIO_NSUB:
    case SIO_FSUB:
    case SIO_SLEN:
    case SIO_SOFF:
    case SIO_PFOF:
    case SIO_GRNR:
    case SIO_PTYP:
      typekey = 0;		/* native int */
      break;

    default:
#ifdef OUTPUT_IN_DOUBLEPRECISION
      typekey = 3;
#else
      typekey = 1;		/* native MyOutputFloat */
#endif
      break;
    }

  return typekey;
}

int block_in_sub(enum siofields blocknr)
{
  int present = 0;

  switch (blocknr)
    {
    case SIO_GLEN:
    case SIO_GOFF:
    case SIO_MTOT:
    case SIO_GPOS:
    case SIO_DELTA_MSUB:
    case SIO_DELTA_RSUB:
#ifdef SUBFIND_ADDIO_VELDISP
    case SIO_DELTA_DISPSUB:
#endif
#ifdef SUBFIND_ADDIO_BARYONS
    case SIO_DELTA_MGASSUB:
    case SIO_DELTA_MSTSUB:
    case SIO_DELTA_TEMPSUB:
    case SIO_DELTA_LXSUB:
#endif
    case SIO_NCON:
    case SIO_MCON:
    case SIO_NSUB:
    case SIO_FSUB:
    case SIO_SLEN:
    case SIO_SOFF:
    case SIO_PFOF:
    case SIO_MSUB:
    case SIO_SPOS:
    case SIO_SVEL:
    case SIO_SCM:
    case SIO_SPIN:
    case SIO_DSUB:
    case SIO_VMAX:
    case SIO_RVMAX:
    case SIO_RHMS:
    case SIO_MBID:
    case SIO_GRNR:
    case SIO_SMST:
#ifdef SUBFIND_SAVE_PARTICLEDATA
    case SIO_PPOS:
    case SIO_PVEL:
    case SIO_PTYP:
    case SIO_PMAS:
#endif
    case SIO_PID:
      present = 1;
      break;
    default:
      present = 0;
      break;
    }

  return present;
}


void subfind_save_local_catalogue(int num)
{
  FILE *fd;
  char buf[500], fname[500], label[] = "--------";

  void *IOBuffer;
  MyOutputFloat *fp;
  int *fp_int, nwrite, ndim, bytes_per_blockelement, bnr, type, datatype, bytes;
  MyIDType *fp_id;
  enum siofields blocknr;
  unsigned int blksize;

  int i, j;

#ifdef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
  unsigned int nextblock;
#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[3], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank = 0;
#endif

#endif

#ifdef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
  if(NTask == 1)
    sprintf(fname, "%s/groups_%03d/%s_%03d", All.OutputDir, num, "sub", num);
  else
    sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "sub", num, ThisTask);
#else
  sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, ThisTask);
#endif
  strcpy(buf, fname);

#if defined(HAVE_HDF5)&& defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
  if(All.SnapFormat == 3)
    {

      sprintf(buf, "%s.hdf5", fname);
      hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

      hdf5_grp[0] = H5Gcreate(hdf5_file, "/Group", 0);
      hdf5_grp[1] = H5Gcreate(hdf5_file, "/Subhalo", 0);
      hdf5_grp[2] = H5Gcreate(hdf5_file, "/IDs", 0);



    }
  else
    {
#endif
      if(!(fd = fopen(buf, "w")))

	{
	  printf("can't open file `%s`\n", buf);
	  endrun(1183);
	}
#if defined(HAVE_HDF5)&& defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
    }
#endif

#ifdef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
  for(i = 0; i < 6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.npartTotalHighWord[i] = 0;
      header.mass[i] = 0;
    }
  header.npart[0] = Ngroups;
  header.npartTotal[0] = TotNgroups;
  header.npart[1] = Nsubgroups;
  header.npartTotal[1] = TotNsubgroups;
  header.npart[2] = Nids;
  header.npartTotal[2] = (unsigned int) TotNids;
  header.npartTotalHighWord[2] = (unsigned int) (TotNids >> 32);

  header.time = All.Time;
  if(All.ComovingIntegrationOn) {header.redshift = 1.0 / All.Time - 1;}
    else {header.redshift = 0;}

#ifdef COOLING
  header.flag_cooling = 1;
#endif
#ifdef GALSF
  header.flag_sfr = 1;
  header.flag_feedback = 1;
  header.flag_stellarage = 1;
#ifdef METALS
  header.flag_metals = 1;
#endif
#endif

  header.num_files = NTask;
  header.BoxSize = All.BoxSize;
  header.OmegaMatter = All.OmegaMatter;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

#ifdef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
  if(All.SnapFormat == 2)
    {
      blksize = sizeof(int) + 4 * sizeof(char);
      SKIP;
      my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
      nextblock = sizeof(header) + 2 * sizeof(int);
      my_fwrite(&nextblock, sizeof(int), 1, fd);
      SKIP;
    }
#ifdef HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      write_header_attributes_in_hdf5(hdf5_headergrp);
    }
  else
    {
#endif
#endif
      blksize = sizeof(header);
      SKIP;
      my_fwrite(&header, sizeof(header), 1, fd);
      SKIP;
#if defined(HAVE_HDF5) &&  defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
    }
#endif


#ifdef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
  if(All.SnapFormat == 2)
    {
      if(!(InfoBlock = (struct info_block *) mymalloc("InfoBlock", bytes = sizeof(struct info_block) * 1000)))
	{
	  printf("failed to allocate memory for `InfoBlock' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}

      int n_info;
      for(bnr = 0, n_info = 0; bnr < 1000; bnr++)
	{
	  blocknr = (enum siofields) bnr;

	  if(blocknr == SIO_LASTENTRY)
	    break;

	  if(block_in_sub(blocknr))
	    {
	      type = get_sub_entrytype_of_block(blocknr);

	      for(i = 0; i < 6; i++)
		InfoBlock[n_info].is_present[i] = 0;
	      InfoBlock[n_info].is_present[type] = 1;

	      InfoBlock[n_info].ndim = get_values_per_sub(blocknr);

	      get_IO_Label_sub(blocknr, label);
	      for(i = 0; i < 4; i++)
		InfoBlock[n_info].label[i] = label[i];

	      datatype = get_datatype_in_sub(blocknr);
	      switch (datatype)
		{
		case 0:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "LONG    ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "LONGN   ", 8);
		  break;
		case 1:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
		  break;
		case 2:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "LLONG   ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "LLONGN  ", 8);
		  break;
		case 3:
		  if(InfoBlock[n_info].ndim <= 1)
		    strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
		  else
		    strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
		  break;
		default:
		  endrun(8712519);
		  break;
		}
	      n_info++;
	    }
	}

      blksize = sizeof(int) + 4 * sizeof(char);
      SKIP;
      my_fwrite((void *) "INFO", sizeof(char), 4, fd);
      nextblock = n_info * sizeof(struct info_block) + 2 * sizeof(int);
      my_fwrite(&nextblock, sizeof(int), 1, fd);
      SKIP;
      blksize = n_info * sizeof(struct info_block);
      SKIP;
      my_fwrite(InfoBlock, sizeof(struct info_block), n_info, fd);
      SKIP;

      myfree(InfoBlock);
    }
#endif

#else
  my_fwrite(&Ngroups, sizeof(int), 1, fd);
  my_fwrite(&TotNgroups, sizeof(int), 1, fd);
  my_fwrite(&Nids, sizeof(int), 1, fd);
  my_fwrite(&TotNids, sizeof(long long), 1, fd);
  my_fwrite(&NTask, sizeof(int), 1, fd);
  my_fwrite(&Nsubgroups, sizeof(int), 1, fd);
  my_fwrite(&TotNsubgroups, sizeof(int), 1, fd);
#endif


  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum siofields) bnr;

      if(blocknr == SIO_LASTENTRY)
	break;

      if(block_in_sub(blocknr))
	{
	  type = get_sub_entrytype_of_block(blocknr);
	  switch (type)
	    {
	    case 0:
	      nwrite = Ngroups;
	      break;
	    case 1:
	      nwrite = Nsubgroups;
	      break;
	    case 2:
	      nwrite = Nids;
	      break;
	    }

	  ndim = get_values_per_sub(blocknr);

#if defined(HAVE_HDF5) &&  defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
	  if(All.SnapFormat == 3)
	    {
	      dims[0] = nwrite;
	      dims[1] = ndim;
	      if(dims[1] == 1) {rank = 1;} else {rank = 2;}

	      switch (get_datatype_in_sub(blocknr))
		{
		case 0:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
		  break;
		case 1:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
		  break;
		case 2:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
		  break;
		case 3:
		  hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
		  break;
		default:
		  endrun(8712519);
		  break;
		}
	      get_IO_Label_HDF5_sub(blocknr, buf);
	      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
	      hdf5_dataset =
		H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
	    }
#endif

	  datatype = get_datatype_in_sub(blocknr);
	  switch (datatype)
	    {
	    case 0:
	      bytes_per_blockelement = 4;
	      break;
	    case 1:
	      bytes_per_blockelement = 4;
	      break;
	    case 2:
	      bytes_per_blockelement = 8;
	      break;
	    case 3:
	      bytes_per_blockelement = 8;
	      break;
	    default:
	      endrun(8712519);
	      break;
	    }

#ifdef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
	  get_IO_Label_sub(blocknr, label);

	  if(nwrite > 0)
	    {
	      if(All.SnapFormat == 2)
		{
		  blksize = sizeof(int) + 4 * sizeof(char);
		  SKIP;
		  my_fwrite((void *) label, sizeof(char), 4, fd);
		  nextblock = nwrite * bytes_per_blockelement * ndim + 2 * sizeof(int);
		  my_fwrite(&nextblock, sizeof(int), 1, fd);
		  SKIP;
		}
	      blksize = nwrite * bytes_per_blockelement * ndim;
	      if(All.SnapFormat != 3)
		{
		  SKIP;
		}

	    }
	  else
	    blksize = 0;
#else
	  blksize = nwrite * bytes_per_blockelement * ndim;
#endif
	  if(ThisTask == 0)
	    printf("Writing block %d (%c%c%c%c), n=%d, ptype=%d, dtype=%d, ndim=%d, bpb=%d bytes=%ud\n", bnr,
		   label[0], label[1], label[2], label[3], nwrite, type, datatype, ndim,
		   bytes_per_blockelement, blksize);

	  if(!(IOBuffer = mymalloc("IOBuffer", bytes = blksize)))
	    {
	      printf("failed to allocate memory for `IOBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	      endrun(2);
	    }

	  fp = (MyOutputFloat *) IOBuffer;
	  fp_int = (int *) IOBuffer;
	  fp_id = (MyIDType *) IOBuffer;

	  switch (blocknr)
	    {
	    case SIO_GLEN:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = Group[i].Len;
	      break;
	    case SIO_GOFF:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = Group[i].Offset;
	      break;
	    case SIO_MTOT:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) Group[i].Mass;
	      break;
	    case SIO_GPOS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) Group[i].Pos[j];
	      break;
	    case SIO_DELTA_MSUB:
          for(i=0;i<nwrite;i++) {for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {fp[i*SUBFIND_ADDIO_NUMOVERDEN + j]=(MyOutputFloat)Group[i].SubHaloProps_vsDelta[j].M200;}}
	      break;
        case SIO_DELTA_RSUB:
          for(i=0;i<nwrite;i++) {for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {fp[i*SUBFIND_ADDIO_NUMOVERDEN + j]=(MyOutputFloat)Group[i].SubHaloProps_vsDelta[j].R200;}}
          break;
#ifdef SUBFIND_ADDIO_VELDISP
        case SIO_DELTA_DISPSUB:
          for(i=0;i<nwrite;i++) {for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {fp[i*SUBFIND_ADDIO_NUMOVERDEN + j]=(MyOutputFloat)Group[i].SubHaloProps_vsDelta[j].Disp200;}}
          break;
#endif
#ifdef SUBFIND_ADDIO_BARYONS
        case SIO_DELTA_MGASSUB:
          for(i=0;i<nwrite;i++) {for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {fp[i*SUBFIND_ADDIO_NUMOVERDEN + j]=(MyOutputFloat)Group[i].SubHaloProps_vsDelta[j].gas_mass;}}
          break;
        case SIO_DELTA_MSTSUB:
          for(i=0;i<nwrite;i++) {for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {fp[i*SUBFIND_ADDIO_NUMOVERDEN + j]=(MyOutputFloat)Group[i].SubHaloProps_vsDelta[j].star_mass;}}
          break;
        case SIO_DELTA_TEMPSUB:
          for(i=0;i<nwrite;i++) {for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {fp[i*SUBFIND_ADDIO_NUMOVERDEN + j]=(MyOutputFloat)Group[i].SubHaloProps_vsDelta[j].temp;}}
          break;
        case SIO_DELTA_LXSUB:
          for(i=0;i<nwrite;i++) {for(j=0;j<SUBFIND_ADDIO_NUMOVERDEN;j++) {fp[i*SUBFIND_ADDIO_NUMOVERDEN + j]=(MyOutputFloat)Group[i].SubHaloProps_vsDelta[j].xlum;}}
          break;
#endif
	    case SIO_NCON:
	      for(i = 0; i < nwrite; i++) fp_int[i] = Group[i].ContaminationLen;
	      break;
	    case SIO_MCON:
	      for(i = 0; i < nwrite; i++) fp[i] = (MyOutputFloat) Group[i].ContaminationMass;
	      break;
	    case SIO_NSUB:
	      for(i = 0; i < nwrite; i++) fp_int[i] = Group[i].Nsubs;
	      break;
	    case SIO_FSUB:
	      for(i = 0; i < nwrite; i++) fp_int[i] = Group[i].FirstSub;
	      break;
	    case SIO_SLEN:
	      for(i = 0; i < nwrite; i++) fp_int[i] = SubGroup[i].Len;
	      break;
	    case SIO_SOFF:
	      for(i = 0; i < nwrite; i++) fp_int[i] = SubGroup[i].Offset;
	      break;
	    case SIO_PFOF:
	      for(i = 0; i < nwrite; i++) fp_int[i] = SubGroup[i].SubParent;
	      break;
	    case SIO_MSUB:
	      for(i = 0; i < nwrite; i++) fp[i] = (MyOutputFloat) SubGroup[i].Mass;
	      break;
	    case SIO_SPOS:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].Pos[j];
	      break;
	    case SIO_SVEL:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].Vel[j];
	      break;
	    case SIO_SCM:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].CM[j];
	      break;
	    case SIO_SPIN:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) SubGroup[i].Spin[j];
	      break;
	    case SIO_DSUB:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubVelDisp;
	      break;
	    case SIO_VMAX:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubVmax;
	      break;
	    case SIO_RVMAX:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubVmaxRad;
	      break;
	    case SIO_RHMS:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) SubGroup[i].SubHalfMass;
	      break;
	    case SIO_MBID:
	      for(i = 0; i < nwrite; i++)
		fp_id[i] = SubGroup[i].SubMostBoundID;
	      break;
	    case SIO_GRNR:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = SubGroup[i].GrNr;
	      break;
	    case SIO_SMST:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 6; j++)
		  fp[i * 6 + j] = (MyOutputFloat) SubGroup[i].MassTab[j];
	      break;
#ifdef SUBFIND_SAVE_PARTICLEDATA
	    case SIO_PPOS:
#ifndef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT	/* open new file in case of old format */
	      fclose(fd);
	      sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_posvel", num, ThisTask);
	      if(!(fd = fopen(buf, "w")))
		{
		  printf("can't open file `%s`\n", buf);
		  endrun(1184);
		}
	      my_fwrite(&Ngroups, sizeof(int), 1, fd);
	      my_fwrite(&TotNgroups, sizeof(int), 1, fd);
	      my_fwrite(&Nids, sizeof(int), 1, fd);
	      my_fwrite(&TotNids, sizeof(long long), 1, fd);
	      my_fwrite(&NTask, sizeof(int), 1, fd);
	      my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
	      my_fwrite(&All.Time, sizeof(double), 1, fd);
#endif
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) ID_list[i].Pos[j];
	      break;
	    case SIO_PVEL:
	      for(i = 0; i < nwrite; i++)
		for(j = 0; j < 3; j++)
		  fp[i * 3 + j] = (MyOutputFloat) (ID_list[i].Vel[j] * sqrt(All.cf_a3inv));
	      break;
	    case SIO_PTYP:
	      for(i = 0; i < nwrite; i++)
		fp_int[i] = ID_list[i].Type;
	      break;
	    case SIO_PMAS:
	      for(i = 0; i < nwrite; i++)
		fp[i] = (MyOutputFloat) ID_list[i].Mass;
	      break;
#endif
	    case SIO_PID:
#ifndef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT	/* open new file in case of old format */
	      fclose(fd);
	      sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_ids", num, ThisTask);
	      if(!(fd = fopen(buf, "w")))
		{
		  printf("can't open file `%s`\n", buf);
		  endrun(1184);
		}
	      my_fwrite(&Ngroups, sizeof(int), 1, fd);
	      my_fwrite(&TotNgroups, sizeof(int), 1, fd);
	      my_fwrite(&Nids, sizeof(int), 1, fd);
	      my_fwrite(&TotNids, sizeof(long long), 1, fd);
	      my_fwrite(&NTask, sizeof(int), 1, fd);
	      my_fwrite(&Send_offset[ThisTask], sizeof(int), 1, fd);
#endif
	      for(i = 0; i < nwrite; i++)
		fp_id[i] = ID_list[i].ID;
	      break;
	    default:
	      break;

	    }			/* closing switch(blocknr) */

	  if(blksize > 0)
	    {
#if defined(HAVE_HDF5) && defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
	      if(All.SnapFormat != 3)
		{
#endif
		  my_fwrite(IOBuffer, blksize, sizeof(char), fd);
#ifdef SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT
		  SKIP;
#endif
#if defined(HAVE_HDF5) && defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
		}
	      else if(All.SnapFormat == 3)
		{
		  start[0] = 0;
		  start[1] = 0;

		  count[0] = nwrite;
		  count[1] = ndim;

		  H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

		  dims[0] = nwrite;
		  dims[1] = ndim;
		  hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

		  hdf5_status =
		    H5Dwrite(hdf5_dataset, hdf5_datatype,
			     hdf5_dataspace_memory, hdf5_dataspace_in_file, H5P_DEFAULT, IOBuffer);

		  H5Sclose(hdf5_dataspace_memory);
		}
#endif

	    }
	  myfree(IOBuffer);

	}			/* closing if block exist */

    }				/* closing for bnr loop */
#if defined(HAVE_HDF5) && defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
  if(All.SnapFormat == 3)
    {

      for(type = 0; type < 3; type++)
	H5Gclose(hdf5_grp[type]);
      H5Gclose(hdf5_headergrp);
      H5Fclose(hdf5_file);
    }
  else
    {
#endif
      fclose(fd);
#if defined(HAVE_HDF5) && defined(SUBFIND_WRITE_OUTPUTS_IN_SNAPSHOT_FORMAT)
    }
#endif
}



int subfind_compare_ID_list(const void *a, const void *b)
{
  if(((struct subfind_id_list *) a)->GrNr < ((struct subfind_id_list *) b)->GrNr) {return -1;}
  if(((struct subfind_id_list *) a)->GrNr > ((struct subfind_id_list *) b)->GrNr) {return +1;}
  if(((struct subfind_id_list *) a)->SubNr < ((struct subfind_id_list *) b)->SubNr) {return -1;}
  if(((struct subfind_id_list *) a)->SubNr > ((struct subfind_id_list *) b)->SubNr) {return +1;}
  if(((struct subfind_id_list *) a)->BindingEgy < ((struct subfind_id_list *) b)->BindingEgy) {return -1;}
  if(((struct subfind_id_list *) a)->BindingEgy > ((struct subfind_id_list *) b)->BindingEgy) {return +1;}
  return 0;
}

int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b)
{
  if(((struct subgroup_properties *) a)->GrNr < ((struct subgroup_properties *) b)->GrNr) {return -1;}
  if(((struct subgroup_properties *) a)->GrNr > ((struct subgroup_properties *) b)->GrNr) {return +1;}
  if(((struct subgroup_properties *) a)->SubNr < ((struct subgroup_properties *) b)->SubNr) {return -1;}
  if(((struct subgroup_properties *) a)->SubNr > ((struct subgroup_properties *) b)->SubNr) {return +1;}
  return 0;
}


int subfind_compare_P_GrNr_DM_Density(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr)) {return -1;}
  if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr)) {return +1;}
  if(((struct particle_data *) a)->u.DM_Density > (((struct particle_data *) b)->u.DM_Density)) {return -1;}
  if(((struct particle_data *) a)->u.DM_Density < (((struct particle_data *) b)->u.DM_Density)) {return +1;}
  return 0;
}

#endif

#if defined(SUBFIND)  // || defined(WRITE_KEY_FILES)

int subfind_compare_P_origindex(const void *a, const void *b)
{
  if(((struct particle_data *) a)->origindex < (((struct particle_data *) b)->origindex)) {return -1;}
  if(((struct particle_data *) a)->origindex > (((struct particle_data *) b)->origindex)) {return +1;}
  return 0;
}


#endif
