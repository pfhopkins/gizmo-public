#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../allvars.h"
#include "../proto.h"



/*
* This file was originally part of the GADGET3 code developed by
* Volker Springel. The code has been modified
* substantially (condensed, new feedback routines added,
* some optimizatins, and new variable/memory conventions added)
* by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
*/

void force_update_tree(void)
{
    PRINT_STATUS("Kick-subroutine will prepare for dynamic update of tree");
    int i, j; GlobFlag++; DomainNumChanged = 0; DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));
    /* note: the current list of active particles still refers to that synchronized at the previous time. */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        force_kick_node(i, P[i].dp);    /* kick the parent nodes with this momentum difference, also updated maximum velocity, softening and soundspeed, if needed */
        for(j = 0; j < 3; j++) {P[i].dp[j] = 0;}
    }
    force_finish_kick_nodes();
    myfree(DomainList);
    PRINT_STATUS(" ..Tree has been updated dynamically");
}


void force_kick_node(int i, MyDouble * dp)
{
  int j, no; MyFloat v, vmax;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MyFloat rt_source_lum_dp[3];
#endif
#ifdef DM_SCALARFIELD_SCREENING
  MyFloat dp_dm[3];
#endif

  for(j = 0; j < 3; j++)
    {
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        double lum[N_RT_FREQ_BINS];
        int active_check = rt_get_source_luminosity(i,-1,lum);
        if(active_check) {rt_source_lum_dp[j]=dp[j];} else {rt_source_lum_dp[j]=0;}
#endif
#ifdef DM_SCALARFIELD_SCREENING
      if(P[i].Type != 0)
	dp_dm[j] = dp[j];
      else
	dp_dm[j] = 0;
#endif
    }

  for(j = 0, vmax = 0; j < 3; j++)
    if((v = fabs(P[i].Vel[j])) > vmax)
      vmax = v;

  no = Father[i];

  while(no >= 0)
    {
      force_drift_node(no, All.Ti_Current);

      for(j = 0; j < 3; j++)
	{
	  Extnodes[no].dp[j] += dp[j];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Extnodes[no].rt_source_lum_dp[j] += rt_source_lum_dp[j];
#endif
#ifdef DM_SCALARFIELD_SCREENING
	  Extnodes[no].dp_dm[j] += dp_dm[j];
#endif
	}

      if(Extnodes[no].vmax < vmax)
	Extnodes[no].vmax = vmax;

      Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);
      Extnodes[no].Ti_lastkicked = All.Ti_Current;

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* top-level tree-node reached */
	{
	  if(Extnodes[no].Flag != GlobFlag)
	    {
	      Extnodes[no].Flag = GlobFlag;
	      DomainList[DomainNumChanged++] = no;
	    }
	  break;
	}

      no = Nodes[no].u.d.father;
    }
}






void force_finish_kick_nodes(void)
{
  int i, j, no, ta, totDomainNumChanged;
  int *domainList_all;
  int *counts, *counts_dp, *offset_list, *offset_dp, *offset_vmax;
  MyLongDouble *domainDp_loc, *domainDp_all;

#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MyLongDouble *domainDp_stellarlum_loc, *domainDp_stellarlum_all;
#endif
#ifdef DM_SCALARFIELD_SCREENING
  MyLongDouble *domainDp_dm_loc, *domainDp_dm_all;
#endif
  MyFloat *domainVmax_loc, *domainVmax_all;

  /* share the momentum-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  counts_dp = (int *) mymalloc("counts_dp", sizeof(int) * NTask);
  offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
  offset_dp = (int *) mymalloc("offset_dp", sizeof(int) * NTask);
  offset_vmax = (int *) mymalloc("offset_vmax", sizeof(int) * NTask);

  domainDp_loc = (MyLongDouble *) mymalloc("domainDp_loc", DomainNumChanged * 3 * sizeof(MyLongDouble));
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    domainDp_stellarlum_loc = (MyLongDouble *) mymalloc("domainDp_stellarlum_loc", DomainNumChanged * 3 * sizeof(MyLongDouble));
#endif
#ifdef DM_SCALARFIELD_SCREENING
  domainDp_dm_loc = (MyLongDouble *) mymalloc("domainDp_dm_loc", DomainNumChanged * 3 * sizeof(MyLongDouble));
#endif
  domainVmax_loc = (MyFloat *) mymalloc("domainVmax_loc", DomainNumChanged * sizeof(MyFloat));

  for(i = 0; i < DomainNumChanged; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  domainDp_loc[i * 3 + j] = Extnodes[DomainList[i]].dp[j];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        domainDp_stellarlum_loc[i * 3 + j] = Extnodes[DomainList[i]].rt_source_lum_dp[j];
#endif
#ifdef DM_SCALARFIELD_SCREENING
	  domainDp_dm_loc[i * 3 + j] = Extnodes[DomainList[i]].dp_dm[j];
#endif
	}
      domainVmax_loc[i] = Extnodes[DomainList[i]].vmax;
    }

  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_dp[0] = 0, offset_vmax[0] = 0; ta < NTask;
      ta++)
    {
      totDomainNumChanged += counts[ta];
      if(ta > 0)
	{
	  offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
	  offset_dp[ta] = offset_dp[ta - 1] + counts[ta - 1] * 3 * sizeof(MyLongDouble);
	  offset_vmax[ta] = offset_vmax[ta - 1] + counts[ta - 1] * sizeof(MyFloat);
	}
    }

  PRINT_STATUS(" ..exchanged kick momenta for %d top-level nodes out of %d", totDomainNumChanged, NTopleaves);
  domainDp_all = (MyLongDouble *) mymalloc("domainDp_all", totDomainNumChanged * 3 * sizeof(MyLongDouble));
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    domainDp_stellarlum_all = (MyLongDouble *) mymalloc("domainDp_stellarlum_all", totDomainNumChanged * 3 * sizeof(MyLongDouble));
#endif
#ifdef DM_SCALARFIELD_SCREENING
  domainDp_dm_all =
    (MyLongDouble *) mymalloc("domainDp_dm_all", totDomainNumChanged * 3 * sizeof(MyLongDouble));
#endif
  domainVmax_all = (MyFloat *) mymalloc("domainVmax_all", totDomainNumChanged * sizeof(MyFloat));

  domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

  MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
		 domainList_all, counts, offset_list, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    {
      counts_dp[ta] = counts[ta] * 3 * sizeof(MyLongDouble);
      counts[ta] *= sizeof(MyFloat);
    }


  MPI_Allgatherv(domainDp_loc, DomainNumChanged * 3 * sizeof(MyLongDouble), MPI_BYTE,
		 domainDp_all, counts_dp, offset_dp, MPI_BYTE, MPI_COMM_WORLD);

#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MPI_Allgatherv(domainDp_stellarlum_loc, DomainNumChanged * 3 * sizeof(MyLongDouble), MPI_BYTE,
                   domainDp_stellarlum_all, counts_dp, offset_dp, MPI_BYTE, MPI_COMM_WORLD);
#endif
#ifdef DM_SCALARFIELD_SCREENING
  MPI_Allgatherv(domainDp_dm_loc, DomainNumChanged * 3 * sizeof(MyLongDouble), MPI_BYTE,
		 domainDp_dm_all, counts_dp, offset_dp, MPI_BYTE, MPI_COMM_WORLD);
#endif

  MPI_Allgatherv(domainVmax_loc, DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
		 domainVmax_all, counts, offset_vmax, MPI_BYTE, MPI_COMM_WORLD);


  /* construct momentum kicks in top-level tree */
  for(i = 0; i < totDomainNumChanged; i++)
    {
      no = domainList_all[i];

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))	/* to avoid that the local one is kicked twice */
	no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  force_drift_node(no, All.Ti_Current);

	  for(j = 0; j < 3; j++)
	    {
	      Extnodes[no].dp[j] += domainDp_all[3 * i + j];
#ifdef RT_SEPARATELY_TRACK_LUMPOS
            Extnodes[no].rt_source_lum_dp[j] += domainDp_stellarlum_all[3 * i + j];
#endif
#ifdef DM_SCALARFIELD_SCREENING
	      Extnodes[no].dp_dm[j] += domainDp_dm_all[3 * i + j];
#endif
	    }

	  if(Extnodes[no].vmax < domainVmax_all[i])
	    Extnodes[no].vmax = domainVmax_all[i];

	  Nodes[no].u.d.bitflags |= (1 << BITFLAG_NODEHASBEENKICKED);
	  Extnodes[no].Ti_lastkicked = All.Ti_Current;

	  no = Nodes[no].u.d.father;
	}
    }

  myfree(domainList_all);
  myfree(domainVmax_all);
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    myfree(domainDp_stellarlum_all);
#endif
#ifdef DM_SCALARFIELD_SCREENING
  myfree(domainDp_dm_all);
#endif
  myfree(domainDp_all);
  myfree(domainVmax_loc);
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    myfree(domainDp_stellarlum_loc);
#endif
#ifdef DM_SCALARFIELD_SCREENING
  myfree(domainDp_dm_loc);
#endif
  myfree(domainDp_loc);
  myfree(offset_vmax);
  myfree(offset_dp);
  myfree(offset_list);
  myfree(counts_dp);
  myfree(counts);
}



void force_drift_node(int no, integertime time1)
{
  int j;
  integertime time0;
  double dt_drift, dt_drift_hmax, fac;

  if(time1 == Nodes[no].Ti_current)
    return;

  time0 = Extnodes[no].Ti_lastkicked;

  if(Nodes[no].u.d.bitflags & (1 << BITFLAG_NODEHASBEENKICKED))
    {
      if(Extnodes[no].Ti_lastkicked != Nodes[no].Ti_current)
	{
	  printf("Task=%d Extnodes[no].Ti_lastkicked=%lld  Nodes[no].Ti_current=%lld\n",ThisTask, (long long)Extnodes[no].Ti_lastkicked, (long long)Nodes[no].Ti_current);
	  terminate("inconsistency in drift node");
	}

      if(Nodes[no].u.d.mass)
	fac = 1 / Nodes[no].u.d.mass;
      else
	fac = 0;

#ifdef RT_SEPARATELY_TRACK_LUMPOS
        double fac_stellar_lum;
        double l_tot=0; for(j=0;j<N_RT_FREQ_BINS;j++) {l_tot += (Nodes[no].stellar_lum[j]);}
        if(l_tot>0) {fac_stellar_lum = 1 / l_tot;} else {fac_stellar_lum = 0;}
#endif

#ifdef DM_SCALARFIELD_SCREENING
      double fac_dm;

      if(Nodes[no].mass_dm)
	fac_dm = 1 / Nodes[no].mass_dm;
      else
	fac_dm = 0;
#endif

      for(j = 0; j < 3; j++)
	{
	  Extnodes[no].vs[j] += fac * FLT(Extnodes[no].dp[j]);
	  Extnodes[no].dp[j] = 0;
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Extnodes[no].rt_source_lum_vs[j] += fac_stellar_lum * FLT(Extnodes[no].rt_source_lum_dp[j]);
        Extnodes[no].rt_source_lum_dp[j] = 0;
#endif
#ifdef DM_SCALARFIELD_SCREENING
	  Extnodes[no].vs_dm[j] += fac_dm * FLT(Extnodes[no].dp_dm[j]);
	  Extnodes[no].dp_dm[j] = 0;
#endif

	}
      Nodes[no].u.d.bitflags &= (~(1 << BITFLAG_NODEHASBEENKICKED));
    }

  if(All.ComovingIntegrationOn)
    {
      dt_drift_hmax = get_drift_factor(Nodes[no].Ti_current, time1);
      dt_drift = dt_drift_hmax;
    }
  else
    {
      dt_drift_hmax = (time1 - Nodes[no].Ti_current) * All.Timebase_interval;
      dt_drift = dt_drift_hmax;
    }

    for(j = 0; j < 3; j++) {Nodes[no].u.d.s[j] += Extnodes[no].vs[j] * dt_drift;}
  Nodes[no].len += 2 * Extnodes[no].vmax * dt_drift;

#ifdef DM_SCALARFIELD_SCREENING
    for(j = 0; j < 3; j++) {Nodes[no].s_dm[j] += Extnodes[no].vs_dm[j] * dt_drift;}
#endif


#ifdef RT_SEPARATELY_TRACK_LUMPOS
    for(j = 0; j < 3; j++) {Nodes[no].rt_source_lum_s[j] += Extnodes[no].rt_source_lum_vs[j] * dt_drift;}
#endif
    
  Extnodes[no].hmax *= exp(Extnodes[no].divVmax * dt_drift_hmax / NUMDIMS);
  Nodes[no].Ti_current = time1;
}





/*! This function updates the hmax-values in tree nodes that hold SPH
 *  particles. These values are needed to find all neighbors in the
 *  hydro-force computation.  Since the Hsml-values are potentially changed
 *  in the SPH-density computation, force_update_hmax() should be carried
 *  out just before the hydrodynamical SPH forces are computed, i.e. after
 *  density().
 */
void force_update_hmax(void)
{
  int i, no, ta, totDomainNumChanged;
  int *domainList_all;
  int *counts, *offset_list, *offset_hmax;
  MyFloat *domainHmax_loc, *domainHmax_all;
  int OffsetSIZE = 2;
  double divVel;

  GlobFlag++;

  DomainNumChanged = 0;
  DomainList = (int *) mymalloc("DomainList", NTopleaves * sizeof(int));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
  {
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
    if(P[i].Mass > 0)
#else
    if(P[i].Type == 0 && P[i].Mass > 0)
#endif
      {
        no = Father[i];
        divVel = P[i].Particle_DivVel;
        
        while(no >= 0)
        {
            force_drift_node(no, All.Ti_Current);
            
            if(PPP[i].Hsml > Extnodes[no].hmax || divVel > Extnodes[no].divVmax)
            {
                if(PPP[i].Hsml > Extnodes[no].hmax) {Extnodes[no].hmax = PPP[i].Hsml;}
                if(divVel > Extnodes[no].divVmax) {Extnodes[no].divVmax = divVel;}
                
                if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node */
                {
                    if(Extnodes[no].Flag != GlobFlag)
                    {
                        Extnodes[no].Flag = GlobFlag;
                        DomainList[DomainNumChanged++] = no;
                    }
                    break;
                }
            }
            else
                break;
            
            no = Nodes[no].u.d.father;
        }
      }
  } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])

  /* share the hmax-data of the pseudo-particles accross CPUs */

  counts = (int *) mymalloc("counts", sizeof(int) * NTask);
  offset_list = (int *) mymalloc("offset_list", sizeof(int) * NTask);
  offset_hmax = (int *) mymalloc("offset_hmax", sizeof(int) * NTask);

  domainHmax_loc = (MyFloat *) mymalloc("domainHmax_loc", DomainNumChanged * OffsetSIZE * sizeof(MyFloat));

  for(i = 0; i < DomainNumChanged; i++)
    {
      domainHmax_loc[OffsetSIZE * i] = Extnodes[DomainList[i]].hmax;
      domainHmax_loc[OffsetSIZE * i + 1] = Extnodes[DomainList[i]].divVmax;
    }


  MPI_Allgather(&DomainNumChanged, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0, totDomainNumChanged = 0, offset_list[0] = 0, offset_hmax[0] = 0; ta < NTask; ta++)
    {
      totDomainNumChanged += counts[ta];
      if(ta > 0)
	{
	  offset_list[ta] = offset_list[ta - 1] + counts[ta - 1];
	  offset_hmax[ta] = offset_hmax[ta - 1] + counts[ta - 1] * OffsetSIZE * sizeof(MyFloat);
	}
    }

  PRINT_STATUS(" ..Hmax exchange: %d topleaves out of %d", totDomainNumChanged, NTopleaves);
  domainHmax_all = (MyFloat *) mymalloc("domainHmax_all", totDomainNumChanged * OffsetSIZE * sizeof(MyFloat));
  domainList_all = (int *) mymalloc("domainList_all", totDomainNumChanged * sizeof(int));

  MPI_Allgatherv(DomainList, DomainNumChanged, MPI_INT,
		 domainList_all, counts, offset_list, MPI_INT, MPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    counts[ta] *= OffsetSIZE * sizeof(MyFloat);

  MPI_Allgatherv(domainHmax_loc, OffsetSIZE * DomainNumChanged * sizeof(MyFloat), MPI_BYTE,
		 domainHmax_all, counts, offset_hmax, MPI_BYTE, MPI_COMM_WORLD);


  for(i = 0; i < totDomainNumChanged; i++)
    {
      no = domainList_all[i];

      if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))	/* to avoid that the hmax is updated twice */
	no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  force_drift_node(no, All.Ti_Current);


	  if(domainHmax_all[OffsetSIZE * i] > Extnodes[no].hmax || domainHmax_all[OffsetSIZE * i + 1] > Extnodes[no].divVmax)
	    {
	      if(domainHmax_all[OffsetSIZE * i] > Extnodes[no].hmax)
		Extnodes[no].hmax = domainHmax_all[OffsetSIZE * i];

	      if(domainHmax_all[OffsetSIZE * i + 1] > Extnodes[no].divVmax)
		Extnodes[no].divVmax = domainHmax_all[OffsetSIZE * i + 1];
	    }
	  else
	    break;

	  no = Nodes[no].u.d.father;
	}
    }


  myfree(domainList_all);
  myfree(domainHmax_all);
  myfree(domainHmax_loc);
  myfree(offset_hmax);
  myfree(offset_list);
  myfree(counts);
  myfree(DomainList);

  CPU_Step[CPU_TREEHMAXUPDATE] += measure_time();
}
