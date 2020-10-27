#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

#include "subfind.h"
#include "../fof.h"

/* this file processes the local groups in serial mode */


#ifndef MAX_NGB_CHECK
#define MAX_NGB_CHECK 2 /* max numbers of neighbours for saddle-point detection (default = 2) */
#endif

static int *Head, *Next, *Tail, *Len;
static struct cand_dat
{
  int head;
  int len;
  int nsub;
  int rank, subnr, parent;
  int bound_length;
}
 *candidates;




/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here!!! -- */
int subfind_process_group_serial(int gr, int Offs)
{
  int i, j, k, p, len, subnr, totlen, ss, ngbs, ndiff, N, head = 0, head_attach, count_cand, len_non_gas;
  int listofdifferent[2], count, prev;
  int ngb_index, part_index, nsubs, rank;
  double SubMass, SubPos[3], SubVel[3], SubCM[3], SubVelDisp, SubVmax, SubVmaxRad, SubSpin[3], SubHalfMass, SubMassTab[6];
  MyIDType SubMostBoundID;
  static struct unbind_data *ud;

  while(P[Offs].GrNr != Group[gr].GrNr)
    {
      Offs++;
      if(Offs >= NumPart)
	{
	  printf("don't find a particle for groupnr=%d\n", Group[gr].GrNr);
	  endrun(312);
	}
    }

  N = Group[gr].Len;
  GrNr = Group[gr].GrNr;

  for(i = 0; i < N; i++)
    {
      if(P[Offs + i].GrNr != Group[gr].GrNr)
	{
	  printf
	    ("task=%d, gr=%d: don't have the number of particles for GrNr=%d group-len=%d found=%d before=%d\n",
	     ThisTask, gr, Group[gr].GrNr, N, P[Offs + i].GrNr, P[Offs - 1].GrNr);
	  endrun(312);
	}
    }


  candidates = (struct cand_dat *)mymalloc("candidates", N * sizeof(struct cand_dat));

  Head = (int *)mymalloc("Head", N * sizeof(int));
  Next = (int *)mymalloc("Next", N * sizeof(int));
  Tail = (int *)mymalloc("Tail", N * sizeof(int));
  Len = (int *)mymalloc("Len", N * sizeof(int));
  ud = (struct unbind_data *) mymalloc("ud", N * sizeof(struct unbind_data));

  Head -= Offs;
  Next -= Offs;
  Tail -= Offs;
  Len -= Offs;

  for(i = 0; i < N; i++)
    {
      ud[i].index = Offs + i;
    }

  subfind_loctree_findExtent(N, ud);

  subfind_loctree_treebuild(N, ud);	/* build tree for all particles of this group */

  for(i = Offs; i < Offs + N; i++)
    Head[i] = Next[i] = Tail[i] = -1;

  /* note: particles are already ordered in the order of decreasing density */

  for(i = 0, count_cand = 0; i < N; i++)
    {
      part_index = Offs + i;

      subfind_locngb_treefind(P[part_index].Pos, All.DesLinkNgb, P[part_index].DM_Hsml);

      /* note: returned neighbours are already sorted by distance */

      for(k = 0, ndiff = 0, ngbs = 0; k < All.DesLinkNgb && ngbs < MAX_NGB_CHECK && ndiff < 2; k++)
	{
	  ngb_index = R2list[k].index;

	  if(ngb_index != part_index)	/* to exclude the particle itself */
	    {
	      /* we only look at neighbours that are denser */
	      if(P[ngb_index].u.DM_Density > P[part_index].u.DM_Density)
		{
		  ngbs++;

		  if(Head[ngb_index] >= 0)	/* neighbor is attached to a group */
		    {
		      if(ndiff == 1)
			if(listofdifferent[0] == Head[ngb_index])
			  continue;

		      /* a new group has been found */
		      listofdifferent[ndiff++] = Head[ngb_index];
		    }
		  else
		    {
		      printf("this may not occur.\n");
		      printf
			("ThisTask=%d gr=%d k=%d i=%d part_index=%d ngb_index = %d  head[ngb_index]=%d P[part_index].DM_Density=%g %g GrNrs= %d %d \n",
			 ThisTask, gr, k, i, part_index, ngb_index, Head[ngb_index],
			 P[part_index].u.DM_Density, P[ngb_index].u.DM_Density, P[part_index].GrNr,
			 P[ngb_index].GrNr);
		      endrun(2);
		    }
		}
	    }
	}

      switch (ndiff)		/* treat the different possible cases */
	{
	case 0:		/* this appears to be a lonely maximum -> new group */
	  head = part_index;
	  Head[part_index] = Tail[part_index] = part_index;
	  Len[part_index] = 1;
	  Next[part_index] = -1;
	  break;

	case 1:		/* the particle is attached to exactly one group */
	  head = listofdifferent[0];
	  Head[part_index] = head;
	  Next[Tail[head]] = part_index;
	  Tail[head] = part_index;
	  Len[head]++;
	  Next[part_index] = -1;
	  break;

	case 2:		/* the particle merges two groups together */

	  head = listofdifferent[0];
	  head_attach = listofdifferent[1];

	  if(Len[head_attach] > Len[head])	/* other group is longer, swap them */
	    {
	      head = listofdifferent[1];
	      head_attach = listofdifferent[0];
	    }

	  /* only in case the attached group is long enough we bother to register is 
	     as a subhalo candidate */

	  if(Len[head_attach] >= All.DesLinkNgb)
	    {
	      candidates[count_cand].len = Len[head_attach];
	      candidates[count_cand].head = Head[head_attach];
	      count_cand++;
	    }

	  /* now join the two groups */
	  Next[Tail[head]] = head_attach;
	  Tail[head] = Tail[head_attach];
	  Len[head] += Len[head_attach];

	  ss = head_attach;
	  do
	    {
	      Head[ss] = head;
	    }
	  while((ss = Next[ss]) >= 0);

	  /* finally, attach the particle */
	  Head[part_index] = head;
	  Next[Tail[head]] = part_index;
	  Tail[head] = part_index;
	  Len[head]++;
	  Next[part_index] = -1;
	  break;

	default:
	  printf("can't be! (a)\n");
	  endrun(1);
	  break;
	}
    }

  /* add the full thing as a subhalo candidate */
  for(i = 0, prev = -1; i < N; i++)
    {
      if(Head[Offs + i] == Offs + i)
	if(Next[Tail[Offs + i]] == -1)
	  {
	    if(prev < 0)
	      head = Offs + i;
	    if(prev >= 0)
	      Next[prev] = Offs + i;

	    prev = Tail[Offs + i];
	  }
    }

  candidates[count_cand].len = N;
  candidates[count_cand].head = head;
  count_cand++;

  /* go through them once and assign the rank */
  for(i = 0, p = head, rank = 0; i < N; i++)
    {
      Len[p] = rank++;
      p = Next[p];
    }

  /* for each candidate, we now pull out the rank of its head */
  for(k = 0; k < count_cand; k++)
    candidates[k].rank = Len[candidates[k].head];

  for(i = Offs; i < Offs + N; i++)
    Tail[i] = -1;

  for(k = 0, nsubs = 0; k < count_cand; k++)
    {
      for(i = 0, p = candidates[k].head, len = 0; i < candidates[k].len; i++, p = Next[p])
	if(Tail[p] < 0)
	  ud[len++].index = p;

      if(len >= All.DesLinkNgb)
	len = subfind_unbind(ud, len, &len_non_gas);

#ifdef SUBFIND_REMOVE_GAS_STRUCTURES
      if(len_non_gas >= All.DesLinkNgb)
#else
      if(len >= All.DesLinkNgb)
#endif
	{
	  /* ok, we found a substructure */

	  for(i = 0; i < len; i++)
	    Tail[ud[i].index] = nsubs;	/* we use this to flag the substructures */

	  candidates[k].nsub = nsubs;
	  candidates[k].bound_length = len;
	  nsubs++;
	}
      else
	{
	  candidates[k].nsub = -1;
	  candidates[k].bound_length = 0;
	}
    }

  PRINT_STATUS("\nGroupLen=%d  (gr=%d)", N, gr);
  PRINT_STATUS("Number of substructures: %d", nsubs);

  Group[gr].Nsubs = nsubs;
  Group[gr].Pos[0] = Group[gr].CM[0];
  Group[gr].Pos[1] = Group[gr].CM[1];
  Group[gr].Pos[2] = Group[gr].CM[2];

  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_boundlength);
  /* now we determine the parent subhalo for each candidate */
  for(k = 0; k < count_cand; k++)
    {
      candidates[k].subnr = k;
      candidates[k].parent = 0;
    }

  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_rank);

  for(k = 0; k < count_cand; k++)
    {
      for(j = k + 1; j < count_cand; j++)
	{
	  if(candidates[j].rank > candidates[k].rank + candidates[k].len)
	    break;

	  if(candidates[k].rank + candidates[k].len >= candidates[j].rank + candidates[j].len)
	    {
	      if(candidates[k].bound_length >= All.DesLinkNgb)
		candidates[j].parent = candidates[k].subnr;
	    }
	  else
	    {
	      printf("k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n",
		     k, count_cand, (int) candidates[k].rank, candidates[k].len,
		     (int) candidates[k].bound_length, candidates[j].rank,
		     (int) candidates[j].len, candidates[j].bound_length);
	      endrun(121235513);
	    }
	}
    }

  qsort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_subnr);
  /* now determine the properties */

  for(k = 0, subnr = 0, totlen = 0; k < nsubs; k++)
    {
      len = candidates[k].bound_length;

      PRINT_STATUS("subnr=%d  SubLen=%d", subnr, len);

      totlen += len;

      for(i = 0, p = candidates[k].head, count = 0; i < candidates[k].len; i++)
	{
	  if(Tail[p] == candidates[k].nsub)
	    ud[count++].index = p;

	  p = Next[p];
	}

      if(count != len)
	endrun(12);


      subfind_determine_sub_halo_properties(ud, len, &SubMass,
					    &SubPos[0], &SubVel[0], &SubCM[0], &SubVelDisp, &SubVmax,
					    &SubVmaxRad, &SubSpin[0], &SubMostBoundID, &SubHalfMass,
					    &SubMassTab[0]);

      if(Nsubgroups >= MaxNsubgroups)
	endrun(899);

      if(subnr == 0)
	{
	  for(j = 0; j < 3; j++)
	    Group[gr].Pos[j] = SubPos[j];
	}

      SubGroup[Nsubgroups].Len = len;
      if(subnr == 0)
	SubGroup[Nsubgroups].Offset = Group[gr].Offset;
      else
	SubGroup[Nsubgroups].Offset = SubGroup[Nsubgroups - 1].Offset + SubGroup[Nsubgroups - 1].Len;
      SubGroup[Nsubgroups].GrNr = GrNr - 1;
      SubGroup[Nsubgroups].SubNr = subnr;
      SubGroup[Nsubgroups].SubParent = candidates[k].parent;
      SubGroup[Nsubgroups].Mass = SubMass;
      SubGroup[Nsubgroups].SubMostBoundID = SubMostBoundID;
      SubGroup[Nsubgroups].SubVelDisp = SubVelDisp;
      SubGroup[Nsubgroups].SubVmax = SubVmax;
      SubGroup[Nsubgroups].SubVmaxRad = SubVmaxRad;
      SubGroup[Nsubgroups].SubHalfMass = SubHalfMass;

      for(j = 0; j < 3; j++)
	{
	  SubGroup[Nsubgroups].Pos[j] = SubPos[j];
	  SubGroup[Nsubgroups].CM[j] = SubCM[j];
	  SubGroup[Nsubgroups].Vel[j] = SubVel[j];
	  SubGroup[Nsubgroups].Spin[j] = SubSpin[j];
	}

      for(j = 0; j < 6; j++) SubGroup[Nsubgroups].MassTab[j] = SubMassTab[j];

      Nsubgroups++;

      /* Let's now assign the subgroup number */

      for(i = 0; i < len; i++)
	P[ud[i].index].SubNr = subnr;

      subnr++;
    }

  PRINT_STATUS("Fuzz=%d", N - totlen);

  myfree(ud);
  myfree(Len + Offs);
  myfree(Tail + Offs);
  myfree(Next + Offs);
  myfree(Head + Offs);

  myfree(candidates);

  return Offs;
}




int subfind_unbind(struct unbind_data *ud, int len, int *len_non_gas)
{
  double *bnd_energy, energy_limit, weakly_bound_limit = 0;
  int i, j, p, minindex, unbound, phaseflag;
  double s[3], dx[3], v[3], dv[3], pos[3];
  double vel_to_phys, H_of_a, atime, pot, minpot = 0;
  double boxsize;
  double TotMass;

  boxsize = All.BoxSize;

  if(All.ComovingIntegrationOn)
    {
      vel_to_phys = 1.0 / All.Time;
      H_of_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    {
      vel_to_phys = atime = 1;
      H_of_a = 0;
    }

  bnd_energy = (double *) mymalloc("bnd_energy", len * sizeof(double));

  phaseflag = 0;		/* this means we will recompute the potential for all particles */

  do
    {
      subfind_loctree_treebuild(len, ud);

      /* let's compute the potential  */

      if(phaseflag == 0)	/* redo it for all the particles */
	{
	  for(i = 0, minindex = -1, minpot = 1.0e30; i < len; i++)
	    {
	      p = ud[i].index;

	      pot = subfind_loctree_treeevaluate_potential(p);
	      /* note: add self-energy */
        double h_grav = All.ForceSoftening[P[p].Type];
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
        h_grav = P[p].AGS_Hsml;
#elif defined(ADAPTIVE_GRAVSOFT_FORGAS)
        if(P[p].Type == 0) h_grav = PPP[p].Hsml;
#endif
      P[p].u.DM_Potential = pot - P[p].Mass / h_grav * kernel_gravity(0,1,1,-1); // subtracts self-contribution
	      P[p].u.DM_Potential *= All.G / atime;

	      if(All.TotN_gas > 0 && (FOF_PRIMARY_LINK_TYPES & 1) == 0 && (FOF_SECONDARY_LINK_TYPES & 1) == 0 && All.OmegaBaryon > 0) {P[p].u.DM_Potential *= All.OmegaMatter / (All.OmegaMatter - All.OmegaBaryon);}

	      if(P[p].u.DM_Potential < minpot || minindex == -1)
		{
		  minpot = P[p].u.DM_Potential;
		  minindex = p;
		}
	    }

	  for(j = 0; j < 3; j++) {pos[j] = P[minindex].Pos[j];}	/* position of minimum potential */
	}
      else
	{
	  /* we only repeat for those close to the unbinding threshold */
	  for(i = 0; i < len; i++)
	    {
	      p = ud[i].index;

	      if(P[p].v.DM_BindingEnergy >= weakly_bound_limit)
		{
		  pot = subfind_loctree_treeevaluate_potential(p);
		  /* note: add self-energy */
        double h_grav = All.ForceSoftening[P[p].Type];
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
        h_grav = P[p].AGS_Hsml;
#elif defined(ADAPTIVE_GRAVSOFT_FORGAS)
        if(P[p].Type == 0) h_grav = PPP[p].Hsml;
#endif
        P[p].u.DM_Potential = pot - P[p].Mass / h_grav * kernel_gravity(0,1,1,-1); // subtract self-contribution
		  P[p].u.DM_Potential = pot + P[p].Mass / h_grav;
		  P[p].u.DM_Potential *= All.G / atime;

		  if(All.TotN_gas > 0 && (FOF_PRIMARY_LINK_TYPES & 1) == 0 && (FOF_SECONDARY_LINK_TYPES & 1) == 0 && All.OmegaBaryon > 0) {P[p].u.DM_Potential *= All.OmegaMatter / (All.OmegaMatter - All.OmegaBaryon);}
		}
	    }
	}

      /* let's get bulk velocity and the center-of-mass */

      v[0] = v[1] = v[2] = 0;
      s[0] = s[1] = s[2] = 0;

      for(i = 0, TotMass = 0; i < len; i++)
	{
	  p = ud[i].index;

        double dp[3]; for(j=0;j<3;j++) {dp[j]=P[p].Pos[j]-pos[j];}
        NEAREST_XYZ(dp[0],dp[1],dp[2],-1);
      for(j = 0; j < 3; j++)
        {
          s[j] += P[p].Mass * dp[j];
	      v[j] += P[p].Mass * P[p].Vel[j];
	    }
	  TotMass += P[p].Mass;
	}

      for(j = 0; j < 3; j++)
	{
	  v[j] /= TotMass;
	  s[j] /= TotMass;	/* center-of-mass */

	  s[j] += pos[j];

#ifdef BOX_PERIODIC
	  while(s[j] < 0) {s[j] += boxsize;}
	  while(s[j] >= boxsize) {s[j] -= boxsize;}
#endif
	}

      for(i = 0; i < len; i++)
	{
	  p = ud[i].index;

        
          double dp[3]; for(j=0;j<3;j++) {dp[j]=P[p].Pos[j]-s[j];}
          NEAREST_XYZ(dp[0],dp[1],dp[2],-1);

        for(j = 0; j < 3; j++)
          {
            dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
              dx[j] = dp[j] * atime;
            dv[j] += H_of_a * dx[j];
          }


	  P[p].v.DM_BindingEnergy = P[p].u.DM_Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);

#ifdef FOF_DENSITY_SPLIT_TYPES
	  if(P[p].Type == 0) {P[p].v.DM_BindingEnergy += P[p].w.int_energy;}
#endif
	  bnd_energy[i] = P[p].v.DM_BindingEnergy;
	}

      qsort(bnd_energy, len, sizeof(double), subfind_compare_binding_energy);	/* largest comes first! */
      energy_limit = bnd_energy[(int) (0.25 * len)];

      for(i = 0, unbound = 0; i < len - 1; i++)
	{
	  if(bnd_energy[i] > 0)
	    unbound++;
	  else
	    unbound--;

	  if(unbound <= 0)
	    break;
	}
      weakly_bound_limit = bnd_energy[i];

      /* now omit unbound particles,  but at most 1/4 of the original size */

      for(i = 0, unbound = 0, *len_non_gas = 0; i < len; i++)
	{
	  p = ud[i].index;
	  if(P[p].v.DM_BindingEnergy > 0 && P[p].v.DM_BindingEnergy > energy_limit)
	    {
	      unbound++;
	      ud[i] = ud[len - 1];
	      i--;
	      len--;
	    }
	  else if(P[p].Type != 0)
	    (*len_non_gas)++;
	}

      if(len < All.DesLinkNgb)
	break;

      if(phaseflag == 0)
	{
	  if(unbound > 0)
	    phaseflag = 1;
	}
      else
	{
	  if(unbound == 0)
	    {
	      phaseflag = 0;	/* this will make us repeat everything once more for all particles */
	      unbound = 1;
	    }
	}
    }
  while(unbound > 0);

  myfree(bnd_energy);

  return (len);
}



int subfind_compare_grp_particles(const void *a, const void *b)
{
  if(((struct particle_data *) a)->GrNr < ((struct particle_data *) b)->GrNr) {return -1;}
  if(((struct particle_data *) a)->GrNr > ((struct particle_data *) b)->GrNr) {return +1;}
  if(((struct particle_data *) a)->SubNr < ((struct particle_data *) b)->SubNr) {return -1;}
  if(((struct particle_data *) a)->SubNr > ((struct particle_data *) b)->SubNr) {return +1;}
  if(((struct particle_data *) a)->v.DM_BindingEnergy < ((struct particle_data *) b)->v.DM_BindingEnergy) {return -1;}
  if(((struct particle_data *) a)->v.DM_BindingEnergy > ((struct particle_data *) b)->v.DM_BindingEnergy) {return +1;}
  return 0;
}

void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					   double *pos, double *vel, double *cm, double *veldisp,
					   double *vmax, double *vmaxrad, double *spin,
					   MyIDType * mostboundid, double *halfmassrad, double *mass_tab)
{
  int i, j, p, num_use, i_use;
  double s[3], v[3], max, vel_to_phys, H_of_a, atime, minpot;
  double lx, ly, lz, dv[3], dx[3], disp, rr_tmp, disp_tmp;
  double boxsize, ddxx;
  sort_r2list *rr_list = 0;
  int minindex;
  double mass, maxrad;
  int nstar = 0, ndm = 0, ngas = 0, nbh = 0;

  boxsize = All.BoxSize;

  if(All.ComovingIntegrationOn)
    {
      vel_to_phys = 1.0 / All.Time;
      H_of_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    {
      vel_to_phys = atime = 1;
      H_of_a = 0;
    }

  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      if(P[p].u.DM_Potential < minpot || minindex == -1)
	{
	  minpot = P[p].u.DM_Potential;
	  minindex = p;
	}
      switch (P[p].Type)
	{
	case 0:
	  ngas++;
	  break;
	case 1:
	case 2:
	case 3:
	  ndm++;
	  break;
	case 4:
	  nstar++;
	  break;
	case 5:
	  nbh++;
	  break;
	}
    }

  if(minindex == -1)
    endrun(875412);

  for(j = 0; j < 3; j++)
    pos[j] = P[minindex].Pos[j];


  /* pos[] now holds the position of minimum potential */
  /* we take it that as the center */


  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
#ifdef SUBFIND_REMOVE_GAS_STRUCTURES
      if(P[p].Type > 0)
#endif
	if(P[p].v.DM_BindingEnergy < minpot || minindex == -1)
	  {
	    minpot = P[p].v.DM_BindingEnergy;
	    minindex = p;
	  }
    }

  if(minindex == -1)
    endrun(875413);

  *mostboundid = P[minindex].ID;


  /* let's get bulk velocity and the center-of-mass */
  /* here we still can take all particles */

  for(j = 0; j < 3; j++)
    s[j] = v[j] = 0;

  for(j = 0; j < 6; j++)
    mass_tab[j] = 0;

  for(i = 0, mass = 0; i < num; i++)
    {
      p = d[i].index;
          double dp[3]; for(j=0;j<3;j++) {dp[j]=P[p].Pos[j]-pos[j];}
          NEAREST_XYZ(dp[0],dp[1],dp[2],-1);
        for(j = 0; j < 3; j++)
          {
            s[j] += P[p].Mass * dp[j];
	  v[j] += P[p].Mass * P[p].Vel[j];
	}
      mass += P[p].Mass;

      mass_tab[P[p].Type] += P[p].Mass;
    }

  *totmass = mass;

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass;		/* center of mass */
      v[j] /= mass;

      vel[j] = vel_to_phys * v[j];
    }

  for(j = 0; j < 3; j++)
    {
      s[j] += pos[j];

#ifdef BOX_PERIODIC
      while(s[j] < 0) {s[j] += boxsize;}
      while(s[j] >= boxsize) {s[j] -= boxsize;}
#endif

      cm[j] = s[j];
    }


  disp = lx = ly = lz = 0;

  /* Here we have to perform only on the dm particles for consistency */
  mass = 0;
  num_use = num;
#ifdef FOF_DENSITY_SPLIT_TYPES
  num_use = ndm;
#endif

  if(num_use > 0)
    rr_list = (sort_r2list *)mymalloc("rr_list", sizeof(sort_r2list) * num_use);

  for(i = 0, i_use = 0; i < num; i++)
    {
      p = d[i].index;

          double dp_s[3], dp_p[3]; for(j=0;j<3;j++) {dp_s[j]=P[p].Pos[j]-s[j]; dp_p[j]=P[p].Pos[j]-pos[j];}
          NEAREST_XYZ(dp_s[0],dp_s[1],dp_s[2],-1); NEAREST_XYZ(dp_p[0],dp_p[1],dp_p[2],-1);

      for(j = 0, rr_tmp = 0, disp_tmp = 0; j < 3; j++)
	{
	  dx[j] = atime * dp_s[j];
	  dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
	  dv[j] += H_of_a * dx[j];

	  disp_tmp += P[p].Mass * dv[j] * dv[j];
	  /* for rotation curve computation, take minimum of potential as center */
	  ddxx = atime * dp_p[j];
	  rr_tmp += ddxx * ddxx;
	}

      rr_tmp = sqrt(rr_tmp);
#ifdef FOF_DENSITY_SPLIT_TYPES
      if(P[p].Type >= 1 && P[p].Type <= 3)  /*-- only for dm part --*/
#endif
      {
	  rr_list[i_use].mass = P[p].Mass;
	  rr_list[i_use].r = rr_tmp;
	  disp += disp_tmp;
	  mass += P[p].Mass;
	  lx += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
	  ly += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
	  lz += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

	  i_use++;
	}
    }

  if(i_use != num_use)
    endrun(564321);

  if(num_use > 0)
    {
      *veldisp = sqrt(disp / (3 * mass));	/* convert to 1d velocity dispersion */

      spin[0] = lx / mass;
      spin[1] = ly / mass;
      spin[2] = lz / mass;

      qsort(rr_list, num_use, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

/*--- Here we still have to fix for individual masses, 
        maybe we even want the total mass for this ? 
        Note however that even within multi mass DM simulation
        the DM mass within a clean halo should be the same ... ------*/
      *halfmassrad = rr_list[num_use / 2].r;

      /* compute cumulative mass */
      for(i = 1; i < num_use; i++)
	rr_list[i].mass = rr_list[i - 1].mass + rr_list[i].mass;

/*--- Note that here we might want to correct for the baryon fraction ? ---*/
      for(i = num_use - 1, max = 0, maxrad = 0; i > 5; i--)
	if(rr_list[i].mass / rr_list[i].r > max)
	  {
	    max = rr_list[i].mass / rr_list[i].r;
	    maxrad = rr_list[i].r;
	  }

      *vmax = sqrt(All.G * max);
      *vmaxrad = maxrad;

      myfree(rr_list);
    }
  else
    *veldisp = *halfmassrad = *vmax = *vmaxrad = spin[0] = spin[1] = spin[2] = 0;

}


void subfind_col_determine_sub_halo_properties(struct unbind_data *d, int num, double *totmass,
					       double *pos, double *vel, double *cm, double *veldisp,
					       double *vmax, double *vmaxrad, double *spin,
					       MyIDType * mostboundid, double *halfmassrad, double *mass_tab)
{
  int i, j, part_index, *npart, numtot, nhalf, offset, num_use, i_use;
  MyIDType mbid;
  double s[3], sloc[3], v[3], vloc[3], max, vel_to_phys, H_of_a, atime;
  double lx, ly, lz, dv[3], dx[3], disp;
  double loclx, locly, loclz, locdisp;
  double boxsize, ddxx;
  sort_r2list *loc_rr_list;
  int minindex, mincpu;
  double mass, mass_tab_loc[6], maxrad, massloc, *masslist;
  MyFloat minpot, *potlist;
  int nstarloc = 0, ndmloc = 0, ngasloc = 0, nbhloc = 0;
  int nstar = 0, ndm = 0, ngas = 0, nbh = 0;
  double rr_tmp, disp_tmp;
  boxsize = All.BoxSize;

  if(All.ComovingIntegrationOn)
    {
      vel_to_phys = 1.0 / All.Time;
      H_of_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    {
      H_of_a = 0;
      atime = vel_to_phys = 1;
    }

  potlist = (MyFloat *) mymalloc("potlist", NTask * sizeof(MyFloat));

  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      if(P[d[i].index].u.DM_Potential < minpot || minindex == -1)
	{
	  minpot = P[d[i].index].u.DM_Potential;
	  minindex = d[i].index;
	}
      switch (P[d[i].index].Type)
	{
	case 0:
	  ngasloc++;
	  break;
	case 1:
	case 2:
	case 3:
	  ndmloc++;
	  break;
	case 4:
	  nstarloc++;
	  break;
	case 5:
	  nbhloc++;
	  break;
	}
    }

  MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE, MPI_COMM_WORLD);
  MPI_Allreduce(&ngasloc, &ngas, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ndmloc, &ndm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nstarloc, &nstar, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nbhloc, &nbh, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0, mincpu = -1, minpot = 1.0e30; i < NTask; i++)
    if(potlist[i] < minpot)
      {
	mincpu = i;
	minpot = potlist[i];
      }

  if(mincpu < 0)
    {
      printf("ta=%d num=%d\n", ThisTask, num);
      endrun(121);
    }

  if(ThisTask == mincpu)
    {
      for(j = 0; j < 3; j++)
	s[j] = P[minindex].Pos[j];
    }

  MPI_Bcast(&s[0], 3, MPI_DOUBLE, mincpu, MPI_COMM_WORLD);

  /* s[] now holds the position of minimum potential */
  /* we take that as the center */
  for(j = 0; j < 3; j++)
    pos[j] = s[j];


  /* the ID of the most bound particle, we take the minimum binding energy */
  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
#ifdef SUBFIND_REMOVE_GAS_STRUCTURES
      if(P[d[i].index].Type > 0)
#endif
	if(P[d[i].index].v.DM_BindingEnergy < minpot || minindex == -1)
	  {
	    minpot = P[d[i].index].v.DM_BindingEnergy;
	    minindex = d[i].index;
	  }
    }

  MPI_Allgather(&minpot, sizeof(MyFloat), MPI_BYTE, potlist, sizeof(MyFloat), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0, minpot = 1.0e30; i < NTask; i++)
    if(potlist[i] < minpot)
      {
	mincpu = i;
	minpot = potlist[i];
      }

  if(ThisTask == mincpu)
    {
      if(minindex == -1)	/* This is to cover the quasi impossible case that one cpu have had only gas */
	endrun(875417);		/* particles and the potential in the simulation was larger than 1e30 !      */
      mbid = P[minindex].ID;
    }

  MPI_Bcast(&mbid, sizeof(mbid), MPI_BYTE, mincpu, MPI_COMM_WORLD);

  myfree(potlist);

  *mostboundid = mbid;

  /* let's get bulk velocity and the center-of-mass */

  for(j = 0; j < 3; j++)
    sloc[j] = vloc[j] = 0;

  for(j = 0; j < 6; j++)
    mass_tab_loc[j] = 0;

  for(i = 0, massloc = 0; i < num; i++)
    {
      part_index = d[i].index;

        double dp[3]; for(j=0;j<3;j++) {dp[j]=P[part_index].Pos[j]-pos[j];}
        NEAREST_XYZ(dp[0],dp[1],dp[2],-1);
      for(j = 0; j < 3; j++)
        {
          sloc[j] += P[part_index].Mass * dp[j];
	  vloc[j] += P[part_index].Mass * P[part_index].Vel[j];
	}
      massloc += P[part_index].Mass;

      mass_tab_loc[P[part_index].Type] += P[part_index].Mass;
    }

  MPI_Allreduce(sloc, s, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vloc, v, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(mass_tab_loc, mass_tab, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *totmass = mass;

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass;		/* center of mass */
      v[j] /= mass;

      vel[j] = vel_to_phys * v[j];

      s[j] += pos[j];

#ifdef BOX_PERIODIC
      while(s[j] < 0) {s[j] += boxsize;}
      while(s[j] >= boxsize) {s[j] -= boxsize;}
#endif

      cm[j] = s[j];
    }


  locdisp = loclx = locly = loclz = 0;


  /* Here we have to perform only on the dm particles for consistency */
  num_use = num;
#ifdef FOF_DENSITY_SPLIT_TYPES
  num_use = ndmloc;
#endif

  if(num_use > 0)
    loc_rr_list = (sort_r2list *)mymalloc("loc_rr_list", sizeof(sort_r2list) * num_use);

  for(i = 0, massloc = 0, i_use = 0; i < num; i++)
    {
      part_index = d[i].index;

        double dp_s[3], dp_p[3]; for(j=0;j<3;j++) {dp_s[j]=P[part_index].Pos[j]-s[j]; dp_p[j]=P[part_index].Pos[j]-pos[j];}
        NEAREST_XYZ(dp_s[0],dp_s[1],dp_s[2],-1); NEAREST_XYZ(dp_p[0],dp_p[1],dp_p[2],-1);

      for(j = 0, rr_tmp = 0, disp_tmp = 0; j < 3; j++)
	{
	  dx[j] = atime * dp_s[j];
	  dv[j] = vel_to_phys * (P[part_index].Vel[j] - v[j]);
	  dv[j] += H_of_a * dx[j];

	  disp_tmp += P[part_index].Mass * dv[j] * dv[j];
	  /* for rotation curve computation, take minimum of potential as center */
	  ddxx = atime * dp_p[j];
	  rr_tmp += ddxx * ddxx;
	}

      rr_tmp = sqrt(rr_tmp);

#ifdef FOF_DENSITY_SPLIT_TYPES
      if(P[part_index].Type >= 1 && P[part_index].Type <= 3)  /*-- only for dm part --*/
#endif
      {
	  locdisp += disp_tmp;
	  massloc += P[part_index].Mass;

	  loclx += P[part_index].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
	  locly += P[part_index].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
	  loclz += P[part_index].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

	  loc_rr_list[i_use].mass = P[part_index].Mass;
	  loc_rr_list[i_use].r = rr_tmp;

	  i_use++;
}
    }

  if(i_use != num_use)
    endrun(664321);

  MPI_Allreduce(&loclx, &lx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&locly, &ly, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&loclz, &lz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&locdisp, &disp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&massloc, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  npart = (int *) mymalloc("npart", NTask * sizeof(int));

  MPI_Allgather(&num_use, 1, MPI_INT, npart, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, numtot = 0; i < NTask; i++)
    numtot += npart[i];

  if(mass > 0)
    {
      *veldisp = sqrt(disp / (3 * mass));	/* convert to 1d velocity dispersion */

      spin[0] = lx / mass;
      spin[1] = ly / mass;
      spin[2] = lz / mass;

      parallel_sort(loc_rr_list, num_use, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

      nhalf = numtot / 2;
      mincpu = 0;

      while(nhalf >= npart[mincpu])
	{
	  nhalf -= npart[mincpu];
	  mincpu++;
	}

/*--- Here we still have to fix for individual masses, 
        maybe we even want the total mass for this ? 
        Note however that even within multi mass DM simulation
        the DM mass within a clean halo should be the same ... ------*/
      if(ThisTask == mincpu)
	*halfmassrad = loc_rr_list[nhalf].r;

      MPI_Bcast(halfmassrad, 1, MPI_DOUBLE, mincpu, MPI_COMM_WORLD);
    }
  else
    *veldisp = *halfmassrad = spin[0] = spin[1] = spin[2] = 0;

  /* compute cumulative mass */

  masslist = (double *) mymalloc("masslist", NTask * sizeof(double));

  for(i = 0, massloc = 0; i < num_use; i++)
    massloc += loc_rr_list[i].mass;

  MPI_Allgather(&massloc, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  for(i = 1; i < NTask; i++)
    masslist[i] += masslist[i - 1];

  for(i = 1; i < num_use; i++)
    loc_rr_list[i].mass += loc_rr_list[i - 1].mass;

  if(ThisTask > 0)
    for(i = 0; i < num_use; i++)
      loc_rr_list[i].mass += masslist[ThisTask - 1];

  for(i = 0, offset = 0; i < ThisTask; i++)
    offset += npart[i];

/*--- Note that here we might want to correct for the baryon fraction ? ---*/
  for(i = num_use - 1, max = 0, maxrad = 0; i + offset > 5 && i >= 0; i--)
    {
      if(loc_rr_list[i].r <= 0)
	endrun(124523);
      if(loc_rr_list[i].mass / loc_rr_list[i].r > max)
	{
	  max = loc_rr_list[i].mass / loc_rr_list[i].r;
	  maxrad = loc_rr_list[i].r;
	}
    }
  MPI_Allreduce(&max, vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if(max < *vmax)
    maxrad = 0;

  *vmax = sqrt(All.G * (*vmax));

  MPI_Allreduce(&maxrad, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  *vmaxrad = max;


  myfree(masslist);
  myfree(npart);
  if(num_use > 0)
    myfree(loc_rr_list);

}



int subfind_compare_serial_candidates_boundlength(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->bound_length > ((struct cand_dat *) b)->bound_length) {return -1;}
  if(((struct cand_dat *) a)->bound_length < ((struct cand_dat *) b)->bound_length) {return +1;}
  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank) {return -1;}
  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank) {return +1;}
  return 0;
}

int subfind_compare_serial_candidates_rank(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->rank < ((struct cand_dat *) b)->rank) {return -1;}
  if(((struct cand_dat *) a)->rank > ((struct cand_dat *) b)->rank) {return +1;}
  if(((struct cand_dat *) a)->len > ((struct cand_dat *) b)->len) {return -1;}
  if(((struct cand_dat *) a)->len < ((struct cand_dat *) b)->len) {return +1;}
  return 0;
}

int subfind_compare_serial_candidates_subnr(const void *a, const void *b)
{
  if(((struct cand_dat *) a)->subnr < ((struct cand_dat *) b)->subnr) {return -1;}
  if(((struct cand_dat *) a)->subnr > ((struct cand_dat *) b)->subnr) {return +1;}
  return 0;
}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float ran1(long *idum)
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if(*idum <= 0 || !iy)
    {
      if(-(*idum) < 1)
	*idum = 1;
      else
	*idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--)
	{
	  k = (*idum) / IQ;
	  *idum = IA * (*idum - k * IQ) - IR * k;
	  if(*idum < 0)
	    *idum += IM;
	  if(j < NTAB)
	    iv[j] = *idum;
	}
      iy = iv[0];
    }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if(*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}



float gasdev(long *idum)
{
  static int iset = 0;
  static float gset;
  float fac, rsq, v1, v2;

  if(iset == 0)
    {
      do
	{
	  v1 = 2.0 * ran1(idum) - 1.0;
	  v2 = 2.0 * ran1(idum) - 1.0;
	  rsq = v1 * v1 + v2 * v2;
	}
      while(rsq >= 1.0 || rsq == 0.0);
      fac = sqrt(-2.0 * log(rsq) / rsq);
      gset = v1 * fac;
      iset = 1;
      return v2 * fac;
    }
  else
    {
      iset = 0;
      return gset;
    }
}

#endif
