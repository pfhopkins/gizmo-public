#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "../declarations/allvars.h"
#include "../core/proto.h"

#if defined(OUTPUT_TWOPOINT_ENABLED)
#include <vector>
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "../system/gpu_particles_arena.h"
#endif

/*! \file twopoint.c
 *  \brief computes the two-point mass correlation function on the fly
 */
/*!
* This file was originally part of the GADGET3 code developed by Volker Springel.
* It has been updated by PFH for basic compatibility with GIZMO.
*/

/* Note: This routine will only work correctly for particles of equal mass ! */


#ifdef OUTPUT_TWOPOINT_ENABLED

#define BINS_TP  40		/* number of bins used */
#define ALPHA  -1.0		/* slope used in randomly selecting radii around target particles */

#ifndef FRACTION_TP
#define FRACTION_TP  0.2
#endif
/* above sets fraction of particles selected for sphere placement. Will be scaled with total particle number so that a fixed value should give roughly the same noise level in the meaurement, indpendent of simulation size */

struct twopointdata_in
{
  Vec3<MyDouble> Pos;
  MyFloat Rs;
  int NodeList[NODELISTLENGTH];
}
 *TwoPointDataIn, *TwoPointDataGet;


#define SQUARE_IT(x) ((x)*(x))


static long long Count[BINS_TP], Count_bak[BINS_TP];
static long long CountSpheres[BINS_TP];
static double Xi[BINS_TP];
static double Rbin[BINS_TP];

static double R0, R1;		/* inner and outer radius for correlation function determination */

static double logR0;
static double binfac;
static double PartMass;

static MyFloat *RsList;



/*  This function computes the two-point function.
 */
void twopoint(void)
{
    int i, j, k, bin, n, ndone, ndone_flag, dummy, nexport, nimport, place, recvTask, ngrp;
    double p, rs, vol, scaled_frac, tstart, tend, mass, masstot; long long *countbuf;
    PRINT_STATUS("begin two-point correlation function..."); tstart = my_second();
    /* set inner and outer radius for the bins that are used for the correlation function estimate */
    R0 = All.ForceSoftening[1] / 3.; R1 = All.BoxSize / 2; /* we assume that type=1 is the primary type */
    scaled_frac = FRACTION_TP * 1.0e7 / All.TotNumPart; logR0 = log(R0); binfac = BINS_TP / (log(R1) - log(R0));
    for(i = 0, mass = 0; i < NumPart; i++) {mass += P[i].Mass;}
    MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); PartMass = masstot / All.TotNumPart;
    for(i = 0; i < BINS_TP; i++) {Count[i] = 0; CountSpheres[i] = 0;}
    /* allocate buffers to arrange communication */
    RsList = (MyFloat *) mymalloc("RsList", NumPart * sizeof(MyFloat));

    gizmo_rng_t saved_rng = random_generator;
    gizmo_rng_init(&random_generator, (uint64_t)P[0].ID + (uint64_t)ThisTask);
    /* Modern path: prebuilt CSR NL, per-source rs, all-types search pool.
     *
     * INTERIM IMPLEMENTATION: drops the legacy tree-mass-aggregation
     * optimization. The legacy used the gravity tree's node-level mass
     * moments — when an entire tree node fit cleanly in a single distance
     * bin, it counted (node->mass / PartMass) particles in one shot without
     * enumerating individuals. This makes large-radius bins effectively free.
     * The modern symmetric NL walks individual pair lists with no node-level
     * aggregation; correct but much slower at large radii.
     *
     * FUTURE PROPER PORT: this routine should NOT use the normal symmetric
     * NL infrastructure at all. Two-point correlation is fundamentally an
     * "all-against-all long-range" problem with binned distance accumulators
     * — exactly what the gravity tree already supports. The proper port
     * walks the modern GRAVITY TREE (Nodes_base / SoA in gravity/forcetree.cc)
     * with a custom op that, instead of accumulating gravitational potential
     * via mass moments, accumulates BINNED PAIR COUNT via node->mass for
     * nodes that fall cleanly within a single radial bin. The walk pattern
     * is identical to gravtree; only the accumulated quantity differs.
     * That port preserves the O(N log N) scaling at large radii.
     *
     * The user signed off on shipping the interim NL-based implementation
     * because OUTPUT_TWOPOINT_ENABLED runs ~1x per 1e4-1e5 timesteps, so the
     * O(N) per-pair cost is tolerable for now. The proper gravity-tree port
     * is needed before any production run with TotNumPart > ~1e8 enables this
     * module.
     *
     * Pair-counting in the interim implementation is global with no
     * double-count: each pair (i, j) is counted once on the rank that owns i
     * (j may be local home gas or imported ghost from another rank's home;
     * the OTHER rank doesn't have i as a sampled source, so the pair appears
     * exactly once). Self pair (j == i) skipped. */
    {
        std::vector<int> active_idx;
        std::vector<double> active_radii;
        active_idx.reserve(NumPart);
        active_radii.reserve(NumPart);

        for(i = 0; i < NumPart; i++) {
            if(gizmo_rng_uniform(&random_generator) < scaled_frac) {
                p = gizmo_rng_uniform(&random_generator);
                rs = pow(pow(R0, ALPHA) + p * (pow(R1, ALPHA) - pow(R0, ALPHA)), 1 / ALPHA);
                bin = (int) ((log(rs) - logR0) * binfac);
                rs = exp((bin + 1) / binfac + logR0);
                RsList[i] = rs;
                for(j = 0; j <= bin; j++) {CountSpheres[j]++;}
                active_idx.push_back(i);
                active_radii.push_back((double)rs);
            }
        }

        int num_src = (int)active_idx.size();
        int num_src_global = 0;
        MPI_Allreduce(&num_src, &num_src_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        gpu_neighbor_list_t gnl = {};
        std::vector<int> gnl_neighbors_host;
        int local_count = ghost_get_num_local();
        if(local_count <= 0) local_count = NumPart;
        int ghost_imported = 0;
        if(num_src_global > 0) {
            int need_import_local = (ghost_get_num_ghosts() <= 0) ? 1 : 0;
            int need_import = 0;
            MPI_Allreduce(&need_import_local, &need_import, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if(need_import) {
                if(ghost_get_num_ghosts() > 0) ghost_exchange_cleanup();
                gizmo_density_prep_ghosts(gizmo_ghost_safety_factor());
                ghost_imported = 1;
            }
            int num_all = local_count + ghost_get_num_ghosts();
            if(num_all <= 0) num_all = NumPart;
            gpu_particles_arena_acquire(num_all, P, CellP);
            struct particle_data *P_gpu = gpu_particles_arena_P();
            /* All-types search pool: bitmask = 0xFF covers all 6 standard types. */
            gpu_ngb_list_build(P_gpu, num_all,
                               active_idx.data(), num_src,
                               NGB_SEARCH_ONEWAY, 0xFF /* all types */,
                               &gnl, NULL, 1.0, active_radii.data());
            /* gnl.neighbors is DEVICE_SPACE; host loop below indexes it. */
            if(gnl.total_pairs > 0) {
                gnl_neighbors_host.resize((size_t)gnl.total_pairs);
                gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
            }
        }
        const int *gnl_neighbors = gnl_neighbors_host.empty() ? NULL : gnl_neighbors_host.data();

        for(int aa = 0; aa < num_src; aa++) {
            int isrc = active_idx[aa];
            double rs_local = active_radii[aa];
            double rs2 = rs_local * rs_local;
            Vec3<double> pos_i = P[isrc].Pos;
            int64_t n_off = gnl.offsets[aa], n_off_end = gnl.offsets[aa+1];
            for(int64_t nn = n_off; nn < n_off_end; nn++) {
                int jp = gnl_neighbors[nn];
                if(jp == isrc) continue; /* skip self */
                double dx_raw = P[jp].Pos[0] - pos_i[0];
                double dy_raw = P[jp].Pos[1] - pos_i[1];
                double dz_raw = P[jp].Pos[2] - pos_i[2];
                MyDouble xtmp = 0;
                double dx = NGB_PERIODIC_BOX_LONG_X(dx_raw, dy_raw, dz_raw, -1);
                double dy = NGB_PERIODIC_BOX_LONG_Y(dx_raw, dy_raw, dz_raw, -1);
                double dz = NGB_PERIODIC_BOX_LONG_Z(dx_raw, dy_raw, dz_raw, -1);
                double r2 = dx*dx + dy*dy + dz*dz;
                if(r2 >= R0 * R0 && r2 < R1 * R1 && r2 < rs2) {
                    int bin_local = (int)((log(sqrt(r2)) - logR0) * binfac);
                    if(bin_local < BINS_TP) Count[bin_local]++;
                }
            }
        }

        if(num_src_global > 0) {
            gpu_ngb_list_free(&gnl, NULL);
            gpu_particles_arena_invalidate();
        }
        if(ghost_imported) ghost_exchange_cleanup();
    }
    random_generator = saved_rng;

    myfree(RsList);

    /* Now compute the actual correlation function */
    countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_TP * sizeof(long long));
    MPI_Allgather(Count, BINS_TP * sizeof(long long), MPI_BYTE, countbuf, BINS_TP * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

    for(i = 0; i < BINS_TP; i++) {Count[i] = 0; for(n = 0; n < NTask; n++) {Count[i] += countbuf[n * BINS_TP + i];}}
    MPI_Allgather(CountSpheres, BINS_TP * sizeof(long long), MPI_BYTE, countbuf, BINS_TP * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

    for(i = 0; i < BINS_TP; i++) {CountSpheres[i] = 0; for(n = 0; n < NTask; n++) {CountSpheres[i] += countbuf[n * BINS_TP + i];}}
    myfree(countbuf);

    for(i = 0; i < BINS_TP; i++)
      {
        vol = 4 * M_PI / 3.0 * (pow(exp((i + 1.0) / binfac + logR0), 3) - pow(exp((i + 0.0) / binfac + logR0), 3));
        if(CountSpheres[i] > 0) {Xi[i] = -1 + Count[i] / ((double) CountSpheres[i]) / (All.TotNumPart / pow(All.BoxSize, 3)) / vol;} else {Xi[i] = 0;}
        Rbin[i] = exp((i + 0.5) / binfac + logR0);
      }
    twopoint_save();
    tend = my_second(); PRINT_STATUS(" ..end two-point: Took=%g seconds", timediff(tstart, tend));
}




void twopoint_save(void)
{
  FILE *fd; char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE]; int i;
  if(ThisTask == 0)
    {
      snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/correl_%03d.txt", All.OutputDir, RestartSnapNum);
      if(!(fd = fopen(buf, "w"))) {printf("can't open file `%s`\n", buf); endrun(1323);}
      fprintf(fd, "%.16g\n", All.Time); i = BINS_TP; fprintf(fd, "%d\n", i);
      for(i = 0; i < BINS_TP; i++) {fprintf(fd, "%g %g %g %g\n", Rbin[i], Xi[i], (double) Count[i], (double) CountSpheres[i]);}
      fclose(fd);
    }
}




/*! This function counts the pairs in a sphere
 */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here! -- */
int twopoint_count_local(int target, int mode, int *nexport, int *nsend_local)
{
  int startnode, listindex = 0;
  MyDouble *pos;
  MyFloat rs;

  if(mode == 0)
    {
      pos = P[target].Pos.data_ptr();
      rs = RsList[target];
      memcpy(Count_bak, Count, sizeof(long long) * BINS_TP);
    }
  else
    {
      pos = TwoPointDataGet[target].Pos.data_ptr();
      rs = TwoPointDataGet[target].Rs;
    }


  /* Now start the actual tree-walk for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = TwoPointDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  if(twopoint_ngb_treefind_variable(pos, rs, target, &startnode, mode, nexport, nsend_local) < 0)
	    {
	      /* in this case restore the count-table */
	      memcpy(Count, Count_bak, sizeof(long long) * BINS_TP); {return -1;} /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = TwoPointDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  return 0;
}





/*! This function finds all particles within the radius "rsearch", and counts them in the bins used for the two-point correlation function.
 *    this is a custom version of "ngb_treefind_variable", hard-coded for a square box (no shearing!) and variable search threshold radii, 
 *    bin-dumping, etc. as a result, updates to the core neighbor search routine will not alter this subroutine
 */
/*!   -- this subroutine is not openmp parallelized at present, so there's not any issue about conflicts over shared memory. if you make it openmp, make sure you protect the writes to shared memory here! -- */
int twopoint_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat rsearch, int target, int *startnode, int mode, int *nexport, int *nsend_local)
{
  double r2, r, ri, ro;
  int no, p, bin, task, bin2, nexport_save;
  struct NODE *current;
  MyDouble dx, dy, dz, dist;
  MyDouble xtmp; xtmp=0;

  nexport_save = *nexport;

  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  dx = NGB_PERIODIC_BOX_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  dy = NGB_PERIODIC_BOX_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  dz = NGB_PERIODIC_BOX_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);

	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 >= R0 * R0 && r2 < R1 * R1)
	    {
	      if(r2 < rsearch * rsearch)
		{
		  bin = (int) ((log(sqrt(r2)) - logR0) * binfac);
		  if(bin < BINS_TP)
		    Count[bin]++;
		}
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes + MaxForeignNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(123127);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes + MaxForeignNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= All.BunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0) {endrun(13004);}	/* in this case, the buffer is too small to process even a single particle */
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
		    DomainNodeIndex[no - (All.MaxPart + MaxNodes + MaxForeignNodes)];

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
		  return 0;
		}
	    }

	  no = current->u.d.sibling;	/* make skipping the branch the default */

	  dist = rsearch + 0.5 * current->len;
	  dx = NGB_PERIODIC_BOX_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dx > dist) continue;
	  dy = NGB_PERIODIC_BOX_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dy > dist) continue;
	  dz = NGB_PERIODIC_BOX_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dz > dist) continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += CUBE_EDGEFACTOR_1 * current->len;
	  if((r2 = dx * dx + dy * dy + dz * dz) > dist * dist) continue;

	  r = sqrt(r2);

	  ri = r - CUBE_EDGEFACTOR_2 * current->len;
	  ro = r + CUBE_EDGEFACTOR_2 * current->len;

	  if(ri >= R0 && ro < R1)
	    {
	      if(ro < rsearch)
		{
		  bin = (int) ((log(ri) - logR0) * binfac);
		  bin2 = (int) ((log(ro) - logR0) * binfac);
		  if(bin == bin2)
		    {
		      if(mode == 1)
			{
			  if((current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL)))
			    continue;
			}
		      Count[bin] += (long long)(current->u.d.mass / PartMass);
		      continue;
		    }
		}
	    }

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return 0;
}




#endif
