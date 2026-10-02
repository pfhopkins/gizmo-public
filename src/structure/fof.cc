#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <inttypes.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../core/wakeup_sidecar.h"
#include "../declarations/gpu_rng.h"
#include "../mesh/neighbor_list.h"
#include "group_search.h"

/*! \file fof.c
 *  \brief parallel FoF group finder
 */
/*!
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
  Vec3<MyDouble> Pos;
  MyFloat KernelRadius;
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
static float *fof_nearest_rkern;


struct fof_label_t
{
  MyIDType id;
  MyIDType task;
};

struct fof_remote_edge_t
{
  int local_index;
  int remote_task;
  int remote_index;
};

struct fof_label_update_t
{
  int index;
  MyIDType id;
  MyIDType task;
};

static inline int fof_label_less(const fof_label_t &a, const fof_label_t &b)
{
  if(a.id < b.id) return 1;
  if(a.id > b.id) return 0;
  return a.task < b.task;
}

static int fof_dsu_find(std::vector<int> &parent, int x)
{
  while(parent[x] != x)
    {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
  return x;
}

static void fof_dsu_union(std::vector<int> &parent, std::vector<int> &rank, std::vector<fof_label_t> &root_label, int a, int b)
{
  int ra = fof_dsu_find(parent, a);
  int rb = fof_dsu_find(parent, b);
  if(ra == rb) return;
  if(rank[ra] < rank[rb]) std::swap(ra, rb);
  parent[rb] = ra;
  if(rank[ra] == rank[rb]) rank[ra]++;
  if(fof_label_less(root_label[rb], root_label[ra])) root_label[ra] = root_label[rb];
}

static void fof_exchange_label_updates(const std::vector<fof_label_update_t> &send_updates,
                                       const std::vector<int> &send_count,
                                       std::vector<fof_label_update_t> &recv_updates)
{
  std::vector<int> recv_count(NTask, 0), send_disp(NTask, 0), recv_disp(NTask, 0);
  MPI_Alltoall((void *)send_count.data(), 1, MPI_INT, recv_count.data(), 1, MPI_INT, MPI_COMM_WORLD);
  int total_send = 0, total_recv = 0;
  for(int task = 0; task < NTask; task++)
    {
      if(task > 0)
        {
          send_disp[task] = send_disp[task - 1] + send_count[task - 1];
          recv_disp[task] = recv_disp[task - 1] + recv_count[task - 1];
        }
      total_send += send_count[task];
      total_recv += recv_count[task];
    }
  recv_updates.resize(total_recv);

  std::vector<int> send_bytes(NTask), recv_bytes(NTask), send_bdisp(NTask), recv_bdisp(NTask);
  for(int task = 0; task < NTask; task++)
    {
      send_bytes[task] = send_count[task] * (int)sizeof(fof_label_update_t);
      recv_bytes[task] = recv_count[task] * (int)sizeof(fof_label_update_t);
      send_bdisp[task] = send_disp[task] * (int)sizeof(fof_label_update_t);
      recv_bdisp[task] = recv_disp[task] * (int)sizeof(fof_label_update_t);
    }
  MPI_Alltoallv((void *)(total_send > 0 ? send_updates.data() : NULL), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                (void *)(total_recv > 0 ? recv_updates.data() : NULL), recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
}

static void fof_find_groups_modern(void)
{
  std::vector<int> primary_indices;
  primary_indices.reserve(NumPart);
  for(int i = 0; i < NumPart; i++)
    if(P[i].Mass > 0 && ((1 << P[i].Type) & MyFOF_PRIMARY_LINK_TYPES))
      primary_indices.push_back(i);

  if(primary_indices.empty()) return;

  std::vector<int> local_to_primary(NumPart, -1);
  std::vector<int> parent(primary_indices.size()), rank(primary_indices.size(), 0);
  std::vector<fof_label_t> root_label(primary_indices.size());
  for(size_t a = 0; a < primary_indices.size(); a++)
    {
      int i = primary_indices[a];
      local_to_primary[i] = a;
      parent[a] = a;
      root_label[a].id = P[i].ID;
      root_label[a].task = ThisTask;
    }

  group_search_import_pool_t import_pool;
  group_search_import_particles(MyFOF_PRIMARY_LINK_TYPES, MyFOF_PRIMARY_LINK_TYPES, LinkL, NULL, NULL, import_pool);

  std::vector<double> radii(primary_indices.size(), LinkL);
  neighbor_list_t nl = {NULL, NULL, 0, 0};
  group_search_build_cross_type_nl(import_pool.particles, primary_indices.data(), primary_indices.size(),
                                   radii.data(), MyFOF_PRIMARY_LINK_TYPES, &nl);

  std::vector<fof_remote_edge_t> remote_edges;
  for(size_t a = 0; a < primary_indices.size(); a++)
    {
      int i = primary_indices[a];
      for(int n = nl.offsets[a]; n < nl.offsets[a + 1]; n++)
        {
          int j = nl.neighbors[n];
          if(j == i) continue;
          if(group_search_distance2_particles(import_pool.particles[i], import_pool.particles[j]) >= LinkL * LinkL) continue;
          if(j < NumPart)
            {
              int b = (j >= 0 && j < NumPart) ? local_to_primary[j] : -1;
              if(b >= 0) fof_dsu_union(parent, rank, root_label, a, b);
            }
          else
            {
              int g = j - NumPart;
              fof_remote_edge_t edge;
              edge.local_index = i;
              edge.remote_task = import_pool.ghost_home_task[g];
              edge.remote_index = import_pool.ghost_home_index[g];
              remote_edges.push_back(edge);
            }
        }
    }
  free_neighbor_list(&nl);

  long long n_remote_edges_local = remote_edges.size(), n_remote_edges_global = 0;
  MPI_Allreduce(&n_remote_edges_local, &n_remote_edges_global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  PRINT_STATUS("FOF modern: local primary components linked; remote boundary edges=%lld", n_remote_edges_global);

  int changed_global = 0;
  int iter = 0;
  do
    {
      std::vector<int> send_count(NTask, 0);
      for(size_t e = 0; e < remote_edges.size(); e++) send_count[remote_edges[e].remote_task]++;
      std::vector<int> cursor(NTask, 0), send_disp(NTask, 0);
      int total_send = 0;
      for(int task = 0; task < NTask; task++)
        {
          if(task > 0) send_disp[task] = send_disp[task - 1] + send_count[task - 1];
          total_send += send_count[task];
        }
      cursor = send_disp;
      std::vector<fof_label_update_t> send_updates(total_send);
      for(size_t e = 0; e < remote_edges.size(); e++)
        {
          const fof_remote_edge_t &edge = remote_edges[e];
          int a = local_to_primary[edge.local_index];
          int root = fof_dsu_find(parent, a);
          int off = cursor[edge.remote_task]++;
          send_updates[off].index = edge.remote_index;
          send_updates[off].id = root_label[root].id;
          send_updates[off].task = root_label[root].task;
        }

      std::vector<fof_label_update_t> recv_updates;
      fof_exchange_label_updates(send_updates, send_count, recv_updates);

      int changed_local = 0;
      for(size_t u = 0; u < recv_updates.size(); u++)
        {
          int idx = recv_updates[u].index;
          if(idx < 0 || idx >= NumPart) continue;
          int a = local_to_primary[idx];
          if(a < 0) continue;
          int root = fof_dsu_find(parent, a);
          fof_label_t incoming = {recv_updates[u].id, recv_updates[u].task};
          if(fof_label_less(incoming, root_label[root]))
            {
              root_label[root] = incoming;
              changed_local = 1;
            }
        }
      MPI_Allreduce(&changed_local, &changed_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      iter++;
      if(iter > MAXITER)
        {
          printf("FOF modern distributed label propagation failed to converge after %d iterations\n", iter);
          fflush(stdout);
          endrun(990503);   /* FOF label-propagation non-convergence: symmetric (iter is lock-stepped by the shared changed_global Allreduce), so all ranks reach the poll together */
          gizmo_exit_bad_stop_if_requested("fof:label_prop_nonconverged");
        }
    }
  while(changed_global);

  PRINT_STATUS("FOF modern: distributed labels converged in %d iterations", iter);

  for(size_t a = 0; a < primary_indices.size(); a++)
    {
      int i = primary_indices[a];
      int root = fof_dsu_find(parent, a);
      MinID[i] = root_label[root].id;
      MinIDTask[i] = root_label[root].task;
      Head[i] = i;
    }
}

static void fof_find_nearest_dmparticle_modern(void)
{
  std::vector<int> secondary_indices;
  secondary_indices.reserve(NumPart);
  for(int i = 0; i < NumPart; i++)
    if(P[i].Mass > 0 && ((1 << P[i].Type) & MyFOF_SECONDARY_LINK_TYPES))
      secondary_indices.push_back(i);

  if(secondary_indices.empty()) return;

  group_search_import_pool_t import_pool;
  group_search_import_particles(MyFOF_SECONDARY_LINK_TYPES, MyFOF_PRIMARY_LINK_TYPES, 4.0 * LinkL,
                                MinID, MinIDTask, import_pool);

  std::vector<double> nearest_dist2(secondary_indices.size(), MAX_REAL_NUMBER);
  std::vector<unsigned char> done(secondary_indices.size(), 0);
  std::vector<double> radius(secondary_indices.size(), 0.1 * LinkL);
  std::vector<int> secondary_pos_by_index(NumPart, -1);
  for(size_t a = 0; a < secondary_indices.size(); a++) secondary_pos_by_index[secondary_indices[a]] = a;

  int nleft_global = 0;
  int iter = 0;
  do
    {
      std::vector<int> active_sources;
      std::vector<double> active_radii;
      active_sources.reserve(secondary_indices.size());
      active_radii.reserve(secondary_indices.size());
      for(size_t a = 0; a < secondary_indices.size(); a++)
        if(!done[a])
          {
            active_sources.push_back(secondary_indices[a]);
            active_radii.push_back(radius[a]);
          }

      if(active_sources.empty()) break;

      neighbor_list_t nl = {NULL, NULL, 0, 0};
      group_search_build_cross_type_nl(import_pool.particles, active_sources.data(), active_sources.size(),
                                       active_radii.data(), MyFOF_PRIMARY_LINK_TYPES, &nl);

      for(size_t aa = 0; aa < active_sources.size(); aa++)
        {
          int i = active_sources[aa];
          int sec_pos = secondary_pos_by_index[i];
          if(sec_pos < 0) continue;
          double best_r2 = nearest_dist2[sec_pos];
          int best_j = -1;
          for(int n = nl.offsets[aa]; n < nl.offsets[aa + 1]; n++)
            {
              int j = nl.neighbors[n];
              double r2 = group_search_distance2_particles(import_pool.particles[i], import_pool.particles[j]);
              if(r2 < best_r2 && r2 < active_radii[aa] * active_radii[aa])
                {
                  best_r2 = r2;
                  best_j = j;
                }
            }
          if(best_j >= 0)
            {
              nearest_dist2[sec_pos] = best_r2;
              done[sec_pos] = 1;
              if(best_j < NumPart)
                {
                  MinID[i] = MinID[best_j];
                  MinIDTask[i] = MinIDTask[best_j];
                }
              else
                {
                  int g = best_j - NumPart;
                  MinID[i] = import_pool.ghost_label_id[g];
                  MinIDTask[i] = import_pool.ghost_label_task[g];
                }
              Head[i] = i;
            }
        }
      free_neighbor_list(&nl);

      int nleft_local = 0;
      for(size_t a = 0; a < secondary_indices.size(); a++)
        {
          if(done[a]) continue;
          if(radius[a] < 4.0 * LinkL)
            {
              radius[a] *= 2.0;
              if(radius[a] > 4.0 * LinkL) radius[a] = 4.0 * LinkL;
              nleft_local++;
            }
          else
            {
              done[a] = 1;
            }
        }
      MPI_Allreduce(&nleft_local, &nleft_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      iter++;
      if(iter > MAXITER)
        {
          printf("FOF modern nearest-primary attachment failed to converge after %d iterations\n", iter);
          fflush(stdout);
          endrun(990504);   /* FOF nearest-primary non-convergence: symmetric (iter lock-stepped by the shared nleft_global Allreduce), all ranks reach the poll together */
          gizmo_exit_bad_stop_if_requested("fof:nearest_primary_nonconverged");
        }
    }
  while(nleft_global > 0);
}



void fof_fof(int num)
{
  int i, ndm, start, lenloc, largestgroup, n = 0;
  double mass, masstot, rhodm, t0, t1;
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

#ifdef SUBFIND
  /* SUBFIND builds its own trees, and it runs after the allocations below. If the gravity tree from
     the last domain decomposition is still standing at that point, its DomainNodeIndex and
     TopNodeNodeIndex sit underneath those allocations: SUBFIND's own force_treeallocate overwrites
     the pointers to them, so nothing can release them and they strand above the domain blocks, where
     they block the next domain_free(). Release that tree here instead, before anything is allocated
     on top of it -- the same state group finding already runs in when it is invoked on a snapshot
     rather than mid-step. Matched to the call condition below, so a build that skips SUBFIND keeps
     its tree; force_treefree() does nothing when no tree is allocated, and SUBFIND discards the tree
     before returning in any case.
     TreeReconstructFlag says the tree is gone: the timestepper otherwise reuses it, updating it in
     place with the last step's kicks (run.cc) and skipping construction before the next walk
     (gravtree.cc), both of which would read storage that no longer exists.  DomainReconstructFlag
     says the same thing to the decomposition, which is what holds the tree's storage: a rebuild on
     its own would build into arrays this just released. */
  if(num >= 0 && DumpFlag != 2) {force_treefree(); TreeReconstructFlag = 1; DomainReconstructFlag = 1;}
#endif

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
    printf("Modern FOF neighbor-list construction.\n");

  (void)n;

  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i] = 1;
      Next[i] = -1;
      MinID[i] = P[i].ID;
      MinIDTask[i] = ThisTask;
    }


  t0 = my_second();

  fof_find_groups_modern();

  t1 = my_second();
  if(ThisTask == 0)
    printf("group finding took = %g sec\n", timediff(t0, t1));


  t0 = my_second();

  fof_find_nearest_dmparticle_modern();

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
	  endrun(8812);		/* corrupt label (MinIDTask out of range): soft bad-stop */
	  Next[i] = P[i].ID;	/* reset this entry to a consistent LOCAL singleton: own ID as MinID
				   (Next feeds FOF_PList[i].MinID below) + own task (Tail, next line), matching
				   the singleton convention set at init, so the later FOF_PList fill and the
				   Send_count[FOF_GList[].MinIDTask]++ in fof_compile_catalogue stay in-bounds
				   before the run drains at the caller phase poll. NOT a meaningful clamp --
				   the run is already stopping via the flag. */
	  Tail[i] = ThisTask;
	}
    }

  for(i = 0; i < NumPart; i++)
    {
      FOF_PList[i].MinID = Next[i];
      FOF_PList[i].MinIDTask = Tail[i];
      FOF_PList[i].Pindex = i;
    }

  
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
      int bad_match = 0;
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NumPart)
	    { endrun(78); bad_match = 1; break; }	/* FOF_PList overrun: soft bad-stop, stop indexing past the array */
	}
      if(bad_match) continue;

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
	{ endrun(123); continue; }	/* group MinID mismatch: soft bad-stop, skip this group */

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

#ifdef SINK_SEED_FROM_FOF
  if(num < 0){   // Make BHs in every call to fof_fof (including the group finding for each snapshot)
      if(All.Time < 1.0/(1.0+All.SeedSinkMinRedshift)) { fof_make_sink_particles(); } else {  printf("skipping sink particle seeding at a = %g \n", All.Time); }
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

}



/* fof_find_groups, fof_find_dmparticles_evaluate, fof_find_nearest_dmparticle,
   fof_find_nearest_dmparticle_evaluate — retired in D2.5.  Modern replacements:
   fof_find_groups_modern() + fof_find_nearest_dmparticle_modern() above. */


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
      int bad_match = 0;
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    { endrun(7973); bad_match = 1; break; }	/* FOF_GList overrun: soft bad-stop */
	}
      if(bad_match) continue;

      if(get_FOF_GList[i].LocCount != 0)
	{ endrun(123); continue; }	/* imported group should have LocCount==0: soft bad-stop, skip */

      if(FOF_GList[start].MinIDTask != (MyIDType)ThisTask)
	{ endrun(124); continue; }	/* MinIDTask mismatch: soft bad-stop, skip */

      FOF_GList[start].ExtCount += get_FOF_GList[i].ExtCount;
#ifdef FOF_DENSITY_SPLIT_TYPES
      FOF_GList[start].ExtDMCount += get_FOF_GList[i].ExtDMCount;
#endif
    }

  /* copy the size information back into the list, to inform the others */
  for(i = 0, start = 0; i < nimport; i++)
    {
      int bad_match = 0;
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    { endrun(797831); bad_match = 1; break; }	/* FOF_GList overrun: soft bad-stop */
	}
      if(bad_match) continue;

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
#ifdef SINK_PARTICLES
  Group[gr].Sink_Mass = 0;
  Group[gr].Sink_Mdot = 0;
#ifdef SINK_SEED_FROM_FOF
  Group[gr].MinPot = SINK_MINPOTVALUE_INIT;
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
	Group[gr].Sfr += CellP[index].Sfr;
#endif

#ifdef SINK_PARTICLES
      if(P[index].Type == 5)
	{
	  Group[gr].Sink_Mdot += P[index].Sink_Mdot;
	  Group[gr].Sink_Mass += P[index].Sink_Mass;
	}

#ifdef SINK_SEED_FROM_FOF
#if (SINK_SEED_FROM_FOF==0)
      if(P[index].Type==0)
#elif (SINK_SEED_FROM_FOF==1)
      if(P[index].Type==4)
#endif
         if(P[index].Potential < Group[gr].MinPot)
        {
          Group[gr].MinPot = P[index].Potential;
          Group[gr].index_maxdens = index;
          Group[gr].task_maxdens = ThisTask;
        }
#endif
#endif // SINK_PARTICLES

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
      int bad_match = 0;
      while(Group[start].MinID < get_Group[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    { endrun(797890); bad_match = 1; break; }	/* Group overrun: soft bad-stop */
	}
      if(bad_match) continue;

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
#ifdef SINK_PARTICLES
      Group[start].Sink_Mdot += get_Group[i].Sink_Mdot;
      Group[start].Sink_Mass += get_Group[i].Sink_Mass;
#ifdef SINK_SEED_FROM_FOF
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
    endrun(876889);	/* group-count invariant off: soft bad-stop; Ngroups<=NgroupsExt so the qsort below stays in-bounds, drains at the caller phase poll */

  qsort(Group, Ngroups, sizeof(group_properties), fof_compare_Group_MinID);
}



void fof_save_groups(int num)
{
  int i, j, start, lenloc, nprocgroup, primaryTask, groupTask, ngr, totlen;
  long long totNids;
  char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
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
      endrun(123123);	/* symmetric (i is the post-Allreduce(SUM) global count, TotNgroups is global): all ranks reach the poll together */
      gizmo_exit_bad_stop_if_requested("fof:save_groupcount_mismatch");
    }

  /* bring the group list back into the original order */
  parallel_sort(FOF_GList, NgroupsExt, sizeof(fof_group_list), fof_compare_FOF_GList_ExtCountMinID);

  /* Assign the group numbers to the group properties array */
  for(i = 0, start = 0; i < Ngroups; i++)
    {
      int bad_match = 0;
      while(FOF_GList[start].MinID < Group[i].MinID)
	{
	  start++;
	  if(start >= NgroupsExt)
	    { endrun(7297890); bad_match = 1; break; }	/* FOF_GList overrun: soft bad-stop */
	}
      if(bad_match) continue;
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
      int bad_match = 0;
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
	{
	  start++;
	  if(start >= NumPart)
	    { endrun(78); bad_match = 1; break; }	/* FOF_PList overrun: soft bad-stop */
	}
      if(bad_match) continue;

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
	{ endrun(1313); continue; }	/* group MinID mismatch: soft bad-stop, skip this group */

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
      endrun(12);	/* symmetric (totNids is the post-sumup_large_ints global sum, TotNids is global): all ranks reach the poll together */
      gizmo_exit_bad_stop_if_requested("fof:save_totnids_mismatch");
    }

  /* sort the particle IDs according to group-number */
  parallel_sort(ID_list, Nids, sizeof(fof_id_list), fof_compare_ID_list_GrNrID);

  t1 = my_second();
  PRINT_STATUS("Group catalogues globally sorted. took = %g sec. Started saving of group catalogue", timediff(t0, t1));
  t0 = my_second();

  if(ThisTask == 0)
    {
      snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d", All.OutputDir, num);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);


  if(NTask < All.NumFilesWrittenInParallel) {if(ThisTask==0) {printf("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");} endrun(241931); gizmo_exit_bad_stop_if_requested("fof:save_numfiles_too_large");}	/* symmetric (NTask + NumFilesWrittenInParallel global), just after the all-rank barrier above */

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  primaryTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (primaryTask + groupTask))	/* ok, it's this processor's turn */
	fof_save_local_catalogue(num);
      gizmo_exit_bad_stop_if_requested("fof:save_turn");	/* all-rank drain (replaces the in-group barrier): catches per-rank save-IO faults */
    }

  myfree(ID_list);

  t1 = my_second();
  PRINT_STATUS("Group catalogues saved. took = %g sec", timediff(t0, t1));
}



void fof_save_local_catalogue(int num)
{
  FILE *fd;
  float *mass, *cm, *vel;
  char fname[DEFAULT_PATH_BUFFERSIZE_TOUSE];
  int i, j, *len;
  MyIDType *ids;

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);
  if(!(fd = fopen(fname, "w")))
    {
      printf("can't open file `%s`\n", fname);
      endrun(1183); return;	/* per-rank group_tab write-open fail: soft bad-stop, skip this rank's writes; the per-turn poll in fof_save_groups drains all ranks */
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

#ifdef SINK_PARTICLES
  /* Sink_Mass */
  mass = (float *)mymalloc("mass", Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].Sink_Mass;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);

  /* Sink_Mdot */
  mass = (float *)mymalloc("mass", Ngroups * sizeof(float));
  for(i = 0; i < Ngroups; i++)
    mass[i] = Group[i].Sink_Mdot;
  my_fwrite(mass, Ngroups, sizeof(float), fd);
  myfree(mass);
#endif

  fclose(fd);


  ids = (MyIDType *) ID_list;
  for(i = 0; i < Nids; i++)
    ids[i] = ID_list[i].ID;

  snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, ThisTask);
  if(!(fd = fopen(fname, "w")))
    {
      printf("can't open file `%s`\n", fname);
      endrun(1184); return;	/* per-rank group_ids write-open fail: soft bad-stop, skip; group_tab already closed; per-turn poll in fof_save_groups drains */
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




#ifdef SINK_SEED_FROM_FOF

void fof_make_sink_particles(void)
{
  int i, j, n, ntot;
  long nexport, nimport;
  int recvTask, level;
  int *import_indices, *export_indices;
  double unitmass_in_msun;

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
#if (SINK_SEED_FROM_FOF==0)
    if(Group[i].MassType[1] >= (All.OmegaMatter - All.OmegaBaryon) / All.OmegaMatter * All.MinFoFMassForNewSeed)
#elif (SINK_SEED_FROM_FOF==1)
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

  import_indices = (int *)mymalloc("import_indices", nimport * sizeof(int));
  export_indices = (int *)mymalloc("export_indices", nexport * sizeof(int));

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
#if (SINK_SEED_FROM_FOF==0)
        if(Group[i].MassType[1] >= (All.OmegaMatter - All.OmegaBaryon) / All.OmegaMatter * All.MinFoFMassForNewSeed)
#elif (SINK_SEED_FROM_FOF==1)
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
  PRINT_STATUS("Making %d new sink particle particles", ntot);
  All.TotSinks += ntot;

  for(n = 0; n < nimport; n++)
    {
#if (SINK_SEED_FROM_FOF==0)
        if(P[import_indices[n]].Type != 0)
#elif (SINK_SEED_FROM_FOF==1)
            if(P[import_indices[n]].Type != 4)
#endif
                { endrun(7772); continue; }	/* imported particle has wrong Type for a FOF sink: soft bad-stop, skip converting it (avoids bad P[] mutation); drains at caller phase poll */
        
        P[import_indices[n]].Mass = CellP[import_indices[n]].Mass; /* sync mass before type conversion */
        P[import_indices[n]].Type = 5;    /* make it a sink particle particle */
        force_tree_note_type_presence(import_indices[n]);  /* the element's type is final here, so tell the standing tree before anything walks it looking for that type */
#ifdef GALSF
        P[import_indices[n]].StellarAge = All.Time; /* reset formation time to match BH formation */
#endif
        /* generate BH mass */
        if(All.SeedSinkMassSigma > 0)
        {
            /* compute gaussian random number: mean=0, sigma=All.SeedSinkMassSigma */
            P[import_indices[n]].Sink_Mass = pow( 10., log10(All.SeedSinkMass) + gizmo_gpu_rand_gaussian((uint64_t)P[import_indices[n]].ID, (uint64_t)(17 + All.NumCurrentTiStep)) * All.SeedSinkMassSigma );
            unitmass_in_msun = UNIT_MASS_IN_SOLAR;
            if( P[import_indices[n]].Sink_Mass < 100./unitmass_in_msun )
                P[import_indices[n]].Sink_Mass = 100./unitmass_in_msun;      // enforce lower limit of Mseed = 100 x Msun
        } else {
            P[import_indices[n]].Sink_Mass = All.SeedSinkMass;
        }
        P[import_indices[n]].Sink_Mdot = 0;
        /* set hydro-ish variables */
        if(P[import_indices[n]].Type == 0){
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            P[import_indices[n]].Mass = CellP[import_indices[n]].MassTrue + CellP[import_indices[n]].dMass;
#endif
            P[import_indices[n]].DensityAroundParticle = CellP[import_indices[n]].Density;
        }
        /* set some specific BH variables that are needed below */
#ifdef SINK_INCREASE_DYNAMIC_MASS
        P[import_indices[n]].Mass *= SINK_INCREASE_DYNAMIC_MASS;
#endif
#ifdef SINK_ALPHADISK_ACCRETION
        P[import_indices[n]].Sink_Mass_Reservoir = All.SeedReservoirMass;
#endif
#ifdef SINK_WIND_SPAWN
        P[import_indices[n]].unspawned_wind_mass = 0;
#ifdef SINGLE_STAR_FB_JETS
        P[import_indices[n]].unspawned_jet_mass = 0; /* slots are reused, so zeroing at allocation is not enough */
#endif
#ifdef SINGLE_STAR_FB_WINDS
        P[import_indices[n]].wind_mode = 0;         P[import_indices[n]].wind_mode_time = All.Time; /* the mode is persistent state that the spawning hysteresis reads, so a recycled slot must not inherit the previous star's channel */
#endif
#endif
#ifdef SINK_COUNTPROGS
        P[import_indices[n]].Sink_CountProgs = 1;
#endif
        /* record that we actually made a BH, count numbers for book-keeping in domains */
#if (SINK_SEED_FROM_FOF != 1)
        Stars_converted++;
        TimeBinCountGas[P[import_indices[n]].TimeBin]--;
#endif

    }

  All.TotN_gas -= ntot;

  myfree(export_indices);
  myfree(import_indices);
}

#endif // SINK_SEED_FROM_FOF



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

void fof_assign_HostHaloMass(void)	/* assigns mass of host FoF group to CellP[].HostHaloMass for fluid cells */
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
    {  /* rank-local: groups_to_export may carry junk mass payload. Keep the matched Alltoallv below; drain at the poll after datatype cleanup, before any HostHaloMass write. */
      printf("fof_assign_HostHaloMass: start != nimport (start=%d nimport=%lld task=%d)\n", start, (long long)nimport, ThisTask); fflush(stdout);
      endrun(90001016);
    }

  /* send group masses to requesting tasks */
  MPI_Alltoallv(groups_to_export, Recv_count, Recv_offset, mpi_groups_mass_MinID,
		required_groups, Send_count, Send_offset, mpi_groups_mass_MinID, MPI_COMM_WORLD);

  myfree(groups_to_export);
  MPI_Type_free(&mpi_groups_mass_MinID);

  gizmo_exit_bad_stop_if_requested("fof_assign_HostHaloMass:start_nimport");  /* drain the start!=nimport bad-stop after the matched Alltoallv, before any HostHaloMass write */

  qsort(required_groups, NgroupsExt, sizeof(struct group_mass_MinID), compare_group_mass_ID);

    for(i = 0; i < N_gas; i++) {if(P[i].Type==0) {CellP[i].HostHaloMass = 0;}}

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      int bad_group = 0;
      while(FOF_PList[start].MinID < required_groups[i].MinID)
	{
	  start++;
	  if(start >= NumPart)
	    {  /* OOB guard before FOF_PList[start] deref */
	      printf("fof_assign_HostHaloMass: start >= NumPart (start=%d NumPart=%d task=%d)\n", start, NumPart, ThisTask); fflush(stdout);
	      endrun(90001017);
	      bad_group = 1; break;
	    }
	}
      if(bad_group) break;   /* halting: exit the group loop (no collective inside) rather than re-enter the while with start==NumPart (OOB) */

      if(FOF_PList[start].MinID != required_groups[i].MinID)
	{  /* group not found: skip before the per-cell HostHaloMass writes */
	  printf("fof_assign_HostHaloMass: FOF_PList[start].MinID != required (i=%d start=%d task=%d)\n", i, start, ThisTask); fflush(stdout);
	  endrun(90001018);
	  continue;
	}

      for(lenloc = 0; start + lenloc < NumPart;)
	if(FOF_PList[start + lenloc].MinID == required_groups[i].MinID)
	  {
	    if(P[FOF_PList[start + lenloc].Pindex].Type == 0)
	      CellP[FOF_PList[start + lenloc].Pindex].HostHaloMass = required_groups[i].mass;

	    lenloc++;
	  }
	else
	  break;

      start += lenloc;
    }

  gizmo_exit_bad_stop_if_requested("fof_assign_HostHaloMass:group_scan");  /* drain a start>=NumPart / MinID-mismatch bad-stop before HostHaloMass is consumed downstream */

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
  char fname[DEFAULT_PATH_BUFFERSIZE_TOUSE];
  float *mass, *cm;
  int i, j, ntask, *len, count;
  MyIDType *ids;
  int *list_of_ngroups, *list_of_nids, *list_of_allgrouplen;
  int *recvoffset;
  int grnr, ngrp, sendTask, recvTask;
  int nprocgroup, primaryTask, groupTask, nid_previous;
  int fof_compare_P_SubNr(const void *a, const void *b);
    PRINT_STATUS("Trying to read preexisting FoF group catalogues...  (presently allocated=%g MB)",AllocatedBytes / (1024.0 * 1024.0));
  domain_Decomposition(1, 0, 0, 0);

  force_treefree();


  /* start reading of group catalogue */

  if(ThisTask == 0)
    {
      snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, 0);
      if(!(fd = fopen(fname, "r")))
	{
	  printf("can't read file `%s`\n", fname);
	  endrun(11831);	/* rank-0 header read fail: soft bad-stop; still participate in the Bcasts below so peers do not block, then all-rank poll */
	}
      else
	{
	  my_fread(&Ngroups, sizeof(int), 1, fd);
	  my_fread(&TotNgroups, sizeof(int), 1, fd);
	  my_fread(&Nids, sizeof(int), 1, fd);
	  my_fread(&TotNids, sizeof(long long), 1, fd);
	  my_fread(&ntask, sizeof(int), 1, fd);
	  fclose(fd);
	}
    }

  MPI_Bcast(&ntask, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotNgroups, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotNids, sizeof(long long), MPI_BYTE, 0, MPI_COMM_WORLD);
  gizmo_exit_bad_stop_if_requested("fof:read_header");	/* all-rank drain: rank-0 header-read fault propagates here after the Bcasts */

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

      ID_list = (fof_id_list *) mymalloc("ID_list", (TotNids / NTask + NTask) * sizeof(fof_id_list));


      int filenr, target, ngroups, nids, nsend, stored;

      int *ngroup_to_get = (int *) mymalloc("ngroup_to_get", NTask * sizeof(NTask));
      int *nids_to_get = (int *) mymalloc("nids_to_get", NTask * sizeof(NTask));
      int *ngroup_obtained = (int *) mymalloc("ngroup_obtained", NTask * sizeof(NTask));
      int *nids_obtained = (int *) mymalloc("nids_obtained", NTask * sizeof(NTask));

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



      /* === GROUP-TAB distribution phase (rank 0 reads + sends; peers recv) === */
      if(ThisTask == 0)
	{
	  for(filenr = 0; filenr < ntask; filenr++)
	    {
	      snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, filenr);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(11831);		/* rank-0 group_tab read fail: soft bad-stop */
		  for(int t = 1; t < NTask; t++)	/* drain peers still in the GROUP recv loop with an abort sentinel (nsend<0) */
		    if(ngroup_to_get[t] > 0)
		      { int sentinel = -1; MPI_Send(&sentinel, 1, MPI_INT, t, TAG_N, MPI_COMM_WORLD); }
		  break;			/* stop distributing; all ranks meet at the GROUP-phase poll below */
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
	      len = (int *) mymalloc("len", ngroups * sizeof(int));
	      my_fread(len, ngroups, sizeof(int), fd);
	      for(i = 0; i < ngroups; i++)
		tmpGroup[i].Len = len[i];
	      myfree(len);

	      /* offset into id-list */
	      len = (int *) mymalloc("len", ngroups * sizeof(int));
	      my_fread(len, ngroups, sizeof(int), fd);
	      for(i = 0; i < ngroups; i++)
		tmpGroup[i].Offset = len[i];
	      myfree(len);

	      /* mass */
	      mass = (float *) mymalloc("mass", ngroups * sizeof(float));
	      my_fread(mass, ngroups, sizeof(float), fd);
	      for(i = 0; i < ngroups; i++)
		tmpGroup[i].Mass = mass[i];
	      myfree(mass);

	      /* CM */
	      cm = (float *) mymalloc("cm", ngroups * 3 * sizeof(float));
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
	}
      else
	{
	  while(ngroup_to_get[ThisTask])
	    {
	      MPI_Recv(&nsend, 1, MPI_INT, 0, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      if(nsend < 0) { endrun(11831); break; }	/* abort sentinel: rank 0 failed a group_tab read */
	      MPI_Recv(&Group[ngroup_obtained[ThisTask]], nsend * sizeof(group_properties), MPI_BYTE,
		       0, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      ngroup_to_get[ThisTask] -= nsend;
	      ngroup_obtained[ThisTask] += nsend;
	    }
	}
      /* all-rank drain BEFORE the ids phase: never enter the ids receive loop / touch partial catalog on a group_tab fault */
      gizmo_exit_bad_stop_if_requested("fof:readback_group_phase");

      /* === GROUP-IDS distribution phase (rank 0 reads + sends; peers recv) === */
      if(ThisTask == 0)
	{
	  for(filenr = 0; filenr < ntask; filenr++)
	    {
	      snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, filenr);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(1184132);		/* rank-0 group_ids read fail: soft bad-stop */
		  for(int t = 1; t < NTask; t++)	/* drain peers still in the IDS recv loop with an abort sentinel (nsend<0) */
		    if(nids_to_get[t] > 0)
		      { int sentinel = -1; MPI_Send(&sentinel, 1, MPI_INT, t, TAG_HEADER, MPI_COMM_WORLD); }
		  break;			/* stop distributing; all ranks meet at the IDS-phase poll below */
		}
	      my_fread(&ngroups, sizeof(int), 1, fd);
	      my_fread(&TotNgroups, sizeof(int), 1, fd);
	      my_fread(&nids, sizeof(int), 1, fd);
	      my_fread(&TotNids, sizeof(long long), 1, fd);
	      my_fread(&ntask, sizeof(int), 1, fd);
	      my_fread(&nid_previous, sizeof(int), 1, fd);	/* this is the number of IDs in previous files */


	      fof_id_list *tmpID_list = (fof_id_list *) mymalloc("tmpID_list", nids * sizeof(fof_id_list));

	      ids = (MyIDType *) mymalloc("ids", nids * sizeof(MyIDType));

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
			       target, TAG_GASDATA, MPI_COMM_WORLD);
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
	  while(nids_to_get[ThisTask])
	    {
	      MPI_Recv(&nsend, 1, MPI_INT, 0, TAG_HEADER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      if(nsend < 0) { endrun(1184132); break; }	/* abort sentinel: rank 0 failed a group_ids read */
	      MPI_Recv(&ID_list[nids_obtained[ThisTask]], nsend * sizeof(fof_id_list), MPI_BYTE,
		       0, TAG_GASDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      nids_to_get[ThisTask] -= nsend;
	      nids_obtained[ThisTask] += nsend;
	    }
	}
      /* all-rank drain BEFORE any use of the (possibly partial) catalog in the common readback path below */
      gizmo_exit_bad_stop_if_requested("fof:readback_ids_phase");

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

	      snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(11831); goto fof_readback_turn_drain;	/* per-rank group_tab read-open fail: soft bad-stop, skip this rank's reads; per-turn poll drains all ranks */
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
	      len = (int *) mymalloc("len", Ngroups * sizeof(int));
	      my_fread(len, Ngroups, sizeof(int), fd);
	      for(i = 0; i < Ngroups; i++)
		Group[i].Len = len[i];
	      myfree(len);

	      /* offset into id-list */
	      len = (int *) mymalloc("len", Ngroups * sizeof(int));
	      my_fread(len, Ngroups, sizeof(int), fd);
	      for(i = 0; i < Ngroups; i++)
		Group[i].Offset = len[i];
	      myfree(len);

	      /* mass */
	      mass = (float *) mymalloc("mass", Ngroups * sizeof(float));
	      my_fread(mass, Ngroups, sizeof(float), fd);
	      for(i = 0; i < Ngroups; i++)
		Group[i].Mass = mass[i];
	      myfree(mass);

	      /* CM */
	      cm = (float *) mymalloc("cm", Ngroups * 3 * sizeof(float));
	      my_fread(cm, Ngroups, 3 * sizeof(float), fd);
	      for(i = 0; i < Ngroups; i++)
		for(j = 0; j < 3; j++)
		  Group[i].CM[j] = cm[i * 3 + j];
	      myfree(cm);

	      fclose(fd);

	      snprintf(fname, DEFAULT_PATH_BUFFERSIZE_TOUSE, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, ThisTask);
	      if(!(fd = fopen(fname, "r")))
		{
		  printf("can't read file `%s`\n", fname);
		  endrun(1184132); goto fof_readback_turn_drain;	/* per-rank group_ids read-open fail: soft bad-stop, skip; per-turn poll drains all ranks */
		}
	      my_fread(&Ngroups, sizeof(int), 1, fd);
	      my_fread(&TotNgroups, sizeof(int), 1, fd);
	      my_fread(&Nids, sizeof(int), 1, fd);
	      my_fread(&TotNids, sizeof(long long), 1, fd);
	      my_fread(&ntask, sizeof(int), 1, fd);
	      my_fread(&nid_previous, sizeof(int), 1, fd);	/* this is the number of IDs in previous files */

	      ID_list = (fof_id_list *) mymalloc("ID_list", Nids * sizeof(fof_id_list));
	      ids = (MyIDType *) mymalloc("ids", Nids * sizeof(MyIDType));

	      my_fread(ids, sizeof(MyIDType), Nids, fd);

	      for(i = 0; i < Nids; i++)
		ID_list[i].ID = ids[i];

	      myfree(ids);

	      fclose(fd);
	    }

	fof_readback_turn_drain:
	  gizmo_exit_bad_stop_if_requested("fof:readback_turn");	/* all-rank drain (replaces in-group barrier): catches per-rank read-IO faults */
	}
    }

  t1 = my_second();
  if(ThisTask == 0)
    printf("reading  took %g sec\n", timediff(t0, t1));



  t0 = my_second();

  /* now need to assign group numbers */


  list_of_ngroups = (int *) mymalloc("list_of_ngroups", NTask * sizeof(int));
  list_of_nids = (int *) mymalloc("list_of_nids", NTask * sizeof(int));

  MPI_Allgather(&Ngroups, 1, MPI_INT, list_of_ngroups, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&Nids, 1, MPI_INT, list_of_nids, 1, MPI_INT, MPI_COMM_WORLD);

  list_of_allgrouplen = (int *) mymalloc("list_of_allgrouplen", TotNgroups * sizeof(int));

  recvoffset = (int *) mymalloc("recvoffset", NTask * sizeof(int));
  for(i = 1, recvoffset[0] = 0; i < NTask; i++)
    recvoffset[i] = recvoffset[i - 1] + list_of_ngroups[i - 1];
  len = (int *) mymalloc("len", Ngroups * sizeof(int));
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
    { endrun(8881); gizmo_exit_bad_stop_if_requested("fof:read_totlen_mismatch"); }	/* symmetric (totlen from the Allgathered list_of_allgrouplen, TotNids global): all ranks reach the poll together */


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


  long long *idoffset = (long long *) mymalloc("idoffset", NTask * sizeof(long long));

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
    endrun(1231);	/* per-rank ID-assignment count off: soft bad-stop; bounded (no OOB), drains at the readback tail poll below */

  myfree(idoffset);

  t1 = my_second();
  if(ThisTask == 0)
    printf("assigning took %g sec\n", timediff(t0, t1));



  t0 = my_second();

  for(i = 0; i < NumPart; i++)
    P[i].SubNr = i;

  qsort(P, NumPart, sizeof(struct particle_data), io_compare_P_ID);
  wakeup_sidecar_invalidate();   /* P[] reordered in place → rebuild WakeupDirty from P[] next scan */
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
		  recv_ID_list = (fof_id_list *) mymalloc("recv_ID_list", list_of_nids[recvTask] * sizeof(fof_id_list));

		  /* get the particles */
		  MPI_Sendrecv(ID_list, Nids * sizeof(fof_id_list), MPI_BYTE, recvTask, TAG_FOF_H,
			       recv_ID_list, list_of_nids[recvTask] * sizeof(fof_id_list), MPI_BYTE,
			       recvTask, TAG_FOF_H, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

	      for(i = 0, j = 0; i < list_of_nids[recvTask]; i++)
		{
		  /* Matching is by identifier alone, so FOF/SUBFIND REQUIRES UNIQUE PARTICLE IDs.
		     Where IDs repeat this scan stops at the first particle carrying the identifier and
		     hands it the group number, so duplicates are assigned arbitrarily. Carrying the child
		     number and generation here would fix the match but not the problem: the group_ids
		     catalog on disk stores bare IDs, so a full-identity key cannot survive a write/read
		     round trip without a versioned format change. Configurations that stamp duplicate IDs
		     (SINK_WIND_SPAWN) are therefore not supported with group catalogs. */
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
    { endrun(543); gizmo_exit_bad_stop_if_requested("fof:read_matches_mismatch"); }	/* symmetric (totlen is the post-sumup_large_ints global match count, TotNids global): all ranks reach the poll together */

  t1 = my_second();
  if(ThisTask == 0)
    printf("assigning GrNr to P[] took %g sec\n", timediff(t0, t1));

  gizmo_exit_bad_stop_if_requested("fof:readback_tail");	/* all-rank drain (replaces tail barrier): catches the per-rank count!=Nids fault above */

  myfree(list_of_allgrouplen);
  myfree(list_of_nids);
  myfree(list_of_ngroups);

  /* restore peano-hilbert order */
  qsort(P, NumPart, sizeof(struct particle_data), fof_compare_P_SubNr);
  wakeup_sidecar_invalidate();   /* P[] reordered in place → rebuild WakeupDirty from P[] next scan */
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
