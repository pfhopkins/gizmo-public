#include <algorithm>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <vector>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/sfc_tiles.h"
#include "group_search.h"
#include "../mesh/gpu_neighbor_list.h"

struct group_search_tile_meta_t
{
  double lo[3], hi[3];
  int first, count;
};

static double group_search_aabb_distance2(const double alo[3], const double ahi[3], const double blo[3], const double bhi[3],
                                          double max_gap)
{
  double dist2 = 0;
  for(int k = 0; k < 3; k++)
    {
      int is_periodic = (k == 0) ? TILE_PERIODIC_X : ((k == 1) ? TILE_PERIODIC_Y : TILE_PERIODIC_Z);
      double bsize = (k == 0) ? boxSize_X : ((k == 1) ? boxSize_Y : boxSize_Z);
      double ca = 0.5 * (alo[k] + ahi[k]);
      double cb = 0.5 * (blo[k] + bhi[k]);
      double ha = 0.5 * (ahi[k] - alo[k]);
      double hb = 0.5 * (bhi[k] - blo[k]);
      double dx = fabs(ca - cb);
      if(is_periodic && dx > 0.5 * bsize) dx = bsize - dx;
      double gap = dx - ha - hb;
      if(gap <= 0) continue;
      if(gap > max_gap) return max_gap * max_gap + 1;
      dist2 += gap * gap;
    }
  return dist2;
}

static void group_search_build_tiles_from_indices(const std::vector<int> &indices, std::vector<group_search_tile_meta_t> &tiles)
{
  tiles.clear();
  if(indices.empty()) return;
  int ntiles = (indices.size() + TILE_TARGET_SIZE - 1) / TILE_TARGET_SIZE;
  tiles.resize(ntiles);
  for(int t = 0; t < ntiles; t++)
    {
      int first = t * TILE_TARGET_SIZE;
      int count = TILE_TARGET_SIZE;
      if(first + count > (int)indices.size()) count = indices.size() - first;
      group_search_tile_meta_t &tile = tiles[t];
      tile.first = first;
      tile.count = count;
      int j0 = indices[first];
      for(int k = 0; k < 3; k++) tile.lo[k] = tile.hi[k] = P[j0].Pos[k];
      for(int s = 1; s < count; s++)
        {
          int j = indices[first + s];
          for(int k = 0; k < 3; k++)
            {
              if(P[j].Pos[k] < tile.lo[k]) tile.lo[k] = P[j].Pos[k];
              if(P[j].Pos[k] > tile.hi[k]) tile.hi[k] = P[j].Pos[k];
            }
        }
    }
}

static void group_search_build_tiles_from_points(const std::vector<double> &centers, const std::vector<double> &radii,
                                                 std::vector<group_search_tile_meta_t> &tiles)
{
  tiles.clear();
  int npoints = radii.size();
  if(npoints <= 0) return;
  int ntiles = (npoints + TILE_TARGET_SIZE - 1) / TILE_TARGET_SIZE;
  tiles.resize(ntiles);
  for(int t = 0; t < ntiles; t++)
    {
      int first = t * TILE_TARGET_SIZE;
      int count = TILE_TARGET_SIZE;
      if(first + count > npoints) count = npoints - first;
      group_search_tile_meta_t &tile = tiles[t];
      tile.first = first;
      tile.count = count;
      for(int k = 0; k < 3; k++)
        {
          tile.lo[k] = centers[3 * first + k] - radii[first];
          tile.hi[k] = centers[3 * first + k] + radii[first];
        }
      for(int s = 1; s < count; s++)
        {
          int p = first + s;
          for(int k = 0; k < 3; k++)
            {
              double lo = centers[3 * p + k] - radii[p];
              double hi = centers[3 * p + k] + radii[p];
              if(lo < tile.lo[k]) tile.lo[k] = lo;
              if(hi > tile.hi[k]) tile.hi[k] = hi;
            }
        }
    }
}

static void group_search_allgather_tiles(const std::vector<group_search_tile_meta_t> &local_tiles,
                                         std::vector<group_search_tile_meta_t> &all_tiles,
                                         std::vector<int> &tile_counts, std::vector<int> &tile_offsets)
{
  tile_counts.assign(NTask, 0);
  tile_offsets.assign(NTask, 0);
  int local_count = local_tiles.size();
  MPI_Allgather(&local_count, 1, MPI_INT, tile_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  int total = 0;
  for(int task = 0; task < NTask; task++)
    {
      tile_offsets[task] = total;
      total += tile_counts[task];
    }
  all_tiles.resize(total);
  std::vector<int> byte_counts(NTask), byte_offsets(NTask);
  for(int task = 0; task < NTask; task++)
    {
      byte_counts[task] = tile_counts[task] * (int)sizeof(group_search_tile_meta_t);
      byte_offsets[task] = tile_offsets[task] * (int)sizeof(group_search_tile_meta_t);
    }
  MPI_Allgatherv((void *)(local_tiles.empty() ? NULL : local_tiles.data()), local_count * (int)sizeof(group_search_tile_meta_t), MPI_BYTE,
                 (void *)(all_tiles.empty() ? NULL : all_tiles.data()), byte_counts.data(), byte_offsets.data(), MPI_BYTE, MPI_COMM_WORLD);
}

static int group_search_particle_matches_group(const int i, const int group_filter)
{
  if(group_filter < 0) return 1;
#ifdef SUBFIND
  return P[i].GrNr == group_filter;
#else
  if(ThisTask == 0)
    printf("group_search_import_particles received group_filter=%d but SUBFIND is not enabled.\n", group_filter);
  fflush(stdout);
  endrun(990511);
  return 0;
#endif
}

void group_search_import_particles(int source_type_mask, int import_type_mask, double search_radius,
                                   const MyIDType *label_id, const MyIDType *label_task,
                                   group_search_import_pool_t &pool, int group_filter, bool import_cells)
{
  pool.particles.assign(P, P + NumPart);
  pool.cells.clear();
  if(import_cells) pool.cells.assign(CellP, CellP + NumPart);
  pool.ghost_home_task.clear();
  pool.ghost_home_index.clear();
  pool.ghost_label_id.clear();
  pool.ghost_label_task.clear();

  if(NTask <= 1) return;
  double global_search_radius = search_radius;
  MPI_Allreduce(&search_radius, &global_search_radius, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(global_search_radius <= 0) return;

  std::vector<int> source_indices, import_indices;
  source_indices.reserve(NumPart);
  import_indices.reserve(NumPart);
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Mass <= 0) continue;
      if(!group_search_particle_matches_group(i, group_filter)) continue;
      if(((1 << P[i].Type) & source_type_mask)) source_indices.push_back(i);
      if(((1 << P[i].Type) & import_type_mask)) import_indices.push_back(i);
    }
  std::vector<group_search_tile_meta_t> local_source_tiles, local_import_tiles;
  std::vector<group_search_tile_meta_t> all_source_tiles, all_import_tiles;
  std::vector<int> source_counts, source_offsets, import_counts, import_offsets;
  group_search_build_tiles_from_indices(source_indices, local_source_tiles);
  group_search_build_tiles_from_indices(import_indices, local_import_tiles);
  group_search_allgather_tiles(local_source_tiles, all_source_tiles, source_counts, source_offsets);
  group_search_allgather_tiles(local_import_tiles, all_import_tiles, import_counts, import_offsets);
  if(all_import_tiles.empty()) return;

  std::vector<unsigned char> send_tile(local_import_tiles.size() * NTask, 0);
  double r2 = global_search_radius * global_search_radius;
  for(int task = 0; task < NTask; task++)
    {
      if(task == ThisTask || source_counts[task] == 0) continue;
      int s0 = source_offsets[task];
      for(size_t lt = 0; lt < local_import_tiles.size(); lt++)
        {
          for(int st = 0; st < source_counts[task]; st++)
            {
              if(group_search_aabb_distance2(local_import_tiles[lt].lo, local_import_tiles[lt].hi,
                                             all_source_tiles[s0 + st].lo, all_source_tiles[s0 + st].hi, global_search_radius) < r2)
                {
                  send_tile[lt * NTask + task] = 1;
                  break;
                }
            }
        }
    }

  std::vector<unsigned char> need_tile(all_import_tiles.size(), 0);
  for(int task = 0; task < NTask; task++)
    {
      if(task == ThisTask || import_counts[task] == 0) continue;
      int p0 = import_offsets[task];
      for(size_t ls = 0; ls < local_source_tiles.size(); ls++)
        {
          for(int pt = 0; pt < import_counts[task]; pt++)
            {
              if(group_search_aabb_distance2(local_source_tiles[ls].lo, local_source_tiles[ls].hi,
                                             all_import_tiles[p0 + pt].lo, all_import_tiles[p0 + pt].hi, global_search_radius) < r2)
                {
                  need_tile[p0 + pt] = 1;
                }
            }
        }
    }

  std::vector<int> send_count(NTask, 0), recv_count(NTask, 0), send_disp(NTask, 0), recv_disp(NTask, 0);
  for(int task = 0; task < NTask; task++)
    {
      if(task == ThisTask) continue;
      for(size_t lt = 0; lt < local_import_tiles.size(); lt++)
        if(send_tile[lt * NTask + task])
          send_count[task] += local_import_tiles[lt].count;
      for(int rt = 0; rt < import_counts[task]; rt++)
        if(need_tile[import_offsets[task] + rt])
          recv_count[task] += all_import_tiles[import_offsets[task] + rt].count;
    }
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

  std::vector<particle_data> send_particles(total_send > 0 ? total_send : 1);
  std::vector<gas_cell_data> send_cells(import_cells ? (total_send > 0 ? total_send : 1) : 0);
  std::vector<int> send_home_index(total_send > 0 ? total_send : 1);
  std::vector<MyIDType> send_label_id(total_send > 0 ? total_send : 1);
  std::vector<MyIDType> send_label_task(total_send > 0 ? total_send : 1);
  std::vector<int> cursor = send_disp;
  for(int task = 0; task < NTask; task++)
    {
      if(task == ThisTask || send_count[task] == 0) continue;
      for(size_t lt = 0; lt < local_import_tiles.size(); lt++)
        {
          if(!send_tile[lt * NTask + task]) continue;
          const group_search_tile_meta_t &tile = local_import_tiles[lt];
          for(int s = 0; s < tile.count; s++)
            {
              int idx = import_indices[tile.first + s];
              int off = cursor[task]++;
              send_particles[off] = P[idx];
              if(import_cells) send_cells[off] = CellP[idx];
              send_home_index[off] = idx;
              send_label_id[off] = label_id ? label_id[idx] : P[idx].ID;
              send_label_task[off] = label_task ? label_task[idx] : ThisTask;
            }
        }
    }

  std::vector<int> send_bytes(NTask), recv_bytes(NTask), send_bdisp(NTask), recv_bdisp(NTask);
  for(int task = 0; task < NTask; task++)
    {
      send_bytes[task] = send_count[task] * (int)sizeof(particle_data);
      recv_bytes[task] = recv_count[task] * (int)sizeof(particle_data);
      send_bdisp[task] = send_disp[task] * (int)sizeof(particle_data);
      recv_bdisp[task] = recv_disp[task] * (int)sizeof(particle_data);
    }
  size_t old_size = pool.particles.size();
  pool.particles.resize(old_size + total_recv);
  MPI_Alltoallv(send_particles.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                total_recv > 0 ? &pool.particles[old_size] : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
  if(import_cells)
    {
      for(int task = 0; task < NTask; task++)
        {
          send_bytes[task] = send_count[task] * (int)sizeof(gas_cell_data);
          recv_bytes[task] = recv_count[task] * (int)sizeof(gas_cell_data);
          send_bdisp[task] = send_disp[task] * (int)sizeof(gas_cell_data);
          recv_bdisp[task] = recv_disp[task] * (int)sizeof(gas_cell_data);
        }
      pool.cells.resize(old_size + total_recv);
      MPI_Alltoallv(send_cells.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                    total_recv > 0 ? &pool.cells[old_size] : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
    }

  pool.ghost_home_index.resize(total_recv > 0 ? total_recv : 0);
  pool.ghost_home_task.resize(total_recv > 0 ? total_recv : 0);
  for(int task = 0; task < NTask; task++)
    {
      send_bytes[task] = send_count[task] * (int)sizeof(int);
      recv_bytes[task] = recv_count[task] * (int)sizeof(int);
      send_bdisp[task] = send_disp[task] * (int)sizeof(int);
      recv_bdisp[task] = recv_disp[task] * (int)sizeof(int);
    }
  MPI_Alltoallv(send_home_index.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                total_recv > 0 ? pool.ghost_home_index.data() : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
  for(int task = 0; task < NTask; task++)
    for(int g = 0; g < recv_count[task]; g++)
      pool.ghost_home_task[recv_disp[task] + g] = task;

  pool.ghost_label_id.resize(total_recv > 0 ? total_recv : 0);
  pool.ghost_label_task.resize(total_recv > 0 ? total_recv : 0);
  for(int task = 0; task < NTask; task++)
    {
      send_bytes[task] = send_count[task] * (int)sizeof(MyIDType);
      recv_bytes[task] = recv_count[task] * (int)sizeof(MyIDType);
      send_bdisp[task] = send_disp[task] * (int)sizeof(MyIDType);
      recv_bdisp[task] = recv_disp[task] * (int)sizeof(MyIDType);
    }
  MPI_Alltoallv(send_label_id.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                total_recv > 0 ? pool.ghost_label_id.data() : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
  MPI_Alltoallv(send_label_task.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                total_recv > 0 ? pool.ghost_label_task.data() : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
}

void group_search_import_particles_around_points(int import_type_mask, const std::vector<double> &centers,
                                                 const std::vector<double> &radii,
                                                 group_search_import_pool_t &pool, bool import_cells)
{
  pool.particles.assign(P, P + NumPart);
  pool.cells.clear();
  if(import_cells) pool.cells.assign(CellP, CellP + NumPart);
  pool.ghost_home_task.clear();
  pool.ghost_home_index.clear();
  pool.ghost_label_id.clear();
  pool.ghost_label_task.clear();

  int npoints = radii.size();
  if((int)centers.size() != 3 * npoints)
    {
      printf("group_search_import_particles_around_points received inconsistent center/radius arrays: centers=%d radii=%d task=%d\n",
             (int)centers.size(), npoints, ThisTask);
      fflush(stdout);
      endrun(990512);
    }
  if(NTask <= 1) return;

  std::vector<int> import_indices;
  import_indices.reserve(NumPart);
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Mass <= 0) continue;
      if(((1 << P[i].Type) & import_type_mask)) import_indices.push_back(i);
    }

  std::vector<group_search_tile_meta_t> local_query_tiles, local_import_tiles, all_query_tiles;
  std::vector<int> query_counts, query_offsets;
  group_search_build_tiles_from_points(centers, radii, local_query_tiles);
  group_search_build_tiles_from_indices(import_indices, local_import_tiles);
  group_search_allgather_tiles(local_query_tiles, all_query_tiles, query_counts, query_offsets);
  if(all_query_tiles.empty()) return;

  std::vector<unsigned char> send_tile(local_import_tiles.size() * NTask, 0);
  for(int task = 0; task < NTask; task++)
    {
      if(task == ThisTask || query_counts[task] == 0) continue;
      int q0 = query_offsets[task];
      for(size_t lt = 0; lt < local_import_tiles.size(); lt++)
        {
          for(int qt = 0; qt < query_counts[task]; qt++)
            {
              if(group_search_aabb_distance2(local_import_tiles[lt].lo, local_import_tiles[lt].hi,
                                             all_query_tiles[q0 + qt].lo, all_query_tiles[q0 + qt].hi, 0.0) == 0)
                {
                  send_tile[lt * NTask + task] = 1;
                  break;
                }
            }
        }
    }

  std::vector<int> send_count(NTask, 0), recv_count(NTask, 0), send_disp(NTask, 0), recv_disp(NTask, 0);
  for(int task = 0; task < NTask; task++)
    {
      if(task == ThisTask) continue;
      for(size_t lt = 0; lt < local_import_tiles.size(); lt++)
        if(send_tile[lt * NTask + task])
          send_count[task] += local_import_tiles[lt].count;
    }
  MPI_Alltoall(send_count.data(), 1, MPI_INT, recv_count.data(), 1, MPI_INT, MPI_COMM_WORLD);

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

  std::vector<particle_data> send_particles(total_send > 0 ? total_send : 1);
  std::vector<gas_cell_data> send_cells(import_cells ? (total_send > 0 ? total_send : 1) : 0);
  std::vector<int> send_home_index(total_send > 0 ? total_send : 1);
  std::vector<int> cursor = send_disp;
  for(int task = 0; task < NTask; task++)
    {
      if(task == ThisTask || send_count[task] == 0) continue;
      for(size_t lt = 0; lt < local_import_tiles.size(); lt++)
        {
          if(!send_tile[lt * NTask + task]) continue;
          const group_search_tile_meta_t &tile = local_import_tiles[lt];
          for(int s = 0; s < tile.count; s++)
            {
              int idx = import_indices[tile.first + s];
              int off = cursor[task]++;
              send_particles[off] = P[idx];
              if(import_cells) send_cells[off] = CellP[idx];
              send_home_index[off] = idx;
            }
        }
    }

  std::vector<int> send_bytes(NTask), recv_bytes(NTask), send_bdisp(NTask), recv_bdisp(NTask);
  for(int task = 0; task < NTask; task++)
    {
      send_bytes[task] = send_count[task] * (int)sizeof(particle_data);
      recv_bytes[task] = recv_count[task] * (int)sizeof(particle_data);
      send_bdisp[task] = send_disp[task] * (int)sizeof(particle_data);
      recv_bdisp[task] = recv_disp[task] * (int)sizeof(particle_data);
    }
  size_t old_size = pool.particles.size();
  pool.particles.resize(old_size + total_recv);
  MPI_Alltoallv(send_particles.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                total_recv > 0 ? &pool.particles[old_size] : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);

  if(import_cells)
    {
      for(int task = 0; task < NTask; task++)
        {
          send_bytes[task] = send_count[task] * (int)sizeof(gas_cell_data);
          recv_bytes[task] = recv_count[task] * (int)sizeof(gas_cell_data);
          send_bdisp[task] = send_disp[task] * (int)sizeof(gas_cell_data);
          recv_bdisp[task] = recv_disp[task] * (int)sizeof(gas_cell_data);
        }
      pool.cells.resize(old_size + total_recv);
      MPI_Alltoallv(send_cells.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                    total_recv > 0 ? &pool.cells[old_size] : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
    }

  pool.ghost_home_index.resize(total_recv > 0 ? total_recv : 0);
  pool.ghost_home_task.resize(total_recv > 0 ? total_recv : 0);
  for(int task = 0; task < NTask; task++)
    {
      send_bytes[task] = send_count[task] * (int)sizeof(int);
      recv_bytes[task] = recv_count[task] * (int)sizeof(int);
      send_bdisp[task] = send_disp[task] * (int)sizeof(int);
      recv_bdisp[task] = recv_disp[task] * (int)sizeof(int);
    }
  MPI_Alltoallv(send_home_index.data(), send_bytes.data(), send_bdisp.data(), MPI_BYTE,
                total_recv > 0 ? pool.ghost_home_index.data() : NULL, recv_bytes.data(), recv_bdisp.data(), MPI_BYTE, MPI_COMM_WORLD);
  for(int task = 0; task < NTask; task++)
    for(int g = 0; g < recv_count[task]; g++)
      pool.ghost_home_task[recv_disp[task] + g] = task;
}

void group_search_build_cross_type_nl(std::vector<particle_data> &particles, int *sources, int nsources,
                                      double *radii, int j_type_mask, neighbor_list_t *nl)
{
  int npool = 0;
  for(size_t p = 0; p < particles.size(); p++)
    if(particles[p].Mass > 0 && ((1 << particles[p].Type) & j_type_mask))
      npool++;
  if(nsources <= 0 || npool <= 0)
    {
      nl->num_active = nsources;
      nl->total_pairs = 0;
      nl->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(nsources + 1) * sizeof(int64_t));
      nl->neighbors = (int *) mymalloc("ngb_neighbors", sizeof(int));
      for(int a = 0; a <= nsources; a++) nl->offsets[a] = 0;
      return;
    }
  gpu_build_cross_type_neighbor_list(particles.data(), particles.size(), sources, nsources, radii,
                                     j_type_mask, NGB_SEARCH_ONEWAY, nl);
  std::vector<int64_t> offsets(nsources + 1, 0);
  std::vector<int> neighbors;
  neighbors.reserve((size_t)nl->total_pairs);
  int invalid_entries = 0;
  for(int a = 0; a < nsources; a++)
    {
      int i = sources[a];
      if(i < 0 || i >= (int)particles.size())
        {
          printf("group_search_build_cross_type_nl received invalid source index=%d n_particles=%d row=%d task=%d\n",
                 i, (int)particles.size(), a, ThisTask);
          fflush(stdout);
          endrun(990510);
        }
      double h2 = radii[a] * radii[a];
      for(int64_t n = nl->offsets[a]; n < nl->offsets[a + 1]; n++)
        {
          int j = nl->neighbors[n];
          if(j < 0 || j >= (int)particles.size())
            {
              invalid_entries++;
              continue;
            }
          const particle_data &pj = particles[j];
          if(pj.Mass <= 0) continue;
          if(!((1 << pj.Type) & j_type_mask)) continue;
          if(group_search_distance2_particles(particles[i], pj) >= h2) continue;
          neighbors.push_back(j);
        }
      offsets[a + 1] = neighbors.size();
    }
  if(invalid_entries > 0 && ThisTask == 0)
    printf("group_search_build_cross_type_nl pruned %d invalid raw neighbor-list entries on task=%d\n", invalid_entries, ThisTask);
  free_neighbor_list(nl);
  nl->num_active = nsources;
  nl->total_pairs = (int64_t)neighbors.size();
  nl->offsets = (int64_t *) mymalloc("ngb_offsets", (size_t)(nsources + 1) * sizeof(int64_t));
  nl->neighbors = (int *) mymalloc("ngb_neighbors", (neighbors.size() > 0 ? neighbors.size() : 1) * sizeof(int));
  for(int a = 0; a <= nsources; a++) nl->offsets[a] = offsets[a];
  for(size_t n = 0; n < neighbors.size(); n++) nl->neighbors[n] = neighbors[n];
}

double group_search_distance2_particles(const particle_data &a, const particle_data &b)
{
  Vec3<double> dr = {a.Pos[0] - b.Pos[0], a.Pos[1] - b.Pos[1], a.Pos[2] - b.Pos[2]};
  nearest_xyz(dr, 1);
  return dr.norm_sq();
}
