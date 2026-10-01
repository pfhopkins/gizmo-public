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

static long long Count[BINS_TP];
static long long CountSpheres[BINS_TP];
static double Xi[BINS_TP];
static double Rbin[BINS_TP];

static double R0, R1;		/* inner and outer radius for correlation function determination */

static double logR0;
static double binfac;



/*  This function computes the two-point function.
 */
void twopoint(void)
{
    int i, j, bin, n;
    double p, rs, vol, scaled_frac, tstart, tend; long long *countbuf;
    PRINT_STATUS("begin two-point correlation function..."); tstart = my_second();
    /* set inner and outer radius for the bins that are used for the correlation function estimate */
    R0 = All.ForceSoftening[1] / 3.; R1 = All.BoxSize / 2; /* we assume that type=1 is the primary type */
    scaled_frac = FRACTION_TP * 1.0e7 / All.TotNumPart; logR0 = log(R0); binfac = BINS_TP / (log(R1) - log(R0));
    for(i = 0; i < BINS_TP; i++) {Count[i] = 0; CountSpheres[i] = 0;}

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
     * PROPER PORT (not yet done): this routine should NOT use the normal
     * symmetric NL infrastructure at all. Two-point correlation is
     * fundamentally an "all-against-all long-range" problem with binned
     * distance accumulators — exactly what the gravity tree already
     * supports. The proper port walks the modern GRAVITY TREE (Nodes_base /
     * SoA in gravity/forcetree.cc) with a custom op that, instead of
     * accumulating gravitational potential via mass moments, accumulates
     * BINNED PAIR COUNT via node->mass for nodes that fall cleanly within a
     * single radial bin. The walk pattern is identical to gravtree; only the
     * accumulated quantity differs. That port preserves the O(N log N)
     * scaling at large radii.
     *
     * The interim NL-based implementation is acceptable because
     * OUTPUT_TWOPOINT_ENABLED runs ~1x per 1e4-1e5 timesteps, so the O(N)
     * per-pair cost is tolerable. The gravity-tree port is needed before any
     * production run with TotNumPart > ~1e8 enables this module.
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
            gpu_ngb_list_free(&gnl);
            gpu_particles_arena_invalidate();
        }
        if(ghost_imported) ghost_exchange_cleanup();
    }
    random_generator = saved_rng;

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
      if(!(fd = fopen(buf, "w"))) {printf("can't open file `%s`\n", buf); endrun(1323); return;}	/* rank-0 only, no collective after: skip writes on open failure; drains downstream */
      fprintf(fd, "%.16g\n", All.Time); i = BINS_TP; fprintf(fd, "%d\n", i);
      for(i = 0; i < BINS_TP; i++) {fprintf(fd, "%g %g %g %g\n", Rbin[i], Xi[i], (double) Count[i], (double) CountSpheres[i]);}
      fclose(fd);
    }
}
#endif
