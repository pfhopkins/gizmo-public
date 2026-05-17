#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "../mesh/mesh_motion.h"
#include "../mesh/neighbor_list.h"
#include "../mesh/sfc_tiles.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "../core/step_phases.h"
#include "../system/gpu_particles_arena.h"
#if defined(HYDRO_VOLUME_CORRECTIONS)
#include <vector>
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/ghost_writeback.h"
#endif

/*! \file density.c
 *  \brief hydro kernel size and neighbor determination, volumetric quantities calculated
 *
 *  This file contains the "first hydro loop", where the gas densities and some
 *  auxiliary quantities are computed.  There is also functionality that corrects the kernel length if needed.
 */
/*!
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * The code has been modified substantially (condensed, different criteria for kernel lengths, optimizatins,
 * rewritten parallelism, new physics included, new variable/memory conventions added, fundamentally different
 * criteria and conditioning and calcuilations actually being done for the modular hydro solvers)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/* density() and density_isactive() live in hydro/density_loop.cc. The
 * cellcorrections routines remain here as peer step-phase entries. */


/* Routines for a loop after the iterative density loop needed to find neighbors, etc, once all have converged, to apply additional correction terms to the cell volumes and faces (for those needed -before- the gradients loop because they alter primitive quantities needed for gradients, such as particle densities, pressures, etc.)
    This was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */
#ifdef HYDRO_VOLUME_CORRECTIONS

/* The old CPU exchange scaffolding for cell corrections has been retired.
 * The modern path walks the prebuilt symmetric CSR neighbor list directly. */

/* final operations for after the updates are computed */
void cellcorrections_final_operations_and_cleanup(void)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int _apl = 0; _apl < (int)ActiveParticleList.size(); _apl++) { int i = ActiveParticleList[_apl]; /* check all active elements */
        if(GasGrad_isactive(i, P, CellP)) /* only cells eligible for gradients and hydro */
        {
            if(CellP[i].Volume_1 > 0) {CellP[i].Density = P[i].Mass / CellP[i].Volume_1;} else {CellP[i].Volume_1 = CellP[i].Volume_0;} // set the updated density. other variables that need volumes will all scale off this, so we can rely on it to inform everything else [if bad value here, revert to the 0th-order volume quadrature]
            set_eos_pressure(i, P, CellP);
        }}
}

/* parent routine which calls the work loop above */
void cellcorrections_calc(void)
{
    CPU_Step[CPU_DENSMISC] += measure_time(); double t00_truestart = my_second();
    double timeall = 0, timecomp = 0, timewait = 0, timecomm = 0;
    PRINT_STATUS(" ..calculating first-order corrections to cell sizes/faces");
    /* Modern path: prebuilt symmetric CSR NL. Walks neighbors per active gas
     * cell and accumulates Volume_1 += V_j^2 * wk(r, h_j). Symmetric search
     * (NGB_SEARCH_SYMMETRIC) ensures r < max(h_i, h_j), then per-pair filter
     * r < h_j matches legacy semantic that uses j's kernel for weighting.
     * Ghost cells are valid neighbors (their Volume_0 set during density);
     * j-side write is not done here (i-only accumulation). */
    double t_kern_start = my_second();
    {
        std::vector<int> active_idx;
        std::vector<double> radii;
        active_idx.reserve(N_gas);
        radii.reserve(N_gas);
        for (int aa = 0; aa < (int)ActiveParticleList.size(); aa++) {
            int i = ActiveParticleList[aa];
            if (P[i].Type != 0 || P[i].Mass <= 0) continue;
            if (!GasGrad_isactive(i, P, CellP)) continue;
            if (P[i].KernelRadius <= 0) continue;
            active_idx.push_back(i);
            radii.push_back(P[i].KernelRadius);
        }

        int num_src = (int)active_idx.size();
        int num_src_global = 0;
        MPI_Allreduce(&num_src, &num_src_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        gpu_neighbor_list_t gnl = {};
        std::vector<int> gnl_neighbors_host;
        int imported_ghosts = 0;
        if (num_src_global > 0) {
            /* Defensive ghost prep: legacy used code_block_xchange MPI export to
             * pull in cross-rank j-neighbors. Modern path replaces that with
             * symmetric ghost particles. cellcorrections_calc is called between
             * density and gradients so ghosts are typically already alive — but
             * if a future caller invokes this with no ghosts (e.g. standalone
             * diagnostic), import them here so cross-rank contributions to
             * Volume_1 are NOT silently dropped. */
            int need_import_local = (ghost_get_num_ghosts() <= 0) ? 1 : 0;
            int need_import = 0;
            MPI_Allreduce(&need_import_local, &need_import, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (need_import) {
                if (ghost_get_num_ghosts() > 0) ghost_exchange_cleanup();
                gizmo_density_prep_ghosts(gizmo_ghost_safety_factor());
                imported_ghosts = 1;
            }
            int local_count = ghost_get_num_local();
            if (local_count <= 0) local_count = NumPart;
            int num_all = local_count + ghost_get_num_ghosts();
            if (num_all <= 0) num_all = NumPart;
            if (num_src > 0) {
                gpu_particles_arena_acquire(num_all, P, CellP);
                struct particle_data *P_gpu = gpu_particles_arena_P();
                gpu_ngb_list_build(P_gpu, num_all,
                                   active_idx.data(), num_src,
                                   NGB_SEARCH_SYMMETRIC, 1 /* gas only */,
                                   &gnl, NULL, 1.0, radii.data(), NULL, "dens-vol1");
                /* gnl.neighbors is DEVICE_SPACE; host loop below indexes it. */
                if (gnl.total_pairs > 0) {
                    gnl_neighbors_host.resize((size_t)gnl.total_pairs);
                    gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
                }
            }
        }
        const int *gnl_neighbors = gnl_neighbors_host.empty() ? NULL : gnl_neighbors_host.data();

        for (int aa = 0; aa < num_src; aa++) {
            int i = active_idx[aa];
            Vec3<MyDouble> pos_i = P[i].Pos;
            int64_t n_off = gnl.offsets[aa], n_off_end = gnl.offsets[aa+1];
            double accum_V1 = 0;
            for (int64_t nn = n_off; nn < n_off_end; nn++) {
                int j = gnl_neighbors[nn];
                Vec3<double> dp = pos_i - P[j].Pos;
                nearest_xyz(dp);
                double r2 = dp.norm_sq();
                double h_j = P[j].KernelRadius;
                if (r2 >= h_j * h_j) continue; /* legacy filter: only contribute when in j's kernel */
                double u, hinv, hinv3, hinv4, wk = 0, dwk = 0;
                kernel_hinv(h_j, &hinv, &hinv3, &hinv4);
                u = sqrt(r2) * hinv;
                kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
                accum_V1 += CellP[j].Volume_0 * CellP[j].Volume_0 * wk;
            }
            CellP[i].Volume_1 += accum_V1;
        }

        if (num_src > 0) {
            gpu_ngb_list_free(&gnl, NULL);
            gpu_particles_arena_invalidate();
        }
        if (imported_ghosts) ghost_exchange_cleanup();
    }
    timecomp = timediff(t_kern_start, my_second());
    cellcorrections_final_operations_and_cleanup(); /* do final operations on results */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_DENSCOMPUTE] += timecomp; CPU_Step[CPU_DENSWAIT] += timewait; CPU_Step[CPU_DENSCOMM] += timecomm;
    CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm); /* collect timings and reset clock for next timing */
}

#endif // parent if statement for all code in the HYDRO_VOLUME_CORRECTIONS block
