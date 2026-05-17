/* grain_physics_gpu.cc — GPU neighbor-loop port for grain_backrx_evaluate (B7a)
 * and interpolate_fluxes_opacities_gasgrains_evaluate (B7b).
 *
 * B7a: grain→gas momentum backreaction with atomic j-writes + ghost writeback.
 * B7b: gas↔grain RT opacity coupling, per-source output only, two kernels.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../system/gpu_particles_arena.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"

#include "../declarations/gpu_numeric_macros.h"
#include "../declarations/gpu_error_check.h"
#include "../declarations/gpu_dispatch_templates.h"
#include "../declarations/macros.h"
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/ghost_writeback.h"
#include "../mesh/ghost_symlist_lifecycle.h"
#include "grain_physics_gpu.h"

#if defined(GRAIN_FLUID)

#include "grain_physics_functions.h"


/* ========================================================================
 *  B7a — grain_backrx (GRAIN_BACKREACTION)
 * ======================================================================== */
#if defined(GRAIN_BACKREACTION)
static void grain_backrx_local_fill(int i,
                                     struct particle_data *P_host,
                                     struct GrainBackrxLocalIn *loc)
{
    loc->Pos = P_host[i].Pos;
    loc->KernelRadius = P_host[i].KernelRadius;
    loc->Grain_DeltaMomentum = P_host[i].Grain_DeltaMomentum;
    loc->Gas_Density = P_host[i].Gas_Density;
    loc->Grain_AccelTimeMin = P_host[i].Grain_AccelTimeMin;
#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & (32|64))
    for(int kv = 0; kv < GRAIN_NUM_VOLATILE_SPECIES; kv++) { loc->Grain_DeltaVolatileMass[kv] = P_host[i].Grain_DeltaVolatileMass[kv]; }
    loc->Grain_DeltaInternalEnergyHeating = P_host[i].Grain_DeltaInternalEnergyHeating;
#endif
}
#endif

/* ================================================================
   GPU grain-physics evaluators (LATENT for fire_m11i — GRAIN_FLUID-only)
   ----------------------------------------------------------------
   Kernel writes (host-visible, by index):
     grain_backrx_pair_kernel (Case 2 — sparse, 4 fields):
       P_gpu[j].Vel[k]                — Kokkos::atomic_add
       P_gpu[j].dp[k]                 — Kokkos::atomic_add
       P_gpu[j].Grain_AccelTimeMin    — Kokkos::atomic_min
       CellP_gpu[j].VelPred[k]        — Kokkos::atomic_add

     gasgrain_rt_gas_search_pair_kernel  — PURE READ-ONLY (writes only
                                            per-source out_arr[aa]).
     gasgrain_rt_grain_search_pair_kernel — PURE READ-ONLY (same).

   Sparse-scatter targets:
     - grain_backrx evaluator: walk neighbors[] CSR, scatter the 4 fields
       above (Case 2).
     - interpolate_fluxes_opacities_gasgrains evaluator: kernels are
       read-only, full P/CellP scatter is cargo-cult — DELETE (Case 1).
   ================================================================ */
void grain_backrx_evaluate_gpu(struct particle_data *P_host,
                                struct gas_cell_data *CellP_host,
                                int num_total,
                                int *i_active_host, int num_active,
                                const double *src_radii_host)
{
#if defined(GRAIN_BACKREACTION)
    GIZMO_GPU_ENSURE_ALL_FRESH();

    int num_src = num_active;

    int imported_ghosts = 0;
    {
        int need_import_local = (ghost_get_num_ghosts() <= 0) ? 1 : 0;
        int need_import = 0;
        MPI_Allreduce(&need_import_local, &need_import, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(need_import) {
            if(ghost_get_num_ghosts() > 0) ghost_exchange_cleanup();
            gizmo_density_prep_ghosts(gizmo_ghost_safety_factor());
            imported_ghosts = 1;
        }
    }

    int num_all = ghost_get_num_local() + ghost_get_num_ghosts();
    if(num_all <= 0) num_all = num_total;

    /* Wrapper fast-path: no source grains → skip arena/allocs/kernel/scatter.
     * Preserve ghost_writeback collectives (each self-guards). */
    if(num_src == 0) {
        ghost_writeback_zero_grainbackrx();
        ghost_writeback_grainbackrx();
        if(imported_ghosts) ghost_exchange_cleanup();
        return;
    }

    /* Step 13 Phase 1 arena. */
    gpu_particles_arena_acquire(num_all, P_host, CellP_host);
    struct particle_data *P_gpu = gpu_particles_arena_P();
    struct gas_cell_data *CellP_gpu = gpu_particles_arena_CellP();

    int alloc_n = (num_src > 0) ? num_src : 1;
    struct GrainBackrxLocalIn *d_local = (struct GrainBackrxLocalIn *)
        Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(struct GrainBackrxLocalIn));
    for(int a = 0; a < num_src; a++) {
        grain_backrx_local_fill(i_active_host[a], P_host, &d_local[a]);
    }

    /* Neighbor list: grain sources (i) → gas neighbors (j_type_bitmask = 1) */
    gpu_neighbor_list_t gnl;
    gpu_ngb_list_build(P_gpu, num_all,
                       i_active_host, num_src,
                       NGB_SEARCH_ONEWAY, 1 /* gas only */,
                       &gnl, NULL, 1.0, src_radii_host);

    PRINT_STATUS("  GPU grain_backrx: %d sources, %lld pairs", num_src, (long long)gnl.total_pairs);

    /* Snapshot ghost fields we will atomically modify. */
    ghost_writeback_zero_grainbackrx();

    {
        int64_t *offsets = gnl.offsets;
        int *neighbors = gnl.neighbors;
        struct GrainBackrxLocalIn *local_arr = d_local;
        struct particle_data *kp = P_gpu;
        struct gas_cell_data *kc = CellP_gpu;

        Kokkos::parallel_for("grain_backrx", num_src, KOKKOS_LAMBDA(int aa) {
            const struct GrainBackrxLocalIn& loc = local_arr[aa];
            if(loc.KernelRadius <= 0) return;
            double h2 = (double)loc.KernelRadius * (double)loc.KernelRadius;

            int64_t start = offsets[aa], end = offsets[aa+1];
            for(int64_t nn = start; nn < end; nn++) {
                int j = neighbors[nn];
                if(kp[j].Mass <= 0) continue;
                Vec3<double> dp = loc.Pos - kp[j].Pos;
                nearest_xyz(dp);
                double r2 = dp.norm_sq();
                if(r2 <= 0 || r2 >= h2) continue;
                grain_backrx_pair_kernel(loc, j, kp, kc, dp, r2);
            }
        });
    }
    Kokkos::fence();
    gizmo_gpu_check_last_error("grain_backrx", num_src);

    /* Row 6d.i of arena-scope sweep: sparse scatter for grain_backrx. Per
     * audit (a06e30ca), grain_backrx_pair_kernel writes only:
     *   P_gpu[j].Vel[k]                 — atomic_add
     *   CellP_gpu[j].VelPred[k]         — atomic_add
     *   P_gpu[j].dp[k]                  — atomic_add
     *   P_gpu[j].Grain_AccelTimeMin     — atomic_min
     * Sparse scatter over gnl.neighbors[] (deep-copied to host). */
    if(gnl.total_pairs > 0) {
        std::vector<int> gnl_neighbors_host((size_t)gnl.total_pairs);
        gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
        for(int64_t idx = 0; idx < gnl.total_pairs; idx++) {
            int j = gnl_neighbors_host[idx];
            P_host[j].Vel               = P_gpu[j].Vel;
            P_host[j].dp                = P_gpu[j].dp;
            P_host[j].Grain_AccelTimeMin = P_gpu[j].Grain_AccelTimeMin;
            CellP_host[j].VelPred       = CellP_gpu[j].VelPred;
        }
    }

    ghost_writeback_grainbackrx();

    /* Full scatter syncs arena→host; ghost_writeback_grainbackrx invalidates. */
    gpu_ngb_list_free(&gnl, NULL);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_local);

    if(imported_ghosts) ghost_exchange_cleanup();
#else
    (void)P_host; (void)CellP_host; (void)num_total;
    (void)i_active_host; (void)num_active; (void)src_radii_host;
#endif /* GRAIN_BACKREACTION */
}


/* ========================================================================
 *  B7b — interpolate_fluxes_opacities_gasgrains (RT_OPACITY_FROM_EXPLICIT_GRAINS)
 * ======================================================================== */
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
static void gasgrain_rt_local_fill(int i,
                                    struct particle_data *P_host,
                                    struct gas_cell_data *CellP_host,
                                    struct GasGrainRTLocalIn *loc)
{
    loc->Type = P_host[i].Type;
    loc->Pos = P_host[i].Pos;
    loc->Vel = P_host[i].Vel;
    loc->KernelRadius = P_host[i].KernelRadius;
    loc->Mass = P_host[i].Mass;
    loc->Grain_Size = 0;
    for(int k = 0; k < N_RT_FREQ_BINS; k++) loc->Grain_Abs_Coeff[k] = 0;
    if(((1 << P_host[i].Type) & GRAIN_PTYPES) && P_host[i].Gas_Density > 0) {
        loc->Grain_Size = P_host[i].Grain_Size;
        double R_grain_code  = (double)P_host[i].Grain_Size / UNIT_LENGTH_IN_CGS;
        double rho_grain_code = All.Grain_Internal_Density / UNIT_DENSITY_IN_CGS;
        double rho_gas_code  = (double)P_host[i].Gas_Density * All.cf_a3inv;
        for(int k_freq = 0; k_freq < N_RT_FREQ_BINS; k_freq++) {
            double Q_abs = grain_extinction_Q_inline(i, k_freq, P_host);
            loc->Grain_Abs_Coeff[k_freq] = (MyFloat)(Q_abs * 3.0 /
                (4.0 * C_LIGHT_CODE_REDUCED * rho_grain_code * R_grain_code * rho_gas_code));
        }
    }
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
    /* For gas sources: pack electron fraction + density for J_dust Z_grain estimate. */
    loc->Ne_gas      = (P_host[i].Type == 0) ? CellP_host[i].Ne : 0.0f;
    loc->Density_gas = (P_host[i].Type == 0) ? CellP_host[i].Density : 0.0f;
#endif
}
#endif

void interpolate_fluxes_opacities_gasgrains_evaluate_gpu(struct particle_data *P_host,
                                                          struct gas_cell_data *CellP_host,
                                                          int num_total,
                                                          int *i_active_gas_host, int num_active_gas,
                                                          const double *src_radii_gas_host,
                                                          int *i_active_grain_host, int num_active_grain,
                                                          const double *src_radii_grain_host)
{
#if defined(RT_OPACITY_FROM_EXPLICIT_GRAINS)
    GIZMO_GPU_ENSURE_ALL_FRESH();

    int imported_ghosts = 0;
    {
        int need_import_local = (ghost_get_num_ghosts() <= 0) ? 1 : 0;
        int need_import = 0;
        MPI_Allreduce(&need_import_local, &need_import, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(need_import) {
            if(ghost_get_num_ghosts() > 0) ghost_exchange_cleanup();
            gizmo_density_prep_ghosts(gizmo_ghost_safety_factor());
            imported_ghosts = 1;
        }
    }

    int num_all = ghost_get_num_local() + ghost_get_num_ghosts();
    if(num_all <= 0) num_all = num_total;

    /* Wrapper fast-path: both directions empty → skip arena+allocs+kernels.
     * No multi-rank writeback collective in this function (per-source updates
     * scatter directly via the per-direction loops below). */
    if(num_active_gas == 0 && num_active_grain == 0) {
        if(imported_ghosts) ghost_exchange_cleanup();
        return;
    }

    /* Step 13 Phase 1 arena. */
    gpu_particles_arena_acquire(num_all, P_host, CellP_host);
    struct particle_data *P_gpu = gpu_particles_arena_P();
    struct gas_cell_data *CellP_gpu = gpu_particles_arena_CellP();

    /* -------------------- Direction 1: gas → grains -------------------- */
    {
        int num_src = num_active_gas;
        int alloc_n = (num_src > 0) ? num_src : 1;
        struct GasGrainRTLocalIn *d_local = (struct GasGrainRTLocalIn *)
            Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(struct GasGrainRTLocalIn));
        struct GasGrainRTOut *d_out = (struct GasGrainRTOut *)
            Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(struct GasGrainRTOut));
        memset(d_out, 0, alloc_n * sizeof(struct GasGrainRTOut));
        for(int a = 0; a < num_src; a++) {
            gasgrain_rt_local_fill(i_active_gas_host[a], P_host, CellP_host, &d_local[a]);
        }

        gpu_neighbor_list_t gnl;
        /* Kernel at line 303 rejects on r >= h_j (not h_i). With ONEWAY the BVH
         * walks at h_i and would miss valid pairs where h_j > h_i; per-Phil's
         * "when in doubt sym" rule + the kernel's h_j filter means the right
         * mode is SYMMETRIC, kernel keeps the h_j accept. */
        gpu_ngb_list_build(P_gpu, num_all,
                           i_active_gas_host, num_src,
                           NGB_SEARCH_SYMMETRIC, GRAIN_PTYPES,
                           &gnl, NULL, 1.0, src_radii_gas_host);

        PRINT_STATUS("  GPU gasgrain_rt (gas->grains): %d sources, %lld pairs", num_src, (long long)gnl.total_pairs);

        {
            int64_t *offsets = gnl.offsets;
            int *neighbors = gnl.neighbors;
            struct GasGrainRTLocalIn *local_arr = d_local;
            struct GasGrainRTOut *out_arr = d_out;
            struct particle_data *kp = P_gpu;
            struct gas_cell_data *kc = CellP_gpu;

            Kokkos::parallel_for("gasgrain_rt_gas_src", num_src, KOKKOS_LAMBDA(int aa) {
                const struct GasGrainRTLocalIn& loc = local_arr[aa];
                if(loc.KernelRadius <= 0) return;
                struct GasGrainRTOut myout;
                memset(&myout, 0, sizeof(myout));

                int64_t start = offsets[aa], end = offsets[aa+1];
                for(int64_t nn = start; nn < end; nn++) {
                    int j = neighbors[nn];
                    if(kp[j].Mass <= 0) continue;
                    Vec3<double> dp = loc.Pos - kp[j].Pos;
                    nearest_xyz(dp);
                    double r2 = dp.norm_sq();
                    if(r2 <= 0) continue;
                    double hj = (double)kp[j].KernelRadius;
                    if(hj <= 0 || r2 >= hj * hj) continue;
                    gasgrain_rt_gas_search_pair_kernel(loc, j, kp, kc, dp, r2, myout);
                }
                out_arr[aa] = myout;
            });
        }
        Kokkos::fence();
        gizmo_gpu_check_last_error("gasgrain_rt_gas_src", num_src);

        /* Scatter per-gas output back to CellP[i] */
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
        /* Reset J_dust_cell once per pass (active sources only) so the
           accumulation below replaces any stale value from a previous step. */
        for(int a = 0; a < num_src; a++) {
            int i = i_active_gas_host[a];
            CellP_host[i].J_dust_cell = {};
        }
#endif
        for(int a = 0; a < num_src; a++) {
            int i = i_active_gas_host[a];
            CellP_host[i].InterpolatedGeometricDustCrossSection += d_out[a].InterpolatedGeometricDustCrossSection;
            for(int k_freq = 0; k_freq < N_RT_FREQ_BINS; k_freq++) {
                CellP_host[i].Interpolated_Opacity[k_freq] += d_out[a].Interpolated_Opacity[k_freq];
            }
#if defined(MHD_BATTERY_MECHANISMS) && (MHD_BATTERY_MECHANISMS & 8)
            for(int kd = 0; kd < 3; kd++) {
                CellP_host[i].J_dust_cell[kd] += d_out[a].J_dust_contribution[kd];
            }
#endif
        }

        gpu_ngb_list_free(&gnl, NULL);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_out);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_local);
    }

    /* -------------------- Direction 2: grains → gas -------------------- */
    {
        int num_src = num_active_grain;
        int alloc_n = (num_src > 0) ? num_src : 1;
        struct GasGrainRTLocalIn *d_local = (struct GasGrainRTLocalIn *)
            Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(struct GasGrainRTLocalIn));
        struct GasGrainRTOut *d_out = (struct GasGrainRTOut *)
            Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(alloc_n * sizeof(struct GasGrainRTOut));
        memset(d_out, 0, alloc_n * sizeof(struct GasGrainRTOut));
        for(int a = 0; a < num_src; a++) {
            gasgrain_rt_local_fill(i_active_grain_host[a], P_host, CellP_host, &d_local[a]);
        }

        gpu_neighbor_list_t gnl;
        gpu_ngb_list_build(P_gpu, num_all,
                           i_active_grain_host, num_src,
                           NGB_SEARCH_ONEWAY, 1 /* gas only */,
                           &gnl, NULL, 1.0, src_radii_grain_host);

        PRINT_STATUS("  GPU gasgrain_rt (grains->gas): %d sources, %lld pairs", num_src, (long long)gnl.total_pairs);

        {
            int64_t *offsets = gnl.offsets;
            int *neighbors = gnl.neighbors;
            struct GasGrainRTLocalIn *local_arr = d_local;
            struct GasGrainRTOut *out_arr = d_out;
            struct particle_data *kp = P_gpu;
            struct gas_cell_data *kc = CellP_gpu;

            Kokkos::parallel_for("gasgrain_rt_grain_src", num_src, KOKKOS_LAMBDA(int aa) {
                const struct GasGrainRTLocalIn& loc = local_arr[aa];
                if(loc.KernelRadius <= 0) return;
                double h2 = (double)loc.KernelRadius * (double)loc.KernelRadius;
                struct GasGrainRTOut myout;
                memset(&myout, 0, sizeof(myout));

                int64_t start = offsets[aa], end = offsets[aa+1];
                for(int64_t nn = start; nn < end; nn++) {
                    int j = neighbors[nn];
                    if(kp[j].Mass <= 0) continue;
                    Vec3<double> dp = loc.Pos - kp[j].Pos;
                    nearest_xyz(dp);
                    double r2 = dp.norm_sq();
                    if(r2 <= 0 || r2 >= h2) continue;
                    gasgrain_rt_grain_search_pair_kernel(loc, j, kp, kc, dp, r2, myout);
                }
                out_arr[aa] = myout;
            });
        }
        Kokkos::fence();
        gizmo_gpu_check_last_error("gasgrain_rt_grain_src", num_src);

        /* Scatter per-grain output back to P[i].GravAccel */
        for(int a = 0; a < num_src; a++) {
            int i = i_active_grain_host[a];
            P_host[i].GravAccel += d_out[a].Interpolated_Radiation_Acceleration / All.cf_a2inv;
        }

        gpu_ngb_list_free(&gnl, NULL);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_out);
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_local);
    }

    /* Row 6d.ii of arena-scope sweep: kernel-writes audit (a06e30ca) shows
     * gasgrain_rt_gas_search_pair_kernel and gasgrain_rt_grain_search_pair_kernel
     * are PURE READ-ONLY on P/CellP — they only write to per-source out_arr[aa].
     * The former full memcpy(P_host, P_gpu, num_all*...) +
     * memcpy(CellP_host, CellP_gpu, num_all*...) here was cargo-cult — DELETED. */

    /* Defensive invalidate for any caller post-mods. */
    gpu_particles_arena_invalidate();

    if(imported_ghosts) ghost_exchange_cleanup();
#else
    (void)P_host; (void)CellP_host; (void)num_total;
    (void)i_active_gas_host; (void)num_active_gas; (void)src_radii_gas_host;
    (void)i_active_grain_host; (void)num_active_grain; (void)src_radii_grain_host;
#endif /* RT_OPACITY_FROM_EXPLICIT_GRAINS */
}



#else /* fallback stubs */

void grain_backrx_evaluate_gpu(struct particle_data *p, struct gas_cell_data *cp, int num_total,
                                int *i_active_host, int num_active, const double *src_radii_host)
{ (void)p; (void)cp; (void)num_total; (void)i_active_host; (void)num_active; (void)src_radii_host; }
void interpolate_fluxes_opacities_gasgrains_evaluate_gpu(struct particle_data *p, struct gas_cell_data *cp, int num_total,
                                                          int *ig, int ng, const double *rg,
                                                          int *ir, int nr, const double *rr)
{ (void)p; (void)cp; (void)num_total; (void)ig; (void)ng; (void)rg; (void)ir; (void)nr; (void)rr; }

#endif /* GRAIN_FLUID */
