/* dm_fuzzy_gpu.cc — GPU-accelerated DMGrad (density-gradient estimator for DM_FUZZY).
 *
 * GPU translation unit: (Kokkos/nvcc_wrapper). Called once per iteration from DMGrad_gradient_calc() in
 * dm_fuzzy.cc. Builds a cross-type CSR (i = DM active, j = DM via
 * bitmask group), runs the iteration-specific accumulation kernel.
 *
 * No j-writes, no RNG, no fragment dependencies — simplest of the
 * SIDM/fuzzy ports. Mirrors the B1 AGS density pattern.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../system/gpu_particles_arena.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/neighbor_list.h"

#include "../declarations/gpu_numeric_macros.h"
#include "../declarations/gpu_error_check.h"
#include "../declarations/gpu_dispatch_templates.h"

#include "sidm_gpu_decls.h"


#if defined(DM_FUZZY)

/* Device-side min/max helpers that update in place (MINMAX_CHECK equivalent). */
KOKKOS_INLINE_FUNCTION void dmgrad_minmax(double v, double &mn, double &mx)
{
    if(v < mn) mn = v;
    if(v > mx) mx = v;
}


void dmgrad_evaluate_gpu(struct particle_data *P_host, int num_total,
                         int *i_active_host, int num_active,
                         const double *i_radii_host,
                         const struct dmgrad_gpu_in *i_inputs_host,
                         int j_type_bitmask,
                         int loop_iteration,
                         void *out_host_void)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();
    struct dmgrad_gpu_out *out_host = (struct dmgrad_gpu_out *)out_host_void;

    /* Wrapper fast-path: no active sources → skip arena/allocs/kernel/scatter. */
    if(num_active == 0) { (void)out_host; (void)loop_iteration; return; }

    /* Step 13 Phase 1 arena. P-only kernel (DM-only loop). */
    gpu_particles_arena_acquire(num_total, P_host, NULL);
    struct particle_data *P_gpu = gpu_particles_arena_P();

    /* Build cross-type CSR: i = caller's group, j_type_bitmask filters j-side */
    gpu_neighbor_list_t gnl;
    gpu_ngb_list_build(P_gpu, num_total, i_active_host, num_active,
                       NGB_SEARCH_ONEWAY, j_type_bitmask, &gnl, NULL,
                       1.0 /* search_radius_factor */, i_radii_host, NULL, "sidm");

    /* SharedSpace mirrors of per-active inputs + radii + outputs */
    struct dmgrad_gpu_in *d_in = (struct dmgrad_gpu_in *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
        ((num_active > 0) ? num_active : 1) * sizeof(struct dmgrad_gpu_in));
    double *d_radii = (double *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
        ((num_active > 0) ? num_active : 1) * sizeof(double));
    struct dmgrad_gpu_out *d_out = (struct dmgrad_gpu_out *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
        ((num_active > 0) ? num_active : 1) * sizeof(struct dmgrad_gpu_out));
    memcpy(d_in, i_inputs_host, num_active * sizeof(struct dmgrad_gpu_in));
    memcpy(d_radii, i_radii_host, num_active * sizeof(double));

    PRINT_STATUS("  GPU DMGrad iter=%d: %d active, %lld DM pairs",
                 loop_iteration, num_active, (long long)gnl.total_pairs);

    {
        int64_t *offsets = gnl.offsets;
        int *neighbors = gnl.neighbors;
        int *active = gnl.d_active;
        struct particle_data *kp = P_gpu;
        double *radii = d_radii;
        struct dmgrad_gpu_in *kin = d_in;
        struct dmgrad_gpu_out *kout = d_out;
        int iter = loop_iteration;

        gizmo_gpu_kernel_launch("dmgrad_kernel", num_active, KOKKOS_LAMBDA(int aa) {
            int ii = active[aa];
            double h_i = radii[aa];
            /* Zero the full output struct each call — unused fields stay at zero */
            memset(&kout[aa], 0, sizeof(struct dmgrad_gpu_out));
            if(!(h_i > 0) || !(kin[aa].AGS_Density > 0)) return;

            double h2_i = h_i * h_i;
            double hinv, hinv3, hinv4;
            kernel_hinv(h_i, &hinv, &hinv3, &hinv4);

            Vec3<double> pos_i = kp[ii].Pos;
            double rho_i = kin[aa].AGS_Density;
#if (DM_FUZZY > 0)
            double psi_re_i = kin[aa].AGS_Psi_Re;
            double psi_im_i = kin[aa].AGS_Psi_Im;
#endif

            /* Local accumulators mirror the CPU per-pair body */
            double g_rho[3] = {0,0,0};
            double max_rho = 0, min_rho = 0;
#if (DM_FUZZY > 0)
            double g_psi_re[3] = {0,0,0};
            double g_psi_im[3] = {0,0,0};
#endif
            double g2_rho[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
            double max_grho[3] = {0,0,0}, min_grho[3] = {0,0,0};
#if (DM_FUZZY > 0)
            double g2_psi_re[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
            double g2_psi_im[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
#endif

            int64_t start = offsets[aa], end = offsets[aa + 1];
            for(int64_t nn = start; nn < end; nn++) {
                int j = neighbors[nn];
                if((kp[j].Mass <= 0) || (kp[j].AGS_Density <= 0)) continue;
                Vec3<double> dp;
                for(int kk = 0; kk < 3; kk++) dp[kk] = pos_i[kk] - kp[j].Pos[kk];
                nearest_xyz(dp);
                double r2 = dp.norm_sq();
                if(!(r2 > 0) || !(r2 < h2_i)) continue;
                double r = sqrt(r2);
                double u = r * hinv;
                double wk_i, dwk_i;
                kernel_main(u, hinv3, hinv4, &wk_i, &dwk_i, -1);

                if(iter <= 0) {
                    /* Iter 0: gradient of density-like scalars */
                    double d_rho = kp[j].AGS_Density - rho_i;
                    dmgrad_minmax(d_rho, min_rho, max_rho);
                    for(int k = 0; k < 3; k++) g_rho[k] += -wk_i * dp[k] * d_rho;
#if (DM_FUZZY > 0)
                    double d_re = kp[j].AGS_Psi_Re_Pred * kp[j].AGS_Density / kp[j].Mass - psi_re_i;
                    for(int k = 0; k < 3; k++) g_psi_re[k] += -wk_i * dp[k] * d_re;
                    double d_im = kp[j].AGS_Psi_Im_Pred * kp[j].AGS_Density / kp[j].Mass - psi_im_i;
                    for(int k = 0; k < 3; k++) g_psi_im[k] += -wk_i * dp[k] * d_im;
#endif
                } else {
                    /* Iter 1: gradient-of-gradient tensor */
                    for(int k = 0; k < 3; k++) {
                        double d_grad_rho = kp[j].AGS_Gradients_Density[k] - kin[aa].AGS_Gradients_Density[k];
                        dmgrad_minmax(d_grad_rho, min_grho[k], max_grho[k]);
                        for(int k2 = 0; k2 < 3; k2++) g2_rho[k2][k] += -wk_i * dp[k2] * d_grad_rho;
#if (DM_FUZZY > 0)
                        double d_gre = kp[j].AGS_Gradients_Psi_Re[k] - kin[aa].AGS_Gradients_Psi_Re[k];
                        for(int k2 = 0; k2 < 3; k2++) g2_psi_re[k2][k] += -wk_i * dp[k2] * d_gre;
                        double d_gim = kp[j].AGS_Gradients_Psi_Im[k] - kin[aa].AGS_Gradients_Psi_Im[k];
                        for(int k2 = 0; k2 < 3; k2++) g2_psi_im[k2][k] += -wk_i * dp[k2] * d_gim;
#endif
                    }
                }
            } /* end neighbor loop */

            /* Scatter locals into output struct */
            if(iter <= 0) {
                for(int k = 0; k < 3; k++) kout[aa].grad_rho[k] = g_rho[k];
                kout[aa].max_rho = max_rho; kout[aa].min_rho = min_rho;
#if (DM_FUZZY > 0)
                for(int k = 0; k < 3; k++) {kout[aa].grad_psi_re[k] = g_psi_re[k]; kout[aa].grad_psi_im[k] = g_psi_im[k];}
#endif
            } else {
                for(int k = 0; k < 3; k++) {
                    kout[aa].max_grho[k] = max_grho[k];
                    kout[aa].min_grho[k] = min_grho[k];
                    for(int k2 = 0; k2 < 3; k2++) kout[aa].grad2_rho[k2][k] = g2_rho[k2][k];
#if (DM_FUZZY > 0)
                    for(int k2 = 0; k2 < 3; k2++) {
                        kout[aa].grad2_psi_re[k2][k] = g2_psi_re[k2][k];
                        kout[aa].grad2_psi_im[k2][k] = g2_psi_im[k2][k];
                    }
#endif
                }
            }
        });
    }

    memcpy(out_host, d_out, num_active * sizeof(struct dmgrad_gpu_out));

    /* Read-only on arena; out post-scattered by caller — invalidate. */
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_out);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_in);
    gpu_ngb_list_free(&gnl, NULL);
    gpu_particles_arena_invalidate();
}


#else /* stubs */

void dmgrad_evaluate_gpu(struct particle_data *, int, int *, int, const double *,
                         const struct dmgrad_gpu_in *, int, int, void *) {}

#endif /* DM_FUZZY */
