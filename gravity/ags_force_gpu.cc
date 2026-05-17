/* ags_force_gpu.cc — GPU-accelerated AGSForce loop (B2 of the GIZMO GPU
 * porting plan).
 *
 * Pairs every AGSForce-active i against its AGS-kernel neighbours, runs the
 * three flux-compute_pair functions (cbe_integrator_flux_compute_pair,
 * dm_fuzzy_flux_compute_pair, sidm_core_flux_compute_pair), and writes both
 * i-side accumulators and j-side deltas in one pass.
 *
 * i-side results are scattered back on the host; j-side deltas (Vel, dp,
 * NInteractions, wakeup) are applied via Kokkos atomics inside the kernel,
 * and the host driver wraps this dispatch with ghost_writeback_{zero_,}agsforce
 * so ghost-side modifications are reverse-communicated to their home ranks.
 *
 * Search geometry matches the CPU tree-walk:
 *   - DM_SIDM:      3x inflated i-side search radius, pair filter r < h_i + h_j
 *   - else (CBE,    1x search radius, pair filter r < max(h_i, h_j).  This
 *     DM_FUZZY):    matches the CPU ngb_treefind_pairs_threads_targeted call
 *                   for non-gas types where P[j].KernelRadius is typically 0,
 *                   reducing SEARCHBOTHWAYS to a ONEWAY walk with h_i.
 *
 * RNG: SIDM scatter uses the counter-based generator in declarations/gpu_rng.h
 * (see sidm_core_flux_functions.h for keying). No GSL on device.
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
#include "../mesh/gpu_neighbor_list.h"
#include "../mesh/neighbor_list.h"

#include "../declarations/gpu_numeric_macros.h"
#include "../declarations/gpu_error_check.h"
#include "../declarations/gpu_dispatch_templates.h"

#include "../core/timestep_functions.h"
#include "ags_gpu_decls.h"
#include "ags_functions.h"

#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)

/* Must include AFTER gpu_all_mirror.h redirects All -> All_dev and AFTER the
   KOKKOS_INLINE_FUNCTION macro definition so the helpers expand with GPU
   device annotations. */
#include "../sidm/cbe_integrator_flux_functions.h"
#include "../sidm/dm_fuzzy_flux_functions.h"
#include "../sidm/sidm_core_flux_functions.h"

/* File-scope named struct for TimeBinActive device-capture. CUDA nvcc rejects
   local/unnamed types in device-lambda captures, so we lift the capture
   wrapper out to file scope here. */
struct ags_force_tba_cap_t { int v[TIMEBINS]; };


/* GPU-kernel local struct that mirrors the AGSForce INPUT_STRUCT_NAME field set,
   filled on-device from kp[ii]. Keeps the flux-pair templates happy without
   copying the full host struct across TUs. */
struct ags_force_local_t {
    double Mass;
    double AGS_KernelRadius;
    Vec3<double> Pos;
    Vec3<double> Vel;
    int Type;
    double dtime;
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    MyDouble NV_T[3][3];
    double V_i;
#endif
#if defined(DM_FUZZY)
    Vec3<double> AGS_Gradients_Density;
    double AGS_Gradients2_Density[3][3];
    double AGS_Numerical_QuantumPotential;
#if (DM_FUZZY > 0)
    double AGS_Psi_Re; Vec3<double> AGS_Gradients_Psi_Re; Mat3<double> AGS_Gradients2_Psi_Re;
    double AGS_Psi_Im; Vec3<double> AGS_Gradients_Psi_Im; Mat3<double> AGS_Gradients2_Psi_Im;
#endif
#endif
#if defined(CBE_INTEGRATOR)
    double CBE_basis_moments[CBE_INTEGRATOR_NBASIS][CBE_INTEGRATOR_NMOMENTS];
#endif
#if defined(DM_SIDM)
    double dtime_sidm;
    MyIDType ID;
#ifdef GRAIN_COLLISIONS
    double Grain_CrossSection_PerUnitMass;
#endif
#endif
#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & 7)
    /* Phase 17b pairwise resolver inputs. Composition[] needed for
     * composition-mixing on the absorber side; Grain_Size needed for the
     * Dominik-Tielens v_coag threshold and for size-update bookkeeping. */
    double Grain_Size;
    double Composition[GRAIN_NUM_SPECIES];
#endif
};


/* Kernel struct matching ags_rkern.cc's kernel_AGSForce, but instantiated
   per-pair inside the GPU kernel. Fields filled by the caller's per-pair code. */
struct ags_force_kernel_t {
    Vec3<double> dp, dv;
    double r, wk_i, wk_j, dwk_i, dwk_j;
    double h_i, hinv_i, hinv3_i, hinv4_i;
    double h_j, hinv_j, hinv3_j, hinv4_j;
};


/* Per-i output struct used on-device; matches ags_force_gpu_out layout. */
typedef struct ags_force_gpu_out ags_force_dev_out_t;


/* ================================================================
   GPU AGS-force evaluator
   ----------------------------------------------------------------
   Kernel writes (host-visible, by index):
     i-side: kout[aa] (per-active output: AGS-corrections, gravitational
                       softening sums)
     j-side (neighbors, indexed by j = neighbors[nn]):
       #if defined(CBE_INTEGRATOR)
       P_gpu[j].wakeup                        — atomic_max(local.TimeBin+1)
                                                  (CBE wakeup; hydro convention)
       #endif
       #if defined(DM_SIDM)
       P_gpu[j].wakeup                        — atomic_max(local.TimeBin+1)
                                                  (SIDM wakeup; hydro convention)
       P_gpu[j].Vel[kv]                       — atomic_add (k = 0..2)
       P_gpu[j].dp[kv]                        — atomic_add (k = 0..2)
       P_gpu[j].NInteractions                 — atomic_add(1)
       #endif
     Without DM_SIDM and without CBE_INTEGRATOR: NO j-side writes.

   Sparse-scatter target (legacy, pre-UVM):
     - Under UVM-canonical particles (commit 0d9e74b4), P_gpu/P_host are pointer
       aliases of the same SharedSpace allocation; device atomic writes ARE the
       host writes. The DM_SIDM-gated sparse scatter below is a no-op under UVM
       but preserved for documentation / discrete-memory fallback.
   ================================================================ */
void ags_force_evaluate_gpu(struct particle_data *P_host,
                            int num_total,
                            int *i_active_host, int num_active,
                            const double *i_radii_host,
                            int j_type_bitmask,
                            struct ags_force_gpu_out *out_host)
{
    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* Wrapper fast-path: caller-side ags_density_isactive() filter has already
     * produced num_active.  When 0, skip arena/allocs/kernel/scatter. */
    if(num_active == 0) { (void)out_host; return; }

    /* Step 13 Phase 1 arena. P-only kernel; pass NULL CellP. j-writes go via
     * Kokkos atomics to the arena P, then full P scatter back at the end. */
    gpu_particles_arena_set_site("ags_force_gpu");
    gpu_particles_arena_acquire(num_total, P_host, NULL);
    struct particle_data *P_gpu = gpu_particles_arena_P();

    /* SIDM inflates search radius 3x (matches CPU AGSForce tree-walk). */
    double sr_fac = 1.0;
#if defined(DM_SIDM)
    sr_fac = 3.0;
#endif

    gpu_neighbor_list_t gnl;
    /* Kernel at lines 292,294 rejects only when both r > h_i AND r > h_j (or
     * r > h_i+h_j); that's the symmetric pair predicate, so neighbor pool must
     * include j with h_j > h_i and r > h_i. Was ONEWAY which silently dropped
     * those pairs. Same class of bug Codex caught for sink/mech_fb. */
    gpu_ngb_list_build(P_gpu, num_total, i_active_host, num_active,
                       NGB_SEARCH_SYMMETRIC, j_type_bitmask, &gnl, NULL,
                       sr_fac, i_radii_host, NULL, "ags_force");

    double *d_radii = (double *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
        ((num_active > 0) ? num_active : 1) * sizeof(double));
    memcpy(d_radii, i_radii_host, num_active * sizeof(double));

    ags_force_dev_out_t *d_out = (ags_force_dev_out_t *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
        ((num_active > 0) ? num_active : 1) * sizeof(ags_force_dev_out_t));

#if defined(DM_SIDM)
    /* GeoFactorTable mirror for g_geo_tab. 1000 doubles (~8KB) is trivial. */
    MyDouble *d_geofactor = (MyDouble *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(
        GEOFACTOR_TABLE_LENGTH * sizeof(MyDouble));
    memcpy(d_geofactor, GeoFactorTable, GEOFACTOR_TABLE_LENGTH * sizeof(MyDouble));
#endif

    /* TimeBinActive is needed by sidm_core_flux_compute_pair and by the
       CBE wakeup test. Capture-by-value via the file-scope wrapper struct
       (nvcc requires a named, non-local type on device-lambda captures). */
    ags_force_tba_cap_t tba_cap;
    for(int k = 0; k < TIMEBINS; k++) tba_cap.v[k] = TimeBinActive[k];

    /* Sticky device flag for NeedToWakeupParticles_local (reduced on host). */
    int *d_need_wakeup = (int *) Kokkos::kokkos_malloc<GIZMO_KOKKOS_SHARED_SPACE>(sizeof(int));
    *d_need_wakeup = 0;

    PRINT_STATUS("  GPU AGS-force: %d active, j_bitmask=%d, %lld pairs",
                 num_active, j_type_bitmask, (long long)gnl.total_pairs);

    {
        int64_t *offsets = gnl.offsets;
        int *neighbors = gnl.neighbors;
        int *active = gnl.d_active;
        struct particle_data *kp = P_gpu;
        double *radii = d_radii;
        ags_force_dev_out_t *kout = d_out;
        int *need_wakeup = d_need_wakeup;
#if defined(DM_SIDM)
        const MyDouble *geofactor = d_geofactor;
#endif

        gizmo_gpu_kernel_launch("ags_force_kernel", num_active, KOKKOS_LAMBDA(int aa) {
            int ii = active[aa];
            /* Zero the full output struct — unused fields stay at zero. */
            memset(&kout[aa], 0, sizeof(ags_force_dev_out_t));
            if(!(kp[ii].Mass > 0) || !(radii[aa] > 0)) return;

            /* Build the per-i local struct (mirrors the CPU INPUTFUNCTION_NAME). */
            ags_force_local_t local;
            local.Mass = kp[ii].Mass;
            local.AGS_KernelRadius = radii[aa];
            local.Pos = kp[ii].Pos;
            local.Vel = kp[ii].Vel;
            local.Type = kp[ii].Type;
            local.dtime = get_particle_timestep_in_physical(ii, kp);
            short int local_TimeBin = kp[ii].TimeBin;  /* active i's bin; for hydro-convention wakeup write */
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
            local.V_i = get_particle_volume_ags_P(ii, kp);
            for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) local.NV_T[a][b] = kp[ii].NV_T[a][b];
#endif
#if defined(DM_FUZZY)
            local.AGS_Gradients_Density = kp[ii].AGS_Gradients_Density;
            for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) local.AGS_Gradients2_Density[a][b] = kp[ii].AGS_Gradients2_Density[a][b];
            local.AGS_Numerical_QuantumPotential = kp[ii].AGS_Numerical_QuantumPotential;
#if (DM_FUZZY > 0)
            local.AGS_Psi_Re = kp[ii].AGS_Psi_Re_Pred * kp[ii].AGS_Density / kp[ii].Mass;
            local.AGS_Gradients_Psi_Re = kp[ii].AGS_Gradients_Psi_Re;
            for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) local.AGS_Gradients2_Psi_Re[a][b] = kp[ii].AGS_Gradients2_Psi_Re[a][b];
            local.AGS_Psi_Im = kp[ii].AGS_Psi_Im_Pred * kp[ii].AGS_Density / kp[ii].Mass;
            local.AGS_Gradients_Psi_Im = kp[ii].AGS_Gradients_Psi_Im;
            for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) local.AGS_Gradients2_Psi_Im[a][b] = kp[ii].AGS_Gradients2_Psi_Im[a][b];
#endif
#endif
#if defined(CBE_INTEGRATOR)
            for(int a = 0; a < CBE_INTEGRATOR_NBASIS; a++)
                for(int b = 0; b < CBE_INTEGRATOR_NMOMENTS; b++)
                    local.CBE_basis_moments[a][b] = kp[ii].CBE_basis_moments[a][b];
#endif
#if defined(DM_SIDM)
            local.dtime_sidm = kp[ii].dtime_sidm;
            local.ID = kp[ii].ID;
#ifdef GRAIN_COLLISIONS
            local.Grain_CrossSection_PerUnitMass = return_grain_cross_section_per_unit_mass_P(ii, kp);
#endif
#endif
#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & 7)
            local.Grain_Size = (double)kp[ii].Grain_Size;
            for(int gs = 0; gs < GRAIN_NUM_SPECIES; gs++) { local.Composition[gs] = (double)kp[ii].Composition[gs]; }
#endif

            /* Per-i kernel invariants. */
            ags_force_kernel_t kernel;
            kernel.h_i = local.AGS_KernelRadius;
            kernel_hinv(kernel.h_i, &kernel.hinv_i, &kernel.hinv3_i, &kernel.hinv4_i);

            /* Per-i out (mirrors CPU OUTPUT_STRUCT_NAME — the flux templates
               accumulate into this and we copy it into kout[aa] at the end). */
            ags_force_dev_out_t out;
            memset(&out, 0, sizeof(ags_force_dev_out_t));
#if defined(DM_SIDM)
            out.dtime_sidm = local.dtime_sidm;
#endif
#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & 7)
            out.Grain_DeltaErosionFrac = 1.0; /* multiplicative -- 1.0 = no FRAG/SHAT loss */
#endif

            int64_t start = offsets[aa], end = offsets[aa + 1];
            for(int64_t nn = start; nn < end; nn++) {
                int j = neighbors[nn];
                if((kp[j].Mass <= 0) || (kp[j].AGS_KernelRadius <= 0)) continue;

                /* Periodic-wrapped separation. */
                kernel.dp = local.Pos - kp[j].Pos;
                nearest_xyz(kernel.dp);
                double r2 = kernel.dp.norm_sq();
                if(r2 <= 0) continue;
                kernel.r = sqrt(r2);
                kernel.h_j = kp[j].AGS_KernelRadius;

                /* Pair overlap filter — mirrors CPU AGSForce_evaluate. */
#if defined(DM_SIDM)
                if(kernel.r > kernel.h_i + kernel.h_j) continue;
#else
                if(kernel.r > kernel.h_i && kernel.r > kernel.h_j) continue;
#endif
                kernel_hinv(kernel.h_j, &kernel.hinv_j, &kernel.hinv3_j, &kernel.hinv4_j);
                double u_i = kernel.r * kernel.hinv_i;
                double u_j = kernel.r * kernel.hinv_j;
                if(u_i < 1) kernel_main(u_i, kernel.hinv3_i, kernel.hinv4_i, &kernel.wk_i, &kernel.dwk_i, 0);
                else { kernel.wk_i = 0; kernel.dwk_i = 0; }
                if(u_j < 1) kernel_main(u_j, kernel.hinv3_j, kernel.hinv4_j, &kernel.wk_j, &kernel.dwk_j, 0);
                else { kernel.wk_j = 0; kernel.dwk_j = 0; }

                /* Atomic read of P[j].Vel before building dv (SIDM may modify it). */
                for(int k = 0; k < 3; k++) {
                    double Vel_j_k = Kokkos::atomic_load(&kp[j].Vel[k]);
                    kernel.dv[k] = local.Vel[k] - Vel_j_k;
                    if(All.ComovingIntegrationOn) {
                        kernel.dv[k] += All.cf_hubble_a * kernel.dp[k] / All.cf_a2inv;
                    }
                }

#if defined(CBE_INTEGRATOR)
                {
                    CbeFluxResult cbe_r = cbe_integrator_flux_compute_pair(local, j, kp, kernel, out, tba_cap.v);
                    if(cbe_r.set_wakeup_j) {
                        /* Hydro-convention wakeup (was -1; MAX-reverse-comm silently dropped). */
                        short int wakeup_val = (short int)(local_TimeBin + 1);
                        Kokkos::atomic_max(&kp[j].wakeup, wakeup_val);
                        Kokkos::atomic_store(need_wakeup, 1);
                    }
                }
#endif

#if defined(DM_FUZZY)
                dm_fuzzy_flux_compute_pair(local, j, kp, kernel, out);
#endif

#if defined(DM_SIDM)
                {
                    SidmScatterResult sidm_r = sidm_core_flux_compute_pair(local, j, kp, kernel, out, geofactor, tba_cap.v);
                    if(sidm_r.scattered) {
                        if(sidm_r.set_wakeup_j) {
                            /* Hydro-convention wakeup (was -1; MAX-reverse-comm silently dropped). */
                            short int wakeup_val = (short int)(local_TimeBin + 1);
                            Kokkos::atomic_max(&kp[j].wakeup, wakeup_val);
                            Kokkos::atomic_store(need_wakeup, 1);
                        }
                        for(int kv = 0; kv < 3; kv++) {
                            Kokkos::atomic_add(&kp[j].Vel[kv], sidm_r.dv_sidm[kv]);
                            Kokkos::atomic_add(&kp[j].dp[kv],  sidm_r.dv_sidm[kv] * kp[j].Mass);
                        }
                        Kokkos::atomic_add(&kp[j].NInteractions, (long unsigned int)1);
                    }
                }
#endif
            } /* end neighbor loop */

            kout[aa] = out;
        });
    }

    /* Row 6b of arena-scope sweep (kernel-writes audit a06e30ca, updated
     * post-047ed629 to reflect the CBE+hydro-convention wakeup fix):
     *   j-side kernel writes:
     *     - CBE_INTEGRATOR: P_gpu[j].wakeup            atomic_max(TimeBin+1)
     *     - DM_SIDM:        P_gpu[j].wakeup            atomic_max(TimeBin+1)
     *                       P_gpu[j].Vel[kv]           atomic_add
     *                       P_gpu[j].dp[kv]            atomic_add
     *                       P_gpu[j].NInteractions     atomic_add(1)
     *     - Without DM_SIDM and without CBE_INTEGRATOR: NO j-side writes
     *       (kernel writes only per-active kout[aa]).
     *   The former full memcpy(P_host, P_gpu, num_total*...) was cargo-cult
     *   and is DELETED.
     *
     *   Host visibility of j-side writes:
     *     - Under UVM-canonical particles (current contract, commit 0d9e74b4):
     *       P_host[j] and P_gpu[j] are aliases of the same SharedSpace
     *       allocation. Device atomic writes ARE the host writes; no scatter
     *       needed for ANY j-side fields. The DM_SIDM-gated sparse scatter
     *       below is therefore a no-op under UVM (memcpy of identical
     *       pointers); it is preserved as a documented fallback for any
     *       future discrete-memory build path.
     *     - If a non-UVM mode is ever resurrected, the DM_SIDM-gated scatter
     *       would NOT cover CBE-only j-side wakeup writes — the gate would
     *       need to be widened to `#if defined(DM_SIDM) || defined(CBE_INTEGRATOR)`.
     *       Flagged here so the audit is explicit, not implicit.
     */
    memcpy(out_host, d_out, num_active * sizeof(struct ags_force_gpu_out));
#if defined(DM_SIDM)
    if(gnl.total_pairs > 0) {
        std::vector<int> gnl_neighbors_host((size_t)gnl.total_pairs);
        gpu_ngb_copy_neighbors_to_host(&gnl, gnl_neighbors_host.data());
        for(int64_t idx = 0; idx < gnl.total_pairs; idx++) {
            int j = gnl_neighbors_host[idx];
            P_host[j].wakeup        = P_gpu[j].wakeup;
            P_host[j].Vel           = P_gpu[j].Vel;
            P_host[j].dp            = P_gpu[j].dp;
            P_host[j].NInteractions = P_gpu[j].NInteractions;
        }
    }
#endif
    if(*d_need_wakeup) NeedToWakeupParticles_local = 1;

    /* Full P scatter syncs arena→host (P-side); CellP unmodified by this kernel. */
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_need_wakeup);
#if defined(DM_SIDM)
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_geofactor);
#endif
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_out);
    Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(d_radii);
    gpu_ngb_list_free(&gnl, NULL);
}



#else /* stubs when disabled */

void ags_force_evaluate_gpu(struct particle_data *, int, int *, int,
                            const double *, int, struct ags_force_gpu_out *) {}

#endif /* AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE */
