/* transport_subcycle.cc - implementation of transport (RT + CR) subcycling
 *
 * Flux exchange uses the standard code_block_xchange pattern with a transport-only
 * core function. This correctly handles cross-domain interactions via MPI.
 * The tree walk is repeated each sub-step (future optimization: cache neighbor lists).
 *
 * Written by Mike Grudic with Claude for GIZMO.
 */

#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "../mesh/neighbor_list.h"
#include "../hydro/hydro_structs.h"
#include "../hydro/compute_finitevol_faces_functions.h"
extern void hydro_evaluate_gpu(struct particle_data *, struct gas_cell_data *,
                               int, int *, int, int64_t *, int *, int64_t, void *);

#ifdef TRANSPORT_SUBCYCLE

#include "transport_subcycle.h"
#include "../radiation/rt_direct_ray_transport_functions.h"
#include "../radiation/rt_diffusion_explicit_functions.h"

/* ======================================================================================== */
/* GPU dispatcher (transport_subcycle_exchange_fluxes below) consumes hydro_data_out
   directly. Legacy CPU-tree scaffolding (code_block_xchange_initialize.h, INPUT_STRUCT_NAME,
   OUTPUT_STRUCT_NAME, particle2in_transport, out2particle_transport, transport_flux_evaluate)
   was retired in Step 5 Phase D1. */


/* ===== INITIAL OPERATIONS (zero rate arrays) ===== */
void transport_flux_initial_operations(void)
{
    int k;
    for(int i : ActiveParticleList) {
        if(P[i].Type != 0 || P[i].Mass <= 0) continue;
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_EVOLVE_ENERGY)
        for(k=0; k<N_RT_FREQ_BINS; k++) CellP[i].Dt_Rad_E_gamma[k] = 0;
#endif
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_EVOLVE_FLUX)
        for(k=0; k<N_RT_FREQ_BINS; k++) CellP[i].Dt_Rad_Flux[k] = {};
#endif
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_INFRARED)
        CellP[i].Dt_Rad_E_gamma_T_weighted_IR = 0;
#endif
    }
}


/* ===== PUBLIC ENTRY POINT ===== */
void transport_subcycle_exchange_fluxes(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    transport_flux_initial_operations();

    /* GPU/neighbor-list path: reuse the hydro GPU kernel with the symmetric CSR list.
       The hydro kernel computes face geometry + RT fluxes (via rt_diffusion_explicit.h);
       we scatter only the RT output fields and discard the hydro force outputs.
       This eliminates the tree-walk + MPI exchange per subcycle step. */
    {
        struct hydro_data_out *hydro_out = (struct hydro_data_out *) mymalloc("transport_hydro_out",
            (gizmo_sym_num_active > 0 ? gizmo_sym_num_active : 1) * sizeof(struct hydro_data_out));

        hydro_evaluate_gpu(P, CellP, NumPart,
                           gizmo_sym_active_indices, gizmo_sym_num_active,
                           gizmo_sym_neighbor_list.offsets,
                           gizmo_sym_neighbor_list.neighbors,
                           gizmo_sym_neighbor_list.total_pairs,
                           (void *)hydro_out);

        /* Scatter only RT flux outputs (discard hydro forces) */
        for(int aa = 0; aa < gizmo_sym_num_active; aa++) {
            int ii = gizmo_sym_active_indices[aa];
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_EVOLVE_ENERGY)
            for(int kf=0; kf<N_RT_FREQ_BINS; kf++)
                CellP[ii].Dt_Rad_E_gamma[kf] += hydro_out[aa].Dt_Rad_E_gamma[kf];
#endif
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_EVOLVE_FLUX)
            for(int kf=0; kf<N_RT_FREQ_BINS; kf++)
                CellP[ii].Dt_Rad_Flux[kf] += hydro_out[aa].Dt_Rad_Flux[kf];
#endif
#if defined(RT_SOLVER_EXPLICIT) && defined(RT_INFRARED)
            CellP[ii].Dt_Rad_E_gamma_T_weighted_IR += hydro_out[aa].Dt_Rad_E_gamma_T_weighted_IR;
#endif
        }

        myfree(hydro_out);
    }

    /* add back the radiation pressure work terms saved during the hydro pass.
       The flux exchange zeroed Dt_Rad_E_gamma and recomputed only transport fluxes;
       the work terms must be re-added so each sub-step's kick includes them. */
#if defined(RT_RAD_PRESSURE_FORCES) && defined(RT_EVOLVE_ENERGY)
    for(int i : ActiveParticleList) {
        if(P[i].Type != 0 || P[i].Mass <= 0) continue;
        for(int kf=0; kf<N_RT_FREQ_BINS; kf++)
            CellP[i].Dt_Rad_E_gamma[kf] += CellP[i].Dt_Rad_E_gamma_Work[kf];
    }
#endif

    CPU_Step[CPU_RTNONFLUXOPS] += measure_time();
}


/* ======================================================================================== */
/* ===== KICK ============================================================================= */
/* ======================================================================================== */

void transport_subcycle_kick(void)
{
    double sub_frac = All.Transport_Subcycle_dt_fraction;
    for(int idx : ActiveParticleList) {
        if(P[idx].Type != 0 || P[idx].Mass <= 0) continue;
        int i = idx;
        double sub_dt = get_particle_timestep_in_physical(i, P) * sub_frac;
        if(sub_dt <= 0) continue;

#ifdef RADTRANSFER
        /* rt_update_driftkick adds IR gas heating to DtInternalEnergy as a rate.
           Since we call it N times, the rate would accumulate N-fold. Save/restore
           DtInternalEnergy each sub-step to prevent accumulation, but keep the final
           sub-step's contribution so the cooling function sees the correct IR heating rate. */
        rt_update_driftkick(i, sub_dt, 0, P, CellP); /* kick: advance conserved variables */
        rt_eddington_update_calculation(i, CellP);
        /* update opacities — Rad_Kappa (especially IR band) depends on Dust_Temperature and
           Radiation_Temperature which were just modified by rt_update_driftkick. Stale opacities
           cause the transport flux computation in subsequent sub-steps to use the wrong diffusion rate. */
        for(int kf = 0; kf < N_RT_FREQ_BINS; kf++) {CellP[i].Rad_Kappa[kf] = rt_kappa(i, kf, P, CellP);}
        /* sync Pred = conserved for the next sub-step's flux exchange */
        for(int kf = 0; kf < N_RT_FREQ_BINS; kf++) {
            CellP[i].Rad_E_gamma_Pred[kf] = CellP[i].Rad_E_gamma[kf];
#ifdef RT_EVOLVE_FLUX
            CellP[i].Rad_Flux_Pred[kf] = CellP[i].Rad_Flux[kf];
#endif
#ifdef RT_EVOLVE_INTENSITIES
            for(int ka = 0; ka < N_RT_INTENSITY_BINS; ka++)
                CellP[i].Rad_Intensity_Pred[kf][ka] = CellP[i].Rad_Intensity[kf][ka];
#endif
        }
#endif

#ifdef COSMIC_RAY_FLUID
        CosmicRay_Update_DriftKick(i, sub_dt, 0, P, CellP);
        for(int kb = 0; kb < N_CR_PARTICLE_BINS; kb++) {
            CellP[i].CosmicRayEnergyPred[kb] = CellP[i].CosmicRayEnergy[kb];
            CellP[i].CosmicRayFluxPred[kb] = CellP[i].CosmicRayFlux[kb];
        }
#endif
    }
}


#endif // TRANSPORT_SUBCYCLE
