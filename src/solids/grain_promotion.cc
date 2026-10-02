/* grain_promotion.cc — grain super-particle → solid body promotion.
 * Called from core/run.cc alongside the other source terms. What it promotes is counted in
 * Grains_promoted and folded into the gas block by the next rearrange_particle_sequence(), which
 * is not part of this call.
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */

#include "../declarations/allvars.h"
#include "../core/wakeup_sidecar.h"
#include "grain_promotion.h"

#if defined(GRAIN_FLUID) && defined(GRAIN_FLUID_PROMOTION)

void grain_promotion_parent_routine(void)
{
    /* Iterate strictly over the active list (already TimeBinActive-filtered),
     * skipping gas particles. Never loop N_gas..NumPart in a per-step routine
     * — see Phil's full-loop rule. */
    for(int i : ActiveParticleList)
    {
        if(i < N_gas) continue;
        if(!((1 << P[i].Type) & GRAIN_PTYPES)) continue;

        /* Mass criterion (always active when threshold > 0) */
        if(P[i].Mass < All.GrainPromotion_MassThresh) continue;

        /* Dust/gas density ratio criterion (0 = disabled).
         * Grain local mass density proxy: NumNgb*Mass / (4π/3 * KernelRadius^3). */
        if(All.GrainPromotion_DustGasRatioThresh > 0.0) {
            double V_kernel = 4.18879020479 * P[i].KernelRadius * P[i].KernelRadius * P[i].KernelRadius;
            double grain_rho = (V_kernel > 0.0) ? ((double)P[i].NumNgb * P[i].Mass / V_kernel) : 0.0;
            double gas_rho   = (P[i].Gas_Density > 0.0) ? P[i].Gas_Density : 1.0e-30;
            if(grain_rho / gas_rho < All.GrainPromotion_DustGasRatioThresh) continue;
        }

        /* Safety: don't overflow gas block */
        if(N_gas + Grains_promoted >= All.MaxPartGas) {
            printf("Task %d: GRAIN_FLUID_PROMOTION safety limit reached (N_gas+Grains_promoted=%d >= MaxPartGas=%d), skipping grain %d\n",
                   ThisTask, N_gas + Grains_promoted, All.MaxPartGas, i);
            continue;
        }

        /* Initialize CellP BEFORE type change */
        grain_promotion_init_cellp(i);

        P[i].Type = 0;
        force_tree_note_type_presence(i);  /* the element's type is final here, so tell the standing
                                            * tree before anything walks it looking for that type */
        TimeBinCountGas[P[i].TimeBin]++;
        P[i].wakeup = -1; wakeup_sidecar_mark(i);
        NeedToWakeupParticles_local = 1;
        Grains_promoted++;
    }
}

#endif /* GRAIN_FLUID && GRAIN_FLUID_PROMOTION */
