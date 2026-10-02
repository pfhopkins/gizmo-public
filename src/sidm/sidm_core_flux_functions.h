/* sidm_core_flux_functions.h -- per-pair SIDM scattering decision + kick.
 *
 * KOKKOS_INLINE_FUNCTION; guarded by DM_SIDM so the sole caller
 * (AgsForceSpec::pair_kernel in gravity/ags_force_loop.h) can invoke
 * unconditionally.
 *
 * SIDM scatter is a genuine two-sided physics event: when a scatter occurs,
 * both i and j receive momentum kicks. The i-side update goes into the caller
 * out struct (merged via AccumData); the j-side delta is returned in
 * SidmScatterResult and applied by the caller via Kokkos::atomic_add on
 * P[j].Vel / P[j].dp (and atomic_max on P[j].wakeup).
 *
 * RNG: counter-based gizmo_gpu_rand_double keyed symmetrically on the pair -- the identifiers XORed
 * plus a mix of each end's position, since identifiers repeat and a same-ID pair would otherwise key
 * every such pair to zero -- with a per-loop salt XOR-mixed into the counter (see gpu_rng.h
 * gizmo_loop_rng_salt). Identical statistics and identical streams on CPU
 * and GPU; no GSL.
 *
 * Helpers (prob_of_interaction_tab, prob_of_grain_interaction_tab,
 * calculate_interact_kick_rng) live in sidm/sidm_helper_functions.h and
 * solids/grain_helper_functions.h.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef SIDM_CORE_FLUX_FUNCTIONS_H
#define SIDM_CORE_FLUX_FUNCTIONS_H

#include "../declarations/allvars.h"
#include "../declarations/gpu_rng.h"
#include "sidm_helper_functions.h"
#ifdef GRAIN_COLLISIONS
#include "../solids/grain_helper_functions.h"
#endif
#ifdef GRAIN_EVOLUTION
#include "../solids/grain_evolution_functions.h"
#endif

struct SidmScatterResult {
    int scattered;          /* 1 if a scatter occurred this pair */
    Vec3<double> dv_sidm;   /* to add to P[j].Vel[k] (atomic); P[j].dp[k] += dv_sidm[k]*P[j].Mass */
    int set_wakeup_j;       /* 1 if P[j].wakeup should be set to -1 */
};

/* `rng_salt` is REQUIRED — a per-loop FNV-1a hash from gizmo_loop_rng_salt()
 * (see declarations/gpu_rng.h). It is XOR-mixed into the counter so two
 * different loops drawing SIDM scatter on the same (Ti_Current, ID_i^ID_j)
 * pair pull independent streams. No default value — every caller must
 * supply a named loop salt. Today's sole caller is the ags_force runner
 * (AGS_FORCE_RNG_SALT in gravity/ags_force_loop.h). */
template <typename LocalT, typename KernelT, typename OutT>
KOKKOS_INLINE_FUNCTION
SidmScatterResult sidm_core_flux_compute_pair(
    const LocalT &local,
    int j,
    struct particle_data *P,
    const KernelT &kernel,
    OutT &out,
    const int *timebin_active,
    uint64_t rng_salt)
{
    SidmScatterResult r;
    r.scattered = 0;
    r.dv_sidm = {0, 0, 0};
    r.set_wakeup_j = 0;
#ifdef DM_SIDM
    double Pj_dtime = get_particle_timestep_in_physical(j, P);
    /* Identity is a separation of zero, not an identifier: particle IDs are not unique -- every
     * spawned wind cell carries one stamped ID and an input can hold duplicates of its own -- so
     * testing them here excluded genuine distinct partners that happened to share an ID, and those
     * pairs silently never scattered.  A pair at zero separation is either the same particle or a
     * degenerate one the pair terms cannot evaluate anyway. */
    if(!( ((1 << local.Type) & (DM_SIDM)) && ((1 << P[j].Type) & (DM_SIDM))
          && (kernel.r > 0) && (local.dtime <= Pj_dtime) )) {
        return r;
    }
    /* Ensure each pair is computed only once.  The separation is antisymmetric, so the first axis on
     * which it is nonzero is positive for exactly one of the two evaluations -- which is the whole
     * requirement.  Ordering on the identifier cannot do this when identifiers repeat: it returned
     * for neither side, and the pair scattered twice. All three axes zero is a separation of zero,
     * already rejected above. */
    if(local.dtime == Pj_dtime) {
        if(kernel.dp[0] != 0)      {if(kernel.dp[0] > 0) {return r;}}
        else if(kernel.dp[1] != 0) {if(kernel.dp[1] > 0) {return r;}}
        else if(kernel.dp[2] > 0)  {return r;}
    }

    /* Pairwise-mean mass: used by the momentum-conserving kick split below.
       The scattering probability no longer needs it -- each side of the pair
       now carries its own mass and kernel. */
    double m_si = 0.5 * (local.Mass + P[j].Mass);
    Vec3<double> dv_local = kernel.dv;
#ifdef GRAIN_COLLISIONS
    double prob = prob_of_grain_interaction(local.Mass, local.Grain_Size, kernel.r, kernel.h_i, kernel.h_j, dv_local, local.dtime, j, P);
#else
    double prob = prob_of_interaction(local.Mass, P[j].Mass, kernel.r, kernel.h_i, kernel.h_j, dv_local, local.dtime);
#endif
    if(prob > 0.2) { out.dtime_sidm = DMIN(out.dtime_sidm, local.dtime * (0.2 / prob)); }

    /* Counter-based RNG, symmetric (i,j) key: both sides of the pair see the
       same stream. Timestep, loop-domain salt, and draw tag are mixed into
       one counter; the tag only needs to distinguish the two SIDM draws
       (0 = threshold, 1 = scatter direction). Ti_Current is left-shifted
       to keep timestep in the high half of the counter; the salt and tag
       occupy the rest. */
    /* The pair key must be SYMMETRIC (both evaluations of a pair draw the same stream) and DISTINCT
     * between pairs.  Identifiers alone give neither once they repeat: a same-ID pair XORs to zero,
     * so every such pair in the run shared one stream and made the identical scatter decision in the
     * identical direction.  Adding a mix of each end's position restores distinctness -- positions
     * are unique at any one moment, a pair at zero separation is already rejected above -- while a
     * SUM over the two ends keeps it symmetric under exchange. */
    uint64_t mix_i = gizmo_position_mix(local.Pos[0], local.Pos[1], local.Pos[2]);
    uint64_t mix_j = gizmo_position_mix(P[j].Pos[0],   P[j].Pos[1],   P[j].Pos[2]);
    uint64_t mix_lo = (mix_i < mix_j) ? mix_i : mix_j;   /* order the two ends so the key is the */
    uint64_t mix_hi = (mix_i < mix_j) ? mix_j : mix_i;   /* same from either side of the pair     */
    uint64_t rng_key = ((uint64_t)local.ID ^ (uint64_t)P[j].ID) + mix_lo * 0x9E3779B97F4A7C15ULL;
    rng_key ^= (mix_hi + 0xBF58476D1CE4E5B9ULL + (rng_key << 6) + (rng_key >> 2));
    rng_key ^= rng_key >> 31;
    uint64_t ti      = ((uint64_t)All.Ti_Current) << 32;
    uint64_t rng_ctr_threshold = ti ^ rng_salt ^ UINT64_C(0);
    uint64_t rng_ctr_direction = ti ^ rng_salt ^ UINT64_C(1);
    if(gizmo_gpu_rand_double(rng_key, rng_ctr_threshold) >= prob) { return r; }

    /* scatter happens */
    r.scattered = 1;
    if(!(timebin_active[P[j].TimeBin])) {
        if(WAKEUP * local.dtime < Pj_dtime) { r.set_wakeup_j = 1; }
    }
    Vec3<double> kick;
    calculate_interact_kick_rng(dv_local, kick, m_si, rng_key, rng_ctr_direction);
    for(int k=0; k<3; k++) {
        double dv_sidm = (local.Mass / m_si) * kick[k];
        out.sidm_kick[k] -= (P[j].Mass / m_si) * kick[k];
        r.dv_sidm[k] = dv_sidm;
    }
    out.si_count++;

#if defined(GRAIN_EVOLUTION) && (GRAIN_EVOLUTION & 7)
    /* Pairwise grain-evolution outcomes (COAG/FRAG/SHAT). The SIDM scatter
     * machinery already handles "does a collision fire" via the prob
     * gate; this resolver decides what HAPPENS to grain mass + size when
     * a collision fires, by comparing |dv| to per-species thresholds. */
    grain_evolution_resolve_pairwise(local, j, P, dv_local, out);
#endif
#else
    (void)local; (void)j; (void)P; (void)kernel; (void)out;
    (void)timebin_active; (void)rng_salt;
#endif /* DM_SIDM */
    return r;
}

#endif /* SIDM_CORE_FLUX_FUNCTIONS_H */
