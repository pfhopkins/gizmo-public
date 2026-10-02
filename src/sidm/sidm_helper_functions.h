/* sidm_helper_functions.h — GPU-callable SIDM helper math:
 *   - prob_of_interaction       per-pair SIDM scattering probability
 *   - calculate_interact_kick_rng  isotropic post-scatter kick using counter-RNG
 *
 * Each function is a pure KOKKOS_INLINE_FUNCTION: no internal references to
 * global state beyond All (available via the All_dev mirror on GPU), and no
 * GSL RNG.
 *
 * RNG migration: calculate_interact_kick used to call gsl_rng_uniform twice
 * (for cos_theta and phi). The new `_rng` variant takes a 64-bit key and a
 * 64-bit counter and uses the counter-based RNG in declarations/gpu_rng.h
 * (tag nibble 1 = SIDM scatter direction). The caller chooses key/counter to
 * guarantee (i,j)-symmetric, timestep-reproducible draws on both CPU and GPU.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef SIDM_HELPER_FUNCTIONS_H
#define SIDM_HELPER_FUNCTIONS_H

#include "../declarations/allvars.h"
#include "../mesh/kernel.h"          /* kernel_main */
#include "../declarations/gpu_rng.h"

/* KOKKOS_INLINE_FUNCTION falls back to plain inline outside a GPU TU. */
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif


#ifdef DM_SIDM

/* Pair scattering probability.  Vec3 dV is the i-minus-j velocity difference
   (code units, as the caller passes to the CPU tree-walk).

   A particle scatters off the mass its PARTNER represents, so one side's rate
   is m_j * W(r,h_i) times the cross-section per unit mass: the kernel share of
   the partner's density at this separation.  Summed over a particle's
   neighbours that returns the density it moves through, so the rate is
   normalised by construction rather than through an overlap integral.

   One event covers both macro-particles, and each stands for a micro-particle
   count proportional to its mass, so the pair rate is the count-weighted mean
   of the two one-sided rates.  That is the reduced mass m_i*m_j/(m_i+m_j)
   times the sum of the two kernels -- symmetric, so both members compute the
   same probability and agree on whether they scatter, and for equal masses it
   is the plain average.  Weighting by each side's OWN mass instead would halve
   the rate a light tracer sees moving through heavy neighbours.

   The kernel vanishes at r >= h, so the probability is non-zero exactly when
   r < max(h_i,h_j): the ordinary symmetric-search criterion, with no widened
   search radius and no separate overlap acceptance. */
KOKKOS_INLINE_FUNCTION
double prob_of_interaction(double m_i, double m_j, double r,
                           double h_i, double h_j,
                           const Vec3<double> &dV, double dt)
{
    double dVmag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]) / All.cf_atime;
    double hinv_i = 1.0 / h_i, hinv3_i = hinv_i*hinv_i*hinv_i;
    double hinv_j = 1.0 / h_j, hinv3_j = hinv_j*hinv_j*hinv_j;
    double wk_i = 0, wk_j = 0, dwk_dummy = 0;
    kernel_main(r * hinv_i, hinv3_i, hinv_i*hinv3_i, &wk_i, &dwk_dummy, -1);
    kernel_main(r * hinv_j, hinv3_j, hinv_j*hinv3_j, &wk_j, &dwk_dummy, -1);
    /* Comoving kernel (h is comoving), so the same a^-3 the density carries. */
    double m_sum = m_i + m_j;
    double m_red = (m_sum > 0) ? (m_i * m_j / m_sum) : 0;
    double rho_eff = m_red * (wk_i + wk_j) * All.cf_a3inv;
    double cx_eff = All.DM_InteractionCrossSection;
    double units = UNIT_SURFDEN_IN_CGS;
    if(All.DM_InteractionVelocityScale > 0) {
        double x = dVmag / All.DM_InteractionVelocityScale;
        cx_eff /= 1 + x*x*x*x;
    }
    return rho_eff * cx_eff * dVmag * dt * units;
}


/* Momentum-conserving isotropic scatter kick. Consumes two draws from the
   counter-based RNG, keyed by (rng_key, 2*rng_counter_base) and
   (rng_key, 2*rng_counter_base + 1). Writes to `kick`. `dV` is i-minus-j
   velocity difference (code units). `m_mean` is the pairwise-mean mass
   (unused in the kick formula itself — kept only so the signature matches
   the prior CPU prototype and makes the math path obvious to readers). */
KOKKOS_INLINE_FUNCTION
void calculate_interact_kick_rng(const Vec3<double> &dV, Vec3<double> &kick,
                                 double /*m_mean*/,
                                 uint64_t rng_key, uint64_t rng_counter_base)
{
    double dVmag = (1.0 - All.DM_DissipationFactor) *
                   sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    if(dVmag < 0) dVmag = 0;
    if(All.DM_KickPerCollision > 0) {
        double v0 = All.DM_KickPerCollision;
        dVmag = sqrt(dVmag*dVmag + v0*v0);
    }
    double u0 = gizmo_gpu_rand_double(rng_key, 2 * rng_counter_base);
    double u1 = gizmo_gpu_rand_double(rng_key, 2 * rng_counter_base + 1);
    double cos_theta = 2.0 * u0 - 1.0;
    double sin_theta_sq = 1.0 - cos_theta * cos_theta;
    double sin_theta = (sin_theta_sq > 0) ? sqrt(sin_theta_sq) : 0.0;
    double phi = u1 * 2.0 * M_PI;
    kick[0] = 0.5 * (dV[0] + dVmag * sin_theta * cos(phi));
    kick[1] = 0.5 * (dV[1] + dVmag * sin_theta * sin(phi));
    kick[2] = 0.5 * (dV[2] + dVmag * cos_theta);
}

#endif /* DM_SIDM */

#endif /* SIDM_HELPER_FUNCTIONS_H */
