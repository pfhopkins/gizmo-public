/* grain_helper_functions.h — GPU-callable helpers for GRAIN_COLLISIONS:
 *   - return_grain_cross_section_per_unit_mass_P
 *   - prob_of_grain_interaction_tab
 *
 * Pure KOKKOS_INLINE_FUNCTION forms that take particle_data and the
 * explicit arguments, so the same code path works on
 * CPU tree-walk and the AGSForce GPU kernel. The kernel is pulled in
 * from sidm/sidm_helper_functions.h.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef GRAIN_HELPER_FUNCTIONS_H
#define GRAIN_HELPER_FUNCTIONS_H

#include "../declarations/allvars.h"
#include "../mesh/kernel.h"          /* kernel_main */
#include "../sidm/sidm_helper_functions.h"

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif


#if defined(DM_SIDM) && defined(GRAIN_COLLISIONS)

/* Mass of ONE grain of the given radius -- the micro-particle a macro-particle
   stands for.  SSOT for the count N = M / grain_particle_mass(R). */
KOKKOS_INLINE_FUNCTION
double grain_particle_mass(double grain_radius)
{
    return (4.0*M_PI/3.0) * grain_radius*grain_radius*grain_radius * All.Grain_Internal_Density;
}

/* Per-particle cross-section per unit mass. Reads P[i].Grain_Size. */
KOKKOS_INLINE_FUNCTION
double return_grain_cross_section_per_unit_mass_P(int i, const struct particle_data *P_arr)
{
    return All.DM_InteractionCrossSection * 0.75 / (P_arr[i].Grain_Size * All.Grain_Internal_Density);
}


/* Pairwise grain-grain collision probability.

   Hard spheres: a grain of radius R_i meets a grain of radius R_j across the
   pair cross-section pi*(R_i+R_j)^2, scaled by the collision efficiency the
   cross-section parameter carries.  What differs from the dark-matter case is
   the counting: a macro-particle stands for M/grain_particle_mass(R) grains,
   so two macro-particles of the SAME mass hold very different numbers when
   their grains differ in size, and the small-grain side dominates.

   One event covers both macro-particles, so the rate is the count-weighted
   mean of the two one-sided rates, which is the reduced COUNT
   N_i*N_j/(N_i+N_j) -- written here over masses as M_i*M_j/(M_i*mu_j +
   M_j*mu_i) -- times the pair cross-section and the sum of the two kernels.
   It preserves the expected number of grain collisions the pair takes part in;
   it cannot also reproduce the two sides' separate rates, which would need a
   velocity distribution inside each macro-particle rather than one velocity.

   The kernels vanish at r >= h, so the support is the ordinary symmetric
   criterion r < max(h_i,h_j).

   NOTE this is four times the previous rate for equal grains: the old form
   used pi*R^2 where hard spheres of equal size present pi*(2R)^2. */
KOKKOS_INLINE_FUNCTION
double prob_of_grain_interaction(double m_i, double grain_radius_i,
                                 double r, double h_i, double h_j,
                                 const Vec3<double> &dV, double dt,
                                 int j_ngb,
                                 const struct particle_data *P_arr)
{
    double dVmag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]) / All.cf_atime;
    double grain_radius_j = (double)P_arr[j_ngb].Grain_Size;
    double mu_i = grain_particle_mass(grain_radius_i);
    double mu_j = grain_particle_mass(grain_radius_j);
    double m_j  = (double)P_arr[j_ngb].Mass;
    /* Reduced grain COUNT, over masses so no division by a zero grain mass. */
    double denom = m_i * mu_j + m_j * mu_i;
    if(!(denom > 0)) return 0;
    double n_red = m_i * m_j / denom;
    double r_sum = grain_radius_i + grain_radius_j;
    double sigma_ij = All.DM_InteractionCrossSection * M_PI * r_sum * r_sum;
    double hinv_i = 1.0 / h_i, hinv3_i = hinv_i*hinv_i*hinv_i;
    double hinv_j = 1.0 / h_j, hinv3_j = hinv_j*hinv_j*hinv_j;
    double wk_i = 0, wk_j = 0, dwk_dummy = 0;
    kernel_main(r * hinv_i, hinv3_i, hinv_i*hinv3_i, &wk_i, &dwk_dummy, -1);
    kernel_main(r * hinv_j, hinv3_j, hinv_j*hinv3_j, &wk_j, &dwk_dummy, -1);
    double rate = n_red * sigma_ij * (wk_i + wk_j) * All.cf_a3inv;
    double units = UNIT_SURFDEN_IN_CGS;
    if(All.DM_InteractionVelocityScale > 0) {
        double x = dVmag / All.DM_InteractionVelocityScale;
        rate /= 1 + x*x*x*x;
    }
    return rate * dVmag * dt * units;
}

#endif /* DM_SIDM && GRAIN_COLLISIONS */

#endif /* GRAIN_HELPER_FUNCTIONS_H */
