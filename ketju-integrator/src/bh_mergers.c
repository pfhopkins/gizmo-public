// Ketju-integrator (MSTAR)
// Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

// The merger properties function is made available to external code for
// testing, but is not really a part of the public interface. For testing
// purposes the input is also passed as a struct instead of reading from the
// ketju_system fields directly
struct bh_merger_input {
    double G, c;
    double m1, m2, x1[3], x2[3], v1[3], v2[3], S1[3], S2[3];
};

struct bh_merger_remnant_properties {
    double m, vkick[3], S[3];
};

// Compute the fitted merger remnant properties based on
// Zlochower and Lousto 2015, PhRvD 92, 024022 (arXiv:1503.07536)
// The published version of the paper has been corrected with an additional
// erratum, and the corrected values are also in the arXiv version.
// It's unclear how accurate the results are for generic spinning binaries,
// and the actual resulting kick values appear to depend sensitively on the
// arbitrary point where the merger is triggered, probably due to the dependence
// on the inplane total spin and difference in spins. The values used in the
// fits are measured just before merger, whereas we need to use the values at
// the separation where the merger is triggered.
// The model also assumes a quasi-circular merger, so in the rare case that the
// orbit is still eccentric at merger the results will be incorrect by
// <~25% for the kick velocity (e.g. Radia et al. 2021 arxiv:2101.11015).
//
// Note that we don't have access to the exact rest frame of the binary,
// as the PN CoM velocity formulae implemented in pn.c are only accurate to
// a few tens of km/s in the worst cases, depending on the configuration of the
// binary (mostly spins).
// As a result, the kick that is calculated by this function and reported in the
// merger remnant properties is not the actual velocity the remnant gets
// relative to the actual binary rest frame, but the actual velocity is offset
// by a small error that can't be determined from the instantaneous
// configuration.
void ketju_calculate_bh_merger_remnant_properties(
    struct bh_merger_input *in, struct bh_merger_remnant_properties *out) {

    const double G = in->G, c = in->c;

    const double m1 = in->m1, m2 = in->m2;
    const double m = m1 + m2;

    // relative position and velocity
    double rv[3], vv[3];
    for (int i = 0; i < 3; i++) {
        rv[i] = in->x1[i] - in->x2[i];
        vv[i] = in->v1[i] - in->v2[i];
    }

    // orbital angular momentum (Newtonian) and corresponding unit vector
    double L[3], Lhat[3];
    vector_cross(rv, vv, L);
    for (int i = 0; i < 3; ++i) {
        L[i] *= m1 * m2 / m;
    }

    // Check that L is non-zero to handle the corner case of a head-on merger
    // separately. This shouldn't really ever happen, but don't want it to
    // crash then anyway.
    if (L[0] == 0 && L[1] == 0 && L[2] == 0) {
        // In this case there won't be significant GW emission,
        // so just use Newtonian values and assume spin conservation.
        out->m = m;
        for (int i = 0; i < 3; ++i) {
            out->vkick[i] = 0;
            out->S[i] = in->S1[i] + in->S2[i];
        }
        return;
    }

    unit_vector(L, Lhat);

    // set n0 equal to separation vector. vertical direction is given
    // by L_hat. right handed orthogonal basis given by n0, m0 and Lhat
    double n0[3], m59[3], m0[3];

    unit_vector(rv, n0);
    vector_cross(Lhat, n0, m0);

    // rotate n0 by -59 degrees to get m59
    {
        const double angle = -59 * M_PI / 180.0;

        for (int i = 0; i < 3; i++) {
            m59[i] = n0[i] * cos(angle) + m0[i] * sin(angle);
        }
    }

    // we need direction of total angular momentum for the direction of remnant
    // spin
    double J[3], Jhat[3];
    for (int i = 0; i < 3; i++) {
        J[i] = in->S1[i] + in->S2[i] + L[i];
    }
    unit_vector(J, Jhat);

    // the paper uses spins in units of M^2, multiply G and c out
    double S1[3], S2[3];
    for (int i = 0; i < 3; i++) {
        S1[i] = in->S1[i] * c / G;
        S2[i] = in->S2[i] * c / G;
    }

    // need a bunch of derived quantities
    double eta, delta_m, Delta_tilde_para, Delta_tilde_ortho, S_tilde_para,
        S_tilde_ortho, S0_tilde_para;
    double Delta_tilde[3], S_tilde[3], S0_tilde[3];

    delta_m = (m1 - m2) / m;
    eta = (1 - delta_m * delta_m) / 4;
    for (int i = 0; i < 3; i++) {
        S_tilde[i] = (S1[i] + S2[i]) / (m * m);
        Delta_tilde[i] = (S2[i] / m2 - S1[i] / m1) / m; // Delta / m^2
        S0_tilde[i] = S_tilde[i] + 0.5 * delta_m * Delta_tilde[i]; // S0 / m^2
    }

    // quantities are decomposed into para (parallel to orbital angular
    // momentum) and ortho components
    Delta_tilde_para = vector_dot(Delta_tilde, Lhat);
    {
        double Deltax = vector_dot(Delta_tilde, n0);
        double Deltay = vector_dot(Delta_tilde, m0);
        Delta_tilde_ortho = hypot(Deltax, Deltay);
    }

    S_tilde_para = vector_dot(S_tilde, Lhat);
    {
        double Sx = vector_dot(S_tilde, n0);
        double Sy = vector_dot(S_tilde, m0);
        S_tilde_ortho = hypot(Sx, Sy);
    }

    S0_tilde_para = vector_dot(S0_tilde, Lhat);

    // Calculate the components of recoil velocity, in km/s
    // Eq. (37)
    // Note that this equation appears to have a typo in the leading eta term
    // in the paper. There it reads 4*pow(eta,2), but using that form gives
    // far too low maximum kick values compared to the simulation results
    // in e.g. table XI. All the other formulae use pow(4*eta,2) as well,
    // including the closely related eq 19, so it seems that this is the correct
    // form. The maximal kick values from table XI are also reproduced within
    // the quoted tolerance using this formula.
    const double v_para =
        pow(4 * eta, 2)
        * (vector_dot(Delta_tilde, n0)
               * (3678. - 2475. * pow(delta_m, 2) + 4962. * S0_tilde_para
                  + 7170. * pow(S0_tilde_para, 2)
                  + 12050. * pow(S0_tilde_para, 3))
           + vector_dot(S0_tilde, m59)
                 * (Delta_tilde_para
                        * (4315. - 1262. * pow(delta_m, 2)
                           + 15970. * S0_tilde_para)
                    - 2256. * delta_m - 2231. * delta_m * S0_tilde_para));

    // Eq. (39)
    double v_ortho2 =
        pow(4 * eta, 4)
            * (2.106e5 * pow(Delta_tilde_para, 2)
               + 4.967e5 * Delta_tilde_para * delta_m
               - 2.116e5 * pow(Delta_tilde_para, 3) * delta_m
               - 5.037e5 * pow(Delta_tilde_para, 2) * pow(delta_m, 2)
               - 1.269e5 * Delta_tilde_para * pow(delta_m, 3)
               - 3.384e5 * pow(Delta_tilde_para, 2) * S0_tilde_para
               - 6.440e5 * pow(delta_m, 2) * S0_tilde_para
               + 2.138e6 * pow(Delta_tilde_para, 2) * pow(S0_tilde_para, 2)
               - 4.905e6 * Delta_tilde_para * delta_m * pow(S0_tilde_para, 2)
               - 1.100e6 * pow(delta_m, 2) * pow(S0_tilde_para, 2)
               - 1.024e7 * pow(Delta_tilde_para, 2) * pow(S0_tilde_para, 3))
        + pow(1.2e4 * eta * eta * delta_m * (1 - 0.93 * eta), 2);
    if (v_ortho2 < 0) { v_ortho2 = 0; }
    const double v_ortho = sqrt(v_ortho2);

    // Conversion between the internal unit system defined by the given value of
    // c and the SI units. The model gives kick velocities in km/s.
    const double km_per_s = c / 2.99792458e5;

    // XXX: we will set the parallel component along the direction of
    // original r, since the PN dynamics can't follow the binary until the
    // initial separation of the simulations.
    // This makes the direction essentially random.
    for (int i = 0; i < 3; i++) {
        out->vkick[i] = km_per_s * (v_para * Lhat[i] + v_ortho * n0[i]);
    }

    // For mass loss and spin, we need (dimensionless) energy and angular
    // momentum evaluated at ISCO of Kerr metric, with spin = S_tilde_para.
    // Formulas are standard, mostly given after Eq. (35) in the corrected arXiv
    // version.
    const double chi = S_tilde_para;
    const double sign_chi = chi == 0 ? 0 : (chi > 0 ? 1 : -1);
    const double Z1 = 1
                      + pow(1 - chi * chi, 1. / 3)
                            * (pow(1 + chi, 1. / 3) + pow(1 - chi, 1. / 3));
    const double Z2 = sqrt(3 * chi * chi + Z1 * Z1);
    // also dimensionless / in units of system mass
    const double r_isco =
        3 + Z2 - sign_chi * sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
    const double L_isco = 2. / sqrt(3 * r_isco) * (3 * sqrt(r_isco) - 2 * chi);
    // E_isco not given in the paper, but is a standard result.
    const double E_isco = sqrt(1 - 2. / (3 * r_isco));

    // Calculate mass loss.
    // Only use the E_c term in Eq. (30), the other term should be small and
    // can't be modelled so that it would improve accuracy.
    {
        const double k2a = -0.024;
        const double k2d = 0;
        const double k3a = -0.055;
        const double k3b = 0.0;
        const double k3d = -0.019;
        const double k4a = -0.119;
        const double k4b = 0.005;
        const double k4e = 0;
        const double k4f = 0.035;
        const double k4g = 0.022;
        const double E_A = 0;
        const double E_B = 0.59;
        const double E_D = 0;
        const double E_E = -0.51;
        const double E_F = 0.056;
        const double E_G = -0.073;
        const double E_H = 0;
        const double e1 = 0.0356;
        const double e2 = 0.096;
        const double e3 = 0.122;
        const double eps1 = 0.0043;
        const double eps2 = 0.0050;
        const double eps3 = -0.009;

        // Eq. (33)
        const double E_HU =
            0.0025829 - 0.0773079 / (2 * S0_tilde_para - 1.693959);
        // Eq. (32)
        const double Ec_para =
            pow(4 * eta, 2)
                * (E_HU + k2a * delta_m * Delta_tilde_para
                   + 0.000743 * pow(Delta_tilde_para, 2) + k2d * pow(delta_m, 2)
                   + k3a * delta_m * Delta_tilde_para * S_tilde_para
                   + k3b * S_tilde_para * pow(Delta_tilde_para, 2)
                   + k3d * pow(delta_m, 2) * S_tilde_para
                   + k4a * delta_m * Delta_tilde_para * pow(S_tilde_para, 2)
                   + k4b * delta_m * pow(Delta_tilde_para, 3)
                   + 0.000124 * pow(Delta_tilde_para, 4)
                   + k4e * pow(Delta_tilde_para, 2) * pow(S_tilde_para, 2)
                   + k4f * pow(delta_m, 4)
                   + k4g * pow(delta_m, 3) * Delta_tilde_para)
            + pow(delta_m, 6) * eta * (1 - E_isco);

        // Eq. (34)
        const double Ec =
            Ec_para
            + pow(4 * eta, 2)
                  * (pow(S_tilde_ortho, 2)
                         * (e1 + e2 * S_tilde_para + e3 * pow(S_tilde_para, 2))
                     + pow(Delta_tilde_ortho, 2)
                           * (eps1 + eps2 * S_tilde_para
                              + eps3 * pow(S_tilde_para, 2))
                     + pow(delta_m, 2) * pow(S_tilde_ortho, 2)
                           * (E_A + S_tilde_para * E_B)
                     + pow(delta_m, 2) * pow(Delta_tilde_ortho, 2)
                           * (E_D + S_tilde_para * E_E)
                     + E_F * delta_m * Delta_tilde_ortho * S_tilde_ortho
                     + E_G * pow(Delta_tilde_para, 2)
                           * pow(Delta_tilde_ortho, 2)
                     + E_H * pow(Delta_tilde_para, 2) * pow(S_tilde_ortho, 2));

        // take delta M = (M1 + M2 - M_remnant)/(M1 + M2) = Ec
        out->m = m * (1 - Ec);
    }

    // Calculate square of remnant spin magnitude, using only the A_c term in
    // Eq. (31). Spin direction is  assumed to be equal to total angular
    // momentum, which the paper assures is accurate to within 20 degrees.
    {

        const double L0 = 0.686710;
        const double L1 = 0.613247;
        const double L2a = -0.145427;
        const double L2b = -0.115689;
        const double L2c = -0.005254;
        const double L2d = 0.801838;
        const double L3a = -0.073839;
        const double L3b = 0.004759;
        const double L3c = -0.078377;
        const double L3d = 1.585809;
        const double L4a = -0.003050;
        const double L4b = -0.002968;
        const double L4c = 0.004364;
        const double L4d = -0.047204;
        const double L4e = -0.053099;
        const double L4f = 0.953458;
        const double L4g = -0.067998;
        const double L4h = 0.001629;
        const double L4i = -0.066693;

        // The A_* coefficients have large uncertainties, so no need to bother
        // with the last decimals.
        const double A_A = 3.0;
        const double A_B = -7.3;
        const double A_D = -2.0;
        const double A_E = 5.1;
        const double A_F = 0;
        const double A_G = -2.8;
        const double A_H = 5.1;
        const double a1 = 0.8401;
        const double a2 = -0.328;
        const double a3 = -0.61;
        const double zeta1 = -0.0209;
        const double zeta2 = -0.038;
        const double zeta3 = 0.04;

        // Eq. (35)
        const double alpha_align =
            pow(4 * eta, 2)
                * (L0 + L1 * S_tilde_para + L2a * delta_m * Delta_tilde_para
                   + L2b * pow(S_tilde_para, 2) + L2c * pow(Delta_tilde_para, 2)
                   + L2d * pow(delta_m, 2)
                   + L3a * Delta_tilde_para * S_tilde_para * delta_m
                   + L3b * S_tilde_para * pow(Delta_tilde_para, 2)
                   + L3c * pow(S_tilde_para, 3)
                   + L3d * S_tilde_para * pow(delta_m, 2)
                   + L4a * Delta_tilde_para * pow(S_tilde_para, 2) * delta_m
                   + L4b * pow(Delta_tilde_para, 3) * delta_m
                   + L4c * pow(Delta_tilde_para, 4) + L4d * pow(S_tilde_para, 4)
                   + L4e * pow(Delta_tilde_para, 2) * pow(S_tilde_para, 2)
                   + L4f * pow(delta_m, 4)
                   + L4g * Delta_tilde_para * pow(delta_m, 3)
                   + L4h * pow(Delta_tilde_para * delta_m, 2)
                   + L4i * pow(S_tilde_para * delta_m, 2))
            + S_tilde_para * (1 + 8 * eta) * pow(delta_m, 4)
            + eta * L_isco * pow(delta_m, 6);

        // Eq. (36)
        double Ac =
            pow(alpha_align, 2)
            + pow(4 * eta, 2)
                  * (pow(S_tilde_ortho, 2)
                         * (a1 + a2 * S_tilde_para + a3 * pow(S_tilde_para, 2))
                     + pow(Delta_tilde_ortho, 2)
                           * (zeta1 + zeta2 * S_tilde_para
                              + zeta3 * pow(S_tilde_para, 2))
                     + pow(delta_m, 2) * pow(S_tilde_ortho, 2)
                           * (A_A + A_B * S_tilde_para)
                     + pow(delta_m, 2) * pow(Delta_tilde_ortho, 2)
                           * (A_D + A_E * S_tilde_para)
                     + A_F * delta_m * Delta_tilde_ortho * S_tilde_ortho
                     + A_G * pow(Delta_tilde_para, 2)
                           * pow(Delta_tilde_ortho, 2)
                     + A_H * pow(Delta_tilde_para, 2) * pow(S_tilde_ortho, 2))
            + pow(delta_m, 6) * (1 + 12 * eta) * pow(S_tilde_ortho, 2);

        // square of dimensionless spin magnitude now given by Ac. multiply in
        // G/c*remnant_mass^2 to get simulation units
        double spin_mag = sqrt(Ac);

        // Safety check if the fits don't extrapolate to extreme cases nicely.
        if (spin_mag > 1.) spin_mag = 1.;

        const double rm2 = out->m * out->m;
        for (int i = 0; i < 3; i++) {
            out->S[i] = (G / c) * rm2 * spin_mag * Jhat[i];
        }
    }
}

// Merge BH (PN particle) j into i
static void merge_bhs(struct ketju_system *R, int i, int j) {
    assert(j < R->num_pn_particles);
    assert(j != i);
    if (j < i) {
        int temp = i;
        i = j;
        j = temp;
    }
    struct ketju_system_physical_state *ps = R->physical_state;

    struct bh_merger_input bh_input = {
        .c = R->constants->c,
        .G = R->constants->G,
        .m1 = ps->mass[i],
        .m2 = ps->mass[j],
        .x1 = {ps->pos[i][0], ps->pos[i][1], ps->pos[i][2]},
        .x2 = {ps->pos[j][0], ps->pos[j][1], ps->pos[j][2]},
        .v1 = {ps->vel[i][0], ps->vel[i][1], ps->vel[i][2]},
        .v2 = {ps->vel[j][0], ps->vel[j][1], ps->vel[j][2]},
        .S1 = {ps->spin[i][0], ps->spin[i][1], ps->spin[i][2]},
        .S2 = {ps->spin[j][0], ps->spin[j][1], ps->spin[j][2]},
    };

    // CoM quantities.
    // Position using Newtonian formulae, since it's close enough (somewhere
    // between the BHs), but velocity needs PN corrections.
    // Errors in CoM position insignificant compared to typical
    // uncertainties/errors in velocity after just a few years even for the most
    // massive BHs.
    double xcom[3], vcom[3];
    for (int k = 0; k < 3; k++) {
        xcom[k] = (bh_input.m1 * bh_input.x1[k] + bh_input.m2 * bh_input.x2[k])
                  / (bh_input.m1 + bh_input.m2);
    }
    ketju_binary_PN_CoM_velocity(R, i, j, vcom);

    struct bh_merger_remnant_properties res;
    ketju_calculate_bh_merger_remnant_properties(&bh_input, &res);

    for (int k = 0; k < 3; ++k) {
        ps->pos[i][k] = xcom[k];
        ps->vel[i][k] =
            vcom[k] + (R->options->enable_bh_merger_kicks ? res.vkick[k] : 0);
        ps->spin[i][k] = res.S[k];
    }
    ps->mass[i] = res.m;

    // Store the info on the merger
    R->merger_infos =
        realloc(R->merger_infos,
                sizeof(struct ketju_bh_merger_info[R->num_mergers + 1]));
    struct ketju_bh_merger_info *info = &R->merger_infos[R->num_mergers];
    info->t_merger = ps->time;
    info->m1 = bh_input.m1;
    info->m2 = bh_input.m2;
    info->m_remnant = res.m;
    info->chi1 = vector_norm(bh_input.S1)
                 / (bh_input.G / bh_input.c * bh_input.m1 * bh_input.m1);
    info->chi2 = vector_norm(bh_input.S2)
                 / (bh_input.G / bh_input.c * bh_input.m2 * bh_input.m2);
    info->chi_remnant =
        vector_norm(res.S) / (bh_input.G / bh_input.c * res.m * res.m);
    info->v_kick = vector_norm(res.vkick);
    info->extra_index1 = i;
    info->extra_index2 = R->num_particles - 1;

    // Update extra_index1 for previous mergers, index2 won't change (see below)
    for (int k = 0; k < R->num_mergers; ++k) {
        struct ketju_bh_merger_info *old_info = &R->merger_infos[k];
        if (old_info->extra_index1 == j) {
            old_info->extra_index1 = info->extra_index2;
        } else if (old_info->extra_index1 > j
                   && old_info->extra_index1 < R->num_particles) {
            old_info->extra_index1 -= 1;
        }
    }

    R->num_mergers += 1;

    // Shift everything right of particle j left by one place to overwrite
    // j's data.
    memmove(&ps->mass[j], &ps->mass[j + 1],
            sizeof(double) * (R->num_particles - (j + 1)));
    memmove(&ps->pos[j], &ps->pos[j + 1],
            sizeof(double[3]) * (R->num_particles - (j + 1)));
    memmove(&ps->vel[j], &ps->vel[j + 1],
            sizeof(double[3]) * (R->num_particles - (j + 1)));
    memmove(&ps->spin[j], &ps->spin[j + 1],
            sizeof(double[3]) * (R->num_pn_particles - (j + 1)));

    // The extra data of the merged particle is moved to the end of the array
    // after the particles that still exist.
    // If there are multiple mergers, the removed particles will stay where they
    // are so that their merger_info extra_indices don't need to be updated.
    size_t elem_size = R->particle_extra_data_elem_size;
    char *extra_data = R->particle_extra_data;
    if (extra_data) {
        char rem_extra_data[elem_size];
        memmove(rem_extra_data, extra_data + j * elem_size, elem_size);
        memmove(extra_data + j * elem_size, extra_data + (j + 1) * elem_size,
                elem_size * (R->num_particles - (j + 1)));
        memmove(extra_data + (R->num_particles - 1) * elem_size, rem_extra_data,
                elem_size);
    }

    R->num_particles -= 1;
    R->num_pn_particles -= 1;
}

int ketju_do_bh_mergers(struct ketju_system *R) {
    int num_merged = 0;

    for (int i = 0; i < R->num_pn_particles; ++i) {
        for (int j = i + 1; j < R->num_pn_particles; ++j) {
            double r2 = 0;
            for (int k = 0; k < 3; ++k) {
                double dr =
                    R->physical_state->pos[i][k] - R->physical_state->pos[j][k];
                r2 += dr * dr;
            }
            // Simple condition: merge if particles within a user-specified
            // factor times their combined Schwarzschild radii.
            double fac = R->options->bh_merger_distance_in_Rs * 2
                         * R->constants->G
                         / (R->constants->c * R->constants->c);
            if (sqrt(r2) < fac
                               * (R->physical_state->mass[i]
                                  + R->physical_state->mass[j])) {
                merge_bhs(R, i, j);
                num_merged += 1;
            }
        }
    }

    return num_merged;
}
