// Ketju-integrator (MSTAR)
// Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.

// Functions related to the post-Newtonian corrections for BHs.

#include <string.h>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

inline static void particle_separation(const struct ketju_system *R, int i,
                                       int j, double rv[3]) {

    struct ketju_system_physical_state *ps = R->physical_state;
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    const int Nd = R->options->max_tree_distance;
    int path[Nd];
    int sign[Nd];
    int d;

    // Use MST or Cartesian coordinates depending on distance in the tree
    if (ketju_check_relative_proximity(i, j, Nd, R, &d, path, sign)) {
        for (int k = 0; k < 3; k++) {
            rv[k] = 0;
            for (int l = 0; l < d; l++) {
                int index = path[l];
                rv[k] += sign[l] * mst->edge_rel_pos[index][k];
            }
        }
    } else {
        for (int k = 0; k < 3; k++) {
            rv[k] = ps->pos[i][k] - ps->pos[j][k];
        }
    }
}

inline static void
particle_separation_and_relative_velocity(const struct ketju_system *R, int i,
                                          int j, double rv[3], double vv[3],
                                          double (*vel)[3]) {

    struct ketju_system_physical_state *ps = R->physical_state;
    const struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    const int Nd = R->options->max_tree_distance;
    int path[Nd];
    int sign[Nd];
    int d;

    if (ketju_check_relative_proximity(i, j, Nd, R, &d, path, sign)) {
        double(*rel_vel)[3] =
            vel == ps->vel ? mst->edge_rel_vel : mst->edge_rel_aux_vel;

        for (int k = 0; k < 3; k++) {
            rv[k] = 0;
            vv[k] = 0;
            for (int l = 0; l < d; l++) {
                int index = path[l];
                rv[k] += sign[l] * mst->edge_rel_pos[index][k];
                vv[k] += sign[l] * rel_vel[index][k];
            }
        }
    } else {
        for (int k = 0; k < 3; k++) {
            rv[k] = ps->pos[i][k] - ps->pos[j][k];
            vv[k] = vel[i][k] - vel[j][k];
        }
    }
}

// The vel and spin are taken as arguments to allow switching between physical
// and auxiliary values.
// Equations from Thorne & Hartle (1985), Phys. Rev. D, 31, 1815
// for the 1PN and spin parts, valid for N-body systems,
// and from Blanchet (2014), Living Rev. Relativity, 17, 2,
// doi:10.12942/lrr-2014-2, for the 2-3.5PN terms, which are only for binaries.
// The higher order terms therefore miss some cross terms in systems where there
// are more than two PN particles (black holes).
// The equations are in (modified) harmonic gauge, as was the case for the
// relative binary equations used in the original implementation (Mora & Will,
// 2004, Phys. Rev. D 69, 104021).
// Note that the Newtonian centre of mass is not conserved with these forces.
// Full three-body corrections are available in
// the literature to at least 2.5PN (e.g. Galaviz & Brügmann (2011)
// doi: 10.1103/PhysRevD.83.084013 and references therein) with also 3.5PN
// radiative corrections available. But these are in ADM hamiltonian formalism,
// which appears hard to apply with the current integrator design, as we need to
// work  with velocities as the primary variable instead of momenta which can't
// be used in a relative coordinate system and the relation between velocity and
// PN momentum is highly non-trivial.
static void compute_PN_acc_single(struct ketju_system *R, int i,
                                  double (*vel)[3], double (*spin)[3]) {

    struct ketju_system_physical_state *ps = R->physical_state;
    // these have been zeroed by the caller
    double *acc = R->internal_state->acc_vd[i];
    double *spin_derivative = R->internal_state->spin_derivative[i];

    const int N = R->num_pn_particles;

    const double G = R->constants->G;
    const double G2 = G * G;
    const double G3 = G2 * G;
    const double G4 = G2 * G2;
    const double c = R->constants->c;
    const double c2 = c * c;
    const double c4 = c2 * c2;
    const double c5 = c * c4;
    const double c6 = c4 * c2;
    const double c7 = c * c6;

    const double pi2 = M_PI * M_PI;

    const double mi = ps->mass[i];
    const double mi2 = mi * mi;
    const double Si2 = vector_dot(spin[i], spin[i]);
    const double vi2 = vector_dot(vel[i], vel[i]);

    // Compute the enabled terms
    const unsigned int flags = R->options->PN_flags;

    // The sum over pairs is common to all terms
    for (int j = 0; j < N; ++j) {
        if (j == i) continue;

        double rv_ij[3], nv_ij[3], vv_ij[3];
        particle_separation_and_relative_velocity(R, i, j, rv_ij, vv_ij, vel);
        const double rij = norm_and_unit_vector(rv_ij, nv_ij);
        const double rij2 = rij * rij;
        const double rij3 = rij2 * rij;
        const double rij4 = rij3 * rij;

        const double mj = ps->mass[j];
        const double mj2 = mj * mj;
        const double mj3 = mj2 * mj;

        const double vivj = vector_dot(vel[i], vel[j]);
        const double vj2 = vector_dot(vel[j], vel[j]);
        const double vj4 = vj2 * vj2;
        const double vij2 = vector_dot(vv_ij, vv_ij);

        const double vjnij = vector_dot(vel[j], nv_ij);
        const double vjnij2 = vjnij * vjnij;
        const double vjnij4 = vjnij2 * vjnij2;

        const double vinij = vector_dot(vel[i], nv_ij);
        const double vinij2 = vinij * vinij;

        const double vijnij = vector_dot(vv_ij, nv_ij);

        if (flags & KETJU_PN_1_0_ACC) {
            // Eq. 4.11a from Thorne & Hartle (1985) with the Newtonian part
            // removed and factors of c2, G added. The index substitutions are
            // K->i, A->j, B,C->k. Sign changes have been made to account for
            // the different directions of nv_ij vs n_AK etc.

            double nij_factor = vi2 + 2. * vj2 - 4. * vivj - 1.5 * vjnij2;

            if (flags & KETJU_PN_THREEBODY) {
                for (int k = 0; k < N; ++k) {
                    double mk = ps->mass[k];

                    if (k != i) {
                        double rv_ik[3];
                        particle_separation(R, i, k, rv_ik);
                        const double rik = vector_norm(rv_ik);

                        nij_factor -= 4. * G * mk / rik;
                    }

                    if (k != j) {
                        double rv_jk[3], nv_jk[3];
                        particle_separation(R, j, k, rv_jk);
                        double rjk = norm_and_unit_vector(rv_jk, nv_jk);

                        nij_factor -=
                            G * mk / rjk
                            * (1. - .5 * rij / rjk * vector_dot(nv_ij, nv_jk));

                        for (int l = 0; l < 3; ++l) {
                            acc[l] -= 3.5 * nv_jk[l] * G2 * mj * mk
                                      / (rij * rjk * rjk * c2);
                        }
                    }
                }
            } else {
                // The two-body contributions from the above expressions,
                // equivalent to the non-velocity dependent part in Eq. 4.18a.
                // The k=j part from above
                nij_factor -= 4. * G * mj / rij;
                // The k=i part from above
                nij_factor -= 5. * G * mi / rij;
            }

            nij_factor *= -G * mj / (rij2 * c2);

            // Here 4*vijnij + vjnij = 4*vinij - 3*vjnij
            const double vij_factor =
                G * mj / (rij2 * c2) * (4. * vijnij + vjnij);

            for (int l = 0; l < 3; ++l) {
                acc[l] += nij_factor * nv_ij[l] + vij_factor * vv_ij[l];
            }
        }

        // Spin terms
        if (flags & KETJU_PN_SPIN_ALL) {
            double nijSj = vector_dot(nv_ij, spin[j]);
            double nijSi = vector_dot(nv_ij, spin[i]);
            if (flags & KETJU_PN_1_5_SPIN_ACC) {
                // F_SO, eq. 4.11c in T&H
                double Sjxnij[3];
                vector_cross(spin[j], nv_ij, Sjxnij);
                double Sixnij[3];
                vector_cross(spin[i], nv_ij, Sixnij);
                double Sjxvij[3];
                vector_cross(spin[j], vv_ij, Sjxvij);
                double Sixvij[3];
                vector_cross(spin[i], vv_ij, Sixvij);
                double factor = G / (c2 * mi * rij3);
                for (int l = 0; l < 3; ++l) {
                    acc[l] +=
                        factor
                        * (6. * nv_ij[l]
                               * (mi * vector_dot(Sjxnij, vv_ij)
                                  + mj * vector_dot(Sixnij, vv_ij))
                           + 4. * mi * Sjxvij[l] + 3. * mj * Sixvij[l]
                           - vijnij
                                 * (6. * mi * Sjxnij[l] + 3. * mj * Sixnij[l]));
                }
            }

            if (flags & KETJU_PN_2_0_SPIN_ACC) {
                double SiSj = vector_dot(spin[i], spin[j]);
                double nijSi = vector_dot(nv_ij, spin[i]);
                double Sj2 = vector_dot(spin[j], spin[j]);
                double factor = G / (c2 * mi * rij4);
                for (int l = 0; l < 3; ++l) {
                    // F_SS (4.11d)
                    acc[l] +=
                        factor
                        * (nv_ij[l] * (-3. * SiSj + 15. * nijSi * nijSj)
                           - 3. * spin[i][l] * nijSj - 3. * spin[j][l] * nijSi);
                    // F_Q (4.11b) specialized for BHs using (1.7)
                    acc[l] +=
                        factor
                        * (mi / mj
                               * ((-1.5 * Sj2 + 7.5 * nijSj * nijSj) * nv_ij[l]
                                  - 3. * nijSj * spin[j][l])
                           + mj / mi
                                 * ((-1.5 * Si2 + 7.5 * nijSi * nijSi)
                                        * nv_ij[l]
                                    - 3. * nijSi * spin[i][l]));
                }
            }

            if (flags & KETJU_PN_SPIN_EVOL) {
                // first Omega_GM (4.17a)
                double buffer[3];
                double factor = G / (c2 * rij3);
                for (int l = 0; l < 3; ++l) {
                    buffer[l] = factor * (-spin[j][l] + 3. * nv_ij[l] * nijSj);
                }
                // then Omega_geod (4.17b)
                double v_diff[3];
                double v_diffxnij[3];
                for (int l = 0; l < 3; ++l) {
                    v_diff[l] = 2. * vel[j][l] - 1.5 * vel[i][l];
                }
                vector_cross(v_diff, nv_ij, v_diffxnij);
                factor = G * mj / (c2 * rij2);
                for (int l = 0; l < 3; ++l) {
                    buffer[l] += factor * v_diffxnij[l];
                }

                // Then N_Q (4.17c) using (1.7) and factoring out the cross
                // product
                factor = 3 * G * mj / (rij3 * mi * c2) * nijSi;
                for (int l = 0; l < 3; ++l) {
                    buffer[l] += factor * nv_ij[l];
                }

                // Finally the cross product with the spin
                double sder[3];
                vector_cross(buffer, spin[i], sder);
                for (int l = 0; l < 3; ++l) {
                    spin_derivative[l] += sder[l];
                }
            }
        }

        // The rest is from Blanchet (2014) eq. 203
        // Everything is decomposed into components along n_ij, v_ij, so
        // calculate the total factors first
        double nij_factor = 0.;
        double vij_factor = 0.;
        if (flags & KETJU_PN_2_0_ACC) {
            nij_factor +=
                (-57. * G3 * mi2 * mj / (4. * rij4)
                 - 69. * G3 * mi * mj2 / (2. * rij4) - 9. * G3 * mj3 / rij4
                 + G * mj / rij2
                       * (vjnij2
                              * (-15. / 8. * vjnij2 + 1.5 * vi2 - 6. * vivj
                                 + 9. / 2. * vj2)
                          - 2. * vivj * vivj + 4. * vivj * vj2 - 2. * vj4)
                 + G2 * mi * mj / rij3
                       * (39. / 2. * vinij2 - 39. * vinij * vjnij
                          + 17. / 2. * vjnij2 - 15. / 4. * vi2 - 5. / 2. * vivj
                          + 5. / 4. * vj2)
                 + G2 * mj2 / rij3
                       * (2. * vinij2 - 4. * vinij * vjnij - 6. * vjnij2
                          - 8. * vivj + 4. * vj2))
                / c4;
            vij_factor +=
                (G2 * mj2 / rij3 * (-2. * vinij - 2. * vjnij)
                 + G2 * mi * mj / rij3 * (-63. / 4. * vinij + 55. / 4. * vjnij)
                 + G * mj / rij2
                       * (-6. * vinij * vjnij2 + 9. / 2. * vjnij2 * vjnij
                          + vjnij * vi2 - 4. * vinij * vivj + 4. * vjnij * vivj
                          + 4. * vinij * vj2 - 5. * vjnij * vj2))
                / c4;
        }

        if (flags & KETJU_PN_2_5_ACC) {
            nij_factor += (208. * G3 * mi * mj2 / (15. * rij4)
                           - 24. * G3 * mi2 * mj / (5. * rij4)
                           + 12. * G2 * mi * mj * vij2 / (5. * rij3))
                          * vijnij / c5;
            vij_factor += (8. * G3 * mi2 * mj / (5. * rij4)
                           - 32. * G3 * mi * mj2 / (5. * rij4)
                           - 4. * G2 * mi * mj * vij2 / (5. * rij3))
                          / c5;
        }

        if (flags & KETJU_PN_3_0_ACC) {
            // Here the logarithmic terms have been removed using the gauge
            // transformation eq. (204), which is explained in more detail in
            // section 4 of Arun et al. 2008 doi: 10.1103/PhysRevD.77.064035,
            // arXiv: 0711.0302. The resulting additional terms are indicated
            // with comments below.

            nij_factor +=
                (G * mj / rij2
                     * (35. / 16. * vjnij4 * vjnij2 - 15. / 8. * vjnij4 * vi2
                        + 15. / 2. * vjnij4 * vivj + 3. * vjnij2 * vivj * vivj
                        - 15. / 2. * vjnij4 * vj2 + 1.5 * vjnij2 * vi2 * vj2
                        - 12. * vjnij2 * vivj * vj2 - 2. * vivj * vivj * vj2
                        + 15. / 2. * vjnij2 * vj4 + 4. * vivj * vj4
                        - 2. * vj4 * vj2)
                 + G2 * mi * mj / rij3
                       * (-171. / 8. * vinij2 * vinij2
                          + 171. / 2. * vinij2 * vinij * vjnij
                          - 723. / 4. * vinij2 * vjnij2
                          + 383. / 2. * vinij * vjnij2 * vjnij
                          - 455. / 8. * vjnij4 + 229. / 4. * vinij2 * vi2
                          - 205. / 2. * vinij * vjnij * vi2
                          + 191. / 4. * vjnij2 * vi2 - 91. / 8. * vi2 * vi2
                          - 229. / 2. * vinij2 * vivj
                          + 244. * vinij * vjnij * vivj
                          - 225. / 2. * vjnij2 * vivj + 91. / 2. * vi2 * vivj
                          - 177. / 4. * vivj * vivj + 229. / 4. * vinij2 * vj2
                          - 283. / 2. * vinij * vjnij * vj2
                          + 259. / 4. * vjnij2 * vj2 - 91. / 4. * vi2 * vj2
                          + 43. * vivj * vj2 - 81. / 8. * vj4)
                 + G2 * mj2 / rij3
                       * (-6. * vinij2 * vjnij2 + 12. * vinij * vjnij2 * vjnij
                          + 6. * vjnij4 + 4. * vinij * vjnij * vivj
                          + 12. * vjnij2 * vivj + 4. * vivj * vivj
                          - 4. * vinij * vjnij * vj2 - 12. * vjnij2 * vj2
                          - 8. * vivj * vj2 + 4. * vj4)
                 + G3 * mj3 / rij4
                       * (-vinij2 + 2. * vinij * vjnij + 43. / 2. * vjnij2
                          + 18. * vivj - 9. * vj2)
                 + G3 * mi * mj2 / rij4
                       * (415. / 8. * vinij2 - 375. / 4. * vinij * vjnij
                          + 1113. / 8. * vjnij2
                          - 615. / 64. * vijnij * vijnij * pi2 + 18. * vi2
                          + 123. / 64. * pi2 * vij2 + 33. * vivj
                          - 33. / 2. * vj2)
                 + G3 * mi2 * mj / rij4
                       * (-45887. / 168. * vinij2 + 24025. / 42. * vinij * vjnij
                          - 10469. / 42. * vjnij2 + 48197. / 840. * vi2
                          - 36227. / 420. * vivj + 36227. / 840. * vj2)
                 + G4 / (rij4 * rij)
                       * (16. * mj2 * mj2 + mi2 * mj2 * (175. - 41. / 16. * pi2)
                          - 3187. / 1260. * mi2 * mi * mj
                          + mi * mj3 * (110741. / 630. - 41. / 16. * pi2))
                 // gauge transform terms:
                 + 22. / 3. * G3 * mi2 * mj
                       * (8 * vijnij * vijnij / rij4 - vij2 / rij4
                          + G * (mi + mj) / (rij4 * rij)))
                / c6;

            vij_factor +=
                (G * mj / rij2
                     * (15. / 2. * vinij * vjnij4 - 45. / 8. * vjnij4 * vjnij
                        - 1.5 * vjnij2 * vjnij * vi2
                        + 6. * vinij * vjnij2 * vivj
                        - 6. * vjnij2 * vjnij * vivj - 2. * vjnij * vivj * vivj
                        - 12. * vinij * vjnij2 * vj2
                        + 12. * vjnij2 * vjnij * vj2 + vjnij * vi2 * vj2
                        - 4. * vinij * vivj * vj2 + 8. * vjnij * vivj * vj2
                        + 4. * vinij * vj4 - 7. * vjnij * vj4)
                 + G2 * mj2 / rij3
                       * (-2. * vinij2 * vjnij + 8. * vinij * vjnij2
                          + 2. * vjnij2 * vjnij + 2. * vinij * vivj
                          + 4. * vjnij * vivj - 2. * vinij * vj2
                          - 4. * vjnij * vj2)
                 + G2 * mi * mj / rij3
                       * (-243. / 4. * vinij2 * vinij
                          + 565. / 4. * vinij2 * vjnij
                          - 269. / 4. * vinij * vjnij2
                          - 95. / 12. * vjnij2 * vjnij + 207. / 8. * vinij * vi2
                          - 137. / 8. * vjnij * vi2 - 36. * vinij * vivj
                          + 27. / 4. * vjnij * vivj + 81. / 8. * vinij * vj2
                          + 83. / 8. * vjnij * vj2)
                 + G3 * mj3 / rij4 * (4. * vinij + 5. * vjnij)
                 + G3 * mi * mj2 / rij4
                       * (-307. / 8. * vinij + 479. / 8. * vjnij
                          + 123. / 32. * vijnij * pi2)
                 + G3 * mi2 * mj / rij4
                       * (31397. / 420. * vinij - 36227. / 420. * vjnij)
                 // gauge transform terms:
                 - 44. / 3. * G3 * mi2 * mj * vijnij / rij4)
                / c6;
        }

        if (flags & KETJU_PN_3_5_ACC) {
            const double vijnij4 = vijnij * vijnij * vijnij * vijnij;
            // Note that the ~1/r^6 terms in the published formula are
            // actually supposed to be ~1/r^5 (dimensional analysis shows this).
            // This error goes back to Nissanke & Blanchet (2005)
            // doi:10.1088/0264-9381/22/6/008.
            // The relative acceleration expressions in both are not affected by
            // this, so seems like it's just a misprint.
            nij_factor +=
                (G4 * mi2 * mi * mj / (rij4 * rij)
                     * (3992. / 105. * vinij - 4328. / 105. * vjnij)
                 + G4 * mi2 * mj2 / (rij4 * rij)
                       * (-13576. / 105. * vinij + 2872. / 21. * vjnij)
                 - 3172. / 21. * G4 * mi * mj3 / (rij4 * rij) * vijnij
                 + G3 * mi2 * mj / rij4
                       * (48. * vinij2 * vinij - 696. / 5. * vinij2 * vjnij
                          + 744. / 5. * vinij * vjnij2
                          - 288. / 5. * vjnij2 * vjnij
                          - 4888. / 105. * vinij * vi2
                          + 5056. / 105. * vjnij * vi2
                          + 2056. / 21. * vinij * vivj
                          - 2224. / 21. * vjnij * vivj
                          - 1028. / 21. * vinij * vj2
                          + 5812. / 105. * vjnij * vj2)
                 + G3 * mi * mj2 / rij4
                       * (-582. / 5. * vinij2 * vinij
                          + 1746. / 5. * vinij2 * vjnij
                          - 1954. / 5. * vinij * vjnij2 + 158. * vjnij2 * vjnij
                          + 3568. / 105. * vijnij * vi2
                          - 2864. / 35. * vinij * vivj
                          + 10048. / 105. * vjnij * vivj
                          + 1432. / 35. * vinij * vj2
                          - 5752. / 105. * vjnij * vj2)
                 + G2 * mi * mj / rij3
                       * (-56. * vijnij4 * vijnij + 60. * vinij2 * vinij * vij2
                          - 180. * vinij2 * vjnij * vij2
                          + 174. * vinij * vjnij2 * vij2
                          - 54. * vjnij2 * vjnij * vij2
                          - 246. / 35. * vijnij * vi2 * vi2
                          + 1068. / 35. * vinij * vi2 * vivj
                          - 984. / 35. * vjnij * vi2 * vivj
                          - 1068. / 35. * vinij * vivj * vivj
                          + 180. / 7. * vjnij * vivj * vivj
                          - 534. / 35. * vinij * vi2 * vj2
                          + 90. / 7. * vjnij * vi2 * vj2
                          + 984. / 35. * vinij * vivj * vj2
                          - 732. / 35. * vjnij * vivj * vj2
                          - 204. / 35. * vinij * vj4 + 24. / 7. * vjnij * vj4))
                / c7;

            vij_factor +=
                (-184. / 21. * G4 * mi2 * mi * mj / (rij4 * rij)
                 + 6224. / 105. * G4 * mi2 * mj2 / (rij4 * rij)
                 + 6388. / 105. * G4 * mi * mj3 / (rij4 * rij)
                 + G3 * mi2 * mj / rij4
                       * (52. / 15. * vinij2 - 56. / 15. * vinij * vjnij
                          - 44. / 15. * vjnij2 - 132. / 35. * vi2
                          + 152. / 35. * vivj - 48. / 35. * vj2)
                 + G3 * mi * mj2 / rij4
                       * (454. / 15. * vinij2 - 372. / 5. * vinij * vjnij
                          + 854. / 15. * vjnij2 - 152. / 21. * vi2
                          + 2864. / 105. * vivj - 1768. / 105. * vj2)
                 + G2 * mi * mj / rij3
                       * (60. * vijnij4 - 348. / 5. * vinij2 * vij2
                          + 684. / 5. * vinij * vjnij * vij2
                          - 66. * vjnij2 * vij2 + 334. / 35. * vi2 * vi2
                          - 1336. / 35. * vi2 * vivj + 1308. / 35. * vivj * vivj
                          + 654. / 35. * vi2 * vj2 - 1252. / 35. * vivj * vj2
                          + 292. / 35. * vj4))
                / c7;
        }
        // Finally add to the acceleration
        for (int l = 0; l < 3; ++l) {
            acc[l] += nij_factor * nv_ij[l] + vij_factor * vv_ij[l];
        }
    }
}

void ketju_compute_PN_acc(struct ketju_system *R, double (*vel_to_use)[3],
                          double (*spin_to_use)[3]) {
    if (!R->options->PN_flags) return;
    timer_start(&R->perf->time_force);

    double(*acc_pn)[3] = R->internal_state->acc_vd;
    memset(acc_pn, 0, sizeof(double[R->num_particles][3]));
    double(*spin_derivative)[3] = R->internal_state->spin_derivative;
    memset(spin_derivative, 0, sizeof(double[R->num_pn_particles][3]));

    // Compute the PN accelerations for SMBHs
    for (int i = 0; i < R->num_pn_particles; ++i) {
        compute_PN_acc_single(R, i, vel_to_use, spin_to_use);
    }

    // Then the chained PN accelerations
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    int hi, lo;
    for (int i = 0; i < R->num_particles - 1; ++i) {
        int v1, v2;
        get_v1_v2_fast(R, i, &v1, &v2);
        if (mst->vertex[v1].level > mst->vertex[v2].level) {
            hi = v1;
            lo = v2;
        } else {
            hi = v2;
            lo = v1;
        }
        for (int k = 0; k < 3; k++) {
            mst->edge_rel_acc_vd[i][k] = acc_pn[hi][k] - acc_pn[lo][k];
        }
    }
    timer_stop(&R->perf->time_force);
}

// PN CoM calculations, needed for merger remnants.
// The PN terms shift the CoM of the system, so that near the merger the
// Newtonian CoM oscillates with a velocity of several hundred km/s.
// This corrects for that, but the formulae aren't perfect and have some
// remaining oscillations.

static double binary_PN_CoM_momentum_single_contribution(struct ketju_system *R,
                                                         int i, int j,
                                                         double CoM_vel[3]);

// Calculate the PN corrected CoM velocity of a binary at indices i, j <
// num_pn_particles (returns immediatedly otherwise).
// The CoM velocity is calculated from the momentum based on the total mass +
// Newtonian energy of the system (~1PN accuracy), which should be accurate
// enough for the main purpose of calculating the binary rest frame for
// mergers. The CoM calculated with these formulae isn't completely
// conserved anyway, likely due to missing higher order terms (e.g. the actual
// conserved CoM for equations of motion at 1PN level contains 2PN and maybe
// higher terms, but now we have only 1PN for that case).
void ketju_binary_PN_CoM_velocity(struct ketju_system *R, int i, int j,
                                  double CoM_vel[3]) {
    if (i > R->num_pn_particles || j > R->num_pn_particles || i < 0 || j < 0)
        return;

    // The contribution is the of the same form for both particles, so use a
    // helper function
    double CoM_p1[3], CoM_p2[3];
    double m1 = binary_PN_CoM_momentum_single_contribution(R, i, j, CoM_p1);
    double m2 = binary_PN_CoM_momentum_single_contribution(R, j, i, CoM_p2);

    for (int k = 0; k < 3; ++k) {
        CoM_vel[k] = (CoM_p1[k] + CoM_p2[k]) / (m1 + m2);
    }
}

// Helper used above.
// Calculates the momentum contribution from a single particle
// into CoM_p and returns the effective mass of the particle corrected with
// binding energy.
static double binary_PN_CoM_momentum_single_contribution(struct ketju_system *R,
                                                         int i, int j,
                                                         double CoM_p[3]) {

    struct ketju_system_physical_state *ps = R->physical_state;
    double(*vel)[3] = ps->vel;
    // Compute the enabled terms
    const unsigned int flags = R->options->PN_flags;

    // Helper values
    const double G = R->constants->G;
    const double G2 = G * G;
    const double G3 = G2 * G;
    const double c = R->constants->c;
    const double c2 = c * c;
    const double c4 = c2 * c2;
    const double c6 = c4 * c2;

    const double mi = ps->mass[i];
    const double mi2 = mi * mi;
    const double mi3 = mi2 * mi;
    const double vi2 = vector_dot(vel[i], vel[i]);
    const double vi4 = vi2 * vi2;

    double rv_ij[3], nv_ij[3], vv_ij[3];
    particle_separation_and_relative_velocity(R, i, j, rv_ij, vv_ij, vel);
    const double rij = norm_and_unit_vector(rv_ij, nv_ij);
    const double rij2 = rij * rij;
    const double rij3 = rij2 * rij;

    const double mj = ps->mass[j];
    const double mj2 = mj * mj;
    const double mj3 = mj2 * mj;

    const double vivj = vector_dot(vel[i], vel[j]);
    const double vj2 = vector_dot(vel[j], vel[j]);
    const double vj4 = vj2 * vj2;

    const double vjnij = vector_dot(vel[j], nv_ij);
    const double vjnij2 = vjnij * vjnij;
    const double vjnij3 = vjnij2 * vjnij;
    const double vjnij4 = vjnij2 * vjnij2;

    const double vinij = vector_dot(vel[i], nv_ij);
    const double vinij2 = vinij * vinij;
    const double vinij3 = vinij2 * vinij;
    const double vinij4 = vinij2 * vinij2;

    // Newtonian part
    double vi_fac = mi;
    double nij_fac = 0;

    // Spin independent part.
    // From de Andrade et al. 2001 Class. Quantum Grav. 18 753
    // doi:10.1088/0264-9381/18/5/301
    // eq. 4.3 with the gauge transform applied to remove the log terms at 3PN
    // level.

    if (flags & KETJU_PN_1_0_ACC) {
        nij_fac += -G * mi * mj / (2. * rij) * vinij / c2;
        vi_fac += (.5 * mi * vi2 - G * mi * mj / (2. * rij)) / c2;
    }

    if (flags & KETJU_PN_2_0_ACC) {
        nij_fac += (G2 * mi2 * mj / rij2 * (29. / 4. * vinij - 9. / 4. * vjnij)
                    + G * mi * mj / rij
                          * (3. / 8. * vinij3 + 3. / 8. * vinij2 * vjnij
                             - 9. / 8. * vinij * vi2 - 7. / 8. * vjnij * vi2
                             + 7. / 4. * vinij * vivj))
                   / c4;

        vi_fac += (-3. * G2 * mi2 * mj / rij2 + 7. * G2 * mi * mj2 / (2. * rij2)
                   + 3. / 8. * mi * vi4
                   + G * mi * mj / rij
                         * (13. / 8. * vinij2 - .25 * vinij * vjnij
                            - 13. / 8. * vjnij2 + 5. / 8. * vi2 - 7. / 4. * vivj
                            + 7. / 8. * vj2))
                  / c4;
    }

    // The contribution from the 2.5PN momentum correction (4.6b)
    if (flags & KETJU_PN_2_5_ACC) {
        const double vij2 = vector_dot(vv_ij, vv_ij);
        nij_fac += 4. * G2 * mi2 * mj / (5. * c4 * c * rij2)
                   * (vij2 - 2. * G * mi / rij);
    }

    if (flags & KETJU_PN_3_0_ACC) {
        nij_fac +=
            (G2 * mi2 * mj / rij2
                 * (-45. / 8. * vinij3 + 59. / 8. * vinij2 * vjnij
                    - 179. / 8. * vinij * vjnij2 + 247. / 24. * vjnij3
                    + 135. / 8. * vinij * vi2 - 87. / 8. * vjnij * vi2
                    - 53. / 2. * vinij * vivj + 87. / 4. * vjnij * vivj
                    + 53. / 4. * vinij * vj2 - 12. * vjnij * vj2)
             + G * mi * mj / rij
                   * (-5. / 16. * vinij3 * vinij2 - 5. / 16. * vinij4 * vjnij
                      - 5. / 16. * vinij3 * vjnij2 + 19. / 16. * vinij3 * vi2
                      + 15. / 16. * vinij2 * vjnij * vi2
                      + 3. / 4. * vinij * vjnij2 * vi2 + 5. / 8. * vjnij3 * vi2
                      - 21. / 16. * vinij * vi4 - 11. / 16. * vjnij * vi4
                      - 5. / 4. * vinij3 * vivj
                      - 9. / 8. * vinij2 * vjnij * vivj
                      + 15. / 8. * vinij * vi2 * vivj + vjnij * vi2 * vivj
                      - 1. / 8. * vinij * vivj * vivj
                      - 15. / 16. * vinij * vi2 * vj2)
             - 175. * G3 * mi2 * mj2 / (8. * rij3) * vinij
             + G3 * mi3 * mj / rij3
                   * (-46517. / 840. * vinij
                      + 34547. / 840. * vjnij
                      // This following term is what comes out of the gauge
                      // transform in addition to the terms that cancel out the
                      // log() terms.
                      + 22. / 3. * (vinij - vjnij)))
            / c6;

        vi_fac +=
            (5. * mi * vi2 * vi4 / 16.
             + G2 * mi2 * mj / rij2
                   * (-39. / 4. * vinij2 + 22. * vinij * vjnij
                      - 21. / 4. * vjnij2 - 2. * vi2 + vivj - .5 * vj2)
             + G2 * mi * mj2 / rij2
                   * (23. / 4. * vinij2 - 101. / 4. * vinij * vjnij
                      + 17. * vjnij2 + 9. / 4. * vi2 - vivj + .5 * vj2)
             + G * mi * mj / rij
                   * (-19. / 16. * vinij4 + 1. / 4. * vinij3 * vjnij
                      + 3. / 16. * vinij2 * vjnij2 + 1. / 8. * vinij * vjnij3
                      + 19. / 16. * vjnij4 + 45. / 16. * vinij2 * vi2
                      - 5. / 8. * vinij * vjnij * vi2 - 2. * vjnij2 * vi2
                      + 23. / 16. * vi4 - 3. / 4. * vinij2 * vivj
                      + 3. / 4. * vinij * vjnij * vivj
                      + 19. / 8. * vjnij2 * vivj - 31. / 8. * vi2 * vivj
                      + 17. / 8. * vivj * vivj + 3. / 8. * vinij2 * vj2
                      - .5 * vinij * vjnij * vj2 - 45. / 16. * vjnij2 * vj2
                      + 31. / 16. * vi2 * vj2 - 3. * vivj * vj2
                      + 19. / 16. * vj4)
             - 19. * G3 * mi2 * mj2 / (8. * rij3)
             + G3 * mi3 * mj / rij3 * 20407. / 2520.
             - G3 * mi * mj3 / rij3 * 21667. / 2520.)
            / c6;
    }

    // Spin dependent correction. This is just the leading correction discussed
    // in Keppel et al. 2009 Phys. Rev. D 80, 124015
    // doi:10.1103/PhysRevD.80.124015, i.e. the part in the canonical momentum
    // of eq. 2.22. Higher order corrections for the spin-spin and quadrupole
    // terms would be useful, but I haven't found them in the literature. The
    // current formula helps quite a lot with the CoM oscillations, but they can
    // still be nearly 100 km/s when the binary is at < 10 Rs separations and
    // the spins are high.
    // Damour et al. 2008 Phys. Rev. D 77, 064032 gives the center of mass
    // position with higher corrections, from which the velocity could be
    // derived, but they use ADM coordinates and a different spin supplementary
    // condition, so using their results isn't straightforward.

    double spin_term[3];
    if (flags & KETJU_PN_1_5_SPIN_ACC) {
        // The spin term depends on the acceleration of the particles.
        // The paper only uses the Newtonian part, and that seems to work ok.
        // Ignores the acceleration from the environment, since the
        // contribution from the other binary component should dominate.
        double acc[3];
        for (int k = 0; k < 3; ++k) {
            acc[k] = -G * mj / rij2 * nv_ij[k];
        }
        vector_cross(acc, ps->spin[i], spin_term);
    } else {
        for (int k = 0; k < 3; ++k) {
            spin_term[k] = 0;
        }
    }

    for (int k = 0; k < 3; ++k) {
        CoM_p[k] = vi_fac * vel[i][k] + nij_fac * nv_ij[k] + spin_term[k] / c2;
    }

    // Return the mass corrected with Newtonian binding energy
    double m_corr = mi + (.5 * mi * vi2 - G * mi * mj / (2 * rij)) / c2;
    return m_corr;
}
