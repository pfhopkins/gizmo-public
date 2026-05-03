// Ketju-integrator (MSTAR)
// Copyright (C) 2020-2022 Antti Rantala, Matias Mannerkoski, and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.


#include <string.h>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

void ketju_compute_CoM(const struct ketju_system *R, double CoM_pos[3],
                       double CoM_vel[3]) {
    memset(CoM_pos, 0, sizeof(double[3]));
    memset(CoM_vel, 0, sizeof(double[3]));

    struct ketju_system_physical_state *ps = R->physical_state;
    double M = 0;

    for (int i = 0; i < R->num_particles; i++) {
        M += ps->mass[i];
        for (int k = 0; k < 3; k++) {
            CoM_pos[k] += ps->mass[i] * ps->pos[i][k];
            CoM_vel[k] += ps->mass[i] * ps->vel[i][k];
        }
    }
    for (int k = 0; k < 3; k++) {
        CoM_pos[k] /= M;
        CoM_vel[k] /= M;
    }
}

// Puts the initial data to the CoM frame and writes the original CoM pos and
// vel into the output arrays.
// The tree coordinate system is not updated to match.
void ketju_into_CoM_frame(struct ketju_system *R, double CoM_pos[3],
                          double CoM_vel[3]) {

    struct ketju_system_physical_state *ps = R->physical_state;
    ketju_compute_CoM(R, CoM_pos, CoM_vel);

    for (int i = 0; i < R->num_particles; i++) {
        for (int k = 0; k < 3; k++) {
            ps->pos[i][k] -= CoM_pos[k];
            ps->vel[i][k] -= CoM_vel[k];
        }
    }
}
