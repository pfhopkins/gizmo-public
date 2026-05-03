// Ketju-integrator Python bindings.
// Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.

#include <mpi.h>
#include <stdio.h>

#include "ketju_integrator/ketju_integrator.h"

int init_mpi() {
    MPI_Init(NULL, NULL);
    int task;
    MPI_Comm_rank(MPI_COMM_WORLD, &task);
    return task;
}

void init_system(struct ketju_system *R, int num_pn_particles,
                 int num_other_particles) {
    int task;
    MPI_Comm_rank(MPI_COMM_WORLD, &task);
    int counts[2];
    if (task == 0) {
        counts[0] = num_pn_particles;
        counts[1] = num_other_particles;
    }
    MPI_Bcast(counts, 2, MPI_INT, 0, MPI_COMM_WORLD);
    ketju_create_system(R, counts[0], counts[1], MPI_COMM_WORLD);
}

int comm_initial_data(struct ketju_system *R, int num_ts) {
    int task;
    MPI_Comm_rank(MPI_COMM_WORLD, &task);
    MPI_Bcast(&num_ts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&R->physical_state->time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(R->physical_state->mass, R->num_particles, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
    MPI_Bcast(R->physical_state->pos, 3 * R->num_particles, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
    MPI_Bcast(R->physical_state->vel, 3 * R->num_particles, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
    MPI_Bcast(R->physical_state->spin, 3 * R->num_pn_particles, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

    MPI_Bcast(R->options, sizeof *R->options, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(R->constants, sizeof *R->constants, MPI_BYTE, 0, MPI_COMM_WORLD);

    return num_ts;
}

void run_integrator(struct ketju_system *R, double tmax) {
    MPI_Bcast(&tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    ketju_run_integrator(R, tmax);
}

void load_state_dump(struct ketju_system *R, const char *fname) {
    ketju_debug_read_integrator_state(R, fname, MPI_COMM_WORLD);
}

void free_system(struct ketju_system *R) {
    ketju_free_system(R);
}

void finalize_mpi() {
    MPI_Finalize();
}



