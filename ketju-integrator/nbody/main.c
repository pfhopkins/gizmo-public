// Ketju-integrator sample N-body driver program.
// Copyright (C) 2020-2022 Antti Rantala, Matias Mannerkoski, and contributors.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.



#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "ketju_integrator/ketju_integrator.h"

int ThisTask;

void init_from_ic_file(struct ketju_system *R, char *ic_file_name) {
    FILE *fp;
    int num_pn_part, num_other_part;
    fp = fopen(ic_file_name, "r");
    if (!fp) {
        printf("Error opening file %s on task %d\n", ic_file_name, ThisTask);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    int status = fscanf(fp, "%d %d\n", &num_pn_part, &num_other_part);
    if (!status) {
        printf("Error reading file %s on task %d\n", ic_file_name, ThisTask);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    ketju_create_system(R, num_pn_part, num_other_part, MPI_COMM_WORLD);
    R->physical_state->time = 0;

    for (int i = 0; i < num_pn_part; ++i) {
        status = fscanf(
            fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            &R->physical_state->mass[i], &R->physical_state->pos[i][0],
            &R->physical_state->pos[i][1], &R->physical_state->pos[i][2],
            &R->physical_state->vel[i][0], &R->physical_state->vel[i][1],
            &R->physical_state->vel[i][2], &R->physical_state->spin[i][0],
            &R->physical_state->spin[i][1], &R->physical_state->spin[i][2]);
        if (!status) {
            printf("Error reading file %s on task %d\n", ic_file_name,
                   ThisTask);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    for (int i = num_pn_part; i < R->num_particles; ++i) {
        status = fscanf(
            fp, "%lf %lf %lf %lf %lf %lf %lf\n", &R->physical_state->mass[i],
            &R->physical_state->pos[i][0], &R->physical_state->pos[i][1],
            &R->physical_state->pos[i][2], &R->physical_state->vel[i][0],
            &R->physical_state->vel[i][1], &R->physical_state->vel[i][2]);
        if (!status) {
            printf("Error reading file %s on task %d\n", ic_file_name,
                   ThisTask);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
}

// A super simple output file format
void write_output_file(FILE *fp, struct ketju_system *R) {
    if (ThisTask != 0) return;

    fprintf(fp, "%g %d %d\n", R->physical_state->time, R->num_particles, R->num_pn_particles);
    for (int i = 0; i<R->num_particles; ++i) {
        fprintf(fp, "%g ", R->physical_state->mass[i]);
    }
    fprintf(fp, "\n");
    for (int i = 0; i<R->num_particles; ++i) {
        for (int k=0; k<3; ++k) {
           fprintf(fp, "%g ", R->physical_state->pos[i][k]);
        }
    }
    fprintf(fp, "\n");
    for (int i = 0; i<R->num_particles; ++i) {
        for (int k=0; k<3; ++k) {
           fprintf(fp, "%g ", R->physical_state->vel[i][k]);
        }
    }
    fprintf(fp, "\n");
    for (int i = 0; i<R->num_pn_particles; ++i) {
        for (int k=0; k<3; ++k) {
           fprintf(fp, "%g ", R->physical_state->spin[i][k]);
        }
    }
    fprintf(fp, "\n");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    if (argc != 5) {
        if (ThisTask == 0) {
            printf("Usage:\nnbody ic_file integration_time substeps tolerance\n");
        }
        MPI_Finalize();
        return -2;
    }

    FILE *outfile;
    if (ThisTask == 0) {
        const char* out_file_name = "nbody_output.txt";
        outfile = fopen(out_file_name, "w");
        if (!outfile) {
            printf("Error opening file %s on task %d\n", out_file_name, ThisTask);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    struct ketju_system R;
    init_from_ic_file(&R, argv[1]);
    // Do the integration in the CoM frame
    double CoM_pos[3], CoM_vel[3];
    ketju_into_CoM_frame(&R, CoM_pos, CoM_vel);

    R.options->gbs_relative_tolerance = strtod(argv[4], NULL);
    R.options->mst_algorithm = KETJU_MST_PRIM;
    R.options->max_tree_distance = 2;
    R.options->PN_flags = KETJU_PN_ALL;

    double time_interval = strtod(argv[2], NULL);
    long substeps = strtol(argv[3], NULL, 0);
    for (long i=0; i<substeps; ++i) {
        ketju_run_integrator(&R, time_interval/substeps);
        write_output_file(outfile, &R);
    }

    if (ThisTask == 0) {
        printf("Timers:\n");
        printf(
            "time_force:       %f (%2.1f%%)\n"
            "time_force_comm:  %f (%2.1f%%)\n"
            "time_coord:       %f (%2.1f%%)\n"
            "time_mst:         %f (%2.1f%%)\n"
            "time_gbscomm:     %f (%2.1f%%)\n"
            "time_gbsextr:     %f (%2.1f%%)\n"
            "time_total:       %f\n"
            "successful_steps: %d\n"
            "failed_steps: %d\n"
            "total_work: %g\n"
            "relative_energy_error: %e\n"
            ,
            R.perf->time_force,
            R.perf->time_force/R.perf->time_total*100,
            R.perf->time_force_comm,
            R.perf->time_force_comm/R.perf->time_total*100,
            R.perf->time_coord,
            R.perf->time_coord/R.perf->time_total*100,
            R.perf->time_mst,
            R.perf->time_mst/R.perf->time_total*100,
            R.perf->time_gbscomm,
            R.perf->time_gbscomm/R.perf->time_total*100,
            R.perf->time_gbsextr,
            R.perf->time_gbsextr/R.perf->time_total*100,
            R.perf->time_total, R.perf->successful_steps, R.perf->failed_steps,
            R.perf->total_work,
            R.perf->relative_energy_error);

        printf("Physical time %e\n", R.physical_state->time);
        int n = R.num_particles > 10 ? 10 : R.num_particles;
        for (int i = 0; i < n; ++i) {
            printf(
                "Particle %d\n"
                " pos = [%e, %e, %e]\n"
                " vel = [%e, %e, %e]\n",
                i, R.physical_state->pos[i][0] + CoM_pos[0],
                R.physical_state->pos[i][1] + CoM_pos[1],
                R.physical_state->pos[i][2] + CoM_pos[2],
                R.physical_state->vel[i][0] + CoM_vel[0],
                R.physical_state->vel[i][1] + CoM_vel[1],
                R.physical_state->vel[i][2] + CoM_vel[2]);
        }
        fclose(outfile);
    }

    ketju_free_system(&R);

    MPI_Finalize();
    return 0;
}
