// Ketju-integrator (MSTAR)
// Copyright (C) 2020-2022 Antti Rantala, Matias Mannerkoski, and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.


// Some utilities to help with debugging

#include <stdio.h>
#include <string.h>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

// create a binary dump of the integrator state
void ketju_debug_dump_integrator_state(const struct ketju_system *R,
                                       const char *fname) {
    FILE *f = fopen(fname, "wb");
    if (!f) {
        char msg[100 + strlen(fname)];
        sprintf(msg,
                "ketju_debug_dump_integrator_state: Error opening file %s\n",
                fname);
        terminate(msg);
    }

    // particle counts
    fwrite(&R->num_particles, sizeof(R->num_particles), 1, f);
    fwrite(&R->num_pn_particles, sizeof(R->num_pn_particles), 1, f);

    // physical state
    fwrite(&R->physical_state->time, sizeof(R->physical_state->time), 1, f);
    fwrite(R->physical_state->mass, sizeof(*R->physical_state->mass),
           R->num_particles, f);
    fwrite(R->physical_state->pos, sizeof(*R->physical_state->pos),
           R->num_particles, f);
    fwrite(R->physical_state->vel, sizeof(*R->physical_state->vel),
           R->num_particles, f);
    fwrite(R->physical_state->spin, sizeof(*R->physical_state->spin),
           R->num_pn_particles, f);

    // internal state
    struct ketju_system_internal_state *is = R->internal_state;
    // Can't load group/task info since that's process specific, so don't bother
    // writing them

    fwrite(is->mst->vertex, sizeof(*is->mst->vertex), R->num_particles, f);

    // local_edge_* are used only during tree build, don't write them

    fwrite(is->mst->edge_in_mst, sizeof(*is->mst->edge_in_mst),
           R->num_particles - 1, f);
    fwrite(is->mst->edge_rel_pos, sizeof(*is->mst->edge_rel_pos),
           R->num_particles - 1, f);
    fwrite(is->mst->edge_rel_vel, sizeof(*is->mst->edge_rel_vel),
           R->num_particles - 1, f);
    fwrite(is->mst->edge_rel_aux_vel, sizeof(*is->mst->edge_rel_aux_vel),
           R->num_particles - 1, f);
    fwrite(is->mst->edge_rel_acc, sizeof(*is->mst->edge_rel_acc),
           R->num_particles - 1, f);
    fwrite(is->mst->edge_rel_acc_vd, sizeof *is->mst->edge_rel_acc_vd,
           R->num_particles - 1, f);

    fwrite(&is->U, sizeof is->U, 1, f);
    fwrite(&is->T, sizeof is->T, 1, f);
    fwrite(&is->B, sizeof is->B, 1, f);

    // octree stuff is temporary, don't write it out for now

    // other stuff
    fwrite(R->perf, sizeof(*R->perf), 1, f);
    fwrite(R->options, sizeof(*R->options), 1, f);
    fwrite(R->constants, sizeof(*R->constants), 1, f);

    // don't care about the extra data

    fclose(f);
}

bool ketju_debug_read_integrator_state(struct ketju_system *R,
                                       const char *fname, MPI_Comm comm) {

    FILE *f = fopen(fname, "rb");
    if (!f) {
        char msg[100 + strlen(fname)];
        sprintf(msg,
                "ketju_debug_read_integrator_state: Error opening file %s\n",
                fname);
        terminate(msg);
    }

    int count;
    bool success = true;
    count = fread(&R->num_particles, sizeof(R->num_particles), 1, f);
    success = success && count == 1;
    count = fread(&R->num_pn_particles, sizeof(R->num_pn_particles), 1, f);
    success = success && count == 1;

    if (comm == MPI_COMM_NULL) { comm = MPI_COMM_SELF; }

    ketju_create_system(R, R->num_pn_particles,
                        R->num_particles - R->num_pn_particles, comm);

    // physical state
    count =
        fread(&R->physical_state->time, sizeof(R->physical_state->time), 1, f);
    success = success && count == 1;
    count = fread(R->physical_state->mass, sizeof(*R->physical_state->mass),
                  R->num_particles, f);
    success = success && count == R->num_particles;
    count = fread(R->physical_state->pos, sizeof(*R->physical_state->pos),
                  R->num_particles, f);
    success = success && count == R->num_particles;
    count = fread(R->physical_state->vel, sizeof(*R->physical_state->vel),
                  R->num_particles, f);
    success = success && count == R->num_particles;
    count = fread(R->physical_state->spin, sizeof(*R->physical_state->spin),
                  R->num_pn_particles, f);
    success = success && count == R->num_pn_particles;

    // internal state
    struct ketju_system_internal_state *is = R->internal_state;
    // Can't load group/task info since that's process specific, so don't bother
    // writing them

    count =
        fread(is->mst->vertex, sizeof(*is->mst->vertex), R->num_particles, f);
    success = success && count == R->num_particles;

    // local_edge_* are used only during tree build, don't read them

    count = fread(is->mst->edge_in_mst, sizeof(*is->mst->edge_in_mst),
                  R->num_particles - 1, f);
    success = success && count == R->num_particles - 1;
    count = fread(is->mst->edge_rel_pos, sizeof(*is->mst->edge_rel_pos),
                  R->num_particles - 1, f);
    success = success && count == R->num_particles - 1;
    count = fread(is->mst->edge_rel_vel, sizeof(*is->mst->edge_rel_vel),
                  R->num_particles - 1, f);
    success = success && count == R->num_particles - 1;
    count = fread(is->mst->edge_rel_aux_vel, sizeof(*is->mst->edge_rel_aux_vel),
                  R->num_particles - 1, f);
    success = success && count == R->num_particles - 1;
    count = fread(is->mst->edge_rel_acc, sizeof(*is->mst->edge_rel_acc),
                  R->num_particles - 1, f);
    success = success && count == R->num_particles - 1;
    count = fread(is->mst->edge_rel_acc_vd, sizeof(*is->mst->edge_rel_acc_vd),
                  R->num_particles - 1, f);
    success = success && count == R->num_particles - 1;

    count = fread(&is->U, sizeof is->U, 1, f);
    success = success && count == 1;
    count = fread(&is->T, sizeof is->T, 1, f);
    count = fread(&is->B, sizeof is->B, 1, f);

    // other stuff
    count = fread(R->perf, sizeof(*R->perf), 1, f);
    success = success && count == 1;
    count = fread(R->options, sizeof(*R->options), 1, f);
    success = success && count == 1;
    count = fread(R->constants, sizeof(*R->constants), 1, f);
    success = success && count == 1;

    fclose(f);

    return success;
}
