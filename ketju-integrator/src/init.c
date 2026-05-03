// Ketju-integrator (MSTAR)
// Copyright (C) 2020-2022 Antti Rantala, Matias Mannerkoski, and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.


// Functions for initialization of the integrator

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

static void allocate_internal_state(struct ketju_system *R, int Ntask) {
    R->internal_state = calloc(1, sizeof(struct ketju_system_internal_state));
    R->internal_state->task_info = calloc(1, sizeof(struct task_info));
    R->internal_state->group_info = calloc(1, sizeof(struct group_info));
    // allocated when gbs comms inited
    R->internal_state->group_info->gbs_group_sizes = NULL;

    R->internal_state->mst =
        calloc(1, sizeof(struct spanning_tree_coordinate_system));
    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    mst->vertex = calloc(R->num_particles, sizeof(struct graph_vertex));

    // This is the maximum number of edges per task allocated by
    // loop_scheduling_N2
    int num_local_edge =
        (R->num_particles * (R->num_particles - 1)) / 2 / Ntask + 1;

    mst->local_edge = calloc(num_local_edge, sizeof(struct graph_edge));
    mst->edge_in_mst = calloc(R->num_particles - 1, sizeof(struct graph_edge));

    mst->edge_rel_pos = calloc(1, sizeof(double[R->num_particles - 1][3]));
    mst->edge_rel_vel = calloc(1, sizeof(double[R->num_particles - 1][3]));
    mst->edge_rel_aux_vel = calloc(1, sizeof(double[R->num_particles - 1][3]));
    mst->edge_rel_acc = calloc(1, sizeof(double[R->num_particles - 1][3]));
    mst->edge_rel_acc_vd = calloc(1, sizeof(double[R->num_particles - 1][3]));

    R->internal_state->octree_data = calloc(1, sizeof(struct octree_data));
    R->internal_state->octree_data->indexlist = calloc(8, sizeof(int *));
    for (int i = 0; i < 8; i++) {
        R->internal_state->octree_data->indexlist[i] =
            calloc(R->num_particles, sizeof(int));
    }
    R->internal_state->octree_data->particle_in_node =
        calloc(R->num_particles, sizeof(struct octree_node *));

    R->internal_state->aux_vel = calloc(1, sizeof(double[R->num_particles][3]));
    R->internal_state->acc = calloc(1, sizeof(double[R->num_particles][3]));
    R->internal_state->acc_vd = calloc(1, sizeof(double[R->num_particles][3]));
    R->internal_state->aux_spin =
        calloc(1, sizeof(double[R->num_pn_particles][3]));
    R->internal_state->spin_derivative =
        calloc(1, sizeof(double[R->num_pn_particles][3]));
}

static void free_internal_state(struct ketju_system *R) {
    free(R->internal_state->aux_vel);
    free(R->internal_state->acc);
    free(R->internal_state->acc_vd);
    free(R->internal_state->aux_spin);
    free(R->internal_state->spin_derivative);

    struct spanning_tree_coordinate_system *mst = R->internal_state->mst;
    free(mst->local_edge);
    free(mst->edge_in_mst);
    free(mst->edge_rel_pos);
    free(mst->edge_rel_vel);
    free(mst->edge_rel_aux_vel);
    free(mst->edge_rel_acc);
    free(mst->edge_rel_acc_vd);
    free(mst->vertex);
    free(mst);

    for (int i = 0; i < 8; i++) {
        free(R->internal_state->octree_data->indexlist[i]);
    }
    free(R->internal_state->octree_data->indexlist);
    free(R->internal_state->octree_data->particle_in_node);
    free(R->internal_state->octree_data);

    // Note: newly created MPI objects need to be freed before this
    free(R->internal_state->task_info);
    free(R->internal_state->group_info);
    free(R->internal_state);
    R->internal_state = NULL;
}

static void allocate_physical_state(struct ketju_system *R) {
    R->physical_state = calloc(1, sizeof(struct ketju_system_physical_state));
    R->physical_state->mass = calloc(1, sizeof(double[R->num_particles]));
    R->physical_state->pos = calloc(1, sizeof(double[R->num_particles][3]));
    R->physical_state->vel = calloc(1, sizeof(double[R->num_particles][3]));
    R->physical_state->spin = calloc(1, sizeof(double[R->num_pn_particles][3]));
}

static void free_physical_state(struct ketju_system *R) {
    free(R->physical_state->mass);
    free(R->physical_state->pos);
    free(R->physical_state->vel);
    free(R->physical_state->spin);
    free(R->physical_state);
    R->physical_state = NULL;
}

void ketju_create_system(struct ketju_system *R, int num_pn_particles,
                         int num_other_particles, MPI_Comm comm) {
    R->num_particles = num_other_particles + num_pn_particles;
    R->num_pn_particles = num_pn_particles;

    int Ntask;
    MPI_Comm_size(comm, &Ntask);

    allocate_internal_state(R, Ntask);
    R->internal_state->task_info->main_group_size = Ntask;
    MPI_Comm_rank(comm, &R->internal_state->task_info->main_group_rank);
    MPI_Comm_dup(comm, &R->internal_state->task_info->main_comm);
    MPI_Comm_rank(MPI_COMM_WORLD, &R->internal_state->task_info->global_rank);
    // The rest of the communication setup is done in ketju_init_gbs_comms,
    // since it depends on the configuration options (at least gbs_kmax);

    allocate_physical_state(R);

    R->perf = calloc(1, sizeof(struct ketju_performance_counters));
    R->options = calloc(1, sizeof(struct ketju_integration_options));
    *R->options = ketju_default_integration_options;
    R->constants = calloc(1, sizeof(struct ketju_natural_constants));
    *R->constants = ketju_default_natural_constants;

    R->particle_extra_data = NULL;
    R->particle_extra_data_elem_size = 0;

    R->num_mergers = 0;
    R->merger_infos = NULL;
}

void ketju_init_gbs_comms(struct ketju_system *R) {
    struct task_info *ti = R->internal_state->task_info;
    struct group_info *gi = R->internal_state->group_info;

    if (ti->gbs_comms_initialized) {
        terminate(
            "ketju_init_gbs_comms: Error! Trying to initialize gbs comms a "
            "second time\n");
    }

    // For now just split into even chunks, giving the extra tasks to the
    // first groups since they get slightly more work in some cases.
    // For the 2*k step sequence we're using, we can get perfectly even
    // work distribution for k/2 groups (or k/2^n, n>=1, but the largest
    // number of groups is likely most efficient), so should aim for that.
    const int max_num_gbs_groups = R->internal_state->gbs_k / 2;
    const int num_gbs_groups = ti->main_group_size > max_num_gbs_groups
                                   ? max_num_gbs_groups
                                   : ti->main_group_size;
    int tasks_per_group = ti->main_group_size / num_gbs_groups;
    int remaining_tasks = ti->main_group_size % num_gbs_groups;
    gi->num_gbs_groups = num_gbs_groups;
    R->internal_state->group_info->gbs_group_sizes =
        calloc(1, sizeof(int[num_gbs_groups]));

    for (int i = 0; i < num_gbs_groups; ++i) {
        gi->gbs_group_sizes[i] =
            tasks_per_group + (i < remaining_tasks ? 1 : 0);
    }
    int group_start_rank = 0, group_end_rank = gi->gbs_group_sizes[0] - 1;
    ti->gbs_group_index = 0; // The starting state is not calloc'd necessarily
    for (int i = 1; ti->main_group_rank > group_end_rank; ++i) {
        group_start_rank += gi->gbs_group_sizes[i - 1];
        group_end_rank = group_start_rank + gi->gbs_group_sizes[i] - 1;
        ti->gbs_group_index = i;
    }

    int range[1][3] = {{group_start_rank, group_end_rank, 1}};
    MPI_Group main_group, gbs_group;
    MPI_Comm_group(ti->main_comm, &main_group);
    MPI_Group_range_incl(main_group, 1, range, &gbs_group);
    MPI_Group_rank(gbs_group, &ti->gbs_group_rank);
    MPI_Group_size(gbs_group, &ti->gbs_group_size);
    MPI_Comm_create_group(ti->main_comm, gbs_group, ti->gbs_group_index,
                          &ti->gbs_comm);
    MPI_Group_free(&gbs_group);

    // Task ranges for each group for the next step
    int group_ranges[num_gbs_groups][2];
    for (int i = 0; i < num_gbs_groups; ++i) {
        group_ranges[i][0] = i == 0 ? 0 : group_ranges[i - 1][1] + 1;
        group_ranges[i][1] = group_ranges[i][0] + gi->gbs_group_sizes[i] - 1;
    }
    // Create the communicators for sharing the computed step results
    ti->gbs_data_comms = malloc(sizeof(MPI_Comm[num_gbs_groups]));
    for (int i = 0; i < num_gbs_groups; ++i) {
        // Each group contains the first task of the gbs group as the root
        // (sender of data) and the tasks in the other groups as receivers
        int ranges[num_gbs_groups][3];
        ranges[0][0] = group_ranges[i][0];
        ranges[0][1] = group_ranges[i][0];
        ranges[0][2] = 1;
        int k = 1;
        for (int j = 0; j < num_gbs_groups; ++j) {
            if (i == j) continue;
            ranges[k][0] = group_ranges[j][0];
            ranges[k][1] = group_ranges[j][1];
            ranges[k][2] = 1;
            ++k;
        }
        MPI_Group data_group;
        MPI_Group_range_incl(main_group, num_gbs_groups, ranges, &data_group);
        MPI_Comm_create_group(ti->main_comm, data_group, i,
                              &ti->gbs_data_comms[i]);

        MPI_Group_free(&data_group);
    }

    MPI_Group_free(&main_group);

    ti->gbs_comms_initialized = true;
}

void ketju_free_gbs_comms(struct ketju_system *R) {
    if (!R->internal_state->task_info->gbs_comms_initialized) {
        // not allocated
        return;
    }

    MPI_Comm_free(&R->internal_state->task_info->gbs_comm);
    free(R->internal_state->group_info->gbs_group_sizes);
    for (int i = 0; i < R->internal_state->group_info->num_gbs_groups; ++i) {
        if (R->internal_state->task_info->gbs_data_comms[i] != MPI_COMM_NULL)
            MPI_Comm_free(&R->internal_state->task_info->gbs_data_comms[i]);
    }
    free(R->internal_state->task_info->gbs_data_comms);

    R->internal_state->task_info->gbs_comms_initialized = false;
}

void ketju_free_system(struct ketju_system *R) {
    ketju_free_gbs_comms(R);
    MPI_Comm_free(&R->internal_state->task_info->main_comm);
    free_internal_state(R);

    free_physical_state(R);

    free(R->perf);
    free(R->options);
    free(R->constants);

    if (R->merger_infos) { free(R->merger_infos); }
}
