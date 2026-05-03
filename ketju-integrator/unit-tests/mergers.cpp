// Ketju-integrator tests.
// Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.


#include <cmath>
#include <gtest/gtest.h>
#include <mpi.h>
#include <random>

#include "ketju_integrator/ketju_integrator.h"

// Some simple tests where the BHs merge after integration
TEST(Integration_with_mergers, head_on_merger) {

    struct ketju_system R;
    ketju_create_system(&R, 2, 0, MPI_COMM_WORLD);

    R.constants->c = 1;
    R.constants->G = 1;
    R.options->PN_flags = KETJU_PN_ALL;

    R.physical_state->mass[0] = 1;
    R.physical_state->pos[0][0] = 40;

    R.physical_state->mass[1] = 1;
    R.physical_state->pos[1][0] = -40;

    R.options->gbs_relative_tolerance = 1e-12;
    R.options->output_time_relative_tolerance = 1e-12;

    const double tend = 1e6;

    ketju_run_integrator(&R, tend);
    const double tol = 1e-8;

    ASSERT_EQ(R.num_pn_particles, 1);
    ASSERT_NEAR(R.physical_state->time, tend,
                tend * R.options->output_time_relative_tolerance);

    ASSERT_LE(R.physical_state->mass[0], 2.);
    ASSERT_GE(R.physical_state->mass[0], 1.99);

    struct expected_vals {
        double pos[3];
        double vel[3];
        double spin[3];
    } expected[] = {{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}};
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(R.physical_state->pos[0][i], expected[0].pos[i], tol);
        EXPECT_NEAR(R.physical_state->vel[0][i], expected[0].vel[i], tol);
        EXPECT_NEAR(R.physical_state->spin[0][i], expected[0].spin[i], tol);
    }
    ketju_free_system(&R);
}

TEST(Integration_with_mergers, merger_from_orbit) {

    struct ketju_system R;
    ketju_create_system(&R, 2, 0, MPI_COMM_WORLD);

    R.constants->c = 1;
    R.constants->G = 1;
    R.options->PN_flags = KETJU_PN_ALL;

    R.physical_state->mass[0] = 1;
    R.physical_state->pos[0][0] = 40;
    const double v0[3] = {0., 0.05, 0.};
    for (int i = 0; i < 3; ++i) {
        R.physical_state->vel[0][i] = v0[i];
    }

    R.physical_state->mass[1] = 1;
    R.physical_state->pos[1][0] = -40;

    R.options->gbs_relative_tolerance = 1e-12;
    R.options->output_time_relative_tolerance = 1e-12;

    const double tend = 1e6;

    ketju_run_integrator(&R, tend);
    const double reltol =
        40 * R.options->gbs_relative_tolerance * R.perf->successful_steps;

    ASSERT_EQ(R.num_pn_particles, 1);
    ASSERT_NEAR(R.physical_state->time, tend,
                tend * R.options->output_time_relative_tolerance);

    // No point in checking exact value which is determined by the GR fits,
    // just check that the result is sensible.
    ASSERT_LT(R.physical_state->mass[0], 2.01);
    ASSERT_GT(R.physical_state->mass[0], 1.8);

    // Don't know what to check for in the coordinates

    ketju_free_system(&R);
}

TEST(Instant_mergers, binary_merger_with_extra_data) {

    struct ketju_system R;
    ketju_create_system(&R, 2, 1, MPI_COMM_WORLD);
    int extra_data[3] = {1, 2, 3};
    R.particle_extra_data = extra_data;
    R.particle_extra_data_elem_size = sizeof(*extra_data);

    R.constants->c = 1;
    R.constants->G = 1;
    R.options->PN_flags = KETJU_PN_ALL;

    R.physical_state->mass[0] = 1;
    R.physical_state->pos[0][0] = 4;
    const double v0[3] = {0., 0.05, 0.};
    for (int i = 0; i < 3; ++i) {
        R.physical_state->vel[0][i] = v0[i];
    }

    R.physical_state->mass[1] = 1;
    R.physical_state->pos[1][0] = -4;

    R.physical_state->mass[2] = 1e-3;
    R.physical_state->pos[2][0] = 10000;

    R.options->gbs_relative_tolerance = 1e-12;
    R.options->output_time_relative_tolerance = 1e-12;
    R.options->bh_merger_distance_in_Rs = 6;

    const double tend = 1;

    ketju_run_integrator(&R, tend);

    ASSERT_EQ(R.num_pn_particles, 1);
    ASSERT_EQ(R.num_particles, 2);
    ASSERT_NEAR(R.physical_state->time, tend,
                tend * R.options->output_time_relative_tolerance);

    ASSERT_EQ(R.num_mergers, 1);
    ASSERT_EQ(extra_data[0], 1);
    ASSERT_EQ(extra_data[1], 3);
    ASSERT_EQ(extra_data[2], 2);

    ASSERT_EQ(R.merger_infos[0].extra_index1, 0);
    ASSERT_EQ(R.merger_infos[0].extra_index2, 2);
    ASSERT_EQ(R.merger_infos[0].m1, 1);
    ASSERT_EQ(R.merger_infos[0].m2, 1);
    ASSERT_EQ(R.merger_infos[0].m_remnant, R.physical_state->mass[0]);

    ketju_free_system(&R);
}

// A bunch of different merger orders to test the indexing of the extra_data
// after mergers.
TEST(Instant_mergers, multiple_merger_with_extra_data1) {

    struct ketju_system R;
    ketju_create_system(&R, 3, 1, MPI_COMM_WORLD);
    int extra_data[4] = {1, 2, 3, 4}; // IDs
    R.particle_extra_data = extra_data;
    R.particle_extra_data_elem_size = sizeof(*extra_data);

    R.constants->c = 1;
    R.constants->G = 1;
    R.options->PN_flags = KETJU_PN_ALL;

    R.physical_state->mass[0] = 1;
    R.physical_state->pos[0][0] = 4;

    // Set position so that on first step IDs 1 and 3 merge,
    // then the remaining pair on the next step
    R.physical_state->mass[1] = 1;
    R.physical_state->pos[1][1] = -24;

    R.physical_state->mass[2] = 1;
    R.physical_state->pos[2][0] = -4;

    R.physical_state->mass[3] = 1e-3;
    R.physical_state->pos[3][0] = 10000;

    R.options->gbs_relative_tolerance = 1e-12;
    R.options->output_time_relative_tolerance = 1e-12;
    R.options->bh_merger_distance_in_Rs = 6;

    const double tend = 1;

    ketju_run_integrator(&R, tend);

    ASSERT_EQ(R.num_pn_particles, 1);
    ASSERT_EQ(R.num_particles, 2);
    ASSERT_NEAR(R.physical_state->time, tend,
                tend * R.options->output_time_relative_tolerance);

    ASSERT_EQ(R.num_mergers, 2);
    ASSERT_EQ(extra_data[0], 1);
    ASSERT_EQ(extra_data[1], 4);
    ASSERT_EQ(extra_data[2], 2);
    ASSERT_EQ(extra_data[3], 3);

    ASSERT_EQ(R.merger_infos[0].extra_index1, 0);
    ASSERT_EQ(R.merger_infos[0].extra_index2, 3);
    ASSERT_EQ(R.merger_infos[0].m1, 1);
    ASSERT_EQ(R.merger_infos[0].m2, 1);
    ASSERT_EQ(R.merger_infos[1].extra_index1, 0);
    ASSERT_EQ(R.merger_infos[1].extra_index2, 2);
    ASSERT_EQ(R.merger_infos[1].m_remnant, R.physical_state->mass[0]);

    ketju_free_system(&R);
}

TEST(Instant_mergers, multiple_merger_with_extra_data2) {

    struct ketju_system R;
    ketju_create_system(&R, 3, 1, MPI_COMM_WORLD);
    int extra_data[4] = {1, 2, 3, 4}; // IDs
    R.particle_extra_data = extra_data;
    R.particle_extra_data_elem_size = sizeof(*extra_data);

    R.constants->c = 1;
    R.constants->G = 1;
    R.options->PN_flags = KETJU_PN_ALL;

    R.physical_state->mass[0] = 1;
    R.physical_state->pos[0][0] = 4;

    // Set position so that on first step IDs 2 and 3 merge,
    // then the remaining pair on the next step
    R.physical_state->mass[1] = 0.9;
    R.physical_state->pos[1][0] = -22;

    R.physical_state->mass[2] = 1;
    R.physical_state->pos[2][0] = -23;

    R.physical_state->mass[3] = 1e-3;
    R.physical_state->pos[3][0] = 10000;

    R.options->gbs_relative_tolerance = 1e-12;
    R.options->output_time_relative_tolerance = 1e-12;
    R.options->bh_merger_distance_in_Rs = 6;

    const double tend = 1;

    ketju_run_integrator(&R, tend);

    ASSERT_EQ(R.num_pn_particles, 1);
    ASSERT_EQ(R.num_particles, 2);
    ASSERT_NEAR(R.physical_state->time, tend,
                tend * R.options->output_time_relative_tolerance);

    ASSERT_EQ(R.num_mergers, 2);
    ASSERT_EQ(extra_data[0], 1);
    ASSERT_EQ(extra_data[1], 4);
    ASSERT_EQ(extra_data[2], 2);
    ASSERT_EQ(extra_data[3], 3);

    ASSERT_EQ(R.merger_infos[0].extra_index1, 2);
    ASSERT_EQ(R.merger_infos[0].extra_index2, 3);
    ASSERT_EQ(R.merger_infos[0].m1, 0.9);
    ASSERT_EQ(R.merger_infos[0].m2, 1);
    ASSERT_EQ(R.merger_infos[1].extra_index1, 0);
    ASSERT_EQ(R.merger_infos[1].extra_index2, 2);
    ASSERT_EQ(R.merger_infos[1].m_remnant, R.physical_state->mass[0]);

    ketju_free_system(&R);
}

TEST(Instant_mergers, multiple_merger_with_extra_data3) {

    struct ketju_system R;
    ketju_create_system(&R, 4, 1, MPI_COMM_WORLD);
    int extra_data[5] = {1, 2, 3, 4, 5}; // IDs
    R.particle_extra_data = extra_data;
    R.particle_extra_data_elem_size = sizeof(*extra_data);

    R.constants->c = 1;
    R.constants->G = 1;
    R.options->PN_flags = KETJU_PN_ALL;

    R.physical_state->mass[0] = 1;
    R.physical_state->pos[0][0] = 10;

    // Set position so that on first step IDs 3 and 4 merge,
    // then 1 and 2 on the next step
    R.physical_state->mass[1] = 0.9;
    R.physical_state->pos[1][0] = 10 + (12 + 1e-12) * 1.9;

    R.physical_state->mass[2] = 1;
    R.physical_state->pos[2][0] = -100;
    R.physical_state->mass[3] = 2;
    R.physical_state->pos[3][0] = -105;

    R.physical_state->mass[4] = 1e-3;
    R.physical_state->pos[4][0] = 10000;

    R.options->gbs_relative_tolerance = 1e-12;
    R.options->output_time_relative_tolerance = 1e-12;
    R.options->bh_merger_distance_in_Rs = 6;

    const double tend = 1;

    ketju_run_integrator(&R, tend);

    ASSERT_EQ(R.num_pn_particles, 2);
    ASSERT_EQ(R.num_particles, 3);
    ASSERT_NEAR(R.physical_state->time, tend,
                tend * R.options->output_time_relative_tolerance);

    ASSERT_EQ(R.num_mergers, 2);
    ASSERT_EQ(extra_data[0], 1);
    ASSERT_EQ(extra_data[1], 3);
    ASSERT_EQ(extra_data[2], 5);
    ASSERT_EQ(extra_data[3], 2);
    ASSERT_EQ(extra_data[4], 4);

    ASSERT_EQ(R.merger_infos[0].extra_index1, 1);
    ASSERT_EQ(R.merger_infos[0].extra_index2, 4);
    ASSERT_EQ(R.merger_infos[0].m1, 1);
    ASSERT_EQ(R.merger_infos[0].m2, 2);
    ASSERT_EQ(R.merger_infos[0].m_remnant, R.physical_state->mass[1]);

    ASSERT_EQ(R.merger_infos[1].extra_index1, 0);
    ASSERT_EQ(R.merger_infos[1].extra_index2, 3);
    ASSERT_EQ(R.merger_infos[1].m_remnant, R.physical_state->mass[0]);
    ASSERT_EQ(R.merger_infos[1].m1, 1);
    ASSERT_EQ(R.merger_infos[1].m2, 0.9);

    ketju_free_system(&R);
}
