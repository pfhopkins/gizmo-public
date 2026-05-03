// Ketju-integrator tests.
// Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.


#include <gtest/gtest.h>
#include <mpi.h>
#include <random>
#include <sstream>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

// Since this is C++, the struct keyword before ketju_system is not necessary,
// but feels nicer to have it for consistency with the rest of the codebase.

template <typename RNG>
static void init_random_positions(struct ketju_system *R, int N, RNG &&gen) {
    // init N random particle positions. masses etc don't matter since we're
    // only interested in the tree construction here.
    ketju_create_system(R, 0, N, MPI_COMM_WORLD);
    std::uniform_real_distribution<> dist(-1, 1);
    for (int i = 0; i < N; ++i) {
        R->physical_state->mass[i] = 1; // need to be non-zero
        for (int k = 0; k < 3; ++k) {
            R->physical_state->pos[i][k] = dist(gen);
        }
    }
}

static bool connected_to_root(struct ketju_system *R, int i) {
    auto *mst = R->internal_state->mst;
    int n = 0;
    while (mst->vertex[i].level > 0 && n < R->num_particles) {
        i = mst->vertex[i].parent;
        ++n;
        if (i < 0 || i >= R->num_particles) return false;
    }
    return mst->vertex[i].level == 0;
}

class ConstructionTest
    : public testing::TestWithParam<std::tuple<enum ketju_mst_algorithm, int>> {
  public:
    static std::mt19937 gen;
};
std::mt19937 ConstructionTest::gen(0x5eed); // Define and init the static member

TEST_P(ConstructionTest, random_stars) {
    struct ketju_system R;
    auto [algo, N] = GetParam();
    init_random_positions(&R, N, gen);
    R.options->mst_algorithm = KETJU_MST_PRIM;
    ketju_construct_MST(&R);

    auto *mst = R.internal_state->mst;
    int root_count = 0;
    for (int i = 0; i < N; ++i) {
        ASSERT_GE(mst->vertex[i].level, 0);
        ASSERT_LT(mst->vertex[i].level, N);
        ASSERT_TRUE(connected_to_root(&R, i));

        if (mst->vertex[i].level == 0) root_count += 1;
    }
    ASSERT_EQ(root_count, 1);

    // Check edges
    for (int i = 0; i < N - 1; ++i) {
        auto &edge = mst->edge_in_mst[i];
        ASSERT_NE(edge.vertex1, edge.vertex2);
        ASSERT_GE(edge.vertex1, 0);
        ASSERT_LT(edge.vertex1, N);
        ASSERT_GE(edge.vertex2, 0);
        ASSERT_LT(edge.vertex2, N);
        if (mst->vertex[edge.vertex1].level < mst->vertex[edge.vertex2].level) {
            ASSERT_EQ(mst->vertex[edge.vertex2].parent, edge.vertex1);
        } else {
            ASSERT_EQ(mst->vertex[edge.vertex1].parent, edge.vertex2);
        }
    }

    ketju_free_system(&R);
}

INSTANTIATE_TEST_SUITE_P(
    MST, ConstructionTest,
    testing::Combine(testing::Values(KETJU_MST_PRIM,
                                     KETJU_MST_DIVIDE_AND_CONQUER_PRIM),
                     testing::Values(2, 4, 10, 20, 50, 100, 500, 1000)),
    [](const auto &info) {
        enum ketju_mst_algorithm algo = std::get<0>(info.param);
        auto N = std::get<1>(info.param);
        std::stringstream ss;
        switch (algo) {
        case KETJU_MST_PRIM:
            ss << "Prim_";
            break;
        case KETJU_MST_DIVIDE_AND_CONQUER_PRIM:
            ss << "DnC_";
            break;
        default:
            ss << "Unknown_algo_";
        }
        ss << N << "_particles";
        return ss.str();
    });
