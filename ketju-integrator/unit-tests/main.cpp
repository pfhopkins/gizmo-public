// Ketju-integrator tests.
// Copyright (C) 2020-2022 Matias Mannerkoski and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.

#include "ketju_integrator/ketju_integrator.h"
#include <gtest/gtest.h>

#include <mpi.h>
// This handy header provides the basic environment necessary for MPI,
// and a basic MPI aware reporter.
// The output is nicer for non-MPI runs, but at least testing is possible.
#include <gtest-mpi-listener.hpp>

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new MPIEnvironment);

    ::testing::TestEventListeners &listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    delete listeners.Release(listeners.default_result_printer());

    listeners.Append(new MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
