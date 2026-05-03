# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${MAINFOLDER}/thirdparty/googletest
                 googletest-build
                 EXCLUDE_FROM_ALL)

# List sources explicitly to allow proper regeneration of the build when new
# files are added
set(UNIT_TEST_SRCS
    main.cpp
    mst_construction.cpp
    integration.cpp
    mergers.cpp
    )

add_executable(ketju-integrator-unit-tests ${UNIT_TEST_SRCS})
target_link_libraries(ketju-integrator-unit-tests  
                      PUBLIC ketju-integrator gtest)

target_include_directories(ketju-integrator-unit-tests PUBLIC
                     ${MAINFOLDER}/thirdparty/gtest-mpi-listener/)
target_include_directories(ketju-integrator-unit-tests
                           PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(ketju-integrator-unit-tests
                      PUBLIC ${MPI_CXX_LIBRARIES})
target_link_libraries(ketju-integrator-unit-tests
                      PUBLIC ${MPI_CXX_LINK_FLAGS})

set(tests_max_num_mpi_tasks "4" CACHE STRING
    "Maximum number of tasks to use when running unit tests")

# unit tests are not run automatically, only through the run-unit-tests target
add_custom_target(run-unit-tests)
foreach(ntasks 1 2 ${tests_max_num_mpi_tasks})
    add_custom_target(run-unit-tests-${ntasks}
                      COMMAND ${MPIEXEC} -np ${ntasks} --oversubscribe
                      $<TARGET_FILE:ketju-integrator-unit-tests>
                      DEPENDS ketju-integrator-unit-tests
                      WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
                      COMMENT "Running unit tests with ${ntasks} task(s)..."
                      )
    add_dependencies(run-unit-tests run-unit-tests-${ntasks})
endforeach()
