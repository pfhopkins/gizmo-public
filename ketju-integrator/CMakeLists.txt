cmake_minimum_required(VERSION 3.0 FATAL_ERROR)


# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Set C++ standard (used for tests)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build (Release, RelWithDebInfo, Debug or MinSizeRel)."
      FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Define project
#
project(ketju-integrator C CXX)

set(PROJECT_NAME "Ketju-integrator")
set(PROJECT_COPYRIGHT_YEARS "2019-2022")
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_VERSION_TWEAK -dev)
set(PROJECT_VERSION_STRING "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}${PROJECT_VERSION_TWEAK}")

# Project output paths
#
set(MAINFOLDER ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH "${MAINFOLDER}/bin")
set(LIBRARY_OUTPUT_PATH "${MAINFOLDER}/lib")

# Project Search Paths
#
LIST (APPEND CMAKE_PREFIX_PATH "${MAINFOLDER}")
LIST (APPEND CMAKE_PREFIX_PATH "${MAINFOLDER}/thirdparty")


# get the library sources
set(lib_headers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ketju_integrator/ketju_integrator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ketju_integrator/internal.h
    )
set(lib_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mst.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integration.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pn.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bh_mergers.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/debug.c
    )
source_group("" FILES ${lib_headers} ${lib_sources})



# Create library target and add flags we want to enable always

#
# Generate an exportable library target with the required paths etc.
#
add_library(ketju-integrator ${lib_headers} ${lib_sources})

target_include_directories(ketju-integrator
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

target_compile_options(ketju-integrator
    PRIVATE
        $<$<C_COMPILER_ID:GNU>:
            -Wall -Wextra -Wpedantic -Wno-unknown-pragmas -fPIC
        >
        $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>>:
            -Wall -Wextra -Wpedantic -Wno-missing-braces -fPIC
        >
        $<$<C_COMPILER_ID:Intel>:
            -w3 -diag-disable:remark -fPIC
        >
)

target_link_libraries(ketju-integrator PUBLIC m) # math lib


# Find packages

find_package(MPI REQUIRED)
target_include_directories(ketju-integrator PUBLIC ${MPI_C_INCLUDE_PATH})
target_link_libraries(ketju-integrator PUBLIC ${MPI_C_LIBRARIES})
# Need the linker flags also, target_link_options is only available in very
# new cmake versions, so use this.
# First need to strip extra whitespace in the MPI variable which causes errors
string(STRIP "${MPI_C_LINK_FLAGS}" MPI_C_LINK_FLAGS)
target_link_libraries(ketju-integrator PUBLIC ${MPI_C_LINK_FLAGS})

# standalone executable

add_executable(nbody "nbody/main.c")
target_link_libraries(nbody ketju-integrator)

# Python shim
add_subdirectory("python_bindings")

# Generate install cmake files
#
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${MAINFOLDER}/ketju-integrator-config-version.cmake"
    VERSION ${PROJECT_VERSION_STRING}
    COMPATIBILITY AnyNewerVersion
)

export(TARGETS ketju-integrator
    FILE "${MAINFOLDER}/ketju-integrator-config.cmake"
)

## unit tests if enabled
option(UNIT_TESTS "Build and run unit-tests" OFF)
if(UNIT_TESTS)
    add_subdirectory("unit-tests")
endif()


