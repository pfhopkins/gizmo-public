// Ketju-integrator (MSTAR)
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

// The public interface header with all the user-facing types and function
// prototypes.

#ifndef KETJU_INTEGRATOR_H
#define KETJU_INTEGRATOR_H

#include <mpi.h>
#include <stdbool.h>

// For compatibility with C++ programs
#ifdef __cplusplus
extern "C" {
#endif

// Types and enums

// Using bitflags for the PN terms allows more fine-grained control of the used
// terms (e.g. only PN2.5 radiation reaction or only conservative terms).
// Combine flags with the bitwise or operator |, remove flags from the
// predefined sets using bitwise and with negation, e.g. KETJU_PN_CONSERVATIVE &
// ~KETJU_PN_THREEBODY for conservative terms without three-body cross-terms.
enum ketju_PN_flags {
    KETJU_PN_NONE = 0,
    KETJU_PN_1_0_ACC = 1 << 0,
    KETJU_PN_THREEBODY = 1 << 1,
    KETJU_PN_2_0_ACC = 1 << 2,
    KETJU_PN_2_5_ACC = 1 << 3,
    KETJU_PN_3_0_ACC = 1 << 4,
    KETJU_PN_3_5_ACC = 1 << 5,
    KETJU_PN_CONSERVATIVE = KETJU_PN_1_0_ACC | KETJU_PN_2_0_ACC
                            | KETJU_PN_3_0_ACC | KETJU_PN_THREEBODY,
    // All the spin independent terms.
    // Name might not be the best, but should do.
    KETJU_PN_DYNAMIC_ALL =
        KETJU_PN_CONSERVATIVE | KETJU_PN_2_5_ACC | KETJU_PN_3_5_ACC,

    KETJU_PN_1_5_SPIN_ACC = 1 << 8,
    KETJU_PN_2_0_SPIN_ACC = 1 << 9,
    KETJU_PN_SPIN_EVOL = 1 << 10,
    KETJU_PN_SPIN_ALL =
        KETJU_PN_1_5_SPIN_ACC | KETJU_PN_2_0_SPIN_ACC | KETJU_PN_SPIN_EVOL,

    KETJU_PN_ALL = KETJU_PN_DYNAMIC_ALL | KETJU_PN_SPIN_ALL,
};

enum ketju_mst_algorithm {
    KETJU_MST_PRIM,
    //< Exact MST construction using Prim's algorithm, the default mode suitable
    // for most systems.
    KETJU_MST_DIVIDE_AND_CONQUER_PRIM,
    //< Approximate MST construction by splitting the domain into smaller pieces
    // and then combining the sub-trees.
    // May be faster for very large systems.
};

struct ketju_integration_options {
    int gbs_kmax;
    //< Maximum k (number of rows) the integrator can use.
    // The actually used value is adjusted adaptively during the integration.
    // Values larger than the default of 10 are not generally helpful as double
    // precision arithmetic can't make use of the extremely high order of the
    // integrator (2*k) anyway.
    // The minimum allowed value is 4.

    double gbs_relative_tolerance;
    //< Relative error tolerance per step.
    // This is the main option for controlling the integrator accuracy, with
    // values 1e-7 - 1e-12 working well depending on the problem.
    // Larger values may lead to unacceptable errors in orbits, while smaller
    // can lead to the integrator stalling as floating point errors prevent
    // reaching the desired accuracy, particularly when close to
    // the double precision epsilon at ~1e-16.
    double gbs_absolute_tolerance;
    //< Minimum absolute size of error per step.
    // Can be used to prevent the integrator doing excess work when a variable
    // is near zero, where the relative tolerance condition can result in
    // unnecessarily short timesteps.

    double output_time_relative_tolerance;
    //< Maximum size of relative error in iteration to requested integration end
    // time
    double output_time_absolute_tolerance;
    //< Minimum absolute size of error in time iteration.

    unsigned int PN_flags;
    //< Bitflags for enabling different PN terms, using the values defined in
    // enum ketju_PN_flags.

    int steps_between_mst_reconstruction;
    //< How often to reconstruct the tree coordinates.
    // Constructing the tree is fairly costly, and often the possible accuracy
    // improvements are not big enough to justify doing this that often.
    int max_tree_distance;
    //< Maximum separation in the tree structure before using CoM coordinate
    // distances.
    int max_step_count;
    //< Maximum number of attempted (successful or not) steps before abandoning
    //  the integration and dumping state.
    //  Set to 0 or a negative value to disable.

    enum ketju_mst_algorithm mst_algorithm;

    double star_star_softening;
    //< force softening length between non-pn particles, disabled if <= 0
    bool no_star_star_gravity;
    //< Disable gravity between stars (non-pn particles), only interactions
    // including pn-particles are computed.
    bool enable_bh_mergers;
    //< Check for mergers between BHs (pn particles).
    double bh_merger_distance_in_Rs;
    //< If BH mergers are enabled, the particles are merged when they come
    // within (this factor) * 2G/c^2 * (m1 + m2).
    // The PN forces break down below about 5, and for suppressing PN CoM
    // oscillations values of ~>10 are recommended.
    bool enable_bh_merger_kicks;
    //< Include the merger kick in BH mergers.
    // If disabled the merger remnant moves at the CoM velocity of the merged
    // BHs.
};

// Default values applied when creating the system
// Note that when including this header in C++ code, this only works if the
// compiler extends the standard to support designated initializers,
// which aren't standard C++ before C++20.
// Both GCC (>=8.1) and Clang support this.
// But if necessary, these definitions could be guarded with ifdef __cplusplus
// to hide them from C++ code while allowing internal use.
static const struct ketju_integration_options
    ketju_default_integration_options = {
        .gbs_kmax = 10,
        .gbs_relative_tolerance = 1e-9,
        .gbs_absolute_tolerance = 1e-100,
        .output_time_relative_tolerance = 1e-6,
        .output_time_absolute_tolerance = 1e-300,
        .PN_flags = KETJU_PN_ALL,
        .steps_between_mst_reconstruction = 10,
        .max_tree_distance = 2,
        .max_step_count = -1,
        .mst_algorithm = KETJU_MST_PRIM,
        .star_star_softening = 0.,
        .no_star_star_gravity = false,
        .enable_bh_mergers = true,
        .bh_merger_distance_in_Rs = 12,
        .enable_bh_merger_kicks = true,
};

// The unit system used by the integrator is specified by setting the values of
// these constants.
struct ketju_natural_constants {
    double G; // Gravitational constant
    double c; // Speed of light
};

// Default unit system based on the Gadget-2/3 default unit system:
// Mass in 10^10 Msol (= 1.989e43 g), length in kpc, velocity in cm/s.
// The resulting unit of time is ~ 0.978 Gyr.
// The value of G is based on the value of 6.678e-8 in cgs units.
static const struct ketju_natural_constants ketju_default_natural_constants = {
    .G = 43007.1, .c = 299792.0};

// The physical state of the system.
// The values for the ith particle are stored in the ith element of each array,
// with the num_pn_particles PN-enabled particles (BHs) coming first.
// The order of the particles stays the same during the integration,
// but mergers cause the data to be shifted in order to remove the merged
// particles. If you want to track individual particles over mergers, add IDs
// using the particle_extra_data field of ketju_system.
//
// The origin of the coordinate system can be set freely,
// but it is numerically favorable to have the CoM of the system stay near the
// origin. Note that the PN corrections depend slightly on the coordinate
// velocities so that different choices of the reference frame lead to slightly
// different results. The PN corrected equations of motion are approximately
// Lorentz-invariant, while the standard Newtonian equations are
// Galilei-invariant, so for general systems mixing these particle types there
// is no real 'correct' way to choose the coordinates. The differences caused by
// these effects are very small for reasonable choices of rest frame.
//
// Spins for PN particles (BHs) are stored as the spin angular momentum vector,
// with units of (length)*(mass)*(velocity),
// not to be confused with the spin length scale a or the dimensionless spin
// parameter chi that are often used when discussing black holes.
struct ketju_system_physical_state {
    double time;       // Physical time of the system.
    double *mass;      // array of length num_particles
    double (*pos)[3];  // particle positions: array of shape [num_particles][3]
    double (*vel)[3];  // particle velocities: shape [num_particles][3]
    double (*spin)[3]; // PN particle spins: shape [num_pn_particles][3]
};

struct ketju_performance_counters {
    // Timers for different parts of the code
    double time_force;
    double time_force_comm;
    double time_coord;
    double time_mst;
    double time_gbscomm;
    double time_gbsextr;
    double time_total;

    // Sub-timers for DnC Prim, time included in time_mst also.
    double tprim1, tprim2, tprim3, tprim4, tprim5;

    int successful_steps; // Number of accepted steps
    int failed_steps;     // Number of failed steps
    double total_work;    //< ~ Number of force evaluations done
    double relative_energy_error;
    //< Relative change in the energy of the system over the last
    //  ketju_run_integrator call.
    //  Gives information on accuracy only when PN terms are disabled.
};

struct ketju_bh_merger_info {
    double t_merger;                // time of the merger
    double m1, m2, m_remnant;       // masses
    double chi1, chi2, chi_remnant; // dimensionless spin magnitudes
    double v_kick;                  // magnitude of kick velocity
    int extra_index1,
        extra_index2; // Indices of the components in particle_extra_data
};

struct ketju_system {
    int num_particles;    //< Total number of particles
    int num_pn_particles; //< Number of PN enabled particles

    // The struct members apart from internal_state could also be direct members
    // instead of pointers, but the performance implications aren't likely that
    // big either way so this is mostly a stylistic choice.

    struct ketju_system_physical_state *physical_state;
    //< Physical state of the system

    struct ketju_system_internal_state *internal_state;
    //< Internal state (implementation detail)

    struct ketju_performance_counters *perf;
    struct ketju_integration_options *options;

    struct ketju_natural_constants *constants;
    //< Values of natural constants which define the unit system

    void *particle_extra_data;
    //< Optional extra data (e.g. IDs) associated with the particles, which is
    // reordered in case of mergers so that the removed particle data are found
    // at the end of the array.
    // If you set this pointer, you must also set particle_extra_data_elem_size
    // so that the reordering of the data during mergers is done correctly.
    // The caller is responsible for managing this memory.

    size_t particle_extra_data_elem_size;
    //< Size of the extra data for a single particle.

    int num_mergers; //< How many mergers there have been
    struct ketju_bh_merger_info *merger_infos; //< contains num_mergers elements
};

// Function prototypes

// Allocate and set default options and physical constants.
// Parameters:
//  R: pointer to an unitialized ketju_system struct,
//     which this function initializes.
//  num_pn_particles, num_other_particles: counts of PN enabled particles (BHs)
//     and purely Newtonian particles (stars) to allocate space for in
//     physical_state.
//  comm: An MPI communicator from which the integrator communicators are
//        created.
// All tasks in the specified comm must call this function with the same
// particle counts, then set the physical_state, options, and constants contents
// to the same values, and finally call ketju_run_integrator with the same
// arguments.
void ketju_create_system(struct ketju_system *R, int num_pn_particles,
                         int num_other_particles, MPI_Comm comm);
// Free the allocated system storage
void ketju_free_system(struct ketju_system *R);

// Compute the centre of mass of the system, puts results in CoM_pos and CoM_vel
void ketju_compute_CoM(const struct ketju_system *R, double CoM_pos[3],
                       double CoM_vel[3]);
// Shift the reference frame origin of the system to the CoM.
// Outputs the CoM coordinates in the original frame in CoM_pos and CoM_vel.
void ketju_into_CoM_frame(struct ketju_system *R, double CoM_pos[3],
                          double CoM_vel[3]);
// Runs the integrator for a given time interval.
// The time interval can be negative for backwards-in-time integration.
void ketju_run_integrator(struct ketju_system *R, double time_interval);

// Debug helpers for dumping the state to a file and reading it later for
// analysis. The state is dumped in case of some fatal errors, such as when the
// potential energy of the system is not finite.
void ketju_debug_dump_integrator_state(const struct ketju_system *R,
                                       const char *fname);
bool ketju_debug_read_integrator_state(struct ketju_system *R,
                                       const char *fname, MPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif // KETJU_INTEGRATOR_H
