// Ketju-integrator (MSTAR)
// Copyright (C) 2020-2022 Antti Rantala, Matias Mannerkoski, and contributors.
// Licensed under the GPLv3 license.
// See LICENSE.md or <https://www.gnu.org/licenses/> for details.

// Main integration and force calculation functions

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "ketju_integrator/internal.h"
#include "ketju_integrator/ketju_integrator.h"

// This is needed only in this file
enum velocity_type { PHYSICAL, AUXILIARY };

/////////////////////////////////////////
// Acceleration and energy calculation //
/////////////////////////////////////////

// Gadget softening kernels
static inline double softened_inverse_r(double r, double h) {
    if (h <= 0) return 1. / r;

    const double hinv = 1. / h;
    const double u = r * hinv;

    if (u >= 1.) {
        return 1. / r;
    } else if (u < 0.5) {
        return -hinv * (-2.8 + u * u * (16. / 3. + u * u * (6.4 * u - 9.6)));
    } else {
        return -hinv
               * (-3.2 + 2. / 30. / u
                  + u * u
                        * (32. / 3. + u * (-16.0 + u * (9.6 - 32. / 15. * u))));
    }
}

static inline double softened_inverse_r3(double r, double h) {
    if (h <= 0) return 1. / (r * r * r);

    const double hinv = 1. / h;
    const double u = r * hinv;

    if (u >= 1) {
        return 1.0 / (r * r * r);
    } else if (u < 0.5) {
        return hinv * hinv * hinv * (32. / 3. + u * u * (32.0 * u - 38.4));
    } else {
        return hinv * hinv * hinv
               * (64. / 3. - 48.0 * u + 38.4 * u * u - 32. / 3. * u * u * u
                  - 2. / 30. / (u * u * u));
    }
}

// This function wraps the distance calculation that would be otherwise repeated
// in the below functions.
// It takes all of the necessary working space variables as arguments,
// so that they don't get created within the loop where this is called from,
// which affects performance quite a bit.
// After the compiler inlines this, the performance is the same as with repeated
// code within the loops (used earlier in the implementation).
static inline double separation_and_r2(const struct ketju_system *R, int i,
                                       int j, double dr[restrict 3], int Nd,
                                       int path[restrict Nd],
                                       int sign[restrict Nd], int *restrict d) {
    double r2 = 0;
    for (int k = 0; k < 3; k++) {
        dr[k] = 0;
    }

    const struct graph_vertex *const vertex = R->internal_state->mst->vertex;

    if (abs(vertex[i].level - vertex[j].level) <= Nd
        && ketju_check_relative_proximity(j, i, Nd, R, d, path, sign)) {
        for (int l = 0; l < *d; l++) {
            const int s = sign[l];
            const double *const rel_pos =
                R->internal_state->mst->edge_rel_pos[path[l]];
            for (int k = 0; k < 3; k++) {
                dr[k] += s * rel_pos[k];
            }
        }
    } else {
        for (int k = 0; k < 3; k++) {
            dr[k] = R->physical_state->pos[j][k] - R->physical_state->pos[i][k];
        }
    }
    for (int k = 0; k < 3; k++) {
        r2 += dr[k] * dr[k];
    }
    return r2;
}

static void compute_U(struct ketju_system *R) {

    timer_start(&R->perf->time_force);

    struct ketju_system_internal_state *const is = R->internal_state;
    struct ketju_system_physical_state *const ps = R->physical_state;

    is->U = 0;
    int istart = 0;
    int cum_num, block_size, jstart;
    loop_scheduling_N2(R->num_particles, is->task_info->gbs_group_size,
                       is->task_info->gbs_group_rank, &istart, &jstart,
                       &block_size, &cum_num);
    int c = 0;

    const int Nd = R->options->max_tree_distance;
    const double ss_soft = R->options->star_star_softening;

    double dr[3];
    int path[Nd];
    int sign[Nd];
    int d;

    for (int i = istart; i < R->num_particles; i++) {
        double mi = ps->mass[i];
        if (i > istart) { jstart = i + 1; }

        for (int j = jstart; j < R->num_particles; j++) {
            if (c >= block_size) { goto after_loops; }
            if (R->options->no_star_star_gravity && i >= R->num_pn_particles
                && j >= R->num_pn_particles) {
                continue;
            }

            const double r2 =
                separation_and_r2(R, i, j, dr, Nd, path, sign, &d);

            const double soft =
                i < R->num_pn_particles || j < R->num_pn_particles ? 0.
                                                                   : ss_soft;
            const double invr = softened_inverse_r(sqrt(r2), soft);
            const double mj = ps->mass[j];

            is->U += mi * mj * invr;
            c++;
        }
    }
after_loops:
    timer_stop(&R->perf->time_force);

    timer_start(&R->perf->time_force_comm);
    if (is->task_info->gbs_group_size > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &is->U, 1, MPI_DOUBLE, MPI_SUM,
                      is->task_info->gbs_comm);
    }
    timer_stop(&R->perf->time_force_comm);

    if (!isfinite(is->U)) {
        if (R->internal_state->task_info->main_group_rank == 0) {
            ketju_debug_dump_integrator_state(
                R, "ketju_integrator_infinite_U_state_dump.bin");
        }
        terminate("U is not finite! State dumped to file.");
    }

    is->U *= R->constants->G;
}

// To support optional softening and ignoring of star-star interactions we need
// to have slightly different force loops depending on the options.
// The force loop is the most time consuming part of the integration, and any
// extra logic in the loop leads to significant performance penalties (~10-20%).
// So instead have a different function for each case.
static void comp_U_acc_loop_standard(struct ketju_system *R) {
    struct ketju_system_internal_state *const is = R->internal_state;
    struct ketju_system_physical_state *const ps = R->physical_state;

    const double G = R->constants->G;

    const int istop = R->num_particles;
    const int jstop = R->num_particles;
    int block_size, cum_num, jstart, istart;
    loop_scheduling_N2(R->num_particles, is->task_info->gbs_group_size,
                       is->task_info->gbs_group_rank, &istart, &jstart,
                       &block_size, &cum_num);
    int c = 0;
    const int Nd = R->options->max_tree_distance;

    double dr[3];
    int path[Nd];
    int sign[Nd];
    int d;

    for (int i = istart; i < istop; i++) {
        const double mi = ps->mass[i];
        if (i > istart) jstart = i + 1;

        for (int j = jstart; j < jstop; j++) {
            if (c >= block_size) { return; }

            const double r2 =
                separation_and_r2(R, i, j, dr, Nd, path, sign, &d);
            const double invr = 1. / sqrt(r2);
            const double Ginvr3 = G * invr * invr * invr;
            const double mj = ps->mass[j];

            is->U += G * mi * mj * invr;

            for (int k = 0; k < 3; k++) {
                double drGinvr3 = dr[k] * Ginvr3;
                is->acc[i][k] += mj * drGinvr3;
                is->acc[j][k] += -mi * drGinvr3;
            }
            c++;
        }
    }
}

// Force without star-star interactions
static void comp_U_acc_loop_no_star_star(struct ketju_system *R) {
    struct ketju_system_internal_state *const is = R->internal_state;
    struct ketju_system_physical_state *const ps = R->physical_state;

    const double G = R->constants->G;

    // Assume there are only a few pn particles/BHs and parallelize only
    // over the inner loop.
    const int istart = 0, istop = R->num_pn_particles;

    const int rank = is->task_info->gbs_group_rank;
    int group_size = is->task_info->gbs_group_size;

    if (group_size > R->num_particles) {
        // in case we for some reason have more tasks than particles
        group_size = R->num_particles;
        if (rank >= group_size) return;
    }

    const int part_per_task = R->num_particles / group_size;
    const int part_per_task_remainder = R->num_particles % group_size;
    int jstart = rank * part_per_task;
    if (rank < part_per_task_remainder) {
        jstart += rank; // the non-evenly split particles go to the first tasks
    } else {
        jstart += part_per_task_remainder;
    }

    const int jstop =
        jstart + part_per_task + (rank < part_per_task_remainder ? 1 : 0);

    const int Nd = R->options->max_tree_distance;

    double dr[3];
    int path[Nd];
    int sign[Nd];
    int d;

    for (int i = istart; i < istop; i++) {
        const double mi = ps->mass[i];

        if (jstart <= i) jstart = i + 1; // avoid double counting interactions
        for (int j = jstart; j < jstop; j++) {
            const double r2 =
                separation_and_r2(R, i, j, dr, Nd, path, sign, &d);
            const double invr = 1. / sqrt(r2);
            const double Ginvr3 = G * invr * invr * invr;
            const double mj = ps->mass[j];

            is->U += G * mi * mj * invr;

            for (int k = 0; k < 3; k++) {
                double drGinvr3 = dr[k] * Ginvr3;
                is->acc[i][k] += mj * drGinvr3;
                is->acc[j][k] += -mi * drGinvr3;
            }
        }
    }
}

// Force with softened star-star interactions, doesn't include the pn-particles
static void comp_U_acc_loop_softened_star_star(struct ketju_system *R) {
    struct ketju_system_internal_state *const is = R->internal_state;
    struct ketju_system_physical_state *const ps = R->physical_state;

    const int istop = R->num_particles;
    const int jstop = R->num_particles;
    int block_size, cum_num, jstart, istart;
    loop_scheduling_N2(
        R->num_particles - R->num_pn_particles, is->task_info->gbs_group_size,
        is->task_info->gbs_group_rank, &istart, &jstart, &block_size, &cum_num);
    istart += R->num_pn_particles;
    jstart += R->num_pn_particles;
    int c = 0;

    const int Nd = R->options->max_tree_distance;
    const double ss_soft = R->options->star_star_softening;
    const double G = R->constants->G;

    double dr[3];
    int path[Nd];
    int sign[Nd];
    int d;

    for (int i = istart; i < istop; i++) {
        const double mi = ps->mass[i];
        if (i > istart) jstart = i + 1;

        for (int j = jstart; j < jstop; j++) {
            if (c >= block_size) { return; }

            const double r2 =
                separation_and_r2(R, i, j, dr, Nd, path, sign, &d);
            const double r = sqrt(r2);
            const double invr = softened_inverse_r(r, ss_soft);
            const double Ginvr3 = G * softened_inverse_r3(r, ss_soft);
            const double mj = ps->mass[j];

            is->U += G * mi * mj * invr;

            for (int k = 0; k < 3; k++) {
                double drGinvr3 = dr[k] * Ginvr3;
                is->acc[i][k] += mj * drGinvr3;
                is->acc[j][k] += -mi * drGinvr3;
            }
            c++;
        }
    }
}

// computes the force function and newtonian acceleration
static void compute_U_and_Newtonian_Acc(struct ketju_system *R) {

    timer_start(&R->perf->time_force);

    struct ketju_system_internal_state *const is = R->internal_state;
    struct spanning_tree_coordinate_system *const mst = is->mst;

    memset(is->acc, 0, sizeof(double[R->num_particles][3]));
    is->U = 0;

    // select the force loop based on options
    if (R->options->no_star_star_gravity) {
        comp_U_acc_loop_no_star_star(R);
    } else if (R->options->star_star_softening > 0) {
        comp_U_acc_loop_no_star_star(R);       // non-softened part
        comp_U_acc_loop_softened_star_star(R); // softened part
    } else {
        comp_U_acc_loop_standard(R);
    }

    timer_stop(&R->perf->time_force);

    timer_start(&R->perf->time_force_comm);
    if (is->task_info->gbs_group_size > 1) {
        MPI_Allreduce(MPI_IN_PLACE, is->acc, 3 * R->num_particles, MPI_DOUBLE,
                      MPI_SUM, is->task_info->gbs_comm);
        MPI_Allreduce(MPI_IN_PLACE, &is->U, 1, MPI_DOUBLE, MPI_SUM,
                      is->task_info->gbs_comm);
    }
    timer_stop(&R->perf->time_force_comm);

    timer_start(&R->perf->time_force);
    int hi, lo;
    for (int i = 0; i < R->num_particles - 1; i++) {

        int v1, v2;
        get_v1_v2_fast(R, i, &v1, &v2);

        if (mst->vertex[v1].level > mst->vertex[v2].level) {
            hi = v1;
            lo = v2;
        } else {
            hi = v2;
            lo = v1;
        }
        for (int k = 0; k < 3; k++) {
            mst->edge_rel_acc[i][k] = is->acc[hi][k] - is->acc[lo][k];
        }
    }

    timer_stop(&R->perf->time_force);
}

static void compute_T(struct ketju_system *R) {
    struct ketju_system_internal_state *const is = R->internal_state;
    struct ketju_system_physical_state *const ps = R->physical_state;
    is->T = 0;

    for (int i = 0; i < R->num_particles; i++) {
        double dv, v2 = 0;
        for (int k = 0; k < 3; k++) {
            dv = ps->vel[i][k];
            v2 += dv * dv;
        }
        is->T += 0.5 * ps->mass[i] * v2;
    }
}

//////////////////////////////
// Some MST related helpers //
//////////////////////////////

// The following few functions could maybe go in mst.c, having them here reduces
// the number of externally visible symbols, and should allow for better
// inlining.
static int get_root_index(const struct ketju_system *R) {
    const struct spanning_tree_coordinate_system *const mst =
        R->internal_state->mst;

    int v1, v2;
    get_v1_v2_fast(R, 0, &v1, &v2);

    return mst->vertex[v1].level == 0 ? v1 : v2;
}

static void into_Cartesian_from_MST(struct ketju_system *R) {

    timer_start(&R->perf->time_coord);

    const struct spanning_tree_coordinate_system *const mst =
        R->internal_state->mst;
    struct ketju_system_physical_state *const ps = R->physical_state;

    int first = get_root_index(R);
    for (int k = 0; k < 3; k++) {
        ps->pos[first][k] = mst->root_pos[k];
        ps->vel[first][k] = mst->root_vel[k];
    }
    for (int i = 0; i < R->num_particles - 1; i++) {

        int v1, v2;
        get_v1_v2_fast(R, i, &v1, &v2);

        int lo, hi;
        if (mst->vertex[v1].level > mst->vertex[v2].level) {
            lo = v2;
            hi = v1;
        } else {
            lo = v1;
            hi = v2;
        }
        for (int k = 0; k < 3; k++) {
            ps->pos[hi][k] = ps->pos[lo][k] + mst->edge_rel_pos[i][k];
            ps->vel[hi][k] = ps->vel[lo][k] + mst->edge_rel_vel[i][k];
        }
    }

    timer_stop(&R->perf->time_coord);
}

static void update_pos(struct ketju_system *R) {

    timer_start(&R->perf->time_coord);

    const struct spanning_tree_coordinate_system *const mst =
        R->internal_state->mst;
    struct ketju_system_physical_state *const ps = R->physical_state;

    int first = get_root_index(R);
    for (int k = 0; k < 3; k++) {
        ps->pos[first][k] = mst->root_pos[k];
    }
    for (int i = 0; i < R->num_particles - 1; i++) {

        int v1, v2;
        get_v1_v2_fast(R, i, &v1, &v2);

        int lo, hi;
        if (mst->vertex[v1].level > mst->vertex[v2].level) {
            lo = v2;
            hi = v1;
        } else {
            lo = v1;
            hi = v2;
        }
        for (int k = 0; k < 3; k++) {
            ps->pos[hi][k] = ps->pos[lo][k] + mst->edge_rel_pos[i][k];
        }
    }

    timer_stop(&R->perf->time_coord);
}

static void update_vel(struct ketju_system *R,
                       enum velocity_type propagated_velocity) {

    timer_start(&R->perf->time_coord);
    const struct spanning_tree_coordinate_system *const mst =
        R->internal_state->mst;

    double(*vel_to_update)[3] = propagated_velocity == PHYSICAL
                                    ? R->physical_state->vel
                                    : R->internal_state->aux_vel;

    int first = get_root_index(R);
    for (int k = 0; k < 3; k++) {
        vel_to_update[first][k] = propagated_velocity == PHYSICAL
                                      ? mst->root_vel[k]
                                      : mst->root_aux_vel[k];
    }
    for (int i = 0; i < R->num_particles - 1; i++) {

        int v1, v2;
        get_v1_v2_fast(R, i, &v1, &v2);

        int lo, hi;
        if (mst->vertex[v1].level > mst->vertex[v2].level) {
            lo = v2;
            hi = v1;
        } else {
            lo = v1;
            hi = v2;
        }
        for (int k = 0; k < 3; k++) {
            if (propagated_velocity == PHYSICAL) {
                vel_to_update[hi][k] =
                    vel_to_update[lo][k] + mst->edge_rel_vel[i][k];
            } else {
                vel_to_update[hi][k] =
                    vel_to_update[lo][k] + mst->edge_rel_aux_vel[i][k];
            }
        }
    }

    timer_stop(&R->perf->time_coord);
}

//////////////////////////
// Leapfrog integration //
//////////////////////////

static void drift(double ds, struct ketju_system *R) {

    struct ketju_system_internal_state *const is = R->internal_state;
    struct spanning_tree_coordinate_system *const mst = is->mst;
    struct ketju_system_physical_state *const ps = R->physical_state;

    compute_T(R);
    // T + B = U, but depends only on velocities
    double dt = ds / (is->T + is->B);
    ps->time += dt;

    for (int i = 0; i < R->num_particles - 1; i++) {
        for (int k = 0; k < 3; k++) {
            mst->edge_rel_pos[i][k] += dt * mst->edge_rel_vel[i][k];
        }
    }
    for (int k = 0; k < 3; k++) {
        mst->root_pos[k] += dt * mst->root_vel[k];
    }
    update_pos(R);
}

static void set_auxiliary_variables(struct ketju_system *R) {

    struct ketju_system_internal_state *const is = R->internal_state;
    struct spanning_tree_coordinate_system *const mst = is->mst;
    struct ketju_system_physical_state *const ps = R->physical_state;

    memcpy(is->aux_vel, ps->vel, sizeof(double[R->num_particles][3]));
    memcpy(mst->edge_rel_aux_vel, mst->edge_rel_vel,
           sizeof(double[R->num_particles - 1][3]));
    memcpy(mst->root_aux_vel, mst->root_vel, sizeof(double[3]));
    memcpy(is->aux_spin, ps->spin, sizeof(double[R->num_pn_particles][3]));
}

static void physical_kick(double dt, struct ketju_system *R) {

    struct ketju_system_internal_state *const is = R->internal_state;
    struct spanning_tree_coordinate_system *const mst = is->mst;
    struct ketju_system_physical_state *const ps = R->physical_state;

    update_vel(R, AUXILIARY);
    ketju_compute_PN_acc(R, is->aux_vel, is->aux_spin);

    double dB = 0;
    for (int i = 0; i < R->num_particles; i++) {
        for (int k = 0; k < 3; k++) {
            dB += ps->mass[i] * is->aux_vel[i][k] * is->acc_vd[i][k];
        }
    }
    is->B -= dB * dt;

    for (int i = 0; i < R->num_particles - 1; i++) {
        for (int k = 0; k < 3; k++) {
            mst->edge_rel_vel[i][k] +=
                dt * (mst->edge_rel_acc[i][k] + mst->edge_rel_acc_vd[i][k]);
        }
    }

    int root = get_root_index(R);
    for (int k = 0; k < 3; k++) {
        mst->root_vel[k] += dt * (is->acc[root][k] + is->acc_vd[root][k]);
    }

    for (int i = 0; i < R->num_pn_particles; ++i) {
        for (int k = 0; k < 3; k++) {
            ps->spin[i][k] += dt * is->spin_derivative[i][k];
        }
    }
}

static void auxiliary_kick(double dt, struct ketju_system *R) {

    struct ketju_system_internal_state *const is = R->internal_state;
    struct spanning_tree_coordinate_system *const mst = is->mst;
    struct ketju_system_physical_state *const ps = R->physical_state;

    update_vel(R, PHYSICAL);
    ketju_compute_PN_acc(R, ps->vel, ps->spin);

    for (int i = 0; i < R->num_particles - 1; i++) {
        for (int k = 0; k < 3; k++) {
            mst->edge_rel_aux_vel[i][k] +=
                dt * (mst->edge_rel_acc[i][k] + mst->edge_rel_acc_vd[i][k]);
        }
    }

    int root = get_root_index(R);
    for (int k = 0; k < 3; k++) {
        mst->root_aux_vel[k] += dt * (is->acc[root][k] + is->acc_vd[root][k]);
    }

    for (int i = 0; i < R->num_pn_particles; ++i) {
        for (int k = 0; k < 3; k++) {
            is->aux_spin[i][k] += dt * is->spin_derivative[i][k];
        }
    }
}

static void kick(double ds, struct ketju_system *R) {

    struct ketju_system_internal_state *const is = R->internal_state;
    struct spanning_tree_coordinate_system *const mst = is->mst;

    compute_U_and_Newtonian_Acc(R);
    double dt = ds / is->U;

    if (!R->options->PN_flags) {
        // no velocity dependent forces (only PN for now), just the newtonian
        // part needed
        for (int i = 0; i < R->num_particles - 1; i++) {
            for (int k = 0; k < 3; k++) {
                mst->edge_rel_vel[i][k] += dt * mst->edge_rel_acc[i][k];
            }
        }
        int root = get_root_index(R);
        for (int k = 0; k < 3; k++) {
            mst->root_vel[k] += dt * is->acc[root][k];
        }
    } else {
        auxiliary_kick(dt / 2, R);
        physical_kick(dt, R);
        auxiliary_kick(dt / 2, R);
    }
    update_vel(R, PHYSICAL);
}

static void compute_total_energy(struct ketju_system *R) {

    struct ketju_system_internal_state *const is = R->internal_state;

    compute_T(R);
    compute_U(R);
    is->B = is->U - is->T;
}

static void mst_leapfrog(struct ketju_system *R, double Hstep,
                         int NumberOfSubsteps) {

    if (R->options->PN_flags) { set_auxiliary_variables(R); }

    double h = Hstep / NumberOfSubsteps;

    drift(h / 2, R);
    for (int step = 0; step < NumberOfSubsteps - 1; step++) {
        kick(h, R);
        drift(h, R);
    }
    kick(h, R);
    drift(h / 2, R);
}

/////////////////////////////
// GBS style extrapolation //
/////////////////////////////

// Many of the theoretical details of to the extrapolation algorithm are
// explained in e.g. Hairer et al. "Solving Ordinary differential Equations I"
// (DOI 10.1007/978-3-540-78862-1), particularly chapter II.9.

static int num_gbs_substeps(int row) { return 2 * row; }
// sum of the number of function evalutations for a given number of rows k
static int gbs_total_cost(int k) {
    // Cache values that are likely to be needed
    static int values[16];
    static bool inited = false;
    if (!inited) {
        int sum = 1;
        for (int i = 1; i <= 16; ++i) {
            sum = sum + num_gbs_substeps(i);
            values[i - 1] = sum;
        }
        inited = true;
    }

    if (k < 1) { return 0; }
    if (k <= 16) { return values[k - 1]; }

    // Handle the unlikely larger values
    int sum = values[15];
    for (int i = 17; i <= k; ++i) {
        sum = sum + num_gbs_substeps(i);
    }
    return sum;
}

static int argmin_int(int *arr, int size) {
    int min = arr[0];
    int ind = 0;
    for (int i = 1; i < size; ++i) {
        if (arr[i] < min) {
            min = arr[i];
            ind = i;
        }
    }
    return ind;
}

// Get the values of k that are done by the gbs groups.
// The values are in descending order
// Using 1-based counting for rows, so that empty fields
// in the table are just 0 and are directly falsy.
// This is also the convention in num_gbs_substeps(row) above.
static void get_gbs_rows_for_groups(struct ketju_system *R,
                                    int *gbs_rows_for_groups) {
    const int k = R->internal_state->gbs_k;
    const int num_groups = R->internal_state->group_info->num_gbs_groups;
    memset(gbs_rows_for_groups, 0, sizeof(int[num_groups][k]));
    int(*gbs_rows)[k] = (int(*)[k])gbs_rows_for_groups;

    int group_workloads[num_groups];
    memset(group_workloads, 0, sizeof group_workloads);
    int row_ind[num_groups];
    memset(row_ind, 0, sizeof row_ind);

    for (int row = k; row > 0; --row) {
        int cost = num_gbs_substeps(row);
        int group_ind = argmin_int(group_workloads, num_groups);
        group_workloads[group_ind] += cost;
        gbs_rows[group_ind][row_ind[group_ind]++] = row;
    }
}

struct state_backup {
    double time, B, U, T;
    double root_pos[3];
    double root_vel[3];
    double (*pos)[3];
    double (*vel)[3];
    double (*spin)[3];
    double (*edge_rel_pos)[3];
    double (*edge_rel_vel)[3];
};

static void allocate_state_backup(struct state_backup *sb, int num_particles,
                                  int num_pn_particles) {
    sb->pos = malloc(sizeof(double[num_particles][3]));
    sb->vel = malloc(sizeof(double[num_particles][3]));
    sb->spin = malloc(sizeof(double[num_pn_particles][3]));
    sb->edge_rel_pos = malloc(sizeof(double[num_particles - 1][3]));
    sb->edge_rel_vel = malloc(sizeof(double[num_particles - 1][3]));
}

static void store_state_backup(struct state_backup *sb,
                               const struct ketju_system *R) {
    sb->time = R->physical_state->time;
    sb->B = R->internal_state->B;
    sb->U = R->internal_state->U;
    sb->T = R->internal_state->T;

    memcpy(sb->root_pos, R->internal_state->mst->root_pos, sizeof(double[3]));
    memcpy(sb->root_vel, R->internal_state->mst->root_vel, sizeof(double[3]));

    const int num_particles = R->num_particles;
    memcpy(sb->pos, R->physical_state->pos, sizeof(double[num_particles][3]));
    memcpy(sb->vel, R->physical_state->vel, sizeof(double[num_particles][3]));
    memcpy(sb->spin, R->physical_state->spin,
           sizeof(double[R->num_pn_particles][3]));
    memcpy(sb->edge_rel_pos, R->internal_state->mst->edge_rel_pos,
           sizeof(double[num_particles - 1][3]));
    memcpy(sb->edge_rel_vel, R->internal_state->mst->edge_rel_vel,
           sizeof(double[num_particles - 1][3]));
}

static void restore_state_backup(struct state_backup *sb,
                                 const struct ketju_system *R) {
    R->physical_state->time = sb->time;
    R->internal_state->B = sb->B;
    R->internal_state->U = sb->U;
    R->internal_state->T = sb->T;

    memcpy(R->internal_state->mst->root_pos, sb->root_pos, sizeof(double[3]));
    memcpy(R->internal_state->mst->root_vel, sb->root_vel, sizeof(double[3]));

    const int num_particles = R->num_particles;
    memcpy(R->physical_state->pos, sb->pos, sizeof(double[num_particles][3]));
    memcpy(R->physical_state->vel, sb->vel, sizeof(double[num_particles][3]));
    memcpy(R->physical_state->spin, sb->spin,
           sizeof(double[R->num_pn_particles][3]));
    memcpy(R->internal_state->mst->edge_rel_pos, sb->edge_rel_pos,
           sizeof(double[num_particles - 1][3]));
    memcpy(R->internal_state->mst->edge_rel_vel, sb->edge_rel_vel,
           sizeof(double[num_particles - 1][3]));
}

static void free_state_backup(struct state_backup *sb) {
    free(sb->pos);
    free(sb->vel);
    free(sb->spin);
    free(sb->edge_rel_pos);
    free(sb->edge_rel_vel);
}

struct gbs_storage {
    int num_vars;
    int num_rows;
    double *data;
};

static void allocate_gbs_storage(struct gbs_storage *gs, int num_particles,
                                 int num_pn_particles, int k) {
    gs->num_vars = (1 /*time*/ + 3 /*root_pos*/ + 3 /*root_vel*/
                    + (num_particles - 1) * (3 /*edge_pos*/ + 3 /*edge_vel*/)
                    + num_pn_particles * (3 /*spin*/));
    gs->num_rows = k + 1; // row for extrapolation result also
    gs->data = calloc(1, sizeof(double[gs->num_rows][gs->num_vars]));
}

static void free_gbs_storage(struct gbs_storage *gs) { free(gs->data); }

// store the data for this row for extrapolation
static void store_gbs_data(struct gbs_storage *gs, const struct ketju_system *R,
                           int row) {
    double(*data)[gs->num_vars] = (double(*)[gs->num_vars])gs->data;
    const int row_ind = row - 1;
    data[row_ind][0] = R->physical_state->time;
    int offset = 1;
    const struct spanning_tree_coordinate_system *const mst =
        R->internal_state->mst;
    memcpy(&data[row_ind][offset], mst->root_pos, sizeof(double[3]));
    offset += 3;
    memcpy(&data[row_ind][offset], mst->root_vel, sizeof(double[3]));
    offset += 3;
    const size_t edge_vec_size = sizeof(double[R->num_particles - 1][3]);
    memcpy(&data[row_ind][offset], mst->edge_rel_pos, edge_vec_size);
    offset += 3 * (R->num_particles - 1);
    memcpy(&data[row_ind][offset], mst->edge_rel_vel, edge_vec_size);
    offset += 3 * (R->num_particles - 1);
    memcpy(&data[row_ind][offset], R->physical_state->spin,
           sizeof(double[R->num_pn_particles][3]));
}

// Make the non-blocking receiving calls for the rows of data that need to be
// received. Called for all tasks
static void post_gbs_receives(MPI_Request *requests, int *request_index,
                              struct gbs_storage *gs, int *gbs_rows_for_groups,
                              const struct ketju_system *R) {
    timer_start(&R->perf->time_gbscomm);

    const struct task_info *ti = R->internal_state->task_info;
    const int k = R->internal_state->gbs_k;
    double(*data)[gs->num_vars] = (double(*)[gs->num_vars])gs->data;
    int(*gbs_rows)[k] = (int(*)[k])gbs_rows_for_groups;
    for (int i = 0; i < R->internal_state->group_info->num_gbs_groups; ++i) {
        // This task will be computing the data, so doesn't receive
        if (i == ti->gbs_group_index) continue;

        for (int j = 0; j < k && gbs_rows[i][j] > 0; ++j) {
            int data_row = gbs_rows[i][j] - 1;
            MPI_Ibcast(data[data_row], gs->num_vars, MPI_DOUBLE, 0,
                       ti->gbs_data_comms[i], &requests[(*request_index)++]);
        }
    }
    timer_stop(&R->perf->time_gbscomm);
}

// Make the non-blocking send call for a row of data corresponding to the given
// gbs k. Only does anything on the root task of the gbs group,
// but safe to call for all tasks in the gbs group handling the gbs k substeps.
static void post_gbs_send(MPI_Request *requests, int *request_index,
                          struct gbs_storage *gs, const struct ketju_system *R,
                          int k) {

    const struct task_info *ti = R->internal_state->task_info;
    double(*data)[gs->num_vars] = (double(*)[gs->num_vars])gs->data;
    const int i = ti->gbs_group_index;
    if (ti->gbs_data_comms[i] == MPI_COMM_NULL) return;

    timer_start(&R->perf->time_gbscomm);
    int data_row = k - 1;
    MPI_Ibcast(data[data_row], gs->num_vars, MPI_DOUBLE, 0,
               ti->gbs_data_comms[i], &requests[(*request_index)++]);
    timer_stop(&R->perf->time_gbscomm);
}

// Wait for the communications to finish
static void finish_gbs_comm(MPI_Request *requests, int request_index,
                            const struct ketju_system *R) {
    timer_start(&R->perf->time_gbscomm);
    if (request_index > 0) {
        MPI_Waitall(request_index, requests, MPI_STATUSES_IGNORE);
    }
    timer_stop(&R->perf->time_gbscomm);
}

static void extrapolate_single_variable(double *y, double yscal, const int N,
                                        double *result, double *error,
                                        double *error_lower_order) {
    // Aitken--Neville extrapolation, e.g. ch. II.9 eq. (9.10) of the Hairer et
    // al. book. Index labels and their order differ here, and obviously
    // indexing is 0-based.
    double T[N][N];
    for (int i = 0; i < N; i++) {
        T[0][i] = y[i];
    }
    for (int j = 0; j < N - 1; j++) {
        for (int i = j + 1; i < N; i++) {
            const double nterm = pow((double)num_gbs_substeps(i + 1)
                                         / (double)num_gbs_substeps(i - j),
                                     2);
            // Simplified form of (9.10) which seems to behave better
            // numerically.
            T[j + 1][i] = (nterm * T[j][i] - T[j][i - 1]) / (nterm - 1.0);
        }
    }
    *result = T[N - 1][N - 1];
    *error = fabs(T[N - 1][N - 1] - T[N - 2][N - 2]) / yscal;
    *error_lower_order = fabs(T[N - 2][N - 2] - T[N - 3][N - 3]) / yscal;
}

// Factor for next stepsize, from the Hairer et al. book.
static double gbs_new_stepsize_fac(double scaled_err, int k) {
    return 0.94 * pow(0.65 / scaled_err, 1.0 / (2 * k - 1));
}

// Indices for the variables each task handles
struct gbs_extrapolation_indices {
    int *start_inds;
    int *var_counts;
};

void init_gbs_extrapolation_indices(struct gbs_extrapolation_indices *inds,
                                    int num_tasks, int num_vars) {
    inds->start_inds = calloc(1, sizeof(int[num_tasks]));
    inds->var_counts = calloc(1, sizeof(int[num_tasks]));
    // figure out which variables each task handles
    const int vars_per_task = num_vars / num_tasks;
    const int remainder = num_vars % num_tasks;

    for (int i = 0; i < num_tasks; ++i) {
        inds->start_inds[i] =
            i > 0 ? inds->start_inds[i - 1] + inds->var_counts[i - 1] : 0;
        inds->var_counts[i] = vars_per_task;
        // handle the remainder on the last (remainder) tasks, one extra var
        // each
        if (i >= num_tasks - remainder) { inds->var_counts[i] += 1; }
    }
}

void free_gbs_extrapolation_inds(struct gbs_extrapolation_indices *inds) {
    free(inds->start_inds);
    free(inds->var_counts);
}

// Extrapolate a subset of variables on each task, comm the results.
// Updates the fictitious timestep Hstep (data in/out) with the new optimal
// timestep estimate for the current k value,
// and new_k and Hstep_new_k (only out) for the values at the new proposed k.
// Whether the k is actually changed is decided in the main loop.
static bool gbs_extrapolate(struct gbs_storage *gs,
                            struct gbs_extrapolation_indices *ei, double *Hstep,
                            int *new_k, double *Hstep_new_k,
                            const struct ketju_system *R) {
    timer_start(&R->perf->time_gbsextr);

    const int rank = R->internal_state->task_info->main_group_rank;
    int start_ind = ei->start_inds[rank];
    int end_ind = start_ind + ei->var_counts[rank];

    // run the extrapolation
    const int k = R->internal_state->gbs_k;
    double y[k];
    double yscal;
    double error, error_lower_order;
    double max_error[2] = {0., 0.};
    double result;
    double(*data)[gs->num_vars] = (double(*)[gs->num_vars])gs->data;

    const double abs_tol = R->options->gbs_absolute_tolerance;

    for (int i = start_ind; i < end_ind; ++i) {
        for (int j = 0; j < k; ++j) {
            y[j] = data[j][i];
        }
        yscal = fabs(y[k - 1]) + abs_tol;

        if (i == 0) { yscal += R->options->output_time_absolute_tolerance; }

        extrapolate_single_variable(y, yscal, k, &result, &error,
                                    &error_lower_order);

        // Explicit nan/inf handling in case they are produced.
        // The code after this can then assume finite values.
        if (!isfinite(result)) {
            // Set the max error so that the step size gets cut by a factor of
            // ~10, and then bail as the step has failed anyway.
            max_error[0] =
                R->options->gbs_relative_tolerance * 0.65 / pow(0.1, 2 * k - 1);
            break;
        }

        data[k][i] = result;
        if (error > max_error[0]) { max_error[0] = error; }
        if (error_lower_order > max_error[1]) {
            max_error[1] = error_lower_order;
        }
    }
    timer_stop(&R->perf->time_gbsextr);

    timer_start(&R->perf->time_gbscomm);
    // check convergence across all task
    MPI_Allreduce(MPI_IN_PLACE, &max_error, 2, MPI_DOUBLE, MPI_MAX,
                  R->internal_state->task_info->main_comm);

    const double scaled_err = max_error[0] / R->options->gbs_relative_tolerance;
    const double scaled_err_lower =
        max_error[1] / R->options->gbs_relative_tolerance;
    const bool converged = max_error[0] < R->options->gbs_relative_tolerance;
    if (converged) {
        // success, communicate data and set new timestep
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, data[k],
                       ei->var_counts, ei->start_inds, MPI_DOUBLE,
                       R->internal_state->task_info->main_comm);

        // Should we change the integration order?
        // Check the conditions like in the Hairer et al. book case when
        // convergence occurs on line k. We can't do the other cases stopping at
        // k-1 or k+1 due to the parallelization strategy.
        // Odd k values cause slightly uneven workload balance for the
        // parallelization, but changing the order 2 at a time would likely
        // also lead to extra work in some cases, so it's not clear which is
        // better in general. Doing the standard change by one for now.
        double H_k = *Hstep * gbs_new_stepsize_fac(scaled_err, k);
        double H_km1 = *Hstep * gbs_new_stepsize_fac(scaled_err_lower, k - 1);
        int A_k = gbs_total_cost(k);
        int A_km1 = gbs_total_cost(k - 1);
        double W_k = A_k / H_k;
        double W_km1 = A_km1 / H_km1;
        double Hnext = H_k;

        // Do division instead of W_km1 * 0.9 to handle negative stepsizes
        // The k is kept in the range 4,kmax since very low orders are poor
        // for parallelization.
        if (W_k / W_km1 < 0.9 && k + 1 <= R->options->gbs_kmax) {
            *new_k = k + 1;
            Hnext = (H_k * gbs_total_cost(k + 1)) / A_k;
        } else if (W_km1 / W_k < 0.9 && k > 4) {
            *new_k = k - 1;
            Hnext = H_km1;
        } else {
            *new_k = k;
        }

        // Avoid increasing the stepsize too aggressively
        if (Hnext / *Hstep > 5.0) { Hnext = 5.0 * *Hstep; }
        if (H_k / *Hstep > 5.0) { H_k = 5.0 * *Hstep; }
        *Hstep_new_k = Hnext;
        *Hstep = H_k;
    } else {
        *new_k = k; // keep the same order
        double Hnextfac = gbs_new_stepsize_fac(scaled_err, k);
        // Ensure that the stepsize is decreased enough if the step failed,
        if (Hnextfac > 0.7) { Hnextfac = 0.7; }
        // but not excessively.
        if (Hnextfac < 1e-3) { Hnextfac = 1e-3; }
        *Hstep *= Hnextfac;
        *Hstep_new_k = *Hstep;
    }

    timer_stop(&R->perf->time_gbscomm);

    return converged;
}

static void update_system_from_gbs_result(struct ketju_system *R,
                                          struct gbs_storage *gs) {
    double *data = ((double(*)[gs->num_vars])gs->data)[gs->num_rows - 1];
    R->physical_state->time = data[0];
    int offset = 1;
    struct spanning_tree_coordinate_system *const mst = R->internal_state->mst;
    memcpy(mst->root_pos, &data[offset], sizeof(double[3]));
    offset += 3;
    memcpy(mst->root_vel, &data[offset], sizeof(double[3]));
    offset += 3;
    const size_t edge_vec_size = sizeof(double[R->num_particles - 1][3]);
    memcpy(mst->edge_rel_pos, &data[offset], edge_vec_size);
    offset += 3 * (R->num_particles - 1);
    memcpy(mst->edge_rel_vel, &data[offset], edge_vec_size);
    offset += 3 * (R->num_particles - 1);
    memcpy(R->physical_state->spin, &data[offset],
           sizeof(double[R->num_pn_particles][3]));

    into_Cartesian_from_MST(R);
}

///////////////////////////////////////
// Main integration loop and helpers //
///////////////////////////////////////

// Estimates a safe initial timestep that is sufficiently small compared to the
// smallest two-body timescale in the system.
static double estimate_initial_timestep(struct ketju_system *R) {
    struct ketju_system_internal_state *const is = R->internal_state;
    struct spanning_tree_coordinate_system *const mst = is->mst;
    struct ketju_system_physical_state *const ps = R->physical_state;

    int istart, jstart;
    int cum_num, block_size;
    loop_scheduling_N2(R->num_particles, is->task_info->main_group_size,
                       is->task_info->main_group_rank, &istart, &jstart,
                       &block_size, &cum_num);

    const int Nd = R->options->max_tree_distance;
    int d;
    int path[Nd];
    int sign[Nd];

    const double ss_soft = R->options->star_star_softening;
    double t2min = DBL_MAX;

    int c = 0;
    for (int i = istart; i < R->num_particles; i++) {

        double mi = ps->mass[i];
        int Li = mst->vertex[i].level;
        if (i > istart) { jstart = i + 1; }
        for (int j = jstart; j < R->num_particles; j++) {
            if (c >= block_size) { goto after_loops; }
            if (R->options->no_star_star_gravity && i >= R->num_pn_particles
                && j >= R->num_pn_particles) {
                continue;
            }

            double mj = ps->mass[j];
            int Lj = mst->vertex[j].level;
            double r2 = 0;
            double v2 = 0;
            bool proximity = false;

            // Use MST coordinates or global coordinates depending on
            // tree proximity. Like in separation_and_r2, but here velocity is
            // needed as well.
            if (abs(Li - Lj) < Nd) {
                proximity =
                    ketju_check_relative_proximity(j, i, Nd, R, &d, path, sign);
            }

            if (proximity) {
                for (int k = 0; k < 3; k++) {
                    double dr = 0;
                    double dv = 0;
                    for (int l = 0; l < d; l++) {
                        int index = path[l];
                        dr += sign[l] * mst->edge_rel_pos[index][k];
                        dv += sign[l] * mst->edge_rel_vel[index][k];
                    }
                    r2 += dr * dr;
                    v2 += dv * dv;
                }
            } else {
                for (int k = 0; k < 3; k++) {
                    double dr = ps->pos[j][k] - ps->pos[i][k];
                    r2 += dr * dr;
                    double dv = ps->vel[j][k] - ps->vel[i][k];
                    v2 += dv * dv;
                }
            }
            double soft = i < R->num_pn_particles || j < R->num_pn_particles
                              ? 0.
                              : ss_soft;
            const double invr = softened_inverse_r(sqrt(r2), soft);

            // Estimate the square of the approximate two-body time scale from a
            // combination of free propagation and keplerian timescales.
            double t2 = 1. / (invr * (invr * v2 + R->constants->G * (mi + mj)));
            t2min = fmin(t2min, t2);

            c++;
        }
    }
after_loops:

    MPI_Allreduce(MPI_IN_PLACE, &t2min, 1, MPI_DOUBLE, MPI_MIN,
                  is->task_info->main_comm);

    const double safetyfac = 0.05;
    return safetyfac * sqrt(t2min);
}

static void trivial_propagation(struct ketju_system *R, double time_interval) {
    // Only makes sense when there's just a single particle, but this guards
    // against the zero particle case (altohugh it shouldn't really happen)
    for (int i = 0; i < R->num_particles; ++i) {
        for (int k = 0; k < 3; ++k) {
            R->physical_state->pos[i][k] +=
                R->physical_state->vel[i][k] * time_interval;
        }
    }
    R->physical_state->time += time_interval;
}

static double step_total_work(const struct ketju_system *R) {
    double substeps = gbs_total_cost(R->internal_state->gbs_k);
    double interactions_per_step =
        .5 * R->num_particles
            * (R->options->no_star_star_gravity ? R->num_pn_particles
                                                : R->num_particles - 1)
        // PN terms are more expensive than other force terms, but
        // that's hard to quantify so just use the number of
        // interactions here as well.
        + ((double)R->num_pn_particles * (R->num_pn_particles - 1)
           * (R->options->PN_flags & KETJU_PN_THREEBODY
                  ? R->num_pn_particles - 2
                  : 1));
    return substeps * interactions_per_step;
}

void ketju_run_integrator(struct ketju_system *R, double time_interval) {

    // special case for when no actual integration is needed
    if (R->num_particles < 2) {
        timer_start(&R->perf->time_total);
        trivial_propagation(R, time_interval);
        timer_stop(&R->perf->time_total);
        return;
    }

    timer_start(&R->perf->time_total);

    // Init GBS stuff
    R->internal_state->gbs_k =
        R->options->gbs_kmax >= 4 ? 4 : R->options->gbs_kmax;
    ketju_init_gbs_comms(R);

    int *gbs_rows_for_groups =
        malloc(sizeof(int[R->internal_state->group_info->num_gbs_groups]
                         [R->internal_state->gbs_k]));
    get_gbs_rows_for_groups(R, gbs_rows_for_groups);
    int *gbs_rows_for_this_group =
        &gbs_rows_for_groups[R->internal_state->task_info->gbs_group_index
                             * R->internal_state->gbs_k];

    struct gbs_storage gbs_storage;
    allocate_gbs_storage(&gbs_storage, R->num_particles, R->num_pn_particles,
                         R->internal_state->gbs_k);
    struct gbs_extrapolation_indices gbs_indices;
    init_gbs_extrapolation_indices(
        &gbs_indices, R->internal_state->task_info->main_group_size,
        gbs_storage.num_vars);

    // Other init stuff

    ketju_construct_MST(R);
    int steps_since_reconstruction = 0;

    compute_total_energy(R);

    // Let's not try to step over whole interval in one go for safety,
    // even if the timescale estimate would allow it.
    double dt = copysign(1, time_interval)
                * fmin(estimate_initial_timestep(R), fabs(time_interval) / 2);
    double Hstep = dt * R->internal_state->U;
    // set the start point to 0 for better accuracy when the t_start is large
    const double t_start = R->physical_state->time;
    R->physical_state->time = 0;
    const double initial_B = R->internal_state->B;

    struct state_backup old_state;
    allocate_state_backup(&old_state, R->num_particles, R->num_pn_particles);
    store_state_backup(&old_state, R);

    bool finished = false;
    do {
        // Check mergers between BHs (PN particles)
        // Need to be propagating forwards in time for mergers to make sense.
        if (R->options->enable_bh_mergers && Hstep > 0) {
            int num_merged = ketju_do_bh_mergers(R);
            if (num_merged) {
                // We had mergers, need to redo tree and other things
                if (R->num_particles < 2) {
                    // Particle motion becomes trivial.
                    // Need to only propagate the CoM
                    trivial_propagation(R, time_interval
                                               - R->physical_state->time);
                    break;
                }
                ketju_construct_MST(R);
                steps_since_reconstruction = 0;
                compute_total_energy(R);
                store_state_backup(&old_state, R);
                gbs_storage.num_vars -= num_merged * 9; // pos, vel, spin
                free_gbs_extrapolation_inds(&gbs_indices);
                init_gbs_extrapolation_indices(
                    &gbs_indices, R->internal_state->task_info->main_group_size,
                    gbs_storage.num_vars);
            }
        }
        // Run the gbs substeps, communicating results in the background.
        const int k = R->internal_state->gbs_k;
        int request_index = 0;
        MPI_Request requests[k];
        post_gbs_receives(requests, &request_index, &gbs_storage,
                          gbs_rows_for_groups, R);
        for (int i = 0; i < k && gbs_rows_for_this_group[i]; ++i) {
            if (i) {
                // not needed on the first round
                restore_state_backup(&old_state, R);
            }
            int row = gbs_rows_for_this_group[i];
            mst_leapfrog(R, Hstep, num_gbs_substeps(row));
            store_gbs_data(&gbs_storage, R, row);
            post_gbs_send(requests, &request_index, &gbs_storage, R, row);
        }
        finish_gbs_comm(requests, request_index, R);

        int new_k;
        double Hstep_new_k;
        bool converged = gbs_extrapolate(&gbs_storage, &gbs_indices, &Hstep,
                                         &new_k, &Hstep_new_k, R);

        R->perf->total_work += step_total_work(R);
        if (converged) {
            update_system_from_gbs_result(R, &gbs_storage);
            compute_total_energy(R);
            R->perf->successful_steps += 1;
            steps_since_reconstruction += 1;

            double remaining_time = time_interval - R->physical_state->time;
            if (fabs(remaining_time / time_interval)
                    < R->options->output_time_relative_tolerance
                || fabs(remaining_time)
                       < R->options->output_time_absolute_tolerance) {
                finished = true;
            } else {
                bool change_order = false;
                // Check that we're not overstepping or already past the end
                // and alter timestep if required.
                // Should work also for negative time intervals
                if (fabs(remaining_time) < fabs(Hstep / R->internal_state->U)
                    || time_interval * remaining_time < 0) {

                    Hstep = R->internal_state->U * remaining_time;

                    if (new_k < k) {
                        // Order can be reduced even if near the end, since it
                        // should always reduce the workload.
                        // But the new proposed timestep can be smaller than the
                        // step to the end, and have the wrong sign.
                        change_order = true;
                        if (fabs(Hstep_new_k) < fabs(Hstep)) {
                            Hstep = copysign(Hstep_new_k, Hstep);
                        }
                    }
                } else if (new_k != k
                           && fabs(remaining_time)
                                  > fabs(Hstep_new_k / R->internal_state->U)) {
                    // Change order only if it doesn't result in overstepping,
                    // otherwise the work/step calculation doesn't give the
                    // correct result.
                    change_order = true;
                    Hstep = Hstep_new_k;
                }

                if (change_order) {
                    R->internal_state->gbs_k = new_k;
                    // Reinit the gbs structures
                    ketju_free_gbs_comms(R);
                    ketju_init_gbs_comms(R);
                    gbs_rows_for_groups = realloc(
                        gbs_rows_for_groups,
                        sizeof(
                            int[R->internal_state->group_info->num_gbs_groups]
                               [R->internal_state->gbs_k]));
                    get_gbs_rows_for_groups(R, gbs_rows_for_groups);
                    gbs_rows_for_this_group =
                        &gbs_rows_for_groups[R->internal_state->task_info
                                                 ->gbs_group_index
                                             * R->internal_state->gbs_k];

                    free_gbs_storage(&gbs_storage);
                    allocate_gbs_storage(&gbs_storage, R->num_particles,
                                         R->num_pn_particles,
                                         R->internal_state->gbs_k);
                    // extrapolation indices not changed here
                }

                if (steps_since_reconstruction
                    >= R->options->steps_between_mst_reconstruction) {
                    ketju_construct_MST(R);
                    steps_since_reconstruction = 0;
                }

                store_state_backup(&old_state, R);
            }
        } else {
            restore_state_backup(&old_state, R);
            // stepsize is updated by the extrapolation function

            R->perf->failed_steps += 1;
        }

        if (R->options->max_step_count > 0
            && R->perf->failed_steps + R->perf->successful_steps
                   > R->options->max_step_count) {
            if (R->internal_state->task_info->main_group_rank == 0) {
                ketju_debug_dump_integrator_state(
                    R, "ketju_integrator_max_steps_state_dump.bin");
            }
            terminate("Exceeded options->max_step_count steps! State dumped.");
        }

    } while (!finished);

    // Restore the initial physical time offset
    R->physical_state->time += t_start;

    // free temp storage
    free_state_backup(&old_state);
    free_gbs_storage(&gbs_storage);
    free_gbs_extrapolation_inds(&gbs_indices);
    free(gbs_rows_for_groups);
    ketju_free_gbs_comms(R);

    // Note that this is not really the error when dissipative PN forces are
    // enabled
    R->perf->relative_energy_error =
        (R->internal_state->B - initial_B) / initial_B;

    timer_stop(&R->perf->time_total);
}
