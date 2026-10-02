/* sfc_tiles.h — SFC-ordered tile spatial index for neighbor finding and ghost exchange.
 *
 * Tiles are runs of TILE_TARGET_SIZE members in Morton order (built in gpu_neighbor_list.cc).
 * Each tile stores a bounding box and max kernel radius (hmax).
 * This provides finer spatial granularity than the top-level tree leaves
 * (~100-1000 leaves) while being cheaper than per-particle checks.
 *
 * Used for:
 *   1. Neighbor finding: tile overlap check replaces uniform grid cell-list
 *   2. Ghost exchange: tile-level overlap criterion replaces leaf-level
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef SFC_TILES_H
#define SFC_TILES_H

#include "neighbor_list.h"
#include "nlr_radius_policy.h"  /* mode_b_radius_policy_t + MODE_B_RADIUS_* + SSOT wrappers */

#define TILE_TARGET_SIZE 64  /* particles per tile (tunable) */
#define TILE_BVH_STACK_SIZE 64  /* traversal stack depth (log2(ntiles) + margin) */

/* Per-type hmax: stored as fixed-size vector indexed by P[].Type ∈ [0,5].
 * Dimension matches the GHOST_TYPE_* bitmask (bit k ↔ Type k). hmax_by_type[t]
 * is the max KernelRadius across particles of Type t in the subtree; 0 if
 * none. Opener computes hmax_eff = max over t in supply_mask of hmax_by_type[t]
 * — this avoids the scalar-hmax contamination where a single DM particle
 * with KernelRadius=582 poisons every search whose supply mask doesn't even
 * include DM. (See ghost_exchange_roadmap_2026-05-05.md "per-type hmax".) */
#define TILE_NUM_PTYPES 6

#define SIDX_ROW_WIDTH 5   /* x, y, z at the reference time; reach now; reach once drifted */

/* Axis periodicity flags for tile-based neighbor search and ghost exchange */
#if defined(BOX_PERIODIC) && !defined(BOX_REFLECT_X) && !defined(BOX_OUTFLOW_X)
#define TILE_PERIODIC_X 1
#else
#define TILE_PERIODIC_X 0
#endif
#if defined(BOX_PERIODIC) && !defined(BOX_REFLECT_Y) && !defined(BOX_OUTFLOW_Y)
#define TILE_PERIODIC_Y 1
#else
#define TILE_PERIODIC_Y 0
#endif
#if defined(BOX_PERIODIC) && !defined(BOX_REFLECT_Z) && !defined(BOX_OUTFLOW_Z)
#define TILE_PERIODIC_Z 1
#else
#define TILE_PERIODIC_Z 0
#endif

/* The index describes its members as of one REFERENCE TIME, the time it was built.  Each member has a row
 * of SIDX_ROW_WIDTH doubles: where the drift puts it at that time (particle_motion_envelope, x,y,z,
 * unwrapped -- the position the pair test reads; the tile and BVH boxes hold its primary-box image), its
 * reach at its own clock, and the most that reach can be once it is drifted, from whatever clock
 * (nlr_particle_symmetric_radius_after_drift, never below the reach itself).  A tile or BVH node bounds its
 * members by a box around those positions, the range of velocities they can advance at from there, a
 * residual speed for motion that is not a straight line, and the largest drifted reach.  Positions are never
 * rewritten while the index is kept: at a later time, an undilated drift interval D on, the box is shifted
 * by the velocity range times D and widened by the residual speed times D, which holds every member
 * wherever the drift has put it or will.  A kick or a direct write of a member's velocity RAISES the
 * velocity range (sfc_member_motion_range); the reach bands are raised when a reach grows.  Both only ever
 * widen, and the walk never writes them. */
struct sfc_tile_t {
    int first;                        /* first pool slot of this tile (tile t covers slots t*TILE_TARGET_SIZE on) */
    int count;                        /* number of particles in this tile */
    int bvh_leaf;                     /* the BVH node that holds this tile */
    double lo[3];                     /* members' positions at the reference time, each widened by its half-width */
    double hi[3];
    double hmax;                      /* max drifted reach in tile (any type) */
    double hmax_by_type[TILE_NUM_PTYPES]; /* max drifted reach PER TYPE */
    double u_min[3], u_max[3];        /* range of the members' velocities per unit undilated drift interval */
    double rho;                       /* largest residual speed among the members */
    double hw;                        /* largest member half-width at the reference time */
};

/* BVH node over the tiles: a midpoint split of the tile range, bounds the union of the children.
 * Enables O(log ntiles) spatial pruning for neighbor search, critical for zoom-in sims with
 * h/box ~ 10^-6.  Nodes are numbered children first, so every node's index is above its
 * children's and the root is the last node. */
struct tile_bvh_node_t {
    double lo[3], hi[3];                  /* bounding box of subtree at the reference time */
    double hmax;                          /* max drifted reach in subtree (any type) */
    double hmax_by_type[TILE_NUM_PTYPES]; /* max drifted reach PER TYPE in subtree */
    double u_min[3], u_max[3];            /* union of the subtree's velocity ranges */
    double rho;                           /* max over the subtree's residual speeds */
    int left, right;                      /* children: >= 0 = internal node index, < 0 = -(tile_index+1) for leaf */
    int parent;                           /* -1 at the root */
};

/* The range of velocities member j can advance at from now on, per unit undilated drift interval, and its
 * residual speed: its transport velocity (particle_transport_velocity), both ends of the range.  Along an
 * axis with a reflecting or outflow side, where the drift can turn the particle round, that axis takes
 * instead the symmetric range of its speed bound, as the tree does; every other axis keeps its own.  The
 * build folds this and a raise reads it, so the two cannot differ.  Returns 0, or 1 when the velocity or
 * speed is not finite: nothing can bound that member. */
KOKKOS_INLINE_FUNCTION
int sfc_member_motion_range(int j, const struct particle_data *P, const struct gas_cell_data *cells,
                            double u_lo[3], double u_hi[3], double *rho)
{
    double u[3];
    particle_transport_velocity(j, P, cells, u, rho);
    for(int k = 0; k < 3; k++) {u_lo[k] = u[k]; u_hi[k] = u[k];}
#if BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE
    {   /* a side the drift acts on: reflect or outflow, lower (code 0 or -1) or upper (0 or 1) */
        const double s = particle_motion_speed_bound(j, P, cells);
        for(int k = 0; k < NUMDIMS; k++) {
            const int rf = special_boundary_condition_xyz_def_reflect[k], of = special_boundary_condition_xyz_def_outflow[k];
            if(rf == 0 || rf == -1 || rf == 1 || of == 0 || of == -1 || of == 1) {u_lo[k] = -s; u_hi[k] = s;}
        }
    }
#endif
    for(int k = 0; k < 3; k++) {if(!(fabs(u_lo[k]) < 1.0e30 && fabs(u_hi[k]) < 1.0e30)) {return 1;}}
    return (*rho >= 0.0 && *rho < 1.0e30) ? 0 : 1;
}

/* How a walk reads the index at the time of the search.  `D` is the undilated drift interval from the
 * index's reference time to now.  `reach_current` says every member is current, so a member's reach is its
 * stored one; otherwise it is its drifted reach.  `exact` says the rows ARE the members' current positions
 * and reaches -- no member has moved or changed since they were written -- so the leaf applies the exact
 * test.  A node opens on its drifted reach either way: over-opening a node costs a few leaf tests, never a
 * neighbour. */
struct sfc_walk_frame {
    double D;
    int reach_current;
    int exact;
};

/* Index membership: the particles of the types in type_bitmask with positive mass.  Membership has a
 * single definition, shared by the spatial index build (gpu_neighbor_list.cc) and the ghost supply pool;
 * it reads no position and no radius, so a pool built once stays valid while particles move. */
KOKKOS_INLINE_FUNCTION
int sfc_pool_member(const struct particle_data *p, int type_bitmask)
{
    if(!((1 << p->Type) & type_bitmask)) return 0;
    if(p->Mass <= 0) return 0;
    return 1;
}

/* The supply pool: the members of P[0..num_total), in P[] order.  Returns num_pool and, when
 * pool_indices_out is non-NULL, a mymalloc'd index array the caller owns. */
int build_sfc_supply_pool(struct particle_data *P, int num_total,
                          int type_bitmask, int **pool_indices_out);

#endif /* SFC_TILES_H */
