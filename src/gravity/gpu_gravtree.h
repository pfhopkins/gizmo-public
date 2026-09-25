/* gpu_gravtree.h
 *
 * Core gravity tree walk on GPU (no optional payloads).
 *
 * Runs BEFORE the OpenMP-parallelized CPU primary loop (gravity_primary_loop)
 * as an opportunistic accelerator. For each active particle:
 *   - GPU thread walks the local tree using the SoA mirror.
 *   - If it encounters a pseudo-particle (remote node owned by another rank),
 *     it sets failed[i]=1 and exits. The host then leaves ProcessedFlag[i]
 *     unset so the CPU primary loop handles it (with the existing MPI export
 *     machinery, unchanged).
 *   - If it completes successfully, it writes P[i].GravAccel and sets
 *     ProcessedFlag[i]=1 so the CPU loop skips it.
 *
 * GPU gravity tree (always active on Kokkos builds). PMGRID, ADAPTIVE_GRAVSOFT_*,
 * EVALPOTENTIAL, RT/SINK/SINGLE_STAR/CR/TIDAL/JERK payloads, periodic Ewald,
 * HERMITE/ATFU, and FIRE_BHS MencInRcrit are all supported via the
 * LET work.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifndef GIZMO_GPU_GRAVTREE_H
#define GIZMO_GPU_GRAVTREE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Run the GPU pre-pass over all active particles.
 *
 * Called from gravity_tree() before the primary_loop OpenMP fan-out. Updates
 * P[i].GravAccel and ProcessedFlag[i] for successful particles. Particles
 * that hit a pseudo-particle (remote node) are left untouched, so the
 * existing CPU primary loop handles them as usual.
 *
 * Returns the number of particles that completed successfully on GPU. When
 * the build lacks Kokkos, this is an unconditional no-op
 * (returns 0) so callers in gravity_tree() stay simple.
 * host_candidates_left (may be NULL) receives how many walk candidates this
 * pass leaves to the host loop -- an upper bound on any early return. */
int gpu_gravtree_walk_primary(int *host_candidates_left);
/* The device schedule the last primary walk used on this rank, written into the per-call timings
 * record. It carries the row the call ASKED for and the row it actually ran, because the two
 * differ exactly when the backend could not launch the requested one -- and a pricing arm that
 * cannot see that difference is reading a measurement of a shape nobody chose. Mode FLAT with
 * zeros is a call that took the ordinary one-lane-per-target walk.
 */
/* Which device schedule ran. 0 = one lane per target, the target's own serial walk (the ordinary
 * schedule, and what a call takes whenever the device has a target for every lane); 1 = one walker
 * per team with the members sharing its traversal; 2 = one target per team with the team's lanes
 * sharing that target's traversal. -1 when no device walk ran at all.
 *
 * None of these is a routing decision: every one of them is on the device, and which of them ran
 * says nothing about whether the call could have gone to the host. */
#define GRAV_PACKET_MODE_NONE        (-1)
#define GRAV_PACKET_MODE_FLAT         0
#define GRAV_PACKET_MODE_PACKET       1
#define GRAV_PACKET_MODE_COOPERATIVE  2

struct gpu_grav_packet_shape_t {
    int mode;              /* one of the above */
    int team;              /* threads in the team */
    int q_dev;             /* members (targets) per packet */
    int n_walkers;         /* threads that traverse; 1 is the depth-first single-walker traversal */
    int frontier;          /* items the shared frontier holds */
    int chunk;             /* records held between flushes */
    int steps_per_round;   /* node steps a walker takes before the round boundary */
    int row_requested;
    int row_effective;
    long long scratch_bytes;
};
void gpu_gravtree_packet_shape(struct gpu_grav_packet_shape_t *out);

/* Packets the device engine gave up on the last primary walk, by reason, so that a traversal
 * exhausting the engine's continuation budget is visible instead of being a silent slow path.
 * Slot 0 is unused; the rest follow the engine's own reason order (malformed index, stale
 * source, pseudo-particle, no continuation, unusable record, no progress in a round). Zero on a
 * host-routed call. */
#define GRAV_PACKET_FAIL_REASON_SLOTS 7
void gpu_gravtree_packet_failures(long long *out, int n);
int  gpu_gravtree_packet_failure_reasons(void);
/* Run totals on this rank: device gravity calls whose sources were brought current by the
 * discovery subset, and calls that fell back to drifting everything. */
void gpu_gravtree_subset_drift_counts(long long *taken, long long *declined);

/* GPU Ewald-correction walk. Called from gravity_tree() when Ewald_iter==1
 * (pure-tree periodic, BOX_PERIODIC && !GRAVITY_NOT_PERIODIC && !PMGRID).
 * Mirrors force_treeevaluate_ewald_correction mode=0: walks the local tree a
 * second time, accumulates the periodic-image correction via trilinear
 * interpolation of the fcorrx/y/z look-up tables, and adds the result to
 * P[i].GravAccel.
 *
 * Only the local tree walk runs on GPU. Pseudo-particle hits leave
 * ProcessedFlag unset, so the CPU Ewald secondary loop finishes those via
 * MPI export.  No-op when the build lacks Kokkos. Returns
 * number of successfully walked targets. */
int gpu_ewald_walk_primary(void);

#ifdef __cplusplus
}
#endif

#endif /* GIZMO_GPU_GRAVTREE_H */
