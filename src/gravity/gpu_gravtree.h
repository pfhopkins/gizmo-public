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
/* the packet shape the last primary walk used on this rank (team size, members per packet), 0/0 when
   every candidate took the single-target walk; written into the per-call timings record */
void gpu_gravtree_packet_shape(int *team, int *q_dev);

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
