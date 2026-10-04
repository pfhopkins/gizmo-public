/* declarations/gpu_recorder_claim.h -- how the particle touched set is claimed from a kernel.
 *
 * The touched set records the local particles a walk reached, so only those are drifted: a
 * stamp array as long as the local particles, a compacted list, a cursor, and a generation that
 * invalidates every stamp in O(1).  It has two claimers -- a host phase and a device kernel --
 * and exactly ONE claim belongs to it, or the two drift apart.
 *
 * The claim lives here rather than beside its recorder because its body needs Kokkos atomics,
 * and mesh/neighbor_list.h is plain data read by host units no device compiler ever sees.
 * (The gravity tree's node dirty set follows the same idiom but is claimed and answered on the
 * host only, so its claim lives with its storage in gravity/gpu_gravity_tree.cc.)
 *
 * Include after Kokkos and allvars.h; every includer already has both.
 */
#ifndef GPU_RECORDER_CLAIM_H
#define GPU_RECORDER_CLAIM_H

#include "../mesh/neighbor_list.h"   /* struct GxTouchedSet, gx_touched_owner_t, the GX_WALK_ANOMALY_* codes -- plain data, no Kokkos */


/* ============================================================================
 * THE PARTICLE TOUCHED SET
 * ========================================================================== */

/* THE claim, and the only one -- the fused walk's recording visitor comes through here.
 *
 * `owner` is the phase the caller believes it is in, checked against the owner the host wrote
 * when it opened the epoch and carried by value in the view.  It is what stops one caller
 * appending into another's live generation: the touched set deliberately keeps ONE generation
 * across several passes of a fused call while resetting the cursor at each consume, so a second
 * claimer arriving mid-call would be indistinguishable from that call's own later pass.
 *
 * Overflow is structurally impossible rather than handled -- the stamp admits each owned slot
 * once per generation and the list is as long as there are owned slots -- so the full case is
 * reported, not recovered from.  There is no `unsafe` word in device memory for this recorder
 * for the same reason; the `anomaly` word both callers already carry is the whole channel. */
KOKKOS_INLINE_FUNCTION void
gx_touched_set_claim_in(const struct GxTouchedSet &ts, int j, int owner, int *anomaly)
{
    if(!ts.seen || !ts.list || !ts.counter) {
        if(anomaly) {Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_TOUCHED_SET_FULL);}
        return;
    }
    if(ts.owner != owner) {
        if(anomaly) {Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_RECORDER_OUT_OF_PHASE);}
        return;
    }
    if(j < 0 || j >= ts.capacity) {
        if(anomaly) {Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_TOUCHED_SET_FULL);}
        return;
    }
    /* Claim the slot for this generation.  The exchange returns what was there, so exactly one
     * work item sees a value other than the current generation and exactly one item appends --
     * several actives reaching the same particle is the ordinary case, not the exception. */
    if(Kokkos::atomic_exchange(&ts.seen[j], ts.gen) == ts.gen) {return;}
    const int slot = Kokkos::atomic_fetch_add(ts.counter, 1);
    if(slot < ts.capacity) {ts.list[slot] = j;}
    else if(anomaly)       {Kokkos::atomic_store(anomaly, GX_WALK_ANOMALY_TOUCHED_SET_FULL);}
}

#endif /* GPU_RECORDER_CLAIM_H */
