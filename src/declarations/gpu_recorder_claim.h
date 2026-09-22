/* declarations/gpu_recorder_claim.h -- how a generation-stamped recorder is claimed from a kernel.
 *
 * Two recorders in this code follow the same idiom: a stamp array as long as the things that
 * can be recorded, a compacted list, a cursor, and a generation that invalidates every stamp in
 * O(1).  The node dirty set records nodes whose mirror needs attention before the next device
 * gravity walk; the touched set records the local particles a walk reached, so only those are
 * drifted.  Each has two claimers -- a host phase and a device kernel -- and exactly ONE claim
 * belongs to each, or the two claimers drift apart.
 *
 * Both claims live here rather than beside their recorders because the bodies need Kokkos
 * atomics, and neither recorder's own public header can carry Kokkos: gpu_gravity_tree.h is a
 * declarations header pulled in by sixteen translation units, several of them host-only, and
 * mesh/neighbor_list.h is plain data read by host units no device compiler ever sees.  One
 * header for the claims keeps the shared idiom visible and stops a third copy growing.
 *
 * Include after Kokkos and allvars.h; every includer already has both.
 */
#ifndef GPU_RECORDER_CLAIM_H
#define GPU_RECORDER_CLAIM_H

#include "../mesh/neighbor_list.h"   /* struct GxTouchedSet, gx_touched_owner_t -- plain data, no Kokkos */


/* What a recording walk reports through `anomaly`.  Distinct values because the states are
 * distinct: one says the tree cannot be walked, one says a caller's own bookkeeping broke, one
 * says a claim arrived outside its owner's phase.  Zero means nothing was reported; callers test
 * against it and must not assume 1.  How fatal each is belongs to the caller, not to the code:
 * the fused walk stops the run on any of them, while the gravity pre-walk answers an out-of-phase
 * claim by declining to the full drift it already has. */
#define GX_WALK_ANOMALY_MALFORMED_TREE       1  /* index in no class, or an unfilled view */
#define GX_WALK_ANOMALY_TOUCHED_SET_FULL     2  /* touched-set list shorter than the set it recorded */
#define GX_WALK_ANOMALY_RECORDER_OUT_OF_PHASE 3 /* a claim in an epoch its owner does not hold */


/* ============================================================================
 * THE NODE DIRTY SET
 * ========================================================================== */

struct gpu_node_dirty_ctl_t {
    unsigned int generation;     /* stamps equal to this are claimed in the current epoch */
    int          count;          /* append cursor into list[] */
    int          unsafe;         /* sticky: out-of-range, overflow, or a claim out of phase */
    int          owner;          /* which phase may claim right now (see gpu_node_dirty_owner_t) */
    long long    unsafe_events;  /* how often the fail-safe fired, run-total */
};

/* Everything a claimer needs, and nothing it does not: copied BY VALUE into a kernel, so
 * the claim never dereferences a host pointer from device code. */
struct gpu_node_dirty_view_t {
    unsigned int                *seen;   /* [cap] generation stamps, never cleared */
    int                         *list;   /* [cap] compacted node indices */
    struct gpu_node_dirty_ctl_t *ctl;
    int                          cap;
    int                          base;   /* All.TreeNodeIndexBase at the epoch's start */
};


/* Ordering.  The claim publishes geometry the claimer wrote just before it, and the
 * consumer must observe that geometry once it observes the claim.  The shipped host-only
 * version rode release/acquire on the stamp itself; one claim now serves host and device,
 * where a per-object memory order is not portably expressible, so the pairing is carried
 * by explicit fences instead -- a full fence before the stamp exchange, and one in the
 * consumer before it reads the list.  This is sufficient under the epoch's phase
 * separation -- every producer phase complete before the consuming one begins, host
 * writers joined, the device kernel fenced -- and it is that separation, not the fences
 * alone, that the cross-boundary runtime gate validates. */
/* Device-callable: the claim below runs in a kernel, and a host-only helper reached from a
 * KOKKOS_INLINE_FUNCTION is a device-annotation defect no host compiler can see. */
KOKKOS_INLINE_FUNCTION void nd_mark_unsafe_ctl_(struct gpu_node_dirty_ctl_t *ctl)
{
    if(!ctl) {return;}
    Kokkos::atomic_store(&ctl->unsafe, 1);
    Kokkos::atomic_fetch_add(&ctl->unsafe_events, 1LL);
}

/* THE claim, and the only one.  `owner` is the phase the caller believes it is in; a
 * mismatch is a stopped invariant rather than a tolerated race, because the whole point of
 * naming an epoch owner is that host and device claims never interleave. */
KOKKOS_INLINE_FUNCTION void
gpu_node_dirty_claim_in(const struct gpu_node_dirty_view_t &v, int no, int owner)
{
    if(!v.seen || !v.list || !v.ctl) {nd_mark_unsafe_ctl_(v.ctl); return;}
    if(v.ctl->owner != owner)       {nd_mark_unsafe_ctl_(v.ctl); return;}
    const int k = no - v.base;
    if(k < 0 || k >= v.cap)         {nd_mark_unsafe_ctl_(v.ctl); return;}
    Kokkos::memory_fence();   /* the geometry this claim refers to is published before the claim is */
    const unsigned int gen  = v.ctl->generation;   /* constant within an epoch */
    const unsigned int prev = Kokkos::atomic_exchange(&v.seen[k], gen);
    if(prev == gen) {return;}                      /* one claim per node per epoch */
    const int slot = Kokkos::atomic_fetch_add(&v.ctl->count, 1);
    if(slot < v.cap) {v.list[slot] = no;}
    else             {nd_mark_unsafe_ctl_(v.ctl);}
}

/* Fill from the recorder's own storage; defined beside that storage, in gpu_gravity_tree.cc. */
struct gpu_node_dirty_view_t gpu_node_dirty_view(void);


/* ============================================================================
 * THE PARTICLE TOUCHED SET
 * ========================================================================== */

/* THE claim, and the only one -- the fused walk's recording visitor and the gravity discovery
 * pre-walk both come through here.
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
