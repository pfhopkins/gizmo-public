/* gravity/gpu_node_dirty_claim.h -- the node dirty set's claim, and only that.
 *
 * The recorder itself (its storage, its epochs, its consumers) lives in gpu_gravity_tree.cc.
 * What has to be visible in a SECOND translation unit is just enough to claim a node from a
 * kernel: the shared control block, the by-value view of the recorder, and the claim.
 *
 * This is a private header with exactly two includers -- the recorder's own implementation and
 * the gravity discovery pre-walk. It is deliberately NOT part of gpu_gravity_tree.h: that is a
 * declarations header pulled in by sixteen translation units, several of them host-only, and a
 * Kokkos-atomic body there would drag Kokkos into every one of them. Anything that needs the
 * recorder's API rather than its claim keeps using gpu_gravity_tree.h.
 *
 * Requires Kokkos and allvars.h already in scope; both includers have them.
 */
#ifndef GPU_NODE_DIRTY_CLAIM_H
#define GPU_NODE_DIRTY_CLAIM_H

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
 * consumer before it reads the list.  Strictly stronger than the ordering it replaces. */
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

#endif /* GPU_NODE_DIRTY_CLAIM_H */
