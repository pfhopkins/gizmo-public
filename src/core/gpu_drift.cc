/* gpu_drift.cc
 *
 * Batched device drift. The bulk drift sites hand an index list to
 * drift_particles_batch, which compacts it to the particles that actually need
 * advancing and, when there are enough of them, runs the drift body on the device
 * over compact staged copies instead of walking P[] and CellP[] on the host.
 *
 * The body is the same one the host runs: drift_particle_impl in
 * core/drift_particle_functions.h. There is no second copy of the physics, and no
 * device-only approximation of any part of it.
 *
 * Routing is on the number of particles needing a drift, against the same
 * GPU_MIN_PARTICLES_FOR_OFFLOAD the other batched device loops use. Below it the
 * caller's work is done by the ordinary host loop, so a step with few active
 * particles pays nothing: no staging buffers, no device memory, no copies.
 *
 * A drift that does not happen is a wrong position, not a missing improvement, so
 * every path that cannot reach the device -- no staging memory, no table mirror --
 * falls back to the host loop and still drifts every particle it was given.
 *
 * There are two device routes.  In-place runs the body over the canonical arrays
 * themselves and is used wherever those arrays are device-visible; the staged route
 * copies compact batches across and answers where they are not.  Both run the same
 * body on the same particles, and the choice is a property of the allocation model,
 * not a tuning knob.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../core/timestep_functions.h"
#include "../declarations/gpu_dispatch_templates.h"
#include "../system/gpu_particles_arena.h"
#include "drift_particle_functions.h"

/* Device-visible copy of the drift and gravkick tables for this translation unit.
   Per-TU because a device symbol cannot be shared across TUs without relocatable
   device code; the allocate-and-fill policy is shared, in the arena. */
static double *drift_kick_table_dev_ = NULL;

/* Batch size, the same cap the cooling loop uses. What it buys is a bounded staging
   footprint: the helper holds this many compact structs on the host and again on the
   device while a call runs. Device stack depth is not the binding reason here -- the
   drift's equation-of-state call reads a cached composition rather than running
   cooling's iterative solve, so this chain is the shallower of the two, and a cap that
   is safe for cooling is safe for it. */
static const int GPU_DRIFT_BATCH_SIZE = 32768;

/* Drift idx[0..n_idx) to time1 on the host, the way the callers did before this
   entry existed. Threaded because all four bulk sites already threaded their loops
   and routing through here must not quietly serialise them; the body acts on one
   particle with no pair coupling and no random draws, so the indices being distinct
   is the whole of the argument. */
static void drift_particles_host_(const int *idx, int n_idx, integertime time1)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) if(n_idx >= 16)
#endif
    for(int k = 0; k < n_idx; k++) {drift_particle(idx[k], time1);}
}

/* The status every exit below reports. A controlled stop is first-set-wins and is
   never cleared, so a request raised anywhere in this call -- a trapped or failed
   kernel consumed by gizmo_gpu_check_last_error, a staging buffer or table mirror
   that could not be served -- is still standing here, and so is one that was
   already pending when the call began. Both mean the same thing to a caller: the
   run is draining and nothing about this particle set may be published.
   Deliberately conservative in the one harmless direction -- an unrelated pending
   stop reports failure, which costs a fallback on a run that is already ending. */
static int drift_batch_status_(void)
{
    return (gizmo_controlled_stop_local_reason() != NULL) ? 1 : 0;
}

/* A null idx means the contiguous range [0, n_idx), which is what the full-drift
   site hands over: materialising an identity array there would be an extra
   allocation and an extra pass to say nothing. */
int drift_particles_batch(const int *idx, int n_idx, integertime time1,
                          int *out_drifted, int *out_n_drifted)
{
    if(out_n_drifted) {*out_n_drifted = 0;}
    if(n_idx <= 0) {return drift_batch_status_();}

    /* Compact to the particles that are not already at time1. The drift body returns
       immediately for the rest, so staging them would be copying a whole struct each
       way to do nothing. At the full-drift site this pass is O(NumPart) at a site that
       is already O(NumPart) by definition; the other sites filter a list they hold. */
    /* The one allocation is made before the scan and caught: a caller may be inside a window with no
       collectives, where an allocation failure has to become a controlled stop taken by every rank
       rather than an exception that ends this one.  Each index writes its own slot -- the particle, or
       -1 when it is already current -- so nothing grows inside the parallel region, whatever team size
       OpenMP provides; a stable compaction then keeps the list in index order at any thread count. */
    std::vector<int> needs_drift;
    try {needs_drift.resize((size_t) n_idx);}
    catch(const std::bad_alloc &) {
        printf("drift_particles_batch: task %d could not reserve the work list for %d particles\n", ThisTask, n_idx);
        fflush(stdout);
        gizmo_request_controlled_stop(7739, "drift_particles_batch: could not allocate its work list",
                                      __FILE__, __LINE__, __FUNCTION__);
        return drift_batch_status_();
    }
    int *slot = needs_drift.data();
    /* Threaded above a handful of indices: the full-drift site hands over every particle, and reading
       Ti_current out of the particle array is a strided pass over the whole of it, which is a large part
       of what running the drift in bulk is meant to remove. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n_idx >= 16)
#endif
    for(int k = 0; k < n_idx; k++)
    {
        const int i = idx ? idx[k] : k;
        slot[k] = (P[i].Ti_current != time1) ? i : -1;
    }
    int n_kept = 0;
    for(int k = 0; k < n_idx; k++) {if(slot[k] >= 0) {slot[n_kept++] = slot[k];}}
    needs_drift.resize((size_t) n_kept);
    const int n_need = (int) needs_drift.size();
    /* Hand back the compaction this routine already had to perform, so a caller
       with per-advanced-particle follow-up does not recompute it. memmove, not
       memcpy: the caller is allowed to point this at its own input array, and the
       compaction is a leftward move within it. */
    if(out_n_drifted) {*out_n_drifted = n_need;}
    if(out_drifted && n_need > 0) {
        memmove(out_drifted, needs_drift.data(), (size_t) n_need * sizeof(int));
    }
    if(n_need <= 0) {return drift_batch_status_();}

    /* Tiny-N and everything below the offload threshold stays exactly as it was. */
    if(n_need < GPU_MIN_PARTICLES_FOR_DRIFT_OFFLOAD) {
        drift_particles_host_(needs_drift.data(), n_need, time1);
        return drift_batch_status_();
    }

    GIZMO_GPU_ENSURE_ALL_FRESH();

    /* Built on the host, captured by value: the owners are host globals and a
       kernel that reached for one would read host memory silently. */
    struct EosTableView eos_tables = eos_tables_view();

    struct DriftKickTableView tables;
    if(drift_kick_table_mirror_refresh(&drift_kick_table_dev_, &tables) != 0) {
        drift_particles_host_(needs_drift.data(), n_need, time1);
        return drift_batch_status_();
    }

    /* Drift the canonical arrays where they already live.  Same body, same particles,
       same physics as the staged route below -- the difference is that nothing is
       copied: the index list crosses at 4 bytes per particle instead of a full
       particle_data + gas_cell_data round trip each way.  Measured on Frontier, the
       staged round trip is dominated by the host-side gather and scatter of whole
       structs, which this removes outright rather than making cheaper.

       Legal only under the allocation model that makes the canonical arrays
       device-visible: allocate.cc serves P and CellP from gpu_particles_uvm_alloc,
       which allocates in GIZMO_KOKKOS_SHARED_SPACE, and the arena holds an alias of
       them rather than a copy.  The test is on that model rather than on a pointer,
       because a pointer carries no record of the space it came from.  Where the model
       does not hold this block is not entered and the staged route answers, out of
       buffers it owns. */
    if(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                  GIZMO_KOKKOS_SHARED_SPACE>::accessible)
    {
        int *idx_dev = (int *) gizmo_gpu_alloc_shared((size_t)n_need * sizeof(int), "drift_inplace_idx");
        if(!idx_dev) {   /* no shared memory for the index list: the particles still
                            have to be drifted, so the host loop takes them. */
            drift_particles_host_(needs_drift.data(), n_need, time1);
            return drift_batch_status_();
        }
        memcpy(idx_dev, needs_drift.data(), (size_t)n_need * sizeof(int));
        /* Captured by value: a lambda reaching for the P/CellP globals would read
           host memory silently on device. */
        struct particle_data *kp = P;
        struct gas_cell_data *kc = CellP;
        const int *kidx = idx_dev;
        const double dc0 = DomainCorner[0], dc1 = DomainCorner[1], dc2 = DomainCorner[2], dlen = DomainLen;
        /* Counts the particles this drift moved outside the domain extent, the same test the host
           drift makes; returned by the reduction, so no device flag has to be read back. */
        const int n_outside_extent = gizmo_gpu_kernel_launch_count("drift_particles_inplace", n_need, KOKKOS_LAMBDA(int k, int &n_outside) {
            drift_particle_impl(kidx[k], time1, kp, kc, &tables, &eos_tables);
            const int i = kidx[k];
            if(position_outside_domain_extent(kp[i].Pos[0], kp[i].Pos[1], kp[i].Pos[2], dc0, dc1, dc2, dlen)) {n_outside++;}
        });
        if(n_outside_extent > 0) {DomainExtentOutgrownLocal = 1;}
        Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(idx_dev);
        gpu_particles_arena_invalidate();   /* P/CellP mutated in place; arena stale */
        return drift_batch_status_();
    }

    const int batch_cap = (n_need < GPU_DRIFT_BATCH_SIZE) ? n_need : GPU_DRIFT_BATCH_SIZE;
    struct ParticleStagingBatch batch = {};
    if(!particle_staging_acquire(&batch, batch_cap)) {
        drift_particles_host_(needs_drift.data(), n_need, time1);
        return drift_batch_status_();
    }

    for(int batch_start = 0; batch_start < n_need; batch_start += GPU_DRIFT_BATCH_SIZE)
    {
        int batch_n = n_need - batch_start;
        if(batch_n > GPU_DRIFT_BATCH_SIZE) {batch_n = GPU_DRIFT_BATCH_SIZE;}

        if(!particle_staging_gather(&batch, needs_drift.data() + batch_start, batch_n, P, CellP)) {
            /* The gather staged nothing. These particles still have to be drifted:
               leaving them behind is a wrong position, and the sites that early-return
               on Ti_current would then never revisit them. */
            drift_particles_host_(needs_drift.data() + batch_start, n_need - batch_start, time1);
            break;
        }

        /* Slot j holds a particle, and holds a cell exactly while j < gas_count --
           the same condition the drift body's own Type==0 guard already expresses,
           so the kernel needs no extra gating for the non-gas slots. */
        struct particle_data *kp = batch.dev_P;
        struct gas_cell_data *kc = batch.dev_Cell;
        const double dc0 = DomainCorner[0], dc1 = DomainCorner[1], dc2 = DomainCorner[2], dlen = DomainLen;
        const int n_outside_extent = gizmo_gpu_kernel_launch_count("drift_particles", batch.count, KOKKOS_LAMBDA(int j, int &n_outside) {
            drift_particle_impl(j, time1, kp, kc, &tables, &eos_tables);
            if(position_outside_domain_extent(kp[j].Pos[0], kp[j].Pos[1], kp[j].Pos[2], dc0, dc1, dc2, dlen)) {n_outside++;}
        }, batch_start);
        if(n_outside_extent > 0) {DomainExtentOutgrownLocal = 1;}

        /* Synchronous by construction: the results are home before this returns, so
           the lazy-drift sites that early-return on Ti_current can never observe a
           particle whose Ti_current has advanced while its fields have not. */
        particle_staging_scatter(&batch, P, CellP);

        /* Assert the work product, not a return code: the launch wrapper reports a
           device failure by requesting a controlled stop, which drains at the next
           phase boundary, and returns nothing to its caller; a failed kernel would
           otherwise have its pre-drift staging scattered straight back. Ti_current is what says whether the body ran on a particle, and any
           that is still short of time1 is drifted here on the host. Without this the
           full-drift site would go on to certify the neighbour pool as current over
           positions that were never advanced, and the sweep that would have caught it
           is the one that certification suppresses. Costs one integer compare per
           staged slot on the path where nothing went wrong. */
        int n_undrifted = 0;
        for(int j = 0; j < batch.count; j++) {
            if(P[batch.index[j]].Ti_current != time1) {n_undrifted++;}
        }
        if(n_undrifted > 0) {
            std::vector<int> undrifted;
            undrifted.reserve((size_t) n_undrifted);
            for(int j = 0; j < batch.count; j++) {
                if(P[batch.index[j]].Ti_current != time1) {undrifted.push_back(batch.index[j]);}
            }
            drift_particles_host_(undrifted.data(), (int) undrifted.size(), time1);
        }
    }

    gpu_particles_arena_invalidate();   /* host P/CellP scattered; arena stale */
    particle_staging_release(&batch);
    return drift_batch_status_();
}
