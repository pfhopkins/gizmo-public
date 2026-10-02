/* gpu_moment_refresh.cc
 *
 * GPU implementation of the local-tree moment refresh that mirrors the
 * CPU code in force_refresh_node_moments (gravity/forcetree.cc:3826).
 * The CPU body is a 4-step pass (zero / particle-accumulate / bottom-up
 * propagate / normalize) over Nodes[tree_base .. tree_base+Numnodestree).
 *
 * On the GPU we use 5 device kernels with dependency-counter atomics so
 * the bottom-up walk parallelises:
 *   1) zero  : zero all moment scratch + init pending[]=0.
 *   2) count : each child node atomically increments its parent's
 *              pending counter (counts internal-node children).
 *   3) parts : each particle atomically accumulates into Father[i]'s
 *              scratch slot. Uses precomputed per-particle RT / sink /
 *              CR arrays (mirrors the gpu_gravtree.cc precomputed-input
 *              pattern).
 *   4) walk  : each thread starts at node k. If pending[k] != 0 it
 *              aborts (a descendant will reach this node). Otherwise it
 *              walks up via father, atomic-adding moments to the parent
 *              and atomic-decrementing the parent's pending counter; if
 *              the previous value was 1, this thread now owns the
 *              parent and continues with curr=parent.
 *   5) norm  : divide mass-weighted sums by mass (or fall back to node
 *              centre if mass==0); set/clear BITFLAG_MULTIPLEPARTICLES.
 *
 * After the kernels finish, gpu_moment_writeback_to_aos copies the SoA
 * back into Nodes[]/Extnodes[] (mechanical inverse of seed_node_), so
 * the CPU pseudo path (force_exchange_pseudodata +
 * force_treeupdate_pseudos) reads identical AoS values.
 *
 * All conditional payloads from force_refresh_node_moments are handled:
 *   GRAVTREE_CALCULATE_GAS_MASS_IN_NODE / RT_USE_GRAVTREE (+CHIMES) /
 *   RT_SEPARATELY_TRACK_LUMPOS / SINK_PHOTONMOMENTUM /
 *   COSMIC_RAY_SUBGRID_LEBRON / SINK_CALC_DISTANCES (with the nested
 *   N_SINK / sink_vel / sink_acc / MaxFeedbackVel) /
 *   ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION / DM_SCALARFIELD_SCREENING.
 *
 * Payload-list expansions are factored into KOKKOS_INLINE_FUNCTION
 * helpers so each payload appears at most a handful of times; the X-style
 * directives `MOMENT_FOR_EACH_*` keep zero / accum / norm / writeback in
 * lock-step.
 *
 * Memory model notes
 *   * Moment accumulators live in `Kokkos::View<...>` allocations on
 *     `Kokkos::DefaultExecutionSpace` (= device-local scratch). This avoids
 *     the MI250X HIPManaged atomic penalty.
 *     `Kokkos::deep_copy` flushes results into the SharedSpace SoA before
 *     the host-side AoS writeback.
 *   * Per-particle precomputed inputs (RT lum, sink lum, CR injection)
 *     live in SharedSpace, mirroring gpu_gravtree.cc:1216+.
 *   * Father[i] is mirrored into a SharedSpace int* once per call.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../declarations/gpu_error_check.h"
#include "../system/gpu_particles_arena.h"
#include "gpu_gravity_tree.h"
#include "forcetree.h"
#include "../core/timestep_functions.h"   /* dilation, for the motion bound */
#include "gravtree_moment_kernel.h"
#include "gravtree_moment_sources.h" /* shared host-only per-particle source-input physics gates (SSOT) */

#ifdef RT_USE_GRAVTREE
#include "../radiation/rt_functions.h"
#endif


/* ------------------------------------------------------------------ */
/* Atomic-max helper for floating-point types using a CAS loop.       */
/* Kokkos::atomic_max for floats/doubles is backend-dependent; CAS is */
/* portable across CUDA/HIP/OpenMP-host. Inlined as KOKKOS_INLINE.    */
/* ------------------------------------------------------------------ */
template <typename T>
KOKKOS_INLINE_FUNCTION static void
atomic_max_fp(T *addr, T val)
{
    T old = *addr;
    while(val > old) {
        T prev = Kokkos::atomic_compare_exchange(addr, old, val);
        if(prev == old) {return;}
        old = prev;
    }
}

/* ------------------------------------------------------------------ */
/* Vec3 atomic add (component-wise).                                  */
/* ------------------------------------------------------------------ */
template <typename T, typename U>
KOKKOS_INLINE_FUNCTION static void
atomic_add_vec3(Vec3<T> *addr, const Vec3<U>& val)
{
    Kokkos::atomic_add(&(*addr)[0], (T)val[0]);
    Kokkos::atomic_add(&(*addr)[1], (T)val[1]);
    Kokkos::atomic_add(&(*addr)[2], (T)val[2]);
}

/* Atomic write policy for the shared moment construction kernel (gravtree_moment_kernel.h). Kept in
 * this TU because the header carries no Kokkos dependency; gpu_moment_refresh injects it as the Ops
 * template argument while the topnode re-sum / force_add_element_to_tree use the header's plain policy. */
struct moment_atomic_ops {
    template <class T> KOKKOS_INLINE_FUNCTION static void add(T *dst, T v)  { Kokkos::atomic_add(dst, v); }
    template <class T> KOKKOS_INLINE_FUNCTION static void add_vec3(Vec3<T> *dst, const Vec3<T>& v) { atomic_add_vec3<T, T>(dst, v); }
    template <class T> KOKKOS_INLINE_FUNCTION static void fmax(T *dst, T v) { atomic_max_fp<T>(dst, v); }
    KOKKOS_INLINE_FUNCTION static void add_long(long *dst, long v) { Kokkos::atomic_add(dst, v); }
    KOKKOS_INLINE_FUNCTION static void add_int (int  *dst, int  v) { Kokkos::atomic_add(dst, v); }
    KOKKOS_INLINE_FUNCTION static void or_uint(unsigned int *dst, unsigned int v) { Kokkos::atomic_or(dst, v); }
};

/* ------------------------------------------------------------------ */
/* GPU-callable accessor for the cached force-softening kernel radius. */
/* Single source of truth: compute_force_softening_kernel_radius() in  */
/* forcetree.cc; populated by compute_all_force_softening() at startup */
/* and at the start of every gravity_tree() call.                      */
/* ------------------------------------------------------------------ */
KOKKOS_INLINE_FUNCTION static double
gpu_force_softening_kernelradius(const struct particle_data *Pp, int p)
{
    return Pp[p].ForceSoftening;
}

/* =================================================================== */
/* Per-particle precomputed inputs (RT, sink, CR). Same pattern as     */
/* gpu_gravtree.cc precomputed-input pattern. Allocated in SharedSpace so the device */
/* kernel can read them, but populated on host because the CPU helpers  */
/* (rt_get_source_luminosity, sink_lum_bol, cr_get_source_injection_rate)*/
/* are not GPU-callable.                                                */
/* =================================================================== */
/* Exhaustion is reported by returning NULL, which the pool guards below -- free
 * the pool, hand the caller a soft bad-stop -- are written against. The space is
 * a template parameter here because the raw mirror pool picks its space by
 * config; on builds where the spaces coincide the branches collapse. */
template <class Space>
static void *mr_alloc(const char *label, size_t bytes)
{
    if constexpr (std::is_same<Space, GIZMO_KOKKOS_SHARED_SPACE>::value) {
        return gizmo_gpu_alloc_shared(bytes, label);
    } else if constexpr (std::is_same<Space, Kokkos::HostSpace>::value) {
        return gizmo_gpu_alloc_host(bytes, label);
    } else {
        return gizmo_gpu_alloc_device(bytes, label);
    }
}

namespace {
struct precomputed_t {
#ifdef RT_USE_GRAVTREE
    MyFloat *src_lum;             /* [NumPart * N_RT_FREQ_BINS] */
#ifdef CHIMES_STELLAR_FLUXES
    double  *src_lum_G0;          /* [NumPart * CHIMES_LOCAL_UV_NBINS] */
    double  *src_lum_ion;         /* [NumPart * CHIMES_LOCAL_UV_NBINS] */
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    MyFloat       *bh_lum;        /* [NumPart] */
    Vec3<MyFloat> *bh_angle;      /* [NumPart] */
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyFloat *cr_inject;           /* [NumPart] */
#endif
};

/* Persistent grow-and-keep pool for the per-particle source-input buffers,
 * mirroring the gpu_gravity_tree.cc SoA idiom: reused across refreshes,
 * reallocated only when NumPart grows, freed when the tree is freed
 * (gpu_moment_refresh_release).  Removes per-refresh allocator churn. */
static precomputed_t pre_persist_ = {};
static int           pre_cap_     = 0;   /* current capacity in particles */
static void precompute_free_(precomputed_t& pre);   /* defined just below */

/* Grow the persistent source-input pool to hold at least N particles (no shrink).
 * On allocation failure frees the whole pool and leaves capacity invalid (0). */
static int precompute_ensure_(int N)
{
    if(pre_cap_ >= N) {return 0;}        /* capacity already sufficient: reuse */
    precompute_free_(pre_persist_);      /* drop the smaller pool before growing */
#if defined(RT_USE_GRAVTREE) || defined(SINK_PHOTONMOMENTUM) || defined(COSMIC_RAY_SUBGRID_LEBRON)
    precomputed_t& pre = pre_persist_;
#ifdef RT_USE_GRAVTREE
    pre.src_lum = (MyFloat *) mr_alloc<GIZMO_KOKKOS_SHARED_SPACE>("treescratch_moment_srclum", (long)N * N_RT_FREQ_BINS * sizeof(MyFloat));
    if(!pre.src_lum) {printf("gpu_moment_refresh: src_lum alloc failed\n"); endrun(913301); precompute_free_(pre); pre_cap_ = 0; return 1;}
#ifdef CHIMES_STELLAR_FLUXES
    pre.src_lum_G0  = (double *) mr_alloc<GIZMO_KOKKOS_SHARED_SPACE>("treescratch_moment_srclum_g0", (long)N * CHIMES_LOCAL_UV_NBINS * sizeof(double));
    pre.src_lum_ion = (double *) mr_alloc<GIZMO_KOKKOS_SHARED_SPACE>("treescratch_moment_srclum_ion", (long)N * CHIMES_LOCAL_UV_NBINS * sizeof(double));
    if(!pre.src_lum_G0 || !pre.src_lum_ion) {printf("gpu_moment_refresh: CHIMES alloc failed\n"); endrun(913302); precompute_free_(pre); pre_cap_ = 0; return 1;}
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    pre.bh_lum   = (MyFloat *)       mr_alloc<GIZMO_KOKKOS_SHARED_SPACE>("treescratch_moment_bhlum", (long)N * sizeof(MyFloat));
    pre.bh_angle = (Vec3<MyFloat> *) mr_alloc<GIZMO_KOKKOS_SHARED_SPACE>("treescratch_moment_bhangle", (long)N * sizeof(Vec3<MyFloat>));
    if(!pre.bh_lum || !pre.bh_angle) {printf("gpu_moment_refresh: bh_lum alloc failed\n"); endrun(913303); precompute_free_(pre); pre_cap_ = 0; return 1;}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    pre.cr_inject = (MyFloat *) mr_alloc<GIZMO_KOKKOS_SHARED_SPACE>("treescratch_moment_crinject", (long)N * sizeof(MyFloat));
    if(!pre.cr_inject) {printf("gpu_moment_refresh: cr_inject alloc failed\n"); endrun(913304); precompute_free_(pre); pre_cap_ = 0; return 1;}
#endif
#endif
    pre_cap_ = N;
    return 0;
}

/* Reinitialize the active [0,N) range every refresh: bulk-zero the buffers (the
 * pool is reused, so prior contents must be cleared), then write the active
 * source payloads via the shared SSOT helper.  The full-NumPart scan is
 * unchanged (source sparsification is a separate, profiling-gated change). */
static void precompute_fill_(int N)
{
#if defined(RT_USE_GRAVTREE) || defined(SINK_PHOTONMOMENTUM) || defined(COSMIC_RAY_SUBGRID_LEBRON)
    precomputed_t& pre = pre_persist_;
#ifdef RT_USE_GRAVTREE
    memset(pre.src_lum, 0, (long)N * N_RT_FREQ_BINS * sizeof(MyFloat));
#ifdef CHIMES_STELLAR_FLUXES
    memset(pre.src_lum_G0,  0, (long)N * CHIMES_LOCAL_UV_NBINS * sizeof(double));
    memset(pre.src_lum_ion, 0, (long)N * CHIMES_LOCAL_UV_NBINS * sizeof(double));
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    memset(pre.bh_lum,   0, (long)N * sizeof(MyFloat));
    memset(pre.bh_angle, 0, (long)N * sizeof(Vec3<MyFloat>));
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    memset(pre.cr_inject, 0, (long)N * sizeof(MyFloat));
#endif

    /* Single host pass: gated physics in the shared SSOT helper, then copy ONLY
     * the active entries into this venue's SharedSpace arrays (already bulk-zeroed
     * above), matching the legacy per-section write pattern. */
    for(int p = 0; p < N; p++) {
        struct gravtree_source_inputs_t in;
        gravtree_fill_particle_source_inputs(p, P, CellP, &in);
#ifdef RT_USE_GRAVTREE
        if(in.rt_active) {
            int kf;
            for(kf = 0; kf < N_RT_FREQ_BINS; kf++) {pre.src_lum[(long)p * N_RT_FREQ_BINS + kf] = in.src_lum[kf];}
#ifdef CHIMES_STELLAR_FLUXES
            for(kf = 0; kf < CHIMES_LOCAL_UV_NBINS; kf++) {
                pre.src_lum_G0[(long)p * CHIMES_LOCAL_UV_NBINS + kf]  = in.src_lum_G0[kf];
                pre.src_lum_ion[(long)p * CHIMES_LOCAL_UV_NBINS + kf] = in.src_lum_ion[kf];
            }
#endif
        }
#endif
#ifdef SINK_PHOTONMOMENTUM
        if(in.bh_active) {
            pre.bh_lum[p]   = in.bh_lum;
            pre.bh_angle[p] = in.bh_angle;
        }
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        if(in.cr_inject != 0) {pre.cr_inject[p] = in.cr_inject;}
#endif
    }
#else
    (void) N;
#endif
}

static void precompute_free_(precomputed_t& pre)
{
#ifdef RT_USE_GRAVTREE
    if(pre.src_lum) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(pre.src_lum); pre.src_lum = NULL;}
#ifdef CHIMES_STELLAR_FLUXES
    if(pre.src_lum_G0)  {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(pre.src_lum_G0);  pre.src_lum_G0  = NULL;}
    if(pre.src_lum_ion) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(pre.src_lum_ion); pre.src_lum_ion = NULL;}
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
    if(pre.bh_lum)   {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(pre.bh_lum);   pre.bh_lum   = NULL;}
    if(pre.bh_angle) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(pre.bh_angle); pre.bh_angle = NULL;}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    if(pre.cr_inject) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(pre.cr_inject); pre.cr_inject = NULL;}
#endif
}
} /* anonymous namespace */

/* ------------------------------------------------------------------ */
/* Father[i] mirror in SharedSpace.  Persistent grow-and-keep pool      */
/* (same idiom as the source-input pool above): grown only when NumPart */
/* grows, refilled from Father[] every refresh, freed when the tree is  */
/* freed (gpu_moment_refresh_release).                                   */
/* ------------------------------------------------------------------ */
static int *fmirror_persist_ = NULL;
static int  fmirror_cap_     = 0;   /* current capacity in particles */

static int *father_mirror_ensure_(int N)
{
    if(fmirror_cap_ < N) {
        if(fmirror_persist_) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(fmirror_persist_); fmirror_persist_ = NULL;}
        fmirror_persist_ = (int *) mr_alloc<GIZMO_KOKKOS_SHARED_SPACE>("treescratch_moment_fmirror", (long)N * sizeof(int));
        if(!fmirror_persist_) {printf("gpu_moment_refresh: Father mirror alloc failed\n"); endrun(913305); fmirror_cap_ = 0; return NULL;}
        fmirror_cap_ = N;
    }
    /* Father[] only has All.TreeParticleSlots entries, while N counts every particle P[] holds, and
     * that can be the larger of the two once the particle capacity is raised mid-step.  Slots beyond
     * the tree's own are filled with the value the tree already uses for a particle it does not
     * contain, so consumers treat them as absent instead of reading past the array.  Imported ghosts
     * need no special handling here: their slots were not in the tree at its last build, so they
     * already carry that same value. */
    const int ncopy = (N < All.TreeParticleSlots) ? N : All.TreeParticleSlots;
    for(int i = 0; i < ncopy; i++) {fmirror_persist_[i] = Father[i];}
    for(int i = ncopy; i < N; i++) {fmirror_persist_[i] = -1;}
    return fmirror_persist_;
}

/* =================================================================== */
/* Bundle scratch Views + factor per-payload reductions into            */
/* KOKKOS_INLINE_FUNCTION helpers.  gpu_moment_refresh_topnodes()       */
/* reuses the exact same helpers for the topnode-subtree ancestor       */
/* re-sum, avoiding ~250 lines of duplicated #ifdef-guarded payload     */
/* logic.  All behavior here is bit-identical to the prior inline       */
/* implementation; this is purely a refactor.                          */
/* =================================================================== */
namespace {

using MrExSpace   = Kokkos::DefaultExecutionSpace;
using MrMemSpace  = MrExSpace::memory_space;
using MrUnmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

/* The scratch is held as UNMANAGED Views over a persistent raw kokkos_malloc
 * pool (mr_rawpool_t below): the pool gives grow-and-keep lifetime + finalize
 * safety (raw pointers, freed via explicit release, never a tracked View
 * outliving Kokkos::finalize), while the Views keep the kernel/accessor code
 * unchanged.  The Views own nothing; they are rebound (pointer + extent n) each
 * refresh by mr_scratch_bind_. */
struct mr_scratch_t {
    Kokkos::View<int*,           MrMemSpace, MrUnmanaged> pending;
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> mass;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> s;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> vs;
    Kokkos::View<long*,          MrMemSpace, MrUnmanaged> Npart;
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> hmax;
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> vmax;
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> divVmax;
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> maxsoft;
    Kokkos::View<unsigned int*,  MrMemSpace, MrUnmanaged> bitflags;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> gasmass;
#endif
#ifdef RT_USE_GRAVTREE
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> stellar_lum;  /* [n * N_RT_FREQ_BINS] */
#ifdef CHIMES_STELLAR_FLUXES
    Kokkos::View<double*,        MrMemSpace, MrUnmanaged> chimes_G0;    /* [n * CHIMES_LOCAL_UV_NBINS] */
    Kokkos::View<double*,        MrMemSpace, MrUnmanaged> chimes_ion;
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> rt_s;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> rt_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> sink_lum;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> sink_lum_grad;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> cr_inject;
#endif
#ifdef SINK_CALC_DISTANCES
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> sink_mass;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
    Kokkos::View<int*,           MrMemSpace, MrUnmanaged> N_SINK;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> sink_vel;
#endif
#if defined(SPECIAL_POINT_MOTION)
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> sink_acc;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> max_fbvel;
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> tidal;        /* [n * 6] */
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Kokkos::View<MyGravFloat*,   MrMemSpace, MrUnmanaged> mass_dm;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> s_dm;
    Kokkos::View<Vec3<MyGravFloat>*, MrMemSpace, MrUnmanaged> vs_dm;
#endif
};

/* Persistent raw scratch pool (MrMemSpace) backing the unmanaged scratch Views.
 * Grow-and-keep; freed via gpu_moment_refresh_release().  Raw pointers, not
 * Views, so nothing tracked outlives Kokkos::finalize. */
struct mr_rawpool_t {
    int*               pending;
    MyGravFloat*       mass;
    Vec3<MyGravFloat>* s;
    Vec3<MyGravFloat>* vs;
    long*              Npart;
    MyGravFloat*       hmax;
    MyGravFloat*       vmax;
    MyGravFloat*       divVmax;
    MyGravFloat*       maxsoft;
    unsigned int*      bitflags;
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    MyGravFloat*       gasmass;
#endif
#ifdef RT_USE_GRAVTREE
    MyGravFloat*       stellar_lum;
#ifdef CHIMES_STELLAR_FLUXES
    double*            chimes_G0;
    double*            chimes_ion;
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Vec3<MyGravFloat>* rt_s;
    Vec3<MyGravFloat>* rt_vs;
#endif
#ifdef SINK_PHOTONMOMENTUM
    MyGravFloat*       sink_lum;
    Vec3<MyGravFloat>* sink_lum_grad;
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyGravFloat*       cr_inject;
#endif
#ifdef SINK_CALC_DISTANCES
    MyGravFloat*       sink_mass;
    Vec3<MyGravFloat>* sink_pos;
#if defined(SINK_NODE_MOTION_TRACKED)
    int*               N_SINK;
    Vec3<MyGravFloat>* sink_vel;
#endif
#if defined(SPECIAL_POINT_MOTION)
    Vec3<MyGravFloat>* sink_acc;
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    MyGravFloat*       max_fbvel;
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    MyGravFloat*       tidal;
#endif
#ifdef DM_SCALARFIELD_SCREENING
    MyGravFloat*       mass_dm;
    Vec3<MyGravFloat>* s_dm;
    Vec3<MyGravFloat>* vs_dm;
#endif
};

static mr_rawpool_t raw_     = {};
static int          raw_cap_ = 0;   /* capacity in nodes */

/* Free every allocated pool buffer.  Visits the same #ifdef field set as the
 * struct; idempotent (NULL-checked), used both on grow and at teardown. */
static void mr_rawpool_free_(void)
{
#define MRFREE(p) do { if(raw_.p) { Kokkos::kokkos_free<MrMemSpace>(raw_.p); raw_.p = NULL; } } while(0)
    MRFREE(pending); MRFREE(mass);    MRFREE(s);       MRFREE(vs);     MRFREE(Npart);
    MRFREE(hmax);    MRFREE(vmax);    MRFREE(divVmax); MRFREE(maxsoft); MRFREE(bitflags);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    MRFREE(gasmass);
#endif
#ifdef RT_USE_GRAVTREE
    MRFREE(stellar_lum);
#ifdef CHIMES_STELLAR_FLUXES
    MRFREE(chimes_G0); MRFREE(chimes_ion);
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MRFREE(rt_s); MRFREE(rt_vs);
#endif
#ifdef SINK_PHOTONMOMENTUM
    MRFREE(sink_lum); MRFREE(sink_lum_grad);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MRFREE(cr_inject);
#endif
#ifdef SINK_CALC_DISTANCES
    MRFREE(sink_mass); MRFREE(sink_pos);
#if defined(SINK_NODE_MOTION_TRACKED)
    MRFREE(N_SINK); MRFREE(sink_vel);
#endif
#if defined(SPECIAL_POINT_MOTION)
    MRFREE(sink_acc);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    MRFREE(max_fbvel);
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    MRFREE(tidal);
#endif
#ifdef DM_SCALARFIELD_SCREENING
    MRFREE(mass_dm); MRFREE(s_dm); MRFREE(vs_dm);
#endif
#undef MRFREE
}

/* Grow the raw pool to hold at least n nodes (no shrink).  On any allocation
 * failure frees the whole pool and leaves capacity invalid (0). */
static int mr_rawpool_ensure_(int n)
{
    if(raw_cap_ >= n) {return 0;}        /* capacity already sufficient: reuse */
    mr_rawpool_free_();                  /* drop the smaller pool before growing */
    int fail = 0;
/* Stop asking once one request has failed: the node is out of memory at the size this pool
 * needs, and the remaining fields would only add to the pressure before the caller stops. */
#define MRALLOC(p, count) do { \
        if(!fail) { \
            raw_.p = (decltype(raw_.p)) mr_alloc<MrMemSpace>("treescratch_moment_raw_" #p, (long)(count) * sizeof(*raw_.p)); \
            if(!raw_.p) { fail = 1; } } } while(0)
    MRALLOC(pending, n); MRALLOC(mass, n);    MRALLOC(s, n);       MRALLOC(vs, n);     MRALLOC(Npart, n);
    MRALLOC(hmax, n);    MRALLOC(vmax, n);    MRALLOC(divVmax, n); MRALLOC(maxsoft, n); MRALLOC(bitflags, n);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    MRALLOC(gasmass, n);
#endif
#ifdef RT_USE_GRAVTREE
    MRALLOC(stellar_lum, (long)n * N_RT_FREQ_BINS);
#ifdef CHIMES_STELLAR_FLUXES
    MRALLOC(chimes_G0, (long)n * CHIMES_LOCAL_UV_NBINS); MRALLOC(chimes_ion, (long)n * CHIMES_LOCAL_UV_NBINS);
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MRALLOC(rt_s, n); MRALLOC(rt_vs, n);
#endif
#ifdef SINK_PHOTONMOMENTUM
    MRALLOC(sink_lum, n); MRALLOC(sink_lum_grad, n);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MRALLOC(cr_inject, n);
#endif
#ifdef SINK_CALC_DISTANCES
    MRALLOC(sink_mass, n); MRALLOC(sink_pos, n);
#if defined(SINK_NODE_MOTION_TRACKED)
    MRALLOC(N_SINK, n); MRALLOC(sink_vel, n);
#endif
#if defined(SPECIAL_POINT_MOTION)
    MRALLOC(sink_acc, n);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    MRALLOC(max_fbvel, n);
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    MRALLOC(tidal, (long)n * 6);
#endif
#ifdef DM_SCALARFIELD_SCREENING
    MRALLOC(mass_dm, n); MRALLOC(s_dm, n); MRALLOC(vs_dm, n);
#endif
#undef MRALLOC
    if(fail) {printf("gpu_moment_refresh: scratch pool alloc failed\n"); endrun(913306); mr_rawpool_free_(); raw_cap_ = 0; return 1;}
    raw_cap_ = n;
    return 0;
}

/* Bind the scratch Views (extent n) over the persistent raw pool, growing the
 * pool first.  Returns 1 on allocation failure (pool left freed/invalid). */
static int mr_scratch_bind_(mr_scratch_t& scr, int n)
{
    if(mr_rawpool_ensure_(n) != 0) {return 1;}
#define MRBIND(p, ext) scr.p = decltype(scr.p)(raw_.p, (ext))
    MRBIND(pending, n); MRBIND(mass, n);    MRBIND(s, n);       MRBIND(vs, n);     MRBIND(Npart, n);
    MRBIND(hmax, n);    MRBIND(vmax, n);    MRBIND(divVmax, n); MRBIND(maxsoft, n); MRBIND(bitflags, n);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    MRBIND(gasmass, n);
#endif
#ifdef RT_USE_GRAVTREE
    MRBIND(stellar_lum, (long)n * N_RT_FREQ_BINS);
#ifdef CHIMES_STELLAR_FLUXES
    MRBIND(chimes_G0, (long)n * CHIMES_LOCAL_UV_NBINS); MRBIND(chimes_ion, (long)n * CHIMES_LOCAL_UV_NBINS);
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    MRBIND(rt_s, n); MRBIND(rt_vs, n);
#endif
#ifdef SINK_PHOTONMOMENTUM
    MRBIND(sink_lum, n); MRBIND(sink_lum_grad, n);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MRBIND(cr_inject, n);
#endif
#ifdef SINK_CALC_DISTANCES
    MRBIND(sink_mass, n); MRBIND(sink_pos, n);
#if defined(SINK_NODE_MOTION_TRACKED)
    MRBIND(N_SINK, n); MRBIND(sink_vel, n);
#endif
#if defined(SPECIAL_POINT_MOTION)
    MRBIND(sink_acc, n);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    MRBIND(max_fbvel, n);
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    MRBIND(tidal, (long)n * 6);
#endif
#ifdef DM_SCALARFIELD_SCREENING
    MRBIND(mass_dm, n); MRBIND(s_dm, n); MRBIND(vs_dm, n);
#endif
#undef MRBIND
    return 0;
}

/* Adapters between this venue's scratch Views and the shared moment kernel's POD interface
 * (gravtree_moment_kernel.h).  mr_ref_ points a moment_node_ref at one node's scratch slots (the
 * array payloads as their per-node base); mr_child_accum_ reads a settled child slot into a by-value
 * accumulator for the bottom-up propagate.  Pure pointer/loads — bit-identical to the prior inline
 * scratch accesses. */
KOKKOS_INLINE_FUNCTION
static moment_node_ref<MyGravFloat> mr_ref_(const mr_scratch_t& scr, int k)
{
    moment_node_ref<MyGravFloat> r;
    r.mass     = &scr.mass(k);
    r.s        = &scr.s(k);
    r.vs       = &scr.vs(k);
    r.Npart    = &scr.Npart(k);
    r.hmax     = &scr.hmax(k);
    r.vmax     = &scr.vmax(k);
    r.divVmax  = &scr.divVmax(k);
    r.maxsoft  = &scr.maxsoft(k);
    r.bitflags = &scr.bitflags(k);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    r.gasmass  = &scr.gasmass(k);
#endif
#ifdef RT_USE_GRAVTREE
    r.stellar_lum = &scr.stellar_lum((long)k * N_RT_FREQ_BINS);
#ifdef CHIMES_STELLAR_FLUXES
    r.chimes_G0  = &scr.chimes_G0 ((long)k * CHIMES_LOCAL_UV_NBINS);
    r.chimes_ion = &scr.chimes_ion((long)k * CHIMES_LOCAL_UV_NBINS);
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    r.rt_s  = &scr.rt_s(k);
    r.rt_vs = &scr.rt_vs(k);
#endif
#ifdef SINK_PHOTONMOMENTUM
    r.sink_lum      = &scr.sink_lum(k);
    r.sink_lum_grad = &scr.sink_lum_grad(k);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    r.cr_inject = &scr.cr_inject(k);
#endif
#ifdef SINK_CALC_DISTANCES
    r.sink_mass = &scr.sink_mass(k);
    r.sink_pos  = &scr.sink_pos(k);
#if defined(SINK_NODE_MOTION_TRACKED)
    r.N_SINK   = &scr.N_SINK(k);
    r.sink_vel = &scr.sink_vel(k);
#endif
#if defined(SPECIAL_POINT_MOTION)
    r.sink_acc = &scr.sink_acc(k);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    r.max_fbvel = &scr.max_fbvel(k);
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    r.tidal = &scr.tidal((long)k * 6);
#endif
#ifdef DM_SCALARFIELD_SCREENING
    r.mass_dm = &scr.mass_dm(k);
    r.s_dm    = &scr.s_dm(k);
    r.vs_dm   = &scr.vs_dm(k);
#endif
    return r;
}

KOKKOS_INLINE_FUNCTION
static moment_node_accum<MyGravFloat> mr_child_accum_(const mr_scratch_t& scr, int curr)
{
    moment_node_accum<MyGravFloat> c = {};
    /* The types this child holds.  They live in its bitflags rather than in a payload of their
     * own, so unlike every other field they have to be read back from there -- without this the
     * union stops one level above the leaves and every node higher up claims to hold nothing. */
    c.type_mask = NODE_TYPE_PRESENCE(scr.bitflags(curr));
    c.mass    = scr.mass(curr);
    c.s       = scr.s(curr);
    c.vs      = scr.vs(curr);
    c.Npart   = scr.Npart(curr);
    c.hmax    = scr.hmax(curr);
    c.vmax    = scr.vmax(curr);
    c.divVmax = scr.divVmax(curr);
    c.maxsoft = scr.maxsoft(curr);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    c.gasmass = scr.gasmass(curr);
#endif
#ifdef RT_USE_GRAVTREE
    for(int b = 0; b < N_RT_FREQ_BINS; b++) {c.stellar_lum[b] = scr.stellar_lum((long)curr * N_RT_FREQ_BINS + b);}
#ifdef CHIMES_STELLAR_FLUXES
    for(int b = 0; b < CHIMES_LOCAL_UV_NBINS; b++) {
        c.chimes_G0 [b] = scr.chimes_G0 ((long)curr * CHIMES_LOCAL_UV_NBINS + b);
        c.chimes_ion[b] = scr.chimes_ion((long)curr * CHIMES_LOCAL_UV_NBINS + b);
    }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    c.rt_s  = scr.rt_s(curr);
    c.rt_vs = scr.rt_vs(curr);
#endif
#ifdef SINK_PHOTONMOMENTUM
    c.sink_lum      = scr.sink_lum(curr);
    c.sink_lum_grad = scr.sink_lum_grad(curr);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    c.cr_inject = scr.cr_inject(curr);
#endif
#ifdef SINK_CALC_DISTANCES
    c.sink_mass = scr.sink_mass(curr);
    c.sink_pos  = scr.sink_pos(curr);
#if defined(SINK_NODE_MOTION_TRACKED)
    c.N_SINK   = scr.N_SINK(curr);
    c.sink_vel = scr.sink_vel(curr);
#endif
#if defined(SPECIAL_POINT_MOTION)
    c.sink_acc = scr.sink_acc(curr);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    c.max_fbvel = scr.max_fbvel(curr);
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    for(int kk = 0; kk < 6; kk++) {c.tidal[kk] = scr.tidal((long)curr * 6 + kk);}
#endif
#ifdef DM_SCALARFIELD_SCREENING
    c.mass_dm = scr.mass_dm(curr);
    c.s_dm    = scr.s_dm(curr);
    c.vs_dm   = scr.vs_dm(curr);
#endif
    return c;
}

/* Helper: zero per-node payload accumulators at slot k.  Bitflags
 * preserves the topology bits (TOPLEVEL/INTERNAL_TOPLEVEL/DEPENDS_
 * ON_LOCAL_ELEMENT) and clears the rest, matching the CPU saved_
 * bitflags mask in force_refresh_node_moments (forcetree.cc:3841). */
KOKKOS_INLINE_FUNCTION
static void mr_zero_payloads_(const mr_scratch_t& scr, int k, unsigned int saved_bitflags)
{
    moment_accum_zero<MyGravFloat>(mr_ref_(scr, k), saved_bitflags);
}

/* Helper: bottom-up child→parent propagation.  Atomically adds curr's
 * accumulated moments + max-reduces hmax/vmax/divVmax/maxsoft into kp.
 * Used by Kernel 4 (walk-up); also reused for the topnode-subtree
 * ancestor re-sum after MPI pseudodata exchange. */
KOKKOS_INLINE_FUNCTION
static void mr_propagate_to_parent_(const mr_scratch_t& scr, int curr, int kp)
{
    moment_accum_add_child_raw<moment_atomic_ops, MyGravFloat>(mr_ref_(scr, kp), mr_child_accum_(scr, curr));
}

/* Finalize a node: turn the mass-weighted sums into stored normalized payloads + patch
 * BITFLAG_MULTIPLEPARTICLES.  `center` is the node's geometric center (COM fallback when the weight
 * is zero).  Delegates to the shared moment_finalize (gravtree_moment_kernel.h); the scratch
 * bitflags carry the topology bits, so the multiple-particles bit is patched in place via the ref. */
KOKKOS_INLINE_FUNCTION
static void mr_normalize_payloads_(const mr_scratch_t& scr, int k, const Vec3<MyFloat>& center)
{
    moment_finalize<MyGravFloat>(mr_ref_(scr, k),
                                 Vec3<MyGravFloat>{(MyGravFloat) center[0],
                                                   (MyGravFloat) center[1],
                                                   (MyGravFloat) center[2]});
}

/* Pending-counter init launcher.  For each k in [0..n_iter), atomically
 * increments pending[father(k) - tree_base] when father falls within
 * [range_lo, range_hi).  The topnode-subtree variant passes a
 * restricted iteration range + filter. */
static void mr_pending_init_launch_(Kokkos::View<int*, MrMemSpace, MrUnmanaged> pending,
                                     int *father_soa, int tree_base,
                                     int n_iter, int range_lo, int range_hi)
{
    Kokkos::parallel_for("mr_pending_init", n_iter, KOKKOS_LAMBDA(int k) {
        int f = father_soa[k];
        if(f >= range_lo && f < range_hi) {
            Kokkos::atomic_inc(&pending(f - tree_base));
        }
    });
}

} /* anonymous namespace (moment-refresh helpers) */

/* =================================================================== */
/*  Main entry: gpu_moment_refresh                                      */
/* =================================================================== */
extern "C" int gpu_moment_refresh(int active_root_node)
{
    (void) active_root_node; /* reserved hook, currently unused */
    if(Numnodestree <= 0) {return 0;}
    /* Defensive idempotent sync of AllDevice before any device-side All.* read. */
    GIZMO_GPU_ENSURE_ALL_FRESH();

    int n          = Numnodestree;       /* number of internal nodes [0..n) in SoA */
    int tree_base    = All.TreeNodeIndexBase;
    int N          = NumPart;

    /* ---------------- 1. arenas / SoA / per-particle precompute -------- */
    gpu_particles_arena_set_site("gpu_moment_refresh(pre-precompute)");
    gpu_particles_arena_acquire(N, P, CellP);
    struct particle_data *P_dev = gpu_particles_arena_P();
    struct gas_cell_data *CellP_dev = gpu_particles_arena_CellP();
    if(!P_dev) {printf("gpu_moment_refresh: arena P_dev null\n"); return 1;}

    /* Tree SoA storage must already be allocated; the GPU build pipeline
     * (gpu_nextnode_backup_suns + topology emit/finalize kernels) populated
     * topology fields (sibling, nextnode, father, center, len, bitflags)
     * directly in the SoA before this kernel runs.  We overwrite the
     * moments here. */
    int min_nodes = MaxNodes + 1;
    gpu_gravity_tree_acquire(min_nodes, Nodes_base, Extnodes_base);
    struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa) {printf("gpu_moment_refresh: SoA null\n"); return 1;}

    if(precompute_ensure_(N) != 0) {return 1;}   /* soft bad-stop: ensure freed the pool on failure; caller (forcetree) drains at next poll */
    precompute_fill_(N);
    int *fmirror = father_mirror_ensure_(N);
    if(!fmirror) {return 1;}   /* soft bad-stop: persistent pools kept; caller drains at next poll */

    /* ---------------- 2. device-local scratch --------------------------- */
    /* Bundled into mr_scratch_t so the per-payload helpers (zero /
     * propagate / normalize / pending_init) above can be reused unchanged
     * by gpu_moment_refresh_topnodes(). */
    mr_scratch_t scr;
    if(mr_scratch_bind_(scr, n) != 0) {return 1;}   /* soft bad-stop: pool freed by ensure on failure; caller drains at next poll */

    /* SoA pointers captured into the lambdas (raw pointers; SharedSpace is
     * accessible from device). The SoA `father`/`bitflags` carry topology
     * already seeded by gpu_gravity_tree_acquire above. */
    int          *father_soa   = soa->father;
    unsigned int *bitflags_soa = soa->bitflags;
    Vec3<MyFloat> *center_soa  = soa->center;

    /* ---------------- Kernel 1: zero scratch + pending ----------------- */
    /* pending is initialized to 1 (resident-thread slot); Kernel 2 adds
     * one more for each internal-node child.  In Kernel 4 each thread
     * claims its starting node with atomic_fetch_sub; only the last to
     * arrive (prev==1) proceeds upward.  The +1 prevents the race where
     * a late resident thread and the bottom-up walk both see pending==0
     * and double-count the node. */
    Kokkos::deep_copy(scr.pending, 1);
    Kokkos::parallel_for("mr_zero_payloads", n, KOKKOS_LAMBDA(int k) {
        mr_zero_payloads_(scr, k, bitflags_soa[k]);
    });

    /* ---------------- Kernel 2: pending counter init ------------------- */
    mr_pending_init_launch_(scr.pending, father_soa, tree_base, n,
                             tree_base, tree_base + n);

    /* ---------------- Kernel 3: particle pass -------------------------- */
    double maxKR = All.MaxKernelRadius;
#ifdef SINK_PHOTONMOMENTUM
    MyFloat       *bh_lum_   = pre_persist_.bh_lum;
    Vec3<MyFloat> *bh_angle_ = pre_persist_.bh_angle;
#endif
#ifdef RT_USE_GRAVTREE
    MyFloat *src_lum_ = pre_persist_.src_lum;
#ifdef CHIMES_STELLAR_FLUXES
    double  *src_lum_G0_  = pre_persist_.src_lum_G0;
    double  *src_lum_ion_ = pre_persist_.src_lum_ion;
#endif
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    MyFloat *cr_inj_ = pre_persist_.cr_inject;
#endif
    Kokkos::parallel_for("mr_part_accum", N, KOKKOS_LAMBDA(int i) {
        int f = fmirror[i];
        if(f < tree_base || f >= tree_base + n) {return;}
        int k = f - tree_base;
        const struct particle_data *pa = &P_dev[i];

        /* Load this leaf's fields + per-particle precomputes into the shared kernel's POD (native
         * double; weighted products are formed at double precision inside the helper).  The gates and
         * formulas live in moment_source_from_particle. */
        moment_particle_src<MyGravFloat> ps = {};
        ps.mass = pa->Mass;
        for(int kk = 0; kk < 3; kk++) {ps.pos[kk] = pa->Pos[kk]; ps.vel[kk] = pa->Vel[kk];}
        ps.motion_bound      = particle_motion_speed_bound(i, P_dev, CellP_dev);
        ps.type              = (int) pa->Type;
        ps.kernel_radius     = pa->KernelRadius;
        ps.max_kernel_radius = maxKR;
        ps.force_softening   = gpu_force_softening_kernelradius(P_dev, i);
        ps.particle_divvel   = pa->Particle_DivVel;
#if defined(SINK_ALPHADISK_ACCRETION) && defined(RT_USE_TREECOL_FOR_NH)
        ps.sink_mass_reservoir = pa->Sink_Mass_Reservoir;
#endif
#ifdef RT_USE_GRAVTREE
        for(int b = 0; b < N_RT_FREQ_BINS; b++) {ps.src_lum[b] = src_lum_[(long)i * N_RT_FREQ_BINS + b];}
#ifdef CHIMES_STELLAR_FLUXES
        for(int b = 0; b < CHIMES_LOCAL_UV_NBINS; b++) {
            ps.src_lum_G0 [b] = src_lum_G0_ [(long)i * CHIMES_LOCAL_UV_NBINS + b];
            ps.src_lum_ion[b] = src_lum_ion_[(long)i * CHIMES_LOCAL_UV_NBINS + b];
        }
#endif
#endif
#ifdef SINK_PHOTONMOMENTUM
        ps.bh_lum = bh_lum_[i];
        for(int kk = 0; kk < 3; kk++) {ps.bh_angle[kk] = bh_angle_[i][kk];}
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        ps.cr_inject = cr_inj_[i];
#endif
#if defined(SPECIAL_POINT_MOTION)
        for(int kk = 0; kk < 3; kk++) {ps.acc_prevstep[kk] = pa->Acc_Total_PrevStep[kk];}
#endif
#if defined(SINK_CALC_DISTANCES) && defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
        ps.max_feedback_vel = pa->MaxFeedbackVel;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        for(int kk = 0; kk < 6; kk++) {ps.tidal_prevstep[kk] = pa->tidal_tensorps_prevstep.data[kk];}
#endif

        moment_accum_add_particle<moment_atomic_ops, MyGravFloat>(mr_ref_(scr, k), ps);
    });

    /* ---------------- Kernel 4: bottom-up walk via father chain -------- */
    /* The counter hands a node to whichever thread arrives last, and that thread then reads its
     * children's accumulated payloads with plain loads. The payload adds are atomic, but Kokkos
     * atomics are relaxed: the decrement orders nothing, so on a device whose caches are not
     * coherent the winner can read values that predate the losers' contributions. x86 and OpenMP
     * get the ordering for free from the hardware model, which is why only the device backend
     * miscounts. Make it explicit: release before signalling a parent, acquire after winning one. */
    Kokkos::parallel_for("mr_walk_up", n, KOKKOS_LAMBDA(int k0) {
        if(Kokkos::atomic_fetch_sub(&scr.pending(k0), 1) != 1) {return;}
        Kokkos::memory_fence();   /* acquire: this node's children are visible before we read them */
        int curr = k0;
        while(true) {
            int f = father_soa[curr];
            if(f < tree_base || f >= tree_base + n) {return;}
            int kp = f - tree_base;
            mr_propagate_to_parent_(scr, curr, kp);
            Kokkos::memory_fence();   /* release: our contribution lands before the count can hand kp on */
            int prev = Kokkos::atomic_fetch_sub(&scr.pending(kp), 1);
            if(prev != 1) {return;}
            Kokkos::memory_fence();   /* acquire: all of kp's children are visible before the next lap */
            curr = kp;
        }
    });

    /* ---------------- Kernel 5: normalize ------------------------------ */
    Kokkos::parallel_for("mr_normalize", n, KOKKOS_LAMBDA(int k) {
        mr_normalize_payloads_(scr, k, center_soa[k]);
    });

    Kokkos::fence();
    gizmo_gpu_check_last_error("mr_moment_refresh", n);

    /* ---------------- Bulk deep-copy scratch → SharedSpace SoA -------- */
    /* Wrap raw SoA pointers in unmanaged Views and deep_copy. */
    using UV_unmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
    using ShSpace = Kokkos::View<MyGravFloat*, GIZMO_KOKKOS_SHARED_SPACE, UV_unmanaged>;
    using ShVec3  = Kokkos::View<Vec3<MyGravFloat>*, GIZMO_KOKKOS_SHARED_SPACE, UV_unmanaged>;
    using ShU32   = Kokkos::View<unsigned int*, GIZMO_KOKKOS_SHARED_SPACE, UV_unmanaged>;
    using ShLong  = Kokkos::View<long*, GIZMO_KOKKOS_SHARED_SPACE, UV_unmanaged>;
    using ShDbl   = Kokkos::View<double*, GIZMO_KOKKOS_SHARED_SPACE, UV_unmanaged>;
    using ShInt   = Kokkos::View<int*, GIZMO_KOKKOS_SHARED_SPACE, UV_unmanaged>;

    Kokkos::deep_copy(ShSpace(soa->mass,    n), scr.mass);
    Kokkos::deep_copy(ShVec3 (soa->s,       n), scr.s);
    Kokkos::deep_copy(ShVec3 (soa->node_vs, n), scr.vs);
    Kokkos::deep_copy(ShLong (soa->N_part,  n), scr.Npart);
    Kokkos::deep_copy(ShSpace(soa->hmax,    n), scr.hmax);
    Kokkos::deep_copy(ShSpace(soa->vmax,    n), scr.vmax);
    Kokkos::deep_copy(ShSpace(soa->divVmax, n), scr.divVmax);
    Kokkos::deep_copy(ShSpace(soa->maxsoft, n), scr.maxsoft);
    Kokkos::deep_copy(ShU32  (soa->bitflags,n), scr.bitflags);
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    Kokkos::deep_copy(ShSpace(soa->gasmass, n), scr.gasmass);
#endif
#ifdef RT_USE_GRAVTREE
    Kokkos::deep_copy(ShSpace(soa->stellar_lum, (long)n * N_RT_FREQ_BINS), scr.stellar_lum);
#ifdef CHIMES_STELLAR_FLUXES
    Kokkos::deep_copy(ShDbl(soa->chimes_stellar_lum_G0,  (long)n * CHIMES_LOCAL_UV_NBINS), scr.chimes_G0);
    Kokkos::deep_copy(ShDbl(soa->chimes_stellar_lum_ion, (long)n * CHIMES_LOCAL_UV_NBINS), scr.chimes_ion);
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
    Kokkos::deep_copy(ShVec3(soa->rt_source_lum_s,  n), scr.rt_s);
    Kokkos::deep_copy(ShVec3(soa->rt_source_lum_vs, n), scr.rt_vs);
#endif
#ifdef SINK_PHOTONMOMENTUM
    Kokkos::deep_copy(ShSpace(soa->sink_lum,      n), scr.sink_lum);
    Kokkos::deep_copy(ShVec3 (soa->sink_lum_grad, n), scr.sink_lum_grad);
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
    Kokkos::deep_copy(ShSpace(soa->cr_injection, n), scr.cr_inject);
#endif
#ifdef SINK_CALC_DISTANCES
    Kokkos::deep_copy(ShSpace(soa->sink_mass, n), scr.sink_mass);
    Kokkos::deep_copy(ShVec3 (soa->sink_pos,  n), scr.sink_pos);
#if defined(SINK_NODE_MOTION_TRACKED)
    Kokkos::deep_copy(ShInt  (soa->N_SINK,    n), scr.N_SINK);
    Kokkos::deep_copy(ShVec3 (soa->sink_vel,  n), scr.sink_vel);
#endif
#if defined(SPECIAL_POINT_MOTION)
    Kokkos::deep_copy(ShVec3 (soa->sink_acc,  n), scr.sink_acc);
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
    Kokkos::deep_copy(ShSpace(soa->MaxFeedbackVel, n), scr.max_fbvel);
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    Kokkos::deep_copy(ShSpace(soa->tidal_tensorps, (long)n * 6), scr.tidal);
#endif
#ifdef DM_SCALARFIELD_SCREENING
    Kokkos::deep_copy(ShSpace(soa->mass_dm, n), scr.mass_dm);
    Kokkos::deep_copy(ShVec3 (soa->s_dm,    n), scr.s_dm);
    Kokkos::deep_copy(ShVec3 (soa->vs_dm,   n), scr.vs_dm);
#endif

    Kokkos::fence();

    /* ---------------- AoS write-back ---------------------------------- */
    gpu_moment_writeback_to_aos(n);

    /* ---------------- cleanup ----------------------------------------- */
    /* The source-input + Father-mirror pools are persistent (grow-and-keep);
     * they are reused next refresh and freed when the tree is freed
     * (gpu_moment_refresh_release).  No per-call free here. */
    (void) fmirror;

    return 0;
}

/* =================================================================== */
/* Free the persistent gpu_moment_refresh scratch pools.  Reached through */
/* gpu_gravity_tree_release() (alongside the SoA + force-drift pools),    */
/* which force_treefree() calls, so these pools are dropped and regrown   */
/* once per tree epoch rather than held for the whole run.                */
/* =================================================================== */
extern "C" void gpu_moment_refresh_release(void)
{
    precompute_free_(pre_persist_);
    pre_cap_ = 0;
    if(fmirror_persist_) {Kokkos::kokkos_free<GIZMO_KOKKOS_SHARED_SPACE>(fmirror_persist_); fmirror_persist_ = NULL;}
    fmirror_cap_ = 0;
    mr_rawpool_free_();
    raw_cap_ = 0;
}

/* =================================================================== */
/* Bulk SoA → Nodes[]/Extnodes[] AoS write-back.  Runs on host         */
/* (host-side mutation of the AoS arrays so the CPU pseudo path picks  */
/* up the new moments).  SharedSpace pages may need a Kokkos::fence    */
/* before this; callers are responsible.                                */
/* =================================================================== */
extern "C" void gpu_moment_writeback_to_aos(int n)
{
    struct gpu_gravity_tree_soa_t *soa = gpu_gravity_tree_soa();
    if(!soa) {return;}
    int tree_base = All.TreeNodeIndexBase;

    for(int k = 0; k < n; k++) {
        int no = tree_base + k;
        /* The host stamped every node current before this refresh (forcetree.cc, the reset
         * that mirrors the build's moment pass), and the kernel above just recomputed vmax.
         * The device walk widens a node's stored length from the time that length describes,
         * using vmax; so the length and its time go to the mirror together with the new
         * moments, or the mirror is left pairing a fresh vmax and an old time with whatever
         * length it last saw, and a walk on it can open too little. */
        soa->len[k] = Nodes[no].len;
        if(soa->node_ti) {soa->node_ti[k] = Nodes[no].Ti_current;}
        Nodes[no].u.d.mass     = (MyFloat) soa->mass[k];
        Nodes[no].u.d.s        = Vec3<MyFloat>{(MyFloat) soa->s[k][0],
                                                (MyFloat) soa->s[k][1],
                                                (MyFloat) soa->s[k][2]};
        Nodes[no].u.d.bitflags = soa->bitflags[k];
        Nodes[no].N_part       = soa->N_part[k];
        Nodes[no].maxsoft      = (MyFloat) soa->maxsoft[k];
        /* GravCost / Ti_current: matches CPU step-1 reset, but we only
         * touch them in the same way force_refresh_node_moments did
         * before our dispatch was wired in. They are reset on the host
         * before this routine is called (see forcetree.cc step 1). */
        Extnodes[no].vs       = Vec3<MyFloat>{(MyFloat) soa->node_vs[k][0],
                                               (MyFloat) soa->node_vs[k][1],
                                               (MyFloat) soa->node_vs[k][2]};
        Extnodes[no].hmax     = (MyFloat) soa->hmax[k];
        Extnodes[no].vmax     = (MyFloat) soa->vmax[k];
        Extnodes[no].divVmax  = (MyFloat) soa->divVmax[k];
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
        Nodes[no].gasmass     = (MyFloat) soa->gasmass[k];
#endif
#ifdef RT_USE_GRAVTREE
        for(int b = 0; b < N_RT_FREQ_BINS; b++) {
            Nodes[no].stellar_lum[b] = (MyFloat) soa->stellar_lum[(long)k * N_RT_FREQ_BINS + b];
        }
#ifdef CHIMES_STELLAR_FLUXES
        for(int b = 0; b < CHIMES_LOCAL_UV_NBINS; b++) {
            Nodes[no].chimes_stellar_lum_G0 [b] = soa->chimes_stellar_lum_G0 [(long)k * CHIMES_LOCAL_UV_NBINS + b];
            Nodes[no].chimes_stellar_lum_ion[b] = soa->chimes_stellar_lum_ion[(long)k * CHIMES_LOCAL_UV_NBINS + b];
        }
#endif
#endif
#ifdef RT_SEPARATELY_TRACK_LUMPOS
        Nodes[no].rt_source_lum_s   = Vec3<MyFloat>{(MyFloat) soa->rt_source_lum_s [k][0],
                                                     (MyFloat) soa->rt_source_lum_s [k][1],
                                                     (MyFloat) soa->rt_source_lum_s [k][2]};
        Extnodes[no].rt_source_lum_vs = Vec3<MyFloat>{(MyFloat) soa->rt_source_lum_vs[k][0],
                                                       (MyFloat) soa->rt_source_lum_vs[k][1],
                                                       (MyFloat) soa->rt_source_lum_vs[k][2]};
#endif
#ifdef SINK_PHOTONMOMENTUM
        Nodes[no].sink_lum      = (MyFloat) soa->sink_lum[k];
        Nodes[no].sink_lum_grad = Vec3<MyFloat>{(MyFloat) soa->sink_lum_grad[k][0],
                                                 (MyFloat) soa->sink_lum_grad[k][1],
                                                 (MyFloat) soa->sink_lum_grad[k][2]};
#endif
#ifdef COSMIC_RAY_SUBGRID_LEBRON
        Nodes[no].cr_injection = (MyFloat) soa->cr_injection[k];
#endif
#ifdef SINK_CALC_DISTANCES
        Nodes[no].sink_mass = (MyFloat) soa->sink_mass[k];
        Nodes[no].sink_pos  = Vec3<MyFloat>{(MyFloat) soa->sink_pos[k][0],
                                             (MyFloat) soa->sink_pos[k][1],
                                             (MyFloat) soa->sink_pos[k][2]};
#if defined(SINK_NODE_MOTION_TRACKED)
        Nodes[no].sink_vel  = Vec3<MyFloat>{(MyFloat) soa->sink_vel[k][0],
                                             (MyFloat) soa->sink_vel[k][1],
                                             (MyFloat) soa->sink_vel[k][2]};
#endif
#if defined(SPECIAL_POINT_MOTION)
        Nodes[no].sink_acc  = Vec3<MyFloat>{(MyFloat) soa->sink_acc[k][0],
                                             (MyFloat) soa->sink_acc[k][1],
                                             (MyFloat) soa->sink_acc[k][2]};
#endif
#if defined(SINK_NODE_MOTION_TRACKED)
        Nodes[no].N_SINK    = soa->N_SINK[k];
#endif
#if defined(SINGLE_STAR_TIMESTEPPING) && defined(SINGLE_STAR_FB_TIMESTEPLIMIT)
        Nodes[no].MaxFeedbackVel = (MyFloat) soa->MaxFeedbackVel[k];
#endif
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        for(int kk = 0; kk < 6; kk++) {
            Nodes[no].tidal_tensorps_prevstep.data[kk] = (MyFloat) soa->tidal_tensorps[(long)k * 6 + kk];
        }
#endif
#ifdef DM_SCALARFIELD_SCREENING
        Nodes[no].mass_dm = (MyFloat) soa->mass_dm[k];
        Nodes[no].s_dm    = Vec3<MyFloat>{(MyFloat) soa->s_dm[k][0],
                                            (MyFloat) soa->s_dm[k][1],
                                            (MyFloat) soa->s_dm[k][2]};
        Extnodes[no].vs_dm = Vec3<MyFloat>{(MyFloat) soa->vs_dm[k][0],
                                            (MyFloat) soa->vs_dm[k][1],
                                            (MyFloat) soa->vs_dm[k][2]};
#endif
    }
}


