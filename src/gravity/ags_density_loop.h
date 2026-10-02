/* gravity/ags_density_loop.h — AGS adaptive-softening kernel-radius / density
 * neighbor loop module.
 *
 * Defines AgsDensitySpec for the iterative neighbor-loop runner
 * (mesh/neighbor_loop_runner.h). Uses the `IterControl = Iterative` +
 * `uses_ghost_writeback = true` path, extended with the partition-by-
 * subgroup contract + rebuild-every-iter correctness fallback.
 *
 * Replaces the now-deleted gravity/ags_density_gpu.cc (retired after
 * two-binary parity validation passed).
 *
 * Per-pair physics: ports the body at ags_density_gpu.cc:118-205 verbatim,
 * with NlrCommonScalars / per-active host-fill routing instead of direct
 * All.* reads (TRAP 1) and the wakeup write going through the
 * hydro-convention `atomic_max(TimeBin+1)` instead of the buggy
 * `atomic_store(-1)` sentinel.
 *
 * Three hard gates: partition-key assertion (debug-gated, validates `bm
 * key drift` is impossible across iters); CSR fallback rule
 * (`mode_a_rebuild_csr_every_iter` defaults false; flip to true only if a
 * real PA-A divergence appears); and cross-rank positive-wakeup
 * propagation (validated by the two-binary parity route, since retired).
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef AGS_DENSITY_LOOP_H
#define AGS_DENSITY_LOOP_H

/* Kokkos_Core must precede allvars.h (its macros may conflict with stdlib
 * names). */
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE

#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"   /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
/* NOTE: caller TUs must include "../mesh/kernel.h" before this header.
 * kernel.h has no include guards; defines kernel_main / kernel_gravity /
 * kernel_hinv used by the inline pair body below. */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* Forward decls — defined caller-side in gravity/ags_rkern.cc. */
int    ags_density_isactive(int i);
int    ags_gravity_kernel_shared_BITFLAG(short int particle_type_primary);
double ags_return_minsoft(int i);
double ags_return_maxsoft(int i);

/* ============================================================================
 * Per-pair physics types (file scope; PascalCase). AgsDensitySpec re-exports
 * via `using CallScalars = ...; using ActiveData = ...; using AccumData = ...`.
 * ========================================================================== */

/* Per-call cosmology + AGS scalars captured once from All.* on the host by
 * populate_call_scalars. Passed to the inline pair body as its own argument,
 * once per call, so the body never reads All.* directly (TRAP 1). TimeBinActive
 * is captured by value (TIMEBINS=60 ints, ~240 bytes) — which is affordable
 * precisely because this is per CALL, not per active particle. */
struct AgsDensityCallScalars {
    NlrCommonScalars common;
    int              TimeBinActive[TIMEBINS];
    double           Time;
    double           TimeBegin;
    double           AGS_DesNumNgb;
    double           AGS_MaxNumNgbDeviation;
    double           fac_mu;                 /* -3.0 / All.cf_atime — vsig prefactor */
};

/* (Per-active state previously held in AgsDensityLocalIn + ctx.per_active_local
 * UVM is now read directly from ctx.P[i] in load_active. The old UVM array
 * was indexed by subgroup-local active_slot, which collides across bm
 * subgroups. ctx.P is path-correct on both Mode A (arena P_gpu) and Mode B
 * (request-driven slab), and the fields we read (Pos, Vel, TimeBin) are
 * read-only.) */

/* AccumData. Per-iter i-side accumulators. Matches legacy ags_density_gpu.cc
 * out struct (Ngb, DrkernNgb, AGS_zeta, Particle_DivVel, AGS_vsig) plus
 * AGS_FACE NV_T upper-triangle. AGS_vsig is MAX-merged across peers; the
 * others are ADD-merged. NV_T lower triangle stays zero (mirrors legacy
 * memset behavior).
 *
 * Zeroed per outer iter (runner provides zero_accum); each bm subgroup
 * contributes into the same per-active accumulator (each active lives in
 * exactly one bm subgroup — `actives_partition_by_subgroup = true`). */
struct AgsDensityAccumData {
    double Ngb;
    double DrkernNgb;
    double AGS_zeta;
    double Particle_DivVel;
    double AGS_vsig_max;
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    /* Upper triangle of the raw 3x3 NV_T moment accumulator. after_iter
     * symmetrises + inverts these into the full face operator P[i].NV_T via
     * ags_invert_nvt_for_faces (the lower triangle is implied symmetric). */
    double NV_T_00;
    double NV_T_01;
    double NV_T_02;
    double NV_T_11;
    double NV_T_12;
    double NV_T_22;
#endif
};

/* IterScratch — per-active state living across iterations.
 *
 * AGS_KernelRadius is NOT here — the runner tracks it via radii_uvm and
 * exposes it via AfterIterContext::h_search_current; mutation is via the
 * IterResult { AdjustRadius, new_h_search } return from after_iter.
 *
 * `converged` is NOT here either — the runner removes converged actives
 * from the active set when after_iter returns IterStatus::Converged.
 *
 * What stays here: legacy bisection state (Left, Right) + the entry-time
 * radius snapshot (AGS_Prev, never updated) + the legacy max/min-rkern
 * marker flags (used by after_iter to drive the redo logic). */
struct AgsDensityIterScratch {
    double Left;
    double Right;
    double AGS_Prev;
    int    set_to_maxrkern;
    int    set_to_minrkern;
    /* Monotonic per-call max condition number of the AGS face NV_T inversion.
     * Zero-memset at iter 0 by the runner; a high value at any iteration
     * ratchets the desired neighbor count up (never down) so the search
     * expands to improve conditioning without oscillating. AGS_FACE only. */
    double condition_number_max;
};

/* DeviceContext extension. Single sticky call-scope flag for "any wakeup
 * written across all iters of this call".
 * Trivially copyable: the runner captures by value into Kokkos device
 * lambdas.
 *
 * need_wakeup_uvm is allocated + zeroed ONCE in populate_device_context
 * and propagated to NeedToWakeupParticles_local ONCE in
 * cleanup_device_context (at end of the whole iterative call). It is NOT
 * reset per iter — a per-iter reset would lose non-final-iter wakeups
 * (e.g. iter-1 sets need_wakeup, iter-N doesn't,
 * cleanup sees 0). Sticky matches the legacy GPU evaluator's per-CALL
 * d_need_wakeup semantic (legacy ags_density() ran one evaluator-call
 * per outer iter so the distinction didn't surface, but the runner does
 * the iter loop inside one call so we must accumulate). */
struct AgsDensityDeviceContext : NeighborLoopDeviceContextBase {
    int                     *need_wakeup_uvm;    /* UVM, single int.
                                                  * Sticky across all iters
                                                  * of one runner call. */
    unsigned char           *wakeup_dirty_base;  /* WakeupDirty sidecar base (global UVM); populate sets from WakeupDirty */
};

/* Active-particle state passed into the pair body. `pos` and `h_search`
 * top-level so the runner's Mode B walker can read them directly;
 * `Vel` and `TimeBin` are the only per-active fields the kernel needs
 * (Vel for the dv computation, TimeBin for the hydro-convention wakeup
 * write). All read from ctx.P[i] in load_active — no per_active_local
 * UVM staging needed (was prone to subgroup-slot collision under
 * multi-bm partitioning). */
struct AgsDensityActiveState {
    Vec3<double>           pos;
    double                 h_search;
    Vec3<double>           Vel;
    short int              TimeBin;
    int                    origin_local_idx;
    int                    origin_rank;
};

/* ============================================================================
 * Inline pair body — single source of truth for ags_density per-pair physics.
 * Called from AgsDensitySpec::pair_kernel; used by Mode A GPU lambda, Mode B
 * local host walker, and Mode B remote peer-to-peer.
 *
 * Verbatim port of ags_density_gpu.cc:118-205 with:
 *
 *   1. NlrCommonScalars / per-active host-fill replacing direct All.* /
 *      P[i].* reads (TRAP 1). All gates preserved: ADAPTIVE_GRAVSOFT_FORALL,
 *      DM_FUZZY, CBE_INTEGRATOR (gate the wakeup test); GALSF (suppress
 *      wakeup for type-4 cosmo, type-2/3 non-cosmo); AGS_FACE_CALCULATION_IS_ACTIVE
 *      (NV_T accumulation).
 *
 *   2. Wakeup write uses hydro-convention atomic_max(local.TimeBin + 1)
 *      instead of legacy atomic_store(-1); the legacy `-1` sentinel would
 *      silently drop under MAX reverse-comm.
 *
 *   3. (removed)
 *      pattern).
 *
 *   4. *need_wakeup is atomic_or-d to 1 on any wakeup write (local or
 *      ghost-side). Caller copies to NeedToWakeupParticles_local after
 *      runner returns. Replaces legacy ags_density_gpu.cc:210
 *      `if(*d_need_wakeup) NeedToWakeupParticles_local = 1`.
 *
 *   5. Mass / kernel-radius guards (legacy lines 121-127) move to the
 *      caller's active-set filter (ags_density_isactive); we still do the
 *      Mass <= 0 guard here as a defensive check (cheap, matches legacy).
 * ========================================================================== */
KOKKOS_INLINE_FUNCTION
static void ags_density_pair_kernel_body(const AgsDensityActiveState& active,
                                          const AgsDensityCallScalars& scalars,
                                          struct particle_data&        neighbor_particle,
                                          struct gas_cell_data*        neighbor_cell,
                                          AgsDensityAccumData&         accum,
                                          int*                         need_wakeup,
                                          unsigned char*               wakeup_dirty_slot)
{
    /* neighbor_cell is read for gas neighbors (CellP[j].VelPred); never
     * written. nullptr for non-gas neighbors and skipped via the
     * neighbor_particle.Type == 0 gate below. */

    if(neighbor_particle.Mass <= 0) return;

    const double h_i  = active.h_search;          /* runner-provided per-iter radius */
    const double h_i2 = h_i * h_i;

    /* dpos with periodic wrap. Sign convention matches legacy
     * (pos_i - kp[j].Pos => dp from j to i). */
    Vec3<double> dp;
    dp[0] = active.pos[0] - (double)neighbor_particle.Pos[0];
    dp[1] = active.pos[1] - (double)neighbor_particle.Pos[1];
    dp[2] = active.pos[2] - (double)neighbor_particle.Pos[2];
    nearest_xyz(dp, -1);
    const double r2 = dp.norm_sq();
    if(!(r2 < h_i2)) return;

    double hinv, hinv3, hinv4;
    kernel_hinv(h_i, &hinv, &hinv3, &hinv4);

    const double r = sqrt(r2);
    const double u = r * hinv;
    double wk, dwk;
    kernel_main(u, hinv3, hinv4, &wk, &dwk, 0);

    accum.Ngb       += wk;
    accum.DrkernNgb += -(NUMDIMS * hinv * wk + u * dwk);
    accum.AGS_zeta  += (double)neighbor_particle.Mass * kernel_gravity(u, hinv, hinv3, 0);

    if(r > 0) {
        Vec3<double> dv;
        /* j-side velocity: gas uses CellP[j].VelPred; non-gas uses P[j].Vel.
         * Matches legacy ags_density_gpu.cc:156-157. */
        if(neighbor_particle.Type == 0 && neighbor_cell) {
            dv[0] = active.Vel[0] - (double)neighbor_cell->VelPred[0];
            dv[1] = active.Vel[1] - (double)neighbor_cell->VelPred[1];
            dv[2] = active.Vel[2] - (double)neighbor_cell->VelPred[2];
        } else {
            dv[0] = active.Vel[0] - (double)neighbor_particle.Vel[0];
            dv[1] = active.Vel[1] - (double)neighbor_particle.Vel[1];
            dv[2] = active.Vel[2] - (double)neighbor_particle.Vel[2];
        }
        NGB_SHEARBOX_BOUNDARY_VELCORR_(active.pos, neighbor_particle.Pos, dv, 1);

        double v_dot_r = dot(dp, dv);
        if(v_dot_r > 0) v_dot_r *= 0.333333;
        const double vsig = 0.5 * fabs(scalars.fac_mu * v_dot_r / r);

        short int TimeBin_j = neighbor_particle.TimeBin;
        if(TimeBin_j < 0) TimeBin_j = -TimeBin_j - 1;
        if(vsig > accum.AGS_vsig_max) accum.AGS_vsig_max = vsig;

#if defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
        /* Wakeup test (legacy ags_density_gpu.cc:165-179). Gated on the AGS
         * kernel-radius SSOT so SIDM/fuzzy/CBE (which all use AGS_vsig for
         * their timestep) are covered, not just ADAPTIVE_GRAVSOFT_FORALL.
         *
         * volatile retained from legacy: nvc++ folds the predicate to the
         * initial value otherwise (per legacy comment). */
        volatile int wakeup_condition = 0;
        if(!(scalars.TimeBinActive[TimeBin_j]) &&
           (scalars.Time > scalars.TimeBegin) &&
           (vsig > WAKEUP * (double)neighbor_particle.AGS_vsig)) {
            wakeup_condition = 1;
        }
#if defined(GALSF)
        /* Suppress wakeup for stars / non-cosmo middle types
         * (legacy ags_density_gpu.cc:170-174). */
        if(is_galsf_stellar_candidate_type(neighbor_particle.Type, scalars.common.comoving_integration_on)) {
            wakeup_condition = 0;
        }
#endif
        if(wakeup_condition) {
            /* Hydro-convention wakeup write — `active.TimeBin + 1` (positive).
             * Cross-rank propagation works through ghost_writeback's MAX
             * reverse-comm. Legacy wrote (short int)-1, which the MAX path
             * silently dropped against home=0; fixed before this port
             * landed. */
            const short int wakeup_val = (short int)(active.TimeBin + 1);
            Kokkos::atomic_max(&neighbor_particle.wakeup, wakeup_val);
            if(need_wakeup) { Kokkos::atomic_or(need_wakeup, 1); }
            if(wakeup_dirty_slot) { *wakeup_dirty_slot = 1; }   /* dirty-sidecar mark (byte store; race-benign) */
        }
#endif

        accum.Particle_DivVel -= dwk * dot(dp, dv) / r;

#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
        accum.NV_T_00 += wk * dp[0] * dp[0];
        accum.NV_T_01 += wk * dp[0] * dp[1];
        accum.NV_T_02 += wk * dp[0] * dp[2];
        accum.NV_T_11 += wk * dp[1] * dp[1];
        accum.NV_T_12 += wk * dp[1] * dp[2];
        accum.NV_T_22 += wk * dp[2] * dp[2];
#endif
    }
}

/* ============================================================================
 * AgsDensitySpec — the NeighborLoopSpec contract for ags_density.
 *
 * Iterative + ghost_writeback + multi-subgroup (one bm-group per j-type
 * bitmask). First production user of:
 *   - actives_partition_by_subgroup = true (each active fixed to one bm
 *     subgroup by P[i].Type — see ags_gravity_kernel_shared_BITFLAG);
 *   - active_subgroup_key returning the bm key (= subgroup's j_type_bitmask);
 *   - hydro-convention `atomic_max(TimeBin+1)` wakeup writes via the
 *     `PARTICLE_MAX(wakeup)` ghost-writeback bundle (vs legacy `-1` sentinel).
 * ========================================================================== */
struct AgsDensitySpec {
    /* ====================================================================
     * PHYSICS BLOCK
     * ==================================================================== */

    static constexpr const char *loop_name = "ags_density";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::EpsilonAtomic; /* EpsilonAtomic: only j-writes are gated atomic_max(wakeup)+atomic_or(need_wakeup), both order-invariant, no read-back */

    /* search_mode is ONEWAY: the legacy ags_density CSR builder used
     * NGB_SEARCH_ONEWAY, and the pair predicate is one-way (r < h_i only —
     * no h_j check). Symmetric search would over-include neighbors via
     * h_j logic and inflate NumNgb. AGSForce_calc is symmetric (its pair
     * terms vanish at r >= max(h_i,h_j)) — that loop is a separate port.
     * (Previously wrongly set to MODE_B_SEARCH_SYMMETRIC in the scaffold.)
     *
     * neighbor_type_mask is overridden per-subgroup by the runner from
     * subgroups[sg].j_type_bitmask (multi-bm support — ags_density
     * partitions actives by Type → bm via ags_gravity_kernel_shared_BITFLAG).
     */
    static constexpr int                     search_mode        = MODE_B_SEARCH_ONEWAY;
    static constexpr unsigned int            neighbor_type_mask = 0xFFFFFFFFu;  /* per-subgroup override */
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;

    static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
    /* AGS partitions actives into subgroups with different j_type_bitmask
     * values. A step-persistent SIDX is built for exactly one type mask, so
     * sharing AllTypes/GasOnly across AGS subgroups or later sink calls is
     * invalid. Build a local SIDX per CSR until/unless a per-mask cache lands. */
    static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::None;

    static bool is_active(int particle_index) { return ags_density_isactive(particle_index) != 0; }

    /* Iterative + multi-subgroup traits (step 0 prep-commit additions) */
    using IterControl       = Iterative;
    using IterScratch       = AgsDensityIterScratch;
    using SupportsSubgroups = std::true_type;

    static constexpr int    max_iters                       = MAXITER;
    static constexpr double mode_a_csr_buffer_factor        = 1.3;   /* same structure as density, where 1.3 measured optimal */
    /* CSR correctness fallback. A regression that motivated flipping this
     * to true was traced to an unrelated bug: per-TU All_dev was
     * unsynced, so AGS_DesNumNgb came through as 0 and bisection
     * inverted. Root cause
     * fixed in populate_call_scalars via gizmo_gpu_host_all_ptr(); leaving
     * this at the optimized default (false). Flip to true only if a NEW
     * divergence is observed that cannot be otherwise explained. */
    static constexpr bool   mode_a_rebuild_csr_every_iter   = false;

    /* Partition-by-subgroup contract per design v0.4.2 §5a. Each active i
     * lives in exactly one bm subgroup determined by P[i].Type (fixed for
     * lifetime). Runner asserts via active_subgroup_key under
     * GIZMO_NLR_ASSERT_PARTITION / DEBUG. */
    static constexpr bool actives_partition_by_subgroup     = true;

    using CallScalars   = AgsDensityCallScalars;
    using ActiveData    = AgsDensityActiveState;
    using AccumData     = AgsDensityAccumData;
    using DeviceContext = AgsDensityDeviceContext;

    /* NeighborData carries NON-CONST pointer to particle_data — pair_kernel
     * does atomic_max into neighbor_particle.wakeup. neighbor_cell may be
     * nullptr for non-gas neighbors. */
    struct NeighborData {
        struct particle_data *neighbor_particle;
        struct gas_cell_data *neighbor_cell;     /* nullptr for non-gas / when no CellP */
        int                  *need_wakeup;       /* shared scratch in ctx */
        unsigned char        *wakeup_dirty_slot; /* &WakeupDirty[j] */
    };

    /* Empty Aux — the Spec contract requires the typedef but ags_density
     * has no per-call host-only aux buffers. Active state is read directly
     * from ctx.P[i] in load_active; per-iter convergence + final P[i]
     * scatter happens inside after_iter on host (where P[] is accessible).
     * apply_active_writeback is a no-op. */
    struct Aux { };

    static constexpr bool uses_ghost_write_detector = true;
    static constexpr bool uses_ghost_writeback      = true;

    /* ghost_write_detector_begin/end: runner default (loop_name = "ags_density"). */
    static void ghost_writeback_begin      (const struct neighbor_loop_args&,
                                             const struct NeighborLoopPlan&);
    static void ghost_writeback_end        (const struct neighbor_loop_args&,
                                             const struct NeighborLoopPlan&);

    /* Per-active host hooks. */
    static double      search_radius(const neighbor_loop_args& args,
                                      int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);

    /* DeviceContext lifecycle. populate allocates + zeros need_wakeup_uvm;
     * cleanup frees it AND propagates the sticky-across-iters flag into
     * NeedToWakeupParticles_local. NO per-iter reset — a per-iter reset
     * would lose non-final-iter wakeups. */
    static void populate_device_context(const neighbor_loop_args& args, DeviceContext& ctx);
    static void cleanup_device_context (const neighbor_loop_args& args, DeviceContext& ctx);

    /* Partition-key hook. Returns the bm key for the
     * subgroup this active belongs to.
     *
     * Inlined verbatim from ags_gravity_kernel_shared_BITFLAG
     * (gravity/ags_rkern.cc:45-79) WITHOUT the legacy helper's All.*
     * read: legacy reads All.ComovingIntegrationOn
     * under GALSF; here we route through scalars.common.comoving_integration_on
     * which is host-snapshotted at populate_call_scalars time. Type comes
     * via ctx.particle_type(i) (path-correct: arena P on Mode A,
     * request-driven slab on Mode B). NO globals; the host-side debug
     * partition assertion in runner.cc and any future device-side
     * caller both see consistent values.
     *
     * Decision tree must match legacy bit-for-bit so the runner's per-iter
     * assertion compares against `subgroups[sg].j_type_bitmask` (built by
     * the caller using the same legacy helper at iter-0). Any drift here
     * triggers the runner's "actives_partition_by_subgroup violation"
     * controlled stop. */
    KOKKOS_INLINE_FUNCTION
    static int active_subgroup_key(const DeviceContext& ctx,
                                    int                  i,
                                    const CallScalars&   scalars)
    {
        const short int t = ctx.particle_type(i);

#ifdef ADAPTIVE_GRAVSOFT_FORALL
        if(!((1 << t) & (ADAPTIVE_GRAVSOFT_FORALL))) { return 0; }
#endif

#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        if(!((1 << t) & (ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION))) {
            return ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION;
        }
#endif

        if(t == 0) { return 1; }

#if (ADAPTIVE_GRAVSOFT_FORALL & 32) && defined(SINK_PARTICLES)
        if(t == 5) { return 1; }
#endif

#if defined(GALSF) && ( (ADAPTIVE_GRAVSOFT_FORALL & 16) || (ADAPTIVE_GRAVSOFT_FORALL & 8) || (ADAPTIVE_GRAVSOFT_FORALL & 4) )
        /* GALSF branch reads comoving_integration_on through scalars,
         * not the All global. */
        if(scalars.common.comoving_integration_on) {
            if(t == 4) { return 17; }       /* 2^0 + 2^4 */
        } else {
            if((t == 4) || (t == 2) || (t == 3)) { return 29; }  /* 2^0 + 2^2 + 2^3 + 2^4 */
        }
#endif

#ifdef DM_SIDM
        if((1 << t) & (DM_SIDM)) { return DM_SIDM; }
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
        return (1 << t);
#endif

        return 0;
    }

    /* Build ActiveData[slot] on device. Reads per-active state directly
     * from ctx.P[i] — path-correct (Mode A: arena P_gpu; Mode B: request-
     * driven slab). No per_active_local UVM staging — the legacy UVM was
     * indexed by subgroup-local active_slot, which collides across bm
     * subgroups in the multi-subgroup port.
     *
     * Read fields are: Pos, Vel, TimeBin. All read-only — pair_kernel
     * never writes to P[i], only to neighbor P[j] (wakeup write). */
    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& dctx,
                                   int                  active_slot,
                                   int                  i,
                                   double               h_search,
                                   const CallScalars&   scalars)
    {
        ActiveData active;
        active.pos[0]           = (double)dctx.P[i].Pos[0];
        active.pos[1]           = (double)dctx.P[i].Pos[1];
        active.pos[2]           = (double)dctx.P[i].Pos[2];
        active.Vel[0]           = (double)dctx.P[i].Vel[0];
        active.Vel[1]           = (double)dctx.P[i].Vel[1];
        active.Vel[2]           = (double)dctx.P[i].Vel[2];
        active.TimeBin          = dctx.P[i].TimeBin;
        active.h_search         = h_search;
        (void)scalars;
        active.origin_local_idx = active_slot;
        active.origin_rank      = -1;
        return active;
    }

    /* Build NeighborData. ctx.P / ctx.CellP are path-correct (Mode A arena
     * or Mode B request-driven slab). need_wakeup comes
     * from the ctx. */
    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& dctx,
                                       int                  j,
                                       const IdentitySidecar& /*id*/,
                                       const ActiveData&    /*active*/)
    {
        NeighborData neighbor;
        neighbor.neighbor_particle = &dctx.P[j];
        neighbor.neighbor_cell     = (dctx.CellP && dctx.P[j].Type == 0) ? &dctx.CellP[j] : nullptr;
        neighbor.need_wakeup       = dctx.need_wakeup_uvm;
        neighbor.wakeup_dirty_slot = dctx.wakeup_dirty_base ? &dctx.wakeup_dirty_base[j] : nullptr;
        return neighbor;
    }

    /* Zero accumulator. Byte-zero of the POD AccumData; AGS_vsig_max starts
     * at 0 (legacy initializes to 0 and uses if(vsig > AGS_vsig) max-merge). */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum)
    {
        for(size_t b = 0; b < sizeof(accum); b++) ((char*)&accum)[b] = 0;
    }

    /* The physics — forward to the inline pair body. */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData&   active,
                             const NeighborData& neighbor,
                             AccumData&          accum,
                             NoScatter&          /*scatter*/,
                            const CallScalars& cs)
    {
        ags_density_pair_kernel_body(active, cs,
                                      *neighbor.neighbor_particle,
                                      neighbor.neighbor_cell,
                                      accum,
                                      neighbor.need_wakeup,
                                      neighbor.wakeup_dirty_slot);
    }

    /* NO-OP for AgsDensitySpec. after_iter writes post-processed per-iter
     * values directly to P[i] on host (legacy semantics — see
     * ags_density_loop.cc::after_iter). No aux buffer; caller's post-runner
     * finalize pass reads P[i] directly. Declared because the runner
     * contract requires the hook to exist. */
    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int                       active_slot,
                                        int                       i,
                                        const AccumData&          accum);

    /* Per-field merge for Mode B remote peer accumulators. */
    /* ============================================================================
     * MERGE_ACCUM — Mode B remote peer accumulator merge.
     *
     * Manifest pattern matching sink_feed_loop.cc::merge_accum. Per-field op
     * MUST match what pair_kernel writes (field parity was validated by the
     * retired two-binary route).
     *
     * AGS_vsig_max is MAX-merged (legacy ASSIGN_MAX in scatter); everything
     * else is ADD-merged. AGS_FACE NV_T fields are ADD-merged (each peer
     * contributes neighbor wk * dpos pairs that sum cleanly).
     * ========================================================================== */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local_accum, const AccumData& peer_accum)
    {
#define ACCUM_ADD(field) local_accum.field += peer_accum.field;
#define ACCUM_MAX(field) if(peer_accum.field > local_accum.field) local_accum.field = peer_accum.field;

        ACCUM_ADD(Ngb)
        ACCUM_ADD(DrkernNgb)
        ACCUM_ADD(AGS_zeta)
        ACCUM_ADD(Particle_DivVel)
        ACCUM_MAX(AGS_vsig_max)
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
        ACCUM_ADD(NV_T_00)
        ACCUM_ADD(NV_T_01)
        ACCUM_ADD(NV_T_02)
        ACCUM_ADD(NV_T_11)
        ACCUM_ADD(NV_T_12)
        ACCUM_ADD(NV_T_22)
#endif

#undef ACCUM_ADD
#undef ACCUM_MAX
    }

    /* Iterative hooks. */
    static IterResult after_iter      (const AfterIterContext<AgsDensitySpec>& ctx,
                                        const AccumData&                       accum);
    static void       after_iter_global(const neighbor_loop_args&              args,
                                         const struct NlrIterDriver<AgsDensitySpec>& drv);

    /* ====================================================================
     * ENGINE APPARATUS
     * ==================================================================== */
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;

};

#endif /* AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE */

#endif /* AGS_DENSITY_LOOP_H */
