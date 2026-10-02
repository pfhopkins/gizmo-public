/* rt_source_injection_loop.h — RtSrcInjectionSpec for the runner-template port
 * of rt_source_injection.
 *
 * Non-iterative scatter loop: active non-gas sources with luminosity and
 * positive `KernelSum_Around_RT_Source` deposit radiation energy / flux /
 * intensities (set varies with RT_* feature flags) into surrounding gas
 * neighbors. Pattern parallels ThermalFBSpec (non-iterative + ghost_writeback
 * bundle for j-side writes), with two distinctions:
 *
 *   (1) Toplevel-built active list. The host orchestration TU
 *       (radiation/rt_source_injection.cc) pre-builds the active list with
 *       per-source host pre-fill (which also clears
 *       P[i].Sink_accreted_photon_energy under RT_REINJECT_ACCRETED_PHOTONS)
 *       and hands the runner a ready-made (active_src, RtSrcLocalIn[]) pair
 *       via Aux. Spec::is_active is a defensive Type re-check only;
 *       nlr_build_active_list is NOT used (it would re-include sources for
 *       which the pre-fill has already mutated host state).
 *
 *   (2) Search mode = MODE_B_SEARCH_SYMMETRIC matches the legacy GPU
 *       evaluator's NGB_SEARCH_SYMMETRIC (rt_source_injection_gpu.cc:178).
 *       This is correctness-required under RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION
 *       (the kernel needs h_j-sided neighbors when All.TimeStep>0); benign
 *       overhead in the other config (the body rejects r2 >= h_i^2 either way).
 *
 * SSOT for rt_source_injection physics types and per-pair kernel — supersedes
 * the pre-port `rt_source_injection_functions.h` and the legacy
 * `rt_source_injection_gpu.{cc,h}`, retired in the cleanup commit.
 *
 * RtSrcLocalIn lives in `rt_source_injection_local_in.h` so the host TU can
 * include it without pulling Kokkos / runner-template headers.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifndef RT_SOURCE_INJECTION_LOOP_H
#define RT_SOURCE_INJECTION_LOOP_H

/* Kokkos_Core must precede allvars.h (its macros may conflict with stdlib names). */
#include <Kokkos_Core.hpp>

#include "../declarations/gpu_all_mirror.h"   /* per-TU AllDeviceMirror; safe in host pass too */
#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"

#ifdef RT_SOURCE_INJECTION

/* NOTE: caller TUs must `#include "../mesh/kernel.h"` BEFORE this header.
 * kernel.h has no include guards — pulling it here would multiply-define
 * static helpers in any caller that already had kernel.h on its include list
 * (runner / _loop.cc / .cc all do). Same pattern as thermal_fb_loop.h. */
#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"      /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
#include "../sinks/sink_functions.h"           /* sink_fb_angleweight_localcoupling_gpu */
#include "rt_functions.h"                      /* rt_get_donation_target_bin, RT_BAND_IS_IONIZING */
#include "rt_source_injection_local_in.h"      /* RtSrcLocalIn (host-clean) */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* ============================================================================
 * Pade-coefficient slab-averaged optical-depth factor — relocated verbatim from
 * the retired rt_source_injection_functions.h. Static-inline + KOKKOS_INLINE
 * keeps it private to this header's instantiations (no duplicate symbol with
 * any orphan _functions.h still on disk during commit 1).
 * ========================================================================== */
KOKKOS_INLINE_FUNCTION
static double rt_slab_avg_loop(double x)
{
    return (1.00000000000 + x*(0.21772719088733913 + x*(0.047076512011644776 + x*0.005068307557496351))) /
           (1.00000000000 + x*(0.71772719088733920 + x*(0.239273440788647680 + x*(0.046750496137263675 + x*0.005068307557496351))));
}

/* ============================================================================
 * Per-call scalars. Minimal: the pair kernel reads All.* directly via the
 * per-TU AllDeviceMirror (#define All AllDeviceMirror is in scope on the
 * device pass; benign on the host pass), so no UNIT_* / All.* pre-aggregation
 * needed here. NlrCommonScalars carries the runner-required cf_*, TimeStep,
 * NumCurrentTiStep snapshot.
 * ========================================================================== */
struct RtSrcInjCallScalars {
    NlrCommonScalars common;
};

/* ============================================================================
 * AccumData — EMPTY for this Spec, and that is correct.
 *
 * This loop accumulates nothing per active. All of its output is the
 * atomic_adds into Pj / Cj inside the pair kernel (gas-side: Rad_Je /
 * Rad_E_gamma / Rad_E_gamma_Pred / Rad_Intensity[/Pred] / Rad_Flux[/Pred] /
 * VelPred; particle-side: P[j].Vel / P[j].dp), carried across ranks by the
 * ghost-writeback bundle. The type remains because the Spec contract requires
 * the alias and the runner's per-active reduce is written in terms of it.
 *
 * Order-independence: every j-side write this loop performs is a fresh
 * atomic_add into fields this loop neither reads nor accumulates from within
 * the pair body, so the result does not depend on source order. Contrast with
 * thermal_fb, which reads Pj.Mass and atomic_adds to it in the same loop.
 * ========================================================================== */
struct RtSrcInjAccum {
};

/* Runner ActiveData. `pos` and `h_search` are top-level (runner's Mode B
 * walker reads them directly). */
struct RtSrcInjActiveState {
    Vec3<double>           pos;
    double                 h_search;
    RtSrcLocalIn           local;
};

/* DeviceContext extension: UVM-resident per-active host-fill array.
 * Trivially copyable; runner captures by value into device lambdas. */
struct RtSrcInjDeviceContext : NeighborLoopDeviceContextBase {
    const RtSrcLocalIn *per_active_local;   /* UVM, [num_active]; nullptr when num_active==0 */
};

/* ============================================================================
 * Inline pair body — SSOT for rt_source_injection per-pair physics.
 *
 * Caller pre-filters Pj.Type==0 / Pj.Mass>0 via the gas-only neighbor mask in
 * load_neighbor (defensive checks still inside as safety). dp = source.Pos -
 * Pj.Pos with nearest_xyz applied. r2 pre-checked by the wrapper.
 *
 * RELOCATED from the legacy rt_source_injection_functions.h verbatim, MINUS
 * the GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION branch — which used
 * atomic_exchange (assignment) semantics that are incompatible with the
 * additive ghost-writeback bundle. The testproblem-specific boundary
 * condition is now imposed by an owner-side postpass in
 * radiation/rt_source_injection.cc, called after the runner returns; see
 * rt_source_injection_apply_grain_rdi_boundary_fixup() there.
 * ========================================================================== */
KOKKOS_INLINE_FUNCTION
static void rt_source_injection_pair_body(
    const RtSrcLocalIn& local,
    struct particle_data& Pj,
    struct gas_cell_data& Cj,
    double r2,
    const Vec3<double>& dp)
{
    double r = sqrt(r2);
    double hinv, hinv3, hinv4;
    kernel_hinv(local.KernelRadius, &hinv, &hinv3, &hinv4);

    double wk;
#ifdef RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION
    if (All.TimeStep > 0) {
        wk = sink_fb_angleweight_localcoupling_gpu(Pj, Cj, 0., r, local.KernelRadius)
             / local.KernelSum_Around_RT_Source;
    } else {
        wk = (1. - r2 * hinv * hinv) / local.KernelSum_Around_RT_Source;
    }
#else
    wk = (1. - r2 * hinv * hinv) / local.KernelSum_Around_RT_Source;
#endif
    if (!(wk >= 0)) { wk = 0; } /* kernel weight must never be negative */

#ifdef RT_EVOLVE_INTENSITIES
    int kx;
    double angle_wt_Inu_sum = 0, angle_wt_Inu[N_RT_INTENSITY_BINS];
    for (kx = 0; kx < N_RT_INTENSITY_BINS; kx++) {
        double cos_t = 0;
        for (int kq = 0; kq < 3; kq++) {
            cos_t += All.Rad_Intensity_Direction[kx][kq] * dp[kq] / r;
        }
        double wt = cos_t * cos_t * cos_t * cos_t;
        if (cos_t < 0) wt = 0;
        angle_wt_Inu[kx] = wt;
        angle_wt_Inu_sum += wt;
    }
#endif

    for (int k = 0; k < N_RT_FREQ_BINS; k++) {
        double dE = wk * local.Luminosity[k];
        Vec3<double> dfluxes = {};

#if !defined(RT_INJECT_PHOTONS_DISCRETELY)
        Kokkos::atomic_add(&Cj.Rad_Je[k], dE);
#endif

#if defined(RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION) || defined(RT_REPROCESS_INJECTED_PHOTONS)
        double x_abs = 2. * Cj.Rad_Kappa[k] * (Cj.Density * All.cf_a3inv)
                       * (DMAX(2. * Pj.Get_Particle_Size(), DMAX(local.KernelRadius, r))) * All.cf_atime;
        double slabfac_x = x_abs * rt_slab_avg_loop(x_abs);
        if (slabfac_x != slabfac_x || slabfac_x <= 0) { slabfac_x = 0; } else if (slabfac_x > 1) { slabfac_x = 1; }
#if !defined(RT_DISABLE_RAD_PRESSURE) && defined(RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION)
        double dv = -slabfac_x * dE / (C_LIGHT_CODE_REDUCED * Pj.Mass);
        for (int kv = 0; kv < 3; kv++) {
            double dv_tmp = dv * (dp[kv] / r) * All.cf_atime;
            Kokkos::atomic_add(&Pj.Vel[kv], dv_tmp);
            Kokkos::atomic_add(&Cj.VelPred[kv], dv_tmp);
            Kokkos::atomic_add(&Pj.dp[kv], dv_tmp * Pj.Mass);
        }
#endif
#endif

#ifdef RT_REPROCESS_INJECTED_PHOTONS
        double dE_donation = 0;
        int donation_bin = rt_get_donation_target_bin(k), do_donation = (donation_bin > -1) ? 1 : 0;
#ifdef RT_CHEM_PHOTOION
        if (RT_BAND_IS_IONIZING(k)) {
            double stellum = 0;
            if (local.Dt > 0) {
#if (RT_CHEM_PHOTOION == 1)
                stellum = local.Luminosity[RT_FREQ_BIN_H0];
#else
                for (int k2 = RT_FREQ_BIN_H0; k2 < RT_FREQ_BIN_H0 + 4; k2++) { stellum += local.Luminosity[k2]; }
#endif
                stellum *= 1. / (C_LIGHT_CODE_REDUCED / C_LIGHT_CODE) / local.Dt * UNIT_LUM_IN_CGS;
            }
            double RHII = 4.01e-9 * pow(stellum, 0.333)
                          * pow(local.Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS, -0.66667)
                          / UNIT_LENGTH_IN_CGS;
            if (DMAX(r, Pj.Get_Particle_Size()) * All.cf_atime < RHII) { do_donation = 0; }
        }
#endif
#ifdef RT_INFRARED
        if (k == RT_FREQ_BIN_INFRARED) { do_donation = 0; }
#endif
        if (do_donation) { dE_donation = slabfac_x * dE; dE *= fabs(1. - slabfac_x); }
#endif /* RT_REPROCESS_INJECTED_PHOTONS */

#if defined(RT_EVOLVE_FLUX) && defined(RT_INJECT_PHOTONS_DISCRETELY_ADD_MOMENTUM_FOR_LOCAL_EXTINCTION)
        {
            double dflux = -dE * C_LIGHT_CODE_REDUCED / r;
            dfluxes += dflux * dp;
        }
#endif

#if defined(RT_INJECT_PHOTONS_DISCRETELY)
        Kokkos::atomic_add(&Cj.Rad_E_gamma[k], dE);
#ifdef RT_EVOLVE_ENERGY
        Kokkos::atomic_add(&Cj.Rad_E_gamma_Pred[k], dE);
#endif
#ifdef RT_REPROCESS_INJECTED_PHOTONS
        if (donation_bin > -1) {
            Kokkos::atomic_add(&Cj.Rad_E_gamma[donation_bin], dE_donation);
#ifdef RT_EVOLVE_ENERGY
            Kokkos::atomic_add(&Cj.Rad_E_gamma_Pred[donation_bin], dE_donation);
#endif
        }
#endif
#ifdef RT_EVOLVE_INTENSITIES
        {
            double dflux_int = dE / angle_wt_Inu_sum;
            for (int kv = 0; kv < N_RT_INTENSITY_BINS; kv++) {
                double dI_temp = dflux_int * angle_wt_Inu[kv];
                Kokkos::atomic_add(&Cj.Rad_Intensity[k][kv], dI_temp);
                Kokkos::atomic_add(&Cj.Rad_Intensity_Pred[k][kv], dI_temp);
            }
        }
#endif
#if defined(RT_EVOLVE_FLUX)
        {
            for (int kv = 0; kv < 3; kv++) {
                dfluxes[kv] += dE * (RSOL_CORRECTION_FACTOR_FOR_VELOCITY_TERMS * local.Vel[kv] / All.cf_atime);
            }
        }
        for (int kv = 0; kv < 3; kv++) {
            Kokkos::atomic_add(&Cj.Rad_Flux[k][kv],      (MyFloat)dfluxes[kv]);
            Kokkos::atomic_add(&Cj.Rad_Flux_Pred[k][kv], (MyFloat)dfluxes[kv]);
        }
#endif /* RT_EVOLVE_FLUX */
#endif /* RT_INJECT_PHOTONS_DISCRETELY */
    } /* for k < N_RT_FREQ_BINS */
}

/* ============================================================================
 * RtSrcInjectionSpec — non-iterative scatter Spec (thermal_fb pattern).
 * ========================================================================== */
struct RtSrcInjectionSpec {
    static constexpr const char *loop_name = "rtsrcinjection";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::EpsilonAtomic; /* EpsilonAtomic: radiation scatter atomic_add to Cj.Rad fields, Pj.Vel, Pj.dp; source-snapshot deltas never read back -> ulp */

    /* Search policy. SYMMETRIC matches legacy rt_source_injection_gpu.cc:178
     * (NGB_SEARCH_SYMMETRIC unconditional). Correctness-required under
     * RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION; benign otherwise. */
    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int            neighbor_type_mask = (1u << 0);   /* gas only */
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;

    /* Write policy. */
    static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::GasOnly;
    static constexpr bool mode_a_active_sources_in_sidx_pool = false; /* non-pool active sources (sink/star/grain) -> runner stages explicit P[].Pos */
    static constexpr bool           uses_ghost_writeback      = true;
    static constexpr bool writes_neighbour_motion = true;   /* the injected radiation momentum kicks the neighbour */
    static constexpr bool           uses_ghost_write_detector = true;

    /* Real precision bound (not algorithmic-parity floor): the pair kernel
     * writes only into Pj/Cj fields that this loop never reads back during
     * the loop. No thermal_fb-style read-then-write self-coupling. */

    /* Type aliases. */
    using CallScalars    = RtSrcInjCallScalars;
    using ActiveData     = RtSrcInjActiveState;
    using AccumData      = RtSrcInjAccum;
    using DeviceContext  = RtSrcInjDeviceContext;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

    /* NeighborData: non-const pointers (kernel writes to *neighbor_*). */
    struct NeighborData {
        struct particle_data *neighbor_particle;
        struct gas_cell_data *neighbor_cell;
    };

    /* Aux — toplevel-owned active source pack. host_locals[active_slot] is
     * the per-source pre-fill (clears P[i].Sink_accreted_photon_energy as a
     * side effect at fill time; happens once per source in the toplevel
     * pre-runner build, NOT in any Spec hook). */
    struct Aux {
        const RtSrcLocalIn *host_locals;     /* [num_active]; toplevel-owned */
        const int          *host_active_src; /* [num_active]; particle indices */
    };

    /* ====================================================================
     * Host hooks (bodies in rt_source_injection_loop.cc).
     * ==================================================================== */

    /* Defensive predicate. NEVER called by nlr_build_active_list (the
     * toplevel in rt_source_injection.cc builds the active list directly via
     * the existing rt_sourceinjection_active_check chain). This is a
     * belt-and-suspenders re-check. */
    static bool        is_active(int particle_index);

    /* Per-active search radius = the source's un-scaled KernelRadius (reads from
     * aux->host_locals[active_slot]). The SYMMETRIC search discovers the j-side
     * (r < h_j) reach independently, so the source radius is NOT inflated to
     * find large-h_j gas. */
    static double      search_radius(const neighbor_loop_args& args,
                                      int active_slot, int i);

    static CallScalars populate_call_scalars(const neighbor_loop_args& args);

    static void populate_device_context(const neighbor_loop_args& args, DeviceContext& ctx);
    static void cleanup_device_context (const neighbor_loop_args& args, DeviceContext& ctx);

    /* Source-side host write — no-op for rt_source_injection (the legacy code
     * had no post-runner source-side mutation; Sink_accreted_photon_energy is
     * cleared during pre-runner host fill instead). */
    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int active_slot, int i,
                                        const AccumData& accum);

    /* Additive merge — manifest in rt_source_injection_loop.cc. */
    /* ============================================================================
     * MERGE — nothing to merge; AccumData is empty for this Spec (see its
     * definition in rt_source_injection_loop.h). Declared because the Spec
     * contract requires it.
     * ========================================================================== */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& /*local_accum*/, const AccumData& /*peer_accum*/)
    {
    }

    /* Ghost-writeback + write-detector bookkeeping.
     * Detector uses runner default (loop_name = "rtsrcinjection"). */
    static void ghost_writeback_begin     (const neighbor_loop_args&, const NeighborLoopPlan&);
    static void ghost_writeback_end       (const neighbor_loop_args&, const NeighborLoopPlan&);

    /* Diagnostics — env-gated. */

    /* ====================================================================
     * Device hooks (header-inline).
     * ==================================================================== */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& /*accum*/) {
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& dctx,
                                  int active_slot, int i,
                                  double h_search,
                                  const CallScalars& scalars) {
        ActiveData a{};
        if (dctx.per_active_local != nullptr) {
            a.local = dctx.per_active_local[active_slot];
            a.pos[0] = (double)a.local.Pos[0];
            a.pos[1] = (double)a.local.Pos[1];
            a.pos[2] = (double)a.local.Pos[2];
            /* h_search comes from the runner (the un-scaled source radius; the
             * SYMMETRIC search adds the j-side reach). */
            a.h_search = h_search;
        } else {
            a.pos[0] = (double)dctx.P[i].Pos[0];
            a.pos[1] = (double)dctx.P[i].Pos[1];
            a.pos[2] = (double)dctx.P[i].Pos[2];
            a.h_search = h_search;
        }
        (void)scalars;
        return a;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& dctx, int j,
                                       const IdentitySidecar& /*id*/,
                                       const ActiveData& /*active*/) {
        NeighborData n{};
        n.neighbor_particle = &dctx.P[j];
        n.neighbor_cell     = (dctx.CellP != nullptr && dctx.P[j].Type == 0)
                              ? &dctx.CellP[j] : nullptr;
        return n;
    }

    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                             const NeighborData& neighbor,
                             AccumData& /*accum*/,
                             NoScatter& /*scatter*/,
                            const CallScalars& /*cs*/) {
        /* Source-level early-outs — mirror the legacy lambda guard at
         * rt_source_injection_gpu.cc:193. wk = (1 - r2/h2) / KernelSum_Around_RT_Source
         * would divide by zero on degenerate sources. */
        if (active.local.KernelRadius              <= 0) return;
        if (active.local.KernelSum_Around_RT_Source <= 0) return;

        if (neighbor.neighbor_particle == nullptr) return;
        if (neighbor.neighbor_cell     == nullptr) return;   /* gas-only safety */
        struct particle_data &Pj = *neighbor.neighbor_particle;
        struct gas_cell_data &Cj = *neighbor.neighbor_cell;

        /* Defensive Type / Mass re-checks (mask should have ensured Type==0). */
        if (Pj.Type != 0)  return;
#ifdef HYDRO_MULTIFLUID_DM
        if (Pj.FluidType == FLUID_DM) return; /* skip dark-fluid neighbors */
#endif
        if (Pj.Mass <= 0)  return;

        Vec3<double> dp;
        dp[0] = (double)active.local.Pos[0] - (double)Pj.Pos[0];
        dp[1] = (double)active.local.Pos[1] - (double)Pj.Pos[1];
        dp[2] = (double)active.local.Pos[2] - (double)Pj.Pos[2];
        nearest_xyz(dp);
        const double r2 = dp.norm_sq();
        if (r2 <= 0) return;

        /* Rejection mirrors legacy lambda at rt_source_injection_gpu.cc:205-208.
         * KernelRadius (loc.KernelRadius) is the source's un-scaled radius, the
         * same radius the search uses; the extended-reach branch below accepts
         * the j-side (r < h_j) pairs the SYMMETRIC search discovers.
         * The r2>=h2 exclusion must be unconditional: All.TimeStep reads 0 at
         * the first step after any sync-point, and gating the exclusion on it
         * let far-outside-kernel pairs through with a negative kernel weight. */
        const double h2 = (double)active.local.KernelRadius * (double)active.local.KernelRadius;
#ifdef RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION
        if (r2 >= h2) {
            bool extended_reach = (All.TimeStep > 0) && (r2 < (double)Pj.KernelRadius * (double)Pj.KernelRadius);
            if (!extended_reach) return;
        }
#else
        if (r2 >= h2) return;
#endif
#ifdef SINK_WIND_SPAWN
        if (Pj.StellarAge == All.Time) return;
#endif

        rt_source_injection_pair_body(active.local, Pj, Cj, r2, dp);
    }
};

#endif /* RT_SOURCE_INJECTION */

#endif /* RT_SOURCE_INJECTION_LOOP_H */
