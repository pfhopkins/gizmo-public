/* galaxy_sf/thermal_fb_loop.h — ThermalFBSpec for the runner-template port of
 * thermal_fb_evaluate_gpu.
 *
 * Non-iterative scatter loop: active Type-4 (stellar) particles with
 * SNe_ThisTimeStep>0 deposit thermal energy / ejecta mass into surrounding gas
 * neighbors. Pattern parallels SinkFeedSpec (ActiveReduceOnly source-side
 * output + ghost_writeback bundle for j-side writes). Single CSR; no mode
 * machine; no per-iter scratch.
 *
 * SSOT for thermal_fb physics types and per-pair kernel — supersedes the
 * pre-port `thermal_fb_functions.h` and the legacy `thermal_fb_gpu.cc`,
 * both retired.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */
#ifndef THERMAL_FB_LOOP_H
#define THERMAL_FB_LOOP_H

/* Kokkos_Core must precede allvars.h (its macros may conflict with stdlib
 * names). */
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"

#ifdef GALSF_FB_THERMAL

/* NOTE: caller TUs must `#include "../mesh/kernel.h"` BEFORE this header.
 * kernel.h has no include guards (defines static inline kernel_main / kernel_hinv
 * used by the inline thermal_fb_pair_kernel body below) — pulling it from inside
 * this header would multiply-define those symbols in any caller that already
 * has kernel.h on its include list (runner / _gpu.cc / _loop.cc / .cc all do).
 * Same pattern as sinks/sink_feed_loop.h. */
#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"   /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* Host-side helpers (defined in thermal_fb_loop.cc). Forward-declared on the
 * incomplete struct ThermalFBLocalIn (definition follows below). SSOT —
 * single owner for thermal_fb host-side fill + per-call scalar build. */
struct ThermalFBLocalIn;
struct ThermalFBCallScalars;
void thermal_fb_local_fill(int i,
                            struct particle_data *P_host,
                            struct gas_cell_data *CellP_host,
                            struct ThermalFBLocalIn *loc);
/* Per-call cosmology + run-invariant unit conversion factors. Routes All.*
 * reads through nlr_host_all_ptr() for explicit host-snapshot intent;
 * correct from any TU under the device-pass-only redirect. SSOT —
 * Spec::populate_call_scalars forwards here. */
ThermalFBCallScalars thermal_fb_build_call_scalars(void);

/* ============================================================================
 * Per-pair physics types.
 *
 * ThermalFBCallScalars carries the cosmology + run-invariant unit conversion
 * factors needed by the inline pair kernel. Passed to it as its own argument,
 * once per call (TRAP 1 compliance — pair kernel reads NO All.* / UNIT_* macros).
 *
 * thermalfb does not use RNG (legacy kernel uses no random numbers); no
 * per-loop XOR shift on NumCurrentTiStep needed. Documented for future
 * reviewers — see feedback_rng_loop_uniqueness.md.
 * ========================================================================== */
struct ThermalFBCallScalars {
    NlrCommonScalars common;             /* cf_atime, cf_a2inv, cf_a3inv, ... */

    /* Host-precomputed unit conversion factors. UNIT_*_IN_* macros in
     * declarations/constants.h expand to All.UnitMass_in_g / All.HubbleParam
     * / All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm references. We
     * evaluate them ONCE host-side in populate_call_scalars (via the
     * canonical nlr_host_all_ptr()) and ship the doubles to the kernel.
     * Only used under GALSF_FB_TURNOFF_COOLING; populated unconditionally
     * for stable layout. */
    double unit_energy_in_cgs;
    double unit_density_in_nhcgs;
    double u_to_temp_units;
    double unit_time_in_myr;
};

/* Per-source input pack — host-filled, UVM-staged, read on device by load_active.
 * Trivially copyable. Field set matches the legacy ThermalFBLocalIn (which
 * lived in `thermal_fb_functions.h` pre-port).
 *
 * `Vel` is added (vs the legacy struct) so apply_active_writeback can compute
 * the source-side `dp -= M_coupled * Vel` without re-reading P_host. Avoids
 * an extra P_host indirection in the per-active post-runner loop. Type is
 * Vec3<MyDouble> (NOT Vec3<MyFloat>) to match P[i].Vel's declared precision
 * — narrowing to float here would silently downgrade source-side momentum
 * book-keeping vs legacy. */
struct ThermalFBLocalIn {
    Vec3<MyDouble> Pos;
    Vec3<MyDouble> Vel;
    MyFloat KernelRadius;
    MyFloat wt_sum;       /* Sum_j W_j (number-weighted kernel sum), reconstructed
                           * from NumNgb + KernelRadius; mass-invariant normalizer */
    MyFloat Msne;         /* ejecta mass */
    MyFloat Esne;         /* ejecta kinetic energy = 0.5*Msne*v_ej^2 */
    MyFloat kernel_zero;  /* kernel_main(0.0, 1.0, 1.0, ...) at u=0 */
    short int TimeBin;    /* source.TimeBin -- encoded as wakeup_val
                           * (TimeBin+1) on the receiver via atomic_max in the
                           * pair kernel, so direct thermal-energy receivers
                           * are marked for positive wakeup at deposition. */
#ifdef METALS
    MyFloat yields[NUM_METAL_SPECIES];
#endif
};

/* Per-source output — additive over neighbors. Applied to source P[i] post-
 * runner by apply_active_writeback. */
struct ThermalFBOut {
    MyDouble M_coupled;
};

/* Runner ActiveData. `pos` and `h_search` are top-level (runner's Mode B
 * walker reads them directly: see neighbor_loop_runner.cc:748,753). */
struct ThermalFBActiveState {
    Vec3<double>         pos;             /* P[i].Pos in double — runner reads directly */
    double               h_search;        /* runner-supplied per-active radius */
    ThermalFBLocalIn     local;
};

/* DeviceContext extension: UVM-resident per-active host-fill array.
 * Trivially copyable; runner captures by value into device lambdas. */
struct ThermalFBDeviceContext : NeighborLoopDeviceContextBase {
    const ThermalFBLocalIn *per_active_local;   /* UVM, [num_active]; nullptr when num_active==0 */
    unsigned char          *wakeup_dirty_base;  /* WakeupDirty sidecar base (global UVM); populate sets from WakeupDirty */
};

/* ============================================================================
 * Inline pair body — SSOT for thermal_fb per-pair physics.
 *
 * Caller pre-filters Pj.Type==0 / Pj.Mass>0 via the gas-only neighbor mask in
 * load_neighbor (defensive checks still inside as safety). dp = source.Pos -
 * Pj.Pos with nearest_xyz applied. r2 > 0 and r2 < KernelRadius^2 pre-checked
 * by the lambda.
 *
 * All `All.*` and `UNIT_*` macro reads go through `scalars` — Mode B no-
 * globals rule (project directive 3).
 * ========================================================================== */
KOKKOS_INLINE_FUNCTION
static void thermal_fb_pair_kernel(
    const ThermalFBLocalIn& local,
    const ThermalFBCallScalars& scalars,
    struct particle_data& Pj,
    struct gas_cell_data& Cj,
    double r2,
    unsigned char* wakeup_dirty_slot,
    ThermalFBOut& out)
{
    if (Pj.Type != 0) return;
#ifdef HYDRO_MULTIFLUID_DM
    if (Pj.FluidType == FLUID_DM) return; /* skip dark-fluid neighbors */
#endif
    /* Other sources may be depositing into this cell at the same time, through
     * the atomic adds below, so its mass, density and metallicity are read
     * atomically: the value seen is one of the deposits' results, never a torn
     * mixture. Which one is seen is the accepted order dependence. */
    double Mass_j = (double)Kokkos::atomic_load(&Pj.Mass);
    if (Mass_j <= 0) return;
    /* Belt-and-suspenders guards — Spec::pair_kernel already short-circuits
     * the source if any of (KernelRadius, wt_sum, Msne) is non-positive, but
     * the legacy thermal_fb_gpu.cc:183 + per-pair body both checked, so we
     * keep both layers (and any future direct caller stays safe — wt_sum
     * goes in the denominator below). */
    if (local.KernelRadius <= 0) return;
    if (local.wt_sum       <= 0) return;
    if (local.Msne         <= 0) return;
    double h2 = (double)local.KernelRadius * (double)local.KernelRadius;
    if (r2 >= h2 || r2 <= 0) return;

    double hinv, hinv3, hinv4;
    kernel_hinv((double)local.KernelRadius, &hinv, &hinv3, &hinv4);
    double r = sqrt(r2);
    double u = r * hinv;
    double wk, dwk;
    if (u < 1) { kernel_main(u, hinv3, hinv4, &wk, &dwk, 0); } else { wk = dwk = 0; }
    if (wk <= 0 || wk != wk) return;

    double rho_j_0 = (double)Kokkos::atomic_load(&Cj.Density);

    /* Number-weighted fraction of this neighbor's claim on the ejecta:
     * W_j / Sum_k W_k. Mass-invariant, so robust to the gas-mass growth the
     * deposits cause during the loop (local.wt_sum holds Sum_k W_k). */
    wk = wk / (double)local.wt_sum;

    /* Mass injected into j. */
    double dM = wk * (double)local.Msne;

    /* Delta density. Branches preserved from the legacy kernel (rho_j_0 sign
     * + per-j kernel radius decide the smoothing length). */
    double delta_rho = 0;
    if (Pj.KernelRadius <= 0) {
        if (rho_j_0 > 0) { delta_rho = rho_j_0 * dM / Mass_j; }
        else             { delta_rho = (double)local.kernel_zero * dM * hinv3; }
    } else {
        double kr3 = (double)Pj.KernelRadius * (double)Pj.KernelRadius
                   * (double)Pj.KernelRadius;
        delta_rho = (double)local.kernel_zero * dM / kr3;
    }

    /* Energy: IE += wk * Esne / Mass_j. */
    double dIE = wk * (double)local.Esne / Mass_j;

    out.M_coupled += dM;

    /* j-side atomic writes. */

    Kokkos::atomic_add(&Pj.Mass, (MyDouble)dM);
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    Kokkos::atomic_add(&Cj.MassTrue, (MyDouble)dM);
#endif
    if (delta_rho != 0) {
        Kokkos::atomic_add(&Cj.Density, (MyDouble)delta_rho);
    }
    for (int k = 0; k < 3; k++) {
        Kokkos::atomic_add(&Pj.dp[k], (MyDouble)(dM * (double)Pj.Vel[k]));
    }
    Kokkos::atomic_add(&Cj.InternalEnergy,     (MyDouble)dIE);
    Kokkos::atomic_add(&Cj.InternalEnergyPred, (MyDouble)dIE);

    /* Mark this directly-perturbed gas receiver for positive wakeup using
     * the standard hydro-convention encoding (source.TimeBin + 1). atomic_max
     * preserves any earlier-written more-aggressive wakeup. Coupled with the
     * process_wake_ups floor at lowest_occupied_active_bin, the receiver
     * wakes at the right pace immediately, killing the multi-
     * generation cascade that would otherwise emerge once shell-by-shell
     * pair-body wakeup-checks discover the stale-low MaxSignalVel. */
    Kokkos::atomic_max(&Pj.wakeup, (short int)((int)local.TimeBin + 1));
    if (wakeup_dirty_slot) { *wakeup_dirty_slot = 1; }   /* dirty-sidecar mark (byte store; race-benign) */

#ifdef METALS
    for (int k = 0; k < NUM_METAL_SPECIES; k++) {
        double dMet = (dM / Mass_j)
                    * ((double)local.yields[k] - (double)Kokkos::atomic_load(&Pj.Metallicity[k]));
        Kokkos::atomic_add(&Pj.Metallicity[k], (MyFloat)dMet);
    }
#endif

#ifdef GALSF_FB_TURNOFF_COOLING
    {
        /* Sedov turnoff time. All All.* / UNIT_* legacy macro reads now route
         * through scalars (Mode B no-globals + GPU TU All_dev safety). */
        double Esne51 = (double)local.Esne * scalars.unit_energy_in_cgs / 1.e51;
        double density_to_n   = scalars.common.cf_a3inv * scalars.unit_density_in_nhcgs;
        double pressure_to_p4 = density_to_n * scalars.u_to_temp_units / 1.0e4;
        double dt_ram = 7.08
                      * pow(Esne51 * rho_j_0 * density_to_n, 0.34)
                      * pow((double)Cj.Pressure * pressure_to_p4, -0.70)
                      / scalars.unit_time_in_myr;
        Kokkos::atomic_max(&Cj.DelayTimeCoolingSNe, (MyDouble)dt_ram);
    }
#endif
}

/* ============================================================================
 * ThermalFBSpec — non-iterative scatter Spec (sink_feed pattern).
 * ========================================================================== */
struct ThermalFBSpec {
    /* Identity. loop_name labels this loop in runner diagnostics and the
     * generic _end_bundle print. */
    static constexpr const char *loop_name = "thermalfb";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::EpsilonAtomic; /* several sources may deposit into one cell; each reads the live mass, density and metallicity before its atomic add, so the result depends on their order -- small where a cell takes O(1) small deposits per step, and accepted at that level */
    /* Each deposit's specific energy is the source's share of its energy over the
       cell's mass at that moment, so the cell's live mass is what the kernel must
       see. On a ghost copy that mass lags the owner's -- it misses every other
       rank's deposits -- and the shipped-home delta then lands on a heavier cell
       and injects more than the share. So this loop is never evaluated on ghost
       copies: the host walker or Mode D answer it where the cell lives. */
    static constexpr bool needs_live_neighbours = true;

    /* Search policy. Legacy thermal_fb_gpu.cc:166 used NGB_SEARCH_ONEWAY +
     * j_type_bitmask=1 (gas only). */
    static constexpr int                     search_mode        = MODE_B_SEARCH_ONEWAY;
    static constexpr unsigned int            neighbor_type_mask = (1u << 0);   /* gas only */
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;

    /* Write policy. */
    static constexpr WritePattern   write_pattern   = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind = SidxCacheKind::GasOnly;   /* tbm = 1 */
    static constexpr bool           uses_ghost_writeback      = true;
    static constexpr bool           uses_ghost_write_detector = true;

    /* Source-order independence of M_coupled. The ejecta are split by the
     * number-weighted kernel fraction W_j / Sum_k W_k, which depends only on
     * neighbor geometry and the source kernel radius — NOT on the gas mass the
     * deposits mutate during the loop. So each source couples exactly
     * M_coupled = Msne no matter what order the sources are processed in; the
     * only residual is floating-point summation order across CSR / Mode-B
     * partial sums. */

    /* Type aliases. */
    using CallScalars    = ThermalFBCallScalars;
    using ActiveData     = ThermalFBActiveState;
    using AccumData      = ThermalFBOut;
    using DeviceContext  = ThermalFBDeviceContext;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

    /* NeighborData carries non-const pointers (kernel writes to *neighbor_*);
     * mirrors SinkFeedSpec::NeighborData. */
    struct NeighborData {
        struct particle_data *neighbor_particle;
        struct gas_cell_data *neighbor_cell;   /* nullptr for non-gas; gas-only mask should prevent */
        unsigned char        *wakeup_dirty_slot; /* &WakeupDirty[j] */
    };

    /* Aux — empty; thermalfb needs no host-only per-call state beyond what
     * sits in DeviceContext / per-active UVM. Reserved struct (non-empty so
     * runner aux machinery has a well-defined type). */
    struct Aux {
        const ThermalFBLocalIn *host_locals;   /* [num_active]; owned by toplevel */
    };

    /* ====================================================================
     * Host hooks (bodies in thermal_fb_loop.cc).
     * ==================================================================== */
    static bool        is_active(int particle_index);
    static double      search_radius(const neighbor_loop_args& args,
                                      int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);

    static void populate_device_context(const neighbor_loop_args& args, DeviceContext& ctx);
    static void cleanup_device_context (const neighbor_loop_args& args, DeviceContext& ctx);

    /* Source-side host write — applied per-active post-runner. Translates the
     * legacy post-kernel loop in thermal_fb_gpu.cc:236-246 (Mass -= M_coupled;
     * dp -= M_coupled * Vel; clamp non-finite). Reads ActiveData::local.Vel
     * for the source velocity (carried in the host-fill pack so no extra
     * P_host indirection here). */
    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int active_slot, int i,
                                        const AccumData& accum);

    /* Per-field merge — manifest in thermal_fb_loop.cc. */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local_accum, const AccumData& peer_accum)
    {
#define ACCUM_ADD(field)  local_accum.field += peer_accum.field;
        ACCUM_ADD(M_coupled)
#undef ACCUM_ADD
    }

    /* Ghost-writeback + write-detector bookkeeping.
     * Detector uses runner default (loop_name = "thermalfb"). */
    static void ghost_writeback_begin     (const neighbor_loop_args&, const NeighborLoopPlan&);
    static void ghost_writeback_end       (const neighbor_loop_args&, const NeighborLoopPlan&);

    /* Diagnostics — env-gated. */

    /* ====================================================================
     * Device hooks (header-inline).
     * ==================================================================== */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        accum.M_coupled = 0;
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
            a.h_search = (double)a.local.KernelRadius;
        } else {
            /* Fallback when per_active_local wasn't staged (shouldn't fire in
             * production — populate_device_context always allocates when
             * num_active>0). Pull straight from ctx.P. */
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
        n.wakeup_dirty_slot = dctx.wakeup_dirty_base ? &dctx.wakeup_dirty_base[j] : nullptr;
        return n;
    }

    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                             const NeighborData& neighbor,
                             AccumData& accum,
                             NoScatter& /*scatter*/,
                            const CallScalars& cs) {
        /* Source-level early-out — mirrors the legacy thermal_fb_gpu.cc:183
         * `if(loc.KernelRadius<=0 || loc.wt_sum<=0 || loc.Msne<=0) return;`
         * which short-circuited the ENTIRE neighbor loop for the source.
         * Re-checked per-pair in the runner template because the runner
         * dispatches pair_kernel per (source, neighbor); semantics-equivalent.
         * wt_sum is in the denominator of the kernel weight, so missing this
         * guard would division-by-zero / produce inf on degenerate sources. */
        if (active.local.KernelRadius <= 0) return;
        if (active.local.wt_sum       <= 0) return;
        if (active.local.Msne         <= 0) return;

        if (neighbor.neighbor_particle == nullptr) return;
        if (neighbor.neighbor_cell     == nullptr) return;   /* gas-only safety */
        struct particle_data &Pj = *neighbor.neighbor_particle;
        struct gas_cell_data &Cj = *neighbor.neighbor_cell;

        Vec3<double> dp;
        dp[0] = (double)active.local.Pos[0] - (double)Pj.Pos[0];
        dp[1] = (double)active.local.Pos[1] - (double)Pj.Pos[1];
        dp[2] = (double)active.local.Pos[2] - (double)Pj.Pos[2];
        nearest_xyz(dp);
        const double r2 = dp.norm_sq();

        /* Pre-filter mirrors the legacy thermal_fb_gpu.cc:188-196 lambda
         * (Type==0 and Mass>0 already checked above via the gas-only mask +
         * nullptr guard; r2>0 / r2<h2 filtered here so the inline body's
         * guards become safety-only). */
        if (r2 <= 0) return;
        const double h2 = (double)active.local.KernelRadius
                        * (double)active.local.KernelRadius;
        if (r2 >= h2) return;

        thermal_fb_pair_kernel(active.local, cs, Pj, Cj,
                                r2, neighbor.wakeup_dirty_slot, accum);
    }
};

#endif /* GALSF_FB_THERMAL */

#endif /* THERMAL_FB_LOOP_H */
