/* turb/difffilter_loop.h — DiffFilterSpec + DynDiffSpec for the runner-template
 * port of difffilter_evaluate_gpu / dynamicdiff_evaluate_gpu (Wave 5, item 1).
 *
 * Both are non-iterative, scaled-symmetric gas-gas loops, pure i-side reduce
 * (no j-side writes, no ghost writeback). The symmetric search reaches
 * max(fac*h_i, fac*h_j) via the Spec hook symmetric_neighbor_radius_scale()
 * — see OPEN_3d_difffilter_design.md §2-3.
 *
 * SSOT for the DiffFilter / DynamicDiff per-pair physics — supersedes the
 * legacy turb/difffilter_gpu.cc (retired in this commit).
 *
 * Caller TUs must `#include "../mesh/kernel.h"` BEFORE this header (kernel.h
 * has no include guards; the inline pair kernels below call kernel_main /
 * kernel_hinv). Same pattern as sink_feed_loop.h / thermal_fb_loop.h.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) and Claude for GIZMO.
 */
#ifndef DIFFFILTER_LOOP_H
#define DIFFFILTER_LOOP_H

/* Kokkos_Core before allvars.h — allvars pulls macros.h which #defines
 * terminate(...), mangling std::terminate in Kokkos_Core. Same as
 * dm_dispersion_loop.h / thermal_fb_loop.h. */
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"
#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"   /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
#include "difffilter_loop_api.h"           /* temporary_data_dyndiff, toplevels */

#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

#ifdef TURB_DIFF_DYNAMIC

/* ============================================================================
 * DiffFilterSpec — velocity-smoothing filter (legacy difffilter_evaluate_gpu).
 * ========================================================================== */

/* CallScalars: cosmology common block + the dynamic-diffusion widening factor. */
struct DiffFilterCallScalars {
    NlrCommonScalars common;
    double           turb_dyn_diff_fac;   /* All.TurbDynamicDiffFac */
};

/* AccumData — per-active reduction. norm_hat / velocity_bar_delta additive;
 * filter_width_bar / max_dist_for_grad are max reductions (start 0; r>0). */
struct DiffFilterAccum {
    double norm_hat;
    double velocity_bar_delta[3];
    double filter_width_bar;
    double max_dist_for_grad;
};

/* ActiveData — i-side snapshot. pos / h_search are top-level (Mode B walker
 * reads them directly). `valid` carries the legacy density>0 / h>0 kernel
 * guard (difffilter_gpu.cc:101): an invalid source no-ops every pair, leaving
 * accum at zero_accum — bit-identical to the legacy memset(out,0). */
struct DiffFilterActiveState {
    Vec3<double> pos;
    double       h_search;            /* runner radius = fac * h_kernel       */
    double       h_kernel;            /* P[i].KernelRadius (narrow hydro h)   */
    Vec3<double> vel_pred;            /* CellP[i].VelPred                     */
    double       density;             /* CellP[i].Density                     */
    double       turb_dyn_diff_fac;
#ifdef GALSF_SUBGRID_WINDS
    double       delay_time;          /* CellP[i].DelayTime                   */
#endif
    int          valid;
};

/* NeighborData — j-side gas particle + cell (read-only). */
struct DiffFilterNeighborData {
    struct particle_data *neighbor_particle;
    struct gas_cell_data *neighbor_cell;
};

struct DiffFilterSpec {
    static constexpr const char *loop_name = "difffilter";

    /* Scaled-symmetric gas-gas: search reaches max(fac*h_i, fac*h_j). */
    static constexpr int                    search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int           neighbor_type_mask = (1u << 0);  /* gas */
    static constexpr mode_b_radius_policy_t radius_policy       = MODE_B_RADIUS_DEFAULT;

    static constexpr WritePattern  write_pattern             = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind sidx_cache_kind           = SidxCacheKind::GasOnly;
    static constexpr bool          uses_ghost_writeback      = false;
    static constexpr bool          uses_ghost_write_detector = false;

    /* Pure-read pair kernel (no j-writes) → no order dependence → tight oracle. */
    static constexpr double accum_tolerance = 1e-10;

    using CallScalars    = DiffFilterCallScalars;
    using ActiveData     = DiffFilterActiveState;
    using AccumData      = DiffFilterAccum;
    using NeighborData   = DiffFilterNeighborData;
    using DeviceContext  = NeighborLoopDeviceContextBase;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

    /* Aux unused — DiffFilter scatters directly to CellP in apply_active_writeback. */
    struct Aux {};

    static_assert(uses_ghost_writeback == false,
        "DiffFilter is pure i-side accumulation; no j-side writes.");

    /* Host hooks — bodies in difffilter_loop.cc. */
    static bool        is_active(int i);
    static double      search_radius(const neighbor_loop_args& args, int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);
    static void        apply_active_writeback(const neighbor_loop_args& args,
                                              int active_slot, int i,
                                              const AccumData& accum);
    static void        merge_accum(AccumData& local, const AccumData& peer);
    static double      compare_accum(const AccumData& local, const AccumData& oracle);
    static void        set_oracle_brute_pass(DeviceContext& ctx, bool on);

    /* Scaled-symmetric search hook: j-side radius scaled by All.TurbDynamicDiffFac. */
    static double      symmetric_neighbor_radius_scale();

    /* ---- Device hooks (header-inline; runner instantiates from GPU TUs). ---- */
    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        accum.norm_hat            = 0;
        accum.velocity_bar_delta[0] = 0;
        accum.velocity_bar_delta[1] = 0;
        accum.velocity_bar_delta[2] = 0;
        accum.filter_width_bar    = 0;
        accum.max_dist_for_grad   = 0;
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx,
                                  int active_slot, int i,
                                  double h_search, const CallScalars& scalars) {
        (void)active_slot;
        ActiveData a{};
        a.pos[0]   = (double)ctx.P[i].Pos[0];
        a.pos[1]   = (double)ctx.P[i].Pos[1];
        a.pos[2]   = (double)ctx.P[i].Pos[2];
        a.h_search = h_search;
        a.h_kernel = (double)ctx.P[i].KernelRadius;
        a.turb_dyn_diff_fac = scalars.turb_dyn_diff_fac;
        a.vel_pred[0] = (double)ctx.CellP[i].VelPred[0];
        a.vel_pred[1] = (double)ctx.CellP[i].VelPred[1];
        a.vel_pred[2] = (double)ctx.CellP[i].VelPred[2];
        a.density     = (double)ctx.CellP[i].Density;
#ifdef GALSF_SUBGRID_WINDS
        a.delay_time  = (double)ctx.CellP[i].DelayTime;
#endif
        /* Legacy difffilter_gpu.cc:101 guard — invalid source emits zero. */
        a.valid = (a.density > 0 && a.h_kernel > 0) ? 1 : 0;
        return a;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx, int j,
                                      const IdentitySidecar& /*id*/,
                                      const ActiveData& /*active*/) {
        NeighborData n{};
        n.neighbor_particle = &ctx.P[j];
        n.neighbor_cell     = (ctx.CellP != nullptr && ctx.P[j].Type == 0)
                              ? &ctx.CellP[j] : nullptr;
        return n;
    }

    /* pair_kernel — port of difffilter_gpu.cc:117-171 (per-pair form). */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& a,
                            const NeighborData& nb,
                            AccumData& accum, NoScatter& /*scatter*/) {
        if (!a.valid) return;
        if (nb.neighbor_particle == nullptr || nb.neighbor_cell == nullptr) return;
        struct particle_data& Pj = *nb.neighbor_particle;
        struct gas_cell_data& Cj = *nb.neighbor_cell;

#ifdef GALSF_SUBGRID_WINDS
        const double dt_j = (double)Cj.DelayTime;
        if (a.delay_time == 0 && dt_j  > 0) return;
        if (a.delay_time  > 0 && dt_j == 0) return;
#endif
        if (Pj.Mass <= 0) return;
        const double density_j = (double)Cj.Density;
        if (density_j <= 0) return;

        Vec3<double> dp;
        dp[0] = a.pos[0] - (double)Pj.Pos[0];
        dp[1] = a.pos[1] - (double)Pj.Pos[1];
        dp[2] = a.pos[2] - (double)Pj.Pos[2];
        nearest_xyz(dp);
        const double r2 = dp.norm_sq();
        if (r2 <= 0) return;

        const double h_i = a.h_kernel;
        const double h_j = (double)Pj.KernelRadius;
        const double fac = a.turb_dyn_diff_fac;
        const double mean_weight = 0.5 * (a.density + density_j) / (a.density * density_j);
        const double V_j  = (double)Pj.Mass * mean_weight;

        const double hhat   = fac * h_i;
        const double hhat_j = fac * h_j;
        if (r2 >= hhat * hhat && r2 >= hhat_j * hhat_j) return;

        const double r = sqrt(r2);

        /* Norm_hat: wide-kernel contribution (one-way at fac*h_i). */
        if (r < hhat) {
            double hhatinv, hhatinv3, hhatinv4, wkhat, dwkhat;
            kernel_hinv(hhat, &hhatinv, &hhatinv3, &hhatinv4);
            const double uhat = DMIN(r * hhatinv, 1.0);
            kernel_main(uhat, hhatinv3, hhatinv4, &wkhat, &dwkhat, 0);
            accum.norm_hat += (double)Pj.Mass * wkhat;
        }

        /* Beyond the narrow-kernel reach → no velocity smoothing. */
        if (r2 >= h_i * h_i && r2 >= h_j * h_j) return;

        if (r > accum.filter_width_bar)  accum.filter_width_bar  = r;
        if (r > accum.max_dist_for_grad) accum.max_dist_for_grad = r;

        const double h_avg = 0.5 * (h_i + h_j);
        double hinv, hinv3, hinv4, wk_i, dwk_i;
        kernel_hinv(h_avg, &hinv, &hinv3, &hinv4);
        const double u = DMIN(r * hinv, 1.0);
        if (u < 1) { kernel_main(u, hinv3, hinv4, &wk_i, &dwk_i, 0); }
        else       { wk_i = dwk_i = 0; }
        const double weight_i = wk_i * V_j;

        if (r < h_avg) {
            for (int k = 0; k < 3; k++) {
                accum.velocity_bar_delta[k] +=
                    ((double)Cj.VelPred[k] - a.vel_pred[k]) * weight_i;
            }
        }
    }
};

/* ============================================================================
 * DynDiffSpec — dynamic Smagorinsky coefficient loop (legacy
 * dynamicdiff_evaluate_gpu). One runner pass per external dynamic_iteration.
 * ========================================================================== */

struct DynDiffCallScalars {
    NlrCommonScalars common;
    double           turb_dyn_diff_fac;   /* All.TurbDynamicDiffFac */
    int              dynamic_iteration;   /* current external-loop index */
};

/* AccumData — fuses the legacy out0 + out_iter. dynamic_fac[/_const] additive
 * over all iterations; filter_width_hat max; the remaining fields are iter-0
 * only (apply_active_writeback gates them). maxima/minima are max/min
 * reductions, zero-initialized — see OPEN_3d_difffilter_design.md §6 for why
 * this fuses legacy kernel-init (±1e30) + host-scatter-from-0 exactly. */
struct DynDiffAccum {
    double dynamic_fac[3][3];
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
    double dynamic_fac_const[3][3];
#endif
    double filter_width_hat;
    double grad_velocity_hat[3][3];
    double dynamic_numerator_hat;
    double dynamic_denominator_hat;
    double product_velocity_hat[3][3];
    double maxima_velocity_hat[3];
    double minima_velocity_hat[3];
};

struct DynDiffActiveState {
    Vec3<double> pos;
    double       h_search;             /* runner radius = fac * kernel_radius */
    double       kernel_radius;        /* P[i].KernelRadius (raw narrow h)    */
    double       turb_dyn_diff_fac;
    int          dynamic_iteration;
    int          sph_gradients_flag;   /* ConditionNumber > CONDITION_NUMBER_DANGER */
    int          valid;               /* Density>0 && KernelRadius>0          */
    double       vel_shear_bar[3][3];
    double       td_dyn_diff_coeff;
    double       mag_shear_bar;
    double       velocity_bar[3];
    double       velocity_hat[3];
    double       norm_hat;
    double       dynamic_numerator;
    double       dynamic_denominator;
    double       filter_width_bar;
    double       prefactor;            /* FilterWidth_bar^2 * TD_DynDiffCoeff * MagShear_bar */
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
    double       prefactor_const;      /* FilterWidth_bar^2 * MagShear_bar     */
#endif
#ifdef GALSF_SUBGRID_WINDS
    double       delay_time;
#endif
};

struct DynDiffNeighborData {
    struct particle_data *neighbor_particle;
    struct gas_cell_data *neighbor_cell;
};

struct DynDiffSpec {
    static constexpr const char *loop_name = "dyndiff";

    static constexpr int                    search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int           neighbor_type_mask = (1u << 0);  /* gas */
    static constexpr mode_b_radius_policy_t radius_policy       = MODE_B_RADIUS_DEFAULT;

    static constexpr WritePattern  write_pattern             = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind sidx_cache_kind           = SidxCacheKind::GasOnly;
    static constexpr bool          uses_ghost_writeback      = false;
    static constexpr bool          uses_ghost_write_detector = false;

    static constexpr double accum_tolerance = 1e-10;

    using CallScalars    = DynDiffCallScalars;
    using ActiveData     = DynDiffActiveState;
    using AccumData      = DynDiffAccum;
    using NeighborData   = DynDiffNeighborData;
    using DeviceContext  = NeighborLoopDeviceContextBase;
    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;

    /* Aux — non-owning ptr to the caller's DynamicDiffDataPasser + the current
     * external-loop iteration index (apply_active_writeback needs both). */
    struct Aux {
        struct temporary_data_dyndiff *dddp;
        int                            dynamic_iteration;
    };

    static_assert(uses_ghost_writeback == false,
        "DynDiff is pure i-side accumulation; no j-side writes.");

    static bool        is_active(int i);
    static double      search_radius(const neighbor_loop_args& args, int active_slot, int i);
    static CallScalars populate_call_scalars(const neighbor_loop_args& args);
    static void        apply_active_writeback(const neighbor_loop_args& args,
                                              int active_slot, int i,
                                              const AccumData& accum);
    static void        merge_accum(AccumData& local, const AccumData& peer);
    static double      compare_accum(const AccumData& local, const AccumData& oracle);
    static void        set_oracle_brute_pass(DeviceContext& ctx, bool on);
    static double      symmetric_neighbor_radius_scale();

    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum) {
        for (int k = 0; k < 3; k++) {
            for (int v = 0; v < 3; v++) {
                accum.dynamic_fac[k][v]          = 0;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
                accum.dynamic_fac_const[k][v]    = 0;
#endif
                accum.grad_velocity_hat[k][v]    = 0;
                accum.product_velocity_hat[k][v] = 0;
            }
            accum.maxima_velocity_hat[k] = 0;
            accum.minima_velocity_hat[k] = 0;
        }
        accum.filter_width_hat        = 0;
        accum.dynamic_numerator_hat   = 0;
        accum.dynamic_denominator_hat = 0;
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const DeviceContext& ctx,
                                  int active_slot, int i,
                                  double h_search, const CallScalars& scalars) {
        (void)active_slot;
        ActiveData a{};
        a.pos[0] = (double)ctx.P[i].Pos[0];
        a.pos[1] = (double)ctx.P[i].Pos[1];
        a.pos[2] = (double)ctx.P[i].Pos[2];
        a.h_search          = h_search;
        a.kernel_radius     = (double)ctx.P[i].KernelRadius;
        a.turb_dyn_diff_fac = scalars.turb_dyn_diff_fac;
        a.dynamic_iteration = scalars.dynamic_iteration;
        const struct gas_cell_data& Ci = ctx.CellP[i];
        for (int k = 0; k < 3; k++) {
            for (int v = 0; v < 3; v++) a.vel_shear_bar[k][v] = (double)Ci.VelShear_bar[k][v];
            a.velocity_bar[k] = (double)Ci.Velocity_bar[k];
            a.velocity_hat[k] = (double)Ci.Velocity_hat[k];
        }
        a.td_dyn_diff_coeff   = (double)Ci.TD_DynDiffCoeff;
        a.mag_shear_bar       = (double)Ci.MagShear_bar;
        a.norm_hat            = (double)Ci.Norm_hat;
        a.dynamic_numerator   = (double)Ci.Dynamic_numerator;
        a.dynamic_denominator = (double)Ci.Dynamic_denominator;
        a.filter_width_bar    = (double)Ci.FilterWidth_bar;
        /* SHOULD_I_USE_SPH_GRADIENTS inlined (no macro coupling) — codex. */
        a.sph_gradients_flag  = ((double)Ci.ConditionNumber > CONDITION_NUMBER_DANGER) ? 1 : 0;
        a.valid = ((double)ctx.P[i].Mass > 0 && (double)Ci.Density > 0
                   && a.kernel_radius > 0) ? 1 : 0;
        /* prefactor — legacy difffilter_gpu.cc:309 (precompute once per active). */
        a.prefactor = a.filter_width_bar * a.filter_width_bar
                    * a.td_dyn_diff_coeff * a.mag_shear_bar;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        a.prefactor_const = a.filter_width_bar * a.filter_width_bar * a.mag_shear_bar;
#endif
#ifdef GALSF_SUBGRID_WINDS
        a.delay_time = (double)Ci.DelayTime;
#endif
        return a;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const DeviceContext& ctx, int j,
                                      const IdentitySidecar& /*id*/,
                                      const ActiveData& /*active*/) {
        NeighborData n{};
        n.neighbor_particle = &ctx.P[j];
        n.neighbor_cell     = (ctx.CellP != nullptr && ctx.P[j].Type == 0)
                              ? &ctx.CellP[j] : nullptr;
        return n;
    }

    /* pair_kernel — port of difffilter_gpu.cc:295-396 (per-pair form). The dead
     * legacy `mass_i` sign-flip (difffilter_gpu.cc:306) is NOT ported. */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& a,
                            const NeighborData& nb,
                            AccumData& accum, NoScatter& /*scatter*/) {
        if (!a.valid) return;
        if (nb.neighbor_particle == nullptr || nb.neighbor_cell == nullptr) return;
        struct particle_data& Pj = *nb.neighbor_particle;
        struct gas_cell_data& Cj = *nb.neighbor_cell;

#ifdef GALSF_SUBGRID_WINDS
        const double dt_j = (double)Cj.DelayTime;
        if (a.delay_time == 0 && dt_j  > 0) return;
        if (a.delay_time  > 0 && dt_j == 0) return;
#endif
        if (Pj.Mass <= 0 || Cj.Density <= 0) return;

        Vec3<double> dp;
        dp[0] = a.pos[0] - (double)Pj.Pos[0];
        dp[1] = a.pos[1] - (double)Pj.Pos[1];
        dp[2] = a.pos[2] - (double)Pj.Pos[2];
        nearest_xyz(dp);
        const double r2 = dp.norm_sq();
        if (r2 <= 0) return;

        const double h_i = a.turb_dyn_diff_fac * a.kernel_radius;
        const double h_j = a.turb_dyn_diff_fac * (double)Pj.KernelRadius;
        if (r2 >= h_i * h_i && r2 >= h_j * h_j) return;

        const double h_avg = 0.5 * (h_i + h_j);
        /* Faithful port — no Norm_hat==0 guard (legacy has none; latent NaN
         * for a gas particle with no wide-kernel neighbors). See design §6. */
        const double mean_weight =
            0.5 * ((double)Cj.Norm_hat + a.norm_hat) / (a.norm_hat * (double)Cj.Norm_hat);
        const double V_j = (double)Pj.Mass * mean_weight;
        const double r   = sqrt(r2);

        if (r > accum.filter_width_hat) accum.filter_width_hat = r;

        double hinv, hinv3, hinv4, wk_i, dwk_i;
        kernel_hinv(h_avg, &hinv, &hinv3, &hinv4);
        const double u = DMIN(r * hinv, 1.0);
        if (u < 1) { kernel_main(u, hinv3, hinv4, &wk_i, &dwk_i, 0); }
        else       { wk_i = dwk_i = 0; }
        const double weight_i = wk_i * V_j;

        const int dyn_iter = a.dynamic_iteration;

        if (dyn_iter == 0) {
            double dv_hat[3];
            for (int k = 0; k < 3; k++) {
                dv_hat[k] = (double)Cj.Velocity_hat[k] - a.velocity_hat[k];
                if (dv_hat[k] < accum.minima_velocity_hat[k]) accum.minima_velocity_hat[k] = dv_hat[k];
                if (dv_hat[k] > accum.maxima_velocity_hat[k]) accum.maxima_velocity_hat[k] = dv_hat[k];
            }
            if (r < a.kernel_radius) {
                double hinv_fg, hinv3_fg, hinv4_fg, wk_fg, dwk_fg;
                kernel_hinv(a.kernel_radius, &hinv_fg, &hinv3_fg, &hinv4_fg);
                const double u_fg = DMIN(r * hinv_fg, 1.0);
                kernel_main(u_fg, hinv3_fg, hinv4_fg, &wk_fg, &dwk_fg, 0);
                if (a.sph_gradients_flag) { wk_fg = -dwk_fg * (double)Pj.Mass / r; }
                for (int k = 0; k < 3; k++) {
                    const double grad_pf = -wk_fg * dp[k];
                    for (int v = 0; v < 3; v++) {
                        accum.grad_velocity_hat[v][k] += grad_pf * dv_hat[v];
                    }
                }
            }
        }

        if (r < h_avg) {
            if (dyn_iter == 0) {
                const double dnum = (double)Cj.Dynamic_numerator   - a.dynamic_numerator;
                const double dden = (double)Cj.Dynamic_denominator - a.dynamic_denominator;
                accum.dynamic_numerator_hat   += dnum * weight_i;
                accum.dynamic_denominator_hat += dden * weight_i;
            }
            const double prefactor_j = (double)Cj.FilterWidth_bar * (double)Cj.FilterWidth_bar
                                     * (double)Cj.TD_DynDiffCoeff * (double)Cj.MagShear_bar;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
            const double prefactor_const_j = (double)Cj.FilterWidth_bar * (double)Cj.FilterWidth_bar
                                           * (double)Cj.MagShear_bar;
#endif
            for (int k = 0; k < 3; k++) {
                for (int v = 0; v < 3; v++) {
                    const double dfac = prefactor_j   * (double)Cj.VelShear_bar[k][v]
                                      - a.prefactor   * a.vel_shear_bar[k][v];
                    accum.dynamic_fac[k][v] += dfac * weight_i;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
                    const double dfac_c = prefactor_const_j * (double)Cj.VelShear_bar[k][v]
                                        - a.prefactor_const  * a.vel_shear_bar[k][v];
                    accum.dynamic_fac_const[k][v] += dfac_c * weight_i;
#endif
                    if (dyn_iter == 0) {
                        const double pv = (double)Cj.Velocity_bar[k] * (double)Cj.Velocity_bar[v]
                                        - a.velocity_bar[k] * a.velocity_bar[v];
                        accum.product_velocity_hat[k][v] += pv * weight_i;
                    }
                }
            }
        }
    }
};

#endif /* TURB_DIFF_DYNAMIC */
#endif /* DIFFFILTER_LOOP_H */
