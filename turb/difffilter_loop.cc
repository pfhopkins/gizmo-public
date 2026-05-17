/* turb/difffilter_loop.cc — DiffFilterSpec / DynDiffSpec host hooks + the two
 * toplevel callers for the runner-template port of difffilter_evaluate_gpu /
 * dynamicdiff_evaluate_gpu (Wave 5, item 1).
 *
 * Device-callable hooks (pair_kernel, zero_accum, load_active, load_neighbor)
 * live in difffilter_loop.h so the runner instantiates them from GPU TUs.
 * This file owns: is_active, search_radius, populate_call_scalars,
 * apply_active_writeback, merge_accum, compare_accum, set_oracle_brute_pass,
 * symmetric_neighbor_radius_scale, and the toplevels.
 *
 * Host-hooks TU — mirrors dm_dispersion_loop.cc: no gpu_all_mirror.h; All.*
 * reads route through nlr_host_all_ptr(). Replaces turb/difffilter_gpu.cc.
 *
 * Written by Phil Hopkins (phopkins@caltech.edu) and Claude for GIZMO.
 */
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <Kokkos_Core.hpp>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"            /* MUST precede difffilter_loop.h (no include guards) */
#include "difffilter_loop.h"

#ifdef TURB_DIFF_DYNAMIC

/* ============================================================================
 * Shared finite-aware per-field relative-difference accumulator for
 * compare_accum. A non-finite value on EITHER side (NaN/Inf — e.g. the
 * deliberate no-guard Norm_hat division in DynDiff::pair_kernel) forces a
 * huge residual so the oracle flags it loudly rather than silently dropping
 * it (a plain `rel > max_rel` compare is false for NaN). The finite test is
 * written inline — (x==x) rejects NaN, (x-x==0) rejects Inf — to avoid the
 * isfinite-macro-vs-std::isfinite ambiguity in GPU TUs (runner checklist §2).
 * ========================================================================== */
namespace {
inline bool nlr_is_finite(double x) { return (x == x) && (x - x == 0.0); }

inline double nlr_rel_update(double max_rel, double a, double b) {
    if (!nlr_is_finite(a) || !nlr_is_finite(b)) return 1e30;
    const double denom = std::fmax(1.0, std::fmax(std::fabs(a), std::fabs(b)));
    const double rel   = std::fabs(a - b) / denom;
    return (rel > max_rel) ? rel : max_rel;
}
} /* anonymous namespace */

/* ============================================================================
 * DiffFilterSpec host hooks.
 * ========================================================================== */

/* is_active — legacy dynamic_diff_vel_calc active filter (gas, live, massive). */
bool DiffFilterSpec::is_active(int i)
{
    if (P[i].Type != 0)   return false;
    if (P[i].TimeBin < 0) return false;
    if (P[i].Mass <= 0)   return false;
    return true;
}

/* search_radius — wide dynamic-diffusion kernel: fac * P[i].KernelRadius. */
double DiffFilterSpec::search_radius(const neighbor_loop_args& args,
                                     int /*active_slot*/, int i)
{
    return nlr_host_all_ptr()->TurbDynamicDiffFac * (double)args.P[i].KernelRadius;
}

DiffFilterSpec::CallScalars
DiffFilterSpec::populate_call_scalars(const neighbor_loop_args& /*args*/)
{
    CallScalars s;
    s.common            = nlr_common_scalars_from_all();
    s.turb_dyn_diff_fac = nlr_host_all_ptr()->TurbDynamicDiffFac;
    return s;
}

/* apply_active_writeback — scatter directly into CellP[i]. The preloop
 * (dynamic_diff_vel_calc_initial_operations_preloop) pre-zeroes Norm_hat /
 * FilterWidth_bar / MaxDistance_for_grad and seeds Velocity_bar = VelPred /
 * TurbDynamicDiffSmoothing, so += / max-merge here matches the legacy
 * scatter in dynamic_diffusion_velocities.cc:82-88. */
void DiffFilterSpec::apply_active_writeback(const neighbor_loop_args& args,
                                            int /*active_slot*/, int i,
                                            const AccumData& accum)
{
    struct gas_cell_data& C = args.CellP[i];
    C.Norm_hat += accum.norm_hat;
    for (int k = 0; k < 3; k++) C.Velocity_bar[k] += accum.velocity_bar_delta[k];
    if (accum.filter_width_bar  > C.FilterWidth_bar)
        C.FilterWidth_bar = (MyFloat)accum.filter_width_bar;
    if (accum.max_dist_for_grad > C.MaxDistance_for_grad)
        C.MaxDistance_for_grad = (MyFloat)accum.max_dist_for_grad;
}

void DiffFilterSpec::merge_accum(AccumData& local, const AccumData& peer)
{
    local.norm_hat += peer.norm_hat;
    for (int k = 0; k < 3; k++) local.velocity_bar_delta[k] += peer.velocity_bar_delta[k];
    if (peer.filter_width_bar  > local.filter_width_bar)  local.filter_width_bar  = peer.filter_width_bar;
    if (peer.max_dist_for_grad > local.max_dist_for_grad) local.max_dist_for_grad = peer.max_dist_for_grad;
}

double DiffFilterSpec::compare_accum(const AccumData& local, const AccumData& oracle)
{
    double m = 0.0;
    m = nlr_rel_update(m, local.norm_hat,          oracle.norm_hat);
    for (int k = 0; k < 3; k++)
        m = nlr_rel_update(m, local.velocity_bar_delta[k], oracle.velocity_bar_delta[k]);
    m = nlr_rel_update(m, local.filter_width_bar,  oracle.filter_width_bar);
    m = nlr_rel_update(m, local.max_dist_for_grad, oracle.max_dist_for_grad);
    return m;
}

/* No j-side writes → oracle brute pass needs no suppression. */
void DiffFilterSpec::set_oracle_brute_pass(DeviceContext& /*ctx*/, bool /*on*/) {}

double DiffFilterSpec::symmetric_neighbor_radius_scale()
{
    return nlr_host_all_ptr()->TurbDynamicDiffFac;
}

/* ============================================================================
 * DynDiffSpec host hooks.
 * ========================================================================== */

/* is_active — legacy dynamic_diff_calc filter (gas + massive; NO TimeBin
 * check — distinct from DiffFilterSpec; see dynamic_diffusion.cc:133). */
bool DynDiffSpec::is_active(int i)
{
    if (P[i].Type != 0) return false;
    if (P[i].Mass <= 0) return false;
    return true;
}

double DynDiffSpec::search_radius(const neighbor_loop_args& args,
                                  int /*active_slot*/, int i)
{
    return nlr_host_all_ptr()->TurbDynamicDiffFac * (double)args.P[i].KernelRadius;
}

DynDiffSpec::CallScalars
DynDiffSpec::populate_call_scalars(const neighbor_loop_args& args)
{
    CallScalars s;
    s.common            = nlr_common_scalars_from_all();
    s.turb_dyn_diff_fac = nlr_host_all_ptr()->TurbDynamicDiffFac;
    s.dynamic_iteration = static_cast<const Aux*>(args.aux)->dynamic_iteration;
    return s;
}

/* apply_active_writeback — scatter into the caller's DynamicDiffDataPasser.
 * dynamic_fac[/_const] every iteration; the iter-0-only block gated on
 * aux->dynamic_iteration. Matches dynamic_diffusion.cc:184-210. */
void DynDiffSpec::apply_active_writeback(const neighbor_loop_args& args,
                                         int /*active_slot*/, int i,
                                         const AccumData& accum)
{
    Aux* aux = static_cast<Aux*>(args.aux);
    struct temporary_data_dyndiff& d = aux->dddp[i];

    for (int k = 0; k < 3; k++) {
        for (int v = 0; v < 3; v++) {
            d.dynamic_fac[k][v]       += (MyDouble)accum.dynamic_fac[k][v];
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
            d.dynamic_fac_const[k][v] += (MyDouble)accum.dynamic_fac_const[k][v];
#endif
        }
    }

    if (aux->dynamic_iteration == 0) {
        if (accum.filter_width_hat > (double)d.FilterWidth_hat)
            d.FilterWidth_hat = (MyFloat)accum.filter_width_hat;
        d.Dynamic_numerator_hat   += (MyDouble)accum.dynamic_numerator_hat;
        d.Dynamic_denominator_hat += (MyDouble)accum.dynamic_denominator_hat;
        for (int k = 0; k < 3; k++) {
            if (accum.maxima_velocity_hat[k] > d.Maxima.Velocity_hat[k])
                d.Maxima.Velocity_hat[k] = accum.maxima_velocity_hat[k];
            if (accum.minima_velocity_hat[k] < d.Minima.Velocity_hat[k])
                d.Minima.Velocity_hat[k] = accum.minima_velocity_hat[k];
            for (int v = 0; v < 3; v++) {
                d.ProductVelocity_hat[k][v] += (MyDouble)accum.product_velocity_hat[k][v];
                d.GradVelocity_hat[k][v]    += (MyDouble)accum.grad_velocity_hat[k][v];
            }
        }
    }
}

void DynDiffSpec::merge_accum(AccumData& local, const AccumData& peer)
{
    for (int k = 0; k < 3; k++) {
        for (int v = 0; v < 3; v++) {
            local.dynamic_fac[k][v]          += peer.dynamic_fac[k][v];
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
            local.dynamic_fac_const[k][v]    += peer.dynamic_fac_const[k][v];
#endif
            local.grad_velocity_hat[k][v]    += peer.grad_velocity_hat[k][v];
            local.product_velocity_hat[k][v] += peer.product_velocity_hat[k][v];
        }
        if (peer.maxima_velocity_hat[k] > local.maxima_velocity_hat[k])
            local.maxima_velocity_hat[k] = peer.maxima_velocity_hat[k];
        if (peer.minima_velocity_hat[k] < local.minima_velocity_hat[k])
            local.minima_velocity_hat[k] = peer.minima_velocity_hat[k];
    }
    if (peer.filter_width_hat > local.filter_width_hat)
        local.filter_width_hat = peer.filter_width_hat;
    local.dynamic_numerator_hat   += peer.dynamic_numerator_hat;
    local.dynamic_denominator_hat += peer.dynamic_denominator_hat;
}

double DynDiffSpec::compare_accum(const AccumData& local, const AccumData& oracle)
{
    double m = 0.0;
    for (int k = 0; k < 3; k++) {
        for (int v = 0; v < 3; v++) {
            m = nlr_rel_update(m, local.dynamic_fac[k][v],          oracle.dynamic_fac[k][v]);
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
            m = nlr_rel_update(m, local.dynamic_fac_const[k][v],    oracle.dynamic_fac_const[k][v]);
#endif
            m = nlr_rel_update(m, local.grad_velocity_hat[k][v],    oracle.grad_velocity_hat[k][v]);
            m = nlr_rel_update(m, local.product_velocity_hat[k][v], oracle.product_velocity_hat[k][v]);
        }
        m = nlr_rel_update(m, local.maxima_velocity_hat[k], oracle.maxima_velocity_hat[k]);
        m = nlr_rel_update(m, local.minima_velocity_hat[k], oracle.minima_velocity_hat[k]);
    }
    m = nlr_rel_update(m, local.filter_width_hat,        oracle.filter_width_hat);
    m = nlr_rel_update(m, local.dynamic_numerator_hat,   oracle.dynamic_numerator_hat);
    m = nlr_rel_update(m, local.dynamic_denominator_hat, oracle.dynamic_denominator_hat);
    return m;
}

void DynDiffSpec::set_oracle_brute_pass(DeviceContext& /*ctx*/, bool /*on*/) {}

double DynDiffSpec::symmetric_neighbor_radius_scale()
{
    return nlr_host_all_ptr()->TurbDynamicDiffFac;
}

/* ============================================================================
 * Toplevels.
 * ========================================================================== */

void difffilter_vel_calc_gpu_toplevel(void)
{
    int *active_list = nullptr;
    int  num_active = 0, num_global_active = 0;
    if (!nlr_build_active_list(DiffFilterSpec::is_active,
                               &active_list, &num_active, &num_global_active,
                               "difffilter_active_list")) {
        return;   /* no active gas anywhere this step */
    }

    DiffFilterSpec::Aux aux;   /* unused — DiffFilter scatters straight to CellP */
    neighbor_loop_args args = nlr_default_args();
    args.active_list = active_list;
    args.num_active  = num_active;
    args.aux         = &aux;
    run_neighbor_loop<DiffFilterSpec>(args);

    nlr_free_active_list(active_list);
}

void dynamicdiff_gpu_toplevel(int dynamic_iteration,
                              struct temporary_data_dyndiff *dddp)
{
    int *active_list = nullptr;
    int  num_active = 0, num_global_active = 0;
    if (!nlr_build_active_list(DynDiffSpec::is_active,
                               &active_list, &num_active, &num_global_active,
                               "dyndiff_active_list")) {
        return;
    }

    DynDiffSpec::Aux aux;
    aux.dddp              = dddp;
    aux.dynamic_iteration = dynamic_iteration;

    neighbor_loop_args args = nlr_default_args();
    args.active_list = active_list;
    args.num_active  = num_active;
    args.aux         = &aux;
    run_neighbor_loop<DynDiffSpec>(args);

    nlr_free_active_list(active_list);
}

#endif /* TURB_DIFF_DYNAMIC */
