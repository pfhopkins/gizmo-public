/* hydro/cellcorrections_loop.h — first-pass cell volume corrections.
 *
 * Defines CellcorrectionsSpec for the neighbor-loop runner. Computes the
 * second-order volume correction Volume_1[i] = sum_j Volume_0[j]^2 wk(r, h_j)
 * for active gas particles. Pure i-side accumulate: no j-side writes, no
 * ghost-writeback, no iteration. Same topology as gradients / hydro_force
 * (symmetric gas-gas, GasOnly cache, per-active P[i].KernelRadius) — this
 * Spec is the first consumer in the hydro corridor, ahead of GradientsSpec.
 * Matches the legacy host walker in hydro/density.cc:62-162; consumes the
 * corridor's shared external CSR when published, else the runner builds
 * its own CSR from a narrow active list.
 *
 * Written by Philip F. Hopkins (phopkins@caltech.edu) for GIZMO. */

#ifndef CELLCORRECTIONS_LOOP_H
#define CELLCORRECTIONS_LOOP_H

/* Kokkos_Core.hpp must precede allvars.h (its macros may conflict with stdlib
 * names). */
#include <Kokkos_Core.hpp>
#include "../declarations/allvars.h"
#include "../declarations/multifluid_helpers.h"

#ifdef HYDRO_VOLUME_CORRECTIONS

#include "../mesh/neighbor_loop_runner.h"
#include "../mesh/mode_b_local_walker.h"      /* MODE_B_SEARCH_*, MODE_B_RADIUS_* */
#include "gradient_functions.h"               /* GasGrad_isactive_gpu — KOKKOS_INLINE
                                                * (host+device callable narrow
                                                * predicate; used by load_active's
                                                * `enabled` flag) */
/* NOTE: caller translation units must include "../mesh/kernel.h" BEFORE
 * this header. kernel.h has no include guards (defines static inline
 * kernel_main / kernel_hinv used by the pair body below); double-include
 * triggers redefinition errors. The runner and cellcorrections_loop.cc
 * both include kernel.h first per this convention (sink_env1_loop.h
 * follows the same pattern). gradient_functions.h is also include-guard-
 * free for kernel.h reasons — same ordering applies. */

/* No KOKKOS_INLINE_FUNCTION fallback here — Kokkos_Core.hpp is included
 * unconditionally above. This Spec carries device-callable pair-kernel
 * accessors; misordered includes must compile-fail loudly, not silently
 * resolve to host-only `inline`. Same convention as neighbor_loop_runner.h. */

/* Forward-decl from hydro/gradients.cc — the GasGrad_isactive predicate
 * gates gradients-eligible cells; cellcorrections must match it
 * (legacy semantic: the cell must pass the same active-gas filter that
 * gradients applies, since Volume_1 is consumed by gradients). */
int GasGrad_isactive(int i, struct particle_data *pp, struct gas_cell_data *cell);

/* ============================================================================
 * Per-pair physics types (file scope; PascalCase).
 * ========================================================================== */

struct CellcorrectionsActiveData {
    Vec3<double> pos;       /* P[i].Pos — only needed to compute dp = i - j */
    double       h_search;  /* per-active radius (P[i].KernelRadius) — required
                             * by runner's Mode B remote walker contract
                             * (mesh/neighbor_loop_runner.cc:927) */
    bool         enabled;   /* external-CSR row gate: the corridor row list
                             * is broad (Type==0 && Mass>0) but the narrow
                             * filter (KernelRadius>0, Density>0,
                             * GasGrad_isactive, …) is finer. Set in
                             * load_active; pair_kernel early-returns if
                             * false so subset-non-members contribute zero.
                             * On the non-corridor / Mode B path the active
                             * list is already narrow, so enabled=true. */
#ifdef HYDRO_MULTIFLUID
    unsigned char FluidType;  /* packed P[i].FluidType — for same_lagrangian_fluid_id() */
#endif
};

struct CellcorrectionsAccum {
    double volume_1;    /* sum_j Volume_0[j]^2 wk(r, h_j) */
};

struct CellcorrectionsNeighborData {
    Vec3<double> pos;           /* P[j].Pos */
    double       kernel_radius; /* P[j].KernelRadius — pair kernel uses j's radius */
    double       volume_0;      /* CellP[j].Volume_0 — the contribution weight */
#ifdef HYDRO_MULTIFLUID
    unsigned char FluidType;    /* packed P[j].FluidType */
#endif
};

/* ============================================================================
 * Inline pair body — single source of truth for cellcorrections physics.
 * Pure i-side accumulate using j's kernel radius for the contribution
 * filter and weight (matches legacy hydro/density.cc:130-149). The early
 * return on r2 >= h_j*h_j is the legacy filter; we still get called for
 * such pairs because MODE_B_SEARCH_SYMMETRIC reaches max(h_i, h_j).
 * ========================================================================== */
KOKKOS_INLINE_FUNCTION
static void cellcorrections_pair_kernel(const CellcorrectionsActiveData &active,
                                        const CellcorrectionsNeighborData &nb,
                                        CellcorrectionsAccum &out)
{
#if defined(HYDRO_MULTIFLUID)
    if (!same_lagrangian_fluid_id(active.FluidType, nb.FluidType)) return;
#endif
    /* Per-row narrow-filter gate (external-CSR consumption): the corridor
     * row list is broad; rows that fail the narrow GasGrad_isactive
     * predicate contribute zero. apply_active_writeback's += keeps
     * CellP[i].Volume_1 untouched for disabled rows (accum stays 0). */
    if(!active.enabled) return;

    Vec3<double> dp;
    dp[0] = active.pos[0] - nb.pos[0];
    dp[1] = active.pos[1] - nb.pos[1];
    dp[2] = active.pos[2] - nb.pos[2];
    nearest_xyz(dp);
    double r2 = dp.norm_sq();
    double h_j = nb.kernel_radius;
    if(r2 >= h_j * h_j) return;
    double hinv, hinv3, hinv4, wk = 0, dwk = 0;
    kernel_hinv(h_j, &hinv, &hinv3, &hinv4);
    double u = sqrt(r2) * hinv;
    kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
    out.volume_1 += nb.volume_0 * nb.volume_0 * wk;
}

/* ============================================================================
 * CellcorrectionsSpec — the NeighborLoopSpec contract.
 * ========================================================================== */

struct CellcorrectionsSpec {
    /* Identity */
    static constexpr const char *loop_name = "cellcorrections";
    static constexpr ModeBEvalOMP modeb_eval_omp = ModeBEvalOMP::BitwiseReadonly; /* i-side AccumData only (out.volume_1 += vol0^2*wk); no j-write, no atomics -> threaded eval bit-identical */

    /* Search policy. Symmetric so the pair contributes whenever r < h_j
     * (j's kernel) even if r > h_i; gas-only neighbor_type_mask matches
     * legacy gpu_ngb_list_build's gas-only j-type filter. */
    static constexpr int                     search_mode        = MODE_B_SEARCH_SYMMETRIC;
    static constexpr unsigned int            neighbor_type_mask = (1u << 0);
    static constexpr mode_b_radius_policy_t  radius_policy      = MODE_B_RADIUS_DEFAULT;

    /* Writeback policy. Pure i-side accumulate; no ghost-writeback. */
    static constexpr WritePattern   write_pattern    = WritePattern::ActiveReduceOnly;
    static constexpr SidxCacheKind  sidx_cache_kind  = SidxCacheKind::GasOnly;
    static constexpr bool           uses_ghost_writeback       = false;
    static constexpr bool           uses_ghost_write_detector  = false;

    /* Active-particle predicate. Matches the legacy filter
     * hydro/density.cc:81-83 (gas + Mass>0 + GasGrad_isactive +
     * KernelRadius>0). The toplevel passes this to nlr_build_active_list. */
    static bool is_active(int i) {
        if(P[i].Type != 0)              return false;
        if(P[i].Mass <= 0)              return false;
        if(P[i].KernelRadius <= 0)      return false;
        if(!GasGrad_isactive(i, P, CellP)) return false;
        return true;
    }

    /* Per-pair physics types */
    using CallScalars   = NlrCommonScalars;   /* nothing loop-specific needed */
    using ActiveData    = CellcorrectionsActiveData;
    using AccumData     = CellcorrectionsAccum;
    using NeighborData  = CellcorrectionsNeighborData;

    using ScatterData    = NoScatter;
    using IdentityFields = NoIdentity;
    using IterControl    = NotIterative;
    using DeviceContext  = NeighborLoopDeviceContextBase;

    /* Per-active host hook: pre-arena search radius from P[i].KernelRadius. */
    static double search_radius(const neighbor_loop_args& /*args*/,
                                 int /*active_slot*/, int i)
    {
        return (double)P[i].KernelRadius;
    }

    /* Per-call scalars: only the common cosmology block. The pair body
     * doesn't actually read these (geometry only), but the runner contract
     * still requires populate_call_scalars to be present. */
    static CallScalars populate_call_scalars(const neighbor_loop_args& /*args*/)
    {
        return nlr_common_scalars_from_all();
    }

    /* Device hooks (KOKKOS_INLINE_FUNCTION — host+device callable). */

    KOKKOS_INLINE_FUNCTION
    static void zero_accum(AccumData& accum)
    {
        accum.volume_1 = 0.0;
    }

    KOKKOS_INLINE_FUNCTION
    static ActiveData load_active(const NeighborLoopDeviceContextBase& ctx,
                                   int /*active_slot*/, int i,
                                   double h_search,
                                   const CallScalars& /*scalars*/)
    {
        ActiveData active;
        active.pos      = ctx.P[i].Pos;
        active.h_search = h_search;
        /* Narrow-active filter via the existing device-callable
         * GasGrad_isactive_gpu() — same predicate gradients uses, ensures
         * cellcorrections's "enabled" rows match what gradients will see.
         * For the 5a / Mode B path the runner-built active list is already
         * narrow so this is tautologically true; for the corridor-broad-
         * row-list / external-CSR path this flags subset-non-members so
         * the pair_kernel early-returns and accum stays zero. */
        active.enabled = ctx.CellP
                          ? (GasGrad_isactive_gpu(i, ctx.P, ctx.CellP) != 0)
                          : false;
#ifdef HYDRO_MULTIFLUID
        active.FluidType = ctx.P[i].FluidType;
#endif
        return active;
    }

    KOKKOS_INLINE_FUNCTION
    static NeighborData load_neighbor(const NeighborLoopDeviceContextBase& ctx,
                                       int j,
                                       const IdentitySidecar& /*id*/,
                                       const ActiveData& /*active*/)
    {
        NeighborData nb;
        nb.pos           = ctx.P[j].Pos;
        nb.kernel_radius = (double)ctx.P[j].KernelRadius;
        nb.volume_0      = (ctx.CellP && ctx.P[j].Type == 0)
                              ? (double)ctx.CellP[j].Volume_0
                              : 0.0;
#ifdef HYDRO_MULTIFLUID
        nb.FluidType     = ctx.P[j].FluidType;
#endif
        return nb;
    }

    /* The physics — forwards to the inline pair body above. */
    KOKKOS_INLINE_FUNCTION
    static void pair_kernel(const ActiveData& active,
                             const NeighborData& neighbor,
                             AccumData& accum,
                             NoScatter& /*scatter*/,
                            const CallScalars& /*cs*/)
    {
        cellcorrections_pair_kernel(active, neighbor, accum);
    }

    /* Host writebacks (post-dispatch) — body in cellcorrections_loop.cc */
    static void apply_active_writeback(const neighbor_loop_args& args,
                                        int active_slot, int i,
                                        const AccumData& accum);

    /* Per-field merge of a peer-rank accum into local accum (Mode B remote),
     * and the join of the Mode-A within-row lane reduction. */
    KOKKOS_INLINE_FUNCTION
    static void merge_accum(AccumData& local, const AccumData& peer)
    {
        local.volume_1 += peer.volume_1;
    }

};

#endif /* HYDRO_VOLUME_CORRECTIONS */
#endif /* CELLCORRECTIONS_LOOP_H */
